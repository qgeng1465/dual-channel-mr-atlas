#!/usr/bin/env Rscript
# =============================================================================
# M3_transcript_grid.R — 转录组全量 cis-MR 扫描层（PREREG v5 修订）
# =============================================================================
# 目的：eQTLGen 全部有显著 cis-eQTL 的基因（预计 ~16,000+）对 T2D/CAD/FBG 的
#   系统扫描，回答"全血转录本层面哪些基因有 MR 证据"（hypothesis-generating）。
# 设计（预注册 v5 修订锁定）：
#   - 工具 = 每基因最强 cis-eQTL（lead variant：cis ±1Mb、Pvalue<5e-6 内 Pvalue 最小）
#   - 方法 = 单工具 Wald ratio（nsnp=1）
#   - 结局提取 proxies=FALSE（扫描层防代理伪阳性；深度层保持 TRUE 不变）
#   - 多重检验：按结局 BH-FDR q<0.05
#   - 全网格落盘（含无工具/无结局匹配/失败行）
# 命中复核（第二阶段，LD 参考就绪后）：本地 plink（1000G EUR）clump r²<0.01@1000kb + IVW-MRE
# 用法：PATH=$R_ENV/bin:$PATH Rscript scripts/M3_transcript_grid.R
# 输出：results/grid/transcript_grid_mr.csv + results/grid/transcript_grid_hits.csv
#       + results/grid/transcript_grid_funnel.tsv
# =============================================================================
suppressMessages({
  library(data.table)
  suppressPackageStartupMessages(library(TwoSampleMR))
  suppressPackageStartupMessages(library(ieugwasr))
  library(jsonlite)
})
proj <- "<repo-root>"
res  <- file.path(proj, "results")
dir.create(file.path(res, "grid"), recursive = TRUE, showWarnings = FALSE)
log <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n")

prereg <- file.path(proj, "docs/PREREGISTRATION.md")
lock   <- file.path(proj, "docs/PREREGISTRATION.md.sha256")
stopifnot(file.exists(prereg), file.exists(lock), tools::md5sum(prereg) == readLines(lock))
cfg <- read_json(file.path(res, "config.json"), auto_unbox = TRUE)
OUTCOMES <- list(t2d = "ebi-a-GCST006867", cad = "ebi-a-GCST005194", fbg = "ebi-a-GCST005186")
OUT_N <- c(t2d = "GCST006867", cad = "GCST005194", fbg = "GCST005186")
log("预注册哈希校验通过 ✔（含 v5 修订） | 扫描层：每基因 lead cis-eQTL → Wald")

# --- 读取 eQTLGen（暴露）+ AF（真实 eaf）--------------------------------------
log("读取 eQTLGen + AF（~10.5M 显著对）...")
eqtl <- fread(cmd = paste0("zcat ", file.path(proj, "data/eqtlgen/cis-eQTL-significant.txt.gz")),
              sep = "\t", header = TRUE, nThread = 4)
setnames(eqtl, c("Pvalue","SNP","SNPChr","SNPPos","AssessedAllele","OtherAllele",
                 "Zscore","Gene","GeneSymbol","GeneChr","GenePos",
                 "NrCohorts","NrSamples","FDR","BonferroniP"))
af <- fread(cmd = paste0("zcat ", file.path(proj, "data/eqtlgen/SNP_AF.txt.gz")),
            sep = "\t", header = TRUE, nThread = 4)
setnames(af, c("SNP","hg19_chr","hg19_pos","AlleleA","AlleleB",
               "allA_total","allAB_total","allB_total","AlleleB_all"))
ivs <- eqtl[Pvalue < cfg$instrument$pval_thresh]
rm(eqtl)
ivs <- merge(ivs, af[, .(SNP, AlleleA, AlleleB, AlleleB_all)], by = "SNP", all.x = TRUE)
rm(af)
ivs[, eaf := fcase(AssessedAllele == AlleleA, 1 - AlleleB_all,
                   AssessedAllele == AlleleB, AlleleB_all, default = NA_real_)]
# cis ±1Mb（GenePos 为基因 TSS，hg19 同源）
ivs <- ivs[SNPChr == GeneChr & abs(SNPPos - GenePos) <= cfg$instrument$cis_window_kb * 1000 & !is.na(eaf)]
log("cis±1Mb & p<5e-6 & 有真实 eaf 的变异对: ", nrow(ivs))

# --- lead variant 选取：每基因 cis 内 Pvalue 最小的单变异 ----------------------
# 注：.SD[1] 已保留 GeneSymbol 列，不可再 merge（会生成 GeneSymbol.x/.y 导致 symbol 丢失）
setorder(ivs, Gene, Pvalue)
lead <- ivs[, .SD[1], by = Gene]
log("lead cis-eQTL 基因数: ", nrow(lead), "（含 symbol 映射 ", sum(!is.na(lead$GeneSymbol)), "）")
lead[, beta := Zscore / sqrt(NrSamples)]
lead[, se   := 1 / sqrt(NrSamples)]
lead_snps <- unique(lead$SNP)
log("唯一 lead SNP 数: ", length(lead_snps), "（结局按此去重提取，分块 API）")

# --- 结局提取（每结局一次，proxies=FALSE，缓存落盘）----------------------------
get_out <- function(on) {
  cache <- file.path(res, "grid", paste0("_outcome_", OUT_N[[on]], ".rds"))
  if (file.exists(cache)) { log("  读缓存结局 ", on); return(readRDS(cache)) }
  log("  提取结局 ", on, " (", OUT_N[[on]], ", proxies=FALSE, ", length(lead_snps), " SNP)...")
  o <- extract_outcome_data(snps = lead_snps, outcomes = OUTCOMES[[on]], proxies = FALSE)
  if (!is.null(o) && nrow(o) > 0) saveRDS(o, cache)
  o
}
outs <- lapply(names(OUTCOMES), get_out); names(outs) <- names(OUTCOMES)
for (on in names(outs)) log("  ", on, ": 结局匹配 SNP 数 = ", if (is.null(outs[[on]])) 0 else nrow(outs[[on]]))

# --- 逐基因 × 结局：harmonise + Wald ------------------------------------------
run_gene <- function(gene, on) {
  ld <- lead[Gene == gene]
  snp <- ld$SNP
  o <- outs[[on]]
  oo <- if (is.null(o)) NULL else o[o$SNP == snp, ]
  if (is.null(oo) || nrow(oo) == 0)
    return(data.frame(gene = gene, symbol = ld$GeneSymbol, outcome = on, nsnp = 0,
                      b = NA, se = NA, pval = NA, q = NA, method = NA, ok = FALSE,
                      note = "lead SNP 无结局匹配"))
  exp_dat <- data.frame(SNP = snp, effect_allele.exposure = ld$AssessedAllele,
    other_allele.exposure = ld$OtherAllele, eaf.exposure = ld$eaf,
    beta.exposure = ld$beta, se.exposure = ld$se, pval.exposure = ld$Pvalue,
    samplesize.exposure = ld$NrSamples, id.exposure = gene, exposure = ld$GeneSymbol,
    stringsAsFactors = FALSE)
  tryCatch({
    dat <- harmonise_data(exp_dat, oo, action = 2); dat <- dat[dat$mr_keep, ]
    if (nrow(dat) == 0)
      return(data.frame(gene = gene, symbol = ld$GeneSymbol, outcome = on, nsnp = 1,
                        b = NA, se = NA, pval = NA, q = NA, method = NA, ok = FALSE,
                        note = "harmonise 无保留（等位/回文不匹配）"))
    rw <- tryCatch(mr(dat, method_list = "mr_wald_ratio"), error = function(e) NULL)
    if (is.null(rw) || nrow(rw) == 0)
      return(data.frame(gene = gene, symbol = ld$GeneSymbol, outcome = on, nsnp = nrow(dat),
                        b = NA, se = NA, pval = NA, q = NA, method = NA, ok = FALSE, note = "Wald 无输出"))
    data.frame(gene = gene, symbol = ld$GeneSymbol, outcome = on, nsnp = nrow(dat),
               b = rw$b, se = rw$se, pval = rw$pval, q = NA, method = "Wald ratio",
               ok = TRUE, note = paste0("lead ", snp, " (cis p=", format(ld$Pvalue, digits = 2), ")"))
  }, error = function(e) data.frame(gene = gene, symbol = ld$GeneSymbol, outcome = on, nsnp = 1,
                                    b = NA, se = NA, pval = NA, q = NA, method = NA, ok = FALSE,
                                    note = conditionMessage(e)))
}

log("逐基因 × 结局 Wald 扫描（", nrow(lead), " 基因 × 3 结局）...")
out <- list(); k <- 0
for (on in names(OUTCOMES)) {
  for (g in lead$Gene) {
    k <- k + 1
    rr <- run_gene(g, on)
    out[[length(out) + 1]] <- rr
    if (k %% 2000 == 0) log("  已处理 ", k, "/", nrow(lead) * 3, " 基因×结局")
  }
}
res_df <- do.call(rbind, out)
res_df$q <- NA_real_
# 按结局做 BH-FDR
for (on in names(OUTCOMES)) {
  sel <- res_df$outcome == on & res_df$ok & !is.na(res_df$pval)
  p <- res_df$pval[sel]
  if (length(p) > 0) {
    q <- p.adjust(p, method = "BH")
    res_df$q[sel] <- q
  }
}
setDT(res_df)
write.csv(res_df, file.path(res, "grid/transcript_grid_mr.csv"), row.names = FALSE)
hits <- res_df[ok == TRUE & !is.na(q) & q < 0.05]
write.csv(hits, file.path(res, "grid/transcript_grid_hits.csv"), row.names = FALSE)
funnel <- data.frame(stage = c("genes_with_lead_cis", "genes_with_eaf", "gene_outcome_pairs",
                               "mr_ok", "fdr05_hits"),
  count = c(length(unique(lead$Gene)), sum(!is.na(lead$eaf)),
            nrow(res_df), sum(res_df$ok, na.rm = TRUE), nrow(hits)))
write.table(funnel, file.path(res, "grid/transcript_grid_funnel.tsv"), sep = "\t", row.names = FALSE)
log("全量扫描落盘 ✔")
for (on in names(OUTCOMES)) {
  h <- hits[outcome == on]
  log("  [", on, "] FDR q<0.05 命中: ", nrow(h))
  if (nrow(h) > 0) for (i in seq_len(min(10, nrow(h))))
    log("    ", h$symbol[i], " b=", round(h$b[i], 3), " p=", format(h$pval[i], digits = 2),
        " q=", format(h$q[i], digits = 2))
}
log("（命中基因第二阶段：本地 plink clump + IVW-MRE 复核，待 LD 参考就绪）")
