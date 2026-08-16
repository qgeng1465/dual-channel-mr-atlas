#!/usr/bin/env Rscript
# =============================================================================
# M3_transcript_drugtargets.R — 已知 T2D 药物靶点的转录通道 cis-MR（P1-2 四态分类）
# =============================================================================
# 目的：对已知 T2D 药物靶点基因，跑 eQTLGen 全血 cis-eQTL → T2D/CAD/FBG 的 MR，
#   得到"转录通道信号"侧；与 deCODE 蛋白通道（M3_protein_decode.R，含同基因）合并
#   构成四态分类（concordant / protein-only / transcript-only / both_null）。
# 方法（与 M3_transcript_mr.R v2 完全一致，预注册锁定）：
#   cis ±1Mb、p<5e-6、真实 eaf、EUR LD clump（r²<0.01@1000kb API）、主 mr_ivw_mre、
#   nsnp=1 Wald 补算；全网格落盘含空结果。
# 药物靶点基因（有全血 cis-eQTL 的）：GLP1R/SLC5A2/PPARG/KCNJ11/DPP4/INSR/TCF7L2/PRKAA1
#   （ABCC8/GCG/PDX1/SLC5A1/G6PC/PCK1 全血无显著 cis-eQTL → 转录通道不可用，如实记录）
# 用法：PATH=$R_ENV/bin:$PATH Rscript scripts/M3_transcript_drugtargets.R
# 输出：results/grid/transcript_drugtarget_mr.csv + results/grid/transcript_drugtarget_avail.csv
# =============================================================================
suppressMessages({
  library(data.table)
  suppressPackageStartupMessages(library(TwoSampleMR))
  suppressPackageStartupMessages(library(ieugwasr))
  library(jsonlite)
})
proj <- "/data/qiushuogeng/projects/dual-channel-mr-atlas"
res  <- file.path(proj, "results")
dir.create(file.path(res, "grid"), recursive = TRUE, showWarnings = FALSE)
log <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n")

# 预注册校验
prereg <- file.path(proj, "docs/PREREGISTRATION.md")
lock   <- file.path(proj, "docs/PREREGISTRATION.md.sha256")
stopifnot(file.exists(prereg), file.exists(lock),
          tools::md5sum(prereg) == readLines(lock))
cfg <- read_json(file.path(res, "config.json"), auto_unbox = TRUE)
log("预注册哈希校验通过 ✔ | 主方法=", cfg$mr_methods$primary)

OUTCOMES <- list(t2d = "ebi-a-GCST006867", cad = "ebi-a-GCST005194", fbg = "ebi-a-GCST005186")

# 药物靶点（第一优先级清单 + 常用补充；转录通道可用性由数据决定）
DRUGS <- c("GLP1R","SLC5A2","PPARG","KCNJ11","ABCC8","DPP4","INSR","TCF7L2","PRKAA1",
           "GCG","PDX1","SLC5A1","G6PC","PCK1")
DRUG_NOTES <- c(
  GLP1R = "GLP-1RA 靶点（受体激动剂）", SLC5A2 = "SGLT2i 靶点", PPARG = "TZD 靶点（核受体）",
  KCNJ11 = "磺脲类靶点 Kir6.2", ABCC8 = "磺脲类靶点 SUR1", DPP4 = "DPP4i 靶点",
  INSR = "胰岛素/增敏靶点", TCF7L2 = "T2D 最强易感基因", PRKAA1 = "AMPKα1（二甲双胍通路）",
  GCG = "胰高血糖素受体", PDX1 = "β细胞转录因子", SLC5A1 = "SGLT1",
  G6PC = "糖异生 G6Pase", PCK1 = "糖异生 PEPCK")

# --- 读 eQTLGen（暴露）+ AF（真实 eaf）----------------------------------------
log("读取 eQTLGen + AF...")
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
ivs <- merge(ivs, af[, .(SNP, AlleleA, AlleleB, AlleleB_all)], by = "SNP", all.x = TRUE)
ivs[, eaf := fcase(AssessedAllele == AlleleA, 1 - AlleleB_all,
                   AssessedAllele == AlleleB, AlleleB_all, default = NA_real_)]

# 可用性表（药物靶点 × 全血转录通道）：有显著 cis-eQTL / 有 p<5e-6 工具 / 真实 eaf
avail <- data.table(Gene = DRUGS)
avail[, n_eqtl_fdr := eqtl[GeneSymbol == DRUGS[.I], .N], by = .I]
avail[, n_iv_p5e6 := ivs[GeneSymbol == DRUGS[.I], .N], by = .I]
avail[, n_iv_eaf := ivs[GeneSymbol == DRUGS[.I] & !is.na(eaf), .N], by = .I]
avail[, transcript_available := n_iv_eaf > 0]
avail[, note := unname(DRUG_NOTES[Gene])]
write.csv(avail, file.path(res, "grid/transcript_drugtarget_avail.csv"), row.names = FALSE)
log("药物靶点转录可用性已落盘（transcript_drugtarget_avail.csv）")
print(avail)

# --- 对转录可用的靶点跑 MR（复用 M3 v2 逻辑）---------------------------------
build_ivs <- function(gene) {
  d <- ivs[GeneSymbol == gene & !is.na(eaf)]
  if (nrow(d) == 0) return(NULL)
  setorder(d, Pvalue)
  cl <- tryCatch(ld_clump(data.frame(rsid = d$SNP, pval = d$Pvalue, id = gene),
                          clump_kb = cfg$instrument$clump_kb,
                          clump_r2 = cfg$instrument$clump_r2,
                          clump_p = 1, pop = "EUR",
                          opengwas_jwt = ieugwasr::get_opengwas_jwt()),
                 error = function(e) { log("  [clump 失败] ", gene, ": ", conditionMessage(e)); NULL })
  if (is.null(cl) || nrow(cl) == 0) return(NULL)
  d[SNP %in% cl$rsid][order(Pvalue)]
}
run_gene_mr <- function(gene, outcome_id) {
  iv <- build_ivs(gene)
  if (is.null(iv) || nrow(iv) == 0)
    return(data.frame(gene = gene, outcome = outcome_id, nsnp = 0, b = NA, se = NA,
                      pval = NA, method = NA, ok = FALSE,
                      note = if (is.null(iv)) "无工具变量" else "clump 后无独立工具"))
  exp_dat <- data.frame(SNP = iv$SNP, effect_allele.exposure = iv$AssessedAllele,
    other_allele.exposure = iv$OtherAllele, eaf.exposure = iv$eaf,
    beta.exposure = iv$Zscore / sqrt(iv$NrSamples), se.exposure = 1 / sqrt(iv$NrSamples),
    pval.exposure = iv$Pvalue, samplesize.exposure = iv$NrSamples,
    id.exposure = gene, exposure = gene, stringsAsFactors = FALSE)
  tryCatch({
    out_dat <- extract_outcome_data(snps = iv$SNP, outcomes = outcome_id, proxies = TRUE)
    if (is.null(out_dat) || nrow(out_dat) == 0)
      return(data.frame(gene = gene, outcome = outcome_id, nsnp = 0, b = NA, se = NA,
                        pval = NA, method = NA, ok = FALSE, note = "无结局匹配"))
    dat <- harmonise_data(exp_dat, out_dat, action = 2); dat <- dat[dat$mr_keep, ]
    if (nrow(dat) == 0)
      return(data.frame(gene = gene, outcome = outcome_id, nsnp = nrow(iv), b = NA, se = NA,
                        pval = NA, method = NA, ok = FALSE, note = "harmonise 后无保留 SNP"))
    r_main <- mr(dat, method_list = cfg$mr_methods$primary)
    r_sens <- tryCatch(mr(dat, method_list = cfg$mr_methods$sensitivity), error = function(e) NULL)
    r_all <- if (!is.null(r_sens)) rbind(r_main, r_sens) else r_main
    if (nrow(r_all) == 0 && nrow(dat) == 1) {
      r_w <- tryCatch(mr(dat, method_list = "mr_wald_ratio"), error = function(e) NULL)
      if (!is.null(r_w) && nrow(r_w) > 0) r_all <- r_w
    }
    if (nrow(r_all) == 0)
      return(data.frame(gene = gene, outcome = outcome_id, nsnp = nrow(dat), b = NA, se = NA,
                        pval = NA, method = NA, ok = FALSE, note = "MR 无输出"))
    data.frame(gene = gene, outcome = outcome_id, nsnp = nrow(dat),
               b = r_all$b, se = r_all$se, pval = r_all$pval, method = r_all$method,
               ok = TRUE, note = paste0("clump ", nrow(exp_dat), "→harmonise ", nrow(dat)))
  }, error = function(e) data.frame(gene = gene, outcome = outcome_id, nsnp = 0, b = NA,
                                    se = NA, pval = NA, method = NA, ok = FALSE,
                                    note = conditionMessage(e)))
}

genes <- avail[transcript_available == TRUE, Gene]
log("转录可用的药物靶点（跑 MR）: ", paste(genes, collapse = ", "))
out <- list()
for (on in names(OUTCOMES)) {
  log("=== 结局: ", on, " ===")
  for (g in genes) {
    rr <- run_gene_mr(g, OUTCOMES[[on]])
    out[[length(out) + 1]] <- rr
    for (i in seq_len(nrow(rr)))
      if (rr$ok[i]) log("  ", g, "×", on, "[", rr$method[i], "] nsnp=", rr$nsnp[i],
                        " b=", round(rr$b[i], 3), " p=", format(rr$pval[i], digits = 2))
  }
}
res_df <- do.call(rbind, out)
write.csv(res_df, file.path(res, "grid/transcript_drugtarget_mr.csv"), row.names = FALSE)
log("药物靶点转录 MR 落盘 ✔ → results/grid/transcript_drugtarget_mr.csv")
log("（转录不可用靶点 ABCC8/GCG/PDX1/SLC5A1/G6PC/PCK1 在全血无 cis-eQTL，如实记录于 avail 表）")
