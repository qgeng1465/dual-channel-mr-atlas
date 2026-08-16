#!/usr/bin/env Rscript
# =============================================================================
# M3_transcript_grid_stage2.R — 转录扫描命中基因 本地 plink clump + IVW-MRE 复核
# =============================================================================
# 目的：对全量扫描层（M3_transcript_grid.R，lead Wald，hypothesis-generating）
#   的 BH-FDR q<0.05 命中基因做预注册锁定的第二阶段复核——本地 plink（1000G
#   EUR）LD clump r²<0.01@1000kb + IVW-MRE。回答"命中在独立化工具后是否存活"。
# 设计（预注册 v5 修订 + 深度层 M3_transcript_mr.R 一致纪律）：
#   - 工具 = 命中基因 cis±1Mb、p<5e-6、真实 eaf 的全部变异
#   - clump = 本地 plink --bfile data/ldref/1kg.v3/EUR（r²<0.01@1000kb，贪心 index）
#   - 结局提取 proxies=FALSE（与扫描层一致，防代理伪阳性）
#   - harmonise action=2（含回文处理）；主方法 mr_ivw_mre；nsnp≥2 补敏感性
#     ivw_fe / weighted_median / egger；nsnp=1 → Wald 退化（诚实标注）
#   - 多重检验：stage-2 结果集内按结局再算 BH-FDR（诚实双报 nominal + q）
#   - 全网格落盘含空/失败行；不选择性报告
# 用法：cd 项目 && PATH=<conda-root>/r-mr/bin:$PATH \
#       Rscript scripts/M3_transcript_grid_stage2.R
# 输入：results/grid/transcript_grid_hits.csv（982 hits）
#        data/eqtlgen/*.txt.gz（本地）| 结局 OpenGWAS API（proxies=FALSE）
#        data/ldref/1kg.v3/EUR.{bed,bim,fam}（results/ldref_ready）
# 输出：results/grid/transcript_grid_stage2.csv（全量+失败行）
#        + transcript_grid_stage2_hits.csv（存活）
#        + transcript_grid_stage2_sens.csv（nsnp≥2 敏感性）
#        + transcript_grid_stage2_funnel.tsv
# 缓存：results/grid/_stage2_clump_cache.rds + _outcome_stage2_*.rds（幂等续跑）
# =============================================================================
suppressMessages({
  library(data.table)
  suppressPackageStartupMessages(library(TwoSampleMR))
  suppressPackageStartupMessages(library(ieugwasr))
  library(jsonlite)
  library(parallel)
})
proj <- "<repo-root>"
res  <- file.path(proj, "results")
dir.create(file.path(res, "grid"), recursive = TRUE, showWarnings = FALSE)
log <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n")

# --- 预注册完整性校验 ----------------------------------------------------------
prereg <- file.path(proj, "docs/PREREGISTRATION.md")
lock   <- file.path(proj, "docs/PREREGISTRATION.md.sha256")
stopifnot(file.exists(prereg), file.exists(lock),
          tools::md5sum(prereg) == readLines(lock))
cfg <- read_json(file.path(res, "config.json"), auto_unbox = TRUE)
OUTCOMES <- list(t2d = "ebi-a-GCST006867", cad = "ebi-a-GCST005194", fbg = "ebi-a-GCST005186")
OUT_N <- c(t2d = "GCST006867", cad = "GCST005194", fbg = "GCST005186")
log("预注册哈希校验通过 ✔ | stage-2: 本地 plink clump r2<", cfg$instrument$clump_r2,
    "@", cfg$instrument$clump_kb, "kb EUR + ", cfg$mr_methods$primary)

PLINK <- file.path(proj, "tools/plink")
BFILE <- file.path(proj, "data/ldref/1kg.v3/EUR")
stopifnot(file.exists(PLINK), file.exists(paste0(BFILE, ".bim")),
          file.exists(file.path(res, "ldref_ready")))
log("plink=", PLINK, " | LD 参考=", BFILE, "（ldref_ready ✔）")

# --- 载入命中（只处理 q<0.05 网格行）-------------------------------------------
hits <- fread(file.path(res, "grid/transcript_grid_hits.csv"))
hits <- hits[ok == TRUE & !is.na(q) & q < 0.05]
tb <- table(hits$outcome)
log("hits 载入: ", nrow(hits), "（", paste(names(tb), tb, sep = "=", collapse = ", "),
    "）| 唯一基因 ", uniqueN(hits$gene))

# --- 读取 eQTLGen + AF（复用扫描层逻辑）----------------------------------------
log("读取 eQTLGen + AF ...")
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
ivs <- ivs[SNPChr == GeneChr & abs(SNPPos - GenePos) <= cfg$instrument$cis_window_kb * 1000 & !is.na(eaf)]
log("cis±1Mb & p<5e-6 & 真实 eaf 变异对: ", nrow(ivs))

genes <- unique(hits$gene)
log("命中基因数: ", length(genes))

# --- 本地 plink clump（按基因一次，跨结局复用；并行）----------------------------
cache_clump <- file.path(res, "grid/_stage2_clump_cache.rds")
if (file.exists(cache_clump)) {
  clumped <- readRDS(cache_clump)
  log("clump 缓存命中: ", length(clumped), " 基因")
} else {
  tmpdir <- file.path(res, "grid", "_stage2_clump_tmp")
  dir.create(tmpdir, showWarnings = FALSE)
  clump_one <- function(g) {
    d <- ivs[Gene == g]
    if (nrow(d) == 0) return(character(0))
    setorder(d, Pvalue)
    # plink 1.9 --clump 输入格式：两列带表头 SNP P（CHR/BP 由 --bfile 的 .bim 提供，
    #   与 ieugwasr::ld_clump_local 一致）；index SNP = .clumped 每行 SNP 列
    clumptxt <- file.path(tmpdir, paste0("c_", g, ".txt"))
    fwrite(d[, .(SNP, P = Pvalue)], clumptxt, sep = " ", col.names = TRUE)
    out <- file.path(tmpdir, paste0("c_", g))
    cmd <- paste(shQuote(PLINK), "--bfile", shQuote(BFILE),
                 "--clump", shQuote(clumptxt),
                 "--clump-p1 1",
                 "--clump-r2", cfg$instrument$clump_r2,
                 "--clump-kb", cfg$instrument$clump_kb,
                 "--out", shQuote(out), "--silent")
    st <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
    clf <- paste0(out, ".clumped")
    if (st != 0 || !file.exists(clf)) return(character(0))
    cl <- fread(clf, header = TRUE, fill = TRUE)
    if (nrow(cl) == 0) return(character(0))
    cl$SNP                        # 每行 = 一个 clump 的 index SNP = 独立工具
  }
  nc <- min(8, max(1, detectCores() - 2))
  log("本地 plink clump 并行（", nc, " 核，", length(genes), " 基因，参考 ",
      gsub(paste0(proj, "/"), "", BFILE), "）...")
  clres <- mclapply(genes, clump_one, mc.cores = nc)
  names(clres) <- genes
  clumped <- clres
  saveRDS(clumped, cache_clump)
  unlink(tmpdir, recursive = TRUE)
  nok <- sum(sapply(clumped, length) > 0)
  log("clump 完成并缓存: ", nok, "/", length(genes), " 基因有独立工具")
}
n_clump_ge1 <- sum(sapply(clumped, length) > 0)
log("clump 独立工具计数: ", n_clump_ge1, "/", length(genes), " 基因（≥1 工具）")

# --- 结局提取（clump 后独立工具的并集，proxies=FALSE，缓存）---------------------
cache_out <- function(on) {
  cf <- file.path(res, "grid", paste0("_outcome_stage2_", OUT_N[[on]], ".rds"))
  if (file.exists(cf)) { log("  读缓存 stage-2 结局 ", on); return(readRDS(cf)) }
  snps <- unique(unlist(clumped))
  log("  提取 stage-2 结局 ", on, " (", OUT_N[[on]], ", proxies=FALSE, ", length(snps), " SNP)...")
  o <- extract_outcome_data(snps = snps, outcomes = OUTCOMES[[on]], proxies = FALSE)
  if (!is.null(o) && nrow(o) > 0) saveRDS(o, cf)
  o
}
outs2 <- lapply(names(OUTCOMES), cache_out); names(outs2) <- names(OUTCOMES)
for (on in names(outs2))
  log("  ", on, ": 结局匹配 SNP 数 = ", if (is.null(outs2[[on]])) 0 else nrow(outs2[[on]]))

# --- 复核：每命中基因×结局 harmonise + IVW-MRE ---------------------------------
build_dat <- function(gene, on, snps) {
  ivs2 <- ivs[Gene == gene & SNP %in% snps]
  if (nrow(ivs2) == 0) return(NULL)
  oo <- outs2[[on]]
  if (is.null(oo)) return(NULL)
  oo <- oo[oo$SNP %in% snps, ]
  if (nrow(oo) == 0) return(NULL)
  exp_dat <- data.frame(SNP = ivs2$SNP, effect_allele.exposure = ivs2$AssessedAllele,
    other_allele.exposure = ivs2$OtherAllele, eaf.exposure = ivs2$eaf,
    beta.exposure = ivs2$Zscore / sqrt(ivs2$NrSamples), se.exposure = 1 / sqrt(ivs2$NrSamples),
    pval.exposure = ivs2$Pvalue, samplesize.exposure = ivs2$NrSamples,
    id.exposure = gene, exposure = ivs2$GeneSymbol[1], stringsAsFactors = FALSE)
  tryCatch({
    dat <- harmonise_data(exp_dat, oo, action = 2)
    dat[dat$mr_keep, ]
  }, error = function(e) NULL)
}

run_hit <- function(hit) {
  gene <- hit$gene; on <- hit$outcome
  n_cis <- nrow(ivs[Gene == gene])
  snps <- clumped[[gene]]
  fail <- function(note) data.frame(gene = gene, symbol = hit$symbol, outcome = on,
    scan_b = hit$b, scan_se = hit$se, scan_pval = hit$pval, scan_q = hit$q,
    n_cis_iv = n_cis, n_clump_index = length(snps), n_harmonised = 0,
    method = NA, b = NA, se = NA, pval = NA, ci_l = NA, ci_u = NA,
    direction_consistent = NA, q_stage2 = NA, survive_nominal = NA, survive_fdr = NA,
    ok = FALSE, note = note)
  if (length(snps) == 0)
    return(fail(paste0("clump 后 0 独立工具（cis 候选 ", n_cis, "）")))
  dat <- build_dat(gene, on, snps)
  if (is.null(dat) || nrow(dat) == 0)
    return(fail("结局无匹配 / harmonise 无保留（proxies=FALSE）"))
  n_har <- nrow(dat)
  r_main <- tryCatch(mr(dat, method_list = cfg$mr_methods$primary), error = function(e) NULL)
  method_used <- cfg$mr_methods$primary
  if (is.null(r_main) || nrow(r_main) == 0) {
    if (n_har == 1) {   # 单工具：IVW-MRE ≡ Wald，诚实标注退化
      r_main <- tryCatch(mr(dat, method_list = "mr_wald_ratio"), error = function(e) NULL)
      method_used <- "Wald ratio (单工具退化)"
    }
  }
  if (is.null(r_main) || nrow(r_main) == 0)
    return(fail("MR 无输出"))
  b <- r_main$b[1]; se <- r_main$se[1]; p <- r_main$pval[1]
  ci_l <- b - qnorm(0.975) * se; ci_u <- b + qnorm(0.975) * se
  dir <- if (!is.na(hit$b) && !is.na(b)) sign(hit$b) == sign(b) else NA
  data.frame(gene = gene, symbol = hit$symbol, outcome = on,
    scan_b = hit$b, scan_se = hit$se, scan_pval = hit$pval, scan_q = hit$q,
    n_cis_iv = n_cis, n_clump_index = length(snps), n_harmonised = n_har,
    method = method_used, b = b, se = se, pval = p, ci_l = ci_l, ci_u = ci_u,
    direction_consistent = dir, q_stage2 = NA, survive_nominal = NA, survive_fdr = NA,
    ok = TRUE, note = "")
}

log("逐命中基因×结局复核（", nrow(hits), " 对）...")
s2 <- rbindlist(lapply(seq_len(nrow(hits)), function(i) {
  rr <- run_hit(hits[i])
  if (i %% 100 == 0) log("  已处理 ", i, "/", nrow(hits))
  rr
}))
log("MR 完成: ", sum(s2$ok), "/", nrow(s2), " 对成功")

# --- 多重检验（stage-2 结果集内按结局 BH-FDR，诚实双报）-------------------------
for (on in names(OUTCOMES)) {
  sel <- s2$outcome == on & s2$ok & !is.na(s2$pval)
  if (sum(sel) > 0) s2$q_stage2[sel] <- p.adjust(s2$pval[sel], method = "BH")
}
s2[, survive_nominal := ok & !is.na(pval) & pval < 0.05 & direction_consistent == TRUE]
s2[, survive_fdr    := ok & !is.na(q_stage2) & q_stage2 < 0.05 & direction_consistent == TRUE]

# --- 敏感性（nsnp≥2：ivw_fe / weighted_median / egger）--------------------------
sens_rows <- list()
for (i in which(s2$ok & s2$n_harmonised >= 2)) {
  hit <- hits[i]; gene <- hit$gene; on <- hit$outcome
  dat <- build_dat(gene, on, clumped[[gene]])
  if (is.null(dat) || nrow(dat) < 2) next
  r_s <- tryCatch(mr(dat, method_list = cfg$mr_methods$sensitivity), error = function(e) NULL)
  if (!is.null(r_s) && nrow(r_s) > 0)
    sens_rows[[length(sens_rows) + 1]] <-
      data.frame(gene = gene, symbol = hit$symbol, outcome = on, nsnp = nrow(dat),
                 method = r_s$method, b = r_s$b, se = r_s$se, pval = r_s$pval)
}
sens_df <- if (length(sens_rows) > 0) rbindlist(sens_rows) else
  data.frame(gene = character(), symbol = character(), outcome = character(), nsnp = integer(),
             method = character(), b = numeric(), se = numeric(), pval = numeric())

# --- 落盘（全网格，含空/失败行；funnel）----------------------------------------
write.csv(s2, file.path(res, "grid/transcript_grid_stage2.csv"), row.names = FALSE)
sh <- s2[survive_nominal == TRUE]
write.csv(sh, file.path(res, "grid/transcript_grid_stage2_hits.csv"), row.names = FALSE)
write.csv(sens_df, file.path(res, "grid/transcript_grid_stage2_sens.csv"), row.names = FALSE)
funnel <- data.frame(stage = c("scan_fdr05_hits", "unique_genes", "genes_with_cis_iv",
                               "genes_clump_ge1", "pairs_outcome_matched", "pairs_mr_ok",
                               "survive_nominal", "survive_fdr"),
  count = c(nrow(hits), length(genes), sum(sapply(genes, function(g) nrow(ivs[Gene == g]) > 0)),
            n_clump_ge1, sum(s2$n_harmonised > 0), sum(s2$ok),
            sum(s2$survive_nominal), sum(s2$survive_fdr)))
write.table(funnel, file.path(res, "grid/transcript_grid_stage2_funnel.tsv"),
            sep = "\t", row.names = FALSE)
log("全网格落盘 ✔ results/grid/transcript_grid_stage2.csv（", nrow(s2), " 行）")
log("复核存活（nominal p<0.05 且方向一致）: ", sum(s2$survive_nominal),
    " | BH-FDR q<0.05 且方向一致: ", sum(s2$survive_fdr))

# --- 摘要：每结局存活 + top 名次 ------------------------------------------------
for (on in names(OUTCOMES)) {
  x <- s2[outcome == on]
  sn <- x[survive_nominal == TRUE]
  log("  [", on, "] 复核: ", nrow(x), " | 存活 nominal: ", nrow(sn),
      " | 存活 FDR: ", sum(x$survive_fdr))
  if (nrow(sn) > 0) for (i in seq_len(min(10, nrow(sn))))
    log("    ", sn$symbol[i], " scan_b=", round(sn$scan_b[i], 3),
        "→ stage2 b=", round(sn$b[i], 3), " p=", format(sn$pval[i], digits = 2),
        " q=", format(sn$q_stage2[i], digits = 2), " nsnp=", sn$n_harmonised[i])
}
log("=== stage-2 复核完成 ✔ ===")
