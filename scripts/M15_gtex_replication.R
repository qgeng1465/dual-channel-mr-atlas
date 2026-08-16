#!/usr/bin/env Rscript
# =============================================================================
# M15_gtex_replication.R — P1：106 个 strong coloc 命中 × GTEx 独立 eQTL 源复现
# =============================================================================
# 目的：用 GTEx v8（独立于 eQTLGen 的 eQTL 源，n≈838 中位组织）复核 106 个
#   strong coloc 命中的"基因-结局方向"是否可复现（Wald 单工具方向核查）。
# 与 M6 的区别（P1 复现性设计）：
#   - 基因集：106 个 strong coloc 命中（非 M6 的 21 候选）
#   - 阈值：GTEx egenes lead（FDR 显著），非 M6 的 p<5e-6 预注册工具阈值
#     → 这是"复现性方向核查"（探索性），不构成预注册工具变更
#   - 报告：覆盖数 / 结局匹配数 / 方向一致率 / 同变异子集一致率
# 数据：GTEx v8 egenes（每组织每基因 top cis-eQTL，含 rs_id_dbSNP151_GRCh38p7）
#   + OpenGWAS 结局（ebi-a-GCST006867/005194/005186），TwoSampleMR 同 M6
# 纪律：单工具敏感性；与 eQTLGen 主分析锚定不同，不移动门柱。
# 用法：PATH=$R_ENV/bin:$PATH Rscript scripts/M15_gtex_replication.R
# 输出：results/grid/gtex_replication_p1.csv + 终端摘要
# =============================================================================
suppressMessages({
  library(data.table)
  suppressPackageStartupMessages(library(TwoSampleMR))
  suppressPackageStartupMessages(library(ieugwasr))
  library(jsonlite)
})
proj <- "/data/qiushuogeng/projects/dual-channel-mr-atlas"
res  <- file.path(proj, "results")
GTEX <- file.path(proj, "data/gtex")
log <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n")
prereg <- file.path(proj, "docs/PREREGISTRATION.md")
lock   <- file.path(proj, "docs/PREREGISTRATION.md.sha256")
stopifnot(file.exists(prereg), file.exists(lock), tools::md5sum(prereg) == readLines(lock))

OUTCOMES <- list(t2d = "ebi-a-GCST006867", cad = "ebi-a-GCST005194", fbg = "ebi-a-GCST005186")
TISSUES <- c("Liver","Pancreas","Whole_Blood","Adipose_Subcutaneous","Muscle_Skeletal","Artery_Coronary")

hits <- fread(file.path(res, "grid/transcript_coloc_hits.csv"))
hits <- hits[tier == "strong"]
log("strong 命中: ", nrow(hits))

load_tissue <- function(t) {
  f <- file.path(GTEX, paste0(t, ".egenes.txt.gz"))
  if (!file.exists(f)) return(NULL)
  d <- fread(cmd = paste0("zcat ", shQuote(f)), sep = "\t", header = TRUE, nThread = 4)
  if (nrow(d) == 0) return(NULL)
  d[, Tissue := t]
  setnames(d, c("gene_name", "rs_id_dbSNP151_GRCh38p7"), c("symbol", "rsid"), skip_absent = TRUE)
  d[, rsid := fifelse(is.na(rsid) | rsid %in% c("", ".", "-"), variant_id, rsid)]
  d
}
tis <- list()
for (t in TISSUES) { tis[[t]] <- load_tissue(t); if (!is.null(tis[[t]])) log("组织 ", t, ": ", nrow(tis[[t]]), " 基因") }
tis <- tis[!vapply(tis, is.null, logical(1))]
egenes <- rbindlist(tis)
egenes <- egenes[abs(tss_distance) <= 1e6]
log("egenes 总行: ", nrow(egenes))

# 每个 hit 基因取跨组织最小 pval 的 lead 作为最佳组织工具（同基因不同组织 lead 可能不同）
best <- egenes[symbol %in% hits$symbol][order(pval_nominal), .SD[1], by = symbol]
log("106 命中在 GTEx 有 lead eQTL 的基因: ", nrow(best))

run_one <- function(h, rowg) {
  # rowg: GTEx best lead 行
  exp_dat <- data.frame(SNP = rowg$rsid,
    effect_allele.exposure = rowg$alt, other_allele.exposure = rowg$ref,
    eaf.exposure = ifelse(rowg$ref_factor == 1, 1 - rowg$maf, rowg$maf),
    beta.exposure = rowg$slope, se.exposure = rowg$slope_se,
    pval.exposure = rowg$pval_nominal,
    id.exposure = h$symbol, exposure = h$symbol, stringsAsFactors = FALSE)
  out <- tryCatch(extract_outcome_data(snps = rowg$rsid, outcomes = OUTCOMES[[h$outcome]], proxies = TRUE),
                  error = function(e) NULL)
  if (is.null(out) || nrow(out) == 0)
    return(data.table(gene = h$gene, symbol = h$symbol, outcome = h$outcome, pp_h4 = h$PP.H4,
                      gtex_tissue = rowg$Tissue, gtex_lead = rowg$variant_id, gtex_rsid = rowg$rsid,
                      gtex_lead_p = rowg$pval_nominal, gtex_b = NA, gtex_p = NA, concordant = NA, note = "无结局匹配"))
  dat <- harmonise_data(exp_dat, out, action = 2); dat <- dat[dat$mr_keep, ]
  if (nrow(dat) == 0)
    return(data.table(gene = h$gene, symbol = h$symbol, outcome = h$outcome, pp_h4 = h$PP.H4,
                      gtex_tissue = rowg$Tissue, gtex_lead = rowg$variant_id, gtex_rsid = rowg$rsid,
                      gtex_lead_p = rowg$pval_nominal, gtex_b = NA, gtex_p = NA, concordant = NA, note = "harmonise 无保留"))
  rw <- tryCatch(mr(dat, method_list = "mr_wald_ratio"), error = function(e) NULL)
  if (is.null(rw) || nrow(rw) == 0)
    return(data.table(gene = h$gene, symbol = h$symbol, outcome = h$outcome, pp_h4 = h$PP.H4,
                      gtex_tissue = rowg$Tissue, gtex_lead = rowg$variant_id, gtex_rsid = rowg$rsid,
                      gtex_lead_p = rowg$pval_nominal, gtex_b = NA, gtex_p = NA, concordant = NA, note = "MR 无输出"))
  gb <- rw$b; gp <- rw$pval
  same_var <- !is.na(h$top_snp) & h$top_snp == rowg$rsid
  conc <- (gb > 0) == (h$stage2_b > 0)
  data.table(gene = h$gene, symbol = h$symbol, outcome = h$outcome, pp_h4 = h$PP.H4,
             gtex_tissue = rowg$Tissue, gtex_lead = rowg$variant_id, gtex_rsid = rowg$rsid,
             gtex_lead_p = rowg$pval_nominal, gtex_b = gb, gtex_p = gp,
             concordant = conc, same_variant = same_var, note = paste0("eQTLGen top=", h$top_snp))
}

out <- list(); done <- 0L
for (i in seq_len(nrow(hits))) {
  h <- hits[i]
  rg <- best[symbol == h$symbol]
  if (nrow(rg) == 0) {
    out[[length(out) + 1]] <- data.table(gene = h$gene, symbol = h$symbol, outcome = h$outcome,
      pp_h4 = h$PP.H4, gtex_tissue = NA, gtex_lead = NA, gtex_rsid = NA, gtex_lead_p = NA,
      gtex_b = NA, gtex_p = NA, concordant = NA, same_variant = NA, note = "GTEx 无显著 cis-eQTL")
    next
  }
  out[[length(out) + 1]] <- run_one(h, rg[1])
  done <- done + 1L
  if (done %% 10 == 0) log("  完成 ", done, "/", nrow(hits))
}
res_df <- rbindlist(out, fill = TRUE)
fwrite(res_df, file.path(res, "grid/gtex_replication_p1.csv"))
log("已写 results/grid/gtex_replication_p1.csv | 行数 ", nrow(res_df))

# ---- 摘要 ----
n_cover  <- sum(!is.na(res_df$gtex_lead))
n_valid  <- sum(!is.na(res_df$concordant))
conc     <- res_df[concordant == TRUE]
cat("\n=== P1 GTEx 复现摘要 ===\n")
cat("106 strong 命中: ", nrow(res_df), "\n")
cat("GTEx 有显著 cis-eQTL（覆盖）: ", n_cover, "\n")
cat("结局匹配 + MR 有效: ", n_valid, "\n")
if (n_valid > 0) cat("方向一致率: ", nrow(conc), "/", n_valid, " = ", nrow(conc)/n_valid, "\n")
sv <- res_df[same_variant == TRUE & !is.na(same_variant)]
if (nrow(sv) > 0) cat("同变异子集: ", nrow(sv), "，方向一致 ", sum(sv$concordant, na.rm=TRUE), "/", nrow(sv), "\n")
cat("GTEx 名义显著(p<0.05)中方向一致: ", sum(res_df$concordant == TRUE & res_df$gtex_p < 0.05, na.rm=TRUE), "/",
    sum(res_df$gtex_p < 0.05, na.rm=TRUE), "\n")
