#!/usr/bin/env Rscript
# =============================================================================
# M3_wald_fallback.R — nsnp=1 单工具 Wald ratio 补算
# 背景：M3_transcript_mr.R 主方法 mr_ivw_mre 需 ≥2 工具；nsnp=1 时 IVW≡Wald ratio，
#      但 mr_ivw_mre 对 nsnp=1 无输出 → 用 mr_wald_ratio 补算（诚实标注）。
# 触发：v2 结果中 ok=FALSE & nsnp=1 的基因×结局对（HMGCR×t2d、HMGCR×fbg）。
# 用法：PATH=$R_ENV/bin:$PATH Rscript scripts/M3_wald_fallback.R
# 输出：results/grid/transcript_mr_v2_wald.csv（追加到主结果的口径，主方法标注 Wald）
# =============================================================================
suppressMessages({
  library(data.table)
  suppressPackageStartupMessages(library(TwoSampleMR))
  suppressPackageStartupMessages(library(ieugwasr))
  library(jsonlite)
})
proj <- "<repo-root>"
res  <- file.path(proj, "results")
log <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n")

# --- 与 M3 相同的读入与工具构建（只取目标基因）-------------------------------
eqtl <- fread(cmd = paste0("zcat ", file.path(proj, "data/eqtlgen/cis-eQTL-significant.txt.gz")),
              sep = "\t", header = TRUE, nThread = 4)
setnames(eqtl, c("Pvalue","SNP","SNPChr","SNPPos","AssessedAllele","OtherAllele",
                 "Zscore","Gene","GeneSymbol","GeneChr","GenePos",
                 "NrCohorts","NrSamples","FDR","BonferroniP"))
af <- fread(cmd = paste0("zcat ", file.path(proj, "data/eqtlgen/SNP_AF.txt.gz")),
            sep = "\t", header = TRUE, nThread = 4)
setnames(af, c("SNP","hg19_chr","hg19_pos","AlleleA","AlleleB",
               "allA_total","allAB_total","allB_total","AlleleB_all"))
ivs <- eqtl[Pvalue < 5e-6]
ivs <- merge(ivs, af[, .(SNP, AlleleA, AlleleB, AlleleB_all)], by = "SNP", all.x = TRUE)
ivs[, eaf := fcase(AssessedAllele == AlleleA, 1 - AlleleB_all,
                   AssessedAllele == AlleleB, AlleleB_all, default = NA_real_)]

build_ivs <- function(gene) {
  d <- ivs[GeneSymbol == gene & !is.na(eaf)]
  if (nrow(d) == 0) return(NULL)
  setorder(d, Pvalue)
  cl <- tryCatch(ld_clump(data.frame(rsid = d$SNP, pval = d$Pvalue, id = gene),
                          clump_kb = 1000, clump_r2 = 0.01, clump_p = 1,
                          pop = "EUR", opengwas_jwt = ieugwasr::get_opengwas_jwt()),
                 error = function(e) NULL)
  if (is.null(cl) || nrow(cl) == 0) return(NULL)
  d[SNP %in% cl$rsid][order(Pvalue)]
}

wald_for <- function(gene, outcome_id) {
  iv <- build_ivs(gene)
  if (is.null(iv) || nrow(iv) == 0)
    return(data.frame(gene = gene, outcome = outcome_id, nsnp = 0, b = NA, se = NA,
                      pval = NA, method = "Wald ratio", ok = FALSE, note = "无工具变量"))
  exp_dat <- data.frame(
    SNP = iv$SNP, effect_allele.exposure = iv$AssessedAllele,
    other_allele.exposure = iv$OtherAllele, eaf.exposure = iv$eaf,
    beta.exposure = iv$Zscore / sqrt(iv$NrSamples),
    se.exposure = 1 / sqrt(iv$NrSamples), pval.exposure = iv$Pvalue,
    samplesize.exposure = iv$NrSamples, id.exposure = gene, exposure = gene,
    stringsAsFactors = FALSE)
  out_dat <- extract_outcome_data(snps = iv$SNP, outcomes = outcome_id, proxies = TRUE)
  dat <- harmonise_data(exp_dat, out_dat, action = 2)
  dat <- dat[dat$mr_keep, ]
  if (nrow(dat) == 0)
    return(data.frame(gene = gene, outcome = outcome_id, nsnp = nrow(iv), b = NA, se = NA,
                      pval = NA, method = "Wald ratio", ok = FALSE, note = "harmonise 后无保留 SNP"))
  r <- mr(dat, method_list = "mr_wald_ratio")
  if (nrow(r) == 0)
    return(data.frame(gene = gene, outcome = outcome_id, nsnp = nrow(dat), b = NA, se = NA,
                      pval = NA, method = "Wald ratio", ok = FALSE, note = "Wald 无输出"))
  data.frame(gene = gene, outcome = outcome_id, nsnp = nrow(dat),
             b = r$b, se = r$se, pval = r$pval, method = r$method, ok = TRUE, note = "单工具 Wald 补算")
}

OUTCOMES <- c(t2d = "ebi-a-GCST006867", cad = "ebi-a-GCST005194", fbg = "ebi-a-GCST005186")
# 仅补 v2 中 ok=FALSE 且 nsnp=1 的对（现为 HMGCR×t2d、HMGCR×fbg）
targets <- list(c("HMGCR","t2d"), c("HMGCR","fbg"))
out <- lapply(targets, function(t) {
  g <- t[1]; oc <- OUTCOMES[t[2]]
  log("补算 ", g, "×", t[2], " (Wald ratio)")
  wald_for(g, oc)
})
res_df <- rbindlist(out)
print(res_df[, .(gene, outcome = names(OUTCOMES)[match(outcome, OUTCOMES)], nsnp, b, se, pval, note)],
      digits = 4)
write.csv(res_df, file.path(res, "grid/transcript_mr_v2_wald.csv"), row.names = FALSE)
log("Wald 补算落盘: results/grid/transcript_mr_v2_wald.csv")
