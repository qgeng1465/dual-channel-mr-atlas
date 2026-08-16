#!/usr/bin/env Rscript
# =============================================================================
# M3_protein_interval.R — 蛋白通道 cis-MR：INTERVAL SomaScan pQTL via OpenGWAS
# 数据源：OpenGWAS prot-a-*（INTERVAL, n=3,301, SomaScan），JWT 已配置。
# 对照可用性（prota_index.rds 核对）：HMGCR=prot-a-1354 ✔  ANGPTL3=prot-a-98 ✔
#   PCSK9/CETP/NPC1L1/APOC3 → SomaScan 未收录（如实记录）
# 方法：cis ±1Mb, p<5e-6, LD clump r²<0.01@1000kb EUR（API 模式）
#   主 mr_ivw_mre；敏感性 ivw_fe/WM/Egger；预注册哈希校验
# 用法：PATH=$R_ENV/bin:$PATH Rscript scripts/M3_protein_interval.R
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
dir.create(file.path(res, "funnel"), showWarnings = FALSE)
log <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n")

prereg <- file.path(proj, "docs/PREREGISTRATION.md")
lock   <- file.path(proj, "docs/PREREGISTRATION.md.sha256")
stopifnot(file.exists(prereg), file.exists(lock),
          tools::md5sum(prereg) == readLines(lock))
cfg <- read_json(file.path(res, "config.json"), auto_unbox = TRUE)
log("预注册哈希校验通过 ✔  | 主方法=", cfg$mr_methods$primary)

PROTEINS <- list(
  list(id = "prot-a-1354", gene = "HMGCR", note = "HMG-CoA reductase"),
  list(id = "prot-a-98",   gene = "ANGPTL3", note = "Angiopoietin-related protein 3")
)
log("INTERVAL 可用对照: HMGCR(prot-a-1354), ANGPTL3(prot-a-98)；PCSK9/CETP/NPC1L1/APOC3 未收录")

OUTCOMES <- list(t2d = "ebi-a-GCST006867", cad = "ebi-a-GCST005194", fbg = "ebi-a-GCST005186")
# hg19 坐标与 M5_protein_coloc.R 对齐（ENSEMBL GRCh37 REST 实查 2026-08-13）；
# 原 5:74652278 / 1:63051658 与 coloc 侧偏差 ~20kb/~12kb，可能让 cis ±1Mb 窗中心偏移
GENE_POS <- c(HMGCR = "5:74632154", ANGPTL3 = "1:63063158")  # hg19

run_protein_mr <- function(pid, gene, outcome_id) {
  pos <- GENE_POS[[gene]]; chr <- as.integer(sub(":.*", "", pos)); gp <- as.integer(sub(".*:", "", pos))
  inst <- tryCatch(
    extract_instruments(outcomes = pid, p1 = cfg$instrument$pval_thresh, clump = FALSE),
    error = function(e) NULL)
  if (is.null(inst) || nrow(inst) == 0)
    return(data.frame(gene = gene, outcome = outcome_id, nsnp = 0, b = NA, se = NA,
                      pval = NA, method = NA, ok = FALSE, note = "无 cis-pQTL 工具"))
  inst$chrom <- as.integer(gsub("chr", "", inst$chr.exposure))
  inst <- inst[inst$chrom == chr & abs(inst$pos.exposure - gp) <= 1e6, ]
  if (nrow(inst) == 0)
    return(data.frame(gene = gene, outcome = outcome_id, nsnp = 0, b = NA, se = NA,
                      pval = NA, method = NA, ok = FALSE, note = "cis ±1Mb 内无工具"))
  cl <- tryCatch(ld_clump(data.frame(rsid = inst$SNP, pval = inst$pval.exposure, id = gene),
                          clump_kb = cfg$instrument$clump_kb, clump_r2 = cfg$instrument$clump_r2,
                          clump_p = 1, pop = "EUR", opengwas_jwt = ieugwasr::get_opengwas_jwt()),
                 error = function(e) NULL)
  if (!is.null(cl) && nrow(cl) > 0) inst <- inst[inst$SNP %in% cl$rsid, ]
  if (nrow(inst) == 0)
    return(data.frame(gene = gene, outcome = outcome_id, nsnp = 0, b = NA, se = NA,
                      pval = NA, method = NA, ok = FALSE, note = "clump 后无独立工具"))
  tryCatch({
    out <- extract_outcome_data(snps = inst$SNP, outcomes = outcome_id, proxies = TRUE)
    dat <- harmonise_data(inst, out, action = 2)
    dat <- dat[dat$mr_keep, ]
    if (nrow(dat) == 0)
      return(data.frame(gene = gene, outcome = outcome_id, nsnp = 0, b = NA, se = NA,
                        pval = NA, method = NA, ok = FALSE, note = "harmonise 后无保留 SNP"))
    r_all <- mr(dat, method_list = c(cfg$mr_methods$primary, cfg$mr_methods$sensitivity))
    # nsnp=1：IVW 无输出但 IVW≡Wald ratio → 补 Wald（§6.4 口径，与 deCODE/转录通道一致）
    if (nrow(r_all) == 0 && nrow(dat) == 1) {
      r_w <- tryCatch(mr(dat, method_list = "mr_wald_ratio"), error = function(e) NULL)
      if (!is.null(r_w) && nrow(r_w) > 0) r_all <- r_w
    }
    if (nrow(r_all) == 0)
      return(data.frame(gene = gene, outcome = outcome_id, nsnp = nrow(dat), b = NA, se = NA,
                        pval = NA, method = NA, ok = FALSE, note = "MR 无输出"))
    data.frame(gene = gene, outcome = outcome_id, nsnp = nrow(dat),
               b = r_all$b, se = r_all$se, pval = r_all$pval,
               method = r_all$method, ok = TRUE, note = "")
  }, error = function(e) data.frame(gene = gene, outcome = outcome_id, nsnp = 0,
                                    b = NA, se = NA, pval = NA, method = NA,
                                    ok = FALSE, note = conditionMessage(e)))
}

out <- list()
for (p in PROTEINS) for (on in names(OUTCOMES)) {
  log("=== ", p$gene, " (", p$id, ") × ", on, " ===")
  rr <- run_protein_mr(p$id, p$gene, OUTCOMES[[on]])
  out[[length(out) + 1]] <- rr
  for (i in seq_len(nrow(rr)))
    if (rr$ok[i]) log("  ", p$gene, "×", on, "[", rr$method[i], "] nsnp=", rr$nsnp[i],
                      " b=", round(rr$b[i], 3), " p=", format(rr$pval[i], digits = 2))
}
res_df <- do.call(rbind, out)
write.csv(res_df, file.path(res, "grid/protein_interval_mr.csv"), row.names = FALSE)
funnel <- data.frame(stage = c("tested_proteins", "gene_outcome_pairs", "mr_completed_rows"),
                     count = c(length(PROTEINS), length(PROTEINS) * length(OUTCOMES),
                               sum(res_df$ok, na.rm = TRUE)))
write.table(funnel, file.path(res, "funnel/funnel_protein_interval.tsv"), sep = "\t", row.names = FALSE)
log("INTERVAL 蛋白通道 MR 完成 ✔ → results/grid/protein_interval_mr.csv")
log("局限：INTERVAL n=3,301 功效有限；deCODE(n=35,559) 为主源。")
