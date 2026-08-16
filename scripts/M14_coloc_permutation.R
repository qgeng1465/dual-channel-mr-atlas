#!/usr/bin/env Rscript
# =============================================================================
# M14_coloc_permutation.R — E1：coloc PP.H4 阈值负对照/permutation 标定（探索性）
# =============================================================================
# 目的：回答"PP.H4≥0.8 的 strong coloc 调用在零假设下（无共享因果变异、仅 LD 结构
#   + 同一区域独立信号）出现多少次？" → 经验假阳性率 + 期望 FP 计数。
# 方法（与 M5 完全同源的数据构建）：
#   - 对象：transcript_coloc_hits.csv 的 106 个 strong coloc 对
#   - 数据：_coloc_gwas/{OUT_N}_{gene}.rds（GWAS 区域）+ eQTLGen 全量 cis（eQTL 侧）
#   - 对齐：M5 harmonize 逻辑（等位翻转 + palindromic eAF 解析，MAF 0.01-0.99）
#   - permutation：把 eQTL beta 向量跨 SNP 打乱（sample），保持 GWAS 信号与 MAF/
#     varbeta 不变 → 打破 eQTL-GWAS 共享因果对齐，保留两套信号的边缘分布
#   - 对每个打乱数据集重跑 coloc.abf(p12=1e-5)，记录 PP.H4
#   - 经验 FP 率 = 全部 permutation 中 PP.H4≥0.8 的比例；期望 FP = 106 × FP率
# 纪律：探索性敏感性，独立输出，不触碰预注册主流程。B 可参数化（默认 100）。
# 用法：PATH=$R_ENV/bin:$PATH Rscript scripts/M14_coloc_permutation.R [B]
# 输出：results/grid/coloc_permutation_calib.csv + results/coloc_permutation_20260813.md
# =============================================================================
suppressMessages({
  library(data.table)
  suppressPackageStartupMessages(library(coloc))
})
B <- if (length(commandArgs(trailingOnly = TRUE)) > 0) as.integer(commandArgs(trailingOnly = TRUE)[1]) else 100
proj <- "/data/qiushuogeng/projects/dual-channel-mr-atlas"
res  <- file.path(proj, "results")
gdir <- file.path(res, "grid")
FULL <- file.path(proj, "data/eqtlgen/cis-eQTLs_full_20180905.txt.gz")
SIG  <- file.path(proj, "data/eqtlgen/cis-eQTL-significant.txt.gz")
AF   <- file.path(proj, "data/eqtlgen/SNP_AF.txt.gz")
cfg <- jsonlite::fromJSON(file.path(res, "config.json"))
S_CASE <- c(t2d = 61714 / 655666, cad = 34541 / 296525)
TYPE <- c(t2d = "cc", cad = "cc", fbg = "quant")
OUT_N <- c(t2d = "GCST006867", cad = "GCST005194", fbg = "GCST005186")
log <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n", sep = "")
set.seed(20260813)

hits <- fread(file.path(gdir, "transcript_coloc_hits.csv"))
pairs <- unique(hits[, .(gene, outcome)])
log("strong coloc 对: ", nrow(pairs), " | B = ", B)

# 基因坐标 + eQTL cis 数据（同 M5 Phase B/C）
gcoord <- unique(fread(cmd = paste0("zcat ", SIG), sep = "\t", header = TRUE, nThread = 4,
                       select = c("Gene", "GeneChr", "GenePos")))
gcoord <- gcoord[!duplicated(Gene)]
pairs <- merge(pairs, gcoord, by.x = "gene", by.y = "Gene", all.x = TRUE)
pairs <- pairs[!is.na(GeneChr)]
log("坐标就绪: ", nrow(pairs))

genes_file <- file.path(gdir, "_coloc_genes.txt")
fwrite(unique(data.table(g = c(pairs$gene, hits$symbol))), genes_file, col.names = FALSE)
cmd <- paste0("zcat ", FULL, " | awk -F'\\t' 'NR==FNR{keep[$1]=1; next} FNR==1 || $8 in keep || $9 in keep {print}' ",
              genes_file, " -")
eqtl_full <- fread(cmd = cmd, sep = "\t", header = TRUE, nThread = 4)
eqtl_full <- eqtl_full[SNPChr == GeneChr & abs(SNPPos - GenePos) <= cfg$instrument$cis_window_kb * 1000]
af <- fread(cmd = paste0("zcat ", AF), sep = "\t", header = TRUE, nThread = 4)
setnames(af, c("SNP","hg19_chr","hg19_pos","AlleleA","AlleleB","allA_total","allAB_total","allB_total","AlleleB_all"))
eqtl_full <- merge(eqtl_full, af[, .(SNP, AlleleA, AlleleB, AlleleB_all)], by = "SNP", all.x = TRUE)
eqtl_full[, eaf := fcase(AssessedAllele == AlleleA, 1 - AlleleB_all,
                         AssessedAllele == AlleleB, AlleleB_all, default = NA_real_)]
eqtl_full <- eqtl_full[!is.na(eaf)]
eqtl_full[, `:=`(beta = Zscore / sqrt(NrSamples), se = 1 / sqrt(NrSamples))]
eqtl_full <- unique(eqtl_full[, .(Gene, SNP, AssessedAllele, OtherAllele, eaf, beta, se, NrSamples)])
log("eQTL cis 子集: ", nrow(eqtl_full), " 行")

harmonize <- function(e, g) {
  e2 <- copy(e[, .(SNP, AssessedAllele, OtherAllele, eaf_e = eaf, beta_e = beta, se_e = se, N_e = NrSamples)])
  g2 <- copy(g[, .(rsid = rsid, ea = ea, nea = nea, eaf_g = as.numeric(eaf),
                   beta_g = as.numeric(beta), se_g = as.numeric(se), N_g = as.numeric(n))])
  m <- merge(e2, g2, by.x = "SNP", by.y = "rsid")
  if (nrow(m) == 0) return(data.table())
  A <- m$AssessedAllele; O <- m$OtherAllele; ea <- m$ea; nea <- m$nea
  gb <- m$beta_g
  ok <- A == ea; flip <- !ok & A == nea & O == ea; gb[flip] <- -gb[flip]
  pal <- !ok & !flip & (((A == "A" & O == "T") | (A == "T" & O == "A")) |
                        ((A == "C" & O == "G") | (A == "G" & O == "C")))
  p_flip <- pal & ((m$eaf_e > 0.58 & m$eaf_g < 0.42) | (m$eaf_e < 0.42 & m$eaf_g > 0.58))
  gb[p_flip] <- -gb[p_flip]
  keep <- ok | flip | p_flip
  out <- data.table(snp = m$SNP[keep], maf = pmin(m$eaf_e, 1 - m$eaf_e)[keep],
                    beta_e = m$beta_e[keep], varbeta_e = m$se_e[keep]^2, N_e = m$N_e[keep],
                    beta_g = gb[keep], varbeta_g = m$se_g[keep]^2, N_g = m$N_g[keep])
  out[maf >= 0.01 & maf <= 0.99 & !is.na(beta_e) & !is.na(beta_g) &
      varbeta_e > 0 & varbeta_g > 0]
}
ppv <- function(s, nm) {
  for (k in c(nm, paste0(nm, ".abf"))) if (length(s) && !is.na(s[k])) return(unname(as.numeric(s[k])))
  NA_real_
}

out_rows <- list()
for (i in seq_len(nrow(pairs))) {
  gene <- pairs$gene[i]; on <- pairs$outcome[i]
  gf <- file.path(gdir, "_coloc_gwas", paste0(OUT_N[[on]], "_", gene, ".rds"))
  if (!file.exists(gf)) { out_rows[[length(out_rows) + 1]] <- data.table(gene = gene, outcome = on, nsnp = 0, n_perm = 0, fp_ge08 = NA, obs_h4 = NA, note = "无 GWAS 缓存"); next }
  g <- readRDS(gf)
  e <- eqtl_full[Gene == gene]
  if (is.null(g) || nrow(g) == 0 || nrow(e) == 0) { out_rows[[length(out_rows) + 1]] <- data.table(gene = gene, outcome = on, nsnp = 0, n_perm = 0, fp_ge08 = NA, obs_h4 = NA, note = "eQTL/GWAS 空"); next }
  m <- harmonize(e, g)
  if (nrow(m) < 10) { out_rows[[length(out_rows) + 1]] <- data.table(gene = gene, outcome = on, nsnp = nrow(m), n_perm = 0, fp_ge08 = NA, obs_h4 = NA, note = paste0("交集后 ", nrow(m), " 变异")); next }
  # 观测 PP.H4（应 ≈ coloc 结果，作一致性校验）
  d1o <- list(snp = m$snp, type = "quant", N = m$N_e, beta = m$beta_e, varbeta = m$varbeta_e, MAF = m$maf)
  d2o <- if (TYPE[[on]] == "cc")
    list(snp = m$snp, type = "cc", N = m$N_g, s = S_CASE[[on]], beta = m$beta_g, varbeta = m$varbeta_g, MAF = m$maf)
  else list(snp = m$snp, type = "quant", N = m$N_g, beta = m$beta_g, varbeta = m$varbeta_g, MAF = m$maf)
  ro <- tryCatch(coloc.abf(d1o, d2o, p12 = 1e-5), error = function(e2) NULL)
  obs_h4 <- if (!is.null(ro)) ppv(ro$summary, "PP.H4") else NA
  # permutation：打乱 eQTL beta 跨 SNP，保持 GWAS 与 varbeta/MAF 不变
  n_ge05 <- n_ge08 <- n_ge09 <- 0L
  beta_pool <- m$beta_e
  for (b in seq_len(B)) {
    perm <- sample(beta_pool)
    d1 <- list(snp = m$snp, type = "quant", N = m$N_e, beta = perm, varbeta = m$varbeta_e, MAF = m$maf)
    rp <- tryCatch(coloc.abf(d1, d2o, p12 = 1e-5), error = function(e2) NULL)
    if (!is.null(rp)) {
      h <- ppv(rp$summary, "PP.H4")
      if (!is.na(h)) { if (h >= 0.5) n_ge05 <- n_ge05 + 1L; if (h >= 0.8) n_ge08 <- n_ge08 + 1L; if (h >= 0.9) n_ge09 <- n_ge09 + 1L }
    }
  }
  out_rows[[length(out_rows) + 1]] <- data.table(
    gene = gene, outcome = on, nsnp = nrow(m), n_perm = B,
    fp_ge05 = n_ge05, fp_ge08 = n_ge08, fp_ge09 = n_ge09,
    fp_rate08 = n_ge08 / B, obs_h4 = obs_h4,
    note = if (is.na(obs_h4)) "观测 coloc 失败" else "")
  if (i %% 10 == 0) log("  完成 ", i, "/", nrow(pairs), " 对（FP≥0.8: ", n_ge08, "/", B, "）")
}
res_tab <- rbindlist(out_rows, fill = TRUE)
fwrite(res_tab, file.path(gdir, "coloc_permutation_calib.csv"))
log("总表已写 results/grid/coloc_permutation_calib.csv | 行数 ", nrow(res_tab))

# 汇总
run_rows <- res_tab[!is.na(fp_rate08)]
tot_perm <- sum(run_rows$n_perm)
n_strong <- nrow(run_rows)
r <- function(col) sum(run_rows[[col]]) / tot_perm
cat("\n=== E1 结果摘要 ===\n")
cat("strong 调用数:", n_strong, "\n")
cat("permutation 总数:", tot_perm, "\n")
cat("零假设下 PP.H4≥0.5 经验 FP 率:", r("fp_ge05"), "\n")
cat("零假设下 PP.H4≥0.8 经验 FP 率:", r("fp_ge08"), "\n")
cat("零假设下 PP.H4≥0.9 经验 FP 率:", r("fp_ge09"), "\n")
cat("期望 FP≥0.8 在 ", n_strong, " 个调用中:", n_strong * r("fp_ge08"), "\n")
log(sprintf("汇总: %d strong 对, %d perms, FP率 0.5/0.8/0.9 = %.4f/%.4f/%.4f",
            n_strong, tot_perm, r("fp_ge05"), r("fp_ge08"), r("fp_ge09")))
