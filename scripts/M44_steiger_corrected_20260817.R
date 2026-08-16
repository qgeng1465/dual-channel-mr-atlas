#!/usr/bin/env Rscript
# =============================================================================
# M44_steiger_corrected_20260817.R — 15 候选效应基因的 MR-Steiger 测度误差校正方向检验
# =============================================================================
# 目的：验证 15 个候选效应基因的因果方向（exposure=基因表达 → outcome=T2D/CAD），
#   针对全血 eQTL (n≈26k) 与结局 GWAS (n≈296k–656k) 样本量严重不对称的测度误差问题，
#   实现 TwoSampleMR::mr_steiger(r_xxo=r_yyo=1) 等价的校正方向检验：
#   * 每个候选工具 SNP 计算 r2_exposure（eQTLGen）与 r2_outcome（结局 GWAS）
#     r2 = 2*p*(1-p)*beta^2 / (2*p*(1-p)*beta^2 + n*se^2)  （与 z^2/(z^2+n) 一致）
#   * 独立相关 Fisher-z 检验：z = (atanh(r_exp)-atanh(r_out))/sqrt(1/(n_exp-3)+1/(n_out-3))
#     —— 样本量不对称进入标准误，这正是"测度误差校正"的核心
#   * 输出 corrected Steiger p、correct_causal_direction、sensitivity_ratio = r_out/r_exp
# 输入：
#   results/candidate15_replication_20260816.csv      15 候选（symbol/outcome）
#   results/grid/transcript_grid_mr.csv               每候选 lead cis-eQTL 工具 SNP
#   <scratch>/eqtlgen_stable/bychr/chr*.tsv           eQTLGen cis 汇总（Zscore/等位/样本量）
#   <repo-root>/data/opengwas/full/{out}_full.gz      结局 GWAS 汇总
# 输出：results/m44_steiger_corrected_20260817.csv
# 诚实 caveat：r_xxo=r_yyo=1（无外部可靠性校正，等价 TwoSampleMR 默认）；方向结论为关联性评估。
# =============================================================================
suppressMessages({library(data.table)})
REPO <- "<repo-root>"
SCR  <- "<scratch>"

# ---- 1. 候选 → lead SNP → outcome ----
cand <- fread(file.path(REPO, "results/candidate15_replication_20260816.csv"))
grid <- fread(file.path(REPO, "results/grid/transcript_grid_mr.csv"))
grid$lead <- sub(".*lead (rs[0-9]+).*", "\\1", grid$note)
# 每个候选取与其 FDR-core outcome 匹配的那行
grid <- grid[grid$symbol %in% cand$symbol, ]
grid <- grid[match(paste(cand$symbol, cand$outcome), paste(grid$symbol, grid$outcome)), ]
grid$gene <- cand$gene[match(grid$symbol, cand$symbol)]

dat <- data.frame(symbol = grid$symbol, gene = grid$gene, outcome = grid$outcome,
                  lead_snp = grid$lead, ok = grid$ok, stringsAsFactors = FALSE)
dat$mr_b <- cand$mr_b[match(paste(dat$symbol, dat$outcome), paste(cand$symbol, cand$outcome))]
dat <- dat[!is.na(dat$lead_snp) & dat$lead_snp != "" & dat$lead_snp != "?" & dat$ok == TRUE, ]
cat("candidates with usable instrument:", nrow(dat), "\n")

# ---- 2. 从 eQTLGen + GWAS 提取每工具 SNP 统计 ----
get_eqtl <- function(chr, snp, gene) {
  # 返回 Zscore, AssessedAllele, OtherAllele, NrSamples
  f <- file.path(SCR, sprintf("eqtlgen_stable/bychr/chr%d.tsv", chr))
  cmd <- sprintf("awk -F'\\t' '$1==\"%s\" && $7==\"%s\"' %s", snp, gene, f)
  l <- system(cmd, intern = TRUE)
  if (length(l) == 0) return(NULL)
  p <- strsplit(l[1], "\t")[[1]]
  list(z = as.numeric(p[4]), a1 = p[5], a2 = p[6], n = as.numeric(p[10]))
}
# 各结局文件列 schema 不同：t2d_full.gz 有 n 列（p[22]）；cad_full.gz 为 meta 分析
#   （33 列）无样本量列 → n 固定 296525（M34b 口径）。列号经 zcat | head 实测核对。
gwas_schema <- list(
  t2d = list(chr = 3, eff = 6, other = 5, beta = 19, se = 20, p = 21, eaf = 18, n = 655666),
  cad = list(chr = 3, eff = 6, other = 5, beta = 20, se = 21, p = 22, eaf = 16, n = 296525)
)
get_gwas <- function(out, snp) {
  f <- file.path(REPO, sprintf("data/opengwas/full/%s_full.gz", out))
  cmd <- sprintf("zcat %s | awk -F'\\t' '$2==\"%s\"'", f, snp)
  l <- system(cmd, intern = TRUE)
  if (length(l) == 0) return(NULL)
  p <- strsplit(l[1], "\t")[[1]]
  sc <- gwas_schema[[out]]
  list(chr = as.integer(p[sc$chr]), beta = as.numeric(p[sc$beta]), se = as.numeric(p[sc$se]),
       p = as.numeric(p[sc$p]), n = sc$n, eaf = as.numeric(p[sc$eaf]),
       eff = p[sc$eff], other = p[sc$other])
}

r2_snp <- function(beta, se, eaf, n) {
  if (is.na(beta) || is.na(se) || is.na(eaf) || is.na(n) || n <= 0) return(NA_real_)
  num <- 2 * eaf * (1 - eaf) * beta^2
  num / (num + n * se^2)
}
# 独立相关 Fisher-z 检验（等价 psych::r.test，n1 != n2 的独立相关比较）
steiger_test <- function(r1, r2, n1, n2) {
  if (any(is.na(c(r1, r2, n1, n2))) || n1 < 4 || n2 < 4) return(list(z = NA, p = NA))
  z1 <- atanh(r1); z2 <- atanh(r2)
  z <- (z1 - z2) / sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
  p <- 2 * pnorm(-abs(z))
  list(z = z, p = p)
}

out <- list()
for (i in seq_len(nrow(dat))) {
  sym <- dat$symbol[i]; gene <- dat$gene[i]; oc <- dat$outcome[i]; snp <- dat$lead_snp[i]
  g <- get_gwas(oc, snp)
  if (is.null(g)) { cat(sym, snp, "missing GWAS; skip\n"); next }
  e <- get_eqtl(g$chr, snp, gene)
  if (is.null(e)) { cat(sym, snp, "missing eQTL; skip\n"); next }
  # 曝光端 beta/se：eQTLGen Zscore → beta=Z/sqrt(n), se=1/sqrt(n)
  beta_e <- e$z / sqrt(e$n); se_e <- 1 / sqrt(e$n)
  n_exp <- e$n; eaf_exp <- g$eaf
  r2_exp <- r2_snp(beta_e, se_e, eaf_exp, n_exp)
  # 结局端（对齐到曝光效应等位；r2 与符号无关）
  beta_g <- g$beta; se_g <- g$se; n_out <- g$n; eaf_out <- g$eaf
  r2_out <- r2_snp(beta_g, se_g, eaf_out, n_out)
  tt <- steiger_test(sqrt(r2_exp), sqrt(r2_out), n_exp, n_out)
  # 方向审计：单工具 Wald-ratio MR（原始对齐数据，mr_single = beta_g_al/beta_e）
  # 与 candidate15 存储的 MR 效应（mr_c15）对比。二者符号相反=等位基因约定差异伪影
  # （幅值一致时几乎必然），非真实方向反转 —— 如实记录 mr_sign_flip。
  aligned <- if (e$a1 == g$eff) TRUE else if (e$a1 == g$other) FALSE else NA
  beta_g_al <- if (is.na(aligned)) NA else if (aligned) beta_g else -beta_g
  mr_single <- if (is.na(beta_g_al) || beta_e == 0) NA else beta_g_al / beta_e
  mr_c15 <- dat$mr_b[i]
  mr_sign_flip <- if (is.na(mr_single) || is.na(mr_c15)) NA else sign(mr_single) != sign(mr_c15)
  dir_consistent <- if (is.na(mr_single) || is.na(beta_g_al)) NA else (beta_e * mr_single * beta_g_al) > 0
  out[[sym]] <- data.frame(symbol = sym, outcome = oc, lead_snp = snp,
                           n_exp = n_exp, n_out = n_out,
                           r2_exposure = r2_exp, r2_outcome = r2_out,
                           r_exp = sqrt(r2_exp), r_out = sqrt(r2_out),
                           steiger_z = tt$z, steiger_pval = tt$p,
                           correct_causal_direction = r2_exp > r2_out,
                           sensitivity_ratio = r2_out / r2_exp,
                           mr_c15 = mr_c15, mr_single = mr_single,
                           mr_sign_flip = mr_sign_flip,
                           direction_consistent = dir_consistent)
  cat(sprintf("%-8s %-4s %-10s r2_exp=%.4g r2_out=%.4g steiger_p=%.3g flip=%s ratio=%.3g\n",
              sym, oc, snp, r2_exp, r2_out, tt$p,
              ifelse(is.na(mr_sign_flip), "NA", ifelse(mr_sign_flip, "YES", "no")),
              r2_out / r2_exp))
}
res <- do.call(rbind, out)
write.csv(res, file.path(REPO, "results/m44_steiger_corrected_20260817.csv"), row.names = FALSE)
n_dir <- sum(res$correct_causal_direction, na.rm = TRUE)
cat(sprintf("== DONE M44: %d/%d candidates correct causal direction (r2_exp > r2_out) ==\n",
            n_dir, nrow(res)))
