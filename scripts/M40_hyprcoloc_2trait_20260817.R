#!/usr/bin/env Rscript
# =============================================================================
# M40_hyprcoloc_2trait_20260817.R - hyprcoloc 2-trait 共定位交叉验证
# =============================================================================
# 目的：对 6 个 top 共定位位点（M34c），用独立的贝叶斯共定位实现 hyprcoloc
#   v1.0 (Foley et al. 2021, PLoS Genet) 做 2-trait（eQTLGen 表达 + 结局 GWAS）
#   交叉验证。hyprcoloc 为联合单因果变体后验（posterior_prob = 区域共享单一
#   因果变体的概率），不假定与 coloc 相同的"比例贡献"模型，因此是 coloc.susie
#   PP.H4 方法学独立性的检验。结果与 coloc.abf PP.H4 / coloc.susie PP.H4 对照。
# 输入：M34b 对齐产物 <scratch>/susie/{tag}_z1.txt, _z2.txt, _R.txt, _snp.txt
# 输出：results/m40_hyprcoloc_2trait_20260817.csv
# 诚实 caveat：hyprcoloc 假设每区域单一共享因果变体；多信号位点（如 LAMC1
#   8 CS vs 3 CS）该假设可能过强，其 posterior_prob 会系统性低于 coloc 的
#   PP.H4 —— 这正是"两种模型对同一证据的独立读数"，差异如实报告而非抹平。
# =============================================================================
suppressMessages(library(hyprcoloc))
TMP <- "/data/qiushuogeng/tmp/susie"
RES <- "/data/qiushuogeng/projects/dual-channel-mr-atlas/results"
N1  <- 31684
N_GWAS <- c(t2d = 655666, cad = 296525)
LOCI <- c(RBM6 = "t2d", CD101 = "t2d", CNNM2 = "cad", PLAUR = "cad", RIC8A = "cad", LAMC1 = "cad")
# 对照：coloc.susie ridge w=0.05 的 PP.H4 与 coloc.abf
ridge <- read.csv(file.path(RES, "m34c_coloc_susie_ridge_20260817.csv"))
ridge <- ridge[ridge$w == 0.05, c("symbol", "susie_pp4")]

out_rows <- list()
for (sym in names(LOCI)) {
  oc <- LOCI[[sym]]; tag <- tolower(sym); n2 <- N_GWAS[[oc]]
  z1 <- as.numeric(readLines(file.path(TMP, paste0(tag, "_z1.txt"))))
  z2 <- as.numeric(readLines(file.path(TMP, paste0(tag, "_z2.txt"))))
  R  <- as.matrix(read.table(file.path(TMP, paste0(tag, "_R.txt"))))
  snp <- readLines(file.path(TMP, paste0(tag, "_snp.txt")))
  n <- length(z1)
  if (n < 20 || nrow(R) != n) { cat(sym, "NS too small\n"); next }
  beta1 <- z1 / sqrt(N1); se1 <- rep(1 / sqrt(N1), n)
  beta2 <- z2 / sqrt(n2); se2 <- rep(1 / sqrt(n2), n)
  betas <- cbind(beta1, beta2); ses <- cbind(se1, se2)
  colnames(betas) <- colnames(ses) <- c("eQTLGen", toupper(oc))
  fit <- tryCatch(
    hyprcoloc(effect.est = betas, effect.se = ses,
              trait.names = c("eQTLGen", toupper(oc)),
              snp.id = snp, ld.matrix = R,
              prior.1 = 1e-4, prior.c = 0.02),
    error = function(e) { cat(sym, "ERROR:", conditionMessage(e), "\n"); NULL })
  if (is.null(fit) || is.null(fit$results)) next
  r <- fit$results[1, ]
  pp_susie <- ridge$susie_pp4[ridge$symbol == sym]
  abf <- read.csv(file.path(RES, paste0("coloc_full_", oc, "_20260815.csv")))
  pp_abf <- max(abf$pp4[abf$symbol == sym], na.rm = TRUE)
  out_rows[[sym]] <- data.frame(symbol = sym, outcome = oc, n_snps = n,
    hypr_pp_shared = r$posterior_prob, hypr_regional_prob = r$regional_prob,
    hypr_candidate_snp = r$candidate_snp,
    coloc_abf_pp4 = pp_abf, coloc_susie_pp4 = pp_susie,
    pp_gap_susie = pp_susie - r$posterior_prob)
  cat(sprintf("%-8s %-4s n=%4d hypr_pp=%.4f (snp=%s) | abf=%.4f susie=%.4f gap=%.4f\n",
              sym, oc, n, r$posterior_prob, r$candidate_snp, pp_abf, pp_susie,
              pp_susie - r$posterior_prob))
}
res <- do.call(rbind, out_rows)
write.csv(res, file.path(RES, "m40_hyprcoloc_2trait_20260817.csv"), row.names = FALSE)
cat(sprintf("\n== DONE M40: %d/%d loci hyprcoloc 2-trait cross-check ==\n", nrow(res), length(LOCI)))
