#!/bin/bash
# =============================================================================
# M34c_coloc_susie_ridge.sh — coloc.susie 的 ridge-shrunk LD 敏感性扫描
# =============================================================================
# 目的：M34b（max_iter=1000）在外样本 LD（1000G EUR）下 susie_rss 的 IBSS 仍可能
#   不收敛。susieR 官方对"外部 LD 与汇总统计不一致"的救济是 ridge 收缩
#   （susierss_diagnostic / susie_rss 文档）：R_w = (1-w)R + w*I。
#   本脚本对每个 top 位点扫描 w ∈ {0, 0.05, 0.10, 0.20}，记录收敛状态与 PP.H4，
#   如实评估 ridge 收缩是否能恢复收敛、PP.H4 是否与 coloc.abf 一致。
# 输入：M34b 的对齐产物 <scratch>/susie/{tag}_z1.txt, _z2.txt, _R.txt, _snp.txt
#   （由 scripts/M34b_coloc_susie_maxiter1000.sh 生成；每个位点 926–3235 SNPs，
#   1000G EUR 区域 LD，PSD 投影后与 eQTLGen/OpenGWAS z-score 对齐）
# 输出：results/m34c_coloc_susie_ridge_20260817.csv
# 用法：bash scripts/M34c_coloc_susie_ridge.sh
# =============================================================================
set -e
cd <repo-root>
R_ENV=<conda-root>/r-mr
TMP=<scratch>/susie
RES=results
OUT=$RES/m34c_coloc_susie_ridge_20260817.csv
echo "gene,symbol,outcome,mr_p,abf_pp4,w,susie_pp4,n_snps,conv_eqtl,conv_gwas,cs_eqtl,cs_gwas,note" > $OUT

# 6 个 top 候选：symbol|ensg|outcome（与 M34b 相同；对齐数据文件 tag=小写 symbol）
LOCI="RBM6|ENSG00000004534|t2d
CNNM2|ENSG00000148842|cad
PLAUR|ENSG00000011422|cad
CD101|ENSG00000134256|t2d
RIC8A|ENSG00000177963|cad
LAMC1|ENSG00000135862|cad"

for L in $LOCI; do
  IFS='|' read SYM ENSG OUTCOME <<< "$L"
  TAG="${SYM,,}"
  if [ ! -f $TMP/${TAG}_R.txt ]; then
    echo "跳过 $SYM：缺对齐数据 $TMP/${TAG}_R.txt（先跑 M34b）"
    continue
  fi
  echo "=== $SYM × $OUTCOME ==="
  PATH=$R_ENV/bin:$PATH Rscript - "$SYM" "$ENSG" "$OUTCOME" "$TAG" >> $OUT 2> $TMP/${TAG}_ridge_Rerr.log <<'REOF'
args <- commandArgs(trailingOnly=TRUE)
sym <- args[1]; ensg <- args[2]; out <- args[3]; tag <- args[4]
suppressMessages({library(coloc); library(susieR)})
TMP <- "<scratch>/susie"
RES <- "<repo-root>/results"
z1 <- as.numeric(readLines(file.path(TMP, paste0(tag, "_z1.txt"))))
z2 <- as.numeric(readLines(file.path(TMP, paste0(tag, "_z2.txt"))))
R  <- as.matrix(read.table(file.path(TMP, paste0(tag, "_R.txt"))))
snp <- readLines(file.path(TMP, paste0(tag, "_snp.txt")))
n <- length(z1)
if (n < 20 || nrow(R) != n) {
  cat(sprintf("%s,%s,%s,NA,NA,NA,NA,%d,NA,NA,NA,NA,SKIP NS=%d\n", ensg, sym, out, n, n)); quit(save="no")
}
colnames(R) <- rownames(R) <- snp
# coloc.abf 对照：读主表 coloc_full 的 PP.H4（真实 beta/se 口径，独立重算验证）
cf <- read.csv(file.path(RES, paste0("coloc_full_", out, "_20260815.csv")))
rowm <- cf[cf$gene == ensg & cf$ok == TRUE, ]
abf <- if (nrow(rowm) >= 1) max(rowm$pp4, na.rm = TRUE) else NA
mr_p_main <- if (nrow(rowm) >= 1) min(rowm$mr_p, na.rm = TRUE) else NA
n_gwas <- if (out == "t2d") 655666 else 296525
for (w in c(0, 0.05, 0.10, 0.20)) {
  Rw <- (1 - w) * R + w * diag(n)
  f1 <- tryCatch(susie_rss(z1, Rw, n = 31684, L = 10, estimate_residual_variance = FALSE,
                           coverage = 0.9, max_iter = 1000),
                 error = function(e) NULL)
  f2 <- tryCatch(susie_rss(z2, Rw, n = n_gwas, L = 10, estimate_residual_variance = FALSE,
                           coverage = 0.9, max_iter = 1000),
                 error = function(e) NULL)
  if (is.null(f1) || is.null(f2)) {
    cat(sprintf("%s,%s,%s,%.4g,%.4f,%.3f,NA,%d,NA,NA,NA,NA,FAILED susie_rss\n",
                ensg, sym, out, mr_p_main, abf, w, n)); next
  }
  su <- tryCatch(coloc.susie(f1, f2, p12 = 1e-5), error = function(e) NULL)
  pp_susie <- if (!is.null(su) && "summary" %in% names(su))
    max(as.numeric(su$summary[["PP.H4.abf"]]), na.rm = TRUE) else NA
  cat(sprintf("%s,%s,%s,%.4g,%.4f,%.3f,%.4f,%d,%s,%s,%d,%d,%s\n",
              ensg, sym, out, mr_p_main, abf, w, pp_susie, n,
              f1$converged, f2$converged, length(f1$sets$cs), length(f2$sets$cs),
              ifelse(is.null(su), "FAILED coloc.susie", "OK")))
}
REOF
done
echo "完成 -> $OUT"
cat $OUT
