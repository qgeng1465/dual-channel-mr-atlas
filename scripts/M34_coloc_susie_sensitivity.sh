#!/bin/bash
# =============================================================================
# M34_coloc_susie_sensitivity.sh — 多因果变异（coloc.susie）敏感性分析
# =============================================================================
# 目的：AJHG 底线（关联需精细定位/机制跟进）——对 top 候选用 coloc.susie 允许
#   多个因果变异，验证 PP.H4（coloc.abf 单因果假设）稳健。若 coloc.susie 仍支持
#   共享因果信号，则候选的"共因变异"证据在更宽松假设下成立。
# 流程（每基因）：
#   1. 区域（基因 pos ±600kb, hg19）1kg EUR --r square → LD 矩阵
#   2. eQTLGen 该基因 cis eQTL（Zscore）+ OpenGWAS 区域 beta/se → harmonize → z1,z2
#   3. R: susie_rss(eQTL) + susie_rss(GWAS) + coloc.susie → PP.H4.susie
# 输出：results/m34_coloc_susie_20260816.csv
# 用法：bash scripts/M34_coloc_susie_sensitivity.sh
# =============================================================================
set -e
cd /data/qiushuogeng/projects/dual-channel-mr-atlas
R_ENV=/data/gengqiushuo/home/miniconda3/envs/r-mr
PLINK=tools/plink
EUR=data/ldref/1kg.v3/EUR
TMP=/data/qiushuogeng/tmp/susie
mkdir -p $TMP
OUT=results/m34_coloc_susie_20260816.csv
echo "gene,symbol,outcome,mr_p,abf_pp4,susie_pp4,n_snps,note" > $OUT

# 6 个 top 候选：symbol|ensg|outcome|chr|pos(hg19)
LOCI="RBM6|ENSG00000004534|t2d|3|50057459
CNNM2|ENSG00000148842|cad|10|104758209
PLAUR|ENSG00000011422|cad|19|44162473
CD101|ENSG00000134256|t2d|1|117561774
RIC8A|ENSG00000177963|cad|11|211312
LAMC1|ENSG00000135862|cad|1|183053661"

for L in $LOCI; do
  IFS='|' read SYM ENSG OUTCOME CHR POS <<< "$L"
  TAG="${SYM,,}"
  echo "=== $SYM × $OUT (chr$CHR:$POS) ==="
  LO=$((POS-600000)); HI=$((POS+600000))
  [ $LO -lt 1 ] && LO=1
  # 1) LD 矩阵
  $PLINK --bfile $EUR --chr $CHR --from-bp $LO --to-bp $HI --r square --keep-allele-order \
    --out $TMP/${TAG}_R 2>/dev/null
  # 2) Python：提取 eQTL+GWAS、harmonize、输出 z1/z2/R 子矩阵 + 状态
  python3 - "$ENSG" "$SYM" "$OUTCOME" "$CHR" "$POS" "$TAG" > $TMP/${TAG}_aligned.log 2>&1 <<'PYEOF'
import sys, os, numpy as np
ensg, sym, out, chrS, pos, tag = sys.argv[1:7]
chrS, pos = int(chrS), int(pos)
TMP = "/data/qiushuogeng/tmp/susie"
LO, HI = max(pos-600000, 1), pos+600000
bim = {}
with open("/data/qiushuogeng/projects/dual-channel-mr-atlas/data/ldref/1kg.v3/EUR.bim") as f:
    for line in f:
        p = line.split()
        if int(p[0])==chrS and LO <= int(p[3]) <= HI:
            bim[p[1]] = int(p[3])
rs_order = list(bim.keys())
rs2idx = {rs: i for i, rs in enumerate(rs_order)}
R = np.fromfile(f"{TMP}/{tag}_R.ld", sep=" ")
n = len(rs_order)
if R.size != n*n:
    print(f"ALIGN_OK=0 LD_SIZE_MISMATCH {R.size}!={n*n}"); sys.exit(0)
R = R.reshape(n, n)
# eQTL：该基因 cis 行
eq = {}
with os.popen(f"awk -F'\\t' '$7==\"{ensg}\"' /data/qiushuogeng/tmp/eqtlgen_stable/bychr/chr{chrS}.tsv") as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < 6: continue
        snp, z = p[0], float(p[3])
        if snp in rs2idx:
            eq[snp] = (p[4], p[5], z)
# GWAS：区域 rsid 匹配（t2d/cad 用 hm_ 列）
g = {}
with os.popen(f"zcat /data/qiushuogeng/projects/dual-channel-mr-atlas/data/opengwas/full/{out}_full.gz") as f:
    hdr = f.readline().rstrip("\n").split("\t")
    ci = {c: i for i, c in enumerate(hdr)}
    ri, ei, oi, bi, sei = ci["hm_rsid"], ci["hm_effect_allele"], ci["hm_other_allele"], ci["hm_beta"], ci["standard_error"]
    for line in f:
        p = line.rstrip("\n").split("\t")
        if p[ri] in rs2idx:
            try: beta, se = float(p[bi]), float(p[sei])
            except ValueError: continue
            if se > 0 and beta == beta:
                g[p[ri]] = (p[ei], p[oi], beta, se)
# harmonize：对齐到 eQTL 的 AssessedAllele
snps, z1, z2 = [], [], []
n_align = n_flip = n_pal = n_miss = 0
for snp in rs_order:
    if snp not in eq or snp not in g: n_miss += 1; continue
    a1, a2, ze = eq[snp]
    ge, go, bg, se = g[snp]
    A1, A2 = a1.upper(), a2.upper(); GE, GO = ge.upper(), go.upper()
    z = bg / se
    if A1 == GE and A2 == GO: z2v = z; n_align += 1
    elif A1 == GO and A2 == GE: z2v = -z; n_flip += 1
    else: n_pal += 1; continue
    snps.append(snp); z1.append(ze); z2.append(z2v)
if len(snps) < 20:
    print(f"ALIGN_OK=0 NS={len(snps)} align={n_align} flip={n_flip} pal={n_pal} miss={n_miss}"); sys.exit(0)
idx = [rs2idx[s] for s in snps]
Rsub = R[np.ix_(idx, idx)]
Rsub = (Rsub + Rsub.T) / 2
for i in range(len(idx)): Rsub[i, i] = 1.0
# PSD 投影：负特征值截断（保证 susie_rss 可用的半正定相关阵）
w, V = np.linalg.eigh(Rsub)
w = np.clip(w, 0, None)
Rsub = V @ np.diag(w) @ V.T
Rsub = (Rsub + Rsub.T) / 2
for i in range(len(idx)): Rsub[i, i] = 1.0
np.savetxt(f"{TMP}/{tag}_z1.txt", np.array(z1)); np.savetxt(f"{TMP}/{tag}_z2.txt", np.array(z2))
np.savetxt(f"{TMP}/{tag}_R.txt", Rsub)
with open(f"{TMP}/{tag}_snp.txt", "w") as f: f.write("\n".join(snps))
print(f"ALIGN_OK=1 NS={len(snps)} align={n_align} flip={n_flip} pal={n_pal} miss={n_miss}")
PYEOF
  if ! grep -q "ALIGN_OK=1" $TMP/${TAG}_aligned.log; then
    echo "  跳过：$(cat $TMP/${TAG}_aligned.log | tail -1)"
    continue
  fi
  echo "  $(grep 'ALIGN_OK=1' $TMP/${TAG}_aligned.log)"
  # 3) R: susie_rss + coloc.susie + coloc.abf 对照
  PATH=$R_ENV/bin:$PATH Rscript - "$SYM" "$ENSG" "$OUTCOME" "$TAG" >> $OUT 2> $TMP/${TAG}_Rerr.log <<'REOF'
args <- commandArgs(trailingOnly=TRUE)
sym <- args[1]; ensg <- args[2]; out <- args[3]; tag <- args[4]
suppressMessages({library(coloc); library(susieR)})
TMP <- "/data/qiushuogeng/tmp/susie"
RES <- "/data/qiushuogeng/projects/dual-channel-mr-atlas/results"
z1 <- as.numeric(readLines(file.path(TMP, paste0(tag, "_z1.txt"))))
z2 <- as.numeric(readLines(file.path(TMP, paste0(tag, "_z2.txt"))))
R  <- as.matrix(read.table(file.path(TMP, paste0(tag, "_R.txt"))))
snp <- readLines(file.path(TMP, paste0(tag, "_snp.txt")))
n <- length(z1)
colnames(R) <- rownames(R) <- snp
# EPV=FALSE + 显式 n：外样本 LD（1000G EUR vs eQTLGen/OpenGWAS 汇总统计）标准做法。
# 不传 n 时 susie_rss 假设 n=Inf → "prior variance unreasonably large"（M6.8/实测）。
# 外样本 LD（1000G EUR vs eQTLGen/OpenGWAS）下 susie_rss 的 IBSS 常不收敛（实测 RBM6
# 500 迭代仍 conv=FALSE）。取 max_iter=200 为统一标准，收敛标志如实写入 note；论文以
# "支持性敏感性证据 + 未完全收敛披露"呈现，不作为主表结论。
f1 <- tryCatch(susie_rss(z1, R, n = 31684, L = 10, estimate_residual_variance = FALSE, coverage = 0.9, max_iter = 200),
               error = function(e) NULL)
f2 <- tryCatch(susie_rss(z2, R, n = if (out == "t2d") 655666 else 296525, L = 10,
                         estimate_residual_variance = FALSE, coverage = 0.9, max_iter = 200),
               error = function(e) NULL)
if (is.null(f1) || is.null(f2)) {
  cat(sprintf("%s,%s,%s,NA,NA,NA,%d,FAILED susie_rss\n", ensg, sym, out, n)); quit(save="no")
}
su <- tryCatch(coloc.susie(f1, f2, p12 = 1e-5), error = function(e) NULL)
if (is.null(su) || !("summary" %in% names(su))) {
  cat(sprintf("%s,%s,%s,NA,NA,NA,%d,FAILED coloc.susie\n", ensg, sym, out, n)); quit(save="no")
}
# coloc.susie 的 summary 是 data.table：每个 (eQTL 可信集 × GWAS 可信集) 对一行。
# 汇总口径 = 所有配对的最大 PP.H4（任一对因果信号共享 → 共定位成立）。
pp_susie <- max(as.numeric(su$summary[["PP.H4.abf"]]), na.rm = TRUE)
# coloc.abf 对照：直接读主表 coloc_full 的 PP.H4（真实 beta/se 口径，已由独立重算验证），
# 不用 z-score 近似——z 近似（beta=z, varbeta=1, sdY=1）会给出错误量级的 PP.H4（实测 RBM6
# 0.1011 vs 主表 0.945），因 coloc 的 Bayes factor 依赖真实效应尺度（se=1/√N 远小于 1）。
cf <- read.csv(file.path(RES, paste0("coloc_full_", out, "_20260815.csv")))
rowm <- cf[cf$gene == ensg & cf$ok == TRUE, ]
abf <- if (nrow(rowm) >= 1) max(rowm$pp4, na.rm = TRUE) else NA
mr_p_main <- if (nrow(rowm) >= 1) min(rowm$mr_p, na.rm = TRUE) else NA
nsig1 <- length(f1$sets$cs); nsig2 <- length(f2$sets$cs)
# 诚实披露：外样本 LD（1000G EUR）下 susie_rss 的 IBSS 可能不收敛，收敛标志如实写入 note。
conv <- paste0("conv_eqtl=", f1$converged, "_gwas=", f2$converged)
cat(sprintf("%s,%s,%s,%.4g,%.4f,%.4f,%d,OK susie_cs_eqtl=%d susie_cs_gwas=%d %s\n",
            ensg, sym, out, mr_p_main, abf, pp_susie, n, nsig1, nsig2, conv))
REOF
done
echo "完成 -> $OUT"
cat $OUT
