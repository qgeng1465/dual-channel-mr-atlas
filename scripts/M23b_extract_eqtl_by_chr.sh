#!/bin/bash
# =============================================================================
# M23b_extract_eqtl_by_chr.sh — eQTLGen full 一次性抽取，按染色体分文件落盘
# =============================================================================
# 背景（2026-08-15）：M23 v1 按基因块 awk 抽取，3 进程并发扫同一 4.3G gz，
#   每块 12-14 分钟（14 块全扫）→ 太慢。改为**一次全扫 + 按染色体分文件明文落盘**
#   （<scratch>/eqtlgen_stable/bychr/chr{1..22}.tsv），
#   之后 M23 v2 每结局按染色体循环 fread 明文块（无重复解压、内存可控）。
# 列序（eQTLGen full，2026-08-15 逐列核实）：
#   1 Pvalue 2 SNP 3 SNPChr 4 SNPPos 5 Zscore 6 AssessedAllele 7 OtherAllele
#   8 Gene 9 GeneSymbol 10 GeneChr 11 GenePos 12 NrCohorts 13 NrSamples 14 FDR
# 输出列：SNP SNPChr SNPPos Zscore AssessedAllele OtherAllele Gene GeneChr GenePos NrSamples
# 用法：bash scripts/M23b_extract_eqtl_by_chr.sh   （走仲裁后台）
# =============================================================================
set -e
FULL=<scratch>/eqtlgen_stable/cis-eQTLs_full_20180905.txt.gz
OUT=<scratch>/eqtlgen_stable/bychr
mkdir -p "$OUT"
cd "$OUT"
[ -f DONE ] && { echo "already done"; exit 0; }

t0=$(date +%s)
echo "[$(date +%H:%M:%S)] 开始按染色体抽取 full -> $OUT"
zcat "$FULL" | awk -F'\t' '
NR==1 {
  hdr = $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$10"\t"$11"\t"$13
  next
}
{
  c = $3
  if (c ~ /^[0-9]+$/ && c >= 1 && c <= 22) {
    if (!(c in seen)) { print hdr > ("chr"c".tsv"); seen[c] = 1 }
    print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$10"\t"$11"\t"$13 >> ("chr"c".tsv")
  }
}'
echo "[$(date +%H:%M:%S)] 抽取完成，耗时 $(( $(date +%s) - t0 ))s"
echo "chr files:"
ls -la "$OUT"/chr*.tsv | awk '{print $5, $9}' | head -25
echo "done" > "$OUT/DONE"
