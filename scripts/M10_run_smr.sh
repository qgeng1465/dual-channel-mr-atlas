#!/bin/bash
# =============================================================================
# M10_run_smr.sh — 跑转录通道 SMR/HEIDI（106 探针 × 3 结局）
# 用法: bash scripts/M10_run_smr.sh
# 产物: results/grid/transcript_smr_heidi.csv（含 outcome 列）
# =============================================================================
set -euo pipefail
PROJ=/data/qiushuogeng/projects/dual-channel-mr-atlas
SMR=$PROJ/tools/smr
LDREF=$PROJ/data/ldref/1kg.v3/EUR
MA=$PROJ/data/smr
BESD=$MA/eqtlgen_trans
OUT=$PROJ/results/grid/transcript_smr_heidi.csv
TMP=$MA/run
mkdir -p "$TMP"

> "$OUT"
header_done=0
for on in t2d cad fbg; do
  echo "[$(date +%H:%M:%S)] running SMR trans_$on ..."
  $SMR --bfile "$LDREF" --gwas-summary "$MA/trans_$on.ma" --beqtl-summary "$BESD" \
       --peqtl-smr 5e-6 --disable-freq-ck --thread-num 6 \
       --out "$TMP/tsmr_$on" >/dev/null 2>&1
  if [ -s "$TMP/tsmr_$on.smr" ]; then
    if [ "$header_done" -eq 0 ]; then
      # 表头只写一次：awk NR==1 按每个输入文件单独计，若逐文件拼接会把表头重复 3 次；
      # 且原 " outcome" 列名带前导空格。改用首文件表头 + 无空格列名。
      awk -F'\t' -v oc="$on" 'NR==1{print $0 "\toutcome"; next} {print $0 "\t" oc}' \
          "$TMP/tsmr_$on.smr" > "$OUT"
      header_done=1
    else
      awk -F'\t' -v oc="$on" 'NR>1{print $0 "\t" oc}' "$TMP/tsmr_$on.smr" >> "$OUT"
    fi
  else
    echo "  WARNING: trans_$on 无输出" >&2
  fi
done
echo "[$(date +%H:%M:%S)] done -> $OUT"
