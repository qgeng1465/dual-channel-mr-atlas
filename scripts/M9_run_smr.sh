#!/bin/bash
# =============================================================================
# M9_run_smr.sh — 跑蛋白通道 SMR/HEIDI（5 蛋白 × 3 结局）
# 用法: bash scripts/M9_run_smr.sh
# 产物: results/grid/protein_smr_heidi.csv
# =============================================================================
set -euo pipefail
PROJ=<repo-root>
SMR=$PROJ/tools/smr
LDREF=$PROJ/data/ldref/1kg.v3/EUR
MA=$PROJ/data/smr
BESD=$MA/decode_pqtl
OUT=$PROJ/results/grid/protein_smr_heidi.csv
TMP=$PROJ/data/smr/run
mkdir -p "$TMP"

> "$OUT"
for on in t2d cad fbg; do
  echo "[$(date +%H:%M:%S)] running SMR $on ..."
  $SMR --bfile "$LDREF" --gwas-summary "$MA/$on.ma" --beqtl-summary "$BESD" \
       --peqtl-smr 5e-6 --disable-freq-ck --thread-num 6 \
       --out "$TMP/smr_$on" >/dev/null 2>&1
  # .smr 有表头，直接附加；按结局加列
  if [ -s "$TMP/smr_$on.smr" ]; then
    awk -F'\t' -v oc="$on" 'NR==1{$0=$0"\t outcome"; print} NR>1{$0=$0"\t"oc; print}' \
        "$TMP/smr_$on.smr" >> "$OUT"
  else
    echo "WARNING: $on 无输出" >&2
  fi
done
echo "[$(date +%H:%M:%S)] done -> $OUT"
