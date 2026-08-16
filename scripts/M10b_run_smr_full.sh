#!/bin/bash
# =============================================================================
# M10b_run_smr_full.sh — 跑转录通道 SMR/HEIDI 全量（759 探针 × 3 结局）
# 前提：M10b_build_heidi_full_inputs.R 已构建 trans_esd_full/ + trans_flist_full.txt
#       + trans_{t2d,cad,fbg}_full.ma
# 用法: bash scripts/M10b_run_smr_full.sh
# 产物: results/grid/transcript_smr_heidi_full.csv（含 outcome 列）
# 与 M10_run_smr.sh 完全一致（SMR 参数、awk 追加 outcome 列），仅输入换 _full。
# =============================================================================
set -euo pipefail
PROJ=/data/qiushuogeng/projects/dual-channel-mr-atlas
SMR=$PROJ/tools/smr
LDREF=$PROJ/data/ldref/1kg.v3/EUR
MA=$PROJ/data/smr
BESD=$MA/eqtlgen_trans_full
OUT=$PROJ/results/grid/transcript_smr_heidi_full.csv
TMP=$MA/run
mkdir -p "$TMP"
THREADS=${THREADS:-8}   # 可被 arbitrate 环境覆盖（小任务并跑时降核）

echo "[$(date +%H:%M:%S)] Step1: 构建全量 BESD (759 探针) ..."
$SMR --eqtl-flist "$MA/trans_flist_full.txt" --make-besd \
     --out "$BESD" --thread-num $THREADS >/dev/null 2>&1 || {
  echo "  ERROR: --make-besd 失败，查看输出。" >&2
  $SMR --eqtl-flist "$MA/trans_flist_full.txt" --make-besd --out "$BESD" --thread-num $THREADS 2>&1 | tail -20 >&2
  exit 1
}
ls -la "$BESD".besd "$BESD".esi "$BESD".epi 2>/dev/null
echo "[$(date +%H:%M:%S)] BESD 构建完成"

> "$OUT"
header_done=0
for on in t2d cad fbg; do
  echo "[$(date +%H:%M:%S)] running SMR trans_${on}_full ..."
  $SMR --bfile "$LDREF" --gwas-summary "$MA/trans_${on}_full.ma" --beqtl-summary "$BESD" \
       --peqtl-smr 5e-6 --disable-freq-ck --thread-num $THREADS \
       --out "$TMP/tsmr_full_$on" >/dev/null 2>&1 || echo "  WARNING: SMR trans_${on}_full 返回非零" >&2
  if [ -s "$TMP/tsmr_full_$on.smr" ]; then
    if [ "$header_done" -eq 0 ]; then
      # 表头只写一次：awk NR==1 按每个输入文件单独计，若逐文件拼接会把表头重复 3 次；
      # 且原 " outcome" 列名带前导空格。改用首文件表头 + 无空格列名。
      awk -F'\t' -v oc="$on" 'NR==1{print $0 "\toutcome"; next} {print $0 "\t" oc}' \
          "$TMP/tsmr_full_$on.smr" > "$OUT"
      header_done=1
    else
      awk -F'\t' -v oc="$on" 'NR>1{print $0 "\t" oc}' "$TMP/tsmr_full_$on.smr" >> "$OUT"
    fi
    echo "  $(wc -l < "$TMP/tsmr_full_$on.smr") 行 SMR 结果"
  else
    echo "  WARNING: trans_${on}_full 无输出" >&2
  fi
done
echo "[$(date +%H:%M:%S)] done -> $OUT"
wc -l "$OUT"
