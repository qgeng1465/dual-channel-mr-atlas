#!/bin/bash
# =============================================================================
# M1_decode_subset.sh — deCODE 血浆 pQTL 分析基因集 cis 窗流式裁剪（实现版）
# =============================================================================
# 原理：deCODE 逐蛋白文件（~953MB，hg38，全基因组变异）只保留分析基因集
#       cis 窗（TSS±1Mb）内变异，写入 data/decode/sub/，随后按 --delete-original 删原文件。
# 用法：bash scripts/M1_decode_subset.sh [--delete-original]
#       无参数：只提取已完成的 deCODE 文件，保留原文件；
#       --delete-original：提取成功且 cis 输出非空后删除原文件（README 体积红线）。
# 坐标：hg38 GRCh38，ENSEMBL REST 实查（2026-08-06），TSS=strand 起点/终点。
# =============================================================================
set -uo pipefail
PROJ=<repo-root>
DECODE=$PROJ/data/decode
SUB=$DECODE/sub
DEL_ORIG=0; [ "${1:-}" = "--delete-original" ] && DEL_ORIG=1
mkdir -p "$SUB"

# 分析基因集 → (chr, TSS, ±Mb)。strand - 的 TSS=end。
# 格式: gene|chr|tss|window_kb
GENES=(
  "PCSK9|chr1|55039445|1000"
  "HMGCR|chr5|75336329|1000"
  "ANGPTL3|chr1|62597464|1000"
  "APOC3|chr11|116827019|1000"
  "APOB|chr2|21044075|1000"
  "LDLR|chr19|11089418|1000"
  "CETP|chr16|56961922|1000"
  "NPC1L1|chr7|44541330|1000"
  "GLP1R|chr6|39048562|1000"
  "DPP4|chr2|162074639|1000"
  "INSR|chr19|7294443|1000"
  "PCK1|chr20|57546220|1000"
  "GCG|chr2|162152404|1000"
)

# 文件名: SeqId_Chr_GeneName_ProteinName.txt.gz → 按 GeneName 匹配
get_coord() {  # $1=gene; 输出 "chr tss win" 或空
  local g=$1 e
  for e in "${GENES[@]}"; do
    local name="${e%%|*}"; [ "$name" = "$g" ] || continue
    echo "$e" | awk -F'|' '{print $2, $3, $4}'
    return 0
  done
  return 1
}

extract_one() {  # $1=file path
  local f="$1" base gene chr tss win start end
  base=$(basename "$f")
  # 未完成（gzip 不完整）→ 跳过
  gzip -t "$f" 2>/dev/null || { echo "[M1_decode] ⏳ $base 下载未完成，跳过"; return; }
  # 基因名 = 第 3 字段（SeqId_Chr_GeneName_ProteinName）
  gene=$(echo "$base" | cut -d'_' -f3)
  coord=$(get_coord "$gene") || { echo "[M1_decode] ⏭ $gene 不在本轮分析基因集，跳过 $base"; return; }
  chr=$(echo "$coord" | awk '{print $1}'); tss=$(echo "$coord" | awk '{print $2}')
  win=$(echo "$coord" | awk '{print $3}')
  start=$((tss - win * 1000)); end=$((tss + win * 1000))
  if [ "$start" -lt 1 ]; then start=1; fi
  echo "[M1_decode] ✂ $base → $gene cis($chr:$start-$end, ±${win}kb)"
  # Chrom 列兼容 "chr1" 与 "1"；保留表头（下游 fread header=TRUE 必需——此前跳表头会导致
  # M3_protein_decode.R 把第一行变异当列名、Pval 列找不到而误判"无工具"）
  zcat "$f" | awk -v C="$chr" -v C2="${chr#chr}" -v S="$start" -v E="$end" '
    (NR==1 && $1=="Chrom") || (($1==C || $1==C2) && $2>=S && $2<=E) {print}
  ' | gzip > "$SUB/${base%.gz}_cis.txt.gz"
  # 校验：表头存在 + 至少 1 行数据（按行数判断，不依赖染色体是否带 chr 前缀）
  hdr=$(gzip -dc "$SUB/${base%.gz}_cis.txt.gz" 2>/dev/null | head -1)
  n=$(gzip -dc "$SUB/${base%.gz}_cis.txt.gz" 2>/dev/null | wc -l)
  if [ "$n" -ge 2 ] && echo "$hdr" | grep -q "^Chrom"; then
    echo "[M1_decode] ✔ $base cis 窗提取 $n 行（含表头）→ $SUB/${base%.gz}_cis.txt.gz"
    if [ "$DEL_ORIG" = "1" ]; then rm -f "$f"; echo "[M1_decode] 🗑 已删原文件 $base"; fi
  else
    echo "[M1_decode] ✘ $base 提取为空/异常（$n 行），保留原文件人工核查"
  fi
}

# 主循环：处理 data/decode/ 下所有 *.txt.gz 顶层文件（不含 sub/）
shopt -s nullglob
for f in "$DECODE"/*.txt.gz; do
  extract_one "$f"
done
echo "[M1_decode] 完成。cis 子集目录: $SUB"
