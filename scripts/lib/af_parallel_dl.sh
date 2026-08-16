#!/bin/bash
# =============================================================================
# af_parallel_dl.sh — eQTLGen SNP 频率文件并行分块下载（curl range, 可断点续传）
# 背景：单连接 ~26KB/s（北京→荷兰国际线路），并行分块可达 ~150-200KB/s
# 用法：nohup scripts/lib/af_parallel_dl.sh > /tmp/af_par.log 2>&1 &
# 输出：data/eqtlgen/SNP_AF.txt.gz（与单线程 wget 同名，覆盖旧残档）
# =============================================================================
URL="https://download.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz"
F=/data/qiushuogeng/projects/dual-channel-mr-atlas/data/eqtlgen/SNP_AF.txt.gz
TMP=/tmp/afpar
SIZE=$(curl -sI --max-time 30 "$URL" | awk -F': ' 'tolower($1)=="content-length"{gsub("\r","");print $2}')
[ -z "$SIZE" ] && SIZE=240045342
N=16
CHUNK=$(( (SIZE + N - 1) / N ))
rm -rf "$TMP"; mkdir -p "$TMP"

echo "$(date +%H:%M:%S) 开始并行下载 AF: $SIZE bytes, $N 连接 × $(($CHUNK/1024/1024))MB"
dl_seg() {
  local i=$1 start=$2 end=$3
  local want=$((end - start + 1)) have
  while :; do
    have=$(stat -c%s "$TMP/c_$i" 2>/dev/null || echo 0)
    [ "$have" -ge "$want" ] && break
    timeout 120 curl -sf -r $((start+have))-$end "$URL" >> "$TMP/c_$i" 2>/dev/null
    have=$(stat -c%s "$TMP/c_$i" 2>/dev/null || echo 0)
  done
}
PIDS=()
for i in $(seq 0 $((N-1))); do
  start=$((i * CHUNK)); end=$((start + CHUNK - 1))
  [ "$i" -eq $((N-1)) ] && end=$((SIZE - 1))
  dl_seg "$i" "$start" "$end" &
  PIDS+=($!)
done

# 进度监视（后台每 20s 报告一次聚合字节数）
{
  while :; do
    tot=0
    for f in "$TMP"/c_*; do tot=$((tot + $(stat -c%s "$f" 2>/dev/null || echo 0))); done
    echo "$(date +%H:%M:%S) 进度: $((tot/1024/1024))MB / $((SIZE/1024/1024))MB ($((tot*100/SIZE))%)"
    sleep 20
  done
} > /tmp/af_par_progress.log 2>&1 &
MON=$!

wait "${PIDS[@]}"   # 只等分块段，不等无限监视循环（否则永久阻塞）
kill $MON 2>/dev/null

# 校验并合并（按数字顺序，勿用 c_* 通配——词法序 c_10 会排在 c_2 前！）
tot=0
for f in "$TMP"/c_*; do tot=$((tot + $(stat -c%s "$f" 2>/dev/null || echo 0))); done
echo "$(date +%H:%M:%S) 下载完成: $tot bytes (期望 $SIZE)"
if [ "$tot" -eq "$SIZE" ]; then
  for i in $(seq 0 $((N-1))); do cat "$TMP/c_$i"; done > "$F.new"
  mv "$F.new" "$F"
  gzip -t "$F" && echo "$(date +%H:%M:%S) ✔ AF 文件完整且 gzip 校验通过: $F" || echo "✘ gzip 校验失败！"
else
  echo "✘ 总字节数不匹配，未合并；重启脚本继续。"
fi
