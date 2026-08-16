#!/bin/bash
# =============================================================================
# gtex_parallel_dl.sh — GTEx v8 eQTL 并行分块下载（Zenodo GV-Rep 镜像，支持 Range/206）
# =============================================================================
# 背景：GTEx 官方 tar 实测 404；导师 agent 核实 OpenGWAS 无 GTEx 组织 eQTL；
#       Zenodo GV-Rep 复传 GTEx_Analysis_v8_eQTL.zip（1.49GB）Range 206 可用。
# 镜像：https://zenodo.org/records/11502840  (DOI 10.5281/zenodo.11502840)
# 用法：nohup bash scripts/lib/gtex_parallel_dl.sh > /tmp/gtex_par.log 2>&1 &
# 输出：data/gtex/GTEx_Analysis_v8_eQTL.zip（16 连接并行，约 5-6h）
# 校验：python zipfile.testzip()（CRC 全检）；Zenodo 真支持 Range，不会混入坏字节。
# =============================================================================
URL="https://zenodo.org/records/11502840/files/GTEx_Analysis_v8_eQTL.zip?download=1"
D=/data/qiushuogeng/projects/dual-channel-mr-atlas/data/gtex
F=$D/GTEx_Analysis_v8_eQTL.zip
TMP=/tmp/gtexpar
SIZE=1563245056
N=16
mkdir -p "$D"
CHUNK=$(( (SIZE + N - 1) / N ))

if [ -s "$F" ] && python3 -c "import zipfile,sys; sys.exit(0 if zipfile.ZipFile('$F').testzip() is None else 1)" 2>/dev/null; then
  echo "$(date +%H:%M:%S) OK $F 已存在且 CRC 通过，跳过下载"
  exit 0
fi
resume=0
if [ -d "$TMP" ]; then
  for i in $(seq 0 $((N-1))); do
    w=$(( i == N-1 ? SIZE - (N-1)*CHUNK : CHUNK ))
    s="$TMP/c_$i"
    [ -f "$s" ] && [ "$(stat -c%s "$s")" -lt "$w" ] && resume=1
  done
fi
if [ "$resume" = "0" ]; then rm -rf "$TMP"; fi
mkdir -p "$TMP"
echo "$(date +%H:%M:%S) 开始 GTEx zip 并行下载: $SIZE bytes, $N 连接 x $(($CHUNK/1024/1024))MB (resume=$resume)"

dl_seg() {
  local i=$1 start=$2 end=$3
  local want=$((end - start + 1)) have tries=0
  while :; do
    have=$(stat -c%s "$TMP/c_$i" 2>/dev/null || echo 0)
    [ "$have" -ge "$want" ] && break
    tries=$((tries+1)); [ $tries -gt 300 ] && { echo "  seg $i 超时放弃"; break; }
    timeout 120 curl -sf -r $((start+have))-$end "$URL" >> "$TMP/c_$i" 2>/dev/null
    have=$(stat -c%s "$TMP/c_$i" 2>/dev/null || echo 0)
    [ "$have" -gt "$want" ] && truncate -s "$want" "$TMP/c_$i" 2>/dev/null
  done
}
PIDS=()
for i in $(seq 0 $((N-1))); do
  start=$((i * CHUNK)); end=$((start + CHUNK - 1))
  [ "$i" -eq $((N-1)) ] && end=$((SIZE - 1))
  dl_seg "$i" "$start" "$end" &
  PIDS+=($!)
done
{
  while :; do
    tot=0; for f in "$TMP"/c_*; do tot=$((tot + $(stat -c%s "$f" 2>/dev/null || echo 0))); done
    echo "$(date +%H:%M:%S) 进度: $((tot/1024/1024))MB / $((SIZE/1024/1024))MB ($((tot*100/SIZE))%)"
    sleep 120
  done
} > /tmp/gtex_par_progress.log 2>&1 &
MON=$!
wait "${PIDS[@]}"
kill $MON 2>/dev/null

for i in $(seq 0 $((N-1))); do
  w=$(( i == N-1 ? SIZE - (N-1)*CHUNK : CHUNK ))
  s="$TMP/c_$i"
  [ "$(stat -c%s "$s" 2>/dev/null || echo 0)" -gt "$w" ] && truncate -s "$w" "$s" 2>/dev/null
done
tot=0; for s in "$TMP"/c_*; do tot=$((tot + $(stat -c%s "$s" 2>/dev/null || echo 0))); done
echo "$(date +%H:%M:%S) 聚合: $tot / $SIZE"
if [ "$tot" -eq "$SIZE" ]; then
  for i in $(seq 0 $((N-1))); do cat "$TMP/c_$i"; done > "$F.new"
  mv "$F.new" "$F"
  if python3 -c "import zipfile,sys; sys.exit(0 if zipfile.ZipFile('$F').testzip() is None else 1)" 2>/dev/null; then
    echo "$(date +%H:%M:%S) OK GTEx zip 完整且 CRC 通过: $F"
  else
    echo "FAIL $F CRC 失败，重启续传"
  fi
else
  echo "FAIL 字节不匹配($tot/$SIZE)，保留分块；重启脚本续传"
fi
