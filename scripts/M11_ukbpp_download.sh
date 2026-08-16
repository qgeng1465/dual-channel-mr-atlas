#!/bin/bash
# =============================================================================
# M11_ukbpp_download.sh — 下载 UKB-PPP（European discovery）8 个目标蛋白 pGWAS tar
# =============================================================================
# 说明：走 mihomo 代理（HK 01，用户惯用）。S3 直连被地域限速(17KB/s)，代理 ~600KB/s。
#   每个 tar ~545MB，8 个共 ~4.4GB，预计 1.5-2.5 小时。
# 用法：bash scripts/M11_ukbpp_download.sh   （建议 nohup 后台）
# 产物：data/ukbpp/{GENE}_*.tar（解压到 data/ukbpp/{GENE}/）
# =============================================================================
set -uo pipefail
PROJ=/data/qiushuogeng/projects/dual-channel-mr-atlas
OUT=$PROJ/data/ukbpp
mkdir -p "$OUT"
export HTTPS_PROXY="http://127.0.0.1:7890" HTTP_PROXY="http://127.0.0.1:7890"
TOKEN=$(awk -F' = ' '/authToken/{print $2}' ~/.synapseConfig)
SYNAPSE="https://repo-prod.prod.sagebase.org"
FENDP="https://file-prod.prod.sagebase.org"

# synid|名字（文件名从 Synapse 实体取）
FILES=(
  "syn51468818"
  "syn52362097"
  "syn51468817"
  "syn51469717"
  "syn52361466"
  "syn51469224"
  "syn52362042"
  "syn51470441"
)

get_psurl() {  # $1=synid -> stdout presigned URL
  local FH=$(curl -s --max-time 20 "$SYNAPSE/repo/v1/entity/$1" \
    -H "Authorization: Bearer $TOKEN" | python3 -c "import json,sys;print(json.load(sys.stdin).get('dataFileHandleId',''))")
  curl -s --max-time 20 "$FENDP/file/v1/fileHandle/batch" -H "Authorization: Bearer $TOKEN" \
    -H "Content-Type: application/json" \
    -d "{\"includeFileHandles\":true,\"includePreSignedURLs\":true,\"requestedFiles\":[{\"fileHandleId\":$FH,\"associateObjectId\":\"$1\",\"associateObjectType\":\"FileEntity\"}]}" \
    | python3 -c "import json,sys;print(json.load(sys.stdin)['requestedFiles'][0]['preSignedURL'])"
}

for SID in "${FILES[@]}"; do
  NAME=$(curl -s --max-time 20 "$SYNAPSE/repo/v1/entity/$SID" -H "Authorization: Bearer $TOKEN" \
         | python3 -c "import json,sys;print(json.load(sys.stdin)['name'])")
  # 避免重下
  if [ -f "$OUT/$NAME" ] && [ -s "$OUT/$NAME" ]; then
    echo "[$(date +%H:%M)] SKIP $NAME (已存在)"; continue
  fi
  URL=$(get_psurl "$SID")
  echo "[$(date +%H:%M)] 下载 $NAME ($SID) ..."
  if curl -sfL --retry 3 --max-time 7200 -x "$HTTPS_PROXY" -o "$OUT/$NAME" "$URL" \
     && tar -tf "$OUT/$NAME" >/dev/null 2>&1; then
    echo "  -> $OUT/$NAME ($(du -h "$OUT/$NAME" | cut -f1))"
  else
    echo "  !! 下载失败/不完整 $NAME (exit $?)"
    rm -f "$OUT/$NAME"   # 不完整则删除，避免下次 -s 误判已下载
  fi
done

# 解压
echo "[$(date +%H:%M)] 解压全部 tar..."
for t in "$OUT"/*.tar; do
  [ -e "$t" ] || continue
  d="${t%.tar}"
  if [ -d "$d" ] || tar -tf "$t" 2>/dev/null | head -1 | grep -q .; then
    mkdir -p "$d"
    tar -xf "$t" -C "$d" && echo "  解压 $t -> $d/ 内容 $(ls "$d" | wc -l) 项"
  fi
done
echo "[$(date +%H:%M)] UKB-PPP 下载+解压完成"
