#!/bin/bash
# =============================================================================
# M19_download_full_gwas.sh — 下载三个结局的全量 GWAS（coloc 全扫描的本地数据前提）
# =============================================================================
# 源：GWAS Catalog FTP（harmonised .h.tsv.gz，hg19/GRCh37；FBG 为原始 MAGIC 文件）
# 输出：data/opengwas/full/<outcome>.gz + 完整性校验 + DONE 标记
# 用法：nohup bash scripts/M19_download_full_gwas.sh > tmp/M19.log 2>&1 &
# 注意：系统盘满，所有输出必须落 /data（本脚本在项目 data/ 下）。
# =============================================================================
set -u
PROJ="/data/qiushuogeng/projects/dual-channel-mr-atlas"
OUT="$PROJ/data/opengwas/full"
mkdir -p "$OUT"
PROXY="http://127.0.0.1:7890"
log(){ echo "[$(date '+%H:%M:%S')] $*"; }

# 已有完成标记则跳过
if [ -f "$OUT/DONE" ]; then log "已下载完成，跳过"; exit 0; fi

dl(){ # dl <name> <url>
  local name="$1" url="$2"
  local dst="$OUT/$name"
  if [ -f "$dst" ] && gzip -t "$dst" 2>/dev/null; then log "$name 已存在且完整"; return 0; fi
  log "下载 $name <- $url"
  # -C - 断点续传；-f 令 HTTP 错误直接失败；重试 5 次
  # 注意：校验失败的 .part 必须删除从头重下（对损坏/错误页续传会永久污染文件）
  local try=0 ok=0
  while [ $try -lt 5 ]; do
    curl -sSfL --retry 3 -C - --max-time 1800 -x "$PROXY" -o "$dst.part" "$url"
    if [ -f "$dst.part" ] && gzip -t "$dst.part" 2>/dev/null; then ok=1; break; fi
    rm -f "$dst.part"
    try=$((try+1)); log "  retry $try（不完整/失败，从头重下）"; sleep 3
  done
  if [ $ok -eq 1 ]; then mv "$dst.part" "$dst"; log "$name OK: $(stat -c%s "$dst") bytes"; else
    log "$name 下载失败"; rm -f "$dst.part"; return 1
  fi
}

T2D="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006867/harmonised/30054458-GCST006867-EFO_0001360.h.tsv.gz"
CAD="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST005001-GCST006000/GCST005194/harmonised/29212778-GCST005194-EFO_0000378.h.tsv.gz"
FBG="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST005001-GCST006000/GCST005186/MAGIC_Manning_et_al_FastingGlucose_MainEffect.txt.gz"

rc=0
dl t2d_full.gz "$T2D" || rc=1
dl cad_full.gz "$CAD" || rc=1
dl fbg_full.gz "$FBG" || rc=1

if [ $rc -eq 0 ]; then
  md5sum "$OUT"/*.gz > "$OUT/MD5SUM.txt"
  echo "$(date '+%Y-%m-%d %H:%M:%S') ALL_DONE" > "$OUT/DONE"
  log "全部下载完成"
else
  log "部分下载失败（rc=$rc），不写 DONE"
  exit $rc
fi
