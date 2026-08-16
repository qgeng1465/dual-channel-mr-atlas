#!/bin/bash
# run_gtex_pipeline.sh — GTEx 组织 eQTL 管线（下载完成 → 解压 → 组织 MR → README）
# 用法：nohup bash scripts/lib/run_gtex_pipeline.sh > /tmp/gtex_pipeline.log 2>&1 &
set -uo pipefail
PROJ=/data/qiushuogeng/projects/dual-channel-mr-atlas
R_ENV=/data/gengqiushuo/home/miniconda3/envs/r-mr
R_BIN=$R_ENV/bin/Rscript
ZIP=$PROJ/data/gtex/GTEx_Analysis_v8_eQTL.zip
log(){ echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

log "=== GTEx 管线启动：等待 zip 下载完成（最长 12h）==="
for i in $(seq 1 144); do
  if [ -s "$ZIP" ] && python3 -c "import zipfile,sys; sys.exit(0 if zipfile.ZipFile('$ZIP').testzip() is None else 1)" 2>/dev/null; then
    break
  fi
  # 下载器死亡则重启
  if ! pgrep -f "gtex_parallel_dl[.]sh" >/dev/null 2>&1; then
    log "  ⚠ 下载器不在运行，重启"
    nohup bash "$PROJ/scripts/lib/gtex_parallel_dl.sh" >> /tmp/gtex_par.log 2>&1 &
    sleep 5
  fi
  sleep 300
done
if ! python3 -c "import zipfile,sys; sys.exit(0 if zipfile.ZipFile('$ZIP').testzip() is None else 1)" 2>/dev/null; then
  log "✘ zip 未在窗口内完成，退出"
  exit 1
fi
log "✔ GTEx zip 完整"

log "跑 GTEx 组织 eQTL MR（M6_gtex_mr.R，含解压 + 肝/胰/全血组织）..."
export PATH="$R_ENV/bin:$PATH"
if "$R_BIN" "$PROJ/scripts/M6_gtex_mr.R" >> /tmp/gtex_mr.log 2>&1; then
  log "✔ GTEx 组织 MR 完成 → results/grid/gtex_mr.csv"
else
  log "✘ GTEx 组织 MR 失败（exit=$?），见 /tmp/gtex_mr.log"
fi
log "=== GTEx 管线结束 ==="
