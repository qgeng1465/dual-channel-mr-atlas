#!/bin/bash
# =============================================================================
# run_feasibility.sh — 破局方案最小可行性验证管道（14:00 启动，落 /data）
# =============================================================================
# 链：M19 下载 3 全量 GWAS → M20 试点 coloc(400 对) → M21 UKB-PPP 覆盖 → 写摘要
# 用法：setsid nohup bash scripts/run_feasibility.sh > tmp/feasibility.log 2>&1 &
# 输出：results/feasibility_20260815.md（人读摘要）
# =============================================================================
set -u
PROJ="/data/qiushuogeng/projects/dual-channel-mr-atlas"
RES="$PROJ/results"
RENV="/data/gengqiushuo/home/miniconda3/envs/r-mr/bin"
log(){ echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
mkdir -p "$PROJ/tmp"   # 重定向 M20/M21 日志前必须先建目录，否则 `> tmp/M20.log` 直接失败

# 等到 14:00（若已过则立即开始）
target=$(date -d 'today 14:00' +%s 2>/dev/null || date -d '14:00' +%s)
now=$(date +%s)
if [ "$now" -lt "$target" ]; then
  secs=$((target-now)); log "等待至 14:00（还有 $((secs/60)) 分钟）..."
  sleep "$secs"
fi
log "=== 开始可行性验证管道 ==="

# 1. 下载（约 30-70 分钟，网络绑定；若已被提前启动，等待其结束）
while pgrep -f 'M19_download_full_gwas.sh' >/dev/null; do log "等待已启动的 M19 下载完成..."; sleep 60; done
bash "$PROJ/scripts/M19_download_full_gwas.sh"
if [ ! -f "$RES/../data/opengwas/full/DONE" ]; then log "下载未完成，中止管道"; exit 1; fi
log "M19 下载完成"

# 2. 试点 coloc（CPU 绑定）
# 注意：R/fread 临时文件走 tempdir()，认 TMPDIR（不认 R_TMPDIR）→ 必须指向 /data
#（根盘 50G 满，/tmp/Rtmp* 曾致 fread 管道缓冲 No space left on device）
export TMPDIR=/data/qiushuogeng/tmp/rtmp
PATH="$RENV/bin:$PATH" "$RENV/Rscript" "$PROJ/scripts/M20_feasibility_pilot.R" \
  > "$PROJ/tmp/M20.log" 2>&1
log "M20 试点 coloc 结束（rc=$?）→ tmp/M20.log"

# 3. UKB-PPP 覆盖
python3 "$PROJ/scripts/M21_ukbpp_coverage.py" > "$PROJ/tmp/M21.log" 2>&1
log "M21 UKB-PPP 覆盖结束 → tmp/M21.log"

# 4. 写摘要
{
  echo "# 破局方案最小可行性验证结果（2026-08-15）"
  echo ""
  echo "> 验证 1-2 小时，先验后投。管道：M19 下载 → M20 试点 coloc → M21 UKB-PPP 覆盖。"
  echo ""
  echo "## 1. 全量 GWAS 下载（M19）"
  echo '```'
  ls -la "$PROJ/data/opengwas/full/" 2>/dev/null | grep -E 'gz|DONE|MD5'
  echo '```'
  echo ""
  echo "## 2. 试点 coloc：MR 显著集之外的 PP.H4≥0.8 分布 + 分层（M20）"
  echo "> 依据 2026-08-15 对抗性评审（verification_20260815.md）："
  echo "> 关键判据不是裸的 PP.H4≥0.8 比例，而是把 coloc-only 命中按 GWAS 区域峰 p 分层："
  echo ">   - GWAS 峰 p<5e-8 且 eQTL 弱 → MR 功率问题（召回叙事成立）"
  echo ">   - GWAS 峰 p≥5e-8 → coloc 伪阳候选（召回叙事不成立，转向）"
  echo '```'
  cat "$PROJ/tmp/M20.log" 2>/dev/null | grep -A 40 '=== 试点 coloc 摘要'
  echo '```'
  echo ""
  echo "## 3. UKB-PPP 覆盖（M21，方案 B 路线判定）"
  echo '```'
  cat "$PROJ/tmp/M21.log" 2>/dev/null
  echo '```'
  echo ""
  echo "## 4. 判定与下一步"
  echo "- 下载可行 → 全量 coloc 本地前提成立（是/否，见 §1）"
  echo "- 召回叙事成立 → 全量 49,866 对 coloc 值得跑（是/否，见 §2）"
  echo "- UKB-PPP 蛋白层路线（见 §3）"
  echo "- 完整评估 agent 评审：见 \`verification_20260815.md\`（若已产出）"
} > "$RES/feasibility_20260815.md"
log "摘要已写 results/feasibility_20260815.md"
log "=== 管道结束 ==="
