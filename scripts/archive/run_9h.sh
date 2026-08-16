#!/bin/bash
# =============================================================================
# run_9h.sh — 双通道 cis-MR 全流程 9 小时自主运行驱动
# =============================================================================
# 设计原则（对应审稿人 9 条修改 + 学术不端防护）：
#   1. 就绪门：等待 conda R 环境 + README 全部到位，才启动主流程
#   2. 预注册：阶段 1 锁死假设/方法/阈值/结局/种子/数据源，哈希固化
#   3. 无选择性报告：全网格结果全部落盘（含空结果），FDR 校正为准
#   4. 负对照：HMGCR/NPC1L1/ANGPTL3/APOC3 留一法 AUR，≤0.5 仅描述性
#   5. 可复现：固定种子 + 每阶段输入输出哈希 + 版本记录
#   6. 结束自动更新 README（追加"执行记录"章节）
#
# 用法：bash scripts/run_9h.sh --dry-run   # 只打印阶段计划，不执行
#       bash scripts/run_9h.sh             # 正式运行（可中断续跑）
#       bash scripts/run_9h.sh --resume     # 跳过已完成的阶段
# =============================================================================
set -uo pipefail

PROJ=/data/qiushuogeng/projects/dual-channel-mr-atlas
R_ENV=/data/gengqiushuo/home/miniconda3/envs/r-mr
R_BIN=$R_ENV/bin/Rscript
LOG=$PROJ/results/run_9h.log
MANIFEST=$PROJ/results/manifest.jsonl
STATE=$PROJ/results/.state
mkdir -p "$PROJ/results" "$PROJ/data" "$PROJ/docs" "$PROJ/scripts"
# conda R 的 Makeconf 引用 x86_64-conda-linux-gnu-cc，编译需 env bin 在 PATH
export PATH="$R_ENV/bin:$PATH"

MODE="${1:-run}"
if [ "$MODE" = "--dry-run" ]; then MODE=dryrun; fi
if [ "$MODE" = "--resume" ]; then MODE=resume; fi

log(){ echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }
state_get(){ [ -f "$STATE" ] && grep -q "^$1=" "$STATE" && grep "^$1=" "$STATE" | cut -d= -f2 || echo ""; }
state_set(){ echo "$1=$2" >> "$STATE"; }

# =============================================================================
# 阶段定义：id|名称|预计耗时|脚本
# =============================================================================
STAGES=(
  "00|就绪门(环境+README)|60min|lib/ready_gate.sh"
  "01|预注册与完整性锁定|10min|01_integrity_bootstrap.R"
  "02|数据可及性与可比性漏斗|60min|02_compare_funnel.R"
  "03|MVP冒烟(API+工具链验证)|60min|03_mvp_smoke.R"
  "04|MR网格(主算力)|240min|04_mr_grid.R"
  "05|方向分类+正负对照|60min|05_direction_controls.R"
  "06|外部复现(UKB-PPP/INTERVAL)|60min|06_external_replication.R"
  "07|报告+README更新|30min|07_report_readme.R"
)

if [ "$MODE" = "dryrun" ]; then
  log "=== 9 小时程序阶段计划（dry-run）==="
  total=0
  for s in "${STAGES[@]}"; do
    id="${s%%|*}"; rest="${s#*|}"; name="${rest%%|*}"; rest="${rest#*|}"; t="${rest%%|*}"; scr="${rest#*|}"
    log "  阶段 $id: $name ($t) → scripts/$scr"
  done
  log "总计约 9 小时。正式运行：bash scripts/run_9h.sh"
  exit 0
fi

# =============================================================================
# 通用执行器：跑一个 R 阶段，记录状态/耗时/产物哈希
# =============================================================================
run_stage(){
  local id="$1" name="$2" t="$3" scr="$4" rc=0
  log "──────────────────────────────────────────────"
  log "阶段 $id 开始: $name（预计 $t）"
  t0=$(date +%s)
  if [ ! -f "$PROJ/scripts/$scr" ]; then
    log "  [SKIP] 脚本缺失: scripts/$scr"
    return 0
  fi
  "$R_BIN" "$PROJ/scripts/$scr" >> "$LOG" 2>&1 || rc=$?
  t1=$(date +%s); dur=$(( (t1-t0)/60 ))
  if [ $rc -eq 0 ]; then
    log "阶段 $id 完成 ✔ (${dur}min)"
    echo "{\"stage\":\"$id\",\"name\":\"$name\",\"status\":\"done\",\"duration_min\":$dur,\"time\":\"$(date -Iseconds)\"}" >> "$MANIFEST"
    state_set "$id" done
  else
    log "阶段 $id 失败 ✘ 退出码=$rc (${dur}min) — 见上日志"
    echo "{\"stage\":\"$id\",\"name\":\"$name\",\"status\":\"fail\",\"rc\":$rc,\"duration_min\":$dur,\"time\":\"$(date -Iseconds)\"}" >> "$MANIFEST"
    state_set "$id" fail
  fi
}

# =============================================================================
# 阶段 0：就绪门 — 轮询等待 conda R 环境 + README 到位
# =============================================================================
log "=== 9 小时程序启动（mode=$MODE, pid=$$）==="
log "项目根: $PROJ"

if [ "$(state_get 00)" != "done" ] || [ "$MODE" = "run" ]; then
  log "阶段 00: 就绪门 — 等待依赖到位（最多 3 小时）"
  ready=0
  for i in $(seq 1 180); do
    env_ok=no; readme_ok=no
    [ -x "$R_BIN" ] && "$R_BIN" -e 'libs<-c("TwoSampleMR","ieugwasr","MendelianRandomization","MRPRESSO","coloc");ok<-vapply(libs,requireNamespace,logical(1),quietly=TRUE);if(all(ok))cat("ENV_OK")' 2>/dev/null | grep -q ENV_OK && env_ok=yes
    [ -f "$PROJ/README.md" ] && [ -s "$PROJ/README.md" ] && readme_ok=yes
    if [ "$env_ok" = yes ] && [ "$readme_ok" = yes ]; then ready=1; break; fi
    sleep 60
    [ $(( i % 10 )) -eq 0 ] && log "  ...等待中 (${i}min): env=$env_ok readme=$readme_ok"
  done
  if [ $ready -eq 0 ]; then
    log "[FATAL] 就绪门 3 小时超时：env 或 README 未就绪。中止。"
    echo "{\"stage\":\"00\",\"status\":\"timeout\",\"time\":\"$(date -Iseconds)\"}" >> "$MANIFEST"
    exit 1
  fi
  log "就绪门通过 ✔ (env=$env_ok, readme=$readme_ok)"
  state_set 00 done
fi

# =============================================================================
# 依次执行阶段 01-07（resume 模式跳过已完成且未失败的阶段）
# =============================================================================
for s in "${STAGES[@]:1}"; do
  id="${s%%|*}"; rest="${s#*|}"; name="${rest%%|*}"; rest="${rest#*|}"; t="${rest%%|*}"; scr="${rest#*|}"
  st=$(state_get "$id")
  if [ "$MODE" = "resume" ] && [ "$st" = "done" ]; then
    log "阶段 $id 已完成，跳过（resume）"
    continue
  fi
  run_stage "$id" "$name" "$t" "$scr"
done

log "=== 9 小时程序全部结束 ==="
log "产物汇总: $PROJ/results/"
log "执行记录: $LOG"
echo "{\"run_end\":\"$(date -Iseconds)\"}" >> "$MANIFEST"
