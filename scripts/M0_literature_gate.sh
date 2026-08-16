#!/bin/bash
# =============================================================================
# M0_literature_gate.sh — 文献门禁（README §6-M0）
# =============================================================================
# 目的：防止审稿意见 1/2/9 所禁表述与错误引用进入项目。
# 实现要点（导师修订后要求）：
#   - 禁词正则用【真词】『蛋白主导→单抗』（[a] 写法永远匹配不到真实禁词，禁用）
#   - 行号白名单豁免：§0 对照表、§0 红线声明、§6-M0 本注释块、§3/§8 否定性/引用性口径
#   - 命中禁词且不在豁免范围 → 非 0 退出，阻断后续步骤
#   - 逐条核对 §9 引用核验表字段非空
# 存根状态：基础实现已就绪，运行前人工复核豁免行号范围与正则。
# =============================================================================
set -uo pipefail
PROJ=/data/qiushuogeng/projects/dual-channel-mr-atlas
PATTERN='同腔室|隔离翻译|蛋白主导→单抗|成药模态|大概率阳性'

# --- 豁免（§0 对照表/红线、§6-M0 注释块、§3/§8 否定性引用、M0 存根说明行）----
# 说明：这些区间的禁词属"引用性/否定性"口径，不阻断。
# 健壮性（2026-08-06 修订）：行号硬编码会随 README 编辑漂移（§4.2 扩写后 M0 块
#   已从 257 行漂到 260 行）。除行号区间外，增加【内容豁免】：任何含
#   "M0_literature_gate.sh"（M0 存根说明行，必然引用禁词示例）或"本注释块"的行
#   一律豁免。他人若在正文写禁词作为己方结论（不含这些标记），仍会被命中阻断。
# 当前豁免区间：
#   §0 对照表       : 13-30 （## 0. 至 ## 1. 前）
#   §6-M0 注释块    : 239-265（### M0 至 ### M1 前，上限放宽容忍漂移）
exempt() {
  local ln="$1"; local content="$2"
  if [ "$ln" -ge 13 ] && [ "$ln" -le 30 ]; then return 0; fi    # §0
  if [ "$ln" -ge 239 ] && [ "$ln" -le 265 ]; then return 0; fi  # M0 块
  case "$content" in
    *M0_literature_gate.sh*) return 0;;  # M0 存根说明行（必引禁词示例）
    *本注释块*)              return 0;;
  esac
  return 1
}
grep_scan() {
  local f="$1"
  grep -nE "$PATTERN" "$f" 2>/dev/null | while IFS=: read -r ln rest; do
    if exempt "$ln" "$rest"; then
      continue
    fi
    echo "  [违规] ${f}:${ln}: $rest"
  done
}

echo "[M0] 文献门禁 — 禁词扫描（真实禁词正则）"
echo "[M0] 目标文件: README.md + docs/"
violations=0
for f in "$PROJ/README.md" "$PROJ"/docs/*.md; do
  [ -f "$f" ] || continue
  hits=$(grep_scan "$f")
  if [ -n "$hits" ]; then
    echo "$hits"
    violations=$((violations+1))
  fi
done
echo "[M0] 禁词扫描完成，违规文件数=$violations"

# --- 引用核验表字段非空检查 ---------------------------------------------------
echo "[M0] 引用核验表检查 (docs/citation_checklist.md)"
CL="$PROJ/docs/citation_checklist.md"
if [ ! -f "$CL" ]; then
  echo "[M0] 警告: citation_checklist.md 不存在（M0 放行，投稿前须补）"
else
  empty=$(grep -cE '\|\s*\|' "$CL" 2>/dev/null || true)
  echo "[M0] citation_checklist 空字段行数=$empty（投稿前须为 0）"
fi

if [ "$violations" -gt 0 ]; then
  echo "[M0] 门禁未通过：存在未豁免禁词。阻断后续步骤。"
  exit 1
else
  echo "[M0] 门禁通过 ✔（禁词扫描全绿；引用表待人工补全）"
  exit 0
fi
