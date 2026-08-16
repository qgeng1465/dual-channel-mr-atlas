#!/bin/bash
# ready_gate.sh — 就绪门：校验环境/依赖/预注册是否到位（幂等，可重复调用）
# 由 run_9h.sh 阶段 00 调用。所有检查通过输出 READY，否则输出缺失项并退出 1。
set -uo pipefail

PROJ=/data/qiushuogeng/projects/dual-channel-mr-atlas
R_ENV=/data/gengqiushuo/home/miniconda3/envs/r-mr
R_BIN=$R_ENV/bin/Rscript
export PATH="$R_ENV/bin:$PATH"
fail=0

echo "[ready_gate] 检查项目目录..."
for d in data/eqtlgen data/decode data/opengwas data/ukb_ppp data/interval data/gtex scripts results docs; do
  [ -d "$PROJ/$d" ] || { echo "  [MISS] 目录 $d 不存在"; fail=1; }
done

echo "[ready_gate] 检查 R 环境..."
if [ ! -x "$R_BIN" ]; then
  echo "  [MISS] R 可执行文件不存在: $R_BIN"; fail=1
else
  "$R_BIN" -e 'libs<-c("TwoSampleMR","ieugwasr","MendelianRandomization","MRPRESSO","coloc","susieR","qvalue","data.table","dplyr","ggplot2","jsonlite");ok<-vapply(libs,requireNamespace,logical(1),quietly=TRUE);if(all(ok))cat("PACKAGES_OK") else cat("MISSING:",paste(libs[!ok],collapse=","))' 2>/dev/null | grep -q PACKAGES_OK \
    && echo "  [OK] 关键 R 包全部就绪（MVP 阶段必需）" \
    || { echo "  [MISS] R 包不齐（TwoSampleMR/ieugwasr/MRPRESSO/coloc 等）"; fail=1; }
  # biomaRt 为可延期包（仅 M2 分泌型注释需要；后台安装中，不阻塞 MVP）
fi

echo "[ready_gate] 检查预注册与配置..."
[ -f "$PROJ/docs/PREREGISTRATION.md" ] || { echo "  [MISS] 预注册文档缺失"; fail=1; }
[ -f "$PROJ/docs/PREREGISTRATION.md.sha256" ] || { echo "  [MISS] 预注册哈希缺失"; fail=1; }
[ -f "$PROJ/results/config.json" ] || { echo "  [MISS] config.json 缺失"; fail=1; }

echo "[ready_gate] 检查 README..."
[ -f "$PROJ/README.md" ] && [ -s "$PROJ/README.md" ] \
  && echo "  [OK] README 就绪" \
  || { echo "  [MISS] README.md 不存在或为空"; fail=1; }

if [ $fail -eq 0 ]; then
  echo "[ready_gate] 全部就绪 ✔ READY"
  exit 0
else
  echo "[ready_gate] 存在 $fail 项缺失"
  exit 1
fi
