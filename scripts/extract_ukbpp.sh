#!/bin/bash
set -e
cd /data/qiushuogeng/projects/dual-channel-mr-atlas/data/ukbpp
shopt -s nullglob   # 无 *.tar 时不进入循环（避免 tar -xf "*.tar" 报错退出）
for t in *.tar; do
  echo "[$(date +%H:%M:%S)] 提取 $t"
  tar -xf "$t"
done
echo "[$(date +%H:%M:%S)] 全部提取完成"
