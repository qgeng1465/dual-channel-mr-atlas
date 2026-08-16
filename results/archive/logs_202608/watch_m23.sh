#!/bin/bash
# M23 完成监控：任一结局 CSV 出现即记录时间，全部出现后输出摘要
RES=/data/qiushuogeng/projects/dual-channel-mr-atlas/results
LOG=/data/qiushuogeng/projects/dual-channel-mr-atlas/logs
while true; do
  done_count=0
  for on in t2d cad fbg; do
    if [ -f "$RES/coloc_full_${on}_20260815.csv" ]; then done_count=$((done_count+1)); fi
  done
  if [ "$done_count" -eq 3 ]; then
    echo "$(date +%H:%M:%S) ALL-3-DONE" >> "$LOG/watch_m23.log"
    for on in t2d cad fbg; do
      tail -6 "$LOG/m23_${on}.log" >> "$LOG/watch_m23.log"
    done
    exit 0
  fi
  # 检查进程是否全部死亡（异常退出）
  if ! pgrep -f "M23_full_scan.R" >/dev/null 2>&1; then
    echo "$(date +%H:%M:%S) NO-M23-PROC-ALIVE (done=$done_count)" >> "$LOG/watch_m23.log"
    for on in t2d cad fbg; do
      echo "--- $on ---" >> "$LOG/watch_m23.log"
      tail -4 "$LOG/m23_${on}.log" >> "$LOG/watch_m23.log"
    done
    exit 1
  fi
  sleep 30
done
