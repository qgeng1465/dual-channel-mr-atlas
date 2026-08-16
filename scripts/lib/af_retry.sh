#!/bin/bash
# eQTLGen SNP 频率文件循环续传（240MB，服务器带宽受限）
URL="https://download.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz"
F="/data/qiushuogeng/projects/dual-channel-mr-atlas/data/eqtlgen/SNP_AF.txt.gz"
for i in $(seq 1 200); do
  gzip -t "$F" 2>/dev/null && { echo "AF 完整 ✔ (轮$i, $(stat -c%s "$F") bytes)"; exit 0; }
  wget -c -q --timeout=100 --tries=2 -O "$F" "$URL" 2>/dev/null
  echo "$(date +%H:%M:%S) 轮$i: $(stat -c%s "$F" 2>/dev/null || echo 0) bytes"
  sleep 3
done
echo "200 轮未完成，退出"
