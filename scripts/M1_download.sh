#!/bin/bash
# =============================================================================
# M1_download.sh — 数据获取（README §4.1/§4.3/§4.4）
# =============================================================================
# 存根状态：URL 已按 README §4 列出，断点续传与 md5 记录逻辑已就绪。
# 运行前人工核验各 URL 当前可用性；deCODE 走 M1_decode_subset.sh（表单周转）。
# 学术不端/可复现要求：每个文件记录 md5 至 docs/CHANGELOG.md。
# =============================================================================
set -uo pipefail
PROJ=/data/qiushuogeng/projects/dual-channel-mr-atlas
CHG=$PROJ/docs/CHANGELOG.md
mkdir -p "$PROJ/data/eqtlgen" "$PROJ/data/gtex" "$PROJ/data/opengwas" "$PROJ/data/ldref" "$PROJ/tools"

dl() {
  local url="$1" out="$2"
  echo "[M1] 下载: $(basename "$out")"
  if ! wget -c -q --tries=5 --timeout=60 "$url" -O "$out"; then
    echo "[M1] ✘ 下载失败: $(basename "$out")" >&2
    return 1
  fi
  if [ ! -s "$out" ]; then
    echo "[M1] ✘ 下载文件为空: $(basename "$out")" >&2
    return 1
  fi
  echo "  → md5: $(md5sum "$out" | cut -d' ' -f1)  ($(basename "$out"))" | tee -a "$CHG"
}

echo "[M1] ===== 通道 1: eQTLGen 全血 cis-eQTL (data/eqtlgen/) ====="
# URL 以 README §4.1 为准（下载前人工核验）
echo "[M1] ⚠ 存根：以下 URL 需人工核验后取消注释执行"
# dl "https://download.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz" "$PROJ/data/eqtlgen/cis-eQTL-significant.txt.gz"
# dl "https://download.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/cis-eQTLs_full_20180905.txt.gz" "$PROJ/data/eqtlgen/cis-eQTL-full.txt.gz"
# dl "https://download.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz" "$PROJ/data/eqtlgen/SNP_AF.txt.gz"

echo "[M1] ===== GTEx v8 组织 eQTL (data/gtex/) ====="
# dl "https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL.tar" "$PROJ/data/gtex/GTEx_Analysis_v8_eQTL.tar"

echo "[M1] ===== 工具: plink 1.9 / SMR 1.3.1 (tools/) ====="
# 存根：文件名以官方实际为准（README §4.4）
# wget -c "https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20220402.zip" -O "$PROJ/tools/plink_linux_x86_64_20220402.zip"
# wget -c "https://yanglab.westlake.edu.cn/software/smr/download/smr-1.3.1-linux-x86_64.zip" -O "$PROJ/tools/smr-1.3.1-linux-x86_64.zip"

echo "[M1] ===== 结局数据 (OpenGWAS, 经 ieugwasr M3 前拉取) ====="
echo "[M1] 结局经 ieugwasr::gwasinfo + extract_outcome_data 拉取（见 README §4.3），不在此下载。"

echo "[M1] 存根完成。deCODE 人工表单后运行: bash scripts/M1_decode_subset.sh decode_links.txt"
