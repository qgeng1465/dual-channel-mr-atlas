#!/bin/bash
# 环境搭建：用户自有 conda-forge R 4.4 + MR 管线包
# 目标位置：/data/gengqiushuo/home/miniconda3（/data 盘，不占系统盘）
# 用法：bash scripts/env_setup.sh  （可重复执行，幂等）
set -euo pipefail

MINI=/data/gengqiushuo/home/miniconda3
ENV_NAME=r-mr
log(){ echo "[$(date '+%H:%M:%S')] $*"; }

# 1. 安装 miniconda（若不存在）
if [ ! -x "$MINI/bin/conda" ]; then
  log "安装 miniconda 到 $MINI ..."
  curl -fsSL -o /data/gengqiushuo/home/miniconda_install.sh \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  bash /data/gengqiushuo/home/miniconda_install.sh -b -p "$MINI"
  rm -f /data/gengqiushuo/home/miniconda_install.sh
fi
export PATH="$MINI/bin:$PATH"
conda config --add channels conda-forge >/dev/null 2>&1 || true
conda config --set channel_priority strict >/dev/null 2>&1 || true

# 接受 Anaconda 默认通道 ToS（新版 conda 要求；幂等）
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main 2>/dev/null || true
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r    2>/dev/null || true
log "Anaconda ToS 已接受（pkgs/main, pkgs/r）"

# 2. 创建 R 环境
if ! conda env list | grep -q "$ENV_NAME"; then
  log "创建 conda R 环境 $ENV_NAME（R 4.4，conda-forge 二进制，免编译）..."
  conda create -y -n "$ENV_NAME" -c conda-forge \
    r-base=4.4 r-data.table r-dplyr r-tidyr r-ggplot2 r-ggraph r-httr r-jsonlite r-readr \
    r-coloc r-susier r-metafor r-curl r-openssl r-remotes r-foreach r-doParallel r-matrix
  # 注意：r-qvalue 不在 conda-forge（README §11：qvalue 走 Bioconductor），
  #       故不在 conda create 列表，由下方 BiocManager::install 安装。
fi
R="$MINI/envs/$ENV_NAME/bin/Rscript"
log "R 版本: $("$R" -e 'cat(R.version.string)')"

# 3. 安装 MR 管线包（GitHub + CRAN + Bioconductor）
# 安装路径（网络实测 2026-08-06，本机环境限定）：
#   TwoSampleMR / ieugwasr → GitHub MRCIEU/*（master 分支可用；mrcieu r-universe 的
#     CDN r2.ropensci.org 被本网络 403 拦截，CRAN 无 TwoSampleMR，故走 GitHub 源码安装）
#   MR-PRESSO → GitHub rondolab/MR-PRESSO
#   biomaRt / qvalue → Bioconductor（qvalue 不在 conda-forge/CRAN）
"$R" -e '
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", quiet=TRUE)
if (!requireNamespace("remotes", quietly=TRUE)) install.packages("remotes", quiet=TRUE)

need <- c("MendelianRandomization","MRPRESSO","TwoSampleMR","ieugwasr","biomaRt","qvalue")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly=TRUE)]
if (length(miss)) {
  cat("待安装缺失包:", miss, "\n")
  if (any(c("TwoSampleMR","ieugwasr") %in% miss)) {
    # GitHub 源码安装（master 分支；若已存在包则跳过）
    if ("TwoSampleMR" %in% miss)
      remotes::install_github("MRCIEU/TwoSampleMR", ref="master",
                              upgrade="never", quiet=TRUE, dependencies=TRUE)
    if ("ieugwasr" %in% miss)
      remotes::install_github("MRCIEU/ieugwasr", ref="master",
                              upgrade="never", quiet=TRUE, dependencies=TRUE)
  }
  cr <- setdiff(miss, c("TwoSampleMR","ieugwasr","biomaRt","qvalue"))
  if ("MRPRESSO" %in% cr) {
    remotes::install_github("rondolab/MR-PRESSO", upgrade="never", quiet=TRUE)
    cr <- setdiff(cr, "MRPRESSO")
  }
  if (length(cr)) install.packages(cr, quiet=TRUE)
  bioc <- intersect(c("biomaRt","qvalue"), miss)  # qvalue 走 Bioconductor（README §11）
  if (length(bioc)) BiocManager::install(bioc, update=FALSE, ask=FALSE, quiet=TRUE)
}

final <- c("TwoSampleMR","ieugwasr","MendelianRandomization","MRPRESSO","coloc","susieR",
           "biomaRt","data.table","dplyr","ggplot2","ggraph","qvalue","metafor","jsonlite","httr")
ok <- vapply(final, requireNamespace, logical(1), quietly=TRUE)
cat("\n=== 包安装结果 ===\n")
print(data.frame(pkg=final, installed=ok, row.names=NULL))
if (!all(ok)) quit(status=1)
'
log "环境搭建完成 ✅"
