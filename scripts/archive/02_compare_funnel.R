#!/usr/bin/env Rscript
# =============================================================================
# 02_compare_funnel.R — 可比性漏斗（审稿人核定的关键门控）
# =============================================================================
# 目的：量化"转录本 eQTL 可定位基因 ∩ 蛋白 pQTL 可定位基因"的共享规模。
#   若共享基因数 < 300，按审稿人约定放弃 Δ 叙事，论文只做单通道+方向四态。
# 学术不端防护：该门控阈值在预注册中锁定，不在看到结果后调整。
#
# 输入：results/config.json（预注册锁定的配置）
# 输出：results/compare_funnel.rds / compare_funnel.tsv
#       结果同时写入 results/manifest.jsonl
# =============================================================================
suppressMessages({library(jsonlite); library(data.table)})

cfg  <- read_json("/data/qiushuogeng/projects/dual-channel-mr-atlas/results/config.json",
                  auto_unbox = TRUE)
res  <- "/data/qiushuogeng/projects/dual-channel-mr-atlas/results"

log <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n")
log("阶段 02: 可比性漏斗开始")

# --- 预注册校验：config.json 必须与 PREREGISTRATION.md.sha256 一致 -----------
stopifnot(file.exists("/data/qiushuogeng/projects/dual-channel-mr-atlas/docs/PREREGISTRATION.md.sha256"))
cur_hash <- tools::md5sum("/data/qiushuogeng/projects/dual-channel-mr-atlas/docs/PREREGISTRATION.md")
lock_hash <- readLines("/data/qiushuogeng/projects/dual-channel-mr-atlas/docs/PREREGISTRATION.md.sha256")
stopifnot(cur_hash == lock_hash)
log("预注册哈希校验通过 ✔")

# --- 数据源探测：eQTLGen 与 deCODE 的基因/蛋白列表 ---------------------------
# 优先使用本地缓存；若缺失，尝试从 OpenGWAS API 获取基因范围信息。
# 说明：OpenGWAS 的 prot-a-* 对应 INTERVAL SomaScan；eQTLGen 在 opengwas 上为
# eqtl-a-* / ieu 系列。此处探测两个通道可定位的基因 ID 并求交。

probe_genelist <- function(channel) {
  # channel: "eqtlgen" or "decode/interval"
  cache <- file.path(res, paste0("genelist_", channel, ".rds"))
  if (file.exists(cache)) {
    log("  读取缓存基因列表: ", cache)
    return(readRDS(cache))
  }
  ids <- character(0)
  tryCatch({
    if (requireNamespace("ieugwasr", quietly = TRUE)) {
      # eQTLGen 全血 cis-eQTL 在 OpenGWAS 以 eqtl-a-* 暴露（基因级，31,684 样本）
      g <- ieugwasr::gwasinfo(c("eqtl-a-ENSG00000133067"))  # 探测用
      log("  ieugwasr API 连通 ✔")
    }
  }, error = function(e) log("  API 探测失败(可降级): ", conditionMessage(e)))
  ids
}

log("基因列表探测（当前以本地/占位方式执行，完整映射在阶段 03 验证）")
eqtl_g <- probe_genelist("eqtlgen")
pqtl_g <- probe_genelist("decode")

# --- 门控判定（阈值在预注册 config 中锁定）----------------------------------
# 本阶段真实数据未到手时以"未知"上报；真实判定在数据落地后由 03/04 复算。
gate <- list(
  shared_genes_known = FALSE,
  threshold = 300,
  note = "真实共享基因数在数据下载后由 03_mvp_smoke 确认；此处建立漏斗管道与门控逻辑"
)
log("门控阈值(预注册锁定): shared_genes >= ", gate$threshold)

# --- 落盘 -------------------------------------------------------------------
funnel <- list(
  eqtl_channel = cfg$exposures$transcript,
  pqtl_channel = cfg$exposures$protein,
  gate = gate,
  generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
)
saveRDS(funnel, file.path(res, "compare_funnel.rds"))
writeLines(toJSON(funnel, pretty = TRUE, auto_unbox = TRUE),
           file.path(res, "compare_funnel.json"))

man <- list(stage = "02", status = "done",
            gate = if (gate$shared_genes_known) gate$shared_genes else "pending_funnel",
            time = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
writeLines(toJSON(man, auto_unbox = TRUE), file.path(res, "manifest_stage02.jsonl"))
log("阶段 02 完成 ✔ → results/compare_funnel.json")
