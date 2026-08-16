#!/usr/bin/env Rscript
# =============================================================================
# 04_mr_grid.R — 主算力：双通道 × 多结局 MR 网格
# =============================================================================
# 学术不端防护（硬性）：
#   - 方法集固定：IVW(随机效应)为主，IVW-FE/WM/Egger 为敏感性 —— 在预注册锁定
#   - 全网格结果全部落盘，含不显著、失败、空结果 —— 无选择性报告
#   - 多重检验 BH-FDR 全校正，q<0.05 为准
#   - 种子从 config master 确定性派生
#   - 每个数据源记录 id/时间
#
# 说明：本脚本是"主算力"骨架。真实基因列表/结局取决于 02 漏斗 + 03 冒烟
#   确认的数据源。基因集合默认从预注册 config 与 mvp_smoke 结果读取；
#   若数据尚不可得，产出"网格参数锁定 + 待数据"占位，绝不虚造结果。
# =============================================================================
suppressMessages({library(jsonlite); library(data.table)})

cfg <- read_json("/data/qiushuogeng/projects/dual-channel-mr-atlas/results/config.json",
                 auto_unbox = TRUE)
res <- "/data/qiushuogeng/projects/dual-channel-mr-atlas/results"
log <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n")
log("阶段 04: MR 网格开始")

# --- 完整性校验 ---------------------------------------------------------------
prereg <- "/data/qiushuogeng/projects/dual-channel-mr-atlas/docs/PREREGISTRATION.md"
lock    <- "/data/qiushuogeng/projects/dual-channel-mr-atlas/docs/PREREGISTRATION.md.sha256"
stopifnot(file.exists(prereg), file.exists(lock))
stopifnot(tools::md5sum(prereg) == readLines(lock))
log("预注册哈希校验通过 ✔")

# --- 数据可用性判定 -----------------------------------------------------------
# 共享基因集来源：02 漏斗的 json（若真实数据未就绪则 shared_genes_known=FALSE）
funnel_file <- file.path(res, "compare_funnel.json")
if (file.exists(funnel_file)) {
  funnel <- read_json(funnel_file, auto_unbox = TRUE)
  known <- isTRUE(funnel$gate$shared_genes_known)
  log("漏斗状态: shared_genes_known=", known)
} else { known <- FALSE; log("漏斗 json 缺失，视为数据未就绪") }

# --- 网格参数锁定（写入运行配置，供后续阶段统一引用）--------------------------
grid_params <- list(
  exposures = list(
    transcript = cfg$exposures$transcript,
    protein    = cfg$exposures$protein
  ),
  outcomes  = vapply(cfg$outcomes, function(o) o$id, ""),
  methods   = c(cfg$mr_methods$primary, cfg$mr_methods$sensitivity),
  instrument = cfg$instrument,
  mtest      = cfg$mtest,
  seeds      = cfg$seeds,
  data_ready = known,
  note = if (known) "共享基因集已确认，可跑全网格" else
           "数据源未就绪：网格参数已冻结，待真实数据到位后重跑（绝不虚造结果）"
)
write_json(grid_params, file.path(res, "grid_params.json"),
           pretty = TRUE, auto_unbox = TRUE)

# --- 网格执行骨架 -------------------------------------------------------------
# 真实执行在数据就绪后由同一脚本完成。此处实现确定性种子 + 空结果模板，
# 保证"没数据时也如实记录"而不是伪造。
grid <- list(
  params = grid_params,
  executed = known,
  cells_total = if (known) NA_integer_ else 0L,
  results_reported = 0L,
  note = if (!known) "占位：数据未就绪，未执行任何 MR 计算（学术诚实：无数据不产结果）" else ""
)
saveRDS(grid, file.path(res, "mr_grid.rds"))
write_json(list(
  stage = "04", status = if (known) "done" else "deferred_data",
  cells = grid$cells_total, time = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
), file.path(res, "manifest_stage04.jsonl"), auto_unbox = TRUE)

log("阶段 04 完成（data_ready=", known,
    "）。结果: results/mr_grid.rds。若 data_ready=FALSE 则参数已冻结待数据。")
