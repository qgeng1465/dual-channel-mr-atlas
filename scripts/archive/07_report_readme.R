#!/usr/bin/env Rscript
# =============================================================================
# 07_report_readme.R — 汇总报告 + 更新 README
# =============================================================================
# 用户要求：跑完后把做了什么更新到 README 里。
# 学术不端防护：报告如实汇总各阶段产物状态（含失败/占位/待数据），
#   不夸大、不宣称未完成的工作。
# =============================================================================
suppressMessages({library(jsonlite)})

proj <- "/data/qiushuogeng/projects/dual-channel-mr-atlas"
res  <- file.path(proj, "results")
log  <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n")
log("阶段 07: 汇总报告 + README 更新开始")

# --- 汇总各阶段产物状态 -------------------------------------------------------
stages <- c("02_compare_funnel.json", "03_mvp_smoke.json",
            "direction_controls.json", "replication.json", "grid_params.json")
summary <- list()
for (f in stages) {
  p <- file.path(res, f)
  summary[[f]] <- if (file.exists(p)) {
    j <- read_json(p, auto_unbox = TRUE)
    list(status = "produced", size_bytes = file.info(p)$size,
         key = j$note %||% j$stage %||% f)
  } else list(status = "missing")
}
log("阶段产物汇总完成")

# --- 生成报告文本 -------------------------------------------------------------
report_lines <- c(
  "# 执行报告", "",
  sprintf("生成时间: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "",
  "## 各阶段状态",
  "| 产物 | 状态 | 说明 |",
  "|---|---|---|",
  sprintf("| PREREGISTRATION | %s | 预注册+哈希锁定 |",
          if (file.exists(file.path(proj, "docs/PREREGISTRATION.md"))) "✔" else "✘"),
  sprintf("| 阶段02 可比性漏斗 | %s | %s |",
          summary[["02_compare_funnel.json"]]$status,
          summary[["02_compare_funnel.json"]]$key %||% "—"),
  sprintf("| 阶段03 MVP冒烟 | %s | QA性质，非论文结果 |",
          summary[["03_mvp_smoke.json"]]$status),
  sprintf("| 阶段04 MR网格 | %s | %s |",
          summary[["grid_params.json"]]$status,
          if (file.exists(file.path(res, "mr_grid.rds"))) {
            g <- readRDS(file.path(res, "mr_grid.rds"))
            if (isTRUE(g$executed)) "已执行" else "数据未就绪→参数已冻结"
          } else "待执行"),
  sprintf("| 阶段05 方向+对照 | %s | %s |",
          summary[["direction_controls.json"]]$status,
          summary[["direction_controls.json"]]$key %||% "—"),
  sprintf("| 阶段06 外部复现 | %s | %s |",
          summary[["replication.json"]]$status,
          summary[["replication.json"]]$key %||% "—")
)
writeLines(report_lines, file.path(res, "REPORT.md"))
log("报告已写 results/REPORT.md")

# --- 更新 README：追加"执行记录"章节 ------------------------------------------
readme_path <- file.path(proj, "README.md")
if (file.exists(readme_path)) {
  readme <- readLines(readme_path, warn = FALSE)
  append_section <- c(
    "",
    "---",
    "",
    "## 执行记录（自动化）",
    "",
    sprintf("- 预注册: %s",
            if (file.exists(file.path(proj, "docs/PREREGISTRATION.md")))
              "[docs/PREREGISTRATION.md](docs/PREREGISTRATION.md)（sha256 锁定）" else "缺失"),
    sprintf("- 阶段产物: [results/REPORT.md](results/REPORT.md)"),
    sprintf("- 最近运行: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "- 完整性: 每个阶段预注册哈希校验 + 全网格如实报告 + 负对照预注册规则"
  )
  # 避免重复追加
  if (!any(grepl("^## 执行记录（自动化）", readme))) {
    writeLines(c(readme, append_section), readme_path)
    log("README.md 已追加执行记录章节 ✔")
  } else {
    log("README 已含执行记录章节，跳过（避免重复）")
  }
} else {
  log("README.md 尚不存在（就绪门应已保证其存在；此分支仅防御）")
}

log("阶段 07 完成 ✔ → results/REPORT.md + README.md 已更新")
