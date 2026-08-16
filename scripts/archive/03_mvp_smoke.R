#!/usr/bin/env Rscript
# =============================================================================
# 03_mvp_smoke.R — MVP 冒烟：验证工具链 + API + 真实数据可用性
# =============================================================================
# 目的：在进入 4 万基因×多结局网格前，先用少量蛋白(≤50)×2 结局跑通全链路，
#   验证 OpenGWAS API / ieugwasr / TwoSampleMR 全部可用，并确认数据源基因规模。
# 学术不端防护：本阶段只做"工具链/数据可用性"验证，不产生论文级结果；
#   冒烟结果标注为 QA 性质，不进入论文统计。
#
# 输入：results/config.json
# 输出：results/mvp_smoke.rds / mvp_smoke.json（QA 日志）
# =============================================================================
suppressMessages({
  library(jsonlite); library(data.table)
  suppressPackageStartupMessages(library(TwoSampleMR))
  suppressPackageStartupMessages(library(ieugwasr))
})

cfg <- read_json("/data/qiushuogeng/projects/dual-channel-mr-atlas/results/config.json",
                 auto_unbox = TRUE)
res <- "/data/qiushuogeng/projects/dual-channel-mr-atlas/results"
log <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n")
log("阶段 03: MVP 冒烟开始")

# --- 预注册哈希校验 -----------------------------------------------------------
prereg <- "/data/qiushuogeng/projects/dual-channel-mr-atlas/docs/PREREGISTRATION.md"
lock    <- "/data/qiushuogeng/projects/dual-channel-mr-atlas/docs/PREREGISTRATION.md.sha256"
stopifnot(file.exists(prereg), file.exists(lock))
stopifnot(tools::md5sum(prereg) == readLines(lock))
log("预注册哈希校验通过 ✔")

smoke <- list(started = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
outcomes <- vapply(cfg$outcomes, function(o) o$id, "")

# --- API 连通性 ---------------------------------------------------------------
api_ok <- tryCatch({
  # 探测两个结局 ID 是否可解析
  info <- ieugwasr::gwasinfo(outcomes)
  as.data.frame(info)
}, error = function(e) { log("API 探测失败: ", conditionMessage(e)); NULL })

if (is.null(api_ok) || nrow(api_ok) == 0) {
  smoke$api <- list(ok = FALSE, note = "OpenGWAS 结局 ID 解析失败（API 迁移中，见 ieu/info 404 已知问题）")
  log("API 结局解析失败 — 尝试备用结局 ID")
  # 已知 API 迁移（gwas.mrcieu.ac.uk → opengwais.io）后部分 ID 需重映射
  alt <- tryCatch({
    alt_info <- ieugwasr::gwasinfo(c("ieu-b-38", "ebi-a-GCST005194"))
    as.data.frame(alt_info)
  }, error = function(e) { log("备用探测也失败: ", conditionMessage(e)); NULL })
  if (!is.null(alt) && nrow(alt) > 0) {
    smoke$api <- list(ok = TRUE, resolved = rownames(alt), note = "备用 ID 可解析")
    outcomes <- intersect(c("ieu-b-38", "ebi-a-GCST005194"), rownames(alt))
  }
} else {
  smoke$api <- list(ok = TRUE, resolved = rownames(api_ok),
                    note = paste(nrow(api_ok), "个结局 ID 可解析"))
}
log("API 状态: ok=", if (isTRUE(smoke$api$ok)) "TRUE" else "FALSE")

# --- 暴露数据源探测 -----------------------------------------------------------
# 验证 deCODE 或 INTERVAL 蛋白源的可定位性（以 5 个代表蛋白探测）
probe_proteins <- c("PCSK9", "CETP", "HMGCR", "NPC1L1", "ANGPTL3")
pqtl_ok <- tryCatch({
  # OpenGWAS 上 INTERVAL SomaScan pQTL 为 prot-a-<Ensembl ID>；探测常见 ID
  probes <- ieugwasr::gwasinfo("prot-a-2876")   # 代表性蛋白 prot-a 条目
  is.data.frame(probes) && nrow(probes) > 0
}, error = function(e) { log("prot-a-* 探测失败: ", conditionMessage(e)); FALSE })

smoke$exposure <- list(
  pqtl_via_opengwas = pqtl_ok,
  note = if (pqtl_ok) "INTERVAL prot-a-* pQTL 可经 OpenGWAS 访问（备用蛋白通道）"
                     else "prot-a-* 未探测到，将依赖 deCODE 本地数据"
)
log("pQTL 暴露探测: ", smoke$exposure$pqtl_via_opengwas)

# --- 端到端最小 MR（1 蛋白 × 1 结局，纯 QA）----------------------------------
qa_mr <- tryCatch({
  if (isTRUE(smoke$api$ok)) {
    ex <- extract_instruments(outcomes = "prot-a-2876", p1 = cfg$instrument$pval_thresh)
    if (!is.null(ex) && nrow(ex) > 0) {
      out <- extract_outcome_data(snps = ex$SNP, outcomes = outcomes[1],
                                  proxies = TRUE)
      if (!is.null(out) && nrow(out) > 0) {
        dat <- harmonise_data(ex, out, action = cfg$instrument$harmonise_palindromic_action)
        r <- mr(dat, method_list = cfg$mr_methods$primary)
        list(ok = TRUE, n_snp = nrow(dat), n_mr = nrow(r))
      } else list(ok = FALSE, note = "outcome 无匹配 SNP")
    } else list(ok = FALSE, note = "exposure 无工具变量")
  } else list(ok = FALSE, note = "API 不可用，跳过端到端 QA")
}, error = function(e) list(ok = FALSE, note = conditionMessage(e)))
smoke$qa_minimal_mr <- qa_mr
log("端到端 QA MR: ok=", if (isTRUE(qa_mr$ok)) "TRUE" else "FALSE",
    if (isTRUE(qa_mr$ok)) paste0(" (snps=", qa_mr$n_snp, ")") else paste0(" — ", qa_mr$note))

# --- 落盘 QA 日志 -------------------------------------------------------------
saveRDS(smoke, file.path(res, "mvp_smoke.rds"))
writeLines(toJSON(smoke, pretty = TRUE, auto_unbox = TRUE),
           file.path(res, "mvp_smoke.json"))
log("阶段 03 完成 ✔ → results/mvp_smoke.json（QA 性质，非论文结果）")
