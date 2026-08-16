#!/usr/bin/env Rscript
# =============================================================================
# report_decode_results.R — deCODE 蛋白通道完成后自动更新 README/CHANGELOG + 独立摘要
# =============================================================================
# 由 run_decode_pipeline.sh 调用（下载完成后自动执行）。
# 输入：results/grid/protein_decode_mr.csv + results/funnel/funnel_protein_decode.tsv
# 输出：docs/CHANGELOG.md 顶部插入新条目
#       results/DECODE_PIPELINE_SUMMARY.md（独立可读结果摘要，供用户直接查看）
#       README（v0.9 多管线叙事）不写回运行记录，保持纯叙事（2026-08-07 P2）
# =============================================================================
suppressMessages({library(data.table)})
proj <- "/data/qiushuogeng/projects/dual-channel-mr-atlas"
res  <- file.path(proj, "results")
README <- file.path(proj, "README.md")
CHANGELOG <- file.path(proj, "docs/CHANGELOG.md")
SUMMARY <- file.path(res, "DECODE_PIPELINE_SUMMARY.md")
log <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n")

stopifnot(file.exists(README), file.exists(CHANGELOG))
g <- file.path(res, "grid/protein_decode_mr.csv")
f <- file.path(res, "funnel/funnel_protein_decode.tsv")
if (!file.exists(g)) { log("✘ 缺少 protein_decode_mr.csv，跳过报告"); quit(save = "no", status = 1) }

pm <- fread(g)
# 2026-08-07 P2：funnel 文件缺失/损坏时不能崩掉整个报告（fread 缺文件即抛错）
fn <- NULL
if (file.exists(f)) fn <- tryCatch(fread(f), error = function(e) NULL)
outmap <- c(`ebi-a-GCST006867` = "T2D", `ebi-a-GCST005194` = "CAD", `ebi-a-GCST005186` = "FBG")
pm[, out := fifelse(outcome %in% names(outmap), unname(outmap[outcome]), outcome)]

# 主方法行（IVW-MRE / Wald）用于结果表；失败行带 note
main <- pm[ok == TRUE & method %in% c("Inverse variance weighted (multiplicative random effects)",
                                      "Wald ratio")]
n_ok   <- sum(pm$ok == TRUE)
n_pair <- nrow(pm)
n_conc <- if (!is.null(fn) && "mr_completed_rows" %in% fn$stage) fn[stage == "mr_completed_rows", count] else NA_integer_
n_prot <- if (!is.null(fn) && "tested_proteins" %in% fn$stage) fn[stage == "tested_proteins", count] else NA_integer_

fmt_p <- function(p) ifelse(is.na(p), "—", formatC(p, format = "g", digits = 2))
fmt_b <- function(x) ifelse(is.na(x), "—", formatC(x, format = "f", digits = 3))

# --- 结果表 markdown -----------------------------------------------------------
rows <- ""
if (nrow(main) > 0) {
  setorder(main, pval, na.last = TRUE)
  for (i in seq_len(nrow(main))) {
    rows <- paste0(rows, "| ", main$gene[i], " | ", main$out[i], " | ", main$nsnp[i],
                   " | ", fmt_b(main$b[i]), " | ", fmt_b(main$se[i]), " | ",
                   fmt_p(main$pval[i]), " | ", main$method[i], " |\n")
  }
} else {
  rows <- "| （全部 18 对均未产出主方法结果，如实记录） | | | | | | |\n"
}

# 失败/空结果说明（去重，简写）
fail <- unique(pm[ok == FALSE | is.na(method), .(gene, out, note = substr(note, 1, 40))])
fail_txt <- ""
if (nrow(fail) > 0) {
  fl <- paste(sprintf("  - %s×%s: %s", fail$gene, fail$out, fail$note), collapse = "\n")
  fail_txt <- paste0("**未产出结果的对（含空结果，如实报告）**\n", fl, "\n")
}

stamp <- format(Sys.time(), "%Y-%m-%d %H:%M")

# --- README（v0.9 多管线叙事）不做运行记录插入 --------------------------------
# 2026-08-07 P2：README v0.9 已删"- 最近运行"行与执行记录章节——README 只承载纯叙事、
# 不承载运行历史（用户"多管线叙事"要求；真实改动轨迹留在 CHANGELOG，是诚实审计轨迹）。
# 运行详情一律落 docs/CHANGELOG.md（顶部条目）+ results/DECODE_PIPELINE_SUMMARY.md。
# 若未来需刷新 README §3 通道 B 的 deCODE 状态行（⏳→✅），在此加 marker 处理。
log("README 保持纯叙事：运行记录写 CHANGELOG + SUMMARY，不插入 README ✔")

# --- CHANGELOG 顶部插入 --------------------------------------------------------
cl <- readLines(CHANGELOG)
cl <- c(cl[1],
        "",
        paste0("## ", format(Sys.time(), "%Y-%m-%d"), " — deCODE 蛋白通道真实运行 + 代码审查修复"),
        "",
        paste0("- deCODE 4 个对照蛋白文件下载完成（PCSK9/HMGCR/ANGPTL3/APOC3，各 ~953MB，16 连接并行），cis 窗自动提取，蛋白通道 cis-MR 真实运行。结果见 `results/grid/protein_decode_mr.csv`、`results/DECODE_PIPELINE_SUMMARY.md`。"),
        paste0("- 代码审查修复（技术 bug fix）：① `M1_decode_subset.sh` cis 子集此前丢弃表头 → 下游 `fread(header=TRUE)` 把第一行变异当列名、Pval 列丢失、误判“无工具”（已改为保留表头，校验改按行数）；② `M3_protein_decode.R` 排除 deCODE readme 注明的多等位基因 bug 行（effectAllele==otherAllele / otherAllele=\"!\"），rsids 为 NA 时回落 Name。"),
        paste0("- 核验为“非 bug”项（诚实记录，防误改）：eQTLGen eaf=1−AlleleB_all（AA/AB/BB 计数公式验证一致）；Z→β=Z/√N、se=1/√N 标准转化；LD clump r²<0.01@1000kb EUR 逻辑正确。"),
        "",
        cl[2:length(cl)])
writeLines(cl, CHANGELOG)
log("CHANGELOG 已插入条目 ✔")

# --- 独立摘要文件 --------------------------------------------------------------
sum_txt <- paste0(
"# deCODE 蛋白通道 MR 结果摘要（", stamp, "）\n\n",
"- 蛋白源：deCODE 血浆 pQTL（Ferkingstad 2021, n=35,559）；结局：T2D/CAD/FBG（OpenGWAS）。\n",
"- 方法：cis ±1Mb（hg38, ENSEMBL 坐标）、p<5e-6、EUR LD clump r²<0.01@1000kb、主 mr_ivw_mre / nsnp=1 Wald。\n\n",
"## 结果表（主方法，按 p 排序）\n\n",
"| 蛋白 | 结局 | nsnp | b | SE | p | 方法 |\n|---|---|---|---|---|---|---|\n",
rows, "\n",
fail_txt, "\n",
"## 漏斗\n\n",
"- 测试蛋白：", n_prot, " 个\n",
"- 蛋白×结局对：", n_pair, " 对\n",
"- 产出 MR 结果的对：", n_conc, " 对\n\n",
"## 产物文件\n\n",
"- `results/grid/protein_decode_mr.csv` — 全网格（含空/失败）\n",
"- `results/figures/F9_funnel_protein.png` — 蛋白通道漏斗\n",
"- `results/figures/F10_protein_mr_forest.png` — 蛋白 MR 森林图\n",
"- 转录通道图：`F1`（漏斗）、`F2`（v1→v2 坍缩）、`F3`（极端 p 消失）\n\n",
"## 说明（诚实报告）\n\n",
"- EAF：assocvariants.annotated 未就绪 → palindromic 保守排除 + ImpMAF 近似（ImpMAF 不总是 effect allele 频率）。\n",
"- 多等位基因 bug 行（effectAllele==otherAllele / \"!\"）已排除。\n",
"- 转录通道“结果一般”主因：全血表达对肝/肠脂质基因是弱代理；cis-eQTL 严格 clump 降功效；非编码问题。详见 CHANGELOG。\n"
)
writeLines(sum_txt, SUMMARY)
log("独立摘要已写: ", SUMMARY)
