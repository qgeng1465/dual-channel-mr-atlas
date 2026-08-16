#!/usr/bin/env Rscript
# =============================================================================
# M7_finngen_replication.R — 结局侧外部复现（FinnGen 独立队列，方向一致率）
# =============================================================================
# 目的：对转录 strong 共定位命中（top SNP）+ 蛋白通道 MR ok 对，在同一 top 变异的
#   FinnGen 结局（finn-b-E4_DM2 T2D / finn-b-I9_CHD CAD）下提取关联，
#   与原始结局（ebi-a-GCST006867 / GCST005194）比方向一致率 + β 相关。
#   这是预注册 §9 "外部复现" 的结局侧承诺（跨结局，非跨 eQTL/pQTL）。
# 口径：等位基因对齐（效应等位匹配/翻转）；只比方向不比幅度；FBG 无 FinnGen 对应
#   端点 → 如实跳过并注明。探索性补充，不移动主门柱。
# 用法：PATH=$R_ENV/bin:$PATH Rscript scripts/M7_finngen_replication.R
# 产物：results/grid/finngen_replication.csv
# =============================================================================
suppressMessages({library(data.table); suppressPackageStartupMessages(library(TwoSampleMR))})
proj <- "/data/qiushuogeng/projects/dual-channel-mr-atlas"
gdir <- file.path(proj, "results", "grid")
log <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n", sep = "")
OUT_FULL <- c(t2d = "ebi-a-GCST006867", cad = "ebi-a-GCST005194", fbg = "ebi-a-GCST005186")
FINN <- c("ebi-a-GCST006867" = "finn-b-E4_DM2", "ebi-a-GCST005194" = "finn-b-I9_CHD")
log("M7 FinnGen 结局侧复现")

c <- fread(file.path(gdir, "transcript_coloc.csv"))
strong <- c[tier == "strong" & !is.na(top_snp) & top_snp != ""]
p <- fread(file.path(gdir, "protein_coloc.csv"))
p_ok <- p[ok == TRUE & !is.na(top_snp) & top_snp != ""]
log("转录 strong:", nrow(strong), " | 蛋白 ok:", nrow(p_ok))

tr <- strong[, .(symbol, gene, outcome = OUT_FULL[outcome], top_snp, channel = "transcript")]
id2short <- c("ebi-a-GCST006867" = "t2d", "ebi-a-GCST005194" = "cad", "ebi-a-GCST005186" = "fbg")
pr <- p_ok[, .(symbol = gene, gene, outcome, top_snp, channel = "protein")]
all <- rbindlist(list(tr, pr), use.names = TRUE, fill = TRUE)
all <- all[!is.na(outcome)]

out_rows <- list()
for (o in unique(all$outcome)) {
  finn <- FINN[o]   # 单括号：缺失名返回 NA（[[ 会下标越界）
  if (is.na(finn)) { log("结局", o, "无 FinnGen 对应端点，跳过（FBG）"); next }
  sub <- all[outcome == o]
  snps <- unique(sub$top_snp)
  log("结局", o, "→", finn, ": ", length(snps), " top 变异")
  o1 <- tryCatch(extract_outcome_data(snps, o), error = function(e) NULL)
  o2 <- tryCatch(extract_outcome_data(snps, finn), error = function(e) NULL)
  if (is.null(o1) || is.null(o2)) { log("  提取失败"); next }
  g1 <- as.data.table(o1)[, .(SNP, b1 = as.numeric(beta.outcome), ea1 = as.character(effect_allele.outcome),
                              oa1 = as.character(other_allele.outcome), p1 = as.numeric(pval.outcome))]
  g2 <- as.data.table(o2)[, .(SNP, b2 = as.numeric(beta.outcome), ea2 = as.character(effect_allele.outcome),
                              oa2 = as.character(other_allele.outcome), p2 = as.numeric(pval.outcome))]
  m <- merge(g1, g2, by = "SNP")
  m <- m[!is.na(b1) & !is.na(b2)]
  if (nrow(m) == 0) { log("  无共同变异（FinnGen 未覆盖）"); next }
  # 等位基因对齐：ea2==ea1 同向；ea2==oa1 翻转（均须双向等位完全一致，防多等位位点仅单等位命中的错配）
  same <- m$ea2 == m$ea1 & m$oa2 == m$oa1
  flip <- m$ea2 == m$oa1 & m$oa2 == m$ea1
  m[, b2a := ifelse(flip, -b2, b2)]
  # 2026-08-13 修正：原 `%in%` 是整列向量判断（ea2 命中表内任意其他变异的等位即算可对齐），
  # 把 rs17716350 的 T>A 次等位编码错判为 aligned。改为本行逐元素 `==`，仅 ea2==ea1 或 ea2==oa1 才可对齐。
  # 2026-08-16 补强：aligned 须 same|flip 双等位匹配——原逻辑只检查 ea2∈{ea1,oa1}，
  #   若 ea2==oa1 但 oa2!=ea1（多等位位点），flip=FALSE 会按未翻转 b2 计入，方向比对错。
  ambig <- !(same | flip)  # 无法对齐的记入但不计入一致率（生效行，2026-08-13 核查确认）
  m[, aligned := !ambig]
  cc <- m[aligned == TRUE]
  concord <- if (nrow(cc)) sum(sign(cc$b1) == sign(cc$b2a)) / nrow(cc) else NA_real_
  rho <- if (nrow(cc) >= 3) suppressWarnings(cor(cc$b1, cc$b2a, method = "spearman")) else NA_real_
  sig_o1 <- if (nrow(cc)) sum(cc$p1 < 0.05) else 0
  sig_both <- if (nrow(cc)) sum(cc$p1 < 0.05 & cc$p2 < 0.05) else 0
  log(sprintf("  n=%d 可对齐 %d | 方向一致 %.1f%% | Spearman ρ=%.3f | 原结局显著 %d / 双显著 %d",
              nrow(m), nrow(cc), 100 * concord, rho, sig_o1, sig_both))
  m[, `:=`(outcome = o, finngen_outcome = finn)]
  out_rows[[length(out_rows) + 1]] <- m
}

if (length(out_rows)) {
  res <- rbindlist(out_rows, fill = TRUE)
  fwrite(res, file.path(gdir, "finngen_replication.csv"))
  log("已写出 results/grid/finngen_replication.csv（逐变异方向对比，n=", nrow(res), "）")
} else {
  log("无任何可复现对")
}
