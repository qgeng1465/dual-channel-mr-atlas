#!/usr/bin/env Rscript
# =============================================================================
# M8_tissue_triangulation.R — 组织三角验证：106 个转录 strong 共定位命中 × GTEx 6 组织
# =============================================================================
# 目的：全血 eQTL 共定位命中的因果信号，是否在疾病相关组织（肝/胰/冠脉等）中
#   同样驱动该基因表达，且该组织 eQTL 的 MR 方向与全血一致 → "组织三角验证"。
#   这把描述性 atlas 升级为"跨组织强化"的 atlas，是冲 eBioMedicine 需要的验证层。
# 口径：
#   - 命中源：results/grid/transcript_coloc.csv tier=="strong"（106 个 gene×outcome）
#   - 组织源：GTEx v8 egenes（data/gtex/*.egenes.txt.gz），lead cis-eQTL = 该基因
#     tss±1Mb 内 p_nominal 最小（预注册工具阈值 p<5e-6 才纳入 → 与主通道同调）
#   - MR：单工具 Wald（与 M6 同口径）；结局 API 按 outcome 批量提取（每结局一次）
#   - 方向一致性：sign(组织 MR b) == sign(全血 coloc stage2_b)
#   - 强化判定：组织 MR p<0.05 且方向与全血一致（"reinforced"）；仅方向一致但
#     p≥0.05 记 "consistent(ns)"（弱信号，诚实不夸大）；无组织工具记 "no_instrument"
#   - 探索性补充，不移动四态主分析门柱（同 M6 设计纪律）
# 用法：PATH=$R_ENV/bin:$PATH Rscript scripts/M8_tissue_triangulation.R
# 产物：results/grid/tissue_triangulation.csv（逐 命中×组织）+ summary 打印
# =============================================================================
suppressMessages({
  library(data.table)
  suppressPackageStartupMessages(library(TwoSampleMR))
})
proj <- "/data/qiushuogeng/projects/dual-channel-mr-atlas"
res  <- file.path(proj, "results")
GTEX <- file.path(proj, "data/gtex")
dir.create(file.path(res, "grid"), recursive = TRUE, showWarnings = FALSE)
log <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n")
OUT_FULL <- c(t2d = "ebi-a-GCST006867", cad = "ebi-a-GCST005194", fbg = "ebi-a-GCST005186")
TYPE <- c(t2d = "cc", cad = "cc", fbg = "quant")
S_CASE <- c(t2d = 61714 / 655666, cad = 34541 / 296525)

TISSUES <- c("Liver","Pancreas","Whole_Blood","Adipose_Subcutaneous","Muscle_Skeletal","Artery_Coronary")

# --- 命中源 ------------------------------------------------------------------
hits <- fread(file.path(res, "grid", "transcript_coloc.csv"))[tier == "strong" & ok == TRUE]
hits <- hits[!is.na(symbol) & symbol != ""]
log("转录 strong 命中:", nrow(hits), "（", length(unique(hits$symbol)), " 基因）")

# --- 读组织 egenes，仅留命中基因 -------------------------------------------------
load_tissue <- function(t) {
  f <- file.path(GTEX, paste0(t, ".egenes.txt.gz"))
  if (!file.exists(f)) return(NULL)
  d <- fread(cmd = paste0("zcat ", shQuote(f)), sep = "\t", header = TRUE, nThread = 4)
  if (nrow(d) == 0) return(NULL)
  setnames(d, c("gene_name", "rs_id_dbSNP151_GRCh38p7"), c("symbol", "rsid"), skip_absent = TRUE)
  d[, rsid := fifelse(is.na(rsid) | rsid %in% c("", ".", "-"), variant_id, rsid)]
  d[, Tissue := t]
  d
}
tis <- lapply(TISSUES, load_tissue)
names(tis) <- TISSUES
tis <- tis[!vapply(tis, is.null, logical(1))]
log("组织就绪:", paste(names(tis), collapse = ", "))

# --- 每组织每命中基因取 lead cis-eQTL（p<5e-6）-----------------------------------
lead_all <- rbindlist(lapply(names(tis), function(t) {
  d <- tis[[t]][symbol %in% unique(hits$symbol) & abs(tss_distance) <= 1e6 & pval_nominal < 5e-6]
  if (nrow(d) == 0) return(data.table())
  d <- d[order(pval_nominal), .SD[1], by = symbol]
  d[, .(tissue = t, symbol, rsid, variant_id, ref, alt, ref_factor, maf,
        p_nominal = pval_nominal, slope, slope_se)]
}))
cat("\n=== 命中基因 × 组织 lead cis-eQTL（p<5e-6）可用性 ===\n")
print(dcast(lead_all[, .(n = .N), by = .(tissue, symbol)], symbol ~ tissue,
            value.var = "n", fill = 0)[, total := Liver + Pancreas + Whole_Blood +
              Adipose_Subcutaneous + Muscle_Skeletal + Artery_Coronary][order(-total)])
log("lead eQTL 总行:", nrow(lead_all))

# --- outcome GWAS 区域：离线复用 coloc 缓存（_coloc_gwas/{GCST}_{gene}.rds）-----
#   （M5 同款缓存：每 结局×基因 cis±1Mb 区域，含 rsid/ea/nea/beta/se/p/n；106/106 命中全覆盖）
OUT_N <- c(t2d = "GCST006867", cad = "GCST005194", fbg = "GCST005186")
gwas_cache <- list()
for (i in seq_len(nrow(hits))) {
  h <- hits[i]
  f <- file.path(res, "grid", "_coloc_gwas", paste0(OUT_N[[h$outcome]], "_", h$gene, ".rds"))
  if (!file.exists(f)) next
  g <- readRDS(f)
  if (!is.data.table(g)) g <- as.data.table(g)
  g <- g[!is.na(rsid) & !is.na(beta)][, .(SNP = rsid, beta.outcome = beta, se.outcome = se,
                                          pval.outcome = p, effect_allele.outcome = ea,
                                          other_allele.outcome = nea)]
  g <- g[!duplicated(SNP)]
  gwas_cache[[paste0(h$gene, "_", h$outcome)]] <- g
}
log("离线 outcome 区域缓存（基因×结局）:", length(gwas_cache), " 个")

# --- 逐 命中×组织 跑 Wald MR ---------------------------------------------------
res_rows <- list()
for (i in seq_len(nrow(hits))) {
  h <- hits[i]
  sym <- h$symbol; o <- h$outcome
  eo <- gwas_cache[[paste0(h$gene, "_", o)]]   # 该命中自己的 结局×基因 区域缓存
  for (t in unique(lead_all$tissue)) {
    ld <- lead_all[symbol == sym & tissue == t]
    if (nrow(ld) == 0) next
    row <- data.table(gene = h$gene, symbol = sym, outcome = o,
                      tissue = t, lead_rsid = ld$rsid,
                      lead_p_eqtl = ld$p_nominal, lead_slope = ld$slope,
                      wb_b = h$stage2_b, wb_p = h$stage2_pval, wb_top_snp = h$top_snp)
    ex <- if (is.null(eo)) data.table() else eo[SNP == ld$rsid]
    if (nrow(ex) == 0 || any(is.na(ex$beta.outcome))) {
      row[, mr_b := NA]; row[, mr_p := NA]; row[, status := "no_outcome_match"]
      res_rows[[length(res_rows) + 1]] <- row; next
    }
    # 等位对齐：GTEx slope 以 alt 为效应；outcome 以 effect_allele 为效应
    A <- ld$alt; O <- ld$ref; ea <- ex$effect_allele.outcome; nea <- ex$other_allele.outcome
    flip <- (A == nea & O == ea)
    ok   <- (A == ea)
    pal  <- (A == "A" & O == "T") | (A == "T" & O == "A") | (A == "C" & O == "G") | (A == "G" & O == "C")
    keep <- (ok | flip) & !pal
    if (!keep) { row[, mr_b := NA]; row[, mr_p := NA]; row[, status := "no_alignment"]; res_rows[[length(res_rows) + 1]] <- row; next }
    # 仅翻 outcome 侧（flip ⇒ outcome 效应等位=GTEx ref，需翻成 alt 方向）；eQTL slope 保持 alt 方向
    b_t <- ld$slope; se_t <- ld$slope_se
    b_o <- if (flip) -ex$beta.outcome else ex$beta.outcome
    se_o <- ex$se.outcome
    # Wald ratio（单工具）
    b_wald <- b_o / b_t
    se_wald <- se_o / abs(b_t)
    p_wald <- 2 * pnorm(-abs(b_wald / se_wald))
    row[, mr_b := b_wald]; row[, mr_se := se_wald]; row[, mr_p := p_wald]
    row[, dir_consistent := sign(mr_b) == sign(wb_b)]
    row[, status := fifelse(!is.na(mr_p) & mr_p < 0.05 & dir_consistent == TRUE, "reinforced",
                     fifelse(dir_consistent == TRUE, "consistent_ns", "discordant"))]
    res_rows[[length(res_rows) + 1]] <- row
  }
  if (i %% 10 == 0) log("  命中进度 ", i, "/", nrow(hits))
}
out_tri <- rbindlist(res_rows, fill = TRUE)
fwrite(out_tri, file.path(res, "grid", "tissue_triangulation.csv"))
log("已写出 results/grid/tissue_triangulation.csv（行数 ", nrow(out_tri), "）")

# --- 摘要 ---------------------------------------------------------------------
cat("\n=== 组织三角验证摘要 ===\n")
cat("命中×组织 对数（有组织 lead eQTL）:", nrow(out_tri), "\n")
cat("\n状态分布:\n"); print(table(out_tri$status, useNA = "ifany"))
cat("\n按组织 reinforced 数:\n"); print(out_tri[status == "reinforced", .N, by = tissue][order(-N)])
cat("\n按结局 reinforced 数:\n"); print(out_tri[status == "reinforced", .N, by = outcome])
reinf <- out_tri[status == "reinforced"]
cat("\n强化命中（组织 MR p<0.05 且方向与全血一致）:", nrow(reinf), " 对\n")
if (nrow(reinf)) print(reinf[, .(symbol, outcome, tissue, wb_b, mr_b, mr_p, lead_rsid)][order(symbol, tissue)])
# 每命中是否有≥1 组织强化
per_hit <- out_tri[, .(n_tissue_any = .N,
                   n_tissue_lead = sum(!is.na(mr_b)),
                   n_reinforced = sum(status == "reinforced", na.rm = TRUE),
                   n_consistent_ns = sum(status == "consistent_ns", na.rm = TRUE)),
               by = .(gene, symbol, outcome)]
cat("\n命中被≥1 组织强化:", sum(per_hit$n_reinforced > 0), "/", nrow(per_hit),
    "（", sprintf("%.0f%%", 100 * sum(per_hit$n_reinforced > 0) / nrow(per_hit)), "）\n")
cat("命中被≥1 组织方向一致（含 ns）:", sum(per_hit$n_consistent_ns + per_hit$n_reinforced > 0), "/", nrow(per_hit), "\n")
fwrite(per_hit, file.path(res, "grid", "tissue_triangulation_hits.csv"))
log("=== M8 组织三角验证完成 ✔ ===")
