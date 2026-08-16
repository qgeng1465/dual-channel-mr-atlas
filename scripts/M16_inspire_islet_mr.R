#!/usr/bin/env Rscript
# =============================================================================
# M16_inspire_islet_mr.R — P5：InsPIRE 胰岛 eQTL × 结局 单工具 MR（正交组织通道）
# =============================================================================
# 目的：用第二个**组织特异性** eQTL 源（InsPIRE 人胰岛 eQTL，~420 donors，Varshney 2019）
# 复核 76 优先基因的基因-结局方向，重点回应"KCNJ11 不在全血表达"红线：
#   血中失明的位点，胰岛通道是否"活过来"。
# 数据：data/gtex/islet/InsPIRE_islets_independent_gene_eQTLs.txt
#   （每基因独立 lead eQTL：SNPid=rsID, A1/A2, MAF, Slope, Nominal_Pval）
# 逻辑：对每个有胰岛 eQTL 的 76 优先基因，取 lead 变异 → OpenGWAS 取结局效应 →
#   harmonise 对齐 → 单工具 Wald MR → 与 eQTLGen 通道方向对比。
# 方向一致 = sign(GWAS β at islet lead) == sign(islet Slope × eQTLGen MR b)。
# 纪律：单工具敏感性；探索性；不构成预注册变更。
# 用法：PATH=$R_ENV/bin:$PATH Rscript scripts/M16_inspire_islet_mr.R
# 输出：results/grid/inspire_islet_mr.csv + 终端摘要
# =============================================================================
suppressMessages({
  library(data.table)
  suppressPackageStartupMessages(library(TwoSampleMR))
  library(jsonlite)
})
proj <- "/data/qiushuogeng/projects/dual-channel-mr-atlas"
res  <- file.path(proj, "results")
log <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n")
prereg <- file.path(proj, "docs/PREREGISTRATION.md")
lock   <- file.path(proj, "docs/PREREGISTRATION.md.sha256")
stopifnot(file.exists(prereg), file.exists(lock), tools::md5sum(prereg) == readLines(lock))

OUTCOMES <- list(t2d = "ebi-a-GCST006867", cad = "ebi-a-GCST005194", fbg = "ebi-a-GCST005186")

insp <- fread(file.path(proj, "data/gtex/islet/InsPIRE_islets_independent_gene_eQTLs.txt"), sep = "\t")
# 2026-08-16 路径修复：原相对路径 tmp/ 已归档至 results/archive/tmp_202608/（且 tmp/ 目录已删除）
ppf <- file.path(proj, "results/archive/tmp_202608/priority_genes_heidi_pass.txt")
if (!file.exists(ppf)) ppf <- file.path(proj, "tmp/priority_genes_heidi_pass.txt")
genes76 <- fread(ppf, sep = "\t", header = FALSE)
setnames(genes76, c("ENSG", "symbol", "outcome", "PPH4", "pHEIDI"))
# 与 eQTLGen MR 方向（stage2_b）对照
s2 <- fread(file.path(res, "grid/transcript_grid_stage2.csv"))
coloc <- fread(file.path(res, "grid/transcript_coloc_hits.csv"))[tier == "strong"]

insp76 <- insp[GeneName %in% unique(genes76$symbol)]
# 该文件无 Slope_se → 由两尾 p 反推 z 再算 se（方向核查足够；标准近似）
# 2026-08-16 修复（与 M6_islet 同款）：p<~2.2e-16 时 1-p 舍入为 1 → qchisq=Inf → z=Inf → se=0
#   → Wald b/se=Inf p=0（假性过显著）。|z| 封顶 qnorm(1-1e-16)=8.209536 防止 se=0。
insp76[, z := sign(Slope) * pmin(sqrt(qchisq(1 - Nominal_Pval, 1)), 8.209536)]
insp76[, Slope_se := abs(Slope) / abs(z)]
log("InsPIRE 有胰岛 eQTL 的优先基因: ", uniqueN(insp76$GeneName))

run_one <- function(sym, on) {
  d <- insp76[GeneName == sym][order(Nominal_Pval)][1]
  if (nrow(d) == 0) return(NULL)
  exp_dat <- data.frame(SNP = d$SNPid,
    effect_allele.exposure = d$A1, other_allele.exposure = d$A2,
    eaf.exposure = d$MAF,
    beta.exposure = d$Slope, se.exposure = d$Slope_se,
    pval.exposure = d$Nominal_Pval,
    id.exposure = sym, exposure = sym, stringsAsFactors = FALSE)
  out <- tryCatch(extract_outcome_data(snps = d$SNPid, outcomes = OUTCOMES[[on]], proxies = TRUE),
                  error = function(e) NULL)
  if (is.null(out) || nrow(out) == 0)
    return(data.table(symbol = sym, outcome = on, islet_lead = d$SNPid,
                      islet_slope = d$Slope, islet_p = d$Nominal_Pval,
                      mr_b = NA, mr_p = NA, concordant = NA, note = "无结局匹配"))
  dat <- harmonise_data(exp_dat, out, action = 2); dat <- dat[dat$mr_keep, ]
  if (nrow(dat) == 0)
    return(data.table(symbol = sym, outcome = on, islet_lead = d$SNPid,
                      islet_slope = d$Slope, islet_p = d$Nominal_Pval,
                      mr_b = NA, mr_p = NA, concordant = NA, note = "harmonise 无保留"))
  rw <- tryCatch(mr(dat, method_list = "mr_wald_ratio"), error = function(e) NULL)
  if (is.null(rw) || nrow(rw) == 0)
    return(data.table(symbol = sym, outcome = on, islet_lead = d$SNPid,
                      islet_slope = d$Slope, islet_p = d$Nominal_Pval,
                      mr_b = NA, mr_p = NA, concordant = NA, note = "MR 无输出"))
  # 两通道方向对比：胰岛通道 MR b（rw$b = β_gwas/β_islet_eqtl）vs 全血通道 eQTLGen MR b
  eg <- s2[gene %in% genes76[symbol == sym, ENSG] & outcome == on, b]
  eg_b <- if (length(eg) && !is.na(eg[1])) eg[1] else NA
  conc <- if (!is.na(eg_b) && rw$b != 0) sign(rw$b) == sign(eg_b) else NA
  # 2026-08-13 修（第二轮核查）：`gwas_b`/`gwas_p` 原误命名——实际存的是 Wald 比率 MR 效应
  #（= β_GWAS/β_islet-eQTL，与 `wald_b` 同值），并非原始 GWAS beta → 改名 `mr_b`/`mr_p`。
  data.table(symbol = sym, outcome = on, islet_lead = d$SNPid,
             islet_slope = d$Slope, islet_p = d$Nominal_Pval,
             mr_b = rw$b, mr_p = rw$pval,
             concordant = conc, eqtlgen_b = eg_b, note = "")
}

out <- list()
syms <- unique(insp76$GeneName)
for (s in syms) {
  oss <- unique(genes76[symbol == s, outcome])
  for (on in oss) {
    r <- run_one(s, on)
    if (!is.null(r)) { out[[length(out) + 1]] <- r; log("  ", s, "×", on,
        if (is.na(r$mr_b)) "无结局" else paste0("wald b=", round(r$mr_b,3), " p=", format(r$mr_p, digits=2))) }
  }
}
outDT <- rbindlist(out, fill = TRUE)   # 2026-08-13 修：`res` 已被 rbindlist 覆盖为表 → 改 outDT，保留 res 目录路径
fwrite(outDT, file.path(res, "grid/inspire_islet_mr.csv"))
valid <- outDT[!is.na(mr_b)]
conc  <- outDT[concordant == TRUE]
cat("\n=== InsPIRE 胰岛通道摘要 ===\n")
cat("优先基因 × 结局对: ", nrow(outDT), "（胰岛 eQTL 覆盖）\n")
cat("结局匹配有效: ", nrow(valid), "\n")
if (nrow(valid)) cat("方向一致率 (vs eQTLGen): ", nrow(conc), "/", nrow(valid), " = ", nrow(conc)/nrow(valid), "\n")
k <- outDT[symbol == "KCNJ11"]
if (nrow(k)) { cat("\nKCNJ11 胰岛通道:\n"); print(k[, .(outcome, islet_lead, islet_slope, mr_b, mr_p, concordant)]) }
