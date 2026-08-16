#!/usr/bin/env Rscript
# =============================================================================
# M6_islet_mr.R — 胰岛 eQTL 通道 cis-MR（P1-1b：InsPIRE 人类胰岛 eQTL）
# =============================================================================
# 数据：InsPIRE（Viñuela 2020 Nat Commun, Zenodo 3408356）420 供体胰岛，
#       PacreaticIslets_independent_gene_eQTLs.txt（已 clump 的独立 eQTL，4639 行）
# 用途：对"全血看不见"的 β细胞/胰岛基因（ABCC8/PDX1/GCG/KCNJ11/GLP1R/TCF7L2 等）
#       用胰岛 eQTL 重跑 MR → 若位点"活过来"，证明是全血代理错了，非基因无信号。
# 方法（对齐预注册）：cis ±1Mb、工具为 InsPIRE 已 clump 独立 eQTL；
#       主 mr_ivw_mre / 单工具 Wald；se 由 Slope 与 Nominal_Pval 反推（z=slope/se）。
# eaf/方向诚实处理：MAF 列语义不保证是效应等位频率 → palindromic 一律排除（保守），
#       MAF 作 eaf 近似并 note 标注（与 deCODE 兜底路径同口径）。
# 用法：PATH=$R_ENV/bin:$PATH Rscript scripts/M6_islet_mr.R
# 输出：results/grid/islet_mr.csv + results/grid/islet_avail.csv
# =============================================================================
suppressMessages({
  library(data.table)
  suppressPackageStartupMessages(library(TwoSampleMR))
  suppressPackageStartupMessages(library(ieugwasr))
  library(jsonlite)
})
proj <- "/data/qiushuogeng/projects/dual-channel-mr-atlas"
res  <- file.path(proj, "results")
dir.create(file.path(res, "grid"), recursive = TRUE, showWarnings = FALSE)
log <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n")

prereg <- file.path(proj, "docs/PREREGISTRATION.md")
lock   <- file.path(proj, "docs/PREREGISTRATION.md.sha256")
stopifnot(file.exists(prereg), file.exists(lock), tools::md5sum(prereg) == readLines(lock))
cfg <- read_json(file.path(res, "config.json"), auto_unbox = TRUE)

OUTCOMES <- list(t2d = "ebi-a-GCST006867", cad = "ebi-a-GCST005194", fbg = "ebi-a-GCST005186")

# 关注基因：T2D 药物靶点 + 脂质对照 + 代表性（胰岛/β细胞相关优先）
GENES <- c("GLP1R","ABCC8","KCNJ11","PDX1","GCG","TCF7L2","INSR","SLC5A2","DPP4",
           "PPARG","G6PC","PCK1","SLC5A1","PCSK9","HMGCR","ANGPTL3","APOC3","CETP","NPC1L1")

is <- fread(file.path(proj, "data/gtex/islet/InsPIRE_islets_independent_gene_eQTLs.txt"),
            sep = "\t", header = TRUE)
log("InsPIRE 独立基因 eQTL: ", nrow(is), " 行")
# 可用性
av <- is[, .(n_eqtl = .N, min_p = min(Nominal_Pval)), by = GeneName]
av[, islet_available := n_eqtl > 0]
write.csv(av, file.path(res, "grid/islet_avail.csv"), row.names = FALSE)
cat("=== 关注基因在胰岛 eQTL 的可用性 ===\n")
for (g in GENES) {
  r <- av[GeneName == g]
  cat(sprintf("  %-8s %s\n", g, if(nrow(r)) sprintf("eQTL=%d minP=%s", r$n_eqtl, formatC(r$min_p,2,format="g")) else "（胰岛无独立 eQTL）"))
}

is_pal <- function(a,b) (a=="A"&b=="T")|(a=="T"&b=="A")|(a=="C"&b=="G")|(a=="G"&b=="C")

run_islet_mr <- function(gene, outcome_id) {
  d <- is[GeneName == gene]
  if (nrow(d) == 0)
    return(data.frame(gene = gene, outcome = outcome_id, nsnp = 0, b = NA, se = NA,
                      pval = NA, method = NA, ok = FALSE, note = "胰岛无独立 eQTL"))
  d[, rsid := SNPid]
  d <- d[!duplicated(rsid)]
  # 2026-08-07 P1 修复：p<~4.4e-16 时 1-p/2 舍入为 1 → qnorm=Inf → se=0 → Wald b/se=Inf p=0（假性过显著）。
  # 2026-08-16 复核修正：封顶必须用 pmin（取 z 与 8.209536 的较小值），原 pmax 在 z=Inf 时
  #   仍得 Inf → se=0 未真正封住；且对中等 z（如 5.5）pmax 会把 z 抬高到 8.2，se 低估 ~33%、
  #   显著性被夸大。pmin 才使 se=|Slope|/min(z,8.2) ≥ |Slope|/8.2 > 0，且小 z 不误伤。
  d[, se := abs(Slope) / pmin(qnorm(1 - Nominal_Pval/2), 8.209536)]   # z = slope/se
  d[, eaf := MAF]
  npal <- d[, sum(is_pal(A1, A2))]
  d <- d[!is_pal(A1, A2)]                              # palindromic 保守排除（MAF 非效应等位频率）
  if (nrow(d) == 0)
    return(data.frame(gene = gene, outcome = outcome_id, nsnp = 0, b = NA, se = NA,
                      pval = NA, method = NA, ok = FALSE, note = paste0("palindromic 排除 ", npal, " 后无工具")))
  exp_dat <- data.frame(SNP = d$rsid, effect_allele.exposure = d$A1,
    other_allele.exposure = d$A2, eaf.exposure = d$eaf,
    beta.exposure = d$Slope, se.exposure = d$se,
    pval.exposure = d$Nominal_Pval, samplesize.exposure = 420,
    id.exposure = gene, exposure = gene, stringsAsFactors = FALSE)
  tryCatch({
    out <- extract_outcome_data(snps = d$rsid, outcomes = outcome_id, proxies = TRUE)
    if (is.null(out) || nrow(out) == 0)
      return(data.frame(gene = gene, outcome = outcome_id, nsnp = nrow(d), b = NA, se = NA,
                        pval = NA, method = NA, ok = FALSE, note = "无结局匹配"))
    dat <- harmonise_data(exp_dat, out, action = 2); dat <- dat[dat$mr_keep, ]
    if (nrow(dat) == 0)
      return(data.frame(gene = gene, outcome = outcome_id, nsnp = nrow(d), b = NA, se = NA,
                        pval = NA, method = NA, ok = FALSE, note = "harmonise 后无保留 SNP"))
    r_main <- mr(dat, method_list = cfg$mr_methods$primary)
    r_sens <- tryCatch(mr(dat, method_list = cfg$mr_methods$sensitivity), error = function(e) NULL)
    r_all <- if (!is.null(r_sens)) rbind(r_main, r_sens) else r_main
    if (nrow(r_all) == 0 && nrow(dat) == 1) {
      r_w <- tryCatch(mr(dat, method_list = "mr_wald_ratio"), error = function(e) NULL)
      if (!is.null(r_w) && nrow(r_w) > 0) r_all <- r_w
    }
    if (nrow(r_all) == 0)
      return(data.frame(gene = gene, outcome = outcome_id, nsnp = nrow(dat), b = NA, se = NA,
                        pval = NA, method = NA, ok = FALSE, note = "MR 无输出"))
    data.frame(gene = gene, outcome = outcome_id, nsnp = nrow(dat),
               b = r_all$b, se = r_all$se, pval = r_all$pval, method = r_all$method,
               ok = TRUE, note = paste0("palindromic 排除 ", npal, "; harmonise ", nrow(dat)))
  }, error = function(e) data.frame(gene = gene, outcome = outcome_id, nsnp = 0, b = NA,
                                    se = NA, pval = NA, method = NA, ok = FALSE,
                                    note = conditionMessage(e)))
}

out <- list()
for (on in names(OUTCOMES)) {
  log("=== 结局: ", on, " ===")
  for (g in GENES) {
    rr <- run_islet_mr(g, OUTCOMES[[on]])
    out[[length(out) + 1]] <- rr
    for (i in seq_len(nrow(rr)))
      if (rr$ok[i]) log("  ", g, "×", on, "[", rr$method[i], "] nsnp=", rr$nsnp[i],
                        " b=", round(rr$b[i], 3), " p=", format(rr$pval[i], digits = 2))
  }
}
res_df <- do.call(rbind, out)
write.csv(res_df, file.path(res, "grid/islet_mr.csv"), row.names = FALSE)
log("胰岛 eQTL MR 落盘 ✔ → results/grid/islet_mr.csv")
