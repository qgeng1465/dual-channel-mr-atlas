#!/usr/bin/env Rscript
# =============================================================================
# M6_gtex_mr.R (v2) — GTEx 组织 eQTL 通道 cis-MR（P1-1a：肝/胰/全血等组织）
# =============================================================================
# 数据：GTEx v8 egenes（Zenodo GV-Rep 镜像 zip；各组织每基因 top cis-eQTL 行）
#   列：gene_id, gene_name, variant_id, tss_distance, chr, variant_pos, ref, alt,
#       rs_id_dbSNP151_GRCh38p7, maf, pval_nominal, slope, slope_se
# v2 修正（相对 signif_variant_gene_pairs 版）：
#   - gene_id 是 ENSG 无 symbol → 改用 gene_name 直接匹配候选基因（v1 全空根因）
#   - signif 对 variant_id(chr_pos_ref_alt_b38) 无 rsid，OpenGWAS 结局按 rsid 检索
#     → 改用 egenes 自带 rs_id_dbSNP151_GRCh38p7（可查结局）；nsnp=1 → Wald
# 设计要点：
#   - cis 过滤直接用 tss_distance（±1Mb），不依赖坐标 build
#   - 工具 = 每基因该组织 top cis-eQTL（p_nominal<5e-6 才纳入），单工具 Wald 敏感性
#   - 组织独立报告：肝 Liver、胰 Pancreas、全血 Whole_Blood 为核心（+Adipose/Muscle/Artery）
# 用途：若肝/胰 eQTL 让 PCSK9/ANGPTL3/APOC3 等全血失明位点"活过来" → 组织特异性叙事救场。
#   注意：单工具敏感性，主分析仍锚定全血 eQTL×血浆 pQTL 四态（不移动门柱）。
# 用法：PATH=$R_ENV/bin:$PATH Rscript scripts/M6_gtex_mr.R
# 输出：results/grid/gtex_mr.csv + results/grid/gtex_avail.csv
# =============================================================================
suppressMessages({
  library(data.table)
  suppressPackageStartupMessages(library(TwoSampleMR))
  suppressPackageStartupMessages(library(ieugwasr))
  library(jsonlite)
})
proj <- "/data/qiushuogeng/projects/dual-channel-mr-atlas"
res  <- file.path(proj, "results")
GTEX <- file.path(proj, "data/gtex")
dir.create(file.path(res, "grid"), recursive = TRUE, showWarnings = FALSE)
log <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n")

prereg <- file.path(proj, "docs/PREREGISTRATION.md")
lock   <- file.path(proj, "docs/PREREGISTRATION.md.sha256")
stopifnot(file.exists(prereg), file.exists(lock), tools::md5sum(prereg) == readLines(lock))
cfg <- read_json(file.path(res, "config.json"), auto_unbox = TRUE)
OUTCOMES <- list(t2d = "ebi-a-GCST006867", cad = "ebi-a-GCST005194", fbg = "ebi-a-GCST005186")

TISSUES <- c("Liver","Pancreas","Whole_Blood","Adipose_Subcutaneous","Muscle_Skeletal","Artery_Coronary")
GENES <- c("PCSK9","HMGCR","ANGPTL3","APOC3","APOB","LDLR","CETP","NPC1L1",
           "GLP1R","ABCC8","KCNJ11","PDX1","GCG","TCF7L2","INSR","SLC5A2","DPP4",
           "PPARG","G6PC","PCK1","SLC5A1")

# --- 从 zip 解出 egenes（若未解出）--------------------------------------------
zipf <- file.path(GTEX, "GTEx_Analysis_v8_eQTL.zip")
for (t in TISSUES) {
  out <- file.path(GTEX, paste0(t, ".egenes.txt.gz"))
  if (!file.exists(out) && file.exists(zipf)) {
    member <- paste0("GTEx_Analysis_v8_eQTL/", t, ".v8.egenes.txt.gz")
    ok <- system2("unzip", c("-o", "-j", zipf, member, "-d", file.path(GTEX, ".zipx")),
                  stdout = TRUE, stderr = TRUE)
    tmp <- file.path(GTEX, ".zipx", paste0(t, ".v8.egenes.txt.gz"))
    if (file.exists(tmp)) { file.rename(tmp, out); log("解压 egenes  ", t) }
  }
  if (!file.exists(out)) log("✘ ", t, " egenes 缺失：", out)
}
unlink(file.path(GTEX, ".zipx"), recursive = TRUE)

# --- 读取组织 egenes + 可用性 + 跑 MR ------------------------------------------
load_tissue <- function(t) {
  f <- file.path(GTEX, paste0(t, ".egenes.txt.gz"))
  if (!file.exists(f)) return(NULL)
  d <- fread(cmd = paste0("zcat ", shQuote(f)), sep = "\t", header = TRUE, nThread = 4)
  if (nrow(d) == 0) return(NULL)
  d[, Tissue := t]
  setnames(d, c("gene_name", "rs_id_dbSNP151_GRCh38p7"), c("symbol", "rsid"), skip_absent = TRUE)
  d[, rsid := fifelse(is.na(rsid) | rsid %in% c("", ".", "-"), variant_id, rsid)]
  d
}

# 单工具（每基因 top cis-eQTL）Wald 敏感性；主终点 T2D 优先打印。
run_gtex_mr <- function(gene, tissue, outcome_id) {
  d <- tis[[tissue]][symbol == gene & abs(tss_distance) <= 1e6]
  if (nrow(d) == 0)
    return(data.frame(tissue = tissue, gene = gene, outcome = outcome_id, nsnp = 0,
                      b = NA, se = NA, pval = NA, method = NA, ok = FALSE,
                      note = "该组织无该基因 egenes 记录（无 cis 信号）"))
  d <- d[pval_nominal < 5e-6]
  if (nrow(d) == 0)
    return(data.frame(tissue = tissue, gene = gene, outcome = outcome_id, nsnp = 0,
                      b = NA, se = NA, pval = NA, method = NA, ok = FALSE,
                      note = "cis eQTL p≥5e-6（未达预注册工具阈值）"))
  d <- d[order(pval_nominal)][1]                      # top cis-eQTL（egenes 即 lead）
  exp_dat <- data.frame(SNP = d$rsid,
    effect_allele.exposure = d$alt, other_allele.exposure = d$ref,
    # 2026-08-07 P1 修复：maf 是次等位频率，slope 是 alt 等位效应；ref_factor==1 表示 ref 为次等位
    # → eaf(alt) = ifelse(ref_factor==1, 1-maf, maf)。此前用 maf 当 eaf，回文位点 harmonise 会错误翻号。
    # NA ref_factor → eaf=NA，回文位点被 action=2 保守丢弃（与胰岛通道一致）。
    eaf.exposure = ifelse(d$ref_factor == 1, 1 - d$maf, d$maf),
    beta.exposure = d$slope, se.exposure = d$slope_se,
    pval.exposure = d$pval_nominal,
    samplesize.exposure = if ("minor_allele_samples" %in% names(d)) d$minor_allele_samples else NA,
    id.exposure = gene, exposure = gene, stringsAsFactors = FALSE)
  tryCatch({
    out <- extract_outcome_data(snps = d$rsid, outcomes = outcome_id, proxies = TRUE)
    if (is.null(out) || nrow(out) == 0)
      return(data.frame(tissue = tissue, gene = gene, outcome = outcome_id, nsnp = 1,
                        b = NA, se = NA, pval = NA, method = NA, ok = FALSE, note = "无结局匹配"))
    dat <- harmonise_data(exp_dat, out, action = 2); dat <- dat[dat$mr_keep, ]
    if (nrow(dat) == 0)
      return(data.frame(tissue = tissue, gene = gene, outcome = outcome_id, nsnp = 1,
                        b = NA, se = NA, pval = NA, method = NA, ok = FALSE,
                        note = "harmonise 无保留（等位不匹配/回文无 eaf 判向）"))
    rw <- tryCatch(mr(dat, method_list = "mr_wald_ratio"), error = function(e) NULL)
    if (is.null(rw) || nrow(rw) == 0)
      return(data.frame(tissue = tissue, gene = gene, outcome = outcome_id, nsnp = nrow(dat),
                        b = NA, se = NA, pval = NA, method = NA, ok = FALSE, note = "MR 无输出"))
    data.frame(tissue = tissue, gene = gene, outcome = outcome_id, nsnp = nrow(dat),
               b = rw$b, se = rw$se, pval = rw$pval, method = "Wald ratio",
               ok = TRUE, note = paste0("lead eQTL ", d$variant_id, " (p=",
                                        format(d$pval_nominal, digits = 3), ")"))
  }, error = function(e) data.frame(tissue = tissue, gene = gene, outcome = outcome_id, nsnp = 1,
                                    b = NA, se = NA, pval = NA, method = NA, ok = FALSE,
                                    note = conditionMessage(e)))
}

tis <- list()
for (t in TISSUES) { tis[[t]] <- load_tissue(t); if (!is.null(tis[[t]])) log("组织 ", t, ": ", nrow(tis[[t]]), " 基因 egenes 行") }
tis <- tis[!vapply(tis, is.null, logical(1))]
if (length(tis) == 0) { log("✘ 无组织 egenes 文件就绪。退出。"); quit(save = "no") }

# 可用性（基因 × 组织：该组织有 p<5e-6 cis eQTL）
av <- rbindlist(lapply(names(tis), function(t) {
  dd <- tis[[t]][symbol %in% GENES & abs(tss_distance) <= 1e6]
  dd[, .(tissue = t, n_cis_p5e6 = sum(pval_nominal < 5e-6), top_pval = min(pval_nominal)), by = symbol]
}))
write.csv(av, file.path(res, "grid/gtex_avail.csv"), row.names = FALSE)
cat("\n=== 基因 × 组织 cis eQTL 可用性（p<5e-6）===\n")
print(dcast(av, symbol ~ tissue, value.var = "n_cis_p5e6", fill = 0))

out <- list()
for (t in names(tis)) {
  genes_t <- intersect(GENES, unique(tis[[t]][symbol %in% GENES, symbol]))
  for (g in genes_t) {
    for (on in names(OUTCOMES)) {
      rr <- run_gtex_mr(g, t, OUTCOMES[[on]])
      out[[length(out) + 1]] <- rr
      for (i in seq_len(nrow(rr)))
        if (rr$ok[i]) log("  [", t, "] ", g, "×", on, " nsnp=", rr$nsnp[i],
                          " b=", round(rr$b[i], 3), " p=", format(rr$pval[i], digits = 2))
    }
  }
}
res_df <- do.call(rbind, out)
write.csv(res_df, file.path(res, "grid/gtex_mr.csv"), row.names = FALSE)
log("GTEx 组织 eQTL MR 落盘 ✔ → results/grid/gtex_mr.csv（单工具 Wald 敏感性；组织通道不移动四态主分析门柱）")
