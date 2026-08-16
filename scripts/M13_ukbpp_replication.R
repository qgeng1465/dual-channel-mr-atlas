#!/usr/bin/env Rscript
# =============================================================================
# M13_ukbpp_replication.R — UKB-PPP 跨平台蛋白复现（探索性，v2）
# =============================================================================
# 数据：Sun 2023 UKB-PPP（Olink，discovery n≈34k，英国生物库），Synapse 下载。
#   REGENIE：CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P
#   GENPOS=hg38（实测：PCSK9 top GENPOS 55,039,974 = hg38，ID 列 hg19）。Beta=每等位效应（蛋白 SD 单位）。
# 面板 8/11：PCSK9 INSR ANGPTL3 LDLR APOB DPP4 GLP1R GCG（无 APOC3/PCK1/HMGCR）
# 复现两层：
#   A. pQTL 全区域方向一致率 + beta 相关：每蛋白 deCODE cis ∩ UKB-PPP cis 共享变异
#      （hg38 位置+等位集合匹配，palindromic 保守跳过）→ 方向一致率 + Pearson r。
#      deCODE/UKB Beta 同为每等位蛋白 SD 单位，尺度直接可比。
#   B. 单工具 MR 复现：UKB-PPP 自己 cis 区最强变异（max LOG10P）→ T2D/CAD/FBG，
#      与 deCODE MR 方向/显著性对比。rsID 由 deCODE 位置匹配获得（UKB-PPP 无 rsID）。
# 纪律：探索性，独立输出 CSV；匹配失败/弱信号如实报告；只做"工具与方向"一致性，不声称新发现。
# 用法：PATH=$R_ENV/bin:$PATH Rscript scripts/M13_ukbpp_replication.R
# 输出：results/ukbpp/ukbpp_pqtl_concordance.csv + ukbpp_mr_replication.csv
# =============================================================================
suppressMessages({
  library(data.table)
  suppressPackageStartupMessages(library(TwoSampleMR))
  suppressPackageStartupMessages(library(ieugwasr))
})
proj <- "/data/qiushuogeng/projects/dual-channel-mr-atlas"
res  <- file.path(proj, "results")
UKB  <- file.path(proj, "data/ukbpp")
SUB  <- file.path(proj, "data/decode/sub")
dir.create(file.path(res, "ukbpp"), showWarnings = FALSE)
log <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n", flush = TRUE)

PROTEINS <- list(
  list(gene = "PCSK9",   chr = "1",  tss = 55039445, tar = "PCSK9_Q8NBP7_OID20235_v1_Cardiometabolic"),
  list(gene = "INSR",    chr = "19", tss = 7294443,  tar = "INSR_P06213_OID30566_v1_Inflammation_II"),
  list(gene = "ANGPTL3", chr = "1",  tss = 62597464, tar = "ANGPTL3_Q9Y5C1_OID20407_v1_Cardiometabolic"),
  list(gene = "LDLR",    chr = "19", tss = 11089418, tar = "LDLR_P01130_OID20240_v1_Cardiometabolic"),
  list(gene = "APOB",    chr = "2",  tss = 21044075, tar = "APOB_P04114_OID30673_v1_Inflammation_II"),
  list(gene = "DPP4",    chr = "2",  tss = 162074639,tar = "DPP4_P27487_OID20406_v1_Cardiometabolic"),
  list(gene = "GLP1R",   chr = "6",  tss = 39048562, tar = "GLP1R_P43220_OID30809_v1_Neurology_II"),
  list(gene = "GCG",     chr = "2",  tss = 162152404,tar = "GCG_P01275_OID21263_v1_Oncology")
)
OUTCOMES <- list(t2d = "ebi-a-GCST006867", cad = "ebi-a-GCST005194", fbg = "ebi-a-GCST005186")
is_pal <- function(a, b) (a == "A" & b == "T") | (a == "T" & b == "A") |
  (a == "C" & b == "G") | (a == "G" & b == "C")

# --- deCODE 全部 cis 变异（hg38 pos + 等位 + beta，SD 单位）--------------------
load_decode_cis <- function(gene) {
  f <- list.files(SUB, pattern = paste0("_", gene, "_"), full.names = TRUE)
  f <- f[grepl("_cis.txt.gz$", f)]
  if (length(f) == 0) return(NULL)
  d <- fread(cmd = paste0("zcat ", shQuote(f[1])), sep = "\t", header = TRUE, nThread = 2)
  if (!"Pval" %in% names(d)) d[, Pval := 10^(-minus_log10_pval)]
  d[, pos := as.numeric(Pos)]
  d[, .(pos, rsid = rsids, eff = effectAllele, other = otherAllele, beta = Beta, P = Pval)]
}
# --- UKB-PPP cis（hg38 GENPOS，TSS±1Mb）--------------------------------------
load_ukbpp_cis <- function(p) {
  dirp <- file.path(UKB, p$tar)
  if (!dir.exists(dirp)) return(NULL)
  f <- list.files(dirp, pattern = paste0("discovery_chr", p$chr, "_"), full.names = TRUE)
  if (length(f) == 0) return(NULL)
  d <- fread(cmd = paste0("zcat ", shQuote(f[1])), sep = " ", header = TRUE, nThread = 4)
  d[CHROM == as.integer(p$chr) & GENPOS >= p$tss - 1e6 & GENPOS <= p$tss + 1e6,
    .(pos = GENPOS, eff = ALLELE1, other = ALLELE0, beta = BETA, P = 10^(-LOG10P),
      se = SE, n = N, eaf = A1FREQ)]
}

# ===== A. pQTL 全区域方向一致率 + beta 相关 ===================================
# deCODE/UKB-PPP Beta 同为每等位蛋白 SD 单位 → 方向一致率（sign match）+ Pearson r 可比。
log("=== A. pQTL 共享变异方向一致率 + beta 相关 ===")
conc <- rbindlist(lapply(PROTEINS, function(p) {
  dc <- load_decode_cis(p$gene)
  uk <- load_ukbpp_cis(p)
  if (is.null(dc) || is.null(uk)) {
    return(data.table(gene = p$gene, n_decode_cis = if (is.null(dc)) NA else nrow(dc),
                      n_ukbpp_cis = if (is.null(uk)) NA else nrow(uk),
                      n_shared = 0, concordant_n = NA, concordant_pct = NA,
                      pearson_r = NA, pearson_p = NA, ukbpp_top_pos = NA, ukbpp_top_P = NA))
  }
  m <- merge(dc, uk, by = "pos", suffixes = c(".dc", ".uk"))
  m <- m[!(is_pal(eff.dc, other.dc))]            # 回文（A/T、C/G）方向歧义，保守排除
  m[, beta_uk_aligned := ifelse(eff.dc == eff.uk & other.dc == other.uk, beta.uk,
                         ifelse(eff.dc == other.uk & other.dc == eff.uk, -beta.uk, NA))]
  m <- m[!is.na(beta_uk_aligned)]
  uk_top <- uk[which.min(P)]                     # 最显著 cis 变异（min P）
  if (nrow(m) == 0) {
    return(data.table(gene = p$gene, n_decode_cis = nrow(dc), n_ukbpp_cis = nrow(uk),
                      n_shared = 0, concordant_n = 0, concordant_pct = NA,
                      pearson_r = NA, pearson_p = NA,
                      ukbpp_top_pos = uk_top$pos, ukbpp_top_P = uk_top$P))
  }
  concordant_n <- m[, sum(sign(beta.dc) == sign(beta_uk_aligned))]
  rp <- tryCatch(cor.test(m$beta.dc, m$beta_uk_aligned), error = function(e) NULL)
  data.table(gene = p$gene, n_decode_cis = nrow(dc), n_ukbpp_cis = nrow(uk),
             n_shared = nrow(m), concordant_n = concordant_n,
             concordant_pct = concordant_n / nrow(m),
             pearson_r = if (!is.null(rp)) unname(rp$estimate) else NA,
             pearson_p = if (!is.null(rp)) rp$p.value else NA,
             ukbpp_top_pos = uk_top$pos, ukbpp_top_P = uk_top$P)
}))
write.csv(conc, file.path(res, "ukbpp/ukbpp_pqtl_concordance.csv"), row.names = FALSE)
log("  A 完成：共享变异方向一致率 + Pearson r 已写")

# ===== B. 单工具 MR 复现（UKB-PPP cis top → T2D/CAD/FBG）=====================
log("=== B. UKB-PPP 单工具 MR → 结局 ===")
mr_rows <- list()
for (p in PROTEINS) {
  uk <- load_ukbpp_cis(p)
  if (is.null(uk) || nrow(uk) == 0) {
    for (on in names(OUTCOMES))
      mr_rows[[length(mr_rows) + 1]] <- data.table(gene = p$gene, outcome = on, rsid = NA,
                                                  ukbpp_b = NA, ukbpp_p = NA, nsnp = 0,
                                                  ukbpp_top_P = NA, decode_b = NA, decode_p = NA,
                                                  note = "UKB-PPP cis 不可读")
    next
  }
  top <- uk[which.min(P)]                       # 最显著 cis 变异（修正：min P）
  # rsID：UKB-PPP top → deCODE 位置匹配（须同位点，防误配远处非同源变异；
  #   2026-08-16 复核：旧版仅取最近邻，曾把 PCSK9 配到 400kb 外 rs4463718）
  dc <- load_decode_cis(p$gene)
  dc_match <- if (!is.null(dc) && nrow(dc) > 0) dc[which.min(abs(pos - top$pos))] else NULL
  top_rsid <- if (!is.null(dc_match) && abs(dc_match$pos - top$pos) <= 1) dc_match$rsid else NA
  if (is.na(top_rsid) || top_rsid %in% c("", ".", "-")) {
    for (on in names(OUTCOMES))
      mr_rows[[length(mr_rows) + 1]] <- data.table(gene = p$gene, outcome = on,
                                                  rsid = top_rsid, ukbpp_b = NA, ukbpp_p = NA,
                                                  nsnp = 0, ukbpp_top_P = top$P, decode_b = NA, decode_p = NA,
                                                  note = paste0("UKB-PPP top (", top$pos, ") 无 deCODE rsID 匹配"))
    next
  }
  # 对齐 UKB top 到 deCODE 同等位（beta 符号统一），供 deCODE MR 方向对照
  drow_all <- fread(file.path(res, "grid/protein_decode_mr.csv"))
  dec_res <- drow_all[gene == p$gene & ok == TRUE]
  exp_dat <- data.frame(
    SNP = top_rsid, effect_allele.exposure = top$eff, other_allele.exposure = top$other,
    eaf.exposure = top$eaf, beta.exposure = top$beta, se.exposure = top$se,
    pval.exposure = top$P, samplesize.exposure = top$n,
    id.exposure = p$gene, exposure = p$gene, stringsAsFactors = FALSE)
  for (on in names(OUTCOMES)) {
    rr <- tryCatch({
      out <- extract_outcome_data(snps = exp_dat$SNP, outcomes = OUTCOMES[[on]], proxies = TRUE)
      dat <- harmonise_data(exp_dat, out, action = 2)
      dat <- dat[dat$mr_keep, ]
      if (nrow(dat) == 0) NULL else { r <- mr(dat, method_list = "mr_wald_ratio"); data.table(b = r$b, p = r$pval) }
    }, error = function(e) NULL)
    drow <- dec_res[outcome == OUTCOMES[[on]]][1]
    mr_rows[[length(mr_rows) + 1]] <- data.table(
      gene = p$gene, outcome = on, rsid = top_rsid,
      ukbpp_b = if (!is.null(rr)) rr$b else NA, ukbpp_p = if (!is.null(rr)) rr$p else NA,
      nsnp = if (!is.null(rr)) 1 else 0, ukbpp_top_P = top$P,
      decode_b = if (nrow(drow) > 0) drow$b else NA, decode_p = if (nrow(drow) > 0) drow$pval else NA,
      note = if (is.null(rr)) "outcome 提取/调和失败" else
        paste0("方向一致=", if (!is.na(drow$b) && !is.na(rr$b)) sign(drow$b) == sign(rr$b) else NA))
  }
}
mr_df <- rbindlist(mr_rows, fill = TRUE)
write.csv(mr_df, file.path(res, "ukbpp/ukbpp_mr_replication.csv"), row.names = FALSE)
log("=== 完成 ✔ → results/ukbpp/ukbpp_pqtl_concordance.csv + ukbpp_mr_replication.csv")
