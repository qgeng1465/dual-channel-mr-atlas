#!/usr/bin/env Rscript
# =============================================================================
# M23_full_scan.R — 全量 coloc 扫描 v2（按染色体循环，图谱资源精确计数）
# =============================================================================
# v2 改动（2026-08-15 16:10）：v1 按基因块 awk 抽取，3 进程并发扫同一 4.3G gz，
#   每块 12-14 分钟太慢。改为：eQTLGen full 已由 M23b 一次抽取为按染色体明文
#   （<scratch>/eqtlgen_stable/bychr/chr{1..22}.tsv），本脚本按染色体
#   fread 明文块处理，无重复解压、内存 ~1-2GB/进程。
# 健全性：MR sig 集内 strong（PP.H4≥0.8）应与 transcript_coloc_hits.csv 的 106 一致。
# 用法：PATH=$R_ENV/bin:$PATH TMPDIR=<scratch>/rtmp \
#         Rscript scripts/M23_full_scan.R t2d|cad|fbg
# 输出：results/coloc_full_{outcome}_20260815.csv
# =============================================================================
suppressMessages({
  library(data.table)
  suppressPackageStartupMessages(library(coloc))
})
args <- commandArgs(trailingOnly = TRUE)
on <- args[1]
if (is.na(on) || !on %in% c("t2d", "cad", "fbg")) stop("用法: M23_full_scan.R t2d|cad|fbg")
proj <- "<repo-root>"
res  <- file.path(proj, "results")
gdir <- file.path(res, "grid")
GWAS <- file.path(proj, "data/opengwas/full")
STABLE <- "<scratch>/eqtlgen_stable"
SIG   <- file.path(STABLE, "cis-EQTL-significant.txt.gz")
AF    <- file.path(STABLE, "SNP_AF.txt.gz")
BYCHR <- file.path(STABLE, "bychr")
stopifnot(file.exists(file.path(BYCHR, "DONE")), file.exists(file.path(GWAS, "DONE")))
log <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n")
WINDOW  <- 1e6
TYPE    <- c(t2d = "cc", cad = "cc", fbg = "quant")
S_CASE  <- c(t2d = 61714/655666, cad = 34541/296525, fbg = NA)
GWAS_N  <- c(t2d = 655666, cad = 296525, fbg = 58074)

# ---- 1. 全量对（该结局）+ 坐标 ----
grid <- fread(file.path(gdir, "transcript_grid_mr.csv"))
grid <- grid[ok == TRUE & !is.na(pval) & outcome == on]
log("结局 ", on, " 对: ", nrow(grid))
gcoord <- unique(fread(cmd = paste0("zcat ", SIG), sep = "\t", header = TRUE, nThread = 4,
                       select = c("Gene", "GeneChr", "GenePos")))
grid <- merge(grid, gcoord, by.x = "gene", by.y = "Gene", all.x = TRUE)[!is.na(GeneChr)]
log("有坐标: ", nrow(grid))

# ---- 2. GWAS 侧 ----
load_gwas <- function(on) {
  f <- file.path(GWAS, paste0(on, "_full.gz"))
  if (on == "fbg") {
    t <- fread(cmd = paste0("zcat ", f), sep = "\t", header = TRUE, nThread = 4,
               select = c("Snp", "effect_allele", "other_allele", "maf", "MainEffects", "MainSE", "MainP"))
    t[, `:=`(rsid = Snp, ea = effect_allele, nea = other_allele, eaf = as.numeric(maf),
             beta = as.numeric(MainEffects), se = as.numeric(MainSE), p = as.numeric(MainP),
             chr = NA_integer_, pos = NA_integer_)]
    setkey(t, rsid)
    return(t[, .(rsid, ea, nea, eaf, beta, se, p, chr, pos)])
  }
  sel <- c("hm_rsid","hm_chrom","hm_pos","hm_effect_allele","hm_other_allele","hm_beta",
           "hm_effect_allele_frequency","standard_error","p_value")
  t <- fread(cmd = paste0("zcat ", f), sep = "\t", header = TRUE, nThread = 4, select = sel)
  t[, `:=`(rsid = as.character(hm_rsid), ea = as.character(hm_effect_allele),
           nea = as.character(hm_other_allele), eaf = as.numeric(hm_effect_allele_frequency),
           beta = as.numeric(hm_beta), se = as.numeric(standard_error),
           p = as.numeric(p_value), chr = as.integer(hm_chrom), pos = as.integer(hm_pos))]
  t <- t[!is.na(beta) & !is.na(se) & se > 0 & !is.na(p) & !is.na(chr) & !is.na(pos) & !is.na(rsid)]
  setkey(t, chr)
  t[, .(rsid, ea, nea, eaf, beta, se, p, chr, pos)]
}
gwas <- load_gwas(on)
log(on, " GWAS 行: ", nrow(gwas))

# ---- 3. harmonize + coloc ----
harmonize <- function(e, g) {
  e2 <- copy(e[, .(SNP, AssessedAllele, OtherAllele, eaf_e = eaf, beta_e = beta,
                   se_e = se, N_e = NrSamples)])
  g2 <- copy(g[, .(rsid, ea, nea, eaf_g = as.numeric(eaf),
                   beta_g = as.numeric(beta), se_g = as.numeric(se))])
  setkey(g2, rsid); setkey(e2, SNP)
  m <- merge(e2, g2, by.x = "SNP", by.y = "rsid")
  if (nrow(m) == 0) return(data.table())
  A <- toupper(m$AssessedAllele); O <- toupper(m$OtherAllele)
  ea <- toupper(m$ea); nea <- toupper(m$nea)
  gb <- m$beta_g
  ok   <- A == ea; ok[is.na(ok)] <- FALSE
  flip <- !ok & A == nea & O == ea; flip[is.na(flip)] <- FALSE
  gb[flip] <- -gb[flip]
  pal <- !ok & !flip & (((A == "A" & O == "T") | (A == "T" & O == "A")) |
                        ((A == "C" & O == "G") | (A == "G" & O == "C")))
  pal[is.na(pal)] <- FALSE
  p_flip <- pal & ((m$eaf_e > 0.58 & m$eaf_g < 0.42) | (m$eaf_e < 0.42 & m$eaf_g > 0.58))
  p_flip[is.na(p_flip)] <- FALSE
  gb[p_flip] <- -gb[p_flip]
  keep <- ok | flip | p_flip
  out <- data.table(snp = m$SNP[keep],
                    maf = pmin(m$eaf_e, 1 - m$eaf_e)[keep],
                    beta_e = m$beta_e[keep], varbeta_e = m$se_e[keep]^2, N_e = m$N_e[keep],
                    beta_g = gb[keep], varbeta_g = m$se_g[keep]^2)
  out <- out[maf >= 0.01 & maf <= 0.99]
  out[!is.na(beta_e) & !is.na(beta_g) & !is.na(varbeta_e) & !is.na(varbeta_g) &
        varbeta_e > 0 & varbeta_g > 0]
}
run_coloc <- function(g, e) {
  if (nrow(e) == 0) return(list(ok = FALSE, note = "eQTL cis 无变异"))
  m <- harmonize(e, g)
  if (nrow(m) < 10) return(list(ok = FALSE, note = paste0("对齐后仅 ", nrow(m), " 变异(<10)")))
  d1 <- list(snp = m$snp, type = "quant", N = m$N_e, beta = m$beta_e, varbeta = m$varbeta_e, MAF = m$maf)
  if (TYPE[[on]] == "cc")
    d2 <- list(snp = m$snp, type = "cc", N = GWAS_N[[on]], s = S_CASE[[on]],
               beta = m$beta_g, varbeta = m$varbeta_g, MAF = m$maf)
  else
    d2 <- list(snp = m$snp, type = "quant", N = GWAS_N[[on]], beta = m$beta_g, varbeta = m$varbeta_g, MAF = m$maf)
  r <- tryCatch({ capture.output(res <- coloc.abf(d1, d2, p12 = 1e-5), type = "output"); res },
                error = function(e) NULL)
  if (is.null(r)) return(list(ok = FALSE, note = "coloc.abf 失败"))
  pp4 <- unname(as.numeric(r$summary[["PP.H4.abf"]]))
  list(ok = TRUE, note = "", nsnp = nrow(m), pp4 = pp4)
}

# ---- 4. AF 表（全局）----
af <- fread(cmd = paste0("zcat ", AF), sep = "\t", header = TRUE, nThread = 4,
            select = c("SNP", "AlleleA", "AlleleB", "AlleleB_all"))

# ---- 5. 按染色体循环（明文块，无重复解压）----
out_all <- vector("list", 22)
t0 <- proc.time()
for (cc in 1:22) {
  fchr <- file.path(BYCHR, paste0("chr", cc, ".tsv"))
  if (!file.exists(fchr)) { log("chr", cc, " 无文件，跳过"); next }
  eqtl <- fread(fchr, sep = "\t", header = TRUE, nThread = 4)
  eqtl <- eqtl[abs(SNPPos - GenePos) <= WINDOW]          # 染色体已由文件限定
  eqtl <- merge(eqtl, af, by = "SNP", all.x = TRUE)
  eqtl[, eaf := fcase(AssessedAllele == AlleleA, 1 - AlleleB_all,
                      AssessedAllele == AlleleB, AlleleB_all, default = NA_real_)][!is.na(eaf)]
  eqtl[, `:=`(beta = Zscore / sqrt(NrSamples), se = 1 / sqrt(NrSamples))]
  eqtl <- unique(eqtl[, .(Gene, SNP, AssessedAllele, OtherAllele, eaf, beta, se, NrSamples)])
  setkey(eqtl, Gene)
  sub <- grid[GeneChr == cc]
  if (nrow(sub) == 0) { rm(eqtl); gc(); next }
  rows <- vector("list", nrow(sub))
  for (i in seq_len(nrow(sub))) {
    pr <- sub[i]
    e <- eqtl[Gene == pr$gene]
    # 2026-08-15 修复（build 错位）：OpenGWAS full 的 hm_pos 是 hg38，eQTLGen 是 hg19，
    #   位置窗口提取错位致 rsid 交集骤减（GIT1: 应 1920 实际 190）。全结局统一改用
    #   纯 rsid 匹配（不依赖坐标 build），与 fbg 原逻辑一致。
    g <- gwas[rsid %in% e$SNP]
    g_min_p <- if (nrow(g)) min(g$p, na.rm = TRUE) else NA
    cl <- run_coloc(g, e)
    e_fmax <- if (nrow(e)) max(abs(e$beta) / e$se, na.rm = TRUE) else NA
    rows[[i]] <- data.table(gene = pr$gene, symbol = pr$symbol, outcome = pr$outcome,
      mr_b = pr$b, mr_p = pr$pval, gwas_min_p = g_min_p, eqtl_F_max = e_fmax,
      nsnp = if (cl$ok) cl$nsnp else NA, pp4 = if (cl$ok) cl$pp4 else NA,
      ok = cl$ok, note = cl$note)
  }
  out_all[[cc]] <- rbindlist(rows, fill = TRUE)
  rm(eqtl, rows, sub); gc()
  log("  chr", cc, " 完成 (", nrow(out_all[[cc]]), " 行)")
}
res_df <- rbindlist(out_all, fill = TRUE)
res_df[!is.na(pp4), pp4 := as.numeric(pp4)]
fwrite(res_df, file.path(res, paste0("coloc_full_", on, "_20260815.csv")))
log("已写 results/coloc_full_", on, "_20260815.csv | ", nrow(res_df), " 行 | ",
    round((proc.time()-t0)[3]/60, 1), " min")

# ---- 6. 摘要 ----
cat("\n=== 全量 coloc 摘要 ", on, " ===\n")
d <- res_df
cat("全部: ", nrow(d), " | QC 通过: ", nrow(d[ok==TRUE]), "\n")
for (g0 in c("sig", "grey", "null")) {
  dd <- d[if (g0=="sig") mr_p < 0.05 else if (g0=="grey") mr_p >= 0.05 & mr_p < 0.5 else mr_p >= 0.5]
  h80 <- dd[ok==TRUE & pp4 >= 0.8]
  h50 <- dd[ok==TRUE & pp4 >= 0.5]
  cat(g0, ": ", nrow(dd), " 对 | PP.H4≥0.8: ", nrow(h80),
      " (", round(100*nrow(h80)/max(1,nrow(dd[ok==TRUE])),2), "% of QC)",
      " | ≥0.5: ", nrow(h50), "\n", sep="")
  if (g0 == "sig" && nrow(h80)) {
    cat("   sig 内 strong 按 GWAS 峰显著: ", sum(h80$gwas_min_p<5e-8, na.rm=TRUE), "/", nrow(h80), "\n")
  }
}
