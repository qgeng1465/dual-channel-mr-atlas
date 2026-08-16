#!/usr/bin/env Rscript
# =============================================================================
# M5_transcript_coloc.R — 转录通道 共定位闸门（coloc.abf，PP.H4 三档）
# =============================================================================
# 目的：对 stage-2 复核存活的 MR 命中（nominal p<0.05 且方向一致）做共定位闸门，
#   回答"MR 信号是否由共享因果变异驱动（PP.H4）而非 LD 混杂"。转录通道（eQTLGen
#   全血 cis-eQTL × OpenGWAS 结局）不依赖 deCODE，可立即运行。
# 数据（与 MR 同一数据源，保证一致性）：
#   - eQTL 侧：eQTLGen 全量 cis-eQTL 汇总（data/eqtlgen/cis-eQTLs_full_20180905.txt.gz，
#     14 列：Pvalue SNP SNPChr SNPPos Zscore AssessedAllele OtherAllele Gene
#     GeneSymbol GeneChr GenePos NrCohorts NrSamples FDR）→ beta=Z/sqrt(N), se=1/sqrt(N)
#   - GWAS 侧：OpenGWAS associations API **区域查询**（variant=chr:start-end，proxies=0），
#     返回该结局在该区域的全部变异（rsid/ea/nea/eaf/beta/se/p/n）——与 MR 用的同一套
#     调和后数据，无需下载整文件。per 对一次请求，逐对缓存（_coloc_gwas/）。
#   - 变异交集 + 等位基因对齐（palindromic 用 eAF 解析，MAF 双高歧义则剔除）。
# 设计（对齐预注册 §3 + README §0.4 M5）：
#   - 区域 = 基因 TSS（eQTLGen GenePos，hg19）± 1Mb（与 cis 窗口一致）
#   - coloc.abf 默认先验 p1=1e-4 p2=1e-4 p12=1e-5；敏感性 p12=1e-6
#   - type：eQTL=quant（N=NrSamples）；T2D/CAD=cc（N=n，s=gwasinfo 病例比例）；
#     FBG=quant（N=n）
#   - PP.H4 三档：≥0.8 强共定位 / 0.5–0.8 中等 / <0.5 无证据；连续 PP.H4 全量报告
#   - 全网格落盘（含失败行）；与 stage-2 MR 的 b 方向一致性和 q 并列报告
# 注意：样本重叠（eQTLGen 与 GWAS 部分队列重叠）为 coloc 已知局限，输出注明不作修正。
# 用法：cd 项目 && PATH=... Rscript scripts/M5_transcript_coloc.R
# 产物：results/grid/transcript_coloc.csv + _hits.csv + _funnel.tsv
# 缓存：results/grid/_coloc_gwas/<outcome>_<gene>.rds（逐对）、_coloc_cis_<gene>.rds
# =============================================================================
suppressMessages({
  library(data.table)
  suppressPackageStartupMessages(library(coloc))
  library(jsonlite); library(httr)
})
proj <- "/data/qiushuogeng/projects/dual-channel-mr-atlas"
res  <- file.path(proj, "results")
gdir <- file.path(res, "grid")
dir.create(file.path(gdir, "_coloc_gwas"), recursive = TRUE, showWarnings = FALSE)
log <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n", sep = "")

# --- 预注册完整性校验 ----------------------------------------------------------
prereg <- file.path(proj, "docs/PREREGISTRATION.md")
lock   <- file.path(proj, "docs/PREREGISTRATION.md.sha256")
stopifnot(file.exists(prereg), file.exists(lock),
          tools::md5sum(prereg) == readLines(lock))
cfg <- read_json(file.path(res, "config.json"), auto_unbox = TRUE)
OUTCOMES <- list(t2d = "ebi-a-GCST006867", cad = "ebi-a-GCST005194", fbg = "ebi-a-GCST005186")
OUT_N <- c(t2d = "GCST006867", cad = "GCST005194", fbg = "GCST005186")
TYPE <- c(t2d = "cc", cad = "cc", fbg = "quant")
# s = 病例比例（OpenGWAS gwasinfo；T2D ncontrol 字段损坏 → ncase/sample_size）
S_CASE <- c(t2d = 61714 / 655666, cad = 34541 / 296525)
GWAS_N <- c(t2d = 655666, cad = 296525, fbg = 58074)
log("预注册哈希校验通过 ✔ | M5 转录共定位闸门 | 区域 ±", cfg$instrument$cis_window_kb, "kb | PP.H4 三档")

FULL <- file.path(proj, "data/eqtlgen/cis-eQTLs_full_20180905.txt.gz")
SIG  <- file.path(proj, "data/eqtlgen/cis-eQTL-significant.txt.gz")
AF   <- file.path(proj, "data/eqtlgen/SNP_AF.txt.gz")
stopifnot(file.exists(SIG), file.exists(AF))

JWT <- Sys.getenv("OPENGWAS_JWT")
if (!nzchar(JWT)) stop("OPENGWAS_JWT 未设置（~/.Renviron）")
H <- httr::add_headers(Authorization = paste("Bearer", JWT))
API <- "https://api.opengwas.io/api/associations"

# --- 载入 stage-2 命中（闸门范围 = nominal 存活且方向一致）----------------------
s2 <- fread(file.path(gdir, "transcript_grid_stage2.csv"))
pairs <- s2[survive_nominal == TRUE & ok == TRUE]
log("闸门对: ", nrow(pairs), "（nominal 存活，唯一基因 ", uniqueN(pairs$gene), "）")

# --- 基因坐标（hg19，来自 eQTLGen 显著文件；全量 cis 未就绪时也能定位区域）-----
gcoord <- unique(fread(cmd = paste0("zcat ", SIG), sep = "\t", header = TRUE, nThread = 4,
                       select = c("Gene", "GeneChr", "GenePos")))
gcoord <- gcoord[!duplicated(Gene)]
pairs <- merge(pairs, gcoord, by.x = "gene", by.y = "Gene", all.x = TRUE)
pairs <- pairs[!is.na(GeneChr)]
log("基因坐标就绪: ", nrow(pairs), " 对")

# --- Phase A：GWAS 区域数据（逐对 range query，缓存，可断点续跑）----------------
fetch_region <- function(on, gene, gchr, gpos) {
  cf <- file.path(gdir, "_coloc_gwas", paste0(OUT_N[[on]], "_", gene, ".rds"))
  if (file.exists(cf)) return(readRDS(cf))
  start <- max(1, gpos - cfg$instrument$cis_window_kb * 1000)
  end   <- gpos + cfg$instrument$cis_window_kb * 1000
  v <- paste0(gchr, ":", start, "-", end)
  for (try in 1:3) {
    r <- tryCatch(httr::POST(API, body = list(variant = v, id = OUTCOMES[[on]], proxies = 0),
                             encode = "form", H, httr::timeout(240)),
                  error = function(e) NULL)
    if (is.null(r)) { Sys.sleep(5); next }
    if (r$status_code == 200) {
      arr <- tryCatch(jsonlite::fromJSON(httr::content(r, as = "text", encoding = "UTF-8")),
                      error = function(e) NULL)
      if (!is.null(arr) && is.data.frame(arr) && nrow(arr) > 0) {
        dt <- as.data.table(arr)
        dt[, `:=`(ea = as.character(ea), nea = as.character(nea), rsid = as.character(rsid),
                   eaf = suppressWarnings(as.numeric(eaf)), beta = suppressWarnings(as.numeric(beta)),
                   se = suppressWarnings(as.numeric(se)), p = suppressWarnings(as.numeric(p)),
                   n = suppressWarnings(as.numeric(n)))]
        saveRDS(dt, cf)
        return(dt)
      }
      # 200 但空/异常 → 记空
      saveRDS(data.table(), cf)
      return(data.table())
    }
    Sys.sleep(5 * try)
  }
  stop(paste0("API 失败 ", on, " ", gene))
}
done <- 0
pairs[, `:=`(gwas_ok = FALSE, n_gwas = 0L)]
for (i in seq_len(nrow(pairs))) {
  p <- pairs[i]
  dt <- tryCatch(fetch_region(p$outcome, p$gene, p$GeneChr, p$GenePos), error = function(e) NULL)
  pairs$n_gwas[i] <- if (is.null(dt)) -1 else nrow(dt)
  pairs$gwas_ok[i] <- !is.null(dt)
  done <- done + 1
  if (done %% 25 == 0)
    log("  GWAS 区域 ", done, "/", nrow(pairs), " | 最近 ", p$symbol, "×", p$outcome,
        " n=", pairs$n_gwas[i])
}
saveRDS(pairs[, .(gene, symbol, outcome)], file.path(gdir, "_coloc_pairs_key.rds"))
log("Phase A 完成: ", sum(pairs$gwas_ok), "/", nrow(pairs), " 对取得 GWAS 区域数据")

# --- Phase B：eQTL 侧 cis 全变异（全量文件就绪后执行）----------------------------
if (!file.exists(FULL)) {
  log("全量 cis 文件未就绪，跳过 Phase B/C（Phase A 已缓存，全量下载完成后再跑一次本脚本）")
  quit(save = "no", status = 0)
}
log("读取全量 cis-eQTL 汇总（awk 预过滤至闸门基因）...")
genes_file <- file.path(gdir, "_coloc_genes.txt")
fwrite(unique(data.table(g = c(pairs$gene, pairs$symbol))), genes_file, col.names = FALSE)
cmd <- paste0("zcat ", FULL, " | awk -F'\\t' 'NR==FNR{keep[$1]=1; next} FNR==1 || $8 in keep || $9 in keep {print}' ",
              genes_file, " -")
eqtl_full <- fread(cmd = cmd, sep = "\t", header = TRUE, nThread = 4)
rm_file <- file.exists(file.path(gdir, "_coloc_genes.txt"))
log("全量 cis 子集行: ", nrow(eqtl_full))
eqtl_full <- eqtl_full[SNPChr == GeneChr & abs(SNPPos - GenePos) <= cfg$instrument$cis_window_kb * 1000]
log("cis±1Mb 过滤后: ", nrow(eqtl_full))
af <- fread(cmd = paste0("zcat ", AF), sep = "\t", header = TRUE, nThread = 4)
setnames(af, c("SNP","hg19_chr","hg19_pos","AlleleA","AlleleB","allA_total","allAB_total","allB_total","AlleleB_all"))
eqtl_full <- merge(eqtl_full, af[, .(SNP, AlleleA, AlleleB, AlleleB_all)], by = "SNP", all.x = TRUE)
rm(af)
eqtl_full[, eaf := fcase(AssessedAllele == AlleleA, 1 - AlleleB_all,
                         AssessedAllele == AlleleB, AlleleB_all, default = NA_real_)]
eqtl_full <- eqtl_full[!is.na(eaf)]
eqtl_full[, `:=`(beta = Zscore / sqrt(NrSamples), se = 1 / sqrt(NrSamples))]
eqtl_full <- unique(eqtl_full[, .(Gene, GeneChr, GenePos, SNP, AssessedAllele, OtherAllele, eaf, beta, se, NrSamples)])
log("eQTL 侧 cis 变异（含 eAF）: ", nrow(eqtl_full), " 对基因×SNP")

# --- Phase C：逐对 coloc -------------------------------------------------------
harmonize <- function(e, g) {
  # e = data.table(SNP, AssessedAllele, OtherAllele, eaf, beta, se, NrSamples)
  # g = data.table(rsid, ea, nea, eaf, beta, se, n)  —— 合并前先显式改名，避免同名列冲突
  e2 <- copy(e[, .(SNP, AssessedAllele, OtherAllele, eaf_e = eaf, beta_e = beta,
                   se_e = se, N_e = NrSamples)])
  g2 <- copy(g[, .(rsid = rsid, ea = ea, nea = nea, eaf_g = as.numeric(eaf),
                   beta_g = as.numeric(beta), se_g = as.numeric(se), N_g = as.numeric(n))])
  setkey(g2, rsid); setkey(e2, SNP)
  m <- merge(e2, g2, by.x = "SNP", by.y = "rsid")
  if (nrow(m) == 0) return(data.table())
  A <- m$AssessedAllele; O <- m$OtherAllele; ea <- m$ea; nea <- m$nea
  gb <- m$beta_g
  ok   <- A == ea
  flip <- !ok & A == nea & O == ea
  gb[flip] <- -gb[flip]
  # palindromic 用 eAF 解析：双方 eAF 一侧 >0.58 另一侧 <0.42 → 翻转；其余歧义剔除
  pal <- !ok & !flip & (((A == "A" & O == "T") | (A == "T" & O == "A")) |
                        ((A == "C" & O == "G") | (A == "G" & O == "C")))
  p_flip <- pal & ((m$eaf_e > 0.58 & m$eaf_g < 0.42) | (m$eaf_e < 0.42 & m$eaf_g > 0.58))
  gb[p_flip] <- -gb[p_flip]
  keep <- ok | flip | p_flip
  out <- data.table(snp = m$SNP[keep],
                    maf = pmin(m$eaf_e, 1 - m$eaf_e)[keep],
                    beta_e = m$beta_e[keep], varbeta_e = m$se_e[keep]^2, N_e = m$N_e[keep],
                    beta_g = gb[keep], varbeta_g = m$se_g[keep]^2, N_g = m$N_g[keep])
  out <- out[maf >= 0.01 & maf <= 0.99]
  out <- out[!is.na(beta_e) & !is.na(beta_g) & !is.na(varbeta_e) & !is.na(varbeta_g) &
               varbeta_e > 0 & varbeta_g > 0]
  out
}
run_coloc <- function(gene, on, e, g) {
  if (is.null(g) || nrow(g) == 0) return(list(ok = FALSE, note = "GWAS 区域无数据"))
  if (nrow(e) == 0) return(list(ok = FALSE, note = "eQTL cis 无变异"))
  m <- harmonize(e, g)
  if (nrow(m) < 10) return(list(ok = FALSE, note = paste0("交集/对齐后仅 ", nrow(m), " 变异（<10）")))
  d1 <- list(snp = m$snp, type = "quant", N = m$N_e, beta = m$beta_e, varbeta = m$varbeta_e, MAF = m$maf)
  if (TYPE[[on]] == "cc")
    d2 <- list(snp = m$snp, type = "cc", N = m$N_g, s = S_CASE[[on]], beta = m$beta_g, varbeta = m$varbeta_g, MAF = m$maf)
  else
    d2 <- list(snp = m$snp, type = "quant", N = m$N_g, beta = m$beta_g, varbeta = m$varbeta_g, MAF = m$maf)
  r <- tryCatch(coloc.abf(d1, d2, p12 = 1e-5), error = function(e) NULL)
  if (is.null(r)) return(list(ok = FALSE, note = "coloc.abf 失败"))
  r2 <- tryCatch(coloc.abf(d1, d2, p12 = 1e-6), error = function(e) NULL)
  top_snp <- if (!is.null(r$results) && nrow(r$results)) r$results$snp[which.max(r$results$SNP.PP.H4)] else NA
  list(ok = TRUE, note = "", nsnp = nrow(m), r = r, r_sens = r2, top_snp = top_snp)
}
# coloc.abf 的 summary 键名在不同版本是 "PP.H4" 或 "PP.H4.abf"——统一兼容
ppv <- function(s, nm) {
  for (k in c(nm, paste0(nm, ".abf"))) if (length(s) && !is.na(s[k])) return(unname(as.numeric(s[k])))
  NA_real_
}
log("逐对 coloc（", nrow(pairs), " 对）...")
cols <- data.frame(gene = character(), symbol = character(), outcome = character(),
  stage2_b = numeric(), stage2_pval = numeric(), stage2_q = numeric(),
  n_eqtl_cis = integer(), n_gwas_region = integer(), n_coloc = integer(),
  PP.H0 = numeric(), PP.H1 = numeric(), PP.H2 = numeric(), PP.H3 = numeric(), PP.H4 = numeric(),
  PP.H4_p12e6 = numeric(), top_snp = character(), tier = character(), ok = logical(), note = character())
out_rows <- vector("list", nrow(pairs))
for (i in seq_len(nrow(pairs))) {
  p <- pairs[i]
  g <- tryCatch(readRDS(file.path(gdir, "_coloc_gwas", paste0(OUT_N[[p$outcome]], "_", p$gene, ".rds"))),
                error = function(e) data.table())
  e <- eqtl_full[Gene == p$gene]
  cl <- run_coloc(p$gene, p$outcome, e, g)
  if (cl$ok) {
    pp <- cl$r$summary
    pp4_s <- if (!is.null(cl$r_sens)) ppv(cl$r_sens$summary, "PP.H4") else NA
    pp4 <- ppv(pp, "PP.H4")
    tier <- ifelse(pp4 >= 0.8, "strong", ifelse(pp4 >= 0.5, "moderate", "none"))
    out_rows[[i]] <- data.frame(gene = p$gene, symbol = p$symbol, outcome = p$outcome,
      stage2_b = p$b, stage2_pval = p$pval, stage2_q = p$q_stage2,
      n_eqtl_cis = nrow(e), n_gwas_region = nrow(g), n_coloc = cl$nsnp,
      PP.H0 = ppv(pp, "PP.H0"), PP.H1 = ppv(pp, "PP.H1"), PP.H2 = ppv(pp, "PP.H2"),
      PP.H3 = ppv(pp, "PP.H3"), PP.H4 = pp4, PP.H4_p12e6 = pp4_s,
      top_snp = cl$top_snp, tier = tier, ok = TRUE, note = "")
  } else {
    out_rows[[i]] <- data.frame(gene = p$gene, symbol = p$symbol, outcome = p$outcome,
      stage2_b = p$b, stage2_pval = p$pval, stage2_q = p$q_stage2,
      n_eqtl_cis = nrow(e), n_gwas_region = nrow(g), n_coloc = NA,
      PP.H0 = NA, PP.H1 = NA, PP.H2 = NA, PP.H3 = NA, PP.H4 = NA,
      PP.H4_p12e6 = NA, top_snp = NA, tier = NA, ok = FALSE, note = cl$note)
  }
  if (i %% 25 == 0) log("  coloc ", i, "/", nrow(pairs), " | ", p$symbol, "×", p$outcome,
                        " PP.H4=", if (cl$ok) round(ppv(cl$r$summary, "PP.H4"), 3) else cl$note)
}
co <- rbindlist(out_rows, fill = TRUE)
write.csv(co, file.path(gdir, "transcript_coloc.csv"), row.names = FALSE)
hits <- co[ok == TRUE & tier == "strong"]
write.csv(hits, file.path(gdir, "transcript_coloc_hits.csv"), row.names = FALSE)
scan_hits <- nrow(fread(file.path(gdir, "transcript_grid_hits.csv")))
funnel <- data.frame(stage = c("scan_fdr05_hits", "stage2_nominal", "coloc_attempted",
                               "coloc_ok_nsnp10", "pp4_strong", "pp4_moderate", "pp4_none"),
  count = c(scan_hits, nrow(pairs), nrow(pairs), sum(co$ok), sum(co$tier == "strong", na.rm = TRUE),
            sum(co$tier == "moderate", na.rm = TRUE), sum(co$tier == "none", na.rm = TRUE)))
write.table(funnel, file.path(gdir, "transcript_coloc_funnel.tsv"), sep = "\t", row.names = FALSE)
log("=== M5 coloc 完成 ✔ | PP.H4≥0.8: ", sum(co$tier == "strong", na.rm = TRUE),
    " | 0.5–0.8: ", sum(co$tier == "moderate", na.rm = TRUE),
    " | <0.5: ", sum(co$tier == "none", na.rm = TRUE), " ===")
