#!/usr/bin/env Rscript
# =============================================================================
# M20_feasibility_pilot.R — 破局方案可行性验证：MR 显著集之外 coloc 分布 + 分层
# =============================================================================
# 依据 2026-08-15 对抗性评审（docs/archive/verification_20260815.md §1/§7）修订：
#   - 400 对 → 3000 灰区 + 3000 阴性（按结局分层），紧化 q 估计
#   - 关键分层：coloc-only 命中按 GWAS 区域峰 p<5e-8 拆为
#       "真共定位（MR 功率问题）" vs "coloc 伪阳（GWAS 侧不显著）"
#   - 顺带按 eQTL 强度（|beta|/se）分箱 → 区分"弱 eQTL 功效陈述"与"工具选择问题"
# 判据：若 coloc-only 中 GWAS 峰显著占多数 → "MR 召回低"叙事成立；
#        若 mostly GWAS 峰不显著 → 召回叙事是 coloc 伪阳，转向
# 数据：data/opengwas/full/*.gz（M19）+ eQTLGen 全量；输出 results/feasibility_pilot_20260815.csv
# 用法：PATH=$R_ENV/bin:$PATH Rscript scripts/M20_feasibility_pilot.R
# =============================================================================
suppressMessages({
  library(data.table)
  suppressPackageStartupMessages(library(coloc))
})
proj <- "<repo-root>"
res  <- file.path(proj, "results")
gdir <- file.path(res, "grid")
GWAS <- file.path(proj, "data/opengwas/full")
FULL <- file.path(proj, "data/eqtlgen/cis-eQTLs_full_20180905.txt.gz")
SIG  <- file.path(proj, "data/eqtlgen/cis-eQTL-significant.txt.gz")
AF   <- file.path(proj, "data/eqtlgen/SNP_AF.txt.gz")
log <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n")
stopifnot(file.exists(file.path(GWAS, "DONE")))

TYPE   <- c(t2d = "cc", cad = "cc", fbg = "quant")
S_CASE <- c(t2d = 61714/655666, cad = 34541/296525, fbg = NA)
GWAS_N <- c(t2d = 655666, cad = 296525, fbg = 58074)
WINDOW <- 1e6
N_PER  <- as.integer(Sys.getenv("M20_NPER", "1000"))   # 每结局×每分组抽样数（6K 对；冒烟测试可覆盖）

# ---- 1. 试点抽样：灰区（0.05≤p<0.5）+ 阴性（p≥0.5），按结局分层 ----
grid <- fread(file.path(gdir, "transcript_grid_mr.csv"))
grid <- grid[ok == TRUE & !is.na(pval)]
grid[, grp := fifelse(pval < 0.05, "sig", fifelse(pval < 0.5, "grey", "null"))]
set.seed(42)
sel <- rbind(
  grid[grp == "grey"][, .SD[sample(.N, min(.N, N_PER))], by = outcome],
  grid[grp == "null"][, .SD[sample(.N, min(.N, N_PER))], by = outcome]
)
cat("试点对：灰区", nrow(sel[grp=="grey"]), " 阴性", nrow(sel[grp=="null"]), "\n")

# ---- 2. 基因坐标 + eQTL cis 子集（awk 落盘 /data，绕开根盘管道临时文件） ----
# 坑（2026-08-15 实测）：fread(cmd=) 读管道需临时文件，R 用 tempdir() 落 /tmp（根盘 50G 满）
# → awk 扫完 153M 行后 fread 写 /tmp/Rtmp* 时 No space left on device。改两步：
#   ① awk 抽取瘦列写 gz 落盘到 /data；② fread 直接读文件（无管道临时文件）。
gcoord <- unique(fread(cmd = paste0("zcat ", SIG), sep = "\t", header = TRUE, nThread = 4,
                       select = c("Gene", "GeneChr", "GenePos")))
sel <- merge(sel, gcoord, by.x = "gene", by.y = "Gene", all.x = TRUE)[!is.na(GeneChr)]
genes_file <- tempfile(); fwrite(unique(data.table(g = c(sel$gene, sel$symbol))), genes_file, col.names = FALSE)
# eQTLGen full 列序（2026-08-15 逐列核实）：$8=Gene $9=GeneSymbol（keep 含 ENSG+symbol，双列匹配）。
# 用 fread(cmd=) 单引号内联 awk（$ 不经 shell 展开，语法正确）+ select 按列名（不按位置）。
# 历史坑记录：① system2("bash", c("-c",cmd)) 内联 awk 的 $ 会被 bash 当位置参数展开成空
#  → awk 错位/语法错，gzip "unexpected end of file" 落盘 0MB；② awk 位置抽列错位
#  （$7 当 Zscore 实为 OtherAllele）→ 必须按列名；③ fread 从管道读大输入需临时文件
#  → 必须 TMPDIR=/data（根盘 50G 满，/tmp 曾 No space left on device）。
eqtl <- fread(cmd = paste0("zcat ", FULL, " | awk -F'\\t' 'NR==FNR{keep[$1]=1; next} FNR==1 || $8 in keep || $9 in keep {print}' ",
                           genes_file, " -"),
              sep = "\t", header = TRUE, nThread = 4,
              select = c("SNP", "SNPChr", "SNPPos", "Zscore", "AssessedAllele", "OtherAllele",
                         "Gene", "GeneChr", "GenePos", "NrSamples"))
eqtl <- eqtl[SNPChr == GeneChr & abs(SNPPos - GenePos) <= WINDOW][, .(SNP, Zscore, AssessedAllele, OtherAllele, Gene, NrSamples)]
af <- fread(cmd = paste0("zcat ", AF), sep = "\t", header = TRUE, nThread = 4,
            select = c("SNP", "AlleleA", "AlleleB", "AlleleB_all"))
eqtl <- merge(eqtl, af, by = "SNP", all.x = TRUE)
eqtl[, eaf := fcase(AssessedAllele == AlleleA, 1 - AlleleB_all,
                    AssessedAllele == AlleleB, AlleleB_all, default = NA_real_)][!is.na(eaf)]
eqtl[, `:=`(beta = Zscore / sqrt(NrSamples), se = 1 / sqrt(NrSamples))]
eqtl <- unique(eqtl[, .(Gene, SNP, AssessedAllele, OtherAllele, eaf, beta, se, NrSamples)])
setkey(eqtl, Gene)   # 6K 次按 Gene 查找必须走键二分，否则扫描 50M 行太慢
gc()
log("eQTL cis 变异（pilot 基因）: ", nrow(eqtl), " | 唯一基因: ", length(unique(eqtl$Gene)))

# ---- 3. GWAS 侧：加载到内存（fread 流式 zcat，瘦列） ----
load_gwas <- function(on) {
  f <- file.path(GWAS, paste0(on, "_full.gz"))
  if (!file.exists(f)) return(NULL)
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
  # 注：GWAS Catalog harmonised 的 hm_ci_lower/upper 全空；se 用原始 standard_error
  #（se 等位取向不变，harmonise 只翻 beta 符号，不翻 se）
  t[, `:=`(rsid = as.character(hm_rsid), ea = as.character(hm_effect_allele),
           nea = as.character(hm_other_allele), eaf = as.numeric(hm_effect_allele_frequency),
           beta = as.numeric(hm_beta), se = as.numeric(standard_error),
           p = as.numeric(p_value), chr = as.integer(hm_chrom), pos = as.integer(hm_pos))]
  t <- t[!is.na(beta) & !is.na(se) & se > 0 & !is.na(p) & !is.na(chr) & !is.na(pos) & !is.na(rsid)]
  setkey(t, chr)
  t[, .(rsid, ea, nea, eaf, beta, se, p, chr, pos)]
}
gwas <- list()
for (on in c("t2d", "cad", "fbg")) { gwas[[on]] <- load_gwas(on); log(on, " 行: ", nrow(gwas[[on]])) }

# ---- 4. harmonise + coloc（复用 M5 逻辑，harmonise 不变）----
harmonize <- function(e, g) {
  e2 <- copy(e[, .(SNP, AssessedAllele, OtherAllele, eaf_e = eaf, beta_e = beta,
                   se_e = se, N_e = NrSamples)])
  g2 <- copy(g[, .(rsid, ea, nea, eaf_g = as.numeric(eaf),
                   beta_g = as.numeric(beta), se_g = as.numeric(se))])
  setkey(g2, rsid); setkey(e2, SNP)
  m <- merge(e2, g2, by.x = "SNP", by.y = "rsid")
  if (nrow(m) == 0) return(data.table())
  A <- toupper(m$AssessedAllele); O <- toupper(m$OtherAllele)
  ea <- toupper(m$ea); nea <- toupper(m$nea)   # FBG 等位基因小写，统一大小写
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
run_coloc <- function(on, e, g) {
  if (is.null(g) || nrow(g) == 0) return(list(ok = FALSE, note = "GWAS 区域无数据"))
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

out_rows <- vector("list", nrow(sel))
for (i in seq_len(nrow(sel))) {
  pr <- sel[i]
  e <- eqtl[Gene == pr$gene]
  # 2026-08-15 修复（build 错位）：OpenGWAS full 的 hm_pos 是 hg38，eQTLGen 是 hg19，
  #   位置窗口提取错位致 rsid 交集骤减。全结局统一纯 rsid 匹配（与 fbg 原逻辑一致）。
  g <- gwas[[pr$outcome]][rsid %in% e$SNP]
  g_min_p <- if (nrow(g)) min(g$p, na.rm = TRUE) else NA
  cl <- run_coloc(pr$outcome, e, g)
  # eQTL 强度：区域最强 cis-eQTL 的 |beta|/se（lead 功效代理）
  e_fmax <- if (nrow(e)) max(abs(e$beta) / e$se, na.rm = TRUE) else NA
  out_rows[[i]] <- data.table(gene = pr$gene, symbol = pr$symbol, outcome = pr$outcome,
    grp = pr$grp, mr_b = pr$b, mr_p = pr$pval,
    gwas_min_p = g_min_p, eqtl_F_max = e_fmax,
    nsnp = if (cl$ok) cl$nsnp else NA, pp4 = if (cl$ok) cl$pp4 else NA,
    ok = cl$ok, note = cl$note)
  if (i %% 500 == 0) { log("  试点 ", i, "/", nrow(sel)); gc() }
}
res_df <- rbindlist(out_rows, fill = TRUE)
fwrite(res_df, file.path(res, "feasibility_pilot_20260815.csv"))
log("已写 results/feasibility_pilot_20260815.csv | ", nrow(res_df), " 行")

# ---- 5. 分层摘要（对抗性评审 §1/§7 要求）----
cat("\n=== 试点 coloc 摘要（MR 显著集之外）===\n")
for (g0 in c("grey", "null")) {
  d <- res_df[grp == g0]
  n_qc <- d[ok == TRUE]
  hit  <- n_qc[pp4 >= 0.8]
  cat(g0, ": ", nrow(d), " 对 | QC 通过 ", nrow(n_qc), " (", round(100*nrow(n_qc)/nrow(d),1), "%)",
      " | PP.H4≥0.8 ", nrow(hit), " (", round(100*nrow(hit)/max(1, nrow(n_qc)),2), "% of QC, ",
      round(100*nrow(hit)/nrow(d),2), "% of 全部)\n", sep="")
  if (nrow(hit)) {
    # 关键分层：coloc-only 命中里 GWAS 峰是否显著
    gsig <- sum(hit$gwas_min_p < 5e-8, na.rm = TRUE)
    cat("   → coloc-only 中 GWAS 峰 p<5e-8（真共定位，MR 功率问题）: ", gsig, "/", nrow(hit),
        " (", round(100*gsig/nrow(hit),1), "%)\n", sep="")
    cat("   → GWAS 峰不显著（coloc 伪阳候选）: ", nrow(hit) - gsig, "/", nrow(hit), "\n", sep="")
    cat("   → 按 eQTL 强度中位数分箱: ",
        nrow(hit[eqtl_F_max >= median(res_df$eqtl_F_max, na.rm=TRUE)]), " 强 eQTL / ",
        nrow(hit[eqtl_F_max < median(res_df$eqtl_F_max, na.rm=TRUE)]), " 弱 eQTL\n", sep="")
  }
  if (nrow(n_qc)) cat("   PP.H4≥0.5: ", sum(n_qc$pp4>=0.5), " | ≥0.6: ", sum(n_qc$pp4>=0.6),
                      " | ≥0.9: ", sum(n_qc$pp4>=0.9), "\n", sep="")
}
cat("coloc QC（nsnp≥10）总通过率: ", round(100*sum(res_df$ok)/nrow(res_df),1), "%\n")
cat("\n判定：若 coloc-only 中 GWAS 峰显著占比高（>50%）→ 'MR 召回低'叙事成立；\n")
cat("      若 GWAS 峰不显著占多数 → 召回叙事多为 coloc 伪阳，需转向。\n")
