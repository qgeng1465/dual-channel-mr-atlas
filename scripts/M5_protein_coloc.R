#!/usr/bin/env Rscript
# =============================================================================
# M5_protein_coloc.R — 蛋白通道 共定位闸门（coloc.abf，PP.H4 三档）
# =============================================================================
# 目的：对 deCODE 血浆 pQTL（cis 窗内全部变异）× OpenGWAS 结局做共定位，
#   回答"蛋白通道 MR 信号是否由共享因果变异驱动（PP.H4）而非 LD 混杂"。
#   与 M5_transcript_coloc.R（转录通道）同口径、同先验，结果可比。
# 数据：
#   - pQTL 侧：deCODE cis 子集（data/decode/sub/*_cis.txt.gz，M1 裁剪，hg38
#     TSS±1Mb）。列：Chrom Pos Name rsids effectAllele otherAllele Beta Pval
#     minus_log10_pval SE N ImpMAF。beta=Beta(per-SD)，se=SE，maf=ImpMAF，
#     type=quant，N=变异级中位数。
#   - GWAS 侧：OpenGWAS associations API 区域查询（variant=chr:start-end，
#     proxies=0），缓存 _coloc_gwas_prot/。基因 hg19 坐标 = ENSEMBL GRCh37 REST
#     实查（2026-08-13），与 deCODE hg38 cis 窗通过 rsid 交集对齐。
# 设计（对齐预注册 §3 + README §0.4 M5，转录/蛋白同口径）：
#   - 范围 = 蛋白通道 MR 有工具的对（protein_decode_mr.csv ok=TRUE），全网格含空
#   - 区域 = 基因 hg19 TSS ± 1Mb（与 cis 窗一致）
#   - coloc.abf 默认先验 p1=1e-4 p2=1e-4 p12=1e-5；敏感性 p12=1e-6
#   - type：pQTL=quant（N=deCODE 变异级中位数）；T2D/CAD=cc（s=gwasinfo 病例比例）；
#     FBG=quant（N=n）
#   - PP.H4 三档：≥0.8 强共定位 / 0.5–0.8 中等 / <0.5 无证据；连续 PP.H4 全量报告
#   - palindromic 处理（诚实声明）：deCODE 未提供 effectAlleleFreq（annotated 未就绪），
#     ImpMAF 不总是 effect allele 频率 → 回文位点一律保守排除（与 M3_protein_decode.R
#     的 MR 同口径）；非回文位点按等位基因匹配/翻转对齐。
# 注意：样本重叠（deCODE 与 GWAS 部分队列重叠）为 coloc 已知局限，输出注明不作修正。
# 用法：PATH=$R_ENV/bin:$PATH Rscript scripts/M5_protein_coloc.R
# 产物：results/grid/protein_coloc.csv + _hits.csv + _funnel.tsv
# 缓存：results/grid/_coloc_gwas_prot/<outcome>_<gene>.rds（逐对）
# =============================================================================
suppressMessages({
  library(data.table)
  suppressPackageStartupMessages(library(coloc))
  library(jsonlite); library(httr)
})
proj <- "/data/qiushuogeng/projects/dual-channel-mr-atlas"
res  <- file.path(proj, "results")
gdir <- file.path(res, "grid")
SUB  <- file.path(proj, "data/decode/sub")
dir.create(file.path(gdir, "_coloc_gwas_prot"), recursive = TRUE, showWarnings = FALSE)
log <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n", sep = "")

# --- 预注册完整性校验 ----------------------------------------------------------
prereg <- file.path(proj, "docs/PREREGISTRATION.md")
lock   <- file.path(proj, "docs/PREREGISTRATION.md.sha256")
stopifnot(file.exists(prereg), file.exists(lock),
          tools::md5sum(prereg) == readLines(lock))
cfg <- read_json(file.path(res, "config.json"), auto_unbox = TRUE)
OUTCOMES <- list(t2d = "ebi-a-GCST006867", cad = "ebi-a-GCST005194", fbg = "ebi-a-GCST005186")
OUT_N <- c(t2d = "GCST006867", cad = "GCST005194", fbg = "GCST005186")
# 2026-08-13 P5：pairs$outcome 是全文 id（ebi-a-GCST006867），而 OUTCOMES/OUT_N 以短名索引。
#   原代码直接 OUTCOMES[[on]] → NULL → API id=NULL 全部失败（n=-1）+ 缓存文件名错。
ID2SHORT <- setNames(names(OUTCOMES), unname(unlist(OUTCOMES)))  # 全 id → 短名
TYPE <- c(t2d = "cc", cad = "cc", fbg = "quant")
S_CASE <- c(t2d = 61714 / 655666, cad = 34541 / 296525)
GWAS_N <- c(t2d = 655666, cad = 296525, fbg = 58074)
log("预注册哈希校验通过 ✔ | M5 蛋白共定位闸门 | deCODE pQTL × OpenGWAS | PP.H4 三档")

JWT <- Sys.getenv("OPENGWAS_JWT")
if (!nzchar(JWT)) stop("OPENGWAS_JWT 未设置（~/.Renviron）")
H <- httr::add_headers(Authorization = paste("Bearer", JWT))
API <- "https://api.opengwas.io/api/associations"

# --- 蛋白文件 + hg38/hg19 坐标 --------------------------------------------------
# hg38 TSS 与 M3 一致（strand− 取 end）；hg19 TSS = ENSEMBL GRCh37 REST 实查（2026-08-13）
PROTEINS <- list(
  list(file = "5231_79_PCSK9_PCSK9.txt.gz",          gene = "PCSK9",   chr = "1",  hg38 = 55039445, hg19 = 55505221, note = "positive control"),
  list(file = "5230_99_HMGCR_HMGR.txt.gz",           gene = "HMGCR",   chr = "5",  hg38 = 75336329, hg19 = 74632154, note = "negative control"),
  list(file = "10391_1_ANGPTL3_ANGL3.txt.gz",        gene = "ANGPTL3", chr = "1",  hg38 = 62597464, hg19 = 63063158, note = "negative control"),
  list(file = "6461_54_APOC3_Apo_C_III.txt.gz",      gene = "APOC3",   chr = "11", hg38 = 116827019, hg19 = 116700422, note = "negative control"),
  list(file = "2797_56_APOB_Apo_B.txt.gz",           gene = "APOB",    chr = "2",  hg38 = 21044075, hg19 = 21266945, note = "optional"),
  list(file = "13129_40_LDLR_LDLR.txt.gz",           gene = "LDLR",    chr = "19", hg38 = 11089418, hg19 = 11200038, note = "optional"),
  list(file = "13085_18_GLP1R_GLP1R.txt.gz",         gene = "GLP1R",   chr = "6",  hg38 = 39048562, hg19 = 39016574, note = "drug target GLP-1RA"),
  list(file = "15460_9_DPP4_CD26.txt.gz",            gene = "DPP4",    chr = "2",  hg38 = 162074639, hg19 = 162931052, note = "drug target DPP4i"),
  list(file = "3448_13_INSR_IR.txt.gz",              gene = "INSR",    chr = "19", hg38 = 7294443, hg19 = 7294045, note = "drug target insulin/insulin sens."),
  list(file = "18182_24_PCK1_PCKGC.txt.gz",          gene = "PCK1",    chr = "20", hg38 = 57546220, hg19 = 56136136, note = "drug target gluconeogenesis"),
  list(file = "4891_50_GCG_Glucagon.txt.gz",         gene = "GCG",     chr = "2",  hg38 = 162152404, hg19 = 163008914, note = "drug target glucagon")
)

# --- 载入蛋白通道 MR 结果（闸门范围 = 有工具的蛋白×结局对）----------------------
mr <- fread(file.path(gdir, "protein_decode_mr.csv"))
mr_prim <- mr[method == "Inverse variance weighted (multiplicative random effects)" |
                method == "Wald ratio"]               # 主方法行（含 nsnp=1 Wald）
pairs <- mr_prim[ok == TRUE, .(gene, outcome, nsnp, b, se, pval)]
pairs[, on := ID2SHORT[outcome]]                      # 短名，供 OUTCOMES/OUT_N 索引
pairs <- merge(pairs, rbindlist(lapply(PROTEINS, as.data.table)), by = "gene", all.x = TRUE)
pairs <- pairs[!is.na(hg19)]
log("闸门对: ", nrow(pairs), "（蛋白 MR 有工具；基因 ", uniqueN(pairs$gene),
    "；结局 x", length(OUTCOMES), "）")

# --- Phase A：GWAS 区域数据（逐对 range query，缓存，可断点续跑）----------------
fetch_region <- function(on, gene, gchr, gpos) {
  cf <- file.path(gdir, "_coloc_gwas_prot", paste0(OUT_N[[on]], "_", gene, ".rds"))
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
      saveRDS(data.table(), cf)                     # 200 但空 → 记空
      return(data.table())
    }
    Sys.sleep(5 * try)
  }
  stop(paste0("API 失败 ", on, " ", gene))
}
for (i in seq_len(nrow(pairs))) {
  p <- pairs[i]
  dt <- tryCatch(fetch_region(p$on, p$gene, p$chr, p$hg19), error = function(e) NULL)
  pairs$gwas_n[i] <- if (is.null(dt)) -1 else nrow(dt)
  if (i %% 3 == 0) log("  GWAS 区域 ", i, "/", nrow(pairs), " | ", p$gene, "×", p$outcome,
                       " n=", pairs$gwas_n[i])
}
saveRDS(pairs[, .(gene, chr, hg19, outcome)], file.path(gdir, "_coloc_prot_pairs_key.rds"))
log("Phase A 完成: 全 ", nrow(pairs), " 对取得 GWAS 区域数据（含空）")

# --- Phase B：pQTL 侧 cis 全变异 + eAF（deCODE）---------------------------------
load_pqtl <- function(prot) {
  # 2026-08-13 P4：M1 命名 = ${base%.gz}_cis.txt.gz（去 .gz 再加 _cis.txt.gz），
  #   原 paste0(file,"_cis.txt.gz") 拼错文件名 → 全部蛋白 0 变异。
  f <- file.path(SUB, paste0(sub("\\.gz$", "", prot$file), "_cis.txt.gz"))
  if (!file.exists(f)) return(data.table())
  d <- fread(cmd = paste0("zcat ", shQuote(f)), sep = "\t", header = TRUE, nThread = 4)
  d[, Chrom2 := as.character(gsub("^chr", "", Chrom))]
  d[, Pos := as.numeric(Pos)]
  # 排除 readme 注明的多等位 bug 行（与 M3 同口径）
  d <- d[effectAllele != otherAllele & !otherAllele %in% c("!", "", NA)]
  # rsid 解析：rsids 列可多 rsid 逗号分隔；NA 回落 Name（chr_pos_ref_alt）
  d[, rsid2 := fifelse(is.na(rsids) | rsids %in% c("", ".", "-"),
                       Name, sub("[,; ].*$", "", rsids))]
  d <- d[!duplicated(rsid2)]
  d[, `:=`(beta = as.numeric(Beta), varbeta = as.numeric(SE)^2,
           maf = as.numeric(ImpMAF), Np = as.numeric(N))]
  # deCODE readme：ImpMAF 即次等位频率（不总是 effect allele 频率），直接作 coloc MAF。
  # 原 pmin(1-ImpMAF, ImpMAF) 兜底对 NA 无效（NA 入 NA 出）且对非 NA 恒等 → 移除，NA 由下方过滤剔除
  d <- d[!is.na(beta) & !is.na(varbeta) & !is.na(maf) & varbeta > 0 &
           maf >= 0.01 & maf <= 0.99]
  d
}
pqtl_list <- lapply(PROTEINS, function(p) {
  d <- load_pqtl(p)
  list(gene = p$gene, d = d)
})
names(pqtl_list) <- sapply(PROTEINS, `[[`, "gene")
log("pQTL 侧 cis 变异载入: ", paste(sapply(pqtl_list, function(x) paste0(x$gene, "=", nrow(x$d))),
                                    collapse = " "))

# --- Phase C：逐对 coloc --------------------------------------------------------
is_pal <- function(a, b) (a == "A" & b == "T") | (a == "T" & b == "A") |
                          (a == "C" & b == "G") | (a == "G" & b == "C")
harmonize <- function(e, g) {
  # e = deCODE pQTL（snp=rsid2, A=effectAllele, O=otherAllele, beta, varbeta, maf, Np）
  # g = GWAS region（rsid, ea, nea, eaf, beta, se, n）
  e2 <- copy(e[, .(snp = rsid2, A = effectAllele, O = otherAllele, beta_e = beta,
                   varbeta_e = varbeta, maf_e = maf, N_e = Np)])
  g2 <- copy(g[, .(snp = rsid, ea = ea, nea = nea, eaf_g = as.numeric(eaf),
                   beta_g = as.numeric(beta), se_g = as.numeric(se), N_g = as.numeric(n))])
  setkey(g2, snp); setkey(e2, snp)
  m <- merge(e2, g2, by = "snp")
  if (nrow(m) == 0) return(data.table())
  A <- m$A; O <- m$O; ea <- m$ea; nea <- m$nea
  gb <- m$beta_g
  ok   <- A == ea
  flip <- !ok & A == nea & O == ea
  gb[flip] <- -gb[flip]
  pal <- !ok & !flip & is_pal(A, O)
  # palindromic 无可靠 effectAlleleFreq（ImpMAF 不总是 eaf）→ 保守排除（诚实声明）
  keep <- (ok | flip) & !pal
  out <- data.table(snp = m$snp[keep], maf = m$maf_e[keep],
                    beta_e = m$beta_e[keep], varbeta_e = m$varbeta_e[keep], N_e = m$N_e[keep],
                    beta_g = gb[keep], varbeta_g = m$se_g[keep]^2, N_g = m$N_g[keep])
  out <- out[maf >= 0.01 & maf <= 0.99 & N_e > 0 & N_g > 0 &
               !is.na(beta_e) & !is.na(beta_g) & !is.na(varbeta_e) & !is.na(varbeta_g) &
               varbeta_e > 0 & varbeta_g > 0]
  out
}
run_coloc <- function(gene, on, e, g) {
  if (is.null(g) || nrow(g) == 0) return(list(ok = FALSE, note = "GWAS 区域无数据"))
  if (nrow(e) == 0) return(list(ok = FALSE, note = "pQTL cis 无变异"))
  m <- harmonize(e, g)
  if (nrow(m) < 10) return(list(ok = FALSE, note = paste0("交集/对齐后仅 ", nrow(m), " 变异（<10）")))
  d1 <- list(snp = m$snp, type = "quant", N = round(median(m$N_e)),
             beta = m$beta_e, varbeta = m$varbeta_e, MAF = m$maf)
  if (TYPE[[on]] == "cc")
    d2 <- list(snp = m$snp, type = "cc", N = m$N_g, s = S_CASE[[on]],
               beta = m$beta_g, varbeta = m$varbeta_g, MAF = m$maf)
  else
    d2 <- list(snp = m$snp, type = "quant", N = m$N_g,
               beta = m$beta_g, varbeta = m$varbeta_g, MAF = m$maf)
  r <- tryCatch(coloc.abf(d1, d2, p12 = 1e-5), error = function(e) NULL)
  if (is.null(r)) return(list(ok = FALSE, note = "coloc.abf 失败"))
  r2 <- tryCatch(coloc.abf(d1, d2, p12 = 1e-6), error = function(e) NULL)
  top_snp <- if (!is.null(r$results) && nrow(r$results)) r$results$snp[which.max(r$results$SNP.PP.H4)] else NA
  list(ok = TRUE, note = "", nsnp = nrow(m), r = r, r_sens = r2, top_snp = top_snp)
}
ppv <- function(s, nm) {
  for (k in c(nm, paste0(nm, ".abf"))) if (length(s) && !is.na(s[k])) return(unname(as.numeric(s[k])))
  NA_real_
}
log("逐对 coloc（", nrow(pairs), " 对）...")
out_rows <- vector("list", nrow(pairs))
for (i in seq_len(nrow(pairs))) {
  p <- pairs[i]
  g <- tryCatch(readRDS(file.path(gdir, "_coloc_gwas_prot",
                                  paste0(OUT_N[[p$on]], "_", p$gene, ".rds"))),
                error = function(e) data.table())
  e <- pqtl_list[[p$gene]]$d
  cl <- run_coloc(p$gene, p$on, e, g)   # P6: TYPE/S_CASE 以短名索引，传 p$on
  if (cl$ok) {
    pp <- cl$r$summary
    pp4 <- ppv(pp, "PP.H4")
    pp4_s <- if (!is.null(cl$r_sens)) ppv(cl$r_sens$summary, "PP.H4") else NA
    tier <- ifelse(pp4 >= 0.8, "strong", ifelse(pp4 >= 0.5, "moderate", "none"))
    out_rows[[i]] <- data.frame(gene = p$gene, outcome = p$outcome,
      mr_nsnp = p$nsnp, mr_b = p$b, mr_pval = p$pval,
      n_pqtl_cis = nrow(e), n_gwas_region = nrow(g), n_coloc = cl$nsnp,
      PP.H0 = ppv(pp, "PP.H0"), PP.H1 = ppv(pp, "PP.H1"), PP.H2 = ppv(pp, "PP.H2"),
      PP.H3 = ppv(pp, "PP.H3"), PP.H4 = pp4, PP.H4_p12e6 = pp4_s,
      top_snp = cl$top_snp, tier = tier, ok = TRUE, note = "")
  } else {
    out_rows[[i]] <- data.frame(gene = p$gene, outcome = p$outcome,
      mr_nsnp = p$nsnp, mr_b = p$b, mr_pval = p$pval,
      n_pqtl_cis = nrow(e), n_gwas_region = nrow(g), n_coloc = NA,
      PP.H0 = NA, PP.H1 = NA, PP.H2 = NA, PP.H3 = NA, PP.H4 = NA,
      PP.H4_p12e6 = NA, top_snp = NA, tier = NA, ok = FALSE, note = cl$note)
  }
  if (i %% 3 == 0 || i == nrow(pairs))
    log("  coloc ", i, "/", nrow(pairs), " | ", p$gene, "×", p$outcome,
        " PP.H4=", if (cl$ok) round(ppv(cl$r$summary, "PP.H4"), 3) else cl$note)
}
co <- rbindlist(out_rows, fill = TRUE)
write.csv(co, file.path(gdir, "protein_coloc.csv"), row.names = FALSE)
hits <- co[ok == TRUE & tier == "strong"]
write.csv(hits, file.path(gdir, "protein_coloc_hits.csv"), row.names = FALSE)
funnel <- data.frame(stage = c("mr_ok_pairs", "coloc_attempted", "coloc_ok_nsnp10",
                               "pp4_strong", "pp4_moderate", "pp4_none"),
  count = c(nrow(pairs), nrow(pairs), sum(co$ok), sum(co$tier == "strong", na.rm = TRUE),
            sum(co$tier == "moderate", na.rm = TRUE), sum(co$tier == "none", na.rm = TRUE)))
write.table(funnel, file.path(gdir, "protein_coloc_funnel.tsv"), sep = "\t", row.names = FALSE)
log("=== M5 蛋白 coloc 完成 ✔ | PP.H4≥0.8: ", sum(co$tier == "strong", na.rm = TRUE),
    " | 0.5–0.8: ", sum(co$tier == "moderate", na.rm = TRUE),
    " | <0.5: ", sum(co$tier == "none", na.rm = TRUE),
    " | 失败: ", sum(co$ok == FALSE, na.rm = TRUE), " ===")
log("注意：deCODE ImpMAF 无 effectAlleleFreq → 回文位点保守排除（如实报告）")
