#!/usr/bin/env Rscript
# =============================================================================
# M5_protein_coloc_susie.R — 蛋白通道 coloc-SuSiE 复核（多因果变异，探索性稳健分析）
# =============================================================================
# 目的：coloc.abf 假设每个性状单一因果变异。APOC3（chr11 APOA1/C3/A4/A5 脂质基因簇，
#   top rs964184 典型多信号区）与 INSR（chr19 多信号）最可能违背单信号假设。
#   对 3 个 "MR 显著但 coloc.abf none" 对 + 2 个 strong 对跑 coloc-SuSiE，
#   SuSiE PP.H4 与 abf PP.H4 双列诚实对照（不替代主闸门，探索性稳健分析）。
#
# 2026-08-13 三修（收敛 + 因果变异保留 + 输出规范化）：
#   (1) susie_rss 完整 cis LD（~2000-4000 变异）estimate_prior_variance 不收敛。
#       修复：前 300 信号变异 + EPV=FALSE + max_iter=10000 直接调用（绕过 runsusie ×100 跳档）。
#   (2) 发现 top-300 按 min(pqtl,gwas) z² 截断会丢弃"一侧极强另一侧中等"的因果变异
#       （PCSK9 rs11591147 R46L：pQTL z=23 但 GWAS z=9.7 → min-z² 排名被压出前 300）。
#       改为 方案B：top-300 by min-z² ∪ 各性状单侧 top-20（保证双侧最强信号必在输入）。
#   (3) coloc.susie 输出多行（每 CS 对一行），旧版只存 PP.H4 向量导致 CSV 混乱。
#       新版输出：每对单行汇总（max PP.H4_susie + 对应 hit1/hit2）+ 详细行表。
#
# 数据：复用 M5_protein_coloc.R 的 pQTL cis + GWAS 区域缓存 + harmonize。
# LD：1000G Phase3 EUR（data/ldref/1kg.v3/EUR）plink --r square，交集变异。
# 用法：PATH=$R_ENV/bin:$PATH Rscript scripts/M5_protein_coloc_susie.R
# 产物：results/grid/protein_coloc_susie_summary.csv（每对 1 行）
#       results/grid/protein_coloc_susie_detail.csv（每 CS 对 1 行，审稿人可核）
# =============================================================================
suppressMessages({
  library(data.table)
  suppressPackageStartupMessages(library(coloc))
  suppressPackageStartupMessages(library(susieR))
})
proj <- "/data/qiushuogeng/projects/dual-channel-mr-atlas"
res  <- file.path(proj, "results")
gdir <- file.path(res, "grid")
SUB  <- file.path(proj, "data/decode/sub")
PLINK <- file.path(proj, "tools/plink")
LDREF <- file.path(proj, "data/ldref/1kg.v3/EUR")
stopifnot(file.exists(paste0(LDREF, ".bed")))
log <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n", sep = "")

OUTCOMES <- list(t2d = "ebi-a-GCST006867", cad = "ebi-a-GCST005194", fbg = "ebi-a-GCST005186")
OUT_N <- c(t2d = "GCST006867", cad = "GCST005194", fbg = "GCST005186")
ID2SHORT <- setNames(names(OUTCOMES), unname(unlist(OUTCOMES)))
TYPE <- c(t2d = "cc", cad = "cc", fbg = "quant")
S_CASE <- c(t2d = 61714 / 655666, cad = 34541 / 296525)

PROTEINS <- list(
  list(file = "5231_79_PCSK9_PCSK9.txt.gz",          gene = "PCSK9",   chr = "1",  hg19 = 55505221),
  list(file = "6461_54_APOC3_Apo_C_III.txt.gz",      gene = "APOC3",   chr = "11", hg19 = 116700422),
  list(file = "3448_13_INSR_IR.txt.gz",              gene = "INSR",    chr = "19", hg19 = 7294045)
)
TARGETS <- data.table(
  gene = c("APOC3", "INSR", "APOC3", "PCSK9", "APOC3"),
  outcome = c("ebi-a-GCST005186", "ebi-a-GCST005186", "ebi-a-GCST006867",
              "ebi-a-GCST005194", "ebi-a-GCST005194"),
  label = c("APOC3×FBG(争议)", "INSR×FBG(争议)", "APOC3×T2D(争议)",
            "PCSK9×CAD(校准)", "APOC3×CAD(校准)"))

load_pqtl <- function(prot) {
  f <- file.path(SUB, paste0(sub("\\.gz$", "", prot$file), "_cis.txt.gz"))
  if (!file.exists(f)) return(data.table())
  d <- fread(cmd = paste0("zcat ", shQuote(f)), sep = "\t", header = TRUE, nThread = 4)
  d <- d[effectAllele != otherAllele & !otherAllele %in% c("!", "", NA)]
  d[, rsid2 := fifelse(is.na(rsids) | rsids %in% c("", ".", "-"),
                       Name, sub("[,; ].*$", "", rsids))]
  d <- d[!duplicated(rsid2)]
  d[, `:=`(beta = as.numeric(Beta), varbeta = as.numeric(SE)^2,
           maf = as.numeric(ImpMAF), Np = as.numeric(N))]
  # deCODE readme：ImpMAF 即次等位频率；pmin(1-ImpMAF,ImpMAF) 对 NA 无效 → 移除，NA 由下方过滤剔除
  d <- d[!is.na(beta) & !is.na(varbeta) & !is.na(maf) & varbeta > 0 &
           maf >= 0.01 & maf <= 0.99]
  d
}
pqtl_list <- lapply(PROTEINS, function(p) list(gene = p$gene, d = load_pqtl(p)))
names(pqtl_list) <- sapply(PROTEINS, `[[`, "gene")

is_pal <- function(a, b) (a == "A" & b == "T") | (a == "T" & b == "A") |
                          (a == "C" & b == "G") | (a == "G" & b == "C")
harmonize <- function(e, g) {
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
  keep <- (ok | flip) & !pal
  out <- data.table(snp = m$snp[keep], maf = m$maf_e[keep],
                    beta_e = m$beta_e[keep], varbeta_e = m$varbeta_e[keep], N_e = m$N_e[keep],
                    beta_g = gb[keep], varbeta_g = m$se_g[keep]^2, N_g = m$N_g[keep])
  out <- out[maf >= 0.01 & maf <= 0.99 & N_e > 0 & N_g > 0 &
               !is.na(beta_e) & !is.na(beta_g) & !is.na(varbeta_e) & !is.na(varbeta_g) &
               varbeta_e > 0 & varbeta_g > 0]
  out
}

ld_from_plink <- function(m) {
  snplist <- tempfile(fileext = ".txt")
  outpref <- tempfile(fileext = "")
  writeLines(m$snp, snplist)
  ok <- system(sprintf("%s --bfile %s --extract %s --r square --write-snplist --out %s 2>/dev/null",
                       PLINK, LDREF, snplist, outpref)) == 0
  if (!ok) return(list(order = character(), LD = NULL))
  ord <- readLines(paste0(outpref, ".snplist"))
  if (length(ord) < 10) return(list(order = ord, LD = NULL))
  ld <- as.matrix(read.table(paste0(outpref, ".ld")))
  colnames(ld) <- rownames(ld) <- ord
  list(order = ord, LD = ld)
}

MAX_SUSIE <- 300
# 变异选择候选集（收敛回退：方案B 优先，不收敛退方案A）：
#   方案A：top-300 by min(pqtl,gwas) z² —— 收敛最稳，但会丢"一侧极强另一侧中等"的因果变异
#          （例：PCSK9 rs11591147 R46L pQTL z=23 / GWAS z=9.7，min-z² 排名被挤出前 300）。
#   方案B：方案A ∪ 各性状单侧 top-20 —— 保双侧最强信号，但多变异可能破坏 susie_rss 收敛。
select_candidates <- function(m2, z_e, z_g) {
  A <- m2[order(-pmin(z_e, z_g))]$snp[seq_len(MAX_SUSIE)]
  B <- unique(c(A, m2[order(-z_e)]$snp[seq_len(20)], m2[order(-z_g)]$snp[seq_len(20)]))
  list(A = A, B = B)
}

susie_coloc_pair <- function(gene, on, e, g) {
  m <- harmonize(e, g)
  if (nrow(m) < 10) return(list(ok = FALSE, note = paste0("交集/对齐后仅 ", nrow(m), " 变异")))
  lk <- ld_from_plink(m)
  if (is.null(lk$LD)) return(list(ok = FALSE, note = paste0("LD 面板可用交集 ", length(lk$order), " 变异（<10）")))
  m2 <- m[snp %in% lk$order]
  ord <- lk$order[lk$order %in% m2$snp]
  m2 <- m2[match(ord, snp)]
  z_e <- abs(m2$beta_e) / sqrt(m2$varbeta_e)
  z_g <- abs(m2$beta_g) / sqrt(m2$varbeta_g)
  # 收敛回退：方案B（保因果变异）优先，不收敛退方案A（纯 top-300 min-z²，最稳）
  cands <- select_candidates(m2, z_e, z_g)
  fit <- NULL; used_cand <- NA_character_
  for (cn in c("B", "A")) {
    keep <- cands[[cn]]
    m3 <- m2[snp %in% keep]
    ord3 <- keep[keep %in% m3$snp]
    m3 <- m3[match(ord3, snp)]
    LD <- lk$LD[ord3, ord3]
    d1 <- list(snp = m3$snp, type = "quant", N = round(median(m3$N_e)),
               beta = m3$beta_e, varbeta = m3$varbeta_e, MAF = m3$maf, LD = LD)
    d2 <- if (TYPE[[on]] == "cc")
      list(snp = m3$snp, type = "cc", N = round(median(m3$N_g)), s = S_CASE[[on]],
           beta = m3$beta_g, varbeta = m3$varbeta_g, MAF = m3$maf, LD = LD) else
      list(snp = m3$snp, type = "quant", N = round(median(m3$N_g)),
           beta = m3$beta_g, varbeta = m3$varbeta_g, MAF = m3$maf, LD = LD)
    f1 <- tryCatch(susie_rss(z = d1$beta / sqrt(d1$varbeta), R = d1$LD, n = d1$N, L = 5,
                             estimate_prior_variance = FALSE, estimate_residual_variance = FALSE,
                             max_iter = 10000), error = function(e) NULL)
    f2 <- tryCatch(susie_rss(z = d2$beta / sqrt(d2$varbeta), R = d2$LD, n = d2$N, L = 5,
                             estimate_prior_variance = FALSE, estimate_residual_variance = FALSE,
                             max_iter = 10000), error = function(e) NULL)
    if (!is.null(f1) && !is.null(f2) && !is.null(f1$converged) && !is.null(f2$converged) &&
        f1$converged && f2$converged) { fit <- list(f1 = f1, f2 = f2, m3 = m3, LD = LD); used_cand <- cn; break }
  }
  if (is.null(fit))
    return(list(ok = FALSE, note = "susie_rss 10000 次迭代内不收敛（方案A/B 均试，LD 面板与人群不匹配）"))
  fit1 <- fit$f1; fit2 <- fit$f2; m2 <- fit$m3
  su <- tryCatch(coloc.susie(fit1, fit2), error = function(e) NULL)
  if (is.null(su)) return(list(ok = FALSE, note = "coloc.susie 失败"))
  s <- su$summary
  if (is.null(s) || nrow(s) == 0 || !"PP.H4.abf" %in% names(s))
    return(list(ok = FALSE, note = paste0("coloc.susie 无可信集对（pQTL 信号=", length(fit1$sets$cs),
                                          "/ GWAS 信号=", length(fit2$sets$cs), "），无法输出 PP.H4")))
  # 每 CS 对一行（coloc.susie 原始输出），保留 hit1/hit2 代表变异 + 各 PP
  detail <- data.table(hit1 = s$hit1, hit2 = s$hit2,
                       PP.H0 = s$PP.H0.abf, PP.H1 = s$PP.H1.abf,
                       PP.H2 = s$PP.H2.abf, PP.H3 = s$PP.H3.abf,
                       PP.H4 = s$PP.H4.abf)
  # 单行汇总：max PP.H4 + 对应 hit 对
  best <- detail[which.max(PP.H4)]
  list(ok = TRUE, note = "", nsnp = nrow(m2), detail = detail,
       PP.H4_susie_max = best$PP.H4, hit1 = best$hit1, hit2 = best$hit2,
       n_signal1 = length(fit1$sets$cs), n_signal2 = length(fit2$sets$cs),
       cs1 = paste(sapply(fit1$sets$cs, function(x) m2$snp[x[1]]), collapse = "|"),
       cs2 = paste(sapply(fit2$sets$cs, function(x) m2$snp[x[1]]), collapse = "|"),
       used_cand = used_cand)
}

sum_rows <- vector("list", nrow(TARGETS)); det_list <- vector("list", nrow(TARGETS))
for (i in seq_len(nrow(TARGETS))) {
  t <- TARGETS[i]
  on <- ID2SHORT[[t$outcome]]
  prot <- PROTEINS[[match(t$gene, sapply(PROTEINS, `[[`, "gene"))]]
  g <- readRDS(file.path(gdir, "_coloc_gwas_prot", paste0(OUT_N[[on]], "_", t$gene, ".rds")))
  e <- pqtl_list[[t$gene]]$d
  abf <- tryCatch({
    mm <- harmonize(e, g)
    if (nrow(mm) < 10) NA_real_ else {
      dd1 <- list(snp = mm$snp, type = "quant", N = round(median(mm$N_e)),
                  beta = mm$beta_e, varbeta = mm$varbeta_e, MAF = mm$maf)
      dd2 <- if (TYPE[[on]] == "cc")
        list(snp = mm$snp, type = "cc", N = mm$N_g, s = S_CASE[[on]],
             beta = mm$beta_g, varbeta = mm$varbeta_g, MAF = mm$maf) else
        list(snp = mm$snp, type = "quant", N = mm$N_g,
             beta = mm$beta_g, varbeta = mm$varbeta_g, MAF = mm$maf)
      r <- coloc.abf(dd1, dd2, p12 = 1e-5)
      as.numeric(r$summary["PP.H4.abf"])
    }
  }, error = function(e) NA_real_)
  su <- susie_coloc_pair(t$gene, on, e, g)
  sum_rows[[i]] <- data.frame(
    pair = t$label, gene = t$gene, outcome_short = toupper(on),
    PP.H4_abf = round(abf, 4),
    PP.H4_susie_max = if (su$ok) round(su$PP.H4_susie_max, 4) else NA,
    hit1 = if (su$ok) su$hit1 else NA, hit2 = if (su$ok) su$hit2 else NA,
    nsnp_susie = if (su$ok) su$nsnp else NA,
    n_signal_pqtl = if (su$ok) su$n_signal1 else NA,
    n_signal_gwas = if (su$ok) su$n_signal2 else NA,
    cand = if (su$ok) su$used_cand else NA,
    ok = su$ok, note = if (su$ok) "" else su$note)
  if (su$ok) det_list[[i]] <- data.table(pair = t$label, su$detail)
  log("  SuSiE ", i, "/", nrow(TARGETS), " | ", t$label,
      " | abf=", round(abf, 3), " susie_max=",
      if (su$ok) round(su$PP.H4_susie_max, 3) else su$note,
      " | 信号 pQTL=", if (su$ok) su$n_signal1 else "-",
      "/GWAS=", if (su$ok) su$n_signal2 else "-",
      if (su$ok) paste0(" | top:", su$hit1, "×", su$hit2) else "")
}
sum_out <- rbindlist(sum_rows, fill = TRUE)
write.csv(sum_out, file.path(gdir, "protein_coloc_susie_summary.csv"), row.names = FALSE)
det_out <- rbindlist(det_list, fill = TRUE)
if (nrow(det_out)) write.csv(det_out, file.path(gdir, "protein_coloc_susie_detail.csv"), row.names = FALSE)
print(sum_out[, .(pair, PP.H4_abf, PP.H4_susie_max, hit1, hit2, n_signal_pqtl, n_signal_gwas, note)])
log("=== M5 protein coloc-SuSiE 完成 ✔ | 产物 summary.csv + detail.csv ===")
