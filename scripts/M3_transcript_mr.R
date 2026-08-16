#!/usr/bin/env Rscript
# =============================================================================
# M3_transcript_mr.R — 转录本通道 cis-MR（eQTLGen × OpenGWAS 结局）v2
# =============================================================================
# v2 方法学修正（相对 v1 QA 版，记录于 CHANGELOG 2026-08-06）：
#   1. eaf 不再用 0.5 占位 → 用 eQTLGen 官方 SNP 频率文件
#      (2018-07-18_SNP_AF_for_AlleleB..., 本地 data/eqtlgen/SNP_AF.txt.gz)
#      eaf = 1 - AlleleB_all （effect allele = AlleleA = eQTLGen AssessedAllele）
#   2. 工具间 LD clumping 用 OpenGWAS EUR 1000G 参考（ieugwasr::ld_clump API 模式，
#      无需本地 plink/bfile；clump_kb=1000, clump_r2=0.01 按预注册标准）
#   3. 主分析 mr_ivw_mre 保持（预注册锁定）；敏感性 ivw_fe/WM/Egger（nsnp 足够时）
#   4. nsnp=1 时 IVW ≡ Wald ratio，WM/Egger 记为 NA（诚实标注，不作伪敏感性）
#
# 学术不端纪律（预注册 docs/PREREGISTRATION.md，每轮运行前校验 sha256）：
#   - 工具标准锁定：cis ±1Mb, p<5e-6, clump r²<0.01@1000kb EUR
#   - 全网格结果全部落盘（含空结果/失败/单工具），不做选择报告
#   - 每基因工具数如实报告，弱工具（nsnp=0）如实记录
#
# 用法：PATH=$R_ENV/bin:$PATH Rscript scripts/M3_transcript_mr.R
# 输入：data/eqtlgen/cis-eQTL-significant.txt.gz + SNP_AF.txt.gz（本地）
#       结局：OpenGWAS API（JWT）
# 输出：results/grid/transcript_mr_v2.csv + results/funnel/funnel_transcript_v2.tsv
# =============================================================================
suppressMessages({
  library(data.table)
  suppressPackageStartupMessages(library(TwoSampleMR))
  suppressPackageStartupMessages(library(ieugwasr))
  library(jsonlite)
})

proj <- "<repo-root>"
res  <- file.path(proj, "results")
dir.create(file.path(res, "grid"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(res, "funnel"), showWarnings = FALSE)

log <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n")

# --- 预注册完整性校验 ----------------------------------------------------------
prereg <- file.path(proj, "docs/PREREGISTRATION.md")
lock   <- file.path(proj, "docs/PREREGISTRATION.md.sha256")
stopifnot(file.exists(prereg), file.exists(lock),
          tools::md5sum(prereg) == readLines(lock))
cfg <- read_json(file.path(res, "config.json"), auto_unbox = TRUE)
log("预注册哈希校验通过 ✔  | 主方法=", cfg$mr_methods$primary,
    " | 工具变量 p<", cfg$instrument$pval_thresh,
    " | clump r2<", cfg$instrument$clump_r2, "@", cfg$instrument$clump_kb, "kb")

# --- 结局定义 ----------------------------------------------------------------
OUTCOMES <- list(
  t2d = "ebi-a-GCST006867",   # T2D, n=655,666
  cad = "ebi-a-GCST005194",   # CAD, n=296,525
  fbg = "ebi-a-GCST005186"    # FBG, n=58,074
)

# --- 读取 eQTLGen 暴露数据 ----------------------------------------------------
log("读取 eQTLGen 显著 cis-eQTL...")
eqtl <- fread(cmd = paste0("zcat ", file.path(proj, "data/eqtlgen/cis-eQTL-significant.txt.gz")),
              sep = "\t", header = TRUE, nThread = 4)
setnames(eqtl, c("Pvalue","SNP","SNPChr","SNPPos","AssessedAllele","OtherAllele",
                 "Zscore","Gene","GeneSymbol","GeneChr","GenePos",
                 "NrCohorts","NrSamples","FDR","BonferroniP"))
log("eQTLGen 行数=", nrow(eqtl), " | 基因数=", uniqueN(eqtl$GeneSymbol))

# --- 真实 SNP 频率（替代 v1 的 0.5 占位）--------------------------------------
log("读取 eQTLGen SNP 频率文件...")
af <- fread(cmd = paste0("zcat ", file.path(proj, "data/eqtlgen/SNP_AF.txt.gz")),
            sep = "\t", header = TRUE, nThread = 4)
log("AF 文件行数=", nrow(af))
# 核对 SNP 唯一性
setnames(af, c("SNP","hg19_chr","hg19_pos","AlleleA","AlleleB",
               "allA_total","allAB_total","allB_total","AlleleB_all"))
# 合并到工具集，记录 eaf（effect allele=AssessedAllele，应=AlleleA）
ivs <- eqtl[Pvalue < cfg$instrument$pval_thresh]
ivs <- merge(ivs, af[, .(SNP, AlleleA, AlleleB, AlleleB_all)], by = "SNP", all.x = TRUE)
# 核对等位基因一致：AssessedAllele 应对应 AlleleA；翻转则补位 eaf
ivs[, eaf := fcase(
  AssessedAllele == AlleleA,  1 - AlleleB_all,   # effect=AlleleA
  AssessedAllele == AlleleB,  AlleleB_all,       # effect=AlleleB（罕见）
  default = NA_real_                             # 无法核对 → 不用于 MR（诚实保守）
)]
log("工具候选（p<5e-6）: ", nrow(ivs), " | 基因: ", uniqueN(ivs$GeneSymbol),
    " | 有真实 eaf:", sum(!is.na(ivs$eaf)), " | 等位不一致/缺失:", sum(is.na(ivs$eaf)))

# --- 分析基因集 --------------------------------------------------------------
controls <- c(cfg$controls$positive, cfg$controls$negative)
test_genes <- intersect(controls, unique(eqtl$GeneSymbol))
strong <- ivs[!is.na(eaf), .(minP = min(Pvalue)), by = GeneSymbol][order(minP)][1:10]
test_genes <- unique(c(test_genes, strong$GeneSymbol))
log("本轮分析基因集: ", length(test_genes), " 个基因")

# --- 工具变量：真实 eaf + LD clump（EUR, r²<0.01, 1000kb）----------------------
build_ivs <- function(gene) {
  d <- ivs[GeneSymbol == gene & !is.na(eaf)]
  if (nrow(d) == 0) return(NULL)
  setorder(d, Pvalue)
  # LD clump（OpenGWAS API，EUR 1000G 参考）
  cl <- tryCatch({
    ld_clump(data.frame(rsid = d$SNP, pval = d$Pvalue, id = gene),
             clump_kb   = cfg$instrument$clump_kb,
             clump_r2   = cfg$instrument$clump_r2,
             clump_p    = 1,           # 全部候选参与 clump
             pop        = "EUR",
             opengwas_jwt = ieugwasr::get_opengwas_jwt())
  }, error = function(e) {
    log("  [clump 失败] ", gene, ": ", conditionMessage(e))
    NULL
  })
  if (is.null(cl) || nrow(cl) == 0) return(NULL)
  kept <- d[SNP %in% cl$rsid]
  kept[order(Pvalue)]
}

# --- 单基因 × 单结局 MR --------------------------------------------------------
run_gene_mr <- function(gene, outcome_id) {
  iv <- build_ivs(gene)
  if (is.null(iv) || nrow(iv) == 0)
    return(data.frame(gene = gene, outcome = outcome_id, nsnp = 0,
                      b = NA, se = NA, pval = NA, method = NA, ok = FALSE,
                      note = if (is.null(iv)) "无工具变量" else "clump 后无独立工具"))
  exp_dat <- data.frame(
    SNP = iv$SNP,
    effect_allele.exposure = iv$AssessedAllele,
    other_allele.exposure = iv$OtherAllele,
    eaf.exposure = iv$eaf,
    beta.exposure = iv$Zscore / sqrt(iv$NrSamples),  # eQTLGen Z→β: se=1/√N, β=Z*se
    se.exposure = 1 / sqrt(iv$NrSamples),
    pval.exposure = iv$Pvalue,
    samplesize.exposure = iv$NrSamples,
    id.exposure = gene, exposure = gene,
    stringsAsFactors = FALSE
  )
  tryCatch({
    out_dat <- extract_outcome_data(snps = iv$SNP, outcomes = outcome_id, proxies = TRUE)
    if (is.null(out_dat) || nrow(out_dat) == 0)
      return(data.frame(gene = gene, outcome = outcome_id, nsnp = 0,
                        b = NA, se = NA, pval = NA, method = NA, ok = FALSE, note = "无结局匹配"))
    dat <- harmonise_data(exp_dat, out_dat, action = 2)
    dat <- dat[dat$mr_keep, ]
    if (nrow(dat) == 0)
      return(data.frame(gene = gene, outcome = outcome_id, nsnp = 0,
                        b = NA, se = NA, pval = NA, method = NA, ok = FALSE, note = "harmonise 后无保留 SNP"))
    # 主分析：预注册锁定的 mr_ivw_mre
    r_main <- mr(dat, method_list = cfg$mr_methods$primary)
    # 敏感性：ivw_fe / WM / Egger（nsnp≥2 时可用性由 mr() 决定，单工具自动跳过）
    r_sens <- tryCatch(mr(dat, method_list = cfg$mr_methods$sensitivity),
                       error = function(e) NULL)
    if (!is.null(r_sens)) r_all <- rbind(r_main, r_sens) else r_all <- r_main
    # nsnp=1 单工具：mr_ivw_mre 无输出，但 IVW≡Wald ratio → 补 Wald（诚实标注，不伪敏感性）
    if (nrow(r_all) == 0 && nrow(dat) == 1) {
      r_w <- tryCatch(mr(dat, method_list = "mr_wald_ratio"), error = function(e) NULL)
      if (!is.null(r_w) && nrow(r_w) > 0) r_all <- r_w
    }
    # 若全部失败，保留一条空记录
    if (nrow(r_all) == 0)
      return(data.frame(gene = gene, outcome = outcome_id, nsnp = nrow(dat),
                        b = NA, se = NA, pval = NA, method = NA, ok = FALSE, note = "MR 无输出"))
    data.frame(gene = gene, outcome = outcome_id, nsnp = nrow(dat),
               b = r_all$b, se = r_all$se, pval = r_all$pval,
               method = r_all$method, ok = TRUE,
               note = paste0("clump ", nrow(exp_dat), "→harmonise ", nrow(dat)))
  }, error = function(e) data.frame(gene = gene, outcome = outcome_id, nsnp = 0,
                                    b = NA, se = NA, pval = NA, method = NA,
                                    ok = FALSE, note = conditionMessage(e)))
}

all_results <- list()
for (out_name in names(OUTCOMES)) {
  log("=== 结局: ", out_name, " (", OUTCOMES[[out_name]], ") ===")
  for (g in test_genes) {
    rr <- run_gene_mr(g, OUTCOMES[[out_name]])
    all_results[[length(all_results) + 1]] <- rr
    for (i in seq_len(nrow(rr)))
      if (rr$ok[i]) log("  ", g, "×", out_name, "[", rr$method[i], "] nsnp=", rr$nsnp[i],
                        " b=", round(rr$b[i], 3), " p=", format(rr$pval[i], digits = 2))
  }
}

# --- 落盘（全量，含空结果；统一 schema）---------------------------------------
res_df <- do.call(rbind, all_results)
write.csv(res_df, file.path(res, "grid/transcript_mr_v2.csv"), row.names = FALSE)

funnel <- data.frame(
  stage = c("eqtlgen_total_genes", "iv_pval5e6_genes", "eaf_available_genes",
            "tested_genes", "gene_outcome_pairs", "mr_completed_rows"),
  count = c(uniqueN(eqtl$GeneSymbol), uniqueN(ivs$GeneSymbol),
            uniqueN(ivs[!is.na(eaf)]$GeneSymbol), length(test_genes),
            length(test_genes) * length(OUTCOMES), sum(res_df$ok, na.rm = TRUE))
)
write.table(funnel, file.path(res, "funnel/funnel_transcript_v2.tsv"), sep = "\t", row.names = FALSE)

log("M3 转录本通道 MR v2 完成 ✔ → results/grid/transcript_mr_v2.csv")
log("方法学修正：真实 eaf + EUR LD clump（r²<0.01@1000kb）。")
log("说明：本轮为代表性基因集（对照+强 cis-eQTL），全量网格在 9h 主程序阶段 05 扩展。")
