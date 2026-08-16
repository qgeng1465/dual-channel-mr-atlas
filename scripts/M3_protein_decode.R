#!/usr/bin/env Rscript
# =============================================================================
# M3_protein_decode.R — 蛋白通道 cis-MR：deCODE 血浆 pQTL（n=35,559, Ferkingstad 2021）
# =============================================================================
# 数据源：deCODE summary data（官网 download.decode.is，用户申请表下载）
#   主输入：data/decode/sub/<SeqId>_<Chr>_<Gene>_<Protein>.txt_cis.txt.gz
#           （M1_decode_subset.sh 裁剪的 TSS±1Mb cis 窗，hg38，小文件）
#   兜底输入：data/decode/<SeqId>_<Gene>_<Protein>.txt.gz（整文件，zcat 流式过滤）
#   列：Chrom, Pos(hg38), Name, rsids, effectAllele, otherAllele, Beta(per-SD),
#       Pval, minus_log10_pval, SE, N, ImpMAF
# 方法（预注册锁定）：cis ±1Mb(TSS)、p<5e-6、LD clump r²<0.01@1000kb EUR(API)
#   主 mr_ivw_mre；敏感性 ivw_fe/WM/Egger；nsnp=1 → mr_wald_ratio
# eaf 处理（诚实声明）：deCODE readme 注明 ImpMAF 不总是 effect allele 频率。
#   ① assocvariants.annotated（effectAlleleFreq，Name 连接）若完整则增强；
#   ② 否则 palindromic（A/T、C/G）一律排除（不依赖 eaf 判向，最保守），
#      ImpMAF 仅作为非 palindromic 的 eaf 近似并在 note 标注局限。
# 坐标：hg38 ENSEMBL REST 实查（2026-08-06），strand− 取 end 为 TSS。
# 纪律：预注册哈希校验；全网格落盘含空结果；工具数/排除数如实报告。
#
# 用法：PATH=$R_ENV/bin:$PATH Rscript scripts/M3_protein_decode.R
# 输出：results/grid/protein_decode_mr.csv + results/funnel/funnel_protein_decode.tsv
# =============================================================================
suppressMessages({
  library(data.table)
  suppressPackageStartupMessages(library(TwoSampleMR))
  suppressPackageStartupMessages(library(ieugwasr))
  library(jsonlite)
})

proj <- "/data/qiushuogeng/projects/dual-channel-mr-atlas"
res  <- file.path(proj, "results")
SUB  <- file.path(proj, "data/decode/sub")
dir.create(file.path(res, "grid"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(res, "funnel"), showWarnings = FALSE)
log <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n")

prereg <- file.path(proj, "docs/PREREGISTRATION.md")
lock   <- file.path(proj, "docs/PREREGISTRATION.md.sha256")
stopifnot(file.exists(prereg), file.exists(lock),
          tools::md5sum(prereg) == readLines(lock))
cfg <- read_json(file.path(res, "config.json"), auto_unbox = TRUE)
log("预注册哈希校验通过 ✔  | 主方法=", cfg$mr_methods$primary)
# --only GENE1,GENE2：只重跑指定蛋白（2026-08-13 用于 DPP4 修复 + INSR×fbg 401 重试，
#   避免扰动已复现的其它蛋白×结局对）。未提供则全量跑。
ONLY <- NULL
av <- commandArgs(trailingOnly = TRUE)
if (length(av) > 0 && startsWith(av[1], "--only="))
  ONLY <- unlist(strsplit(sub("^--only=", "", av[1]), ","))
OUT_CSV <- file.path(res, "grid/protein_decode_mr.csv")

# --- 蛋白文件 + hg38 坐标（TSS±1Mb，ENSEMBL 实查）-----------------------------
# 文件名与 gene/染色体/TSS(GRCh38)：strand− 的 TSS=end
PROTEINS <- list(
  list(file = "5231_79_PCSK9_PCSK9.txt.gz",   gene = "PCSK9",   chr = "1",  tss = 55039445, note = "positive control"),
  list(file = "5230_99_HMGCR_HMGR.txt.gz",    gene = "HMGCR",   chr = "5",  tss = 75336329, note = "negative control"),
  list(file = "10391_1_ANGPTL3_ANGL3.txt.gz", gene = "ANGPTL3", chr = "1",  tss = 62597464, note = "negative control"),
  list(file = "6461_54_APOC3_Apo_C_III.txt.gz", gene = "APOC3", chr = "11", tss = 116827019, note = "negative control"),
  list(file = "2797_56_APOB_Apo_B.txt.gz",    gene = "APOB",    chr = "2",  tss = 21044075, note = "optional"),
  list(file = "13129_40_LDLR_LDLR.txt.gz",    gene = "LDLR",    chr = "19", tss = 11089418, note = "optional"),
  # T2D 药物靶点（P1-2 四态分类，deCODE 收录的 5 个）
  list(file = "13085_18_GLP1R_GLP1R.txt.gz",  gene = "GLP1R",   chr = "6",  tss = 39048562, note = "drug target GLP-1RA"),
  list(file = "15460_9_DPP4_CD26.txt.gz",     gene = "DPP4",    chr = "2",  tss = 162074639, note = "drug target DPP4i"),
  list(file = "3448_13_INSR_IR.txt.gz",       gene = "INSR",    chr = "19", tss = 7294443, note = "drug target insulin/insulin sens."),
  list(file = "18182_24_PCK1_PCKGC.txt.gz",   gene = "PCK1",    chr = "20", tss = 57546220, note = "drug target gluconeogenesis"),
  list(file = "4891_50_GCG_Glucagon.txt.gz",  gene = "GCG",     chr = "2",  tss = 162152404, note = "drug target glucagon")
)
OUTCOMES <- list(t2d = "ebi-a-GCST006867", cad = "ebi-a-GCST005194", fbg = "ebi-a-GCST005186")

# --- annotated eaf 增强（仅当 gzip 完整时加载）--------------------------------
ann_path <- file.path(proj, "data/decode/assocvariants.annotated.txt.gz")
ANN <- NULL
if (file.exists(ann_path) &&
    system(paste("gzip -t", shQuote(ann_path), "2>/dev/null"),
           ignore.stdout = TRUE, ignore.stderr = TRUE) == 0) {
  ANN <- fread(cmd = paste0("zcat ", shQuote(ann_path)), sep = "\t", header = TRUE, nThread = 4)
  setnames(ANN, c("Chrom","Pos","Name","rsids","effectAllele","otherAllele","effectAlleleFreq"))
  log("annotated 已加载: ", nrow(ANN), " 变异（提供 effectAlleleFreq）")
} else {
  log("⚠ annotated 未就绪/不完整 → palindromic 一律排除（保守），ImpMAF 作近似")
}

is_pal <- function(a, b) {
  (a == "A" & b == "T") | (a == "T" & b == "A") |
  (a == "C" & b == "G") | (a == "G" & b == "C")
}

# --- 数据入口：优先 sub/ cis 窗文件，否则整文件 zcat 流式过滤 -------------------
load_cis <- function(prot) {
  # 2026-08-13 P4：M1 命名 = ${base%.gz}_cis.txt.gz（去 .gz 再加 _cis.txt.gz）。
  #   原 paste0(file,"_cis.txt.gz") 拼错文件名 → 一直回退整文件 zcat（慢，且绕开 M1 裁剪）。
  sub_f <- file.path(SUB, paste0(sub("\\.gz$", "", prot$file), "_cis.txt.gz"))
  full_f <- file.path(proj, "data/decode", prot$file)
  f <- if (file.exists(sub_f) && file.size(sub_f) > 1000) sub_f
       else if (file.exists(full_f) && file.size(full_f) > 1e7) full_f else NULL
  if (is.null(f)) return(NULL)
  d <- fread(cmd = paste0("zcat ", shQuote(f)), sep = "\t", header = TRUE, nThread = 4)
  if (!"Pval" %in% names(d) && "minus_log10_pval" %in% names(d)) d[, Pval := 10^(-minus_log10_pval)]  # 2026-08-07 P2：实列名 minus_log10_pval
  if (!"Pval" %in% names(d)) return(NULL)
  # 2026-08-07 P2：deCODE readme 要求排除 quality-excluded 变异（assocvariants.excluded.txt.gz）
  excl_path <- file.path(proj, "data/decode/assocvariants.excluded.txt.gz")
  if (file.exists(excl_path) && file.size(excl_path) > 0) {
    ex <- tryCatch(fread(cmd = paste0("zcat ", shQuote(excl_path)), sep = "\t", header = TRUE, nThread = 2),
                   error = function(e) NULL)
    if (!is.null(ex) && nrow(ex) > 0 && "Name" %in% names(ex)) {
      n0 <- nrow(d); d <- d[!Name %in% ex$Name]
      if (nrow(d) < n0) log("  [quality-excluded 排除] ", n0 - nrow(d), " 行")
    }
  }
  d[, Chrom2 := as.character(gsub("^chr", "", Chrom))]
  d[, Pos := as.numeric(Pos)]
  d
}

run_decode_mr <- function(prot, outcome_id) {
  # 2026-08-07 P2：load_cis 包 tryCatch——损坏/半截文件通过 size 启发式时不能崩掉整个蛋白通道
  d <- tryCatch(load_cis(prot), error = function(e) { log("  [读取失败] ", prot$gene, ": ", conditionMessage(e)); NULL })
  if (is.null(d))
    return(data.frame(gene = prot$gene, outcome = outcome_id, nsnp = 0, b = NA, se = NA,
                      pval = NA, method = NA, ok = FALSE, note = "蛋白文件未下载/不完整/读取失败"))
  d <- d[Chrom2 == prot$chr & abs(Pos - prot$tss) <= 1e6 & Pval < cfg$instrument$pval_thresh]
  if (nrow(d) == 0)
    return(data.frame(gene = prot$gene, outcome = outcome_id, nsnp = 0, b = NA, se = NA,
                      pval = NA, method = NA, ok = FALSE, note = "cis ±1Mb、p<5e-6 内无变异"))
  # 排除 readme 注明的多等位基因 bug 行：summary 文件有时 effectAllele==otherAllele
  # （正确 otherAllele 应为 "!"，即"非 effectAllele 的多等位集合"）。这类变异无法可靠
  # harmonise，一律排除（annotated 文件才对这类行做了修正，而 annotated 未就绪）。
  n_ma <- d[, sum(effectAllele == otherAllele | otherAllele %in% c("!", "", NA), na.rm = TRUE)]
  d <- d[effectAllele != otherAllele & !otherAllele %in% c("!", "", NA)]
  # 工具标识：rsid 优先，缺失/NA 用 Name
  d[, rsid2 := fifelse(is.na(rsids) | rsids %in% c("", ".", "-"), Name, rsids)]
  d <- d[!duplicated(rsid2)]
  if (n_ma > 0) log("  [多等位/等位缺失排除] ", prot$gene, ": ", n_ma, " 行")
  # eaf：annotated 增强 → 否则 ImpMAF
  d[, eaf := ImpMAF]; d[, eaf_src := "ImpMAF 近似"]
  if (!is.null(ANN)) {
    setkey(ANN, Name)
    d[ANN, c("eaf", "eaf_src") := .(i.effectAlleleFreq, "annotated")]
    # 仅非回文变异用 ImpMAF 兜底：ImpMAF 是次等位频率（未必 = effect allele 频率），
    # 回文位点若缺少 annotated 效应等位频率而用 ImpMAF 判向会翻转错误 → 交下方过滤排除
    d[is.na(eaf) & !is_pal(effectAllele, otherAllele), eaf := ImpMAF]
  }
  # palindromic：有 annotated 频率则保留（可判向），否则排除
  npal <- sum(is_pal(d$effectAllele, d$otherAllele))
  if (is.null(ANN)) {
    d <- d[!is_pal(effectAllele, otherAllele)]
    note_eaf <- paste0("palindromic 排除 ", npal, "（无 annotated）")
  } else {
    d <- d[!(is_pal(effectAllele, otherAllele) & is.na(eaf))]
    note_eaf <- paste0("palindromic 保留（annotated 频率），排除无频率 ", sum(is.na(d$eaf)), " 个")
  }
  if (nrow(d) == 0)
    return(data.frame(gene = prot$gene, outcome = outcome_id, nsnp = 0, b = NA, se = NA,
                      pval = NA, method = NA, ok = FALSE, note = paste0(note_eaf, " → 无工具")))
  # LD clump（EUR 1000G API）
  cl <- tryCatch(ld_clump(data.frame(rsid = d$rsid2, pval = d$Pval, id = prot$gene),
                          clump_kb = cfg$instrument$clump_kb,
                          clump_r2 = cfg$instrument$clump_r2,
                          clump_p = 1, pop = "EUR",
                          opengwas_jwt = ieugwasr::get_opengwas_jwt()),
                 error = function(e) NULL)
  # 2026-08-13 P3（DPP4 修复）：clump 无输出/全被移除 → 无独立工具，如实报告而非继续抽取。
  if (is.null(cl) || nrow(cl) == 0)
    return(data.frame(gene = prot$gene, outcome = outcome_id, nsnp = 0, b = NA, se = NA,
                      pval = NA, method = NA, ok = FALSE,
                      note = "clump 后无独立工具（变异不在 LD 参考面板或全被 LD 吸收）"))
  # 2026-08-13 P3：OpenGWAS ld_clump 对"不在 1000G EUR 参考面板的罕见位点"行为不确定
  #   （DPP4 的 cis p<5e-6 变异全为罕见位点，曾出现 clump 全保留→outcome 抽取报 API 错的
  #   假错误行）。仅当 clump 保留全部候选且候选数≤30（罕见位点独占的小集场景）时，
  #   用 ld_reflookup 判参考面板存在性；正常蛋白 clump 会削减数量、不触发，结果不受影响。
  if (nrow(cl) >= nrow(d) && nrow(d) <= 30) {
    lk <- tryCatch(ieugwasr::ld_reflookup(d$rsid2, pop = "EUR",
                                          opengwas_jwt = ieugwasr::get_opengwas_jwt()),
                   error = function(e) NULL)
    if (!is.null(lk) && !any(d$rsid2 %in% lk))
      return(data.frame(gene = prot$gene, outcome = outcome_id, nsnp = 0, b = NA, se = NA,
                        pval = NA, method = NA, ok = FALSE,
                        note = "cis ±1Mb、p<5e-6 变异均不在 1000G EUR LD 参考面板（罕见位点）→ 无独立工具"))
  }
  d <- d[d$rsid2 %in% cl$rsid, ]
  exp_dat <- data.frame(
    SNP = d$rsid2,
    effect_allele.exposure = d$effectAllele,
    other_allele.exposure = d$otherAllele,
    eaf.exposure = d$eaf,
    beta.exposure = d$Beta, se.exposure = d$SE,
    pval.exposure = d$Pval, samplesize.exposure = d$N,
    id.exposure = prot$gene, exposure = prot$gene, stringsAsFactors = FALSE)
  tryCatch({
    out <- extract_outcome_data(snps = exp_dat$SNP, outcomes = outcome_id, proxies = TRUE)
    dat <- harmonise_data(exp_dat, out, action = 2)
    dat <- dat[dat$mr_keep, ]
    if (nrow(dat) == 0)
      return(data.frame(gene = prot$gene, outcome = outcome_id, nsnp = nrow(d),
                        b = NA, se = NA, pval = NA, method = NA, ok = FALSE,
                        note = "harmonise 后无保留 SNP"))
    r_main <- mr(dat, method_list = cfg$mr_methods$primary)
    r_sens <- tryCatch(mr(dat, method_list = cfg$mr_methods$sensitivity),
                       error = function(e) NULL)
    r_all <- if (!is.null(r_sens)) rbind(r_main, r_sens) else r_main
    if (nrow(r_all) == 0 && nrow(dat) == 1) {   # nsnp=1 → Wald
      r_w <- tryCatch(mr(dat, method_list = "mr_wald_ratio"), error = function(e) NULL)
      if (!is.null(r_w) && nrow(r_w) > 0) r_all <- r_w
    }
    if (nrow(r_all) == 0)
      return(data.frame(gene = prot$gene, outcome = outcome_id, nsnp = nrow(dat),
                        b = NA, se = NA, pval = NA, method = NA, ok = FALSE, note = "MR 无输出"))
    note_final <- paste0(ifelse(any(d$eaf_src == "ImpMAF 近似"), note_eaf, "eaf=annotated"),
                         "; clump ", nrow(exp_dat), "→harmonise ", nrow(dat))
    data.frame(gene = prot$gene, outcome = outcome_id, nsnp = nrow(dat),
               b = r_all$b, se = r_all$se, pval = r_all$pval,
               method = r_all$method, ok = TRUE, note = note_final)
  }, error = function(e) data.frame(gene = prot$gene, outcome = outcome_id, nsnp = 0,
                                    b = NA, se = NA, pval = NA, method = NA,
                                    ok = FALSE, note = conditionMessage(e)))
}

out <- list()
for (p in PROTEINS) {
  if (!is.null(ONLY) && !p$gene %in% ONLY) next
  # 2026-08-13 P4 一致性：M1 命名 = ${base%.gz}_cis.txt.gz（与 load_cis 相同），
  #   原 paste0(p$file, "_cis.txt.gz") 拼错 → 仅当整文件缺失时会误判"未就绪"跳过
  sub_f <- file.path(SUB, paste0(sub("\\.gz$", "", p$file), "_cis.txt.gz"))
  if (!file.exists(file.path(proj, "data/decode", p$file)) && !file.exists(sub_f)) {
    log("跳过 ", p$gene, "：文件未就绪")
    for (on in names(OUTCOMES))
      out[[length(out) + 1]] <- data.frame(gene = p$gene, outcome = OUTCOMES[[on]], nsnp = 0,
                                           b = NA, se = NA, pval = NA, method = NA, ok = FALSE,
                                           note = "蛋白文件未下载/不完整")
    next
  }
  for (on in names(OUTCOMES)) {
    log("=== ", p$gene, " × ", on, " ===")
    rr <- run_decode_mr(p, OUTCOMES[[on]])
    out[[length(out) + 1]] <- rr
    for (i in seq_len(nrow(rr)))
      if (rr$ok[i]) log("  ", p$gene, "×", on, "[", rr$method[i], "] nsnp=", rr$nsnp[i],
                        " b=", round(rr$b[i], 3), " p=", format(rr$pval[i], digits = 2))
  }
}
res_df <- as.data.table(do.call(rbind, out))
# 2026-08-13 resume 合并：重跑过的 gene×outcome 以本次为准，未重跑的旧行保留
if (file.exists(OUT_CSV)) {
  old <- fread(OUT_CSV)
  new_keys <- unique(res_df[, .(gene, outcome)])
  old <- old[!paste(gene, outcome) %in% new_keys[, paste(gene, outcome)]]
  res_df <- rbindlist(list(res_df, old), fill = TRUE)
}
write.csv(res_df, OUT_CSV, row.names = FALSE)
funnel <- data.frame(stage = c("tested_proteins", "gene_outcome_pairs", "mr_completed_rows"),
                     count = c(length(PROTEINS), length(PROTEINS) * length(OUTCOMES),
                               sum(res_df$ok, na.rm = TRUE)))
write.table(funnel, file.path(res, "funnel/funnel_protein_decode.tsv"), sep = "\t", row.names = FALSE)
log("deCODE 蛋白通道 MR 完成 ✔ → results/grid/protein_decode_mr.csv")
log("坐标修正：PCSK9 chr1:55.04Mb / HMGCR chr5:75.34Mb（hg38 ENSEMBL 实查，原脚本 hg19 坐标已修正）")
