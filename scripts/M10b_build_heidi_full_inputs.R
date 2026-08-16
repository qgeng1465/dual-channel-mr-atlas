#!/usr/bin/env Rscript
# =============================================================================
# M10b_build_heidi_full_inputs.R — HEIDI 全量输入构建（P3：818 转录 MR 显著测试）
# =============================================================================
# 目标：对 results/grid/transcript_coloc.csv 的全部唯一基因（819 行 MR 显著测试、
#   759 个唯一 ENSG）构建 SMR/HEIDI 输入，把已跑的 107 探针扩展到全量：
#     * data/smr/trans_esd_full/{ENSG}.esd     — 每基因一个探针文件
#     * data/smr/trans_flist_full.txt          — SMR --eqtl-flist 输入
#     * data/smr/trans_{t2d,cad,fbg}_full.ma   — 各结局 GWAS 汇总（复用 _coloc_gwas 缓存）
#
# 本脚本不运行 SMR 二进制（--make-besd 与 HEIDI 是独立步骤，需经资源仲裁跑，
#   见 M10_run_smr.sh 的用法）。
#
# 与 M10_transcript_smr_build_inputs.R 完全一致的约定（2026-08-13 复核）：
#   * eQTLGen cis-eQTL-significant 为 hg19（SNPPos 与 1kg bim 逐位点一致，已验证）
#     → .esd 直接用文件内 SNPChr/SNPPos，无需坐标转换。
#   * eQTLGen 只有 Zscore → Beta = Z/sqrt(NrSamples)，se = 1/sqrt(NrSamples)。
#   * 频率：SNP_AF 的 AlleleB_all 按 AssessedAllele 取向算 eaf；缺失则剔除该 SNP。
#   * LD 参考面板：1kg EUR (hg19)，.esd 只保留参考面板内变异（SMR 需参考算 LD）。
#   * flist 列：Chr ProbeID GeneticDistance ProbeBp Gene Orientation PathOfEsd；
#     Chr = 该基因显著 cis-SNP 首行的 SNPChr；ProbeBp = round(mean(range(SNPPos)))；
#     Orientation = "+"。.esd 列：Chr SNP Bp A1 A2 Freq Beta se p（空格分隔，含表头）。
#   * GWAS .ma：离线复用 _coloc_gwas/{OUTCOME}_{ENSG}.rds 缓存（与 coloc 完全同源）。
#
# 纪律：所有产物均为新建文件（trans_esd_full/、trans_flist_full.txt、
#   trans_*_full.ma），绝不覆盖已跑完的 107 探针输入
#   （data/smr/trans_esd/、data/smr/trans_flist.txt、data/smr/eqtlgen_trans.*、
#   data/smr/trans_{t2d,cad,fbg}.ma）。
#
# 用法：PATH=<conda-root>/r-mr/bin:$PATH \
#       Rscript scripts/M10b_build_heidi_full_inputs.R
# =============================================================================
suppressMessages(library(data.table))
proj   <- "<repo-root>"
gdir   <- file.path(proj, "results/grid")
LDREF  <- file.path(proj, "data/ldref/1kg.v3/EUR")
OUT    <- file.path(proj, "data/smr")
ODIR   <- file.path(OUT, "trans_esd_full")
dir.create(ODIR, recursive = TRUE, showWarnings = FALSE)
EQTL_SIG <- file.path(proj, "data/eqtlgen/cis-eQTL-significant.txt.gz")
AF       <- file.path(proj, "data/eqtlgen/SNP_AF.txt.gz")
stopifnot(file.exists(EQTL_SIG), file.exists(AF), file.exists(paste0(LDREF, ".bim")))

OUT_N <- c(t2d = "GCST006867", cad = "GCST005194", fbg = "GCST005186")

# ---- 目标基因集：transcript_coloc.csv 全部唯一基因 ----
coloc <- fread(file.path(gdir, "transcript_coloc.csv"))
targets <- unique(coloc[, .(gene, symbol)])
cat("目标基因数:", nrow(targets), "(transcript_coloc.csv 唯一基因)\n")
cat("  MR 显著测试行数:", nrow(coloc), "，无 ENSG 缺失:",
    !any(is.na(targets$gene)), "\n")

# ---- 单趟流式抽取目标基因的显著 cis-eQTL ----
gene_ids <- targets$gene
tmp <- file.path(ODIR, "_extract.tsv")
if (file.exists(tmp)) unlink(tmp)
idf <- file.path(ODIR, "_ids.txt")
writeLines(gene_ids, idf)
cmd <- sprintf("zcat %s | awk -F'\\t' 'NR==FNR{a[$1]=1;next} (FNR>1) && a[$8]{print}' %s - > %s",
               shQuote(EQTL_SIG), shQuote(idf), shQuote(tmp))
cat("抽取显著 cis-eQTL（单趟流式）...\n")
system(cmd, wait = TRUE)
eq <- fread(tmp, sep = "\t", header = TRUE, nThread = 4)
setnames(eq, c("Pvalue","SNP","SNPChr","SNPPos","AssessedAllele","OtherAllele",
               "Zscore","Gene","GeneSymbol","GeneChr","GenePos",
               "NrCohorts","NrSamples","FDR","BonferroniP"))
cat("  ", nrow(eq), "SNP-gene 对,", uniqueN(eq$Gene), "基因（显著 cis-eQTL 可用性）\n")

# ---- AF + 1kg bim ----
af  <- fread(cmd = paste0("zcat ", shQuote(AF)), sep = "\t", header = TRUE, nThread = 4)
setnames(af, c("SNP","hg19_chr","hg19_pos","AlleleA","AlleleB",
               "allA_total","allAB_total","allB_total","AlleleB_all"))
bim <- fread(paste0(LDREF, ".bim"), header = FALSE,
             col.names = c("chr", "snp", "cm", "pos", "a1", "a2"))
setkey(bim, snp)

eq2 <- merge(eq, af[, .(SNP, AlleleA, AlleleB, AlleleB_all)], by = "SNP", all.x = TRUE)
eq2[, eaf := fifelse(AssessedAllele == AlleleA, 1 - AlleleB_all,
                 fifelse(AssessedAllele == AlleleB, AlleleB_all, NA_real_))]
# 仅保留参考面板内变异（SMR 需参考算 LD），且频率有效
eq2 <- merge(eq2, bim[, .(chr, snp, pos)], by.x = "SNP", by.y = "snp")
eq2 <- eq2[!is.na(eaf) & eaf > 0 & eaf < 1 & !is.na(Zscore) & NrSamples > 0]
cat("  合并 AF + 1kg 交集后:", nrow(eq2), "行,", uniqueN(eq2$Gene), "基因\n")

# ---- 逐基因写 .esd + flist ----
eq2[, `:=`(Beta = Zscore / sqrt(NrSamples), se = 1 / sqrt(NrSamples))]
genes_esd <- sort(unique(eq2$Gene))
cat("  ", length(genes_esd), "个基因生成 .esd\n")
fl <- lapply(genes_esd, function(g) {
  d <- eq2[Gene == g]
  setorder(d, Pvalue)
  esd <- d[, .(Chr = SNPChr, SNP, Bp = SNPPos,
               A1 = AssessedAllele, A2 = OtherAllele, Freq = eaf,
               Beta, se, p = Pvalue)]
  fn <- file.path(ODIR, paste0(g, ".esd"))
  fwrite(esd, fn, sep = " ", quote = FALSE, na = "NA")
  data.table(Chr = unique(d$SNPChr)[1], ProbeID = g, GeneticDistance = 0,
             ProbeBp = round(mean(range(d$SNPPos))), Gene = d$GeneSymbol[1],
             Orientation = "+", PathOfEsd = fn)
})
fl <- rbindlist(fl)
flf <- file.path(OUT, "trans_flist_full.txt")
fwrite(fl, flf, sep = " ", quote = FALSE, na = "NA")
cat("flist:", nrow(fl), "个探针 ->", flf, "\n")

# 缺 .esd 的基因（显著 cis-eQTL 但在 AF/1kg 过滤后无存活 SNP）
missing <- setdiff(targets$gene, genes_esd)
if (length(missing)) {
  cat("  未生成 .esd 的基因:", length(missing), "个:\n")
  print(targets[gene %in% missing])
} else {
  cat("  全部", length(genes_esd), "个目标基因均生成 .esd\n")
}

# ---- 每结局 .ma（复用 _coloc_gwas 缓存；与 M10 同源同格式）----
gcol <- file.path(gdir, "_coloc_gwas")
for (on in names(OUT_N)) {
  fpat <- file.path(gcol, paste0(OUT_N[[on]], "_", targets$gene, ".rds"))
  ok <- file.exists(fpat)
  cat("  outcome", on, ":", sum(ok), "/", nrow(targets), "基因有 GWAS 缓存\n")
  parts <- lapply(targets$gene[ok], function(g) {
    d <- readRDS(file.path(gcol, paste0(OUT_N[[on]], "_", g, ".rds")))
    d[, .(snp = rsid, A1 = ea, A2 = nea, freq = eaf, b = beta, se = se, p = p, n = n)]
  })
  ma <- rbindlist(parts)
  ma <- ma[!duplicated(snp)]
  ma <- ma[A1 != A2 & !is.na(b) & !is.na(se) & se > 0 & !is.na(p)]
  maf <- file.path(OUT, paste0("trans_", on, "_full.ma"))
  fwrite(ma, maf, sep = " ", quote = FALSE, na = "NA")
  cat("  ->", maf, ":", nrow(ma), "SNPs\n")
}

cat("=== HEIDI 全量输入构建完成 ===\n")
cat("  .esd:", length(list.files(ODIR, pattern = "\\.esd$")),
    "个 ->", ODIR, "\n")
cat("  flist:", nrow(fl), "行 ->", flf, "\n")
