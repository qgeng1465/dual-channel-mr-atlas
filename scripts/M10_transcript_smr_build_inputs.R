#!/usr/bin/env Rscript
# =============================================================================
# M10_transcript_smr_build_inputs.R — 转录通道 SMR/HEIDI 输入构建
# =============================================================================
# 目标：对转录通道 106 个 strong coloc 命中（PP.H4≥0.8）+ KCNJ11（药物靶点星）
#   跑 SMR + HEIDI（Zhu 2016），eQTL 侧用 eQTLGen 全血 cis-eQTL。
#
# 设计与诚实声明（2026-08-13）：
#   * eQTLGen cis-eQTL-significant 文件为 hg19（SNPPos 与 1kg bim 逐位点一致，已验证）
#     → .esd 直接用文件内 SNPChr/SNPPos，无需坐标转换。
#   * eQTLGen 只有 Zscore → Beta = Z/sqrt(NrSamples)，se = 1/sqrt(NrSamples)
#     （标准化表达量 var≈1 的常规转换，与转录 MR 管道一致）。
#   * 频率：SNP_AF（AlleleB_all）按 AssessedAllele 取向计算 eaf；缺失则剔除该 SNP。
#   * LD 参考面板：1kg EUR（hg19），.esd 只保留参考面板内的变异（SMR 需参考算 LD）。
#   * GWAS .ma：离线复用 _coloc_gwas/{OUTCOME}_{ENSG}.rds 缓存（与 coloc 完全同源）。
#
# 产物：
#   data/smr/trans_esd/{ENSG}.esd, data/smr/trans_flist.txt
#   → data/smr/eqtlgen_trans.besd/.esi/.epi
#   data/smr/trans_{t2d,cad,fbg}.ma
# 用法：PATH=$R_ENV/bin:$PATH Rscript scripts/M10_transcript_smr_build_inputs.R
# =============================================================================
suppressMessages(library(data.table))
proj <- "/data/qiushuogeng/projects/dual-channel-mr-atlas"
gdir <- file.path(proj, "results/grid")
LDREF <- file.path(proj, "data/ldref/1kg.v3/EUR")
OUT <- file.path(proj, "data/smr")
dir.create(file.path(OUT, "trans_esd"), recursive = TRUE, showWarnings = FALSE)
EQTL_SIG <- file.path(proj, "data/eqtlgen/cis-eQTL-significant.txt.gz")
AF <- file.path(proj, "data/eqtlgen/SNP_AF.txt.gz")

OUT_N <- c(t2d = "GCST006867", cad = "GCST005194", fbg = "GCST005186")

# ---- 目标基因集：106 strong 命中（唯一基因）+ KCNJ11 ----
hits <- fread(file.path(gdir, "transcript_coloc_hits.csv"))
targets <- unique(hits[, .(gene, symbol)])
targets <- rbind(targets, data.table(gene = "ENSG00000187486", symbol = "KCNJ11"))
targets <- targets[!duplicated(gene)]
cat("目标基因数:", nrow(targets), "(106 strong 命中 + KCNJ11, 去重)\n")

# ---- 单趟流式抽取目标基因的显著 cis-eQTL ----
gene_ids <- targets$gene
tmp <- file.path(OUT, "trans_esd", "_extract.tsv")
if (file.exists(tmp)) unlink(tmp)
# awk 按列 8（Gene, ENSG）匹配；写成临时文件避免 shell 引号地狱
idf <- file.path(OUT, "trans_esd", "_ids.txt")
writeLines(gene_ids, idf)
cmd <- sprintf("zcat %s | awk -F'\\t' 'NR==FNR{a[$1]=1;next} (FNR>1) && a[$8]{print}' %s - > %s",
               shQuote(EQTL_SIG), shQuote(idf), shQuote(tmp))
cat("抽取显著 cis-eQTL（单趟）...\n")
system(cmd, wait = TRUE)
eq <- fread(tmp, sep = "\t", header = TRUE, nThread = 4)
setnames(eq, c("Pvalue","SNP","SNPChr","SNPPos","AssessedAllele","OtherAllele",
               "Zscore","Gene","GeneSymbol","GeneChr","GenePos",
               "NrCohorts","NrSamples","FDR","BonferroniP"))
cat("  ", nrow(eq), "SNP-gene 对,", uniqueN(eq$Gene), "基因\n")

# ---- AF + 1kg bim ----
af <- fread(cmd = paste0("zcat ", shQuote(AF)), sep = "\t", header = TRUE, nThread = 4)
setnames(af, c("SNP","hg19_chr","hg19_pos","AlleleA","AlleleB",
               "allA_total","allAB_total","allB_total","AlleleB_all"))
bim <- fread(paste0(LDREF, ".bim"), header = FALSE,
             col.names = c("chr", "snp", "cm", "pos", "a1", "a2"))
setkey(bim, snp)

eq2 <- merge(eq, af[, .(SNP, AlleleA, AlleleB, AlleleB_all)], by = "SNP", all.x = TRUE)
eq2[, eaf := fifelse(AssessedAllele == AlleleA, 1 - AlleleB_all,
                 fifelse(AssessedAllele == AlleleB, AlleleB_all, NA_real_))]
# 仅保留参考面板内变异（SMR 需参考 LD），且频率有效
eq2 <- merge(eq2, bim[, .(chr, snp, pos)], by.x = "SNP", by.y = "snp")
eq2 <- eq2[!is.na(eaf) & eaf > 0 & eaf < 1 & !is.na(Zscore) & NrSamples > 0]
cat("  合并 AF + 1kg 交集后:", nrow(eq2), "行\n")

# ---- 逐基因写 .esd + flist ----
eq2[, `:=`(Beta = Zscore / sqrt(NrSamples), se = 1 / sqrt(NrSamples))]
fl <- lapply(unique(eq2$Gene), function(g) {
  d <- eq2[Gene == g]
  setorder(d, Pvalue)
  esd <- d[, .(Chr = SNPChr, SNP, Bp = SNPPos,
               A1 = AssessedAllele, A2 = OtherAllele, Freq = eaf,
               Beta, se, p = Pvalue)]
  fn <- file.path(OUT, "trans_esd", paste0(g, ".esd"))
  fwrite(esd, fn, sep = " ", quote = FALSE, na = "NA")
  data.table(Chr = unique(d$SNPChr)[1], ProbeID = g, GeneticDistance = 0,
             ProbeBp = round(mean(range(d$SNPPos))), Gene = d$GeneSymbol[1],
             Orientation = "+", PathOfEsd = fn)
})
fl <- rbindlist(fl)
fwrite(fl, file.path(OUT, "trans_flist.txt"), sep = " ", quote = FALSE, na = "NA")
cat("flist:", nrow(fl), "个探针 ->", file.path(OUT, "trans_flist.txt"), "\n")

# ---- 构建 BESD ----
cmd <- sprintf("%s --eqtl-flist %s --make-besd --out %s --thread-num 4 2>&1",
               file.path(proj, "tools/smr"), file.path(OUT, "trans_flist.txt"),
               file.path(OUT, "eqtlgen_trans"))
cat("运行:", cmd, "\n")
system(cmd)

# ---- 每结局 .ma（复用 _coloc_gwas 缓存）----
gcol <- file.path(gdir, "_coloc_gwas")
for (on in names(OUT_N)) {
  fpat <- file.path(gcol, paste0(OUT_N[[on]], "_", targets$gene, ".rds"))
  ok <- file.exists(fpat)
  cat("  outcome", on, ":", sum(ok), "/", length(ok), "基因有 GWAS 缓存\n")
  parts <- lapply(targets$gene[ok], function(g) {
    d <- readRDS(file.path(gcol, paste0(OUT_N[[on]], "_", g, ".rds")))
    d[, .(snp = rsid, A1 = ea, A2 = nea, freq = eaf, b = beta, se = se, p = p, n = n)]
  })
  ma <- rbindlist(parts)
  ma <- ma[!duplicated(snp)]
  ma <- ma[A1 != A2 & !is.na(b) & !is.na(se) & se > 0 & !is.na(p)]
  fwrite(ma, file.path(OUT, paste0("trans_", on, ".ma")), sep = " ", quote = FALSE, na = "NA")
  cat("  -> trans_", on, ".ma:", nrow(ma), "SNPs\n", sep = "")
}
cat("=== 转录 SMR 输入构建完成 ===\n")
