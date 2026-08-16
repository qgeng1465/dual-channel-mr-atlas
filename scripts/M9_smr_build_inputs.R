#!/usr/bin/env Rscript
# =============================================================================
# M9_smr_build_inputs.R — 构建 SMR/HEIDI 输入（deCODE pQTL BESD + GWAS .ma）
# =============================================================================
# 目标：对蛋白通道 5 个有 cis 工具的蛋白（PCSK9/APOC3/INSR/ANGPTL3/PCK1）× 3 结局
#   （T2D/CAD/FBG）跑 SMR + HEIDI（多效性 vs LD 裁决，Zhu 2016 Nat Genet 方法）。
#
# 设计（2026-08-13）：
#   * deCODE 为 hg38，参考面板 1kg 为 hg19 → 写 .esd 时用 1kg bim 的 hg19 位置（rsid 匹配），
#     探针位置用蛋白侧 GWAS 缓存（hg19 ±1Mb）中心（已与 M5 已知 TSS 验证，偏差 <200bp）。
#   * deCODE 无 effectAlleleFreq（annotated 文件未下）→ .esd Freq 用 ImpMAF 近似，SMR 用
#     --disable-freq-ck 跳过频率 QC，如实记录（探索性敏感性分析）。
#   * GWAS .ma 用 GCTA-COJO 格式（SNP A1 A2 freq b se p n），来自蛋白侧缓存（ea=eaf 真值）。
#   * SMR 用 rsid 在 BESD/GWAS/参考面板之间对齐，build 无关。
#
# 产物：
#   data/smr/esd/{PROT}.esd, data/smr/flist.txt → data/smr/decode_pqtl.{besd,esi,epi}
#   data/smr/{t2d,cad,fbg}.ma
# 用法：PATH=$R_ENV/bin:$PATH Rscript scripts/M9_smr_build_inputs.R
# =============================================================================
suppressMessages(library(data.table))
proj <- "/data/qiushuogeng/projects/dual-channel-mr-atlas"
gdir <- file.path(proj, "results/grid"); GPROT <- file.path(proj, "results/grid/_coloc_gwas_prot")
SUB  <- file.path(proj, "data/decode/sub")
LDREF <- file.path(proj, "data/ldref/1kg.v3/EUR")
SMR <- file.path(proj, "tools/smr")
OUT <- file.path(proj, "data/smr")
dir.create(file.path(OUT, "esd"), recursive = TRUE, showWarnings = FALSE)

# ---- 5 蛋白 × 3 结局定义（探针位置 = 蛋白侧 GWAS 缓存中心，hg19）----
OUT_N <- c(t2d = "GCST006867", cad = "GCST005194", fbg = "GCST005186")
PROTEINS <- data.table(
  gene = c("PCSK9", "APOC3", "INSR", "ANGPTL3", "PCK1"),
  chr  = c("1", "11", "19", "1", "20"),
  file = c("5231_79_PCSK9_PCSK9.txt.gz", "6461_54_APOC3_Apo_C_III.txt.gz",
           "3448_13_INSR_IR.txt.gz", "10391_1_ANGPTL3_ANGL3.txt.gz",
           "18182_24_PCK1_PCKGC.txt.gz"))
for (i in seq_len(nrow(PROTEINS))) {
  g <- readRDS(file.path(GPROT, paste0("GCST005194_", PROTEINS$gene[i], ".rds")))
  PROTEINS$probe_bp[i] <- round(mean(range(g$position)))
}
cat("探针位置（hg19，缓存中心）:\n"); print(PROTEINS[, .(gene, chr, probe_bp)])

# ---- 1kg bim（hg19 位置 + 等位基因）----
bim <- fread(paste0(LDREF, ".bim"), header = FALSE,
             col.names = c("chr", "snp", "cm", "pos", "a1", "a2"))
setkey(bim, snp)

# ---- 1) 构建 GWAS .ma（每结局 1 个，GCTA-COJO 格式）----
for (on in names(OUT_N)) {
  parts <- lapply(PROTEINS$gene, function(g) {
    d <- readRDS(file.path(GPROT, paste0(OUT_N[[on]], "_", g, ".rds")))
    d[, .(snp = rsid, A1 = ea, A2 = nea, freq = eaf, b = beta, se = se, p = p, n = n)]
  })
  ma <- rbindlist(parts)
  ma <- ma[!duplicated(snp)]            # 基因区间不重叠，保险去重
  ma <- ma[A1 != A2 & !is.na(b) & !is.na(se) & se > 0 & !is.na(p)]
  fwrite(ma, file.path(OUT, paste0(on, ".ma")), sep = " ", quote = FALSE, na = "NA")
  cat("GWAS .ma", on, ":", nrow(ma), "SNPs ->", file.path(OUT, paste0(on, ".ma")), "\n")
}

# ---- 2) 构建 per-protein .esd + .flist ----
load_pqtl <- function(prot) {
  f <- file.path(SUB, paste0(sub("\\.gz$", "", prot$file), "_cis.txt.gz"))
  d <- fread(cmd = paste0("zcat ", shQuote(f)), sep = "\t", header = TRUE, nThread = 4)
  d <- d[effectAllele != otherAllele & !otherAllele %in% c("!", "", NA)]
  d[, rsid2 := fifelse(is.na(rsids) | rsids %in% c("", ".", "-"),
                       Name, sub("[,; ].*$", "", rsids))]
  d <- d[!duplicated(rsid2)]
  d
}
flist <- lapply(seq_len(nrow(PROTEINS)), function(i) {
  g <- PROTEINS$gene[i]
  d <- load_pqtl(PROTEINS[i])
  # rsid 匹配 1kg（hg19 位置）；仅保留参考面板有的变异（SMR 需参考算 LD）
  d2 <- merge(d, bim[, .(chr, snp, pos)], by.x = "rsid2", by.y = "snp")
  d2 <- d2[!is.na(Beta) & !is.na(SE) & SE > 0 & !is.na(Pval)]
  d2[, `:=`(A1 = effectAllele, A2 = otherAllele, Freq = ImpMAF)]  # Freq 近似（ImpMAF 不总是 effect allele 频率），用 --disable-freq-ck
  esd <- d2[, .(Chr = sub("^chr", "", Chrom), SNP = rsid2, Bp = pos,
                A1, A2, Freq, Beta, se = SE, p = Pval)]
  esd <- esd[Freq > 0 & Freq < 1]
  fn <- file.path(OUT, "esd", paste0(g, ".esd"))
  fwrite(esd, fn, sep = " ", quote = FALSE, na = "NA")
  cat("  .esd", g, ":", nrow(esd), "SNPs (1kg 交集)\n")
  data.table(Chr = PROTEINS$chr[i], ProbeID = g, GeneticDistance = 0,
             ProbeBp = PROTEINS$probe_bp[i], Gene = g, Orientation = "+",
             PathOfEsd = fn)
})
fl <- rbindlist(flist)
fwrite(fl, file.path(OUT, "flist.txt"), sep = " ", quote = FALSE, na = "NA")
cat("flist 写好了:", file.path(OUT, "flist.txt"), "\n")

# ---- 3) 构建 BESD ----
cmd <- sprintf("%s --eqtl-flist %s --make-besd --out %s --thread-num 4 2>&1",
               SMR, file.path(OUT, "flist.txt"), file.path(OUT, "decode_pqtl"))
cat("运行:", cmd, "\n")
rc <- system(cmd)
if (rc != 0) stop("SMR --make-besd 失败，exit=", rc)   # 失败即停，勿谎报完成
cat("=== BESD 构建完成 →", file.path(OUT, "decode_pqtl.besd"), "===\n")
