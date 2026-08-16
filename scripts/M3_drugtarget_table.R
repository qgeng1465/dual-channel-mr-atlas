#!/usr/bin/env Rscript
# =============================================================================
# M3_drugtarget_table.R — 已知 T2D 药物靶点 × 通道可用性/四态表（P1-2，描述性）
# =============================================================================
# 输入：
#   results/grid/transcript_drugtarget_avail.csv（eQTLGen 全血转录可用性）
#   results/grid/transcript_drugtarget_mr.csv（全血转录 MR）
#   deCODE 蛋白可用性：硬编码自 data/decode 文件夹清单核验（2026-08-06）
# 输出：results/grid/drugtarget_fourstate.csv + 控制台表
# 纪律：严格描述性；不结论"成药模态"（M0 红线）；显式声明"血浆 pQTL 轴 ≠ 药物靶点蛋白轴"。
# =============================================================================
suppressMessages({library(data.table)})
proj <- "<repo-root>"
res  <- file.path(proj, "results")

a <- fread(file.path(res, "grid/transcript_drugtarget_avail.csv"))
m <- fread(file.path(res, "grid/transcript_drugtarget_mr.csv"))
outmap <- c(`ebi-a-GCST006867`="T2D",`ebi-a-GCST005194`="CAD",`ebi-a-GCST005186`="FBG")
m[, out := outmap[outcome]]

# deCODE 蛋白可用性（2026-08-06 文件夹清单核验：GLP1R/DPP4/INSR/PCK1/GCG 收录，其余否）
decode_avail <- c(GLP1R=TRUE, SLC5A2=FALSE, PPARG=FALSE, KCNJ11=FALSE, ABCC8=FALSE,
                  DPP4=TRUE, INSR=TRUE, TCF7L2=FALSE, PRKAA1=FALSE, GCG=TRUE,
                  PDX1=FALSE, SLC5A1=FALSE, G6PC=FALSE, PCK1=TRUE)
DRUG_NOTES <- c(GLP1R="GLP-1RA 靶点（受体激动剂）", SLC5A2="SGLT2i 靶点", PPARG="TZD 靶点（核受体）",
  KCNJ11="磺脲类靶点 Kir6.2", ABCC8="磺脲类靶点 SUR1", DPP4="DPP4i 靶点",
  INSR="胰岛素/增敏靶点", TCF7L2="T2D 最强易感基因", PRKAA1="AMPKα1（二甲双胍通路）",
  GCG="胰高血糖素受体", PDX1="β细胞转录因子", SLC5A1="SGLT1", G6PC="糖异生 G6Pase", PCK1="糖异生 PEPCK")

# 主方法结果（IVW-MRE/Wald），取每基因×结局最显著一行
main <- m[ok==TRUE & method %in% c("Inverse variance weighted (multiplicative random effects)","Wald ratio")]
main[, keep := pval == min(pval), by = .(gene, out)]
sig <- main[keep==TRUE, .(gene, out, nsnp, b, se, pval, method)]

# 每基因：转录通道最显著主方法（全 3 结局里 min p，方向）
best_tr <- sig[, .(n_sig = sum(pval<0.05),
                   min_p = min(pval), min_out = out[which.min(pval)],
                   min_b = b[which.min(pval)],
                   tr_signal = ifelse(any(pval<0.05), "transcript-signal", "transcript-null")), by = gene]

tab <- data.table(Gene = a$Gene)
tab[, drug := unname(DRUG_NOTES[Gene])]
tab[, transcript_cis := a$n_iv_eaf[match(Gene, a$Gene)] > 0]
tab[, decode_protein := unname(decode_avail[Gene])]
tab[, tr_signal := best_tr$tr_signal[match(Gene, best_tr$gene)]]
tab[is.na(tr_signal), tr_signal := "transcript-null"]
tab[, tr_min_p := best_tr$min_p[match(Gene, best_tr$gene)]]
tab[, tr_min_pair := best_tr$min_out[match(Gene, best_tr$gene)]]
tab[, tr_min_b := best_tr$min_b[match(Gene, best_tr$gene)]]
# 通道分类（描述性，非"成药模态"结论）：
#   both（转录+蛋白都可用）/ protein-only / transcript-only / neither
tab[, channel_class := fifelse(transcript_cis & decode_protein, "both",
                        fifelse(decode_protein, "protein-only",
                         fifelse(transcript_cis, "transcript-only", "neither")))]
setorder(tab, -transcript_cis, -decode_protein, Gene)

cat("===== 已知 T2D 药物靶点 × 通道可用性/信号（描述性）=====\n")
for (i in seq_len(nrow(tab))) {
  cat(sprintf("  %-7s %-16s 转录=%-3s 蛋白=%-3s | 通道=%s | 转录信号=%s%s\n",
    tab$Gene[i], tab$drug[i],
    ifelse(tab$transcript_cis[i],"YES","no"), ifelse(tab$decode_protein[i],"YES","no"),
    tab$channel_class[i], tab$tr_signal[i],
    ifelse(!is.na(tab$tr_min_p[i]), sprintf(" (min p=%s, %s, b=%s)", formatC(tab$tr_min_p[i],2,format="g"), tab$tr_min_pair[i], formatC(tab$tr_min_b[i],3,format="f")), "")))
}
cat("\n注：转录信号来自 eQTLGen 全血 cis-eQTL（p<5e-6, LD clump）→ T2D/CAD/FBG 主方法 IVW-MRE/Wald。\n")
cat("    蛋白通道 MR 待 deCODE 数据（当前下载受服务器 Range 限制阻塞）。\n")
cat("    血浆 pQTL 轴 ≠ 药物靶点蛋白轴；本表为观察性通道可用性，非成药建议（M0 红线）。\n")

write.csv(tab, file.path(res, "grid/drugtarget_fourstate.csv"), row.names = FALSE)
cat("✔ 落盘: results/grid/drugtarget_fourstate.csv\n")
