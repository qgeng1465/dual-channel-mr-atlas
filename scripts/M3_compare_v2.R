#!/usr/bin/env Rscript
# =============================================================================
# M3_compare_v2.R — 转录本通道 MR v1(QA) vs v2(修正) 对比
# v1: results/grid/transcript_mr_qa.csv  (eaf=0.5 占位, 无 LD clump)
# v2: results/grid/transcript_mr_v2.csv   (真实 eaf, EUR LD clump r²<0.01@1000kb)
# 输出: results/grid/compare_transcript_v1v2.csv + 控制台摘要
# =============================================================================
suppressMessages({library(data.table)})
proj <- "/data/qiushuogeng/projects/dual-channel-mr-atlas"
res  <- file.path(proj, "results")

f1 <- file.path(res, "grid/transcript_mr_qa.csv")
f2 <- file.path(res, "grid/transcript_mr_v2.csv")
stopifnot(file.exists(f1), file.exists(f2))

# 仅取主方法 IVW MRE（两版一致）
ivw1 <- fread(f1)[method == "Inverse variance weighted (multiplicative random effects)"]
ivw2 <- fread(f2)[method == "Inverse variance weighted (multiplicative random effects)"]
out_name <- setNames(c("t2d","cad","fbg"),
                     c("ebi-a-GCST006867","ebi-a-GCST005194","ebi-a-GCST005186"))
ivw1[, outcome := fifelse(outcome %in% names(out_name), unname(out_name[outcome]), outcome)]
ivw2[, outcome := fifelse(outcome %in% names(out_name), unname(out_name[outcome]), outcome)]
setnames(ivw1, c("b","se","pval"), c("b1","se1","p1"))
setnames(ivw2, c("b","se","pval"), c("b2","se2","p2"))

cmp <- merge(ivw1[, .(gene, outcome, nsnp1 = nsnp, b1, se1, p1)],
             ivw2[, .(gene, outcome, nsnp2 = nsnp, b2, se2, p2)],
             by = c("gene","outcome"), all = TRUE)
cmp[, sign_flip := (sign(b1) != sign(b2) & !is.na(b1) & !is.na(b2))]
cmp[, p_change := fifelse(is.na(p1) | is.na(p2), NA_real_, log10(p1) - log10(p2))]

cat("===== v1 vs v2 转录本通道 MR 对比（主方法 IVW-MRE）=====\n")
cat("v1 行数:", nrow(ivw1), " | v2 行数:", nrow(ivw2),
    " | 成功对数:", sum(!is.na(cmp$b1) & !is.na(cmp$b2)), "\n\n")

key <- c("CETP","HMGCR","PCSK9","ANGPTL3","APOC3","NPC1L1","LDLR","APOB","A2MP1")
tab <- cmp[gene %in% key & outcome == "t2d", .(gene, nsnp1, nsnp2, b1, b2, p1, p2)]
setorder(tab, gene)
print(tab, digits = 3)
cat("\n符号翻转对数:", cmp[sign_flip == TRUE, .N],
    " | nsnp 缩减（clump 生效）:", cmp[nsnp2 < nsnp1 & !is.na(nsnp1), .N], "\n")
cat("极端 p（p2<1e-50）v1:", ivw1[p1 < 1e-50, .N], " → v2:", ivw2[p2 < 1e-50, .N], "\n")

# HMGCR×CAD 关键核对（v1 阴性 p=0.93 与已知 LDL-C→CAD 不一致）
cat("\n—— HMGCR×CAD ——\n")
print(cmp[gene=="HMGCR" & outcome=="cad", .(gene, outcome, nsnp1, nsnp2, b1, b2, p1, p2)], digits = 4)

write.csv(cmp, file.path(res, "grid/compare_transcript_v1v2.csv"), row.names = FALSE)
cat("\n✔ 对比表已落盘: results/grid/compare_transcript_v1v2.csv\n")
