#!/usr/bin/env Rscript
# =============================================================================
# M12_lipid_validation.R — 杠杆 B：脂质结局 MR 验证（探索性，非预注册主流程）
# =============================================================================
# 动机（2026-08-13，用户选定 A+B 并进）：
#   蛋白通道（deCODE pQTL）在 T2D/CAD/FBG 下仅 2 对 coloc strong，且阳性对照 PCSK9×CAD
#   过不了 HEIDI → 蛋白侧未校准。但 PCSK9/APOC3/ANGPTL3/LDLR/APOB/HMGCR 都是脂质生物学
#   教科书因果基因：若 deCODE 蛋白工具在 LDL/HDL/TG/apoB 结局下能恢复已知脂质效应
#   （方向正确 + 强显著），则证明蛋白工具本身可用，"蛋白通道在心血管结局未达标"是
#   结局侧而非工具侧问题 → 决定蛋白通道能否被救活/如何表述。
# 结局（Richardson 2020 UKB，欧洲，与 deCODE 冰岛 n=35,559 无样本重叠 → 干净双样本）：
#   LDL  ieu-b-110 (n=440,546) / HDL ieu-b-109 (n=403,943)
#   TG   ieu-b-111 (n=441,016) / apoB ieu-b-108 (n=439,214)
# 方法：与 M3 完全一致 —— cis ±1Mb(TSS)、p<5e-6、LD clump r²<0.01@1000kb EUR(API)、
#   主 mr_ivw_mre、FE/WM/Egger 敏感性、nsnp=1 → Wald；deCODE 列名/多等位排除/annotated
#   eaf 增强/ImpMAF 近似照搬。
# 纪律：探索性扩展，独立输出 CSV（不触碰 protein_decode_mr.csv 主结果）；
#   全网格落盘含空结果；显著命中按已知脂质生物学方向标注校准（非新发现声称）。
# 用法：PATH=$R_ENV/bin:$PATH Rscript scripts/M12_lipid_validation.R
# 输出：results/grid/protein_decode_lipid_mr.csv
# =============================================================================
suppressMessages({
  library(data.table)
  suppressPackageStartupMessages(library(TwoSampleMR))
  suppressPackageStartupMessages(library(ieugwasr))
  library(jsonlite)
})

proj <- "<repo-root>"
res  <- file.path(proj, "results")
SUB  <- file.path(proj, "data/decode/sub")
dir.create(file.path(res, "grid"), recursive = TRUE, showWarnings = FALSE)
log <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n")

cfg <- read_json(file.path(res, "config.json"), auto_unbox = TRUE)
log("M12 脂质验证（探索性）| 主方法=", cfg$mr_methods$primary)

# --- 蛋白列表 + hg38 坐标（与 M3 一致）-----------------------------------------
PROTEINS <- list(
  list(file = "5231_79_PCSK9_PCSK9.txt.gz",   gene = "PCSK9",   chr = "1",  tss = 55039445, note = "positive control"),
  list(file = "5230_99_HMGCR_HMGR.txt.gz",    gene = "HMGCR",   chr = "5",  tss = 75336329, note = "negative control"),
  list(file = "10391_1_ANGPTL3_ANGL3.txt.gz", gene = "ANGPTL3", chr = "1",  tss = 62597464, note = "negative control"),
  list(file = "6461_54_APOC3_Apo_C_III.txt.gz", gene = "APOC3", chr = "11", tss = 116827019, note = "negative control"),
  list(file = "2797_56_APOB_Apo_B.txt.gz",    gene = "APOB",    chr = "2",  tss = 21044075, note = "optional"),
  list(file = "13129_40_LDLR_LDLR.txt.gz",    gene = "LDLR",    chr = "19", tss = 11089418, note = "optional"),
  list(file = "13085_18_GLP1R_GLP1R.txt.gz",  gene = "GLP1R",   chr = "6",  tss = 39048562, note = "drug target GLP-1RA"),
  list(file = "15460_9_DPP4_CD26.txt.gz",     gene = "DPP4",    chr = "2",  tss = 162074639, note = "drug target DPP4i"),
  list(file = "3448_13_INSR_IR.txt.gz",       gene = "INSR",    chr = "19", tss = 7294443, note = "drug target insulin/insulin sens."),
  list(file = "18182_24_PCK1_PCKGC.txt.gz",   gene = "PCK1",    chr = "20", tss = 57546220, note = "drug target gluconeogenesis"),
  list(file = "4891_50_GCG_Glucagon.txt.gz",  gene = "GCG",     chr = "2",  tss = 162152404, note = "drug target glucagon")
)
LIPIDS <- list(ldl = "ieu-b-110", hdl = "ieu-b-109", tg = "ieu-b-111", apob = "ieu-b-108")
LIPID_LABEL <- c(ldl = "LDL-C", hdl = "HDL-C", tg = "Triglycerides", apob = "apoB")

# --- 已知脂质生物学校准方向（文献确定，用于标注"工具正确回收"，非新发现）------
#   MR b 的方向 = 蛋白水平↑(per-SD) 对结局的因果效应方向（不依赖工具是 LOF/GOF）。
#   "up"=蛋白↑→结局↑（b 应>0）、"down"=蛋白↑→结局↓（b 应<0）。
#   来源：PCSK9i/volanesorsen/evinacumab/他汀/FH 等既有干预表型。
EXPECTED_DIR <- list(
  PCSK9   = c(ldl = "up", apob = "up", tg = "up"),                      # PCSK9↑→LDL/apoB/TG↑（LOF 全降）
  APOC3   = c(tg = "up", ldl = "up", hdl = "down"),                     # APOC3↑→TG/LDL↑、HDL↓（volanesorsen 表型）
  ANGPTL3 = c(tg = "up", ldl = "up", hdl = "down"),                     # ANGPTL3↑→TG/LDL↑、HDL↓（evinacumab 表型）
  LDLR    = c(ldl = "down", apob = "down"),                             # LDLR↑→LDL/apoB↓（受体清除）
  APOB    = c(ldl = "up", apob = "up"),                                 # APOB↑→LDL/apoB↑（结构蛋白）
  HMGCR   = c(ldl = "up", apob = "up"),                                 # HMGCR↑→LDL/apoB↑（胆固醇合成）
  GLP1R   = c(tg = "down", ldl = "down")                                # GLP1R↑(激活)→TG/LDL↓（GLP-1RA 表型）
)
# 注：PCSK9→HDL 临床试验基本中性，不列入校准；APOC3→HDL 仅作参考（volanesorsen 升 HDL）。

# --- annotated eaf 增强 + palindromic 处理（照搬 M3）---------------------------
ann_path <- file.path(proj, "data/decode/assocvariants.annotated.txt.gz")
ANN <- NULL
if (file.exists(ann_path) &&
    system(paste("gzip -t", shQuote(ann_path), "2>/dev/null"),
           ignore.stdout = TRUE, ignore.stderr = TRUE) == 0) {
  ANN <- fread(cmd = paste0("zcat ", shQuote(ann_path)), sep = "\t", header = TRUE, nThread = 4)
  setnames(ANN, c("Chrom","Pos","Name","rsids","effectAllele","otherAllele","effectAlleleFreq"))
  log("annotated 已加载（effectAlleleFreq）")
} else log("⚠ annotated 未就绪 → palindromic 排除（保守），ImpMAF 近似")

is_pal <- function(a, b) {
  (a == "A" & b == "T") | (a == "T" & b == "A") |
  (a == "C" & b == "G") | (a == "G" & b == "C")
}

load_cis <- function(prot) {
  sub_f <- file.path(SUB, paste0(sub("\\.gz$", "", prot$file), "_cis.txt.gz"))
  full_f <- file.path(proj, "data/decode", prot$file)
  f <- if (file.exists(sub_f) && file.size(sub_f) > 1000) sub_f
       else if (file.exists(full_f) && file.size(full_f) > 1e7) full_f else NULL
  if (is.null(f)) return(NULL)
  d <- fread(cmd = paste0("zcat ", shQuote(f)), sep = "\t", header = TRUE, nThread = 4)
  if (!"Pval" %in% names(d) && "minus_log10_pval" %in% names(d)) d[, Pval := 10^(-minus_log10_pval)]
  if (!"Pval" %in% names(d)) return(NULL)
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

run_lipid_mr <- function(prot, out_id, out_label) {
  d <- tryCatch(load_cis(prot), error = function(e) { log("  [读取失败] ", prot$gene, ": ", conditionMessage(e)); NULL })
  if (is.null(d))
    return(data.frame(gene = prot$gene, outcome = out_label, outid = out_id, nsnp = 0,
                      b = NA, se = NA, pval = NA, method = NA, ok = FALSE, note = "蛋白文件未下载/不完整/读取失败"))
  d <- d[Chrom2 == prot$chr & abs(Pos - prot$tss) <= 1e6 & Pval < cfg$instrument$pval_thresh]
  if (nrow(d) == 0)
    return(data.frame(gene = prot$gene, outcome = out_label, outid = out_id, nsnp = 0,
                      b = NA, se = NA, pval = NA, method = NA, ok = FALSE, note = "cis ±1Mb、p<5e-6 内无变异"))
  n_ma <- d[, sum(effectAllele == otherAllele | otherAllele %in% c("!", "", NA), na.rm = TRUE)]
  d <- d[effectAllele != otherAllele & !otherAllele %in% c("!", "", NA)]
  d[, rsid2 := fifelse(is.na(rsids) | rsids %in% c("", ".", "-"), Name, rsids)]
  d <- d[!duplicated(rsid2)]
  d[, eaf := ImpMAF]; d[, eaf_src := "ImpMAF 近似"]
  if (!is.null(ANN)) {
    setkey(ANN, Name); d[ANN, c("eaf", "eaf_src") := .(i.effectAlleleFreq, "annotated")]
    d[is.na(eaf), eaf := ImpMAF]
  }
  npal <- sum(is_pal(d$effectAllele, d$otherAllele))
  if (is.null(ANN)) {
    d <- d[!is_pal(effectAllele, otherAllele)]
    note_eaf <- paste0("palindromic 排除 ", npal, "（无 annotated）")
  } else {
    d <- d[!(is_pal(effectAllele, otherAllele) & is.na(eaf))]
    note_eaf <- paste0("palindromic 保留（annotated），排除无频率 ", sum(is.na(d$eaf)), " 个")
  }
  if (nrow(d) == 0)
    return(data.frame(gene = prot$gene, outcome = out_label, outid = out_id, nsnp = 0,
                      b = NA, se = NA, pval = NA, method = NA, ok = FALSE, note = paste0(note_eaf, " → 无工具")))
  cl <- tryCatch(ld_clump(data.frame(rsid = d$rsid2, pval = d$Pval, id = prot$gene),
                          clump_kb = cfg$instrument$clump_kb, clump_r2 = cfg$instrument$clump_r2,
                          clump_p = 1, pop = "EUR",
                          opengwas_jwt = ieugwasr::get_opengwas_jwt()),
                 error = function(e) NULL)
  if (is.null(cl) || nrow(cl) == 0)
    return(data.frame(gene = prot$gene, outcome = out_label, outid = out_id, nsnp = 0,
                      b = NA, se = NA, pval = NA, method = NA, ok = FALSE,
                      note = "clump 后无独立工具"))
  if (nrow(cl) >= nrow(d) && nrow(d) <= 30) {
    lk <- tryCatch(ieugwasr::ld_reflookup(d$rsid2, pop = "EUR",
                                          opengwas_jwt = ieugwasr::get_opengwas_jwt()),
                   error = function(e) NULL)
    if (!is.null(lk) && !any(d$rsid2 %in% lk))
      return(data.frame(gene = prot$gene, outcome = out_label, outid = out_id, nsnp = 0,
                        b = NA, se = NA, pval = NA, method = NA, ok = FALSE,
                        note = "cis 变异均不在 1000G EUR LD 面板（罕见位点）→ 无工具"))
  }
  d <- d[d$rsid2 %in% cl$rsid, ]
  exp_dat <- data.frame(
    SNP = d$rsid2, effect_allele.exposure = d$effectAllele, other_allele.exposure = d$otherAllele,
    eaf.exposure = d$eaf, beta.exposure = d$Beta, se.exposure = d$SE,
    pval.exposure = d$Pval, samplesize.exposure = d$N,
    id.exposure = prot$gene, exposure = prot$gene, stringsAsFactors = FALSE)
  tryCatch({
    out <- extract_outcome_data(snps = exp_dat$SNP, outcomes = out_id, proxies = TRUE)
    dat <- harmonise_data(exp_dat, out, action = 2)
    dat <- dat[dat$mr_keep, ]
    if (nrow(dat) == 0)
      return(data.frame(gene = prot$gene, outcome = out_label, outid = out_id, nsnp = nrow(d),
                        b = NA, se = NA, pval = NA, method = NA, ok = FALSE, note = "harmonise 后无保留 SNP"))
    r_main <- mr(dat, method_list = cfg$mr_methods$primary)
    r_sens <- tryCatch(mr(dat, method_list = cfg$mr_methods$sensitivity), error = function(e) NULL)
    r_all <- if (!is.null(r_sens)) rbind(r_main, r_sens) else r_main
    if (nrow(r_all) == 0 && nrow(dat) == 1) {
      r_w <- tryCatch(mr(dat, method_list = "mr_wald_ratio"), error = function(e) NULL)
      if (!is.null(r_w) && nrow(r_w) > 0) r_all <- r_w
    }
    if (nrow(r_all) == 0)
      return(data.frame(gene = prot$gene, outcome = out_label, outid = out_id, nsnp = nrow(dat),
                        b = NA, se = NA, pval = NA, method = NA, ok = FALSE, note = "MR 无输出"))
    note_final <- paste0(ifelse(any(d$eaf_src == "ImpMAF 近似"), note_eaf, "eaf=annotated"),
                         "; clump ", nrow(exp_dat), "→harmonise ", nrow(dat))
    data.frame(gene = prot$gene, outcome = out_label, outid = out_id, nsnp = nrow(dat),
               b = r_all$b, se = r_all$se, pval = r_all$pval, method = r_all$method,
               ok = TRUE, note = note_final)
  }, error = function(e) data.frame(gene = prot$gene, outcome = out_label, outid = out_id, nsnp = 0,
                                    b = NA, se = NA, pval = NA, method = NA, ok = FALSE,
                                    note = conditionMessage(e)))
}

out <- list()
for (p in PROTEINS) {
  if (!file.exists(file.path(proj, "data/decode", p$file)) && !file.exists(file.path(SUB, paste0(p$file, "_cis.txt.gz")))) {
    log("跳过 ", p$gene, "：文件未就绪")
    for (on in names(LIPIDS))
      out[[length(out) + 1]] <- data.frame(gene = p$gene, outcome = LIPID_LABEL[[on]], outid = LIPIDS[[on]],
                                           nsnp = 0, b = NA, se = NA, pval = NA, method = NA, ok = FALSE,
                                           note = "蛋白文件未下载/不完整")
    next
  }
  for (on in names(LIPIDS)) {
    log("=== ", p$gene, " × ", LIPID_LABEL[[on]], " (", LIPIDS[[on]], ") ===")
    rr <- run_lipid_mr(p, LIPIDS[[on]], LIPID_LABEL[[on]])
    out[[length(out) + 1]] <- rr
    for (i in seq_len(nrow(rr)))
      if (rr$ok[i]) log("  ", p$gene, "×", LIPID_LABEL[[on]], "[", rr$method[i], "] nsnp=", rr$nsnp[i],
                        " b=", round(rr$b[i], 3), " p=", format(rr$pval[i], digits = 2))
  }
}
res_df <- as.data.table(do.call(rbind, out))
# 校准标注：按 EXPECTED_DIR 标"方向符合已知生物/非预期"
res_df[, calib := ""]
for (g in names(EXPECTED_DIR)) {
  for (on in names(EXPECTED_DIR[[g]])) {
    lbl <- LIPID_LABEL[[on]]; want <- EXPECTED_DIR[[g]][[on]]
    sig <- res_df[gene == g & outcome == lbl & ok == TRUE]
    for (i in seq_len(nrow(sig))) {
      if (is.na(sig$b[i])) next
      got <- ifelse(sig$b[i] > 0, "up", "down")
      res_df[gene == g & outcome == lbl & method == sig$method[i] & !is.na(b),
             calib := ifelse(got == want, paste0("方向 ✓ (蛋白", want, "→结局 ", lbl, "，符合已知)"),
                             paste0("方向 ✗ (蛋白", want, "→结局 ", lbl, "，与已知 ", want, " 相反)"))]
    }
  }
}
write.csv(res_df, file.path(res, "grid/protein_decode_lipid_mr.csv"), row.names = FALSE)
log("脂质验证完成 ✔ → results/grid/protein_decode_lipid_mr.csv (", nrow(res_df), " 行，探索性)")
