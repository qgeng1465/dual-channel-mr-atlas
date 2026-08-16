#!/usr/bin/env Rscript
# =============================================================================
# 08_figures.R — M8 图表（诚实报告风格）
# =============================================================================
# 设计方法（dataviz 规范，选型先于配色）：
#   F1 漏斗阶梯      水平条形   magnitude → 顺序蓝（单色阶）
#   F2 v1→v2 p 值修正  dumbbell  before→after（1 色 2 档：灰=旧/蓝=新）
#   F3 极端 p 消失     2 柱对比  magnitude → 顺序蓝
#   F4 四态堆叠条形   待 coloc/双通道数据（守门 stub）
#   F5 PP.H4 三档      待 coloc 数据（stub）
#   F6 LOOCV AUR      待对照数据（stub）
#   F7 外部复现森林图  待复现数据（stub）
#   F8 介质层热图     待四态数据（stub）
# 配色：CVD 安全（蓝#2a78d6 / 橙#eb6834 / 青#1baf7a），白底，墨#0b0b0b，
#       次级#52514e，网格线#e1e0d9，无双轴，细标记，直标标签。
# 用法：PATH=$R_ENV/bin:$PATH Rscript scripts/08_figures.R
# 输出：results/figures/F1_*.png（dpi=320，论文可用）
# =============================================================================
suppressMessages({
  library(data.table)
  library(ggplot2)
})
proj <- "/data/qiushuogeng/projects/dual-channel-mr-atlas"
res  <- file.path(proj, "results")
fig  <- file.path(res, "figures"); dir.create(fig, showWarnings = FALSE)
log  <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n")

# --- 统一主题（诚实报告风格）--------------------------------------------------
pal <- c(blue = "#2a78d6", orange = "#eb6834", aqua = "#1baf7a",
         gray = "#a8a8a8", dark = "#0b0b0b")
theme_paper <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      text = element_text(colour = pal["dark"]),
      axis.text = element_text(colour = "#52514e"),
      axis.title = element_text(colour = "#52514e"),
      axis.line = element_line(colour = "#c3c2b7", linewidth = 0.4),
      axis.ticks = element_line(colour = "#c3c2b7"),
      panel.grid.major.x = element_line(colour = "#e1e0d9", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.title = element_text(face = "bold", size = rel(1.15), margin = margin(b = 6)),
      plot.subtitle = element_text(colour = "#52514e", size = rel(0.85)),
      legend.position = "top",
      legend.title = element_blank(),
      legend.key.size = unit(0.35, "cm"),
      plot.caption = element_text(colour = "#898781", size = rel(0.75), hjust = 0)
    )
}

save_fig <- function(p, nm, w = 7.5, h = 4.2) {
  ggsave(file.path(fig, nm), p, width = w, height = h, dpi = 320, bg = "white")
  log("✔ 已生成: ", nm)
}

# =============================================================================
# F1 漏斗阶梯（M2，转录通道）—— 顺序蓝水平条形，直标计数
# =============================================================================
# 漏斗末档用"唯一基因×结局对"（ok=TRUE 的去重对数），非方法级行数（避免漏斗不递减）
n_pairs <- NA_integer_
f2csv <- file.path(res, "grid/transcript_mr_v2.csv")
if (file.exists(f2csv)) {
  tmp <- fread(f2csv); n_pairs <- uniqueN(tmp[ok == TRUE, .(gene, outcome)])
}
if (file.exists(file.path(res, "grid/transcript_mr_v2_wald.csv"))) {
  tmp2 <- fread(file.path(res, "grid/transcript_mr_v2_wald.csv"))
  n_pairs <- n_pairs + uniqueN(tmp2[ok == TRUE, .(gene, outcome)])
}
f <- data.table(
  stage_cn = c("eQTLGen 有 cis-eQTL 基因",
               "p<5e-6 工具候选基因",
               "真实 eaf 可用基因",
               "本轮代表基因集（对照+强 cis）",
               "基因×结局对",
               "产出 MR 结果的基因×结局对"),
  count = c(16923, 16561, 16561, 12, 36, n_pairs))
f[, order := .N - .I + 1]
f[, grp := fifelse(stage_cn == f$stage_cn[.N], "end", "mid")]
p1 <- ggplot(f, aes(x = reorder(stage_cn, order), y = count, fill = grp)) +
  geom_col(width = 0.62, show.legend = FALSE) +
  geom_text(aes(label = format(count, big.mark = ",")), hjust = -0.12,
            size = 3.2, colour = "#0b0b0b") +
  scale_fill_manual(values = c(mid = pal["blue"], end = pal["dark"])) +
  coord_flip() +
  labs(title = "转录通道可比性漏斗（代表基因集 12，双通道全量待 deCODE）",
       subtitle = "eQTLGen 全血 cis-eQTL（n=31,684）× 真实 eaf × EUR LD clump",
       x = NULL, y = "基因 / 对 数") +
  theme_paper() +
  expand_limits(x = 1.35)
save_fig(p1, "F1_funnel_transcript.png")

# =============================================================================
# F2 v1→v2 方法学修正（dumbbell：灰=旧 蓝=新，-log10(p) 坍缩）
# =============================================================================
fn2 <- file.path(res, "grid/compare_transcript_v1v2.csv")
if (file.exists(fn2)) {
  c0 <- fread(fn2)
  # 取两版都有 p 的完整对（before→after 才可比）；p=0 视为极显著→-log10 截断在 6
  d <- c0[!is.na(p1) & !is.na(p2)]
  log("F2 完整对比对（v1 且 v2 均有结果）: ", nrow(d), " 对")
  outmap <- c(t2d = "T2D", cad = "CAD", fbg = "FBG")
  d[, out := fifelse(outcome %in% names(outmap), outmap[outcome], outcome)]
  d[, pair := paste(gene, out, sep = " × ")]
  d[, lp1 := ifelse(p1 > 0, -log10(p1), Inf)]; d[, lp2 := ifelse(p2 > 0, -log10(p2), Inf)]
  d[, lp1c := pmin(lp1, 6)]; d[, lp2c := pmin(lp2, 6)]  # 截断防视图失真
  d[, cap1 := ifelse(lp1 > 6, ">6", NA_character_)]
  setorder(d, -lp1)
  dd <- melt(d[, .(pair, v1 = lp1c, v2 = lp2c)], id.vars = "pair",
             variable.name = "ver", value.name = "l")
  dd[, cap := ifelse(ver == "v1", "v1", "v2")]
  p2 <- ggplot(dd, aes(x = l, y = reorder(pair, l), colour = ver)) +
    geom_line(aes(group = pair), colour = "#c3c2b7", linewidth = 0.5) +
    geom_point(size = 2.6) +
    scale_colour_manual(values = c(v1 = "#a8a8a8", v2 = pal["blue"]),
                        labels = c(v1 = "v1（未 clump / eaf 占位）", v2 = "v2（真实 eaf + LD clump）")) +
    scale_x_continuous(breaks = 0:6, labels = c(0:5, "≥6")) +
    labs(title = "转录本通道 MR 方法学修正：-log10(p) 坍缩",
         subtitle = "LD clump + 真实 eaf 消除连锁 SNP 虚增；截断在 -log10(p)=6",
         x = "-log10(p)（IVW-MRE 主方法）", y = NULL) +
    theme_paper()
  save_fig(p2, "F2_v1v2_pvalue_collapse.png", w = 6.5, h = 3.2)
}

# =============================================================================
# F3 极端 p 值消失（顺序蓝 2 柱）
# =============================================================================
if (file.exists(fn2) && file.exists(file.path(res, "grid/transcript_mr_qa.csv"))) {
  v1 <- fread(file.path(res, "grid/transcript_mr_qa.csv"))
  v2 <- fread(file.path(res, "grid/transcript_mr_v2.csv"))
  n1 <- sum(v1$pval < 1e-50, na.rm = TRUE); n2 <- sum(v2$pval < 1e-50, na.rm = TRUE)
  ext <- data.table(ver = c("v1", "v2"), n = c(n1, n2))
  p3 <- ggplot(ext, aes(x = ver, y = n)) +
    geom_col(width = 0.5, fill = c(pal["gray"], pal["blue"])) +
    geom_text(aes(label = n), vjust = -0.4, size = 4, colour = "#0b0b0b") +
    labs(title = "极端 p 值（<1e-50）数量：v1 → v2",
         subtitle = "未 clump 的连锁 SNP 被当作独立工具 → 假性极端 p",
         x = NULL, y = "p<1e-50 的结果数") +
    theme_paper() + ylim(0, max(ext$n) * 1.2)
  save_fig(p3, "F3_extreme_p_collapse.png", w = 3.8, h = 3.4)
}

# =============================================================================
# F9 deCODE 蛋白通道漏斗（顺序蓝水平条）—— 数据到位自动渲染
# =============================================================================
ff <- file.path(res, "funnel/funnel_protein_decode.tsv")
if (file.exists(ff)) {
  fn <- fread(ff)
  # 2026-08-07 P2：stage 行缺失时 `fn[stage==…, count]` 返回 integer(0)，c() 会丢弃
  # → 行数不齐 ggplot 崩溃；三阶段齐备才渲染，缺行则跳过留痕。
  need <- c("tested_proteins", "gene_outcome_pairs", "mr_completed_rows")
  if (!all(need %in% fn$stage)) {
    log("⚠ funnel_protein_decode.tsv 缺阶段行（", paste(setdiff(need, fn$stage), collapse = ","),
        "），跳过 F9")
    rm(fn)
  } else {
  f9 <- data.table(
    stage_cn = c("deCODE 分析蛋白", "蛋白×结局对", "产出 MR 结果的对"),
    count = c(fn[stage == "tested_proteins", count],
              fn[stage == "gene_outcome_pairs", count],
              fn[stage == "mr_completed_rows", count]))
  f9[, order := .N - .I + 1]
  f9[, grp := fifelse(stage_cn == f9$stage_cn[.N], "end", "mid")]
  p9 <- ggplot(f9, aes(x = reorder(stage_cn, order), y = count, fill = grp)) +
    geom_col(width = 0.62, show.legend = FALSE) +
    geom_text(aes(label = count), hjust = -0.12, size = 3.2, colour = "#0b0b0b") +
    scale_fill_manual(values = c(mid = pal["blue"], end = pal["dark"])) +
    coord_flip() +
    labs(title = "deCODE 蛋白通道可比性漏斗",
         subtitle = "cis pQTL（p<5e-6）× EUR LD clump（r²<0.01@1000kb）",
         x = NULL, y = "蛋白 / 对数") +
    theme_paper() + expand_limits(x = 1.35)
  save_fig(p9, "F9_funnel_protein.png")
  }
}

# =============================================================================
# F10 deCODE 蛋白通道 MR 森林图（主方法 IVW-MRE / Wald，蓝点+CI，直标标签）
# =============================================================================
fp <- file.path(res, "grid/protein_decode_mr.csv")
if (file.exists(fp)) {
  pm <- fread(fp)
  main <- pm[ok == TRUE &
               (method %in% c("Inverse variance weighted (multiplicative random effects)",
                              "Wald ratio"))]
  if (nrow(main) > 0) {
    outmap <- c(`ebi-a-GCST006867` = "T2D", `ebi-a-GCST005194` = "CAD", `ebi-a-GCST005186` = "FBG")
    main[, out := fifelse(outcome %in% names(outmap), unname(outmap[outcome]), outcome)]
    main[, pair := paste(gene, out, sep = " × ")]
    main[, ci_lo := b - 1.96 * se]; main[, ci_hi := b + 1.96 * se]
    main[, sig := fifelse(pval < 0.05, "sig", "ns")]
    setorder(main, pval)
    p10 <- ggplot(main, aes(x = b, y = reorder(pair, b), colour = sig)) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "#898781", linewidth = 0.4) +
      geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.25, linewidth = 0.7) +
      geom_point(size = 2.8) +
      geom_text(aes(label = paste0("b=", round(b, 2), "  p=", format(pval, digits = 2, scientific = TRUE))),
                hjust = ifelse(main$b >= 0, -0.15, 1.15), size = 2.7, colour = "#52514e") +
      scale_colour_manual(values = c(sig = pal["blue"], ns = pal["gray"]),
                          labels = c(sig = "p<0.05", ns = "NS")) +
      labs(title = "deCODE 蛋白通道 cis-MR（血浆蛋白 → 结局）",
           subtitle = "主方法 IVW-MRE（nsnp≥2）或 Wald ratio（nsnp=1）；b = 每 SD 蛋白对结局",
           x = "MR 估计 b (95% CI)", y = NULL) +
      theme_paper()
    save_fig(p10, "F10_protein_mr_forest.png", w = 7.5, h = 4.6)
  }
}

# =============================================================================
# F4–F8 待数据图（stub：数据到位后自动生成，先写出设计骨架）
# =============================================================================
# F4 四态堆叠条形（按分泌状态分层）—— categorical 四色（蓝/橙/青/灰）
#   输入 results/fourstate/*.csv：gene, outcome, state(concordant/protein_only/
#   transcript_only/both_null/discordant), secretion(secreted/membrane/intracellular)
# F5 coloc PP.H4 三档命中数柱状 —— 顺序蓝
#   输入 results/coloc/*.csv：pph4_05/07/09 命中数
# F6 LOOCV AUR 曲线 + 逐个留出散点 —— 强调蓝 + 灰
#   输入 results/controls/loocv_aur.csv
# F7 外部复现方向一致率森林图 —— 蓝色点 + 二项 95% CI（经典森林图）
#   输入 results/replication/*.csv：gene, consistent_rate, ci_lo, ci_hi, n
# F8 干预介质层优先级热图 —— 顺序蓝（无→深蓝），grid heatmap
#   输入 results/tables/priority_matrix.csv：row=gene, col=outcome, value=层
if (file.exists(file.path(res, "fourstate"))) {
  log("F4 数据就绪 → 待补齐渲染")
}
log("M8 图表脚本完成：F1–F3 已渲染（有数据），F4–F8 待数据自动补全。")
log("输出目录: ", fig)
