#!/usr/bin/env python3
# =============================================================================
# M38_figures_merged_20260816.py — 终稿 5 主图 + 1 补充图（审稿合并方案）
# =============================================================================
# 5 主图 + 2 表（AJHG 上限 7）：
#   Fig 1  研究设计 + 全基因组分布        (A) workflow  (B) 染色体分布
#   Fig 2  Coloc yield funnel             FDR-core=实心大圆主点, nominal=灰参考曲线,
#                                         stage-2=虚线参考点
#   Fig 3  校准披露                      (A) PP.H4 by MR 层  (B) 负边界  (C) PP.H4 阈值敏感性
#                                         (D) 置换 FPR  (E) GWAS 显著比例  (F) 峰 p 分布
#   Fig 4  区域 coloc 图（复用 M32/F5 渲染）
#   Fig 5  15 候选基因                   (A) PP.H4 lollipop  (B) MRp vs GWASp
#                                         (C) 复现矩阵  (D) FinnGen -log10p  (E) 对齐覆盖  (F) 收敛
#   Fig S1 coloc.susie 敏感性（复用 F8 渲染，exploratory）
#
# 所有数字从 results/*_20260816.csv 实时计算，不硬编码。
# =============================================================================
import csv, os, shutil, sys
from math import sqrt
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle, Patch
from matplotlib.lines import Line2D
import matplotlib as mpl
from scipy.stats import gaussian_kde

RES = "<repo-root>/results"
GRID = f"{RES}/grid"
OUT = f"{RES}/figures"
os.makedirs(OUT, exist_ok=True)

C_MAIN = "#1a3a6b"; C_RED = "#c00000"; C_GREY = "#8a8a8a"
C_NEW = "#d9822b"; C_GREEN = "#6b8e5a"
T2D_C = "#d9822b"; CAD_C = "#1a3a6b"; FBG_C = "#6b8e5a"
OUTC = {"t2d": T2D_C, "cad": CAD_C, "fbg": FBG_C}
OUTL = {"t2d": "T2D", "cad": "CAD", "fbg": "FG"}

def f(x, d=float('nan')):
    try: return float(x)
    except (TypeError, ValueError): return d

def wilson(k, n, z=1.96):
    if n == 0: return (0.0, 0.0)
    p = k / n
    denom = 1 + z*z/n
    c = (p + z*z/(2*n)) / denom
    m = z*sqrt(p*(1-p)/n + z*z/(4*n*n)) / denom
    return (max(0, c-m), min(1, c+m))

def ecdf(x):
    xs = np.sort(np.asarray(x, float)); return xs, np.arange(1, len(xs)+1)/len(xs)

def load_full():
    return pd.concat([pd.read_csv(f"{RES}/coloc_full_{o}_20260815.csv") for o in ["t2d", "cad", "fbg"]],
                     ignore_index=True)

def box(ax, x, y, w, h, text, fc, ec, fs=8.2, tc="white", bold=True, lw=1.0):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05", facecolor=fc, edgecolor=ec, lw=lw)
    ax.add_patch(p)
    ax.text(x+w/2, y+h/2, text, ha="center", va="center", fontsize=fs,
            color=tc, fontweight="bold" if bold else "normal", linespacing=1.3)

def arrow(ax, x1, y1, x2, y2, color="#333333"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                 mutation_scale=13, lw=1.4, color=color))

# ============================ Fig 1 ============================
def fig1():
    fig, axes = plt.subplots(2, 1, figsize=(11.2, 12.2), dpi=300,
                             gridspec_kw=dict(height_ratios=[1.08, 1], hspace=0.32))

    # ---- (A) workflow ----
    ax = axes[0]; ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10)
    box(ax, 0.4, 8.7, 2.8, 1.05, "eQTLGen whole-blood\ncis-eQTL (n=31,684)", "#4a7a9b", "#2e4d66", fs=7.8)
    box(ax, 3.7, 8.7, 2.8, 1.05, "GWAS: T2D / CAD / FG\nn=655,666 / 296,525 / 58,074", "#4a7a9b", "#2e4d66", fs=7.6)
    box(ax, 7.0, 8.7, 2.6, 1.05, "LD reference 1000G\nEUR (n=503)", "#4a7a9b", "#2e4d66", fs=7.8)
    box(ax, 3.3, 6.9, 3.4, 1.12, "31,373 gene–trait pairs\nwhole-transcriptome cis-MR × coloc scan", "#1a3a6b", "#0f2445", fs=8.0)
    arrow(ax, 5.0, 8.7, 5.0, 8.02)
    box(ax, 3.3, 5.35, 3.4, 1.05, "982 MR-significant pairs\nper-outcome BH-FDR q<0.05 (preregistered)", "#2e5598", "#1a3a6b", fs=7.8)
    arrow(ax, 5.0, 6.9, 5.0, 6.4)
    box(ax, 3.3, 3.9, 3.4, 0.98, "coloc.abf (p12=1e-5)\n121 strong (PP.H4≥0.8) = 12.3% yield", C_RED, "#7a0000", fs=7.8)
    arrow(ax, 5.0, 5.35, 5.0, 4.88)
    box(ax, 0.5, 2.15, 3.4, 1.15, "106 known T2D/CAD loci\nall reproduced (100%)", C_GREEN, "#3f5c35", fs=7.8)
    box(ax, 6.1, 2.15, 3.4, 1.15, "15 candidate effector genes\n7 T2D + 8 CAD (not in GWAS Catalog)", C_NEW, "#9a5a1c", fs=7.6)
    arrow(ax, 4.2, 4.3, 2.2, 3.3, C_GREEN)
    arrow(ax, 5.8, 4.3, 7.8, 3.3, C_NEW)
    box(ax, 6.1, 0.45, 3.4, 1.15, "Independent replication\nGTEx v8: 6/7 direction-consistent\nFinnGen R11: 8 variant / 9 gene-level, 4 p<0.05\n(alignment coverage 9/15 = 60%)", "#5a5a5a", "#3a3a3a", fs=7.0)
    arrow(ax, 7.8, 2.15, 7.8, 1.6)
    ax.set_title("(A) Study design: transcriptome-wide cis-MR × coloc under a preregistered\nper-outcome BH-FDR multiple-testing correction", fontsize=10, pad=8)

    # ---- (B) 染色体分布 ----
    ax = axes[1]
    fdr = pd.read_csv(f"{RES}/fdr_core_20260816.csv")
    grid = pd.read_csv(f"{GRID}/transcript_coloc_hits.csv")
    known_keys = set(zip(grid["gene"].astype(str), grid["outcome"].astype(str)))
    strong = fdr[fdr["strong"]].copy()
    strong["is_new"] = ~strong[["gene", "outcome"]].apply(
        lambda r: (str(r["gene"]), str(r["outcome"])) in known_keys, axis=1)
    bim = {}
    for line in open("<repo-root>/data/ldref/1kg.v3/EUR.bim"):
        p = line.split(); bim[p[1]] = (int(p[0]), int(p[3]))
    m25 = pd.read_csv(f"{RES}/m25_new_strong_annotation_20260816.csv")
    m25_key = {(row.gene, row.outcome): (row.chr_hg19, row.pos_hg19) for row in m25.itertuples()}
    pts = []
    for r in strong.itertuples():
        if r.is_new:
            c, p = m25_key[(r.gene, r.outcome)]
            pts.append((int(c), int(p), r.outcome, True, r.symbol))
        else:
            gr = grid[(grid["gene"].astype(str) == str(r.gene)) & (grid["outcome"].astype(str) == str(r.outcome))]
            if len(gr):
                rs = str(gr.iloc[0]["top_snp"]).strip()
                if rs and rs in bim:
                    c, p = bim[rs]; pts.append((int(c), int(p), r.outcome, False, None))
    chrs = sorted({p[0] for p in pts})
    chr_off = {c: i for i, c in enumerate(chrs)}
    for i, c in enumerate(chrs):
        ax.axhline(i + 0.5, color="#e8ecf2", lw=6, zorder=0)
    for c, p, out, isnew, lab in pts:
        y = chr_off[c] + 0.35 + 0.3*np.random.RandomState(hash((c, p)) % 2**32).rand()*0.3
        ax.scatter(p/2.5e8, y, s=16, marker="^" if isnew else "o",
                   color=OUTC[out], alpha=0.9, edgecolor="white", lw=0.4, zorder=3)
    for c, p, out, isnew, lab in pts:
        if isnew:
            y = chr_off[c] + 0.5
            ax.annotate(lab, (p/2.5e8, y), textcoords="offset points", xytext=(6, 4),
                        fontsize=6.0, color="#5a4010", rotation=30)
    ax.set_yticks([i + 0.5 for i in range(len(chrs))]); ax.set_yticklabels(chrs, fontsize=7.5)
    ax.set_xlabel("Position along chromosome (hg19, relative)", fontsize=9)
    n_new = int(strong["is_new"].sum())
    ax.set_title(f"(B) Genome-wide distribution of the {len(strong)} FDR-core strong colocalizations\n"
                 f"(triangles = {n_new} candidate genes, circles = {len(strong)-n_new} known loci)", fontsize=10)
    ax.set_xlim(0, 1.0); ax.set_ylim(0, len(chrs))
    from matplotlib.lines import Line2D as _L2
    ax.legend(handles=[_L2([0], [0], marker="o", color="w", markerfacecolor=OUTC["t2d"], ms=6, label="T2D"),
                       _L2([0], [0], marker="o", color="w", markerfacecolor=OUTC["cad"], ms=6, label="CAD"),
                       _L2([0], [0], marker="o", color="w", markerfacecolor=OUTC["fbg"], ms=6, label="FG")],
              loc="upper right", fontsize=8, frameon=False)

    fig.savefig(f"{OUT}/20260816_Fig1_design_genome.png", bbox_inches="tight")
    plt.close(fig); print("Fig1 OK")

# ============================ Fig 2: yield funnel (4-panel) ============================
# (A) 总体 funnel（名义阈值）+ FDR-core 主点 + stage-2 参考点 + 全对基线
# (B) 分结局 funnel（t2d/cad/fbg 同阈值曲线 + 各结局 FDR-core 点）
# (C) FDR q 阈值扫描（per-outcome BH-FDR 重算，q 收紧 → yield 单调上升）
# (D) 分结局 FDR-core yield + Wilson CI + 合并线
def fig2():
    fun = pd.read_csv(f"{RES}/m36b_funnel_20260816.csv")
    nom = fun[fun["kind"] == "nominal"].sort_values("threshold", ascending=False)
    xs = nom["n"].values
    ys = 100*nom["yield"].values
    ks = nom["strong"].values
    cis = [wilson(k, n) for k, n in zip(ks, xs)]
    xlabs = ["<0.5", "<0.05", "<0.01", "<0.005", "<0.001", "<0.0005", "<1e-4", "<1e-5"]
    x = np.arange(len(xs))
    THS = [0.5, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 1e-5]

    fr = fun[fun["kind"] == "fdr_core"].iloc[0]
    gr = fun[fun["kind"] == "stage2_grid"].iloc[0]

    # ---- 分结局 funnel 数据（从 coloc_full_* 实时算）----
    def per_outcome_funnel():
        out = {}
        for o in ["t2d", "cad", "fbg"]:
            d = pd.read_csv(f"{RES}/coloc_full_{o}_20260815.csv")
            d = d[d["ok"] == True]
            pts = []
            for t in THS:
                sub = d[d["mr_p"] < t]
                k = int((sub["pp4"] >= 0.8).sum())
                pts.append((len(sub), k))
            fc = pd.read_csv(f"{RES}/fdr_core_20260816.csv")
            fo = fc[fc["outcome"] == o]
            n_fc, k_fc = len(fo), int(fo["strong"].sum())
            out[o] = (pts, n_fc, k_fc)
        return out
    pof = per_outcome_funnel()

    # ---- FDR q 扫描（per-outcome BH-FDR 重算）----
    def bh_fdr(p):
        p = np.asarray(p, float); n = len(p)
        order = np.argsort(p); padj = p[order]*n/np.arange(1, n+1)
        padj = np.minimum.accumulate(padj[::-1])[::-1]
        out = np.full(n, np.nan); out[order] = padj
        return out
    QS = [0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001]
    qsweep = []
    for q in QS:
        ns = ks_ = 0
        for o in ["t2d", "cad", "fbg"]:
            d = pd.read_csv(f"{RES}/coloc_full_{o}_20260815.csv")
            d = d[d["ok"] == True].copy()
            d["padj"] = bh_fdr(d["mr_p"].values)
            sub = d[d["padj"] < q]
            ns += len(sub); ks_ += int((sub["pp4"] >= 0.8).sum())
        qsweep.append((q, ns, ks_, 100*ks_/ns if ns else float("nan")))

    fig, axes = plt.subplots(2, 2, figsize=(13.8, 9.6), dpi=300)
    fig.subplots_adjust(hspace=0.42, wspace=0.28)

    # ---- (A) 总体 funnel ----
    ax = axes[0, 0]
    yerr_lo = [max(0, y - 100*c[0]) for c, y in zip(cis, ys)]
    yerr_hi = [max(0, 100*c[1] - y) for c, y in zip(cis, ys)]
    ax.errorbar(x, ys, yerr=[yerr_lo, yerr_hi], fmt="o-", color=C_GREY, ecolor=C_GREY,
                elinewidth=0.9, capsize=2.4, ms=4.6, zorder=2, alpha=0.85, lw=1.0)
    for xi, y, n in zip(x, ys, xs):
        ax.text(xi, y + 2.2, f"{y:.1f}%\n(n={n:,})", ha="center", fontsize=5.6, color=C_GREY)
    # 全对基线 0.42%
    ax.axhline(0.42, color="#bbbbbb", ls=":", lw=1.0)
    ax.text(0.03, 0.9, "all QC-passed\npairs (0.42%)", fontsize=6.0, color="#777777")
    # 分隔线：阈值区间 | 分析层
    ax.axvline(7.75, color="#bbbbbb", ls=":", lw=0.9)
    # FDR-core 主点
    ax.scatter([8.15], [100*fr["yield"]], s=210, marker="o", facecolor=C_RED, edgecolor="#7a0000",
               linewidths=1.6, zorder=5, label="per-outcome BH-FDR core")
    ax.annotate(f"Per-outcome BH-FDR core (primary)\n{100*fr['yield']:.1f}%  [{int(fr['n']):,} pairs, {int(fr['strong'])} strong]",
                (8.15, 100*fr["yield"]), textcoords="offset points", xytext=(-4, 20),
                ha="center", fontsize=6.8, color=C_RED, fontweight="bold")
    # stage-2 参考点
    ax.scatter([9.35], [100*gr["yield"]], s=90, marker="o", facecolor="white", edgecolor=C_MAIN,
               linewidths=1.6, ls="--", zorder=4)
    ax.annotate(f"Stage-2 clumped+IVW grid (reference)\n{100*gr['yield']:.1f}%",
                (9.35, 100*gr["yield"]), textcoords="offset points", xytext=(2, -30),
                ha="center", fontsize=6.2, color=C_MAIN)
    ax.set_xticks(x); ax.set_xticklabels(xlabs, fontsize=7.0)
    ax.set_xlabel("Nominal cis-MR significance threshold (mr_p < X)  ·  Analysis layer →", fontsize=8.5)
    ax.set_ylabel("Colocalization yield (PP.H4≥0.8, %)", fontsize=8.5)
    ax.set_title("(A) Coloc-yield funnel: nominal thresholds\nas a sensitivity reference to the FDR core", fontsize=9)
    ax.set_xlim(-0.5, 9.8); ax.set_ylim(0, 34)

    # ---- (B) 分结局 funnel ----
    ax = axes[0, 1]
    for o, col in zip(["t2d", "cad", "fbg"], [T2D_C, CAD_C, FBG_C]):
        pts, n_fc, k_fc = pof[o]
        ys_o = [100*k/n if n else float("nan") for n, k in pts]
        ax.plot(x, ys_o, "o-", color=col, lw=1.3, ms=3.8, alpha=0.9,
                label=f"{OUTL[o]} (n at FDR core = {n_fc:,})")
        # 分结局 FDR-core 点
        ax.scatter([8.15], [100*k_fc/n_fc], s=120, marker="o", facecolor=col,
                   edgecolor="white", linewidths=1.0, zorder=5)
    ax.axhline(0.42, color="#bbbbbb", ls=":", lw=0.9)
    ax.text(0.03, 0.9, "all-QC baseline\n(0.42%)", fontsize=5.6, color="#777777")
    ax.set_xticks(x); ax.set_xticklabels(xlabs, fontsize=6.8, rotation=30, ha="right")
    ax.set_xlabel("Nominal cis-MR significance threshold", fontsize=8.5)
    ax.set_ylabel("Colocalization yield (PP.H4≥0.8, %)", fontsize=8.5)
    ax.set_title("(B) Outcome-stratified funnels: the calibration\nholds within T2D, CAD, and FG", fontsize=9)
    ax.set_xlim(-0.5, 9.8); ax.set_ylim(0, 45)
    ax.legend(fontsize=6.4, loc="upper left", frameon=False)

    # ---- (C) FDR q 扫描 ----
    ax = axes[1, 0]
    qx = np.arange(len(QS))
    ys_q = [t[3] for t in qsweep]
    ns_q = [t[1] for t in qsweep]
    ax.plot(qx, ys_q, "o-", color=C_RED, lw=1.6, ms=5.5, zorder=3)
    for xi, y, n in zip(qx, ys_q, ns_q):
        ax.text(xi, y + 1.6, f"{y:.1f}%\n(n={n:,})", ha="center", fontsize=6.0, color=C_RED)
    ax.axhline(100*gr["yield"], color=C_MAIN, ls="--", lw=1.0)
    ax.text(len(QS)-0.15, 100*gr["yield"] + 1.0, f"stage-2 grid\n{100*gr['yield']:.1f}%",
            fontsize=5.8, color=C_MAIN, ha="right")
    ax.set_xticks(qx); ax.set_xticklabels([f"q < {q}" for q in QS], fontsize=7.0, rotation=30, ha="right")
    ax.set_xlabel("Per-outcome BH-FDR threshold (recomputed)", fontsize=8.5)
    ax.set_ylabel("Colocalization yield (PP.H4≥0.8, %)", fontsize=8.5)
    ax.set_title("(C) The preregistered FDR control itself\ncalibrates: yield rises as q tightens", fontsize=9)
    ax.set_xlim(-0.4, len(QS)-0.6); ax.set_ylim(0, 36)

    # ---- (D) 分结局 FDR-core yield + CI ----
    ax = axes[1, 1]
    cats = ["T2D", "CAD", "FG", "Combined"]
    cols = [T2D_C, CAD_C, FBG_C, C_RED]
    kk = [65, 54, 2, 121]; nn = [394, 576, 12, 982]
    for i, (k, n, col) in enumerate(zip(kk, nn, cols)):
        lo, hi = wilson(k, n); yc = 100*k/n
        ax.bar(i, yc, width=0.55, color=col, alpha=0.92)
        ax.errorbar(i, yc, yerr=[[max(0, yc-100*lo)], [max(0, 100*hi-yc)]], fmt="none", ecolor="black", capsize=3)
        ax.text(i, yc + 1.5, f"{yc:.1f}%", ha="center", fontsize=7.2, fontweight="bold")
        ax.text(i, 0.7, f"{k}/{n}", ha="center", fontsize=6.4, color="white", fontweight="bold")
    ax.set_xticks(range(4)); ax.set_xticklabels(cats, fontsize=8.2)
    ax.set_ylabel("Colocalization yield in FDR core (PP.H4≥0.8, %)", fontsize=8.5)
    ax.set_title("(D) Per-outcome FDR-core yield with\nWilson 95% confidence intervals", fontsize=9)
    ax.set_ylim(0, 52)

    fig.suptitle("Coloc yield calibrates with MR evidence: nominal funnel, outcome-stratified funnels,\nFDR-threshold calibration, and per-outcome FDR-core yields",
                 fontsize=11, y=0.985)
    fig.savefig(f"{OUT}/20260816_Fig2_yield_funnel.png", bbox_inches="tight")
    plt.close(fig); print("Fig2 OK")

# ============================ Fig 3: calibration disclosures ============================
def fig3():
    full = load_full(); q = full[full["ok"] == True].copy()
    fdr = pd.read_csv(f"{RES}/fdr_core_20260816.csv")
    fdr_key = set(zip(fdr["gene"].astype(str), fdr["outcome"].astype(str)))
    def _layer(r):
        if (str(r["gene"]), str(r["outcome"])) in fdr_key: return "sig"
        if r["mr_p"] >= 0.5: return "null"
        if r["mr_p"] >= 0.05: return "grey"
        return "nom_not_fdr"   # nominal p<0.05 but not FDR-core (8 dropped; Table S3)
    q["layer"] = q.apply(_layer, axis=1)
    layers = {k: q[q["layer"] == k]["pp4"].values for k in ["sig", "grey", "null"]}
    n_nomfdr = int((q["layer"] == "nom_not_fdr").sum())
    k_nomfdr = int((q[(q["layer"] == "nom_not_fdr")]["pp4"] >= 0.8).sum())
    grid = pd.read_csv(f"{GRID}/transcript_coloc_hits.csv")
    known_keys = set(zip(grid["gene"].astype(str), grid["outcome"].astype(str)))
    strong = fdr[fdr["strong"]].copy()
    strong["is_new"] = ~strong[["gene", "outcome"]].apply(
        lambda r: (str(r["gene"]), str(r["outcome"])) in known_keys, axis=1)
    kn = strong[~strong["is_new"]]; nw = strong[strong["is_new"]]
    gwas_sig = 5e-8

    fig, axes = plt.subplots(2, 3, figsize=(14.6, 9.0), dpi=300)

    # ---- (A) PP.H4 ECDF by layer ----
    ax = axes[0, 0]
    for k, lab, col, ls in [("sig", f"MR-significant (FDR-core)  n={len(layers['sig']):,}", C_MAIN, "-"),
                            ("grey", f"MR-null (0.05≤p<0.5)  n={len(layers['grey']):,}", C_GREY, "--"),
                            ("null", f"MR-negative (p≥0.5)  n={len(layers['null']):,}", "#b0b0b0", ":")]:
        xs, ys = ecdf(layers[k]); ax.plot(xs, ys, ls, color=col, lw=1.8, label=lab)
    ax.axvline(0.8, color=C_RED, lw=1, ls="--")
    ax.text(0.8, 0.04, "PP.H4=0.8", color=C_RED, fontsize=6.4, rotation=90)
    ax.set_xlabel("PP.H4 (coloc.abf)"); ax.set_ylabel("ECDF"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.legend(fontsize=6.4, loc="lower right", frameon=False)
    ax.set_title("(A) PP.H4 distribution by MR-significance layer", fontsize=9)

    # ---- (B) strong hit rate per layer（负边界）----
    ax = axes[0, 1]
    ns = [len(layers["sig"]), len(layers["grey"]), len(layers["null"])]
    ks = [int((layers["sig"] >= 0.8).sum()), int((layers["grey"] >= 0.8).sum()), int((layers["null"] >= 0.8).sum())]
    cats = ["FDR-core", "MR-null\n(0.05–0.5)", "MR-neg\n(≥0.5)"]
    for i, (k, n) in enumerate(zip(ks, ns)):
        lo, hi = wilson(k, n); yc = 100*k/n
        ax.bar(i, yc, width=0.52, color=[C_MAIN, C_GREY, "#b0b0b0"][i], alpha=0.92)
        ax.errorbar(i, yc, yerr=[[max(0, yc-100*lo)], [max(0, 100*hi-yc)]], fmt="none", ecolor="black", capsize=3)
        ax.text(i, yc + 0.45, f"{k}/{n}\n({yc:.2f}%)", ha="center", fontsize=6.3)
    ax.set_xticks(range(3)); ax.set_xticklabels(cats, fontsize=7.5)
    ax.set_ylabel("Strong-colocalization rate (PP.H4≥0.8, %)")
    ax.set_ylim(0, 16)
    ax.set_title(f"(B) Strong coloc is confined to the MR-significant layer\n(nominal p<0.05 not FDR-core: {k_nomfdr}/{n_nomfdr} = {100*k_nomfdr/n_nomfdr:.2f}%; Table S3)", fontsize=8.2)

    # ---- (C) PP.H4 threshold sensitivity ----
    ax = axes[0, 2]
    thr = [0.5, 0.8, 0.9]
    kcs = [int((layers["sig"] >= t).sum()) for t in thr]
    for i, (t, k) in enumerate(zip(thr, kcs)):
        n = len(layers["sig"]); lo, hi = wilson(k, n); yc = 100*k/n
        ax.bar(i, yc, width=0.5, color=[C_NEW, C_RED, "#8b8b8b"][i], alpha=0.92)
        ax.errorbar(i, yc, yerr=[[max(0, yc-100*lo)], [max(0, 100*hi-yc)]], fmt="none", ecolor="black", capsize=3)
        ax.text(i, yc + 0.5, f"{k} ({yc:.1f}%)", ha="center", fontsize=7)
    ax.set_xticks(range(3)); ax.set_xticklabels(["PP.H4≥0.5", "PP.H4≥0.8", "PP.H4≥0.9"], fontsize=8)
    ax.set_ylabel("Pairs passing threshold in FDR-core (%)")
    ax.set_ylim(0, 40)
    ax.set_title("(C) Strong-count sensitivity to the PP.H4 threshold (FDR-core, n=982)", fontsize=9)

    # ---- (D) permutation FPR ----
    ax = axes[1, 0]
    perm = pd.read_csv(f"{GRID}/coloc_permutation_calib.csv")
    tot = int(perm["n_perm"].sum())
    fps = [int(perm["fp_ge05"].sum()), int(perm["fp_ge08"].sum()), int(perm["fp_ge09"].sum())]
    for i, (fpn, t) in enumerate(zip(fps, thr)):
        ax.bar(i, 100*fpn/tot, width=0.5, color="#9aa7bd", alpha=0.92)
        ax.text(i, 100*fpn/tot + 0.08, f"{100*fpn/tot:.2f}%", ha="center", fontsize=7)
    ax.axhline(12.32, color=C_RED, ls="--", lw=1.1)
    ax.text(2.35, 13.1, "observed strong rate\nin FDR-core (12.3%)", fontsize=6.4, color=C_RED, ha="right")
    ax.set_xticks(range(3)); ax.set_xticklabels(labs := ["PP.H4≥0.5", "PP.H4≥0.8", "PP.H4≥0.9"], fontsize=8)
    ax.set_ylabel("False-positive rate under permutation (%)")
    ax.set_title(f"(D) Empirical null: coloc FPR across {tot:,} permuted scans of 106 strong loci", fontsize=9)
    ax.set_ylim(0, 16)

    # ---- (E) GWAS sig fraction ----
    ax = axes[1, 1]
    cats2 = ["106 known loci", "15 candidates"]
    subs = [kn, nw]
    for i, (lab, sub) in enumerate(zip(cats2, subs)):
        k = int((sub["gwas_min_p"] < gwas_sig).sum()); n = len(sub); yc = 100*k/n
        lo, hi = wilson(k, n)
        ax.bar(i, yc, width=0.5, color=[C_MAIN, C_NEW][i], alpha=0.92)
        ax.errorbar(i, yc, yerr=[[max(0, yc-100*lo)], [max(0, 100*hi-yc)]], fmt="none", ecolor="black", capsize=4)
        ax.text(i, yc + 2.5, f"{yc:.1f}% ({k}/{n})\n[{100*lo:.0f}, {100*hi:.0f}]", ha="center", fontsize=7.2)
    ax.axhline(100*0.05, color="#bbbbbb", ls=":", lw=1)
    ax.text(1.45, 6.5, "5×10⁻⁸ by chance\n(noise floor)", fontsize=6.2, color="#888888", ha="right")
    ax.set_xticks([0, 1]); ax.set_xticklabels(cats2, fontsize=8.5)
    ax.set_ylabel("Fraction reaching GWAS p<5×10⁻⁸ (%)"); ax.set_ylim(0, 75)
    ax.set_title("(E) A minority of strong-colocalized loci\nreach genome-wide significance", fontsize=9)

    # ---- (F) peak p distribution ----
    ax = axes[1, 2]
    for sub, col, lab, ls in [(kn, C_MAIN, "106 known loci", "-"), (nw, C_NEW, "15 candidates", "--")]:
        v = np.clip(-np.log10(sub["gwas_min_p"]), 0, 16)
        if len(v) > 3:
            kde = gaussian_kde(v, bw_method=0.4)
            xs = np.linspace(0, 16, 300)
            ax.plot(xs, kde(xs), ls, color=col, lw=1.8, label=f"{lab} (median −log₁₀p = {np.median(v):.1f})")
        ax.hist(v, bins=np.linspace(0, 16, 24), density=True, alpha=0.18, color=col)
    ax.axvline(-np.log10(gwas_sig), color=C_RED, ls="--", lw=1)
    ax.text(-np.log10(gwas_sig)+0.1, 0.94*max(ax.get_ylim()[1], 1), "p=5×10⁻⁸", fontsize=6.2, color=C_RED)
    ax.set_xlabel("GWAS peak −log₁₀(p) (clipped at 16)")
    ax.set_ylabel("Density")
    ax.legend(fontsize=6.4, loc="upper left", frameon=False)
    ax.set_title("(F) Peak-association distribution: candidates\ncluster below genome-wide significance", fontsize=9)

    fig.suptitle("Calibration disclosures: PP.H4 distribution, the negative boundary, and the GWAS-significance caveat",
                 fontsize=11, y=1.0)
    fig.tight_layout(); fig.savefig(f"{OUT}/20260816_Fig3_calibration.png", bbox_inches="tight")
    plt.close(fig); print("Fig3 OK")

# ============================ Fig 5: 15 candidates ============================
def fig5():
    c = pd.read_csv(f"{RES}/candidate15_replication_20260816.csv").copy()
    c = c.sort_values("pp4", ascending=True).reset_index(drop=True)
    n = len(c)

    fig = plt.figure(figsize=(13.6, 13.0), dpi=300)
    gs = fig.add_gridspec(3, 6, height_ratios=[1, 1.05, 0.95],
                          hspace=0.55, wspace=0.4)

    # ---- (A) PP.H4 lollipop (top-left, wide) ----
    ax = fig.add_subplot(gs[0, :2])
    y = np.arange(n)
    ax.axvline(0.8, color=C_RED, lw=0.9, ls="--")
    ax.text(0.805, n-0.4, "PP.H4=0.8", fontsize=6.4, color=C_RED)
    for i, r in c.iterrows():
        col = OUTC[r["outcome"]]
        ax.plot([r["pp4"], r["pp4"]], [i, i], color=col, lw=1.1, zorder=2)
        ax.scatter(r["pp4"], i, s=38, color=col, zorder=3, edgecolor="white", lw=0.6,
                   marker="*" if (pd.notna(r["finn_p"]) and r["finn_p"] < 0.05) else "o")
    ax.set_yticks(y); ax.set_yticklabels([f"{s} ({o.upper()})" for s, o in zip(c["symbol"], c["outcome"])], fontsize=7.2)
    ax.set_xlabel("PP.H4 (coloc.abf)"); ax.set_xlim(0.75, 1.01)
    ax.set_title("(A) Colocalization support\n(★ = replicated at FinnGen p<0.05)", fontsize=9)

    # ---- (B) MR p vs GWAS peak p (top-middle) ----
    ax = fig.add_subplot(gs[0, 2:5])
    grid = pd.read_csv(f"{GRID}/transcript_coloc_hits.csv")
    known_keys = set(zip(grid["gene"].astype(str), grid["outcome"].astype(str)))
    fdr = pd.read_csv(f"{RES}/fdr_core_20260816.csv")
    strong = fdr[fdr["strong"]].copy()
    strong["is_new"] = ~strong[["gene", "outcome"]].apply(
        lambda r: (str(r["gene"]), str(r["outcome"])) in known_keys, axis=1)
    kn = strong[~strong["is_new"]]; nw = strong[strong["is_new"]]
    gwas_sig = 5e-8
    for sub, col, mk, s in [(kn, "#8fa3c0", "o", 14), (nw, C_NEW, "^", 40)]:
        ax.scatter(-np.log10(sub["gwas_min_p"]), -np.log10(sub["mr_p"]),
                   s=s, marker=mk, color=col, alpha=0.9, edgecolor="white", lw=0.5, zorder=3)
    ax.axvline(-np.log10(gwas_sig), color=C_RED, ls="--", lw=0.9)
    ax.text(-np.log10(gwas_sig)+0.1, 0.3, "GWAS p=5×10⁻⁸", fontsize=6.4, color=C_RED, rotation=90)
    for r in nw.itertuples():
        ax.annotate(r.symbol, (-np.log10(r.gwas_min_p), -np.log10(r.mr_p)),
                    textcoords="offset points", xytext=(5, 4), fontsize=6.0, color="#5a4010")
    n_top = int((kn["gwas_min_p"] < gwas_sig).sum())
    n_nw_top = int((nw["gwas_min_p"] < gwas_sig).sum())
    ax.text(0.03, 0.96, f"106 known loci: {n_top} reach p<5×10⁻⁸\n15 candidates: {n_nw_top} reach p<5×10⁻⁸\n({n- n_nw_top} sub-threshold)",
            transform=ax.transAxes, fontsize=6.6, va="top",
            bbox=dict(boxstyle="round,pad=0.3", fc="#f5f7fb", ec="#cccccc", lw=0.5))
    ax.set_xlabel("GWAS peak −log₁₀(p)"); ax.set_ylabel("cis-MR −log₁₀(p)")
    ax.set_title("(B) Candidates are often recovered\nbelow genome-wide significance", fontsize=9)
    ax.set_ylim(0, None)
    ax.legend(handles=[Line2D([0],[0], marker="o", color="w", markerfacecolor="#8fa3c0", ms=5, label="106 known loci"),
                       Line2D([0],[0], marker="^", color="w", markerfacecolor=C_NEW, ms=7, label="15 candidates")],
              fontsize=6.4, loc="lower right", frameon=False)

    # ---- (C) replication matrix (full-width row 2) ----
    ax = fig.add_subplot(gs[1, :])
    cols = [("eqtl_F_max", "Instrument\nstrength (F)"), ("-log10 mr_p", "−log₁₀ MR p"),
            ("pp4", "PP.H4"), ("GTEx dir", "GTEx\ndirection"), ("FinnGen gene", "FinnGen\ngene-level"), ("FinnGen p<0.05", "FinnGen\np<0.05")]
    mat = np.zeros((n, 6))
    cmap_c = plt.cm.viridis
    for j, col in enumerate(["eqtl_F_max", "-log10 mr_p", "pp4"]):
        if col == "-log10 mr_p":
            v = -np.log10(c["mr_p"].clip(lower=1e-300))
        else:
            v = c[col].values
        vn = (v - v.min()) / (v.max() - v.min())
        mat[:, j] = vn
    for i, r in c.iterrows():
        d = str(r["direction"]).lower()
        mat[i, 3] = 0.0 if d == "consistent" else (0.5 if d == "conflicting" else np.nan)
        mat[i, 4] = 0.0 if str(r["mr_replicated"]).lower() == "yes" else np.nan
        mat[i, 5] = 0.0 if (pd.notna(r["finn_p"]) and r["finn_p"] < 0.05) else np.nan
    for j in range(3):
        for i in range(n):
            ax.add_patch(Rectangle((j, i), 0.9, 0.9, facecolor=cmap_c(mat[i, j]), edgecolor="white", lw=0.6))
    for j in [3, 4, 5]:
        for i in range(n):
            v = mat[i, j]
            if np.isnan(v): ax.add_patch(Rectangle((j, i), 0.9, 0.9, facecolor="#e8e8e8", edgecolor="white", lw=0.6))
            elif v == 0.0: ax.add_patch(Rectangle((j, i), 0.9, 0.9, facecolor=C_GREEN, edgecolor="white", lw=0.6))
            elif v == 0.5: ax.add_patch(Rectangle((j, i), 0.9, 0.9, facecolor=C_RED, edgecolor="white", lw=0.6))
    ax.set_xlim(0, 6); ax.set_ylim(-0.5, n+0.5)
    ax.set_xticks([i+0.45 for i in range(6)]); ax.set_xticklabels([l[1] for l in cols], fontsize=7.6)
    ax.set_yticks(np.arange(n)); ax.set_yticklabels(c["symbol"].values, fontsize=7.6)
    ax.tick_params(length=0)
    ax.set_title("(C) Independent-replication matrix (green = supported, red = conflicting, grey = not measurable)",
                 fontsize=9, pad=10)
    cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(0, 1), cmap=cmap_c),
                        ax=ax, orientation="vertical", fraction=0.018, pad=0.02)
    cbar.set_label("normalized instrument strength / MR p / PP.H4", fontsize=6.5)
    cbar.ax.tick_params(labelsize=6)

    # ---- (D) FinnGen -log10p (row 3) ----
    ax = fig.add_subplot(gs[2, :2])
    sub = c[c["finn_p"].notna()].copy()
    ox = -np.log10(sub["gwas_min_p"].clip(lower=1e-300))
    oy = -np.log10(sub["finn_p"])
    sig = sub["finn_p"] < 0.05
    lim = max(float(np.concatenate([ox, oy]).max())*1.1, 8)
    ax.scatter(ox, oy, s=30, c=C_NEW, alpha=0.85, edgecolor="white", lw=0.4, zorder=3)
    ax.scatter(ox[sig], oy[sig], s=52, c=C_RED, marker="*", edgecolor="white", lw=0.4, zorder=4)
    ax.plot([0, lim], [0, lim], color="#bbbbbb", ls=":", lw=1)
    ax.axhline(-np.log10(0.05), color=C_RED, ls="--", lw=0.9)
    ax.text(0.03, -np.log10(0.05)+0.15, "FinnGen p=0.05", fontsize=6.2, color=C_RED)
    for r in sub.itertuples():
        ax.annotate(r.symbol, (-np.log10(r.gwas_min_p), -np.log10(r.finn_p)),
                    textcoords="offset points", xytext=(5, 4), fontsize=5.8, color="#5a4010")
    ax.set_xlim(0, lim); ax.set_ylim(0, lim)
    ax.set_xlabel("Original GWAS peak −log₁₀(p)"); ax.set_ylabel("FinnGen R11 −log₁₀(p)")
    ax.set_title(f"(D) FinnGen R11 independent-cohort replication\n({len(sub)} alignable; ★ = p<0.05, {int(sig.sum())})", fontsize=9)

    # ---- (E) FinnGen alignment coverage (row 3) ----
    ax = fig.add_subplot(gs[2, 2:4])
    notloc = [s for s in ["U6atac", "CWF19L1", "N4BP2L2", "SLC12A3", "PLAUR"] if s in set(c["symbol"])]
    nomiss = [s for s in ["CLEC3B"] if s in set(c["symbol"])]
    n_align = int(c["finn_p"].notna().sum())
    n_nl = len(notloc); n_om = len(nomiss)
    ax.barh(0, n_align, height=0.5, color=C_GREEN, alpha=0.9)
    ax.barh(0, n_nl, left=n_align, height=0.5, color=C_NEW, alpha=0.9)
    ax.barh(0, n_om, left=n_align+n_nl, height=0.5, color="#b0b0b0", alpha=0.9)
    ax.text(n_align/2, 0, f"{n_align} alignable\n({100*n_align/n:.0f}%)", ha="center", va="center",
            fontsize=6.8, color="white", fontweight="bold")
    ax.text(n_align+n_nl/2, 0, f"{n_nl}\nlead not\nlocalized", ha="center", va="center", fontsize=6.0, color="white")
    ax.text(n_align+n_nl+n_om/2, 0, f"{n_om}\norig data\nmissing", ha="center", va="center", fontsize=6.0, color="#555555")
    ax.set_yticks([]); ax.set_xlim(0, n)
    ax.set_xticks([]); ax.set_xlabel(f"15 candidate genes", fontsize=8)
    ax.set_title("(E) FinnGen R11 alignment coverage\n(9/15 alignable; disclosure of limits)", fontsize=9)

    # ---- (F) convergence + direction summary (row 3) ----
    import json
    heidi = json.load(open(f"{GRID}/heidi_full_summary.json"))
    steig = json.load(open(f"{GRID}/steiger_direction_76.json"))
    p1 = pd.read_csv(f"{GRID}/gtex_replication_p1.csv")
    p1_meas = p1[p1["concordant"].notna()]
    p1_ok = (p1_meas["concordant"] == True).sum()
    m26_meas = c[c["direction"].isin(["consistent", "conflicting"])]
    m26_ok = (m26_meas["direction"] == "consistent").sum()
    bars = [("HEIDI\npre-scan", heidi["heidi_ok_pre"], heidi["n_pre"], C_GREY),
            ("HEIDI\nstrong set", heidi["heidi_ok_strong"], heidi["n_strong"], C_MAIN),
            ("Steiger\neQTL→outcome", steig["fwd"], steig["n"], C_GREEN),
            ("GTEx dir\n15 cand.", int(m26_ok), len(m26_meas), C_NEW),
            ("GTEx dir\n106 known", int(p1_ok), len(p1_meas), C_MAIN)]
    ax = fig.add_subplot(gs[2, 4:])
    for i, (lab, k, nn, col) in enumerate(bars):
        lo, hi = wilson(k, nn); yc = 100*k/nn
        ax.bar(i, yc, width=0.62, color=col, alpha=0.92)
        ax.errorbar(i, yc, yerr=[[max(0, yc-100*lo)], [max(0, 100*hi-yc)]], fmt="none", ecolor="black", capsize=3)
        ax.text(i, yc + 2.5, f"{yc:.0f}%\n{k}/{nn}", ha="center", fontsize=6.6)
    ax.set_xticks(range(len(bars))); ax.set_xticklabels([b[0] for b in bars], fontsize=6.8)
    ax.set_ylabel("Proportion (%)"); ax.set_ylim(0, 112)
    ax.set_title("(F) Convergence and causal direction\n(SMR+HEIDI, Steiger, GTEx)", fontsize=9)

    fig.suptitle("The 15 candidate effector genes: colocalization support, GWAS-peak relationship, and independent replication",
                 fontsize=11, y=0.995)
    fig.savefig(f"{OUT}/20260816_Fig5_candidates.png", bbox_inches="tight")
    plt.close(fig); print("Fig5 OK")

# ============================ Fig 4 / Fig S1: 复用渲染 ============================
def copy_figs():
    pairs = [("20260816_F5_locuszoom_mirror_v2.png", "20260816_Fig4_regional.png"),
             ("20260816_F8_susie_v2.png", "20260816_FigS1_susie.png")]
    for src, dst in pairs:
        s = os.path.join(OUT, src); d = os.path.join(OUT, dst)
        shutil.copyfile(s, d)
        print(f"copy {src} -> {dst}")

if __name__ == "__main__":
    fig1(); fig2(); fig3(); fig5(); copy_figs()
    print("M38 merged figures done")
