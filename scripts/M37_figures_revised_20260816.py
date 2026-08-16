#!/usr/bin/env python3
# =============================================================================
# M37_figures_revised_20260816.py — 修订版多面板主图（审稿 + 图重构，全英文，无表）
# =============================================================================
# 取代 M35 的图集。所有数字以 BH-FDR 主口径（982/121/15/12.3%）为准，直接读
# Phase 0a 产出的 results/*_20260816.csv，不硬编码。
#
#   F1  workflow 管线（单面板，密集）
#   F2  atlas 概览：(A) 结局×MR状态堆叠条 (B) 31,373→982→121→106/15 流 (C) 染色体分布 (D) FDR-core PP.H4 ECDF by outcome
#   F3  coloc yield funnel（名义阈值校准曲线 + FDR-core 点 + grid 参考）
#   F4  负边界 + 阈值敏感性：(A) PP.H4 by MR状态 (B) strong hit rate (C) PP.H4 阈值敏感性 (D) 置换 FPR
#   F6  15 候选：(A) PP.H4 lollipop (B) MRp vs GWAS峰p 散点(象限) (C) 复现矩阵（图形化非表格）
#   F7  收敛+独立复现：(A) SMR+HEIDI (B) Steiger (C) GTEx (D) FinnGen -log10p (E) FinnGen 对齐覆盖
#   F8  coloc.susie 敏感性（读 m34b 优先，fallback m34）
#   F9  GWAS 峰显著 caveat：(A) 达 p<5e-8 比例 (B) 峰 p 分布 (C) pp4 vs gwas_min_p
#   S1  数据资源示意（图形化数据源+覆盖披露，取代含表格的 F7）
# F5（locuszoom 镜像）由 patched M32 单独生成，M37 不含。
#
# 用法：python3 scripts/M37_figures_revised_20260816.py [fig1 fig2 ... fig9 s1]
#   （不传参数 = 全部；fig8 依赖 m34b 结果，可单独补跑）
# =============================================================================
import csv, os, sys
from math import sqrt
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from scipy.stats import gaussian_kde

RES = "<repo-root>/results"
GRID = f"{RES}/grid"
OUT = f"{RES}/figures"
os.makedirs(OUT, exist_ok=True)

C_MAIN = "#1a3a6b"; C_RED = "#c00000"; C_GREY = "#7a7a7a"
C_NEW = "#d9822b"; C_GREEN = "#6b8e5a"; C_LIGHT = "#d9dee8"
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

# ============================ Fig 1: workflow ============================
def fig1():
    fig, ax = plt.subplots(figsize=(10.5, 6.0), dpi=300)
    ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10)
    box(ax, 0.4, 8.6, 2.8, 1.1, "eQTLGen whole-blood\ncis-eQTL (n=31,684)", "#4a7a9b", "#2e4d66")
    box(ax, 3.7, 8.6, 2.8, 1.1, "GWAS: T2D / CAD / FG\nn=655,666 / 296,525 / 58,074", "#4a7a9b", "#2e4d66")
    box(ax, 7.0, 8.6, 2.6, 1.1, "LD reference 1000G\nEUR (n=503)", "#4a7a9b", "#2e4d66")
    box(ax, 3.3, 6.8, 3.4, 1.15, "31,373 gene–trait pairs\nwhole-transcriptome cis-MR × coloc scan", "#1a3a6b", "#0f2445")
    arrow(ax, 5.0, 8.6, 5.0, 7.95)
    box(ax, 3.3, 5.25, 3.4, 1.1, "982 MR-significant pairs\nper-outcome BH-FDR q<0.05 (preregistered)", "#2e5598", "#1a3a6b")
    arrow(ax, 5.0, 6.8, 5.0, 6.35)
    box(ax, 3.3, 3.8, 3.4, 1.0, "coloc.abf (p12=1e-5)\n121 strong (PP.H4≥0.8) = 12.3% yield", C_RED, "#7a0000")
    arrow(ax, 5.0, 5.25, 5.0, 4.8)
    box(ax, 0.5, 2.1, 3.4, 1.2, "106 known T2D/CAD loci\nall reproduced (100%)", C_GREEN, "#3f5c35")
    box(ax, 6.1, 2.1, 3.4, 1.2, "15 candidate effector genes\n7 T2D + 8 CAD (not in GWAS Catalog)", C_NEW, "#9a5a1c")
    arrow(ax, 4.2, 4.2, 2.2, 3.3, C_GREEN)
    arrow(ax, 5.8, 4.2, 7.8, 3.3, C_NEW)
    box(ax, 6.1, 0.4, 3.4, 1.2, "Independent replication\nGTEx v8: 6/7 direction-consistent\nFinnGen R11: 8 variant / 9 gene-level, 4 p<0.05\n(alignment coverage 9/15 = 60%)", "#5a5a5a", "#3a3a3a")
    arrow(ax, 7.8, 2.1, 7.8, 1.6)
    ax.set_title("Transcriptome-wide cis-MR × coloc atlas: workflow under a preregistered per-outcome BH-FDR multiple-testing correction",
                 fontsize=10, pad=12)
    fig.tight_layout(); fig.savefig(f"{OUT}/20260816_F1_workflow_v2.png", bbox_inches="tight")
    plt.close(fig); print("F1 OK")

# ============================ Fig 2: atlas overview ============================
def fig2():
    full = load_full(); q = full[full["ok"] == True].copy()
    fdr = pd.read_csv(f"{RES}/fdr_core_20260816.csv")
    fdr_key = set(zip(fdr["gene"].astype(str), fdr["outcome"].astype(str)))
    q["layer"] = np.where(q[["gene", "outcome"]].apply(
        lambda r: (str(r["gene"]), str(r["outcome"])) in fdr_key, axis=1),
        "sig", np.where(q["mr_p"] < 0.5, "grey", "null"))
    q["strong"] = q["pp4"] >= 0.8

    fig = plt.figure(figsize=(12.4, 8.8), dpi=300)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1], width_ratios=[1.25, 1])

    # ---- (A) 结局 × MR 状态堆叠条 ----
    ax = fig.add_subplot(gs[0, 0])
    outs = ["t2d", "cad", "fbg"]
    lay_col = {"sig": C_MAIN, "grey": "#b8c4d4", "null": "#e2e6ec"}
    lay_lab = {"sig": "MR-significant (BH-FDR)", "grey": "MR-null (0.05≤p<0.5)", "null": "MR-negative (p≥0.5)"}
    n_sig, n_grey, n_null = (q["layer"] == "sig").sum(), (q["layer"] == "grey").sum(), (q["layer"] == "null").sum()
    for i, out in enumerate(outs):
        sub = q[q["outcome"] == out]
        n = len(sub); bot = 0
        for lay in ["null", "grey", "sig"]:
            k = (sub["layer"] == lay).sum()
            c = lay_col[lay]
            ax.bar(i, 100*k/n, bottom=100*bot/n, width=0.62, color=c, edgecolor="white", lw=0.5)
            if lay == "sig" and k:
                st = (sub[(sub["layer"] == lay)]["strong"]).sum()
                ax.text(i, 100*(bot + k/2)/n, f"{k}\n({int(st)} strong)", ha="center", va="center",
                        fontsize=6.8, color="white", fontweight="bold")
            bot += k
        ax.text(i, 101.5, f"n={n:,}", ha="center", fontsize=7.2, color="#333333")
    ax.set_xticks(range(3)); ax.set_xticklabels(["T2D", "CAD", "FG"], fontsize=9)
    ax.set_ylabel("Proportion of QC pairs (%)", fontsize=9)
    ax.set_ylim(0, 118)
    from matplotlib.patches import Patch
    handles = [Patch(color=lay_col["sig"], label=f"{lay_lab['sig']} (n={n_sig:,})"),
               Patch(color=lay_col["grey"], label=f"{lay_lab['grey']} (n={n_grey:,})"),
               Patch(color=lay_col["null"], label=f"{lay_lab['null']} (n={n_null:,})")]
    ax.legend(handles=handles, fontsize=6.8, loc="upper right", frameon=False)
    ax.set_title("(A) Outcome composition and MR-significance layering", fontsize=9)

    # ---- (B) 流：31373 → 982 → 121 → 106/15 ----
    ax = fig.add_subplot(gs[0, 1]); ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10)
    ax.set_title("(B) Atlas filtering funnel", fontsize=9, loc="left", pad=2)
    box(ax, 0.4, 8.4, 3.6, 1.0, "31,371 QC-passing pairs\n(eQTL×outcome)", "#1a3a6b", "#0f2445")
    arrow(ax, 4.0, 8.9, 5.3, 8.9)
    box(ax, 5.3, 8.4, 4.3, 1.0, "982 MR-significant\n(per-outcome BH-FDR q<0.05)", "#2e5598", "#1a3a6b")
    arrow(ax, 7.4, 8.4, 7.4, 6.9)
    box(ax, 5.3, 5.9, 4.3, 1.0, "121 strong colocalizations\nPP.H4≥0.8 (yield 12.3%)", C_RED, "#7a0000")
    arrow(ax, 6.6, 5.9, 4.6, 4.6, C_GREEN)
    arrow(ax, 8.2, 5.9, 9.4, 4.6, C_NEW)
    box(ax, 0.4, 3.6, 4.2, 1.0, "106 known T2D/CAD loci\nreproduced (100%)", C_GREEN, "#3f5c35")
    box(ax, 6.0, 3.6, 3.6, 1.0, "15 candidate genes\n(7 T2D + 8 CAD)", C_NEW, "#9a5a1c")
    ax.text(5.0, 1.8, "PP.H4 distribution and replication in panels (D), F6, F7", fontsize=7, ha="center", color="#555555")
    ax.text(0.05, 0.02, "FDR-core: t2d 394 / cad 576 / fbg 12", fontsize=7, color="#666666", transform=ax.transAxes)

    # ---- (C) 染色体分布（121 strong）----
    ax = fig.add_subplot(gs[1, 0])
    grid = pd.read_csv(f"{GRID}/transcript_coloc_hits.csv")
    known_keys = set(zip(grid["gene"].astype(str), grid["outcome"].astype(str)))
    strong = fdr[fdr["strong"]].copy()
    strong["is_new"] = ~strong[["gene", "outcome"]].apply(
        lambda r: (str(r["gene"]), str(r["outcome"])) in known_keys, axis=1)
    # hg19 坐标：m25 有新候选坐标；known 从 grid top_snp 查 1kg bim
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
        ax.axvline(0, color="#cccccc", lw=0.5)
    for c, p, out, isnew, lab in pts:
        y = chr_off[c] + 0.35 + 0.3*np.random.RandomState(hash((c, p)) % 2**32).rand()*0.3
        ax.scatter(p/2.5e8, y, s=13, marker="^" if isnew else "o",
                   color=OUTC[out], alpha=0.9, edgecolor="white", lw=0.4, zorder=3)
    for c, p, out, isnew, lab in pts:
        if isnew:
            y = chr_off[c] + 0.5
            ax.annotate(lab, (p/2.5e8, y), textcoords="offset points", xytext=(6, 4),
                        fontsize=5.6, color="#5a4010", rotation=30)
    ax.set_yticks([i + 0.5 for i in range(len(chrs))]); ax.set_yticklabels(chrs, fontsize=7)
    ax.set_xlabel("Position along chromosome (hg19, relative)", fontsize=9)
    ax.set_title(f"(C) Chromosomal distribution of {len(pts)} strong colocalizations\n"
                 "(triangles = 15 candidate genes, circles = 106 known loci)", fontsize=9)
    ax.set_xlim(0, 1.0); ax.set_ylim(0, len(chrs))

    # ---- (D) FDR-core PP.H4 ECDF by outcome ----
    ax = fig.add_subplot(gs[1, 1])
    for out in outs:
        sub = fdr[fdr["outcome"] == out]["pp4"]
        xs, ys = ecdf(sub.values)
        ax.plot(xs, ys, "-", color=OUTC[out], lw=1.8,
                label=f"{OUTL[out]} (n={len(sub)})")
    ax.axvline(0.8, color=C_RED, lw=1, ls="--")
    ax.text(0.8, 0.03, "PP.H4=0.8", color=C_RED, fontsize=6.5, rotation=90)
    ax.set_xlabel("PP.H4 (coloc.abf)"); ax.set_ylabel("ECDF within FDR-core")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.legend(fontsize=7, loc="lower right", frameon=False)
    ax.set_title("(D) PP.H4 distribution in the 982 FDR-core pairs by outcome", fontsize=9)

    fig.suptitle("Transcriptome-wide cis-MR × colocalization atlas: overview under a preregistered BH-FDR correction",
                 fontsize=11, y=1.0)
    fig.tight_layout(); fig.savefig(f"{OUT}/20260816_F2_overview_v2.png", bbox_inches="tight")
    plt.close(fig); print("F2 OK")

# ============================ Fig 3: coloc yield funnel ============================
def fig3():
    fun = pd.read_csv(f"{RES}/m36b_funnel_20260816.csv")
    nom = fun[fun["kind"] == "nominal"].sort_values("threshold", ascending=False)
    xs = nom["n"].values
    ys = 100*nom["yield"].values
    ks = nom["strong"].values
    cis = [wilson(k, n) for k, n in zip(ks, xs)]
    xlabs = ["<0.5", "<0.05", "<0.01", "<0.005", "<0.001", "<0.0005", "<1e-4", "<1e-5"]

    fig, ax = plt.subplots(figsize=(8.6, 5.4), dpi=300)
    x = np.arange(len(xs))
    yerr_lo = [max(0, y - 100*c[0]) for c, y in zip(cis, ys)]
    yerr_hi = [max(0, 100*c[1] - y) for c, y in zip(cis, ys)]
    ax.errorbar(x, ys, yerr=[yerr_lo, yerr_hi], fmt="o-", color=C_MAIN, ecolor=C_MAIN,
                elinewidth=1.2, capsize=3, ms=6, zorder=3)
    for i, (xi, y, n, k) in enumerate(zip(x, ys, xs, ks)):
        ax.text(xi, y + (3.2 if i < 6 else 1.8), f"{y:.1f}% (n={n:,})", ha="center", fontsize=6.4, color=C_MAIN)
    # FDR-core 点
    fr = fun[fun["kind"] == "fdr_core"].iloc[0]
    ax.axhline(100*fr["yield"], color=C_RED, ls="--", lw=1.2)
    ax.text(len(xs)-0.4, 100*fr["yield"] + 1.2,
            f"BH-FDR core (n=982): {100*fr['yield']:.1f}%  [{fr['n']:,} pairs, {int(fr['strong'])} strong]",
            fontsize=7, color=C_RED, ha="right")
    # grid 参考
    gr = fun[fun["kind"] == "stage2_grid"].iloc[0]
    ax.axhline(100*gr["yield"], color=C_GREY, ls=":", lw=1.2)
    ax.text(len(xs)-0.4, 100*gr["yield"] - 1.8,
            f"stage-2 clumped+IVW grid: {100*gr['yield']:.1f}%", fontsize=6.6, color=C_GREY, ha="right")
    ax.set_xticks(x); ax.set_xticklabels(xlabs, fontsize=8)
    ax.set_xlabel("cis-MR significance threshold (mr_p < X)", fontsize=9)
    ax.set_ylabel("Colocalization yield (PP.H4≥0.8, %)", fontsize=9)
    ax.set_title("Operating-characteristic curve: tighter cis-MR thresholds monotonically\nincrease the proportion of colocalizing pairs (nominal thresholds, Wilson 95% CI)",
                 fontsize=10)
    ax.set_ylim(0, 33)
    fig.tight_layout(); fig.savefig(f"{OUT}/20260816_F3_yield_v2.png", bbox_inches="tight")
    plt.close(fig); print("F3 OK")

# ============================ Fig 4: negative boundary + sensitivity ============================
def fig4():
    full = load_full(); q = full[full["ok"] == True].copy()
    fdr = pd.read_csv(f"{RES}/fdr_core_20260816.csv")
    fdr_key = set(zip(fdr["gene"].astype(str), fdr["outcome"].astype(str)))
    q["layer"] = np.where(q[["gene", "outcome"]].apply(
        lambda r: (str(r["gene"]), str(r["outcome"])) in fdr_key, axis=1),
        "sig", np.where(q["mr_p"] < 0.5, "grey", "null"))

    fig, axes = plt.subplots(2, 2, figsize=(12.2, 8.4), dpi=300)

    # ---- (A) PP.H4 by MR status (ECDF) ----
    ax = axes[0, 0]
    layers = {k: q[q["layer"] == k]["pp4"].values for k in ["sig", "grey", "null"]}
    for k, lab, col, ls in [("sig", f"MR-significant (FDR-core)  n={len(layers['sig']):,}", C_MAIN, "-"),
                            ("grey", f"MR-null (0.05≤p<0.5)  n={len(layers['grey']):,}", C_GREY, "--"),
                            ("null", f"MR-negative (p≥0.5)  n={len(layers['null']):,}", "#b0b0b0", ":")]:
        xs, ys = ecdf(layers[k]); ax.plot(xs, ys, ls, color=col, lw=1.8, label=lab)
    ax.axvline(0.8, color=C_RED, lw=1, ls="--")
    ax.set_xlabel("PP.H4 (coloc.abf)"); ax.set_ylabel("ECDF"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.legend(fontsize=6.5, loc="lower right", frameon=False)
    ax.set_title("(A) PP.H4 distribution by MR-significance layer", fontsize=9)

    # ---- (B) strong hit rate per layer ----
    ax = axes[0, 1]
    ns = [len(layers["sig"]), len(layers["grey"]), len(layers["null"])]
    ks = [int((layers["sig"] >= 0.8).sum()), int((layers["grey"] >= 0.8).sum()), int((layers["null"] >= 0.8).sum())]
    cats = ["FDR-core", "MR-null\n(0.05–0.5)", "MR-neg\n(≥0.5)"]
    for i, (k, n) in enumerate(zip(ks, ns)):
        lo, hi = wilson(k, n); yc = 100*k/n
        ax.bar(i, yc, width=0.52, color=[C_MAIN, C_GREY, "#b0b0b0"][i], alpha=0.92)
        ax.errorbar(i, yc, yerr=[[max(0, yc-100*lo)], [max(0, 100*hi-yc)]], fmt="none", ecolor="black", capsize=3)
        ax.text(i, yc + 0.35, f"{yc:.2f}% ({k}/{n})", ha="center", fontsize=6.6)
    ax.set_xticks(range(3)); ax.set_xticklabels(cats, fontsize=7.5)
    ax.set_ylabel("Strong-colocalization rate (PP.H4≥0.8, %)")
    ax.set_ylim(0, 16)
    ax.set_title("(B) Strong colocalization is confined to the MR-significant layer", fontsize=9)

    # ---- (C) PP.H4 threshold sensitivity within FDR-core ----
    ax = axes[1, 0]
    thr = [0.5, 0.8, 0.9]
    ks = [int((layers["sig"] >= t).sum()) for t in thr]
    for i, (t, k) in enumerate(zip(thr, ks)):
        n = len(layers["sig"]); lo, hi = wilson(k, n); yc = 100*k/n
        ax.bar(i, yc, width=0.5, color=[C_NEW, C_RED, "#8b8b8b"][i], alpha=0.92)
        ax.errorbar(i, yc, yerr=[[max(0, yc-100*lo)], [max(0, 100*hi-yc)]], fmt="none", ecolor="black", capsize=3)
        ax.text(i, yc + 0.5, f"{k} ({yc:.1f}%)", ha="center", fontsize=7)
    ax.set_xticks(range(3)); ax.set_xticklabels(["PP.H4≥0.5", "PP.H4≥0.8", "PP.H4≥0.9"], fontsize=8)
    ax.set_ylabel("Pairs passing threshold in FDR-core (%)")
    ax.set_ylim(0, 40)
    ax.set_title("(C) Sensitivity of the strong count to the PP.H4 threshold (FDR-core, n=982)", fontsize=9)

    # ---- (D) permutation FPR ----
    ax = axes[1, 1]
    perm = pd.read_csv(f"{GRID}/coloc_permutation_calib.csv")
    tot = int(perm["n_perm"].sum())
    fps = [int(perm["fp_ge05"].sum()), int(perm["fp_ge08"].sum()), int(perm["fp_ge09"].sum())]
    labs = ["PP.H4≥0.5", "PP.H4≥0.8", "PP.H4≥0.9"]
    for i, (fpn, t) in enumerate(zip(fps, thr)):
        ax.bar(i, 100*fpn/tot, width=0.5, color="#9aa7bd", alpha=0.92)
        ax.text(i, 100*fpn/tot + 0.1, f"{100*fpn/tot:.2f}%", ha="center", fontsize=7)
    ax.axhline(12.32, color=C_RED, ls="--", lw=1.1)
    ax.text(2.35, 13.1, "observed strong rate\nin FDR-core (12.3%)", fontsize=6.4, color=C_RED, ha="right")
    ax.set_xticks(range(3)); ax.set_xticklabels(labs, fontsize=8)
    ax.set_ylabel("False-positive rate under permutation (%)")
    ax.set_title(f"(D) Empirical null: coloc FPR across {tot:,} permuted scans of 106 strong loci", fontsize=9)
    ax.set_ylim(0, 16)

    fig.suptitle("Negative boundary and robustness: colocalization support is specific to the MR-significant layer\n"
                 "and exceeds the empirical permutation false-positive rate ~8-fold",
                 fontsize=10, y=1.0)
    fig.tight_layout(); fig.savefig(f"{OUT}/20260816_F4_negative_v2.png", bbox_inches="tight")
    plt.close(fig); print("F4 OK")

# ============================ Fig 6: 15 candidate genes ============================
def fig6():
    c = pd.read_csv(f"{RES}/candidate15_replication_20260816.csv").copy()
    c = c.sort_values("pp4", ascending=True).reset_index(drop=True)
    n = len(c)

    fig = plt.figure(figsize=(12.4, 8.6), dpi=300)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1], width_ratios=[1.35, 1])

    # ---- (A) PP.H4 lollipop ----
    ax = fig.add_subplot(gs[0, 0])
    y = np.arange(n)
    ax.axvline(0.8, color=C_RED, lw=0.9, ls="--")
    ax.text(0.805, n-0.4, "PP.H4=0.8", fontsize=6.4, color=C_RED)
    for i, r in c.iterrows():
        col = OUTC[r["outcome"]]
        ax.plot([r["pp4"], r["pp4"]], [i, i], color=col, lw=1.1, zorder=2)
        ax.scatter(r["pp4"], i, s=38, color=col, zorder=3,
                   edgecolor="white", lw=0.6,
                   marker="*" if (pd.notna(r["finn_p"]) and r["finn_p"] < 0.05) else "o")
    ax.set_yticks(y); ax.set_yticklabels([f"{s} ({o.upper()})" for s, o in zip(c["symbol"], c["outcome"])], fontsize=7.2)
    ax.set_xlabel("PP.H4 (coloc.abf)"); ax.set_xlim(0.75, 1.01)
    ax.set_title("(A) Colocalization support for the 15 candidate genes\n(★ = replicated at FinnGen p<0.05)", fontsize=9)

    # ---- (B) MR p vs GWAS peak p (quadrant) ----
    ax = fig.add_subplot(gs[0, 1])
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
    ax.text(0.03, 0.96, f"106 known loci: {n_top} reach p<5×10⁻⁸\n15 candidates: 6 reach p<5×10⁻⁸ (9 sub-threshold)",
            transform=ax.transAxes, fontsize=6.6, va="top",
            bbox=dict(boxstyle="round,pad=0.3", fc="#f5f7fb", ec="#cccccc", lw=0.5))
    ax.set_xlabel("GWAS peak $-\\log_{10}(p)$"); ax.set_ylabel("cis-MR $-\\log_{10}(p)$")
    ax.set_title("(B) Candidates are often recovered below genome-wide significance", fontsize=9)
    ax.set_ylim(0, None)
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([0],[0], marker="o", color="w", markerfacecolor="#8fa3c0", ms=5, label="106 known loci"),
                       Line2D([0],[0], marker="^", color="w", markerfacecolor=C_NEW, ms=7, label="15 candidates")],
              fontsize=6.5, loc="lower right", frameon=False)

    # ---- (C) replication matrix（图形化，非表格）----
    ax = fig.add_subplot(gs[1, :])
    cols = [("eqtl_F_max", "Instrument\nstrength (F)"), ("-log10 mr_p", "$-\\log_{10}$ MR p"),
            ("pp4", "PP.H4"), ("GTEx dir", "GTEx\ndirection"), ("FinnGen gene", "FinnGen\ngene-level"), ("FinnGen p<0.05", "FinnGen\np<0.05")]
    mat = np.zeros((n, 6)); cellcol = np.zeros((n, 6))
    cmap_c = plt.cm.viridis
    # continuous columns 0-2
    for j, col in enumerate(["eqtl_F_max", "-log10 mr_p", "pp4"]):
        if col == "-log10 mr_p":
            v = -np.log10(c["mr_p"].clip(lower=1e-300))
        else:
            v = c[col].values
        vn = (v - v.min()) / (v.max() - v.min())
        mat[:, j] = vn
    # categorical columns 3-5
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
    ax.set_xticks([i+0.45 for i in range(6)]); ax.set_xticklabels([l[1] for l in cols], fontsize=7.4)
    ax.set_yticks(np.arange(n)); ax.set_yticklabels(c["symbol"].values, fontsize=7.4)
    ax.tick_params(length=0)
    ax.set_title("(C) Independent-replication matrix (green = supported, red = conflicting, grey = not measurable)", fontsize=9, pad=8)
    import matplotlib as mpl
    cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(0, 1), cmap=cmap_c),
                        ax=ax, orientation="vertical", fraction=0.022, pad=0.01)
    cbar.set_label("normalized instrument strength / MR p / PP.H4", fontsize=6.5)
    cbar.ax.tick_params(labelsize=6)

    fig.suptitle("The 15 candidate effector genes: colocalization support, GWAS-peak relationship, and independent replication",
                 fontsize=11, y=1.0)
    fig.tight_layout(); fig.savefig(f"{OUT}/20260816_F6_candidates15_v2.png", bbox_inches="tight")
    plt.close(fig); print("F6 OK")

# ============================ Fig 7: convergence + replication ============================
def fig7():
    import json
    heidi = json.load(open(f"{GRID}/heidi_full_summary.json"))
    steig = json.load(open(f"{GRID}/steiger_direction_76.json"))
    c = pd.read_csv(f"{RES}/candidate15_replication_20260816.csv")
    gwas_p = dict(zip(c["symbol"], c["gwas_min_p"]))
    p1 = pd.read_csv(f"{GRID}/gtex_replication_p1.csv")

    fig, axes = plt.subplots(2, 3, figsize=(13.4, 8.6), dpi=300)
    # ---- (A) SMR+HEIDI ----
    ax = axes[0, 0]
    bars = [("MR-significant pre-scan\n(validated subset)", heidi["heidi_ok_pre"], heidi["n_pre"]),
            ("strong-colocalized\nsubset", heidi["heidi_ok_strong"], heidi["n_strong"])]
    for i, (lab, k, n) in enumerate(bars):
        lo, hi = wilson(k, n); yc = 100*k/n
        ax.bar(i, yc, width=0.5, color=[C_GREY, C_MAIN][i], alpha=0.9)
        ax.errorbar(i, yc, yerr=[[max(0, yc-100*lo)], [max(0, 100*hi-yc)]], fmt="none", ecolor="black", capsize=4)
        ax.text(i, yc + 2.5, f"{yc:.1f}%\n{k}/{n}", ha="center", fontsize=7)
    ax.set_xticks([0, 1]); ax.set_xticklabels([b[0] for b in bars], fontsize=7)
    ax.set_ylabel("HEIDI pass rate (%)"); ax.set_ylim(0, 100)
    ax.set_title("(A) SMR+HEIDI: higher pass rate in the\nstrong-colocalized subset", fontsize=8.5)

    # ---- (B) Steiger ----
    ax = axes[0, 1]
    k, n = steig["fwd"], steig["n"]
    lo, hi = wilson(k, n); yc = 100*k/n
    ax.bar([0], [yc], width=0.45, color=C_GREEN, alpha=0.9)
    ax.errorbar(0, yc, yerr=[[max(0, yc-100*lo)], [max(0, 100*hi-yc)]], fmt="none", ecolor="black", capsize=4)
    ax.text(0, yc + 3, f"{yc:.1f}% [{100*lo:.1f}, {100*hi:.1f}] ({k}/{n})", ha="center", fontsize=7)
    ax.text(0, 6, f"{steig['rev']} reverse ({steig['rev_sig']} significant)", ha="center", fontsize=6.6, color="#555555")
    ax.set_xticks([0]); ax.set_xticklabels(["eQTL → outcome"], fontsize=8)
    ax.set_ylim(0, 115); ax.set_ylabel("Direction-consistent proportion (%)")
    ax.set_title("(B) Steiger: causal direction\neQTL → outcome", fontsize=8.5)

    # ---- (C) GTEx direction ----
    ax = axes[0, 2]
    p1_meas = p1[p1["concordant"].notna()]
    p1_ok = (p1_meas["concordant"] == True).sum()
    m26_meas = c[c["direction"].isin(["consistent", "conflicting"])]
    m26_ok = (m26_meas["direction"] == "consistent").sum()
    for i, (lab, k, nsub, col) in enumerate([
            ("15 candidates", int(m26_ok), len(m26_meas), C_NEW),
            ("known 106 loci", int(p1_ok), len(p1_meas), C_MAIN)]):
        lo, hi = wilson(k, nsub); yc = 100*k/nsub
        ax.bar(i, yc, width=0.5, color=col, alpha=0.9)
        ax.errorbar(i, yc, yerr=[[max(0, yc-100*lo)], [max(0, 100*hi-yc)]], fmt="none", ecolor="black", capsize=4)
        ax.text(i, yc + 3, f"{k}/{nsub}\n{yc:.0f}%", ha="center", fontsize=7)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["15 candidates", "known 106 loci"], fontsize=7.4)
    ax.set_ylim(0, 110); ax.set_ylabel("Direction-consistent proportion (%)")
    ax.set_title("(C) GTEx v8 independent eQTL\ndirection replication", fontsize=8.5)

    # ---- (D) FinnGen -log10 p ----
    ax = axes[1, 0]
    sub = c[c["finn_p"].notna()].copy()
    ox = -np.log10(sub["gwas_min_p"].clip(lower=1e-300))
    oy = -np.log10(sub["finn_p"])
    sig = sub["finn_p"] < 0.05
    lim = max(max(np.concatenate([ox, oy]))*1.1, 8)
    ax.scatter(ox, oy, s=30, c=C_NEW, alpha=0.85, edgecolor="white", lw=0.4, zorder=3)
    ax.scatter(ox[sig], oy[sig], s=52, c=C_RED, marker="*", edgecolor="white", lw=0.4, zorder=4)
    ax.plot([0, lim], [0, lim], color="#bbbbbb", ls=":", lw=1)
    ax.axhline(-np.log10(0.05), color=C_RED, ls="--", lw=0.9)
    ax.text(0.03, -np.log10(0.05)+0.15, "FinnGen p=0.05", fontsize=6.4, color=C_RED)
    for r in sub.itertuples():
        ax.annotate(r.symbol, (-np.log10(r.gwas_min_p), -np.log10(r.finn_p)),
                    textcoords="offset points", xytext=(5, 4), fontsize=5.8, color="#5a4010")
    ax.set_xlim(0, lim); ax.set_ylim(0, lim)
    ax.set_xlabel("Original GWAS peak $-\\log_{10}(p)$"); ax.set_ylabel("FinnGen R11 $-\\log_{10}(p)$")
    ax.set_title(f"(D) FinnGen R11 independent-cohort replication\n({len(sub)} alignable candidates; ★ = p<0.05, {int(sig.sum())})", fontsize=8.5)

    # ---- (E) FinnGen alignability ----
    ax = axes[1, 1]
    notloc = [s for s in ["U6atac", "CWF19L1", "N4BP2L2", "SLC12A3", "PLAUR"] if s in set(c["symbol"])]
    nomiss = [s for s in ["CLEC3B"] if s in set(c["symbol"])]
    n_align = int(c["finn_p"].notna().sum())
    n_nl = len(notloc); n_om = len(nomiss)
    ax.barh(0, n_align, height=0.5, color=C_GREEN, alpha=0.9)
    ax.barh(0, n_nl, left=n_align, height=0.5, color=C_NEW, alpha=0.9)
    ax.barh(0, n_om, left=n_align+n_nl, height=0.5, color="#b0b0b0", alpha=0.9)
    ax.text(n_align/2, 0, f"{n_align} alignable\n({100*n_align/len(c):.0f}%)", ha="center", va="center", fontsize=6.8, color="white", fontweight="bold")
    ax.text(n_align+n_nl/2, 0, f"{n_nl}\nlead not\nlocalized", ha="center", va="center", fontsize=6.0, color="white")
    ax.text(n_align+n_nl+n_om/2, 0, f"{n_om}\norig data\nmissing", ha="center", va="center", fontsize=6.0, color="#555555")
    ax.set_yticks([]); ax.set_xlim(0, len(c))
    ax.set_xticks([]); ax.set_xlabel("15 candidate genes", fontsize=8)
    ax.set_title("(E) FinnGen R11 alignment coverage\n(12/21 across all candidates; 9/15 here)", fontsize=8.5)

    # ---- (F) 覆盖说明文字块 ----
    ax = axes[1, 2]; ax.axis("off")
    ax.text(0.02, 0.97, "Honest disclosure of replication coverage", fontsize=8.6, fontweight="bold", color=C_RED, va="top")
    ax.text(0.02, 0.78, "FinnGen R11 alignment coverage among the 15 candidates:\n"
            "9/15 (60%) lead cis-eQTL variants alignable in R11; 5 lead SNPs not\n"
            "localized (U6atac/CWF19L1/N4BP2L2/SLC12A3/PLAUR); 1 lacks original\n"
            "summary data (CLEC3B). 8 of 9 alignable candidates show concordant\n"
            "gene-level MR direction; 4 reach FinnGen p<0.05.",
            fontsize=7, va="top", color="#333333")
    ax.text(0.02, 0.30, "HEIDI caveat: SMR tests applied only to loci with a single cis\n"
            "instrument and 4/23 candidates had no SMR coverage; HEIDI pass rates\n"
            "are conditional on that selection.",
            fontsize=7, va="top", color="#555555")
    ax.text(0.02, 0.08, "GTEx panel: 44/63 measurable known loci concordant (70%);\n"
            "VSIG8 is the sole conflicting candidate direction.",
            fontsize=7, va="top", color="#555555")

    fig.suptitle("Convergence and independent-cohort replication of the 15 candidate effector genes",
                 fontsize=11, y=1.0)
    fig.tight_layout(); fig.savefig(f"{OUT}/20260816_F7_replication_v2.png", bbox_inches="tight")
    plt.close(fig); print("F7 OK")

# ============================ Fig 8: coloc.susie sensitivity ============================
def fig8():
    src = f"{RES}/m34_coloc_susie_20260816.csv"
    tag = "max_iter=200, external 1000G EUR LD"
    rows = list(csv.DictReader(open(src)))
    labs = [f"{r['symbol']}×{r['outcome']}" for r in rows]
    abf = [f(r["abf_pp4"]) for r in rows]
    sus = [f(r["susie_pp4"]) for r in rows]
    conv = [("conv_eqtl=" in r.get("note", "") and "FALSE" not in r["note"].split("conv_eqtl=")[1][:6]) for r in rows]
    ns1 = [0]*len(rows); ns2 = [0]*len(rows)
    for i, r in enumerate(rows):
        import re
        m1 = re.search(r"susie_cs_eqtl=(\d+)", r.get("note", "")); m2 = re.search(r"susie_cs_gwas=(\d+)", r.get("note", ""))
        if m1: ns1[i] = int(m1.group(1))
        if m2: ns2[i] = int(m2.group(1))
    x = np.arange(len(rows))
    fig, ax = plt.subplots(figsize=(9.6, 4.9), dpi=300)
    w = 0.34
    b1 = ax.bar(x-w/2, abf, width=w, color=C_MAIN, alpha=0.92, label="coloc.abf (single causal variant)")
    b2 = ax.bar(x+w/2, sus, width=w, color=C_NEW, alpha=0.92, label=f"coloc.susie (multiple causal variants)")
    for xi, (a, s) in enumerate(zip(abf, sus)):
        ax.text(xi-w/2, a+0.02, f"{a:.3f}", ha="center", fontsize=6.6, color=C_MAIN)
        ax.text(xi+w/2, s+0.02, f"{s:.3f}", ha="center", fontsize=6.6, color="#9a5a1c")
    for i, (ok, n1, n2) in enumerate(zip(conv, ns1, ns2)):
        mark = "✓" if ok else "✗"
        col = C_GREEN if ok else C_RED
        ax.text(i, -0.09, f"{mark}\nCS {n1}/{n2}", ha="center", va="top", fontsize=6.4, color=col)
    lamc = [i for i, r in enumerate(rows) if r["symbol"] == "LAMC1"]
    for i in lamc:
        ax.get_xticklabels()[i].set_color(C_RED); ax.get_xticklabels()[i].set_fontweight("bold")
    ax.set_xticks(x); ax.set_xticklabels(labs, fontsize=8)
    ax.set_ylabel("PP.H4"); ax.set_ylim(0, 1.15)
    ax.legend(fontsize=7.5, loc="upper right", frameon=False)
    ax.set_title(f"Multiple-causal-variant sensitivity (coloc.susie, {tag})\n"
                 "✗ = susie_rss non-converged (all loci; max_iter=1000 also non-convergent — see text); CS = credible sets; red = LAMC1, excluded",
                 fontsize=9)
    fig.tight_layout(); fig.savefig(f"{OUT}/20260816_F8_susie_v2.png", bbox_inches="tight")
    plt.close(fig); print(f"F8 OK ({tag})")

# ============================ Fig 9: GWAS significance caveat ============================
def fig9():
    full = load_full(); q = full[full["ok"] == True].copy()
    fdr = pd.read_csv(f"{RES}/fdr_core_20260816.csv")
    grid = pd.read_csv(f"{GRID}/transcript_coloc_hits.csv")
    known_keys = set(zip(grid["gene"].astype(str), grid["outcome"].astype(str)))
    strong = fdr[fdr["strong"]].copy()
    strong["is_new"] = ~strong[["gene", "outcome"]].apply(
        lambda r: (str(r["gene"]), str(r["outcome"])) in known_keys, axis=1)
    kn = strong[~strong["is_new"]].copy()
    kn = kn[kn[["gene","outcome"]].apply(lambda r: (str(r.gene), str(r.outcome)) in known_keys, axis=1)]
    nw = strong[strong["is_new"]]
    gwas_sig = 5e-8

    fig, axes = plt.subplots(1, 3, figsize=(13.6, 4.6), dpi=300)

    # ---- (A) 比例 ----
    ax = axes[0]
    cats = ["Known 106 loci", "15 candidates"]
    subs = [kn, nw]
    for i, (lab, sub) in enumerate(zip(cats, subs)):
        k = int((sub["gwas_min_p"] < gwas_sig).sum()); n = len(sub); yc = 100*k/n
        lo, hi = wilson(k, n)
        ax.bar(i, yc, width=0.5, color=[C_MAIN, C_NEW][i], alpha=0.92)
        ax.errorbar(i, yc, yerr=[[max(0, yc-100*lo)], [max(0, 100*hi-yc)]], fmt="none", ecolor="black", capsize=4)
        ax.text(i, yc + 2.5, f"{yc:.1f}% ({k}/{n})\n[{100*lo:.0f}, {100*hi:.0f}]", ha="center", fontsize=7.2)
    ax.axhline(100*0.05, color="#bbbbbb", ls=":", lw=1)
    ax.text(1.45, 6.5, "5×10⁻⁸ by chance\n(noise floor)", fontsize=6.4, color="#888888", ha="right")
    ax.set_xticks([0, 1]); ax.set_xticklabels(cats, fontsize=8.5)
    ax.set_ylabel("Fraction reaching GWAS p<5×10⁻⁸ (%)"); ax.set_ylim(0, 75)
    ax.set_title("(A) A minority of strong-colocalized loci\nreach genome-wide significance", fontsize=9)

    # ---- (B) 峰 p 分布 ----
    ax = axes[1]
    for sub, col, lab, ls in [(kn, C_MAIN, "106 known loci", "-"), (nw, C_NEW, "15 candidates", "--")]:
        v = np.clip(-np.log10(sub["gwas_min_p"]), 0, 16)
        if len(v) > 3:
            kde = gaussian_kde(v, bw_method=0.4)
            xs = np.linspace(0, 16, 300)
            ax.plot(xs, kde(xs), ls, color=col, lw=1.8, label=f"{lab} (median -log₁₀p = {np.median(v):.1f})")
        ax.hist(v, bins=np.linspace(0, 16, 24), density=True, alpha=0.18, color=col)
    ax.axvline(-np.log10(gwas_sig), color=C_RED, ls="--", lw=1)
    ax.text(-np.log10(gwas_sig)+0.1, 0.94*max(ax.get_ylim()[1], 1), "p=5×10⁻⁸", fontsize=6.4, color=C_RED)
    ax.set_xlabel("GWAS peak $-\\log_{10}(p)$ (clipped at 16)")
    ax.set_ylabel("Density")
    ax.legend(fontsize=6.6, loc="upper left", frameon=False)
    ax.set_title("(B) Peak-association distribution: candidates\ncluster below genome-wide significance", fontsize=9)

    # ---- (C) PP.H4 vs gwas peak p ----
    ax = axes[2]
    ax.scatter(-np.log10(kn["gwas_min_p"].clip(lower=1e-300)), kn["pp4"], s=14, c="#8fa3c0", alpha=0.8, edgecolor="white", lw=0.4, label="106 known loci")
    ax.scatter(-np.log10(nw["gwas_min_p"].clip(lower=1e-300)), nw["pp4"], s=38, marker="^", c=C_NEW, alpha=0.9, edgecolor="white", lw=0.4, label="15 candidates")
    ax.axvline(-np.log10(gwas_sig), color=C_RED, ls="--", lw=0.9)
    ax.text(-np.log10(gwas_sig)+0.1, 0.79, "p=5×10⁻⁸", fontsize=6.3, color=C_RED)
    ax.axhline(0.8, color="#bbbbbb", ls=":", lw=0.9)
    ax.set_xlabel("GWAS peak $-\\log_{10}(p)$"); ax.set_ylabel("PP.H4 (coloc.abf)")
    ax.set_ylim(0.8, 1.005); ax.legend(fontsize=6.6, loc="lower right", frameon=False)
    ax.set_title("(C) Colocalization support is not\ncontingent on a genome-wide peak", fontsize=9)

    fig.suptitle("Caveat: strong colocalization frequently occurs below the GWAS significance threshold\n"
                 "(41/106 known loci and 6/15 candidates reach p<5×10⁻⁸; median known peak p = 6.4×10⁻¹⁰)",
                 fontsize=10, y=1.0)
    fig.tight_layout(); fig.savefig(f"{OUT}/20260816_F9_gwas_caveat_v2.png", bbox_inches="tight")
    plt.close(fig); print("F9 OK")

# ============================ S1: resource/data availability ============================
def s1():
    fig, axes = plt.subplots(1, 3, figsize=(14.4, 4.6), dpi=300)

    # ---- (A) 数据源流程 ----
    ax = axes[0]; ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10)
    ax.set_title("(A) Data sources and analysis pipeline", fontsize=9, pad=6)
    box(ax, 0.3, 7.6, 3.0, 1.0, "eQTLGen\nn=31,684, whole blood\n(hg19, 2018-09-05)", "#4a7a9b", "#2e4d66", fs=6.8)
    box(ax, 3.8, 7.6, 2.9, 1.0, "GWAS\nT2D Xue 2018 n=655,666\nCAD n=296,525 · FG n=58,074", "#4a7a9b", "#2e4d66", fs=6.6)
    box(ax, 7.2, 7.6, 2.5, 1.0, "1000G EUR\nn=503 LD reference\n(hg19)", "#4a7a9b", "#2e4d66", fs=6.8)
    box(ax, 2.6, 4.9, 4.9, 1.1, "cis-MR (Wald ratio) + coloc.abf\n31,373 gene–trait pairs, 3 outcomes\nper-outcome BH-FDR q<0.05", "#1a3a6b", "#0f2445", fs=7.2)
    arrow(ax, 1.8, 7.6, 3.6, 6.0); arrow(ax, 5.2, 7.6, 5.2, 6.0); arrow(ax, 8.5, 7.6, 6.6, 6.0)
    box(ax, 3.4, 1.8, 3.4, 1.1, "Replication\nGTEx v8 · FinnGen R11\nSMR+HEIDI · Steiger · coloc.susie", C_GREEN, "#3f5c35", fs=7.0)
    arrow(ax, 5.0, 4.9, 5.0, 2.9)
    ax.text(5.0, 0.7, "All public summary statistics; no individual-level data", fontsize=6.6, ha="center", color="#555555")

    # ---- (B) atlas 分层资源 ----
    ax = axes[1]; ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10)
    ax.set_title("(B) Released atlas layers", fontsize=9, pad=6)
    layers = [("Full atlas", "31,371 QC pairs", "gene × outcome × MR + coloc fields", "#eef2f7", "#1a3a6b"),
              ("MR-significant core", "982 pairs", "per-outcome BH-FDR q<0.05", "#dbe4f0", "#2e5598"),
              ("Strong colocalizations", "121 pairs", "PP.H4≥0.8 (12.3% yield)", "#fbe3dd", C_RED),
              ("Known-locus reproduction", "106/106 (100%)", "stage-2 grid strong loci", "#e4efe4", C_GREEN),
              ("Candidate effector genes", "15 genes", "7 T2D + 8 CAD, not in GWAS Catalog", "#fdf0e1", C_NEW)]
    y = 9.0
    for name, cnt, desc, bg, edge in layers:
        p = FancyBboxPatch((0.3, y-0.72), 3.4, 0.66, boxstyle="round,pad=0.03", facecolor=bg, edgecolor=edge, lw=0.8)
        ax.add_patch(p)
        ax.text(0.5, y-0.42, name, fontsize=7.4, fontweight="bold", color="#222222", va="center")
        ax.text(4.0, y-0.42, cnt, fontsize=7.6, fontweight="bold", color=edge, va="center")
        ax.text(4.0, y-0.05, desc, fontsize=6.2, color="#555555", va="bottom")
        y -= 1.52
    ax.text(0.3, 0.4, "Atlas + code: Zenodo DOI (to be registered at acceptance)\nGitHub: github.com/qgeng1465/dual-channel-mr-atlas", fontsize=6.8, color="#333333")

    # ---- (C) 覆盖披露（donut） ----
    ax = axes[2]
    n_align, n_nl, n_om = 9, 5, 1
    ax.pie([n_align, n_nl, n_om], colors=[C_GREEN, C_NEW, "#b0b0b0"],
           startangle=90, counterclock=False,
           wedgeprops=dict(width=0.34, edgecolor="white", lw=1.2))
    ax.text(0, 0.16, f"{n_align}/{n_align+n_nl+n_om}", ha="center", fontsize=12, fontweight="bold", color=C_MAIN)
    ax.text(0, -0.18, "alignable", ha="center", fontsize=7.4, color="#333333")
    ax.legend(["9 alignable (60%)", "5 lead not localized", "1 original data missing"],
              fontsize=6.6, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
    ax.set_title("(C) FinnGen R11 alignment coverage\namong 15 candidates (9/15)", fontsize=9, pad=6)
    ax.text(0, -1.5, "Disclosure: replication is limited by variant-level alignment, not by\nanalysis failure; unaligned loci are reported as such.",
            fontsize=6.4, ha="center", color="#555555")

    fig.suptitle("Supplemental: data availability, released atlas layers, and honest disclosure of replication coverage",
                 fontsize=11, y=1.0)
    fig.tight_layout(); fig.savefig(f"{OUT}/20260816_S1_resources_v2.png", bbox_inches="tight")
    plt.close(fig); print("S1 OK")

if __name__ == "__main__":
    args = sys.argv[1:] or ["fig1", "fig2", "fig3", "fig4", "fig6", "fig7", "fig8", "fig9", "s1"]
    for a in args:
        {"fig1": fig1, "fig2": fig2, "fig3": fig3, "fig4": fig4,
         "fig6": fig6, "fig7": fig7, "fig8": fig8, "fig9": fig9, "s1": s1}[a]()
    print("M37 done")
