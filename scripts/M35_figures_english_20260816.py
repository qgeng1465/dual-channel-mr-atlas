#!/usr/bin/env python3
# =============================================================================
# M35_figures_english_20260816.py — English-label versions of the main figures
# for the AJHG manuscript submission. Same data, same layout, English text.
#   EN-Fig1 = F1_workflow_schematic_EN
#   EN-Fig2 = F4_fuji_genome_EN            (manuscript Fig 2)
#   EN-Fig3 = F2_precision_funnel_EN       (manuscript Fig 3)
#   EN-Fig4 = F3_pph4_panel_EN             (manuscript Fig 4)
#   EN-Fig6 = F6_convergence_EN            (manuscript Fig 7)
#   EN-S1   = F7_resource_EN               (supplemental Fig S1)
# F5 (locuszoom mirror, M32) and the candidate lollipop (M33b) are already English.
# =============================================================================
import csv, json, os
from math import sqrt
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

RES = "<repo-root>/results"
GRID = f"{RES}/grid"
OUT = f"{RES}/figures"
os.makedirs(OUT, exist_ok=True)

C_MAIN = "#1a3a6b"; C_RED = "#c00000"; C_GREY = "#7a7a7a"
C_NEW = "#d9822b"; C_GREEN = "#6b8e5a"
T2D_C = "#d9822b"; CAD_C = "#1a3a6b"; FBG_C = "#6b8e5a"

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

# ============================ Fig 1: workflow ============================
def fig1():
    fig, ax = plt.subplots(figsize=(10.5, 6.4), dpi=300)
    ax.axis("off"); ax.set_xlim(0, 10); ax.set_ylim(0, 10)

    def box(x, y, w, h, text, fc, ec, fs=8.2, tc="white", bold=True):
        p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05",
                           facecolor=fc, edgecolor=ec, lw=1.0)
        ax.add_patch(p)
        ax.text(x+w/2, y+h/2, text, ha="center", va="center", fontsize=fs,
                color=tc, fontweight="bold" if bold else "normal", linespacing=1.3)

    def arrow(x1, y1, x2, y2, color="#333333"):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                     mutation_scale=13, lw=1.4, color=color))

    box(0.4, 8.7, 2.8, 1.1, "eQTLGen whole-blood\ncis-eQTL (n=31,684)", "#4a7a9b", "#2e4d66")
    box(3.7, 8.7, 2.8, 1.1, "GWAS: T2D / CAD / FG\nn=655,666 / 296,525 / 58,074", "#4a7a9b", "#2e4d66")
    box(7.0, 8.7, 2.6, 1.1, "LD reference 1000G\nEUR (n=503)", "#4a7a9b", "#2e4d66")
    box(3.3, 6.9, 3.4, 1.15, "31,373 gene–trait pairs\nwhole-transcriptome cis-MR × coloc scan", "#1a3a6b", "#0f2445")
    arrow(5.0, 8.7, 5.0, 8.05)
    box(3.3, 5.35, 3.4, 1.1, "4,248 MR-significant pairs\n(mr_p < 0.05, Wald ratio)", "#2e5598", "#1a3a6b")
    arrow(5.0, 6.9, 5.0, 6.45)
    box(3.3, 3.85, 3.4, 1.1, "coloc.abf (p12=1e-5)\n131 strong colocalizations (PP.H4≥0.8)", C_RED, "#7a0000")
    arrow(5.0, 5.35, 5.0, 4.95)
    box(0.5, 2.2, 3.4, 1.2, "106 known T2D/CAD loci\nall reproduced (100%)", C_GREEN, "#3f5c35")
    box(6.1, 2.2, 3.4, 1.2, "23 novel candidate effector genes\n15 known-locus + 8 weak-locus", C_NEW, "#9a5a1c")
    arrow(4.2, 4.25, 2.2, 3.4, C_GREEN)
    arrow(5.8, 4.25, 7.8, 3.4, C_NEW)
    box(6.1, 0.5, 3.4, 1.2, "Independent replication\nGTEx v8: 11/12 direction-consistent\nFinnGen R11: 12/12 gene-level, 5 p<0.05\n(alignment coverage 12/21 = 57%)", "#5a5a5a", "#3a3a3a")
    arrow(7.8, 2.2, 7.8, 1.7)

    ax.set_title("Transcriptome-wide cis-MR × coloc atlas: workflow and operating-characteristic calibration",
                 fontsize=10, pad=12)
    fig.tight_layout(); fig.savefig(f"{OUT}/20260816_F1_workflow_EN.png", bbox_inches="tight")
    plt.close(fig)
    print("EN Fig 1 OK")

# ============================ Fig 2: precision funnel ============================
def fig2():
    rows = list(csv.DictReader(open(f"{RES}/m27_precision_funnel_20260816.csv")))
    m = rows[:8]
    xs = [f(r["n"]) for r in m]
    ys = [f(r["strong_rate"])*100 for r in m]
    ks = [int(f(r["strong"])) for r in m]
    cis = [wilson(k, n) for k, n in zip(ks, xs)]
    xlabs = ["<0.5", "<0.05", "<0.01", "<0.005", "<0.001", "<0.0005", "<1e-4", "<1e-5"]

    fig, ax = plt.subplots(figsize=(8.2, 5.2), dpi=300)
    x = np.arange(len(xs))
    yerr_lo = [max(0, y - 100*lo) for lo, hi, y in zip([c[0] for c in cis], [c[1] for c in cis], ys)]
    yerr_hi = [max(0, 100*hi - y) for lo, hi, y in zip([c[0] for c in cis], [c[1] for c in cis], ys)]
    ax.errorbar(x, ys, yerr=[yerr_lo, yerr_hi],
                fmt="o-", color=C_MAIN, ecolor=C_MAIN, elinewidth=1.2, capsize=3, ms=6, zorder=3)
    for i, (xi, y, n) in enumerate(zip(x, ys, xs)):
        ax.text(xi, y + (4.0 if i < 6 else 2.2), f"{y:.1f}%", ha="center", fontsize=7, color=C_MAIN)
    ax.set_xticks(x); ax.set_xticklabels(xlabs, fontsize=8)
    ax.set_xlabel("cis-MR significance threshold (mr_p < X)", fontsize=9)
    ax.set_ylabel("Strong-colocalization rate (PP.H4≥0.8, %)", fontsize=9)
    ax.set_title("Operating-characteristic curve of single-instrument cis-MR:\n"
                 "tightening the threshold monotonically calibrates colocalization support", fontsize=10)
    ax.axhline(0.418, color="#bbbbbb", ls=":", lw=1)
    ax.text(7.6, 1.0, "All 31,371 pairs: 0.42%", fontsize=7, color="#777777", ha="right")
    ax.axhline(12.96, color=C_RED, ls="--", lw=1.2)
    ax.text(7.6, 13.6, "stage-2 clumped + IVW grid: 12.96%\n[10.83, 15.43]", fontsize=7, color=C_RED, ha="right")
    ax.text(0.02, 0.97, "Within MR-significant set: 3.04% [2.56, 3.60]\nOutside MR-significant (14,294 pairs): 0.014% (only 2)",
            transform=ax.transAxes, fontsize=7, va="top", color="#333333",
            bbox=dict(boxstyle="round,pad=0.3", fc="#f0f4fa", ec="#cccccc", lw=0.5))
    ax.set_ylim(0, 30)
    fig.tight_layout(); fig.savefig(f"{OUT}/20260816_F2_precision_funnel_EN.png", bbox_inches="tight")
    plt.close(fig)
    print("EN Fig 2 OK")

# ============================ Fig 3: PP.H4 panel ============================
def fig3():
    layers = {"sig": [], "grey": [], "null": []}
    for out in ["t2d", "cad", "fbg"]:
        for r in csv.DictReader(open(f"{RES}/coloc_full_{out}_20260815.csv")):
            p4 = f(r["pp4"]); mp = f(r["mr_p"])
            if np.isnan(p4) or np.isnan(mp): continue
            if mp < 0.05: layers["sig"].append(p4)
            elif mp < 0.5: layers["grey"].append(p4)
            else: layers["null"].append(p4)

    def ecdf(x):
        xs = np.sort(x); ys = np.arange(1, len(xs)+1)/len(xs)
        return xs, ys

    fig, axes = plt.subplots(1, 3, figsize=(12.2, 3.9), dpi=300, gridspec_kw={"width_ratios": [1.35, 1, 1]})
    ax = axes[0]
    for key, lab, col, ls in [("sig", f"MR-significant (mr_p<0.05)  n={len(layers['sig']):,}", C_MAIN, "-"),
                              ("grey", f"MR-null (0.05≤p<0.5)  n={len(layers['grey']):,}", C_GREY, "--"),
                              ("null", f"MR-negative (p≥0.5)  n={len(layers['null']):,}", "#b0b0b0", ":")]:
        xs, ys = ecdf(np.array(layers[key]))
        ax.plot(xs, ys, ls, color=col, lw=1.8, label=lab)
    for t in (0.75, 0.8):
        ax.axvline(t, color=C_RED, lw=1, ls="-." if t == 0.75 else "-")
        ax.text(t, 0.94, f"PP.H4={t}", color=C_RED, fontsize=6.5, rotation=90, ha="right")
    ax.set_xlabel("PP.H4 (coloc.abf)"); ax.set_ylabel("ECDF (cumulative proportion)")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.legend(fontsize=6.5, loc="lower right", frameon=False)
    ax.set_title("(A) PP.H4 distribution by MR status", fontsize=9)
    ax = axes[1]
    n_sig, n_grey, n_null = map(len, (layers["sig"], layers["grey"], layers["null"]))
    for t, name in [(0.8, "≥0.8"), (0.75, "≥0.75")]:
        ks = [sum(1 for v in layers["sig"] if v >= t),
              sum(1 for v in layers["grey"] if v >= t),
              sum(1 for v in layers["null"] if v >= t)]
        off = -0.17 if t == 0.75 else 0.03
        for i, (k, n) in enumerate(zip(ks, [n_sig, n_grey, n_null])):
            lo, hi = wilson(k, n)
            yc = 100*k/n
            ax.errorbar(i+off, yc, yerr=[[max(0, yc-100*lo)], [max(0, 100*hi-yc)]],
                        fmt="o", ms=5, capsize=2.5, elinewidth=1,
                        color=C_RED if t == 0.8 else "#8b8b8b")
            ax.text(i+off, 100*k/n + 0.7, f"{k} ({100*k/n:.1f}%)",
                    ha="center", fontsize=6.2, color=C_RED if t == 0.8 else "#666666")
    ax.set_xticks(range(3)); ax.set_xticklabels(["MR-sig", "MR-null", "MR-neg"], fontsize=8)
    ax.set_ylabel("PP.H4 hit rate (%)"); ax.set_ylim(0, 8)
    ax.set_title("(B) Threshold hit rates (Wilson 95% CI)", fontsize=9)
    ax.legend(["PP.H4≥0.8", "PP.H4≥0.75"], fontsize=6.5, frameon=False, loc="upper left")
    ax = axes[2]
    k_sig = sum(1 for v in layers["sig"] if v >= 0.8)
    k_grey = sum(1 for v in layers["grey"] if v >= 0.8)
    k_null = sum(1 for v in layers["null"] if v >= 0.8)
    cats = ["MR-sig", "MR-null", "MR-neg"]
    total = [n_sig, n_grey, n_null]
    strong = [k_sig, k_grey, k_null]
    xw = 0.55
    for i, (tot, st) in enumerate(zip(total, strong)):
        ax.bar(i, tot, width=xw, color="#d9dee8", label="All pairs" if i == 0 else None)
        ax.bar(i, st, width=xw, color=C_RED if i == 0 else C_NEW, label="strong (PP.H4≥0.8)" if i == 0 else None)
        ax.text(i, st + tot*0.02, f"{st} strong\n{100*st/tot:.2f}% of {tot}", ha="center", fontsize=6.5, color="#333333")
    ax.set_xticks(range(3)); ax.set_xticklabels(cats, fontsize=8)
    ax.set_ylabel("Number of pairs"); ax.set_ylim(0, max(total)*1.25)
    ax.set_title("(C) MR status → strong conversion", fontsize=9)
    ax.legend(fontsize=6.5, frameon=False, loc="upper right")
    fig.suptitle("Full-atlas PP.H4 distribution: strong colocalization concentrates in the MR-significant layer",
                 fontsize=10, y=1.0)
    fig.tight_layout(); fig.savefig(f"{OUT}/20260816_F3_pph4_panel_EN.png", bbox_inches="tight")
    plt.close(fig)
    print("EN Fig 3 OK")

# ============================ Fig 4: Fuji genome distribution ============================
def fig4():
    bim = []
    for line in open("<repo-root>/data/ldref/1kg.v3/EUR.bim"):
        p = line.split(); bim.append((p[1], int(p[0]), int(p[3])))
    rs2pos = {rs: (ch, pos) for rs, ch, pos in bim}
    RS_FALLBACK_POS = {"rs147526786": (1, 47071399),   # LRRC41
                       "rs17716350":  (15, 41259910)}  # C15orf62
    multi = {}
    gwas_p = {}
    for out in ["t2d", "cad", "fbg"]:
        for r in csv.DictReader(open(f"{RES}/coloc_full_{out}_20260815.csv")):
            if r["ok"] == "TRUE" and f(r["pp4"]) >= 0.8:
                multi.setdefault(r["gene"], set()).add(out)
            if f(r["gwas_min_p"]) == f(r["gwas_min_p"]):
                gwas_p.setdefault((r["gene"], out), f(r["gwas_min_p"]))
    pts = []
    for r in csv.DictReader(open(f"{GRID}/transcript_coloc_hits.csv")):
        if r["tier"] != "strong" or not r["top_snp"].strip(): continue
        rs = r["top_snp"].strip()
        if rs in rs2pos:
            ch, pos = rs2pos[rs]
        elif rs in RS_FALLBACK_POS:
            ch, pos = RS_FALLBACK_POS[rs]
        else:
            continue
        gp = gwas_p.get((r["gene"], r["outcome"]))
        if gp is None or gp <= 0: gp = 1e-300
        pts.append((ch, pos, gp, "known", r["outcome"], None, len(multi.get(r["gene"], {out for out in [r["outcome"]]}))))
    for r in csv.DictReader(open(f"{RES}/m25_new_strong_annotation_20260816.csv")):
        ch, pos = f(r["chr_hg19"]), f(r["pos_hg19"])
        if np.isnan(ch) or np.isnan(pos): continue
        pts.append((int(ch), int(pos), f(r["gwas_min_p"]), "new", r["outcome"], r["symbol"], len(multi.get(r["gene"], {r["outcome"]}))))
    pts.sort(key=lambda x: (x[0], x[1]))
    chrs = sorted({p[0] for p in pts})
    chr_off = {}; cum = 0
    for c in chrs: chr_off[c] = cum; cum += 1
    fig, ax = plt.subplots(figsize=(11.5, 4.6), dpi=300)
    outc_map = {"t2d": (T2D_C, "T2D"), "cad": (CAD_C, "CAD"), "fbg": (FBG_C, "FG")}
    shape_map = {"known": "o", "new": "^"}
    for ch, pos, pval, cat, outcome, label, neff in pts:
        x = chr_off[ch] + pos/2.5e8
        y = -np.log10(max(pval, 1e-300))
        col, _ = outc_map[outcome]
        ax.scatter(x, y, s=16 + 9*neff, marker=shape_map[cat], color=col, alpha=0.88,
                   edgecolor="white", linewidth=0.4, zorder=3)
    new_pts = [p for p in pts if p[3] == "new"]
    for i, (ch, pos, pval, cat, outcome, label, neff) in enumerate(new_pts):
        x = chr_off[ch] + pos/2.5e8; y = -np.log10(max(pval, 1e-300))
        dy = 10 if i % 2 == 0 else 18
        ax.annotate(label, (x, y), textcoords="offset points", xytext=(0, dy),
                    fontsize=6.4, ha="center", color="#5a4010")
    handles = [plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=C_MAIN, ms=6, label="known 106 (replicated loci)"),
               plt.Line2D([0], [0], marker="^", color="w", markerfacecolor=C_NEW, ms=7, label="23 novel candidates")]
    for out, (col, lab) in outc_map.items():
        handles.append(plt.Line2D([0], [0], marker="s", color="w", markerfacecolor=col, ms=5, label=lab))
    handles.append(plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#999999", ms=9, label="size ∝ n strong outcomes"))
    ax.legend(handles=handles, fontsize=7, loc="upper right", frameon=False, ncol=2)
    ax.set_xticks([chr_off[c] + 0.5 for c in chrs])
    ax.set_xticklabels(chrs, fontsize=7)
    ax.set_xlabel("Chromosome (hg19)", fontsize=9); ax.set_ylabel("GWAS peak $-\\log_{10}(p)$", fontsize=9)
    ax.set_title("Genome-wide distribution of 129 MR-significant strong-colocalization pairs (Fuji-style)\n"
                 "Gray-zone pairs AP3S2×T2D and ZNF19×CAD not plotted (see text)", fontsize=10)
    ax.set_ylim(0, None)
    fig.tight_layout(); fig.savefig(f"{OUT}/20260816_F4_fuji_genome_EN.png", bbox_inches="tight")
    plt.close(fig)
    print("EN Fig 4 OK, points:", len(pts))

# ============================ Fig 6: convergence panel ============================
def fig6():
    heidi = json.load(open(f"{GRID}/heidi_full_summary.json"))
    steig = json.load(open(f"{GRID}/steiger_direction_76.json"))
    m28 = list(csv.DictReader(open(f"{RES}/m28_finngen_replication_new23_20260816.csv")))
    gwas_p = {r["gene"]: f(r["gwas_min_p"]) for r in csv.DictReader(open(f"{RES}/m25_new_strong_annotation_20260816.csv"))}
    aligned = [r for r in m28 if (r["mr_replicated"] or "").lower() in ("yes", "true")]
    sig5 = [r for r in aligned if r["finn_p"] and float(r["finn_p"]) < 0.05]

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.6), dpi=300)
    ax = axes[0, 0]
    bars = [("stage-2 full set", heidi["heidi_ok_pre"], heidi["n_pre"]),
            ("strong subset", heidi["heidi_ok_strong"], heidi["n_strong"])]
    for i, (lab, k, n) in enumerate(bars):
        lo, hi = wilson(k, n)
        ax.bar(i, 100*k/n, width=0.5, color=[C_GREY, C_MAIN][i], alpha=0.9)
        ax.errorbar(i, 100*k/n, yerr=[[100*lo], [100*hi]], fmt="none", ecolor="black", capsize=4, elinewidth=1)
        ax.text(i, 100*k/n + 3, f"{100*k/n:.1f}%\n{n} pairs", ha="center", fontsize=7.5)
    ax.set_xticks([0, 1]); ax.set_xticklabels([b[0] for b in bars], fontsize=8)
    ax.set_ylabel("HEIDI pass rate (%)"); ax.set_ylim(0, 95)
    ax.set_title("(A) SMR+HEIDI: higher pass rate in the strong subset", fontsize=9)
    ax = axes[0, 1]
    k, n = steig["fwd"], steig["n"]
    lo, hi = wilson(k, n)
    ax.bar([0], [100*k/n], width=0.45, color=C_GREEN, alpha=0.9)
    ax.errorbar(0, 100*k/n, yerr=[[100*lo], [100*hi]], fmt="none", ecolor="black", capsize=4)
    ax.text(0, 100*k/n + 3, f"{100*k/n:.1f}% [{100*lo:.1f}, {100*hi:.1f}]\n{n} pairs", ha="center", fontsize=7.5)
    ax.text(0.0, 6, f"{steig['rev']} reverse (0 reverse-significant)", ha="center", fontsize=7, color="#555555")
    ax.set_xticks([0]); ax.set_xticklabels(["eQTL → outcome"], fontsize=8)
    ax.set_ylim(0, 115); ax.set_ylabel("Direction-consistent proportion (%)")
    ax.set_title("(B) Steiger: causal-direction validation", fontsize=9)
    ax = axes[1, 0]
    m26r = list(csv.DictReader(open(f"{RES}/m26_gtex_replication_new23_20260816.csv")))
    m26_meas = [r for r in m26r if (r.get("gtex_p") or "").strip() and r["direction"] not in ("", "no_gtex_gwas_variant")]
    m26_ok = sum(1 for r in m26_meas if r["direction"] == "consistent")
    p1r = list(csv.DictReader(open(f"{GRID}/gtex_replication_p1.csv")))
    p1_meas = [r for r in p1r if (r["concordant"] or "").strip().upper() in ("TRUE", "FALSE")]
    p1_ok = sum(1 for r in p1_meas if r["concordant"].strip().upper() == "TRUE")
    for i, (lab, k, n) in enumerate([("novel 23-candidate layer", m26_ok, len(m26_meas)),
                                     ("known 106-locus layer", p1_ok, len(p1_meas))]):
        lo, hi = wilson(k, n)
        ax.bar(i, 100*k/n, width=0.5, color=[C_NEW, C_MAIN][i], alpha=0.9)
        ax.errorbar(i, 100*k/n, yerr=[[100*lo], [100*hi]], fmt="none", ecolor="black", capsize=4)
        ax.text(i, 100*k/n + 3, f"{k}/{n}\n{100*k/n:.1f}%", ha="center", fontsize=7.5)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["novel 23-candidate layer", "known 106-locus layer"], fontsize=8)
    ax.set_ylim(0, 100); ax.set_ylabel("Direction-consistent proportion (%)")
    ax.set_title("(C) GTEx v8 independent eQTL direction replication (Wilson 95% CI)", fontsize=9)
    ax = axes[1, 1]
    ox, oy, cx, cy = [], [], [], []
    for r in aligned:
        fp = f(r["finn_p"]); op = gwas_p.get(r["gene"], f(r["orig_p"]))
        if np.isnan(op) or np.isnan(fp) or op <= 0 or fp <= 0: continue
        ox.append(-np.log10(op)); oy.append(-np.log10(fp))
        if fp < 0.05: cx.append(-np.log10(op)); cy.append(-np.log10(fp))
    lim = max(max(ox+oy+[1])*1.15, 6)
    ax.scatter(ox, oy, s=26, c=C_NEW, alpha=0.85, edgecolor="white", lw=0.4, zorder=3)
    if cx:
        ax.scatter(cx, cy, s=42, c=C_RED, marker="*", edgecolor="white", lw=0.4, zorder=4)
    ax.plot([0, lim], [0, lim], color="#bbbbbb", ls=":", lw=1)
    ax.axhline(-np.log10(0.05), color=C_RED, ls="--", lw=0.9)
    ax.text(0.03, -np.log10(0.05)+0.15, "FinnGen p=0.05", fontsize=6.5, color=C_RED)
    ax.set_xlim(0, lim); ax.set_ylim(0, lim)
    ax.set_xlabel("Original GWAS peak $-\\log_{10}(p)$"); ax.set_ylabel("FinnGen R11 $-\\log_{10}(p)$")
    ax.set_title("(D) FinnGen independent-cohort replication (23 novel candidates)", fontsize=9)
    ax.text(0.98, 0.03, f"alignable {len(aligned)}/21 (57%)\n12/12 gene-level direction concordant\n★ = FinnGen p<0.05 ({len(sig5)})",
            transform=ax.transAxes, fontsize=7, ha="right", va="bottom",
            bbox=dict(boxstyle="round,pad=0.3", fc="#fdf6ec", ec="#d9822b", lw=0.6))
    fig.suptitle("Convergence and independent-replication validation (Wilson 95% CI)", fontsize=10, y=1.0)
    fig.tight_layout(); fig.savefig(f"{OUT}/20260816_F6_convergence_EN.png", bbox_inches="tight")
    plt.close(fig)
    print("EN Fig 6 OK")

# ============================ Fig 7: resource snapshot (supplemental S1) ============================
def fig7():
    fig = plt.figure(figsize=(9.5, 5.6), dpi=300)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.25, 1], height_ratios=[1, 0.8])
    ax = fig.add_subplot(gs[0, 0])
    ax.axis("off"); ax.set_title("(A) Data sources and versions", fontsize=9, loc="left")
    items = [("eQTLGen (n=31,684, whole blood, hg19)", "cis-eQTL + AF, 2018-09-05"),
             ("GWAS T2D (Xue 2018, n=655,666)", "GCST006867"),
             ("GWAS CAD (n=296,525)", "GCST005194"),
             ("GWAS FG (n=58,074)", "GCST005186"),
             ("1000 Genomes EUR (n=503)", "LD reference, hg19"),
             ("GTEx v8 (n=838)", "multi-tissue eQTL replication"),
             ("FinnGen R11", "independent outcome replication")]
    y = 0.95
    for name, ver in items:
        ax.text(0.02, y, name, fontsize=7, va="top", color="#1a3a6b")
        if ver: ax.text(0.98, y, ver, fontsize=6.5, va="top", ha="right", color="#555555")
        y -= 0.14
    ax = fig.add_subplot(gs[0, 1])
    ax.axis("off"); ax.set_title("(B) Atlas resource table (released with the Zenodo DOI)", fontsize=9, loc="left")
    data = [
        ["Full atlas", "31,373 pairs", "gene×trait×MR+coloc full fields"],
        ["Strong subset", "131", "PP.H4≥0.8 (T2D 67 / CAD 60 / FG 4)"],
        ["Known-locus reproduction", "106/106", "100%"],
        ["Novel candidate effector genes", "23", "15 known-locus + 8 weak-locus"],
        ["GTEx replication", "11/12 direction-consistent", "1 conflict (VSIG8)"],
        ["FinnGen replication", "12/12 gene-level", "5 p<0.05; coverage 12/21=57%"],
    ]
    tab = ax.table(cellText=data, colLabels=["Layer", "Count", "Description"], loc="upper center",
                   bbox=[0.0, 0.12, 1.0, 0.88], cellLoc="left")
    tab.auto_set_font_size(False); tab.set_fontsize(7.2)
    for (row, col), cell in tab.get_celld().items():
        if row == 0: cell.set_facecolor("#1a3a6b"); cell.set_text_props(color="white", fontweight="bold")
        elif col == 0: cell.set_facecolor("#eef2f7")
        if row > 0: cell.set_edgecolor("#dddddd")
    ax = fig.add_subplot(gs[1, :])
    ax.axis("off")
    cols = ["gene", "symbol", "outcome", "mr_b", "mr_p", "gwas_min_p", "eqtl_F_max", "nsnp", "pp4"]
    example = next(csv.DictReader(open(f"{RES}/coloc_full_t2d_20260815.csv")))
    col_data = [[c, "float" if c in ("mr_b","mr_p","gwas_min_p","eqtl_F_max","nsnp","pp4") else "str",
                 f"{example[c][:22]}" + ("…" if len(str(example[c])) > 22 else "")] for c in cols]
    tab2 = ax.table(cellText=col_data, colLabels=["Column", "Type", "Example value"], loc="upper center",
                    bbox=[0.0, 0.30, 0.55, 0.60], cellLoc="left")
    tab2.auto_set_font_size(False); tab2.set_fontsize(6.8)
    for (row, col), cell in tab2.get_celld().items():
        if row == 0: cell.set_facecolor("#4a7a9b"); cell.set_text_props(color="white", fontweight="bold")
    ax.text(0.58, 0.62, "Disclosure (replication coverage)", fontsize=8.5, fontweight="bold", color="#c00000")
    ax.text(0.58, 0.50, "FinnGen R11 alignment coverage 12/21 (57%):\n7 lead SNPs not localized in R11\n(PLAUR/SLC12A3/CWF19L1/U6atac/\nN4BP2L2/ZBTB46/ZNF100), 2 missing original data\n(RPL13/CLEC3B), 2 FG no-phenotype skipped",
            fontsize=6.8, va="top", color="#333333")
    ax.text(0.58, 0.06, "Use: directly queryable for any cis-MR target screen. Schema published with the bundle\n(Data Availability). Candidates are hypothesis-generating, not causal findings.",
            fontsize=6.8, color="#555555", va="bottom")
    fig.suptitle("Atlas resource snapshot: data sources, table structure, and independent-replication coverage",
                 fontsize=10, y=1.0)
    fig.tight_layout(); fig.savefig(f"{OUT}/20260816_F7_resource_EN.png", bbox_inches="tight")
    plt.close(fig)
    print("EN Fig 7 OK")

if __name__ == "__main__":
    fig1(); fig2(); fig3(); fig4(); fig6(); fig7()
    print("All English figures done")
