#!/usr/bin/env python3
# FigS1 (replacement): method-robustness of strong-colocalization calls
#   Panel A: coloc.susie ridge-shrunk LD sweep (M34c) — PP.H4 vs w
#   Panel B: hyprcoloc 2-trait cross-check vs coloc.abf (M40)
#   Panel C: GCTA-COJO conditional/joint fine-mapping (M46) — index signals + lead pJ
# 输出：results/figures/20260817_FigS1_method_robustness.png
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = "/data/qiushuogeng/projects/dual-channel-mr-atlas"
FIG = os.path.join(REPO, "results", "figures")
os.makedirs(FIG, exist_ok=True)

ridge = pd.read_csv(os.path.join(REPO, "results/m34c_coloc_susie_ridge_20260817.csv"))
hypr  = pd.read_csv(os.path.join(REPO, "results/m40_hyprcoloc_2trait_20260817.csv"))
cojo  = None
if os.path.exists(os.path.join(REPO, "results/m46_cojo_20260817.csv")):
    cojo = pd.read_csv(os.path.join(REPO, "results/m46_cojo_20260817.csv"))

fig, axes = plt.subplots(1, 3, figsize=(17, 5.8), gridspec_kw={"width_ratios": [1.15, 1, 1]})

# ---------- Panel A: ridge sweep ----------
ax = axes[0]
marker = {"RBM6": "o", "CD101": "s", "CNNM2": "^", "RIC8A": "D", "LAMC1": "v"}
for sym, g in ridge.groupby("symbol"):
    g = g.sort_values("w")
    color = "#c00000" if sym == "PLAUR" else "#1a3a6b"
    lw = 2.2 if sym == "PLAUR" else 1.3
    ax.plot(g["w"], g["susie_pp4"], marker=marker.get(sym, "o"), ms=6,
            color=color, lw=lw, label=sym, clip_on=False)
# annotate non-convergence at w=0
for sym, g in ridge.groupby("symbol"):
    g = g.sort_values("w")
    r0 = g.iloc[0]
    if r0["susie_pp4"] == 0.0 or str(r0["note"]).startswith("FAILED"):
        ax.annotate("✗ not\nconverged", (r0["w"], -0.10), ha="center", fontsize=7.5, color="0.4")
    if r0["susie_pp4"] == 1.0 and "OK" in str(r0["note"]):
        # RIC8A converged at w=0 too
        if r0["conv_eqtl"] == "TRUE" and r0["conv_gwas"] == "TRUE":
            ax.annotate("✓", (r0["w"], -0.10), ha="center", fontsize=7.5, color="0.4")
ax.axhline(0.8, color="0.5", lw=0.8, ls=":")
ax.text(0.205, 0.805, "PP.H4 = 0.8", fontsize=8, color="0.4", ha="right")
ax.set_xlabel("ridge weight w  (R_w = (1−w)R + w·I)", fontsize=11)
ax.set_ylabel("PP.H4 (coloc.susie)", fontsize=11)
ax.set_title("Ridge-shrunk LD sweep (6 loci)\ncoloc.susie recovery", fontsize=12)
ax.set_ylim(-0.2, 1.06); ax.set_xlim(-0.03, 0.24)
ax.legend(fontsize=8, loc="lower right", ncol=1)
# PLAUR w=0.20 NA (no credible sets)
ax.text(0.205, 0.13, "PLAUR w=0.20\nNA (0 CS)", fontsize=7.5, color="#c00000", ha="right")

# ---------- Panel B: hyprcoloc vs coloc.abf ----------
ax = axes[1]
for _, r in hypr.iterrows():
    if pd.notna(r["hypr_pp_shared"]):
        ax.scatter(r["coloc_abf_pp4"], r["hypr_pp_shared"], s=80,
                   color="#1a3a6b", zorder=3, edgecolors="white", linewidths=0.5)
        ax.annotate(r["symbol"], (r["coloc_abf_pp4"], r["hypr_pp_shared"]),
                    textcoords="offset points", xytext=(7, 5), fontsize=9)
    else:
        ax.scatter(r["coloc_abf_pp4"], 0.06, marker="x", s=90, color="#c00000", zorder=4)
        ax.annotate(f"{r['symbol']} (dropped eQTLGen;\nregional P = {r['hypr_regional_prob']:.2f})",
                    (r["coloc_abf_pp4"], 0.06), textcoords="offset points",
                    xytext=(-5, -26), ha="right", fontsize=7.5, color="#c00000")
ax.plot([0.5, 1], [0.5, 1], ls="--", color="0.4", lw=0.8)
ax.axhline(0.8, color="0.5", lw=0.8, ls=":")
ax.axvline(0.8, color="0.5", lw=0.8, ls=":")
ax.set_xlabel("PP.H4 (coloc.abf)", fontsize=11)
ax.set_ylabel("PP.shared (hyprcoloc)", fontsize=11)
ax.set_title("Hyprcoloc 2-trait cross-check\n6 loci, shared-variant posterior", fontsize=12)
ax.set_xlim(0.5, 1.02); ax.set_ylim(-0.02, 1.02)
ax.text(0.83, 0.93, "PLAUR hypr 0.95\nsupports abf", fontsize=8, color="#c00000")

# ---------- Panel C: COJO ----------
ax = axes[2]
if cojo is not None:
    loci = ["RBM6", "CD101", "CNNM2", "PLAUR", "RIC8A", "LAMC1"]
    n_slct = cojo.set_index("symbol")["n_index_p5e8"].reindex(loci).fillna(0)
    y = np.arange(len(loci))[::-1]
    w = 0.36
    ax.barh(y + w/2, n_slct.values, height=w, color="#1a3a6b", label="independent signals (p < 5e-8)")
    n_locus = cojo.set_index("symbol")["n_index_p1e4"].reindex(loci).fillna(0)
    ax.barh(y - w/2, n_locus.values, height=w, color="#c8a95c", label="locus-level (p < 1e-4)")
    ax.set_yticks(y); ax.set_yticklabels(loci, fontsize=10)
    ax.set_xlabel("number of conditional independent signals", fontsize=11)
    ax.set_title("GCTA-COJO conditional/joint fine-mapping\n6 loci (1000G EUR LD)", fontsize=12)
    for i, sym in enumerate(loci):
        r = cojo[cojo["symbol"] == sym]
        if len(r):
            pj = r.iloc[0].get("lead_pJ")
            txt = f"lead pJ = {pj:.2g}" if pd.notna(pj) else "lead not index"
            ax.text(0.03, y[i] - w/2, txt, va="center", fontsize=8, color="0.3")
    ax.legend(fontsize=8, loc="lower right")
else:
    ax.text(0.5, 0.5, "COJO panel pending\n(M46 still running)", ha="center", va="center",
            fontsize=11, color="0.4", transform=ax.transAxes)
    ax.axis("off")

fig.tight_layout()
fig.savefig(os.path.join(FIG, "20260817_FigS1_method_robustness.png"), dpi=300, bbox_inches="tight")
plt.close(fig)
print("saved FigS1_method_robustness (COJO panel:", "OK" if cojo is not None else "pending)")
