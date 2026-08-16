#!/usr/bin/env python3
# FigS2: MR-Steiger corrected direction (M44) + empirical-Bayes coloc enrichment (M45)
# 输出：results/figures/20260817_FigS2_steiger_ebayes.png
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = "/data/qiushuogeng/projects/dual-channel-mr-atlas"
FIG = os.path.join(REPO, "results", "figures")
os.makedirs(FIG, exist_ok=True)

st = pd.read_csv(os.path.join(REPO, "results/m44_steiger_corrected_20260817.csv"))
eb = pd.read_csv(os.path.join(REPO, "results/m45_ebayes_20260817.csv"))
summ = pd.read_csv(os.path.join(REPO, "results/m45_ebayes_summary_20260817.csv")).iloc[0]

fig, axes = plt.subplots(1, 2, figsize=(13, 5.6), gridspec_kw={"width_ratios": [1, 1.05]})

# --- Panel A: MR-Steiger r2_exposure vs r2_outcome ---
ax = axes[0]
order = st.sort_values("r2_exposure", ascending=False)["symbol"].tolist()
y = np.arange(len(order))
w = 0.36
for i, sym in enumerate(order):
    r = st[st["symbol"] == sym].iloc[0]
    ax.barh(y[i] + w/2, r["r2_exposure"], height=w, color="#1a3a6b", label="r2_exp" if i == 0 else "")
    ax.barh(y[i] - w/2, r["r2_outcome"], height=w, color="#c00000", label="r2_out" if i == 0 else "")
ax.set_yticks(y); ax.set_yticklabels(order, fontsize=9)
ax.set_xscale("log")
ax.set_xlabel("r² (variance explained)", fontsize=11)
ax.set_title("MR-Steiger corrected direction\n15 candidates: 15/15 r2_exp > r2_out", fontsize=12)
ax.legend(fontsize=9, loc="lower right")
ax.axvline(1, color="0.5", lw=0.6, ls=":")
# annotate sig steiger
for i, sym in enumerate(order):
    r = st[st["symbol"] == sym].iloc[0]
    if r["steiger_pval"] < 0.05:
        ax.text(max(r["r2_exposure"], r["r2_outcome"]) * 1.5, y[i], "*", ha="left", va="center", fontsize=10, color="0.25")
ax.text(0.002, len(order) - 0.4, "* corrected Steiger p < 0.05 (13/15)",
        fontsize=8, color="0.4", va="center")

# --- Panel B: EB enrichment prior PP.H4 vs PP.H4_enr ---
ax = axes[1]
ax.scatter(eb["pp4"], eb["pp4_enr"], s=14, alpha=0.6, color="#1a3a6b",
           edgecolors="none", label="121 strong pairs")
# identity + enr thresholds
lim = [0.6, 1.001]
ax.plot(lim, lim, ls="--", color="0.4", lw=0.8, label="identity")
ax.plot(lim, [0.8, 0.8], ls=":", color="0.4", lw=0.8)
ax.text(0.995, 0.8, "0.8", ha="right", va="bottom", fontsize=8, color="0.4")
ax.set_xlabel("PP.H4 (coloc.abf)", fontsize=11)
ax.set_ylabel("PP.H4_enr (empirical-Bayes prior)", fontsize=11)
ax.set_title(f"Empirical-Bayes coloc enrichment\nρ = {summ['rho_enrichment']:.1f} "
             f"[{summ['rho_ci']}] ; 121/121 retained", fontsize=12)
ax.legend(fontsize=9)
ax.set_xlim(0.6, 1.001); ax.set_ylim(0.6, 1.001)
ax.text(0.63, 0.985, f"background π0 = {summ['pi0_global']:.4f}\nMR-sig rate = {summ['pi_mr']:.3f}",
        fontsize=8, color="0.4")

fig.tight_layout()
fig.savefig(os.path.join(FIG, "20260817_FigS2_steiger_ebayes.png"), dpi=300, bbox_inches="tight")
plt.close(fig)
print("saved FigS2 (MR-Steiger + EB enrichment)")
