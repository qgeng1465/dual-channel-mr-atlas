#!/usr/bin/env python3
# =============================================================================
# M33b_candidate_genes_EN_20260816.py — English version of the 23-candidate
# lollipop figure for the AJHG manuscript (Word submission).
# Same data and layout as M29, labels in English only.
# Output: results/figures/20260816_F6_candidate_genes_23_EN.png
# =============================================================================
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

RES = "/data/qiushuogeng/projects/dual-channel-mr-atlas/results"
OUT = f"{RES}/figures"
os.makedirs(OUT, exist_ok=True)

def f(x, d=float('nan')):
    try: return float(x)
    except (TypeError, ValueError): return d

m25 = pd.read_csv(f"{RES}/m25_new_strong_annotation_20260816.csv", dtype=str).fillna("")
m26 = pd.read_csv(f"{RES}/m26_gtex_replication_new23_20260816.csv", dtype=str).fillna("")
m28 = None
if os.path.exists(f"{RES}/m28_finngen_replication_new23_20260816.csv"):
    m28 = pd.read_csv(f"{RES}/m28_finngen_replication_new23_20260816.csv", dtype=str).fillna("")
    m28 = m28[m28["outcome"] != "fbg"] if "outcome" in m28 else m28

m = m25.merge(m26[["gene", "outcome", "direction"]],
              on=["gene", "outcome"], how="left")
if m28 is not None:
    m = m.merge(m28[["gene", "outcome", "mr_replicated"]],
                on=["gene", "outcome"], how="left")
m = m.sort_values("pp4", ascending=True).reset_index(drop=True)

out_color = {"t2d": "#c00000", "cad": "#1a3a6b", "fbg": "#2e7d32"}
out_label = {"t2d": "T2D", "cad": "CAD", "fbg": "FG"}
tier_face = {"known_locus_100kb_novel_gene": "full",
             "known_locus_250kb_novel_gene": "full",
             "no_catalog_t2dcad": "none"}

fig, ax = plt.subplots(figsize=(11.5, 7.2), dpi=200)
y = np.arange(len(m))
for yi, (_, r) in zip(y, m.iterrows()):
    c = out_color.get(r["outcome"], "#666666")
    mfc = tier_face.get(r["tier"], "full")
    ax.scatter(f(r["pp4"]), yi, s=95, color=c, marker="o", zorder=3,
               facecolor=(c if mfc == "full" else "white"),
               edgecolor=c, linewidth=1.4)

ax.set_yticks(y)
ax.set_yticklabels(m["symbol"], fontsize=9.5)
ax.set_xlim(0.79, 1.18)
ax.set_xlabel("coloc PP.H4", fontsize=11)
ax.set_title("Twenty-three additional strong-colocalization candidates\n"
             "outside the previously reported set "
             "(filled = novel effector gene at a known T2D/CAD risk locus; "
             "open = no T2D/CAD catalog record)",
             fontsize=10)
ax.grid(axis="x", alpha=0.3, lw=0.5)
ax.axvline(0.8, color="#999999", ls=":", lw=1)

sep_x = 1.04
ax.axvline(sep_x, color="#bbbbbb", lw=0.8)
col_gtex, col_finn = 1.055, 1.105
ax.text((sep_x + col_gtex) / 2 - 0.006, len(m) - 0.35, "GTEx\n(independent\nEQL direction)",
        fontsize=7.5, color="#333333", va="bottom", ha="center")
if m28 is not None:
    ax.text((col_gtex + col_finn) / 2 + 0.005, len(m) - 0.35, "FinnGen\n(independent cohort)",
            fontsize=7.5, color="#333333", va="bottom", ha="center")
for yi, (_, r) in zip(y, m.iterrows()):
    d = r["direction"]
    if d == "consistent":
        ax.text(col_gtex, yi, "✓", color="#2e7d32", ha="center", va="center", fontsize=11)
    elif d == "conflicting":
        ax.text(col_gtex, yi, "×", color="#c00000", ha="center", va="center", fontsize=11)
    else:
        ax.text(col_gtex, yi, "·", color="#999999", ha="center", va="center", fontsize=11)
    if m28 is not None:
        rd = r.get("mr_replicated", "")
        if rd == "yes":
            ax.text(col_finn, yi, "✓", color="#2e7d32", ha="center", va="center", fontsize=11)
        elif rd == "no":
            ax.text(col_finn, yi, "×", color="#c00000", ha="center", va="center", fontsize=11)
        else:
            ax.text(col_finn, yi, "·", color="#999999", ha="center", va="center", fontsize=11)

leg = [
    Line2D([0], [0], marker="o", color="w", markerfacecolor="#c00000", markeredgecolor="#c00000", ms=8, label="T2D"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="#1a3a6b", markeredgecolor="#1a3a6b", ms=8, label="CAD"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="#2e7d32", markeredgecolor="#2e7d32", ms=8, label="FG"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="white", markeredgecolor="#666666", ms=8,
           label="No T2D/CAD catalog record"),
    Line2D([0], [0], marker="", color="w", ls="",
           label="✓ consistent direction  · × conflicting  · · not testable"),
]
ax.legend(handles=leg, loc="upper center", bbox_to_anchor=(0.42, -0.03), fontsize=8,
          frameon=False, ncol=3)
fig.tight_layout()
fig.savefig(f"{OUT}/20260816_F6_candidate_genes_23_EN.png", dpi=200)
plt.close(fig)
print("English candidate figure written:", f"{OUT}/20260816_F6_candidate_genes_23_EN.png")
