#!/usr/bin/env python3
# FigS4: cross-tissue lead-eQTL concordance (M40a, GTEx v8 per-gene lead eQTL)
# 输入：results/m40a_cross_tissue_20260817.csv
# 输出：results/figures/20260817_FigS4_cross_tissue.png
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

REPO = "/data/qiushuogeng/projects/dual-channel-mr-atlas"
FIG = os.path.join(REPO, "results", "figures")
os.makedirs(FIG, exist_ok=True)

df = pd.read_csv(os.path.join(REPO, "results/m40a_cross_tissue_20260817.csv"))

TISSUES = ["Whole_Blood", "Adipose_Subcutaneous", "Adipose_Visceral_Omentum",
           "Liver", "Pancreas", "Muscle_Skeletal", "Artery_Coronary", "Artery_Aorta"]
# 排序：T2D 基因在前
order = df[df["outcome"] == "t2d"]["symbol"].unique().tolist() + \
        df[df["outcome"] == "cad"]["symbol"].unique().tolist()

nlog = df.pivot_table(index="symbol", columns="tissue", values="qval", aggfunc="first")
w50  = df.pivot_table(index="symbol", columns="tissue", values="within_50kb", aggfunc="first")

Z = -np.log10(nlog.loc[order, TISSUES]).values   # NaN = 该组织无该基因 lead-eQTL 条目
M = w50.loc[order, TISSUES].values.astype(object)

fig, ax = plt.subplots(figsize=(11.5, 7.6))
im = ax.imshow(Z, aspect="auto", cmap="YlGnBu", vmin=0, vmax=100,
               interpolation="nearest")
ax.set_xticks(range(len(TISSUES)))
ax.set_xticklabels([t.replace("_", "\n") for t in TISSUES], fontsize=8.5)
ax.set_yticks(range(len(order)))
ax.set_yticklabels(order, fontsize=10)
ax.set_xticks(np.arange(-0.5, len(TISSUES), 1), minor=True)
ax.set_yticks(np.arange(-0.5, len(order), 1), minor=True)
ax.grid(which="minor", color="0.9", linewidth=0.6)
ax.tick_params(which="minor", length=0)

# 一致性标记：实心圆 = 该组织 GTEx lead 落在 eQTLGen lead ±50kb
for i in range(len(order)):
    for j in range(len(TISSUES)):
        if np.isnan(Z[i, j]):
            continue
        if M[i, j] is True:
            ax.add_patch(Circle((j, i), 0.24, facecolor="none", edgecolor="#c00000", lw=1.6, zorder=5))
        elif M[i, j] is False:
            ax.add_patch(Circle((j, i), 0.14, facecolor="none", edgecolor="0.5", lw=0.8, zorder=5))

cbar = fig.colorbar(im, ax=ax, shrink=0.9, pad=0.02)
cbar.set_label("-log10 q  (GTEx v8 lead-eQTL)", fontsize=10)

# 图注
ax.text(-1.5, len(order) + 0.9,
        "Cross-tissue lead-eQTL concordance (15 candidate effector genes × 8 GTEx v8 tissues)\n"
        "Fill = -log10(q) of per-tissue lead eQTL; red ring = tissue lead within ±50 kb of eQTLGen blood lead; "
        "grey dot = not concordant; blank = no lead-eQTL test in tissue",
        fontsize=8.5, color="0.25", va="bottom")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "20260817_FigS4_cross_tissue.png"), dpi=300, bbox_inches="tight")
plt.close(fig)

# 摘要（诚实口径）
sub = df[df["qval"] < 0.05]
genes_sig = sub.groupby("symbol")["tissue"].nunique()
genes_w50 = sub[sub["within_50kb"]].groupby("symbol").size()
print(f"FigS4 saved: {len(order)} genes × {len(TISSUES)} tissues")
print(f"  genes with ≥1 sig tissue: {len(genes_sig)}/{len(order)}")
print(f"  genes with ≥1 concordant (within50) sig tissue: {len(genes_w50)}/{len(order)}")
print(f"  U6atac: {'present' if 'U6atac' in order else 'ABSENT from GTEx lead-eQTL (not testable)'}")
