#!/usr/bin/env python3
# =============================================================================
# M43_celltype_enrichment_20260817.py - 细胞类型富集（诚实降级：PanglaoDB 标记库）
# =============================================================================
# 目的：审稿人任务/模块 2（单细胞细胞类型富集）。本仓库无原始 scRNA 矩阵，
#   网络对 HPA/Tabula Sapiens/DESCARTES 均不可达（实测：proteinatlas 反爬、
#   descartes 403、cellxgene 403）。诚实降级为：用 PanglaoDB (Franzén et al. 2019,
#   Nucleic Acids Res) 的人类细胞类型标记库，对 106 known（strong-coloc）与
#   15 候选基因集做 Enrichr 式超几何富集。
#   * 背景 = PanglaoDB 人类标记库的基因全集（标准 marker-enrichment 做法）
#   * 检验 = 超几何（Fisher 精确），BH-FDR 多重校正
#   * 输出每细胞类型的 overlap/p/基因列表，气泡图按 -log10(p) 与 overlap 比例
# 诚实 caveat：这是"细胞类型标记富集"而非单细胞表达特异性的直接测量；PanglaoDB
#   为文献标记库，偏重免疫与已注释细胞；未被富集不等于该细胞类型无关。
# 输入：
#   /data/qiushuogeng/tmp/celltype/panglao.tsv        PanglaoDB 标记库
#   results/grid/transcript_coloc_hits.csv  106 known
#   results/candidate15_replication_20260816.csv  15 candidates
# 输出：
#   results/m43_celltype_20260817.csv      每细胞类型 × 基因集富集
#   results/figures/20260817_FigS_celltype_enrichment.png
# =============================================================================
import os
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = "/data/qiushuogeng/projects/dual-channel-mr-atlas"
SCR = "/data/qiushuogeng/tmp"
FIG = os.path.join(REPO, "results", "figures")
os.makedirs(FIG, exist_ok=True)

pang = pd.read_csv(os.path.join(SCR, "celltype/panglao.tsv"), sep="\t")
hum = pang[pang["species"].str.contains("Hs", na=False)].copy()
print(f"panglao rows={len(pang)} human rows={len(hum)}")
# 背景 = 人类标记库基因全集
universe = set(hum["official gene symbol"].dropna())
print(f"human universe genes={len(universe)}")
# 每个细胞类型的标记集
ct_markers = {ct: set(g["official gene symbol"].dropna())
              for ct, g in hum.groupby("cell type")}
print(f"human cell types={len(ct_markers)}")

known = set(pd.read_csv(os.path.join(REPO, "results/grid/transcript_coloc_hits.csv"))["symbol"].dropna())
cand = set(pd.read_csv(os.path.join(REPO, "results/candidate15_replication_20260816.csv"))["symbol"].dropna())
sets = {"known106": known, "candidate15": cand}

def enrich(gene_set, label):
    rows = []
    for ct, markers in ct_markers.items():
        overlap = gene_set & markers
        if len(overlap) == 0:
            continue
        N = len(universe); m = len(markers); n = len(gene_set); k = len(overlap)
        # 超几何：P(X >= k)
        p = stats.hypergeom.sf(k - 1, N, m, n)
        rows.append(dict(cell_type=ct, overlap=";".join(sorted(overlap)),
                         k=k, n_markers=m, n_genes=n, n_universe=N,
                         frac=k / max(n, 1), pval=p))
    df = pd.DataFrame(rows).sort_values("pval")
    if len(df):
        df["fdr"] = multipletests(df["pval"], method="fdr_bh")[1]
        df["neglog10p"] = -np.log10(df["pval"].clip(lower=1e-300))
    df.insert(0, "set", label)
    return df

res = pd.concat([enrich(g, l) for l, g in sets.items()], ignore_index=True)
res.to_csv(os.path.join(REPO, "results/m43_celltype_20260817.csv"), index=False)

# ---- 摘要 ----
for label in ("known106", "candidate15"):
    sub = res[(res["set"] == label) & (res["fdr"] < 0.05)].sort_values("pval")
    print(f"\n== {label}: {len(sub)}/{res[res['set']==label].shape[0]} cell types FDR<0.05 ==")
    print(sub[["cell_type", "k", "n_markers", "frac", "pval", "fdr"]].head(12).to_string(index=False))
    if len(sub):
        print("top genes:", sub.head(1)["overlap"].iloc[0])

# ---- 气泡图 ----
sel = res[res["fdr"] < 0.05].copy()
if len(sel):
    sel = sel.sort_values(["set", "neglog10p"], ascending=[True, False])
    top = sel.groupby("set").head(10)
    x = np.arange(len(top))
    fig, ax = plt.subplots(figsize=(11, 6.5))
    colors = {"known106": "#c00000", "candidate15": "#1a3a6b"}
    for s in ("known106", "candidate15"):
        d = top[top["set"] == s]
        ax.scatter(d.index, d["neglog10p"], s=80 + 800 * d["frac"],
                   color=colors[s], alpha=0.7, label=s, edgecolors="white", linewidths=0.8)
    ax.set_xticks(np.arange(len(top)))
    ax.set_xticklabels([f"{r['cell_type']}" for _, r in top.iterrows()], rotation=60, ha="right", fontsize=7)
    ax.set_ylabel("-log10(p), BH-FDR<0.05")
    ax.set_title("PanglaoDB cell-type marker enrichment of coloc effector genes\n"
                 "(bubble size = overlap fraction of the gene set)")
    ax.axhline(1.301, color="0.6", ls="--", lw=0.7)
    ax.legend(loc="upper right", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "20260817_FigS_celltype_enrichment.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("\nsaved bubble figure")
else:
    print("\nno FDR<0.05 cell types; no figure")
print("\n== DONE M43 ==")

# ---- 气泡图（始终输出：top 名义，标注 FDR 存活情况）----
top = res.sort_values(["set", "pval"], ascending=[True, True]).groupby("set").head(12)
x = np.arange(len(top))
fig, ax = plt.subplots(figsize=(12, 6.5))
colors = {"known106": "#c00000", "candidate15": "#1a3a6b"}
nfdr = {s: int((res[(res["set"]==s)]["fdr"] < 0.05).sum()) for s in ("known106", "candidate15")}
for s in ("known106", "candidate15"):
    d = top[top["set"] == s]
    sig = d["fdr"] < 0.05
    ax.scatter(d.index, d["neglog10p"], s=60 + 900 * d["frac"],
               color=colors[s], alpha=0.75, edgecolors="white", linewidths=0.8,
               label=f"{s} (FDR<0.05: {nfdr[s]})")
    if sig.any():
        ax.scatter(d.index[sig], d["neglog10p"][sig], s=60 + 900 * d["frac"][sig],
                   facecolors="none", edgecolors="black", linewidths=1.2)
ax.set_xticks(np.arange(len(top)))
ax.set_xticklabels([r["cell_type"] for _, r in top.iterrows()], rotation=60, ha="right", fontsize=7)
ax.set_ylabel("-log10(p) [hypergeometric]")
ax.set_title("PanglaoDB cell-type marker enrichment of coloc effector genes\n"
             "(top nominal; open circles = BH-FDR<0.05; none survive correction)")
ax.legend(loc="upper right", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "20260817_FigS_celltype_enrichment.png"), dpi=300, bbox_inches="tight")
plt.close(fig)
print("saved bubble figure (top nominal, FDR annotated)")
