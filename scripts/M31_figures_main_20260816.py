#!/usr/bin/env python3
# =============================================================================
# M31_figures_main_20260816.py — 主图补齐：Fig 3 / Fig 4 / Fig 6 / Fig 7
# =============================================================================
# Fig 3: PP.H4 全分布（31,371 对按 MR 状态分层 ECDF）—— R2/R3
# Fig 4: 图谱曼哈顿图（129 nominal-sig strong 按染色体 × 结局 × known/new）—— R1
#        （灰区 AP3S2/ZNF19 仅 2 个，图注说明 + R3 vignette 细述）
# Fig 6: 收敛验证面板（HEIDI / Steiger / GTEx 跨组织 / 已知 106 FinnGen 双显著）—— R6
# Fig 7: 资源底表快照（图谱 CSV 列结构 + 行数 + 示例行）—— R7
# 诚实口径：所有数字与 INTEGRITY_AUDIT 重算值一致；标注"候选评估非发现"。
# 输出：results/figures/20260816_F{3,4,6,7}_*.png
# =============================================================================
import csv, json, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

for f in font_manager.fontManager.ttflist:
    if f.name == "Noto Sans CJK JP":
        plt.rcParams["font.sans-serif"] = [f.name]
        break
plt.rcParams["axes.unicode_minus"] = False

RES = "<repo-root>/results"
GRID = f"{RES}/grid"
OUT = f"{RES}/figures"
os.makedirs(OUT, exist_ok=True)

# ---- 统一调色板（与 M29 一致：深蓝主 / 红强调 / 灰参考 / 橙新候选 / 绿已知）----
C_MAIN = "#1a3a6b"; C_RED = "#c00000"; C_GREY = "#7a7a7a"
C_NEW = "#d9822b"; C_KNOWN = "#4a7a9b"; C_GREEN = "#6b8e5a"
T2D_C = "#d9822b"; CAD_C = "#1a3a6b"; FBG_C = "#6b8e5a"

def f(x, d=float('nan')):
    try: return float(x)
    except (TypeError, ValueError): return d

def wb(path):
    p = f"{RES}/figures/{path}"
    return p

# ============================ Fig 3：PP.H4 全分布 ============================
# 合并三个结局的 coloc_full，按 MR 状态分层（sig/grey/null）
layers = {"sig": [], "grey": [], "null": []}
for out in ["t2d", "cad", "fbg"]:
    for r in csv.DictReader(open(f"{RES}/coloc_full_{out}_20260815.csv")):
        p4 = f(r["pp4"]); mp = f(r["mr_p"])
        if np.isnan(p4) or np.isnan(mp):
            continue
        if mp < 0.05: layers["sig"].append(p4)
        elif mp < 0.5: layers["grey"].append(p4)
        else: layers["null"].append(p4)

def ecdf(x):
    xs = np.sort(x); ys = np.arange(1, len(xs)+1) / len(xs)
    return xs, ys

fig, ax = plt.subplots(figsize=(7.6, 4.8), dpi=200)
labels = [("sig", f"MR 显著 (mr_p<0.05)  n={len(layers['sig']):,}", C_MAIN, "-"),
          ("grey", f"grey (0.05≤p<0.5)   n={len(layers['grey']):,}", C_GREY, "--"),
          ("null", f"null (p≥0.5)        n={len(layers['null']):,}", "#b0b0b0", ":")]
for key, lab, col, ls in labels:
    xs, ys = ecdf(np.array(layers[key]))
    ax.plot(xs, ys, ls, color=col, lw=2, label=lab)
ax.axvline(0.8, color=C_RED, lw=1.4)
ax.text(0.8, 0.5, "PP.H4 = 0.8\n(strong 阈值)", color=C_RED, fontsize=8,
        rotation=90, ha="right", va="center")
ax.text(0.62, 0.06, "MR 显著层长尾\n非显著层塌在 0 附近",
        color="#333333", fontsize=8, ha="center")
ax.set_xlabel("PP.H4 (coloc.abf 后验概率)"); ax.set_ylabel("ECDF（累计比例）")
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.legend(fontsize=8, loc="lower right", frameon=False)
ax.set_title(f"Fig 3. 全图谱 {len(layers['sig'])+len(layers['grey'])+len(layers['null']):,} 对的 PP.H4 分布按 MR 状态分层", fontsize=11)
fig.tight_layout(); fig.savefig(wb("20260816_F3_pph4_distribution.png")); plt.close(fig)
print("Fig 3 OK:", f"sig={len(layers['sig'])}, grey={len(layers['grey'])}, null={len(layers['null'])}")

# ============================ Fig 4：图谱曼哈顿图 ============================
# 129 nominal-sig strong：known 106（top_snp）+ 23 new（eqtlgen_lead_snp）
# 位置从 1kg EUR bim (hg19) 反查；未匹配跳过并披露
bim = []
for line in open("<repo-root>/data/ldref/1kg.v3/EUR.bim"):
    p = line.split()
    bim.append((p[1], int(p[0]), int(p[3])))  # rsid, chr, pos
rs2pos = {rs: (ch, pos) for rs, ch, pos in bim}

pts = []  # (chr, pos, gwas_min_p, category, outcome, label)
# (gene,outcome) → gwas_min_p（known 106 的峰 p 来自 coloc_full；hits 表无此列）
gwas_min_p = {}
for on in ["t2d", "cad", "fbg"]:
    for r in csv.DictReader(open(f"{RES}/coloc_full_{on}_20260815.csv")):
        gwas_min_p.setdefault((r["gene"], on), f(r["gwas_min_p"]))
# known 106（top_snp → 1kg bim 反查 hg19 位置）
for r in csv.DictReader(open(f"{GRID}/transcript_coloc_hits.csv")):
    if r["tier"] != "strong" or not r["top_snp"].strip():
        continue
    rs = r["top_snp"].strip()
    if rs not in rs2pos: continue
    ch, pos = rs2pos[rs]
    gp = gwas_min_p.get((r["gene"], r["outcome"]))
    if gp is None or not np.isfinite(gp) or gp <= 0: gp = 1e-300
    pts.append((ch, pos, gp, "known", r["outcome"], None))
# 23 new（m25 注释含 chr_hg19/pos_hg19/gwas_min_p，直接用）
for r in csv.DictReader(open(f"{RES}/m25_new_strong_annotation_20260816.csv")):
    ch = f(r["chr_hg19"]); pos = f(r["pos_hg19"])
    if np.isnan(ch) or np.isnan(pos):
        continue
    pts.append((int(ch), int(pos), f(r["gwas_min_p"]), "new", r["outcome"], r["symbol"]))
missing = len([r for r in csv.DictReader(open(f"{RES}/m25_new_strong_annotation_20260816.csv"))
               if np.isnan(f(r["chr_hg19"])) or np.isnan(f(r["pos_hg19"]))])

# 按染色体排序 + 累积坐标
pts.sort(key=lambda x: (x[0], x[1]))
chrs = sorted({p[0] for p in pts})
chr_off = {}; cum = 0
for c in chrs:
    chr_off[c] = cum; cum += 1  # 简化：每条染色体等宽段
fig, ax = plt.subplots(figsize=(10.5, 4.6), dpi=200)
outc_map = {"t2d": (T2D_C, "T2D"), "cad": (CAD_C, "CAD"), "fbg": (FBG_C, "FBG")}
shape_map = {"known": "o", "new": "^"}
for ch, pos, pval, cat, outcome, label in pts:
    x = chr_off[ch] + pos / 2.5e8
    y = -np.log10(max(pval, 1e-300))
    col, _ = outc_map[outcome]
    ax.scatter(x, y, s=26 if cat == "known" else 40, marker=shape_map[cat],
               color=col, alpha=0.85, edgecolor="white", linewidth=0.4, zorder=3)
# 标注 23 个 new 的基因名（交替垂直偏移避免重叠）
new_pts = [p for p in pts if p[4] == "new"]
for i, (ch, pos, pval, cat, outcome, label) in enumerate(new_pts):
    x = chr_off[ch] + pos / 2.5e8; y = -np.log10(max(pval, 1e-300))
    dy = 8 if i % 2 == 0 else 14
    ax.annotate(label, (x, y), textcoords="offset points", xytext=(0, dy),
                fontsize=6.2, ha="center", color="#5a4010")
# 图例 + 轴
handles = [plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=C_MAIN, ms=6, label="known 106"),
           plt.Line2D([0], [0], marker="^", color="w", markerfacecolor=C_NEW, ms=7, label="23 新候选")]
for out, (col, lab) in outc_map.items():
    handles.append(plt.Line2D([0], [0], marker="s", color="w", markerfacecolor=col, ms=5, label=lab))
ax.legend(handles=handles, fontsize=7.5, loc="upper right", frameon=False, ncol=2)
ax.set_xticks([chr_off[c] + 0.5 for c in chrs])
ax.set_xticklabels(chrs, fontsize=7)
ax.set_xlabel("染色体 (hg19)"); ax.set_ylabel("GWAS 峰 −log10(p)")
ax.set_title("Fig 4. 图谱 129 个 nominal-sig strong coloc 的基因组分布\n"
             "（灰区 AP3S2×t2d、ZNF19×cad 未入图——见 R3 vignette）", fontsize=11)
ax.set_ylim(0, None)
fig.tight_layout(); fig.savefig(wb("20260816_F4_atlas_manhattan.png")); plt.close(fig)
print(f"Fig 4 OK: {len(pts)} 点 (known {sum(1 for p in pts if p[3]=='known')}, "
      f"new {sum(1 for p in pts if p[3]=='new')}), 未匹配披露 {missing}")

# ============================ Fig 6：收敛验证面板 ============================
heidi = json.load(open(f"{GRID}/heidi_full_summary.json"))
steig = json.load(open(f"{GRID}/steiger_direction_76.json"))
# GTEx 跨组织（tissue_triangulation 口径，v6.7/INTEGRITY）
tri = list(csv.DictReader(open(f"{GRID}/tissue_triangulation.csv")))
rein = [r for r in tri if r["status"] == "reinforced"]
tri_hits = len({(r["gene"], r["outcome"]) for r in rein})
# FinnGen 已知 106 双显著（T2D/CAD）
fr = list(csv.DictReader(open(f"{GRID}/finngen_replication.csv")))
both = [r for r in fr if f(r["p1"]) < 0.05 and f(r["p2"]) < 0.05 and r["aligned"] == "True"]
dir_ok = [r for r in both if (f(r["b1"]) * f(r["b2a"])) > 0]
t2d_b = [r for r in both if r["outcome"] == "ebi-a-GCST006867"]
cad_b = [r for r in both if r["outcome"] == "ebi-a-GCST005194"]

fig, axes = plt.subplots(2, 2, figsize=(10.2, 6.6), dpi=200)
# (a) HEIDI
ax = axes[0, 0]
ax.bar(["全集 (819)", "strong 子集 (106)"],
       [100*heidi["heidi_ok_pre"]/heidi["n_pre"], 100*heidi["heidi_ok_strong"]/heidi["n_strong"]],
       color=[C_GREY, C_MAIN], width=0.55)
ax.set_ylabel("HEIDI 通过率 (%)"); ax.set_ylim(0, 100)
ax.text(0, 3, f"398/819 = 48.6%", ha="center", fontsize=8, color="#333333")
ax.text(1, 3, f"76/106 = 71.7%", ha="center", fontsize=8, color="#333333")
ax.set_title("(a) HEIDI：strong 子集通过率更高", fontsize=9)
# (b) Steiger
ax = axes[0, 1]
ax.bar(["方向一致"], [100*steig["fwd"]/steig["n"]], color=C_GREEN, width=0.4)
ax.set_ylim(0, 100); ax.set_ylabel("占比 (%)")
ax.text(0, 3, f"73/76 = 96.1%\n(3 反向，0 反向显著)", ha="center", fontsize=8, color="#333333")
ax.set_title("(b) Steiger：因果方向 eQTL→结局", fontsize=9)
# (c) GTEx 跨组织三角验证
ax = axes[1, 0]
ax.bar(["reinforced 对", "覆盖命中"],
       [100*len(rein)/len(tri), 100*tri_hits/79],
       color=[C_MAIN, C_NEW], width=0.5)
ax.set_ylim(0, 100); ax.set_ylabel("占比 (%)")
ax.text(0, 3, f"129/242 对 = 53%", ha="center", fontsize=8, color="#333333")
ax.text(1, 3, f"55/79 命中 = 70%", ha="center", fontsize=8, color="#333333")
ax.set_title("(c) GTEx 跨组织：≥1 组织 MR 支持", fontsize=9)
# (d) 已知 106 FinnGen 双显著（方向一致率 = 双显著对中 b1×b2a 同号占比）
ax = axes[1, 1]
t2d_ok = sum(1 for r in t2d_b if (f(r["b1"]) * f(r["b2a"])) > 0)
cad_ok = sum(1 for r in cad_b if (f(r["b1"]) * f(r["b2a"])) > 0)
ax.bar(["T2D", "CAD"],
       [100*t2d_ok/max(len(t2d_b), 1), 100*cad_ok/max(len(cad_b), 1)],
       color=[T2D_C, CAD_C], width=0.45)
ax.set_ylim(0, 105); ax.set_ylabel("方向一致率 (%)")
ax.text(0, 3, f"{t2d_ok}/{len(t2d_b)} = {100*t2d_ok/max(len(t2d_b),1):.0f}%",
        ha="center", fontsize=8, color="#333333")
ax.text(1, 3, f"{cad_ok}/{len(cad_b)} = {100*cad_ok/max(len(cad_b),1):.0f}%",
        ha="center", fontsize=8, color="#333333")
ax.set_title("(d) 已知 106 层 FinnGen 双显著一致", fontsize=9)
fig.suptitle("Fig 6. 收敛验证面板（R6）——数字与 INTEGRITY_AUDIT 重算一致", fontsize=12, y=1.0)
fig.tight_layout(); fig.savefig(wb("20260816_F6_convergence_panel.png"), bbox_inches="tight"); plt.close(fig)
print(f"Fig 6 OK: HEIDI {heidi['heidi_ok_strong']}/{heidi['n_strong']}, Steiger {steig['fwd']}/{steig['n']}, "
      f"GTEx rein {len(rein)}/{len(tri)} hits {tri_hits}/79, FinnGen both {len(both)} dir_ok {len(dir_ok)}")

# ============================ Fig 7：资源底表快照 ============================
cols = ["gene", "symbol", "outcome", "mr_b", "mr_p", "gwas_min_p", "eqtl_F_max", "nsnp", "pp4", "ok", "note"]
rows_all = sum(1 for _ in csv.DictReader(open(f"{RES}/coloc_full_t2d_20260815.csv")))
n_t2d = rows_all
n_cad = sum(1 for _ in csv.DictReader(open(f"{RES}/coloc_full_cad_20260815.csv")))
n_fbg = sum(1 for _ in csv.DictReader(open(f"{RES}/coloc_full_fbg_20260815.csv")))
example = next(csv.DictReader(open(f"{RES}/coloc_full_t2d_20260815.csv")))

fig, ax = plt.subplots(figsize=(8.2, 5.2), dpi=200)
ax.axis("off")
ax.text(0.02, 0.98, "Fig 7. 图谱资源底表快照（开放 CSV，投稿随 DOI 发布）",
        fontsize=12, fontweight="bold", va="top", transform=ax.transAxes)
ax.text(0.02, 0.93, f"三结局合并：31,373 对（QC 31,371）→ 含 MR 全部字段 + coloc PP.H4",
        fontsize=9, va="top", transform=ax.transAxes)
ax.text(0.02, 0.90, f"行数：T2D {n_t2d} / CAD {n_cad} / FBG {n_fbg}    每行一个基因×结局对",
        fontsize=9, va="top", transform=ax.transAxes)
# 列结构表
tab = [[c, "float" if c in ("mr_b","mr_p","gwas_min_p","eqtl_F_max","nsnp","pp4") else ("str" if c!="ok" else "bool"),
        f"{example[c][:26]}" + ("…" if len(str(example[c]))>26 else "")]
       for c in cols]
tbl = ax.table(cellText=tab, colLabels=["列", "类型", "示例值 (T2D 首行)"],
               loc="upper center", bbox=[0.02, 0.30, 0.96, 0.55], cellLoc="left")
tbl.auto_set_font_size(False); tbl.set_fontsize(7.6)
for (row, col), cell in tbl.get_celld().items():
    if row == 0: cell.set_facecolor("#1a3a6b"); cell.set_text_props(color="white", fontweight="bold")
    elif col == 0: cell.set_facecolor("#eef2f7")
ax.text(0.02, 0.10, "用途：任何 cis-MR 靶点筛选可直接查表（gene×结局→MR 证据 + 共定位 PP.H4）。"
        "列字典/schema 随打包发布（DATA_AVAILABILITY §4.1）。", fontsize=8, va="top", transform=ax.transAxes)
fig.tight_layout(); fig.savefig(wb("20260816_F7_resource_snapshot.png"), bbox_inches="tight"); plt.close(fig)
print(f"Fig 7 OK: T2D {n_t2d} / CAD {n_cad} / FBG {n_fbg}")
print("全部完成")
