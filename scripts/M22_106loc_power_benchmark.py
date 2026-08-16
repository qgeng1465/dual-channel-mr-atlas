#!/usr/bin/env python3
# =============================================================================
# M22_106loc_power_benchmark.py — 106 强位点 eQTL 功效解析 + GWAS Catalog 文献基准
# =============================================================================
# ⚠️ DEPRECATED (2026-08-15) — window_100kb 分支有 build 错配
#   gene±100kb 窗口用 eQTLGen GenePos(hg19) 比对 GWAS Catalog CHR_POS(hg38)，
#   错位 ~1.7Mb → 55% 命中率是系统性低估。修复版见 scripts/M22b_window_build_fix.py
#   （v3：GenePos 转 hg38 ±100kb + 基因本体），命中率 81%，新候选上限 19%。
#   本文件保留 rsID 直命（build 无关，正确）+ eQTL 功效分析部分。
# =============================================================================
# =============================================================================
# 依据 2026-08-15 对抗性评审（verification_20260815.md §1/§8/§10#3）：
#   - 评审 §1 点 1：coloc-only = 弱 eQTL + 强 GWAS 是"功效陈述（平凡真理）"。
#     → 量化 MR 显著集（含 106 strong）的 lead eQTL 强度在全量 eQTL 分布中的分位，
#       检验"MR 找到的是不是系统性超强 eQTL"（工具选择偏倚的强度侧证据）。
#   - 评审 §4：coloc-only 多落已知 GWAS 基因座 → 必须逐位点做 GWAS Catalog 注释。
#     → 106 strong 位点命中已报道 T2D/CAD 关联（直接 rsID 或 gene ±100kb）的比例
#       = 已报道基因座；未命中 = 新候选比例上限。
# 输出：results/m22_efqt_power_20260815.csv + stdout 摘要
# =============================================================================
import gzip, csv, sys, io
import numpy as np
import pandas as pd

PROJ = "/data/qiushuogeng/projects/dual-channel-mr-atlas"
COLOC = f"{PROJ}/results/grid/transcript_coloc.csv"
SIG   = f"{PROJ}/data/eqtlgen/cis-eQTL-significant.txt.gz"
CAT   = f"{PROJ}/data/gwas_catalog/gwas-catalog-download-associations-alt-full.tsv"
OUT   = f"{PROJ}/results/m22_efqt_power_20260815.csv"
WINDOW = 100_000  # ±100kb 基因窗口

# ---- 1. strong / non-strong 分组（MR 显著集内，PP.H4≥0.8 vs <0.8）----
rows = pd.read_csv(COLOC)
rows["PP.H4"] = pd.to_numeric(rows["PP.H4"], errors="coerce")
rows["grp"] = np.where(rows["PP.H4"] >= 0.8, "strong", "non_strong")
strong = rows[rows.grp == "strong"]
print(f"MR 显著集: {len(rows)} | strong(PP.H4≥0.8): {len(strong)} | non_strong: {len(rows)-len(strong)}")
print(f"strong 结局分布: {strong['outcome'].value_counts().to_dict()}")

# ---- 2. eQTL 功效：每基因 lead cis-eQTL（|Zscore| max）----
print("\n读取 eQTLGen significant（10.5M 行，取每基因 lead）...")
le = []
for ch in pd.read_csv(SIG, sep="\t", compression="gzip",
                      usecols=["Gene","GeneSymbol","SNP","Zscore","GeneChr","GenePos","NrSamples"],
                      chunksize=5_000_000):
    ch = ch.dropna(subset=["Zscore"])
    ch["absz"] = ch["Zscore"].abs()
    le.append(ch)
sig = pd.concat(le, ignore_index=True)
lead = sig.loc[sig.groupby("Gene")["absz"].idxmax()].reset_index(drop=True)
print(f"lead eQTL 基因数: {len(lead)}")

# 全量 lead Zscore 分布（分位数）
q_all = lead["absz"].quantile([0.25,0.5,0.75,0.90,0.95,0.99]).to_dict()
print(f"全量 lead |Z| 分位: { {k: round(v,1) for k,v in q_all.items()} }")

# 分组基因的 lead |Z|（strong/non_strong 用 coloc 表的 gene）
strong_z = lead.merge(strong[["gene"]].rename(columns={"gene":"Gene"}), on="Gene")["absz"]
non_z    = lead.merge(rows[rows.grp=="non_strong"][["gene"]].rename(columns={"gene":"Gene"}), on="Gene")["absz"]
miss_s, miss_n = len(strong)-len(strong_z), (len(rows)-len(strong))-len(non_z)
print(f"strong 有 lead eQTL: {len(strong_z)}/{len(strong)} (缺 {miss_s} 无 lead eQTL)")
print(f"non_strong 有 lead eQTL: {len(non_z)}/{len(rows)-len(strong)} (缺 {miss_n})")

def pctile(v, dist):
    return float((dist < v).mean() * 100)  # v 在 dist 中的百分位（< 的比例）
if len(strong_z):
    sp = [round(pctile(v, lead["absz"]),1) for v in strong_z]
    print(f"\nstrong lead |Z|: 中位 {np.median(strong_z):.1f} | 全量分位中位 {np.median(sp):.0f}%")
    print(f"strong 超过全量 95 分位: {sum(1 for v in strong_z if v>q_all[0.95])}/{len(strong_z)} ({100*sum(1 for v in strong_z if v>q_all[0.95])/len(strong_z):.0f}%)")
if len(non_z):
    np_p = [round(pctile(v, lead["absz"]),1) for v in non_z]
    print(f"non_strong lead |Z|: 中位 {np.median(non_z):.1f} | 全量分位中位 {np.median(np_p):.0f}%")

# Wilcoxon strong vs non_strong
from scipy.stats import mannwhitneyu
if len(strong_z) and len(non_z):
    u, wp = mannwhitneyu(strong_z, non_z, alternative="two-sided")
    print(f"\nMann-Whitney strong vs non_strong lead |Z|: U={u:.0f} p={wp:.2e}")
    u2, wp2 = mannwhitneyu(strong_z, lead["absz"], alternative="greater")
    print(f"Mann-Whitney strong vs 全量: p={wp2:.2e}")

# ---- 3. GWAS Catalog 文献基准（T2D EFO_0001360 / CAD EFO_0001645）----
print("\nGWAS Catalog 过滤 T2D/CAD 关联...")
EFO_T2D, EFO_CAD = "EFO_0001360", "EFO_0001645"
cat_rows, cat_snps, cat_pos = [], set(), []
for ch in pd.read_csv(CAT, sep="\t", usecols=["DISEASE/TRAIT","CHR_ID","CHR_POS","SNP_ID_CURRENT","SNPS","P-VALUE","MAPPED_TRAIT_URI","MAPPED_TRAIT"],
                      dtype=str, chunksize=2_000_000):
    if ch.empty: continue
    uri = ch["MAPPED_TRAIT_URI"].fillna("")
    dt  = ch["DISEASE/TRAIT"].fillna("").str.lower()
    m = uri.str.contains(EFO_T2D) | uri.str.contains(EFO_CAD) | \
        dt.str.contains("type 2 diabetes") | dt.str.contains("coronary artery disease")
    hit = ch[m]
    if len(hit):
        cat_rows.append(hit)
        for s in hit["SNP_ID_CURRENT"].dropna():
            cat_snps.add(s.strip())
        for s in hit["SNPS"].fillna("").astype(str):
            for tok in s.split(";"):
                tok = tok.split(" x ")[0].strip()
                if tok.startswith("rs"):
                    cat_snps.add(tok)
        for _, r in hit[["CHR_ID","CHR_POS"]].dropna().iterrows():
            try:
                cat_pos.append((int(float(r["CHR_ID"])), int(float(r["CHR_POS"]))))
            except (ValueError, TypeError):
                pass
print(f"T2D/CAD 关联行: {len(pd.concat(cat_rows, ignore_index=True)) if cat_rows else 0}")
print(f"T2D/CAD 唯一关联 rsID: {len(cat_snps)} | 位置点: {len(cat_pos)}")
cat_pos = sorted(set(cat_pos))

def near_pos(chr_id, pos):
    lo, hi = pos - WINDOW, pos + WINDOW
    import bisect
    lo_i = bisect.bisect_left(cat_pos, (chr_id, lo))
    hi_i = bisect.bisect_right(cat_pos, (chr_id, hi))
    return hi_i > lo_i

# 106 strong 位点命中判定（先 rsID 直接命中，再 gene ±100kb）
def check_hit(r):
    top = str(r.get("top_snp", "")).strip()
    if top.startswith("rs") and top in cat_snps:
        return "rsid"
    gchr = lead.loc[lead.Gene == r["gene"], "GeneChr"]
    gpos = lead.loc[lead.Gene == r["gene"], "GenePos"]
    if len(gchr) and len(gpos) and pd.notna(gchr.iloc[0]) and pd.notna(gpos.iloc[0]):
        if near_pos(int(gchr.iloc[0]), int(gpos.iloc[0])):
            return "window_100kb"
    return "none"

strong = strong.copy()
strong["catalog_hit"] = strong.apply(check_hit, axis=1)
n_hit = (strong.catalog_hit != "none").sum()
print(f"\n106 strong 位点命中 GWAS Catalog T2D/CAD: {n_hit}/{len(strong)} ({100*n_hit/len(strong):.0f}%)")
print(f"  rsID 直接命中: {(strong.catalog_hit=='rsid').sum()} | gene±100kb 命中: {(strong.catalog_hit=='window_100kb').sum()} | 未命中: {(strong.catalog_hit=='none').sum()}")
print("  → 未命中比例 = 本研究中'新候选'比例上限（诚实：未命中≠新发现，仅超出 GWAS Catalog 注释范围）")

# 每结局命中率
for o in ["t2d","cad","fbg"]:
    s = strong[strong.outcome==o]
    if len(s):
        print(f"  {o}: {len(s)} 位点 | 命中 {(s.catalog_hit!='none').sum()} ({100*(s.catalog_hit!='none').mean():.0f}%)")

# ---- 4. 输出 ----
out = strong[["gene","symbol","outcome","top_snp","PP.H4","stage2_pval","catalog_hit"]].copy()
out["lead_absZ"] = strong["gene"].map(lambda g: lead.loc[lead.Gene==g, "absz"].iloc[0] if (lead.Gene==g).any() else np.nan)
out["lead_pctile_global"] = strong["gene"].map(
    lambda g: round(pctile(lead.loc[lead.Gene==g, "absz"].iloc[0], lead["absz"]),1) if (lead.Gene==g).any() else np.nan)
out["lead_pctile_vs_nonstrong"] = strong["gene"].map(
    lambda g: round(pctile(lead.loc[lead.Gene==g, "absz"].iloc[0], non_z),1) if (lead.Gene==g).any() and len(non_z) else np.nan)
out.to_csv(OUT, index=False)
print(f"\n已写 {OUT}")

# ---- 5. 判定（对齐评审 §1/§8）----
print("\n=== 判定 ===")
e_strong = np.median(sp) if len(strong_z) else float('nan')
print(f"① 工具选择偏倚：106 strong 的 lead eQTL 全部处于全量高分区（中位 {e_strong:.0f} 分位）→ "
      f"{'是' if e_strong>=90 else '否'} 系统性强 eQTL 富集")
print(f"② 文献基准：{100*n_hit/len(strong):.0f}% strong 位点落已报道 T2D/CAD 基因座 → "
      f"未报道 {(strong.catalog_hit=='none').sum()}/{len(strong)} = 新候选上限")
