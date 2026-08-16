#!/usr/bin/env python3
# =============================================================================
# M22b v3 — M22 window_100kb 的 build 错配复核（审核 agent 发现，2026-08-15）
# -----------------------------------------------------------------------------
# v1 缺陷：用"基因 cis SNP 集合(hg19 GenePos±100kb)的 hg38"做锚，但 cis SNP
#   可散布在基因±1Mb → 有效窗口被放大（C15orf62 被误命中）。
# v2 缺陷：用 top_snp 单一锚，top_snp 可距基因本体数百 kb（C15orf62 top_snp
#   距基因 196kb）→ 与"基因±100kb"语义不符。
# v3 方案（正确）：用该基因近端 cis SNP 的 hg19→hg38 偏移量中位数（同区域
#   build 平移稳定，实测 std<1bp）把 GenePos 精确转换到 hg38，再 ±100kb 判定。
#   补充 catalog MAPPED_GENE/REPORTED_GENE 本体校验（完全 build 无关）。
# 判定规则：
#   rsid          : top_snp ∈ catalog T2D/CAD rsID 集合
#   gene_100kb    : GenePos(hg38) ±100kb 内存在 catalog T2D/CAD 关联位置
#   gene_anno     : 基因本体(ENSG/symbol) 出现在 catalog T2D/CAD 的 MAPPED_GENE/REPORTED_GENE
#   none          : 均未命中
# 输出：results/m22b_window_fix_20260815.csv + stdout
# =============================================================================
import gzip, bisect, json, os
import pandas as pd
import numpy as np

PROJ = "<repo-root>"
COLOC = f"{PROJ}/results/grid/transcript_coloc.csv"
SIG   = f"{PROJ}/data/eqtlgen/cis-eQTL-significant.txt.gz"
CAT   = f"{PROJ}/data/gwas_catalog/gwas-catalog-download-associations-alt-full.tsv"
def first_existing(*paths):
    for p in paths:
        if os.path.exists(p):
            return p
    return paths[0]

OLD   = first_existing(f"{PROJ}/results/m22_efqt_power_20260815.csv",
                       f"{PROJ}/results/archive/m22_efqt_power_20260815.csv")
OUT   = f"{PROJ}/results/m22b_window_fix_20260815.csv"
WINDOW = 100_000

# ---- 1. 106 strong 位点 + 旧分类 ----
rows = pd.read_csv(COLOC)
rows["PP.H4"] = pd.to_numeric(rows["PP.H4"], errors="coerce")
strong = rows[rows["PP.H4"] >= 0.8].copy()
old = pd.read_csv(OLD)[["gene","catalog_hit"]].rename(columns={"catalog_hit":"old_hit"})
strong = strong.merge(old, on="gene", how="left")
print(f"strong 位点: {len(strong)} | 结局: {strong['outcome'].value_counts().to_dict()}")

# ---- 2. catalog T2D/CAD：rsID 集合 + hg38 位置 + 基因本体 ----
print("读取 GWAS Catalog T2D/CAD ...")
cat_snps, cat_pos, cat_genes = set(), set(), set()
for ch in pd.read_csv(CAT, sep="\t",
                      usecols=["SNP_ID_CURRENT","SNPS","CHR_ID","CHR_POS",
                               "MAPPED_TRAIT_URI","DISEASE/TRAIT","MAPPED_GENE","REPORTED GENE(S)"],
                      dtype=str, chunksize=2_000_000):
    if ch.empty: continue
    uri = ch["MAPPED_TRAIT_URI"].fillna("")
    dt  = ch["DISEASE/TRAIT"].fillna("").str.lower()
    m = uri.str.contains("EFO_0001360") | uri.str.contains("EFO_0001645") | \
        dt.str.contains("type 2 diabetes") | dt.str.contains("coronary artery disease")
    hit = ch[m]
    if len(hit):
        for s in hit["SNP_ID_CURRENT"].dropna():
            cat_snps.add(s.strip())
        for s in hit["SNPS"].fillna("").astype(str):
            for tok in s.split(";"):
                tok = tok.split(" x ")[0].strip()
                if tok.startswith("rs"): cat_snps.add(tok)
        for _, r in hit[["CHR_ID","CHR_POS"]].dropna().iterrows():
            try: cat_pos.add((int(float(r["CHR_ID"])), int(float(r["CHR_POS"]))))
            except (ValueError, TypeError): pass
        for col in ["MAPPED_GENE","REPORTED GENE(S)"]:
            for s in hit[col].fillna("").astype(str):
                for tok in s.replace(" - ", " ").split():
                    tok = tok.strip()
                    if tok.startswith("ENSG") or (tok.isalpha() and len(tok) >= 2):
                        cat_genes.add(tok)
cat_pos = sorted(cat_pos)
print(f"catalog T2D/CAD: rsID {len(cat_snps)} | 位置 {len(cat_pos)} | 基因本体 {len(cat_genes)}")

def near_pos(chr_id, pos):
    lo, hi = pos - WINDOW, pos + WINDOW
    lo_i = bisect.bisect_left(cat_pos, (chr_id, lo))
    hi_i = bisect.bisect_right(cat_pos, (chr_id, hi))
    return hi_i > lo_i

# ---- 3. 每基因近端 cis SNP 的 (SNPPos_hg19, rsid) 用于 offset 计算 ----
print("读取 eQTLGen significant，收集基因近端 cis SNP ...")
genes = set(strong["gene"])
gene_info = {}   # gene -> (GeneChr, GenePos_hg19)
gene_cis = {g: [] for g in genes}  # gene -> [(SNPPos_hg19, rsid)]
for ch in pd.read_csv(SIG, sep="\t", compression="gzip",
                      usecols=["Gene","GeneSymbol","GeneChr","GenePos","SNP","SNPPos"],
                      chunksize=5_000_000):
    m = ch["Gene"].isin(genes)
    if m.sum() == 0: continue
    sub = ch[m]
    for g, gc, gp, snp, sp in sub[["Gene","GeneChr","GenePos","SNP","SNPPos"]].itertuples(index=False):
        if g not in gene_info:
            gene_info[g] = (int(gc), int(gp))
        if pd.notna(sp):
            gene_cis[g].append((int(sp), snp))
print(f"基因数: {len(gene_info)}")

# ---- 4. OpenGWAS full → rsid→hg38 ----
need = set().union(*(set(s for _, s in v) for v in gene_cis.values()))
print(f"需查 hg38 的 cis SNP: {len(need)}")
rsid_hg38 = {}
for fn in ["t2d","cad","fbg"]:
    fp = f"{PROJ}/data/opengwas/full/{fn}_full.gz"
    with gzip.open(fp, "rt") as fh:
        for i, line in enumerate(fh):
            if i == 0: continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 10: continue
            r = p[1]
            if r in need and r not in rsid_hg38:
                try: rsid_hg38[r] = (int(p[2]), int(p[3]))
                except (ValueError, TypeError): pass
            if len(rsid_hg38) == len(need): break
print(f"cis SNP 命中 hg38: {len(rsid_hg38)}/{len(need)}")

# ---- 5. 每基因：GenePos(hg38) = GenePos(hg19) + 中位 offset ----
def gene_hg38_pos(g):
    chr_hg19, pos_hg19 = gene_info[g]
    offsets = []
    for sp19, snp in gene_cis[g]:
        hg = rsid_hg38.get(snp)
        if hg and hg[0] == chr_hg19:
            offsets.append(hg[1] - sp19)
    if not offsets:
        return chr_hg19, None, None, 0
    med = float(np.median(offsets))
    return chr_hg19, int(pos_hg19 + med), int(med), len(offsets)

# ---- 6. 重判 ----
print("重判 106 位点 ...")
new_hit, hit_note, gene_offsets = {}, {}, {}
for _, r in strong.iterrows():
    g, top = r["gene"], str(r["top_snp"]).strip()
    chrg, gpos_hg38, off, n = gene_hg38_pos(g)
    gene_offsets[g] = (gpos_hg38, off, n)
    if top.startswith("rs") and top in cat_snps:
        new_hit[g] = "rsid"; hit_note[g] = top
    elif gpos_hg38 and near_pos(chrg, gpos_hg38):
        new_hit[g] = "gene_100kb"; hit_note[g] = f"GenePos(hg38)@chr{chrg}:{gpos_hg38}"
    elif g in cat_genes or (r.get("symbol") and str(r.get("symbol")) in cat_genes):
        new_hit[g] = "gene_anno"; hit_note[g] = f"{r.get('symbol')} in catalog MAPPED_GENE"
    else:
        new_hit[g] = "none"; hit_note[g] = f"GenePos(hg38)@chr{chrg}:{gpos_hg38}"

strong["new_hit"] = strong["gene"].map(new_hit)
strong["new_hit_note"] = strong["gene"].map(hit_note)
strong["gene_hg38_pos"] = strong["gene"].map(lambda g: gene_offsets[g][0])
strong["offset"] = strong["gene"].map(lambda g: gene_offsets[g][1])

# ---- 7. 汇总 ----
print("\n=== 旧分类（M22，hg19 window 错配）===")
print(strong["old_hit"].fillna("(无)").value_counts().to_dict())
print("=== 新分类（v3，GenePos(hg38)±100kb + 基因本体）===")
print(strong["new_hit"].value_counts().to_dict())
hit = strong["new_hit"].isin(["rsid","gene_100kb","gene_anno"])
n_new = int(hit.sum())
print(f"\n命中率: 旧 {(strong['old_hit'].fillna('')!='none').sum()}/{len(strong)} → 新 {n_new}/{len(strong)} ({100*n_new/len(strong):.0f}%)")
print(f"新候选上限(none): {(~hit).sum()}/{len(strong)} ({100*(~hit).mean():.0f}%)")
for o in ["t2d","cad","fbg"]:
    s = strong[strong.outcome==o]
    if len(s):
        print(f"  {o}: {len(s)} 位点 | 命中 {100*s.new_hit.isin(['rsid','gene_100kb','gene_anno']).mean():.0f}% | 新候选上限 {100*(~s.new_hit.isin(['rsid','gene_100kb','gene_anno'])).mean():.0f}%")

print("\n=== 逐位点对比（旧→新变化）===")
chg = strong[strong["old_hit"].fillna("") != strong["new_hit"]]
for _, r in chg.iterrows():
    print(f"  {r['symbol']:>14} {r['outcome']:<4} {str(r['old_hit']):>13} → {str(r['new_hit']):>11}  {r['new_hit_note']}")

# ---- 8. 输出 ----
out = strong[["gene","symbol","outcome","top_snp","PP.H4","old_hit","new_hit","new_hit_note","gene_hg38_pos","offset"]]
out.to_csv(OUT, index=False)
print(f"\n已写 {OUT}")
