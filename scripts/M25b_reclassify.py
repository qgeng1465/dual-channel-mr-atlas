#!/usr/bin/env python3
# =============================================================================
# M25b_reclassify.py — 23 new strong 的 catalog 重分类（宽口径 CAD/T2D 过滤）
# =============================================================================
# 动机（诚实性修正）：M25/M22b 的 catalog 过滤只匹配 "coronary artery disease"，
#   漏掉 catalog 里 "coronary heart disease / myocardial infarction / ischaemic
#   heart disease / atherosclerosis" 等表述 → 低估"已知位点"，把该算已知的误判为
#   新位点。本脚本用宽口径 trait 字符串 + EFO/MONDO URI 重读 catalog 重分类，
#   得到保守（不夸大新奇度）的分类。
# 输入：results/m25_new_strong_annotation_20260816.csv（含每基因 hg38 坐标）
# 输出：results/m25b_reclassify_20260816.csv
# =============================================================================
import gzip, bisect, csv
import pandas as pd

PROJ = "<repo-root>"
RES  = f"{PROJ}/results"
CAT  = f"{PROJ}/data/gwas_catalog/gwas-catalog-download-associations-alt-full.tsv"
OUT  = f"{RES}/m25b_reclassify_20260816.csv"
W100 = 100_000
W250 = 250_000

# 宽口径 T2D/CAD 关键词（DISEASE/TRAIT 小写）
T2D_KW = ["type 2 diabetes", "type ii diabetes", "non-insulin dependent",
          "niddm", "maturity onset diabetes", "diabetes mellitus type 2", "t2d"]
CAD_KW = ["coronary", "myocardial", "ischaemic heart", "ischemic heart",
          "coronary heart", "atheroscler", "angina", "infarction"]
# URI（EFO 标准 + MONDO 常见 T2D/CAD/MI）
URI_KW = ["efo_0001360", "efo_0001645", "mondo_0005147", "mondo_0005148",
          "mondo_0004995", "mondo_0005267", "mondo_0005213", "mondo_0023450",
          "mondo_0008577", "mondo_0001360"]

# ---- 1. 读入 M25 坐标 ----
m25 = pd.read_csv(f"{RES}/m25_new_strong_annotation_20260816.csv", dtype=str).fillna("")
print(f"M25 基因: {len(m25)}")

# ---- 2. 宽口径读 catalog ----
cat_pos, cat_genes_t2dcad = set(), set()
for ch in pd.read_csv(CAT, sep="\t", dtype=str,
                      usecols=["CHR_ID","CHR_POS","MAPPED_TRAIT_URI","DISEASE/TRAIT",
                               "MAPPED_GENE","REPORTED GENE(S)"],
                      chunksize=2_000_000):
    if ch.empty: continue
    uri = ch["MAPPED_TRAIT_URI"].fillna("").str.lower()
    dt  = ch["DISEASE/TRAIT"].fillna("").str.lower()
    is_t2dcad = uri.str.contains("|".join(URI_KW)) | \
                dt.str.contains("|".join(T2D_KW)) | \
                dt.str.contains("|".join(CAD_KW))
    hit = ch[is_t2dcad]
    if len(hit):
        for _, r in hit[["CHR_ID","CHR_POS"]].dropna().iterrows():
            try: cat_pos.add((int(float(r["CHR_ID"])), int(float(r["CHR_POS"]))))
            except (ValueError, TypeError): pass
        for col in ["MAPPED_GENE","REPORTED GENE(S)"]:
            for s in hit[col].fillna("").astype(str):
                for tok in s.replace(" - "," ").split():
                    tok = tok.strip()
                    if tok.startswith("ENSG") or (tok.isalpha() and len(tok)>=2):
                        cat_genes_t2dcad.add(tok)
cat_pos = sorted(cat_pos)
print(f"宽口径 catalog T2D/CAD: 位置 {len(cat_pos)} | 基因本体 {len(cat_genes_t2dcad)}")

def near(chr_id, pos, w):
    return bisect.bisect_right(cat_pos, (chr_id, pos+w)) > bisect.bisect_left(cat_pos, (chr_id, pos-w))
def nearest(chr_id, pos, w):
    a = bisect.bisect_right(cat_pos, (chr_id, pos-w)); b = bisect.bisect_right(cat_pos, (chr_id, pos+w))
    best = None
    for i in range(a, b):
        d = abs(cat_pos[i][1]-pos)
        if best is None or d < best: best = d
    return best

# ---- 3. 重分类 ----
out = []
for _, r in m25.iterrows():
    sym = r["symbol"]; g = r["gene"]
    if r["pos_hg38"] == "":
        tier = "no_gene_pos"; nd = None; in100 = in250 = False; gene_anno = False
    else:
        chr38, pos38 = int(float(r["chr_hg38"])), int(r["pos_hg38"])
        gene_anno = sym in cat_genes_t2dcad or g in cat_genes_t2dcad
        in100 = near(chr38, pos38, W100); in250 = near(chr38, pos38, W250)
        nd = nearest(chr38, pos38, W250)
        if gene_anno: tier = "catalog_gene_t2dcad"
        elif in100:   tier = "known_locus_100kb_novel_gene"
        elif in250:   tier = "known_locus_250kb_novel_gene"
        else:         tier = "no_catalog_t2dcad"
    out.append({**{k: r[k] for k in r.keys()}, "gene_anno_broad": gene_anno,
                "cat_pos_100kb_broad": in100, "cat_pos_250kb_broad": in250,
                "nearest_cat_kb_broad": round(nd/1000,1) if nd is not None else "",
                "tier_broad": tier})

df = pd.DataFrame(out)
df.to_csv(OUT, index=False)
print(f"\n=== 宽口径重分类汇总 ===")
print(df["tier_broad"].value_counts().to_string())
print(f"\n=== 新旧 tier 对照 ===")
df2 = df[["symbol","outcome","pp4","gwas_min_p","pos_hg38","tier","tier_broad",
          "cat_pos_100kb","cat_pos_100kb_broad","cat_pos_250kb","cat_pos_250kb_broad"]]
for _, r in df2.iterrows():
    mark = "◀改变" if r["tier"] != r["tier_broad"] else ""
    print(f"{r['symbol']:<10}{r['outcome']:<5} {r['tier']:<28}→ {r['tier_broad']:<28} {mark}")
print(f"\n已写 {OUT}")
