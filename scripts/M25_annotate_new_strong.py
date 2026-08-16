#!/usr/bin/env python3
# =============================================================================
# M25_annotate_new_strong.py — 23 个 sig 内新增 strong coloc 的 GWAS Catalog 注释
# =============================================================================
# 目的（路径 B：把负结果翻成正向基因发现）：
#   M24 发现全量扫描在 MR 显著集内多出 23 个 strong coloc（不在 M5 已知 106）。
#   本脚本为这 23 个基因做注释，判定：
#     - 该基因是否已在 GWAS Catalog T2D/CAD 被报道（MAPPED_GENE/REPORTED_GENE 本体）
#     - 该基因所在位点（GenePos hg38 ±100kb / ±250kb）是否有 catalog T2D/CAD 关联
#       → "已知风险位点内的新效应基因"（novel effector gene）或"全新位点"
#     - 基因本体在 catalog 任意表型被报道过与否（已知度）
#   build 一致：GenePos(hg19)→hg38 用近端 cis-SNP 中位偏移（M22b v3 同法）
# 输出：results/m25_new_strong_annotation_20260816.csv + stdout 摘要
# =============================================================================
import gzip, bisect, sys, csv
import pandas as pd
import numpy as np

PROJ  = "<repo-root>"
RES   = f"{PROJ}/results"
SIG   = f"{PROJ}/data/eqtlgen/cis-eQTL-significant.txt.gz"
CAT   = f"{PROJ}/data/gwas_catalog/gwas-catalog-download-associations-alt-full.tsv"
COLOC = f"{RES}/grid/transcript_coloc.csv"
OUT   = f"{RES}/m25_new_strong_annotation_20260816.csv"

def f(x, d=float('nan')):
    try: return float(x)
    except (TypeError, ValueError): return d

# ---- 1. 23 个 new strong（M24 口径：sig 内 pp4>=0.8 且不在 known 106）----
rows = []
for on in ["t2d","cad","fbg"]:
    with open(f"{RES}/coloc_full_{on}_20260815.csv") as fh:
        for r in csv if False else __import__('csv').DictReader(fh):
            r["outcome"] = on
            rows.append(r)
known = set()
with open(f"{RES}/grid/transcript_coloc_hits.csv") as fh:
    for r in __import__('csv').DictReader(fh):
        known.add((r["gene"], r["outcome"]))
new23 = [r for r in rows if r["ok"].strip()=="TRUE" and f(r["pp4"])>=0.8
         and f(r["mr_p"])<0.05 and (r["gene"], r["outcome"]) not in known]
print(f"new strong: {len(new23)}")
genes = {r["gene"] for r in new23}

# top_snp（若有）：从 M5 coloc 表取同 gene-outcome 的 top_snp
top_of = {}
try:
    for r in __import__('csv').DictReader(open(COLOC)):
        top_of[(r["gene"], r["outcome"])] = r.get("top_snp","")
except FileNotFoundError:
    pass

# ---- 2. eQTLGen：基因位置 + 近端 cis SNP（同 M22b）----
gene_info = {}
gene_cis  = {g: [] for g in genes}
for ch in pd.read_csv(SIG, sep="\t", compression="gzip",
                      usecols=["Gene","GeneChr","GenePos","SNP","SNPPos"],
                      chunksize=5_000_000):
    m = ch["Gene"].isin(genes)
    if m.sum() == 0: continue
    for g, gc, gp, snp, sp in ch[m][["Gene","GeneChr","GenePos","SNP","SNPPos"]].itertuples(index=False):
        if g not in gene_info:
            gene_info[g] = (int(gc), int(gp))
        if pd.notna(sp):
            gene_cis[g].append((int(sp), snp))
print(f"基因定位: {len(gene_info)}/{len(genes)}")

# ---- 3. OpenGWAS full → rsid→hg38（三结局文件，命中即收）----
need = set().union(*(set(s for _, s in v) for v in gene_cis.values()))
print(f"需查 hg38: {len(need)} cis SNP")
rsid_hg38 = {}
for fn in ["t2d","cad","fbg"]:
    with gzip.open(f"{PROJ}/data/opengwas/full/{fn}_full.gz","rt") as fh:
        for i, line in enumerate(fh):
            if i == 0: continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 10: continue
            r = p[1]
            if r in need and r not in rsid_hg38:
                try: rsid_hg38[r] = (int(p[2]), int(p[3]))
                except (ValueError, TypeError): pass
            if len(rsid_hg38) == len(need):
                print(f"  ({fn} 收集完毕)")
                break
print(f"cis SNP 命中 hg38: {len(rsid_hg38)}/{len(need)}")

def gene_hg38_pos(g):
    chr19, pos19 = gene_info[g]
    offs = [rsid_hg38[s][1]-sp for sp, s in gene_cis[g]
            if s in rsid_hg38 and rsid_hg38[s][0]==chr19]
    if not offs: return chr19, None, 0, len(gene_cis[g])
    return chr19, int(pos19 + np.median(offs)), int(np.median(offs)), len(gene_cis[g])

# ---- 4. GWAS Catalog T2D/CAD：位置集合 + 基因本体集合 + 任意表型基因本体 ----
cat_pos, cat_genes_t2dcad, cat_genes_any, cat_snps = set(), set(), set(), set()
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
        for _, r in hit[["CHR_ID","CHR_POS"]].dropna().iterrows():
            try: cat_pos.add((int(float(r["CHR_ID"])), int(float(r["CHR_POS"]))))
            except (ValueError, TypeError): pass
        for col in ["MAPPED_GENE","REPORTED GENE(S)"]:
            for s in hit[col].fillna("").astype(str):
                for tok in s.replace(" - "," ").split():
                    tok = tok.strip()
                    if tok.startswith("ENSG") or (tok.isalpha() and len(tok)>=2):
                        cat_genes_t2dcad.add(tok)
    # 任意表型的基因本体（已知度）
    for col in ["MAPPED_GENE","REPORTED GENE(S)"]:
        for s in ch[col].fillna("").astype(str):
            for tok in s.replace(" - "," ").split():
                tok = tok.strip()
                if tok.startswith("ENSG") or (tok.isalpha() and len(tok)>=2):
                    cat_genes_any.add(tok)
cat_pos = sorted(cat_pos)
print(f"catalog T2D/CAD: 位置 {len(cat_pos)} | 基因本体 {len(cat_genes_t2dcad)} | 任意表型基因 {len(cat_genes_any)}")

def near_pos(chr_id, pos, w):
    lo, hi = pos-w, pos+w
    return bisect.bisect_right(cat_pos, (chr_id, hi)) > bisect.bisect_left(cat_pos, (chr_id, lo))
def nearest_pos(chr_id, pos, w):
    lo, hi = pos-w, pos+w
    a = bisect.bisect_right(cat_pos, (chr_id, lo)); b = bisect.bisect_right(cat_pos, (chr_id, hi))
    best = None
    for i in range(a, b):
        d = abs(cat_pos[i][1]-pos)
        if best is None or d < best: best = d
    return best, (a,b)

# ---- 5. 判分类 + 注释 ----
out_rows = []
for r in sorted(new23, key=lambda x: -f(x["pp4"])):
    g, sym, on = r["gene"], r["symbol"], r["outcome"]
    if g not in gene_info:
        out_rows.append({**{k: r.get(k,"") for k in ["gene","symbol","outcome","mr_p","gwas_min_p","eqtl_F_max","nsnp","pp4"]},
                         "chr_hg19":"","pos_hg19":"","chr_hg38":"","pos_hg38":"","offset":"",
                         "catalog_gene_t2dcad":"?","catalog_gene_any":"?","cat_pos_100kb":"?","cat_pos_250kb":"?",
                         "nearest_cat_kb":"","tier":"no_gene_pos","top_snp":top_of.get((g,on),""),"n_cis_snp":""})
        continue
    chr19, pos19 = gene_info[g]
    chr38, pos38, off, ncis = gene_hg38_pos(g)
    gene_t2dcad = sym in cat_genes_t2dcad or g in cat_genes_t2dcad
    gene_any    = sym in cat_genes_any or g in cat_genes_any
    if pos38 is not None:
        in100 = near_pos(chr38, pos38, 100_000)
        in250 = near_pos(chr38, pos38, 250_000)
        nd, _ = nearest_pos(chr38, pos38, 250_000)
    else:
        in100 = in250 = False; nd = None
    # tier
    if gene_t2dcad: tier = "catalog_gene_t2dcad"          # 基因本体已被 T2D/CAD 报道
    elif in100:     tier = "known_locus_100kb_novel_gene"  # 已知位点±100kb，基因本体未报道
    elif in250:     tier = "known_locus_250kb_novel_gene"
    else:           tier = "no_catalog_t2dcad"             # 该区域无 catalog T2D/CAD
    out_rows.append({
        "gene": g, "symbol": sym, "outcome": on,
        "mr_p": r["mr_p"], "gwas_min_p": r["gwas_min_p"], "eqtl_F_max": r["eqtl_F_max"],
        "nsnp": r["nsnp"], "pp4": r["pp4"],
        "chr_hg19": chr19, "pos_hg19": pos19, "chr_hg38": chr38,
        "pos_hg38": pos38 if pos38 is not None else "", "offset": off,
        "n_cis_snp": ncis,
        "catalog_gene_t2dcad": gene_t2dcad, "catalog_gene_any": gene_any,
        "cat_pos_100kb": in100, "cat_pos_250kb": in250,
        "nearest_cat_kb": round(nd/1000,1) if nd is not None else "",
        "tier": tier, "top_snp": top_of.get((g,on),""),
    })

df = pd.DataFrame(out_rows)
df.to_csv(OUT, index=False)

print(f"\n=== 分类汇总 ===")
print(df["tier"].value_counts().to_string())
print(f"\n=== 明细 ===")
for _, r in df.iterrows():
    print(f"{r['symbol']:<12}{r['outcome']:<5} pp4={f(r['pp4']):.3f} gwas_min_p={f(r['gwas_min_p']):.1e} "
          f"chr{r['chr_hg38']}:{r['pos_hg38']} cat100={r['cat_pos_100kb']} cat250={r['cat_pos_250kb']} "
          f"gene_t2dcad={r['catalog_gene_t2dcad']} gene_any={r['catalog_gene_any']} → {r['tier']}")
print(f"\n已写 {OUT}")
