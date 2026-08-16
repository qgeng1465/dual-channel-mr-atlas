#!/usr/bin/env python3
# =============================================================================
# M40a_cross_tissue_20260817.py — 跨组织 lead-eQTL 一致性（GTEx v8 8 组织）
# =============================================================================
# 诚实降级说明：GTEx v8 全 cis eQTL 文件（~13GB）当前网络不可下载，本地只有
#   per-gene lead-eQTL 文件（每基因 1 行）。Module 1 多组织维度以 15 候选在各
#   组织的 lead-eQTL 锚点呈现：每组织是否有显著 eQTL、lead 变体与 eQTLGen 全血
#   lead 是否同一变体（same_variant，GRCh38 chr_pos 精确匹配）或同区域
#   （within_50kb，同 chr ±50kb）。完整多组织 coloc 需全 cis SNP×基因矩阵
#   （如实报告不可得）。
# 输入：data/gtex/*.egenes.txt.gz；eQTLGen lead → GRCh38 映射（GWAS full 文件实测）
# 输出：results/m40a_cross_tissue_20260817.csv + 摘要
# =============================================================================
import pandas as pd, numpy as np, glob, os, re, json

REPO = "<repo-root>"
GTEX = os.path.join(REPO, "data", "gtex")

TISSUE_N = {"Liver": 208, "Pancreas": 328, "Whole_Blood": 755,
            "Adipose_Subcutaneous": 581, "Adipose_Visceral_Omentum": 469,
            "Muscle_Skeletal": 803, "Artery_Coronary": 213, "Artery_Aorta": 387}

# eQTLGen lead → GRCh38 (variant_id, pos)；13 个唯一 lead，实测自 GWAS full 文件
LEAD_B38 = {"rs10494191": ("1_117001393_A_C", 117001393), "rs4851758": ("2_105335960_C_T", 105335960),
            "rs10049087": ("3_50085984_A_G", 50085984), "rs1004877": ("3_45025925_G_A", 45025925),
            "rs12332382": ("5_244828_C_T", 244828), "rs10430625": ("10_100136849_G_A", 100136849),
            "rs12780155": ("10_12219853_T_A", 12219853), "rs2789422": ("1_159922298_G_A", 159922298),
            "rs11191447": ("10_102892566_C_T", 102892566), "rs6598075": ("11_207275_C_G", 207275),
            "rs1123462": ("13_32571325_C_T", 32571325), "rs56228609": ("16_56953853_C_T", 56953853),
            "rs4760": ("19_43648948_A_G", 43648948)}

cand = pd.read_csv(f"{REPO}/results/candidate15_replication_20260816.csv")
grid = pd.read_csv(f"{REPO}/results/grid/transcript_grid_mr.csv")
grid["lead"] = grid["note"].str.extract(r"lead (rs[0-9]+)")
grid = grid[grid.symbol.isin(cand.symbol)]
gmap = {}
for _, r in grid.iterrows():
    if pd.isna(r.lead) or r.outcome not in cand[cand.symbol == r.symbol].outcome.values:
        continue
    gmap[r.symbol] = r.lead
# 同符号多 outcome 取首个
gmap = {s: v for s, v in gmap.items()}

rows = []
for f in sorted(glob.glob(f"{GTEX}/*.egenes.txt.gz")):
    tissue = os.path.basename(f).replace(".egenes.txt.gz", "").replace(".v8", "").replace("GTEx_Analysis_v8_", "")
    n = TISSUE_N.get(tissue, np.nan)
    d = pd.read_csv(f, sep="\t", compression="gzip", usecols=[
        "gene_id", "gene_name", "variant_id", "slope", "slope_se", "pval_nominal", "qval"])
    d["gene_core"] = d["gene_id"].str.split(".").str[0]
    d = d[d["gene_core"].isin(cand.gene.values)]
    if d.empty:
        continue
    d = d.drop_duplicates(subset="gene_core")
    d = d.set_index("gene_core")
    for _, c in cand.iterrows():
        gene, sym, oc = c.gene, c.symbol, c.outcome
        if gene not in d.index:
            continue
        vid = d.loc[gene, "variant_id"]
        m = re.match(r"chr(\d+)_(\d+)_([A-Z]+)_([A-Z]+)_b38", str(vid))
        if not m:
            continue
        chr_b38, pos_b38 = "chr" + m.group(1), int(m.group(2))
        eqtlgen_lead = gmap.get(sym)
        vid_e, pos_e = LEAD_B38.get(eqtlgen_lead, (None, None)) if eqtlgen_lead else (None, None)
        same_variant = (vid == vid_e) if vid_e else False
        within50 = (pos_e is not None and chr_b38 == "chr" + str(vid_e).split("_")[0] and
                    abs(pos_b38 - pos_e) <= 50000) if vid_e else False
        rows.append(dict(symbol=sym, gene=gene, outcome=oc, tissue=tissue, n_tissue=n,
                         lead_variant=vid, chr=chr_b38, pos=pos_b38,
                         slope=d.loc[gene, "slope"], slope_se=d.loc[gene, "slope_se"],
                         pval=d.loc[gene, "pval_nominal"], qval=d.loc[gene, "qval"],
                         eqtlgen_lead=eqtlgen_lead or "", eqtlgen_variant=vid_e or "",
                         same_variant=same_variant, within_50kb=within50))

df = pd.DataFrame(rows)
df.to_csv(f"{REPO}/results/m40a_cross_tissue_20260817.csv", index=False)

print(f"{'candidate':10s} {'tissues':6s} {'sig':4s} {'sameVar':7s} {'within50':8s}  notes")
for sym in sorted(set(df.symbol)):
    sub = df[df.symbol == sym]
    sig = sub[sub.qval < 0.05]
    same = sig[sig.same_variant]
    w50 = sig[sig.within_50kb]
    print(f"{sym:10s} {len(sub):6d} {len(sig):4d} {len(same):7d} {len(w50):8d}  "
          f"same_tissue={','.join(same.tissue) if len(same) else '-'}  w50_tissue={','.join(w50.tissue) if len(w50) else '-'}")
print("\nrows:", len(df), "| sig rows:", int(df.qval.lt(0.05).sum()))
