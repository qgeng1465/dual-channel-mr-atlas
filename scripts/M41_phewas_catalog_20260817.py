#!/usr/bin/env python3
# =============================================================================
# M41_phewas_catalog_20260817.py - 候选 lead SNP 的 GWAS Catalog PheWAS 多效性屏幕
# =============================================================================
# 目的：审稿人任务 4（系统性 PheWAS 多效性过滤）。对 15 个候选效应基因的
#   每个 lead cis-eQTL 工具 SNP，在 GWAS Catalog 全量关联表（738MB，P<5e-8
#   收录口径）中做全表型扫描，计算跨域多效性指标：
#   * n_traits / n_studies：该 SNP 在全目录中的关联性状/研究数
#   * n_cardio：心血管代谢域（糖尿病/血糖、CAD、血脂、BMI、血压）
#   * n_noncardio：跨域（非心肺代谢）关联数 —— 多效性顾虑的量化
#   * cross_domain：n_noncardio>0 → 该工具的排除性(exclusion restriction)存疑
# 诚实 caveat：GWAS Catalog 是"收录级"数据库（默认 P<5e-8），无统一效应量比较
#   与 LD 调整，性状间不独立；跨域关联计数为多效性的保守屏幕，非因果结论。
# 输入：
#   data/gwas_catalog/gwas-catalog-download-associations-alt-full.tsv  (738MB)
#   results/candidate15_replication_20260816.csv / results/grid/transcript_grid_mr.csv
# 输出：
#   results/m41_phewas_snp_20260817.csv       每 lead SNP 多效性汇总
#   results/m41_phewas_gene_20260817.csv      每候选基因多效性汇总（SNP→基因映射）
#   results/m41_phewas_hits_20260817.csv      每 SNP×性状 关联明细（volcano 数据）
#   results/figures/20260817_FigS_phewas_volcano.png
# =============================================================================
import os, re, subprocess
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = "/data/qiushuogeng/projects/dual-channel-mr-atlas"
CAT = os.path.join(REPO, "data/gwas_catalog/gwas-catalog-download-associations-alt-full.tsv")
FIG = os.path.join(REPO, "results", "figures")
os.makedirs(FIG, exist_ok=True)

# ---- 1. 候选 → lead SNP ----
cand = pd.read_csv(os.path.join(REPO, "results/candidate15_replication_20260816.csv"))
grid = pd.read_csv(os.path.join(REPO, "results/grid/transcript_grid_mr.csv"))
grid["lead"] = grid["note"].str.extract(r"lead (rs[0-9]+)")
sub = grid[grid["symbol"].isin(cand["symbol"])].dropna(subset=["lead"])
sub = sub.merge(cand[["symbol", "outcome", "mr_p", "mr_b", "is_new", "strong"]],
                on=["symbol", "outcome"], how="left")
snps = sorted(sub["lead"].unique())
print(f"candidate lead SNPs: {len(snps)} -> {snps}")

# ---- 2. 单趟 awk 提取 ----
alt_esc = "|".join(snps)
# 匹配列 21（STRONGEST SNP-RISK ALLELE, "rsX-?"）或列 22（SNPS）或列 24（SNP_ID_CURRENT）
awk = (f'awk -F"\\t" \'{{ s=""; f=0; '
       f'if ($21 ~ /^({alt_esc})[-,]?/ || $22 ~ /({alt_esc})(,|$| )/ || $24 ~ /^({alt_esc})$/) {{ '
       f'print $2 "\\t" $8 "\\t" $21 "\\t" $22 "\\t" $24 "\\t" $28 "\\t" $29 "\\t" $31 "\\t" $35 "\\t" $37 "\\t" $13 }} }}\' '
       f'"{CAT}"')
print("running awk extraction...")
res = subprocess.run(awk, shell=True, capture_output=True, text=True)
if res.returncode != 0:
    print("AWK ERROR:", res.stderr[-500:]); raise SystemExit(1)
lines = [l for l in res.stdout.split("\n") if l.strip()]
print(f"raw hit lines: {len(lines)}")
cols = ["PUBMEDID", "DISEASE_TRAIT", "STRONGEST_SNP", "SNPS", "SNP_ID_CURRENT",
        "P_VALUE", "PVALUE_MLOG", "OR_BETA", "MAPPED_TRAIT", "STUDY_ACC", "CHR_POS"]
hits = pd.DataFrame([l.split("\t") for l in lines], columns=cols)
hits["P_VALUE"] = pd.to_numeric(hits["P_VALUE"], errors="coerce")
hits["PVALUE_MLOG"] = pd.to_numeric(hits["PVALUE_MLOG"], errors="coerce")
hits["OR_BETA"] = hits["OR_BETA"].str.extract(r"([-0-9.]+)")[0].astype(float)
hits["CHR_POS"] = pd.to_numeric(hits["CHR_POS"], errors="coerce")
# 归属 SNP：优先 STRONGEST_SNP，其次 SNP_ID_CURRENT / SNPS
def assign_snp(row):
    for c in (row["STRONGEST_SNP"], row["SNP_ID_CURRENT"], row["SNPS"]):
        if pd.isna(c):
            continue
        m = re.search(r"(rs\d+)", str(c))
        if m and m.group(1) in set(snps):
            return m.group(1)
    return None
hits["snp"] = hits.apply(assign_snp, axis=1)
hits = hits.dropna(subset=["snp"])
print(f"assigned hits: {len(hits)}")

# ---- 3. 领域分类 ----
def classify(trait):
    t = str(trait).lower()
    dm = re.search(r"diabet|glucose|hba1c|glycated|insulin|glycem|homa|c-peptide", t)
    cad = re.search(r"coronary|myocardial|heart|angina|ischemic|ischaemic|atherosclero|cardiovascular|chd\b", t)
    lip = re.search(r"cholesterol|ldl|hdl|triglycer|lipid|apolipoprotein|dyslipid|non-hdl", t)
    bmi = re.search(r"bmi|body mass|waist|obes|adiposit|body fat|overweight", t)
    bp = re.search(r"blood pressure|hypertension|systolic|diastolic", t)
    flag = ""
    if dm: flag += "D"
    if cad: flag += "C"
    if lip: flag += "L"
    if bmi: flag += "B"
    if bp: flag += "P"
    if flag:
        return ("cardio", flag)
    return ("noncardio", "")
hits["domain"], hits["sub"] = zip(*hits["MAPPED_TRAIT"].fillna(hits["DISEASE_TRAIT"]).apply(classify))

# ---- 4. SNP 级汇总 ----
rows = []
for s in sorted(hits["snp"].unique()):
    h = hits[hits["snp"] == s]
    hcard = h[h["domain"] == "cardio"]
    hnc = h[h["domain"] == "noncardio"]
    best = h.loc[h["P_VALUE"].idxmin()]
    rows.append(dict(snp=s, n_traits=h["MAPPED_TRAIT"].nunique(),
                     n_studies=h["PUBMEDID"].nunique(), n_hits=len(h),
                     n_cardio=hcard["MAPPED_TRAIT"].nunique(),
                     n_diabetes=int(h["sub"].str.contains("D").sum()),
                     n_cad=int(h["sub"].str.contains("C").sum()),
                     n_lipid=int(h["sub"].str.contains("L").sum()),
                     n_bmi=int(h["sub"].str.contains("B").sum()),
                     n_bp=int(h["sub"].str.contains("P").sum()),
                     n_noncardio=hnc["MAPPED_TRAIT"].nunique(),
                     cross_domain=bool(len(hnc) > 0),
                     min_p=best["P_VALUE"], best_trait=best["MAPPED_TRAIT"],
                     top_cardio_traits="; ".join(hcard["MAPPED_TRAIT"].dropna().unique()[:6]),
                     top_noncardio_traits="; ".join(hnc["MAPPED_TRAIT"].dropna().unique()[:6])))
snp_sum = pd.DataFrame(rows)
snp_sum.to_csv(os.path.join(REPO, "results/m41_phewas_snp_20260817.csv"), index=False)

# ---- 5. 基因级汇总（候选→lead SNP 映射）----
gene_map = sub.drop_duplicates("symbol")[["symbol", "lead"]].rename(columns={"lead": "snp"})
gene_sum = gene_map.merge(snp_sum, on="snp", how="left")
gene_sum = gene_sum.merge(cand[["symbol", "outcome", "mr_p", "is_new"]].drop_duplicates("symbol"),
                          on="symbol", how="left")
gene_sum.to_csv(os.path.join(REPO, "results/m41_phewas_gene_20260817.csv"), index=False)
hits.to_csv(os.path.join(REPO, "results/m41_phewas_hits_20260817.csv"), index=False)

# ---- 6. Volcano ----
fig, ax = plt.subplots(figsize=(10, 7))
cols_dom = {"cardio": "#c00000", "noncardio": "#4a6b8a"}
for d, col in cols_dom.items():
    hd = hits[hits["domain"] == d]
    ax.scatter(hd["OR_BETA"], hd["PVALUE_MLOG"].fillna(0), s=18, alpha=0.55,
               color=col, label=f"{d} (n={len(hd)})", edgecolors="none")
# 标注重叠最高的 cardio 锚点
anchor = hits[hits["domain"] == "cardio"].sort_values("PVALUE_MLOG", ascending=False).drop_duplicates("MAPPED_TRAIT").head(12)
for _, r in anchor.iterrows():
    ax.annotate(str(r["MAPPED_TRAIT"])[:28], (r["OR_BETA"], r["PVALUE_MLOG"]),
                fontsize=7, color="#c00000", alpha=0.85,
                textcoords="offset points", xytext=(3, 3))
ax.axhline(7.3, color="0.6", ls="--", lw=0.7)  # P=5e-8
ax.text(0.01, 7.5, "P = 5×10⁻⁸ - GWAS Catalog inclusion line", fontsize=8, color="0.4", va="bottom")
ax.set_xlabel("OR or β (from GWAS Catalog)", fontsize=11)
ax.set_ylabel("−log₁₀(P)", fontsize=11)
ax.set_title("GWAS Catalog PheWAS of 13 candidate lead cis-eQTL SNPs\n"
             f"{len(hits)} associations across {hits['PUBMEDID'].nunique()} studies; "
             f"{snp_sum['cross_domain'].sum()}/{len(snp_sum)} SNPs with cross-domain hits",
             fontsize=12)
ax.legend(loc="upper right", frameon=True, fontsize=10)
ax.axvline(0, color="0.5", lw=0.6)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "20260817_FigS_phewas_volcano.png"), dpi=300, bbox_inches="tight")
plt.close(fig)
print("saved volcano figure")

# ---- 7. 摘要 ----
print("\n== SNP-level pleiotropy summary ==")
print(snp_sum[["snp", "n_traits", "n_studies", "n_cardio", "n_noncardio",
               "cross_domain", "min_p", "best_trait"]].to_string(index=False))
nc = snp_sum[snp_sum["cross_domain"]]
print(f"\ncross-domain (pleiotropy) SNPs: {len(nc)}")
print("== DONE M41 ==")
