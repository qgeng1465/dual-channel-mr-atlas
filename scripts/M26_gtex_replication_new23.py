#!/usr/bin/env python3
# =============================================================================
# M26_gtex_replication_new23.py — 23 个 new strong 的 GTEx 独立 eQTL 方向复现
# =============================================================================
# 目的：为 M25 注释出的 23 个 sig 内新增 strong coloc（均未在 GWAS Catalog
#   T2D/CAD 报道过）做独立 eQTL 源（GTEx v8）方向核查——评估这些候选基因的
#   可信度。纯离线：GTEx egenes + eQTLGen sig + OpenGWAS full（hg38）。
# 方法（与 M15 对 106 命中的设计一致）：
#   - eQTLGen lead eQTL SNP（max |Z|）为单工具
#   - GTEx：每基因跨 6 组织取 pval_nominal 最小的 egenes lead
#   - 两个 SNP 各自在结局 GWAS full 里的效应方向，与 coloc_full 的 mr_b
#     （eQTLGen-MR 方向）比对 → 方向一致率
# 诚实口径：这是候选评估，不是发现宣称；只报方向一致性 + 各方 p 值。
# 输出：results/m26_gtex_replication_new23_20260816.csv
# =============================================================================
import gzip, csv
import pandas as pd
import numpy as np

PROJ = "/data/qiushuogeng/projects/dual-channel-mr-atlas"
RES  = f"{PROJ}/results"
SIG  = f"{PROJ}/data/eqtlgen/cis-eQTL-significant.txt.gz"
GTEX = f"{PROJ}/data/gtex"
OUT  = f"{RES}/m26_gtex_replication_new23_20260816.csv"

TISSUES = ["Artery_Coronary","Liver","Pancreas","Adipose_Subcutaneous","Whole_Blood","Muscle_Skeletal"]

def f(x, d=float('nan')):
    try: return float(x)
    except (TypeError, ValueError): return d

# ---- 1. 23 个 new strong（复用 M24 口径）----
rows = []
for on in ["t2d","cad","fbg"]:
    with open(f"{RES}/coloc_full_{on}_20260815.csv") as fh:
        for r in csv.DictReader(fh):
            r["outcome"] = on
            rows.append(r)
known = set()
with open(f"{RES}/grid/transcript_coloc_hits.csv") as fh:
    for r in csv.DictReader(fh):
        known.add((r["gene"], r["outcome"]))
new23 = [r for r in rows if r["ok"].strip()=="TRUE" and f(r["pp4"])>=0.8
         and f(r["mr_p"])<0.05 and (r["gene"], r["outcome"]) not in known]
genes = {r["gene"] for r in new23}
print(f"new strong: {len(new23)}")

# ---- 2. eQTLGen lead per gene（max |Z|）----
eqtl_lead = {}
for ch in pd.read_csv(SIG, sep="\t", compression="gzip",
                      usecols=["Gene","SNP","AssessedAllele","OtherAllele","Zscore"],
                      chunksize=5_000_000):
    m = ch["Gene"].isin(genes)
    if m.sum()==0: continue
    for r in ch[m].itertuples(index=False):
        gg = r.Gene
        if not np.isfinite(r.Zscore): continue
        cur = eqtl_lead.get(gg)
        if cur is None or abs(r.Zscore) > abs(cur["z"]):
            eqtl_lead[gg] = {"snp": r.SNP, "a1": r.AssessedAllele, "a2": r.OtherAllele, "z": r.Zscore}
print(f"eQTLGen lead: {len(eqtl_lead)}")

# ---- 3. GTEx egenes：每基因跨组织 best lead ----
gtex_lead = {}
for t in TISSUES:
    fp = f"{GTEX}/{t}.egenes.txt.gz"
    d = pd.read_csv(fp, sep="\t", compression="gzip")
    d = d[d["gene_name"].isin([r["symbol"] for r in new23])]
    if len(d)==0: continue
    for _, r in d.iterrows():
        sym = r["gene_name"]
        if sym not in gtex_lead or f(r["pval_nominal"]) < f(gtex_lead[sym]["p"]):
            gtex_lead[sym] = {"tissue": t, "rsid": r["rs_id_dbSNP151_GRCh38p7"], "vid": r["variant_id"],
                              "alt": r["alt"], "ref": r["ref"], "slope": r["slope"], "se": r["slope_se"],
                              "p": r["pval_nominal"]}
print(f"GTEx lead（跨组织 best）: {len(gtex_lead)}")

# ---- 4. 需要查 hg38 GWAS 的 rsid ----
need = set()
for g in genes:
    need.add(eqtl_lead[g]["snp"])
for r in new23:
    if r["symbol"] in gtex_lead: need.add(gtex_lead[r["symbol"]]["rsid"])
need.discard(""); need.discard("."); need.discard("-")
print(f"需查 GWAS full: {len(need)} rsid")

gwas_of = {"t2d": {}, "cad": {}, "fbg": {}}
for fn in ["t2d","cad","fbg"]:
    with gzip.open(f"{PROJ}/data/opengwas/full/{fn}_full.gz","rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        if "hm_rsid" in hdr:                       # OpenGWAS harmonized（t2d/cad）
            di = {c:i for i,c in enumerate(hdr)}
            mp = {"Snp":"hm_rsid", "beta":"hm_beta", "p":"p_value",
                  "eff":"hm_effect_allele", "oth":"hm_other_allele"}
            def gcol(row, c): return row[di[mp.get(c, c)]]
        else:                                      # MAGIC 原始格式（fbg）
            di = {c:i for i,c in enumerate(hdr)}
            def gcol(row, c):
                mp = {"Snp":"Snp", "beta":"MainEffects", "p":"MainP",
                      "eff":"effect_allele", "oth":"other_allele"}
                return row[di[mp.get(c, c)]]
        for line in fh:
            p = line.rstrip("\n").split("\t")
            rs = gcol(p, "Snp" if "Snp" in hdr else "hm_rsid")
            if rs in need and rs not in gwas_of[fn]:
                gwas_of[fn][rs] = {"beta": f(gcol(p, "beta")), "p": f(gcol(p, "p")),
                                   "eff": gcol(p, "eff"), "oth": gcol(p, "oth")}
        print(f"  {fn}: 命中 {len(gwas_of[fn])}")

# ---- 5. 方向一致性 ----
def gtex_mr_dir(slope, alt, gwas_row):
    """GTEx eQTL（slope 相对 alt 等位基因）× GWAS 效应 → 基因对结局 MR 方向"""
    if gwas_row is None or not np.isfinite(gwas_row["beta"]): return None
    if alt == gwas_row["eff"]:   gw = gwas_row["beta"]
    elif alt == gwas_row["oth"]: gw = -gwas_row["beta"]
    else: return None
    if not np.isfinite(slope) or abs(slope) < 1e-12: return None
    return np.sign(gw) * np.sign(slope)

out_rows = []
for r in sorted(new23, key=lambda x: -f(x["pp4"])):
    g, sym, on = r["gene"], r["symbol"], r["outcome"]
    el = eqtl_lead.get(g, {})
    gl = gtex_lead.get(sym, {})
    mr_b = f(r["mr_b"]); mr_p = f(r["mr_p"])
    # eQTLGen-MR 方向（coloc_full 的 mr_b 已 harmonize）
    eq_dir = None if not np.isfinite(mr_b) else np.sign(mr_b)
    # GWAS 处 eQTLGen lead 的 p
    gw_e = gwas_of[on].get(el.get("snp",""), None)
    # GTEx-MR 方向
    gt_row = gwas_of[on].get(gl.get("rsid",""), None) if gl else None
    gt_dir = gtex_mr_dir(f(gl.get("slope", float('nan'))), gl.get("alt",""), gt_row)
    # 一致性
    if eq_dir is not None and gt_dir is not None:
        consistent = "consistent" if eq_dir == gt_dir else "conflicting"
    elif eq_dir is not None and gt_dir is None:
        consistent = "no_gtex_gwas_variant" if not gt_row else "allele_mismatch"
    else:
        consistent = "no_eqtlgen_dir"
    out_rows.append({
        "gene": g, "symbol": sym, "outcome": on,
        "pp4": f(r["pp4"]), "mr_b": r["mr_b"], "mr_p": r["mr_p"],
        "eqtlgen_lead_snp": el.get("snp",""), "eqtlgen_z": el.get("z",""),
        "eqtlgen_gwas_p": gw_e["p"] if gw_e else "",
        "gtex_tissue": gl.get("tissue",""), "gtex_lead_rsid": gl.get("rsid",""),
        "gtex_slope": gl.get("slope",""), "gtex_p": gl.get("p",""),
        "gtex_lead_gwas_p": gt_row["p"] if gt_row else "",
        "eqtlgen_mr_dir": eq_dir, "gtex_mr_dir": gt_dir,
        "direction": consistent,
    })
df = pd.DataFrame(out_rows)
df.to_csv(OUT, index=False)
print(f"\n=== 方向一致性汇总 ===")
print(df["direction"].value_counts().to_string())
print(f"\n=== 明细 ===")
for _, r in df.iterrows():
    print(f"{r['symbol']:<10}{r['outcome']:<5} pp4={f(r['pp4']):.3f} eQTLGen-MR dir={r['eqtlgen_mr_dir']} "
          f"GTEx {str(r['gtex_tissue'])[:18]:<18} slope={f(r['gtex_slope']):+.3f} p={f(r['gtex_p']):.1e} "
          f"GWASp@GTExlead={f(r['gtex_lead_gwas_p']):.1e} → {r['direction']}")
print(f"\n已写 {OUT}")
