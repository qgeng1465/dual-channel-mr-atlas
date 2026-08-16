#!/usr/bin/env python3
# =============================================================================
# M28_finngen_replication_new23.py — 23 个 new strong 候选的 FinnGen R11 独立队列方向复现
# =============================================================================
# 目的（路径 B 强化）：M25 注释出的 23 个 sig 内新增 strong coloc（未在 GWAS Catalog
#   T2D/CAD 报道）需要独立证据。M26 做了 eQTL 侧独立复现（GTEx，同结局 GWAS）；
#   本脚本做**结局侧独立队列复现**（FinnGen R11 全量 sumstats，芬兰人群）：
#     - 工具不变：eQTLGen cis lead eQTL（max |Z|）
#     - 变异级复现：同一 lead SNP 在原始结局 GWAS（t2d=ebi-a-GCST006867,
#       cad=ebi-a-GCST005194）与 FinnGen R11（T2D / I9_CHD）的 β 方向一致率
#     - 基因级复现：sign(mr_b) vs sign(beta_finn × Z_eQTL)（FinnGen 队列的 MR 方向）
#   FBG 无 FinnGen 端点 → 如实跳过（同 M7）。
# 诚实口径：候选评估非发现宣称；比方向不比幅度；同时报告 FinnGen p 值。
# 输入：
#   results/coloc_full_*_20260815.csv（mr_b/mr_p/pp4）
#   results/grid/transcript_coloc_hits.csv（known 106）
#   data/eqtlgen/cis-eQTL-significant.txt.gz（lead a1/a2/z）
#   data/opengwas/full/{t2d,cad}_full.gz（原结局 hm_*，hg38 坐标）
#   data/finngen/finngen_R11_{T2D,I9_CHD}.gz（独立队列）
# 输出：results/m28_finngen_replication_new23_20260816.csv
# =============================================================================
import gzip, csv, os, sys
import pandas as pd
import numpy as np

PROJ = "<repo-root>"
RES  = f"{PROJ}/results"
SIG  = f"{PROJ}/data/eqtlgen/cis-eQTL-significant.txt.gz"
FINN = f"{PROJ}/data/finngen"
FINN_STUDY = {"t2d": "finngen_R11_T2D", "cad": "finngen_R11_I9_CHD"}
OUT  = f"{RES}/m28_finngen_replication_new23_20260816.csv"

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

# ---- 2. eQTLGen lead per gene（a1/a2/z）----
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
print(f"eQTLGen lead: {len(eqtl_lead)}（非 rsid lead: {sum(1 for v in eqtl_lead.values() if not v['snp'].startswith('rs'))}）")

# 每基因归属结局 → 需要哪些 lead SNP 在哪个文件里查
need_by_on = {"t2d": set(), "cad": set()}
gene_of_on = {r["gene"]: r["outcome"] for r in new23}
for g in genes:
    on = gene_of_on.get(g)
    if on in need_by_on and g in eqtl_lead:
        need_by_on[on].add(eqtl_lead[g]["snp"])
print(f"各结局需查 lead: t2d={len(need_by_on['t2d'])} cad={len(need_by_on['cad'])}")

# ---- 3. 原结局 GWAS full 读 lead SNP（hm_*，含 hg38 坐标）----
def load_gwas_leads(on, need_rs):
    """返回 rsid -> {beta, p, eff, oth, chrom38, pos38}"""
    out = {}
    with gzip.open(f"{PROJ}/data/opengwas/full/{on}_full.gz","rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        di = {c:i for i,c in enumerate(hdr)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            rs = p[di["hm_rsid"]]
            if rs in need_rs and rs not in out:
                out[rs] = {"beta": f(p[di["hm_beta"]]), "p": f(p[di["p_value"]]),
                           "eff": p[di["hm_effect_allele"]], "oth": p[di["hm_other_allele"]],
                           "chrom38": p[di["hm_chrom"]], "pos38": p[di["hm_pos"]]}
            if len(out) >= len(need_rs):
                break
    return out

gwas_lead = {on: load_gwas_leads(on, need_by_on[on]) for on in ["t2d","cad"]}
for on in ["t2d","cad"]:
    print(f"原结局 {on} full 命中 lead: {len(gwas_lead[on])}/{len(need_by_on[on])}")

# 原结局里能拿到 hg38 坐标的 lead → 用于 Finngen (chrom,pos) 精确匹配
pos_of_rs = {}   # rsid -> (chrom, pos)
for on in ["t2d","cad"]:
    for rs, d in gwas_lead[on].items():
        try:
            pos_of_rs[rs] = (int(float(d["chrom38"])), int(float(d["pos38"])))
        except (ValueError, TypeError):
            pass
# 反向：给 Finngen 扫描用的 (chrom,pos)->rsid 集合
pos_need = {v: k for k, v in pos_of_rs.items()}
print(f"可经 hg38 坐标匹配 Finngen 的 lead: {len(pos_need)}")

# ---- 4. Finngen R11 扫描 lead SNP ----
finn_of = {"t2d": {}, "cad": {}}   # on -> rsid -> {beta,p,ref,alt}
for on, fn in FINN_STUDY.items():
    fp = f"{FINN}/{fn}.gz"
    if not os.path.exists(fp):
        print(f"  [SKIP] {fp} 不存在")
        continue
    n = 0
    with gzip.open(fp,"rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ci = {c:i for i,c in enumerate(hdr)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            # (chrom,pos) 精确匹配优先
            try:
                key = (int(p[ci["#chrom"]]), int(p[ci["pos"]]))
            except (ValueError, TypeError, IndexError):
                continue
            rs = pos_need.get(key)
            if rs is None:
                # 兜底：rsids 字段逗号分隔的精确成员匹配（防 rs1234 子串误配 rs12345）
                import re
                rsids = p[ci["rsids"]]
                for cand in need_by_on[on]:
                    if re.search(rf"(^|,){re.escape(cand)}(,|$)", rsids):
                        rs = cand; break
            if rs is None or rs in finn_of[on]:
                continue
            finn_of[on][rs] = {"beta": f(p[ci["beta"]]), "p": f(p[ci["pval"]]),
                               "ref": p[ci["ref"]], "alt": p[ci["alt"]]}
            n += 1
            if n >= len(need_by_on[on]):
                break
    print(f"FinnGen {fn}: 命中 {len(finn_of[on])}/{len(need_by_on[on])} lead")

# ---- 5. 对齐 + 方向判定 ----
def align_gwas(eq, gw):
    """把原结局 hm_beta 对齐到 eQTLGen a1 → 返回 (beta_aligned, ok)"""
    if gw is None or not np.isfinite(gw["beta"]): return None, False
    a1 = eq["a1"]
    if a1 == gw["eff"]:  return gw["beta"], True
    if a1 == gw["oth"]:  return -gw["beta"], True
    return None, False

def align_finn(eq, fn):
    """把 FinnGen beta（相对 alt）对齐到 eQTLGen a1 → (beta_aligned, p, ok)"""
    if fn is None: return None, None, False
    a1 = eq["a1"]
    if a1 == fn["alt"]:  return fn["beta"], fn["p"], True
    if a1 == fn["ref"]:  return -fn["beta"], fn["p"], True
    return None, None, False

out_rows = []
for r in sorted(new23, key=lambda x: -f(x["pp4"])):
    g, sym, on = r["gene"], r["symbol"], r["outcome"]
    if on not in FINN_STUDY:
        out_rows.append({**{k: r.get(k,"") for k in ["gene","symbol","outcome","pp4","mr_b","mr_p"]},
                         "eqtlgen_lead_snp": eqtl_lead.get(g,{}).get("snp",""),
                         "eqtlgen_z": eqtl_lead.get(g,{}).get("z",""),
                         "orig_beta_a1":"", "orig_p":"", "finn_beta_a1":"", "finn_p":"",
                         "variant_dir":"fbg_no_finngen", "mr_dir_finn":"", "variant_replicated":"",
                         "mr_replicated":"", "note":""})
        continue
    eq = eqtl_lead.get(g, {})
    if not eq:
        out_rows.append({**{k: r.get(k,"") for k in ["gene","symbol","outcome","pp4","mr_b","mr_p"]},
                         "eqtlgen_lead_snp":"", "eqtlgen_z":"", "orig_beta_a1":"", "orig_p":"",
                         "finn_beta_a1":"", "finn_p":"", "variant_dir":"no_eqtl_lead",
                         "mr_dir_finn":"", "variant_replicated":"", "mr_replicated":"", "note":""})
        continue
    gw = gwas_lead[on].get(eq["snp"], None)
    fn = finn_of[on].get(eq["snp"], None)
    bo, ok_o = align_gwas(eq, gw)
    bf, pf, ok_f = align_finn(eq, fn)
    mr_b = f(r["mr_b"])
    # 变异级复现：原结局 vs FinnGen 同一 lead 的 β 方向
    if ok_o and ok_f:
        vdir = "consistent" if (np.sign(bo) == np.sign(bf)) else "conflicting"
        vrep = "yes" if vdir == "consistent" else "no"
    elif ok_o and not ok_f:
        vdir, vrep = "finn_unalignable_or_missing", ""
    else:
        vdir, vrep = "orig_missing", ""
    # 基因级复现：FinnGen MR 方向 = sign(z) × sign(beta_finn_a1)
    if ok_f and np.isfinite(mr_b) and eq["z"] != 0 and np.isfinite(eq["z"]):
        mr_dir_finn = int(np.sign(bf) * np.sign(eq["z"]))
        mrep = "yes" if mr_dir_finn == np.sign(mr_b) else "no"
    else:
        mr_dir_finn, mrep = "", ""
    out_rows.append({
        "gene": g, "symbol": sym, "outcome": on,
        "pp4": r["pp4"], "mr_b": r["mr_b"], "mr_p": r["mr_p"],
        "eqtlgen_lead_snp": eq["snp"], "eqtlgen_a1": eq["a1"], "eqtlgen_a2": eq["a2"], "eqtlgen_z": eq["z"],
        "orig_beta_a1": bo if ok_o else "", "orig_p": gw["p"] if gw else "",
        "finn_beta_a1": bf if ok_f else "", "finn_p": pf if ok_f else "",
        "variant_dir": vdir, "variant_replicated": vrep,
        "mr_dir_orig": np.sign(mr_b) if np.isfinite(mr_b) else "",
        "mr_dir_finn": mr_dir_finn, "mr_replicated": mrep,
        "note": "",
    })

df = pd.DataFrame(out_rows)
df.to_csv(OUT, index=False)

# ---- 6. 汇总 ----
print(f"\n=== Finngen R11 结局侧复现汇总 ===")
rep = df[df["variant_replicated"].isin(["yes","no"])]
print(f"变异级（原 vs Finngen 同 lead β 方向）：一致 {sum(rep['variant_replicated']=='yes')}/{len(rep)}")
mr = df[df["mr_replicated"].isin(["yes","no"])]
print(f"基因级（FinnGen MR 方向 vs mr_b）：一致 {sum(mr['mr_replicated']=='yes')}/{len(mr)}")
fp = df[(df["finn_p"]!="") & (df["finn_p"].apply(lambda x: f(x)<0.05))]
print(f"FinnGen lead p<0.05: {len(fp)}（共 {sum(df['finn_p']!='')} 可测）")
print(f"\n=== 明细 ===")
def fmt_num(v, fmt):
    return format(v, fmt) if v != "" else "—"
for _, r in df.iterrows():
    print(f"{r['symbol']:<10}{r['outcome']:<5} pp4={f(r['pp4']):.3f} "
          f"lead={str(r['eqtlgen_lead_snp']):<12} "
          f"origβ={fmt_num(r['orig_beta_a1'], '+.3f'):<8} "
          f"finnβ={fmt_num(r['finn_beta_a1'], '+.3f'):<8} finnP={fmt_num(r['finn_p'], '.1e'):<10} "
          f"变异级={r['variant_dir']:<24} MR级={r['mr_replicated'] or '—'}")
print(f"\n已写 {OUT}")
