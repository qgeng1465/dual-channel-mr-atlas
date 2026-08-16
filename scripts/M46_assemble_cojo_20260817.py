#!/usr/bin/env python3
# =============================================================================
# M46_assemble_cojo_20260817.py — GCTA-COJO 结果组装（供 M49 FigS1 与手稿使用）
# =============================================================================
# 读取 <scratch>/cojo/{tag}_A.log（p < 5e-8）与 {tag}_B.log（p < 1e-4）的 stepwise
# 选中信号数，{tag}_B.jma.cojo 的独立信号清单，{tag}_B.cma.cojo 中 lead SNP 的
# 条件 p（pC；col-linear 时 bC/pC=NA），组装为 results/m46_cojo_20260817.csv。
# 列：symbol,outcome,lead_snp,n_index_p5e8,n_index_p1e4,lead_in_index,lead_pC,
#     lead_pJ,index_snps,notes
# =============================================================================
import os
import re
import pandas as pd

REPO = "/data/qiushuogeng/projects/dual-channel-mr-atlas"
SCR = "/data/qiushuogeng/tmp"

LOCI = [  # symbol|outcome|lead
    ("RBM6", "t2d", "rs10049087"),
    ("CD101", "t2d", "rs10494191"),
    ("CNNM2", "cad", "rs11191447"),
    ("PLAUR", "cad", "rs4760"),
    ("RIC8A", "cad", "rs6598075"),
    ("LAMC1", "cad", "rs10458355"),
]

RE_FINALLY = re.compile(r"Finally,\s+(\d+)\s+associated SNPs?\s+are selected\.")
RE_NONE = re.compile(r"No SNPs have been selected")


def n_selected(log_path):
    if not os.path.exists(log_path):
        return None
    txt = open(log_path).read()
    m = RE_FINALLY.search(txt)
    if m:
        return int(m.group(1))
    if RE_NONE.search(txt):
        return 0
    return None


def read_jma(path):
    """返回 (SNP 列表, {SNP: pJ})；pJ=联合条件 p（jma 文件）。"""
    if not os.path.exists(path):
        return [], {}
    d = pd.read_csv(path, sep="\t")
    pj = dict(zip(d["SNP"], d["pJ"])) if "pJ" in d.columns else {}
    return list(d["SNP"]), pj


def lead_pc(path, lead):
    """从 cma.cojo 取 lead 的条件 p（bC/pC；col-linear 时为 NA）。"""
    if not os.path.exists(path):
        return None
    d = pd.read_csv(path, sep="\t")
    r = d[d["SNP"] == lead]
    return None if r.empty else r["pC"].iloc[0]


rows = []
for sym, out, lead in LOCI:
    tag = sym.lower()
    n5 = n_selected(os.path.join(SCR, "cojo", f"{tag}_A.log"))
    n1 = n_selected(os.path.join(SCR, "cojo", f"{tag}_B.log"))
    idx, pj_map = read_jma(os.path.join(SCR, "cojo", f"{tag}_B.jma.cojo"))
    pC = lead_pc(os.path.join(SCR, "cojo", f"{tag}_B.cma.cojo"), lead)
    lead_in = lead in idx
    pJ = pj_map.get(lead) if lead_in else None
    notes = ""
    if n1 is None:
        notes = "COJO failed/incomplete"
    elif lead_in:
        notes = "lead indexes an independent signal"
    elif pC is None:
        notes = "lead excluded by collinearity/freq filter (not in conditional output)"
    elif pd.isna(pC):
        notes = "lead col-linear with a selected signal (bC/pC=NA)"
    else:
        notes = f"lead not significant conditional on index signal (pC={pC:.2g})"
    rows.append(dict(symbol=sym, outcome=out, lead_snp=lead,
                     n_index_p5e8=n5, n_index_p1e4=n1,
                     lead_in_index=lead_in, lead_pC=pC, lead_pJ=pJ,
                     index_snps=";".join(idx), notes=notes))
    print(f"{sym:8s} p5e8={n5} p1e4={n1} lead_in_idx={lead_in} lead_pC={pC} pJ={pJ} idx={idx}")

df = pd.DataFrame(rows)
out = os.path.join(REPO, "results", "m46_cojo_20260817.csv")
df.to_csv(out, index=False)
print(f"== DONE -> {out} ==")
