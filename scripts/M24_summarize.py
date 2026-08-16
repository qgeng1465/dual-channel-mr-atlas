#!/usr/bin/env python3
# =============================================================================
# M24_summarize.py — 全量 coloc 扫描结果汇总 + 与 106 已知 strong 交叉验证
# =============================================================================
# 输入：results/coloc_full_{t2d,cad,fbg}_20260815.csv（M23 输出）
# 交叉验证：MR sig 集内的 strong（PP.H4≥0.8）应与 transcript_coloc_hits.csv 的
#   106 命中一致（健全性检查，coloc 是确定性的——同输入应同输出）。
# 输出：results/coloc_full_summary_20260815.csv + 终端摘要
# =============================================================================
import csv, os, sys

RES = "/data/qiushuogeng/projects/dual-channel-mr-atlas/results"

def f(x, default=float('nan')):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default

# ---- 1. 读三个结局 ----
rows = []
for on in ["t2d", "cad", "fbg"]:
    p = os.path.join(RES, f"coloc_full_{on}_20260815.csv")
    if not os.path.exists(p):
        print(f"[FATAL] 缺 {p}——等 M23 完成再跑")
        sys.exit(1)
    with open(p) as fh:
        for r in csv.DictReader(fh):
            r["outcome"] = on
            rows.append(r)
n_all = len(rows)
qc = [r for r in rows if r["ok"].strip() == "TRUE"]
print(f"全量对: {n_all} | QC 通过(nsnp≥10): {len(qc)} ({100*len(qc)/max(1,n_all):.1f}%)")

# ---- 2. 与 106 已知 strong 交叉验证 ----
known = {}
with open(os.path.join(RES, "grid/transcript_coloc_hits.csv")) as fh:
    for r in csv.DictReader(fh):
        known[(r["gene"], r["outcome"])] = f(r["PP.H4"])

# 口径修正（2026-08-16）：sig/nsig/null 必须基于 QC 通过子集 qc，与 n_qc 及 strong_* 自洽。
# 原用全量 rows，QC 失败对（RP11-764K9.1/RP11-87H9.2，mr_p=0.143 恰落灰区）被计入，
# 使 n_sig+n_nsig+n_null=31,373 ≠ n_qc=31,371（学术诚信审计指出）。数字不变（sig 4,248 全过 QC）。
sig   = [r for r in qc if r["mr_p"] != "" and f(r["mr_p"]) < 0.05]
nsig  = [r for r in qc if r["mr_p"] != "" and f(r["mr_p"]) >= 0.05 and f(r["mr_p"]) < 0.5]
null  = [r for r in qc if r["mr_p"] != "" and f(r["mr_p"]) >= 0.5]

strong_sig = [r for r in sig if r["ok"].strip() == "TRUE" and f(r["pp4"]) >= 0.8]
strong_ns  = [r for r in nsig if r["ok"].strip() == "TRUE" and f(r["pp4"]) >= 0.8]
strong_null= [r for r in null if r["ok"].strip() == "TRUE" and f(r["pp4"]) >= 0.8]

sk  = {(r["gene"], r["outcome"]) for r in strong_sig}
hk  = set(known.keys())
print(f"\nMR 显著: {len(sig)} | 其中 strong: {len(strong_sig)}")
print(f"灰区(0.05≤p<0.5): {len(nsig)} | strong: {len(strong_ns)}")
print(f"阴性(p≥0.5): {len(null)} | strong: {len(strong_null)}")
print(f"全量 strong 总数: {len(strong_sig)+len(strong_ns)+len(strong_null)}")
print(f"\n健全性交叉验证:")
print(f"  重现已知 106: {len(sk & hk)}/{len(hk)}")
print(f"  sig 内新 strong(不在 106): {len(sk - hk)}")
print(f"  已知 106 未重现: {len(hk - sk)}")

# 精度（操作特性）
prec = len(strong_sig)/max(1, len(sig))
print(f"\n精度(操作特性): MR 显著集内 strong 占比 = {len(strong_sig)}/{len(sig)} = {100*prec:.1f}%")

# ---- 3. 按结局统计 ----
print(f"\n按结局:")
for on in ["t2d", "cad", "fbg"]:
    o   = [r for r in rows if r["outcome"] == on]
    oq  = [r for r in o if r["ok"].strip() == "TRUE"]
    osg = [r for r in o if r["ok"].strip() == "TRUE" and r["mr_p"] != "" and f(r["mr_p"])<0.05 and f(r["pp4"])>=0.8]
    ong = [r for r in o if r["ok"].strip() == "TRUE" and r["mr_p"] != "" and f(r["mr_p"])>=0.05 and f(r["pp4"])>=0.8]
    # GWAS 峰显著率（strong 内）
    gsig = sum(1 for r in osg if r["gwas_min_p"] != "" and f(r["gwas_min_p"]) < 5e-8)
    print(f"  {on}: 全部 {len(o)} | QC {len(oq)} | strong(全部) {len(osg)+len(ong)} "
          f"(sig内 {len(osg)}, sig外 {len(ong)}) | strong 内 GWAS峰显著 {gsig}/{len(osg) if osg else 0}")

# ---- 3.5 M20 vs M23 一致性检查（健全性：同基因同结局对 pp4 应完全一致）----
# M20 抽样了 MR 显著集外的 6000 对；M23 全量含这些对。coloc 是确定性算法，
# 同输入（eQTLGen 同数据源 + 同 GWAS + 同 harmonize/p12）应产出相同 pp4。
m20 = {}
m20p = os.path.join(RES, "feasibility_pilot_20260815.csv")
if not os.path.exists(m20p):
    m20p = os.path.join(RES, "archive/feasibility_pilot_20260815.csv")
with open(m20p) as fh:
    for r in csv.DictReader(fh):
        m20[(r["gene"], r["outcome"])] = r
m23_idx = {(r["gene"], r["outcome"]): r for r in rows if r["ok"].strip() == "TRUE" and r["pp4"] != ""}
common = set(m20.keys()) & set(m23_idx.keys())
match = miss = 0
maxdiff = 0.0
for k in common:
    if m20[k]["ok"].strip() != "TRUE" or m20[k]["pp4"] == "":
        miss += 1
        continue
    d = abs(f(m20[k]["pp4"]) - f(m23_idx[k]["pp4"]))
    maxdiff = max(maxdiff, d)
    if d < 1e-6:
        match += 1
    else:
        miss += 1
print(f"\nM20↔M23 一致性: 共同对 {len(common)} | pp4 一致 {match} | 不一致 {miss} | 最大差 {maxdiff:.2e}")

# ---- 4. 落盘汇总 CSV ----
with open(os.path.join(RES, "coloc_full_summary_20260815.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["metric", "value"])
    w.writerow(["n_all", n_all])
    w.writerow(["n_qc", len(qc)])
    w.writerow(["n_sig", len(sig)])
    w.writerow(["n_nsig", len(nsig)])
    w.writerow(["n_null", len(null)])
    w.writerow(["strong_sig", len(strong_sig)])
    w.writerow(["strong_nsig", len(strong_ns)])
    w.writerow(["strong_null", len(strong_null)])
    w.writerow(["strong_total", len(strong_sig)+len(strong_ns)+len(strong_null)])
    w.writerow(["precision_sig", round(prec, 4)])
    w.writerow(["reproduced_known", len(sk & hk)])
    w.writerow(["known_total", len(hk)])
    w.writerow(["new_in_sig", len(sk - hk)])
    w.writerow(["known_missing", len(hk - sk)])
print(f"\n已写 results/coloc_full_summary_20260815.csv")
