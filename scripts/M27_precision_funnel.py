#!/usr/bin/env python3
# =============================================================================
# M27_precision_funnel.py — cis-MR coloc 操作特性：MR 显著度 / GWAS 峰强 → strong 率
# =============================================================================
# 目的（方法学正向叙事）：全转录组扫描的"精度-召回"校准曲线。
#   - strong 率 = 该层内 PP.H4≥0.8 的占比（coloc 后验，非真值）
#   - 随 MR p 收紧：strong 率单调上升？（校准）——可操作的"该采信哪些 hit"
#   - 对比：预过滤 stage-2 网格（M5, 819 对）的 strong 率 vs 全量 sig（4248 对）
#   → "全量扫描以召回换精度，找回 23 个新候选；预过滤网格精度高 4 倍"
# 输出：results/m27_precision_funnel_20260816.csv + stdout
# =============================================================================
import csv, os
import numpy as np

RES = "<repo-root>/results"

def f(x, d=float('nan')):
    try: return float(x)
    except (TypeError, ValueError): return d

rows = []
for on in ["t2d","cad","fbg"]:
    with open(f"{RES}/coloc_full_{on}_20260815.csv") as fh:
        for r in csv.DictReader(fh):
            r["outcome"] = on
            rows.append(r)
qc = [r for r in rows if r["ok"].strip()=="TRUE"]

def strong_rate(pairs):
    if not pairs: return float('nan'), 0, 0
    s = sum(1 for r in pairs if f(r["pp4"])>=0.8)
    return s/len(pairs), s, len(pairs)

# ---- 1. MR p 阈值扫描（全量 QC） ----
print("=== 1. MR p 阈值 → strong 率（全量 QC 对）===")
thresh = [0.5, 0.05, 0.01, 0.005, 1e-3, 5e-4, 1e-4, 1e-5]
rows_out = []
for t in thresh:
    sub = [r for r in qc if f(r["mr_p"]) < t]
    rate, s, n = strong_rate(sub)
    print(f"  mr_p<{t:<7}: strong {s:>3}/{n:<6} = {100*rate if not np.isnan(rate) else float('nan'):.2f}%")
    rows_out.append({"stratum": f"mr_p<{t}", "n": n, "strong": s,
                     "strong_rate": round(rate,5) if not np.isnan(rate) else ""})
    if n == 0: break

# ---- 2. 在 MR 显著内，GWAS 峰强分层 ----
print("\n=== 2. MR 显著内按 GWAS 峰 p 分层 ===")
sig = [r for r in qc if f(r["mr_p"]) < 0.05]
for g in [5e-8, 1e-6, 1e-5, 1e-4, 1e-3, 0.01, 1.0]:
    sub = [r for r in sig if f(r["gwas_min_p"]) < g]
    rate, s, n = strong_rate(sub)
    print(f"  gwas_min_p<{g:<7}: strong {s:>3}/{n:<6} = {100*rate:.2f}%")
    rows_out.append({"stratum": f"mr_p<0.05 & gwas_min_p<{g}", "n": n, "strong": s,
                     "strong_rate": round(rate,5)})

# ---- 3. 对比 stage-2 网格（M5 819 对） ----
print("\n=== 3. 预过滤 stage-2 网格（M5, 819 对）vs 全量 sig ===")
grid = list(csv.DictReader(open(f"{RES}/grid/transcript_coloc.csv")))
gq = []
for r in grid:
    rr = dict(r)
    if "PP.H4" in rr and "pp4" not in rr:
        rr["pp4"] = rr["PP.H4"]
    gq.append(rr)
gq = [r for r in gq if str(r.get("ok","TRUE")).strip()=="TRUE"]
rate, s, n = strong_rate(gq)
print(f"  stage-2 网格: strong {s}/{n} = {100*rate:.2f}%")
rows_out.append({"stratum": "stage2_grid(M5)", "n": n, "strong": s, "strong_rate": round(rate,5)})

rate, s, n = strong_rate(sig)
grid_rate, gs, gn = strong_rate(gq)
lift = rate / grid_rate if grid_rate else float('nan')
print(f"  全量 MR sig:  strong {s}/{n} = {100*rate:.2f}%  "
      f"(stage-2 网格精度为全量 sig 的 {1/lift if lift else float('nan'):.1f}×)")
rows_out.append({"stratum": "full_scan_sig", "n": n, "strong": s, "strong_rate": round(rate,5)})

# ---- 4. 全量分层计数（漏斗底） ----
for name, sub in [("all_QC", qc), ("nsig(0.05-0.5)", [r for r in qc if f(r["mr_p"])>=0.05 and f(r["mr_p"])<0.5]),
                  ("null(>=0.5)", [r for r in qc if f(r["mr_p"])>=0.5])]:
    rate, s, n = strong_rate(sub)
    print(f"  {name:<16}: strong {s:>3}/{n:<6} = {100*rate if not np.isnan(rate) else 0:.3f}%")
    rows_out.append({"stratum": name, "n": n, "strong": s,
                     "strong_rate": round(rate,5) if not np.isnan(rate) else ""})

with open(f"{RES}/m27_precision_funnel_20260816.csv","w",newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["stratum","n","strong","strong_rate"])
    w.writeheader(); w.writerows(rows_out)
print(f"\n已写 results/m27_precision_funnel_20260816.csv")
