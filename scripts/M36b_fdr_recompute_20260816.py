#!/usr/bin/env python3
# =============================================================================
# M36b_fdr_recompute_20260816.py — 按预注册口径（分结局 BH-FDR q<0.05）重算
# cis-MR 显著性，取代代码中实际使用的原始 p<0.05 口径（4,248 对）。
#
#   FDR-core  = { ok & per-outcome BH-FDR(q)<0.05 }  → 982 对（t2d 394/cad 576/fbg 12）
#   strong    = PP.H4 ≥ 0.8
#   106 已知  = stage-2 grid strong（transcript_coloc_hits.csv）
#   15 候选   = strong ∩ FDR-core − 106 已知（= 23 − 8 个 FDR 掉出的）
#
# 产出：
#   results/fdr_core_20260816.csv                   （982 对 + padj + strong 标记）
#   results/candidate15_replication_20260816.csv    （15 候选 + GTEx/FinnGen 复现）
#   results/m36b_funnel_20260816.csv                （名义阈值校准曲线 + FDR-core + grid）
#   results/m36b_summary_20260816.csv               （关键数字供 FACTS/manuscript 引用）
# =============================================================================
import os, sys
import numpy as np
import pandas as pd

BASE = "/data/qiushuogeng/projects/dual-channel-mr-atlas"
RES = f"{BASE}/results"

OUTCOMES = ["t2d", "cad", "fbg"]

def bh_fdr(pvals):
    """Benjamini–Hochberg FDR on a vector of p-values (per-outcome).
    Must sort before the suffix-min accumulate, then map back to original order.
    out[i] = smallest q such that rank(i)-sorted p <= q."""
    p = np.asarray(pvals, float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = adj
    return out


def load_full(outcome):
    d = pd.read_csv(f"{RES}/coloc_full_{outcome}_20260815.csv")
    d["outcome"] = d["outcome"].astype(str).str.strip()
    assert d["outcome"].eq(outcome).all(), f"outcome mismatch in {outcome}"
    return d


# ---------------------------------------------------------------------------
# 1. 分结局 BH-FDR → FDR-core
# ---------------------------------------------------------------------------
frames = []
for out in OUTCOMES:
    d = load_full(out)
    q = d[d["ok"]].copy()
    q["padj"] = bh_fdr(q["mr_p"].values)
    fdr = q[q["padj"] < 0.05].copy()
    fdr["strong"] = fdr["pp4"] >= 0.8
    frames.append(fdr)
    print(f"{out}: QC {len(q)} | FDR-core {len(fdr)} | strong-in-core {int(fdr['strong'].sum())}")

fdr_core = pd.concat(frames, ignore_index=True)
n_sig = len(fdr_core)
n_strong = int(fdr_core["strong"].sum())
print(f"\nFDR-core total: {n_sig} pairs | strong: {n_strong} | yield {n_strong/n_sig:.4f}")

fdr_core_out = [
    "gene", "symbol", "outcome", "mr_b", "mr_p", "padj", "gwas_min_p",
    "eqtl_F_max", "nsnp", "pp4", "strong",
]
fdr_core[fdr_core_out].to_csv(f"{RES}/fdr_core_20260816.csv", index=False)

# ---------------------------------------------------------------------------
# 2. 与 106 已知 strong 交叉验证
# ---------------------------------------------------------------------------
known = set()
with open(f"{RES}/grid/transcript_coloc_hits.csv") as fh:
    grid = pd.read_csv(fh)
for _, r in grid.iterrows():
    known.add((str(r["gene"]), str(r["outcome"]).strip()))
print(f"\nknown (stage-2 grid strong): {len(known)}")

core_strong = fdr_core[fdr_core["strong"]].copy()
core_strong_keys = set(zip(core_strong["gene"].astype(str), core_strong["outcome"].astype(str)))
reproduced = core_strong_keys & known
missing = known - core_strong_keys
new = core_strong_keys - known
print(f"known reproduced in FDR-core strong: {len(reproduced)}/{len(known)}")
if missing:
    print("  !! known MISSING under FDR:", sorted(missing)[:20])
print(f"new strong in FDR-core (not in 106): {len(new)}")
assert len(new) == 15, f"expected 15 candidate genes, got {len(new)}"

# ---------------------------------------------------------------------------
# 3. 15 候选 vs M25 的 23 新 strong（应 = 23 − 8 个 FDR 掉出的）
# ---------------------------------------------------------------------------
m25 = pd.read_csv(f"{RES}/m25_new_strong_annotation_20260816.csv")
m25_keys = set(zip(m25["gene"].astype(str), m25["outcome"].astype(str)))
assert new <= m25_keys, f"FDR-core new strong outside M25 list: {new - m25_keys}"
dropped = m25_keys - new
print(f"M25 23 候选 → FDR-core 15 候选（掉 {len(dropped)}）：")
for g, o in sorted(dropped):
    row = fdr_core if False else None
    # look up padj for the dropped pair from full scan
    sym = m25[(m25["gene"] == g) & (m25["outcome"] == o)]["symbol"].iloc[0]
    print(f"   {sym}({o})")

# ---------------------------------------------------------------------------
# 4. 15 候选复现统计（GTEx m26 / FinnGen m28）
# ---------------------------------------------------------------------------
cand = core_strong[~core_strong[["gene", "outcome"]].apply(
    lambda r: (str(r["gene"]), str(r["outcome"])) in known, axis=1)].copy()
cand["is_new"] = True
cand = cand.merge(
    m25[["gene", "outcome", "tier", "catalog_gene_t2dcad", "nearest_cat_kb", "top_snp"]],
    on=["gene", "outcome"], how="left",
)

m26 = pd.read_csv(f"{RES}/m26_gtex_replication_new23_20260816.csv")
m28 = pd.read_csv(f"{RES}/m28_finngen_replication_new23_20260816.csv")
for df in (m26, m28):
    df["gene"] = df["gene"].astype(str)

c = cand.merge(m26[["gene", "outcome", "direction", "gtex_p"]],
               on=["gene", "outcome"], how="left")
c = c.merge(m28[["gene", "outcome", "variant_replicated", "mr_replicated", "finn_p", "note"]],
            on=["gene", "outcome"], how="left")

print(f"\n=== 15 候选复现统计 ===")
_dir = c["direction"].fillna("")
n_meas = int(_dir.isin(["consistent", "conflicting"]).sum())
n_consistent = int(_dir.eq("consistent").sum())
n_conflict = int(_dir.eq("conflicting").sum())
print(f"GTEx 可测 {n_meas}（consistent {n_consistent} / conflicting {n_conflict} / 其余不可比）")
n_meas_f = int(c["finn_p"].notna().sum())
n_var_rep = int(c["variant_replicated"].eq("yes").sum())
n_mr_rep = int(c["mr_replicated"].eq("yes").sum())
n_p05 = int((c["finn_p"].fillna(1.0) < 0.05).sum())
print(f"FinnGen 可测 {n_meas_f} | variant-level 复现 {n_var_rep} | gene-level(MR) 复现 {n_mr_rep} | finn p<0.05 {n_p05}")

sub = c[c["finn_p"].notna()].sort_values("finn_p")[["symbol", "outcome", "finn_p"]]
for _, r in sub.iterrows():
    print(f"   {r['symbol']}({r['outcome']}) finn_p={r['finn_p']:.2e}")

c.to_csv(f"{RES}/candidate15_replication_20260816.csv", index=False)

# ---------------------------------------------------------------------------
# 5. 新 funnel：名义阈值校准曲线 + FDR-core 点 + grid 点
# ---------------------------------------------------------------------------
funnel_old = pd.read_csv(f"{RES}/m27_precision_funnel_20260816.csv")
nominal = funnel_old[funnel_old["stratum"].isin([
    "mr_p<0.5", "mr_p<0.05", "mr_p<0.01", "mr_p<0.005", "mr_p<0.001",
    "mr_p<0.0005", "mr_p<0.0001", "mr_p<1e-05"])].copy()
nominal["threshold"] = nominal["stratum"].str.replace("mr_p<", "").astype(float)
funnel = nominal[["threshold", "n", "strong", "strong_rate"]].rename(
    columns={"strong_rate": "yield"}).copy()
funnel["kind"] = "nominal"
fdr_row = pd.DataFrame([{"threshold": np.nan, "n": n_sig, "strong": n_strong,
                         "yield": n_strong / n_sig, "kind": "fdr_core"}])
grid_row = pd.DataFrame([{"threshold": np.nan, "n": 818, "strong": 106,
                          "yield": 106 / 818, "kind": "stage2_grid"}])
funnel = pd.concat([funnel, fdr_row, grid_row], ignore_index=True)
funnel.to_csv(f"{RES}/m36b_funnel_20260816.csv", index=False)
print(f"\nfunnel saved ({len(funnel)} rows incl FDR-core {n_strong/n_sig:.4f} and grid 12.96%)")

# ---------------------------------------------------------------------------
# 6. 汇总数字（供 FACTS / manuscript 直接引用）
# ---------------------------------------------------------------------------
summary = {
    "n_pairs_total": 31373,
    "n_pairs_qc": 31371,
    "n_sig_raw_p05": 4248,
    "n_sig_fdr": n_sig,
    "n_strong_fdr_core": n_strong,
    "yield_fdr": round(n_strong / n_sig, 4),
    "n_known_reproduced": len(reproduced),
    "n_candidates": len(new),
    "n_candidates_dropped_by_fdr": len(dropped),
    "gtex_measured": n_meas,
    "gtex_consistent": n_consistent,
    "gtex_conflict": n_conflict,
    "finngen_measured": n_meas_f,
    "finngen_variant_replicated": n_var_rep,
    "finngen_gene_replicated": n_mr_rep,
    "finngen_p05": n_p05,
}
pd.DataFrame([summary]).to_csv(f"{RES}/m36b_summary_20260816.csv", index=False)
print("\n=== m36b_summary_20260816.csv ===")
for k, v in summary.items():
    print(f"  {k}: {v}")

# ---------------------------------------------------------------------------
# 7. 断言（防回归）
# ---------------------------------------------------------------------------
assert n_sig == 982, n_sig
assert n_strong == 121, n_strong
assert len(new) == 15
assert len(reproduced) == 106
assert len(missing) == 0
drops = {("LAMC1"), ("TPD52"), ("SENP6"), ("HMGN3"), ("MT3"), ("RPL13"), ("ZBTB46"), ("ZNF100")}
dropped_symbols = {m25[(m25["gene"] == g) & (m25["outcome"] == o)]["symbol"].iloc[0]
                   for g, o in dropped}
assert dropped_symbols == drops, f"dropped set mismatch: {dropped_symbols ^ drops}"
print("\nALL ASSERTIONS PASSED ✓")
