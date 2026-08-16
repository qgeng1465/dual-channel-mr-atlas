#!/usr/bin/env python3
"""蛋白通道 MR 首读统计量修正（2026-08-13，多 agent 评审后落地）。

背景：mr_ivw_mre 在 nsnp=2-3（工具间比率近同质）时 SE 被压缩 7-16 倍，
导致 INSR×FBG p=0、APOC3×FBG p=9.1e-46 等头条 p 值为伪影（见 journal 评审）。
修正口径：nsnp=1 → Wald ratio；nsnp<=3 → IVW-FE 首读；nsnp>=4 → IVW-MRE 首读。
同时补齐 F 统计量（deCODE 侧 F≈(b/se)^2，对首读方法；Wald 单独标注）。
本脚本只改报告口径，不动任何门柱（FE 本就是预注册敏感性方法）。
"""
import pandas as pd
import numpy as np
import os

HERE = os.path.dirname(os.path.abspath(__file__))
GRID = os.path.join(HERE, "..", "results", "grid")

MRE = "Inverse variance weighted (multiplicative random effects)"
FE = "Inverse variance weighted (fixed effects)"
WM = "Weighted median"
EGGER = "MR Egger"
WALD = "Wald ratio"

OUTCOME_SHORT = {
    "ebi-a-GCST006867": "T2D",
    "ebi-a-GCST005194": "CAD",
    "ebi-a-GCST005186": "FBG",
}

def primary_method(nsnp: int) -> str:
    """修正后首读方法：nsnp<=3 → FE（避免 MRE 小样本 SE 塌缩），nsnp>=4 → MRE。"""
    if nsnp == 1:
        return WALD
    if nsnp <= 3:
        return FE
    return MRE

def main():
    df = pd.read_csv(os.path.join(GRID, "protein_decode_mr.csv"))
    ok = df[df["ok"] == True].copy()

    rows = []
    for (gene, outcome), g in ok.groupby(["gene", "outcome"]):
        nsnp = int(g["nsnp"].iloc[0])
        pm = primary_method(nsnp)
        row = g[g["method"] == pm]
        if len(row) == 0:  # 例如 nsnp=1 的 Wald 行
            row = g.iloc[[0]]
        r = row.iloc[0]
        f_stat = (float(r["b"]) / float(r["se"])) ** 2 if pd.notna(r["b"]) and float(r["se"]) > 0 else np.nan
        rows.append({
            "gene": gene,
            "outcome": outcome,
            "outcome_short": OUTCOME_SHORT.get(outcome, outcome),
            "nsnp": nsnp,
            "primary_method": pm.replace("Inverse variance weighted ", "IVW-"),
            "b": r["b"], "se": r["se"], "pval": r["pval"],
            "F_stat": f_stat,
            "sig_fe": _sig(g, FE), "sig_mre": _sig(g, MRE), "sig_wm": _sig(g, WM),
            "sig_egger_slope": _sig_egger(g),
            "MRE_FE_SE_ratio": _se_ratio(g),
        })

    out = pd.DataFrame(rows)
    out = out.sort_values(["pval"], ascending=True).reset_index(drop=True)
    outpath = os.path.join(GRID, "protein_decode_mr_primary.csv")
    out.to_csv(outpath, index=False)

    print("修正后首读结果表（按 p 升序）:")
    print(out[["gene", "outcome_short", "nsnp", "primary_method", "b", "se", "pval", "F_stat",
               "sig_fe", "sig_mre", "sig_wm", "sig_egger_slope", "MRE_FE_SE_ratio"]]
          .to_string(index=False, float_format=lambda x: f"{x:.3g}"))

    # H1 计数准备：分泌 vs 胞内（描述性，诚实报告 n）
    print("\n--- 预注册 H1 计数（分泌 vs 胞内显著性率，修正后）---")
    secr = {"PCSK9", "APOC3", "ANGPTL3"}  # 分泌/循环
    intra = {"INSR", "PCK1"}               # 跨膜受体 / 胞内酶
    for grp, genes in (("分泌型", secr), ("胞内/膜", intra)):
        sub = out[out["gene"].isin(genes)]
        sig = sub[sub["pval"] < 0.05]
        print(f"{grp}: {len(sub)} 对, MR 显著 {len(sig)} ({len(sig)}/{len(sub)} = "
              f"{len(sig)/max(len(sub),1):.0%})")
        if len(sig):
            print("   显著对:", ", ".join(f"{r.gene}×{r.outcome_short}" for r in sig.itertuples()))

    print(f"\n已写出: {outpath}")
    return out

def _sig(g: pd.DataFrame, method: str) -> bool:
    r = g[g["method"] == method]
    return bool(len(r) and pd.notna(r["pval"].iloc[0]) and r["pval"].iloc[0] < 0.05)

def _sig_egger(g: pd.DataFrame) -> bool:
    r = g[g["method"] == EGGER]
    return bool(len(r) and pd.notna(r["pval"].iloc[0]) and r["pval"].iloc[0] < 0.05)

def _se_ratio(g: pd.DataFrame) -> float:
    m = g[g["method"] == MRE]; f = g[g["method"] == FE]
    if len(m) and len(f) and f["se"].iloc[0] > 0:
        return float(f["se"].iloc[0]) / float(m["se"].iloc[0])
    return np.nan

if __name__ == "__main__":
    main()
