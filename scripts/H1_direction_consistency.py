#!/usr/bin/env python3
"""预注册 H1：分泌型 vs 胞内/膜蛋白的 pQTL-MR 显著性率 + 双通道方向一致性（2026-08-13）。

数据：results/grid/protein_decode_mr_primary.csv（修正后首读）+ transcript_drugtarget_mr.csv。
诚实口径：工具级双通道重叠基因极少 → 方向一致性为描述性报告，如实声明 n。
"""
import pandas as pd, numpy as np, os
from scipy.stats import fisher_exact

HERE = os.path.dirname(os.path.abspath(__file__))
GRID = os.path.join(HERE, "..", "results", "grid")
OUT_SHORT = {"ebi-a-GCST006867": "T2D", "ebi-a-GCST005194": "CAD", "ebi-a-GCST005186": "FBG"}

SECRETED = {"PCSK9", "APOC3", "ANGPTL3"}   # 分泌/循环蛋白
INTRA = {"INSR", "PCK1"}                    # 跨膜受体 / 胞内酶

def primary(df, gene, outcome):
    g = df[(df.gene == gene) & (df.outcome == outcome)].sort_values("pval")
    # 选首读方法行（p 最小的就是首读，因同一对 FE/MRE/Wald b 相同）
    return g.iloc[0] if len(g) else None

def main():
    prim = pd.read_csv(os.path.join(GRID, "protein_decode_mr_primary.csv"))
    tr = pd.read_csv(os.path.join(GRID, "transcript_drugtarget_mr.csv"))

    # ---- 1) 分泌 vs 胞内：pQTL-MR 显著性率 ----
    print("=" * 70)
    print("H1 分支一：分泌型 vs 胞内/膜蛋白 的 pQTL-MR 显著率（修正后首读 p<0.05）")
    print("=" * 70)
    for grp, genes in (("分泌型", SECRETED), ("胞内/膜", INTRA)):
        sub = prim[prim.gene.isin(genes)]
        n_sig = (sub.pval < 0.05).sum()
        print(f"  {grp} {sorted(genes)}: {n_sig}/{len(sub)} = {n_sig/len(sub):.0%}")
    table = [[int((prim[prim.gene.isin(SECRETED)].pval < 0.05).sum()),
              int((prim[prim.gene.isin(SECRETED)].pval >= 0.05).sum())],
             [int((prim[prim.gene.isin(INTRA)].pval < 0.05).sum()),
              int((prim[prim.gene.isin(INTRA)].pval >= 0.05).sum())]]
    or_, p = fisher_exact(table)
    print(f"  Fisher 2×2 (分泌显著/不显著, 胞内显著/不显著): {table}")
    print(f"  OR={or_:.2f}, p={p:.3f}  → 功效严重不足，仅描述性报告（预注册纪律：标 CI、不排序）")

    # ---- 2) 双通道方向一致性 ----
    print()
    print("=" * 70)
    print("H1 分支二：双通道方向一致性（转录 vs 蛋白，同一基因×同一结局）")
    print("=" * 70)
    pairs = []
    for _, r in prim.iterrows():
        g, o = r.gene, r.outcome
        tp = tr[(tr.gene == g) & (tr.outcome == o) & tr.ok].sort_values("pval")
        if len(tp) == 0:
            continue
        t = tp.iloc[0]
        both_sig = (r.pval < 0.05) and (t.pval < 0.05)
        pair = {
            "gene": g, "outcome": OUT_SHORT.get(o, o),
            "protein_b": round(r.b, 4), "protein_p": r.pval,
            "transcript_b": round(t.b, 4), "transcript_p": t.pval,
            "both_sig": both_sig,
            "direction_consistent": both_sig and (np.sign(r.b) == np.sign(t.b)),
        }
        pairs.append(pair)
    if pairs:
        pd.DataFrame(pairs).to_csv(os.path.join(GRID, "H1_direction_consistency.csv"), index=False)
        for p in pairs:
            dc = "一致" if p["direction_consistent"] else ("—" if not p["both_sig"] else "相反")
            print(f"  {p['gene']}×{p['outcome']}: 蛋白 b={p['protein_b']} p={p['protein_p']:.3g}"
                  f" | 转录 b={p['transcript_b']} p={p['transcript_p']:.3g}"
                  f" | 双显著={'是' if p['both_sig'] else '否'} | 方向{dc}")
        n_dual = sum(p["both_sig"] for p in pairs)
        print(f"\n  工具级双通道重叠对（两通道均有 MR）：{len(pairs)}；其中双通道均显著的：{n_dual}")
        if n_dual:
            cons = sum(p["direction_consistent"] for p in pairs)
            print(f"  方向一致率: {cons}/{n_dual} = {cons/n_dual:.0%}（n 极小，仅描述性）")
        else:
            print("  → 双通道均显著的对为 0，方向一致性无法计算（功效不足），如实声明 n=0")
    else:
        print("  无工具级双通道重叠基因（仅 INSR 有工具级重叠，且其转录显著对×蛋白显著对落在不同结局）")
    print("\n已写出: results/grid/H1_direction_consistency.csv")

if __name__ == "__main__":
    main()
