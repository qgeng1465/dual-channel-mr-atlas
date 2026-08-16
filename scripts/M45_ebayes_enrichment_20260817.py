#!/usr/bin/env python3
# =============================================================================
# M45_ebayes_enrichment_20260817.py — 经验贝叶斯富集先验校准（TORUS/fastENLOC 式）
# =============================================================================
# 目的：回应审稿人任务 2 —— 固定 p12=1e-5 的先验是任意的；用数据的全局共定位背景率
#   校准一个经验富集先验（等价 TORUS 的 π_d 富集参数 / fastENLOC 的区域共定位概率 RCP）：
#   * 全局背景率 π0 = P(H4) 来自全部测试位点（coloc_full_{t2d,cad} 的 pp4 均值）
#   * MR 显著集富集率 π_MR = 121/982（strong 共定位占比）
#   * 富集先验优势比 ρ = π_MR/(1-π_MR) ÷ π0/(1-π0)   （log-odds 尺度，Wilson CI）
#   * 每位点"富集标定 PP.H4（RCP 等价）"：pp4_enr = ρ·pp4/(1+(ρ-1)·pp4)
#     —— 由默认先验下的后验 pp4 与先验优势比缩放解析推出（见下推导）
#   * 检验 106 known + 15 候选在富集先验下是否仍 ≥0.8（稳健性，非固定 p12 伪影）
# 推导：默认先验 posterior pp4 = p12·BF4 / S，S=Σ p_i·BF_i。
#   富集先验 p12' = ρ·p12（其余先验不变）→ S' = S + (ρ-1)p12·BF4 = S[1+(ρ-1)pp4]
#   → pp4_enr = ρ·p12·BF4/S' = ρ·pp4/(1+(ρ-1)·pp4)。仅依赖 pp4 与 ρ，不需重跑 coloc。
# 输入：results/coloc_full_{t2d,cad}_20260815.csv（全部测试位点 + pp4）
#       results/grid/transcript_coloc_hits.csv（106 known）
#       results/candidate15_replication_20260816.csv（15 候选）
# 输出：results/m45_ebayes_20260817.csv + 摘要
# 诚实 caveat：这是"富集先验敏感性校准"而非完整 TORUS 逐 SNP 注释模型；ρ 为全局-显著集
#   边际富集（经验贝叶斯），非基因特异注释先验。方向结论仍为关联性评估。
# =============================================================================
import pandas as pd, numpy as np
from statsmodels.stats.proportion import proportion_confint

REPO = "<repo-root>"

# ---- 1. 全局背景率（两结局合并） ----
t2d = pd.read_csv(f"{REPO}/results/coloc_full_t2d_20260815.csv")
cad = pd.read_csv(f"{REPO}/results/coloc_full_cad_20260815.csv")
t2d = t2d[t2d["ok"] == True]; cad = cad[cad["ok"] == True]
allc = pd.concat([t2d, cad])
n_tot = len(allc)
n_strong_global = int((allc["pp4"] >= 0.8).sum())
pi0 = n_strong_global / n_tot

# ---- 2. MR 显著集富集率（982 sig / 121 strong） ----
n_mr, n_strong_mr = 982, 121
pi_mr = n_strong_mr / n_mr
# 富集先验优势比
logit = lambda p: np.log(p / (1 - p))
rho = np.exp(logit(pi_mr) - logit(pi0))
# Wilson CI on rho（对数尺度 ±1.96*se；se 用 delta 法）
se_logit0 = 1 / np.sqrt(n_tot * pi0 * (1 - pi0))
se_logit_mr = 1 / np.sqrt(n_mr * pi_mr * (1 - pi_mr))
se_logrho = np.sqrt(se_logit0**2 + se_logit_mr**2)
rho_lo, rho_hi = np.exp(np.log(rho) - 1.96 * se_logrho), np.exp(np.log(rho) + 1.96 * se_logrho)
# 直接 Wilson CI（比例差方向，供报告）
ci0 = proportion_confint(n_strong_global, n_tot, method="wilson")
cimr = proportion_confint(n_strong_mr, n_mr, method="wilson")

# ---- 3. 富集标定 pp4（RCP 等价） ----
def enrich(pp4, rho):
    return rho * pp4 / (1 + (rho - 1) * pp4)

res = []
for name, df in [("known", None), ("candidate", None)]:
    if name == "known":
        known = pd.read_csv(f"{REPO}/results/grid/transcript_coloc_hits.csv")
        # 用 known 表自带 PP.H4（其与 coloc_full 的 outcome 行一一对应），不另做 merge
        df = known[["gene", "symbol", "outcome", "PP.H4", "tier", "ok"]].copy()
        df = df.rename(columns={"PP.H4": "pp4"})
    else:
        cand = pd.read_csv(f"{REPO}/results/candidate15_replication_20260816.csv")
        df = cand[["gene", "symbol", "outcome", "mr_b", "mr_p", "pp4", "is_new", "tier"]].copy()
    df["pp4_enr"] = df["pp4"].apply(lambda x: enrich(x, rho))
    df["set"] = name
    res.append(df)

df = pd.concat(res, ignore_index=True)
# 富集先验下仍 ≥0.8 计数（稳健性）
n_stay = int((df["pp4_enr"] >= 0.8).sum())
# 汇总行
summ = pd.DataFrame([{
    "n_tested_all": n_tot, "n_strong_global": n_strong_global, "pi0_global": pi0,
    "n_mr_sig": n_mr, "n_strong_mr": n_strong_mr, "pi_mr": pi_mr,
    "rho_enrichment": rho, "rho_ci": f"[{rho_lo:.3g},{rho_hi:.3g}]",
    "pi0_wilson": f"[{ci0[0]:.4f},{ci0[1]:.4f}]", "pi_mr_wilson": f"[{cimr[0]:.4f},{cimr[1]:.4f}]",
    "n_known": int((df["set"] == "known").sum()),
    "n_cand": int((df["set"] == "candidate").sum()),
    "n_pp4enr_ge08": n_stay,
}])
summ.to_csv(f"{REPO}/results/m45_ebayes_summary_20260817.csv", index=False)
df.to_csv(f"{REPO}/results/m45_ebayes_20260817.csv", index=False)

print(f"global: n={n_tot} strong={n_strong_global} pi0={pi0:.4f} [{ci0[0]:.4f},{ci0[1]:.4f}]")
print(f"mr-sig: n={n_mr} strong={n_strong_mr} pi_mr={pi_mr:.4f} [{cimr[0]:.4f},{cimr[1]:.4f}]")
print(f"enrichment prior odds rho={rho:.3g} 95% CI [{rho_lo:.3g},{rho_hi:.3g}]")
print(f"calibrated pp4>=0.8 retained: {n_stay}/{len(df)}  (known+candidates)")
print(f"  -> known {df[df.set=='known'].pp4_enr.ge(0.8).sum()}/{sum(df.set=='known')}, "
      f"candidates {df[df.set=='candidate'].pp4_enr.ge(0.8).sum()}/{sum(df.set=='candidate')}")
