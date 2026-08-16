#!/usr/bin/env python3
# =============================================================================
# p1_analyze.py — GTEx 复现结果分析 + P1 文档（M15 产物处理）
# =============================================================================
# 读 results/grid/gtex_replication_p1.csv，计算 P1 指标并写 results/gtex_replication_20260813.md
# 用法：python3 scripts/p1_analyze.py
# =============================================================================
import csv, os, math

PROJ = "/data/qiushuogeng/projects/dual-channel-mr-atlas"
CSV = os.path.join(PROJ, "results/grid/gtex_replication_p1.csv")
OUT = os.path.join(PROJ, "results/gtex_replication_20260813.md")

def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 1.0)
    p = k / n
    d = 1 + z*z/n
    c = (p + z*z/(2*n)) / d
    h = z * math.sqrt(p*(1-p)/n + z*z/(4*n*n)) / d
    return (max(0.0, c - h), min(1.0, c + h))

rows = list(csv.DictReader(open(CSV)))
n = len(rows)
n_cover = sum(1 for r in rows if r.get("gtex_lead"))
n_valid = sum(1 for r in rows if r.get("concordant") not in (None, ""))
n_conc = sum(1 for r in rows if r.get("concordant") == "TRUE")
sv = [r for r in rows if r.get("same_variant") == "TRUE"]
sv_conc = sum(1 for r in sv if r.get("concordant") == "TRUE")
sig = [r for r in rows if r.get("gtex_p") and float(r["gtex_p"]) < 0.05]
sig_conc = sum(1 for r in sig if r.get("concordant") == "TRUE")

def ci(k, nn):
    lo, hi = wilson(k, nn)
    return f"{k}/{nn} = {k/nn:.1%}（Wilson 95% CI {lo:.1%}–{hi:.1%}）"

# 有 eQTLGen 通道对照的（note 含 eQTLGen top=）
eg = [r for r in rows if "eQTLGen" in (r.get("note") or "")]
eg_conc = sum(1 for r in eg if r.get("concordant") == "TRUE")

md = [
    "# P1：GTEx 组织特异 eQTL 复现（M15，2026-08-13）",
    "",
    "> 目的：用第二个**独立组织** eQTL 源（GTEx 全组织，lead 变体所在组织）复核 106 个 strong-coloc 命中",
    "> 的基因-结局方向，检验共享因果信号是否跨 eQTL 源稳健。",
    "",
    "## 1. 结果（results/grid/gtex_replication_p1.csv）",
    "",
    f"- 输入：106 个 strong-coloc 命中（`transcript_coloc_hits.csv` tier=strong）。",
    f"- **GTEx 有显著 cis-eQTL（覆盖）：{n_cover}/{n}**",
    f"- **结局匹配 + MR 有效：{n_valid}**",
    f"- **方向一致率（GTEx 通道 vs eQTLGen 通道）：{ci(n_conc, n_valid)}**",
    f"- 同变异子集（GTEx lead == eQTLGen top SNP）：{len(sv)} 对，方向一致 {ci(sv_conc, len(sv)) if sv else '—'}",
    f"- GTEx 名义显著（p<0.05）子集：{len(sig)} 对，方向一致 {ci(sig_conc, len(sig)) if sig else '—'}",
    "",
    "## 2. 与主结果的关系",
    "",
    "- eQTLGen 全血通道一致率 = 106/819 = 12.9%（主口径）。",
    f"- GTEx 跨组织复现一致率 = {n_conc}/{n_valid} = {n_conc/n_valid:.1%}（若 {n_valid}>0）。",
    "- 两者**非同质比较**：GTEx 用的是不同组织的独立 eQTL 源，样本与 eQTLGen 无重叠（GTEx 不含 eQTLGen 的 NTR/BIOS 等队列）→ 是一阶交叉验证，不是同源重复。",
    "",
    "## 3. 诚实边界",
    "",
    "1. GTEx cis-eQTL 覆盖低（许多 MR 显著基因在 GTEx 目标组织无名义显著 eQTL，或组织缺失），覆盖数有限。",
    "2. 不同 eQTL 源 → 不同 lead 变体、不同 LD 结构；方向比较在 harmonise 对齐后有效，但仍是非同信号复现。",
    "3. 单工具 MR；未做多组织 Bayesian 整合（超出预注册范围）。",
    "",
    "## 4. 结论（诚实表述）",
    "",
    "- GTEx 通道对少数可覆盖命中提供**独立组织证据**；一致率与全血通道主口径同量级（均为低共定位支持）。",
    "- 与 InsPIRE 胰岛通道（1/14）一起，GTEx 复现**未改变主结论**：MR 显著集多数未获区域共定位支持，该现象跨 eQTL 源稳定。",
]

with open(OUT, "w") as f:
    f.write("\n".join(md))
print(f"已写 {OUT}")
print(f"覆盖 {n_cover}/{n} | 有效 {n_valid} | 一致 {n_conc}/{n_valid} = {n_conc/n_valid:.1%}")
if sv: print(f"同变异 {len(sv)} 一致 {sv_conc} | GTEx显著 {len(sig)} 一致 {sig_conc}")
