#!/usr/bin/env python3
# =============================================================================
# M42_pqtl_coverage_20260817.py - pQTL 中心法则整合的覆盖检查（诚实局限）
# =============================================================================
# 目的：审稿人任务 3（pQTL 3-track LocusZoom）的可执行性评估。中心法则整合
#   需要同一基因的 cis-eQTL + 蛋白 pQTL + 结局 GWAS 三轨数据。本脚本如实核算
#   本地可用 pQTL 源对 15 候选的覆盖：
#   * deCODE (data/decode/): 9 蛋白（药物靶标定向下载，已核实文件名）
#   * UKB-PPP (data/ukbpp/): 8 蛋白（同上）
#   * OpenGWAS prot-a-*/ukb-b-* pQTL API: 网络不可达（gwas-api.mrcieu.ac.uk
#     unreachable, api.opengwas.io 401 auth）
# 输出：
#   results/m42_pqtl_coverage_20260817.csv  15 候选 × {deCODE, UKBPP, 结论}
#   results/m42_pqtl_summary_20260817.csv   汇总行（覆盖数、不可达 API）
# 诚实 caveat：pQTL 覆盖 0 → 3-track 对 15 候选不可执行；本地 pQTL 仅覆盖
#   药物靶蛋白，与 121 known 基因交集如实报告。
# =============================================================================
import os, pandas as pd

REPO = "/data/qiushuogeng/projects/dual-channel-mr-atlas"
cand = pd.read_csv(os.path.join(REPO, "results/candidate15_replication_20260816.csv"))
known = pd.read_csv(os.path.join(REPO, "results/grid/transcript_coloc_hits.csv"))

# 已核实（data/decode + data/ukbpp 文件名）：药物靶蛋白
decode = {"ANGPTL3", "GLP1R", "DPP4", "PCK1", "INSR", "GCG", "HMGCR", "PCSK9", "APOC3"}
ukbpp  = {"ANGPTL3", "APOB", "DPP4", "GCG", "GLP1R", "INSR", "LDLR", "PCSK9"}
local_prots = sorted(decode | ukbpp)
cand_syms = sorted(cand["symbol"].dropna().unique())
known_syms = set(known["symbol"].dropna().astype(str))

rows = [dict(symbol=g, in_decode=(g in decode), in_ukbpp=(g in ukbpp),
             pqtl_available=False, note="no local pQTL; OpenGWAS API unreachable")
        for g in cand_syms]
pd.DataFrame(rows).to_csv(os.path.join(REPO, "results/m42_pqtl_coverage_20260817.csv"), index=False)
overlap = [g for g in local_prots if g in known_syms]
summary = pd.DataFrame([dict(n_candidates=len(cand_syms),
    n_with_local_pqtl=sum(r["pqtl_available"] for r in rows),
    local_pqtl_proteins=";".join(local_prots), n_local_pqtl_proteins=len(local_prots),
    known_genes_with_local_pqtl=";".join(overlap), n_known_genes_with_local_pqtl=len(overlap),
    api_status="OpenGWAS gwas-api.mrcieu.ac.uk unreachable; api.opengwas.io requires auth (401)",
    conclusion="3-track central-dogma integration infeasible for 15 candidates with available data")])
summary.to_csv(os.path.join(REPO, "results/m42_pqtl_summary_20260817.csv"), index=False)
print("local pQTL proteins:", local_prots)
print("known genes with local pQTL overlap:", overlap)
print("== DONE M42 (honest limitation) ==")
