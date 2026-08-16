#!/usr/bin/env python3
# =============================================================================
# M21_ukbpp_coverage.py — strong coloc 编码基因 × UKB-PPP Olink 面板覆盖核查
# =============================================================================
# 目的（BREAKTHROUGH_PLAN §6 方案 B 验证）：判定 UKB-PPP 蛋白层扩展值不值得做。
#   核 strong coloc 命中里的编码基因（biotype 判定，~76）在 UKB-PPP Olink
#   Explore（European discovery, syn51365303, 2940 蛋白 tar）里被测量几个。
#   ≥25 → UKB-PPP 路线可行；<20 → 覆盖不足，改走 deCODE 定向或放弃蛋白层扩展。
# 2026-08-15 改版：弃 S3 list（ukb-pqtl bucket 不存在），改用 Synapse children
#   （synapseclient 列 syn51365303 的 2940 个 `Gene_UniProt_OID_v#_Panel.tar`）。
# 用法：HTTPS_PROXY=http://127.0.0.1:7890 python3 scripts/M21_ukbpp_coverage.py
# 输出：results/ukbpp_coverage_20260815.csv + 终端摘要
# =============================================================================
import csv, os, re
import synapseclient

PROJ = "/data/qiushuogeng/projects/dual-channel-mr-atlas"
RES  = os.path.join(PROJ, "results")

# 1. strong coloc 命中 + biotype
hits = {}
with open(os.path.join(RES, "grid/transcript_coloc_hits.csv")) as f:
    for row in csv.DictReader(f):
        hits[row["symbol"]] = row   # 106 行 → 105 unique symbol（1 个跨结局重复）
bt = {}
with open(os.path.join(RES, "grid/ensg_biotype_20260813.csv")) as f:
    for row in csv.DictReader(f):
        bt[row["ENSG"]] = row.get("biotype", "")
coding = {s for s, h in hits.items() if bt.get(h.get("gene", ""), "") == "protein_coding"}
# 兜底：若 ENSG 全对不上，用全部 strong 命中当编码
if len(coding) < 50:
    coding = set(hits.keys())
print(f"strong coloc 命中: {len(hits)} unique symbol | 判定编码基因: {len(coding)}")

# 2. UKB-PPP European discovery 蛋白目录（Synapse children）
syn = synapseclient.Synapse()
syn.login(silent=True)
print("列 UKB-PPP European discovery (syn51365303) children...")
kids = list(syn.getChildren("syn51365303"))
gene_from_dir = set()
for k in kids:
    n = k["name"]                       # `Gene_UniProt_OID_v#_Panel.tar`
    g = n.split("_")[0]
    if g:
        gene_from_dir.add(g)
print(f"UKB-PPP 蛋白 tar: {len(kids)} | 可解析 gene 名: {len(gene_from_dir)}")

# 3. 覆盖
overlap = sorted(coding & gene_from_dir)
print(f"\n编码基因 ∩ UKB-PPP 面板: {len(overlap)}/{len(coding)}")
print("命中:", ", ".join(overlap))
missing = sorted(coding - gene_from_dir)
print(f"未命中 {len(missing)}: ", ", ".join(missing))

with open(os.path.join(RES, "ukbpp_coverage_20260815.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["coding_gene", "in_ukbpp", "outcome", "pp_h4"])
    for g in sorted(coding):
        h = hits[g]
        w.writerow([g, g in gene_from_dir, h.get("outcome", ""), h.get("PP.H4", "")])

print("\n判定: ",
      "UKB-PPP 路线可行（≥25）" if len(overlap) >= 25 else
      ("边缘（20-24）" if len(overlap) >= 20 else
       "UKB-PPP 覆盖不足（<20），改走 deCODE 定向或放弃蛋白层扩展"))
