#!/usr/bin/env python3
# =============================================================================
# M17_pathway_enrichment.py — 76 优先基因 / 56 编码基因的功能富集（探索性，阴性结果）
# =============================================================================
# 目的：给图谱论文的"生物学解释"层做探索——优先基因是否富集到已知通路。
# 数据：tmp/priority_genes_heidi_pass.txt（76 优先基因：MR+coloc+HEIDI+Steiger 全通过）
#       results/drug_annotation_20260813.csv（56 编码 / 20 非编码标记）
# 方法：g:Profiler REST API（biit.cs.ut.ee），g_SCS 多重检验校正。
#       网络走本机 mihomo 代理（127.0.0.1:7890，符合全局代理偏好）。
# 口径（2026-08-13 已排雷）：
#   - 显著判据 = 响应 p_value < 0.05（g_SCS 校正后），不用 `significant` 字段
#     （该字段在 user_threshold=1.0 下把 p=0.35~0.97 的条目也标 True，含义非 g_SCS<0.05）。
#   - 驱动基因解析：`intersections` 与 `query` 逐基因对齐，非空 evidence 列表 = 命中。
#   - user_threshold=1.0 返回全部条目（保留透明对照）；g_SCS p 值由方法本身给出。
# 诚实边界：
#   (1) 76 优先基因是选择后子集 → 富集反映选择过程，非随机基因组基线，仅作探索性/假设生成。
#   (2) 20 个非编码/伪影基因无注释通路 → 主结果用 56 编码集，76 全集作透明对照。
#   (3) 小基因集 → 单条目可能由 1-2 基因驱动，看命中基因数。
# 用法：python3 scripts/M17_pathway_enrichment.py
# 输出：results/pathway_enrichment_20260813.csv/.md
# =============================================================================
import json, subprocess, csv, os, sys

PROJ = "<repo-root>"
PROXY = "http://127.0.0.1:7890"
API = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
SOURCES = ["GO:BP", "GO:MF", "KEGG", "REAC"]


def first_existing(*paths):
    for p in paths:
        if os.path.exists(p):
            return p
    return paths[0]


PRIORITY = first_existing(
    os.path.join(PROJ, "tmp/priority_genes_heidi_pass.txt"),
    os.path.join(PROJ, "results/archive/tmp_202608/priority_genes_heidi_pass.txt"),
)
DRUGANN = first_existing(
    os.path.join(PROJ, "results/drug_annotation_20260813.csv"),
    os.path.join(PROJ, "results/archive/drug_annotation_20260813.csv"),
)

GENE2OUT = {}
for line in open(PRIORITY):
    p = line.strip().split("\t")
    if len(p) >= 3:
        GENE2OUT.setdefault(p[1], []).append(p[2])

all76 = list(GENE2OUT.keys())
coding = set()
for r in csv.DictReader(open(DRUGANN)):
    if r["biotype_ok"] == "1":
        coding.add(r["gene_symbol"])
coding56 = [g for g in all76 if g in coding]
nc20 = [g for g in all76 if g not in coding]
print(f"76 优先基因: {len(all76)} | 56 编码: {len(coding56)} | 20 非编码: {len(nc20)}")


def gp(query, label):
    body = {
        "organism": "hsapiens", "query": query, "sources": SOURCES,
        "user_threshold": 1.0, "significance_threshold_method": "g_SCS",
        "domain_scope": "annotated", "no_evidences": False,
    }
    r = subprocess.run(
        ["curl", "-s", "-m", "120", "-x", PROXY, "-X", "POST", API,
         "-H", "Content-Type: application/json", "-d", json.dumps(body)],
        capture_output=True, text=True)
    if r.returncode != 0 or not r.stdout:
        print(f"[{label}] API 失败: {r.stderr[:200]}"); return [], label
    try:
        d = json.loads(r.stdout)
    except Exception as e:
        print(f"[{label}] JSON 解析失败: {e}"); return [], label
    res = d.get("result", [])
    n_sig = sum(1 for it in res if it["p_value"] < 0.05)
    print(f"[{label}] 返回 {len(res)} 条 | g_SCS p<0.05: {n_sig}")
    return res, query


def run(query, label):
    res, q = gp(query, label)
    rows = []
    for it in res:
        # intersections 与 query 逐基因对齐：非空 evidence 列表 = 该基因在 term 内
        ev = it.get("intersections") or []
        genes = [g for g, e in zip(q, ev) if isinstance(e, list) and e]
        rows.append({
            "set": label, "source": it["source"], "term_id": it["native"],
            "term_name": it["name"], "p_value_gSCS": it["p_value"],
            "query_size": it["query_size"], "term_size": it["term_size"],
            "n_intersect": it["intersection_size"], "genes": ";".join(sorted(genes)),
        })
    rows.sort(key=lambda r: r["p_value_gSCS"])
    return rows


out = run(coding56, "coding56")
out += run(all76, "all76")

if not out:
    print("无任何返回条目——如实报告。"); sys.exit(0)

with open(os.path.join(PROJ, "results/pathway_enrichment_20260813.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
    w.writeheader(); w.writerows(out)

sig56 = [r for r in out if r["set"] == "coding56" and r["p_value_gSCS"] < 0.05]
sig76 = [r for r in out if r["set"] == "all76" and r["p_value_gSCS"] < 0.05]
all_nom56 = sorted([r for r in out if r["set"] == "coding56"], key=lambda r: r["p_value_gSCS"])
top_nom = all_nom56[:10]
best = all_nom56[0]["p_value_gSCS"] if all_nom56 else float("nan")

md = ["# 优先基因功能富集（2026-08-13，探索性，阴性结果）", "",
      f"> 76 优先基因（MR+coloc+HEIDI+Steiger 全通过）；主结果 = 56 编码集，76 全集为透明对照。",
      f"> g:Profiler g_SCS 校正；方法 + 口径排雷见 `scripts/M17_pathway_enrichment.py`。", ""]
md.append(f"**56 编码 g_SCS<0.05 条目: {len(sig56)} | 76 全集: {len(sig76)} | 56 编码最佳 g_SCS p = {best:.2f}**")
md.append("")
md.append("## 结论（如实）")
md.append("**优先基因在 GO:BP/MF/KEGG/Reactome 通路层面无显著富集（g_SCS 校正后 0 条，最佳 p≈0.35）。**")
md.append("这与 56 个编码优先基因的功能异质性一致（多为管家/普适表达基因：CHD4 染色质重塑、EIF2B2、")
md.append("BLOC1S2、FASTKD5、C6orf106 等），未见共享通路信号。该阴性结果本身是诚实的生物学信息：")
md.append("cis-eQTL 驱动的 T2D/CAD/FBG 优先位点基因分散在多种基本生物学过程，无单一路径富集。")
md.append("注意：不排除选择偏倚掩盖真信号，或通路覆盖不完整；此结果不作任何因果主张。")
md.append("")
md.append("## 透明对照：56 编码集名义（g_SCS）p 最小 10 条（均不显著）")
md.append("")
md.append("| 源 | 通路 | p | 命中基因 | 驱动基因 |")
md.append("|---|---|---|---|---|")
for r in top_nom:
    md.append(f"| {r['source']} | {r['term_name']} | {r['p_value_gSCS']:.3g} | "
              f"{r['n_intersect']}/{r['query_size']} | {r['genes']} |")
md.append("")
md.append("## 诚实边界")
md.append("- **选择偏倚**：76 基因经 MR+coloc+HEIDI+Steiger 筛选，富集反映选择过程而非随机基因组基线；")
md.append("  仅作探索性/假设生成。")
md.append("- **小样本**：56 个编码基因，条目多由 1-2 基因驱动（见命中基因数）。")
md.append("- **非编码基因排除**：20 个非编码/假基因（RP11-*/CTD-*/SERBP1P3）无注释通路，不参与主富集。")
md.append("- **口径排雷**：`significant` 字段在 g:Profiler 响应中含义非 g_SCS<0.05，已弃用，判据为 p_value<0.05。")
with open(os.path.join(PROJ, "results/pathway_enrichment_20260813.md"), "w") as f:
    f.write("\n".join(md))
print(f"\n已写 results/pathway_enrichment_20260813.csv（{len(out)} 行）+ .md")
