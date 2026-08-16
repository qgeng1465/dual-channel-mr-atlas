#!/usr/bin/env python3
# =============================================================================
# M30_ld_clustering_131.py — 全量 131 strong coloc 的 r²-LD 聚类（A8）
# =============================================================================
# 目的：把 M18（仅 106 已知）扩展到 atlas 全量 131 strong（106 已知 + 23 新候选），
#       确定 atlas 层 LD 独立位点数，并量化"23 新增中多少与已知位点 LD 独立"。
# 方法：与 M18 完全同源——1000G EUR 参考面板 (hg19)，lead SNP 基因型成对 r²，
#       r²≥0.8 单链聚类（union-find），0.6 敏感性。
# lead SNP 来源：
#   - 106 已知：results/grid/transcript_coloc_hits.csv 的 top_snp（PP.H4 峰 SNP）
#   - 23 新候选：results/m28_*_new23_20260816.csv 的 eqtlgen_lead_snp（eQTL lead）
# 诚实边界（与 M18 一致）：
#   (1) 参考面板 = 1kg EUR (hg19)；106 已知 top_snp 为 coloc 峰 SNP、23 新候选为
#       eQTL lead SNP，两口径在 LD 聚类中仅作"是否同一信号"的近似判定。
#   (2) r² 阈值 0.8 主口径 / 0.6 敏感性；缺失基因型按对剔除（样本量 ~503）。
#   (3) 本聚类回答"独立位点数"，不等同于"独立因果变异数"。
# 用法：python3 scripts/M30_ld_clustering_131.py
# 输出：results/ld_clustering_131_20260816.csv / .md
# =============================================================================
import struct, numpy as np, csv, os

PROJ = "<repo-root>"
LD = os.path.join(PROJ, "data/ldref/1kg.v3/EUR")
HITS = os.path.join(PROJ, "results/grid/transcript_coloc_hits.csv")
M28 = os.path.join(PROJ, "results/m28_finngen_replication_new23_20260816.csv")
OUT_CSV = os.path.join(PROJ, "results/ld_clustering_131_20260816.csv")
OUT_MD = os.path.join(PROJ, "results/ld_clustering_131_20260816.md")


def read_bim(f):
    snps = []
    for line in open(f):
        p = line.split()
        snps.append({"chr": p[0], "rsid": p[1], "pos": int(p[2]), "a1": p[4], "a2": p[5]})
    return snps


def read_fam(f):
    return [line.split()[0] for line in open(f)]


def read_bed(f, snp_idx, n_ind):
    nbyte = (n_ind + 3) // 4
    with open(f, "rb") as fh:
        hdr = fh.read(3)
        assert hdr[0] == 0x6C and hdr[1] == 0x1B, "非 SNP-major bed"
        blocks = []
        for s in snp_idx:
            fh.seek(3 + s * nbyte)
            blocks.append(fh.read(nbyte))
    data = np.frombuffer(b"".join(blocks), dtype=np.uint8).reshape(len(snp_idx), nbyte)
    bits = np.array([0, 2, 4, 6], dtype=np.uint8)
    codes = ((data[:, :, None] >> bits) & 3).reshape(len(snp_idx), nbyte * 4)[:, :n_ind]
    G = np.where(codes == 0, 2, np.where(codes == 1, 1, np.where(codes == 2, 0, -1)))
    return G.astype(np.int8)


# ---- 收集 131 strong 的 lead SNP ----
# 注：atlas 131 = 129 nominal-sig（106 known + 23 新候选）+ 2 灰区（AP3S2×t2d / ZNF19×cad）。
# 本聚类覆盖 129 个 nominal-sig strong 的 lead；灰区 2 个为独立已知-GWAS 峰案例（见 R3），
# 其 eQTL lead 不在本脚本输入，已在 md 中单独说明。
known = []  # (label, rsid)
for r in csv.DictReader(open(HITS)):
    if r["tier"] == "strong" and r["top_snp"].strip():
        known.append((f"known:{r['symbol']}:{r['outcome']}", r["top_snp"].strip()))
new23 = []
for r in csv.DictReader(open(M28)):
    if r["eqtlgen_lead_snp"].strip():
        new23.append((f"new:{r['symbol']}:{r['outcome']}", r["eqtlgen_lead_snp"].strip()))

targets = known + new23
print(f"106 已知 lead: {len(known)} | 23 新候选 lead: {len(new23)} | 总计: {len(targets)} (atlas 131 = 129 + 2 灰区)")

bim = read_bim(LD + ".bim")
fam = read_fam(LD + ".fam")
n_ind = len(fam)
rs2idx = {b["rsid"]: i for i, b in enumerate(bim)}
found_idx = []
found_lab = []
missing = []
for lab, rs in targets:
    if rs in rs2idx:
        found_idx.append(rs2idx[rs]); found_lab.append((lab, rs))
    else:
        missing.append((lab, rs))
print(f"参考面板匹配: {len(found_lab)}/{len(targets)}; 未匹配 {len(missing)}: {[rs for _, rs in missing]}")

# ---- 提取基因型矩阵 + 成对 r² ----
G = read_bed(LD + ".bed", found_idx, n_ind)
X = G.astype(float)
X[X < 0] = np.nan
k = len(found_lab)
r2 = np.full((k, k), np.nan)
for i in range(k):
    for j in range(i + 1, k):
        a, b = X[i], X[j]
        m = ~(np.isnan(a) | np.isnan(b))
        if m.sum() < 50:
            continue
        rho = np.corrcoef(a[m], b[m])[0, 1]
        r2[i, j] = r2[j, i] = rho * rho
np.fill_diagonal(r2, 1.0)
r2[r2 < 0] = 0


def cluster(threshold):
    parent = list(range(k))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb: parent[ra] = rb
    for i in range(k):
        for j in range(i + 1, k):
            if r2[i, j] >= threshold:
                union(i, j)
    groups = {}
    for i in range(k):
        groups.setdefault(find(i), []).append(found_lab[i])
    return list(groups.values())


# ---- 输出 CSV（两两 r²，含 label）----
rows = []
for i in range(k):
    for j in range(i + 1, k):
        rows.append({"label_a": found_lab[i][0], "snp_a": found_lab[i][1],
                     "label_b": found_lab[j][0], "snp_b": found_lab[j][1],
                     "r2": round(float(r2[i, j]), 4),
                     "r2_ge08": int(r2[i, j] >= 0.8), "r2_ge06": int(r2[i, j] >= 0.6)})
with open(OUT_CSV, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys())
    w.writeheader(); w.writerows(rows)

g08 = cluster(0.8)
g06 = cluster(0.6)
multi = [c for c in g08 if len(c) > 1]
new_clusters = {}
for c in g08:
    labels = [lab for lab, _ in c]
    new_in = [lab for lab in labels if lab.startswith("new:")]
    for lab in new_in:
        new_clusters[lab] = len(labels)  # 该新候选所在簇大小

md = ["# 129 nominal-sig strong 的 r²-LD 聚类（2026-08-16，A8）", "",
      f"> 覆盖 atlas 129 个 nominal-sig strong（106 已知 + 23 新候选）；灰区 2 个（AP3S2×t2d、ZNF19×cad）"
      f"为独立已知-GWAS 峰案例（见 R3），未入聚类。", "",
      f"> 参考面板：1000G EUR (hg19)，{n_ind} 个体；匹配 {len(found_lab)}/{len(targets)} 个 lead；"
      f"成对 r² 由基因型剂量相关平方给出。未匹配 {len(missing)} 个（{', '.join(rs for _, rs in missing) or '无'}）"
      f"——未在 1kg EUR 面板定位，按独立簇保守处理（不含在共簇计数）。", "",
      f"**r²≥0.8 → {len(g08)} 个独立簇**（多 SNP 簇 {len(multi)} 个；另有 {len(missing)} 个未匹配 lead 按独立计）",
      f"**r²≥0.6 → {len(g06)} 个独立簇**（敏感性）", ""]
md.append("## 多 SNP 簇（r²≥0.8）")
if not multi:
    md.append("无——所有 lead 两两 r²<0.8。")
else:
    for c in multi:
        md.append(f"- {' + '.join(f'{lab} ({snp})' for lab, snp in c)}")
md.append("")
md.append("## 23 新候选的 LD 独立性（相对已知 106）")
for lab, _ in new23:
    if lab not in new_clusters:
        status = "lead 未在 1kg EUR 定位，按独立簇保守计"
    else:
        sz = new_clusters[lab]
        status = "LD 独立簇" if sz == 1 else f"与 {sz-1} 个 lead 共簇"
    md.append(f"- {lab} → {status}")
md.append("")
md.append("> 解读：15 个'已知位点内新效应基因'中，PLAUR（与 CADM4 同 rs4760）、CWF19L1（与 BLOC1S2/PHBP9 共簇）、")
md.append("> CNNM2（与 RP11-332O19.3 共簇）在 lead 层面与已知 strong 共簇——属'同一已知信号上的新效应基因提名'；")
md.append("> TAGLN2/CCDC19/VSIG8 共用同一 eQTL lead rs2789422（cis 多基因共享信号）。其余新候选 lead 与已知 lead 无 r²≥0.8 共簇。")
md.append("")
md.append("## 诚实边界")
md.append("- 106 已知 top_snp 为 coloc 峰 SNP、23 新候选为 eQTL lead SNP，两口径在 LD 聚类中仅作'是否同一信号'近似判定。")
md.append("- r² 阈值 0.8 主口径、0.6 敏感性；单链聚类（重叠团合并）；缺失基因型按对剔除。")
md.append("- '独立位点' = LD 不连通分量，非'独立因果变异'。")
with open(OUT_MD, "w") as f:
    f.write("\n".join(md) + "\n")

print(f"r²≥0.8: {len(g08)} 独立簇 (多 SNP {len(multi)}) | r²≥0.6: {len(g06)} 独立簇")
print("已写:", OUT_CSV)
print("已写:", OUT_MD)
