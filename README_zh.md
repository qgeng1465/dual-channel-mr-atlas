# 2 型糖尿病、冠心病与空腹血糖的全转录组 cis-MR × 共定位图谱

[English](README.md) | **简体中文**

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

一个可复现的全转录组 **cis-MR × 贝叶斯共定位** 图谱，覆盖 **2 型糖尿病（T2D）**、**冠心病（CAD）** 与 **空腹血糖（FG）** 三个结局。本项目量化了全转录组 cis 孟德尔随机化（cis-MR）作为基因优先级筛选工具的**操作特性（operating characteristics）**，提供一张经校准的共定位图谱作为公共资源，并识别出具备正交证据（MR + coloc + 方向复现）支撑的候选效应基因。

> **关于分析版本。** 分结局 BH-FDR 流水线（"FDR-core"，2026-08-16 定稿）是权威分析，也是论文所采用的分析。所有面向投稿的数字见 [`docs/manuscript/manuscript.md`](docs/manuscript/manuscript.md) 与 `results/` 中的结果表。

---

## 摘要

全基因组关联研究（GWAS）已为 2 型糖尿病（T2D）、冠心病（CAD）和空腹血糖（FG）鉴定出数百个位点，但多数效应基因仍属未知。全转录组孟德尔随机化（cis-MR）联合贝叶斯共定位被广泛用于基因优先级排序，但其操作特性尚未被系统量化。我们对 eQTLGen 全血 cis-eQTL（n = 31,684）与三个 GWAS 结局在 31,371 个基因–性状对上进行 cis-MR × coloc 扫描。在分结局 FDR 控制（q < 0.05）下，982 对达到 MR 显著，其中 121 对达到强共定位（PP.H4 ≥ 0.8）：即 12.3% 的共定位产出率。在阈值敏感性分析中，产出率自名义 p < 0.05 时的 3.0% 单调升至 p < 1e-5 时的 25.6%，且 MR 显著集之外强共定位几乎不存在（2/27,123 对）。全部 106 个既往报道位点均被重现，并新鉴定出 15 个候选效应基因，在 GTEx 与 FinnGen 中获得方向性复现（4 个在 FinnGen 中达到名义显著）。我们如实披露校准局限，并提供一张经校准的共定位支持公共图谱。

---

## 核心发现

| 发现 | 数值 |
|---|---|
| 扫描的基因–性状对（QC 通过） | 31,371（原始 31,373） |
| 分结局 BH-FDR 显著（q < 0.05）MR 对 | **982**（T2D 394 / CAD 576 / FG 12） |
| 强共定位（PP.H4 ≥ 0.8） | **121**（T2D 65 / CAD 54 / FG 2）→ **共定位产出率 12.3%** |
| 重现的既往报道位点 | **106 / 106** |
| 新候选效应基因 | **15**（9 个位于已知位点 + 6 个无 GWAS Catalog 记录） |
| 产出率校准（名义 MR p 阈值） | p < 0.05 时 3.0% → p < 1e-5 时 25.6%（单调） |
| 与 stage-2 网格扫描的一致性 | 12.3%（FDR-core） vs 12.96%（grid） |
| MR 显著集之外的强共定位 | 2 / 27,123（决定性负边界） |

**候选效应基因（15 个）：** SLC12A3、CWF19L1、U6atac、CD101、RBM6、CNNM2、N4BP2L2、RIC8A、C2orf49（位于已知 T2D/CAD 风险位点内）；PLAUR、TAGLN2、VSIG8、PDCD6、CLEC3B、CCDC19（无 T2D/CAD GWAS Catalog 记录的区域）。

**独立复现（15 候选口径）：**
- **GTEx v8** eQTL 方向：6/7 一致（如实报告 1 个冲突：VSIG8）。
- **FinnGen R11**（独立队列）：可定位子集内 9/9 基因级方向一致（8/9 变异级），4 个在 FinnGen 中达到名义显著（RBM6、CNNM2、CD101、RIC8A）；对齐覆盖 9/15 = 60%。

**方法学操作特性：**
- SMR + HEIDI 一致性：76/106 = 71.7%。
- Steiger 方向：73/76 = 96.1%。
- 零假设下置换假阳性率：1.45%（154 / 10,600）。

**诚实性 caveat（论文中如实披露）：**
- 41/106 = 38.7% 的已知强共定位区域 GWAS 峰 p < 5e-8 → 多数共定位区域不构成新的 GWAS 位点，被解读为候选效应基因。
- `coloc.susie` 在外样本 LD 估计下未完全收敛（探索性）；LAMC1 同时依据 FDR 与多信号证据被排除。
- 候选效应基因是**候选评估**而非因果发现。

---

## 论文稿件

- [`docs/manuscript/manuscript.md`](docs/manuscript/manuscript.md) — 完整稿件源文件（标题页、摘要、IMRaD、表格、图、参考文献）。
- [`docs/manuscript/AJHG_submission_Qiushuo_Geng_20260816.docx`](docs/manuscript/AJHG_submission_Qiushuo_Geng_20260816.docx) — 格式化投稿文档。
- [`docs/manuscript/cover_letter.md`](docs/manuscript/cover_letter.md) — 投稿信。

---

## 仓库结构

```
scripts/   全部分析脚本（M20–M28、M33–M38；见"复现指南"）
results/   全部结果表（CSV）、汇总统计与图
  coloc_full_{t2d,cad,fbg}_20260815.csv   全量扫描输出（31,373 对）
  fdr_core_20260816.csv                   FDR-core MR 显著集（982）
  strong_all_subset_20260816.csv          121 个强共定位 + 2 个灰区对
  candidate15_replication_20260816.csv    15 个候选及复现统计
  grid/                                    stage-2 网格扫描输出
  figures/                                 论文图
docs/      论文
LICENSE    CC BY 4.0
```

---

## 复现指南

完整流水线（输入 → 脚本 → 输出）见各分析脚本本身。关键步骤：

| 步骤 | 脚本 | 输入 | 输出 |
|---|---|---|---|
| 全转录组 cis-MR × coloc 全量扫描 | `M20*`–`M24` | eQTLGen + 三个 GWAS（rsID 匹配，hg19） | `results/coloc_full_{t2d,cad,fbg}_20260815.csv` |
| 新强共定位候选发现 | `M25` / `M25b` | `coloc_full_*` | `results/m25_new_strong_annotation_20260816.csv` |
| GTEx v8 独立方向复现 | `M26` | GTEx v8 | `results/m26_gtex_replication_new23_20260816.csv` |
| 名义精度漏斗 | `M27` | `coloc_full_*` | `results/m27_precision_funnel_20260816.csv` |
| FinnGen R11 独立队列复现 | `M28` | FinnGen R11 sumstats | `results/m28_finngen_replication_new23_20260816.csv` |
| **分结局 BH-FDR 重算（FDR-core）** | `M36b_fdr_recompute_20260816.py` | `coloc_full_*` + grid | `results/fdr_core_20260816.csv`、`candidate15_replication_20260816.csv`、`m36b_funnel_20260816.csv` |
| coloc.susie 敏感性（探索性） | `M34b` | 6 位点 + 1000G EUR LD | `results/m34_coloc_susie_20260816.csv` |
| 图（5 主图 + S1） | `M37` / `M38` | 各结果表 | `results/figures/20260816_Fig{1..5}_*.png`、`FigS1_susie.png` |
| Word 稿件 | `M36_build_word_ajhg_20260816.py` | `docs/manuscript/*` | `docs/manuscript/AJHG_submission_Qiushuo_Geng_20260816.docx` |

**路径占位符。** 脚本使用可移植占位符而非机器特定路径：`<repo-root>`（仓库根目录）、`<scratch>`（下载的源数据暂存目录）、`<conda-root>`（R 环境前缀）。重跑前请替换为你本地的实际路径；批量源数据（eQTLGen、OpenGWAS、GTEx、FinnGen、1000 Genomes）需另行下载。

---

## 数据可用性

**输入数据（公开）：**
- eQTLGen 全血 cis-eQTL（n = 31,684；Võsa et al., 2021）。
- GWAS：T2D（`ebi-a-GCST006867`）、CAD（`ebi-a-GCST005194`）、FG（`ebi-a-GCST005186`），经 OpenGWAS。
- GTEx v8、FinnGen R11、1000 Genomes Phase 3（EUR LD 参考）。

**输出数据：** 全部分析输出（结果表、汇总统计、图、脚本）均在本仓库中提供。

**存档副本：** 本仓库的版本化存档副本将存放于 Zenodo。

> DOI 待 Zenodo 归档后填入。

---

## 许可

本仓库基于 [Creative Commons Attribution 4.0 International (CC BY 4.0)](LICENSE) 许可发布。

## 引用

Qiushuo Geng. *Transcriptome-wide cis-MR and colocalization atlas for type 2 diabetes, coronary artery disease, and fasting glucose: operating characteristics and candidate effector genes.* Zenodo: DOI 待分配。
