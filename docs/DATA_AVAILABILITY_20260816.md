# Data and Code Availability — dual-channel MR atlas（2026-08-16）

> 本文件是投稿时的 **Data Availability 声明蓝本** + **发布打包清单**。所有承诺必须真实可兑现——凡标注【待打包】的项在投稿提交前必须完成，否则删除对应承诺。

---

## 1. 预注册定位（诚实声明，正文 Methods 引用）

**本分析为探索性全量枚举，未在第三方平台（OSF/AsPredicted）做时间戳前瞻性预注册。**
分析启动前已生成内部设计锁定文档 `docs/PREREGISTRATION.md`（工具变量标准、MR 方法集、
coloc 闸门、多重检验、正/负对照），作为分析纪律与透明化记录；因分析已完成（2026-08-15/16），
不再补作第三方预注册，亦不将本地锁定文档表述为正式预注册。

**投稿措辞（逐字可用）**：
> *"This study is a hypothesis-generating, exhaustive enumeration of a publicly available
> summary-statistics atlas. The analytical protocol (instrument definition, MR methods,
> colocalization threshold, multiple-testing control) was fixed in an internal design-lock
> document before analysis; the analysis was not prospectively registered on a third-party
> platform. All data and code are available (see Data Availability), enabling independent
> verification."*

---

## 2. 数据来源与版本（全部公开 sumstats，无个体数据）

| 来源 | 数据 | ID/URL | 获取日期 |
|---|---|---|---|
| eQTLGen | 全血 cis-eQTL，n=31,684 | eqtlgen.org | 2026-08 |
| T2D GWAS | Mahajan 2018 | GCST006867（n=655,666；ncase=61,714） | 2026-08 |
| CAD GWAS | van der Harst 2018 | GCST005194（n=296,525） | 2026-08 |
| FBG GWAS | MAGIC | GCST005186（n=58,074） | 2026-08 |
| FinnGen R11 | 结局侧独立复现 | FinnGen release 11 sumstats | 2026-08 |
| GTEx v8 | eQTL 方向复现 | gtexportal.org | 2026-08 |
| 1000 Genomes EUR | LD 参考面板 | Phase 3 (hg19) | 2026-08 |

**坐标 build 声明**：eQTLGen 与 1000G 为 hg19；OpenGWAS 结局元数据 hm_pos 为 hg38。所有
coloc 采用**纯 rsID 匹配**（build 无关），规避跨数据源坐标错配；早期窗口错配版本已存档并
内部披露（`results/archive/buildbug_20260815/`）。

### 2b. 数据溯源明细（本地文件与完整性；合并自 `docs/archive/data_sources.md`）

| 数据 | 版本 | 样本量 | 来源 URL | 下载日期 | 本地文件 | 完整性 |
|---|---|---|---|---|---|---|
| eQTLGen 全血 cis-eQTL（显著） | Võsa 2021 | n=31,684 | https://download.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz | 2026-08-06 | data/eqtlgen/cis-eQTL-significant.txt.gz | ✅ GZIP_OK, 10.5M 行, 16,923 基因 |
| eQTLGen SNP 频率（等位基因定向） | Võsa 2021 | — | https://download.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz | 2026-08-06 | data/eqtlgen/SNP_AF.txt.gz | 下载中（240MB，服务器带宽受限） |
| T2D GWAS（Xue 2018） | GCST006867 | n=655,666（ncase=61,714） | OpenGWAS ebi-a-GCST006867 / GWAS Catalog FTP | 2026-08 | data/opengwas/full/t2d_full.gz | ✅ API 解析成功 |
| CAD GWAS | GCST005194 | n=296,525 | OpenGWAS ebi-a-GCST005194 | 2026-08 | data/opengwas/full/cad_full.gz | ✅ API 解析成功 |
| FBG GWAS | GCST005186 | n=58,074 | OpenGWAS ebi-a-GCST005186 | 2026-08 | data/opengwas/full/fbg_full.gz | ✅ API 解析成功 |
| FinnGen R11 | Release 11 | — | FinnGen sumstats（结局侧复现） | 2026-08 | data/finngen/ | ✅ 已核验 |
| GTEx v8 | GTEx 2020 | n=838 | gtexportal.org | 2026-08 | data/gtex/ | ✅ 已核验 |
| 1000 Genomes EUR | Phase 3 | n=503 | 1000G 官网 | 2026-08 | data/ldref/1kg.v3/EUR | ✅ LD 参考面板 |

> **诚实标注**：`docs/archive/data_sources.md` 原件中"暴露层通道 2（循环蛋白 deCODE/INTERVAL
> pQTL）"条目属于论文早期**双通道设计**阶段的产物，最终论文已重定位于转录组 cis-MR × coloc
> 图谱（去双通道），蛋白通道数据**未进入论文分析**，故未合并进本表，仅随原件归档。

---

## 3. 代码可用性

- 全部分析脚本（M20–M30 + 早期网格）存放 `scripts/`，含版本化注释与参数。
- 软件版本：coloc 5.2.3 / R env r-mr / PLINK 1.9 / SMR 1.3.1 / python 3.10+（numpy/pandas/matplotlib）。
- 随机种子：master = 20260805，下游确定性派生（`docs/PREREGISTRATION.md` §8）。

---

## 4. 数据可用性（图谱资源，发布打包清单）

### 4.1 核心交付（投稿时必须可下载，拿 DOI；2026-08-16 已按 FDR-core 更新）
- [x] **FDR-core 显著集**：`results/fdr_core_20260816.csv`（982 对：gene/symbol/outcome/mr_b/mr_p/padj/gwas_min_p/eqtl_F_max/nsnp/pp4/strong）
- [x] **15 候选复现表**：`results/candidate15_replication_20260816.csv`（tier/PP.H4/mr_p/padj/GTEx direction/FinnGen variant+gene/p）
- [x] **全量图谱表**：`results/coloc_full_{t2d,cad,fbg}_20260815.csv`（31,373 对全字段：gene/symbol/outcome/mr_b/mr_p/gwas_min_p/eqtl_F_max/nsnp/pp4/ok/note）
- [x] **操作特性漏斗**：`results/m36b_funnel_20260816.csv`（名义曲线 + FDR-core 12.32% + grid 12.96%）+ `results/m36b_summary_20260816.csv`
- [x] **coloc.susie 敏感性**：`results/m34_coloc_susie_20260816.csv` + 评估文档 `docs/COLOC_SUSIE_ASSESSMENT_20260816.md`
- [x] **全量 strong 子集**：`results/strong_all_subset_20260816.csv`（121 FDR-core strong + 2 灰区 AP3S2×T2D/ZNF19×CAD，set/status 标注；106 known 见 `results/grid/transcript_coloc_hits.csv`）
- [x] **schema 文档**：`docs/SCHEMA_20260816.md`（全部结果表列字典 + 口径定义 + build + 阈值 + 缺失标记）
- [x] **README 复现指南**：`README.md` 摘要 + 新增「复现指南（Reproducibility）」流水线表（数据源→脚本→产物，M20–M36b）

### 4.2 发布渠道（按 REVIEWER_IF_ADVICE 分级）
1. **Zenodo**（免费，期刊普遍要求，拿 DOI）→ 全量表 + 脚本 + schema + README 打包上传（**需作者账号**）
2. **GitHub**（公开仓库，脚本 + 结果 + 图 + manuscript）——2026-08-16 已 push 至 qgeng1465/dual-channel-mr-atlas（**默认私有，需作者一键公开**）
3. OSF 项目页（可选，作为数据/代码链接聚合页，**非预注册**）

### 4.3 版本号
图谱资源以 **v1.0** 发布；后续修正递增并写 CHANGELOG，保证"v1.0"成为可引用句柄。

---

## 5. 伦理与合规

- 全部输入为**公开汇总统计（GWAS/eQTL sumstats）**，不涉及个体级数据或人类受试者，**无需伦理批准**；投稿声明照此填写。
- 无利益冲突（无资助方介入分析/发表）。
- 作者贡献按 CRediT 列出（投稿时填写）。

---

## 6. 学术端正自查清单（投稿前逐项打钩）

| 项 | 状态 | 落点 |
|---|---|---|
| 无捏造/伪造/抄袭 | ✅ | INTEGRITY_AUDIT_20260816.md 正面重算成立 |
| **FDR-core 主口径（982/121/15/12.32%）vs raw 名义漏斗并排披露** | ✅ | FACTS §1b/§2 + manuscript §3.1/§3.2 |
| 15 候选 = 候选评估非因果发现（非新 GWAS 位点） | ✅ | manuscript §3.5 + FACTS §4 |
| "catalog 未报道"≠"文献未报道" | ✅ | 全文措辞纪律 |
| 41/106=38.7% GWAS 峰显著 caveat 披露 | ✅ | manuscript §3.5 + Fig 9 |
| 覆盖 9/15=60% FinnGen 对齐披露 | ✅ | manuscript §3.6 |
| coloc.susie 外样本 LD 不收敛（exploratory）披露 | ✅ | COLOC_SUSIE_ASSESSMENT + manuscript §3.7 |
| 样本重叠（Burgess 2016）作为首要局限 | ✅ | manuscript §2.2 + §4 |
| 全网格落盘（含失败/空结果） | ✅ | 31,373 对含 note 列 |
| p12 主口径 vs 敏感性、Wilson CI | ✅ | FACTS §2 |
| **Data/Code 可复现（DOI + schema）** | ⏳【待打包】 | 本节 4.1 |
| **预注册定位诚实（无第三方前瞻性注册）** | ✅ 本文档 | §1 |
| 软件版本/种子 pinning | ✅ | §3 |
| 利益冲突 / CRediT | ✅ 无冲突 | 投稿表 |
