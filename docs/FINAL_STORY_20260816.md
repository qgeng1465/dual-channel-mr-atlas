# 最终故事：我们能诚实讲一个什么样的正结果（2026-08-16，FDR-core 定稿）

> 目的：回答"现在能将一个怎么样的故事就直接讲了"。本文档把现有全部证据按可投稿口径汇总，
> 明确**能说的**与**不能说的**，列出把故事做实的可选补实验，并直接回答
> "这次改动还有意义吗？之后还要做改动吗？想发更高分"。
> 配套：审计 `docs/INTEGRITY_AUDIT_20260816.md`、写作方案 `docs/PAPER_WRITING_PLAN_20260816.md`、
> IF 建议 `docs/REVIEWER_IF_ADVICE_20260816.md`、图 `results/figures/20260816_F*_v2.png`。

> **⚠️ 口径声明（2026-08-16）**：投稿主口径 = **分结局 BH-FDR q<0.05（FDR-core）**。
> 本文档早期快照中的 raw 口径数字（131/129/23/11/12/12/21）已被 FDR-core 取代，仅作过程参照；
> **唯一数字权威** = `docs/FACTS_20260816.md` + `docs/manuscript/manuscript.md`。

---

## 〇、一句话判定

**可以讲一个"图谱资源 + 方法学操作特性 + 候选效应基因发现"三合一的诚实正结果**：
31,373 对全转录组 cis-MR×coloc 全量扫描，**FDR-core 982 对 MR-significant → 121 个 strong coloc**
（PP.H4≥0.8；另有 2 个灰区 AP3S2×T2D / ZNF19×CAD 在 MR-sig 集外；106 个已知位点 100% 重现），
**新发现 15 个 GWAS Catalog 未报道的候选效应基因**，其中 9 个落在已知 T2D/CAD 风险位点 ±100–250 kb 内
（新效应基因提名）、6 个在 catalog 无记录区域（弱候选）；独立 GTEx v8 eQTL 方向 6/7 一致复现、
FinnGen R11 独立队列基因级 9/9 方向一致（4 个 p<0.05）；方法学上首次给出单工具 cis-MR 的
**操作特性**（MR p 阈值单调校准 strong 率 0.71%→25.6%；FDR-core 12.32% 与 stage-2 网格 12.96% 收敛）。
这比"警示性阴性结论"（08-13 叙事）强得多，且**全部数字经审计重算成立**。

---

## 一、故事的三股线（投稿叙事）

### 线 1 —— 图谱资源（public resource）
全转录组 cis-MR×coloc 图谱：eQTLGen（n=31,684）→ T2D/CAD/FBG 三个结局，**31,373 对**（QC 31,371），
**FDR-core 982 对 MR-significant → 121 个 strong coloc（PP.H4≥0.8）**（另有 2 个灰区 MR-null）。
这是一张药物靶点 MR 社区可直接复用的证据底座表（全量 123 strong 子集 = `results/strong_all_subset_20260816.csv`）。
（注：蛋白通道 deCODE 仅 PCSK9/APOC3 两对，已降级为补充材料，不占主叙事。）

### 线 2 —— 方法学：单工具 cis-MR 的操作特性（methodology）
审稿人关心的核心问题：**"筛出来的有多少是真"**（精度/yield）与**"真的有多少被筛出"**（召回）从未被系统量化。
本项目给了答案（M24/M27/M36b）：
- **MR 召回好，不丢信号**：FDR-core MR-significant 集内 121 个 strong 全部命中，MR-sig 集外 strong 仅
  2 个灰区（AP3S2×T2D、ZNF19×CAD，均 MR-null；AP3S2 eQTL F=55.6/MR p=0.38/PP.H4=0.935、ZNF19 eQTL F=7.13
  弱 eQTL/MR p=0.33/PP.H4=0.903 → "真实位点上单工具 MR 功效不足"，决定维度在 GWAS 侧）。
- **但 yield 低**：名义口径下全量 MR 显著集（4,248 对，raw p<0.05）strong 率仅 **3.04%**；应用预注册的
  **分结局 BH-FDR（q<0.05）** 把核心集收敛到 **982 对 → 12.32% yield**（Wilson 95% CI 10.41–14.53%），
  与标准 stage-2 复核网格（clump+IVW，818 可评估对 → **12.96%**，富集 4.3×）**收敛于 ~12–13%**——两种
  独立过滤策略（FDR 控制 vs LD clump + 工具复现）给出同一共享因果变异率。
- **精度漏斗**（M27）：收紧 MR p 阈值，strong 率**单调校准** 0.71%→25.6%（`20260816_F3_yield_v2.png`）→
  给出了一个**可操作的工具**："用多严的 MR p 阈值，期望多少 coloc 支持率"，直接给全领域"MR 筛药靶"实践贴了一张操作特性曲线。

### 线 3 —— 发现：15 个候选效应基因（discovery）
FDR-core 全量扫描比 stage-2 预过滤网格**多出 15 个 strong coloc**（M25/M25b + M36b，不在已知 106、
不在 GWAS Catalog T2D/CAD）：
- **9 个落在已知 T2D/CAD 风险位点 ±100–250kb 内**（基因本体未被 T2D/CAD 报道）→ **候选效应基因提名**：
  SLC12A3×CAD、CWF19L1×T2D、U6atac×T2D、CD101×T2D、RBM6×T2D、CNNM2×CAD、N4BP2L2×CAD、RIC8A×CAD（100 kb）、C2orf49×T2D（250 kb）
- **6 个所在区域无 T2D/CAD catalog 关联**（PLAUR、TAGLN2、VSIG8、PDCD6、CLEC3B、CCDC19），GWAS 峰多 >5e-8 → **只算低置信假设，不算新位点**
- **8 个 raw strong 被 FDR 掉出**（LAMC1、TPD52、SENP6、HMGN3、MT3、RPL13、ZBTB46、ZNF100，padj≥0.05），**不作为候选**；LAMC1 另经 coloc.susie 多信号证据（8 eQTL CS vs 3 GWAS CS）独立排除
- **独立复现（15 候选口径）**：
  - GTEx v8（M26）：7 个可测基因中 **6 个方向一致**（1 个冲突：VSIG8）；8 个不可测（GTEx lead 不在结局 GWAS）
  - FinnGen R11 独立队列（M28，全量 sumstats）：9 个可检（6 个 lead 未在 FinnGen R11 定位），**基因级 9/9、变异级 8/9 方向一致，0 冲突**；其中 **4 个在 FinnGen 自身达 p<0.05**（RBM6 6.9e-06、CNNM2 6.3e-04、CD101 1.2e-03、RIC8A 1.1e-02）；**诚实披露：对齐覆盖 9/15（60%）**
- **强调**：这是**候选评估**，不是发现宣称；"catalog 未报道"仅指 GWAS Catalog，不排除其他文献。

---

## 二、证据链数字一览（全部经审计重算，FDR-core 主口径）

| 数字 | 值 | 出处 |
|---|---|---|
| 全量对数 | 31,373（QC 31,371） | M24 `coloc_full_summary` |
| MR 显著集（raw p<0.05） | 4,248 对（仅作名义漏斗） | M24 |
| **FDR-core MR-significant** | **982 对**（t2d 394/cad 576/fbg 12） | M36b `fdr_core` |
| **strong coloc（FDR-core）** | **121**（t2d 65/cad 54/fbg 2） | M36b |
| 灰区 strong（MR-null） | 2（AP3S2×T2D、ZNF19×CAD） | M24/M36b |
| 已知 106 重现 | 106/106（100%） | M36b |
| **新增候选（strong ∩ FDR-core − 106）** | **15**（9 已知位点候选 + 6 无 catalog 弱候选） | M25/M25b/M36b |
| FDR 掉出的 raw strong | 8（LAMC1/TPD52/SENP6/HMGN3/MT3/RPL13/ZBTB46/ZNF100） | M36b |
| **FDR-core yield** | **12.32%** [Wilson 10.41–14.53%] | M36b |
| stage-2 网格 yield | **12.96%**（106/818）[10.83–15.43%] | M5/M27 |
| 名义漏斗（mr_p<0.05→<1e-5） | 3.04% → 6.56% → 8.67% → 14.81% → 17.49% → 24.59% → 25.59% | M27 |
| MR-sig 集外 strong | 2/27,123（AP3S2、ZNF19，单侧 95% 上界 0.0232%/对） | M24/M36b |
| null（mr_p≥0.5）strong | 0 | M24 |
| 已知集 HEIDI 通过 | 76/106=71.7% | M22/M26 |
| 已知集 Steiger 通过 | 73/76=96.1% | M22 |
| 已知集 GWAS 峰 p<5e-8 | 41/106=38.7% | M24 |
| GTEx 独立方向（15 候选） | 6 一致 / 1 冲突 / 8 不可测 | M26 |
| FinnGen R11（9/15 可检） | 基因级 9/9、变异级 8/9 方向一致，0 冲突；FinnGen 自身 p<0.05：**4 个**；对齐覆盖 **9/15（60%）** | M28 |
| coloc.susie（外样本 LD，exploratory） | 5/6 稳健（SuSiE PP.H4≥0.999），LAMC1 排除 | M34/M34b |

---

## 三、诚实边界（这些话投稿时绝不能写）

1. **15 个是候选，不是因果发现**。未经 SMR+HEIDI、独立队列全效应复现、或功能验证前，正文措辞只能是
   "candidate effector genes" / "hypothesis-generating"，不是 "new causal genes"。
2. **"catalog 未报道" ≠ "文献未报道"**。GWAS Catalog 是关联条目库，不是文献综述。9 个在已知位点的基因
   只能说"该位点已知，但此基因未在 Catalog T2D/CAD 关联中作为标注基因"，不能写"全新基因-疾病关联"。
3. **6 个弱峰候选**（无 catalog 记录、GWAS 峰多 >5e-8）只能作为低置信列表，最多进补充材料，不进主表头条。
4. **coloc 在弱 GWAS 峰下的行为不可靠**（fbg 2 个 strong 峰均不显著）。正文必须写：41/106（38.7%）
   已知强候选的区域 GWAS 峰 p<5e-8，"MR 显著集内 strong coloc 大部分位于未达全基因组显著的区域"
   （审计 P0-3/P1-3）。
5. **方向复现 ≠ 效应量复现**。GTEx/FinnGen 只比符号，不比幅度。
6. **VSIG8 方向冲突如实报告**（M26 已记），不隐藏。
7. **无家族级显著性声明**：121 是 FDR-core 内 strong 调用；置换标定（FP=1.45%）只覆盖 raw MR 显著子集，
   未做全量尺度置换（审计 P0-3）。正文写"图谱为描述性/假设生成性质"。
8. **coloc.susie 外样本 LD 不收敛**（max_iter=200；max_iter=1000 未解决）→ 仅作 exploratory 敏感性，
   LAMC1 排除依据多信号证据，不依赖非收敛后验。

---

## 四、把故事做实的可选补实验（用户允许"补实验"）

按**性价比从高到低**：

| 可选项 | 成本 | 收益 | 判定 |
|---|---|---|---|
| **FinnGen 独立队列 MR 方向复现**（本机已下 R11 sumstats，M28 跑） | ~30 min | 结局侧独立复现，候选列表可信度大增 | ✅ 已完成（9/15 可检，4 个 p<0.05） |
| GTEx 全组织而非 6 组织方向复现 | 已有 6 组织；扩到 49 组织 + coloc.susie 条件分析 | 更强的 eQTL 侧证据 | 可选，成本低 |
| **SMR+HEIDI 候选层** | 需重建 .esd+.ma | 独立方法学验证 | ⚠️ 已评估：15 候选 trans 探针覆盖不足，跳过，写进 Limitation（已知集 HEIDI 76/106 作为参考） |
| Finngen 在线 API 复现 | API 已 404（证书过期），已改用直下全量 | — | 已完成替代 |
| **TWAS（FUSION/PrediXcan）交叉验证** | 需 LDREF+权重（已有 1kg.v3 EUR 8.8G） | 另一种独立方法学 | 可选，若时间允许 |
| 湿实验（细胞/组织表达敲低验证 1-3 个 top 候选，如 PLAUR/CNNM2 在巨噬细胞/胰岛） | 数千元–数万 + 数月 | 把候选坐实为"可发表的功能证据" | 超出纯干实验范围，作为 future work 写在 Discussion |
| 在线公开数据外延：UK Biobank 表型-基因型关联（M13 已有 ukbpp 覆盖） | 已有部分 | 备选 | 视情况 |

**结论**：纯干实验范围内，分析已收官；其余属锦上添花或已评估不可行。
湿实验是 Discussion "future work" 而非必做项。

---

## 五、对三个直接问题的回答

### Q1：这次做的改动还有意义吗？
**有意义，且是决定性的。** 08-13 叙事是"无新因果发现"的阴性结论（106/818=12.96% 精度，无候选）。
M24 全量扫描 + M25–M28 + M36b FDR-core 把它改成**阳性**：15 个候选 + 操作特性漏斗 + 双独立复现。
对一篇目标 AJHG 的文章，"阴性+方法学"与"阳性发现+方法学+资源"的接受概率天差地别。
**这些改动直接决定了这篇论文是发得出去还是发不出去。**

### Q2：之后还要做改动吗？
**分析层面：不再开新战场。** 完成即止。理由：
(1) 15 个候选的下一步（SMR/功能验证）超出纯干实验数据可得性；
(2) 操作特性漏斗已把方法学故事说满；(3) 继续加分析 → 边际收益递减、过度编辑风险上升。
**收尾层面：已完成** —— 审计 P0/P1 文档修复、README 同步 v0.13、图定稿（9 主图 + S1 v2）、commit。

### Q3：发到哪？怎么发？
> **2026-08-16 定案：定投 AJHG（IF 7.7，Article 类型）**。直接对标样稿 Ray et al. 2025
> "Single-cell transcriptome-wide MR and colocalization... ASCVD"（AJHG 112(7):1597–1609，
> DOI 10.1016/j.ajhg.2025.06.001）——方法学骨架几乎一致（cis-MR + coloc + drug target），
> 方向完全在 AJHG scope 内。投稿准备指南见 `docs/AJHG_SUBMISSION_GUIDE_20260816.md`。
> 其他期刊（eBioMedicine/Diabetologia/HMG/NC）的推测性定位已全部删除。

- **与 Ray et al. 的差异化（Cover Letter 卖点）**：Ray 回答"因果基因作用在哪个细胞类型"（单细胞分辨）；
  我们回答"cis-MR 筛选工具本身有多可靠"（操作特性标定）——互补而非竞争。我们的 31,373 对图谱可供
  单细胞研究（如 Ray）直接查询。
- **差异化卖点**：① 全转录组尺度（31,373 对，同类多为细胞/组织子集）② **操作特性曲线（精度漏斗）是
  任何同类都没有的方法学原创** ③ 跨 3 结局泛化 ④ GTEx+FinnGen 双独立复现。
- **诚实提醒**：AJHG 竞争者靠单细胞分辨或功能验证撑新意（§二）。我们不追这两点，用"全转录组底表 +
  操作特性曲线"讲差异化；候选一律"hypothesis-generating"，不做因果/细胞特异宣称。
- **不要做的提档动作**：不要为数字好看放宽口径（如把 6 个弱候选吹成新位点、把"候选"写成"因果"）——
  审计已盯死这类动作，被审稿人抓到就是撤稿级风险。诚实提档的唯一路径是补证据（FinnGen 已做）而非改措辞。

---

## 六、交付物清单（2026-08-16 定稿）

1. **定稿稿件**：`docs/manuscript/manuscript.md` + `docs/manuscript/AJHG_submission_Qiushuo_Geng_20260816.docx`
2. **数字权威**：`docs/FACTS_20260816.md`；列字典 `docs/SCHEMA_20260816.md`
3. **核心结果**：`results/fdr_core_20260816.csv`（982）、`results/candidate15_replication_20260816.csv`（15）、
   `results/strong_all_subset_20260816.csv`（123）、`results/m36b_funnel_20260816.csv`、`results/m36b_summary_20260816.csv`
4. **图**：`results/figures/20260816_F{1..9}_v2.png` + `S1_resources_v2.png`
5. **数据可用性**：`docs/DATA_AVAILABILITY_20260816.md`（GitHub qgeng1465/dual-channel-mr-atlas + Zenodo DOI 待注册）
6. `docs/INTEGRITY_AUDIT_20260816.md` —— 审计（P0/P1 修复状态见 §七）

## 七、审计 P0/P1/P2 修复状态（2026-08-16 定稿快照）

| 项 | 状态 |
|---|---|
| P0-1 双精度并排（3.04% vs 12.96%） | ✅ 名义漏斗（M27）+ FDR-core 12.32% 并排（M36b 漏斗 CSV + manuscript §3.1 + F3 图） |
| P0-2 README/故事同步 M20-M36b | ✅ README v0.13（摘要、复现指南、§0、§4.7、§6、§8 全同步）+ 本文档 |
| P0-3 弱 GWAS caveat + "多数应峰显著"错误 | ✅ 披露 41/106=38.7%；FACTS/manuscript §3.5/Fig 9 |
| P1-1 候选注释 | ✅ M25/M25b 已注释；FDR-core 15 候选明细见 `candidate15_replication_20260816.csv` |
| P1-2 T2D N | ✅ 全量文件实测 per-variant N mode≈573,704（闭合） |
| P1-3 弱峰 caveat 对 t2d/cad | ✅ 并入 manuscript §3.5 |
| P1-4 README FinnGen"未修"过期描述 | ✅ §5 局限 7 已改"已修" |
| P1-5 多重检验 | ✅ 分结局 BH-FDR 为预注册主口径（M36b），raw p<0.05 仅作名义漏斗 |
| P2-1 灰区 0→1 | ✅ 灰区 2（AP3S2、ZNF19）如实标注 |
| P2-2 ≥0.5 仅 2 | ✅ 负边界 0 strong 于 mr_p≥0.5 |
| P2-3 旧 55%/45% | ✅ verification §10b 已加修正注 |
| P2-4 56+20 | ✅ verification §6 + 表已加修正注 |
