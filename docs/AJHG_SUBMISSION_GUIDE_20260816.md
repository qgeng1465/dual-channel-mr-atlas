# AJHG 投稿准备指南（2026-08-16）

> 定投 AJHG（American Journal of Human Genetics，IF 7.7），Article 类型。
> 依据：用户提供的 AJHG 官方投稿要求 + Ray et al. 2025 对标样稿 + 两个后台调研 agent 结果
> （投稿导师调研 / figure 风格调研，均于 2026-08-16 完成，已整合进本文档）。
> 目标：所有投稿文件直接按本文档产出。
> ⚠️ 调研局限性：WebFetch 对 cell.com/pmc/x-mol/medrxiv/scispace 等域名被网络策略拦截，
> 结论基于 WebSearch 聚合的索引摘要/图注/全文片段交叉验证。**提交前以 cell.com/ajhg/authors
> 复核一遍硬指标**。
>
> **⚠️ 口径声明（2026-08-16 定稿）**：投稿主口径 = **分结局 BH-FDR q<0.05（FDR-core：982 MR-sig → 121 strong → 15 候选）**。
> 本文档正文若干图表建议仍写早期 raw 口径（131/129/23/11/12/12/21），已全部被 FDR-core 取代；
> **数字一律以 `docs/FACTS_20260816.md` 与 `docs/manuscript/manuscript.md` 为准**，下述建议只取"形式/结构"，不取数字。

---

## 一、AJHG 硬性要求（官方 Author Guidelines，两个调研 agent 交叉印证）

| 项 | 要求 | 我们的执行 |
|---|---|---|
| 摘要 | **非结构化（单段、无小标题、禁引用）**；上限调研口径：用户指南 **≤200 词**，投稿导师调研 2026 模板 **≤250 词**（旧版 150 词已过时） | 目标 ≤200 词（更保守，见 §三 ~195 词模板） |
| 正文字数 | Article **≤60,000 字符 ≈ 9,000 词**（含参考文献） | 估算目标 ~7,500 词 |
| 图表 | **≤7 个**（主图+主表合并计数，Reports ≤4） | 精选方案见 §四 |
| 章节结构 | **Introduction → Material and Methods → Results → Discussion**（传统 IMRaD，**不用** Cell 主刊 STAR Methods/Lead Contact） | 按 IMRaD 组织 |
| Limitations | Discussion 内写 **"Limitations of the study" 小节**（Ray et al. 2025 就有独立小节 5 条） | 放 Discussion 内 |
| 参考文献格式 | **数字上标引用**（按出现顺序编号，Science/Nature 风格）——**不是** Cell 的 author–date | 用 AJHG 专用 CSL/EndNote 样式 |
| Data Availability | **必须**，指向 dbGaP/EGA/GEO/**Zenodo**/课题组站 + **联系人**（Mohlke lab 先例：marginal+条件化 stats 公开） | Zenodo DOI + 联系人，见 §五 + DATA_AVAILABILITY_20260816.md |
| 伦理声明 / IRB / 知情同意 | 必须 | §五模板 |
| 性别/祖先报告 | **2023 NASEM 要求**：报告样本性别构成与祖先组成 | §五模板（eQTLGen/FinnGen 均已报告） |
| CRediT Author Contributions | 必须 | §五模板 |
| 推荐审稿人 | **3–5 个**（cover letter 内） | §六建议名单 |
| 代码可用性 | 新程序必须发表时公开，URL 放正文 **Web Resources** 段 | 打包清单 §五 |
| 图形摘要 | 标注可选，**建议提交**（Cell Press 通用规格 1200×1200 方形单面板；scispace 模板写 1200×300 条图，投前确认） | 做一张"漏斗→15 候选→双复现"浓缩图 |
| 投稿系统/格式 | Editorial Manager（editorialmanager.com/ajhg）；接受 **Word/LaTeX，不接受 PDF/RTF**；12pt 双倍行距；图单独文件 | 投前按此组装 |

> ⚠️ 字数硬限：正文+摘要 ≤60,000 字符。Ray et al. 2025 全文 ~58k 字符，可作标尺。
> 编辑初筛（desk rejection ~40–50%）：title/abstract/首图必须立刻让审稿人看到"方法的信任边界
> + 产出了什么"——**Fig 1 必须是浓缩图**（31,373 对→FDR-core 982→121 strong→15 候选→双复现 或 精度漏斗）。

## 二、对标样稿与结构模板（投稿导师调研核实的 AJHG 同类论文）

### 2.1 Ray et al. 2025（AJHG 112(7):1597–1609）——"管线 + 单基因深挖"模板
- DOI 10.1016/j.ajhg.2025.06.001 —— "Single-cell transcriptome-wide MR and colocalization… ASCVD"
- 摘要七步流水：① 问题（GWAS 定位多但连因果基因难）→ ② 缺口（bulk TWAS 无细胞分辨率）→
  ③ "Here, we present a pipeline combining cis-MR with colocalization" → ④ 头号数字（440 关联、88% 被 bulk 漏掉）
  → ⑤ 复现（21 个外部 eQTL 复现）→ ⑥ 单基因深挖（LIPA，PheWAS + scRNA/IHC）→ ⑦ 落点（靶点选择）。
- **Results 按管线流水组织**：Fig1 管线 → discovery MR → replication MR → coloc（PPH4≥0.8）→ 细胞类型 → LIPA 案例 → 验证。
- Discussion 内 "Our study has limitations" 独立小节 5 条（含欧洲血统、缺乏相关组织 sc-eQTL）。
- 腔调：全程 "genetically proxied expression" / "drive associations" / "consistent with" 弱因果措辞，"causal" 只留给方法名。

### 2.2 Wilson et al. 2025（AJHG 112(11):2693–2707）——**我们这篇的最强结构模板**
- 骨骼肌 eQTL 图谱：18,818 条件信号 → 与 26 GWAS coloc → **2,252 coloc → 1,342 候选基因** →
  方法学洞见（仅 37% 对应最近基因、44% 离 TSS >50kb → "物理邻近不足以推断靶基因"）→ 多组织 T2D → 功能验证（luciferase）。
- **这正是"全转录组 MR/coloc 图谱 + 大量候选 + 一条可概括的方法学洞见 + 复现"的 AJHG 成功配方**：
  我们用"MR 精度/召回操作特性（漏斗 0.71%→25.6%）"当那条方法学洞见，与 Wilson 的"邻近性不足"是同类贡献。

### 2.3 Broadaway et al. 2024（AJHG 111(9):1899–1913）——"atlas 资源"模板
- 肝脏 eQTL 图谱：方法学洞见（条件分析先行多 37% coloc）+ 功能验证（MPRA）。
- 数据可用性实况：marginal+条件化 summary stats 放课题组站 + Zenodo + 联系人（Mohlke lab）——我们的范本。

### 2.4 Hemerich et al. 2024（AJHG 111(6):1035–1046）——"诚实披露"模板
- 8 方法基因优先排序 536 BMI 位点；**诚实句**："For most of the high-scoring genes, however, we found limited or no evidence for a role in obesity, including the top-scoring gene BPTF"。
- **证明"大多数候选缺乏先验证据"不是回避点，是诚实卖点**——前提是定位为"待功能验证的优先候选目录"。

### 2.5 mintMR（Lu et al. 2024, AJHG 111(8):1736–1749）——方法学论文先例
- 摘要结构（MR 价值→两挑战→方法→应用 35 性状→落点）；Simulation → Data analysis 组织；Limitations 融 Discussion 无独立小节。

### 2.6 操作特性/方法局限直接先例（决定我们卖点合法性的关键文献）
- **Zuber et al. 2022（AJHG 109(5):767–782，MR×coloc 综述）**：coloc 阳性通常蕴含 MR 非零，反之不成立。
  → 我们"MR sig 集内 coloc 精度仅 12.9%、强 coloc 罕见"正是该框架下的**预期行为**，要用 Zuber 背书。
- **van der Graaf/Auwerx/Kutalik et al. 2025（AJHG 113(2):309–323）**：**MR 操作特性论文直接先例**——
  用 945 条黄金标准酶-底物关系基准，**精度 35–47% 但召回仅 3.2–4.6%** → **AJHG 完全接受"诚实的操作特性（低召回）+ 靶点回报"的方法学评估文章**。我们的"精度漏斗"就是同一文体。

### 2.7 我们的差异化
- 我们**不打** Ray 的单细胞分辨与功能验证；我们打**全转录组尺度（31,373 对）+ 操作特性曲线（方法学原创）+ 双独立复现**。
- 卖点定位一句话：**给全转录组 cis-MR×coloc 一个可验证的信任边界，并交付可下载的优先基因图谱**。

## 三、摘要模板（≤200 词，已压缩）

> 背景句（一句话问题）→ 方法（数据+工具）→ 结果（三个数字组）→ 结论（诚实定位）。
> 数字全部与 FINAL_STORY/INTEGRITY_AUDIT 一致；候选措辞"hypothesis-generating"。

草稿（~195 词）：

```
Cis-Mendelian randomization (cis-MR) and colocalization are widely used to prioritize
blood-expressed genes as drug targets, but the operating characteristics of single-instrument
cis-MR across the whole transcriptome have never been systematically quantified. We performed
an exhaustive cis-MR x colocalization scan of 31,373 gene-outcome pairs (eQTLGen, n=31,684)
against type 2 diabetes, coronary artery disease, and fasting glucose. After per-outcome
Benjamini-Hochberg FDR control (q<0.05), 982 pairs were MR-significant, of which 121 reached
strong colocalization (PP.H4>=0.8): coloc yield 12.3%, and tightening the MR p-value threshold
calibrated the yield monotonically from 0.7% to 25.6% — an actionable operating-characteristic
curve for MR-based target screening. All 106 previously reported loci were reproduced, and 15
additional candidates were identified (9 candidate effector genes at known risk loci and 6
without T2D/CAD catalog records), with direction supported by GTEx v8 (6/7) and FinnGen R11
(9/9 gene-level), 4 reaching nominal significance in FinnGen. These results provide a
genome-wide evidence base for cis-MR target screening, quantify its operating characteristics,
and nominate 15 hypothesis-generating candidate effector genes, without implying causality.
```

（注：投前逐字再核一遍词数 ≤200。）

## 四、图表精选方案（≤7）——按 figure 调研升级版

> **figure 风格调研核心结论（2026-08-16 已整合）**：
> - AJHG 单张主图普遍 **3–6 面板**、每面板只讲一个结论；"一图讲一件事 + 多面板递进"是标准形态。
> - **审稿人对 coloc 类文章最期待的正是"逐基因座 GWAS+QTL 镜像图 + PP.H4 标注"（区域共定位图）**——
>   我们第一版几乎缺这张，必须补。
> - 具体图面规范（面板编号/字体/尺寸/配色/统计标注）见 §九。

主图方案（7 个，图号按引用序重排）：

| # | 主图 | 内容 | 关键升级点（vs 第一版） |
|---|---|---|---|
| Fig 1 | **浓缩图 / 管线示意图**（编辑第一眼） | 31,373 对→FDR-core 982→121 strong→106 已知重现→15 候选→GTEx/FinnGen 双复现 的漏斗式工作流 | **必须能单独回答"这套方法哪里可信、产出了什么"**（AJHG desk rejection 第一关） |
| Fig 2 | **精度漏斗（操作特性曲线）** | MR p 阈值→strong 率 0.71%→25.6%，配 95% Wilson CI；MR sig 集内 3.0% vs 集外 0.014% 双色分层 | 加 Wilson CI、分层漏斗感、Kutalik/Zuber 背书措辞进图注 |
| Fig 3 | **PP.H4 全分布按 MR 分层** | sig/grey/null ECDF + 0.75/0.8 双阈值线 + 阈值处堆叠条形（配 Wilson CI） | **多面板化**：(A) ECDF (B) 阈值命中条形 (C) MR 状态→H4/H3/其他转换计数 |
| Fig 4 | **图谱分布（Fuji-plot 式）** | 121 strong（+2 灰区）按染色体分布，点色=结局、点径=多效性（关联性状数）；known/candidate 形状区分 | 曼哈顿升级为 Fuji 式 + 多效性点径编码 |
| Fig 5 | **代表性基因座区域共定位镜像图**（**必须新增**） | 4–6 个代表性位点（如 RBM6/CNNM2/PLAUR/LAMC1 + 灰区 AP3S2/ZNF19）：上 GWAS −log10p（LD 着色、lead 紫色菱形带 rs 号）、下 eQTL 或基因级 PP.H4 条形、共享基因轨道、PP.H4 直标 | **审稿人对图谱类文章的核心期待**；参照 Abood Fig2/mintMR Fig2c |
| Fig 6 | **收敛验证面板** | HEIDI/Steiger/GTEx/FinnGen | 柱状图升级：(a) SMR b vs P_HEIDI 散点 (b) Steiger 方向 z 镜像或比例+CI (c) GTEx 跨组织效应森林图 (d) FinnGen 标准 MR 森林图（OR+CI，p<0.05 加粗） |
| Fig 7 | **资源一览 + 复现覆盖诚实披露** | 资源工作流示意 + 摘要热图/矩阵 + 关键计数（982/121/106/15）配 Wilson CI + 对齐覆盖 9/15（60%） | 纯底表快照改版为"资源一览图 + 概要表"混合 |

降补充材料（S1–S6）：灰区 AP3S2/ZNF19 完整 vignette、PP.H4 阈值敏感性（A4）、LD 聚类（A8）、
6 弱候选列表、全部 123 strong 子集明细表（`strong_all_subset_20260816.csv`）、方法学细节。

> 图注自足要求（AJHG 审稿惯例）：每个图注独立可读——含 n、阈值（PP.H4≥0.8、P_HEIDI=0.05）、
> 方法（coloc.abf/SMR/Steiger）、LD 参考人群（1000G EUR）、样本量、CI 类型（Wilson）、显著性符号定义。

## 五、必备投稿文件模板

### 5.1 Data Availability 声明
逐字可用蓝本见 `docs/DATA_AVAILABILITY_20260816.md` §1。要点：
- 全量图谱表（31,373 对）+ 123 strong 子集（`strong_all_subset_20260816.csv`）+ 15 候选注释（`candidate15_replication_20260816.csv`）+ schema 文档（`docs/SCHEMA_20260816.md`）→ **Zenodo DOI**（首选）
- 代码：GitHub（投稿时公开，作者 qgeng1465）
- 第三方数据出处表（eQTLGen/FinnGen R11/GTEx v8/1kg.v3）全部带版本与 URL

### 5.2 Ethics 声明（模板）
```
This study used only de-identified summary-level data from publicly available studies
(eQTLGen, FinnGen, GTEx, 1000 Genomes). No individual-level data were accessed. Ethical
approval and informed consent were obtained by the original studies; this analysis of
summary statistics did not require additional institutional review board approval.
```

### 5.3 性别/祖先报告（2023 NASEM 要求）
- **eQTLGen**：31,684 名欧洲祖先志愿者（全血），性别构成已在 eQTLGen 公开说明中（多数队列为男女性别混合，按队列原样，无性别分层）；本分析使用汇总统计，未做性别分层。
- **FinnGen R11**：芬兰人群（欧洲祖先），表型分布按 FinnGen release 报告。
- **GTEx v8**：838 名已故捐赠者，跨 49 组织；性别构成在 GTEx 官方元数据中（~34% 女性 / ~66% 男性，需投前核对）。
- **1000 Genomes EUR**：503 名欧洲祖先样本，作为 LD 参考。
- **诚实声明**：本分析未按性别分层做亚组分析；祖先构成均为欧洲，结论外推到非欧洲人群需谨慎。

### 5.4 Author Contributions（模板，投稿时按真实分工定稿）
```
Q.G. conceived and designed the study, performed the analyses, and wrote the manuscript.
```
（如导师参与，按其真实贡献补充。）

### 5.5 Code Availability
GitHub 仓库 URL + 运行环境（Python/R/PLINK 版本）+ 关键脚本索引（M* 脚本编号映射到仓库内 scripts/ 路径）。

## 六、推荐审稿人（3–5 个，投稿时按真实名单定稿）

原则：**方法论（cis-MR/coloc/操作特性）+ 应用（T2D/CAD 遗传学）+ 资源（eQTL 图谱）** 三领域各 1–2 人，
避开近 3 年共作者/直接竞争者（**重点避开：Ray et al. 2025、Kutalik 组、Wilson/Broadaway 的 Mohlke 组、
Hemerich 组**——皆为潜在竞争者或可能同审同一稿件）。

候选方向（需投前逐一核实无利益冲突；**人名必须以真实、可检索为准，不允许编造**）：
1. **操作特性/方法学方向**（首选）：引用过 Kutalik 2025 或 Zuber 2022 的 MR×coloc 方法学作者。
2. **cis-MR/coloc 方法论专家**（e.g. 近期 ciseQTL-MR 方法论文通讯作者）。
3. **eQTLGen/全血 eQTL 资源作者**（图谱资源评价方）。
4. **2 型糖尿病遗传学**（FinnGen/DIAGRAM 方向）。
5. **冠心病遗传学**（CARDIoGRAMplusC4D 方向）。

> ⚠️ 名单必须以真实、可检索的人名为准，不允许编造姓名。投稿前由用户/导师确认。

## 七、诚实措辞红线（投稿必守）

1. 15 个候选 → "candidate effector genes not currently indexed in GWAS Catalog" / "hypothesis-generating"，**绝不写** "causal"/"novel"。
2. "GWAS Catalog 未报道" ≠ "文献未报道"。
3. **41/106 = 38.7% 的 strong coloc 区域 GWAS 峰 p<5e-8 的 caveat 必须披露**（摘要/正文/图注择一处，正文必写）。
4. 6 个弱峰候选只进补充材料。
5. 方向复现 ≠ 效应量复现；VSIG8 冲突如实报告。
6. 对齐覆盖 9/15（60%）诚实披露（6 个 lead 未在 FinnGen R11 定位）。
7. 不承诺家族级显著性；图谱定位为"假设生成 + 描述性"。

## 八、执行清单（投稿前）

- [x] 后台调研 agent 结果整合进本文档（投稿导师调研 + figure 风格调研，2026-08-16）
- [ ] 图全部按调研结论升级（区域共定位镜像图、Fuji-plot、森林图、配色、图注自足）→ 7 主图定稿
- [ ] 摘要最终压到 ≤200 词（单段、无引用）
- [ ] 正文草稿按 Ray/Wilson 结构成稿（IMRaD），字符数 ≤60,000，参考文献用数字上标
- [ ] 标题定稿（动词+机制+资源，含 "atlas" + "operating characteristics" 关键词）
- [ ] 图形摘要做一张（漏斗→15 候选→双复现 浓缩单面板）
- [ ] Data Availability 落 Zenodo 包并取 DOI
- [ ] Ethics/性别祖先/Contributions/Code 模板定稿
- [ ] 推荐审稿人名单核实（不编造）
- [ ] Cover Letter（差异化卖点：全转录组尺度 + 操作特性曲线 + 双复现，见 FINAL_STORY §五）
- [ ] 方法学审稿雷点堵住（§八.2）：条件分析敏感性 / hg38-hg19 坐标校验 / 欧洲祖先披露 / PP.H4 阈值敏感性 / FinnGen 57% 覆盖解释
- [ ] cover letter <1 页（§十：首段"可验证的信任边界 + 可下载图谱"、为何 AJHG、披露预印本+AI 使用）
- [ ] commit 收口（作者 qgeng1465，无 Co-Authored-By；push 暂缓）

### 八.2 方法学审稿雷点（投稿导师调研预警，必须提前堵住）

1. **eQTL 多信号/邻近性**：Broadaway/Wilson 证明"用 marginal eQTL 直接 coloc 会低估"。→ 我们未做条件分析（APEX/SuSiE），**必须**在 Limitation 写，或补一列敏感性结果。
2. **坐标构建（hg38 vs hg19）**：历史上有 OpenGWAS(hg38)×eQTLGen(hg19) 错位 bug。→ 正文写明 harmonization + liftover 校验（已知位点 sanity check 进方法）。
3. **样本重叠/人群**：eQTLGen 与多数 GWAS 均欧洲人群 → 主动披露 European-ancestry 限制普适性（学 Ray 写法）。
4. **多重检验/阈值预设**：PP.H4 阈值写死并给敏感性（A4 已有）；MR 区段按 Zuber 框架交代"MR 非零但无 coloc"计数（正是 106/818 读数——819 复核存活中 818 可评估、106 获共定位支持——的意义）。
5. **FinnGen 复现覆盖率 9/15=60%**：写成 "independently replicated in a second biobank with 60% signal coverage"，解释未覆盖位点原因（lead 未定位/无显著 eQTL），**不写"15 个全复现"**。

## 九、图面规范（figure 风格调研定稿）

1. **面板结构**：AJHG 单张主图 3–6 面板、每面板一个结论、多面板递进；单面板图不达标。
2. **面板编号**：大写加粗 **A / B / C** 左上角（Cell Press 规范）；图中**不写标题**（标题只放图注）；跨面板字号一致。
3. **字体**：Arial/Helvetica（首选 Avenir）；出版尺寸下 6–8 pt，全图统一，最小 ≥6 pt；字体嵌入。
4. **尺寸**：AJHG 双栏，主图多用双栏 **170 mm** 横向；单栏 85 mm。
5. **分辨率**：位图 300 dpi、半色调 600 dpi、线稿 ≥1000 dpi；RGB；优先矢量 PDF/EPS。
6. **配色**：低饱和 + 色盲安全（CUD），推荐 Okabe-Ito/viridis/ColorBrewer；**全论文配色一致性**（MR 状态、结局、组织同色贯穿）；中性灰 #8C8C8C 用于辅助元素；主图元素 ≤6 色；禁纯红标显著性、禁彩虹渐变、热图避免红蓝+绿紫双极色；可用线型+形状冗余编码。
7. **统计标注**：比例配 **95% Wilson CI** 与 n；ECDF 加 0.75/0.8 阈值线与累计计数；点±CI 或森林图式粗横线+菱形。
8. **区域共定位图规格**（审稿人最期待）：上板 GWAS −log10(p)（x=Mb 标注 build，点按 lead LD r² 着色，lead 紫色菱形带 rs 号，typed/imputed 区分，可叠重组率线）；下板 eQTL −log10p 或基因级 PP.H4 条形；共享基因轨道；**PP.H4 数值直标面板上**；图注注明工具与 LD 参考人群（1000G EUR）。

## 十、Cover Letter 要点（<1 页）

1. 首段直接点出**对广泛人类遗传学读者的可见后果**：给全转录组 cis-MR×coloc 一个可验证的信任边界 + 可下载图谱。
2. 说明**为何是 AJHG**（而非 HGG Advances/Genome Research/Nat Genet）。
3. 披露：预印本（如适用）、受限数据、dbGaP/EGA 登记、ancestry 表述、**AI 工具使用**。
4. 建议审稿人 3–5 名（§六）。
5. 差异化卖点（FINAL_STORY §五）：全转录组尺度 31,373 对 + 操作特性曲线（方法学原创）+ 跨 3 结局 + GTEx/FinnGen 双独立复现；与 Ray 互补（他答细胞类型，我们答工具信任边界）。
