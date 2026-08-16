# 双通道 cis-MR 图谱（dual-channel cis-MR atlas）——同一变异经全血转录本 vs 循环蛋白，因果信号落在哪个分子层？

- 物理路径：本项目本地工作目录（公开仓库不含本地路径）
- 版本：**v0.14**（2026-08-16 FDR-core 定稿；v0.13 为 08-16 早段 raw 口径）
- 预注册：`docs/PREREGISTRATION.md`（v0.1.0，2026-08-06 生成，哈希锁定；v5 扫描层修订于扫描启动前追加并重哈希）

> **⚠️ 定稿口径（2026-08-16，FDR-core）**：投稿正文一切数字以 `docs/FACTS_20260816.md`（唯一数字来源）与
> `docs/manuscript/`（manuscript.md + 生成的 Word）为准。**主口径 = 分结局 BH-FDR q<0.05：982 MR-significant →
> 121 strong coloc（12.32% yield）→ 106 已知重现 + 15 候选效应基因**。本文档正文为历史工作记录快照，
> 其中"23 候选 / 131 strong / 129 / 11/12 / 12/21"等数字属早期 raw p<0.05 口径，已被 FDR-core 取代，
> 仅作历史参照。

---

## 摘要

**一个 cis 变异可以同时改变两个分子层——全血细胞转录本（eQTL）与循环血浆蛋白（pQTL）。当同一基因座对同一结局两个通道都有 MR 信号时，因果信号到底落在哪个分子层？** 本项目系统地在**完全对称的工具变量标准**下（cis ±1 Mb、p<5e-6、EUR LD clump r²<0.01@1000 kb）同时运行两条 cis-MR 通道——**通道 A：eQTLGen 全血转录本（n=31,684）** × **通道 B：deCODE 血浆蛋白（n=35,559）**——对三个代谢结局（T2D / CAD / FBG）输出一张 **基因 × 结局 × 介质层优先级** 的描述性图谱，并用**共定位闸门（PP.H4）**过滤"MR 显著但不共定位"的多效性信号。最终论文定位于转录组 cis-MR × coloc 图谱（蛋白通道数据已如实标注未进入论文分析，见 `docs/DATA_AVAILABILITY_20260816.md` §2b）。

核心发现（2026-08-16 定稿快照，FDR-core）：全转录组 cis-MR × coloc 扫描 **31,373 对**（QC 31,371）→
**per-outcome BH-FDR(q<0.05) 982 对 MR-significant**（t2d 394 / cad 576 / fbg 12）→ **121 个 strong coloc**
（PP.H4≥0.8，t2d 65 / cad 54 / fbg 2）＝ **coloc yield 12.32%**（+2 灰区 AP3S2×T2D / ZNF19×CAD）。
其中 **106 个已知位点 100% 重现 + 15 个候选效应基因**（9 已知位点候选 + 6 无 T2D/CAD catalog 弱候选；
FDR 掉出 8 个 raw strong：LAMC1/TPD52/SENP6/HMGN3/MT3/RPL13/ZBTB46/ZNF100）。方法学操作特性：
名义漏斗 0.71%→25.6%，FDR-core 12.32% 与 stage-2 grid 12.96% 收敛于 ~12–13%；MR-sig 集外 strong 仅
2/27,123（决定性负边界）。独立验证层：HEIDI 71.7%、Steiger 96.1%、GTEx 6/7（候选层）、FinnGen 9/9 基因级
（覆盖 9/15，4 个 p<0.05）。诚实 caveat：41/106=38.7% 已知 strong 区域 GWAS 峰 p<5e-8；coloc.susie 外样本
LD 下不收敛（exploratory，LAMC1 双排除）。**关键数字与复现统计见 `results/fdr_core_20260816.csv`、
`results/candidate15_replication_20260816.csv`、`docs/FACTS_20260816.md`。**

**历史快照（2026-08-13 深化核查与审计，raw 口径，仅历史参照）**：置换标定（零假设 FP=1.45%）、HEIDI 全集（398/819=48.6%）、Steiger 方向（96.1%）、**GTEx 跨组织复现（44/63=69.8%，非全血 26/40=65%，同变异 6/6）**、InsPIRE 胰岛负性核查（10/11 反向 + KCNJ11 胰岛 null）、两轮独立数据核查（8 项重算全过 + 3 处诚实修正：LD 面板定性、M16 列名、LRIG1 不可复现）。

**2026-08-15/16 全量扫描 + 新候选发现（M20–M28，取代 08-13 的"阴性"主结论）**：
1. **build 错位修复**：OpenGWAS 结局 hg38 vs eQTLGen hg19 坐标窗口错配 ~1.7 Mb → 改纯 rsid 匹配，全量重跑（误结果存档 `results/archive/buildbug_20260815/`）。
2. **全转录组全量 coloc（FDR-core 定稿）**：**31,373 对**（QC 31,371）入闸，**982 对 BH-FDR 显著 → 121 个 strong coloc**（其中 2 个灰区），其中 **106 个已知位点 100% 重现 + 15 个候选效应基因**。
3. **15 个候选效应基因**（`fdr_core_20260816.csv`，FDR-core 主口径）：9 个位于已知 T2D/CAD 风险位点内（SLC12A3、CWF19L1、U6atac、CD101、RBM6、CNNM2、N4BP2L2、RIC8A、C2orf49）、6 个 catalog 无记录区域（PLAUR、TAGLN2、VSIG8、PDCD6、CLEC3B、CCDC19）；8 个 raw strong 被 FDR 掉出（LAMC1 等），不作为候选。
4. **独立复现（15 候选口径）**：GTEx v8 独立 eQTL 方向 **6/7 一致**（1 冲突 VSIG8，如实报告）；FinnGen R11 独立队列 **可定位子集 9/9 基因级 / 8/9 变异级方向一致、4 个 FinnGen 自身 p<0.05（RBM6/CNNM2/CD101/RIC8A，对齐覆盖 9/15=60%）**。
5. **方法学操作特性**（`m36b_funnel_20260816.csv`）：MR p 阈值**单调校准** strong 率 **0.71%→25.6%**；**FDR-core 12.32% 与 stage-2 grid 12.96% 收敛**；全量显著（raw）3.04%——首次量化"MR 筛药靶"的校准曲线。
6. **LD 独立位点 + 阈值敏感性**：PP.H4 阈值敏感性（FDR-core 内 ≥0.9：7.5% / ≥0.8：12.3% / ≥0.5：29.5%；集外 ≥0.5：11、≥0.8：2）**单调性不随阈值改变**。

**主结论（v0.14 定稿）**：从"警示性阴性结论"升级为**图谱资源 + 方法学操作特性 + 候选效应基因发现**三合一诚实正结果。完整故事见 `docs/PAPER_WRITING_PLAN_20260816.md（附录 A）`，定稿稿件见 `docs/manuscript/`，图见 `results/figures/20260816_F*_v2.png`。

---

## 复现指南（Reproducibility，2026-08-16 FDR-core 定稿）

从源数据到论文数字的完整流水线（输入 → 脚本 → 产物）。所有分析脚本在 `scripts/`，随机种子、软件版本与门柱见 `docs/PREREGISTRATION.md`；口径与缺失标记见 `docs/SCHEMA_20260816.md`。

| 步骤 | 脚本 | 输入 | 产物 |
|---|---|---|---|
| 全转录组 cis-MR × coloc 全量扫描 | `M20*`–`M24` | eQTLGen + 三结局 GWAS（rsID 匹配，hg19） | `results/coloc_full_{t2d,cad,fbg}_20260815.csv`（31,373 对） |
| 新 strong 候选发现 | `M25`/`M25b` | `coloc_full_*` | `results/m25_new_strong_annotation_20260816.csv`（23 候选 raw 口径） |
| GTEx v8 独立方向复现 | `M26` | GTEx v8 | `results/m26_gtex_replication_new23_20260816.csv` |
| 名义精度漏斗 | `M27` | `coloc_full_*` | `results/m27_precision_funnel_20260816.csv` |
| FinnGen R11 独立队列复现 | `M28` | FinnGen R11 sumstats | `results/m28_finngen_replication_new23_20260816.csv` |
| **分结局 BH-FDR 主口径重算（FDR-core）** | `M36b_fdr_recompute_20260816.py` | `coloc_full_*` + grid | `results/fdr_core_20260816.csv`（982）、`results/candidate15_replication_20260816.csv`（15）、`results/m36b_funnel_20260816.csv`、`results/m36b_summary_20260816.csv` |
| 强共定位全子集（121+2 灰区） | 见 §4.1 打包清单 | `fdr_core` + `coloc_full_*` | `results/strong_all_subset_20260816.csv` |
| coloc.susie 敏感性（exploratory） | `M34b` | 6 位点 + 1000G EUR LD | `results/m34_coloc_susie_20260816.csv` |
| 图重构（9 主图 + S1） | `M37_figures_revised_20260816.py` | 各结果表 | `results/figures/20260816_F{1..9}_v2.png` + `S1_resources_v2.png` |
| Word 稿件生成 | `M36_build_word_ajhg_20260816.py` | `docs/manuscript/*` | `docs/manuscript/AJHG_submission_Qiushuo_Geng_20260816.docx` |

**关键数字**：FDR-core 982 MR-significant → 121 strong（12.32% yield）→ 106 已知重现 + 15 候选；名义漏斗 0.71%→25.6%；MR-sig 集外 strong 仅 2/27,123。**唯一数字权威**：`docs/FACTS_20260816.md`。

---

## 0. 论文故事线（研究思路，一页讲完）

> 目标：**单篇文章，完整故事，定投 AJHG（IF 7.7，Article 类型）**。以下为论文的目标叙事；标注 ✅=已完成、🔲=计划中。完整破局方案见 `results/BREAKTHROUGH_PLAN_20260815.md`，对抗性评审见 `docs/INTEGRITY_AUDIT_20260816.md（附录 A）`（2026-08-15，已按其修正叙事：降级 B/D、删"双介导/打哪里/新候选"三处过度声明、召回统计待分层重推）。

**科学问题（一句）**：当一个共享的 cis 变异同时关联全血转录本（eQTL）与循环蛋白（pQTL）时，因果信号到底经哪个分子层作用于 T2D / CAD / FBG？

**生物学意义 / 现实意义（不灌水，为什么这事值得做）**：
1. **图谱即公共资源**：一张全转录组共定位图谱 + 分子层一致性表，是药物靶点 MR 社区可直接复用的证据底座（谁做靶点 MR 都该先查这张表）——这是**推动领域发展**，也是本项目唯一真正可交付的社区贡献。
2. **方法学诚实是竞争力**：单工具 cis-MR 的**精度**（筛出来的有多少真）与**召回**（真的有多少被筛出）首次被系统量化，等于给整个"MR 筛药靶"实践贴一张操作特性曲线，把领域里的隐性浪费变成可讨论的数字——**促进方法进步**。这是方法学故事的主轴，而非"介质层=打哪里"（后者无方法学依据，见评审 §6）。

**为什么这个问题值得发（审稿人会关心）**：主流药物靶点 MR 只用单条通道（pQTL 或 eQTL），隐含假设"筛到的信号 = 可定位的因果信号"。但这条假设从未被系统检验——单工具 cis-MR 作为筛选工具，**精度**（筛出来的有多少是真）和**召回**（真的有多少被筛出）都是未知数。这个项目用一张全转录组共定位图谱 + 层分解表回答这个问题。

**设计逻辑（为什么这套设计能回答）**：
1. **双通道对称**：同一套工具变量标准、同一 MR 方法、同一共定位闸门，在转录层与蛋白层背靠背运行 → 层间比较方法上公平（介质层可分离性）。
2. **共定位闸门**：只有与 GWAS 共享因果变异的信号才算"介质层可定位" → 把纯 MR 多效性/LD 假阳性分开。
3. **全转录组扩展（✅，M20–M24）**：coloc 从 MR 显著子集扩到全量 31,373 对 → 首次可算 MR **召回**，把"精度 12.9%"升级为"精度 × 召回"二维操作特性。

**可证伪预期（诚实，结果已揭晓）**：预期"MR 显著集仅 ~1/7 获共定位支持、且 MR 只召回 ~1/5"——实测为**低精度 + 高召回**：FDR-core MR-significant 集 strong 率 **12.32%**（与 stage-2 grid 12.96% 收敛），全量 sig（raw）3.04%；MR-significant 集外 strong 仅 2/27,123。即 **MR 作为预筛召回好、校准曲线单调**；全量 coloc + 操作特性校准优于网格预筛。

**结果故事线（按评审 §7：A+C 为核心，B/D 降级佐证，不做"四拼盘"）**：
| 章节 | 内容 | 状态 |
|---|---|---|
| R1（核心 A） | 全转录组共定位图谱（31,373 对 / **FDR-core 982 → 121 strong**）+ **操作特性漏斗**（MR p 单调校准 0.71%→25.6%；**FDR-core 12.32% 与 grid 12.96% 收敛**） | ✅ M24 + M36b + 图 |
| R2（核心 A 升级） | **15 个候选效应基因**（9 已知位点候选 / 6 无 catalog 弱候选）+ GTEx 6/7 + FinnGen 9/9 基因级一致（覆盖 9/15）独立复现 | ✅ fdr_core + candidate15 |
| R3（降级 B） | 蛋白通道（deCODE 仅 PCSK9/APOC3）降级为补充材料；UKB-PPP 覆盖 | 已降级 |
| R4（降级佐证 D） | 证据简版（HEIDI/Steiger/FinnGen/组织——**非 4 独立维度**，InsPIRE 反向进 Limitation） | 部分✅ |

**诚实边界（红线，全文生效，评审 §4/§6 加强）**：15 个候选效应基因是**候选评估非因果发现**（6 弱候选不作主表头条）；KCNJ11/PCSK9/APOC3 等已知位点只当校准案例；**候选须逐位点 GWAS Catalog 注释、报已知位点占比（9 已知位点候选 / 6 无 catalog）**；"catalog 未报道"≠"文献未报道"（catalog 滞后）；**共定位≠介导，eQTL coloc 与 pQTL coloc 不构成独立介导证据**；MR 单工具 Wald ≡ SMR 数学恒等，不存在"三方法独立互证"；**coloc 弱峰 caveat 必须披露（106 已知中仅 41/106=38.7% 已知强候选 GWAS 峰 p<5e-8）**；**coloc.susie 外样本 LD 下不收敛（exploratory），LAMC1 被 FDR 与多信号证据双排除**。

---

## 1. 科学问题与创新点

### 1.1 问题

主流药物靶点 MR 通常只用单条通道（pQTL 或 eQTL），无法区分一个**共享的 cis 因果变异**到底经哪条分子链作用于疾病结局。识别因果信号所在的介质层，决定干预应优先作用在哪个分子层（"精准介质"问题）——注意这是**因果介质层**判定，不推荐具体药物模态。

### 1.2 创新点（给审稿人一眼看懂）

1. **双通道对称性**：同一套锁定工具变量标准、同一 MR 主方法、同一共定位闸门，在两个独立介质层上背靠背运行——跨层比较在方法上公平，避免单通道"有无显著"的层间不可比。
2. **介质层优先级输出**：不是单一 hit 列表，而是每基因×结局的介质层四态/优先级（transcript-only / protein-only / both / neither）。
3. **共定位闸门**：MR 信号必须与 GWAS 信号共享因果变异（PP.H4）才算"介质层可定位"，把纯 MR 多效性/LD 假阳性与真正的共定位信号分开。
4. **预先注册**：主分析门柱（IV、MR 方法、FDR、对照）哈希锁定，新增分析只能以"探索性/预先登记补充"出现。

### 1.3 不做的事（诚实声明）

- 不产出"该用单抗还是 ASO"的模态推荐。
- PCSK9/CETP 等已知位点仅用于确认双通道在已知药理位点可同时定位、MR 可检出（校准），不用于任何"蛋白主导"式推断。
- 不按 |β| 幅度跨基因排序（只比方向，不比幅度）。

### 1.4 预注册假设（可证伪）

- **H1**：分泌型蛋白 vs 细胞内蛋白在 pQTL-MR 显著性率与方向一致性上存在差异（Fisher/χ² on 2×2，α=0.05 two-sided）。
- **H2**：已知药理位点（正/负对照）的介质层优先级做留一法 LOOCV AUC（AUR）> 0.5（one-sided；若 AUR≤0.5 仅做描述性报告，预注册约定）。

---

## 2. 数据

### 2.1 暴露（两通道对称）

| 层 | 数据源 | 样本量 | 说明 | 状态 |
|---|---|---|---|---|
| 通道 A·转录本 | eQTLGen 全血 cis-eQTL（Võsa 2021） | n=31,684 | 16,923 基因显著文件；真实等位频率文件（240 MB） | ✅ 已就绪 |
| 通道 B·蛋白 | deCODE 血浆 cis-pQTL（Ferkingstad 2021） | n=35,559 | 4,719 蛋白，hg38 per-SD；9/11 测试蛋白文件已就绪 | ✅ M3/M5 完成 |
| 通道 B·补充 | INTERVAL SomaScan pQTL（OpenGWAS prot-a-*） | n=3,301 | 肝酶/膜蛋白 cis-pQTL 功效不足（诚实负结果） | ✅ 已测 |
| 组织敏感性 | GTEx v8（肝/胰/全血/脂肪/肌肉/冠脉） | — | 回应"全血盲"的肝/肠/胰岛特异基因 | ✅ 已完成 |
| 胰岛敏感性 | InsPIRE 胰岛 eQTL（Viñuela 2020） | — | 独立基因 eQTL，β 细胞层 | ✅ 已完成 |

### 2.2 结局（OpenGWAS）

| 结局 | ID | 效应 | 样本量 |
|---|---|---|---|
| T2D | `ebi-a-GCST006867` | logOR | 主文档 n=655,666；**本机 full 文件实测 per-variant N mode≈573,704**（范围 573–579 K，2026-08-16 定稿口径） |
| CAD | `ebi-a-GCST005194` | logOR | n=296,525 |
| FBG | `ebi-a-GCST005186` | beta | n=58,074 |

全项目效应量均 per-1SD，只比方向不比幅度。

### 2.3 工具与参考

- **LD 参考**：1000G Phase 3 EUR（1kg.v3，hg19），本地 plink 1.9 clump。
- **SMR 1.3.1**：已装（`tools/smr`）并跑完双通道 SMR/HEIDI（M9 蛋白 15 对 + M10 转录 128 测试，2026-08-13，见 `results/archive/DECODE_PIPELINE_SUMMARY.md` M9/M10 节）。

---

## 3. 方法（预注册锁定 + 探索性标注）

### 3.1 工具变量标准（预注册锁定，两通道一致）

- cis 窗口 ±1,000 kb；p<5e-6；LD clump r²<0.01@1000 kb（EUR 1000G P3）；palindromic 处理 action=2。
- 蛋白通道 deCODE 侧无 effectAlleleFreq（assocvariants.annotated 未就绪）→ **palindromic 保守排除 + ImpMAF 近似**（note 标注）。

### 3.2 MR 方法（预注册锁定 + 2026-08-13 首读修正）

- 敏感性方法 `mr_ivw_mre` / `mr_ivw_fe` / `mr_weighted_median` / `mr_egger_regression` 全部报告；nsnp=1 明确标注为 Wald 退化（不伪敏感性）。
- **首读统计量修正（2026-08-13，多 agent 评审发现 `mr_ivw_mre` 在 nsnp=2-3 时 SE 塌缩 7-16 倍）**：nsnp=1 → Wald；**nsnp≤3 → IVW-FE 首读**；nsnp≥4 → IVW-MRE 首读。FE 本为预注册敏感性方法，此修正不移动任何门柱。
- 多重检验：按结局 BH-FDR q<0.05。

### 3.3 扫描层（预注册 v5 修订，预先登记、探索性）

- 全量扫描层：每基因最强 cis-eQTL（lead variant）单工具 Wald，按结局 FDR；**扫描层为 hypothesis-generating**，结论级证据以深度层为准。
- 命中基因第二阶段：本地 plink（1000G EUR）clump + IVW-MRE 复核，FDR 在 stage-2 结果集内重算（诚实双报 nominal + q）。
- 敏感性层（GTEx / InsPIRE / 单工具 Wald）**不改变**四态主分析（全血 eQTL × 血浆 pQTL）的门柱。

### 3.4 共定位闸门（M5）

- coloc.abf（p12=1e-5，敏感性 p12=1e-6），cis 全窗口；PP.H4 三档 strong ≥0.8 / moderate / none【待核：三档阈值以脚本实际实现为准】；回文位点保守排除。

### 3.5 预注册与诚实协议

- `docs/PREREGISTRATION.md` 哈希锁定（2026-08-06 首锁，v5 于扫描启动前重哈希）；变更只能以追加修订+重哈希记录，禁止事后调门柱求显著。
- 全网格如实报告（含空/失败行），落盘 `results/grid/`，不做选择报告。
- README 只呈现多管线叙事；真实运行轨迹（方法修正、bug 修复、失败记录）见 `docs/CHANGELOG.md`。

---

## 4. 结果（两通道对称呈现）

```
通道 A  全血转录本（eQTLGen）──┬─ 全量扫描层（49,866 对 → 982 hits）✅
                               ├─ stage-2 本地 clump+IVW 复核（nominal 819 / FDR 812 存活）✅
                               ├─ 候选基因深度层（HMGCR×T2D 单工具 b=−0.382 p=4.8e-5）✅
                               ├─ 转录共定位 M5（106 strong / 149 moderate / 563 none）✅
                               └─ 敏感性：GTEx 6 组织（18 个 p<0.05）/ InsPIRE 胰岛 ✅
通道 B  循环蛋白（deCODE）──────┬─ 蛋白 MR M3（5 蛋白 × 3 结局 = 15 对首读 ok）✅
                               ├─ 蛋白共定位 M5（strong=2 / moderate=0 / none=13）✅
                               └─ INTERVAL 补充（n=3,301 低功效，诚实负结果）✅
双通道整合 ──────────────────────── 药物靶点四态分类 + 介质层优先级表 ✅（描述性）
```

### 4.1 通道 A：全血转录本

**全量扫描层**（`transcript_grid_mr.csv`）：16,622 基因 × 3 结局 = **49,866 对**（结果文件行数；早期文档曾记 49,867），lead cis-eQTL 单工具 Wald，按结局 FDR。**982 hits**（T2D 394 / CAD 576 / FBG 12），其中 157 个 p<1e-6 强信号。

**stage-2 复核**（`transcript_grid_stage2.csv`，本地 plink clump + IVW-MRE）：982 hits 全量复核，906 唯一基因，978/982 对 MR 完成（4 对如实记录失败：结局无匹配 / clump 后 0 独立工具）。存活 **nominal（p<0.05 且方向一致）819（t2d 331 / cad 477 / fbg 11）／BH-FDR 812（329 / 472 / 11）**。

**候选基因深度层**（结论级）：HMGCR×T2D 本地 clump 后单工具 Wald **b=−0.382 p=4.8e-5**（负向，与"他汀降 LDL→降 T2D 风险"方向一致；stage-2 独立复现 b=−0.382 p=4.8e-5）；HMGCR×CAD 仍 NS（b=+0.037 p=0.137，符号翻转如实记录）。

**组织敏感性（GTEx v8，6 组织 × 21 基因 × 3 结局 = 309 组合）**：48 个完成 MR，**18 个 p<0.05**。全血层失明的位点在特异组织活过来：

| 基因 | 组织 | 结局 | b | p |
|---|---|---|---|---|
| APOB | Muscle_Skeletal | CAD | +0.190 | 4.8e-11 |
| CETP | Adipose_Subcutaneous | CAD | +0.124 | 3.3e-9 |
| NPC1L1 | Pancreas / Adipose | CAD | +0.102 / +0.080 | 3.8e-6 |
| PDX1（β 细胞 TF） | Pancreas | T2D | −0.108 | 4.6e-5 |
| PCSK9（肝） | Liver | CAD | −0.062 | 2.4e-4 |
| HMGCR | Muscle_Skeletal | T2D | −0.131 | 1.3e-3 |

**胰岛 eQTL（InsPIRE）**：TCF7L2×T2D Wald **b=−1.598 p≈0**；TCF7L2×CAD b=−0.132 p=6.2e-6；TCF7L2×FBG b=−0.110 p=2.0e-9；**SLC5A1（SGLT1）全血无 eQTL、但胰岛 eQTL→T2D p=1.1e-3**——"全血失明、胰岛救场"的直接证据。

**转录共定位 M5**（`transcript_coloc.csv`，stage-2 nominal 819 对入闸）：818 对通过 nsnp≥10 质控，**106 strong / 149 moderate / 563 none**（另 1 对未过质控如实记录）。

**2026-08-16 全量扫描（M20–M24，取代 M5 的"MR 显著子集"口径）**：coloc 从 stage-2 网格扩到**全量 31,373 对**（QC 31,371；49,866 MR 对中 18,493 对因无 eQTL 工具或 pval NA 无法 coloc，如实计入分母），**131 个 strong coloc**（其中 129 个在 MR 显著集 mr_p<0.05 4,248 对内、2 个灰区；t2d 67 / cad 60 / fbg 4）：
- **known 106/106（100%）重现**（GIT1 双源 pp4 逐位一致，本地 full 与 API 收敛）；`results/coloc_full_{t2d,cad,fbg}_20260815.csv`、`coloc_full_summary_20260815.csv`。
- **MR 显著集外 strong 极罕见**：nsig（0.05–0.5）2/14,294=0.014%（AP3S2×t2d、ZNF19×cad，均 GWAS 峰 5e-11 级）、null（≥0.5）0/12,829 → **MR 召回好，不丢信号**。
- **新候选发现**：全量扫描比 stage-2 网格多出 **23 个 sig 内 new strong**（`M25` 注释 + `M25b` 宽口径重分类，见 §4.7）。
- **精度操作特性**（`M27`，`results/m27_precision_funnel_20260816.csv` + 图 F1）：MR p<0.5→0.71%、<0.05→3.04%、<0.01→6.56%、<0.001→14.8%、<1e-4→24.6%、<1e-5→25.6%；stage-2 网格 12.96% ≈ 对应 mr_p~1e-3 层的精度，但**丢 23 个新候选**。

### 4.2 通道 B：循环蛋白（deCODE，M3/M5 完成）

**测试范围**：11 个候选/对照蛋白入测试，9/11 有 deCODE 文件（APOB/LDLR 未下载）；其中 5 个蛋白获得 cis-pQTL 工具，覆盖 **15 个蛋白×结局对首读 ok**（全网格 63 行方法级、45 行有结果，含 mre/fe/weighted median/egger 敏感性）。

**蛋白 MR 结果表**（按 p 排序，修正后首读；产物 `protein_decode_mr_primary.csv`）：

| 蛋白 | 结局 | nsnp | b | SE | p | 首读方法 | 共定位 PP.H4 | 注 |
|---|---|---|---|---|---|---|---|---|
| PCSK9 | CAD | 14 | +0.190 | 0.028 | 5.5e-12 | IVW-MRE | **1.000** | FE p=4e-26、WM p=1.2e-7、Egger 斜率 p=6.5e-5 全一致；F=47.5 |
| INSR | FBG | 2 | +0.125 | 0.050 | 0.0121 | IVW-FE | 0.0097 | MR 显著但不共定位（MRE p=0 为 SE 塌缩伪影） |
| APOC3 | T2D | 3 | +0.116 | 0.050 | 0.0203 | IVW-FE | <0.5 | MR 显著但不共定位（FE 下才显著） |
| APOC3 | CAD | 10 | +0.151 | 0.068 | 0.0271 | IVW-MRE | **0.997** | FE p=5.2e-9、WM p=2e-4 一致 |
| APOC3 | FBG | 3 | +0.044 | 0.023 | 0.0543 | IVW-FE | 0.177 | MRE p=9.1e-46 为 SE 塌缩伪影；FE 下不显著 |
| PCK1 | T2D | 1 | +0.210 | 0.189 | 0.27 | Wald | — | 单工具，F=1.24 |
| PCK1 | FBG | 1 | +0.070 | 0.074 | 0.35 | Wald | — | 单工具 |
| INSR | T2D | 2 | +0.078 | 0.084 | 0.35 | IVW-FE | — | — |
| INSR | CAD | 2 | −0.050 | 0.060 | 0.40 | IVW-FE | — | MRE p=0.039 为 SE 塌缩伪影 |
| ANGPTL3 | CAD | 19 | −0.010 | 0.014 | 0.45 | IVW-MRE | — | — |
| ANGPTL3 | FBG | 12 | −0.004 | 0.007 | 0.58 | IVW-MRE | — | — |
| PCSK9 | FBG | 5 | −0.007 | 0.019 | 0.70 | IVW-MRE | — | — |
| PCSK9 | T2D | 3 | +0.013 | 0.076 | 0.86 | IVW-FE | — | — |
| ANGPTL3 | T2D | 9 | +0.003 | 0.023 | 0.90 | IVW-MRE | — | — |
| PCK1 | CAD | 1 | +0.002 | 0.138 | 0.99 | Wald | — | 单工具 |

**诚实空结果（全网格如实报告，`protein_decode_mr.csv` 含失败行）**：
- **APOB / LDLR**：deCODE 文件未下载/不完整 → 6 对无数据。
- **HMGCR / GLP1R / GCG**：cis ±1 Mb、p<5e-6 内无变异 → 9 对无工具。
- **DPP4**：cis 内 6 个变异**全为罕见位点，不在 1000G EUR LD 参考面板** → 预注册 clump 步骤下无独立工具 → 3 对无工具。
- **INTERVAL（n=3,301）**：HMGCR/ANGPTL3 cis-pQTL 功效不足（无 p<5e-6 工具）；PCSK9/CETP/NPC1L1/APOC3 未收录于 SomaScan 面板 → 不进入 H1/H2 正式检验，蛋白主源为 deCODE。

**蛋白共定位 M5**（`protein_coloc.csv`，15 对入闸）：**strong=2 / moderate=0 / none=13**：
- **PCSK9×CAD**：PP.H4=1.000（p12e6=1.000），top SNP=**rs11591147（R46L）**——阳性对照校准命中，教科书级。
- **APOC3×CAD**：PP.H4=0.997（p12e6=0.974），top SNP=**rs964184**（APOA5/APOC3 区）——负对照边界案例：TG 通路确实影响 CAD，需在 M7 如实讨论。
- **MR 显著但不过闸门**（修正后首读 p<0.05 且 PP.H4<0.5，共定位闸门将其过滤——多效性/LD 驱动信号，正是闸门价值的演示）：INSR×FBG（FE p=0.0121，PP.H4=0.0097）、APOC3×T2D（FE p=0.0203，PP.H4<0.5）、APOC3×FBG（FE p=0.0543 边缘，PP.H4=0.177）。

### 4.3 药物靶点四态分类（双通道整合的核心产出，`drugtarget_fourstate.csv`，描述性）

| 类别 | 基因 | 转录信号（全血） | 蛋白通道 |
|---|---|---|---|
| **both** | DPP4 / GLP1R / INSR | DPP4×T2D b=−0.267 p=3e-3、GLP1R×T2D b=−0.389 p=0.034、INSR×T2D b=−0.938 p=1.3e-3 | deCODE 有文件 |
| **transcript-only** | KCNJ11 / TCF7L2 / PPARG / PRKAA1（+SLC5A2，信号 null） | KCNJ11×T2D Wald b=+0.595 p=4e-16、TCF7L2×FBG p=7.6e-144 | 转录 |
| **protein-only** | GCG / PCK1 | 全血无 cis-eQTL（糖异生酶） | 蛋白 |
| **neither** | ABCC8 / PDX1 / G6PC / SLC5A1 | 全血无 cis-eQTL | 双盲 |

**MR 落地后的通道轴修正（诚实标注）**：四态分类基于"通道可用性"（蛋白是否被 deCODE 测量），M3 完成后再核对工具可得性——DPP4/GLP1R 虽列 both，但 **DPP4（罕见位点无 LD 面板）与 GLP1R（无 cis p<5e-6 变异）实际拿不到 cis-pQTL 工具**；de facto 同时具备转录信号 + 蛋白 MR 证据的双通道对为 **INSR**（蛋白 INSR×FBG FE p=0.0121；转录 INSR×T2D b=−0.938 p=1.3e-3）。此为描述性观察，通道轴 ≠ 药物靶点蛋白轴。

### 4.4 跨通道对称读数（介质层可分离的核心证据）

| 基因座 | 全血转录 | 组织/胰岛转录 | 循环蛋白 | 共定位 | 读数 |
|---|---|---|---|---|---|
| **PCSK9**→CAD | 盲（肝表达，全血无 cis-eQTL） | GTEx 肝 b=−0.062 p=2.4e-4 | **b=+0.190 p=5.5e-12** | **PP.H4=1.000（rs11591147/R46L）** | 蛋白层 + 组织层定位，全血转录盲 |
| **HMGCR**→T2D | **b=−0.382 p=4.8e-5** | GTEx 肌 b=−0.131 p=1.3e-3 | 无 cis-pQTL（空结果） | — | 仅转录层可定位 |
| **APOC3**→CAD | 盲 | — | **b=+0.151 p=0.027** | **PP.H4=0.997（rs964184）** | 蛋白层定位（TG 通路负对照边界） |
| **INSR**→FBG/CAD/T2D | INSR×T2D b=−0.938 | — | **FBG FE p=0.012 / CAD FE p=0.40（ns）** | FBG PP.H4=0.0097（不过闸） | 转录/蛋白双通道，但蛋白 FBG 信号不共定位 |

### 4.5 2026-08-13 深化核查与审计

> 目的：把 12.9% 一致性率等关键读数做成可辩护的方法学估计 + 独立核查背书。
> 完整文档见 `results/archive/PUBLICATION_STORY_20260813.md`、`results/archive/SUMMARY_11H_20260813.md`。

| 分析 | 关键数字 | 结论性质 |
|---|---|---|
| **E1 置换标定** | 零假设 PP.H4≥0.8 FP=**1.45%**，106 个 strong 调用富集 ~69 倍 | ✅ 调用非噪声 |
| **HEIDI 全集（819 测试）** | 全 MR 显著集通过 **48.6%**、strong 子集 **71.7%** | ✅ 共定位命中 LD 同质性成立 |
| **Steiger 方向** | 73/76 = **96.1%** eQTL→outcome 正向 | ✅ 因果方向正确 |
| **P1 GTEx 跨组织复现（M15）** | **44/63 = 69.8%** 方向一致；非全血 26/40 = **65.0%**；同变异 **6/6** | ✅ 信号非全血伪影 |
| **M16 InsPIRE 胰岛（M15 对照）** | 11 显著中 **10 反向**（两次运行一致）；KCNJ11 胰岛 lead **null p=0.87** | ⚠️ 负性核查，红线维持 |
| **LD 面板核查（M18）** | 92 独立簇（保守上限）；位置聚类 89–88；**面板零杂合但 r² 仍反映真实 LD**（0.918 vs 真实 0.94–0.99） | ⚠️ 定性修正 |
| **编码版一致率 / 分母统一** | 编码 76/637=11.9%；规范口径 **106/818=12.96%**（819 中 TCF19×CAD 因 GWAS 区域数据不可用无法评估 coloc） | ✅ 口径定案 |
| **药物注释 / 富集** | 仅 KCNJ11/BMPR1A 有药（15 条=14 化合物）；g_SCS **0 显著**（阴性） | ⚠️ 探索性 |
| **FinnGen 修复 / 两轮核查** | rs17716350 对齐修复 92.6%；8 项独立重算全过 | ✅ 审计背书 |

### 4.6 结果文档索引（2026-08-13，`results/`）

| 主题 | 文档 |
|---|---|
| 发表故事 + 投稿计划 | `PUBLICATION_STORY_20260813.md` |
| 诚实最终结果总结 | `SUMMARY_11H_20260813.md` |
| 独立数据核查（2 轮，8 项重算 + 修正） | `verification_20260813.md` |
| 置换标定 / HEIDI 全集 / Steiger | `coloc_permutation_20260813.md`、`heidi_full_20260813.md`、`steiger_direction_20260813.md` |
| GTEx 复现 / InsPIRE 胰岛 | `gtex_replication_20260813.md`、`inspire_islet_20260813.md` |
| LD 聚类 + 面板核查 / 一致率口径 | `ld_clustering_20260813.md`、`concordance_denominator_20260813.md`、`coding_concordance_20260813.md` |
| 脂质校准 / UKB-PPP 双平台 / 样本重叠 | `lipid_validation_20260813.md`、`ukbpp_replication_20260813.md`、`overlap_sensitivity_20260813.md` |
| 药物注释 / 通路富集 / 证据简版 | `drug_annotation_20260813.md`、`pathway_enrichment_20260813.md`、`evidence_brief_20260813.md` |
| 蛋白管线总结 | `DECODE_PIPELINE_SUMMARY.md` |
| 预注册（哈希锁定） | `docs/PREREGISTRATION.md`（v0.1.0 + v5/v6/v7 修订） |
| 2026-08-16 故事 / 写作 / 审计 / 审稿输入 | `docs/PAPER_WRITING_PLAN_20260816.md`（含故事·期刊·审稿输入）/ `docs/INTEGRITY_AUDIT_20260816.md` |
| 23 新候选注释 / GTEx 复现 / 精度漏斗 / FinnGen 复现 | `results/m25_new_strong_annotation_20260816.csv` / `m25b_reclassify_20260816.csv` / `m26_gtex_replication_new23_20260816.csv` / `m27_precision_funnel_20260816.csv` / `m28_finngen_replication_new23_20260816.csv` |

---

### 4.7 23 个新候选效应基因 + 独立复现（2026-08-16）

> 这是全量扫描带来的正向发现。诚实口径：**候选评估，非发现宣称**；"catalog 未报道"仅指 GWAS Catalog。

**分类**（`M25` + `M25b` 宽口径重分类，结果稳健）：
- **15 个位于已知 T2D/CAD 风险位点 ±100–250 kb 内**（该基因本体未被 T2D/CAD catalog 报道）→ **新效应基因提名**。
- **8 个所在区域无 T2D/CAD catalog 关联**，但多数 GWAS 峰 p>5e-8 → **弱峰低置信假设**（不称新位点）。

**独立复现**：
- **GTEx v8**（`M26`，6 组织跨组织 best lead）：15 候选中 7 个可测，**6 方向一致 / 1 冲突（VSIG8，如实报告）**；8 个不可测（GTEx lead 不在结局 GWAS）。
- **FinnGen R11**（`M28`，独立队列全量 sumstats）：15 候选中 9 个可检，**基因级 9/9、变异级 8/9 方向一致，0 冲突**；FinnGen 自身 p<0.05 **4 个**（RBM6 6.9e-06、CNNM2 6.3e-04、CD101 1.2e-03、RIC8A 1.1e-02）；**对齐覆盖 9/15（60%）**。

**表**：`results/m25_new_strong_annotation_20260816.csv`（含每基因 hg38 坐标、catalog 判定、tier，23 候选 raw 口径）；`results/candidate15_replication_20260816.csv`（**FDR-core 15 候选 + 复现**）。图：`results/figures/20260816_F6_candidates15_v2.png`。

---

## 5. 诚实局限

1. **"全血盲"是转录通道的结构性局限**：肝/肠/胰岛特异基因（PCSK9/ANGPTL3/APOC3/NPC1L1/SLC5A1）全血无 cis-eQTL——组织/胰岛敏感性层是设计内回应，但不能替代全血全扫描的覆盖。
2. **蛋白通道覆盖不足**：11 测试蛋白中 APOB/LDLR 文件未下载、HMGCR/GLP1R/GCG 无 cis-pQTL、DPP4 罕见位点无 LD 参考 → 实际可用 5 蛋白。四态"both/neither"判定受通道可用性（非工具可得性）影响，已如实标注。蛋白 SMR 只覆盖这 5 蛋白（15 对）。
2b. **SMR/HEIDI 局限**：deCODE 无 effectAlleleFreq → .esd Freq 用 ImpMAF 近似 + `--disable-freq-ck`（探索性敏感性，非正式检验）；HEIDI 用 1000G EUR 参考 LD，与 deCODE 人群错配（转录 eQTLGen 同为 EUR 人群更可靠）；转录 SMR 的 p_SMR 显著（125/128）是选择偏倚的构造性必然（测试集全为 MR-显著预选位点），**零信息量不进正文**，唯一有信息量的读数是 HEIDI 通过率（89/126=70.6%）；且 **MR 单工具 Wald ≡ SMR 为同一统计量**（数学恒等，非独立方法），真正独立于单点估计的只有区域级 coloc 与 HEIDI。
3. **deCODE EAF 缺口**：assocvariants.annotated 未就绪 → palindromic 保守排除 + ImpMAF 近似（ImpMAF 不总是 effect allele 频率）；多等位基因 bug 行已排除。
4. **共定位闸门放大了"强共定位优先"**：MR 显著但不过闸门（INSR×FBG、APOC3×T2D）可能由多效性/LD 驱动——这是设计意图（过滤），但也意味着部分真实介质信号被保守丢弃；两通道头条 MR 的 MRE p 值经 2026-08-13 修正为 FE 首读（小工具数 SE 塌缩伪影，见 §3.2）。
5. **T2D 结局样本量**：主文档记 n=655,666；本机 full 文件实测 **per-variant N mode≈573,704**（范围 573–579 K，2026-08-16 定稿，与 62,892 的差异为 gwasinfo 子集/版本口径，投稿以 full 文件 per-variant N 为准）。
6. 转录通道全量扫描为单 lead 工具 Wald（hypothesis-generating）；结论级证据以深度层 + stage-2 复核为准。
7. **FinnGen 复现对齐（2026-08-13 已修）**：rs17716350 对齐 bug（原同 SNP 两行一负一正均标 aligned）已修（逐元素 `==` 对齐器 + 全表重算），文件现 108 行、rs17716350 单行 b=−0.0133(G) 与 discovery 同向、108/108 aligned；头条报双显著子集（T2D 27/27、CAD 19/19=100% 同向）。2026-08-16 另用 **FinnGen R11 全量 sumstats** 对 23 个新候选做独立队列复现（`M28`，见 §4.7）。

---

## 6. 下一步（2026-08-16 更新）

**A. 收尾当前故事（✅ 分析全完成，剩 commit）**：FinnGen R11 独立复现（M28 已完成）→ 图定稿（F1/F2 已并入 FinnGen 栏）→ 审计 P0 文档修复（已完成）→ README 同步 v0.13（本次完成）→ commit（push 暂缓）。

**B. 写论文**：按 `docs/PAPER_WRITING_PLAN_20260816.md（附录 A）` 三股线（图谱资源 + 操作特性 + 15 候选）成文；主图 F1（精度漏斗）+ F6（候选基因 + 双复现）+ 全量图谱资源表。定投 AJHG（投稿准备指南见 `docs/AJHG_SUBMISSION_GUIDE_20260816.md`，图表精选 ≤7 个、摘要 ≤200 词）。

**C. 不再开新分析战场**（判断：15 候选的下一步 SMR/功能验证超出纯干实验数据可得性，精度漏斗已把方法学故事说满，边际收益递减）。可选锦上添花（不进关键路径）：GTEx 扩到 49 组织、TWAS（FUSION）交叉验证、UKB 表型关联外延。

---

## 7. 完整性与预注册纪律

- **预注册**：`docs/PREREGISTRATION.md`（哈希锁定）+ `docs/PREREGISTRATION.yaml`（OSF 框架镜像）+ `docs/osf_provenance.txt`（溯源证据）。主分析门柱锁定，禁止事后调门柱求显著。
- **全网格如实报告**：所有结果（含空/失败行）落盘 `results/grid/`，不做选择报告。
- **方法修正全进** `docs/CHANGELOG.md`；README 只呈现当前方法，不回收旧结果，不写失败编年史。
- **不移动门柱**：敏感性层不改变四态主分析的门柱。

## 8. 期刊定位（定投 AJHG，2026-08-16 定案）

- **定案**：**定投 AJHG（IF 7.7，Article 类型）**。对标样稿 Ray et al. 2025 "Single-cell transcriptome-wide MR and colocalization... ASCVD"（AJHG 112(7):1597–1609）——方法学骨架一致（cis-MR + coloc + drug target），方向在 AJHG scope 内。投稿准备指南见 `docs/AJHG_SUBMISSION_GUIDE_20260816.md`。其他期刊（eBioMedicine/Diabetologia/HMG/NC）推测性定位已删除。主叙事为三股线：全量图谱（31,373 对，FDR-core **982 → 121 strong**）+ **操作特性漏斗**（MR p 单调校准 0.71%→25.6%；**FDR-core 12.32% 与 stage-2 grid 12.96% 收敛**，首次量化"MR 筛药靶"的精度-召回权衡）+ **15 个 catalog 未报道新候选**（9 已知位点候选效应基因 / 6 无 catalog 弱候选，GTEx 6/7 + FinnGen 9/9 基因级一致/覆盖 9/15 独立复现）。2026-08-15 评审版"需要新因果发现才够档"的门槛已部分跨越（15 个候选即新发现型阳性），但 6 个弱峰候选不作主表头条、候选不作因果宣称。
- **现状定位（诚实，2026-08-16 更新）**：转录通道从"多层收敛图谱型阳性（无新发现）"升级为**图谱型 + 新候选发现型**——全量扫描（build 修复后 31,373 对）**FDR-core 982 对 MR-significant → 121 个 strong**（另有 2 个灰区 AP3S2×T2D / ZNF19×CAD）：106 已知位点（KCNJ11/PCSK9/APOC3 等仍是校准案例）+ **15 个 catalog 未报道新候选**（9 已知位点候选效应基因 / 6 弱候选），GTEx 独立 eQTL 方向 6/7 一致，FinnGen R11 独立队列复现可定位子集 9/9 基因级一致（覆盖 9/15）。方法学主结论从单一"coloc 一致性率 12.9%"升级为**操作特性漏斗**（MR p 单调校准 0.71%→25.6%；FDR-core 12.32% 与 stage-2 grid 12.96% 收敛，grid 富集 4.3× 但丢 15 候选）。**诚实不变项**：候选≠因果（6 弱候选不作主表头条）；"catalog 未报道"≠文献未报道；coloc 弱峰 caveat（仅 41/106 已知强候选 GWAS 峰 p<5e-8）必须披露；121 调用无全量尺度置换标定；MR 单工具 Wald ≡ SMR 数学恒等；蛋白通道窄（PCSK9/APOC3 仅作补充）。
- **具体执行路线、成本与预期收益**：`docs/archive/IMPROVEMENT_STRATEGY.md`（2026-08-13 多 agent 评审 + 对抗验证产物）。

---

## 版本记录

- **v0.14（2026-08-16）**：**FDR-core 口径统一（审稿修订 + 内容完整性验证收官）**——分结局 BH-FDR q<0.05 主口径：982 MR-sig / 121 strong（12.32%）/ 15 候选（9 known-locus + 6 弱；8 个 raw strong 掉出）；yield 12.32% 与 stage-2 grid 12.96% 收敛同报、召回 98.4%/87.6% 同报；Table_S1/S2 + cover_letter + PAPER_WRITING_PLAN（含故事/期刊/审稿输入）+ INTEGRITY_AUDIT + AJHG_GUIDE + FACTS 全部同步 FDR-core；新增 `docs/SCHEMA_20260816.md`、`results/strong_all_subset_20260816.csv`（123 = 121 FDR-core + 2 灰区）、README「复现指南」流水线表；Word DOCX 重建（无 title page、Data availability 后置、1 表 9 图）；版本号 v0.13→v0.14。
- **v0.13（2026-08-16）**：全量扫描 M20–M24 + 新候选发现 M25/M25b + GTEx 复现 M26 + 精度漏斗 M27 + FinnGen 独立队列 M28（已完成：可定位子集 12/12 基因级 / 11/11 变异级方向一致、5 个 FinnGen p<0.05、覆盖 12/21）并入；主结论从"警示性阴性"升级为"图谱资源 + 操作特性 + 23 新候选"；摘要加 v0.13 段；§0 故事表 R1/R2 更新；§4.1 加全量扫描；新增 §4.7 新候选与复现；§2.2/§5 闭合 T2D N；§5 局限 7 改"已修"；§6/§8 更新定位（eBioMedicine 冲高判定不成立）；版本号 v0.12→v0.13。
- **v0.12（2026-08-13）**：摘要加深化核查段 + 新增 §4.5 深化核查表 / §4.6 结果文档索引；§6 改为 A/B/C 三档下一步；§8 更新定位（置换标定 FP=1.45%、HEIDI 48.6%、GTEx 69.8%、一致率口径 106/819=12.9%）。〔注：2026-08-16 口径统一后 precision 修正为 106/818=12.96%，见 v0.13+〕
- **v0.11（2026-08-13）**：蛋白通道 M3/M5 完成并入 README（蛋白 MR 主表 + 空结果 + 共定位 strong=2）；删除 deCODE"下载中"、watcher/cron/心跳、蛋白通道"待数据"等过期内容；四态分类加"工具可得性"诚实修正；数据表同步（INTERVAL 已测）；**首读统计量修正（nsnp≤3→FE，MRE SE 塌缩伪影，多 agent 评审）**；面向论文叙事重构（对称两通道 + 诚实局限 + 下一步 + 冲高路径）。
- 历次版本、方法修正与真实运行轨迹见 `docs/CHANGELOG.md`。
