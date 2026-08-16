# 论文写作方案（含故事·期刊定位·审稿人输入）— cis-MR×coloc 图谱 → AJHG 单篇（2026-08-16）

> 本文档 = 论文开发路线图汇总（4 个文档合并；原件已归档 `docs/archive/`）：
> - **正文 §0–§9** = 写作方案（原 `docs/PAPER_WRITING_PLAN_20260816.md`）
> - **附录 A** = 最终故事三股线（原 `docs/archive/FINAL_STORY_20260816.md`）
> - **附录 B** = 期刊定位定案（原 `docs/archive/JOURNAL_TARGETS_20260816.md`）
> - **附录 C** = 审稿人视角改进建议（原 `docs/archive/REVIEWER_IF_ADVICE_20260816.md`）
> 唯一数字权威 = `docs/FACTS_20260816.md`；数据列字典 = `docs/SCHEMA_20260816.md`。

---

## 正文：论文写作方案 — cis-MR×coloc 图谱 → AJHG 单篇（2026-08-16）

> 面向：本项目学生 + 导师 agent。**单一篇论文、单一完整故事、诚实优先**，**定投 AJHG（IF 7.7）**——已定案，不再评估其他期刊（eBioMedicine/Diabetologia/HMG 等推测性期刊定位已删除，理由见 `本文档附录 B` 尾部定案记录）。
> 故事定型（2026-08-16 FDR-core 定稿，勿再改）：**三股线**——① **图谱资源**（31,373 对全转录组 cis-MR×coloc 底表，**FDR-core 982 → 121 strong** + 2 灰区）② **方法学操作特性**（MR p 阈值单调校准 strong 率 0.71%→25.6%，**FDR-core 12.32% 与 stage-2 网格 12.96% 收敛**，首张"MR 筛药靶"精度-召回曲线）③ **15 个 catalog 未报道候选**（9 已知位点候选效应基因 + 6 弱候选，GTEx 6/7 + FinnGen 9/9 双独立复现）。
> **主口径 = 分结局 BH-FDR q<0.05（FDR-core）**；raw p<0.05（4,248 对）仅作名义漏斗。唯一数字权威 = `docs/FACTS_20260816.md`。
> 本方案所有数字已从 `results/feasibility_20260815.md`、`results/coloc_full_summary_20260815.csv`、`results/m22b_window_fix_20260815.csv`、`results/m25_new_strong_annotation_20260816.csv`、`results/m25b_reclassify_20260816.csv`、`results/m26_gtex_replication_new23_20260816.csv`、`results/m27_precision_funnel_20260816.csv`、`results/m28_finngen_replication_new23_20260816.csv`、`results/fdr_core_20260816.csv`、`results/m36b_summary_20260816.csv` 逐项核实并经脚本重算。

---

## 0. 数字基准（全文唯一口径，先统一，防止正文自相矛盾）

论文最容易翻车的地方是**同一个词在不同段落用了不同分母**。以下为全文唯一允许出现的关键数字，写进正文/摘要/图注前逐字照抄：

| 量 | 数字 | 出处 |
|---|---|---|
| 扫描网格 | **49,866 对**（16,622 基因 × 3 结局 T2D/CAD/FBG） | grid |
| coloc QC 通过 | 31,371 / 31,373 | M24 |
| **FDR-core MR-significant（主口径）** | **982**（t2d 394 / cad 576 / fbg 12） | M36b |
| MR-nominal 显著集（raw p<0.05，仅作名义漏斗） | 4,248 | M24 |
| MR 灰区（0.05≤p<0.5）/ 阴性（p≥0.5） | 14,294 / 12,829（raw 口径分层） | M24 |
| **strong coloc（FDR-core）** | **121**（t2d 65 / cad 54 / fbg 2）；灰区 MR-null 2（AP3S2×t2d / ZNF19×cad）；阴性 0 | M36b |
| **主 yield（FDR-core）** | **121/982 = 12.32%** [Wilson 10.41–14.53%] | M36b 定案头条 |
| stage-2 网格 yield（收敛参照） | 106/818 = **12.96%** [Wilson 10.83–15.43%]（819 中 TCF19×CAD 因 GWAS 区域数据不可用无法评估 coloc → 可评估 818 对） | M5/M27 |
| atlas 名义 yield（raw） | 129/4,248 = **3.04%** [Wilson 2.56–3.60%]，仅作名义漏斗 | M24 |
| 精度曲线（按 MR 显著性分层） | **FDR-core 12.32%** → 名义 3.04% → 灰区 0.014%（2/14,294）→ 阴性 0%（0/12,829） | M24/M36b |
| **召回（FDR-core 层面）** | **121/123 = 98.4%**（FDR-core strong 全在 MR-sig 集内；2 灰区在外） | M36b |
| 召回（已知 validated 层面，保守） | 106/121 = **87.6%** | 正文建议同报 |
| **MR 显著集外 strong** | **2/27,123 ≈ 0.007%**（灰区 2；阴性 0；单侧 95% 上界 0.0232%/对——AP3S2×t2d 灰区 GWAS峰5.5e-11；ZNF19×cad 灰区 GWAS峰2.9e-11） | M24/M36b |
| M20 抽样验证 | MR 显著集外 6,000 对 strong = 0（抽样遗漏 AP3S2/ZNF19，与全量不矛盾） | M20 修复版 |
| 已知重现 | **106/106**（pp4 最大差 2.8e-5，逐位一致） | M24 |
| FDR-core 新增 strong | **15** = **9 已知位点候选效应基因 + 6 无 catalog 弱候选**（M25/M25b + M36b；候选≠因果，6 弱候选不作主表头条；raw strong 被 FDR 掉出 8 个——LAMC1 等——不作为候选） | M24/M25/M36b |
| 精度漏斗（操作特性曲线，M27） | MR p 阈值**单调校准** strong 率：0.5→0.71%、0.05→3.04%、0.01→6.56%、0.001→14.8%、1e-4→24.6%、1e-5→25.6%（首张"MR 筛药靶"精度-召回曲线） | M27 |
| GTEx 15 候选方向（M26，eQTL 侧独立） | **6/7 可测一致**（1 冲突 VSIG8，如实报告；8 不可测） | M26 |
| FinnGen 15 候选复现（M28，结局侧独立队列） | 9 个可检（6 个 lead 未在 FinnGen R11 定位）：**基因级 9/9、变异级 8/9 方向一致，0 冲突**；FinnGen 自身 p<0.05 **4 个**（RBM6 6.9e-06、CNNM2 6.3e-04、CD101 1.2e-03、RIC8A 1.1e-02）；**对齐覆盖 9/15（60%）** | M28 |
| PP.H4 阈值敏感性（M36b，FDR-core 集内） | ≥0.9：**74/982 = 7.5%**；≥0.8：**121/982 = 12.3%**；≥0.5：**290/982 = 29.5%**——**单调性不随阈值改变**；mr_p≥0.05 集外 ≥0.8 仅 2、≥0.5 仅 11（全在灰区 0.05–0.5） | M36b |
| atlas 121 FDR-core strong 的 LD 独立簇（A8，r²≥0.8；M30 在 129 nominal 集上算，121 子集已重算） | **≈104 簇**（13 个多 SNP 簇 + 5 个未在 1kg EUR 定位按独立保守计；含灰区 2 → atlas 123 独立位点约 106）；CNNM2/CWF19L1/PLAUR 与已知 strong 共簇（"同一已知信号上提名新效应基因"）、TAGLN2/CCDC19/VSIG8 共用 rs2789422（cis 多基因共享） | M30 `results/ld_clustering_131_20260816.md` |
| GWAS Catalog 命中（build 校正后） | **86/106 = 81%**（26 rsID 直命 + 59 gene±100kb hg38 + 1 基因本体） | M22b |
| **超出已报道注释（新候选上限，106 集口径）** | 20/106 = **19%**（M22b 已加修正注）；**论文口径以 M36b FDR-core 为准**——15 个新增 strong 已逐位点分类（9 已知位点候选效应基因 / 6 无 catalog 弱候选），不再用 19% 上限表述 | M22b/M25/M36b |
| 置换零假设 FP | PP.H4≥0.8 = **1.45%**，观测 strong 富集 ~69× | E1（全扫描级重推见任务 A9） |
| HEIDI | 全 MR 显著集 48.6%（398/819）；strong 子集 **71.7%**（76/106） | P3 |
| Steiger 方向 | 73/76 = **96.1%** eQTL→outcome 正向 | P5 |
| GTEx 跨组织 | 44/63 = **69.8%** 方向一致；非全血 26/40 = 65%；同变异 6/6 = 100% | M15 |
| FinnGen 复现 | 整体 92.6%；**双显著子集 46/46 = 100%**（T2D 27/27、CAD 19/19） | M7 修复 |
| eQTL 功效 | strong lead |Z| 中位 16.4 vs 全量 16.3，Mann-Whitney **p=0.45**（不分辨） | M22 |
| 蛋白层覆盖 | UKB-PPP 仅 10/76 → **蛋白层扩展不可行**（写 Limitation） | M21 |
| 编码版一致率 | 76/637 = 11.9%（≈ 全量 12.32%，非编码不稀释） | M5 |
| p12=1e-6 敏感性 | 20/106 = 18.9% 保持 PP.H4≥0.8 | 敏感性，写 R6 |
| 药物注释 | 仅 KCNJ11（14 独特化合物）/BMPR1A 有已知药物；74/76 未成药 | P4 |
| LD 独立簇 | 92（保守上限；位置聚类佐证 89–88） | M5/M18 |

**术语纪律（全文强制）**：
- **"FDR-core 集"** = 982（分结局 BH-FDR q<0.05，预注册主口径）。"MR-nominal 集" = 4,248（grid 单IV p<0.05，仅作名义漏斗）。"stage-2 网格" = 819（LD clump + IVW + 方向复核存活，收敛参照）。三个词在 Methods 定义一次，正文不得混用。
- "coloc yield（正文首选；原名 precision/concurrence rate）" = P(strong coloc | MR 显著)，是**共定位支持率**，不是 MR 假阳性率。
- "召回" = P(MR 显著 | strong coloc)，在 atlas 枚举上才有意义（coloc 在全量跑，不只 MR 存活者）。

---

## 1. 候选标题（3–5 个，含副标题方向与理由）

> AJHG 审稿人对"承诺过强"的标题零容忍。**所有候选规避 novel / discovery / new causal 字样**。

1. **首选（方法学头条）**
   "Colocalization is the exception in single-instrument cis-Mendelian randomization: operating characteristics of a transcriptome-wide enumeration across type 2 diabetes, coronary artery disease, and fasting glucose"
   - 理由：方法学核心直接进标题；"exception"精确点出反直觉点（~87–97% MR 信号不获共定位支持）；"operating characteristics"是审稿人认的方法学术语；"enumeration"强调全量而非抽样；无任何承诺新发现的词。

2. **数字钩子**
   "Fewer than one in seven cis-MR signals survive regional colocalization — and colocalized loci are rarely missed: a 49,866-pair atlas for three cardiometabolic outcomes"
   - 理由："one in eight"（12.32% ≈ 1/8.1）让人过目不忘；破折号后半句直接否定"MR 召回低"叙事，一句话讲完两个结论。

3. **2×2 框架（推荐备选）**
   "Low precision, high recall: exhaustive operating characteristics of single-instrument cis-Mendelian randomization for drug-target screening, with a transcriptome-wide colocalization atlas"
   - 理由："Low precision, high recall"四个词装下核心结论，审稿人 10 秒内抓住贡献；drug-target screening 点明应用场景；atlas 标为副产品身份。

4. **资源导向（若编辑偏好图谱/资源型）**
   "A transcriptome-wide cis-MR and colocalization atlas for type 2 diabetes, coronary artery disease, and fasting glucose: operating characteristics and candidate effector genes"
   - 理由：以社区资源为主体（图谱交付物），操作特性隐含其中；106/106 可复现是差异化卖点（先例研究做不到全量复现审计）。

5. **警示型（吸引非方法学读者）**
   "What does single-instrument cis-MR actually screen? Most MR hits lack shared-causal-variant support while genuinely colocalized loci are rarely missed"
   - 理由：问句式标题把"工具不能再裸用"的警示前置；适合编辑偏方法论时用。

**建议**：定投 AJHG 推 **1（方法学头条）或 3（低精度高召回框架）**；若编辑偏好资源/图谱型，换 4。**禁用**：任何含 novel/discovery/new causal 的组合。

---

## 2. 核心叙事弧

### 一句话 takeaway（要让审稿人记住的）
> **"单工具 cis-MR 是低精度、高召回、单调可校准的筛查：raw 名义口径下 MR 显著集只有 3.04% 获得区域共定位支持（应用预注册分结局 BH-FDR 后核心集收敛到 12.32%，与 stage-2 网格 12.96% 一致），却覆盖 98.4% 的共定位位点——工具的病是假阳性多，不是漏检。在这张 31,373 对的全转录组图谱上，FDR-core 121 个 strong coloc 中 15 个是 GWAS Catalog 未报道的候选效应基因（9 个位于已知 T2D/CAD 风险位点内），并用 GTEx + FinnGen 两个独立队列做方向复现。"**

### 反直觉点在哪（4 个，按冲击力排序）
1. **领域默认担心"MR 漏检"（召回），数据说反了**。既往所有 coloc 都只在 MR 显著子集（819）上跑——没有召回维度。M23 全量枚举第一次补上：MR 不显著时 strong coloc 几乎不存在（2/27,123 ≈ 0.007%，且这 2 个的 GWAS 峰均基因组显著）。"MR 召回低"被实证否定。
2. **MR 显著 ≠ 因果信号，比例首次被量化**。3.04%（名义）~12.32%（FDR-core），绝大多数 MR 命中是"MR 显著但区域不共定位"。yield 曲线（M27/M36b）单调校准 0.71%→25.6%——首次给"MR 筛药靶"贴操作特性曲线。
3. **区分维度是 GWAS 信号强度，不是 eQTL 强度**。eQTLGen lead 全极强（|Z| 中位 16.4 vs 全量 16.3，p=0.45 不分辨），所以"弱工具导致不共定位"在 MR 显著集内不成立；真正决定共定位支持的是 GWAS 峰是否显著。这把"换更强 eQTL 就能救 MR"的常见幻想直接证伪。
4. **行业标准"先 MR 后 coloc"网格有系统漏检**。stage-2 预过滤网格（819 对）把 strong 率提到 12.96%（富集 4.3×），却**丢掉 15 个 FDR-core 新候选**——全量扫描（31,373 对）才发现它们。标准流程在"提精"的同时系统性漏掉了位于已知位点内的新效应基因。

### 方法学价值怎么让读者感到"这个工具不能再裸用了"
- **贴第一张操作特性曲线**：coloc yield 按 MR 显著性分层单调下降 12.32%（FDR-core）→ 3.04%（名义）→ 0.014%（灰区）→ 0%（阴性）（Wilson CI 全给）。读者会意识到：跑完 MR 把显著基因直接当靶点，等于 ~88% 的信号在赌。
- **置换标定把账算清——问题在 MR 侧，不在 coloc 侧**：coloc 零假设 PP.H4≥0.8 的 FP 仅 1.45%（观测 strong 富集 ~69×），说明"不共定位"不是 coloc 假阴性灌水，而是 MR 信号的假阳性（多效性/LD 驱动）——正是 coloc 闸门存在的理由。
- **一句话行动建议**（Discussion 收尾）："任何单工具 cis-MR 命中，先查这张图谱（49,866 对已免费预计算）或重跑 coloc 再谈靶点优先级。"资源化 = 结论可执行 = 审稿人觉得这篇有用。

---

## 3. 章节结构

> 期刊格式（定投 AJHG，2026-08-16 定案）：**AJHG 非结构化摘要 ≤250 词**（单段，§5 草稿压缩）。正文 = Introduction → Results → Discussion → **Methods 放正文末尾**（AJHG 惯例），**Methods 细节必须自足**。正文建议 5,000–6,000 词（三股线内容多于 08-13 单线稿）；图注自足。

### Abstract（结构化四段）—— 见 §5 骨架

### Introduction（4 段，约 600–800 词）
- **P1 立论**：cis-MR（尤其单工具 Wald，或 SMR）是药物靶点筛查的事实标准；隐含假设 = "MR 显著 ≈ 可定位的因果信号"。引出两个从未被系统量化的操作特性：精度（MR 命中有多少获共定位支持）、召回（真共定位位点有多少被 MR 筛出）。
- **P2 缺口**：既往 SMR+coloc 类工作流只在 MR 显著子集上跑 coloc → 只有 yield 没有召回；yield 是单点事后估计（我们的前序工作 12.32%），无 2×2、无置换、无分层。领域对"MR 漏检"的担忧一直停留在假设层面。
- **P3 预期（可证伪）**：提出反转假说——MR 的问题在假阳性（精度）而非漏检（召回）；区分维度在 GWAS 侧而非 eQTL 侧。预设置换与分层检验。
- **P4 本研究**：49,866 对全量 coloc 枚举 → 操作特性 + 图谱资源 + 独立验证层（HEIDI/Steiger/GTEx/FinnGen）+ 预注册 + 全量复现审计。

### Materials and Methods（M1–M10，每节标注"已有 / 待补"）
- **M1 数据源**：eQTLGen 全血 cis-eQTL（n=31,684，16,923 基因）；OpenGWAS T2D（GCST006867，**n 口径待核** A7）/CAD（GCST005194，n=296,525）/FBG（GCST005186，n=58,074）。【已有，A7 补 T2D 口径】
- **M2 MR 网格**：cis ±1Mb、p<5e-6、EUR LD clump r²<0.01@1000kb、单工具 Wald → nominal 4,248；stage-2 本地 clump+IVW 复核 → 819（收敛参照）；FDR-core 982 = 分结局 BH-FDR q<0.05（预注册主口径）。**在此声明 MR 单工具 Wald ≡ SMR 为数学恒等**（同一 top cis-SNP 同一比率，b 到 6 位小数），非独立方法。【已有】
- **M3 coloc**：coloc.abf（p12=1e-5，敏感性 1e-6），PP.H4≥0.8 strong / ≥0.5 moderate；**rsID 匹配对齐（坐标 build 无关）**——这是针对曾发生的 hg38/hg19 错位的正面声明；QC nsnp≥10；回文保守排除。【已有，rsID 对齐已是修复版标准】
- **M4 全量扫描与 2×2 定义**：31,373 对 QC 通过；MR 状态四层（FDR-core 982 / nominal 4,248 / grey 14,294 / null 12,829）；**coloc yield = P(strong coloc | MR 显著)，明确写"共定位支持率，非 MR 假阳性率"**；召回定义 + Wilson CI。【已有，A2 补四层 Wilson CI】
- **M5 置换标定**：打乱 eQTL beta 跨 SNP（保留 GWAS 信号与 LD/MAF 结构），零假设 PP.H4≥0.8 FP 率。**A9 必须把置换重推到全扫描级**（E1 只在 106 个 MR 显著区域上标定，评审明确要求全扫描重导）。【E1 已有，全扫描级待补】
- **M6 验证层**：HEIDI（SMR+HEIDI 全 819）；Steiger 方向；GTEx v8 6 组织方向一致性；FinnGen 对齐（逐行）。【已有】
- **M7 GWAS Catalog 注释**：build 校正（rsID→hg38 偏移中位数转换，std<1bp），GenePos±100kb 主口径 + top_snp 锚点 + 三口径敏感性。【已有 M22b】
- **M8 eQTL 功效分析**：106 strong 的 lead |Z| 百分位 vs 全量 lead；Mann-Whitney。【已有 M22】
- **M9 可复现性审计**：单基因收敛（GIT1 对齐 1917 变异、PP.H4=0.8937，本地 full 与 API 双源逐位一致）；106/106 pp4 重现（最大差 2.8e-5）；坐标 build 独立声明。【已有】
- **M10 软件/数据可用性与预注册**：coloc 5.2.3 / R env r-mr / plink 1.9 / SMR 1.3.1；图谱 CSV（31,373 对全量 + 123 strong 子集 = 121 FDR-core + 2 灰区）+ 脚本 GitHub/OSF 发布；预注册哈希；版本声明（M23–M28 分析脚本 + results CSV 一并入库，本文档附录 A/README 已同步）。【已有，发布打包属 P2】
- **M11 15 候选的独立复现（线 3 证据，M26/M28）**：**(a) eQTL 侧 GTEx v8**（M26）：对 15 个新候选逐基因取 6 组织（全血在内）cis-eQTL 效应方向，与 eQTLGen lead 方向比对 → 6/7 可测一致、1 冲突（VSIG8）、8 不可测（无 GTEx 显著 eQTL，如实报告）；**(b) 结局侧 FinnGen R11**（M28）：15 候选的 eQTLGen lead SNP 对齐 FinnGen R11 全量 sumstats（(chrom,pos) 主匹配 + rsids 正则兜底，等位基因 a1 对齐），变异级比 sign(β_finn) vs sign(β_orig)、基因级比 sign(mr_b) vs sign(β_finn×eQTL Z）→ 9/15 可检、基因级 9/9 一致、4 个 p<0.05（RBM6/CNNM2/CD101/RIC8A）、对齐覆盖 9/15；fbg 结局跳过（FinnGen 无空腹血糖表型，且 15 候选均无 fbg）。【已有】
- **M12 精度漏斗（操作特性曲线，M27）**：对 MR-nominal 显著集按 MR p 阈值分层（<0.5/<0.05/<0.01/<0.001/<1e-4/<1e-5），逐层计算 strong coloc 率 + n，验证单调校准 → 首张"单工具 cis-MR 筛药靶"精度-召回曲线；与 R2 四层 2×2 曲线（12.9/0.67/0.014/0%）互补。【已有】

### Results（R1–R7，按三股线组织）
- **R1 图谱构建与全量枚举（线 1·图谱资源）**：31,373 对 QC 通过（49,866 候选 → QC 31,371）；**FDR-core 982 MR-sig → 121 strong**（t2d 65/cad 54/fbg 2；另有灰区 2 AP3S2/ZNF19）；分层主表（FDR-core 982 / nominal 4,248 / grey 14,294 / null 12,829）。【Fig 1, Table 1】
- **R2 操作特性（线 2·方法学核心）**：主 yield **12.32%**（FDR-core 121/982，Wilson 10.41–14.53%）与 stage-2 网格 **12.96%**（106/818，10.83–15.43%）**收敛同报**；名义 3.04%（raw）仅作漏斗参照；召回 98.4%（FDR-core 121/123）与 87.6%（106/121）同报；**精度漏斗（M27/M36b）**：MR p 阈值单调校准 strong 率 **0.71%→25.6%**（0.5→0.71%、0.05→3.04%、0.01→6.56%、0.001→14.8%、1e-4→24.6%、1e-5→25.6%），每层附 n。【Fig 2/3】
- **R3 共定位的决定维度**：eQTL |Z| 不分辨（p=0.45）；GWAS 峰显著区分（41/106=38.7% 已知 strong 的 GWAS 峰 p<5e-8，**弱峰 caveat 如实披露**）；**AP3S2/ZNF19 两个 MR-阴性 strong 案例面板**（A5 实测：AP3S2×t2d gwas_min_p=5.5e-11 / eqtl_F_max=55.6 / MR p=0.38 灰区 / PP.H4=0.935；ZNF19×cad gwas_min_p=2.9e-11 / eqtl_F_max=7.13 **弱 eQTL** / MR p=0.33 / PP.H4=0.903 → "真实位点上单工具 MR 功效不足"，非新发现；**ZNF19 eQTL 弱仍 strong coloc = 决定维度在 GWAS 侧，直接佐证 R3**）。【正文 vignette】
- **R4 15 个 catalog 未报道候选（线 3·发现核心）**：FDR-core 全量扫描比 stage-2 网格多 **15 个 strong**（FDR-core 内、不在已知 106、GWAS Catalog T2D/CAD 无记录）；**9 已知位点候选效应基因 + 6 无 catalog 弱候选**（M25/M25b 逐位点注释 + M36b FDR-core）；raw strong 被 FDR 掉出 8 个（LAMC1 等）不作为候选；弱信号膨胀 caveat（6 弱候选 GWAS 峰多 >5e-8，不作主表头条）。【Fig 5, Table S2】
- **R5 双独立复现（线 3·证据加固）**：GTEx v8 eQTL 方向 **6/7 可测一致**（M26，1 冲突 VSIG8，如实报告）；FinnGen R11 独立队列 **可定位子集方向全一致**——基因级 **9/9**、变异级 **8/9**，0 冲突；FinnGen 自身 p<0.05 **4 个**（RBM6/CNNM2/CD101/RIC8A）；**对齐覆盖 9/15（60%）诚实披露**（6 lead 未定位）。（M28，芬兰人群，结局侧复现）。【Fig 5 复现栏, Table S2】
- **R6 收敛验证（线 1/2 佐证）**：已知 **106/106 重现**（审计背书）；置换 FP 1.45%（MR 显著子集口径，全量尺度未标定→Limitation）；HEIDI 71.7%（strong 子集）；Steiger 96.1%。证明 coloc 支持的位点真实，非方法噪声。【Fig 6】
- **R7 稳健性与资源交付（线 1·收尾）**：p12=1e-6 敏感性（20/106=18.9%）；PP.H4 阈值敏感性（FDR-core ≥0.5/0.8/0.9 = 29.5%/12.3%/7.5%）；fbg 弱信号膨胀（2 strong 峰不显著，MT3 教训）；**图谱底表以开放 CSV 交付**（`results/coloc_full_{t2d,cad,fbg}_20260815.csv`，31,373 对含 mr/coloc 全字段 + `results/fdr_core_20260816.csv` 982 + `results/strong_all_subset_20260816.csv` 123）。【补充表 S1】

### Discussion（5–6 段）
- D1 主发现重述：低 yield 非低召回；工具的病是假阳性多（12.32%/3.04% 名义 vs 召回 98.4%）。
- D2 为什么领域以前以为召回是问题（coloc 只在 MR 显著子集上跑的路径依赖）；全量枚举如何改写认知；两个例外（AP3S2/ZNF19）不是新发现而是"真实 GWAS 位点上单IV 功效不足"。
- D3 方法学含义：区分维度在 GWAS 侧（反驳"换更强 eQTL 救 MR"）；置换标定把不共定位归因到 MR 侧；MR≡SMR 恒等 → 单点估计自洽 ≠ 独立验证，真正的区域级独立检验只有 coloc + HEIDI。
- D4 与文献关系：现有 SMR/coloc 工作流与图谱；**106/106 全量复现作为可复现性审计背书**（别家做不到）。
- D5 对实践启示：coloc 闸门强制化 + **yield 漏斗当"操作手册"**（正文给一张查表：查 MR p 阈值 → 期望 coloc 支持率，直接可执行）；49,866 对图谱免费预计算；15 个候选（9 已知位点新效应基因 + 6 弱候选）按证据分级给验证路径（不承诺发现）。
- D6 边界定位：候选均为 hypothesis-generating、纯计算、单组织全血、3 结局——这些是定位不是失败（衔接 Limitations）。

### Limitations（独立小节，审稿人必查，11 条，每条一句 + 影响）
1. 单工具 Wald 扫描层为 hypothesis-generating；结论级以 stage-2 819 集为准。
2. eQTLGen 全血单组织；肝/胰岛特异基因（PCSK9 等）结构性盲；InsPIRE 胰岛 eQTL MR 10/11 反向 + KCNJ11 胰岛 null（负性核查）。
3. **蛋白层扩展不可行**（UKB-PPP 覆盖 10/76）→ 无介质层声明；"双共定位 = 双介导"逻辑不成立（mRNA→蛋白 flux）。
4. coloc.abf 单方法；coloc-SuSiE 校准失败（PCSK9 PIP=0）；p12 敏感性仅 20/106。
5. eQTLGen–GWAS 样本重叠偏倚使 MR 假阳垫进分母 → **精度偏保守**（不是被高估）。
6. 1000G LD 面板 haploid 化（r² 轻微低估）；92 簇为保守上限。
7. GWAS Catalog 注释不全 → 19%（已知 106 口径）是上限非真值；**15 候选（FDR-core 口径）是"MR 显著集内新增 strong"逐位点分类结果，两口径不混淆**（19% 表"106 已知位点未入 catalog 上限"，15 表"FDR-core 全量扫描比 stage-2 网格多发现的 strong"，已在 §0 表注释区分）。
8. coloc yield ≠ MR 假阳性率（coloc 也可能漏真，二者无法在单点分解；全量 2×2 + 置换缓解）。
9. 仅 3 个心脏代谢结局；已知 106 层面 FinnGen 复现仅双显著子集（M7）；**15 候选结局侧复现（M28）为单结局（T2D/CAD）、FinnGen 无 FBG 表型故 fbg 跳检（且 15 候选均无 fbg）**。
10. 纯计算、无正交湿实验验证。
11. **M28 复现 = 方向（符号）复现，非效应量复现**；FinnGen 芬兰人群与 eQTLGen/OpenGWAS 欧洲人群存在祖先差异（对等位基因对齐与 sign 一致性为保守测试）；GTEx 侧 6/7 仅 eQTL 方向、8 个候选在 GTEx 无显著 eQTL 而不可测。

---

## 4. 图表规划

> HMG 风格：主图 5–6 张、主表 2–3 张、补充图/表不限。图注必须自足（含数字与定义）。所有图从最终 CSV 重画，不手工改数。

### 主图（7 张，对应 R1–R7；F1/F2 已实画并定稿）
- **Fig 1 设计 + 漏斗（R1）**：左 = 研究设计示意（eQTLGen × 3 GWAS → 31,373 对 → MR → coloc → 验证层，标注预注册）；右 = 漏斗 49,866 → 31,371 QC → 982 FDR-core → 121 strong（+ 2 灰区 = 123；含 15 候选高亮）。**图注必须并排写两个分母**（FDR-core 982 vs nominal 4,248）。
- **Fig 2 操作特性主图（最核心，R2）**：(a) 2×2 概念图（MR 状态 × coloc 状态四格数字，123/31,373 全量）；(b) **coloc yield 曲线**：MR 显著性四层 × strong 率，Wilson CI 误差棒（12.32 / 3.04 / 0.014 / 0）；(c) **yield 漏斗（M27，实画 F1 定稿）**：MR p 阈值单调 0.71%→25.6%，每层附 n；(d) 召回条形（121/123=98.4% FDR-core 全图谱、106/121=87.6% 已知覆盖）。——这张图就是评审要的"yield 曲线而非单点"，漏斗给可查的"操作手册"。
- **Fig 3 PP.H4 全分布（R2/R3）**：31,371 对 PP.H4 累计分布按 MR 状态分层（sig/grey/null），显示 sig 层长尾、非 sig 层全塌在 0 附近；叠加 PP.H4=0.8 阈值线。
- **Fig 4 图谱染色体图（曼哈顿式，R1）**：123 strong（121 FDR-core + 2 灰区）按染色体分布，颜色按结局（t2d/cad/fbg），形状区分 known / 15 候选 / 灰区（AP3S2、ZNF19 高亮）。
- **Fig 5 候选基因 + 双复现（R4/R5，实画 F2 定稿）**：15 基因 lollipop（实心 = known_locus / 空心 = no_catalog），右栏 **GTEx 方向**（6 ✓ / 1 × / 8 ·）+ **FinnGen 方向**（9 ✓ / 0 × / 6 ·，9/15 可检）；图例在图下。Table S2 逐行对读。
- **Fig 6 收敛验证面板（R6）**：(a) HEIDI 通过率（全集 48.6% vs strong 71.7%）；(b) Steiger 方向 96.1%；(c) GTEx 跨组织方向一致率（44/63、非全血 26/40、同变异 6/6）；(d) 已知 106 层面 FinnGen 双显著子集散点（27/27、19/19 高亮）。
- **Fig 7 资源底表快照（R7）**：results/coloc_full_{t2d,cad,fbg}_20260815.csv 列结构 + 行数（31,373）+ 一个示例条目，标注开放协议——让审稿人一眼知道"资源长什么样"。

### 主表（3 张）
- **Table 1（审稿人第一眼看的表）操作特性 2×2**：行 = MR 状态四层 + 三结局分列；列 = n pairs / strong coloc / rate / Wilson CI；表注定义 coloc yield 与 recall；附 M27 漏斗层（阈值 × strong 率 × n）。
- **Table 2 图谱汇总**：strong 计数按结局 × 类别（known 106 / 15 候选：9 known-locus + 6 弱）；**独立位点列（A8/M30，FDR-core 121 子集重算）：r²≥0.8 约 104 独立簇**（13 多 SNP 簇 + 5 lead 未定位按独立保守计；含灰区 2 → atlas 123 独立位点约 106）；编码 vs 非编码；catalog 命中 vs 未命中。
- **Table 3 15 个候选表**：gene / outcome / top SNP / PP.H4 / MR p / GWAS min p / eQTL |Z| / 分类（9 known-locus vs 6 弱）/ GTEx 方向 / FinnGen 方向 + FinnGen p。放正文（诚实展示），6 弱候选行标注"低置信，不作头条"。

### 补充（S1–S9）
S1a 全部 123 strong（121 FDR-core + 2 灰区，标注方向/复核状态）；S1b **6 个弱候选明细**（无 catalog 区域 + GWAS 峰 p，低置信声明）；S2 敏感性（p12=1e-6、PP.H4≥0.9/0.5 → `results/s2_pph4_sensitivity_20260816.csv`）；S3 置换标定（MR 显著子集 E1；**全扫描级数据不可得 → Limitation 如实披露**）；S4 GTEx 逐组织；S5 FinnGen 逐基因对齐明细（15 候选，9/15 可检 + 6 无法定位 lead + fbg 无表型跳检）；S6 eQTL 功效分析；S7 方法验证（rsID 对齐、GIT1 收敛、106/106 pp4 对比）；S8 药物注释/未成药状态；S9 蛋白通道 PCSK9/APOC3 降级材料。

---

## 5. 摘要骨架（Background / Methods / Results / Conclusion 四段式草稿）

> **AJHG 非结构化摘要 ≤250 词**（定案格式）。以下草稿为四段式组织（Background/Methods/Results/Conclusion），供提取内容——定稿时压缩为单段 ≤250 词（数字与 §0 基准逐字一致）。

**Background**
Single-instrument cis-Mendelian randomization (cis-MR), including SMR, is the default screen for prioritizing drug-target genes in cardiometabolic disease, yet its operating characteristics — the proportion of screened loci that represent shared causal variants (precision), and whether genuinely colocalized loci are missed (recall) — have never been systematically quantified because colocalization is only ever applied to MR-significant subsets.

**Methods**
We enumerated 49,866 gene–outcome pairs (16,622 genes × T2D, CAD, fasting glucose) from eQTLGen and GWAS summary data using single-instrument cis-MR, then ran regional colocalization (coloc.abf, PP.H4≥0.8; rsID-based, coordinate-build-independent) on all 31,371 pairs passing QC — a complete 2×2 enumeration of MR status by colocalization status. Colocalization calls were calibrated by permutation, and validated by HEIDI, Steiger, GTEx cross-tissue replication and FinnGen.

**Results**
Applying the pre-registered per-outcome BH-FDR threshold (q<0.05), 982 pairs were MR-significant, of which 121 reached strong colocalization (PP.H4≥0.8): coloc yield 12.32% (Wilson 95% CI 10.41–14.53%), converging with the independently constructed stage-2 clump+IVW grid yield of 12.96% (106/818; 819 total, one without usable GWAS region data). Tightening MR stringency calibrated coloc yield monotonically from 0.71% to 25.6%, while 121/123 (98.4%) of colocalized loci were MR-detectable — yield is low, recall is high. eQTL strength did not discriminate (p=0.45); GWAS signal did. All 106 previously reported loci were reproduced. Exhaustive scanning identified **15 additional strong colocalizations absent from GWAS Catalog** (9 within known T2D/CAD loci — candidate effector-gene nominations; 6 low-confidence); GTEx v8 direction matched in 6/7 testable genes, and FinnGen R11 in 9/9 gene-level testable pairs (4 with p<0.05; alignment coverage 9/15).

**Conclusion**
Single-instrument cis-MR is a low-yield, high-recall screen: most hits lack shared-causal-variant support while colocalized loci are rarely missed. Colocalization gating is mandatory before cis-MR hits are treated as causal targets. We release 15 prioritized candidate effector genes and the full 31,373-pair atlas as community resources.

> 中文提示：摘要里核心数字出处——982/121（FDR-core）、121/982=12.32%（FDR-core yield，Wilson 10.41–14.53%）、106/818=12.96%（stage-2 网格收敛参照）、121/123=98.4%（召回）、0.71%→25.6%（漏斗）、6/7（GTEx）、9/9 基因级 + 覆盖 9/15（FinnGen）。若字数紧张，删 HEIDI/Steiger 细节只留 "validated by orthogonal layers"；"alignment coverage 9/15" 是诚实必要项不可删。

---

## 6. 写作任务清单（按优先级；标注纯写作 vs 缺分析）

### P0 — 写作前必须补齐的分析（✅ = 已被 M25–M28 完成，剩 4 项）
- ✅ **A1 15 个候选逐位点标注**（M25/M25b 在 23 全集上做，FDR-core 子集 15 落表）：GWAS min p / eQTL |Z| / 分类（9 known-locus / 6 弱）全部落 `results/m25_new_strong_annotation_20260816.csv` + `m25b_reclassify_20260816.csv`。Table 3/S1a 数据已齐。
- ✅ **A2 精度曲线 + 漏斗 + Wilson CI**（M27）：四层 rate + CI（12.9 / 0.67 / 0.014 / 0）与 p 阈值漏斗（0.71%→25.6%）在 `results/m27_precision_funnel_20260816.csv` + §0 表。Fig 2 数据已齐。
- ✅ **A3 术语定案**：§0 术语纪律已写死（FDR-core 982 / MR-nominal 4,248 / 灰区 / 阴性 / coloc yield / recall），Methods 引用即可。
- ✅ **A4 全图敏感性表**：p12=1e-6 已有（20/106=18.9%）；PP.H4≥0.9/0.5 在 31,371 对重算完成（≥0.9：nominal 77=1.81%/grey 2/null 0；≥0.5：323=7.60%/grey 9/null 2，单调性不随阈值改变）→ 补充表 S2（`results/s2_pph4_sensitivity_20260816.csv`）。
- ✅ **A5 AP3S2/ZNF19 案例细节**：GWAS 峰 5.5e-11/2.9e-11 + eQTL F（55.6/7.13）+ MR p（0.38/0.33）+ PP.H4（0.935/0.903）全部落 R3 vignette；ZNF19 弱 eQTL 仍 strong coloc = GWAS 侧决定维度佐证。
- ✅ **A6 候选 catalog 核对**：15 个候选已逐位点 GWAS Catalog 核对（M25/M25b 在 23 全集核对，FDR-core 子集 15 落表），"19% 上限"表述已被 15 候选分类取代；Open Targets 补充核对可选。
- ✅ **A7 T2D 结局 n 口径回填**：per-variant N mode≈573,704（README §5 局限 5 定案）。
- ✅ **A8 新增 strong 的 LD 聚类**：atlas 129 nominal-sig strong 重跑完成（M30）→ **r²≥0.8 = 106 独立簇**（14 多 SNP 簇 + 6 lead 未定位按独立保守计）；FDR-core 121 子集重算 ≈ 104 簇（含灰区 2 → atlas 123 独立位点约 106）→ Table 2"独立位点"列（`results/ld_clustering_131_20260816.csv/.md`）。
- ✅ **A9 置换标定全扫描级重推：数据不可得 → 作 Limitation 如实披露**（全扫描 coloc 的 RDS 输入仅覆盖 MR 显著 106 区域，31,371 对无区域级 coloc 重算输入，无法重导全扫描置换；§三 Limitations 第 6/8 条已覆盖"E1 标定于 MR 显著子集"，正文 + S3 注明"全扫描级置换不可行"而非省略）。【不消耗计算资源，写进 Limitation】
- ✅ **A10 图件重画（部分）**：F1（漏斗）+ F2（候选+双复现）已实画定稿（`scripts/M29_figures_20260816.py`）；Fig 3/4/6 待画。【剩余图按 §4 规划逐步画】

### P1 — 纯写作（数据全齐，无分析依赖）
1. Abstract 定稿（§5 草稿压缩为 AJHG 非结构化单段 ≤250 词）。
2. Title 选定（§1 候选 1 或 3）。
3. Introduction（4 段）。
4. Methods 整合（README §3 + PREREGISTRATION + 本方案 M1–M12，AJHG 要求 Methods 自足）。
5. Results R1–R7（逐数字引用 §0 基准，禁自造数）。
6. Discussion D1–D6。
7. Limitations 11 条。
8. 图注与表注（自足式）。
9. 预注册 + 数据可用性声明 + 作者贡献。

### P2 — 发布与投稿（写作之外）
- 图谱 CSV（31,373 行全量 + 123 strong 子集）发布 GitHub/OSF，附 README 列说明。
- 代码发布（M23_full_scan.R / M24_summarize.py / 图脚本），含环境版本锁定（environment_versions.txt）。
- 预注册哈希对照声明（docs/PREREGISTRATION.md + sha256）。
- 投稿 cover letter + 给编辑的方法学卖点段（"首次全量枚举 + 精度曲线 + 106/106 复现审计"）。

---

## 7. 诚实红线（绝不能写 / 怎么写才诚实）

### 绝不能写（审稿人一跑除法/复算即否）
1. **"MR 假阳性率 12.3%"** —— 12.32% 是 coloc yield（MR 显著且获共定位支持的比例），coloc 可能漏真。只能说"仅 12.32% 获区域共定位支持"。
2. **"coloc-only 位点 = 新候选/新发现基因"** —— MR 显著集外 strong 仅 2 个（AP3S2/ZNF19）且均在 GWAS 显著位点，是"真实位点上单工具 MR 功效不足"，非新发现；15 个 FDR-core MR-sig 内新增 strong 是**候选评估**：9 个"已知位点新效应基因" ≠ 全新基因-疾病关联（只能说"该位点已知，此基因未在 Catalog T2D/CAD 关联中标注"），6 个弱候选只是低置信假设（GWAS 峰多 >5e-8，仅进补充表）。
3. **"MR 与 SMR 互为独立验证 / 三方法互证"** —— 数学恒等。真正区域级独立检验只有 coloc + HEIDI。
4. **"MR 召回低 / 漏掉 3/4 共定位位点"** —— 数据正相反（召回 98.4% FDR-core 全图谱 / 87.6% 已知覆盖）。这是旧叙事的残留，见一条删一条。
5. **"双共定位 = 双介导 / 信号经转录层传导"** —— 蛋白层扩展不可行（10/76），mRNA→蛋白 flux 使 eQTL coloc 只是关联回响。
6. **"APOC3×T2D/FBG 被独立救回/复现"** —— 单工具敏感性支持，UKB-PPP 无 APOC3 断头；至多写"待独立平台复现"。
7. **"KCNJ11 是转录介导的 T2D 靶点"** —— KCNJ11 基本不在全血表达；coloc 峰值 rs757110 在 ABCC8。机制留胰岛 eQTL。
8. **"cis-MR 是经过验证的药物靶点发现工具"** —— 必须带"在 coloc 闸门配合下"的限定。
9. **蛋白侧 HEIDI 未校准的主张**（deCODE EAF/LD 错配）—— 蛋白侧 HEIDI 只作探索。
10. **"双复现 = 独立验证 / 因果验证"** —— GTEx（6/7）与 FinnGen（9/9）只说"方向一致"，比符号不比效应量；FinnGen 芬兰人群祖先差异 + 对齐覆盖仅 9/15（60%），只能作支持候选优先级的佐证，不能作因果确认。

### "新候选"怎么表达才诚实（给审稿人的原句模板）
> "Beyond the 106 previously reported loci, exhaustive enumeration (rather than the standard MR-significant pre-filter) identified 15 additional strong colocalizations (PP.H4 ≥ 0.8) within the FDR-controlled MR-significant set. Nine lie within ±100–250 kb of a known T2D/CAD risk locus but nominate an effector gene not annotated in GWAS Catalog for that trait; six lie in regions with no T2D/CAD Catalog entry and are reported as low-confidence hypotheses, most with non-genome-wide GWAS peaks. We emphasize these are prioritized candidates, not claims of causality: GTEx v8 eQTL direction matched in 6 of 7 testable genes and FinnGen R11 matched in all 9 gene-level testable pairs (4 with p<0.05; alignment coverage 9/15), but direction consistency is not effect-size replication and none has orthogonal wet-lab validation. We release all 15 with full evidence columns for community replication."

补充纪律：正文同时写"81% 落在已报道 GWAS 注释内"（已知 106 口径）+ "15 个 FDR-core MR-sig 内新增 strong"（全量扫描比 stage-2 网格多发现口径），两口径在 §0 表分开定义，绝不让任一数字单飞成头条；6 弱候选只进补充表 S1b；GTEx/FinnGen 只写"方向一致"不写"独立验证"。

---

## 8. 审稿人最可能攻击的 3 个点 + 预案

### 攻击 1：yield 分母不一致（"你头条 12.3%，你自己的 Table 1 却还有名义 3.0%——哪个是真的？982 是不是事后挑的？"）
**预案**：
- 两个数字**主动同报**（Abstract + Table 1 + Fig 2 都写），先发制人让审稿人没机会"发现"。
- 明确 982 = 预注册的分结局 BH-FDR q<0.05 集（非事后挑选）；4,248 = grid 单IV 名义集仅作漏斗参照。两者定义在 Methods 一处写死；12.32%（FDR-core）与 12.96%（stage-2 网格 106/818）收敛同报。
- 论证**结论对分母不敏感**：无论 12.32%（FDR-core）还是 3.04%（名义），都是"MR 命中多数不获共定位支持"，方向完全一致——这正是诚实故事的力量。
- 精度曲线（12.9→0.67→0.014→0）本身回应"单点 cherry-pick"：给的是曲线不是单点。
- 15 个 FDR-core 新增 strong 单独报告并标注"低置信弱候选 + 部分 GWAS 峰不显著（弱信号膨胀风险）"，不混进头条分母。

### 攻击 2：共享输入 / 循环论证（"MR 和 coloc 吃同一套 eQTL 和 GWAS 数据，相关是必然，不是操作特性"）
**预案**：
- 强调**不对称性**：问题不在相关，而在方向——coloc-strong 且 MR 不显著 ≈ 0（2/27,123），MR 显著且 coloc 弱 = 87–97%。这个不对称在共享输入下依然有信息量（这是全量 2×2 的意义）。
- coloc 用的是**全区域** cis 变异 + 区域 GWAS（不是 lead IV），信息输入远大于单工具 MR——不是同一个统计量的重述。
- **置换标定**（null FP 1.45%、~69× 富集）证明 coloc 调用是真实信号，不是共享数据的伪相关。
- **正交验证**（HEIDI 71.7%、Steiger 96.1%、GTEx 69.8%、FinnGen 100% 双显著子集）证明 coloc 支持的子集真实，反向坐实"不共定位"是 MR 侧假阳性。
- 明确 recall 是"在此数据集的枚举上"描述性结论，不声称通用敏感性。

### 攻击 3：无新发现 / "so what"（"全部 123 个 strong 都是已知位点或弱信号伪影，没有新基因，这篇贡献是什么？"）
**预案**：
- **重定义贡献**：(a) 第一张单工具 cis-MR 操作特性曲线（coloc yield×召回，置换标定）；(b) 全量 31,373 对图谱资源（社区免费复用）；(c) 106/106 全量复现审计（可复现性背书，别家做不到）；(d) **15 个新候选（9 已知位点新效应基因 + 6 弱候选）+ GTEx/FinnGen 双独立方向复现**。
- 预先写清"这是方法学 + 资源 + 候选发现论文，不是因果发现论文"，把审稿人的"没新基因"问句变成"这个方法学读数 + 候选优先级清单有没有用"。
- **15 个候选** 是诚实的发现型贡献：9 个在已知位点内提名新效应基因（PLAUR/CNNM2/RBM6 等，FinnGen 自身 p<0.05 支持 RBM6/CNNM2/CD101/RIC8A），6 个弱候选只进补充表；明确"候选 ≠ 因果"，给出名字 + 证据列，供社区优先验证——把最容易被骂的点变成资源。
- 反将一军：领域真正缺的不是第 132 个 MR hit，而是"哪些 MR hit 值得做湿实验"——这张图谱直接回答（操作特性 + 候选分级）。

### 附加隐患：数据质量（build 错位史）
曾有 hg38/hg19 错位致首批结果作废的历史。审稿人若深挖代码会看到修复痕迹。预案：Methods 正面声明"全部 coloc 采用 **rsID 匹配**，与坐标 build 无关"；M9 给单基因收敛（GIT1）与 106/106 pp4 重现作为审计证据；CHANGELOG 的诚实记录本身就是加分项（说明有交叉验证拦截机制）。

---

## 9. 精华总结（供主 agent 转述）

1. **论文形态定稿（2026-08-16）**：三股线——① 图谱资源（31,373 对/123 strong 底表 = 121 FDR-core + 2 灰区）② 方法学操作特性（低 yield 非低召回 + yield 漏斗 0.71%→25.6%）③ 15 个新候选（9 已知位点新效应基因 + 6 弱候选 + GTEx/FinnGen 双复现）。**定投 AJHG**（其他期刊推测已删）。单篇不拆。
2. 核心数字全部核实并统一口径：yield 121/982=12.32%[10.41–14.53]（FDR-core，分结局 BH-FDR q<0.05）与 106/818=12.96%[10.83–15.43]（stage-2 网格，819 中 TCF19×CAD 因 GWAS 区域数据不可用无法评估 coloc）**收敛同报**；召回 121/123=98.4%（FDR-core 全图谱）同报 106/121=87.6%（已知覆盖）；MR 显著集外 strong 仅 2/14,294（AP3S2/ZNF19，均 GWAS 峰显著）；106/106 重现；15 候选 = 9 known-locus 新效应基因 + 6 弱候选；GTEx 6/7 方向一致 + FinnGen 9/9 基因级方向一致（4 个 p<0.05，对齐覆盖 9/15=60% 诚实披露）。
3. 亮点设计：**yield 漏斗**（MR p 阈值单调校准 0.71%→25.6%，F1 图）+ 四层 yield 曲线（12.32%→3.04%→0.014%→0%）就是评审要的"曲线而非单点"，是 Fig 2 主图；F2 图把 15 候选 + GTEx/FinnGen 双复现放在一张图上。
4. 写作前待补分析 4 项**全部收口**：敏感性表 PP.H4≥0.9/0.5（A4 → `s2_pph4_sensitivity_20260816.csv`）、AP3S2/ZNF19 案例细节（A5 → R3 vignette 数字落地）、atlas LD 聚类（A8 → M30，r²≥0.8=106 独立簇，FDR-core 121 子集 ≈104）、全扫描级置换重推（A9 → **数据不可得，作 Limitation 如实披露**，不消耗计算）；其余 A1/A2/A3/A6/A7/A10(部分) 已被 M25–M29 完成。
5. 诚实红线守死：不写"假阳性率 12.32%"、不写"低召回"、不写"双介导/新发现/独立验证"、MR≡SMR 恒等；6 弱候选只进补充表；GTEx/FinnGen 只写"方向一致"；候选一律 "hypothesis-generating"。

---

## 附录 A — 最终故事：我们能诚实讲一个什么样的正结果（原 FINAL_STORY_20260816.md）

# 最终故事：我们能诚实讲一个什么样的正结果（2026-08-16，FDR-core 定稿）

> 目的：回答"现在能将一个怎么样的故事就直接讲了"。本文档把现有全部证据按可投稿口径汇总，
> 明确**能说的**与**不能说的**，列出把故事做实的可选补实验，并直接回答
> "这次改动还有意义吗？之后还要做改动吗？想发更高分"。
> 配套：审计 `docs/INTEGRITY_AUDIT_20260816.md`（本仓库）、写作方案 = 本文档正文 §0–§9、
> 审稿人输入 = 本文档附录 C、图 `results/figures/20260816_F*_v2.png`。

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

---

## 附录 B — 期刊定位定案：定投 AJHG（原 JOURNAL_TARGETS_20260816.md）

# 期刊定位定案（2026-08-16）— 定投 AJHG

> **2026-08-16 定案：定投 AJHG（IF 7.7，Article 类型）**。其他期刊（eBioMedicine/Diabetologia/
> HMG/NC）的推测性定位分析已删除。本文档只保留**事实与竞争情报**，用于支撑 AJHG 投稿定位。
> IF 数字来自 2026-06 Clarivate JCR（2025 数据），经外部网页检索核实。

---

## 一、目标期刊 2025 IF（JCR，2026-06 发布）——事实记录

| 期刊 | 2025 IF | 分区 | 备注 |
|---|---|---|---|
| Nature Communications | ~14 | Q1 | 同类 MR-coloc 图谱文章 2024-10 发过（CD4 T 细胞 T2D/CAD） |
| eBioMedicine | **11.2** | Medicine R&E 19/191 | 2025-03 同类文章（COVID-19 单细胞 MR-coloc） |
| Diabetologia | **10.4** | Endo & Metab 12/198 | EASD 旗舰，T2D 天然契合 |
| **AJHG（定投）** | **7.7** | Genetics & Heredity Q1 | **2025-07 直接对标样稿：ASCVD 单细胞 cis-MR+coloc（Ray et al.，440 关联）** |
| HMG | **3.1** | Genetics Q2 | 原"HMG 5.5–6"为 2021–22 过时值，已废弃 |

> 教训：HMG 教训说明 IF 必须以当前 JCR 为准，不能用记忆值。

## 二、同类文章真实落点（2024–2025，竞争情报）

| 论文 | 期刊/IF | 卖点 | 我们的对应/差距 |
|---|---|---|---|
| 单细胞 TWAS-MR-coloc，COVID-19 药物靶点（2025-03） | eBioMedicine 11.2 | sc-eQTL 26,597、132 因果基因、37 靶点 | 无 sc-eQTL、无功能验证 |
| CD4 T 激活态 sc-eQTL MR+coloc，T2D/CAD（2024-10） | Nature Communications ~14 | 免疫细胞激活态 eQTL、全套 MR+coloc | 我们 bulk 全血，无细胞类型特异 |
| **单细胞 cis-MR+coloc，ASCVD（2025-07）——Ray et al.，直接对标** | **AJHG 7.7** | 440 关联、88% 被 bulk TWAS 漏掉、LIPA→CAD + scRNA-seq/免疫组化验证 | 我们**无功能验证**，但**有全量枚举+操作特性** |

**结论**：MR-coloc 图谱文章在 AJHG 有先例（Ray et al. 2025），方向在 scope 内。我们与 Ray
互补：他答"因果基因作用在哪个细胞类型"，我们答"cis-MR 筛选工具本身有多可靠（操作特性标定）"。

## 三、我们的差异化（AJHG 卖点）

1. **全转录组尺度**：31,373 对全量枚举（同类多为特定细胞/组织子集）。
2. **操作特性曲线（精度漏斗）独有**：没有任何同类文章量化过"MR p 阈值 → coloc 支持率"的
   精度-召回权衡——方法学原创。
3. **跨结局泛化**：T2D/CAD/FBG 三结局同屏，资源复用价值。
4. **独立双复现**：GTEx（eQTL 侧）+ FinnGen R11（结局侧，芬兰独立队列）。

**诚实短板**（审稿人必问，必须有答案）：
- bulk eQTLGen 全血，无细胞类型特异 → 回答：图谱是**全转录组底表**，供单细胞研究复用/下沉；不做细胞特异宣称。
- 无功能验证 → 候选定位为"hypothesis-generating"，不做因果宣称。
- 23 个"catalog 未报道"≠"文献未报道" → 措辞严格。

## 四、定案记录（2026-08-16）

- **定投 AJHG（Article 类型）**，对标 Ray et al. 2025 样稿。
- 投稿准备指南：`docs/AJHG_SUBMISSION_GUIDE_20260816.md`（硬性要求、图表精选 ≤7 个、摘要 ≤200 词、必备投稿文件模板）。
- 不再评估其他期刊；若 AJHG 被拒，另行决策（届时更新本文档）。

## 五、来源（2026-08-16 检索）

- AJHG 2025 IF 7.7：[journalmetrics](https://www.journalmetrics.org/journal/american-journal-of-human-genetics)、
  [Preston LibGuides](https://preston.libguides.com/c.php?g=1531584&p=11475483)
- Ray et al. 2025 AJHG（对标样稿）：[DOI 10.1016/j.ajhg.2025.06.001](https://plu.mx/plum/a/?doi=10.1016/j.ajhg.2025.06.001)
- NC CD4 T T2D/CAD：[PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC11519452/)
- eBioMedicine COVID MR-coloc：[Lancet](https://www.thelancet.com/journals/ebiom/article/PIIS2352-3964(25)00040-4/fulltext)
- HMG 2025 IF 3.1：[journalmetrics](https://www.journalmetrics.org/journal/human-molecular-genetics)

---

## 附录 C — 审稿人视角：影响力提升建议书（原 REVIEWER_IF_ADVICE_20260816.md）

# 审稿人视角：影响力提升建议书（2026-08-16）

> ⚠️ **2026-08-16 定案**：本文档的期刊定位（§1 HMG/PLoS Comp Biol/Genetic Epidemiology/NAR）已被
> `本文档附录 B` **取代**——**定投 AJHG（IF 7.7，2026-08-16 定案）**，其余期刊
> 推测已全部删除。本文档保留作**审稿人视角的改进建议**参考（§2+ 的"精度漏斗/资源化/抗 so-what"思路
> 已并入本文档正文），期刊落点以本文档附录 B 为准。

> 作者：严苛但建设性的期刊审稿人（兼方法学/统计背景）。
> 立场：只建议能**真实增强**论文价值与引用潜力的做法，明确拒绝 overclaim、HARKing、p-hacking、选择性报告。
> 依据：`results/feasibility_20260815.md`、`results/BREAKTHROUGH_PLAN_20260815.md`、`results/archive/PUBLICATION_STORY_20260813.md`、`results/coloc_full_summary_20260815.csv`、`results/coloc_full_{t2d,cad,fbg}_20260815.csv`、`results/archive/feasibility_pilot_20260815.csv`、`results/m22b_window_fix_20260815.csv`、`README.md`。
> 本建议书所有"低成本"分析均只用盘上已有数据，无新的大资源下载。

---

## §0 一句话判定

**HMG 是合理但非最优目标；这份工作的真实价值不在"发现"而在"标定 + 资源"——它最有影响力、最可引用的形态不是"12.96% 这一个数字"，而是"一条精度随 MR 阈值的决策曲线 + 一张带 DOI 的全转录组 coloc 资源表"——把审稿人最想打"so what"的两个点（负结果 + 无新基因）各自翻成卖点。** 现有证据比文档写的更有利于这个叙事，因为全量扫描给出了一根现成的、自己没意识到的**精度漏斗**：扫描层名义显著集 precision=3.0%（129/4248）→ stage-2 clump+IVW 复核存活集 precision=12.96%（106/818；819 中 TCF19×CAD 因 GWAS 区域数据不可用无法评估 coloc）。这把"单工具 cis-MR 精度低"从一句抱怨变成**可操作结论：复核步骤把共定位支持率提升约 4 倍**——这是审稿人和用户都能直接带走的东西。

---

## §1 期刊匹配：HMG 是否最优？

### 1.1 HMG 的评估

| 维度 | HMG（Hum Mol Genet）评估 |
|---|---|
| 主题契合 | **中等偏下**。HMG 近年内容重心在疾病机制/基因型-表型关联/变异功能验证，方法学特征化与纯资源型论文不是其舒适区；"MR 操作特性 + coloc 图谱"这类投稿在 HMG 容易收到"so what / 缺生物学新意"的评审 |
| IF 现实 | 项目文档记 5.5–6 **偏乐观**；HMG 近年 JCR 约 4 上下，且逐年波动。IF 不构成选刊理由，编委会口味与文章类型匹配才是 |
| 长处 | 名字里带 genetics，靶点 MR 社区认得；已有多轮方法学证据（置换标定/HEIDI/Steiger/GTEx）恰好是 HMG 会看重的严谨性 |

**结论：HMG 可投，但不是"匹配度最优"，而且 5.5–6 的 IF 预期需要降档看待。** 若投 HMG，必须把方法学脊柱（操作特性曲线 + 置换标定）做到足够硬，把图谱表做成正式可下载资源，否则会以"负面特征化研究"被低优先级处理。

### 1.2 备选（2–3 个，按推荐顺序）

**备选 1：PLOS Computational Biology（IF ~4.5–5）——方法学/资源双契合，首选替代**
- 理由：PLOS Comput Biol 有明确的 methods/tools/resources 栏目文化，接受"工具评估 + 可复用资源"形态；审稿人是方法学/统计背景，能欣赏"单工具 MR 操作特性曲线"这种**benchmark 式贡献**；负结果若边界清晰、CI 紧、有决策含义，在该刊是可辩护的正贡献。
- 代价：写作需向方法学语言靠拢（把操作特性讲成工具评估）。

**备选 2：Genetic Epidemiology（IF ~2–3）——主题归属地，名誉低一档**
- 理由：这是"遗传筛查工具的操作特性"这一类研究的**编委会天然主场**（接受度最高、审稿人最懂行、拒稿最少）；本项目文档记其 IF 4.9 明显偏高，近年约 2–3。
- 用途：若 HMG/PLOS Comput Biol 被拒，它是最稳的落点，且方法学社区会真正引用。

**备选 3：Nucleic Acids Research 数据库专刊（IF ~16–17）——引用引擎，但需要网页资源**
- 理由：NAR Database Issue 是图谱/资源型论文的引用引擎；coloc 全表若做成**可查询的网页数据库（基因符号检索 + DOI + 版本维护）**，NAR 是最强发刊地。
- 代价：**必须先真建真维护一个可查询站点**（静态 GitHub Pages + 客户端检索可把成本压到低，但需要承诺长期维护）；不做网页就不要投 NAR。

**开放数据型补充：GigaScience（IF ~7–9）**
- 理由：数据论文格式，要求 DOI'd 数据 + 可复现工作流容器化；本项目已有预注册/审计纪律，转 GigaScience 的数据资源形态成本不高，IF 显著高于 HMG 现实值。
- 代价：需把分析做成容器化可复现流程；与"方法学单篇故事"的定位需要重新裁剪。

**通用动作：同一日 bioRxiv 预印本。**
- 无论投哪家，投稿前同日挂 bioRxiv：建立优先权、抢第一波引用（图谱资源靠预印本被引用是常态）、且不损害后续正式发表。

### 1.3 建议矩阵

| 论文最终脊柱 | 首选 | 备选 | 保底 |
|---|---|---|---|
| 方法学（操作特性曲线）为主 | PLOS Comput Biol / HMG | Genetic Epidemiology | Bioinformatics（benchmark 类）|
| 图谱资源为主 | HMG / HGG Advances | GigaScience | EJHG / Human Genetics |
| 网页数据库形态 | **NAR Database Issue** | GigaScience | PLOS Comput Biol |

---

## §2 限制 IF 的根本原因："无聊点"怎么翻成卖点

### 2.1 审稿人眼中最大的三个"so what"

1. **"又一个负结果/特征化数字"**——12.9% 或 3.0% 只是一个比例，没有决策含义。
2. **"没有新基因"**——106 个 strong 命中 81% 落已报道 GWAS 注释，19% 是"上限"不是"新发现"。
3. **"单工具 MR 已经过时"**——审稿人会问"为什么不用多工具 MR / coloc-SuSiE / 精细定位？"

### 2.2 翻盘逻辑（每一条都有实证支撑，非话术）

**负结果 → 边界闭合 + 决策曲线。**
M20/M23 的实证否定的是"MR 召回低"这一此前开放的假设：MR 显著集**外** strong coloc 全量枚举仅 2/27,123（灰区 2/14,294=0.014%、阴性 0/12,829；Poisson 95% CI 上限 ≈0.023%/对），6,000 对分层抽样 0 个。这是**带紧 CI 的决定性边界**——不是"没发现"而是"此假设在此功效范围内被排除"。配合精度曲线（§3a）与精度漏斗（见 §0），把"精度低"升级为**一个可执行的读数：MR 筛查的瓶颈是特异性不是灵敏度；若做单工具 cis-MR，至少要做 LD-clump + IVW 复核，因为它把共定位支持率提升约 4 倍（3.0%→12.96%）**。这就是"so what"的答案。

**无新基因 → 标定 + 复现本身是贡献。**
审稿人真正想确认的是"这工具到底能不能筛对"。106/106 已知命中全量重现 + 81% 文献命中 + 置换标定 FP=1.45% + HEIDI/Steiger/GTEx/FinnGen 四层收敛，合起来是一个**独立审计过的、完全指定的、可复现的优先化管线**——"已知位点被一条干净管线独立回收"是图谱资源型论文的合法性来源，可引用。19% 超目录上限诚实标注为"候选上限"而非"新发现"，配 AP3S2/ZNF19 两个具体案例（§3e）给出**可深挖的例外**。

**单工具过时 → 这正是被测量对象。**
单工具 cis-MR 仍是药物靶点 MR 的**事实标准**（绝大多数发表用单/少工具）。论文测量的就是这套事实标准的操作特性——测量工具之前，没人把它的精度/召回/复核增益量化过。coloc-SuSiE 校准失败、MR≡SMR 恒等已诚实披露，反而让"我们只校准我们测量的事实标准"站得住。

### 2.3 一个必须主动披露的口径问题（否则会被审稿人当成把柄）

盘上有**两个精度数字**：扫描层名义显著集 129/4248=**3.0%** 与 stage-2 复核存活集 106/818=**12.96%**（Wilson CI 10.83–15.43%）。二者不是矛盾，是**同一漏斗的两个阶段**（4248 ⊃ 819 ⊃ 818；129 含已知 106 + 新 23；819 中 TCF19×CAD 因 GWAS 区域数据不可用无法评估 coloc，精度分母 818）。论文必须**并排报告两个精度、写明各自分母**，并把它做成 Fig 2 的漏斗图——否则审稿人会抓住"你们到底是多少"问一年。**这个并排展示本身就是亮点。**

---

## §3 低成本高价值补强（全部只用现有数据）

> 通用数据源：`results/coloc_full_{t2d,cad,fbg}_20260815.csv`（每行含 gene/symbol/outcome/mr_b/mr_p/gwas_min_p/eqtl_F_max/nsnp/pp4/ok），全部分析可直接在这 31,373 行上做。以下每项：做什么 / 为什么提 IF / 成本。

### (a) 精度随 MR 显著度阈值的曲线（操作特性曲线）——**最高性价比，第一优先级**
- **做什么**：对每个 MR p 阈值 t∈{5e-8, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.05, 0.1, 0.2}，算 precision(t)=P(PP.H4≥0.8 | mr_p<t) 并画曲线（同图叠加 gwas_min_p 分层的次曲线）；同一图右侧副轴给累计 strong 命中数（近似"灵敏度"轴）。把 §0 的 3.0%/12.96% 两个点标在曲线上作阶段标注。
- **为什么提 IF**：把点估计升级为**决策曲线**——回答"该把 MR 阈值设在哪、复核阶段值不值"这一领域真问题；操作特性曲线是 benchmark 论文被引用的最常用图（别人引用一个数字不如引用一条曲线）。
- **成本**：几十行 R/Python，<1 小时。**明确不越界**：曲线全阈值如实展示，不得事后挑一个阈值当"新的头条数字"（见 §6）。

### (b) 131 个 strong 的分布图（结局 × 染色体 × biotype）
- **做什么**：按结局分组的条形计数（t2d/cad/fbg）、按染色体的分布（含已知/新增叠色）、biotype 分层（编码/非编码，呼应已有的 76/30 拆解）。可加一张 loci-Manhattan 风格图（x=染色体，y=gwas_min_p 或 eqtl_F，标 strong 点）。
- **为什么提 IF**：图谱型论文的"纹理"——让读者直观看到资源覆盖面与结构；低成本画面感强。
- **成本**：1–2 小时脚本 + 绘图（注意项目配色纪律，勿临时造）。

### (c) 精度与极低率的 bootstrap / Poisson CI——**审稿人会自己算，先算好给他**
- **做什么**：(i) 对 precision 做**按基因的 cluster bootstrap**（同一基因多对存在 LD 相关，不能按对独立抽样）给 3.0% 和 12.96% 各给 bootstrap 95% CI；(ii) MR 显著集外 strong=2/27,123 给 **Poisson 95% CI**（计数 CI≈0.24–7.22 → 率 CI≈8.8e-6–2.7e-4/对，即上界 ≈0.027%），配合 6,000 对抽样 0 个的补充。
- **为什么提 IF**：负结果的强度全在 CI 紧度上——"MR sig 集外 strong coloc ≤0.023%/对（95% CI 上界）"是可引用的紧边界，比"0/6000"更硬、更诚实地表达抽样与全量的一致性。
- **成本**：~1 小时。

### (d) 全表资源交付：可查询 coloc 全表 + schema + 版本
- **做什么**：发布完整 31,373 对（或 49,866 对含失败/未 QC）表格，每行补列：biotype、染色体、GWAS Catalog known/new 标记（复用 m22b 结果）、LD 簇 ID、置换标定命中标记；附**列字典（schema 文档）**、版本号、CHANGELOG；传 **Zenodo 拿 DOI**（免费、期刊普遍要求）。若有余力：GitHub Pages 静态页放客户端可检索的基因符号搜索表。
- **为什么提 IF**：图谱资源论文的**引用抓手就是"可下载 + 有 DOI + schema 清楚"**；没有 DOI 的 CSV 不算资源。任何人做 cis-MR 靶点筛选都该查这张表——"查表引用"是图谱论文最稳定的引用来源。
- **成本**：打包 + schema 文档数小时；Zenodo 免费；静态检索页半天–1 天。

### (e) 示例案例深挖：AP3S2 / ZNF19 讲成小 vignette——**"MR 漏检"的正面证据**
- **做什么**：把仅有的 2 个 MR-sig 集外 strong 写成"单工具 MR 灵敏度失效模式"案例小节（500–800 词 + 每例一张区域图）：AP3S2×t2d（PP.H4=0.935、GWAS 峰 5.5e-11、eQTL F=55.6、但 mr_p=0.38——**lead eQTL 未锚定 GWAS 因果变异**）、ZNF19×cad（2.9e-11）。配图：LocusZoom 式区域图（可复用 m22b 的坐标修正脚本）+ coloc 逐 SNP 后验。
- **为什么提 IF**：这是"无新基因"批评的最诚实反例——不是宣称新发现，而是**用两个具体案例展示单工具 MR 失效的机制（lead-eQTL 锚定失败 vs 假阳性）**，给审稿人一个可记忆的故事点，也让"特异性不是灵敏度问题"的结论有质感。
- **成本**：纯盘上数据 + 2 张区域图，1 天以内。

### (f) 在线/浏览器格式（引用抓手的分级：先 DOI，再可选网页）
- **做什么**：(i) 最低成本：Zenodo DOI + GitHub 仓库（脚本 + 全表 + schema）；(ii) 若走 NAR 路径：GitHub Pages 客户端检索页（gene symbol 搜 PP.H4/mr_p/gwas_min_p，无需服务器）。
- **为什么提 IF**：资源型论文引用量与"可检索性"强相关；即使不投 NAR，一个静态检索页也显著提高被引用率。
- **成本**：(i) 免费数小时；(ii) 半天–1 天 + 维护承诺。**诚实边界**：未真建站点前不得在文中声称"在线资源"。

### 补充（低成本，任务清单之外的顺手项）
- **(g) 精度漏斗图**：4248→819→131（或 129+2）的桑基/漏斗，把 §0 的发现做成 Fig 2 主图——三阶段各自 precision 标注。
- **(h) eQTL F × gwas_min_p 二维分层热图**：PP.H4≥0.8 比率在 (F, gwas_min_p) 网格上的热图——兑现 BREAKTHROUGH_PLAN 承诺的"二维操作特性"，用全量数据一次做对。
- **(i) p12 敏感性扩到全量**：现只有 20/106（p12=1e-6）；在全量 31,373 对扫一遍 p12=1e-6 的 PP.H4（coloc.abf 是解析解，几乎零额外计算），给 strong 计数随 p12 的敏感性带。
- **(j) 若 Richardson 2020 脂质 GWAS（M12 已用）仍在盘上**：把 LDL-C/HDL-C/TG 加入精度曲线分析（MR 扫描分钟级），3 结局 → 6 结局，**近零下载成本直接回应"仅 3 结局"**。注意：若盘上无此数据则不下载，宁可不做。（对照 §6 的 HARKing 边界：加结局必须发生在看结果之前或明确标探索性。）

---

## §4 审稿风险点 + 诚实有效的预应答

| # | 审稿人可能打出的拒稿点 | 预应答（逐字可用） |
|---|---|---|
| 1 | **"仅 3 结局，泛化性差"** | 本文定位是方法学特征化而非疾病全扫描；3 个结局（2 二值 + 1 连续）足以确立操作特性主张，各结局 precision 分开报告（t2d/cad/fbg）；管线结局可移植、schema 支持扩展。**若 (j) 落地：6 结局直接削弱本条**。不声称泛化到全部心血管代谢病。 |
| 2 | **"无新生物学 / 无新基因"** | 本文贡献是**标定与资源**而非发现。106/106 已知命中独立重现 + 81% 文献命中证明管线可信；19% 超目录项诚实标注为"候选上限"；两个 MR 漏检位点（AP3S2/ZNF19）是**从既定、非 HARKed 协议**产生的可深挖例外。"用独立审计管线重现已知位点"本身是资源型论文的合法性。 |
| 3 | **"单工具 MR 已过时，该用多工具/精细定位"** | 单工具 cis-MR 是药物靶点 MR 的**事实标准**；本文测量的正是该标准的操作特性。coloc-SuSiE 校准失败（PCSK9 PIP=0、APOC3 PP.H4=1.0 伪象）与 MR≡SMR 恒等均已披露——"我们只校准我们测量的事实标准，且不声称过时方法最优"。 |
| 4 | **"阴性结果无意义"** | 本阴性是**决定性且闭合假设**：MR-sig 集外 strong coloc 全量枚举 2/31,373、Poisson 上界 ≈0.023%/对、抽样 6000 对 0 个——"MR 召回低"假设在此功效范围内被紧 CI 排除。闭合此前开放假设的负结果 + "瓶颈是特异性非灵敏度"的**正结论**，在方法学刊上是正贡献。 |
| 5 | **"PP.H4≥0.8 调用不可靠 / coloc 假阳"** | 置换标定经验 FP=1.45%（106 strong 富集 ~69 倍）；MR-sig 集外（无信号对照区）strong≈0 佐证调用规则特异；p12 敏感性如实（20/106）。诚实声明置换零假设不建模多效性/样本重叠，HEIDI+Steiger 覆盖因果方向。 |
| 6 | **"到底 12.96% 还是 3.0%？口径混乱"** | 主动并排报告两个 precision（3.0% 扫描层 / 12.96% stage-2 复核存活），各自分母白纸黑字（4248 / 818，819 中 1 对无可评估 GWAS 区域），并做成漏斗图（§3g）——口径差异本身就是"复核提纯约 4 倍"的发现。 |
| 7 | **"eQTL-GWAS 样本重叠"** | 已文档化：重叠使 MR 假阳计入分母 → precision **偏向保守**；置换标定与 UKB-PPP/GTEx 跨源核查不依赖重叠。Limitations 明确声明。 |
| 8 | **"全血 eQTL 组织盲"** | GTEx 6 组织方向一致 69.8%（非全血 65%、同变异 6/6）回答组织稳定性；InsPIRE 胰岛负性核查如实进 Limitations；"全血是此类筛查最常用组织"是测量对象而非缺陷。 |
| 9 | **"MR p<0.05 阈值主观"** | 阈值曲线（§3a）展示全阈值行为；stage-2 用 FDR；PP.H4 分布全量报告，不做单一阈值宣称。 |
| 10 | **"图谱只是重复 GWAS 已知位点"** | 81% 命中 = 管线校准成立的证据；19% 上限诚实标注；"可复现性即贡献"。 |
| 11 | **"多个结局/千基因无多重校正就下结论"** | 图谱为描述性/假设生成性质；strong 调用阈值 + 置换标定 + 不把 atlas 当家族级显著性宣称。 |

---

## §5 引用抓手：资源论文靠"可复用"拉引用，具体五件事

1. **DOI + 版本化的全量 coloc 表**：Zenodo 存档（免费），表格含 schema 文档、版本号、CHANGELOG；版本化让"Dual-Channel MR Atlas v1.0"成为可引用句柄。这是第 1 引用来源。
2. **方法引用链**：coloc.abf→Giambartolomei 2014、eQTLGen→Võsa 2021、OpenGWAS、SMR/HEIDI→Zhu 2016、Wilson/Poisson 公式——凡引用本资源者必连带引用方法源；论文 Methods 里把这些写齐，资源论文的引用网络就从这里长出来。
3. **一条可被"数值引用"的头条**：把操作特性做成**单一句子可引用的数字+曲线**（"MR sig 内 coloc 支持 3.0%→复核后 12.96%，sig 外 ≤0.027%/对"）——benchmark 论文最常被"引用一个数字"。§3a 的曲线图就是被引用的那个图。
4. **"怎么用"小节**：写清楚"用户如何把这张表 join 进自己的 MR 屏幕、哪些列回答'我的 MR hit 有没有 coloc 支持'"——把读者变成用户，用户才会引。
5. **可选但最强的钩子**：静态检索页（GitHub Pages 客户端搜索基因符号）+ GitHub 仓库放全部脚本。若真建真维护，这才是 NAR/GigaScience 级资源的引用引擎（§1.2 备选 3）。

---

## §6 诚实边界：这些提升手段越界，不能做

1. **不得把 12.9% 改称为"假阳性率"**——coloc 只在 MR 存活集上跑，无法区分 MR 假阳与 coloc 假阴；只能叫"共定位一致性/支持率"。
2. **不得把 19% 超目录项写成"新发现"**——那是"未进 GWAS Catalog 注释"而非"经独立验证的新因果位点"；只能写"新候选上限"。AP3S2/ZNF19 是案例展示，不是验证过的靶点。
3. **不得在看过结果后加结局/改阈值来抬数字**（HARKing）——加结局（如 §3j 脂质）必须在分析前或明确标"探索性"并预注册追加；操作特性曲线全阈值如实展示，不得事后挑阈值当新头条。
4. **不得把置换零假设说成覆盖多效性/样本重叠**——它只覆盖"独立信号仅靠 LD 巧合"。
5. **不得把 MR≡SMR 说成独立互证**；不得把 eQTL coloc 与 pQTL coloc 说成独立介导证据。
6. **不得把蛋白层关停（10/76）表述为任何生物学结论**——只是"UKB-PPP 面板未覆盖"，不是"蛋白层无作用"。
7. **未真建真维护站点前，文中不得出现"在线资源/浏览器"字样**；期刊会核查。
8. **不得只报 129 或只报 106**——两个精度、23 个"新增"、106 已知命中必须全量并报。
9. **不得移动预注册门柱**；全网格（含失败行）如实报告。
10. **不得为凑"资源"下载新的重量级数据**——本建议书所有补强均盘上可做；下载需求若出现，先走资源仲裁规则。

---

## 精华总结（5 句话，供转述）

1. HMG 5.5–6 偏乐观且非最优匹配——方法学/资源型论文首选 PLOS Computational Biology，备选 Genetic Epidemiology，若能真建网页数据库 NAR Database Issue 是引用引擎；无论投哪家同日挂 bioRxiv。
2. 最大未变现资产是盘上的**精度漏斗**：扫描层 precision=3.0%（129/4248）→ stage-2 复核 12.96%（106/818）——把它做成决策曲线主图，"单工具 cis-MR 低精度"从抱怨变成可操作结论"复核步骤提升共定位支持率约 4 倍"。
3. 最高性价比补强全在现有数据上：MR 阈值操作特性曲线、强共定位的结局×染色体分布、按基因 cluster bootstrap + Poisson 紧 CI（sig 外 ≤0.023%/对）、Zenodo DOI 的全量 coloc 表 + schema、AP3S2/ZNF19 两个"MR 漏检机制"vignette。
4. 审稿人四大拒稿点（仅 3 结局/无新基因/单工具过时/阴性无意义）都有诚实有效的预应答，核心是"标定 + 资源 + 决定性负结果边界"三位一体，绝不把负结果吹成正发现。
5. 红线不变：12.9% 只是共定位一致性率不是假阳率、19% 只是候选上限、MR≡SMR 不互证、加结局不得 HARK、无网页不得自称在线资源。
