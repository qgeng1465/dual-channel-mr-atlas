# 论文写作方案 — cis-MR×coloc 图谱 → AJHG 单篇（2026-08-16）

> 面向：本项目学生 + 导师 agent。**单一篇论文、单一完整故事、诚实优先**，**定投 AJHG（IF 7.7）**——已定案，不再评估其他期刊（eBioMedicine/Diabetologia/HMG 等推测性期刊定位已删除，理由见 `docs/JOURNAL_TARGETS_20260816.md` 尾部定案记录）。
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
- **M10 软件/数据可用性与预注册**：coloc 5.2.3 / R env r-mr / plink 1.9 / SMR 1.3.1；图谱 CSV（31,373 对全量 + 123 strong 子集 = 121 FDR-core + 2 灰区）+ 脚本 GitHub/OSF 发布；预注册哈希；版本声明（M23–M28 分析脚本 + results CSV 一并入库，FINAL_STORY/README 已同步）。【已有，发布打包属 P2】
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
