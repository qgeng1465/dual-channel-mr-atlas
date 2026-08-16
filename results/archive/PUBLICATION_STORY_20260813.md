# 发表故事与投稿计划（2026-08-13，多 agent 终审固化）

> 来源：`mr-atlas-verify-and-story` 工作流（3 verify + 3 story + 1 synthesize，7 agent 全部完成，
> 0 错误，614k tokens）。所有数字经 agent 逐行复算，且与本日独立数据核查 agent 交叉验证。
> 本文档是用户「怎么串起现有研究、讲什么故事、如何发表」的可执行答案。

---

## 0. 2026-08-13 更新（M12 脂质校准 + M13 UKB-PPP 双平台核查落地）

> 本节是 2026-08-13 第一轮新增，不改写上文结论，只把蛋白通道的证据等级提升记录在案。

**蛋白通道从"未校准"升级为"已校准 + 双平台核查"（但红线 3 的"无新发现"仍成立）：**

- **M12 脂质校准（杠杆 B）**：deCODE 蛋白工具在 Richardson 2020 UKB 脂质结局下**回收教科书因果**
  —— PCSK9↑→LDL-C/apoB↑（p≈1e-15，MRE/FE/WM/Egger 全稳健）、APOC3↑→HDL-C↓（p=2.6e-8 全稳健）；
  PCSK9×TG/HDL 与 INSR 全 null（负对照正确）。→ **蛋白工具本身有效，心血管结局未达标是结局侧信号问题**。
- **M13 UKB-PPP 双平台核查**：PCSK9×CAD 在**两个独立 pQTL 平台**（deCODE 冰岛 n=35,559 / UKB-PPP Olink n≈34k）
  同一 cis 变异 rs11591147（R46L）同向强显著（deCODE b=+0.190 p=5.5e-12 / UKB-PPP b=+0.209 p=2.0e-22）；
  cis 共享变异方向一致率：PCSK9 r=0.338 (p<1e-232)、ANGPTL3 r=0.396 (p<1e-287)（强 cis 信号平台一致；
  弱信号 GLP1R/GCG r≈0 如实报告）。
- **GTEx 全血交叉**：PCSK9→CAD 经 GTEx 全血 eQTL 亦同向（b=+0.13 p=2.7e-5），转录通道内部跨数据源一致。
- **HEIDI 全集（P3 完成）**：818/819 测试全量 SMR+HEIDI → **全 MR 显著集通过 48.6%**（多数不一致立住）、
  strong coloc 子集 71.7%（76/106）、PP.H4 单调 63.6%→79.3% → Fig 5 桑基图数据齐。
- **Steiger 方向（P5 部分）**：76 优先基因中 **73/76（96.1%）eQTL→outcome 正向**、67 达 Steiger p<0.05；
  仅 3 个方向未决（不显著反向，如实标注）→ 三层收敛（coloc+HEIDI+方向）立住。

**叙事影响（诚实）**：
- ✅ 蛋白通道可以写为"**经脂质校准 + 双平台核查的次要层**"，工具可用性有独立证据（回应审稿人"蛋白工具是坏的"）。
- ⚠️ 红线不变：仍**无新因果发现**（PCSK9 是已知降脂靶点）；APOC3 仍无 UKB-PPP 复现（面板缺如）→
  心血管命中仍集中在已知靶点，双通道图谱"双显著 0 对"的事实未变。
- 期刊评估：HMG 5-6.5 路径的 P2（UKB-PPP）已实质完成、**P1（GTEx eQTL 复现）2026-08-13 晚已完整
  落地**（69.8%/非全血 65%/同变异 6/6），比正文 §3 的勾选状态领先。

### 同日第三轮（E1 置换标定 + P4 药物注释落地）

> 2026-08-13 第二轮新增。见 `results/coloc_permutation_20260813.md` + `results/drug_annotation_20260813.md`。

- **E1 负对照/permutation 标定完成（eBioMedicine 闸门第 1 项）**：106 个 strong 对 × B=100 置换
  （打乱 eQTL beta 跨 SNP，保留 GWAS 信号与 LD/MAF 结构）→ **零假设经验 FP 率（PP.H4≥0.8）= 1.45%**，
  期望 FP 仅 1.54/106；观测 106 个 strong 调用富集 ~69 倍。有效性校验：重建 coloc 与 `transcript_coloc_hits`
  逐对 PP.H4 **偏差全为 0.0000**（106/106）→ 标定针对的就是产生命中集的同一流程。
  意义：**coloc 的 PP.H4≥0.8 阈值在"独立信号仅靠 LD 巧合"零假设下特异性极高**，
  "13% 一致性率"可辩护为"MR 显著集约 1/7 获共定位支持、该子集假阳率 ~1.5%"。
  **唯一脆弱对点名**：ZNF34×t2d（置换 FP 22%、区域仅 239 SNP）→ 降级为"脆弱/weak"标注。
  边界（不过度声称）：置换零假设不建模多效性/样本重叠/共享上游混杂，仍需 HEIDI+Steiger 覆盖。
- **P4 药物注释层完成（76 优先基因，Open Targets + ChEMBL 双源）**：**56 编码 / 20 非编码伪影**；
  **仅 2 基因有已知药物**（KCNJ11 15 条药物条目/14 个独特化合物——8 种磺脲类一线药物 +
  repaglinide/nateglinide/mitiglinide + K⁺ 通道开放剂 minoxidil/diazoxide/pinacidil，其中 mitiglinide
  以游离碱与钙盐双条目计入；BMPR1A 2 个已批 BMP 生物药），
  **74 基因无任何已知药物**（ChEMBL 机制交叉核验 54 编码基因 0 临床机制）；疾病关联显著者仅
  KCNJ11→T2D 0.864、LIPA→MI/CAD 0.534/0.453。Pharos 从本机不可达（403）→ TDL 列为 `(est)` 估算。
  → **76 优先基因 97% 是"未成药靶点"**，这是图谱型论文的新贡献点：优先化表里绝大多数是
  尚无药物的可干预候选（配合 P5 Steiger 方向，形成"转录介导、方向确定、未成药"的可验证靶点清单）。
- 期刊影响：E1 让"一致性率"读数具备置换标定背书；P4 让 Fig 3/Table 的优先基因表可加"可成药性"列。
  但 eBioMedicine 闸门 E2（≥1 真新发现）**仍不满足**——无新发现就是无新发现，药表不是发现。
- **独立核查背书（2026-08-13 晚）**：第三方 agent 对 E1 FP 率（154/10600=1.4528%）、FinnGen 对齐
  （108/108 可对齐、92.59%、双显著 46/46=100%）、独立信号（100 唯一 top_snp）、P4 药物表
  （76 行/56+20/仅 KCNJ11+BMPR1A）**四项逐列独立重算全部通过**，记录见
  `results/verification_20260813.md`；两处非阻断细节已修（KCNJ11 药物条目口径 15=14 独特化合物、M7 注释行号）。

### 同日第四轮（口径统一 + LD 聚类 + 编码一致率 + 通路富集阴性）

> 2026-08-13 第三轮新增。分别见 `concordance_denominator_20260813.md`、`ld_clustering_20260813.md`、
> `coding_concordance_20260813.md`、`pathway_enrichment_20260813.md`。

- **一致率分母口径统一（诚实修正）**：规范 MR 显著测试集 = `transcript_coloc.csv` **819 行**（759 唯一基因，
  全 stage2_pval<0.05）；"818"实为 HEIDI 可运行数（819−1 失败测试），"820"含 1 个未入规范集的额外行。
  选定规范口径 **106/819 = 12.9%（Wilson 95% CI 10.82–15.42%）**，四份文档（story/overlap/coloc/evidence_brief）
  已统一并加注；三种口径一致率差 <0.1pp，结论不变（`concordance_denominator_20260813.md`）。
- **r²-LD 聚类（M5 完成，1000G EUR 参考面板）**：100 个唯一 top-SNP → 98 匹配参考面板，
  **r²≥0.8 下 92 个独立 LD 簇**（6 个多 SNP 对，含两个 Metazoa_SRP ENSG 共享 LD）；
  r²≥0.6 敏感性 91 簇（+rs2239540/rs3824535）。**106 命中 → 92 LD 独立位点**，M5 独立位点表升级为 LD 级。
- **⚠️ 参考面板质量核查（第二轮独立核查触发，已闭环）**：本地 1kg.v3 面板**零杂合、~40% 个体缺失**
  （疑为 haploid 化/简化，非标准二倍体）；但亲测对照验证 r² 仍反映真实 LD 结构（K-ATP block
  rs5215×rs757110 = 0.918 vs 真实 0.94–0.99；远距离背景 median 0.002 正常），**数值轻微低估**。
  92 簇结论保留但为**保守上限**（真实可能 ~90–91）；stage2 本地 plink clump（r²<0.01@1000kb）
  受影响很小（阈值对噪声不敏感，实测吸收率 median 0.011 正常）。位置聚类（bim 坐标，不依赖
  基因型内容）独立佐证 **89–88 独立位点**。主口径仍为「106→100 唯一 top_snp」文本去重
  （`ld_clustering_20260813.md` 终版）。
- **仅编码版一致性率（M5 完成，Ensembl biotype 注释 759 基因）**：**编码 76/637 = 11.9% ≈ 全量 12.9%**；
  非编码/未注释 30/182 = **16.5% 反而更高** → **头部读数不被非编码/伪影基因稀释**，且非编码层并非"灌水"。
  主结论"MR 显著集多数未获共定位支持"对纯编码基因同样成立。
- **通路富集（探索性，阴性）**：56 编码优先基因在 GO:BP/MF/KEGG/Reactome 上 **g_SCS 校正后 0 显著**（最佳 p=0.35）。
  优先基因功能异质（管家基因为主），无共享通路信号——阴性本身是诚实的生物学信息；选择偏倚边界已写明。
- 期刊影响：编码版一致率 + LD 独立位点数（92）使 Fig 2 的"一致性率×biotype 分层"与 Fig 3 的独立位点表
  具备精确口径；通路富集阴性进 Discussion/Limitations（不构成发现）。

### 同日第五轮（InsPIRE 胰岛通道 MR 完成，阴性）

> 见 `inspire_islet_20260813.md`。回应"KCNJ11 非转录介导"红线的组织特异核查。

- **M16 InsPIRE 胰岛 eQTL MR（完成，探索性；第二轮核查后修正）**：15/76 优先基因有胰岛 eQTL；
  首跑 14 对、重跑 13 对结局匹配。**跨组织方向一致：首跑 1/14（LRIG1×t2d）但 LRIG1 重跑
  OpenGWAS API 空返回、不可复现 → 当前可复现口径 0/13**；11 个显著胰岛 MR 中 **10 个与全血
  通道方向相反且两次运行全一致**（C18orf8/EIF2B2/MLH3/HMBS/HSD17B12/ZNF34/LPIN3/PGAP3/SNX16/DCAF16）。
  **KCNJ11 胰岛 lead rs2283253 在 T2D GWAS null（p=0.87，两次一致）** → **胰岛通道未救活 KCNJ11**，
  "KCNJ11 的 T2D 信号非转录介导（全血或胰岛皆否）"红线维持。另修 M16 CSV 列名 `gwas_b`→`mr_b`
  （原列名实为 MR Wald 效应，第二轮核查发现，已重跑更正）。
- 诚实限定：两种组织用不同 lead（独立工具）非同信号复现、n 小、探索性——**不得写"1/14"**，
  唯一稳定读数 = **10/11 反向 + KCNJ11 胰岛 null**（均两次运行一致），作负性核查进 Limitations。

### 同日第六轮（P1 GTEx 复现完成 + 第二轮独立核查闭环）

> 见 `gtex_replication_20260813.md` + `verification_20260813.md`（第二轮）。

- **P1 GTEx 跨组织复现（M15 完成，正向）**：106 strong 命中，93/106 在 GTEx（6 组织）有显著
  cis-eQTL；63 对结局匹配，**方向一致 44/63 = 69.8%**；**排除 Whole_Blood 后非全血组织
  26/40 = 65.0%**（脂肪 76%、胰岛 83%）；同变异子集 6/6 = 100%。→ **eQTLGen 全血通道的优先化
  信号跨独立组织方向稳定，非全血组织伪影**。注意：这是**方向一致性**指标，与 coloc 一致率
  （12.9%）是不同指标，不改变主结论。
- **第二轮独立核查（agent 触发）→ 三处诚实修正闭环**（`verification_20260813.md`）：
  ① M16 CSV 列名 `gwas_b` 实为 MR Wald 效应 → 改名 `mr_b` 重跑；② LD 面板核查：本地 1kg.v3
  零杂合/40% 个体缺失，但亲测对照 r² 仍反映真实 LD（rs5215×rs757110=0.918 vs 真实 0.94–0.99）→
  92 簇保留为保守上限，位置聚类 89–88 独立佐证，stage2 clump 受影响很小；③ M16 LRIG1 唯一一致对
  重跑 OpenGWAS API 空返回**不可复现** → 一致率改 0/13 可复现口径（10 反向 + KCNJ11 null 两次全一致，稳健）。

---

## 1. 最终阳性判断（必须诚实，防学术不端）

**现在可发表，但请勿称"新因果发现型阳性"**。真实构成是：

- ✅ **可辩护的阳性 = 转录通道多层收敛图谱型**：106 个 strong coloc 位点（76 编码 / 30 非编码伪影倾向；
  PP.H4≥0.9 共 62）+ HEIDI 通过 89/126（70.6%）+ 组织三角 55/79（70%）+ FinnGen 双显著子集 100% 同向
  （T2D 27/27、CAD 19/19）。四维独立证据收敛。
- ⚠️ **三处措辞红线（审稿人一跑除法/复算即否，投稿前必须改）**：
  1. **MR-Wald ≡ SMR 是数学恒等**（同一 top cis-SNP 的同一比率，b 同到 6 位小数）——不是两个独立方法。
     "三方法独立互证"必须撤为"**单点 Wald 自洽 + 两个区域级独立检验（coloc、HEIDI）**"。
  2. **"MR 假阳性率 13%"是错误标签**——106/819=12.9%（Wilson 95% CI 10.8–15.4%）是"MR 显著且获区域共定位支持的**一致性率
     （concurrence rate）**"；coloc 只在 MR 存活者上跑，无法区分 MR 假阳性 vs coloc 假阴性。
  3. **蛋白通道不构成"双通道图谱"**：工具级重叠仅 3 对（全 INSR）、双显著 0、INSR 跨通道方向相反；
     APOC3×T2D/FBG 是单工具边缘 p（0.011/0.007）无多重校正 + 两 coloc 方法判不共定位 + FinnGen 不显著
     + **UKB-PPP 无 APOC3（复现断头）** → 只写"单工具敏感性支持，待独立平台复现"，**不得写"独立救回/新发现"**。
- ❌ **没有新因果发现**：所有强命中均为已知位点（KCNJ11/PCSK9/APOC3 都是既有药物靶点）或伪探针。
  想上 eBioMedicine 必须先补 ≥1 个带正交验证的新发现。

**一句话**：一篇"转录通道阳性的多层收敛图谱 + coloc 一致性率方法学读数"的论文，现在能投 IF 3-4；
补做清单后 HMG 5-6.5；eBioMedicine ≥9 是路线图不是现状。

---

## 2. 主故事（选定）

**A（方法学脊柱）+ C（图谱框架）融合，B 的校准案例作验证图。**

- C 提供生物学主体：106 基因转录图谱 + 诚实蛋白次要层
- A 提供唯一的新贡献骨架：coloc 一致性率 + MR≡SMR 校准声明
- B 降为校准展品："管线能正确回收已知磺脲类/降脂靶点（KCNJ11/PCSK9/APOC3）" → 方法可信，新候选才可信

### 论文标题（首选 + 2 备选）

1. **首选**："Shared causal loci are the exception in cis-Mendelian randomization: a transcriptome-wide colocalization atlas for type 2 diabetes and coronary artery disease with cross-tissue and cross-cohort convergence"
2. 备选（方法学脊柱）："Only one in seven cis-MR signals survives regional colocalization: a concordance atlas for cardiometabolic disease and the case against counting single-instrument MR and SMR as independent confirmations"
3. 备选（校准/靶点）："Triangulated cis-MR colocalization for drug-target prioritization in T2D and CAD: calibration against established targets (KCNJ11, PCSK9, APOC3) and prioritized transcript-mediated genes"

### 一句话摘要

> In a cis-eQTL Mendelian randomization screen of ~1,500 genes across T2D/CAD/FBG, 106 loci reached
> strong colocalization (PP.H4≥0.8; 62 ≥0.9; 76 protein-coding); 70.6% pass HEIDI homogeneity, 70% of
> tested hits are reinforced by ≥1 tissue eQTL, and both-significant FinnGen subsets are 100% directionally
> concordant (27/27, 19/19). Only 12.9% of cis-MR signals survive regional colocalization, and single-instrument
> MR and SMR are shown to be the same statistic — a concordance framework for reading cis-MR drug-target evidence.

### 核心图 5 张

- **Fig 1（设计 + 诚实双通道漏斗，最重要）**：转录漏斗 1,513 基因×3 结局 → MR 存活 818 测试 → **106 strong**
  （58 T2D/46 CAD/2 FBG），下嵌 biotype 拆解（76 编码 / 30 非编码）；蛋白漏斗 11→5→4→2（标注"均为已知校准位点"）。
  两漏斗**视觉不等高**，图注 "parallel screening; integration limited by instrument overlap (n=3)"。
- **Fig 2（一致性率方法学头条）**：(a) 12.9% [Wilson CI 10.82–15.42%]（106/819，规范口径 = `transcript_coloc.csv` MR 显著测试数；
  HEIDI 上下文的"818"=819−1 个非可运行测试，见 `concordance_denominator_20260813.md`）按结局/通道/biotype 分层；
  (b) 先验敏感性瀑布：主设定 → p12=1e-6 后 20/106（18.9%）；(c) 818 测试 PP.H4 分布；(d) 一致性率 × F 分箱。
- **Fig 3（收敛展品：校准 + 新优先位点）**：KCNJ11/ABCC8×T2D（Wald 0.595 自洽 + coloc + HEIDI 0.128）、
  PCSK9×CAD（R46L PP.H4=1.000，HEIDI 异质性如实标注）、APOC3×CAD（0.997）、LIPA×CAD（0.999/HEIDI 0.457）、
  RASD1×CAD（0.986）。图注"单点自洽 + 两个区域级独立检验"。
- **Fig 4（外部复现层）**：FinnGen β 散点（双显著子集 27/27、19/19 高亮，弱信号行降权，修 rs17716350 后重画）
  + 组织三角热图（55/79 × 6 GTEx 组织）。
- **Fig 5（二选一）**：优先 HEIDI 全集（818 测试）→ 不一致分类桑基（多信号/低功效/伪探针/多效性）；
  否则画 ~60-70 编码优先基因证据得分表。

### 正文结构

Intro（立论：cis-MR 命中大多不能直接当靶点；MR 显著但区域共定位不过的比例未被系统量化 + 单工具 MR 与 SMR 常被双重计数）
→ Results（R1 漏斗与图谱；R2 一致性率与敏感性；R3 收敛展品与校准；R4 跨组织/跨队列复现；R5 蛋白通道如实次要层；R6 不一致分类）
→ Methods（预注册 v6 修订 + 结局 ID 偏差记录；数据/工具/MR/coloc.abf + 敏感性/HEIDI/分类协议/正负对照/诚实协议）
→ Discussion（13% 含义、MR≡SMR、何时单点估计够用）→ **Limitations + 预注册合规（独立小节，审稿人必查）**。

---

## 3. 投稿前差什么（must-do / 加分项）

### must-do（任何档位绕不开，纯改稿零资源争议）
- [ ] **M1 措辞纪律 7 处**（§1 三条红线 + 具体替换句，见 IMPROVEMENT_STRATEGY §2.1 措辞纪律项）
- [x] **M2 修 FinnGen rs17716350 对齐 bug**（**已修，2026-08-13**：根因 = M7 第 54 行 `%in%` 整列向量判断
      而非逐行 `==`，rs17716350 的 T>A 次等位编码被错判 aligned；改逐行对齐后重算——
      移除 T>A 行、108/108 可对齐、整体一致率 91.7%→**92.6%**；头条双显著子集
      **46/46=100% 一致**（T2D+CAD，与原 27/27、19/19 一致））
- [ ] **M3 HEIDI 计数更正 87/128 → 89/126（70.6%）**（分母 126 非 NA）——**2026-08-13 可进一步升级**：
      全量 strong coloc 集 76/106=71.7%（`heidi_full_20260813.md`），正文引用全量数
- [ ] **M4 p12=1e-6 敏感性转正为头条数字**（20/106=18.9%）
- [x] **M5 LD 聚类 → ~100 独立位点表 + biotype 注释 + 仅编码版一致性率**（**2026-08-13 已完成**：
      ① top-SNP 去重：106 strong → **100 唯一 top-SNP 信号**（5 簇多基因共享，全 T2D：rs11257655[CAMK1D/NUDT5]、
      rs223490[KRT8P46/LRRC37A15P/RP11-10L12.2 三假基因]、rs2246012[ARG1/MED23]、
      rs2280018[NTAN1/RRN3]、rs9991574[DKFZP434I0714/FBXW7]）；② **r²-LD 聚类（M18，1000G EUR）：
      r²≥0.8 → 92 个独立 LD 簇**（rs10888385/rs11204675 等 6 个多 SNP 对，含两 Metazoa_SRP ENSG）；
      **⚠️ 2026-08-13 面板核查：本地 1kg.v3 零杂合/40% 个体缺失（非标准二倍体），但亲测对照
      r² 仍反映真实 LD（rs5215×rs757110=0.918 vs 真实 0.94–0.99；远距离背景 0.002 正常），
      92 簇为保守上限（真实可能 ~90–91）；位置聚类（bim 坐标）独立佐证 89–88 独立位点；
      主口径仍为文本去重 100**（`ld_clustering_20260813.md` 终版）；
      ③ biotype 注释（Ensembl 759 基因）+ **仅编码版一致率 76/637=11.9% ≈ 全量 12.9%**
      （`coding_concordance_20260813.md`）；④ 分母口径统一 **106/819=12.9%[10.8–15.4%]**
      （`concordance_denominator_20260813.md`））
- [x] **M6 预注册 v6/v7 修订**（**已完成，先前轮次**：`docs/PREREGISTRATION.md` §6.1–6.9——
      结局 ID 变更记录、FBG 引入、T2D N 口径（655,666 vs 62,892）、MD5/sha256 澄清、
      H1/H2 定位、组织三角、coloc-SuSiE 校准失败记录；哈希锁已重算一致）。关键诚实项：
      **APOC3×T2D coloc PP.H4=1.0 被 coloc-SuSiE 识破为假阳性（GWAS 侧 p=0.054 不显著），
      PCSK9×CAD 校准失败 → coloc-SuSiE 不作证据，主闸门维持 coloc.abf**
- [~] **M7 一致性率加 CI + 分层**（**CI 已加**：106/819 = 12.9% Wilson 95% CI 10.82–15.42%；
      **biotype 分层已做**：编码 11.9% / 非编码 16.5%；待补：按结局/通道/F 分层）
- [ ] **M8 蛋白通道所有主张降为探索性**；APOC3×T2D/FBG 标"待 INTERVAL/ARIC/SCALLOP 或文献复现，否则撤回"

### 加分项（IF 5-7 前提）
- [x] **P1 eQTL 侧独立复现**（GTEx 独立 eQTL 源，106 命中方向/效应复现）——**M6 v2 已产出 6 组织**
      eQTL MR（全血 PCSK9→CAD +0.13 p=2.7e-5 与 eQTLGen 同向）；**2026-08-13 全量核对完成
      （M15）**：106 命中 93/106 有 GTEx 显著 cis-eQTL，63 对结局匹配，**方向一致 44/63 = 69.8%**，
      非全血组织 26/40 = 65.0%，同变异 6/6 = 100%（`gtex_replication_20260813.md`）
- [x] **P2 UKB-PPP 蛋白跨平台**（**M13 已完成**：8/11 面板，PCSK9×CAD 双平台强复现，
      cis 共享变异一致率 PCSK9/ANGPTL3 显著；APOC3 缺如显式写缺口）
- [x] **P3 HEIDI 跑满全 818 测试**（**M10b 已完成**：全量 SMR+HEIDI 落地，见
      `results/heidi_full_20260813.md`）——**全 MR 显著集 HEIDI 通过 398/819=48.6%**（多数不一致叙事立住），
      strong coloc 子集 76/106=71.7%，PP.H4 单调 63.6%→79.3%；Fig 5 桑基图数据齐
- [x] **P4 药物性注释层**（**已完成**，2026-08-13：Open Targets + ChEMBL 双源，56 编码/20 非编码，
      仅 KCNJ11/BMPR1A 有药、74 无药；Pharos 不可达故 TDL 为 `(est)`；见
      `drug_annotation_20260813.md`）——优先基因表加"可成药性"列
- [~] **P5 非 KCNJ11 基因正交验证**——**Steiger 已完成**（76 优先基因 96.1% 正向，`steiger_direction_20260813.md`）；
      **InsPIRE 胰岛 eQTL MR 已完成（M16，2026-08-13，阴性；第二轮核查修正）**：15/76 有胰岛 eQTL、
      首跑 14 对、重跑 13 对结局匹配，**跨组织方向一致首跑 1/14（LRIG1，重跑 API 空返回不可复现）
      → 可复现口径 0/13**；11 显著胰岛 MR 中 **10 个与全血方向相反（两次运行全一致）**；
      **KCNJ11 胰岛 lead rs2283253 在 T2D null（p=0.87）→ 胰岛通道未救活 KCNJ11，红线维持**；
      探索性、不同 lead 工具、n 小，见 `inspire_islet_20260813.md`（§2b 可复现性核查）
- [x] **P6 eQTLGen-GWAS 样本重叠敏感性**（**已写，2026-08-13**：`overlap_sensitivity_20260813.md`——
      定性+方向性：重叠偏倚让 MR 假阳垫进分母 → "13% 一致率"**偏向保守**；E1 置换标定与
      M13 跨平台核查是不依赖重叠的独立证据；定量校正需 cohort 重叠矩阵，Limitations 明确声明）
- [ ] **P7 coloc 对 1000G vs UKB 双 LD 面板敏感性**

### eBioMedicine ≥9 硬闸门（缺一不可）
- [x] **E1 负对照/permutation 标定**（**已完成**，2026-08-13：零假设 PP.H4≥0.8 FP 率 1.45%、
      期望 FP 1.54/106、重建偏差 0.0；见 `coloc_permutation_20260813.md`）——"13% 一致性率"
      升级为可辩护假阳性读数（"该共定位子集假阳率 ~1.5%"）
- [ ] **E2 ≥1 个真新因果发现 + 正交/功能支持**（**现状没有就是没有**；P4 药表 97% 未成药
      靶点是图谱贡献，不是新发现）
- [ ] E3 UKB-PPP 落地；- [ ] E4 方法学正式估计当主图

---

## 4. 期刊排序

| 档位 | 期刊（IF） | 达到条件 |
|---|---|---|
| 现状可投（措辞+bug 修完） | Human Genetics (~3.9) ＞ EJHG (~3.7) ＞ Front Genet (~3.5) | M1–M8 |
| 中等补做后（首选现实目标） | **Hum Mol Genet (~5.5)** ＞ Genetic Epidemiology (~4.9) ＞ Cardiovasc Diabetology (~5.5) ＞ Atherosclerosis (~5) | must-do + P1–P7 |
| 冲高（路线图非现状） | Diabetologia (~8-9) → eBioMedicine (~9-11) | E1–E4 全满足 |

---

## 5. 诚实清单（每条都是审稿人可能抓到、且已实锤）

1. **MR≡SMR 恒等**：全部 143 行精确满足 b_SMR=b_GWAS/b_eQTL；一致是自洽性检验非独立互证；真正独立的只有 coloc + HEIDI。
2. **125/128 SMR 显著是选择偏倚构造性必然**（测试集全来自 MR-显著预选位点）——零信息量，不进正文；唯一读数 HEIDI 89/126。
3. **"13%"是 coloc 一致性率不是假阳性率**；coloc 只在 MR 存活者上跑，无法区分 MR-FP vs coloc-FN。
4. **coloc 先验敏感性**：p12=1e-6 后仅 20/106（18.9%）保持 PP.H4≥0.8；44/106 落 0.80-0.90 边缘；strong 依赖主设定 p12=1e-5。
5. **30/106 命中是非编码/假基因/伪影倾向位点**（RP11-*/CTD-*/Metazoa_SRP/U6/hsa-mir-296/SERBP1P3/KRT8P46）。
6. **蛋白通道 HEIDI 判定未校准**：deCODE（冰岛）配 1000G EUR LD、无 effectAlleleFreq（ImpMAF+--disable-freq-ck）、
   阳性对照 PCSK9×CAD 过不了 HEIDI（p=1.28e-4）→ 蛋白侧 HEIDI 只作探索。**（2026-08-13 补充）**此"未校准"
   仅指蛋白侧 LD/频率近似的 HEIDI 判定；工具级可用性已由 M12 脂质校准（PCSK9→LDL/apoB、APOC3→HDL 教科书回收）
   + M13 UKB-PPP 双平台核查（PCSK9×CAD 同变异同向双 GWS）独立补足，两类主张不可混为一谈。
7. **coloc-SuSiE 校准失败已披露**（PCSK9 低频 PIP=0；APOC3×T2D PP.H4=1.0 为 LD 错配伪可信集）→ 多信号稳健性缺失。
8. **APOC3×T2D/FBG "救回"不成立为独立发现**：15 对零假设期望 0.75 个 p<0.05、无多重校正、两 coloc 方法判不共定位、
   FinnGen rs964184×T2D p=0.075 不显著、rs964184 在多基因脂质簇、UKB-PPP 无 APOC3（复现断头）。
9. **KCNJ11 不是"转录介导"**：基本不在全血表达；coloc 峰值 rs757110 在 ABCC8（S1369A）；"表达中介"主张无据，
   机制留胰岛 eQTL。
10. **FinnGen rs17716350 对齐 bug** 未修前，全样本方向一致率（94.7%/88.5%）不作头条，只报双显著子集。
11. **T2D 结局 N 口径**（655,666 vs 62,892）同一研究不同分母，已核清非矛盾，投稿前回填。
12. **eQTLGen-GWAS 样本重叠**未做敏感性（SMR 在重叠下偏假阳）——投稿前至少讨论。

---

## 6. 执行顺序建议

1. **第一步（本周，纯改稿）**：M1 措辞 + M2 rs17716350 + M3/M4/M5 数字更正与独立位点表。零资源争议。
2. **同步**：P3 HEIDI 全集（性价比最高，818 测试）跑起来；UKB-PPP 下载继续（已 3/8）。
3. **P3 完成后**：Fig 5 与"多数不一致"叙事立住，评估投 HMG 5-6.5。
4. **冲 eBioMedicine 再评估**：负对照标定（E1）+ 一个新发现（E2）到位后才谈，不承诺。

> **在 M1-M5 落地前，不投任何期刊。**
