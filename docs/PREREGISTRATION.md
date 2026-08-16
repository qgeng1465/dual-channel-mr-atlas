# 预注册声明 / 内部设计锁定记录 (Pre-registration / design-lock record)

项目: dual-channel MR atlas（cis-MR×coloc 双通道/全转录组图谱）

> **v7 定位澄清（2026-08-16，本文档当前版本）**：本文档为**分析启动前的内部设计锁定 + 运行期
> 执行纪律记录**，非第三方平台（OSF/AsPredicted）的时间戳前瞻性注册。分析已全部完成（M20–M30，
> 2026-08-15/16），**不再补作第三方预注册**，亦不将本文档表述为正式预注册。投稿的 Data
> Availability 声明见 `docs/DATA_AVAILABILITY_20260816.md` §1（诚实措辞逐字可用）。
> 本次 v7 修订仅为**删除历史重复段落 + 定位澄清**，不改变任何分析门柱；哈希锁重算。
> 修订历史：v1 2026-08-06（双通道核心）→ v5 2026-08-07（转录组全量扫描层）→ v6 2026-08-13
> （执行记录澄清）→ v7 2026-08-16（去重+定位）。

---

## v1 — 原始预注册（2026-08-06 13:35:32 CST，双通道核心）

版本: 0.1.0

**v1 由审稿人 agent 有条件放行(MAJOR_REVISION)后，在运行任何分析之前生成。
  md5 哈希锁定（文件名 .sha256 实为 MD5，见 v6 §6.3 澄清）。分析后任何修改将导致
  后续阶段完整性校验失败（M3/M15/M16 等脚本运行时校验）。**

### 1. 研究问题

全血细胞转录本(cis-eQTL)与循环蛋白(cis-pQTL)两个介质层，哪个更能作为
代谢疾病(T2D/CAD)的干预介质层？——回答的是 **因果介质物种** 而非药物模态。

### 2. 暴露通道与数据

- 通道A(转录本): eQTLGen 全血 cis-eQTL, n=31684
- 通道B(蛋白): deCODE 血浆 cis-pQTL, n=35559
- 备用蛋白源: INTERVAL SomaScan pQTL via OpenGWAS prot-a-*（若 deCODE 不可得）
- 结局: ieu-b-38(T2D), ieu-a-7(CAD), finn-b-E4_DM2(T2D_FINNGEN)

### 3. 工具变量标准（锁定）

cis 窗口 ±1000kb，p<5e-06
LD clumping r²<0.01 @ 1000kb (EUR (1000G Phase 3))，
palindromic 处理 action=2

### 4. MR 方法（锁定）

主方法: mr_ivw_mre
敏感性: mr_ivw_fe, mr_weighted_median, mr_egger_regression
多重检验: BH-FDR q<0.05

### 5. 方向四态分类（核心统计量，取代 Δ 点值）

早期版本字段曾迭代（concordant / protein_only / transcript_only / both_null，为 v1 迭代残留）；
定稿口径以 CI 是否跨 0 与 p 值双重判定，CI/PP 为基础，不做 |Δ| 跨基因排序。

### 6. 正/负对照（审稿人核定）

- 正对照(校准): PCSK9, CETP —— 仅用于确认 pQTL 可定位 + 药物可及性校准，不用于药物模态预测
- 负对照: HMGCR, NPC1L1, ANGPTL3, APOC3 —— 留一法 LOOCV AUR；若 AUR≤0.5 仅做描述性报告（预注册约定）

### 7. 预注册假设（可证伪）

- **H1**: 分泌型蛋白 vs 细胞内蛋白 在 pQTL-MR 显著性率与方向一致性上存在差异.
  检验: Fisher/χ² on 2×2 (secreted vs intracellular × significant-yes/no), α=0.05 two-sided
- **H2**: 负对照蛋白 LOOCV 分类器 AUR > 0.5.
  检验: leave-one-out AUC on 4 negative-control proteins, α=0.05 one-sided (AUR>0.5)

### 8. 种子

master = 20260805（全部下游随机从 master 确定性派生）

### 9. 报告纪律（学术不端防护）

- 全网格结果全部落盘报告，**无论是否显著**，含空结果与失败记录
- 禁止事后调整阈值/方法以追求显著（p-hacking）
- 阴性结果如实报告，不做阳性概率的夸大表述（审稿人否决该表述）
- 每个数据源记录 ID/URL/日期，全部可追溯
- 复现: UKB-PPP pQTL (Sun 2023, 54,306 UKB) / INTERVAL independent split

### 10. 论文定位（审稿人核定）

目标: eBioMedicine(天花板/IF≈9, CAS Q1) → Human Genetics(现实, IF≈3.9) →
J Lipid Res(备选, IF≈3.6)。产品: 干预介质层优先级表。
（注：2026-08-16 期刊定位已定案为**主投 AJHG**，其余期刊推测已删除，见 JOURNAL_TARGETS。）

---

## v5 — 修订（2026-08-07）：转录组全量扫描层（transcriptome-wide cis-MR scan）

目标：在候选基因深度层之外，补全"转录通道全量网格"——eQTLGen 全部有显著 cis-eQTL 的基因
（预计 ~16,000+ 基因）对 T2D/CAD/FBG 的系统扫描，回答"全血转录本层面哪些基因有 MR 证据"，
为论文提供系统扫描层（假设生成）。候选基因深度层（四态分类、药物靶点、组织敏感性）不移动。

设计（对 §3 工具标准的部分偏离，**仅限扫描层**，预先登记）：
- 工具 = 每基因最强 cis-eQTL（lead variant：该基因 cis ±1Mb、Pvalue<5e-6 内 Pvalue 最小的单变异）
- 方法 = 单工具 Wald ratio（nsnp=1）；主分析 mr_ivw_mre 仍用于深度层
- 不逐基因 LD clump（~16k 基因 × API clump 不可行、耗时失控）；命中基因（FDR q<0.05）进入
  第二阶段：本地 plink（1000G EUR，1kg.v3 参考）clump r²<0.01@1000kb + IVW-MRE 复核
- 结局提取 proxies=FALSE（扫描层避免代理 SNP 引入伪阳性；深度层保持 proxies=TRUE 不变）
- 多重检验：按结局对全基因 Wald p 做 BH-FDR q<0.05
- 报告：全网格落盘（含无工具/无结局匹配/失败行）；命中表如实报告数量、方向与 FDR q
- 严格性声明：扫描层是 hypothesis-generating；结论级证据以深度层为准（不移动门柱）

---

## v6 — 修订（2026-08-13）：实际执行记录澄清（结局 ID / 样本量口径 / 首读统计量 / H1 执行状态）

v6 为**运行期执行记录澄清**，不移动任何主分析门柱（IV 标准、MR 方法集、FDR、正/负对照均不变）。
内容为把"预注册文档 §2 列出的结局 ID"与"实际执行所用结局"之间已发生的偏差、以及运行期
诚实报告口径，以追加修订方式如实记录并重新哈希。

### 6.1 结局集实际执行记录（预注册 §2 的偏差，如实记录）

预注册 §2 列：`ieu-b-38`(T2D), `ieu-a-7`(CAD), `finn-b-E4_DM2`(T2D_FINNGEN)。

实际执行结局（全网格与全部 MR/共定位均使用）：

| 结局 | 实际 ID | 效应 | gwasinfo sample_size | gwasinfo ncase/ncontrol |
|---|---|---|---|---|
| T2D | `ebi-a-GCST006867` | logOR | 655,666 | 61,714 / 1,178 |
| CAD | `ebi-a-GCST005194` | logOR | 296,525 | 34,541 / 261,984 |
| FBG | `ebi-a-GCST005186` | beta | 58,074 | —（定量） |

- `ieu-b-38` / `ieu-a-7` 在扫描运行期（2026-08-07）经 OpenGWAS 查询不可用/覆盖不足，
  实际改用 `ebi-a-GCST006867`/`ebi-a-GCST005194`（相同表型，GCST 系列，人群 EUR）。
- `finn-b-E4_DM2`（FinnGen T2D）未进入主网格；保留为**结局侧外部复现**目标（见 6.6）。
- **FBG（`ebi-a-GCST005186`）为预注册 §2 结局集之外的**新增结局，作为连续代谢表型补充；
  属预先登记式新增（探索性补充），不改变预注册主假设 H1/H2。

### 6.2 T2D 样本量口径澄清（CHANGELOG 记录 62,892 vs README 655,666）

非矛盾：同一研究（Mahajan 2018, GCST006867）的两套元数据口径。
- `sample_size=655,666` = 汇总样本量（OpenGWAS metadata，README §2 采用）。
- `ncase=61,714 + ncontrol=1,178 = 62,892` = 直接病例对照计数（gwasinfo）。
- 论文统一采用 sample_size=655,666，并在方法节注明 ncase=61,714。

### 6.3 哈希锁澄清（文件名 .sha256 实为 MD5）

`docs/PREREGISTRATION.md.sha256` 文件名为历史遗留，**内容为 32 位十六进制 MD5**；
脚本校验用 `tools::md5sum(prereg) == readLines(lock)`（如 M5_protein_coloc.R 第 47 行）。
不重命名（避免破坏现有校验脚本）。

### 6.4 首读统计量报告口径（不移动门柱的报告修正）

多 agent 评审（对抗验证可复现）发现 `mr_ivw_mre` 在 nsnp=2-3（工具间比率近同质）时 SE 被
压缩 7-16 倍（如 deCODE INSR×FBG MRE se=0.003 vs FE se=0.050）。预注册 §4 的敏感性方法
`mr_ivw_fe` 即为预设的固定效应对照。报告口径修正（**仅影响报告，不改变方法集与门柱**）：
- nsnp=1 → Wald ratio（不变）
- **nsnp≤3 → 以 IVW-FE 为首读统计量**
- nsnp≥4 → 以 IVW-MRE 为首读统计量
- 全方法（MRE/FE/WM/Egger）仍全量报告于 `results/grid/protein_decode_mr.csv`；
  首读表单独落盘 `protein_decode_mr_primary.csv`，F 统计量一并补齐（F≈(β/se)²）。

### 6.5 H1/H2 执行状态（预注册 §7）

- **H1（分泌型 vs 胞内 pQTL-MR 显著率）已执行**：分泌型 3/9=33% vs 胞内/膜 1/6=17%，
  Fisher 2×2 OR=2.50 p=0.604 → 功效不足，按预注册纪律仅描述性报告（标 CI、不排序）。
  产物 `results/grid/H1_direction_consistency.csv`。
- **H1 第二支（双通道方向一致性）已执行**：工具级双通道重叠对仅 3（全为 INSR）、双通道
  均显著对 0 → 方向一致性无法计算，如实声明 n=0（功效不足的必然结果）。
- **H2（负对照 LOOCV AUR）尚未执行**：负对照 HMGCR 无 cis p<5e-6 工具、NPC1L1 deCODE 文件
  未下载，仅 ANGPTL3/APOC3 可用（且 APOC3×CAD 为边界正信号）。待续下 NPC1L1/CETP 文件或
  按 CHANGELOG:114 设计调整后再执行；不事后改主门柱。

### 6.6 外部复现：结局侧 FinnGen（预注册 §9 承诺执行）

对转录 106 strong 共定位命中 + 蛋白 15 ok 对的 top SNP，在 FinnGen 独立端点
（finn-b-E4_DM2 T2D / finn-b-I9_CHD CAD）提取关联，等位基因对齐（ea2==oa1 → 翻号）后比方向。
- T2D：n=57，方向一致 **94.7%**，Spearman ρ=0.866，双显著 27/57。
- CAD：n=52，方向一致 **88.5%**，Spearman ρ=0.769，双显著 19/52。
- FBG 无 FinnGen 对应端点 → 如实跳过。
- 产物 `results/grid/finngen_replication.csv`（n=109）；脚本 `scripts/M7_finngen_replication.R`。
- 口径：等位对齐失败的行不计入一致率；探索性补充，不移动主门柱。

### 6.7 组织三角验证（GTEx 6 组织，探索性补充）

对 106 个转录 strong 共定位命中 × GTEx v8 6 组织，取每基因该组织 lead cis-eQTL
（tss±1Mb、p<5e-6，同预注册工具阈值），对同一结局跑单工具 Wald MR，与全血 eQTL MR 比方向。
- 强化（组织 MR p<0.05 且方向一致）：**129/242 对（53%），覆盖 55/79 命中（70%）**。
- 疾病相关组织均有支持：肝 12 / 胰 21 / 冠脉 17 / 脂肪 23 / 肌肉 23 / GTEx 全血 33（跨数据集复现）。
- 产物 `results/grid/tissue_triangulation.csv` + `tissue_triangulation_hits.csv`；脚本 `scripts/M8_tissue_triangulation.R`。
- 探索性补充，不移动四态主分析门柱。

### 6.8 coloc-SuSiE 多信号稳健分析（探索性，校准失败，不作为证据）

对 3 个"MR 显著但 coloc.abf none"争议对 + 2 个校准对跑 coloc-SuSiE。收敛问题已解决：
- susie_rss 默认 EPV=TRUE 在完整 cis LD（~2000-4000 变异）上不收敛（LD 面板 1000G EUR 与
  deCODE 人群不匹配，susieR 自带诊断警告）→ 前 300 信号变异 + EPV=FALSE + Rfast +
  max_iter=10000 直接调用（实测 300 变异 864 次迭代 3.7s 收敛）。变异选择方案B
  （top-300 by min-z² ∪ 各性状单侧 top-20）保因果变异，不收敛自动回退方案A。
- 另修 2 个独立 bug：susie_rss 要求 n 为标量（`condition length>1`）、coloc 5.x summary
  列名 PP.H4.abf 非 PP.H4。
- **结果（5 对全跑完）与校准失败**：
  - 阳性对照 **PCSK9×CAD（abf PP.H4=1.0，rs11591147=R46L 低频）未复现**：coloc-SuSiE 全信号对
    PP.H4≈0（rs11591147 MAF≈1.2%，susie_rss 在 1000G EUR LD 上 PIP=0，无法识别低频因果变异）。
  - **APOC3×T2D 的 PP.H4=1.0（rs11216103×rs11216103）为伪阳性**：该变异 T2D GWAS 侧 p=0.054
    不显著却被 susie 放入 GWAS 可信集 → LD 面板错配产生退化单变异可信集。
  - APOC3×CAD 不收敛、INSR×FBG GWAS 无可信集、APOC3×FBG 两方法一致判不共定位。
- **执行记录**：coloc-SuSiE 未通过阳性对照校准，其输出在本数据（deCODE 人群 vs 1000G EUR LD
  错配 + 低频因果变异）下不可信，**不作为证据**。主闸门维持 coloc.abf。多信号共定位留待
  UKB-PPP 独立平台 + 匹配人群 LD 再评估。探索性稳健分析。

### 6.9 声明

以上 v6 修订全部为运行期**执行记录澄清/报告口径说明**，不改变预注册主分析门柱
（工具变量标准 §3、MR 方法集 §4、方向四态 §5、正/负对照 §6、假设 §7）。
重新哈希锁：2026-08-13。

---

## 哈希锁

- 锁文件：`docs/PREREGISTRATION.md.sha256`（内容为 MD5，见 §6.3）
- 当前版本：v7（2026-08-16，去重+定位澄清，不改变门柱）
- 校验方式：`tools::md5sum(PREREGISTRATION.md) == readLines(PREREGISTRATION.md.sha256)`
