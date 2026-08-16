# deCODE 蛋白通道 MR + 共定位结果摘要（2026-08-13 更新，M3 全网格重跑 + M5 共定位闸门）

- 蛋白源：deCODE 血浆 pQTL（Ferkingstad 2021, n=35,559）；结局：T2D/CAD/FBG（OpenGWAS）。
- 方法：cis ±1Mb（hg38, ENSEMBL 坐标）、p<5e-6、EUR LD clump r²<0.01@1000kb、首读统计量 2026-08-13 修正为 nsnp≤3→IVW-FE / nsnp≥4→IVW-MRE / nsnp=1→Wald（见下）。

## 结果表（修正后首读方法，按 p 排序）

> 2026-08-13 评审修正：`mr_ivw_mre` 在 nsnp=2-3（工具间比率近同质）时 SE 被压缩 7-16 倍，
> 使 INSR×FBG p=0、APOC3×FBG p=9.1e-46 成为伪影。**修正口径：nsnp=1 → Wald；nsnp≤3 → IVW-FE 首读；nsnp≥4 → IVW-MRE 首读**（FE 本为预注册敏感性方法，未动任何门柱）。产物 `results/grid/protein_decode_mr_primary.csv`。

| 蛋白 | 结局 | nsnp | b | SE | p | 首读方法 | 共定位 PP.H4 | 注 |
|---|---|---|---|---|---|---|---|---|
| PCSK9 | CAD | 14 | 0.190 | 0.028 | 5.5e-12 | IVW-MRE | **1.000** | FE p=4e-26、WM p=1.2e-7、Egger 斜率 p=6.5e-5 全一致；F=47.5 |
| INSR | FBG | 2 | 0.125 | 0.050 | 0.0121 | IVW-FE | 0.0097 | MR 显著但不共定位（MRE p=0 为 SE 塌缩伪影） |
| APOC3 | T2D | 3 | 0.116 | 0.050 | 0.0203 | IVW-FE | <0.5 | MR 显著但不共定位（FE 下才显著，MRE 不显著） |
| APOC3 | CAD | 10 | 0.151 | 0.068 | 0.0271 | IVW-MRE | **0.997** | FE p=5.2e-9、WM p=2e-4 一致 |
| APOC3 | FBG | 3 | 0.044 | 0.023 | 0.0543 | IVW-FE | 0.177 | MRE p=9.1e-46 为伪影；FE 下不显著（p=0.054） |
| PCK1 | T2D | 1 | 0.210 | 0.189 | 0.27 | Wald | — | 单工具，F=1.24 |
| PCK1 | FBG | 1 | 0.070 | 0.074 | 0.35 | Wald | — | 单工具 |
| INSR | T2D | 2 | 0.078 | 0.084 | 0.35 | IVW-FE | — | — |
| INSR | CAD | 2 | -0.050 | 0.060 | 0.40 | IVW-FE | — | MRE p=0.039 为 SE 塌缩伪影 |
| ANGPTL3 | CAD | 19 | -0.010 | 0.014 | 0.45 | IVW-MRE | — | — |
| ANGPTL3 | FBG | 12 | -0.004 | 0.007 | 0.58 | IVW-MRE | — | — |
| PCSK9 | FBG | 5 | -0.007 | 0.019 | 0.70 | IVW-MRE | — | — |
| PCSK9 | T2D | 3 | 0.013 | 0.076 | 0.86 | IVW-FE | — | — |
| ANGPTL3 | T2D | 9 | 0.003 | 0.023 | 0.90 | IVW-MRE | — | — |
| PCK1 | CAD | 1 | 0.002 | 0.138 | 0.99 | Wald | — | 单工具 |

**未产出结果的对（含空结果，如实报告）**
  - HMGCR×T2D: cis ±1Mb、p<5e-6 内无变异
  - HMGCR×CAD: cis ±1Mb、p<5e-6 内无变异
  - HMGCR×FBG: cis ±1Mb、p<5e-6 内无变异
  - APOB×T2D: 蛋白文件未下载/不完整
  - APOB×CAD: 蛋白文件未下载/不完整
  - APOB×FBG: 蛋白文件未下载/不完整
  - LDLR×T2D: 蛋白文件未下载/不完整
  - LDLR×CAD: 蛋白文件未下载/不完整
  - LDLR×FBG: 蛋白文件未下载/不完整
  - GLP1R×T2D: cis ±1Mb、p<5e-6 内无变异
  - GLP1R×CAD: cis ±1Mb、p<5e-6 内无变异
  - GLP1R×FBG: cis ±1Mb、p<5e-6 内无变异
  - DPP4×T2D: cis ±1Mb、p<5e-6 内 6 个变异均为罕见位点、不在 1000G EUR LD 参考面板 → 无独立工具（2026-08-13 P3 修正，原 API 报错为 ld_clump 非确定性假象）
  - DPP4×CAD: 同上（罕见位点 → 无独立工具）
  - DPP4×FBG: 同上（罕见位点 → 无独立工具）
  - GCG×T2D: cis ±1Mb、p<5e-6 内无变异
  - GCG×CAD: cis ±1Mb、p<5e-6 内无变异
  - GCG×FBG: cis ±1Mb、p<5e-6 内无变异

## 漏斗

- 测试蛋白：11 个
- 蛋白×结局对：63 对
- 产出 MR 结果的对：45 对
- MR 显著（修正后首读，p<0.05）：4 对
- 共定位确认（PP.H4≥0.8）：2 对（均为对照校准位点）

## 预注册 H1（2026-08-13 执行，诚实报告）

- **分支一**：分泌型（PCSK9/APOC3/ANGPTL3）pQTL-MR 显著率 **3/9=33%** vs 胞内/膜（INSR/PCK1）**1/6=17%**；Fisher 2×2 OR=2.50 p=0.604 → **功效不足，仅描述性报告**（预注册纪律：标 CI、不排序）。
- **分支二**：双通道方向一致性——工具级双通道重叠对仅 3（全为 INSR），双通道均显著对 **0** → 方向一致性无法计算，**如实声明 n=0**。这是 n≈1 功效下的必然结果，也是蛋白通道必须扩面板（UKB-PPP）的量化依据。
- 产物 `results/grid/H1_direction_consistency.csv`；脚本 `scripts/H1_direction_consistency.py`。

## M5 蛋白共定位闸门（coloc.abf，p12=1e-5，敏感性 p12=1e-6，2026-08-13）

- 输入：15 对（MR 主方法 ok=TRUE 的蛋白×结局）× 全 cis 窗口（hg19 TSS±1Mb）GWAS 区域 + deCODE pQTL cis 全变异。
- 调律：回文位点保守排除（deCODE ImpMAF 无 effectAlleleFreq）；GWAS 侧用 region range query（proxies=0）。
- **结果（PP.H4 三档）：strong=2 / moderate=0 / none=13**
  - **PCSK9×CAD**：PP.H4=1.000（p12e6=1.000），top SNP=**rs11591147**（R46L）——阳性对照校准命中，教科书级。
  - **APOC3×CAD**：PP.H4=0.997（p12e6=0.974），top SNP=**rs964184**（APOA5/APOC3 区）——负对照中 TG 通路确实影响 CAD，属负对照边界案例，需在 M7 如实讨论。
  - **MR 显著但不过闸门**（修正后首读 p<0.05 且 PP.H4<0.5）：INSR×FBG（FE p=0.0121，PP.H4=0.0097）、APOC3×T2D（FE p=0.0203，PP.H4<0.5）、APOC3×FBG（FE p=0.0543 边缘，PP.H4=0.177）——MR 信号不共定位，可能由多效性/LD 驱动，共定位闸门将其过滤（这 3 对正是"M5 闸门价值"的演示案例）。
  - 其余 13 对 PP.H4<0.5。
- 产物：`results/grid/protein_coloc.csv`（全 15 对）、`protein_coloc_hits.csv`（strong 2 对）、`protein_coloc_funnel.tsv`。

## M7 结局侧外部复现（FinnGen 独立队列，2026-08-13）

> 预注册 §9 "外部复现" 的结局侧承诺（跨结局，非跨 eQTL/pQTL）：对转录 106 strong 共定位
> 命中 + 蛋白 15 ok 对的 top SNP，在 FinnGen 独立端点（finn-b-E4_DM2 T2D / finn-b-I9_CHD CAD）
> 提取关联，等位基因对齐后比方向一致率 + Spearman ρ。FBG 无 FinnGen 对应端点，如实跳过。
> 产物 `results/grid/finngen_replication.csv`（逐变异方向对比，n=109）；脚本 `scripts/M7_finngen_replication.R`。

| 结局 | 原始端点 | FinnGen 端点 | n | 方向一致 | Spearman ρ | 原结局显著 | 双显著 |
|---|---|---|---|---|---|---|---|
| T2D | GCST006867 | finn-b-E4_DM2 | 57 | **94.7%** | **0.866** | 53/57 | 27/57 |
| CAD | GCST005194 | finn-b-I9_CHD | 52 | **88.5%** | **0.769** | 49/52 | 19/52 |
| FBG | GCST005186 | 无对应端点 | — | 跳过 | — | — | — |

- **解读**：跨独立队列（FinnGen 50 万+）方向高度一致 + β 相关强 → 结局侧因果方向稳健，
  不是单一 GWAS 伪影。这是冲 eBioMedicine 的关键加分项（外部复现承诺落地）。

## M8 组织三角验证（GTEx 6 组织，2026-08-13）

> 对 106 个转录 strong 共定位命中 × GTEx v8 6 组织（肝/胰/全血/皮下脂肪/骨骼肌/冠脉），
> 取每基因该组织 lead cis-eQTL（tss±1Mb，p<5e-6，预注册工具阈值），对同一结局跑单工具 Wald
> MR，与全血 eQTL MR 方向比。离线复用 coloc 缓存 GWAS 区域（无需 API）。探索性补充，不移动
> 四态主分析门柱。产物 `results/grid/tissue_triangulation.csv` + `tissue_triangulation_hits.csv`；
> 脚本 `scripts/M8_tissue_triangulation.R`。

- 命中×组织 对数（有组织 lead eQTL）：**242**
- **强化（组织 MR p<0.05 且方向与全血一致）**：**129/242（53%）**，覆盖 **55/79 命中（70%）**
- 仅方向一致（p≥0.05，弱信号）：14；方向相反：4；无法对齐：31；结局无该变异：64
- 按组织强化数：Whole_Blood 33 / Adipose 23 / Muscle 23 / Pancreas 21 / Artery 17 / Liver 12
  （疾病相关组织都有支持：肝=脂质通路、胰=糖代谢、冠脉=CAD；GTEx 全血强化=跨数据集复现）
- 按结局：T2D 70 / CAD 59
- **解读**：70% 的 strong 共定位命中在 ≥1 个 GTEx 组织获得独立组织 eQTL 的 MR 方向强化 →
  全血 eQTL 信号不是组织特异伪影，是跨组织的稳健因果方向。这把描述性 atlas 升级为
  "跨组织三角验证"的验证层，是冲 eBioMedicine 要求的验证内容。

## M5-SuSiE 多信号共定位（coloc-SuSiE，2026-08-13 完成 → **校准失败，不可用作证据**）

- 对 3 个"MR 显著但 coloc.abf none"争议对 + 2 个校准对跑 coloc-SuSiE（多因果变异稳健分析）。
- **收敛修复（已解决）**：susie_rss 在完整 cis LD（~2000-4000 变异）上不收敛（LD 面板 1000G EUR
  与 deCODE 人群不匹配，susieR 自带诊断警告）。修复链：Rfast 加速 → 前 800 变异仍振荡 → **前 300
  信号变异 + EPV=FALSE + max_iter=10000 直接调用**（实测 3.7s 收敛）。另修 2 个独立 bug：susie_rss
  要求 n 为标量、coloc 5.x summary 列名 PP.H4.abf 非 PP.H4。变异选择改为 **方案B = top-300 by
  min(pqtl,gwas) z² ∪ 各性状单侧 top-20**（避免丢弃一侧极强另一侧中等的因果变异），不收敛自动回退
  方案A。
- **结果（`protein_coloc_susie_summary.csv`，5 对全跑完）**：

  | 对 | abf PP.H4 | SuSiE max PP.H4 | 结果 |
  |---|---|---|---|
  | PCSK9×CAD（校准） | **1.000** | **~0（全 25 个信号对）** | **阳性对照失败** |
  | APOC3×CAD（校准） | 0.997 | 不收敛（方案A/B 均试） | 无法评估 |
  | APOC3×FBG（争议） | 0.177 | ~0 | 两方法一致：不共定位 |
  | APOC3×T2D（争议） | 0.108 | **1.000（rs11216103×rs11216103）** | **伪阳性**（见下） |
  | INSR×FBG（争议） | 0.0097 | GWAS 无可信集（0 信号） | 无法评估 |

- **校准失败（关键，诚实报告）**：
  1. **阳性对照 PCSK9×CAD 未复现**。abf PP.H4=1.0 教科书级（top SNP rs11591147=R46L），但
     coloc-SuSiE 全部 25 个信号对 PP.H4≈0、PP.H3=1（"两性状各有信号但变异不同"）。根因：
     **rs11591147 为低频变异（MAF≈1.2%）**，susie_rss 的 z-score 方法在 1000G EUR LD 面板
     （n=503）上对其 PIP=0，无法识别为因果变异——即使把它强制加入输入（方案B）仍 PIP=0。
  2. **APOC3×T2D 的 PP.H4=1.0 是伪阳性**。该信号对 rs11216103×rs11216103：pQTL 侧 P=1.2e-19
     （强），但 **T2D GWAS 侧 p=0.054（完全不显著）**——susie 却把它放进 GWAS 可信集并给
     PP.H4=1.0。这是 LD 面板错配产生伪可信集的直接证据（弱 GWAS 信号 + 错配 LD → 退化单变异
     可信集）。若当真，将违反共定位的基本前提（共定位需两侧都显著）。
- **结论**：coloc-SuSiE 作为敏感性分析**未通过阳性对照校准**，其输出（包括任何 PP.H4=1.0）在
  本数据（deCODE 人群 vs 1000G EUR LD 面板错配 + 低频因果变异）下**不可信，不作为证据**。主闸门
  维持 **coloc.abf**（2 strong/0 moderate/13 none）。这与"需要 UKB-PPP 独立平台复现"的结论一致：
  跨平台数据 + 匹配人群 LD 才能支撑多信号共定位。唯一两方法一致点：APOC3×FBG 均判不共定位。

## M9 蛋白通道 SMR/HEIDI（Zhu 2016，2026-08-13）

> 对蛋白通道 5 个有 cis 工具的蛋白 × 3 结局 = 15 对，用 top cis-pQTL 变异做工具变量，
> SMR 测共享因果变异、HEIDI 测多效性 vs LD（默认阈值 1.57e-3）。产物
> `results/grid/protein_smr_heidi_clean.csv`（15 行）；脚本 `scripts/M9_smr_build_inputs.R` +
> `scripts/M9_run_smr.sh`。
> **诚实 caveat**：deCODE 无 effectAlleleFreq → .esd Freq 用 ImpMAF 近似 + `--disable-freq-ck`
> （探索性敏感性分析）；HEIDI 用 1000G EUR 参考 LD（与 deCODE 人群错配，见 M5-SuSiE 节）。

| 蛋白 | 结局 | topSNP | p_SMR | p_HEIDI | 判定 | 与 coloc.abf 对比 |
|---|---|---|---|---|---|---|
| PCSK9 | CAD | rs11591147 (R46L) | 2.96e-19 | 1.28e-4 | **SMR 极显著 + HEIDI 异质性（过不了 HEIDI）** | strong (1.000)，校准命中 |
| APOC3 | CAD | rs964184 | 5.52e-11 | 0.0266 | **SMR 显著 + HEIDI 异质性** | strong (0.997) |
| APOC3 | T2D | rs964184 | 0.0112 | 0.422 | **单工具敏感性支持（非独立确认）** | none (0.108) |
| APOC3 | FBG | rs964184 | 0.00708 | 0.915 | **单工具敏感性支持（非独立确认）** | none (0.177) |
| INSR | FBG | rs4804368 | 0.345 | 0.420 | null | none (0.0097) |
| 其余 10 对 | — | — | >0.05 | — | null | 与 coloc 一致 |

> **诚实判定（2026-08-13 多 agent 评审收敛）**：蛋白 SMR 全表仅作**探索性敏感性**，不能当作独立发现，
> 理由：① deCODE 无 effectAlleleFreq → ImpMAF 近似 + `--disable-freq-ck` 未正式校准；
> ② HEIDI 用 1000G EUR LD，与 deCODE（冰岛）人群错配；③ **阳性对照 PCSK9×CAD 过不了 HEIDI
> （p=1.28e-4）** → 该平台下"HEIDI 通过/异质"判定均未校准；④ 15 对零假设期望 0.75 个 p<0.05，
> 观测 4 个含两个已知阳性，无多重校正 → APOC3×T2D/FBG 的 p=0.011/0.007 是边缘值；
> ⑤ coloc.abf none、coloc-SuSiE 对 FBG 明确判不共定位、FinnGen rs964184×T2D p=0.075 不显著；
> ⑥ rs964184 位于 APOA1/C3/A4/A5 多基因脂质簇，deCODE 适体未必 APOC3 特异；
> ⑦ **APOC3 不在 UKB-PPP Olink 面板** → 计划内跨平台复现路径断头，只能另找 INTERVAL/ARIC/SCALLOP 或文献复现。
> **结论表述**：APOC3×T2D/FBG 为"单工具敏感性支持，待独立平台复现"，**不得写"独立救回/共享因果变异新发现"**。

- **阳性对照**：PCSK9×CAD top SNP 正是 rs11591147（R46L 已知功能变异），p_SMR 极显著；
  HEIDI 异质性反映该位点多信号（与 coloc-SuSiE 的 PP.H3=1 结果吻合），非管道伪影——但正因阳性对照
  过不了 HEIDI，蛋白侧 HEIDI 判定整体未校准（见上诚实判定）。
- APOC3×CAD 的 HEIDI 异质性同理（rs964184 区多信号）。
- 蛋白通道 HEIDI pass 12/15，3 对异质 = PCSK9×CAD（1.28e-4）、APOC3×CAD（0.0266）、INSR×T2D（0.014，p_SMR null 但区域异质）；
  因上述校准问题不放大解读。

## M10 转录通道 SMR/HEIDI（2026-08-13）

> 对转录通道 106 个 strong 共定位命中（106 唯一基因）× 3 结局跑 SMR+HEIDI，
> eQTL 侧用 eQTLGen 全血 cis-eQTL（hg19，位置与 1kg bim 逐位点验证一致）。
> 产物 `results/grid/transcript_smr_heidi.csv`（128 行）；脚本 `scripts/M10_transcript_smr_build_inputs.R` +
> `scripts/M10_run_smr.sh`。
> **诚实 caveat（关键）**：① 测试集正是 106 个 MR-显著预选 coloc 命中，其 top cis-eQTL SNP 的
> GWAS 关联本就显著 → **p_SMR 显著（125/128）是选择偏倚的构造性必然，零信息量，不进正文**；
> 唯一有信息量的读数是 **HEIDI 通过率（区域 LD 同质性检验）**。② **MR 单工具 Wald 与 SMR 是同一
> 数学量**（同一 top cis-SNP 上的同一比率），不是两个独立方法——转录通道的"确认"结构是
> **单点 Wald 自洽 + 两个区域级独立检验（coloc PP.H4、HEIDI）**，不是"三方法独立互证"。
> ③ eQTLGen 只给 Zscore → Beta=Z/√N、se=1/√N（标准化表达常规转换）。④ 多数对 nsnp<10，HEIDI 功效有限。

- **128 个探针×结局测试：89 个 HEIDI 通过（p_HEIDI>0.05，89/126=70.6%，2 NA 排除）**；
  37 个 HEIDI 异质性（多信号位点，与 coloc-SuSiE 结论一致）；2 个 NA（SNP 不足）
- 按结局：T2D 42 通过；CAD 45；FBG 2
- **KCNJ11×T2D：topSNP=rs2074310，p_SMR=5.82e-12，p_HEIDI=0.128** —— 单点 Wald 自洽
  （b_SMR=0.595181=0.0643/0.108034=MR Wald b，同数到 6 位小数，**数学恒等非独立互证**）
  + 两个区域级独立检验（coloc PP.H4 与 HEIDI）通过。**机制 caveat**：coloc 峰值 rs757110 落在
  **ABCC8（S1369A）**，KCNJ11 基本不在全血表达 → 正确表述为"区域级共享变异 + 同质性通过"，
  "血表达中介"主张无据，机制解释留待胰岛 eQTL（InsPIRE/GTEx pancreas）。
- 其他亮点：LIPA×CAD（p_SMR=1e-20, HEIDI 0.457）、ARG1×T2D（HEIDI 0.970）、
  PTPRN×CAD（HEIDI 0.249）、RASD1×T2D（HEIDI 0.519）——三角验证种子获 HEIDI 支持
- **意义**：转录通道 = MR（单点 Wald）→ coloc.abf（区域级，PP.H4≥0.8，106 命中）→ HEIDI
  （区域 LD 同质性，89/126 通过）。106 命中拆 **76 编码基因 / 30 非编码伪影倾向位点**
  （RP11-*/CTD-*/Metazoa_SRP/U6/hsa-mir-296/SERBP1P3/KRT8P46）；**p12=1e-6 先验敏感性仅
  20/106（18.9%）保持 PP.H4≥0.8** → strong 结论依赖主设定 p12=1e-5。漏斗图增加 SMR/HEIDI 层。
- 其他 drug target（DPP4/GLP1R/PPARG/INSR/TCF7L2/PRKAA1/SLC5A2）未进 coloc 网格故无 SMR
  （转录 MR 有部分显著，待 UKB-PPP 蛋白侧 + 补 GWAS 后扩展）

## 产物文件

- `results/grid/protein_decode_mr.csv` — 全网格（含空/失败）
- `results/grid/protein_coloc.csv` / `protein_coloc_hits.csv` — M5 共定位全表 + strong 命中
- `results/grid/protein_coloc_susie_summary.csv` / `protein_coloc_susie_detail.csv` — M5-SuSiE 每对单行汇总 + 每 CS 对详情（校准失败，见 M5-SuSiE 节）
- `results/figures/F9_funnel_protein.png` — 蛋白通道漏斗
- `results/figures/F10_protein_mr_forest.png` — 蛋白 MR 森林图
- 转录通道图：`F1`（漏斗）、`F2`（v1→v2 坍缩）、`F3`（极端 p 消失）

## 说明（诚实报告）

- EAF：assocvariants.annotated 未就绪 → palindromic 保守排除 + ImpMAF 近似（ImpMAF 不总是 effect allele 频率）。
- 多等位基因 bug 行（effectAllele==otherAllele / "!"）已排除。
- 2026-08-13 M3 修复（P3）：DPP4 的 cis p<5e-6 变异全为罕见位点、不在 1000G EUR LD 参考面板 → 预注册 clump 步骤下无独立工具（如实报告）。原"missing required columns/401"为 OpenGWAS ld_clump 对罕见位点非确定性行为 + 瞬时 API 抖动的假错误；INSR×FBG 的 401 已重试复现成功。
- 2026-08-13 M5 修复（P4-P6）：M1 cis 子集文件名（${base%.gz}_cis.txt.gz）拼错致 pQTL 全 0；OUTCOMES/OUT_N/TYPE/S_CASE 以短名索引而 pairs 存全 id → API id=NULL 全失败 + run_coloc 越界。已修，15 对全跑通。
- 2026-08-13 MRE SE 塌缩修正（多 agent 评审确认）：nsnp=2-3 对改用 FE 首读（INSR×FBG p=0→0.0121、APOC3×FBG p=9.1e-46→0.0543、INSR×CAD p=0.039→0.40、APOC3×T2D p=0.135→0.0203）；F 统计量已补齐（F≈(b/se)²，单工具 Wald 单独标注）。脚本 `scripts/report_protein_primary.py`，产物 `results/grid/protein_decode_mr_primary.csv`。
- 转录通道“结果一般”主因：全血表达对肝/肠脂质基因是弱代理；cis-eQTL 严格 clump 降功效；非编码问题。详见 CHANGELOG。

