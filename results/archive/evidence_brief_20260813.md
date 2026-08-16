# 双通道 MR 图谱：当前证据状态简报（2026-08-13，供 agent 评估）

> 用途：给多 agent 评审工作的输入。目标是回答：(1) 现在能否达到"阳性结果"？
> (2) 如何串起现有研究讲一个可发表的、够 IF≥9 的故事？

## 1. 项目定位

双样本 cis-MR 图谱：测试共享 cis-variant 因果信号是否通过**转录**（eQTLGen 全血）或
**蛋白**（deCODE 血浆 pQTL）层作用于 T2D/CAD/FBG。输出基因×结局×介质层优先级图谱。
预注册（docs/PREREGISTRATION.md，MD5 哈希锁校验通过）。全部结果如实报告（含空/阴性）。

## 2. 转录通道（最强通道）

- 漏斗：scan 16,622 基因×3 结局 → FDR 存活 982 → stage-2 906 基因/812 FDR 存活 → **106 个 strong 共定位（PP.H4≥0.8）**（106 唯一基因，拆 **76 编码 / 30 非编码伪影倾向**；PP.H4≥0.9 共 62）
- 证据结构（2026-08-13 多 agent 核查收敛，**措辞已诚实化**）：
  - **单点 Wald 自洽 + 两个区域级独立检验**：MR（预注册，单工具 Wald）→ coloc.abf（区域级，PP.H4≥0.8）+ HEIDI（区域 LD 同质性）
  - **注意：MR 单工具 Wald ≡ SMR 是数学恒等**（同一 top cis-SNP 的同一比率），不是两个独立方法，不得写"三方法独立互证"
  - SMR/HEIDI 结果：128 个探针×结局测试，**HEIDI 通过 89/126=70.6%**（2 NA）；37 个 HEIDI 异质性；**p_SMR 显著（125/128）是选择偏倚构造性必然，零信息量不进正文**
- **KCNJ11×T2D：MR Wald p=4e-16 + coloc PP.H4 + HEIDI 0.128 通过** = 单点自洽 + 两区域级检验过；
  **机制 caveat**：coloc 峰值 rs757110 在 ABCC8（S1369A），KCNJ11 基本不在全血表达 → 只可称"区域级共享变异"，非"血表达中介"
- 其他 drug target（DPP4/GLP1R/PPARG/INSR/TCF7L2/PRKAA1/SLC5A2）：转录 MR 部分显著（GLP1R×T2D FE p=0.002、
  PPARG×T2D FE p=0.0003、DPP4×T2D p=0.003、INSR×T2D p=0.001），但未做 coloc/SMR（无 GWAS 缓存）
- 外部复现：
  - **FinnGen 结局侧**：T2D 方向一致 94.7%（ρ=0.866, n=57）、CAD 88.5%（ρ=0.769, n=52）。
    **已知 bug 未修**：rs17716350 同 SNP 两行 FinnGen beta 一负一正（-0.0133(G) vs +0.1344(T)）均标 aligned=TRUE
    （`finngen_replication.csv` 30-31 行，恰为 C15orf62 top SNP）→ 去重后 CAD 46/51=90.2%；
    投稿前须修对齐器 + 全表重算，头条只报双显著子集 27/27、19/19
  - **组织三角验证**：106 命中 × GTEx 6 组织，55/79（70%）获 ≥1 组织 MR 方向强化
    （肝/胰/冠脉/脂肪/肌肉/全血全有支持）

## 3. 蛋白通道（弱但诚实）

- 11 候选蛋白 → 5 有 deCODE 工具 → MR 修正后显著 4 对 → coloc.abf strong 2 对（均校准位点）：
  PCSK9×CAD（PP.H4=1.000, rs11591147 R46L，教科书阳性对照）、APOC3×CAD（0.997）
- **SMR/HEIDI 结果（15 对，仅探索性、未校准）**：
  - **APOC3×T2D：p_SMR=0.011, p_HEIDI=0.42 → 单工具敏感性支持（非独立确认）**
  - **APOC3×FBG：p_SMR=0.007, p_HEIDI=0.91 → 单工具敏感性支持（非独立确认）**
  - PCSK9×CAD：p_SMR=2.96e-19 极显著，但 HEIDI=1.28e-4 异质性（多信号位点）
  - APOC3×CAD：p_SMR=5.5e-11，HEIDI=0.027 异质性
  - 其余 11 对 null（与 coloc 一致）
- **诚实解读（不夸大）**：15 对零假设期望 0.75 个 p<0.05、无多重校正；coloc.abf none、coloc-SuSiE 判不共定位、
  FinnGen rs964184×T2D p=0.075 不显著、**APOC3 不在 UKB-PPP 面板（复现断头）**；
  **阳性对照 PCSK9×CAD 过不了 HEIDI（p=1.28e-4）→ 蛋白侧 HEIDI 判定整体未校准**。
  故 APOC3→T2D/FBG 只能作为"待独立平台复现的单工具敏感性"写入，**不得写"独立救回/共享因果变异新发现"**
- H1 分支一（分泌 vs 胞内）功效不足（n≈1）仅描述性；分支二双通道方向一致性 n=0（如实声明）

## 4. 方法学诚实声明（必须保留）

- MRE SE 塌缩已修正（nsnp≤3 → FE 首读；INSR×FBG p=0→0.0121）
- **coloc-SuSiE 校准失败**（PCSK9×CAD 阳性对照未复现：rs11591147 低频 + 1000G LD 错配 → PIP=0；
  APOC3×T2D 的 PP.H4=1.0 为伪可信集）→ 不作为证据，主闸门维持 coloc.abf
- HEIDI 用 1000G EUR 参考 LD，与 deCODE 人群错配的 caveat（转录 eQTLGen 同为 EUR 人群更可靠）
- deCODE 无 effectAlleleFreq → 频率用 ImpMAF 近似 + --disable-freq-ck（探索性）
- T2D N 矛盾已核清（655,666 vs 62,892 同一研究不同口径）

## 5. UKB-PPP 跨平台复现（进行中，最高杠杆）

- 用户已注册 Synapse + PAT 已存 ~/.synapseConfig（gitignored）
- Synapse API 被中国大陆 GEO 屏蔽 → 走 mihomo 代理（用户已被告知无镜像、直连被限速 100 倍）
- **面板覆盖核查**：UKB-PPP 有 PCSK9/INSR/ANGPTL3/LDLR/APOB/DPP4/GLP1R/GCG（8 个）；
  **无 APOC3/PCK1/HMGCR**（Olink 面板不覆盖，HMGCR 的"无工具"仍是数据缺口）
- 8 个蛋白 tar 正后台下载（~4.4GB，~38min），下载后按 cis±1Mb 提取 → MR + coloc + 方向一致率

## 6. 待 agent 裁决的问题

1. 现在是否构成"阳性结果"？（单点自洽 + 两个区域级独立检验 + FinnGen + 组织，是否构成可辩护的多层收敛证据？——多 agent 裁决：转录通道是，但非新发现型，见结论）
2. 主故事应该是什么？候选（2026-08-13 多 agent 已裁决，见下）：
   a. 方法学：单点 Wald 自洽 + 两个区域级独立检验（coloc/HEIDI）的"coloc 一致性率 106/819=12.9%[Wilson CI 10.8–15.4%]" + MR≡SMR 校准声明（唯一可辩护新贡献，需负对照标定后才可谈假阳性）
   b. 药物靶点：KCNJ11/PCSK9/APOC3 作校准案例 + ~60-70 编码优先基因列表（优先化而非发现）
   c. 多组学图谱：转录图谱为主 + 蛋白如实次要层 + 跨组织/跨平台（整合 n=0，不得称"双通道整合"）
3. 每个故事缺什么证据、能投哪档期刊、怎么补？
4. MD 文档哪些过时该删/该更新？（DECODE_PIPELINE_SUMMARY / IMPROVEMENT_STRATEGY / PREREGISTRATION / README）
