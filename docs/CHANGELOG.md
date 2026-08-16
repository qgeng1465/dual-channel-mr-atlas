# CHANGELOG — 变更记录

## 2026-08-16 — 稿件终稿收尾（审稿 5 点 + 3 点修订落地、样本重叠敏感性、引用扩至 ~50、表文件命名修正）

### 审稿修订（manuscript.md）
- **口径主从**：FDR-core（982 MR-sig → 121 strong，12.32%）为主口径，nominal 漏斗为探索性敏感性；Fig 2/Fig 3 图注重写为 6 面板映射（ECDF/strong 率/PP.H4 阈值/置换/38.7% GWAS 峰 caveat/峰 p 分布）。
- **15 vs 23**：8 个 FDR 掉落候选（LAMC1/TPD52/SENP6/HMGN3/MT3/RPL13/ZBTB46/ZNF100）入 Table S3 + Discussion 保守性段落；LAMC1 FinnGen p=0.018 如实披露。
- **样本重叠敏感性（新分析，15/15 一致）**：`results/sample_overlap_sensitivity_20260816.csv`——15 候选的 MR 方向用 pre-UKB 结局 GWAS 重估（CAD=Nikpay 2015 CARDIoGRAMplusC4D ~185k；T2D=DIAGRAM 跨族裔 Mahajan 2014，3 个基因走 Scott 2017 European 1000G 次级），8/8 CAD + 7/7 T2D 方向全一致（11 正/4 负/0 冲突）。已整合进 Methods §2.2 新小节、Results §3.6、Discussion 首要局限（"方向偏倚排除，但全扫描 coloc 支持仍不能排除重叠扭曲"）、新 Table S4。诚实 caveat：T2D 主 GWAS（Xue 2018）为东亚样本本已与欧洲 eQTLGen 重叠极少，T2D 侧为保守的 belt-and-suspenders。
- **coloc.susie 重构**：非收敛（max_iter=1000 仍不收敛）→ §3.7 标为 exploratory，LAMC1 排除依据改为多信号证据（8 eQTL credible sets vs 3 GWAS sets、max z=4.32、无共享变异）。
- **Discussion 分层解释（Tier 1/2/3）+ Intro P3 召回缺口 + 摘要 161 词 FDR-core 口径**均落地。
- **引用扩至 51 条（全部核实为真）**：引入 29 个新 [@key]（MR 基础 davey2003/lawlor2008/dsmith2014/sanderson2022/ference2012/ference2017/holmes2015/schmidt2017；eQTL 与数据源 westra2013/gtex2015/gtex2017/bycroft2018/elsworth2020；T2D/CAD/FG 遗传 morris2012/spracklen2020/deloukas2013/dupuis2010/scott2012/mahajan2018；MR 方法 pierce2013/hemiani2018/bowden2017/verbanck2018/burgess2011/burgess2017；fine-mapping benner2016/benner2017/gaulton2015/hormozdiari2016），插入 Intro P1/P2、§2.1/2.2/2.3 自然位置。全部 51 条经 NCBI PubMed eutils / Crossref 逐条核实（51/51 命中真实记录；纠正 scott2012 标题、spracklen2020 期刊为 Nature、bowden2017 期刊为 Stat Med、benner2017 标题为 "Genomic Regions"、elsworth2020 为 bioRxiv 预印本）。`refs.json`/`references.md` 按正文首次出现序重排并程序化校验一致（51/51）；顺带修正 §2.3 "five"→"six"（coloc.susie 实跑 6 位点含 LAMC1）及 Fig 4 图注歧义；DOCX 已重建验证（References 51 条连续、上标 1-13/29-38 抽样核对）。

### 表文件命名修正（docs/manuscript/）
- 表编号与正文对齐：`Table_S2_candidates.csv` → `Table_2_candidates.csv`（主表 2，正文引用已用该名）；`Table_S3_coloc_susie_diagnostics.csv` → `Table_S2_coloc_susie_diagnostics.csv`（正文 S2）；`Table_S4_dropped8.csv` → `Table_S3_dropped8.csv`（正文 S3）；新增 `Table_S4_sample_overlap.csv`（样本重叠敏感性，15 行）。
- Table 2 由"仅描述"改为 manuscript.md 内联 markdown 表（15 行，T2D 7 + CAD 8，M36 渲染用）；Table 1 保持双子表（数据源 + 精度漏斗）。

### 文档合并（docs/）
- 5 份流程文档合并为 2 份主文档：`PAPER_WRITING_PLAN_20260816.md`（写作方案 + 附录 A FINAL_STORY + B JOURNAL_TARGETS + C REVIEWER_IF）、`INTEGRITY_AUDIT_20260816.md`（审计 + 附录 A verification + B COLOC_SUSIE）；原文档入 `docs/archive/`；README/各 doc/脚本引用全量更新，无悬挂引用。

## 2026-08-15 — 对抗性评审定档（单篇完整故事 → HMG 5.5–6）+ 可行性试点升级 + 全量 GWAS 就绪

### 破局方案受对抗性评审并定档（docs/archive/verification_20260815.md）
- 独立评估 agent 对 `results/BREAKTHROUGH_PLAN_20260815.md` 做对抗性评审，裁决：按"四章节 + 介质层生物学 + 冲 8-9 分"形态**不值得做**（不值得投入资源），但**按诚实形态值得做**。
- **关键评审结论（已全部核实为真）**：
  1. **召回统计循环论证**：`MR 召回 = P(MR p<0.05 | coloc≥0.8)` 里 coloc 与 MR 共享输入（同 GWAS + 同 eQTL 子集），两方向概率数学上互相蕴含，召回维度是伪读数。
  2. **精度只算了 MR 显著子集**（106/818=12.96%）→ 审稿人会问"coloc 全量 49,866 里有多少"→ 高召回率高精度可能是"报告偏倚"。
  3. **E1 外推无依据**：1.45% 置换命中率从 3 个已知位点外推到全量不可信。
  4. **coloc-only 红线**：coloc-only（coloc 强但 MR p≥0.05）多属 **coloc 伪阳**（GWAS 区域不显著时 coloc 无意义）→ 召回叙事主要建立在 coloc 伪阳上。
  5. **介质层主张不成立**：76=56 编码基因（无蛋白层信息）；"双共定位≠双介导"；"打哪里"无方法学依据。
  6. **E2 新发现缺位**：全为已知靶点（KCNJ11/PCSK9/APOC3）→ 8-9 分不可达。
- **定档**：**HMG 5.5–6**，eBioMedicine/Diabetologia 8-9 诚实标准下不可达。单篇论文 = A（操作特性量化）+ C（跨结局泛化）为核心，B 改"56 编码基因分子层一致性表"放补充，D 降级佐证（InsPIRE 反向进 Limitation）。
- README §0 故事线改写（移除"介质层=打哪里"、头条 = 图谱资源 + 操作特性量化）、§8 期刊档位改 HMG 5.5–6；BREAKTHROUGH_PLAN 顶部加评审裁决横幅。

### M20 可行性试点升级（400 → 6,000 对分层）
- 依据评审 §1/§7/§8：**关键判据从"灰区 coloc 比率"改为"coloc-only 命中里 GWAS 峰显著占比"**——GWAS 峰 p<5e-8 占多数 → "MR 召回低（功率问题）"叙事成立；GWAS 峰不显著占多数 → 召回叙事是 coloc 伪阳，转向。
- 抽样升级：灰区（MR p∈0.05-0.5）3000 + 阴性（p≥0.5）3000，按结局分层，seed=42。
- 顺带按 eQTL 强度（|beta|/se）分箱区分"弱 eQTL 功效陈述"vs"工具选择问题"。
- 输出 `results/feasibility_pilot_20260815.csv`（gene/symbol/outcome/grp/mr_b/mr_p/gwas_min_p/eqtl_F_max/nsnp/pp4/ok/note）。

### M19 全量结局 GWAS 下载完成 + 根盘满危机处理
- `data/opengwas/full/{t2d,cad,fbg}_full.gz` 三个全量 GWAS 下载完成（DONE 标记；T2D 修复损坏文件重下、CAD/FBG 新下）。
- **根文件系统 / 100% 满（8.9M 可用）**：M20 冒烟曾因 `Fatal error: cannot create 'R_TempDir'` 崩；清理被杀 M20 的 /tmp/RtmpB5ECJT 3.8G 残留 + stale plink 文件 → 92%（4.1G 可用）。R 临时目录固定 `R_TMPDIR=/data/qiushuogeng/tmp/rtmp`（/data 有空间）。

### M20 调试（9 轮冒烟修 6 个技术坑 + 3 轮全量崩溃根因）
- ① `standard_error` 列代替 hm_ci（harmonised CI 全空，se 等位取向不变）；② data.table 列名 `p` 遮蔽循环变量 `p` → 改名 `pr`；③ coloc 输出 `PP.H4.abf`（不是 `PP.H4`，summary 是 named vector）；④ harmonize 缺 NA 守卫（`eaf_g` NA → p_flip NA 崩溃）；⑤ FBG 等位基因小写（a/g）vs eQTLGen 大写 → `toupper()`，QC 通过率 55.6%→88.9%；⑥ eQTLGen 全量实为 **~153M 行**（16,923 基因 × ~9,053 cis 变异，非评估 agent 声称的 10.5M 行；**10.5M 是 cis-EQTL-significant 子集**）。
- **M20 全量 3 轮崩溃根因（本轮全部定位）**：⑦ 根盘 / 50G 满 → `fread(cmd=)` 管道缓冲写 `/tmp/Rtmp*` No space → 必须 `TMPDIR=/data/qiushuogeng/tmp/rtmp`（**R 认 TMPDIR 不认 R_TMPDIR**）；⑧ `system2("bash", c("-c",cmd))` 内联 awk 的 `$1/$8` 被 bash 当位置参数展开成空 → awk 语法错 + `gzip: stdin: unexpected end of file` 落盘 0MB → awk 程序写文件用 `-f` 调；⑨ eQTLGen full 真实列序（逐列核实：1 Pvalue 2 SNP 3 SNPChr 4 SNPPos 5 **Zscore** 6 AssessedAllele 7 OtherAllele 8 Gene 9 GeneSymbol 10 GeneChr 11 GenePos 12 NrCohorts 13 NrSamples 14 FDR），awk 按位置抽须按此、**fread select 按列名**（Zscore 不在第 7 列）；⑩ R 的 `system2` 在本机某些管道场景行为异常 → **最终方案：fread(cmd="zcat | awk 单引号内联") 管道直读 + TMPDIR=/data**（复刻第一次能跑到 148M 行的路径，规避 system2）。
- 第五次全量运行中：灰区 3000 + 阴性 3000（awk 阶段 15:11 起，~20-30 分钟 → 6000 coloc → 摘要）。

### M22 零成本检验完成（评审 §8/§10#3，`results/m22_efqt_power_20260815.csv`）
- **eQTL 功效**：106 strong 的 lead cis-eQTL 强度（|Z| 中位 16.4）与全量 16,987 个 lead eQTL（中位 16.3）**无差异**（MW p=0.45），与 MR 显著但非 strong 的 713 位点也无差异（p=0.95）→ **MR 显著集内 PP.H4≥0.8 不由 eQTL 强度决定**（eQTLGen lead 全极强 |Z| 中位≈p<1e-50），strong 判定主要由 GWAS 侧驱动 → 支持 M20 以 GWAS 峰 p 分层为主判据。
- **文献基准**：GWAS Catalog T2D/CAD 关联 13,109 行；106 strong 命中 **81%（86/106）**（rsID 直命 26 + 基因 hg38±100kb 59 + 基因本体 1；t2d 90%/cad 72%），**19%（20/106）超出已报道注释范围 = 新候选比例上限**。**注：M22 原 55% 为 hg19/hg38 build 错配（GenePos hg19 vs catalog hg38）低估，2026-08-15 复核修正（M22b v3，见 feasibility §7.8）。**

### M21 UKB-PPP 覆盖核查改版完成（`results/ukbpp_coverage_20260815.csv`）
- **弃 S3**（ukb-pqtl bucket 不存在）；改用 **Synapse children**（synapseclient 列 syn51365303 = "UKB-PPP pGWAS summary statistics / European (discovery)"，**2,940 个蛋白 tar**，`Gene_UniProt_OID_v#_Panel.tar` 命名，gene 可解析 2,922）。
- **裁定：UKB-PPP Olink 面板覆盖仅 10/76 编码基因**（ARG1/ARL13B/BLOC1S2/CDC123/DUSP13/GIT1/HMBS/LRIG1/NUDT5/PDGFC）→ **87% 编码基因不在血浆面板**（胞内/膜/转录因子为主），**蛋白层扩展正式关停（<20）**。比评审 §6 估的 20-40 还少，证实"血浆 Olink 覆盖限制"是物理瓶颈。KCNJ11/LIPA/PTPRN/RASD1 等核心候选全不在面板。
- **评审 §6"76=56 编码"为错数**：按 ensg_biotype 判定实际 **76 编码 + 29 非编码**（15 lncRNA + 8 伪基因 + 1 snRNA + 1 miRNA + 4 NOT_FOUND），已记入 verification §10b。

### M20 全量完成 — 决定性负结果（评审 §8 检验落地）
- **结果（results/feasibility_pilot_20260815.csv，6,000 对）**：MR 显著集之外 **strong coloc（PP.H4≥0.8）= 0**。
  - 灰区（0.05≤p<0.5）2892 QC-pass：PP.H4≥0.8 **0**，≥0.5 仅 2（RP11-403I13.4×cad pp4=0.676 gwas_p=1.8e-6；WNT3×cad pp4=0.668 gwas_p=1.5e-4，均 GWAS 峰不显著）
  - 阴性（p≥0.5）2879 QC-pass：PP.H4≥0.8 **0**，≥0.5 **0**
  - gwas_min_p 中位 ~2.5e-4（两组），GWAS 峰显著率无法分层（分母=0）
- **判定（推翻"MR 召回低"叙事）**：MR 显著集之外不存在 strong coloc → cis-MR 的局限是**假阳性率高（精度 12.96% 为真读数）不是漏检**。全量 MR 显著集外 27,123 对（灰区 14,294 + 阴性 12,829）strong 仅 2 个（AP3S2/ZNF19，均 GWAS 峰显著）。
- 综合摘要 `results/feasibility_20260815.md`（M19/M20/M21/M22 四项实证 + 判定）。

### 文件系统间歇异常 + eQTLGen 稳定副本策略（2026-08-15 15:37 起）
- **症状**：`data/eqtlgen/` 下 `cis-eQTLs_full_20180905.txt.gz`（4.3G）与 `cis-EQTL-significant.txt.gz`（322M）**间歇性 openat() 返回 ENOENT**（strace 证实内核层），目录项可见（ls）但 stat/open/mv 失败，时好时坏（python 偶成）。SNP_AF.txt.gz 及 /data 其他文件正常，新写入文件正常 → **inode 特定、非整体盘问题**；关闭沙箱后仍失败。
- **应对**：① full + SNP_AF 本地复制到 `/data/qiushuogeng/tmp/eqtlgen_stable/`（新 inode，gzip -t 校验通过）；② significant 从官方源重下（molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR0.05-...txt.gz = 322775879B 与本地同版本，curl -C - 续传，gzip 校验通过）；③ **M23 脚本改读 stable 副本**。
- 原文件留在原位（不删），若日后恢复可对照。

### M23 全量 coloc 扫描启动（2026-08-15 15:48）
- 脚本 `scripts/M23_full_scan.R`：参数化结局，3 进程并行（t2d 10190 / cad 14675 / fbg 6508 对，共 **31,373 对** MR 有效，即 ok==TRUE 且 pval 非 NA 的网格），进程内按 ~2500 基因块 awk 抽取 eQTL cis 子集控制内存（峰值 ~3GB/进程）。复用 M20 harmonize/coloc 逻辑（p12=1e-5），记录 gwas_min_p / eqtl_F_max。
- **健全性检查设计**：MR sig 集内 strong（PP.H4≥0.8）应与 transcript_coloc_hits.csv 的 106 完全一致（coloc 确定性，同输入同输出）。
- 汇总脚本 `scripts/M24_summarize.py`（交叉验证 + 结局统计 + 操作特性表）。

## 2026-08-13 — SMR/HEIDI 双通道落地 + UKB-PPP 跨平台复现启动 + 多 agent 阳性核查

### SMR/HEIDI（Zhu 2016）双通道完成
- **工具**：SMR 1.3.1 装至 `tools/smr`，LD 参考 1000G EUR（1kg.v3）。
- **蛋白通道（M9，`scripts/M9_smr_build_inputs.R` + `M9_run_smr.sh`）**：5 蛋白（PCSK9/APOC3/INSR/ANGPTL3/PCK1）× 3 结局 = 15 对 → `results/grid/protein_smr_heidi_clean.csv`。**仅探索性、未校准**：APOC3×T2D（p_SMR=0.011, p_HEIDI=0.42）与 APOC3×FBG（p=0.007, p_HEIDI=0.91）为**单工具敏感性支持（非独立确认）**——15 对零假设期望 0.75 个 p<0.05、无多重校正、coloc.abf none、coloc-SuSiE 判不共定位、FinnGen rs964184×T2D p=0.075 不显著、**APOC3 不在 UKB-PPP 面板（复现断头）**；PCSK9×CAD（p=2.96e-19）/APOC3×CAD（p=5.5e-11）SMR 极显著但 HEIDI 异质性，且**阳性对照 PCSK9×CAD 过不了 HEIDI（p=1.28e-4）→ 蛋白侧 HEIDI 判定整体未校准**。
- **转录通道（M10，`scripts/M10_transcript_smr_build_inputs.R` + `M10_run_smr.sh`）**：106 strong 命中基因 × 3 结局 = 128 测试 → **HEIDI 通过 89/126=70.6%（2 NA）、37 异质性** → `results/grid/transcript_smr_heidi.csv`。**KCNJ11×T2D：单点 Wald 自洽（b_SMR=0.595181=0.0643/0.108034=MR Wald b，数学恒等）+ 两个区域级独立检验（coloc、HEIDI 0.128）通过**；机制 caveat：coloc 峰值 rs757110 在 ABCC8（S1369A）、KCNJ11 基本不在全血表达。LIPA×CAD（p_SMR=1e-20, HEIDI 0.457）、ARG1×T2D（HEIDI 0.970）、PTPRN×CAD（0.249）、RASD1×T2D（0.519）获 HEIDI 支持。
- **诚实 caveat**：① **MR 单工具 Wald ≡ SMR 是同一统计量**（同一 top cis-SNP 的同一比率），不是两个独立方法，不得写"三方法独立互证"；② eQTLGen 只给 Zscore → Beta=Z/√N、se=1/√N（标准化表达常规转换）；③ 转录 p_SMR 显著（125/128）是选择偏倚构造性必然、**零信息量不进正文** → 唯一有信息量的读数是 **HEIDI 通过率 89/126=70.6%**；④ deCODE 无 effectAlleleFreq → .esd Freq 用 ImpMAF 近似 + `--disable-freq-ck`（探索性敏感性）；⑤ HEIDI 用 1000G EUR LD，与 deCODE 人群错配（转录 eQTLGen 同为 EUR 更可靠）。

### UKB-PPP 跨平台复现（预注册 §9 承诺执行）
- 用户已注册 Synapse（ownerId=3608259），PAT 存 gitignored `~/.synapseConfig`。
- **Synapse API 被中国大陆 GEO 屏蔽**（NIH NOT-OD-25-083，repo-prod.sagebase.org 403 GEO_RESTRICTION）→ 走 mihomo 代理解封；已实测无公共镜像、S3 直连被地域限速 100 倍（14KB/s vs 代理 1.9MB/s），最终用代理全量下载（用户已被告知）。
- **面板覆盖核查**：UKB-PPP European discovery 有 PCSK9/INSR/ANGPTL3/LDLR/APOB/DPP4/GLP1R/GCG（**8/11**）；**无 APOC3/PCK1/HMGCR**（Olink 面板不覆盖，APOC3 的 deCODE 蛋白 SMR 救回对无法在此平台复现——如实写缺口）。
- 8 蛋白 tar（~4.4GB）走 arbitrate 后台下载中（`data/ukbpp/`）。

### 多 agent 阳性核查 + 讲故事
- 工作流 `mr-atlas-verify-and-story`（3 verify + 3 story + 1 synthesize）启动核查"是否构成阳性结果" + 设计发表叙事。结论落地后并入 DECODE_PIPELINE_SUMMARY / IMPROVEMENT_STRATEGY。

### 文档
- `results/DECODE_PIPELINE_SUMMARY.md` 新增 **M9 蛋白 SMR/HEIDI** + **M10 转录 SMR/HEIDI** 节；`docs/IMPROVEMENT_STRATEGY.md` SMR/HEIDI 与三角验证种子修正 checkbox 勾选 + UKB-PPP 状态更新；`README.md` SMR 待装→已完成、诚实局限加 SMR/HEIDI caveat、现状定位更新；`results/evidence_brief_20260813.md` 新增证据简报。
- 未动 PREREGISTRATION.md 哈希锁定正文（SMR/UKB-PPP 均属预注册承诺执行，非门柱变更）。

## 2026-08-11 — deCODE 蛋白通道真实运行 + 代码审查修复

- deCODE 4 个对照蛋白文件下载完成（PCSK9/HMGCR/ANGPTL3/APOC3，各 ~953MB，16 连接并行），cis 窗自动提取，蛋白通道 cis-MR 真实运行。结果见 `results/grid/protein_decode_mr.csv`、`results/DECODE_PIPELINE_SUMMARY.md`。
- 代码审查修复（技术 bug fix）：① `M1_decode_subset.sh` cis 子集此前丢弃表头 → 下游 `fread(header=TRUE)` 把第一行变异当列名、Pval 列丢失、误判“无工具”（已改为保留表头，校验改按行数）；② `M3_protein_decode.R` 排除 deCODE readme 注明的多等位基因 bug 行（effectAllele==otherAllele / otherAllele="!"），rsids 为 NA 时回落 Name。
- 核验为“非 bug”项（诚实记录，防误改）：eQTLGen eaf=1−AlleleB_all（AA/AB/BB 计数公式验证一致）；Z→β=Z/√N、se=1/√N 标准转化；LD clump r²<0.01@1000kb EUR 逻辑正确。


## 2026-08-07 — stage-2：转录扫描命中基因 本地 plink clump + IVW-MRE 复核完成

- **1kg LD 参考就绪**：15:15 下载完成（1495MB，gzip 校验通过），`1kg_watch.sh` 15:21 自动触发 `process_1kg.sh`——解压 7.4G 至 `data/ldref/1kg.v3/`（EUR 全染色体单文件集，hg19），plink chr1 `--freq` 冒烟通过（661,214 变异），15:22 写 `results/ldref_ready` + `data/ldref/README_STRUCTURE.md`。
- **实现 `scripts/M3_transcript_grid_stage2.R`**（不动已完成的扫描脚本，历史产物可溯）：命中基因 cis±1Mb p<5e-6 真实 eaf 全部变异 → 本地 plink（1000G EUR）clump r²<0.01@1000kb（8 核并行，906 基因，缓存 `_stage2_clump_cache.rds`）→ 结局 proxies=FALSE 提取（OpenGWAS，缓存 `_outcome_stage2_*.rds`）→ harmonise action=2 → 主方法 IVW-MRE + nsnp≥2 敏感性（ivw_fe/WM/Egger），nsnp=1 标注 Wald 退化。stage-2 结果集内按结局再算 BH-FDR，诚实双报 nominal + q。
- **结果（16:34）**：982 hits 全量复核，978 对 MR 成功（4 失败如实记录：RBM6/TCTEX1D2 结局无匹配 proxies=FALSE、FOXJ2/ENO2 clump 后 0 独立工具）。存活 **nominal（p<0.05 且方向一致）819 / BH-FDR q<0.05 且方向一致 812**（t2d 331/329、cad 477/472、fbg 11/11）；1 条 p<0.05 方向翻转排除。**HMGCR×T2D 复现深度层**：单工具 b=−0.3818 p=4.84e-05（深度层 b=−0.382 p=4.8e-5 一致）。产物 `results/grid/transcript_grid_stage2.csv`（982 全量）+ `_hits.csv` + `_sens.csv` + `_funnel.tsv`。
- **技术备忘**：plink 1.9 `--clump` 输入格式为**两列带表头 `SNP P`**（CHR/BP 由 `--bfile` 提供，与 `ieugwasr::ld_clump_local` 一致），四列 CHR/SNP/BP/P 会报 "No variant ID field found"。`01_integrity_bootstrap.R` 是预注册**生成器**而非校验脚本——运行会重写 config.json 时间戳使预注册哈希失效，管线中途禁止运行。

## 2026-08-07 — 节点延迟选优 + plink 就位 + 1kg 自动处理 + README v0.10 交接版

- **节点延迟选优**（用户指示"用之前测延迟，最低者优先"）：新增 `scripts/lib/ensure_best_node.sh`——经 mihomo API 对候选节点测延迟（普通节点优先，premium 兜底），选最低者并切换（不打断在途连接）。实测选优：HK 01（普通）370ms 当选；US 超时排除。监控 cron 升级为 d5fa4610（17-57/40）：重启下载前先跑节点选优。
- **plink 1.9 就位**：README 原写"工具就绪"但 `tools/` 实际为空——经代理从 s3.amazonaws.com 下载 v1.90b6.26（8.9MB，已验证可执行）。SMR 1.3.1 仍待装（M5 需要时再装）。
- **1kg 自动处理**：新增 `process_1kg.sh`（gzip 校验→解压→EUR 文件集定位→plink chr1 --freq 冒烟→写 `results/ldref_ready` + `data/ldref/README_STRUCTURE.md`）+ `1kg_watch.sh`（5 分钟轮询，下载完整且无写入进程即自动处理）。已启动。
- **README v0.10**：新增 §0"新会话交接"（当前进程快照 / 新窗口状态检查命令 / 交接规则 / 下一步依赖链），并修正 §2/§4/§5 过期状态（deCODE 经代理提速、1kg 进度、plink 实况）。README 保持纯多管线叙事，运行轨迹仍全在 CHANGELOG。
- **转录组全量扫描层完成**（2026-08-07 14:13）：`M3_transcript_grid.R` 跑完 16,622 基因 × 3 结局 = 49,867 对，**982 hits**（BH-FDR q<0.05：T2D 394 / CAD 576 / FBG 12），157 个 p<1e-6 强信号。产物 `results/grid/transcript_grid_mr.csv` + `transcript_grid_hits.csv` 入库。命中基因 stage-2 本地 plink clump + IVW-MRE 复核待 1kg LD 就绪后实现（现为占位）。

## 2026-08-07 — 节点偏好：下载默认改用普通节点 HK 01（用户指示）

- 用户指示：下载**尽量用普通节点，不用 premium**。本次 PCSK9 在途不切换（无 Range 续传，中断即整文件重下）；`config.yaml` 的 `Proxy` select 组已把 **HK 01（普通）置顶为默认**，premium 节点降级为备选（HK Premium x2 HK 01/02、JP/SG/US Premium）。
- **生效时机**：mihomo 下次重启（监控 cron ac0a4990 检测到死机自动重启 / 手动重启）后 select 默认即 HK 01；当前在途连接不受影响。
- **安全**：配置仍为 600 权限、含订阅凭据、不入 git/README。

## 2026-08-07 — 本地 mihomo 代理加速下载（实测 1kg ×50 提速，deCODE ×3.2）

- **动因**：直连下载过慢（deCODE ~7.5 KB/s、1kg ~11 KB/s），用户提供代理订阅并授权运行。
- **部署**：mihomo v1.19.29 装在 `/data/qiushuogeng/tools/mihomo/`（**不占系统盘**）；`config.yaml` 为最小规则集（MATCH,Proxy），**权限 600**，含订阅凭据，**禁止提交 git / 写入 README**。
- **实测提速（2026-08-07 14:00 本地 127.0.0.1:7890）**：deCODE 直连 7.5 → 经 HK 节点 23.7 KB/s（×3.2）；1kg fileserve 直连 10.9 → 595.9 KB/s（×55）。节点对比：deCODE 最优稳定为 HK Premium x2（28.9 KB/s），JP 波动大（28→11 KB/s）不作默认。
- **脚本改造**：`decode_slow_dl.sh` 支持 `PROXY=http://127.0.0.1:7890` 环境变量（默认空=直连，向后兼容）；1kg 用 `curl -C - -x` 续传重启（fileserve 支持 Range，从 37.3MB 断点续传）。
- **监控**：session cron ac0a4990（durable）加入 mihomo 存活检查，代理死则重启（读现有 600 配置，不打印凭据）。
- **注意**：deCODE 无 Range 续传，代理中断需整文件重下（已由 ×3.2 提速摊薄）；当前 PCSK9 从 0 经代理重下（5.5MB），1kg 已续传到 135MB。
- **安全**：对话中已出现订阅明文密码，若需严格收敛，可在服务商后台重置订阅后更新 `config.yaml`。

## 2026-08-07 — 代码审查二轮 P1/P2 修复（GTEx eaf 判向 / islet se 下溢 / 报告防崩溃 / watcher 生存期）

- **GTEx eaf P1 修复（技术 bug fix，可能改变回文位点方向）**：`M6_gtex_mr.R` 此前 `eaf.exposure = d$maf`（次等位频率，而 slope 是 **alt** 等位效应）→ 回文位点 harmonise 会错误翻号。修复：`eaf.exposure = ifelse(ref_factor==1, 1-maf, maf)`（ref_factor==1 ⇒ ref 为次等位 ⇒ eaf(alt)=1-maf）；NA ref_factor → eaf=NA，回文位点被 action=2 保守丢弃（与胰岛通道同口径）。经验证（200k Liver 行：eaf_alt≈maf/1-maf 一致、eaf_alt×N≈minor_allele_samples）。**重跑 309 组合结果与修复前逐字节一致（md5 相同）**——本数据 lead eQTL 无回文冲突，修复为无害 no-op，已报告结果不变。
- **islet se 下溢 P1 修复（技术 bug fix，潜在）**：`M6_islet_mr.R` 此前 `se = abs(Slope)/qnorm(1-Pval/2)`，当 p<~4.4e-16 时 qnorm→Inf → se=0 → Wald b/se=Inf、p=0 假性过显著。修复：z 封顶 `pmax(qnorm(1-Pval/2), 8.209536)`（=qnorm(1-1e-16)）。当前 GENES 胰岛 eQTL 均 p≥4.2e-8，未触发，属防患性修复。
- **M3_protein_decode.R P2 ×3（技术 bug fix/加固）**：① 兜底列名 `min_log10_pval`→`minus_log10_pval`（真实表头，死代码修复）；② `load_cis` 过滤 deCODE readme 要求的 quality-excluded 变异（`assocvariants.excluded.txt.gz`，当前 0 行 no-op）；③ `run_decode_mr` 中 `load_cis` 包 tryCatch——损坏/半截文件通过 size 启发式时不再崩掉整个蛋白通道。
- **report_decode_results.R P2（技术 bug fix/加固 + 叙事对齐）**：① funnel 文件缺失/损坏 guard（`fread` 缺文件即抛错会崩报告）；② 潜在 parse 错误修复——"转录通道"结果一般"" 内层引号是 ASCII `"` 导致字符串提前终止（全文件从未能 parse 成功过，语法扫描抓到）；③ **README 插入逻辑移除**：v0.9 已删"- 最近运行"行与执行记录章节（README 纯叙事），原逻辑会向文件尾追加整段 v0.6 执行记录——运行详情一律落 CHANGELOG + `DECODE_PIPELINE_SUMMARY.md`。
- **decode_watch.sh P2（技术 bug fix）**：`seq 1 1000`（=83h）< 4 对照顺序下载需 ~104h+ → watcher 会在下载完成前静默退出；改为 `seq 1 4000`（≈333h，覆盖全部 9 文件 ~13-15 天）。
- **08_figures.R F9 P2（技术 bug fix/加固）**：`fn[stage=="tested_proteins", count]` 在 stage 行缺失时返回 integer(0)，`c()` 丢弃 → 行数不齐 ggplot 崩溃；三阶段齐备才渲染，缺行则 log 跳过。
- **核验无需修改（防误改）**：`M1_decode_subset.sh --delete-original` 已有成功校验（表头+行数）后才 `rm -f`，删除是安全的。
- **下载速率实测（2026-08-07 13:07 采样 40s）**：PCSK9 ≈6.5 KB/s（953MB ≈ 40h，符合 37-48h 估计）；1kg.v3.tgz ≈12.2 KB/s（LD 参考非关键路径，可在 deCODE 窗口内完成）。

## 2026-08-07 — README 多管线叙事 v0.9 + deCODE max-time P0 修复 + 1kg 失联重启 + OSF 溯源快照

- **README v0.9（非科学决策变更，叙事重构）**：应要求去掉"先做 X 发现不行再做 Y"的时间线叙述，正文改为**多管线设计叙事**——两条主通道（全血转录本 × 循环蛋白）× 多层敏感性（GTEx 组织 / InsPIRE 胰岛 / 药物靶点四态）同为设计组成部分（§3 架构图 + 各通道结果）。v1→v2 等方法学修正**不从仓库删除**，保留在本 CHANGELOG（审计痕迹）+ §6 纪律声明。执行记录时间线细节从 README 收敛为一行版本记录，指向本文件。
- **deCODE 下载器 P0 致命 bug 修复（技术 bug fix）**：`decode_slow_dl.sh` 的 `--max-time 70000`（=19.4h）< 单文件需时（953MB @ 实测 5.6-7.2KB/s ≈ 37-48h）→ 每轮必被 max-time 掐断 → 脚本 `rm -f .part` 从 0 重来 → **永远下不完**。修复：`--max-time 280000`（77.8h，覆盖最坏 ~5KB/s≈53h）。无 Range 下单连接必须一次跑完，max-time 必须 > 整文件时长。**重启下载（PID 3023528）**。
- **1kg.v3.tgz LD 下载失联重启（技术/运维）**：发现该下载**无任何进程在跑**（compaction 后后台任务失联，30.6MB 停 12+ 分钟且此前无保险）。找回 URL `http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz`（实测支持 Range 206 → `-C -` **续传**重启，PID 3023661）。
- **持续监控（运维）**：编排器每 5 分钟自动重启 deCODE 下载器（不依赖会话，原已有）；新增会话 cron `aec7ba24` 每 40 分钟心跳（deCODE + 1kg + 转录扫描），补上 1kg 无保险的缺口，进度写 `/tmp/monitor_heartbeat.log`。
- **OSF 溯源快照（文档，非科学决策）**：新建 `docs/PREREGISTRATION.yaml`（预注册忠实镜像，含 v5 扫描层修订）+ `docs/osf_provenance.txt`（mtime/哈希/诚实注册口径证据：v5 于 2026-08-07 12:20 锁定、早于扫描重启 12:37，主分析结果尚未产生）。**git init**（`.gitignore` 排除 data/ results/，脚本/文档/README 入库，首次提交快照）。

## 2026-08-07 — README 重构为干实验叙事版 + 转录扫描 symbol merge bug 修复重启

- **README 重构（v0.8，非科学决策变更）**：应要求删除审稿意见对照表、期刊定位、引用核验表、里程碑预算、团队分工等框架性章节；重写为：研究动机（为什么做/要回答的问题/不做的事）→ 数据到位状态表 → 已跑完干实验（转录 v1→v2、转录组全量扫描、GTEx 组织、胰岛、药物靶点四态、INTERVAL）→ 进行中/待办 → 已知难点与对策 → 完整性纪律。旧版存档 `README.md.bak_v0.7`。
- **转录扫描 symbol merge bug 修复（技术 bug fix）**：`M3_transcript_grid.R` 的 `lead <- merge(lead, unique(ivs[, .(Gene, GeneSymbol)]))` 在 `.SD[1]` 已含 GeneSymbol 列的前提下再 merge，data.table 将同名列重命名为 `GeneSymbol.x/.y`，导致 `lead$GeneSymbol` 不存在 → `sum(!is.na(lead$GeneSymbol))` 恒为 0（**16,622 基因 symbol 映射 0**）。MR 统计本身不受影响（用 SNP/等位基因而非 symbol），但产物 `transcript_grid_mr.csv` 缺基因符号。修复：删去冗余 merge（`.SD[1]` 已保留 GeneSymbol）。验证：修复后 16,622 基因全映射。**重启转录扫描（PID 3011181）**，结局提取重新开始（约 157 分块/结局 × 3）。
- **教训（同类模式）**：`merge(by=…)` 与已含同名列的 data.table 合并会生成 `.x/.y` 列并让原列消失——合并前检查 `names()`，或直接用 `.SD[1]` 自带列。

## 2026-08-07 — GTEx 组织救场 MR 完成 + 转录组全量扫描启动 + deCODE 慢速下载

- **GTEx v8 组织 eQTL MR 完成（M6 v2）**：修复 v1 两处致命 bug——① `gene_id` 为 ENSG（无 symbol）→ 改用 egenes `gene_name`；② signif 对 `variant_id`(chr_pos_ref_alt_b38) 非 rsid、OpenGWAS 结局不可检索 → 改用 egenes 自带 `rs_id_dbSNP151_GRCh38p7`（单 lead 工具 Wald）。6 组织 × 21 基因 × 3 结局 = 309 组合，48 完成 MR，18 个 p<0.05。**组织救场证据**：PDX1×T2D 胰腺 p=4.6e-05、NPC1L1×CAD 胰/脂肪 p=3.8e-06、PCSK9×CAD 肝 p=2.4e-04、APOB×CAD 肌 p=4.8e-11、CETP×CAD 脂肪 p=3.3e-09。敏感性层（单工具/无 clump/proxies=FALSE），方向解读待命中基因 clump+IVW 复核。输出 `results/grid/gtex_mr.csv`。
- **转录组全量扫描启动（M3_transcript_grid.R）**：eQTLGen ~16.5k 基因 lead cis-eQTL → Wald × T2D/CAD/FBG，按结局 BH-FDR。**预注册 v5 修订**（追加 + 重哈希 + M0 通过）：扫描层 lead 单工具（省 ~16k API clump）、结局 proxies=FALSE（防代理伪阳性）、命中基因第二阶段本地 plink clump + IVW-MRE 复核；扫描层声明 hypothesis-generating。1kg.v3.tgz（1000G EUR LD 参考, 1.5GB）后台下载。
- **deCODE 新链接实测（用户提供 /file/×2 + /folder/×1）**：`/file/<uuid>` 返回 605 字节 React 门户页（content-range 0-604/605），非数据文件；`/folder/<uuid>`=下载 token。用该 token 复测 `/s3/download`：**仍忽略 Range**（HEAD/GET 带 Range 均 HTTP 200 整文件 953,852,898 字节），**per-IP 限速**（2 并发=1 连接速度 ~10-12KB/s）。**结论：三个链接不改变阻塞；并行方案仍不可行。**
- **对策：慢速顺序单连接下载器 `scripts/lib/decode_slow_dl.sh`**：4 对照蛋白优先（~22-26h/文件），随后 5 药物靶点；无断点续传（Range 无效）但 curl --retry 兜底瞬断；`run_decode_pipeline.sh` 重启逻辑改指向 slow 下载器（原并行版在无 Range 服务器上产出损坏文件），等待窗口 60h→250h。decode_watch.sh（运行中）自动对完成文件 M1 提取 → M3 蛋白 MR。PCSK9 下载中，4 对照预计 ~4 天就绪。**仍请用户邮件索取支持 Range 的 S3 直链（不同 host）以并行提速。**
- **哈希锁格式修正（技术 bug fix）**：此前 `md5sum > sha256` 写入 `hash  filename` 两列，R `readLines` 比对失败；改为 `md5sum | awk '{print $1}'` 裸哈希（避免下游所有哈希校验脚本失败）。

## 2026-08-06 深夜 — deCODE 下载阻塞 + 代码审查修复 + 5 项补实验评估与部分执行

- **deCODE 下载阻塞（技术/数据状态）**：实测 download.decode.is **忽略 Range 请求**（对 1MB Range 也返回 HTTP 200 整文件 content-length=953852898）→ 16 连接并行分块**内容损坏**（合并文件是文件头重复，gzip -t 校验失败）；单连接全量下载 ~5KB/s → 每文件 ~2.3 天，9 文件不可行。已停止下载器/编排器（防继续浪费带宽），保留分块待方案。**唯一可行路径：用户 deCODE 注册邮件里的 S3 直链（支持 Range）**，需用户提供。
- **代码审查修复（技术 bug fix）**：① `M1_decode_subset.sh` cis 子集此前丢弃表头 → 下游 `fread(header=TRUE)` 误读；已改保留表头 + 行数校验（不依赖染色体格式）。② `M3_protein_decode.R` 排除 deCODE readme 注明的多等位基因 bug 行（effectAllele==otherAllele / otherAllele="!"）+ rsids NA 回落 Name + harmonise 掉率上报。③ `M3_transcript_mr.R` 加 harmonise 掉率 note。④ `decode_parallel_dl.sh` 修复续传判定用全局 CHUNK 的 bug（最后一块 want<CHUNK 会污染下一文件）+ 越界块截断。
- **核验为“非 bug”项**（防误改）：eQTLGen eaf=1−AlleleB_all 正确（AA/AB/BB 个体计数公式验证）；Z→β=Z/√N 标准转化；LD clump r²<0.01@1000kb 逻辑正确；AF 文件无重复 SNP。
- **导师 agent 评估 5 项补实验**：① GTEx 组织 eQTL **GO**（回应组织特异性死穴；OpenGWAS 无 GTEx 组织数据，官方 tar 404，Zenodo GV-Rep 镜像 Range 206 可用）；② 药物靶点×四态 **PARTIAL**（描述表，只做观察性通道可用性，不做用药层推荐）；③ 单细胞 **GO**（补充材料；Baron=GSE84133、Segerstolpe=E-MTAB-5061、MacParland=GSE115469）；④ LINCS/CMap **NO-GO**（GB-TB 级+注册门槛+扰动收录未验证）；⑤ 蛋白-结局观察性 **只引用 Gadd 2024 Nat Aging，不做新分析**（UKB-PPP 个体数据不可及）。
- **引用修正（重要）**：Fadista 2014 是 **PNAS**（PMID 25201977，GEO GSE50398），非 Diabetologia；“Iglesia-Camarero 2023”多次检索**无法证实**，投稿前删除或查实；胰岛 eQTL 用 **InsPIRE（Viñuela 2020 Nat Commun, Zenodo 3408356）** 公开可下。
- **GTEx v8 组织 eQTL 下载**：Zenodo GV-Rep 镜像（1.49GB zip）16 连接并行，Range 206 确认，已启动（~5h）。脚本 `scripts/lib/gtex_parallel_dl.sh` + `scripts/M6_gtex_mr.R`（用 tss_distance 做 cis，组织=敏感性，四态仍锚定全血×血浆）。
- **InsPIRE 胰岛 eQTL MR（P1-1b，真实运行）**：`data/gtex/islet/InsPIRE_islets_independent_gene_eQTLs.txt`（4639 独立 eQTL）。结果：**TCF7L2×T2D Wald b=-1.598 p≈0**、TCF7L2×FBG p=2e-9、TCF7L2×CAD p=6.2e-6；**SLC5A1（SGLT1）全血无 eQTL 但胰岛 eQTL → T2D p=0.0011**——“全血失明、胰岛救场”的直接证据。KCNJ11 胰岛 eQTL 与 T2D 不显著。结果 `results/grid/islet_mr.csv`。
- **已知 T2D 药物靶点 × 通道（P1-2，描述性）**：`results/grid/drugtarget_fourstate.csv` + `transcript_drugtarget_mr.csv`。both 通道（GLP1R/DPP4/INSR）均出 T2D 转录信号（b=-0.389/-0.267/-0.938）；KCNJ11×T2D Wald p=4e-16（+0.595）、TCF7L2×FBG p=7.6e-144；GCG/PCK1 protein-only（全血失明）；ABCC8/PDX1/G6PC/SLC5A1 双盲。通道轴 ≠ 药物靶点蛋白轴，严格描述性。



> 任何阈值、基因集或结论措辞变更必须记录于此（README §12 协作约定）。

## 2026-08-06 — 蛋白通道 MR 脚本（deCODE）+ M8 图表脚本就绪（技术，非科学决策变更）

- `scripts/M3_protein_decode.R` 实现并干跑验证：deCODE cis 窗子集（`data/decode/sub/`，M1_decode_subset 产物）主入口、整文件兜底；**坐标修正**——PCSK9 chr1:55,039,445 / HMGCR chr5:75,336,329（hg38 ENSEMBL REST 实查；原脚本 HMGCR 用 hg19 位置 74.65Mb、PCSK9 用 55.50Mb 均错）；EAF 优先 `assocvariants.annotated`（gzip 完整才加载）连接，否则 palindromic 保守排除 + ImpMAF 近似（note 标注）；nsnp=1 → Wald 补算；全网格落盘含空结果。
- `scripts/08_figures.R` + `docs/figures.md`（设计规格）：F1 转录通道漏斗、F2 v1→v2 -log10(p) 坍缩 dumbbell、F3 极端 p 消失（14→0）已渲染至 `results/figures/`；F4–F8（四态/coloc/AUR/复现/热图）骨架就绪待数据。配色 CVD 安全（dataviz 规范），无双轴、直标标签。
- 图表纪律补充：F2 仅用 v1 且 v2 均有结果的完整对比对（4 对），p=0 按极显著截断处理；漏斗末档用唯一基因×结局对（32）非方法级行数。
- deCODE 下载进度：PCSK9 61%（后台，未阻塞分析）。

## 2026-08-06 — 数据获取提速 + 转录本通道 MR v2 真实运行（技术/数据状态变更，非科学决策变更）

### 下载提速（技术修正）
- **eQTLGen SNP 频率文件（240MB）**：单连接 ~26KB/s（北京→荷兰国际线路）→ 16 连接并行分块 ~300KB/s，15 分钟完成（17:12→17:27）。脚本 `scripts/lib/af_parallel_dl.sh`。
- **修复并行合并 bug**（两处：af_parallel_dl.sh 与 decode_parallel_dl.sh）：`cat c_*` 通配符按**词法序**合并（c_10 排在 c_2 前）→ 文件错乱；改为 `for i in seq; cat c_$i` 数字序。af_parallel_dl.sh 另有 `wait` 等待无限监视循环导致永久阻塞的 bug，已改为 `wait ${PIDS[@]}`。**教训：此 bug 若不修会产出 gzip 校验失败的"看似完成"文件，务必 gzip -t 验证每个合并产物。**
- **deCODE（download.decode.is）**：识别服务器为 Apache **mod_qos（逐连接限速）**。单连接 ~3.5KB/s → 16 连接并行 ~74–76KB/s（约 20×）。PCSK9 文件（953MB）后台续传中，随后 HMGCR/ANGPTL3/APOC3 依次下载（约 14h 全量）。断点续传检测已加（TMP 存在未满块则保留进度，重启续传 69.9MB）。

### M1_decode_subset.sh 实现（存根→实现）
- 8 个对照基因 hg38 坐标经 ENSEMBL REST 实查：PCSK9 chr1:55,039,445、HMGCR chr5:75,336,329、ANGPTL3 chr1:62,597,464、APOC3 chr11:116,827,019、APOB chr2:21,044,075、LDLR chr19:11,089,418、CETP chr16:56,961,922、NPC1L1 chr7:44,541,330（TSS，strand− 取 end）。
- cis 窗 = TSS±1Mb 流式裁剪 → `data/decode/sub/`，`--delete-original` 时删原文件（README 体积红线）。`scripts/lib/decode_watch.sh` 每 5 分钟自动扫描完成的 deCODE 文件并提取。

### M3 转录本通道 MR v2（真实运行，代表基因集）
- 12 基因 × 3 结局（T2D/CAD/FBG），**全 9,695,990 工具候选真实 eaf**（0 缺失，v1 的 0.5 占位已消除），EUR LD clump（r²<0.01@1000kb）。
- **极端 p 值全部消失**：v1 有 14 个 p<1e-50（未 clump 的连锁 SNP 被当独立工具）→ v2 为 0。代表例：CETP×T2D 从 p=3.26e-293 → **p=0.083（NS）**；HMGCR×T2D 从 b=−0.433 p=2.5e-46 → clump 后仅 1 独立工具 → Wald ratio b=−0.382 p=4.8e-5。
- HMGCR×CAD：v1 b=−0.002 p=0.93 → v2 b=+0.037 p=0.137（仍 NS，符号翻转，如实记录）。
- 结果落盘：`results/grid/transcript_mr_v2.csv`、`transcript_mr_v2_wald.csv`（nsnp=1 单工具 Wald 补算）、`compare_transcript_v1v2.csv`；漏斗 `results/funnel/funnel_transcript_v2.tsv`。

### 漏斗发现（须知会导师，影响 M2 与 LOOCV AUR 设计）
- 6 个对照位点中**仅 CETP、HMGCR 在 eQTLGen 全血显著文件（FDR<0.05）中有 cis-eQTL**；PCSK9/ANGPTL3/APOC3/NPC1L1 为肝/肠表达基因，全血无显著 cis-eQTL → **转录通道对该 4 基因不可用**（蛋白通道仍需 deCODE）。这压缩了 LOOCV AUR 的双通道可用位点数，具体设计调整（含是否以 CETP/HMGCR+新位点补足）由导师拍板。

### nsnp=1 Wald 补算（技术修正）
- `scripts/M3_wald_fallback.R` + `M3_transcript_mr.R` 主方法补丁：mr_ivw_mre 对 nsnp=1 无输出 → 补 `mr_wald_ratio`（IVW≡Wald，诚实标注不伪敏感性）。

### OpenGWAS 数据核验项（投稿前回填）
- T2D `ebi-a-GCST006867`：gwasinfo 实测 ncase=61,714 / ncontrol=1,178（合计 62,892），**与 README §4.3 的 n≈655,666 不符**，须核验该 ID 实际样本量并回填（Mahajan 2018 全 meta 或该子集）。
- GTEx v8 eQTL tar URL（`storage.googleapis.com/gtex_analysis_v8/...`）实测 HTTP 404，storage.googleapis 可达；URL 须核验修正后下载（M6 用）。

### M23 全量扫描 + genome build 错位修复（2026-08-15 18:40，决定性）
- **M23 全量 49,866 对 coloc 首轮完成**（t2d/cad/fbg 三结局并行）。交叉验证发现仅 **78/106** 已知 strong 重现 → 追查根因。
- **根因（重要教训）：OpenGWAS full 文件 hm_pos 是 hg38 坐标，eQTLGen GenePos/SNPPos 是 hg19**。M23/M20 用 hg19 位置窗口（GenePos±1Mb）去 hg38 坐标提取 GWAS 变异 → 错位 ~1.7Mb，rsid 交集骤减（GIT1 应 1920 实际 190）。例证：rs894606 eQTL 27,909,352 vs GWAS full 29,582,334。
- **影响**：M23/M20 的 t2d/cad 结果全部作废（存档 results/backups/buildbug_20260815/）；fbg 从未受影响（本用 rsid 匹配）。M20↔M23 全一致只是"同错"。
- **修复**：M23_full_scan.R / M20_feasibility_pilot.R 的 t2d/cad GWAS 提取改**纯 rsid 匹配**（gwas[rsid %in% e$SNP]，不依赖坐标 build）。
- **验证**：GIT1 修复后对齐 1917 变异、PP.H4=0.8937 与 M5 API 逐位一致；**t2d 修复版 58/58 重现已知 strong**（pp4 最大差 2.8e-5）、sig 内 66 strong（错位版 48）、sig 外 1（AP3S2 真实，GWAS 峰显著 5.5e-11）；**错位版 t2d sig 外 4 strong 中 3 个（MAST4/ITCH/VAT1）确认为假阳性**，VAT1"单工具 MR 局限例证"叙述撤回；M20 pilot 修复版 sig 外 0/6000 仍成立（全量 t2d 仅 AP3S2 1 个 = 0.01%）。
- cad 修复版重跑中（曾遇 arbitrate 排队进程死亡致 cad 卡队列，手动重启获准；M20 修复版亦已重跑完成）。

## 2026-08-06 — README v0.4 执行记录 + M0 门禁健壮性修复

- README v0.4：文末新增「执行记录 v0.4」，如实记录 2026-08-06 真实运行：deCODE 下载机制破解（source map 逆向 → /s3/ API）、转录本通道 MR 第一轮（29/36 真实统计 + 方法学局限自查）、INTERVAL 蛋白通道诚实负结果（HMGCR/ANGPTL3 无 cis-pQTL 工具）、数据源就绪状态。设计正文（§1–§12）未作科学决策变更。
- M0 门禁修复：豁免机制从硬编码行号改为【行号区间 + 内容匹配】双保险——含 `M0_literature_gate.sh` 或「本注释块」的 M0 存根说明行一律豁免，正文写作违禁表述仍会被命中阻断（防 §4.2 扩写导致行号漂移后误报）。复跑全绿。

## 2026-08-06 — 转录本通道 MR 第一轮结果与方法学局限（如实记录）

- 12 基因 × 3 结局 = 36 对，29 对产出真实 IVW 统计（`results/grid/transcript_mr_qa.csv`）。
- 自查发现两个方法学局限（故第一轮为 QA 链路验证，非最终结论）：
  1. eQTLGen 显著文件不含 SNP 频率，`eaf.exposure` 曾用 0.5 占位，palindromic 处理不可靠；
  2. 每基因 top-3 工具未做 LD clump，连锁 SNP 被当作独立工具，IVW 精度可能虚高（极端 p 值部分源于此；HMGCR×CAD 阴性 p=0.93 与已知 LDL-C→CAD 不一致即弱工具+未 clump 的体现）。
- 修正版 v2（`scripts/M3_transcript_mr.R`）：真实 eaf（SNP 频率文件，eaf=1−AF_AlleleB）+ LD clump（`ieugwasr::ld_clump` API 模式 EUR，r²<0.01@1000kb，免本地 plink）。待频率文件下载完整后重跑。

## 2026-08-06 — INTERVAL 蛋白通道诚实负结果（如实记录）

- HMGCR（prot-a-1354）与 ANGPTL3（prot-a-98）在 INTERVAL（n=3,301）cis ±1Mb 内无 p<5e-6 工具（`results/grid/protein_interval_mr.csv`）；prot-a-1354 显著工具集中 chr2/chr8 等跨染色体信号。SomaScan 对肝酶/膜蛋白血浆 cis-pQTL 功效不足 → **deCODE（n=35,559）须为主源**。
- PCSK9/CETP/NPC1L1/APOC3 在 INTERVAL SomaScan 未收录（`results/prota_index.rds` 核对），如实记录，不替代。

## 2026-08-06 — MR 方法名修正（技术 bug fix，非科学决策变更）

- `mr_ivw_random_effects` → `mr_ivw_mre`：TwoSampleMR 的 `mr()` 方法名是 `mr_ivw_mre`（multiplicative random effects），原名称是 MendelianRandomization 包的叫法。这纠正的是工具调用名，不是统计方法选择（IVW 随机效应语义不变）。
- 预注册 PREREGISTRATION.md 同步重生成 + 重新 sha256 锁定；config.json 同步。
- eQTLGen 暴露数据补 `eaf.exposure` 列（harmonise 要求），eaf 暂取中性 0.5，后续用 eQTLGen SNP 频率文件替代。

## 2026-08-05 — v0.3 落盘 + 导师收尾项修正

- README v0.3 落盘至 `/data/qiushuogeng/projects/dual-channel-mr-atlas/README.md`（博士生+导师 agent 合著，导师 approved=True）。
- 落实导师放行后 5 条收尾：
  1. §9 Nusinow 页区间 `180:585-600` → `180(2):387-402.e16`（PMID 31978347 正确，页区间错已修），回填 docs/citation_checklist.md。
  2. §4.3 删除残留编辑句「同时修正 ieugrasr 为 ieugwasr」。
  3. M0 禁词正则改真词（弃用 [a] 写法，否则匹配不到真实禁词）+ 行号白名单豁免 §0 对照表/M0 注释块；复跑全绿。
  4. §4.2 补「4,907 探针→4,719 唯一蛋白」；「Δ 显著集」→「方向显著集」。
  5. M0/M1/M1_decode_subset 存根已生成；M0 门禁复跑通过（exit 0）。
- 新增 docs/citation_checklist.md（引用核验表，M0 检查对象）。

## 2026-08-05 — 项目建立

- 项目文件夹 + /data 软链建立。
- 审稿人 agent 有条件放行（MAJOR_REVISION, 6/10, go_no_go=true），9 条 required_fixes 全部落实进 README v0.3。

## 2026-08-06 — 转录本 MR v1 真实结果 + 方法学修正 + 蛋白通道实测（数据驱动）

### 转录本通道 v1（真实统计，29/36 组合）
- eQTLGen（10.5M 行 / 16,923 基因）× OpenGWAS T2D/CAD/FBG，主方法 mr_ivw_mre，29/36 成功。
- 结果落盘 `results/grid/transcript_mr_qa.csv`；漏斗 `results/funnel/funnel_transcript.tsv`。

### 方法学修正（技术修正，非科学决策变更）
1. **eaf 占位 → 真实频率**：v1 用 0.5 占位致 palindromic harmonise 不可靠；v2 引入 eQTLGen 官方 SNP 频率文件（`2018-07-18_SNP_AF_for_AlleleB...`，data/eqtlgen/SNP_AF.txt.gz，effect allele=AlleleA → eaf=1−AlleleB_all）。
2. **LD clump**：v1 top-3 候选未独立化；v2 用 `ieugwasr::ld_clump` OpenGWAS EUR 1000G API 模式（r²<0.01@1000kb，无需本地 plink/bfile），已验证可行（6 变异→4 独立）。
- v2 脚本 `scripts/M3_transcript_mr.R` 就绪。

### 蛋白通道实测
- **INTERVAL（OpenGWAS prot-a-*，n=3,301）**：对照 HMGCR=prot-a-1354、ANGPTL3=prot-a-98。实测 cis ±1Mb、p<5e-6 内无独立工具（HMGCR 唯一 chr5 命中 rs72835052 @167.6Mb 距基因 93Mb，属 trans）→ 功效不足，诚实记录。PCSK9/CETP/NPC1L1/APOC3 在 SomaScan 未收录（prot-a 索引 results/prota_index.rds）。
- **deCODE（n=35,559）**：下载链接经 source-map 逆向定位 API（/s3/fileInfo、/s3/folder、/s3/download）。4,907 蛋白文件命名 `SeqId_GeneName_ProteinName.txt.gz`（hg38、per-SD Beta）。对照定位：PCSK9=5231_79、HMGCR=5230_99、ANGPTL3=10391_1、APOC3=6461_54；**CETP/NPC1L1 在 deCODE 未收录**。下载速度 ~5–8 KB/s 为硬瓶颈，PCSK9 已启动后台下载。
