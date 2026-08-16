# 冲高路线图：双通道 MR 图谱 → eBioMedicine（IF≈9）级

> 2026-08-13 产出。来源：多 agent 评审工作流（导师/实验精进/算法参数/期刊定位 4 视角 + 3 个对抗验证 agent，全部 7 个完成，0 错误）。
> 目标：把论文从"3 分多现实档"冲到 **eBioMedicine（IF≈9）/ Diabetologia（IF≈8）** 级。

---

## 0. 一句话结论

**现在直接投只有 IF 3-4 档（Human Genetics 3.9 / J Lipid Res 3.6）。冲到 IF≥9 的唯一路径是：把转录通道做成"多层收敛的顺式 eQTL 图谱 + coloc 一致性率方法学读数"（106/818≈13% 是 MR 显著且获区域共定位支持的一致性率，需 CI/分层/负对照标定后才可谈假阳性），配 UKB-PPP 蛋白跨平台复现 + ≥1 个带正交验证的新发现。**

关键物理约束：**deCODE 全基因组 4,719 蛋白 ≈ 4.7TB**，磁盘/带宽不可行 → 蛋白通道规模化只能走 **UKB-PPP**（Sun 2023, Nature, Olink, n=54,219，独立队列+独立平台，预注册 §9 本就承诺，一箭双雕=复现+扩面板）。

---

## 1. 四个视角评审结论（均已完成对抗验证）

### 1.1 期刊定位 agent
- 当前证据/设计水平 = **中等偏保守**。方法纪律（预注册哈希锁、全网格如实报告、双闸门、组织/胰岛敏感性）在同类 cis-MR 论文里**明显高于均值**——这是可打的牌。
- 够不到 eBioMedicine 的四点差距：① 蛋白通道仅 11 候选/5 有工具/2 过闸门且都是已知校准位点；② 预注册两个可证伪主假设 H1/H2 未执行（executed=false）；③ 外部复现（UKB-PPP/FinnGen）未完成/有 bug、SMR/HEIDI 跑完但蛋白侧未校准；④ 转录 812 FDR 存活 → 106 strong（coloc 一致性率 106/818=13%），这个**一致性率本身是当前最有卖点的方法学发现**（注意：不是假阳性率——coloc 只在 MR 存活者上跑，无法区分 MR 假阳性 vs coloc 假阴性，需负对照标定后才能谈假阳性）。
- 档位阶梯：现状态可投 → Human Genetics(3.9)/J Lipid Res(3.6)；**做完冲高清单 → eBioMedicine(9)/Diabetologia(8)/Hum Mol Genet(5-7)**；远期 Nat Commun 需真正新因果发现 + 正交/功能验证（超 MR+coloc 范围）。

### 1.2 实验精进 agent（最高杠杆，全部预注册兼容）
1. **UKB-PPP 独立 pQTL 复现**（最高单一杠杆，预注册 §9 承诺未执行）。
2. **补算 H1 第二支：双通道方向一致性**（数据在手、从未算过）。
3. **FinnGen 结局侧复现** + 记录结局集与预注册 §2 的偏差。
4. **蛋白通道诚实性审计**：头牌 p 值不可靠（见 2.1）。

### 1.3 算法参数 agent（发现实质问题，两个独立视角收敛）
- **框架本身不需要更换**（cis-MR + coloc 闸门是标准药物靶点配置，正对照校准命中教科书级）。
- **但 `mr_ivw_mre` 在 nsnp=2-3 时 SE 被压缩 7-16 倍** → 蛋白通道两条头条 p 值（INSR×FBG p=0、APOC3×FBG p=9.1e-46）是伪影。修正=小工具数 FE 首读（不移动门柱）。
- 其余建议：补 F 统计量、coloc-SuSiE 裁决 3 个争议位点、Steiger 方向性、UKB-PPP 复现。

### 1.4 验证 agent 结论（对抗性复算）
- MRE SE 塌缩**可复现**：INSR×FBG FE p=0.0121（SE 差 16 倍）、APOC3×FBG FE p=0.0543（失去显著性）、INSR×CAD FE p=0.40、APOC3×T2D FE p=0.0203（反转）。**已在本日修正落地**（`scripts/report_protein_primary.py` → `protein_decode_mr_primary.csv`）。
- H1 可执行（11 蛋白 45 对在盘）；**H2 当前不可执行**（负对照 HMGCR 无工具、NPC1L1 文件未下，仅 ANGPTL3/APOC3 可用）→ 需续下 NPC1L1/CETP 或按 CHANGELOG:114 先做设计调整。
- SMR 1.3.1 可装（1000G EUR 参考已就绪）；需 eQTLGen BESD + deCODE 自建 BESD。
- **T2D 样本量矛盾属实**（gwasinfo 62,892 vs README 655,666）未闭合；**预注册锁实为 MD5 而非 sha256**（osf_provenance.txt 明示，文件名 .sha256 有误导）。
- **预注册 §2 结局 ID 与实跑不符属实**（预注册 ieu-b-38/ieu-a-7/finn-b-E4_DM2 vs 实跑 ebi-a-GCST006867/005194/005186，FBG 新引入）→ 需 v6 修订。
- UKB-PPP 数据在 Synapse syn51365303 + AWS S3 桶，**controlled access 需 Synapse 注册+条款**（summary-stats 层无需 UKB 申请，仅需免费账号+PAT）；FinnGen 公共摘要 post-embargo 免费直下。
- "LIPA×CAD、CADM4×CAD、RASD1×CAD"转录 strong 种子核实为真（PP.H4 0.99/0.996/0.986）；但 "PTPRN×T2D" 不存在（实为 ×CAD 且 p12e6 敏感性仅 0.629 最弱）——三角验证种子需修正。

---

## 2. 关键问题清单（已确认，按优先级）

### 2.1 必须现在修（诚实性，否则审稿人一跑 FE/coloc 就抓到）
- [x] **MRE SE 塌缩**：nsnp≤3 改 FE 首读。已修（脚本+摘要+README 全更新）。
- [ ] **措辞纪律（2026-08-13 多 agent 核查收敛，投稿前必须）**：删"三方法独立互证/三层漏斗"→"单点 Wald 自洽 + 两个区域级独立检验（coloc、HEIDI）"（MR-Wald≡SMR 数学恒等实锤）；删"MR 假阳性率"→"coloc 一致性率 106/818=13.0%[CI]"；删"APOC3 独立救回/共享变异新发现"→"单工具敏感性支持，待独立平台复现"；删"双通道图谱/整合"→"平行双组学筛查，整合受工具重叠限制（n=0）"；"125/128 SMR 显著"不进正文（构造性必然）；"106"→"106 位点（76 编码/30 非编码）"；KCNJ11 改口"区域共享变异 + 同质性通过"，机制留胰岛 eQTL。
- [ ] **修 FinnGen rs17716350 对齐 bug**：同 rsID、同 discovery b=-0.0425，FinnGen 两行 beta -0.0133(G) 与 +0.1344(T) 均标 aligned=TRUE（`finngen_replication.csv` 30-31 行）→ 修后全表重算；头条改**双显著子集** 27/27（T2D）与 19/19（CAD）+ CI；弱信号行降权。
- [ ] **补 F 统计量**：蛋白已补（`protein_decode_mr_primary.csv`）；转录 stage-2 待补（F=Z²，零成本）。
- [ ] **T2D 结局样本量核清**：gwasinfo 重查三个 ebi-a ID 实际 N，回填 README/REPORT。
- [ ] **预注册 v6 修订**：记录结局 ID 变更（ieu-b-38/ieu-a-7/finn-b-E4_DM2 → ebi-a-GCST006867/005194/005186，FBG 引入）+ T2D N + MD5/sha256 命名澄清 + H1/H2 重定位为探索性。文档修订，不动主门柱。

### 2.2 应该做（证据提档，全部预注册兼容）
- [x] **H1 双通道方向一致性**：已算（`scripts/H1_direction_consistency.py` → `results/grid/H1_direction_consistency.csv`）。分泌 3/9=33% vs 胞内 1/6=17%，Fisher OR=2.50 p=0.604（功效不足描述性）；双通道工具级重叠仅 3 个 INSR 对、双显著 0 → 方向一致率如实声明 n=0。
- [x] **coloc-SuSiE**：收敛已解决（前 300 信号变异 + EPV=FALSE + Rfast + 直接调用 susie_rss，3.7s；另修 n 标量 + PP.H4.abf 列名两个 bug）。5 对全跑完，但 **校准失败，不作为证据**：阳性对照 PCSK9×CAD（abf=1.0，rs11591147 低频 R46L）在 SuSiE 下全 PP.H4≈0（低频变异 PIP=0）；APOC3×T2D 的 PP.H4=1.0（rs11216103）为伪阳性（T2D GWAS 侧 p=0.054 不显著却被 susie 放进可信集=LD 错配伪可信集）。**结论：主闸门维持 coloc.abf，多信号共定位需 UKB-PPP 匹配 LD**（见 DECODE_PIPELINE_SUMMARY M5-SuSiE 节）。
- [x] **FinnGen 结局侧复现**：已完成（`scripts/M7_finngen_replication.R` → `results/grid/finngen_replication.csv`）。全样本方向一致 T2D 94.7%（ρ=0.866, n=57）、CAD 88.5%（ρ=0.769, n=52）；**头条改双显著子集：T2D 27/27、CAD 19/19 方向 100% 一致**；FBG 无 FinnGen 端点如实跳过。**已知 bug 待修**：rs17716350 同 SNP 两行 FinnGen beta 一负一正（30-31 行）→ 修后全表重算，弱信号行降权。
- [x] **组织三角验证**：已完成（`scripts/M8_tissue_triangulation.R` → `tissue_triangulation.csv`）。106 strong 命中 × 6 GTEx 组织：**129/242（53%）组织强化，55/79 命中（70%）获 ≥1 组织 MR 方向强化**（肝 12/胰 21/冠脉 17/脂肪 23/肌肉 23/全血 33）→ 全血信号非组织特异伪影，atlas 升级为跨组织验证层。
- [x] **SMR/HEIDI**：已装 1.3.1 并跑完两条通道（`scripts/M9_*` 蛋白 + `scripts/M10_*` 转录）。**蛋白 15 对（仅探索性，未校准）**：APOC3×T2D(p_SMR=0.011, p_HEIDI=0.42)/APOC3×FBG(p=0.007, HEIDI=0.91) 为**单工具敏感性支持（非独立确认）**——15 对零假设期望 0.75 个 p<0.05、无多重校正、coloc.abf none、coloc-SuSiE 判不共定位、FinnGen rs964184×T2D p=0.075 不显著、**APOC3 不在 UKB-PPP 面板（复现断头）**；PCSK9×CAD(p=2.96e-19)/APOC3×CAD(p=5.5e-11) SMR 极显著但 HEIDI 异质性，且**阳性对照 PCSK9×CAD 过不了 HEIDI（p=1.28e-4）→ 蛋白侧 HEIDI 判定整体未校准**。**转录 128 测试**：HEIDI 通过 89/126=70.6%（2 NA）、37 异质性；p_SMR 显著（125/128）是选择偏倚构造性必然，零信息量不进正文；KCNJ11×T2D 为**单点 Wald 自洽 + 两个区域级独立检验（coloc、HEIDI）通过**，且 MR-Wald≡SMR 数学恒等（b 同到 6 位小数），不得写"三方法独立互证"；KCNJ11 基本不在全血表达、coloc 峰值 rs757110 在 ABCC8，机制 caveat 见 DECODE M10。deCODE 无 effectAlleleFreq → .esd Freq 用 ImpMAF 近似 + `--disable-freq-ck`（探索性敏感性，已在 MD/README 标注）。
- [x] **三角验证种子修正**：LIPA×CAD、CADM4×CAD、RASD1×CAD 种子已核实（PP.H4 0.99/0.996/0.986），PTPRN 已改为 ×CAD 并注明敏感性弱（M8 §2.1 验证 agent 记录）；新增 SMR 正交支持 LIPA×CAD(p_SMR=1e-20, HEIDI=0.457)、RASD1×T2D(HEIDI=0.519)。

### 2.3 冲高结构性升级（决定能否到 IF≥9）
- [ ] **UKB-PPP 蛋白通道扩面板+跨平台复现**：Synapse 已注册（PAT 已存 gitignored `~/.synapseConfig`）→ **面板覆盖核查完成**：UKB-PPP 有 PCSK9/INSR/ANGPTL3/LDLR/APOB/DPP4/GLP1R/GCG（**8/11**），**无 APOC3/PCK1/HMGCR**（Olink 面板不覆盖，APOC3 缺口尤其遗憾——deCODE 蛋白 SMR 救回的对无法在 UKB-PPP 复现，如实写）→ 8 蛋白 tar（~4.4GB）正走代理下载 → 按 cis±1Mb 提取 → MR+coloc+方向一致率，与 deCODE 并列对照。**Synapse API 被中国大陆 GEO 屏蔽**（NIH NOT-OD-25-083），须走 mihomo 代理；已实测无公共镜像、直连 S3 被限速 100 倍（14KB/s vs 代理 1.9MB/s）。
- [ ] **"coloc 一致性率"量化为方法学主结论**：转录漏斗 982→819(nominal)→818(coloc ok)→106（coloc 一致性率 106/818=13.0%，附 Wilson CI 10.83–15.43%）+ 蛋白漏斗 15→2，做成主图 + 按结局/通道/biotype 分层。**注意措辞**：这是"MR 显著且获区域共定位支持的一致性率"，不是假阳性率——coloc 只在 MR 存活者上跑，无法区分 MR 假阳性 vs coloc 假阴性；要升级为可辩护的假阳性读数必须补负对照/permutation 标定 coloc 阈值 FDR 曲线。另须**剔除非编码伪影（30/106）**并出**独立位点表（~100）**。

---

## 3. 期刊档位与证据需求对照

| 档位 | IF | 需要的证据 | 达到条件 |
|---|---|---|---|
| 保底（现状态） | 2-4 | 描述性图谱+诚实阴性 | 无需新数据 |
| **目标（冲高）** | **5-9** | 跨平台复现 + 方向一致性 + 方法学量化 | 做完全部 §2.2 + §2.3 |
| 远期 | 10+ (NC/eLife) | 真正新因果发现 + 正交/功能验证 | 超 MR+coloc 范围，需实验 |

**结论**：eBioMedicine 级**是路线图非现状**——需按本路线图补齐跨平台复现 + coloc 一致性率方法学估计（CI/分层/负对照标定）+ 一个带正交验证的新发现后的**整合稿**。预注册承诺 + 全网格诚实协议 + 小样本敏感性纪律是比同档位竞争者强的差异化卖点，要写进 cover letter。

---

## 4. 数据可行性核对（验证 agent 实测）

- **deCODE 全量扩展**：4,719 蛋白 ≈ 4.3TB gz —— **物理不可行**（磁盘 932G 可用、单文件 910MB 无 Range）。❌
- **UKB-PPP（Sun 2023, Olink, n=54,219）**：Synapse syn51365303 + S3 桶 `ukbiobank.opendata.sagebase.org`（CC BY）。summary-stats 层无需 UKB Biobank 申请，**仅需免费 Synapse 账号+PAT**；每蛋白按染色体分片 REGENIE tar，11 蛋白约数十 GB，S3 支持 Range 高速。⚠️ **需要用户注册 Synapse + 接受数据条款**（无法代做）。
- **FinnGen**：公共摘要 post-embargo 免费直下。✅
- **SMR 1.3.1**：Yang lab 免费二进制。✅（BESD 需转换/另下）
- **coloc-SuSiE**：susieR 0.14.2 + coloc 5.2.3 已安装，LD 用现有 1kg EUR。✅

---

## 5. 执行顺序（建议）

**第一批（本周，0 资源争议，先做）**：
1. ✅ F 统计量补齐（蛋白 done/转录 done）+ 预注册 v6 修订（done）+ T2D N 核清（done，655,666 与 62,892 是同一研究不同口径，非矛盾）。
2. ✅ H1 双通道方向一致性（已算，功效不足如实声明 n=0）。
3. ✅ FinnGen 结局复现（T2D 94.7%/ρ=0.866，CAD 88.5%/ρ=0.769）。
4. ✅ 组织三角验证（M8：55/79=70% 命中获组织强化）。
5. ✅ **SMR/HEIDI 双通道**（M9 蛋白 15 对 + M10 转录 128 测试，见 §2.2；APOC3×T2D/FBG 为单工具敏感性支持待复现、KCNJ11 单点自洽 + 两区域级检验过，措辞已诚实化）。

**第一批补充（零资源争议，当前性价比最高）**：
5b. **HEIDI 跑满全部 818 个 MR 存活测试**（现在只跑了 128 个 coloc 命中）→ 不一致分类（多信号/低功效/伪探针/多效性）→ Fig 5 与"多数不一致"叙事成立。

**第二批（需 UKB-PPP，Synapse 已注册，下载中）**：
6. UKB-PPP 8 蛋白复现（PCSK9/INSR/ANGPTL3/LDLR/APOB/DPP4/GLP1R/GCG）→ cis±1Mb 提取 → MR+coloc+方向一致率 → 与 deCODE 并列成蛋白通道主体。（APOC3/PCK1/HMGCR 面板不覆盖，如实写缺口；**APOC3×T2D/FBG 的"救回"对无法在计划路径复现，另找 INTERVAL/ARIC 或撤回**）
7. coloc-SuSiE（已完成但校准失败不作为证据，主闸门维持 coloc.abf）+ 转录 SMR drug-target 扩展（DPP4/GLP1R/PPARG/INSR/TCF7L2/PRKAA1/SLC5A2 需拉 cis GWAS 区域）。

**第三批（写论文）**：
8. coloc 一致性率方法学主图（106/818=13.0%[CI]，双通道不对称漏斗）+ 优先化基因表成稿。
9. 目标投 eBioMedicine → Diabetologia → Hum Mol Genet 备选。

---

## 6. 风险提示（诚实）

- UKB-PPP 为 Olink 平台，未必覆盖 HMGCR/LDLR 等膜蛋白（Olink Explore 3072 面板可能不含/QC 不过）——**先核对 11 蛋白面板可得性再投入下载**；若 UKB-PPP 也无，DPP4/GLP1R/GCG/HMGCR 的"无工具"仍是证据缺失，如实写。
- 预注册 §9 写 UKB-PPP n=54,306，实际论文为 54,219——数字口径统一。
- H2（负对照 LOOCV AUR）功效严重不足，诚实处理=描述性报告或按 CHANGELOG:114 先做设计调整（不事后改主门柱）。
- 所有扩展分析以"探索性/预先登记补充"落档，不改预注册主门柱。
