# 对抗性学术诚信审计报告 — dual-channel MR atlas（2026-08-16）

> 审计对象：`/data/qiushuogeng/projects/dual-channel-mr-atlas`（目标投稿 HMG，方法学+图谱资源论文）
> 审计方式：独立复核全部面向文档 + 用 python 对 `results/` 下 CSV 逐项重算（小任务，未申请资源）
> 审计日期：2026-08-16。审计未修改任何文件。
> 结论先行：**未发现捏造/伪造/抄袭**；全部头条数字（131 strong、12.9%、3.0%、81%、19%、0/6000、置换 1.45%）经重算**成立**。但存在**一处口径级选择性报告风险（3.0% vs 12.9% 从未并排披露）**、论文面向文档（README / PUBLICATION_STORY）**未同步 2026-08-15 build 修复与 M20–M24 全量结果**、以及若干内部数字不一致。按"宁可多报"原则，凡有疑点均列入。

---

## 〇、已用 python 重算验证为真的关键声称（正面清单）

审计不仅报负面。以下数字我逐项重算过，**站得住**（详见 §3 表）：

| 声称 | 出处 | 重算值 | 判定 |
|---|---|---|---|
| M24 全量：n_all=31,373 / n_qc=31,371 / n_sig=4,248 / strong=131 | coloc_full_summary_20260815.csv | 31,373 / 31,371 / 4,248（t2d 1557+cad 2242+fbg 449）/ 131 | ✅ |
| strong 131 分解：t2d 67 / cad 60 / fbg 4 | 记忆/摘要 | t2d 67（66 sig + AP3S2）、cad 60（59 sig + ZNF19）、fbg 4（全 sig） | ✅ |
| precision_sig=3.0% | summary CSV | 129/4248=0.0304 | ✅（**但全文档无散文叙述**，见 P0-1） |
| 106/106 已知 strong 全量重现；new_in_sig=23；known_missing=0 | summary CSV | 106/106；129−106=23 | ✅ |
| MR sig 集外 strong 全量仅 2（AP3S2×t2d / ZNF19×cad，均 GWAS 峰显著） | feasibility §7.7 | AP3S2 pp4=0.935、gwas 5.5e-11、F=55.6、mr_p=0.38；ZNF19 pp4=0.903、gwas 2.9e-11 | ✅ |
| M22b 文献命中 81%（rsid 26 + gene_100kb 59 + gene_anno 1）+ none 20=19% | feasibility §7.8 | 86/106=81%（t2d 52/58、cad 33/46、fbg 1/2） | ✅ |
| 旧 M22 55%（58/106）为 build 错配低估 | feasibility §7.8 | 58/106=54.7%；old/new 命中状态 42/106 改变 | ✅ |
| M20 试点 6000 对 sig 外 strong=0 | feasibility §2 | 灰区 3000（QC 2998）+阴性 3000，PP.H4≥0.8 全 0；灰区 ≥0.5 仅 1（RP11-403I13.4 pp4=0.676） | ✅（**但表格列写 0，见 P2-1**） |
| AP3S2/ZNF19 未被 M20 抽样抽到（抽样遗漏，0/6000 与全量 2 不矛盾） | feasibility §7.7 | 两者都不在 pilot 6000 对中 | ✅ |
| 置换标定 E1：FP=1.45%（154/10600），期望 FP 1.54/106，ZNF34×t2d 22% 最脆弱 | PUBLICATION_STORY / coloc_permutation | 154/10600=1.4528%；顶部 ZNF34 位点 fp08=22 | ✅ |
| HEIDI 全集：398/819=48.6%；strong 子集 76/106=71.7% | heidi_full_summary.json | 398/819、76/106 | ✅ |
| GTEx 复现 44/63=69.8%、非全血 26/40=65%、同变异 6/6 | gtex_replication | 同值 | ✅ |
| Steiger 96.1%（73/76） | steiger_direction | 73/76 | ✅ |
| InsPIRE 负性核查：11 显著 10 反向、KCNJ11 胰岛 null p=0.87、可复现口径 0/13 | inspire_islet | 同值 | ✅ |
| FinnGen 双显著子集 46/46（T2D 27/CAD 19）=100% 方向一致；108/108 可对齐 | PUBLICATION_STORY M2 | 46/46（两结局 p<0.05 双显著），rs17716350 bug 行已在文件修复（现单行） | ✅（但 README §5 仍写"未修"，见 P1-4） |
| 转录漏斗 49,866→982（394/576/12）→stage2 nominal 819（331/477/11）/FDR 812；157 p<1e-6；906 唯一基因；978/982 完成 | README §4.1 | 全部吻合 transcript_grid_mr.csv / stage2 | ✅ |
| transcript_coloc 819 行：106 strong/149 moderate/563 none/1 NA | README §4.1 | tier 计数 106/149/563/1 | ✅ |
| 106 已知：58 t2d/46 cad/2 fbg；PP.H4≥0.9 共 62；100 唯一 top-SNP | PUBLICATION_STORY/README | 62；100 唯一 top_snp | ✅ |
| 76 编码 + 30 非编码（biotype） | README §8/ukbpp CSV | ukbpp_coverage 76 行全匹配；106−76=30 | ✅（"56+20"旧口径已废，见 P2-4） |
| ukbpp 覆盖 10/76 | ukbpp_coverage | 10 True | ✅ |
| Wilson CI 106/818 = 10.83–15.43% | PUBLICATION_STORY | 重算 [10.83, 15.43]%（106/818；819 中 TCF19×CAD 因 GWAS 区域数据不可用无法评估 coloc） | ✅ |
| eQTL 功效不分辨：106 strong |Z| 中位 16.4、pctile 中位 ~50 | feasibility §4/M22 CSV | 16.41；lead_pctile_global 中位 50.4 | ✅ |
| fbg 4 个 strong 全 GWAS 峰不显著、MT3 F=17/hsa-mir-296 F=7.4（弱信号，已降级） | feasibility §7.5 | 全部吻合 | ✅ |
| GIT1 修复后 pp4=0.8937 与 M5 API 版逐位一致（本地 full 与 API 双源收敛） | feasibility §7.7 | 106/106 已知全部被全量扫描重现（pp4≥0.8） | ✅ |

**正面结论**：核心数据链（网格 → stage2 → coloc → 全量扫描 → 文献基准）内部一致、可复现、无捏造迹象。build 错位的**内部**披露充分（CHANGELOG:200–205、feasibility §7.7、BREAKTHROUGH_PLAN §8.4，含误结果存档 `results/archive/buildbug_20260815/`）。

---

## 一、P0 — 学术不端风险，投稿前必须改

### P0-1 两个"精度"从未并排披露；文档把 12.9% 标为"唯一真读数"，全量扫描给出的 3.0% 只躺在 CSV 里
- **现象**：
  - `results/coloc_full_summary_20260815.csv` 明确给出 `precision_sig = 129/4248 = 0.0304`（M24_summarize.py 的注释就叫"精度（操作特性）：MR 显著集内 strong 占比"）。
  - 但**没有任何 .md 文档提到 3.0% / 0.030 / 4248**。feasibility_20260815.md §1 与 §5（line 13/48/85）反复称"MR 显著集内 106/819 = 12.9% 精度是**唯一真读数**"。README §8 / PUBLICATION_STORY 主图计划均以 12.9% 为头条。
  - 已重算确认 **819 ⊂ 4248**（819 个 stage-2 存活对的 mr_p 在全量网格里全部 <0.05）。两个数字是**同一个操作定义、同一根数据链、嵌套的两组分母**：扫描层名义显著 4,248 对 → strong 129（3.0%）；stage-2 clump+IVW 复核存活 819 对 → strong 106（12.9%）。分母不同，精度差 ~4 倍。
- **为什么是 P0**：这是论文最核心的定量声称。只报 12.9%（更"好"的数字）而不披露同管线直接算出的 3.0%（更"差"的数字），审稿人一跑数就会被抓成"口径不一致/选择分母"，且属于对不利读数的不披露。文档自己最新一份评审建议（`docs/REVIEWER_IF_ADVICE_20260816.md` §0/§2.3/§6 第 8 条）也明确把它列为"必须主动披露，否则被当把柄"，并建议并排报两个 precision 做成漏斗图——说明项目已意识到，但尚未写进任何论文面向文档。
- **改法**：论文并排报告 `扫描层精度 129/4248=3.0% [2.56–3.60%]` 与 `stage-2 复核层精度 106/818=12.96% [10.83–15.43%]`（2026-08-16 口径统一：819 中 TCF19×CAD 因 GWAS 区域数据不可用无法评估 coloc，精度分母用可评估的 818 对），写清各自分母与漏斗关系（4248 ⊃ 819 ⊃ 818；129 = 已知 106 + 新 23），把"复核步骤把共定位支持率提升 ~4 倍"做成明确结论（这正是 REVIEWER_IF_ADVICE 给的诚实卖点）。删除"唯一真读数"这类措辞。

### P0-2 论文面向文档（README / PUBLICATION_STORY）完全没有 build 错位修复与 M20–M24 全量结果
- **现象**：README 版本仍为 v0.12（2026-08-13），§0/§4/§8 讲的全是"49,866 网格 → 982 → 819 → 106 strong"、且 §0 line 34 把全转录组扩展标为"🔲 计划中"。PUBLICATION_STORY_20260813.md 全文无 build bug、无 M20/M21/M22/M23/M24。而事实上 M23/M24 全量扫描已完成（31,373 对、131 strong、sig 外 2、precision 3.0%、106/106 重现），build 错位已修复。
- **为什么是 P0**：若论文按 M24 全量结果写，现在 README 与 PUBLICATION_STORY 的全部结果章节需要重写；若论文按 08-13 快照写，则必须说明为何不使用已完成的全量扫描，否则属于"做完的分析不进论文"。build 错位的**内部**披露是充分的（CHANGELOG/feasibility/BREAKTHROUGH_PLAN），但**论文/方法学面向文档里没有**——投稿 Methods 必须如实描述坐标处理（GWAS hg38 + eQTL hg19，最终改用纯 rsid 匹配），并说明"曾用坐标窗口匹配出错、经交叉验证发现后全部重跑"这一修复史（可写进 Reproducibility 段落，属加分诚实项而非减分项）。
- **改法**：投稿前把 M20–M24 全量结果 + build 修复史同步进 README 与论文故事线；版本号/日期一并更新。

### P0-3 全量扫描的 131 个 strong 缺全量尺度置换标定，且文档声称"多数应峰显著"被数据否定
- **现象**：E1 置换标定（FP=1.45%）只在 **106 个 MR 显著区域**上做；verification §3 自己承认"不能外推到全扫描"，BREAKTHROUGH_PLAN/feasibility §6 计划"置换零假设对全扫描重导"——但 **M24 未做全量置换标定**。与此同时 131 个 strong 是对 31,373 次 coloc 检验的调用，无全尺度校准。
- 另外我已重算：**106 个已知 strong 中仅 41 个（39%）的区域 GWAS 峰 p<5e-8**（5e-8–1e-6 22 个、1e-6–1e-4 42 个、≥1e-4 1 个；gwas_min_p 中位 5.1e-7）。feasibility §7.5（line 131）写"已知 106 中 t2d 58/cad 46 **多数应峰显著**"——**数据不支持**（t2d 24/58、cad 17/46、fbg 0/2 峰显著）。该句在同一节里把 fbg 的"两侧均无强信号 → PP.H4 不可靠"作为降级理由，却假定 t2d/cad 不适用，从未核对。
- **为什么是 P0/P1 交界**：把 39% 的"GWAS 峰显著率"与自定的弱信号降级标准放在一起，一个严格审稿人可以据此质疑 106/131 中相当比例 strong 的可信度。这是"未披露的 caveat"，按任务要求"宁可多报"，列为 P0-3。
- **改法**：① 补全量置换标定（或明确写"131 个 strong 调用未做全量尺度置换标定"进 Limitation）；② 报告 131/106 的 gwas_min_p 分层分布；③ 删除"多数应峰显著"这句（或改为如实数字）。

---

## 二、P1 — overclaim / 未披露 caveat，建议改

### P1-1 "19%（20/106）新候选上限"只对已知 106 成立；全量扫描新增的 23 个 strong 未做 GWAS Catalog 注释
- m22b 的 81%/19% 只注释了 **106 个已知 strong**。M24 全量又发现 **23 个 sig 内新 strong**（如 CD101、C2orf49、RBM6、CLEC3B、PDCD6、CWF19L1、LAMC1、TAGLN2、CNNM2、SLC12A3、PLAUR、ZNF100 等），**从未跑 m22b 注释**。若论文把"19% 超出已报道注释 = 新候选上限"套到 131/129 上，是错误的。
- **改法**：对全部 131 个 strong 逐位点跑 GWAS Catalog 注释，或明确写"19% 上限仅基于 M5 已知 106 集，全量新增 23 个尚未注释"。

### P1-2 T2D 结局样本量仍未闭合（655,666 vs 62,892）
- README §2.2/§5 limitation 5（line 92/269）标"待核：投稿前回填"；PREREGISTRATION v6 说两数是同一研究（Mahajan 2018 GCST006867）两套元数据口径、非矛盾。投稿前必须用 gwasinfo 实查并写入 Methods/数据表，否则头部结局 GWAS 描述不完整。

### P1-3 弱 GWAS 信号下 coloc 高 PP.H4 的 caveat 未对 t2d/cad 披露（承接 P0-3 数据）
- 39% 峰显著 + 22% 处 5e-8–1e-6 边缘。论文 Discussion/Limitation 应统一写清楚"MR 显著集内的 strong coloc 大部分位于未达全基因组显著的区域，coloc 在此区域的行为需与置换标定、HEIDI、Steiger 一起解读"，而不是只在 fbg 内部降级。

### P1-4 README §5 limitation 7 与事实矛盾：FinnGen rs17716350 对齐 bug 已修，README 仍写"已知未修"
- README line 271 写"`finngen_replication.csv` 第 30–31 行 rs17716350 同 SNP 两行一负一正均 aligned=TRUE…已知未修"；但文件已是修复态（108 行、rs17716350 只剩 1 行 b2=−0.0133(G) 与 discovery 同向、108/108 aligned），PUBLICATION_STORY M2 与 verification 均记为"已修（92.6%、双显著 46/46=100%）"。README 该条是**过期的假阴性陈述**，投稿前须改写为"已修复"。

### P1-5 多重检验问题未在论文口径中讨论
- 131 strong 是 31,373 次 coloc 检验的调用；MR 用 BH-FDR（stage-2 812/819），但 coloc 的 PP.H4≥0.8 调用没有家族级显著性声明，置换标定也仅覆盖 MR 显著子集。论文需明确"图谱为描述性/假设生成性质，不做家族级显著性宣称"（REVIEWER_IF_ADVICE §4 第 11 条已有预案，需落实到正文）。

---

## 三、P2 — 小问题 / 内部不一致

| # | 位置 | 问题 | 数据/证据 |
|---|---|---|---|
| P2-1 | feasibility §2 表格（line 40） | 灰区"PP.H4≥0.5"列写 **0**，但同节 note 与数据都是 **1**（RP11-403I13.4 pp4=0.676） | 重算 grey ≥0.5 = 1 |
| P2-2 | BREAKTHROUGH_PLAN §8.2（line 170） | "≥0.5 仅 2" 是 **build 错位版**数字；修复版为 1。此处未加修正注 | feasibility §2/§7.7 |
| P2-3 | verification_20260815 §10b（line 99） | 报 M22 旧 55%/45% 未标注"已被 build 错配修复推翻"（修正只在 feasibility §7.8 / CHANGELOG） | 55%=58/106 为错配值 |
| P2-4 | verification_20260815 §6（line 55） | "76 基因口径错误：实际 56 编码 + 20 非编码"已被 M21 推翻（正确 76 编码 + 30 非编码）；verification 文档未加更正注 | ukbpp_coverage 76 行全匹配 |
| P2-5 | PUBLICATION_STORY 摘要（line 159–160） | "screen of ~1,500 genes"应为 16,622 基因；"70.6% pass HEIDI"是旧 89/126，现统一口径 398/819=48.6% 或 strong 76/106=71.7% | 网格 49,866 对 = 16,622 基因 × 3 |
| P2-6 | BREAKTHROUGH_PLAN §1/§2（line 43） | "coloc-only 且经多指标复核的基因都是新分析产出的**新候选**"是评审前旧框架，已由后续修正（§8.2 改"上限"）覆盖，但该行仍在正文 | 建议加"已按评审修订"注 |
| P2-7 | 图谱分母口径 | "全转录组 coloc 图谱（49,866 对）" vs 实际 coloc 运行 31,373 对（63%；18,493 对无 eQTL 工具或 pval NA）。CHANGELOG 披露了 31,373，但论文口径需写清 coloc 测试分母 | grid 49,866 = 31,373 有效 + 18,493 NA |
| P2-8 | figures/ 陈旧 | 仅有 5 张 08-11 旧图（F1/F2/F3/F9/F10）；论文计划 Fig 1–5（PUBLICATION_STORY）未渲染；figures.md（M8）F1 引用的 funnel_transcript_v2.tsv 是 12 基因 drugtarget 时代漏斗（tested_genes=12），不是 49,866 全漏斗。投稿前须按当前数据重渲染 | funnel_transcript_v2.tsv |
| P2-9 | README 头部 | 版本 v0.12（2026-08-13）但 §0/§8 已含 08-15 评审修订，版本/日期不同步 | README.md:4,21,299 |
| P2-10 | PREREGISTRATION | `.sha256` 锁文件实为 MD5 裸哈希（osf_provenance + IMPROVEMENT_STRATEGY 已披露）；结局 ID 与实跑不符已在 v6 修订如实记录；"预注册"严格意义上是"设计先于结果 + 可能 post-hoc 提交"，osf_provenance §3 已诚实说明 | 无需改，仅提醒投稿时按披露口径自证 |

---

## 四、统计完整性专项（任务第 5 项）

- **131 是否有多重检验问题未讨论**：有，见 P0-3/P1-5。MR 侧做了 FDR；coloc 侧无全尺度置换标定。
- **6000 抽样 → 全量枚举的不一致是否如实**：是，如实。M20 抽样 0/6000 与全量 t2d 仅 AP3S2（+cad ZNF19）的差异被明确归因为"抽样未抽到 AP3S2"（feasibility §7.7 line 150），并改为"极罕见（全量 t2d 1/10,190=0.01%）"而非"严格 0"。Poisson 上界（REVIEWER_IF_ADVICE：2/31,373 的 95% CI 上界 ≈0.023%/对）口径合理。
- **CI 是否给出**：106/818 有 Wilson CI（10.83–15.43%，已重算吻合）；129/4248 的 3.0% 有 CI [2.56–3.60%]；0/6000 只有 Poisson 上界，无 CI——需在论文里给紧 CI 以支撑"决定性边界"。

## 五、build bug 披露专项（任务第 4 项）

- **内部披露**：充分、诚实。根因（OpenGWAS hm_pos=hg38 vs eQTLGen=hg19）、影响（M23/M20 t2d/cad 全作废、fbg 幸免、M20↔M23"同错"）、修复（纯 rsid 匹配）、验证（GIT1 双源 pp4 逐位一致、58/58 重现）、误结果存档——均在 CHANGELOG:200–205、feasibility §7.7/§7.8、BREAKTHROUGH_PLAN §8.4 完整记录。错位版误判（MAST4/ITCH/VAT1 假阳性撤回、WNT3 假象、C15orf62 误命中）全部如实撤回。
- **论文面向文档**：README 与 PUBLICATION_STORY **零提及**（二者均为 08-13 版本）。**必须补**，见 P0-2。
- 注意一个优点：build bug 是被**健全性交叉验证（M24 只重现 78/106）主动拦截**的，不是被隐瞒后暴露的——这一点在投稿 Reproducibility 段落是加分项，应保留叙述。

## 六、来源/引用专项（任务第 7 项）

- eQTLGen（Võsa 2021, n=31,684）、deCODE（Ferkingstad 2021, n=35,559, 4,719 蛋白）、INTERVAL（Sun 2018, n=3,301）、UKB-PPP（Sun 2023, n=54,219，预注册曾写 54,306 已在 IMPROVEMENT_STRATEGY 记录口径统一）、GTEx v8、GWAS Catalog 出处均有 PMID/DOI/样本量记录（`docs/DATA_AVAILABILITY_20260816.md`（合并自 `docs/archive/data_sources.md`）、`docs/archive/citation_checklist.md`）。`citation_checklist.md` 的"待人工核验"节如实列出未闭合引用（2026 预印本、Yuan 2026、VTE 2026 等缺卷期/DOI）。**未发现编造来源**。
- 唯一的"出处待核"是结局 N：T2D 655,666 vs 62,892 未闭合（P1-2）；OpenGWAS 结局 ID 与预注册不符已在 v6 修订记录。

## 七、图/表声称专项（任务第 6 项）

- 现有 figures/ 为 08-11 旧产物，与论文计划（Fig 1–5）不对应（P2-8）。无"图注与数据矛盾"的实锤——因为论文图尚未渲染。注意漏斗图必须用当前数据（49,866→31,373/982→819→106→131 的新分层），不要沿用 12 基因 funnel。

---

## 一句话总结

**这个项目可以诚实投稿，但现在不能**——核心数据全部重算成立、无捏造、build 错位内部披露充分，但必须先：(1) 并排披露 3.0% 与 12.9% 两个精度及其漏斗关系（P0-1）；(2) 把 M20–M24 全量结果与 build 修复史同步进 README/论文故事线（P0-2）；(3) 报告 106/131 的 GWAS 峰显著分层并补全量置换标定口径（P0-3）；以及闭合 T2D N、修正 README 里已过期/错误的两处描述（FinnGen"未修"、feasibility 表 0/1）。做完这些，它就是一份经审计的、诚实的阴性-方法学-图谱资源论文。
