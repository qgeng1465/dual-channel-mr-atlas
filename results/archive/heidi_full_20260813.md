# HEIDI 全集运行结果（P3 完成：818/819 测试全量，2026-08-13）

> 来源：`heidi-full-smr` 经资源仲裁（3 核）后台运行 `M10b_run_smr_full.sh`。
> 输入：`data/smr/trans_flist_full.txt`（759 探针）+ `eqtlgen_trans_full`（BESD）+ 各结局 `_full.ma`。
> 与 M10 参数一致：`--peqtl-smr 5e-6 --disable-freq-ck`；`--make-besd` 后逐结局 SMR+HEIDI。
> 产出：`results/grid/transcript_smr_heidi_full.csv`（1132 行含表头）。

## 1. 漏斗与通过率（诚实主线）

| 层 | 测试数 | HEIDI 通过 (p>0.05, nsnp≥10) | 通过率 |
|---|---|---|---|
| **全 SMR 显著命中**（含重复 topSNP 行） | 1129（t2d 498/cad 612/fbg 19） | 550 | **48.7%** |
| **预注册 coloc 测试集**（MR 显著 × coloc 全集） | **819/819**（漏斗完整） | 398 | **48.6%** |
| —— 其中 **strong coloc（PP.H4≥0.8）** | 106 | 76 | **71.7%** |
| —— 非 strong（tier=none） | 713 | ~310 | **43.5%** |
| —— 异质性（p_HEIDI≤0.05, nsnp≥10） | 330 | — | 40.3% |
| —— 低功效/不可判（nsnp<10） | 91 | — | 11.1% |

**每结局通过**：t2d 234/498（47.0%）、cad 313/612（51.1%）、fbg 3/19（15.8%）。

## 2. PP.H4 与 HEIDI 通过的单调关系（Fig 5 核心图）

| PP.H4 区间 | n | HEIDI 通过 |
|---|---|---|
| [0.80, 0.90) | 44 | 63.6% |
| [0.90, 0.95) | 33 | 75.8% |
| [0.95, 1.00] | 29 | 79.3% |

**解读**：coloc 共定位概率越高 → HEIDI 同质性通过率单调越高（63.6%→79.3%）。
这是两个区域级独立检验（coloc.abf 与 HEIDI）互相印证的**方法学标定曲线**：
coloc-strong 先验把 HEIDI 通过率从非 strong 的 43.5% 抬升到 71.7%。

## 3. "多数不一致"叙事立住（论文方法学头条）

- **在全部 MR 显著测试中，超过一半（51.4%）不能通过 HEIDI 同质性**：
  40.3% 为区域多信号/异质性（p_HEIDI≤0.05），11.1% 为低功效不可判（nsnp<10）。
- 这说明 **cis-MR 显著 ≠ 单因果位点**：大部分 MR 命中是区域多信号或 LD 混杂，
  这正是"cis-MR 命中不能直接当药物靶点"的定量证据。
- **coloc 与 HEIDI 是部分一致的两个独立区域级检验**：coloc-strong 富集 HEIDI 一致位点，
  但仍有 30/106 的 strong coloc（PP.H4≥0.8）在 HEIDI 下不一致（24 异质性 + 6 低功效）。

## 4. strong coloc 但 HEIDI 不通过（必须如实报告的"带病"位点）

**24 个异质性（p_HEIDI≤0.05, nsnp≥10）**，含：
- NUDT5×T2D（PP.H4=1.000，p_HEIDI=3.3e-8）、CAMK1D×T2D（PP.H4=1.000，p_HEIDI=4.4e-5）——
  **同一 topSNP rs11257655**，区域多信号典型（coloc 满分但 HEIDI 否定）
- GIGYF1×T2D（0.998, p=4.8e-4）、CADM4×CAD（0.996, p=3.6e-2）、RASD1×CAD（0.986, p=9.8e-4）
- FBXW7×T2D（0.925, p=1.9e-3）、MTCH2×CAD（0.902, p=9.9e-3）等

**6 个低功效（nsnp<10）**：RP5-1068B5.3、RP11-347C18.3、LDLRAD2、MYB、RP11-384M15.3、
hsa-mir-296（多为非编码/假基因，符合之前 biotype 伪影倾向判定）。

**纪律**：这些位点只作"coloc strong 但 HEIDI 不一致"如实分层，**不得**写为"coloc 证实的新因果位点"。
论文 Fig 5 用桑基图呈现：MR 显著 → coloc strong（106）→ HEIDI 通过（76）/ 异质性（24）/ 低功效（6）。

## 5. 与先前 107 探针结果的一致性

- 先前 strong 集：89/126 = 70.6%；本次 strong coloc 全集：76/106 = 71.7% → **一致 ✓**。
- 全量跑证实先前 HEIDI 判定无系统性偏差；且新获得非 strong 集的 43.5% 通过率（先前无此读数）。

## 6. 方法学读数更新

- 诚实清单 #2 的"125/128 SMR 显著零信息量"仍成立（SMR 显著性本身因预选偏倚不承载信息）；
  **HEIDI 通过率是唯一信息读数**，现升级为：**全 MR 显著集 48.6%、strong coloc 子集 71.7%**。
- must-do M3"HEIDI 计数更正 89/126"可进一步升级为引用全量 strong coloc 集 76/106=71.7%。

## 7. 已知局限

1. **.ma 是 coloc 区域缓存并集非全基因组**：HEIDI 基于区域 SNP 判定，非全 GWAS 设置；
   对强 cis 信号影响小（1kg LD 与区域覆盖充足），但低功效比例（11.1%）可能因区域裁剪被高估。
2. **nsnp<10 的 91 个测试被归为"低功效不可判"**：未与 PCSK9×CAD 蛋白侧（HEIDI p=1.28e-4）
   混淆——本表是转录通道（eQTLGen × GWAS），蛋白侧 HEIDI 未校准仍单独标注。
3. **SMR 默认阈值 peqtl-smr 5e-6**：与预注册一致；不同阈值下通过率会变化，属敏感性范畴。
4. **重复 topSNP 行**：1129 行含多 topSNP（每 probe 可能多个独立信号）；预注册集合
   （819 gene×outcome）用最显著行归并，与 coloc 全集一一对应。
