# 仅编码版一致性率（M5 剩余项，2026-08-13）

> 目的：检验"13% 一致性率"是否被非编码/伪影基因稀释。结论：**不是**。

## 1. 方法

- 用 Ensembl REST（走本机代理）给 `transcript_coloc.csv` 全部 759 个唯一 ENSG 注释 biotype
  （保存于 `results/grid/ensg_biotype_20260813.csv`）。
- 一致率按 gene×outcome 对计算，分编码 / 非编码分层。

## 2. 结果

| 分层 | MR 显著对 | strong coloc 对 | 一致率 | 95% CI (Wilson) |
|---|---|---|---|---|
| **全量** | 819 | 106 | **12.94%** | 10.82–15.42% |
| **蛋白编码** | 637 | 76 | **11.93%** | 9.55–14.66% |
| 非编码/未注释 | 182 | 30 | 16.48% | 11.64–22.44% |

- 编码基因 biotype 分布：protein_coding 585、lncRNA 85、NOT_FOUND 34（Ensembl 查不到 = 版本过旧）、
  pseudogene 各型 43、miRNA/snRNA/snoRNA/misc_RNA 11、IG_V 1。

## 2b. 按结局分层（M7 补充）

| 结局 | MR 显著对 | strong 对 | 一致率 | Wilson 95% CI |
|---|---|---|---|---|
| **T2D** | 331 | 58 | **17.5%** | 13.8–22.0% |
| **CAD** | 477 | 46 | **9.6%** | 7.3–12.6% |
| FBG | 11 | 2 | 18.2% | 5.1–47.7%（n 太小） |
| 全量 | 819 | 106 | 12.9% | 10.8–15.4% |

- **T2D 共定位率（17.5%）显著高于 CAD（9.6%）**：非重叠的 Wilson CI → 结局间真实差异。
  可能反映 eQTLGen 全血 eQTL 对 T2D 胰岛/代谢信号的相关性高于 CAD 动脉信号，或 CAD 区域更多
  multiple causal variants（coloc 单因果假设被违反 → PP.H4 低）。这是可写入论文的分层事实。
- FBG n=11（历史 coloc 范围所限），不单独解读。

## 3. 诚实解读

1. **编码版一致率（11.9%）≈ 全量（12.9%）**：头部读数不依赖非编码基因；
   "MR 显著集约 1/8 获共定位支持"对纯编码基因同样成立。
2. **非编码层一致率反而更高（16.5%）**：非编码元素（lncRNA/伪基因）并非"灌水"——
   它们在 MR 显著集里获得 strong coloc 支持的比例不低。这与"转录图谱"性质一致
   （非编码 eQTL 是真实中介），但仍是**优先化/探索性**主张，不作因果声称。
3. **34 个 NOT_FOUND**：Ensembl 现行版本查不到（旧 ENSG 版本号），归入非编码层；
   对结论无实质影响（它们 MR-sig 对数 <20）。
4. **分层一致性**：无论编码/非编码，一致率都在 12–17%，远低于 HEIDI-pass 的 71.7%——
   主结论"MR 显著集多数未获区域共定位支持"稳健。

## 4. 论文表述建议

> "Among 819 cis-MR-significant gene–outcome tests, 106 (12.9%) achieved regional colocalization
> (PP.H4≥0.8). The concordance was similar when restricted to protein-coding tests
> (76/637 = 11.9%), indicating the low concordance is not driven by non-coding elements;
> non-coding loci showed if anything a higher rate (30/182 = 16.5%)."
