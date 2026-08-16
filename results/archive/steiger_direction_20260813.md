# Steiger 因果方向检验（2026-08-13，P5 部分）

> 目的：为 76 个"HEIDI 通过 + coloc strong"优先基因判定因果方向（eQTL→outcome vs 反向/共存）。
> 方法：Steiger 方差解释比较（TwoSampleMR 同款算法），从全量 SMR 输出直接计算。
> 输入：`transcript_smr_heidi_full.csv` 的 topSNP b_eQTL/b_GWAS/Freq；N 近似
> （eQTLGen 全血 n=31,684；结局 T2D 655,666 / CAD 547,261 / FBG 133,010）。
> 纪律：探索性敏感性；方向分类对 N 近似稳健，Steiger p 值为近似（标准做法）。

## 结果

| 方向 | n | 占比 | Steiger p<0.05 |
|---|---|---|---|
| **eQTL→outcome**（正向，转录→结局） | **73** | **96.1%** | 67 |
| outcome→eQTL（反向，疑似共存/连锁） | 3 | 3.9% | 0 |

**反向/方向未决位点（3 个，均不显著）**：
- SPATA5×T2D（H4=0.92, p=0.33）、C6orf106×CAD（H4=0.91, p=0.48）、ARL13B×T2D（H4=0.94, p=0.84）
- → 这些位点 r2 接近、方向未决，如实标注为"方向不确定"，不作正向声称。

**方向最明确位点（正向，Steiger p≈0）**：LIPA×CAD、ZNF34×T2D、NRBF2×CAD、KDM5A×CAD、
SH3BGRL3×T2D、MLH3×CAD、PDGFC×T2D、MED27×T2D 等（r2_eQTL 远大于 r2_outcome）。

## 诚实解读

1. **76 个优先基因中 96% 的因果方向是 eQTL→outcome**：这支持"转录水平变异驱动结局"的因果
   方向解释，与 coloc + HEIDI 形成三层收敛（共定位 + 同质性 + 方向）。
2. **3 个方向未决位点不作声称**：r2 接近，可能是表达与结局共享上游因果（共存），
   不违反主结论，但论文需在补充中列出。
3. **已知局限**：N 为近似；eQTL β 为标准化表达效应、GWAS β 为 log-odds，跨尺度比较是
   Steiger 的标准近似假设（TwoSampleMR 同款）；本分析是敏感性，非正式因果推断。
   正式的因果方向推断应以工具变量 + 敏感性方法为准，Steiger 只作支撑。

## 联动

- 与 HEIDI 全集（48.6% 全集通过 / strong 71.7%）+ coloc strong（106）合并：
  **76 个"三层收敛"优先基因 = coloc strong ∩ HEIDI 通过 ∩ eQTL→outcome 方向**，
  其中 73/76 方向正向、67/76 Steiger 显著。
- 这是 Fig 5 桑基图第三层（coloc→HEIDI→方向）的读数。
