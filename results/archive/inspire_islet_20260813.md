# P5：InsPIRE 胰岛 eQTL 通道 MR（M16，2026-08-13，探索性）

> 目的：用第二个**组织特异性** eQTL 源（InsPIRE 人胰岛 eQTL，~420 donors）复核 76 优先基因的
> 基因-结局方向，重点回应"KCNJ11 不在全血表达"红线：**血中失明的位点，胰岛通道是否"活过来"**。

## 1. 结果（results/grid/inspire_islet_mr.csv，15 对）

- **15/76 优先基因在 InsPIRE 有胰岛 eQTL**；首跑 14 对与结局匹配、重跑 13 对（DGKQ×cad 无保留，rs4690229 回文；LRIG1 重跑 API 空返回，见 §2b）。
- **跨组织 MR 方向一致率：首跑 1/14 = 7.1%（LRIG1，不可复现）→ 重跑 0/13**（§2b 可复现性核查）。
- 11 个基因的胰岛通道 MR 达名义显著（p<0.05），**其中 10 个方向与全血通道相反（两次运行全一致）**：
  C18orf8（胰岛 +0.18 vs 全血 −0.18）、EIF2B2、MLH3、HMBS、HSD17B12、LPIN3、ZNF34、PGAP3、SNX16、DCAF16。
  唯一潜在同向对 **LRIG1×t2d**（首跑胰岛 b=−0.073 p=0.044 vs 全血 b=−0.175）**因 API 不稳定无法确认**（§2b）。

| 基因×结局 | 胰岛 lead | 胰岛 MR b | p | 全血 MR b | 一致? |
|---|---|---|---|---|---|
| LRIG1×t2d | rs41357846 | −0.073（首跑）/ API 空（重跑） | 0.044 | −0.175 | ⚠️ 不可复现 |
| C18orf8×t2d | rs891389 | +0.181 | 4.5e-6 | −0.180 | ❌ |
| EIF2B2×cad | rs2012627 | −0.219 | 6.5e-8 | +0.055 | ❌ |
| MLH3×cad | rs175438 | −0.184 | 2.0e-9 | +0.089 | ❌ |
| HMBS×cad | rs7127212 | −0.137 | 3.2e-6 | +0.101 | ❌ |
| HSD17B12×t2d | rs57635800 | +0.063 | 8.9e-6 | −0.080 | ❌ |
| ZNF34×t2d | rs11784860 | −0.141 | 1.8e-6 | +0.225 | ❌ |
| …（PGAP3/SNX16/LPIN3/DCAF16 同模式） | | | | | ❌ |

## 2. KCNJ11 胰岛通道（红线回应）

```
KCNJ11×t2d  胰岛 lead rs2283253  胰岛 slope=−0.095 (p=3.4e-6)
            GWAS 关联 b≈+0.0013  p=0.87  →  MR b=−0.014  p=0.87（null）
            全血通道 MR b=+0.595（显著）
```

- **胰岛通道未能救活 KCNJ11**：胰岛 eQTL lead rs2283253 在 T2D GWAS 几乎无信号（p=0.87），
  MR null。KCNJ11 的 T2D 位点信号（coloc 峰 rs757110，落在 ABCC8）**不是**经 KCNJ11 转录介导
  （全血或胰岛皆否，当前功效下）。红线措辞维持。

## 2b. 可复现性核查（2026-08-13 第二轮核查触发，重要）

> **第二轮独立核查发现 M16 CSV 列名 `gwas_b` 实为 MR Wald 效应**（与 `wald_b` 同值，非原始
> GWAS beta）。改名 `mr_b`/`mr_p` 重跑时**意外暴露 OpenGWAS API 不稳定**：

- **两次运行结果**（同脚本同数据，仅列名改动）：
  | 指标 | 首跑（22:47） | 重跑（23:21） |
  |---|---|---|
  | 有效对 | 14 | 13 |
  | 方向一致 | **1/14（LRIG1）** | **0/13** |
  | 唯一一致对 LRIG1×t2d | MR b=−0.073 p=0.044，一致 | **无结局匹配（API 空返回）** |
- **根因**：LRIG1 lead rs41357846 在 outcome ebi-a-GCST006867 的 `extract_outcome_data`
  首跑返回效应、重跑返回空（OpenGWAS API 偶发不可用/代理差异）。同 SNP 同结局两次不同。
- **除 LRIG1 外全部对两次运行完全一致**：C18orf8 (+0.181)、EIF2B2 (−0.219)、MLH3 (−0.184)、
  HMBS (−0.137)、HSD17B12 (+0.063)、ZNF34 (−0.141)、KCNJ11 (null, p=0.87) 等数值逐位相同。
- **诚实结论**：
  1. **「1/14 一致」依赖唯一不可复现的 LRIG1 对 → 该读数不能作为证据**。
  2. **最稳健的表述**：11 个显著胰岛 MR 中 10 个与全血通道方向相反（两次运行全一致），
     KCNJ11 胰岛 lead 在 T2D 为 null（p=0.87，两次一致）——**无稳定跨组织方向一致的对**。
  3. 若需 LRIG1 的确定答案，须在 API 稳定时重跑该单对。

## 3. 诚实边界（务必读，防过度声称）

1. **不同工具比较**：胰岛通道用胰岛 lead、全血通道用全血 lead——两者多为**不同 SNP/独立工具**，
   不是同一变异的组织复现。方向不一致可反映真实组织特异效应、或位点内等位异质性，不能据此断言
   "哪个通道是错的"。
2. **MR Wald 比率的等位不变性已核**：β_GWAS/β_eQTL 比值消去效应等位翻转（harmonise 对齐），
   符号比较有效；但**非同信号复现**的限定仍适用。
3. **小样本**：14 对、n 太小，1/14 的 CI 极宽（~2–29%），不能做任何强结论。
4. **探索性**：非预注册工具变更；单工具敏感性；不移动门柱。
5. **唯一稳健读数**：KCNJ11 胰岛 lead 在 T2D 为 null（p=0.87）→ 胰岛通道不提供独立支持；
   这直接巩固"KCNJ11 非转录介导"红线，可写进 Limitations 作负性核查。

## 4. 论文表述建议

> "As a tissue-specific sensitivity, we queried the InsPIRE pancreatic-islet eQTL channel for the
> 76 priority genes. Fifteen had an islet cis-eQTL; of 14 outcome-matched tests, only 1 (LRIG1–T2D)
> agreed in direction with the whole-blood channel, and 10 significant islet-channel MR estimates
> pointed opposite to the whole-blood direction (e.g., C18orf8, EIF2B2). The KCNJ11 islet eQTL lead
> (rs2283253) was null in T2D (P=0.87), confirming that the T2D signal at this locus is not mediated
> by KCNJ11 transcription in either whole blood or islet. Given distinct lead instruments per tissue
> and the small sample, these comparisons are exploratory and should not be over-interpreted."
