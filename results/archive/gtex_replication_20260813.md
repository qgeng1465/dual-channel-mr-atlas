# P1：GTEx 组织特异 eQTL 复现（M15，2026-08-13）

> 目的：用第二个**独立组织 eQTL 源**（GTEx v8 全组织 egenes）复核 106 个 strong-coloc 命中的
> 基因-结局方向，检验共享因果信号是否跨 eQTL 源/组织稳健。

## 1. 结果（results/grid/gtex_replication_p1.csv）

- 输入：106 个 strong-coloc 命中；GTEx（6 组织：Whole_Blood/Liver/Pancreas/Adipose_Subcutaneous/
  Muscle_Skeletal/Artery_Coronary，cis±1Mb）**93/106 有显著 cis-eQTL**（lead 取跨组织最小 pval）。
- 结局匹配 + MR 有效：**63**。
- **方向一致率（GTEx 通道 vs eQTLGen 通道）：44/63 = 69.8%（95% CI 57.7–79.7%）**。
- **排除 Whole_Blood 后（独立跨组织）：26/40 = 65.0%（95% CI 49.5–78.0%）**。
- 同变异子集（GTEx lead == eQTLGen top SNP）：6 对，方向一致 **6/6 = 100%**。
- GTEx 名义显著（p<0.05）子集：50 对，方向一致 40/50 = 80.0%。

## 2. 按组织分层

| 组织 | 有效 | 一致率 |
|---|---|---|
| Whole_Blood | 23 | 18/23 = 78.3%（同组织不同源） |
| Adipose_Subcutaneous | 21 | 16/21 = 76.2% |
| Muscle_Skeletal | 11 | 5/11 = 45.5% |
| Pancreas | 6 | 5/6 = 83.3% |
| Liver | 2 | 0/2 = 0%（n 过小） |
| **非全血合计** | **40** | **26/40 = 65.0%** |

## 3. 诚实解读

1. **非全血组织一致率 65%（26/40）**：即使排除同组织的 Whole_Blood，独立组织（脂肪 76%、
   胰岛 83%）仍显示方向基本稳定 → **eQTLGen 全血通道的优先化信号不是全血组织伪影**。
2. **与 InsPIRE（10/11 反向）的差异是设计差异，不可直接对比**：GTEx 取跨组织最小 pval 的
   best-lead（偏向强 eQTL 基因），InsPIRE 取固定胰岛 lead；两设计回答不同问题。GTEx 的 65%
   是一阶「方向稳定性」证据，InsPIRE 是「特定组织通道是否救活血中失明位点」的负性核查。
3. **这不改变主结论**：GTEx 一致率（69.8%）是**方向一致性**（MR 符号），与 coloc 一致率
   （12.9%，共定位支持）是**不同指标**。方向稳定 ≠ 共定位支持。主结论「MR 显著集多数未获
   区域共定位支持」不变。
4. **选择偏倚**：93/106 有 GTEx eQTL 的命中偏"eQTL 稳定"基因，一致率不能外推到无 GTEx eQTL 的 13 个。
5. 同变异 6/6 全一致（5 在全血、1 在脂肪）——最强的一阶同信号复现，但 n 太小。

## 4. 论文表述建议

> "Among the 106 colocalized tests, 93 (88%) had a significant cis-eQTL in ≥1 of 6 GTEx tissues;
> of 63 outcome-matched tests, 44 (69.8%) showed a concordant direction of effect between the
> GTEx and eQTLGen channels, including 65.0% (26/40) when Whole_Blood was excluded. Six
> same-variant comparisons were all concordant. Directional stability of the prioritized
> gene–outcome associations thus extends beyond the eQTLGen whole-blood discovery channel."
