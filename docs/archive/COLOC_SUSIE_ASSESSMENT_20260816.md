# coloc.susie 收敛性评估（2026-08-16，审稿意见 §3 处理记录）

## 结论（方案 B 落地）

**coloc.susie（外部 1000G EUR LD）在本数据集中不收敛，且提高迭代次数无法解决。**

- 6 个代表位点在 max_iter=200 下 **全部 conv_eqtl=FALSE / conv_gwas=FALSE**（`results/m34_coloc_susie_20260816.csv`）。
- 对最强案例 RBM6 提高至 **max_iter=1000** 实测：susie_rss 仍报 "IBSS algorithm did not converge
  in 1000 iterations! Please check consistency between summary statistics and LD matrix"
  （`docs/archive/m34b_partial_maxiter1000_20260816.csv` 记录 RBM6 FAILED susie_rss）；
  单次拟合即超出实用运行时间（>2 min/位点 × 12 次拟合）。此警告是 susieR 对外部 LD 与汇总统计
  不一致的诊断，属**原理性限制**而非迭代不足。

**处理**：coloc.susie 降级为 exploratory 敏感性，不作为主表结论；LAMC1 排除依据改用多信号证据
（非收敛后验单独不采用）。

## 数字（m34，max_iter=200，全部 conv=FALSE）

| symbol | outcome | abf_pp4 | susie_pp4 | CS_eqtl | CS_gwas | 说明 |
|---|---|---|---|---|---|---|
| RBM6 | t2d | 0.9448 | 1.0000 | 2 | 2 | 与 abf 方向一致 |
| CNNM2 | cad | 0.9366 | 0.9996 | 10 | 6 | 与 abf 方向一致 |
| PLAUR | cad | 0.9957 | 0.9990 | 10 | 9 | 与 abf 方向一致 |
| CD101 | t2d | 0.9488 | 1.0000 | 10 | 8 | 与 abf 方向一致 |
| RIC8A | cad | 0.8895 | 1.0000 | 3 | 3 | 与 abf 方向一致 |
| LAMC1 | cad | 0.9139 | **0.0000** | 8 | 3 | 多信号不一致 → 排除 |

**多信号证据（LAMC1 排除依据，独立于非收敛后验）**：
eQTL 侧 8 个可信集 vs GWAS 侧 3 个可信集；区域 max |z| = 4.32（低于多数共定位位点的单峰 z）；
两数据集无可定位的共享因果变异。coloc.abf 的 PP.H4=0.9139 依赖单因果变异假设，多信号框架下
该支持消失 → LAMC1 不作为候选。

## 落点

- F8 图：abf vs susie 成对条形 + 收敛标记（全 ✗）+ CS 数 + LAMC1 高亮；注释注明 max_iter=1000 实测无效。
- §3.7/§4：coloc.susie 为 exploratory，明确"外部 LD 下不收敛（提升迭代不解决）"。
- Table S3（m34 输出）：保留 6 基因 SuSiE 诊断。
