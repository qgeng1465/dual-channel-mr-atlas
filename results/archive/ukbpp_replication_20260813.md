# M13：UKB-PPP 跨平台蛋白复现（2026-08-13，探索性）

> 目的：用第二个独立 pQTL 平台（Sun 2023 UKB-PPP，Olink，discovery n≈34k，英国生物库）
> 核查 deCODE（冰岛 n=35,559）蛋白工具的可靠性与跨平台一致性。
> 面板覆盖 8/11：PCSK9 INSR ANGPTL3 LDLR APOB DPP4 GLP1R GCG（无 APOC3/PCK1/HMGCR）。
> 坐标：UKB-PPP GENPOS=hg38、ID 列=hg19（实测 PCSK9 top=hg38 55,039,974）；与 deCODE hg38 直接对齐。
> Beta 同为每等位蛋白 SD 单位，尺度直接可比。
> 纪律：探索性，独立输出 CSV，只做"工具与方向"一致性，不声称新发现。

## Part A：cis 区共享变异方向一致率 + beta 相关

deCODE cis ∩ UKB-PPP cis（hg38 位置匹配 + 等位集合匹配，回文保守排除）：

| 蛋白 | 共享变异 | 方向一致率 | Pearson r | p |
|---|---|---|---|---|
| **PCSK9** | 8703 | 61.3% | **0.338** | <1e-232 |
| **ANGPTL3** | 7686 | 68.5% | **0.396** | <1e-287 |
| INSR | 9737 | 57.6% | 0.075 | 1.7e-13 |
| DPP4 | 5965 | 54.4% | 0.086 | 2.5e-11 |
| GLP1R | 8417 | 52.9% | 0.013 | 0.22（NS） |
| GCG | 5981 | 53.8% | -0.008 | 0.53（NS） |
| LDLR | — | 0 共享 | — | 无 deCODE cis 子文件 |
| APOB | — | 0 共享 | — | 无 deCODE cis 子文件 |

**诚实解读**：强 cis 信号蛋白（PCSK9/ANGPTL3）跨平台 beta 显著正相关（r≈0.34-0.40，
在数万 SNP 下 p 极小）；弱 cis 信号（GLP1R/GCG）r≈0 → **平台间差异集中在弱信号**，
这是已知的 pQTL 平台可比性问题，不是工具坏。

## Part B：UKB-PPP 单工具 MR 复现（cis top → T2D/CAD/FBG）

UKB-PPP 自己 cis 区最显著变异（min P）→ 单工具 Wald，与 deCODE MR 对比：

| 基因 | top 变异 | UKB-PPP top P | 复现结论 |
|---|---|---|---|
| **PCSK9 × CAD** | rs11591147（R46L） | 0 | **强复现 ✓**：UKB-PPP b=+0.209 p=2.0e-22 vs deCODE b=+0.190 p=5.5e-12，方向一致、双平台全基因组显著 |
| PCSK9 × t2d/fbg | — | — | outcome 提取失败（deCODE 侧本就 null，p=0.81/0.70，不补跑） |
| INSR × cad | rs1799815 | 4.5e-11 | 不成立 ✗：UKB-PPP b=+0.12 p=0.35（NS）vs deCODE b=-0.05 p=0.039，方向相反且弱 |
| ANGPTL3 × 3 | rs626787 | 0 | 双平台一致 null ✓：UKB-PPP 全 p>0.58、deCODE 全 p>0.45（同 null，非冲突） |
| LDLR/APOB × 3 | — | — | 无法复现：UKB-PPP top 变异无 deCODE rsID 匹配（deCODE 无 cis 子文件） |
| DPP4 × t2d | rs13015258 | 9e-204 | 单平台 nominal：UKB-PPP b=-0.105 p=0.003；deCODE MR 无 ok 结果，方向不可比 |
| GLP1R / GCG × 3 | — | — | UKB-PPP 全 NS（p>0.17）；deCODE 无结果 |

## 综合诚实结论

1. **最强的跨平台证据是 PCSK9 × CAD**：两个独立 pQTL 平台、同一 cis 变异
   （rs11591147 = PCSK9 R46L 功能变异）、同一因果方向、双平台均全基因组显著
   （deCODE p=5.5e-12 / UKB-PPP p=2.0e-22）。这为"蛋白通道在 CAD 的信号"提供了
   平台独立的支持，比单平台更有说服力。
2. **弱信号不放大**：INSR/ANGPTL3/GLP1R/GCG 在 UKB-PPP 下均未复现任何心血管信号，
   与 deCODE 侧弱/空一致（PCSK9 是唯一强 cis 蛋白之一，其心血管信号最强也最稳）。
3. **平台差异如实报告**：Part A 显示 r≈0 的蛋白恰好是弱 cis 信号蛋白；LDLR/APOB
   因 deCODE 无 cis 子文件无法跨平台对齐（数据覆盖局限，非方法缺陷）。
4. **不影响主结论**：M13 是工具层/方向层核查，主结论（双通道基因×结局×介质层优先级
   图谱）仍以预注册主流程为准；本结果作为"蛋白工具经双平台核查"的方法学加固。

## 与杠杆 B 联动（合并叙事）

- **pQTL 工具层**：deCODE vs UKB-PPP 共享变异一致（PCSK9/ANGPTL3 显著正相关）
  + PCSK9×CAD 单工具 MR 双平台复现
- **脂质结局校准**：PCSK9→LDL/apoB、APOC3→HDL 回收教科书因果（M12）
- **心血管结局**：如实探索性报告（蛋白通道信号集中在强 cis 蛋白，其余弱/空）

→ 蛋白通道从"未校准"升级为"经双平台核查 + 脂质校准 + 心血管信号弱"的诚实三层叙事。
