# 图表设计规格（M8）

> 生成脚本：`scripts/08_figures.R`（有数据即渲染，缺数据自动跳过）。
> 主题：CVD 安全配色（dataviz 规范），白底，细网格线，无双轴，直标标签，dpi=320。
> 配色：蓝 `#2a78d6`（主） / 橙 `#eb6834` / 青 `#1baf7a` / 灰 `#a8a8a8`（去强调） / 墨 `#0b0b0b`。

## 已渲染（数据就绪）

| 图 | 名称 | 选型 | 数据 | 状态 |
|---|---|---|---|---|
| F1 | 转录通道可比性漏斗 | 水平条形（magnitude → 顺序蓝，末档深色强调） | `results/funnel/funnel_transcript_v2.tsv` + grid CSV | ✔ |
| F2 | v1→v2 方法学修正（-log10(p) 坍缩） | dumbbell（before→after：灰=旧/蓝=新） | `compare_transcript_v1v2.csv`（4 完整对比对） | ✔ |
| F3 | 极端 p 值消失 | 2 柱对比（顺序蓝） | v1/v2 CSV 计数 | ✔ |

**F1 语义说明**：漏斗为"基因生存 → 对展开"混合结构，末档用唯一基因×结局对（32），非方法级行数（避免漏斗不递减的假象）；12→36 是基因→对的单位展开，纵轴已标注"基因 / 对数"。

**F2 关键读数**：CETP×T2D/CAD/FBG 的 v1 p 值（≤1e-30 或 0）全部坍缩到 v2 不显著（p=0.083/0.205/0.025）；HMGCR×CAD v1 p=0.93 → v2 p=0.137（均 NS）。直观呈现"未 clump 的连锁 SNP 当独立工具 → 假性极端 p"的方法学修正。

## 待数据图（骨架已写入脚本，数据到位自动渲染）

| 图 | 名称 | 选型 | 输入（待） | 设计要点 |
|---|---|---|---|---|
| F4 | 四态分类堆叠条形（按分泌状态分层） | 水平堆叠条形（part-to-whole） | `results/fourstate/*.csv` | 四色 categorical（蓝/橙/青/灰）固定顺序，discordant 单列 |
| F5 | coloc PP.H4 三档命中数 | 柱状（顺序蓝） | `results/coloc/*.csv` | 0.5/0.7/0.9 三档并列，直标 |
| F6 | LOOCV AUR 曲线 + 逐个留出散点 | 强调蓝 + 灰（单曲线是重点） | `results/controls/loocv_aur.csv` | 曲线+留出点，AUR≤0.5 仅描述性（诚实降级） |
| F7 | 外部复现方向一致率森林图 | 森林图（蓝点 + 二项 95% CI） | `results/replication/*.csv` | 经典 MR 森林图；0.5 随机线 |
| F8 | 干预介质层优先级热图 | grid heatmap（顺序蓝 无→深蓝） | `results/tables/priority_matrix.csv` | row=基因 × col=结局，值为介质层 |

## 图表纪律（对齐 README §7）
- 所有数字从 `results/` 表格直接引用，不手写。
- 幅度不跨基因/通道排序（审稿意见 3）→ 无 |Δ| 图。
- 阳性定义不依赖命中数 → 无"网格命中数主卖点"图。
- 漏斗/四态/AUR/复现全部如实呈现，负结果不隐藏。
