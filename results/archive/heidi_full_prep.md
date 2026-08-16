# HEIDI 全量输入准备报告（P3：818/819 转录 MR 显著测试）

日期：2026-08-13　作者：qgeng1465　状态：输入已构建完成（未跑 SMR 二进制）

## 1. 背景与目标

转录通道 cis-MR（eQTLGen 全血 eQTL exposure × T2D/CAD/FBG GWAS outcome）此前已跑
107 探针（106 个 strong coloc 命中 + KCNJ11）的 SMR/HEIDI。本任务把输入扩展到
`results/grid/transcript_coloc.csv` 的全部 **819 个 MR 显著测试（759 个唯一基因）**，
为后续经资源仲裁全量跑 HEIDI 做准备。本任务只读 + 构建输入，未运行 SMR 二进制。

## 2. 数据可用性核查（对应任务点 1、2）

### 2.1 目标基因集
- `results/grid/transcript_coloc.csv`：**819 数据行**（19 列），**759 个唯一 ENSG 基因**，
  含 symbol。全部为 `ENSG` 前缀，无缺失、无 NA。

### 2.2 eQTLGen 数据可用性（能建 .esd 的基因数）
本地 eQTLGen 全血 cis-eQTL 源文件齐全：
- `data/eqtlgen/cis-eQTL-significant.txt.gz`（hg19，322 MB）— 显著 cis-eQTL 全量
- `data/eqtlgen/SNP_AF.txt.gz`（hg19，240 MB）— 频率
- LD 参考 `data/ldref/1kg.v3/EUR.{bim,fam,bed}` — 1kg EUR（hg19）

流式核查：**759/759 个基因在 eQTLGen cis-eQTL-significant 中至少有一个显著 SNP**。
抽取得到 **591,093 个 SNP-gene 对**（759 基因）；合并 AF + 1kg EUR 参考面板、
剔除 eaf 缺失/无效、Zscore 缺失后存活 **527,000 行 / 759 基因**。

**结论：全部 759 个基因都能建 .esd（100% 可用）。缺 0 个基因。**

### 2.3 GWAS outcome 缓存可用性（.ma 输入）
`.ma` 复用 `results/grid/_coloc_gwas/{OUTCOME}_{ENSG}.rds` 缓存（与 coloc/SMR 完全同源）：

| outcome | GCST ID | 有缓存的基因数 / 759 | .ma SNP 数 |
|---------|---------|---------------------|-----------|
| t2d | GCST006867 | 331 | 605,005 |
| cad | GCST005194 | 477 | 1,261,496 |
| fbg | GCST005186 | 11  | 9,909 |

每个 coloc 行（gene-outcome 对）都有对应缓存（819 个 .rds），故 331+477+11=819，
与 MR 显著测试行数一致。FBG 只有 11 个基因是历史 coloc 范围使然（与已跑 107 探针相同）。

## 3. 已构建的输入文件

| 文件 | 内容 | 行数 |
|------|------|------|
| `data/smr/trans_esd_full/{ENSG}.esd` | 759 个探针文件（每基因一个） | 759 个文件 |
| `data/smr/trans_flist_full.txt` | SMR `--eqtl-flist` 输入 | 759 数据行 + 1 表头 |
| `data/smr/trans_t2d_full.ma` | 结局 t2d GWAS 汇总 | 605,005 SNP |
| `data/smr/trans_cad_full.ma` | 结局 cad GWAS 汇总 | 1,261,496 SNP |
| `data/smr/trans_fbg_full.ma` | 结局 fbg GWAS 汇总 | 9,909 SNP |

- `trans_flist_full.txt` 列格式与现 `trans_flist.txt` 完全一致：
  `Chr ProbeID GeneticDistance ProbeBp Gene Orientation PathOfEsd`（空格分隔，含表头）。
- 全部 759 行 `PathOfEsd` 均指向 `trans_esd_full/`，且文件逐一验证存在（759/759）。
- 样本抽查：基因 1（PLXND1）flist 的 Chr=3、ProbeBp=129030642 与 .esd 内
  `Chr=3`、`round(mean(range(Bp)))=129030642` 一致；.esd 无 NA/NaN。

### 输入构建约定（与 M10 完全一致）
- 坐标：全程 hg19。eQTLGen cis-eQTL `SNPChr/SNPPos` 与 1kg bim 逐位点一致（M10 已验证），
  .esd 直接用文件内坐标，**无需 liftOver**。
- 效应量：eQTLGen 只有 Zscore → `Beta = Z/sqrt(NrSamples)`，`se = 1/sqrt(NrSamples)`。
- 频率：`eaf` 按 `AssessedAllele` 取向由 SNP_AF `AlleleB_all` 计算；缺失剔除。
- LD 参考：1kg EUR，.esd 只保留参考面板内变异（SMR 需参考算 HEIDI LD）。
- `.esd` 列：`Chr SNP Bp A1 A2 Freq Beta se p`（空格分隔，含表头）。
- flist 的 `ProbeBp = round(mean(range(SNPPos)))`（cis-SNP 范围中点，M10 同约定，非基因 TSS），
  `Chr` = 该基因首行 SNP 的 `SNPChr`，`Orientation = "+"`，`Gene` 列 = 基因符号（来自 GeneSymbol），
  `ProbeID` = ENSG。

## 4. 脚本与用法

`scripts/M10b_build_heidi_full_inputs.R` 已写好。单趟流式（awk 按列 8=Gene 匹配）抽取
显著 cis-eQTL → 合并 AF + 1kg → 逐基因写 .esd → 写 flist_full → 复用 `_coloc_gwas` 缓存
写各结局 `_full.ma`。产物全部为**新建文件**，已确认不覆盖既有 107 探针输入
（`trans_esd/`、`trans_flist.txt`、`eqtlgen_trans.*`、`trans_{t2d,cad,fbg}.ma`，
其 mtime 保持在 2026-08-13 20:20 未变）。

运行（本任务已在资源仲裁下成功跑完一遍）：
```bash
PATH=/data/gengqiushuo/home/miniconda3/envs/r-mr/bin:$PATH \
  Rscript scripts/M10b_build_heidi_full_inputs.R
# 或经本机资源仲裁工具排队运行（命令略）
```
本脚本**不调用 SMR 二进制**。

## 5. 下一步（SMR/HEIDI 全量跑，需经资源仲裁）

BESD 构建 + HEIDI 为独立步骤（仿 `scripts/M10_run_smr.sh`，仅改输入名）：
```bash
# 1) 由 flist_full 建全量 BESD
tools/smr --eqtl-flist data/smr/trans_flist_full.txt --make-besd \
  --out data/smr/eqtlgen_trans_full --thread-num N
# 2) 各结局 SMR + HEIDI（默认含 HEIDI；--peqtl-smr 5e-6 --disable-freq-ck 同 M10）
tools/smr --bfile data/ldref/1kg.v3/EUR --gwas-summary data/smr/trans_t2d_full.ma \
  --beqtl-summary data/smr/eqtlgen_trans_full \
  --peqtl-smr 5e-6 --disable-freq-ck --thread-num 6 --out data/smr/run/tsmr_full_t2d
# cad / fbg 同理
```

## 6. 风险与坑（记录）

1. **坐标全 hg19**：eQTLGen 显著文件、SNP_AF、1kg EUR、_coloc_gwas 缓存均 hg19，一致，
   无需转换。若将来混入 hg38 数据源需先 liftOver。
2. **.ma 是 coloc 区域缓存的并集，非全 GWAS**：沿用 M10 同源策略（与 coloc 严格同源，
   cis 测试只需区域 SNP）。若审稿人要求全基因组 .ma，属额外数据获取，非本管线缺陷。
3. **ProbeBp 为 cis-SNP 范围中点**（M10 约定），非基因 TSS/探针 TSS。SMR 仅用该位置做
   probe 定位，不参与 LD 计算，影响可忽略；但报告 ProbeBp 时应注明此约定。
4. **Orientation 一律 "+"**：SMR 默认；eQTLGen 提供的是 cis-eQTL 效应（非 allele-strand
   定向），与 M10 一致。若个别基因需按链翻转需人工核对，当前不适用。
5. **频率剔除**：591,093 对中有 ~11% 因不在 1kg EUR 或无 AF 被剔除（527,000 保留）。
   全部 759 基因仍有 ≥1 存活 SNP，无基因被整体剔除；HEIDI 的 nsnp 取决于该过滤后的存活 SNP。
6. **FBG 覆盖小**：`trans_fbg_full.ma` 仅 11 基因 / 9,909 SNP（与 coloc 范围一致），
   属预期，不是数据缺失。
7. **SMR --make-besd 依赖系统 GLIBC**：tools/smr 为 2024-03 编译的二进制，当前主机可跑
   （107 探针 BESD 已成功构建过）。全量 759 探针文件更大，跑时给足 --thread-num 与内存。
8. **命名隔离**：全量产物统一 `_full` 后缀 + `trans_esd_full/` 目录，避免与已跑 107 探针
   混淆；BESD 建议输出 `eqtlgen_trans_full`。

## 7. 关键数字小结

- 能建 .esd 的基因：**759 / 759（100%）**；缺 0 个。
- `trans_flist_full.txt`：**已生成，759 数据行**（+表头），位于 `data/smr/trans_flist_full.txt`。
- `.esd` 文件：759 个，位于 `data/smr/trans_esd_full/`。
- 各结局 `_full.ma`：t2d 605,005 SNP、cad 1,261,496 SNP、fbg 9,909 SNP。
- 脚本：`scripts/M10b_build_heidi_full_inputs.R`（已写好并验证）。
- 未覆盖任何既有结果；未运行 SMR 二进制。
