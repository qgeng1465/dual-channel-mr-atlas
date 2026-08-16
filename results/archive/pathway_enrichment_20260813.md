# 优先基因功能富集（2026-08-13，探索性，阴性结果）

> 76 优先基因（MR+coloc+HEIDI+Steiger 全通过）；主结果 = 56 编码集，76 全集为透明对照。
> g:Profiler g_SCS 校正；方法 + 口径排雷见 `scripts/M17_pathway_enrichment.py`。

**56 编码 g_SCS<0.05 条目: 0 | 76 全集: 0 | 56 编码最佳 g_SCS p = 0.35**

## 结论（如实）
**优先基因在 GO:BP/MF/KEGG/Reactome 通路层面无显著富集（g_SCS 校正后 0 条，最佳 p≈0.35）。**
这与 56 个编码优先基因的功能异质性一致（多为管家/普适表达基因：CHD4 染色质重塑、EIF2B2、
BLOC1S2、FASTKD5、C6orf106 等），未见共享通路信号。该阴性结果本身是诚实的生物学信息：
cis-eQTL 驱动的 T2D/CAD/FBG 优先位点基因分散在多种基本生物学过程，无单一路径富集。
注意：不排除选择偏倚掩盖真信号，或通路覆盖不完整；此结果不作任何因果主张。

## 透明对照：56 编码集名义（g_SCS）p 最小 10 条（均不显著）

| 源 | 通路 | p | 命中基因 | 驱动基因 |
|---|---|---|---|---|
| GO:BP | regulation of nucleobase-containing compound metabolic process | 0.351 | 20/48 | EIF2B2;FAM184A;GIT1;HSD17B12;KDM5A;LIPA;LRRC41;MLH3;NRBF2;NTAN1;PGAP3;RC3H2;SAP130;SNX16;SPDYE2;THAP8;YTHDF2;ZBTB6;ZNF268;ZNF34 |
| GO:MF | transcription coactivator activity | 0.36 | 5/50 | EIF2B2;FAM184A;HSD17B12;PGAP3;RC3H2 |
| GO:BP | acylglycerol metabolic process | 0.369 | 4/48 | CDC123;DGKQ;EIF2B2;SAP130 |
| GO:BP | neutral lipid metabolic process | 0.38 | 4/48 | CDC123;DGKQ;EIF2B2;SAP130 |
| GO:BP | regulation of RNA metabolic process | 0.385 | 19/48 | EIF2B2;FAM184A;GIT1;HSD17B12;KDM5A;LIPA;LRRC41;MLH3;NRBF2;PGAP3;RC3H2;SAP130;SNX16;SPDYE2;THAP8;YTHDF2;ZBTB6;ZNF268;ZNF34 |
| GO:BP | T cell receptor signaling pathway | 0.515 | 4/48 | CDC123;SIK2;TMEM136;ZBTB6 |
| GO:BP | regulation of primary metabolic process | 0.728 | 23/48 | DGKQ;EIF2B2;FAM184A;GIT1;HERPUD2;HSD17B12;KDM5A;LIPA;LRRC41;MLH3;NRBF2;NTAN1;PGAP3;RC3H2;SAP130;SIK2;SNX16;SPDYE2;THAP8;YTHDF2;ZBTB6;ZNF268;ZNF34 |
| KEGG | Glycerolipid metabolism | 0.81 | 2/23 | EIF2B2;SAP130 |
| REAC | Effects of PIP2 hydrolysis | 0.926 | 2/32 | DGKQ;SAP130 |
| GO:BP | positive regulation of RNA metabolic process | 0.952 | 12/48 | EIF2B2;FAM184A;HSD17B12;KDM5A;LIPA;LRRC41;NRBF2;PGAP3;RC3H2;SNX16;YTHDF2;ZBTB6 |

## 诚实边界
- **选择偏倚**：76 基因经 MR+coloc+HEIDI+Steiger 筛选，富集反映选择过程而非随机基因组基线；
  仅作探索性/假设生成。
- **小样本**：56 个编码基因，条目多由 1-2 基因驱动（见命中基因数）。
- **非编码基因排除**：20 个非编码/假基因（RP11-*/CTD-*/SERBP1P3）无注释通路，不参与主富集。
- **口径排雷**：`significant` 字段在 g:Profiler 响应中含义非 g_SCS<0.05，已弃用，判据为 p_value<0.05。