# 数据溯源表 (Data Provenance)

> 学术诚信要求：每个数据源记录 ID/URL/下载日期/样本量。本表持续更新。

## 暴露层（通道 1：细胞转录本）

| 数据 | 版本 | 样本量 | 来源 URL | 下载日期 | 本地文件 | 完整性 |
|---|---|---|---|---|---|---|
| eQTLGen 全血 cis-eQTL（显著） | Võsa 2021 | n=31,684 | https://download.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz | 2026-08-06 | data/eqtlgen/cis-eQTL-significant.txt.gz | ✅ GZIP_OK, 10.5M 行, 16,923 基因 |
| eQTLGen SNP 频率（等位基因定向） | Võsa 2021 | — | https://download.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz | 2026-08-06 | data/eqtlgen/SNP_AF.txt.gz | ⏳ 240MB 循环续传中（服务器带宽受限） |

## 暴露层（通道 2：循环蛋白）

| 数据 | 状态 | 样本量 | 来源 | 本地文件 | 完整性 |
|---|---|---|---|---|---|
| deCODE 血浆 pQTL（4907 探针/4719 蛋白，HG38，n=35,559） | ✅ 2026-08-06 已获表单链接并破解下载 API | n=35,559 | https://www.decode.com/summarydata/ → 邮件链接；下载机制：source map 逆向 `GET /s3/fileInfo` / `/s3/folder` / `/s3/download`（见 README §4.2） | data/decode/ | readme ✅；6 个对照蛋白文件（PCSK9/HMGCR/ANGPTL3/APOC3/APOB/LDLR，各 ~950MB）循环续传中；assocvariants.annotated/excluded 下载中 |
| INTERVAL SomaScan pQTL（OpenGWAS prot-a-*） | ✅ JWT 可用 | n=3,301 | OpenGWAS API（JWT） | 经 API 拉取 | ✅ 对照可用性：HMGCR=prot-a-1354、ANGPTL3=prot-a-98；PCSK9/CETP/NPC1L1/APOC3 未收录 |

## 结局层（OpenGWAS，JWT 已配置，样本量已实查）

| 结局 | ID | 样本量（gwasinfo 实查） | 用途 |
|---|---|---|---|
| T2D（Mahajan/meta） | ebi-a-GCST006867 | n=655,666 | 主结局（logOR） |
| CAD | ebi-a-GCST005194 | n=296,525 | 主结局（logOR） |
| FBG（MAGIC Manning） | ebi-a-GCST005186 | n=58,074 | 血糖轴（beta） |

## 结局层（GWAS Catalog 免 token 源，备用/交叉验证）

| 数据 | 来源 | 样本量 | URL | 状态 |
|---|---|---|---|---|
| T2D meta (Xue 2018) | GWAS Catalog FTP | n≈568,947 | https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006867/Xue_et_al_T2D_META_Nat_Commun_2018.gz | ⏳ 下载中 |
| MAGIC FG (Manning) | GWAS Catalog FTP | n≈58,074 | https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST005001-GCST006000/GCST005186/MAGIC_Manning_et_al_FastingGlucose_MainEffect.txt.gz | 待下载 |

## OpenGWAS（标准结局源）

- ✅ JWT 已配置到 `~/.Renviron`（`OPENGWAS_JWT=...`），`gwasinfo()`/`extract_outcome_data()` 实测可用；结局 ID 与样本量已实查回填（上表）。

## 运行状态（2026-08-06）

| 数据 | 状态 | 更新说明 |
|---|---|---|
| eQTLGen SNP 频率（2018-07-18_SNP_AF_for_AlleleB...） | ⏳ 下载中 | 240 MB；URL https://download.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz ；列 SNP/hg19_chr/hg19_pos/AlleleA/AlleleB/.../AlleleB_all；effect=AlleleA → eaf=1−AlleleB_all |
| INTERVAL prot-a-1354 = HMGCR；prot-a-98 = ANGPTL3 | ✅ 已核验 | n=3,301；cis p<5e-6 内无独立工具（功效不足，诚实记录） |
| deCODE 4,907 蛋白 | ⏳ 下载中 | 官网唯一（https://download.decode.is，API: /s3/fileInfo /s3/folder /s3/download）；对照定位 PCSK9=5231_79、HMGCR=5230_99、ANGPTL3=10391_1、APOC3=6461_54（各 ~910 MB）；CETP/NPC1L1 未收录；格式 hg38 per-SD Beta + assocvariants.annotated 补 effectAlleleFreq；速度瓶颈 ~5–8 KB/s |
| 结局 OpenGWAS | ✅ JWT 已验证 | T2D ebi-a-GCST006867 / CAD ebi-a-GCST005194 / FBG ebi-a-GCST005186 全部 API 解析成功 |
