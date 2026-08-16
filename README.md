# Transcriptome-wide cis-MR × colocalization atlas for type 2 diabetes, coronary artery disease, and fasting glucose

**English** | [简体中文](README_zh.md)

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

A reproducible, transcriptome-wide **cis-MR × Bayesian colocalization** atlas for **type 2 diabetes (T2D)**, **coronary artery disease (CAD)**, and **fasting glucose (FG)**. The project quantifies the *operating characteristics* of transcriptome-wide cis-Mendelian randomization (cis-MR) as a gene-prioritization screen, provides a calibrated colocalization atlas as a public resource, and identifies candidate effector genes with orthogonal (MR + coloc + directional-replication) support.

> **Note on the analysis version.** The per-outcome BH-FDR pipeline ("FDR-core", finalized 2026-08-16) is the authoritative analysis and the one used in the manuscript. All submission-facing numbers appear in [`docs/manuscript/manuscript.md`](docs/manuscript/manuscript.md) and the result tables in `results/`.

---

## Abstract

Genome-wide association studies (GWAS) have identified hundreds of loci for type 2 diabetes (T2D), coronary artery disease (CAD), and fasting glucose (FG), yet most effector genes remain unknown. Transcriptome-wide Mendelian randomization (cis-MR) followed by Bayesian colocalization is widely used to prioritize genes, but its operating characteristics have not been systematically quantified. We performed a cis-MR × coloc scan of eQTLGen whole-blood cis-eQTLs (n = 31,684) against three GWAS across 31,371 gene–trait pairs. Under per-outcome FDR control (q < 0.05), 982 pairs were MR-significant, of which 121 reached strong colocalization (PP.H4 ≥ 0.8): a coloc yield of 12.3%. Yield rose monotonically from 3.0% at nominal p < 0.05 to 25.6% at p < 1e-5 in threshold-sensitivity analyses, and strong colocalization was essentially absent outside the MR-significant set (2/27,123 pairs). All 106 previously reported loci were reproduced, and 15 additional candidate effector genes were identified, with directional replication in GTEx and FinnGen (four reaching nominal significance in FinnGen). We disclose calibration caveats and provide a public atlas of calibrated colocalization support.

---

## Key findings

| Finding | Value |
|---|---|
| Gene–trait pairs scanned (QC-passed) | 31,371 (31,373 raw) |
| MR-significant under per-outcome BH-FDR (q < 0.05) | **982** (T2D 394 / CAD 576 / FG 12) |
| Strong colocalization (PP.H4 ≥ 0.8) | **121** (T2D 65 / CAD 54 / FG 2) → **coloc yield 12.3%** |
| Previously reported loci reproduced | **106 / 106** |
| New candidate effector genes | **15** (9 at known loci + 6 not in GWAS Catalog) |
| Yield calibration (nominal MR p threshold) | 3.0% at p < 0.05 → 25.6% at p < 1e-5 (monotonic) |
| Consistency with stage-2 grid scan | 12.3% (FDR-core) vs 12.96% (grid) |
| Strong coloc outside the MR-significant set | 2 / 27,123 (decisive negative boundary) |

**Candidate effector genes (15):** SLC12A3, CWF19L1, U6atac, CD101, RBM6, CNNM2, N4BP2L2, RIC8A, C2orf49 (within known T2D/CAD risk loci); PLAUR, TAGLN2, VSIG8, PDCD6, CLEC3B, CCDC19 (regions without a T2D/CAD GWAS Catalog record).

**Independent replication (15-candidate set):**
- **GTEx v8** eQTL direction: 6/7 consistent (1 conflict reported: VSIG8).
- **FinnGen R11** (independent cohort): 9/9 gene-level direction-consistent in the alignable subset (8/9 variant-level), 4 reaching nominal significance in FinnGen (RBM6, CNNM2, CD101, RIC8A); alignment coverage 9/15 = 60%.

**Methodological operating characteristics:**
- SMR + HEIDI concordance: 76/106 = 71.7%.
- Steiger direction: 73/76 = 96.1%.
- Permutation false-positive rate under the null: 1.45% (154 / 10,600).

**Honesty caveats (disclosed in the manuscript):**
- 41/106 = 38.7% of known strong-colocalization regions have a GWAS peak p < 5e-8 → most colocalized regions do not qualify as new GWAS loci and are interpreted as candidate effector genes.
- `coloc.susie` does not fully converge under external-LD estimation (exploratory); LAMC1 is excluded on both FDR and multi-signal evidence.
- Candidate effector genes are *candidate assessments*, not causal discoveries.

---

## Manuscript

- [`docs/manuscript/manuscript.md`](docs/manuscript/manuscript.md) — full manuscript source (title page, abstract, IMRaD, tables, figures, references).
- [`docs/manuscript/AJHG_submission_Qiushuo_Geng_20260816.docx`](docs/manuscript/AJHG_submission_Qiushuo_Geng_20260816.docx) — formatted submission document.
- [`docs/manuscript/cover_letter.md`](docs/manuscript/cover_letter.md) — cover letter.

---

## Repository structure

```
scripts/   All analysis scripts (M20–M28, M33–M49; see Reproducibility)
results/   All result tables (CSV), summary statistics, and figures
  coloc_full_{t2d,cad,fbg}_20260815.csv   full scan output (31,373 pairs)
  fdr_core_20260816.csv                   FDR-core MR-significant set (982)
  strong_all_subset_20260816.csv          121 strong coloc + 2 gray-zone pairs
  candidate15_replication_20260816.csv    15 candidates with replication stats
  grid/                                    stage-2 grid scan outputs
  figures/                                 publication figures
docs/      Manuscript
LICENSE    CC BY 4.0
```

---

## Reproducibility

The full pipeline (input → script → output) is documented in the analysis scripts themselves. Key steps:

| Step | Script | Input | Output |
|---|---|---|---|
| Full transcriptome cis-MR × coloc scan | `M20*`–`M24` | eQTLGen + 3 GWAS (rsID-matched, hg19) | `results/coloc_full_{t2d,cad,fbg}_20260815.csv` |
| New strong-coloc candidate discovery | `M25` / `M25b` | `coloc_full_*` | `results/m25_new_strong_annotation_20260816.csv` |
| GTEx v8 independent direction replication | `M26` | GTEx v8 | `results/m26_gtex_replication_new23_20260816.csv` |
| Nominal precision funnel | `M27` | `coloc_full_*` | `results/m27_precision_funnel_20260816.csv` |
| FinnGen R11 independent-cohort replication | `M28` | FinnGen R11 sumstats | `results/m28_finngen_replication_new23_20260816.csv` |
| **Per-outcome BH-FDR recomputation (FDR-core)** | `M36b_fdr_recompute_20260816.py` | `coloc_full_*` + grid | `results/fdr_core_20260816.csv`, `candidate15_replication_20260816.csv`, `m36b_funnel_20260816.csv` |
| coloc.susie sensitivity (exploratory) | `M34b` | 6 loci + 1000G EUR LD | `results/m34_coloc_susie_20260816.csv` |
| **Ridge-shrunk LD sweep (coloc.susie recovery)** | `M34c_coloc_susie_ridge.sh` | M34b aligned z/R | `results/m34c_coloc_susie_ridge_20260817.csv` |
| **Hyprcoloc 2-trait cross-check** | `M40_hyprcoloc_2trait_20260817.R` | eQTLGen + GWAS + 1000G LD | `results/m40_hyprcoloc_2trait_20260817.csv` |
| **GCTA-COJO conditional/joint fine-mapping** | `M46_cojo_slct_20260817.sh` | OpenGWAS + 1000G EUR bfile | `results/m46_cojo_20260817.csv` |
| **MR-Steiger corrected direction (15 candidates)** | `M44_steiger_corrected_20260817.R` | candidate15 + eQTLGen + GWAS | `results/m44_steiger_corrected_20260817.csv` |
| **Empirical-Bayes coloc enrichment (TORUS/fastENLOC-style)** | `M45_ebayes_enrichment_20260817.py` | coloc_full + FDR-core | `results/m45_ebayes_summary_20260817.csv` |
| **STRING RWR network proximity** | `M39_rwr_network_20260817.py` | STRING v12 + 106 seeds | `results/m39_rwr_20260817.csv` |
| **Cross-tissue lead-eQTL concordance (GTEx v8)** | `M40a_cross_tissue_20260817.py` | GTEx v8 egenes | `results/m40a_cross_tissue_20260817.csv` |
| **PheWAS pleiotropy screen (GWAS Catalog)** | `M41_phewas_catalog_20260817.py` | GWAS Catalog + 13 leads | `results/m41_phewas_snp_20260817.csv` |
| **Cell-type enrichment (PanglaoDB)** | `M43_celltype_enrichment_20260817.py` | PanglaoDB markers | `results/m43_celltype_20260817.csv` |
| **pQTL coverage (deCODE/UKB-PPP)** | `M42_pqtl_coverage_20260817.py` | OpenGWAS pQTL | `results/m42_pqtl_coverage_20260817.csv`, `m42_pqtl_summary_20260817.csv` |
| Figures (5 main + S1–S6) | `M37` / `M38` / `M47`–`M49` | result tables | `results/figures/20260816_Fig{1..5}_*.png`, `20260817_FigS{1..6}_*.png` |
| Word manuscript | `M36_build_word_ajhg_20260816.py` | `docs/manuscript/*` | `docs/manuscript/AJHG_submission_Qiushuo_Geng_20260816.docx` |

**Path placeholders.** Scripts use portable placeholders instead of machine-specific paths: `<repo-root>` (repository root), `<scratch>` (staging directory for downloaded source data), and `<conda-root>` (R environment prefix). Replace them with your local paths before re-running; bulk source data (eQTLGen, OpenGWAS, GTEx, FinnGen, 1000 Genomes) must be downloaded separately.

---

## Data availability

**Input data (public):**
- eQTLGen whole-blood cis-eQTLs (n = 31,684; Võsa et al., 2021).
- GWAS: T2D (`ebi-a-GCST006867`), CAD (`ebi-a-GCST005194`), FG (`ebi-a-GCST005186`) via OpenGWAS.
- GTEx v8, FinnGen R11, 1000 Genomes Phase 3 (EUR LD reference).

**Output data:** all analysis outputs (result tables, summary statistics, figures, scripts) are provided in this repository.

**Archival copy:** a versioned archival copy of this repository is deposited on Zenodo.

> DOI: [10.5281/zenodo.21967917](https://doi.org/10.5281/zenodo.21967917)

---

## License

This repository is released under the [Creative Commons Attribution 4.0 International (CC BY 4.0)](LICENSE) license.

## Citation

Qiushuo Geng. *Transcriptome-wide cis-MR and colocalization atlas for type 2 diabetes, coronary artery disease, and fasting glucose: operating characteristics and candidate effector genes.* Zenodo: [10.5281/zenodo.21967917](https://doi.org/10.5281/zenodo.21967917).
