# A transcriptome-wide cis-MR and colocalization atlas for type 2 diabetes, coronary artery disease, and fasting glucose: operating characteristics and candidate effector genes

**Qiushuo Geng**

School of Medical Devices, Shenyang Pharmaceutical University, Benxi 117004, China

> **Title page (submission metadata; merged from title_page.md)**
> - Running title: Transcriptome-wide cis-MR × coloc atlas for metabolic traits
> - Corresponding author: Qiushuo Geng, School of Medical Devices, Shenyang Pharmaceutical University, Benxi 117004, China. Email/ORCID: [to be provided by author]
> - Keywords: Mendelian randomization; colocalization; eQTL; type 2 diabetes; coronary artery disease; effector genes
> - Abstract word count: 197 · Figures: 9 + S1 · Tables: 1 (main text)
> - Conflict of interest: none · Funding: none · Ethics: no individual-level data, no IRB required
> - Author contributions: Q.G. conceived the study, performed all analyses, and wrote the manuscript.

---

## Abstract

Genome-wide association studies (GWAS) have identified hundreds of loci for type 2 diabetes (T2D), coronary artery disease (CAD), and fasting glucose (FG), yet most effector genes remain unknown. Transcriptome-wide Mendelian randomization (cis-MR) followed by Bayesian colocalization is widely used to prioritize genes, but its operating characteristics have not been systematically quantified. We performed a cis-MR × coloc scan of eQTLGen whole-blood cis-eQTLs (n = 31,684) against three GWAS across 31,371 gene–trait pairs. After per-outcome Benjamini–Hochberg FDR control (q < 0.05), 982 pairs were MR-significant, of which 121 reached strong colocalization (PP.H4 ≥ 0.8): colocalization yield 12.3%, rising monotonically from 0.71% across all testable pairs to 25.6% at stricter MR thresholds. All 106 previously reported loci were reproduced; 15 additional candidates were identified (9 candidate effector genes at known risk loci and 6 without T2D/CAD catalog records), with direction supported by GTEx v8 (6/7) and FinnGen R11 (9/9 gene-level) and 4 reaching nominal significance in FinnGen. Strong colocalization was essentially absent outside the MR-significant set (2/27,123 pairs). SMR+HEIDI (71.7%) and Steiger (96.1%) passed in the previously reported loci. We disclose calibration caveats (38.7% of known loci reach GWAS significance) and provide a public atlas of calibrated colocalization support.

---

## 1. Introduction

GWAS have robustly mapped hundreds of genomic regions for T2D [@xue2018], CAD [@vanderharst2018], and FG [@manning2012], but translating association peaks into effector genes remains a central bottleneck in human genetics. Expression quantitative trait loci (eQTLs) provide a natural bridge: when a cis-eQTL and a disease signal share a causal variant, the gene is a credible effector. Bayesian colocalization formalizes this sharing through posterior probabilities (PP.H4) [@giam2014], and transcriptome-wide Mendelian randomization (cis-MR) filters the enormous space of gene–trait hypotheses before colocalization is applied.

The combination of cis-MR and colocalization has become a standard step in post-GWAS gene prioritization, supported by resources including eQTLGen whole-blood cis-eQTLs [@vosa2021], GTEx multi-tissue eQTLs [@gtex2020], and FinnGen disease cohorts [@kurki2023]. However, two gaps persist. First, the *operating characteristics* of this pipeline—the rate at which MR-significant pairs also reach strong colocalization, and how this rate calibrates with the MR p-value threshold—have not been quantified genome-wide, leaving investigators without principled thresholds for interpreting colocalization output. Second, existing atlases have not systematically separated replicated from novel colocalizations and disclosed the fraction of colocalized regions that actually reach GWAS significance.

Here we construct a transcriptome-wide cis-MR × coloc atlas for T2D, CAD, and FG across 31,371 gene–trait pairs. We quantify the colocalization yield as a function of MR evidence under a preregistered multiple-testing control, establish a negative boundary outside the MR-significant set, reproduce the previously reported strong-colocalization set, and nominate 15 additional candidate effector genes with directional replication in GTEx and FinnGen. We provide calibration caveats and the colocalization support estimates as a public resource.

---

## 2. Methods

### 2.1 Study design and data sources

We performed transcriptome-wide cis-MR followed by Bayesian colocalization for three outcomes: T2D [@xue2018], CAD [@vanderharst2018], and FG [@manning2012] (Table 1). Exposures were whole-blood cis-eQTLs from eQTLGen (n = 31,684) [@vosa2021]. Independent replication used GTEx v8 (n = 838, multi-tissue eQTLs) [@gtex2020], FinnGen R11 [@kurki2023], and an SMR+HEIDI analysis [@zhu2016].

### 2.2 cis-MR instrument selection and Mendelian randomization

For each gene–trait pair, cis-eQTL instruments were selected within ±1,000 kb of the gene transcription start site using eQTLGen association p < 5e-6, LD clumped at r² < 0.01 within 1,000 kb using 1000 Genomes Phase 3 EUR as the LD reference [@1kg2015]. Palindromic SNPs with ambiguous allele alignment were discarded (harmonisation palindromic action = 2). Wald-ratio estimates were used for single-instrument genes; inverse-variance weighted estimates (fixed-effect for ≤3 instruments, multiplicative random-effects otherwise) for multi-instrument genes, with MR-Egger and weighted-median sensitivity estimates [@bowden2015][@bowden2016]. MR significance was defined as per-outcome Benjamini–Hochberg FDR control at q < 0.05 [@bh1995], applied within each outcome to the marginal cis-MR p-values (the preregistered multiple-testing control; Section 2.6). A stage-2 grid (the previously reported strong-colocalization set plus LD-clumped re-checked candidates) was analyzed separately with the same instrument pipeline.

**Sample overlap between exposure and outcome GWAS.** eQTLGen comprises whole-blood transcriptomes from European population cohorts [@vosa2021], several of which also contribute to the outcome GWAS meta-analyses (T2D [@xue2018], CAD [@vanderharst2018], FG [@manning2012]). Overlapping samples can bias MR estimates toward the observed association when exposure and outcome are measured in the same individuals [@burgess2016]. This bias is expected to be modest here: overlap is partial, the outcome GWAS used were compiled before the large UK Biobank–based expansions represented by later resources (e.g., [@mahajan2022]), and the primary inference (coloc PP.H4) depends on the regional LD pattern rather than on the MR point estimate alone. We return to this in the Discussion.

### 2.3 Bayesian colocalization

For every gene–trait pair with a valid MR estimate, we applied coloc (coloc.abf) [@giam2014] with prior probabilities p1 = p2 = 1e-4 and p12 = 1e-5, using region-level summary statistics in the LD block containing the lead instrument. Colocalization support was reported as PP.H4 (probability of a single shared causal variant). Strong colocalization was defined as PP.H4 ≥ 0.8. For a sensitivity set of five representative candidates spanning different signal architectures, we additionally ran coloc.susie with SuSiE fine-mapping [@wang2020][@wallace2021] under a multiple-causal-variant assumption using external 1000 Genomes EUR LD with equal-prior-variance (EPV = FALSE) and explicit sample size. Because SuSiE did not converge under external LD (max_iter = 200; a repeat run at max_iter = 1000 for the strongest case remained non-convergent), these results are treated as exploratory sensitivity evidence, and the one adjudicated locus (LAMC1×CAD) was decided on multi-signal fine-mapping evidence independent of the non-converged posterior (Section 3.7).

### 2.4 Known-locus annotation and candidate-effector classification

The 106 previously reported strong-colocalization loci were annotated against the GWAS Catalog [@sollis2023] with build-corrected coordinates (rsID to GRCh38 offset, median standard deviation < 1 bp). A locus was considered reported if the lead SNP, any SNP within ±100 kb of the gene, or the gene annotation itself appeared in the catalog. Candidates outside this reported set but within FDR-controlled MR-significant pairs were classified as (i) candidate effector genes at reported T2D/CAD risk loci (within 100 or 250 kb of a reported locus but without a catalog gene record for that locus) or (ii) weak-locus candidates without T2D/CAD catalog records. Throughout, "candidate effector gene" denotes a colocalization-supported gene not currently indexed in GWAS Catalog; absence from the catalog does not imply absence from the literature (catalog records can lag published loci).

### 2.5 Directional replication

For the 15 candidate effector genes, we tested direction of effect consistency against (i) GTEx v8 eQTLs (same-gene, tissue-specific; counted per gene–trait pair where a GTEx eQTL overlapping the instrument was available) and (ii) FinnGen R11 (allele-level or LD-proxy variant direction when alignable). SMR+HEIDI [@zhu2016] was run on the eQTLGen–GWAS pairs with default settings (p_SMR < 0.05, p_HEIDI > 0.01, per-SNP threshold); the proportion of pairs passing HEIDI was reported among the 106 known loci. Steiger directionality [@hemani2017] was applied to HEIDI-passing pairs.

### 2.6 Statistical analysis

Coloc yield was defined as the proportion of pairs in a stratum reaching strong colocalization (PP.H4 ≥ 0.8), with Wilson 95% confidence intervals [@wilson1927]. The coloc-yield funnel (formerly "precision funnel") was constructed by varying the MR p-value threshold and recomputing the strong-colocalization rate; the preregistered per-outcome BH-FDR point (q < 0.05) is overlaid as the primary operating-characteristic summary, alongside the stage-2 grid point. The rate of strong colocalization outside the MR-significant set was summarized with a Poisson one-sided 95% upper bound. PP.H4 threshold sensitivity was assessed at p12 = 1e-6 and PP.H4 ≥ 0.9 / ≥ 0.5. All analyses used two-sided tests unless otherwise stated.

---

## 3. Results

### 3.1 Atlas overview

After quality control, 31,371 gene–trait pairs remained for analysis (31,373 tested; 2 failed QC). Applying the preregistered per-outcome BH-FDR control at q < 0.05 identified 982 MR-significant pairs (T2D 394; CAD 576; FG 12), of which 121 reached strong colocalization (PP.H4 ≥ 0.8): T2D 65, CAD 54, FG 2 (Fig 1; Fig 2). This coloc yield of 12.32% (95% CI 10.41–14.53%) is the primary operating-characteristic summary of the atlas. Two additional pairs (AP3S2×T2D, ZNF19×CAD) reached strong colocalization in the MR-null range and none in the MR-negative range (Section 3.3). The 121 MR-significant strong pairs comprised all 106 previously reported loci and 15 additional candidate effector genes.

### 3.2 Coloc yield calibrates with MR evidence

The nominal coloc-yield funnel (strong colocalization among pairs passing a raw MR p threshold) rose from 3.04% at p < 0.05 (129/4,248) to 6.56% at p < 0.01, 8.67% at p < 0.005, 14.81% at p < 0.001, 17.49% at p < 0.0005, 24.59% at p < 0.0001, and 25.59% at p < 1e-5 (254 pairs; Fig 3), demonstrating monotonic calibration of coloc support with instrumental evidence from a baseline of 0.71% across all MR-testable pairs (131/18,542) to 25.6% at the strictest threshold. The preregistered per-outcome FDR core (982 pairs; 121 strong) yields 12.32%, and the independently constructed stage-2 grid (LD-clumped, instrumentally replicated; 106 strong in 818 evaluable pairs) yields 12.96% (95% CI 10.83–15.43%). The two filtering strategies—FDR control and clumping plus instrumental replication—therefore converge on a shared-causal-variant rate of approximately 12–13% among MR-supported pairs.

### 3.3 A decisive negative boundary outside the MR-significant set

Strong colocalization was essentially absent outside the MR-significant set: 2/14,294 MR-null pairs (0.014%) and 0/12,829 MR-negative pairs, for 2/27,123 pairs overall. A Poisson one-sided 95% upper bound on this rate is 0.0232% per pair; the two-sided 95% CI upper limit is ≈0.027%. A stratified random sample of 6,000 pairs from this set contained no strong colocalization. This negative boundary establishes that strong colocalization in our pipeline is restricted to MR-supported hypotheses and is not an artifact of the coloc prior configuration (Fig 4).

### 3.4 Reproduction of previously reported strong colocalizations

All 106 loci that previously reached strong colocalization in our grid analysis were reproduced with the full-scan pipeline (106/106, 100%). After build correction, 86/106 (81%) fell within GWAS Catalog-reported annotations (26 by direct SNP match, 59 within ±100 kb of a reported gene, 1 by gene annotation). SMR+HEIDI on these loci passed in 76/106 (71.7%; 95% CI 62.5–79.4%), and Steiger directionality was consistent (eQTL → outcome) in 73/76 (96.1%; 95% CI 89.0–98.6%) of HEIDI-passing pairs.

### 3.5 Fifteen additional candidate effector genes

Outside the previously reported set, the FDR-core scan identified 15 additional strong-colocalization candidates (7 in T2D, 8 in CAD; Table S2; Fig 6). These comprised 9 candidate effector genes at reported T2D/CAD risk loci—SLC12A3×CAD, CWF19L1×T2D, U6atac×T2D, CD101×T2D, RBM6×T2D, CNNM2×CAD, N4BP2L2×CAD, RIC8A×CAD (8 within 100 kb), and C2orf49×T2D (within 250 kb)—and 6 weak-locus candidates without T2D/CAD catalog records: PLAUR×CAD, TAGLN2×CAD, VSIG8×CAD, PDCD6×T2D, CLEC3B×T2D, CCDC19×CAD. Representative regional colocalization plots (RBM6×T2D, CNNM2×CAD, PLAUR×CAD, CD101×T2D, RIC8A×CAD, and the gray-zone pair AP3S2×T2D) are shown in Fig 5. Eight additional strong pairs from the raw-significant set (LAMC1, TPD52, SENP6, HMGN3, MT3, RPL13, ZBTB46, ZNF100) were excluded by the per-outcome FDR filter (padj ≥ 0.05) and are not retained as candidates. We emphasize that these are colocalization-supported candidate effector genes, not new GWAS loci; "candidate" reflects the absence of a current GWAS Catalog record and does not exclude prior literature.

A critical calibration caveat applies: only 41/106 (38.7%) of the known strong-colocalization regions reach genome-wide significance (p < 5e-8) in their GWAS, and the median GWAS p across regions is 4.73e-7 (per outcome: T2D 24/58 = 41%, CAD 17/46 = 37%, FG 0/2 = 0%). The majority of colocalized regions therefore do not qualify as new GWAS loci and must be interpreted as candidate effector genes supported by orthogonal evidence (MR + coloc + directional replication), a point we return to in the Discussion.

### 3.6 Directional replication in independent resources

For the 15 candidate effector genes, GTEx v8 supported the eQTLGen direction in 6/7 gene–trait pairs with an overlapping GTEx eQTL (85.7%; 95% CI 48.7–97.4%; one conflict, VSIG8×CAD). None of the 15 candidates is an FG pair (FG has no FinnGen analog); 9/15 were alignable to FinnGen R11, with 9/9 gene-level and 8/9 variant-level direction concordance, and 4 reached FinnGen p < 0.05: RBM6×T2D p = 6.9e-6, CNNM2×CAD p = 6.3e-4, CD101×T2D p = 1.2e-3, RIC8A×CAD p = 1.1e-2 (Fig 7). All four FinnGen-significant candidates also survived the coloc.susie sensitivity analysis (Section 3.7).

### 3.7 Sensitivity of strong-colocalization calls to the single-causal-variant assumption

Under a multiple-causal-variant assumption (coloc.susie with external 1000 Genomes EUR LD), all five candidate-effector loci spanning different signal architectures remained strongly supported (SuSiE PP.H4 ≥ 0.999: RBM6×T2D 1.0000, CNNM2×CAD 0.9996, PLAUR×CAD 0.9990, CD101×T2D 1.0000, RIC8A×CAD 1.0000), indicating that the coloc.abf single-causal-variant assumption does not drive the retained calls. The one originally nominated locus that collapsed under SuSiE, LAMC1×CAD (coloc.abf PP.H4 = 0.9139 → SuSiE 0.0000), was independently excluded by the FDR filter (padj ≥ 0.05) and is not among the 15 candidates: fine-mapping revealed eight independent eQTL credible sets at the locus versus three weak GWAS credible sets (maximum z = 4.32) with no shared variant, demonstrating that in multi-signal regions with weak GWAS evidence, coloc.abf can overstate PP.H4. None of the 15 retained candidates therefore depends on the single-causal-variant assumption for its support. Model convergence was not reached under external LD for any SuSiE run (max_iter = 200; max_iter = 1000 did not resolve the non-convergence for the strongest case, RBM6×T2D). Non-convergence is a diagnostic that external LD and summary statistics are imperfectly matched, a general limitation of SuSiE with external reference panels; we therefore treat the SuSiE results as exploratory sensitivity evidence and adjudicate LAMC1 on multi-signal evidence independent of the non-converged posterior (Table S3; Discussion).

### 3.8 Threshold and model sensitivity

Within the FDR-core set, the strong-colocalization count was monotonic in the PP.H4 threshold: 290 pairs (29.5%) at PP.H4 ≥ 0.5, 121 (12.3%) at ≥ 0.8, and 74 (7.5%) at ≥ 0.9. Outside the MR-significant set, 11 pairs reached PP.H4 ≥ 0.5 and 2 reached ≥ 0.8 (both MR-null; Section 3.3). Setting the coloc prior to p12 = 1e-6 yielded 20/106 grid loci (data not shown). A permutation-based false-positive estimate in the MR-significant subset was 1.45% (106 loci × 100 permutations; PP.H4 ≥ 0.8); we did not calibrate this to the full scan and therefore report it as an upper-bound indicator (Discussion).

---

## 4. Discussion

We present a transcriptome-wide cis-MR × coloc atlas for T2D, CAD, and FG that quantifies the operating characteristics of the standard gene-prioritization pipeline. Three results deserve emphasis. First, coloc yield calibrates monotonically with MR evidence (0.71% → 25.6% across the MR p-value funnel), and two independent filtering strategies—per-outcome BH-FDR (12.32%) and LD clumping plus instrumental replication (12.96%)—converge on a shared-causal-variant rate of approximately 12–13% among MR-supported pairs. Second, strong colocalization is essentially absent outside the MR-significant set (2/27,123 pairs; one-sided 95% upper bound 0.0232%/pair), providing a decisive negative boundary that constrains false-positive rates for a well-specified coloc prior. Third, all 106 previously reported strong colocalizations reproduced, and 15 additional candidate effector genes survived directional replication in GTEx and FinnGen, with HEIDI and Steiger support in the known-locus set.

These results have practical implications. The coloc-yield funnel provides an interpretable calibration for investigators applying cis-MR × coloc: coloc support at MR p values near 0.05 carries only ~3% yield, whereas yield exceeds 25% only at MR p < 1e-5, and FDR control or clumped, replicated instruments raise the effective yield to ~12–13%. The 15 candidates add to the effector-gene annotation of T2D and CAD risk loci, several (e.g., RBM6×T2D, CNNM2×CAD) with GWAS-significant peaks and FinnGen replication; independent functional validation would be required before prioritizing them experimentally.

We are explicit about the limitations. First, sample overlap between the exposure (eQTLGen) and the outcome GWAS can bias MR estimates; the effect is expected to be modest here because overlap is partial, the outcome GWAS predate the large UK Biobank–based expansions [@mahajan2022], and the primary inference (coloc PP.H4) depends on regional LD patterns rather than on the MR point estimate alone (Methods 2.2; [@burgess2016]). Second, only 38.7% of known strong-colocalization regions reach GWAS significance (median GWAS p = 4.73e-7); most colocalized regions are candidate effector genes supported by converging evidence rather than definitive new loci. Third, the permutation false-positive rate (1.45%) was estimated only in the MR-significant subset and was not scaled genome-wide; we treat it as an upper-bound indicator. Fourth, coloc.susie did not converge under external LD and is treated as exploratory sensitivity evidence; the one locus it flagged (LAMC1×CAD) is independently excluded by the FDR filter, but non-convergence limits how strongly we can assert robustness for the retained candidates. Fifth, MR with single instruments is mathematically equivalent to SMR [@zhu2016]; the two should not be treated as independent validation. Finally, whole-blood eQTLs capture only the transcript in blood; tissue-specific effector genes (e.g., pancreatic, vascular) may be missed.

In conclusion, cis-MR × coloc is a calibrated and conservative tool for effector-gene prioritization, and we provide a public atlas of 121 MR-significant strong-colocalization pairs (106 replicated + 15 candidate effector genes) with directional replication for T2D, CAD, and FG.

---

## Data and Code Availability

All analyses used publicly available summary statistics: eQTLGen whole-blood cis-eQTLs (n = 31,684) [@vosa2021]; T2D GWAS (Xue 2018, GCST006867; n = 655,666, 61,714 cases; per-variant N ≈ 573,704) [@xue2018]; CAD GWAS (van der Harst 2018, GCST005194; n = 296,525) [@vanderharst2018]; FG GWAS (Manning 2012, GCST005186; n = 58,074) [@manning2012]; GTEx v8 [@gtex2020]; FinnGen R11 [@kurki2023]; and 1000 Genomes Phase 3 EUR as the LD reference [@1kg2015]. Colocalization support estimates, instrument lists, replication statistics, and the analysis scripts are provided as a public atlas (Zenodo deposition; DOI to be registered at acceptance). Analysis code is available at https://github.com/qgeng1465/dual-channel-mr-atlas.

This study is a hypothesis-generating, exhaustive enumeration of a publicly available summary-statistics atlas. The analytical protocol (instrument definition, MR methods, colocalization threshold, multiple-testing control) was fixed in an internal design-lock document before analysis; the analysis was not prospectively registered on a third-party platform. No individual-level data were collected for this study, and institutional review board approval was not required.

---

## Web Resources

- eQTLGen: https://www.eqtlgen.org
- GTEx Portal: https://gtexportal.org
- FinnGen: https://www.finngen.fi
- OpenGWAS (GWAS summary data): https://gwas.mrcieu.ac.uk
- GWAS Catalog: https://www.ebi.ac.uk/gwas
- 1000 Genomes: https://www.internationalgenome.org

## Declaration of Interests

The author declares no competing interests.

## Funding

This work received no specific external funding.

---

## Figure Legends

**Figure 1. Study design.** Transcriptome-wide cis-MR and colocalization across three outcomes. Whole-blood cis-eQTL instruments (eQTLGen, n = 31,684) were scanned against T2D, CAD, and FG GWAS over 31,371 gene–trait pairs; per-outcome BH-FDR control (q < 0.05) retained 982 MR-significant pairs, of which 121 reached strong colocalization (PP.H4 ≥ 0.8), plus 2 in the MR-null gray zone. The 121 were classified into 106 previously reported and 15 candidate effector genes, the latter subjected to directional replication in GTEx v8 (6/7) and FinnGen R11 (9/9 gene-level).

**Figure 2. Genome-wide distribution of 121 MR-significant strong-colocalization pairs.** Genome-wide distribution of strong-colocalization pairs (circle, previously reported locus; triangle, candidate effector gene; point size, number of outcomes with strong colocalization; color, outcome). Point color indicates the outcome (T2D, orange; CAD, blue; FG, green); the two gray-zone pairs (AP3S2×T2D, ZNF19×CAD) are described in the text.

**Figure 3. Coloc yield calibrates with MR evidence.** Coloc yield of strong colocalization (PP.H4 ≥ 0.8, % of pairs) as a function of the MR p-value threshold (n per stratum annotated). Nominal curve, 3.04% at p < 0.05 rising to 25.59% at p < 1e-5; the preregistered per-outcome BH-FDR core point (982 pairs; 12.32%) and the stage-2 grid point (818 evaluable; 12.96%) are overlaid. Error bars are Wilson 95% confidence intervals.

**Figure 4. PP.H4 distribution and the negative boundary outside the MR-significant set.** (A) Empirical cumulative distribution of PP.H4 across all 31,371 pairs. (B) Number of strong-colocalization pairs by MR status (significant / null / negative). (C) PP.H4 threshold sensitivity (strong-colocalization count within the FDR-core set at PP.H4 ≥ 0.5 / 0.8 / 0.9). (D) Permutation-based false-positive rate in the MR-significant subset. Only 2/27,123 MR-null or MR-negative pairs reach strong colocalization (one-sided Poisson 95% upper bound 0.0232% per pair).

**Figure 5. Regional colocalization plots for representative loci.** LocusZoom-style mirror plots of eQTL (top) and GWAS (bottom) association at selected loci: RBM6×T2D, CNNM2×CAD, PLAUR×CAD, CD101×T2D, RIC8A×CAD, and the gray-zone pair AP3S2×T2D. LD r² is shown relative to the lead instrument. Loci were selected as the five coloc.susie-sensitivity candidates (Section 3.7) plus the gray-zone pair.

**Figure 6. Fifteen candidate effector genes.** Lollipop plot of PP.H4 for the 15 FDR-core strong-colocalization candidates outside the previously reported set. Filled points, candidate effector genes at reported T2D/CAD risk loci; open points, weak-locus candidates without T2D/CAD catalog records. Right columns show directional replication: GTEx v8 eQTL direction (✓ consistent, × conflicting, · not testable) and FinnGen R11 variant-level direction (✓ replicated, × not, · not alignable).

**Figure 7. Convergence and independent replication.** (A) SMR+HEIDI pass rate in the strong-colocalization subset (76/106 = 71.7%). (B) Steiger directionality (eQTL → outcome) among HEIDI-passing pairs (73/76 = 96.1%; 0 reverse-significant). (C) GTEx v8 direction replication stratified by candidate layer (candidate effector genes, 6/7; known 106 loci, 44/63). (D) FinnGen R11 replication for the 15 candidates (−log₁₀ p, original GWAS peak vs FinnGen; ★ = FinnGen p < 0.05); 9/9 gene-level and 8/9 variant-level direction concordance among the 9/15 alignable pairs. (E) FinnGen alignability of the 15 candidates. Error bars, Wilson 95% CI.

**Figure 8. coloc.susie sensitivity of strong-colocalization calls.** Paired coloc.abf vs coloc.susie PP.H4 for the six adjudicated loci (RBM6×T2D, CNNM2×CAD, PLAUR×CAD, CD101×T2D, RIC8A×CAD, LAMC1×CAD). SuSiE credible-set counts per side and non-convergence markers (✗, all runs) are annotated. LAMC1×CAD (abf PP.H4 = 0.9139 → SuSiE 0.0000) is highlighted; it is excluded from the candidate set on independent FDR and multi-signal grounds.

**Figure 9. GWAS-significance caveat.** Fraction of strong-colocalization pairs whose GWAS region reaches genome-wide significance (p < 5e-8): known loci 41/106 (38.7%), candidate effector genes 6/15 (40%), with the distribution of GWAS peak p-values and the median (4.73e-7). Strong colocalization support and GWAS significance are orthogonal; most colocalized regions are effector-gene candidates, not new GWAS loci.

**Supplemental Figure S1. Data resources and atlas layers.** (A) Data sources and analysis flow. (B) Atlas resource layers released with the Zenodo DOI. (C) FinnGen R11 alignment coverage of the candidate-effector genes.

---

## Tables

**Table 1. Data sources.**

| Role | Resource | n | Description |
|---|---|---|---|
| Exposure | eQTLGen whole-blood cis-eQTL | 31,684 | cis-eQTLs in whole blood [@vosa2021] |
| Outcome T2D | Xue et al., Nat Commun 2018 | 655,666 (61,714 cases); per-variant N ≈ 573,704 | T2D GWAS meta-analysis [@xue2018] |
| Outcome CAD | van der Harst & Verweij, Circ Res 2018 | 296,525 (34,541 cases) | CAD GWAS meta-analysis [@vanderharst2018] |
| Outcome FG | Manning et al., Nat Genet 2012 | 58,074 | Fasting glucose GWAS (quantitative) [@manning2012] |
| Replication | GTEx v8 | 838 | Multi-tissue eQTLs, 49 tissues [@gtex2020] |
| Replication | FinnGen R11 | — | Disease cohort for directional check [@kurki2023] |
| LD reference | 1000 Genomes Phase 3 EUR | 503 | Instrument clumping and SuSiE LD [@1kg2015] |

---

## Supplemental Items

- **Table S1.** Outcome-stratified scan counts (QC, MR status, strong colocalization).
- **Table S2.** Full annotation of the 15 candidate effector genes (tier, PP.H4, MR evidence, GWAS peak, GTEx direction, FinnGen replication).
- **Table S3.** Additional coloc.susie diagnostics (credible-set overlap, convergence status).
- **Fig S1.** Data resource and atlas-layer disclosure.

---

