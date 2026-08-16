# A transcriptome-wide cis-MR and colocalization atlas for type 2 diabetes, coronary artery disease, and fasting glucose: operating characteristics and candidate effector genes

**Qiushuo Geng**<sup>1*</sup>,<sup>1#</sup>

Affiliations:

1 School of Medical Device, Shenyang Pharmaceutical University, Benxi, 117004, China.
2 School of Biomedical Engineering, Tsinghua Medicine, Tsinghua University, Beijing, 100084, China.

<sup>#</sup> Corresponding author
<sup>*</sup> Co-first author

> **Title page (submission metadata)**
> - Running title: Transcriptome-wide cis-MR × coloc atlas for metabolic traits
> - Corresponding author: Qiushuo Geng (co-first and corresponding author)
> - Keywords: Mendelian randomization; colocalization; eQTL; type 2 diabetes; coronary artery disease; effector genes
> - Abstract word count: 169 · Figures: 5 · Tables: 2 (main text) + 4 supplementary
> - Conflict of interest: none · Funding: none · Ethics: no individual-level data, no IRB required
> - Author contributions: Q.G. conceived the study, performed all analyses, and wrote the manuscript.

---

## Abstract

Genome-wide association studies (GWAS) have identified hundreds of loci for type 2 diabetes (T2D), coronary artery disease (CAD), and fasting glucose (FG), yet most effector genes remain unknown. Transcriptome-wide Mendelian randomization (cis-MR) followed by Bayesian colocalization is widely used to prioritize genes, but its operating characteristics have not been systematically quantified. We performed a cis-MR × coloc scan of eQTLGen whole-blood cis-eQTLs (n = 31,684) against three GWAS across 31,371 gene–trait pairs. Under per-outcome FDR control (q < 0.05), 982 pairs were MR-significant, of which 121 reached strong colocalization (PP.H4 ≥ 0.8): a coloc yield of 12.3%. Yield rose monotonically from 3.0% at nominal p < 0.05 to 25.6% at p < 1e-5 in threshold-sensitivity analyses, and strong colocalization was essentially absent outside the MR-significant set (2/27,123 pairs). All 106 previously reported loci were reproduced, and 15 additional candidate effector genes were identified, with directional replication in GTEx and FinnGen (four reaching nominal significance in FinnGen). We disclose calibration caveats and provide a public atlas of calibrated colocalization support.

---

## 1. Introduction

GWAS have robustly mapped hundreds of genomic regions for T2D [@xue2018], CAD [@vanderharst2018], and FG [@manning2012], but translating association peaks into effector genes remains a central bottleneck in human genetics. Expression quantitative trait loci (eQTLs) provide a natural bridge: when a cis-eQTL and a disease signal share a causal variant, the gene is a credible effector. Bayesian colocalization formalizes this sharing through posterior probabilities (PP.H4) [@giam2014] and has been used to identify target genes at disease-associated loci [@hormozdiari2016]. Transcriptome-wide Mendelian randomization (cis-MR) uses genotype as an instrument [@davey2003][@lawlor2008][@dsmith2014][@sanderson2022] and filters the enormous space of gene–trait hypotheses before colocalization is applied. Yet the effector genes through which most GWAS signals act remain unresolved.

The combination of cis-MR and colocalization has become a standard step in post-GWAS gene prioritization, supported by resources including eQTLGen whole-blood cis-eQTLs [@vosa2021], GTEx multi-tissue eQTLs [@gtex2020][@gtex2017][@gtex2015], and FinnGen disease cohorts [@kurki2023]. The eQTLGen resource builds on earlier whole-blood eQTL maps [@westra2013]. However, two gaps persist. First, the *operating characteristics* of this pipeline have not been quantified genome-wide. These include the rate at which MR-significant pairs also reach strong colocalization and how that rate calibrates with the MR p-value threshold. Investigators therefore lack principled thresholds for interpreting colocalization output. Second, because colocalization is almost always applied to MR-significant subsets, the field has been unable to estimate the recall of cis-MR screening (that is, whether colocalized loci are ever missed by the MR filter). Thus, precision (yield) can be estimated, but recall cannot, leaving the operating characteristics of the pipeline incomplete.

Here we construct a transcriptome-wide cis-MR × coloc atlas for T2D, CAD, and FG across 31,371 gene–trait pairs. We quantify the colocalization yield as a function of MR evidence under a preregistered multiple-testing control and establish a negative boundary outside the MR-significant set. We reproduce the previously reported strong-colocalization set and nominate 15 additional candidate effector genes with directional replication in GTEx and FinnGen. We provide calibration caveats and the colocalization support estimates as a public resource.

---

## 2. Methods

### 2.1 Study design and data sources

We performed transcriptome-wide cis-MR followed by Bayesian colocalization for three outcomes: T2D [@xue2018], CAD [@vanderharst2018], and FG [@manning2012] (Table 1). The underlying associations trace to successive T2D meta-analyses [@morris2012][@spracklen2020], CAD meta-analyses [@deloukas2013], and fasting-glucose loci mapped by the MAGIC consortium [@dupuis2010][@scott2012], refined by fine-mapping of T2D loci [@mahajan2018]. The CAD GWAS includes UK Biobank participants [@bycroft2018], and summary statistics were retrieved through OpenGWAS [@elsworth2020]. Exposures were whole-blood cis-eQTLs from eQTLGen (n = 31,684) [@vosa2021]. Independent replication used GTEx v8 (n = 838, multi-tissue eQTLs) [@gtex2020], FinnGen R11 [@kurki2023], and an SMR+HEIDI analysis [@zhu2016].

### 2.2 cis-MR instrument selection and Mendelian randomization

For each gene–trait pair, cis-eQTL instruments were selected within ±1,000 kb of the gene transcription start site using eQTLGen association p < 5e-6, LD clumped at r² < 0.01 within 1,000 kb using 1000 Genomes Phase 3 EUR as the LD reference [@1kg2015]. Palindromic SNPs with ambiguous allele alignment were discarded (harmonisation palindromic action = 2). Wald-ratio estimates were used for single-instrument genes. Inverse-variance weighted estimates (fixed-effect for ≤3 instruments, multiplicative random-effects otherwise) were used for multi-instrument genes with two-sample summary-data MR [@pierce2013][@hemiani2018]. MR-Egger and weighted-median estimates provided sensitivity analyses robust to certain forms of horizontal pleiotropy [@bowden2015][@bowden2016][@bowden2017][@verbanck2018]. We report the maximum cis-eQTL Z-statistic among cis variants for each gene–trait pair as a measure of instrument strength (Z = √F) to assess weak-instrument bias [@burgess2011]; sensitivity analyses followed established guidelines [@burgess2017]. Instrument selection at eQTL p < 5e-6 corresponds to an instrument F-statistic above 20 for every retained instrument. Instrument strength is reported in the Zenodo atlas and summarized visually in Figure 5C. MR significance was defined as per-outcome Benjamini–Hochberg FDR control at q < 0.05 [@bh1995], applied within each outcome to the marginal cis-MR p-values (the preregistered multiple-testing control; Section 2.6). A stage-2 grid was analyzed separately with the same instrument pipeline. The grid is an independently constructed cis-MR × coloc screen from our earlier pipeline, built with LD-clumped, instrumentally replicated instruments; among 818 evaluable gene–trait pairs, 106 reached strong colocalization (12.96%; 95% CI 10.83–15.43%). These 106 are the previously reported loci reproduced by the full-scan pipeline (Section 3.4) and serve as an internal positive control that benchmarks the FDR-core yield (Section 3.2; Table 1).

**Sample overlap between exposure and outcome GWAS.** eQTLGen comprises whole-blood transcriptomes from European population cohorts [@vosa2021], several of which also contribute to the outcome GWAS meta-analyses (T2D [@xue2018], CAD [@vanderharst2018], FG [@manning2012]). Overlapping samples can bias MR estimates toward the observed association when exposure and outcome are measured in the same individuals [@burgess2016]. This bias is expected to be modest here. Overlap is partial, and the outcome GWAS used were compiled before the large UK Biobank–based expansions represented by later resources (e.g., [@mahajan2022]). The primary inference (coloc PP.H4) depends on the regional LD pattern rather than on the MR point estimate alone. We test this directly below for the candidate effector genes and return to it in the Discussion.

**Sample-overlap sensitivity.** To test whether sample overlap between the eQTLGen instrument and the outcome GWAS could distort cis-MR directions, we re-estimated the MR direction of each of the 15 candidate effector genes (Section 3.5). Re-estimation used outcome GWAS with minimal or no UK Biobank overlap: the entirely pre-UKB CARDIoGRAMplusC4D meta-analysis for CAD (Nikpay et al. 2015, ~185,000 samples) [@nikpay2015] and the DIAGRAM trans-ethnic meta-analysis for T2D (Mahajan et al. 2014) [@mahajan2014]. Three genes absent from the trans-ethnic meta (U6atac, PDCD6, CLEC3B) were checked in a European 1000 Genomes–based secondary analysis (Scott et al. 2017) [@scott2017] (Table S4). All 15 candidates (8/8 CAD, 7/7 T2D) showed MR directions concordant with the primary analysis (11 positive, 4 negative; 0 discordant). Because the primary T2D GWAS (Xue et al. 2018) is East Asian, it shares minimal participants with the European eQTLGen exposure. The T2D check is therefore conservative.

### 2.3 Bayesian colocalization

For every gene–trait pair with a valid MR estimate, we applied coloc (coloc.abf) [@giam2014] with prior probabilities p1 = p2 = 1e-4 and p12 = 1e-5, using region-level summary statistics in the LD block containing the lead instrument. Colocalization support was reported as PP.H4 (probability of a single shared causal variant). Strong colocalization was defined as PP.H4 ≥ 0.8. For a sensitivity set of six representative candidates spanning different signal architectures, we additionally ran coloc.susie with SuSiE fine-mapping [@wang2020][@wallace2021] under a multiple-causal-variant assumption. The analysis used external 1000 Genomes EUR LD with equal-prior-variance (EPV = FALSE) and explicit sample size. Because SuSiE did not converge under external LD (max_iter = 200; a repeat run at max_iter = 1000 for the strongest case remained non-convergent), these results are treated as exploratory sensitivity evidence. The one adjudicated locus (LAMC1×CAD) was decided on multi-signal fine-mapping evidence independent of the non-converged posterior (Section 3.7).

### 2.4 Known-locus annotation and candidate-effector classification

The 106 previously reported strong-colocalization loci were annotated against the GWAS Catalog [@sollis2023] with build-corrected coordinates (rsID to GRCh38 offset, median standard deviation < 1 bp). A locus was considered reported if the lead SNP, any SNP within ±100 kb of the gene, or the gene annotation itself appeared in the catalog. Candidates outside this reported set but within FDR-controlled MR-significant pairs were classified as (i) candidate effector genes at reported T2D/CAD risk loci (within 100 or 250 kb of a reported locus but without a catalog gene record for that locus) or (ii) weak-locus candidates without T2D/CAD catalog records. Throughout, "candidate effector gene" denotes a colocalization-supported gene not currently indexed in GWAS Catalog; absence from the catalog does not imply absence from the literature (catalog records can lag published loci).

### 2.5 Directional replication

For the 15 candidate effector genes, we tested direction of effect consistency against (i) GTEx v8 eQTLs (same-gene, tissue-specific) and (ii) FinnGen R11. For GTEx v8, a gene–trait pair was counted where a GTEx eQTL overlapping the instrument was available. For FinnGen R11, we used allele-level or LD-proxy variant direction when alignable. SMR+HEIDI [@zhu2016] was run on the eQTLGen–GWAS pairs with default settings (p_SMR < 0.05, p_HEIDI > 0.01, per-SNP threshold); the proportion of pairs passing HEIDI was reported among the 106 known loci. Steiger directionality [@hemani2017] was applied to HEIDI-passing pairs.

### 2.6 Statistical analysis

Coloc yield was defined as the proportion of pairs in a stratum reaching strong colocalization (PP.H4 ≥ 0.8), with Wilson 95% confidence intervals [@wilson1927]. The coloc-yield funnel (formerly "precision funnel") was constructed by varying the MR p-value threshold and recomputing the strong-colocalization rate. The preregistered per-outcome BH-FDR point (q < 0.05) is overlaid as the primary summary, alongside the stage-2 grid point. The funnel was recomputed within each outcome. Sensitivity of the primary summary to the FDR threshold itself was assessed by recomputing the per-outcome BH-FDR control at stricter q values (q < 0.01, 5e-3, 1e-3, 5e-4, 1e-4). The rate of strong colocalization outside the MR-significant set was summarized with a Poisson one-sided 95% upper bound. PP.H4 threshold sensitivity was assessed at p12 = 1e-6 and PP.H4 ≥ 0.9 / ≥ 0.5. All analyses used two-sided tests unless otherwise stated.

---

## 3. Results

### 3.1 Atlas overview

After quality control, 31,371 gene–trait pairs remained for analysis (31,373 tested; 2 failed QC). Applying the preregistered per-outcome BH-FDR control at q < 0.05 identified 982 MR-significant pairs (T2D 394; CAD 576; FG 12). Of these, 121 reached strong colocalization (PP.H4 ≥ 0.8): T2D 65, CAD 54, FG 2 (Fig 1). The preregistered FDR-core yield was 12.32% (95% CI 10.41–14.53%). All retained cis-eQTL instruments exceeded the conventional weak-instrument threshold (F > 10), mitigating weak-instrument bias (Methods 2.2). The FG scan was the least informative: its GWAS (n = 58,074) is substantially smaller than the T2D (n = 655,666) and CAD (n = 296,525) meta-analyses, so few FG gene–trait pairs survived multiple testing (12 MR-significant, 2 strong colocalizations), and FG has no FinnGen analog for replication; replication and candidate nomination therefore focus on T2D and CAD. Two additional pairs (AP3S2×T2D, ZNF19×CAD) reached strong colocalization in the MR-null range and none in the MR-negative range (Section 3.3). The 121 MR-significant strong pairs comprised all 106 previously reported loci and 15 additional candidate effector genes.

### 3.2 Coloc yield calibrates with MR evidence

The preregistered per-outcome BH-FDR core (Methods 2.6) is the primary summary: among the 982 MR-significant pairs, 121 reached strong colocalization, a coloc yield of 12.32% (95% CI 10.41–14.53%) (Fig 2). The independently constructed stage-2 grid (LD-clumped, instrumentally replicated; 106 strong in 818 evaluable pairs) yields 12.96% (95% CI 10.83–15.43%). The two filtering strategies (per-outcome FDR control and clumping plus instrumental replication) therefore converge on a shared-causal-variant rate of approximately 12–13% among MR-supported pairs.

In an exploratory threshold-sensitivity analysis, the nominal coloc-yield funnel (strong colocalization among pairs passing a raw MR p threshold) rose from 0.71% at p < 0.5 (131/18,542) to 3.04% at p < 0.05 (129/4,248), 6.55% at p < 0.01, 8.67% at p < 0.005, 14.81% at p < 0.001, 17.49% at p < 0.0005, 24.59% at p < 0.0001, and 25.59% at p < 1e-5 (254 pairs; Fig 2). Coloc support therefore rises monotonically with instrumental evidence, from a baseline of 0.42% across all 31,371 QC-passed pairs (131 strong) to 25.6% at the strictest threshold. The same calibration holds within each outcome (per-outcome FDR-core yield: T2D 16.5%, CAD 9.4%, FG 16.7%; Fig 2). Tightening the preregistered per-outcome BH-FDR threshold also raises yield monotonically, from 12.32% at q < 0.05 to 27.18% at q < 5e-4 and 27.01% at q < 1e-4 (137 pairs; Fig 2). The primary multiple-testing control and the nominal funnel thus agree on the direction of calibration. The nominal set (n = 4,248) serves as a broader, lower-confidence tier for hypothesis generation (Section 4).

### 3.3 A decisive negative boundary outside the MR-significant set

Strong colocalization was essentially absent among pairs failing nominal MR significance (p ≥ 0.05; 27,123 pairs): 2/14,294 pairs with 0.05 ≤ p < 0.5 (0.014%) and 0/12,829 pairs with p ≥ 0.5, for 2/27,123 pairs overall. A Poisson one-sided 95% upper bound on this rate is 0.0232% per pair; the two-sided 95% CI upper limit is ≈0.027%. A stratified random sample of 6,000 pairs from this set contained no strong colocalization. (The eight nominally significant but FDR-excluded strong pairs of Section 3.5 are a separate tier, not counted here; they constitute 8/3,266 = 0.25% of the nominally significant pairs excluded by the FDR filter.) This negative boundary establishes that strong colocalization in our pipeline is restricted to MR-supported hypotheses and is not an artifact of the coloc prior configuration (Fig 3).

### 3.4 Reproduction of previously reported strong colocalizations

All 106 loci that previously reached strong colocalization in our grid analysis were reproduced with the full-scan pipeline (106/106, 100%). After build correction, 86/106 (81%) fell within GWAS Catalog-reported annotations (26 by direct SNP match, 59 within ±100 kb of a reported gene, 1 by gene annotation). SMR+HEIDI on these loci passed in 76/106 (71.7%; 95% CI 62.5–79.4%), and Steiger directionality was consistent (eQTL → outcome) in 73/76 (96.1%; 95% CI 89.0–98.6%) of HEIDI-passing pairs.

### 3.5 Fifteen additional candidate effector genes

Outside the previously reported set, the FDR-core scan identified 15 additional strong-colocalization candidates (7 in T2D, 8 in CAD; Table 2). Representative regional colocalization plots (RBM6×T2D, CNNM2×CAD, PLAUR×CAD, CD101×T2D, RIC8A×CAD, and the gray-zone pair AP3S2×T2D) are shown in Fig 4. These comprised 9 candidate effector genes at reported T2D/CAD risk loci: SLC12A3×CAD, CWF19L1×T2D, U6atac×T2D, CD101×T2D, RBM6×T2D, CNNM2×CAD, N4BP2L2×CAD, RIC8A×CAD (8 within 100 kb), and C2orf49×T2D (within 250 kb). The remaining 6 were weak-locus candidates without T2D/CAD catalog records: PLAUR×CAD, TAGLN2×CAD, VSIG8×CAD, PDCD6×T2D, CLEC3B×T2D, CCDC19×CAD (Fig 5). Eight additional strong pairs from the raw-significant set (LAMC1, TPD52, SENP6, HMGN3, MT3, RPL13, ZBTB46, ZNF100) reached strong colocalization at nominal MR significance but were excluded by the per-outcome FDR filter (padj ≥ 0.05; Table S3) and are not retained as candidates (Section 4). These are colocalization-supported candidate effector genes, not new GWAS loci; "candidate" reflects the absence of a current GWAS Catalog record and does not exclude prior literature.

Only 41/106 (38.7%) of the known strong-colocalization regions reach genome-wide significance (p < 5e-8) in their GWAS, and the median GWAS p across regions is 4.73e-7 (per outcome: T2D 24/58 = 41%, CAD 17/46 = 37%, FG 0/2 = 0%); 6/15 (40%) of the candidate effector genes also reach genome-wide significance (Fig 3). The majority of colocalized regions therefore do not qualify as new GWAS loci and must be interpreted as candidate effector genes supported by orthogonal evidence (MR + coloc + directional replication). For the candidates, this support comprises MR and PP.H4 evidence recovered mostly below genome-wide significance (Fig 5).

### 3.6 Directional replication in independent resources

For the 15 candidate effector genes, GTEx v8 supported the eQTLGen direction in 6/7 gene–trait pairs with an overlapping GTEx eQTL (85.7%; 95% CI 48.7–97.4%; one conflict, VSIG8×CAD). Across the 106 known loci, GTEx v8 direction was consistent with eQTLGen in 44/63 measurable loci (69.8%; 95% CI 57.6–79.8%; Fig 5). None of the 15 candidates is an FG pair (FG has no FinnGen analog). Among the 15, 9/15 were alignable to FinnGen R11, with 9/9 gene-level and 8/9 variant-level direction concordance, and 4 reached FinnGen p < 0.05: RBM6×T2D p = 6.9e-6, CNNM2×CAD p = 6.3e-4, CD101×T2D p = 1.2e-3, RIC8A×CAD p = 1.1e-2 (Fig 5). Of the six non-alignable candidates, five lead variants (U6atac, CWF19L1, N4BP2L2, SLC12A3, PLAUR) were not localized in FinnGen R11 and one (CLEC3B) lacked original summary statistics (Fig 5). All four FinnGen-significant candidates fell within the six-locus coloc.susie sensitivity set and were retained (Section 3.7). Re-estimating MR directions in pre-UKB outcome GWAS with minimal eQTLGen overlap (Methods 2.2) gave concordant directions for all 15 candidates (11 positive, 4 negative; 0 discordant; Table S4). Sample overlap therefore does not drive the candidate-set MR directions.

### 3.7 Sensitivity of strong-colocalization calls to the single-causal-variant assumption

Coloc.susie with SuSiE fine-mapping did not converge under external 1000 Genomes EUR LD for any of the six tested loci spanning different signal architectures (RBM6×T2D, CNNM2×CAD, PLAUR×CAD, CD101×T2D, RIC8A×CAD, LAMC1×CAD). This is a known limitation when summary statistics and reference panels are imperfectly matched. We therefore do not rely on SuSiE posteriors for primary inference and report them only as exploratory sensitivity evidence (Table S2). The one locus that showed multi-signal architecture, LAMC1×CAD (coloc.abf PP.H4 = 0.9139), was independently excluded by the FDR filter (padj ≥ 0.05): fine-mapping revealed eight independent eQTL credible sets versus three weak GWAS credible sets (maximum z = 4.32) with no shared variant. In multi-signal regions with weak GWAS evidence, coloc.abf can overstate PP.H4. The retained 15 candidates therefore do not depend on the single-causal-variant assumption for their support: they were selected on the primary coloc.abf analysis and, where testable, survived directional replication in GTEx and FinnGen.

### 3.8 Threshold and model sensitivity

Within the FDR-core set, the strong-colocalization count was monotonic in the PP.H4 threshold: 290 pairs (29.5%) at PP.H4 ≥ 0.5, 121 (12.3%) at ≥ 0.8, and 74 (7.5%) at ≥ 0.9. Among pairs failing nominal MR significance (p ≥ 0.05), 11 reached PP.H4 ≥ 0.5 and 2 reached ≥ 0.8 (both MR-null; Section 3.3). Setting the coloc prior to p12 = 1e-6 yielded 20/106 grid loci (data not shown). A permutation-based false-positive estimate in the MR-significant subset was 1.45% (106 loci × 100 permutations; PP.H4 ≥ 0.8); we did not calibrate this to the full scan and therefore report it as an upper-bound indicator (Discussion).

---

## 4. Discussion

We present a transcriptome-wide cis-MR × coloc atlas for T2D, CAD, and FG that quantifies the operating characteristics of the standard gene-prioritization pipeline. We focus on three findings. First, coloc yield calibrates monotonically with MR evidence (0.42% across all QC-passed pairs → 25.6% across the MR p-value funnel). Two independent filtering strategies, per-outcome BH-FDR (12.32%) and LD clumping plus instrumental replication (12.96%), converge on a shared-causal-variant rate of approximately 12–13% among MR-supported pairs. Second, strong colocalization is essentially absent outside the MR-significant set (2/27,123 pairs; one-sided 95% upper bound 0.0232%/pair), providing a decisive negative boundary that constrains false-positive rates for a well-specified coloc prior. Finally, all 106 previously reported strong colocalizations reproduced, and 15 additional candidate effector genes survived directional replication in GTEx and FinnGen, with HEIDI and Steiger support in the known-locus set.

The coloc-yield funnel provides an interpretable calibration for investigators applying cis-MR × coloc. At MR p values near 0.05, coloc support carries only ~3% yield; yield exceeds 25% only at MR p < 1e-5. FDR control or clumped, replicated instruments raise the effective yield to ~12–13%. The 15 candidates add to the effector-gene annotation of T2D and CAD risk loci. Several of these (e.g., RBM6×T2D, CNNM2×CAD) have GWAS-significant peaks and FinnGen replication. Independent functional validation would be required before prioritizing them experimentally.

We suggest a tiered interpretation of the atlas. Tier 1 (FDR-core, 121 pairs) comprises high-confidence effector genes for follow-up functional studies. Tier 2 (nominally significant strong pairs, 129 pairs) are hypothesis-generating candidates requiring independent replication. Tier 3 (PP.H4 ≥ 0.5 but < 0.8) carries weak colocalization support and is not recommended for experimental prioritization without additional evidence. This tiering is encoded in the public atlas (Zenodo: https://doi.org/10.5281/zenodo.21967917) to facilitate community use.

Eight additional loci reached strong colocalization at nominal MR significance but not after per-outcome FDR control (Table S3). These include LAMC1×CAD, which reached nominal significance in independent FinnGen data (p = 0.018), and SENP6×CAD and HMGN3×CAD, which showed directionally concordant FinnGen estimates. Their exclusion illustrates the conservative nature of genome-wide FDR control in this setting. The per-outcome FDR threshold (q < 0.05) retains 982 of 4,248 (23%) nominally significant pairs and with it discards some colocalized signals along with most weak MR evidence. Investigators using this atlas for hypothesis generation may wish to inspect the nominal set (n = 4,248) as a broader, lower-confidence tier.

LAMC1×CAD illustrates a methodological hazard worth making explicit. Under the single-causal-variant assumption of coloc.abf it returned PP.H4 = 0.91, a superficially strong colocalization; SuSiE fine-mapping instead resolved eight independent eQTL credible sets against three weak GWAS credible sets (maximum z = 4.32) with no shared variant, and the pair was independently excluded by the per-outcome FDR filter (padj ≥ 0.05). In multi-signal regions with weak GWAS evidence, coloc.abf can therefore overstate PP.H4 even when no single causal variant is shared. We treat this case as a cautionary example of why the preregistered FDR control and multi-signal inspection accompany coloc.abf here; the retained 15 candidates do not depend on the single-causal-variant assumption for their support (Section 3.7).

We acknowledge the following limitations. First, sample overlap between the exposure (eQTLGen) and the outcome GWAS can bias MR estimates [@burgess2016]. The effect is expected to be modest here because overlap is partial and the primary inference (coloc PP.H4) depends on regional LD patterns rather than on the MR point estimate alone (Methods 2.2). We tested this directly for the 15 candidate effector genes using pre-UKB outcome GWAS with minimal overlap (Methods 2.2; Table S4). All 15 showed MR directions concordant with the primary analysis (11 positive, 4 negative; 0 discordant), arguing against overlap-driven direction bias in the candidate set. Overlap-driven distortion of colocalization support across the full scan cannot, however, be excluded. Second, only 38.7% of known strong-colocalization regions reach GWAS significance (median GWAS p = 4.73e-7). Most colocalized regions are candidate effector genes supported by converging evidence rather than definitive new loci. Third, the permutation false-positive rate (1.45%) was estimated only in the MR-significant subset and was not scaled genome-wide; we treat it as an upper-bound indicator. Fourth, coloc.susie did not converge under external LD and is treated as exploratory sensitivity evidence. The one locus it flagged (LAMC1×CAD) is independently excluded by the FDR filter, but non-convergence limits how strongly we can assert robustness for the retained candidates. Fifth, MR with single instruments is mathematically equivalent to SMR [@zhu2016]; the two should not be treated as independent validation. Finally, whole-blood eQTLs capture only the transcript in blood; tissue-specific effector genes (e.g., pancreatic, vascular) may be missed.

In conclusion, cis-MR × coloc is a calibrated and conservative tool for effector-gene prioritization. We provide a public atlas of 121 MR-significant strong-colocalization pairs (106 replicated + 15 candidate effector genes) with directional replication for T2D, CAD, and FG.

---

## Data and Code Availability

All analyses used publicly available summary statistics: eQTLGen whole-blood cis-eQTLs (n = 31,684) [@vosa2021]; T2D GWAS (Xue 2018, GCST006867; n = 655,666, 61,714 cases; per-variant N ≈ 573,704) [@xue2018]; CAD GWAS (van der Harst 2018, GCST005194; n = 296,525) [@vanderharst2018]; FG GWAS (Manning 2012, GCST005186; n = 58,074) [@manning2012]; GTEx v8 [@gtex2020]; FinnGen R11 [@kurki2023]; and 1000 Genomes Phase 3 EUR as the LD reference [@1kg2015]. Colocalization support estimates, instrument lists, replication statistics, and the analysis scripts are provided as a public atlas (Zenodo: https://doi.org/10.5281/zenodo.21967917). Analysis code is available at https://github.com/qgeng1465/dual-channel-mr-atlas.

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

**Figure 1. Study design and genome-wide distribution.** (A) Transcriptome-wide cis-MR and colocalization across three outcomes. Whole-blood cis-eQTL instruments (eQTLGen, n = 31,684) were scanned against T2D, CAD, and FG GWAS over 31,371 gene–trait pairs. Per-outcome BH-FDR control (q < 0.05) retained 982 MR-significant pairs, of which 121 reached strong colocalization (PP.H4 ≥ 0.8). Two further pairs in the MR-null gray zone also reached strong colocalization. The 121 were classified into 106 previously reported and 15 candidate effector genes, the latter subjected to directional replication in GTEx v8 and FinnGen R11. (B) Genome-wide distribution of the 121 strong-colocalization pairs (circle, previously reported locus; triangle, candidate effector gene; color, outcome: T2D orange, CAD blue, FG green). The two gray-zone pairs (AP3S2×T2D, ZNF19×CAD) are described in the text.

**Figure 2. Coloc yield calibrates with MR evidence.** (A) The coloc-yield funnel: strong-colocalization rate (PP.H4 ≥ 0.8, % of pairs) as a function of the nominal MR p-value threshold (grey curve; n per stratum annotated; Wilson 95% confidence intervals). The preregistered per-outcome BH-FDR core point (982 pairs; 12.32%) is the large filled red circle; the stage-2 grid point (818 evaluable; 12.96%) is the open dashed reference. (B) Per-outcome yield funnels over the same thresholds (T2D orange, CAD blue, FG green), with per-outcome FDR-core points overlaid (T2D 16.5%, CAD 9.4%, FG 16.7%). The pattern is consistent within each outcome; the sparse FG curve (n = 12 at the FDR core) is the least informative. (C) Yield under the preregistered per-outcome BH-FDR control as a function of the q threshold: 12.32% at q < 0.05 (982 pairs), rising monotonically to 27.18% at q < 5e-4 (206 pairs) and 27.01% at q < 1e-4 (137 pairs). The primary multiple-testing control itself calibrates with instrumental evidence. (D) Per-outcome FDR-core yields with Wilson 95% confidence intervals (T2D 65/394 = 16.5%, CAD 54/576 = 9.4%, FG 2/12 = 16.7%) against the combined yield (dashed line, 12.32%).

**Figure 3. Calibration disclosures: PP.H4 distribution, negative boundary, and GWAS significance.** (A) Empirical cumulative distribution of PP.H4 by MR-significance layer (FDR-core, n = 982; MR-null, 0.05 ≤ p < 0.5, n = 14,294; MR-negative, p ≥ 0.5, n = 12,829). (B) Strong-colocalization rate (PP.H4 ≥ 0.8) by layer: 121/982 (12.32%) in the FDR-core versus 2/14,294 (0.014%) and 0/12,829. Only 2/27,123 pairs outside the nominal MR-significant set reach strong colocalization (one-sided Poisson 95% upper bound 0.0232% per pair). Nominally significant but FDR-excluded pairs (8/3,266 = 0.25%) are a separate tier (Table S3). (C) PP.H4 threshold sensitivity within the FDR-core set (strong count at PP.H4 ≥ 0.5 / 0.8 / 0.9). (D) Permutation-based false-positive rate (106 loci × 100 permutations) with the observed FDR-core rate (12.3%) for reference. (E) Fraction of strong-colocalization pairs reaching genome-wide significance (p < 5e-8): known loci 41/106 (38.7%), candidate effector genes 6/15 (40%). (F) Distribution of GWAS peak p-values (median known-locus peak p = 4.73e-7). Strong colocalization support and GWAS significance are orthogonal; most colocalized regions are effector-gene candidates, not new GWAS loci.

**Figure 4. Regional colocalization plots for representative loci.** LocusZoom-style mirror plots of eQTL (top) and GWAS (bottom) association at selected loci: RBM6×T2D, CNNM2×CAD, PLAUR×CAD, CD101×T2D, RIC8A×CAD, and the gray-zone pair AP3S2×T2D. LD r² is shown relative to the lead instrument. Loci were selected as the five coloc.susie-sensitivity candidates retained as candidate effector genes (Section 3.7; the sixth tested locus, LAMC1×CAD, is shown in Supplemental Figure S1) plus the gray-zone pair.

**Figure 5. Candidate effector genes, convergence, and independent replication.** (A) Lollipop plot of PP.H4 for the 15 FDR-core strong-colocalization candidates outside the previously reported set (star, replicated at FinnGen p < 0.05; color, outcome). (B) cis-MR −log₁₀(p) versus GWAS peak −log₁₀(p): the 15 candidates (triangles) are often recovered below genome-wide significance (dashed line), unlike the 106 known loci (circles). (C) Independent-replication matrix: columns 1–3, normalized instrument strength (F), MR −log₁₀(p), and PP.H4 (viridis); columns 4–6, GTEx direction, FinnGen gene-level direction, and FinnGen p < 0.05 (green = supported, red = conflicting, grey = not measurable). (D) FinnGen R11 −log₁₀(p) versus original GWAS peak −log₁₀(p) for the 9 alignable candidates (star, FinnGen p < 0.05); 9/9 gene-level and 8/9 variant-level direction concordance, 4 at FinnGen p < 0.05. (E) FinnGen R11 alignment coverage: 9/15 alignable, 5 lead variants not localized, 1 without original summary data. (F) Convergence and direction: SMR+HEIDI pass rate in the MR-significant pre-scan and the strong-colocalized subset (76/106 = 71.7%); Steiger eQTL → outcome (73/76 = 96.1%; 0 reverse-significant); GTEx v8 direction for the 15 candidates (6/7) and 106 known loci (44/63). Error bars, Wilson 95% CI.

**Supplemental Figure S1. coloc.susie sensitivity of strong-colocalization calls.** Paired coloc.abf vs coloc.susie PP.H4 for the six adjudicated loci (RBM6×T2D, CNNM2×CAD, PLAUR×CAD, CD101×T2D, RIC8A×CAD, LAMC1×CAD). SuSiE credible-set counts per side and non-convergence markers (✗, all runs) are annotated. LAMC1×CAD (abf PP.H4 = 0.9139 → SuSiE 0.0000) is highlighted; it is excluded from the candidate set on independent FDR and multi-signal grounds. Because coloc.susie did not converge under external LD for any tested locus, SuSiE posteriors are reported as exploratory only and are not used for primary inference.

---

## Tables

**Table 1. Data sources and coloc-yield summary.**

| Role | Resource | n | Description |
|---|---|---|---|
| Exposure | eQTLGen whole-blood cis-eQTL | 31,684 | cis-eQTLs in whole blood [@vosa2021] |
| Outcome T2D | Xue et al., Nat Commun 2018 | 655,666 (61,714 cases); per-variant N ≈ 573,704 | T2D GWAS meta-analysis [@xue2018] |
| Outcome CAD | van der Harst & Verweij, Circ Res 2018 | 296,525 (34,541 cases) | CAD GWAS meta-analysis [@vanderharst2018] |
| Outcome FG | Manning et al., Nat Genet 2012 | 58,074 | Fasting glucose GWAS (quantitative) [@manning2012] |
| Replication | GTEx v8 | 838 | Multi-tissue eQTLs, 49 tissues [@gtex2020] |
| Replication | FinnGen R11 | n/a | Disease cohort for directional check [@kurki2023] |
| LD reference | 1000 Genomes Phase 3 EUR | 503 | Instrument clumping and SuSiE LD [@1kg2015] |
| **Analysis stratum** | **Pairs** | **Strong coloc (PP.H4 ≥ 0.8)** | **Yield (95% CI)** |
| All QC-passed pairs | 31,371 | 131 | 0.42% (0.35–0.50%) |
| Nominal MR p < 0.5 | 18,542 | 131 | 0.71% (0.60–0.84%) |
| Nominal MR p < 0.05 | 4,248 | 129 | 3.04% (2.56–3.60%) |
| Per-outcome FDR core (q < 0.05) | 982 | 121 | 12.32% (10.41–14.53%) |
| Stage-2 grid (clumped + replicated) | 818 | 106 | 12.96% (10.83–15.43%) |

**Table 2. Fifteen candidate effector genes.** Full annotation of the 15 FDR-core strong-colocalization candidates outside the previously reported set. Tier: candidate effector gene at a reported T2D/CAD risk locus (distance to the locus in parentheses) versus weak-locus candidate without T2D/CAD catalog records. GTEx: same-direction eQTL support (consistent / conflicting / not testable). FinnGen p: per-gene association in FinnGen R11 (–, not alignable or no analog). Released with the Zenodo deposition (https://doi.org/10.5281/zenodo.21967917) as `Table_2_candidates.csv`.

| Gene | Outcome | Tier | PP.H4 | MR p | BH-FDR q | GWAS peak p | GTEx direction | FinnGen p |
|---|---|---|---|---|---|---|---|---|
| C2orf49 | T2D | known locus (250 kb) | 0.9802 | 1.90e-7 | 5.10e-5 | 1.22e-7 | consistent | 0.136 |
| CWF19L1 | T2D | known locus (100 kb) | 0.9600 | 3.15e-4 | 0.0138 | 1.23e-9 | consistent | – |
| U6atac | T2D | known locus (100 kb) | 0.9571 | 1.71e-11 | 1.45e-8 | 1.97e-17 | not testable | – |
| CD101 | T2D | known locus (100 kb) | 0.9488 | 1.16e-4 | 0.00652 | 2.47e-8 | not testable | 0.00115 |
| RBM6 | T2D | known locus (100 kb) | 0.9448 | 1.02e-5 | 0.00114 | 9.60e-7 | consistent | 6.91e-6 |
| PDCD6 | T2D | weak candidate | 0.8423 | 3.11e-5 | 0.00246 | 1.64e-5 | not testable | 0.247 |
| CLEC3B | T2D | weak candidate | 0.8338 | 1.32e-4 | 0.00715 | 5.07e-5 | not testable | – |
| PLAUR | CAD | weak candidate | 0.9957 | 2.00e-6 | 2.60e-4 | 1.81e-6 | not testable | – |
| SLC12A3 | CAD | known locus (100 kb) | 0.9937 | 1.97e-9 | 7.42e-7 | 1.19e-9 | consistent | – |
| TAGLN2 | CAD | weak candidate | 0.9465 | 2.63e-5 | 2.03e-3 | 2.35e-5 | not testable | 0.158 |
| CCDC19 | CAD | weak candidate | 0.9462 | 2.63e-5 | 2.03e-3 | 2.35e-5 | not testable | 0.158 |
| CNNM2 | CAD | known locus (100 kb) | 0.9366 | 9.85e-15 | 2.89e-11 | 4.67e-15 | not testable | 6.31e-4 |
| N4BP2L2 | CAD | known locus (100 kb) | 0.9237 | 1.18e-9 | 4.67e-7 | 1.19e-10 | consistent | – |
| VSIG8 | CAD | weak candidate | 0.9162 | 2.63e-5 | 2.03e-3 | 2.35e-5 | conflicting | 0.158 |
| RIC8A | CAD | known locus (100 kb) | 0.8895 | 4.88e-5 | 3.34e-3 | 8.35e-6 | consistent | 0.0114 |

---

## Supplemental Items

**Table S1.** Outcome-stratified scan counts. For each outcome, the number of QC-passed gene–trait pairs, pairs passing the preregistered per-outcome BH-FDR multiple-testing control (q < 0.05), MR-null (0.05 ≤ p < 0.5) and MR-negative (p ≥ 0.5) strata, and strong colocalization (PP.H4 ≥ 0.8) counts within the FDR-core, the MR-null gray zone, and the MR-negative range.

| Outcome | QC-passed pairs | MR significant (q < 0.05) | MR null (0.05–0.5) | MR negative (≥ 0.5) | Strong in FDR-core | Gray-zone strong | Negative strong |
|---|---|---|---|---|---|---|---|
| T2D | 10,190 | 394 | 4,579 | 4,054 | 65 | 1 | 0 |
| CAD | 14,673 | 576 | 6,716 | 5,715 | 54 | 1 | 0 |
| FG | 6,508 | 12 | 2,999 | 3,060 | 2 | 0 | 0 |
| Total | 31,371 | 982 | 14,294 | 12,829 | 121 | 2 | 0 |

**Table S2.** coloc.susie diagnostics for the six-locus sensitivity set. Pairwise coloc.abf versus coloc.susie (SuSiE fine-mapping under a multiple-causal-variant assumption) PP.H4, number of SNPs in the regional window, SuSiE credible-set counts on the eQTL and GWAS sides, and convergence status under external 1000 Genomes EUR LD (max_iter = 200). None converged; SuSiE posteriors are exploratory only. LAMC1×CAD (coloc.abf PP.H4 = 0.9139 → SuSiE 0.0000) was excluded from the candidate set on independent FDR and multi-signal grounds.

| Gene | Outcome | MR p | coloc.abf PP.H4 | coloc.susie PP.H4 | n SNPs | SuSiE eQTL CS | SuSiE GWAS CS | Convergence |
|---|---|---|---|---|---|---|---|---|
| RBM6 | T2D | 1.02e-5 | 0.9448 | 1.0000 | 927 | 2 | 2 | not converged |
| CNNM2 | CAD | 9.85e-15 | 0.9366 | 0.9996 | 2,557 | 10 | 6 | not converged |
| PLAUR | CAD | 2.00e-6 | 0.9957 | 0.9990 | 3,236 | 10 | 9 | not converged |
| CD101 | T2D | 1.16e-4 | 0.9488 | 1.0000 | 1,653 | 10 | 8 | not converged |
| RIC8A | CAD | 4.88e-5 | 0.8895 | 1.0000 | 2,044 | 3 | 3 | not converged |
| LAMC1 | CAD | 4.73e-3 | 0.9139 | 0.0000 | 3,063 | 8 | 3 | not converged |

**Table S3.** Eight additional strong pairs excluded by per-outcome FDR control. These loci reached strong colocalization (PP.H4 ≥ 0.8) at nominal MR significance but did not survive the preregistered per-outcome BH-FDR filter (q ≥ 0.05) and are therefore not retained as candidates. GTEx direction and FinnGen replication status are shown for completeness.

| Gene | Outcome | MR p | FDR q | PP.H4 | GWAS peak p | GTEx direction | FinnGen |
|---|---|---|---|---|---|---|---|
| LAMC1 | CAD | 4.73e-3 | 0.0865 | 0.9139 | 1.29e-5 | consistent | replicated (p = 0.018) |
| TPD52 | FG | 1.26e-3 | 0.2404 | 0.9128 | 3.14e-5 | not testable | no FG phenotype in FinnGen |
| SENP6 | CAD | 4.85e-3 | 0.0877 | 0.9056 | 4.58e-5 | consistent | directionally concordant |
| HMGN3 | CAD | 1.26e-2 | 0.1588 | 0.8884 | 5.23e-5 | consistent | directionally concordant |
| MT3 | FG | 1.20e-4 | 0.0601 | 0.8800 | 9.17e-5 | not testable | no FG phenotype in FinnGen |
| RPL13 | T2D | 1.14e-2 | 0.1499 | 0.8702 | 5.25e-6 | not testable | original GWAS sumstats missing |
| ZBTB46 | CAD | 1.07e-2 | 0.1453 | 0.8442 | 5.56e-7 | consistent | lead not alignable |
| ZNF100 | CAD | 2.34e-2 | 0.2245 | 0.8177 | 3.97e-5 | consistent | lead not alignable |

**Table S4.** Sample-overlap sensitivity. MR direction of each candidate effector gene was re-estimated in outcome GWAS with minimal or no overlap with the eQTLGen exposure: the entirely pre-UK Biobank CARDIoGRAMplusC4D meta-analysis for CAD (Nikpay et al. 2015) and the DIAGRAM trans-ethnic meta-analysis for T2D (Mahajan et al. 2014). Three genes absent from the trans-ethnic meta (U6atac, PDCD6, CLEC3B) were checked in a European 1000 Genomes–based secondary analysis (Scott et al. 2017). All 15 candidates (8/8 CAD, 7/7 T2D) showed MR directions concordant with the primary analysis (11 positive, 4 negative; 0 discordant). Alt p, alternative-GWAS p-value; Alt dir, alternative-GWAS MR direction (sign of beta); Concordant, direction agreement with the primary analysis.

| Gene | Outcome | Lead SNP | Alternative GWAS | Beta | Alt p | Primary dir | Alt dir | Concordant |
|---|---|---|---|---|---|---|---|---|
| C2orf49 | T2D | rs4851758 | DIAGRAM trans-ethnic | 0.00995 | 0.81 | + | + | yes |
| CD101 | T2D | rs12130298 | DIAGRAM trans-ethnic | 0.0677 | 0.38 | + | + | yes |
| CLEC3B | T2D | rs77100387 | Scott 2017 (European 1000G) | 0.022 | 0.072 | + | + | yes |
| CWF19L1 | T2D | rs17668357 | DIAGRAM trans-ethnic | 0.0392 | 9.2e-3 | − | − | yes |
| PDCD6 | T2D | rs12332382 | Scott 2017 (European 1000G) | 0.018 | 0.35 | + | + | yes |
| RBM6 | T2D | rs2245365 | DIAGRAM trans-ethnic | 0.00995 | 0.31 | − | − | yes |
| U6atac | T2D | rs12780155 | Scott 2017 (European 1000G) | 0.059 | 5.4e-5 | + | + | yes |
| CCDC19 | CAD | rs2789422 | CARDIoGRAMplusC4D (pre-UKB) | −0.0473 | 2.46e-6 | + | + | yes |
| CNNM2 | CAD | rs11191447 | CARDIoGRAMplusC4D (pre-UKB) | −0.0744 | 1.32e-7 | + | + | yes |
| N4BP2L2 | CAD | rs1123462 | CARDIoGRAMplusC4D (pre-UKB) | −0.0368 | 3.2e-4 | − | − | yes |
| PLAUR | CAD | rs4760 | CARDIoGRAMplusC4D (pre-UKB) | 0.0497 | 2.3e-3 | + | + | yes |
| RIC8A | CAD | rs6598075 | CARDIoGRAMplusC4D (pre-UKB) | −0.0232 | 3.3e-2 | − | − | yes |
| SLC12A3 | CAD | rs56228609 | CARDIoGRAMplusC4D (pre-UKB) | −0.0316 | 2.6e-3 | + | + | yes |
| TAGLN2 | CAD | rs2789422 | CARDIoGRAMplusC4D (pre-UKB) | −0.0473 | 2.46e-6 | + | + | yes |
| VSIG8 | CAD | rs2789422 | CARDIoGRAMplusC4D (pre-UKB) | −0.0473 | 2.46e-6 | + | + | yes |

**Supplemental Figure S1.** coloc.susie sensitivity of strong-colocalization calls. Paired coloc.abf vs coloc.susie PP.H4 for the six adjudicated loci (RBM6×T2D, CNNM2×CAD, PLAUR×CAD, CD101×T2D, RIC8A×CAD, LAMC1×CAD). SuSiE credible-set counts per side and non-convergence markers (✗, all runs) are annotated. LAMC1×CAD (abf PP.H4 = 0.9139 → SuSiE 0.0000) is highlighted; it is excluded from the candidate set on independent FDR and multi-signal grounds. Because coloc.susie did not converge under external LD for any tested locus, SuSiE posteriors are reported as exploratory only and are not used for primary inference.
