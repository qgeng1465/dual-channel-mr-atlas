# Drug Annotation Layer (P4) — 76 HEIDI-pass + coloc-strong priority genes

Date: 2026-08-13
Input: `tmp/priority_genes_heidi_pass.txt` (76 genes, transcription-channel cis-MR, coloc PP.H4>=0.8 and HEIDI p>0.05)
Output CSV: `results/drug_annotation_20260813.csv`
Raw API payloads (audit trail): `tmp/drug_annotate_raw.json` (Open Targets, all 76), `tmp/chembl_crosscheck.json` (ChEMBL cross-check)

## 1. Summary counts

| Metric | Count |
|---|---|
| Total genes | 76 |
| Protein-coding (biotype_ok=1) | **56** |
| Non-coding / pseudogene artifact (biotype_ok=0) | **20** |
| Genes with >=1 known drug (Open Targets) | **2** |
| Genes in any clinical phase (>=1) | **2** |
| Genes with approved/launched drugs (phase 4) | **2** |
| Genes with no known drug | **74** |

Outcome breakdown: t2d 39 (28 coding / 11 non-coding), cad 36 (27 / 9), fbg 1 (ISCA1, coding).

## 2. Genes with known drugs (Open Targets `drugAndClinicalCandidates`)

Only 2 of 76 priority genes are known drug targets:

| Gene | Outcome | PP.H4 | n drugs | max clinical phase | Drugs (API) |
|---|---|---|---|---|---|
| KCNJ11 | t2d | 0.897 | 15 | 4 (approved) | GLIPIZIDE; CHLORPROPAMIDE; REPAGLINIDE; MITIGLINIDE CALCIUM DIHYDRATE; GLIMEPIRIDE; ACETOHEXAMIDE; TOLAZAMIDE; MINOXIDIL; DIAZOXIDE; NATEGLINIDE; MITIGLINIDE; TOLBUTAMIDE; GLICLAZIDE; GLYBURIDE; PINACIDIL |
| BMPR1A | cad | 0.910 | 2 | 4 (approved) | EPTOTERMIN ALFA (recombinant BMP-7, Protein); DIBOTERMIN ALFA (recombinant BMP-2, Protein) |

- KCNJ11 = Kir6.2 K-ATP channel subunit; sulfonylureas/meglitinides are first-line insulin secretagogues. OT T2D association score 0.864 (strongest in the list).
- BMPR1A = BMP type I receptor; the two approved biologics are BMP-2/-7 ortholog agonists used in bone repair. No T2D/CAD association in OT (top associations: juvenile polyposis syndrome 0.826, bone disorder 0.462).

## 3. Genes with T2D / CAD associations (Open Targets `associatedDiseases`, API-returned)

Association scores are Open Targets target-disease association scores (0–1; strong >=~0.5, moderate 0.2–0.5, weak <0.2). All entries below are API-returned.

**T2D-associated (25 genes):** ARG1 (0.444), SH3BGRL3 (0.204), C18orf8/ILRUN (0.385), MED23 (0.464), PDGFC (0.395), PTPRN (0.029), ZNF236 (0.425), ARL13B (0.074), HERPUD2 (0.012), MED27 (0.348), LRIG1 (0.276), CDC123 (0.461 + CAD 0.307 + MI 0.307), SPATA5 (0.262), ZBTB6 (0.072), YTHDF2 (0.210), KCNJ11 (0.864), HMBS (0.329), DUSP13 (0.029 + CAD 0.020), ZNF268 (0.343), ZNF34 (0.086), RC3H2 (0.064), KDM5A (0.321), NTAN1 (0.168), NRBF2 (0.234), HSD17B12 (0.369 + CAD 0.236).

**CAD/MI-associated (8 genes):** SERBP1P3 (MI 0.004, pseudogene), DAGLB (CAD 0.292), TSPAN14 (MI 0.410; CAD 0.379), LIPA (MI 0.534; CAD 0.453), FCHO1 (MI 0.298; CAD 0.095), plus HSD17B12, CDC123, DUSP13 (also in T2D list above).

Most OT association scores are weak-to-moderate; only KCNJ11→T2D (0.864) and LIPA→MI/CAD (0.534/0.453) are notable. LIPA (lysosomal acid lipase) having a CAD/MI association is consistent with its established role in cholesterol handling.

## 4. Genes with no known drug (74)

All 74 remaining genes have zero drug-target entries in Open Targets `drugAndClinicalCandidates`. A ChEMBL cross-check (mechanism-of-action records, clinical-only max_phase>=1) was run for all 54 protein-coding genes without OT drugs and returned **0 clinical mechanisms** — including genes sometimes discussed in early-stage medicinal chemistry (ARG1, KDM5A, CHD4, DAGLB, LOXL4, SIK2, HMBS). These are therefore best treated as undrugged / not-yet-clinically-validated targets.

## 5. Non-coding / pseudogene artifacts (20 genes)

These are intergenic lncRNA, pseudogene or non-coding RNA loci from the cis-eQTL LD block; they have no drug-target meaning and are flagged biotype_ok=0:

`SERBP1P3` (processed_pseudogene), `LRRC37A15P` (processed_pseudogene), `RP11-10L12.2` (processed_pseudogene), `KRT8P46` (processed_pseudogene), `PHBP9` (processed_pseudogene), `RP11-332O19.3` (processed_pseudogene), `MAST4-AS1` (lncRNA), `RP11-464F9.1` (lncRNA), `RP11-856B14.1` (lncRNA), `AC002116.7` (lncRNA), `RP4-756G23.5` (lncRNA), `RP1-159A19.4` (lncRNA), `RP11-98D18.9` (lncRNA), `CTD-2037K23.2` (lncRNA), `AP000487.5` (lncRNA), `RP11-552F3.10` (lncRNA), `U6` (snRNA), `MIR4513` (miRNA), `RP11-464F9.9` (OT returned empty biotype; RP11-* locus, non-coding), `Metazoa_SRP` (OT returned empty biotype; SRP RNA, non-coding).

The two genes with empty OT biotype were classified non-coding by symbol/locus; raw payloads are in `tmp/drug_annotate_raw.json` for verification.

## 6. API / network status (honest report)

- **Open Targets GraphQL** (`api.platform.opentargets.org/api/v4/graphql`): success for all 76/76 genes, 0 failures. Biotype, drugs, disease associations and tractability all returned.
- **Pharos** (`pharos.nih.gov/api*` and GraphQL variants): **unreachable from this host** (HTTP 403 nginx / HTML app responses) — no TDL could be fetched. TDL column therefore uses an **estimation derived from Open Targets** tractability + clinical-stage flags, explicitly suffixed `(est)`, and left blank for non-coding genes. If a real Pharos TDL is required, rerun with network access to pharos.nih.gov.
- **ChEMBL** (`ebi.ac.uk/chembl/api/data`): reachable; used for the 54-gene cross-check (0 clinical mechanisms found) and as a corroborating second source.
- No keys required for any working API. No data was fabricated; drug names, phases and disease scores are verbatim API returns.

## 7. Literature-level caveats (NOT from API — provided for context, treat as common knowledge)

These are well-known public facts, explicitly separated from API output. They were NOT written into the CSV's API columns:

- **KCNJ11**: sulfonylureas (e.g., glibenclamide/glyburide, gliclazide) and meglitinides are licensed insulin secretagogues for T2D — consistent with the API result (API confirms: phase 4).
- **LIPA**: sebelipase alfa (Kanuma) is an approved enzyme-replacement product for lysosomal acid lipase deficiency; it **is** the recombinant LIPA protein itself, so it does not appear as a "drug targeting LIPA" in Open Targets/ChEMBL. Worth noting separately for the paper, but not in `drug_names` (API-sourced only).
- **ARG1**: arginase inhibitors (e.g., INCB001158/CB-1158) have entered early-phase oncology trials; neither Open Targets nor ChEMBL currently annotates a clinical mechanism for ARG1, so this is flagged here as literature knowledge only.
- **KDM5A / CHD4 / DAGLB / SIK2**: chromatin- and kinase-family targets with active tool-compound chemistry, but no approved or late-phase drug specific to these targets per the APIs queried.
