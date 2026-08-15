# CodonSafe Dataset Manifest
**Last verified: 2026-04-13**
**Total on disk: 103 MB across 9 studies**

## Status Legend
- [x] Downloaded, verified as real data (not HTML redirect)
- [~] PDF supplement only — data tables require extraction
- [ ] Not downloaded / not applicable

---

## 1. Napolitano et al. 2016 (PNAS) — AGR recoding + CRAM
**DOI:** 10.1073/pnas.1605856113 | **PMCID:** PMC5035903
**Dir:** `napolitano2016/` (2.1 MB)

| File | Rows | Cols | Content | Status |
|------|------|------|---------|--------|
| sd01.xlsx | 124 | 13 | All 123 AGR positions, recalcitrance | [x] |
| sd02.xlsx | 12 | 6 | Recalcitrant codon doubling costs | [x] |
| sd03.xlsx | 10 | 6 | mRNA folding energy statistics | [x] |
| sd04.xls | 418 | 11 | RBS strength per position (old Excel format) | [x] |
| sd05.xlsx | 4 | 5 | tRNA kinetics (ArgU) | [x] |
| sd06.xlsx | 16 | 5 | CRAM summary (16 positions, convergence data) | [x] |
| sd07.xlsx | **3,222** | **13** | **All AGR positions with mRNA/RBS SRZ covariates** | [x] |
| sd08.xlsx | 123 | 2 | MAGE oligo sequences | [x] |
| SI PDF | — | — | Supplementary methods and figures | [x] |

**Key for CodonSafe:** sd07 is the main table (3,222 AGR positions across 1,676 genes with mRNA structure + RBS strength deviations per CGN alternative). sd06 has CRAM convergence data for 16 positions. The full per-codon CRAM fitness time-series (all 64 alternatives at 14 loci) was deposited in BioProject — may need SRA retrieval for the raw tracking data.

---

## 2. Fredens et al. 2019 (Nature) — Syn61
**DOI:** 10.1038/s41586-019-1192-5 | **GenBank:** CP040347.1
**Dir:** `fredens2019/` (30 MB)

| File | Rows | Content | Status |
|------|------|---------|--------|
| Sup Data 1 - MDS42 genome.txt | — | Parent genome GenBank | [x] |
| Sup Data 2 - Designed MDS42 Fragment 1-37b.txt | — | Designed genome GenBank | [x] |
| **Sup Data 3 - Table Target codons.xlsx** | **10,393** | **All target codon positions (TCA/TCG/TAG)** | [x] |
| Sup Data 4 - Table Overlaps and refactoring.xlsx | — | Gene overlap refactoring | [x] |
| Sup Data 5 - Table MDS42 10kb stretches.xlsx | — | 10kb stretch coordinates | [x] |
| Sup Data 9 - Table BAC construction.xlsx | — | BAC assembly details | [x] |
| Sup Data 10 - Table BAC assemblies.xlsx | — | BAC assembly results | [x] |
| Sup Data 11 - Table REXER.xlsx | — | REXER experiment results | [x] |
| Sup Data 14 - Table oligos for fixing.xlsx | — | Fix oligonucleotides | [x] |
| Sup Data 18 - Syn61 genome.txt | — | Final Syn61 genome GenBank | [x] |
| **Sup Data 19 - Design optimisations and mutations.xlsx** | **21** | **7 design corrections + non-programmed mutations** | [x] |
| Sup Data 20 - Proteomics.xlsx | — | Label-free quantification | [x] |
| Data guide.pdf | — | Column descriptions | [x] |
| + 8 more GenBank/primer files | — | Supporting data | [x] |

**Key for CodonSafe:** Data 3 has all 18,218 target codon positions (split across TCA=7,587, TCG=10,392, TAG=239 columns). Data 19 has the 7 positions that required fixes — these are the topology-relevant failure cases. All Ser swaps cross the disconnection boundary (UCN→AGY).

---

## 3. Ostrov et al. 2016 (Science) — 57-codon design
**DOI:** 10.1126/science.aaf3639
**Dir:** `ostrov2016/` (15 MB)

| File | Content | Status |
|------|---------|--------|
| table_s3.xlsx | Segment-level viability (55 segments) | [x] |
| table_s4.xlsx | 13 lethal design exceptions | [x] |
| table_s5.xlsx | Additional segment data | [x] |
| aaf3639_Ostrov_rE.coli-57.Design.gb | Designed genome GenBank (14 MB) | [x] |
| ostrov_sm.pdf | Supplementary methods | [x] |

**Key for CodonSafe:** Table S4 has the 13 lethal exceptions. The GenBank file has all 62,214 recoded codon positions extractable by comparing to the parent MDS42 genome.

---

## 4. Robertson et al. 2025 (Science) — Syn57
**DOI:** 10.1126/science.ady4368 | **bioRxiv:** 10.1101/2025.05.02.651837
**Dir:** `robertson2025_syn57/` (41 MB)

| File | Rows | Content | Status |
|------|------|---------|--------|
| **Data S1 - Recoding schemes for 20kb REXERs.xlsx** | **37** | **All tested recoding schemes + efficiency** | [x] |
| **Data S2 - vS33A7 genome design.gb** | — | **Designed 57-codon genome (GenBank)** | [x] |
| Data S3 - MDS42 genome annotation updates.xlsx | — | Annotation corrections | [x] |
| Data S4 - uREXER_BAC.gb | — | BAC vector GenBank | [x] |
| **Data S5 - Summary of fixes.xlsx** | **11** | **Fragments requiring troubleshooting** | [x] |
| Data S6 - Summary of conjugations.xlsx | — | Assembly details | [x] |
| Data S7 - 3Chi 100k36 BAC.gb | — | BAC GenBank | [x] |
| **Data S8 - Syn57 final genome.gb** | — | **Final Syn57 genome (GenBank)** | [x] |
| **Data S9 - Mutations in Syn57.xlsx** | **102** | **All deviations from design (programmed + spontaneous)** | [x] |
| **Data S10 - RNA-seq of Syn57.xlsx** | **3,753** | **Differential expression (log2FC, padj)** | [x] |
| robertson_supp.pdf | — | Supp Figs 1-29 (recoding landscapes in fig S20) | [x] |
| ady4368_captions_data_s1_to_s10.docx | — | Data file descriptions | [x] |

**Key for CodonSafe:** Data S2 (designed genome) vs Data S8 (final genome) gives every codon change. Data S5 has the 11 fix events. Data S9 has 102 deviations (101 programmed + 5 non-programmed). Data S10 gives transcriptomic impact. The recoding-fitness linkage map is in Supplementary Fig S20 in the PDF — will need visual/tabular extraction.

---

## 5. Grome/Robertson et al. 2025 (Nature) — Ochre
**DOI:** 10.1038/s41586-024-08501-x | **GenBank:** CP173712
**Dir:** `robertson2025/` (2.5 MB)

| File | Rows | Content | Status |
|------|------|---------|--------|
| MOESM1 (PDF) | — | Supplementary Information | [x] |
| MOESM4.xlsx | 22 | Supplementary Table data | [x] |
| MOESM5.xlsx | 144 | RBS predictions (Fig 1d / Supp Table 2) | [x] |
| MOESM6.xlsx | 250 | **RBS impact for 249 TGA→TAA overlapping genes** | [x] |
| MOESM7.xlsx | 55 | Growth curve replicates (Fig 2e-f) | [x] |
| MOESM8.xlsx | 10 | Fitness data (Fig 3) | [x] |
| MOESM9.xlsx | 47 | Mass spec peptide data (Fig 4c-d) | [x] |
| **MOESM10.xlsx** | **88** | **Growth data biological replicates (Fig 5)** | [x] |
| MOESM11.xlsx | 99 | OTS GFP fluorescence data (Fig 6) | [x] |
| **MOESM12.xlsx** | **88** | **Extended Data Fig 3 — dissected fitness effects** | [x] |
| MOESM13.xlsx | 4 | Liquid complementation (ED Fig 4c) | [x] |
| MOESM14.xlsx | 325 | Genomic data (ED Fig 5b-d) | [x] |
| **MOESM15.xlsx** | **218** | **Kinetic growth curves (LB)** | [x] |
| MOESM16.xlsx | 99 | OTS data (ED Fig 10) | [x] |

**Key for CodonSafe:** MOESM6 has RBS predictions for 249 overlapping genes affected by TGA→TAA. MOESM10/12/15 have the fitness data. The 1,195 individual TGA→TAA positions are derivable from the CP173712 genome vs parent C321.ΔA (CP006698).

---

## 6. Frumkin et al. 2018 (PNAS) — CGG recoding
**DOI:** 10.1073/pnas.1719375115
**Dir:** `frumkin2018/` (48 KB)

| File | Rows | Content | Status |
|------|------|---------|--------|
| sd01.xlsx | — | Gene recoding details (8 genes, 60 mutations) | [x] |
| sd02.xlsx | — | Translation efficiency data | [x] |

**Key for CodonSafe:** 60 CGU/CGC→CGG mutations in 8 highly expressed genes + proteome-wide translation efficiency changes. Ribosome profiling raw data in NCBI SRA.

---

## 7. Lajoie et al. 2013 (Science) — C321.ΔA
**DOI:** 10.1126/science.1241459 | **PMCID:** PMC4924538 | **GenBank:** CP006698
**Dir:** `lajoie2013/` (228 KB)

| File | Content | Status |
|------|---------|--------|
| Table_S16.xlsx | Strain genotypes | [x] |
| Table_S19.xlsx | Codon conversion data | [x] |
| Table_S26.xlsx | Growth/fitness data | [x] |

**Key for CodonSafe:** 321 UAG→UAA conversions. Baseline for Ochre analysis.

---

## 8. Ding et al. 2024 (Science) — Mammalian rare codon recoding
**DOI:** 10.1126/science.adm8143
**Dir:** `ding2024/` (3.6 MB)

| File | Content | Status |
|------|---------|--------|
| science.adm8143_data_s1.xlsx | Recoding position data | [x] |
| science.adm8143_data_s2.xlsx | ncAA incorporation data | [x] |
| science.adm8143_data_s3.xlsx | Additional data | [x] |
| science.adm8143_sm.pdf | Supplementary figs S1-S9 | [x] |

**Key for CodonSafe:** TCG codon recoding in mammalian cells — cross-organism validation.

---

## 9. Nyerges et al. 2024 (bioRxiv) — Ec_Syn57 multi-omics
**DOI:** 10.1101/2024.06.16.599206
**Dir:** `nyerges2024/` (9.5 MB)

| File | Content | Status |
|------|---------|--------|
| nyerges_supp.pdf | **Supplementary figs + methods (9.9 MB)** | [x] |

**Key for CodonSafe:** Multi-omics data (transcriptome, translatome, proteome) at 62,007 recoded positions. Data tables are in the PDF supplement — require extraction. Raw sequencing data in SRA (BioProject accession in paper, pending public release). The differential expression and Ribo-seq data are the key tables to extract.

---

## NOT Downloaded (excluded per user)

| Study | Dir | Reason |
|-------|-----|--------|
| Fang et al. 2026 (Nat. Chem.) | `fang2026/` | Excluded |
| Wannier et al. 2018 (PNAS) | `wannier2018/` | Excluded |

---

## Observation Counts (Estimated)

| Study | Codon-level obs | Key topology features |
|-------|----------------|-----------------------|
| Napolitano 2016 | 3,222 AGR + CRAM | Nonsynonymous alternatives, SRZ covariates |
| Fredens 2019 | 18,218 | Ser boundary crossing (all swaps) |
| Ostrov 2016 | 62,214 (from GenBank) | Ser, Leu, Arg topology |
| Robertson Syn57 2025 | ~101,302 (from GenBank) | **4 Ser + 2 Ala codons removed** |
| Ochre 2025 | 1,195 + 321 | Stop codon topology |
| Nyerges 2024 | 62,007 | Multi-omics at recoded sites |
| Frumkin 2018 | 60 | Arg neighborhood effects |
| Lajoie 2013 | 321 | Stop codon baseline |
| Ding 2024 | Variable | Cross-organism (mammalian) |
| **TOTAL** | **>248,000** | — |
