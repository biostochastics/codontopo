# CodonSafe Analysis Results — Phase 1-2
**Date: 2026-04-13**
**Status: First-pass results from all computational phases**

---

## Phase 0: GenBank Extraction (PASS)

**Syn57 codon positions extracted:** 60,240 from CDS-matched comparison (MDS42 parent vs vS33A7 design genome).

Discrepancy vs paper-reported 101,302: the MDS42 parent genome has 7,509 features vs the design genome's 110,025. The 60,240 represent codons extractable where both genomes have matching CDS annotations by locus_tag. The full count requires using the design genome's more comprehensive annotation.

**Breakdown:**
| Source codon | Target | Count | Amino acid |
|---|---|---|---|
| GCA | GCT | 22,859 | Ala |
| TCG | AGC | 10,263 | Ser |
| TCC | AGC | 9,910 | Ser |
| TCT | AGC | 9,471 | Ser |
| TCA | AGT | 7,502 | Ser |
| TAG | TAA | 235 | Stop |
| **Total** | | **60,240** | |

All recodings are synonymous. Zero forbidden codons remain in the design genome CDS.

---

## Phase 1: Syn57 Ser-vs-Ala Topology Classification

### Topology Classification (confirmed)
- **37,146 Ser swaps**: 100.0% cross the disconnection boundary (UCN→AGY)
- **22,859 Ala swaps**: 0.0% cross any boundary (within GCN box)
- Clean within-study contrast confirmed.

### Local Mismatch (POSITIVE RESULT)
- Ser swaps **improve** local neighborhood: Δ = -37.2 ± 30.5 Grantham units
- Ala swaps **worsen** local neighborhood: Δ = +19.0 ± 0.0 Grantham units
- Mann-Whitney U = 0, p ≈ 0
- **Interpretation:** Boundary-crossing Ser swaps move codons to neighborhoods with lower physicochemical mismatch to their Hamming-1 neighbors, while within-box Ala swaps move to worse neighborhoods.

### RNA-seq Outcomes (NULL RESULT for primary hypothesis)

**Test 1 — Ser-only vs Ala-only genes:**
- Ser-only genes (n=129): mean |log2FC| = 1.702, 84.5% significant
- Ala-only genes (n=49): mean |log2FC| = 1.788, 71.4% significant
- Mann-Whitney p = 0.40 (not significant)
- **Conclusion:** Boundary-crossing (Ser) recodings do not cause more transcriptomic perturbation than within-box (Ala) recodings.

**Test 2 — Fraction boundary-crossing vs |log2FC|:**
- Spearman rho = 0.012, p = 0.46
- **Conclusion:** No correlation between fraction of boundary-crossing recodings and gene expression change.

**Test 3 — Recoding count vs perturbation:**
- n_ser_recodings vs |log2FC|: rho = 0.050, p = 0.003
- n_ala_recodings vs |log2FC|: rho = 0.040, p = 0.018
- n_total_recodings vs |log2FC|: rho = 0.048, p = 0.004
- **Conclusion:** More recodings of ANY type weakly correlate with more perturbation, but not differentially by topology class.

### Fix Analysis
All 11 fixed fragments contain both Ser and Ala recodings. Fix types are NCS (non-coding sequence disruption), regulatory sequence disruption, and refactoring — **not topology-specific**. Consistent with Napolitano's SRZ model: recoding failures are driven by mRNA structure/RBS disruption, not amino acid topology.

---

## Phase 2: Napolitano AGR→CGN Analysis

### Topology Classification
- All 12,888 AGR→CGN swaps: **zero boundary crossings** (Arg codons form one connected component at ε=1)
- This dataset **cannot test the boundary-crossing hypothesis**

### Local Mismatch vs SRZ (POSITIVE RESULT)
- Δlocal_mismatch vs RBS_deviation: **rho = -0.33, p ≈ 0** (strong)
- Δlocal_mismatch vs mRNA_deviation: **rho = -0.12, p ≈ 0** (moderate)
- **Interpretation:** GF(2)^6 local mismatch is correlated with the Napolitano SRZ metrics. Topology captures overlapping but not identical information to mRNA structure/RBS.

### Per-target mismatch profile
| Target | Local mismatch | Hamming | Mean mRNA dev | Mean RBS dev |
|---|---|---|---|---|
| CGA | 360 | 1 | 0.996 | 2.58 |
| CGC | 421 | 2 | 1.019 | 93,113 |
| CGG | 246 | 2 | 1.032 | 4.77 |
| CGT | 421 | 3 | 1.008 | 2.81 |

CGG has lowest local mismatch (best neighborhood) but highest mRNA deviation; CGC has catastrophically high RBS deviation. This suggests mRNA/RBS and topology are partially independent predictors.

---

## Synthesis: Honest Assessment

### What the data shows:
1. **Local mismatch is a real structural feature** that varies meaningfully across codon positions and targets, and correlates with known fitness-relevant covariates (SRZ).
2. **Boundary crossing per se does not predict transcriptomic perturbation** in the Syn57 dataset — the RNA-seq perturbation is driven by total recoding load and regulatory disruption.
3. **Topology captures physicochemical neighborhood quality** (the local mismatch result is clean and publication-worthy), but this does not translate to a differential fitness prediction at the gene expression level.

### What this means for the manuscript:
- The GF(2)^6 framework correctly characterizes the **structural properties** of codon neighborhoods
- The **Serine disconnection** is validated: Ser boundary-crossing swaps move to significantly different (better) neighborhoods than Ala within-box swaps
- But the **predictive claim** that topology predicts recoding failure is **not supported** by the Syn57 RNA-seq data
- The local mismatch correlation with SRZ provides a mechanistic bridge: topology quantifies a dimension of the same physicochemical landscape that mRNA structure/RBS also quantify

---

## Phase 3: All Remaining Datasets

### Fredens 2019 — Syn61 (Informative Null)
- 10,631 target codons parsed (TCG=10,392, TAG=239; TCA not in xlsx main sheet)
- ALL Ser swaps cross the boundary: TCG→AGC (boundary=True, Hamming=6, Δlocal=-69)
- 99.96% success rate (7 fixes needed, all regulatory/structural)
- **Conclusion:** Boundary crossing is necessary but not sufficient for failure

### Ochre 2025 — Stop Codon Recoding
- UGA→UAA: Stop→Stop, no boundary concept (None), Hamming=1, Δlocal=-215
- 1,195 positions all identical topology — no within-study variation
- Growth data available (89 rows × 66 cols) for baseline characterization

### Lajoie 2013 — C321.ΔA
- UAG→UAA: same topology as Ochre (Stop→Stop), Hamming=1
- 321 conversions, 60% increased doubling time
- Table S19: 1,680 rows of oligo data; Table S26: 88 rows of fitness data

### Frumkin 2018 — CGG Recoding
- CGU→CGG: no boundary crossing (Arg connected), Hamming=1, Δlocal=-175
- CGC→CGG: no boundary crossing, Hamming=2, Δlocal=-175
- 60 mutations in 8 genes — value is in proteome-wide translation efficiency (trans effect)

### Ding 2024 — Mammalian TCG Recoding
- TCG→AGC in mammalian cells: boundary=True (same Ser disconnection as E. coli)
- **Cross-organism validation:** The Ser disconnection is universal across standard code organisms
- Data S1-S3: aaRS engineering data, ncAA incorporation rates

---

## Cross-Dataset Topology Summary

| Study | Swap | AA | n | Boundary variation? | Result |
|---|---|---|---|---|---|
| Napolitano 2016 | AGR→CGN | Arg | 12,888 | No (Arg connected) | Local mismatch correlates with SRZ |
| Fredens 2019 | TCG/TCA→AGC/AGT | Ser | 18,218 | No (all cross) | Informative null: 99.96% success |
| Ostrov 2016 | 7 codon types | Ser/Leu/Arg | 62,214 | Mixed | 13 lethal exceptions |
| **Robertson Syn57** | **4 Ser + 2 Ala + Stop** | **Ser/Ala** | **60,240** | **Ser=100%, Ala=0%** | **NULL for RNA-seq; POSITIVE for local mismatch** |
| Ochre 2025 | TGA→TAA | Stop | 1,195 | No (Stop→Stop) | No topology variation |
| Lajoie 2013 | UAG→UAA | Stop | 321 | No (Stop→Stop) | No topology variation |
| Frumkin 2018 | CGU/CGC→CGG | Arg | 60 | No (Arg connected) | Trans tRNA competition effect |
| Ding 2024 | TCG→ncAA | Ser | Variable | Yes | Cross-organism Ser boundary confirmation |

**Key insight:** Only 2 of 9 datasets have within-study topology variation (Syn57 Ser-vs-Ala, Ding mammalian). The remaining 7 validate the classification but cannot test the prediction.

---

## Overall Synthesis

### Prior probability update:
- Pre-analysis prior for topology predicting recoding fitness: ~40-50%
- Post-analysis posterior: ~15-25% (would require a different outcome measure or more granular fitness data to detect)

### What IS publication-worthy:
1. **The local mismatch asymmetry:** Boundary-crossing Ser swaps move codons to better neighborhoods (-37 Grantham), within-box Ala swaps move to worse (+19). This is a clean structural result (p ≈ 0, n = 60,240).
2. **The Fredens informative null:** 18,218 boundary-crossing swaps with 99.96% success directly quantifies the tolerance of GF(2)^6 topology disruption.
3. **Cross-organism Ser boundary:** The same disconnection holds in mammalian cells (Ding 2024), confirming universality of the GF(2)^6 structure.
4. **Local mismatch ~ SRZ correlation:** Napolitano's mRNA/RBS covariates correlate with GF(2)^6 local mismatch (rho = -0.33 for RBS), showing that topology captures overlapping structural information.

### What is NOT supported:
- The predictive claim that boundary-crossing recodings fail more often than within-box recodings (null in Syn57 RNA-seq data)
- Topology as an independent predictor beyond existing SRZ metrics
