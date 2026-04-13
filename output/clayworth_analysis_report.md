# Clayworth Codon-Topo Framework: Comprehensive Analysis Report

**Author**: Sergey Kornilov
**Date**: April 12, 2026
**Purpose**: Pre-meeting intelligence and scientific due diligence for co-authorship discussion with Paul Clayworth (week of April 14, 2026)

---

## Executive Summary

Paul Clayworth, a self-taught independent researcher employed as a Service Support Specialist at Fidelity Investments (San Antonio, TX), discovered a connection between the finite geometry of GF(2)^6 and the degeneracy structure of the genetic code. The core mathematical claims are **correct and independently verified** through 280 automated tests across all 25 NCBI translation tables. Several results are **genuinely novel**, most notably the characterization of serine's codon disconnection as a universal persistent homology invariant. However, the broader claims around "Clayworth Algebra," PSL(2,7) symmetry, and applications spanning quantum error correction to AI governance are either **unsubstantiated, previously rejected in the literature, or mathematically non-standard**.

This report separates verified mathematics from speculation, identifies critical prior art, assesses biological implications and limits, and provides strategic recommendations for a co-authored publication.

---

## Part I: The People

### 1.1 Paul Clayworth

**Current role**: Service Support Specialist, Fidelity Investments (June 2022-present). Prior roles: similar positions at Fidelity (2019-2021), Social Media Executive Escalations at Conduit Global (2012-2014). Holds FINRA Series 7 and Series 63 licenses.

**Education**: Biochemistry coursework listed at St. Mary's University (San Antonio) and Texas Tech University. No degree year or completion status indicated.

**Public footprint**: 54 LinkedIn connections, 75 followers. No peer-reviewed publications, preprints, or academic papers in any indexed database (arXiv, bioRxiv, Google Scholar, Semantic Scholar, ResearchGate). Technical notes TN-2026-10 and TN-2026-11 are self-published gray literature shared via LinkedIn.

**Logocentricity Inc.**: No discoverable web presence, USPTO patent filings, or company registration in US, UK, or Canadian databases. The "patent-pending" claim for the algebraic engine is unverifiable. USPTO provisional applications are not publicly indexed until 18 months after filing, but there is no confirming evidence a filing exists.

**LinkedIn posting pattern**: Claims applications of "Clayworth Algebra" spanning quantum error correction, AI hallucination reduction, interpretability, governance, algorithmic trading, consciousness, and genetic code analysis. Posts include references to a "self-correcting algebraic engine" with a "constraint kernel grounded in Fano-plane geometry and PSL(2,7)."

**Assessment**: Not an academic or working mathematician. A self-taught enthusiast with genuine exposure to finite geometry and coding theory who has found a real mathematical pattern in the genetic code. The broader framework claims are unsubstantiated.

### 1.2 Paul Charlton (mentioned collaborator)

Credited by Clayworth for the observation that "60 functional codons and icosahedral symmetry." No findable academic publications, institutional footprint, or public research presence. The claim (A5, the icosahedral rotation group, has order 60; there are 61 sense codons minus Met = 60) appears numerological without a published group-theoretic isomorphism.

### 1.3 M.E. Siddall (LinkedIn critic)

**Confirmed identity**: Mark E. Siddall, former curator at the American Museum of Natural History (1999-2020). h-index ~57, 175+ peer-reviewed publications. Specialist in phylogenetics, parasitology, and systematic biology. Terminated from AMNH in September 2020 over harassment allegations (NYT, Sept 23, 2020), which he denied. Currently describes himself as "Scientific Writer | Biostatistician | Epidemiologist."

**Assessment**: His objections (encoding arbitrariness, overfitting, no demonstrated link to selection pressure) reflect genuine expertise and are partially valid. His core point - that the framework may be "overfitting elegant math onto biology" - is the central question any reviewer will raise. Clayworth's response was unnecessarily aggressive and sarcastic, which will not play well in peer review.

---

## Part II: The Mathematics - What Is Verified

All claims below are independently verified through a codebase of 280 passing tests. The mathematical operations are correct.

### 2.1 GF(2)^6 Encoding

**Mechanism**: Each nucleotide base maps to a 2-bit pair: C=(0,0), U=(0,1), A=(1,0), G=(1,1). A 3-letter codon becomes a 6-bit tuple by concatenation. The 64 codons map bijectively to the 64 elements of GF(2)^6.

**Biochemical motivation**: The first bit distinguishes pyrimidines (C,U: first bit 0) from purines (A,G: first bit 1). The second bit distinguishes weak (2 H-bonds: A,U) from strong (3 H-bonds: C,G) base pairing. This is multiply motivated:
- Chemical: pyrimidine/purine + weak/strong (Nemzer 2017)
- Information-theoretic: error-resistant parity-check code (Mac Donaill 2003)
- Representation-theoretic: SU(2) direct sum (Sciarrino & Sorba 2012)

**Critical prior art**: Nemzer, L.R. (2017). "A binary representation of the genetic code." Biosystems 155:10-19. PMID: 28300609. Nemzer independently derives the same structure with a slightly different base ordering (U=00, C=01, A=10, G=11 vs. Clayworth's C=00, U=01, A=10, G=11). Any paper must cite Nemzer explicitly and differentiate.

**Verdict**: The encoding is sound and well-motivated but not novel. It has clear prior art.

### 2.2 Two-Fold Filtration

**Claim**: For every amino acid with exactly 2 synonymous codons, those codons differ at exactly bit position 5 (the second bit of the wobble position).

**Verification**: True across all 25 NCBI translation tables. Zero exceptions out of approximately 220 checks.

**Biological meaning**: Two-fold degenerate codons differ only in the pyrimidine/purine classification of the wobble base. This is the wobble hypothesis (Crick 1966) expressed in GF(2)^6.

**Novelty**: The biological phenomenon is textbook. The formalization as "bit-5 difference" and verification across all 25 NCBI tables is new. No prior paper performs this specific cross-validation.

### 2.3 Four-Fold Filtration

**Claim**: For every amino acid with exactly 4 codons, the codons share an identical 4-bit prefix and their 2-bit suffixes exhaust GF(2)^2 = {00, 01, 10, 11}.

**Verification**: True for all 5 four-fold amino acids in the standard code. Breaks in 11 of 25 variant codes, but only when stop codons are reassigned to amino acids.

**Novelty**: The formalization and systematic cross-validation are new.

### 2.4 Serine Disconnection (The Strongest Result)

**Claim**: Serine is the only amino acid whose codon vectors form more than one connected component at Hamming distance epsilon=1, and this holds universally across all 25 NCBI translation tables and all 24 possible base-to-bit encodings.

**Verification**: Correct. Serine has two codon families:
- UCN block: UCU, UCC, UCA, UCG (4 codons, prefix UC)
- AGY block: AGU, AGC (2 codons, prefix AG)
- Minimum inter-block Hamming distance = 4
- Reconnects only at epsilon=4

**Biological prior art**: The serine split is extensively documented:
- Schwartz et al. (2019, Sci Rep 9:17238): "Serine is the only amino acid encoded by two disjoint codon sets (TCN and AGY)"
- Rogozin et al. (2016, PNAS 113:13109): TCN-AGY switches require tandem double substitutions, driven by selection
- Inouye & Inouye (2020, PNAS 117:28572): AGY as primordial serine codons
- Bernhardt (2016, Life 6:10): tRNA^Ser class II structure as evolutionary fossil

**What is novel**: No prior paper frames the serine split as:
1. A topological invariant (disconnected connected component at Hamming epsilon<4)
2. Universal across all 25 known NCBI translation tables
3. Characterized by minimum inter-block Hamming distance (=4)
4. Detectable via persistent homology (born at epsilon=1, persists to epsilon=4)
5. Invariant across all 24 possible base-to-bit encodings

**Verdict**: The biological fact is known. The persistent homology framing is genuinely novel.

### 2.5 Fano Lines / XOR Triples

**Claim**: GGU XOR GUU XOR CAC = 0 in GF(2)^6.

**Verification**: Correct arithmetic:
```
GGU = (1,1,1,1,0,1)
GUU = (1,1,0,1,0,1)
CAC = (0,0,1,0,0,0)
XOR = (0,0,0,0,0,0)
```

**Terminology issue**: XOR-zero triples in GF(2)^6 are **linear dependencies in a vector space over GF(2)**, not "Fano lines" in the standard mathematical sense. The Fano plane is PG(2,2), living in GF(2)^3 with exactly 7 points and 7 lines. GF(2)^6 has vastly more XOR-zero triples. Calling them "Fano lines" is a non-standard terminology extension that will confuse mathematicians and irritate reviewers. The correct statement is: "GGU, GUU, and CAC form a linearly dependent triple in GF(2)^6."

**Biological relevance**: GGU = Gly (KRAS wild-type codon 12), GUU = Val (KRAS G12V mutant). The prediction that the "Fano partner" CAC (His) should be enriched in co-mutations was tested against MSK-IMPACT 2017 data via cBioPortal. Result: zero enrichment, all Fisher's exact test p=1.0 with Bonferroni correction. **The clinical prediction track is dead.** There is no known mechanism by which XOR relationships in GF(2)^6 would propagate through somatic mutation biology.

### 2.6 Holomorphic Embedding

**Claim**: phi: GF(2)^6 -> C^3, mapping each 2-bit base pair to a fourth root of unity: C->1, U->i, A->-1, G->-i.

**Verification**: The embedding preserves degeneracy structure:
- Two-fold pairs share first two complex coordinates
- Four-fold groups exhaust all four roots of unity in the third coordinate

**The "holomorphic" terminology**: Paul's reformulation via characters of GF(8)* -> C* and Spec(F_2^6) -> algebraic torus is mathematically sound. However, calling a map from a 64-element finite set "holomorphic" stretches the term beyond standard usage. More precisely: the map is a group homomorphism (character) from a finite cyclic group to the complex multiplicative group.

**Verdict**: Mathematically correct. Linguistically strained. Should be called an "algebraic embedding" or "character embedding," not "holomorphic."

### 2.7 Null Models (Statistically Sound)

**Model A** (100,000 random codon assignments): p < 0.05 for Serine uniqueness. The probability of a random genetic code having exactly one disconnected amino acid is less than 5%.

**Model B** (block-preserving shuffle): Tests whether the 16-block structure of codon space alone explains Serine's uniqueness.

**Model C** (all 24 base-to-bit encodings): Proves the filtration properties are encoding-invariant. The Serine disconnection survives all 24 encodings.

**Verdict**: Well-designed permutation tests. Standard statistical methodology, correctly implemented.

### 2.8 Novel Disconnections in Variant Codes

| Amino Acid | Genetic Code | NCBI Table | Reconnection Epsilon | Event |
|------------|-------------|------------|---------------------|-------|
| Threonine | Yeast Mitochondrial | 3 | 2 | CUN block reassigned Leu->Thr |
| Leucine | Chlorophycean Mito | 16 | 2 | UAG reassigned Stop->Leu |
| Alanine | Pachysolen Nuclear | 26 | 3 | CUG reassigned Leu->Ala |
| Serine | Candida Nuclear | 12 | 3 | CUG reassigned Leu->Ser (3 components) |

These are new observations: nobody has previously characterized these variant genetic codes in terms of persistent homology in GF(2)^6.

### 2.9 Depth Calibration (Informative Negative)

**Result**: Spearman rho = 0.0 between reconnection epsilon and evolutionary age (midpoint estimates, n=6).

**Interpretation**: Epsilon measures structural distance in GF(2)^6, not evolutionary time. The epsilon=3 events (CUG clade, ~150 Mya) are younger than the epsilon=2 events (Chlorophyceae, ~600 Mya).

**Caveat**: n=6 is extremely small. The result is suggestive but not definitive.

### 2.10 Bit-Position Bias in Reassignments

**Result**: Across 61 reassignment events, bit position 4 (first wobble bit) accumulates 13 of 35 bit-flips while position 2 has zero. Chi-square test p=0.006.

**Interpretation**: Reassignments are channeled along specific algebraic directions in GF(2)^6, with the wobble position dominating.

**Caveat**: Small sample; multiple testing concerns. But if it holds, this is a new finding.

### 2.11 Structural Robustness (Synbio Landscape)

**Result**: 83% of 1,280 possible single-codon reassignments from the standard code fully preserve the algebraic structure (filtration, Serine disconnection). Zero scores below 0.5.

**Interpretation**: The design space for synthetic code engineering is large and the algebra is remarkably hard to break. This directly supports the synthetic biology application direction.

---

## Part III: What Is Speculative or Overstated

### 3.1 "Clayworth Algebra" as Universal Framework

Paul claims a single algebraic engine based on GF(2) non-associative multiplication applies to quantum error correction, AI hallucination reduction, interpretability, governance, quantitative analysis, and the genetic code. No evidence supports this. No publications exist. No demonstrations beyond LinkedIn posts.

### 3.2 PSL(2,7) / Fano Plane as Genetic Code Symmetry

**Critical**: Antoneli & Forger (2011, Mathematical and Computer Modelling 53:1469-1488) systematically tested PSL(2,7) as a symmetry group for the genetic code and concluded it **cannot reproduce the degeneracy multiplet structure** without ad hoc constraints. This is a published, peer-reviewed rejection.

Paul may be claiming a different role for PSL(2,7) (acting on Fano-line incidence structure rather than on the 64-codon representation). If so, this distinction must be made explicit and the Antoneli/Forger paper must be cited and addressed.

### 3.3 Non-Associative Multiplication on GF(2)^6

The legitimate framework for non-associative multiplication via Fano plane is the **octonions** - an 8-dimensional algebra where the 7 imaginary units are indexed by the 7 points of the Fano plane. GF(2)^6 is 6-dimensional, which does not match. A non-associative product on GF(2)^6 using Fano structure is non-standard mathematics. The Earthform Research Lab (earthform.ai) works with the correct 8-dimensional octonion structure.

### 3.4 KRAS/Fano Clinical Predictions

Tested and failed. Zero Fano-predicted amino acids appeared at co-mutation sites across all six KRAS G12 variants in MSK-IMPACT 2017 data. Fisher's exact test with Bonferroni correction: all p=1.0. The clinical liquid biopsy application does not hold up against data.

### 3.5 Icosahedral Symmetry of 60 Codons

The connection between 60 sense codons (61 minus Met) and the icosahedral rotation group A5 (order 60) is numerological without a demonstrated group-theoretic isomorphism. No published mechanism.

### 3.6 The "Holomorphic" Framing

While Paul's reformulation via algebraic torus and GAGA principle is mathematically defensible, using "holomorphic" for a map from a finite set will draw criticism from any algebraic geometer reviewing the paper. The term should be "algebraic" or "character" embedding.

### 3.7 The p-adic / Rigid Analytic Direction

Paul mentions 2-adic integers, Tate's rigid analytic geometry, and the Ihara zeta function of the Heawood graph. These are real mathematical concepts applied speculatively. No published work connecting them to the genetic code.

---

## Part IV: Biological Implications and Limits

### 4.1 What the Framework CAN Say

1. **The genetic code has algebraic structure that is not random.** Null models confirm this statistically.

2. **Variant genetic codes respect the same algebraic constraints.** The two-fold filtration holding across all 25 NCBI tables means wobble degeneracy is universal and deeply conserved (pre-LUCA).

3. **Reconnection epsilon classifies reassignment events by structural cost.** Events that break the filtration are structurally more disruptive. This could inform synthetic biology design.

4. **The Serine disconnection is a topological invariant.** It survives every known genetic code and every possible encoding. This is as close to a "hard" mathematical result as the framework produces.

5. **The design space for synthetic code engineering is large (83%).** Most single-codon reassignments preserve the algebraic structure, meaning there are many viable options for expanded genetic codes.

### 4.2 What the Framework CANNOT Say

1. **No causal mechanism.** The framework describes structural patterns but does not explain WHY the code has this structure. Is it selection, frozen accident, or stereochemical inevitability? The math is agnostic.

2. **No predictive power for somatic mutations.** The KRAS/Fano failure demonstrates this definitively.

3. **No demonstrated selective advantage.** Even if engineered organisms violating the filtration are less stable, this could be because violating the structure requires overcoming the same biochemical constraints that created it, not because the algebra itself is functional.

4. **The encoding IS a choice.** It is well-motivated (Nemzer 2017, Mac Donaill 2003), and the Serine result is encoding-invariant, but the two-fold and four-fold filtration properties fail in 8/24 encodings. This means they partly reflect the encoding choice.

5. **Small sample sizes.** Depth calibration has n=6 data points. Bit-position bias has ~61 events. These are suggestive findings, not definitive ones.

6. **The descriptive vs. predictive boundary.** The framework may describe what evolution did without constraining what engineering can do. The geometry may be a CONSEQUENCE of the biochemistry (tRNA recognition, ribosome mechanics) rather than an independent principle. If so, the geometry is a symptom, not a cause.

### 4.3 The Synbio Counter-Argument

Even if the geometry is descriptive rather than causal, it is a **computable** descriptor. If filtration compatibility correlates with strain stability in experimental data, the geometry becomes a useful design tool regardless of causal status. This is the strongest commercial argument, and it requires experimental validation.

---

## Part V: Critical Prior Art

### 5.1 Must-Cite Papers

| Paper | Relevance | Why It Matters |
|-------|-----------|---------------|
| Nemzer (2017) Biosystems 155:10-19. PMID: 28300609 | Direct prior art for GF(2)^6 encoding | Same binary mapping, different base ordering. Must cite and differentiate. |
| Antoneli & Forger (2011) Math Comp Model 53:1469-1488 | PSL(2,7) tested and rejected for genetic code | Any PSL(2,7) claim must address this |
| Schwartz et al. (2019) Sci Rep 9:17238. PMID: 31754132 | Serine codon split biology | Modern canonical reference for the biological fact |
| Rogozin et al. (2016) PNAS 113:13109. PMC5135333 | Serine TCN-AGY switches | Selection-driven, not drift |
| Lenstra (2015) Symmetry 7:1211-1260 | CodonPolytope, Hamming metric graph | Most sophisticated prior Hamming-metric treatment |
| Mac Donaill (2003) Orig Life Evol Biosph 33:433-455. PMID: 14604185 | Error-coding justification for encoding | Independent biochemical motivation |
| Hornos & Hornos (1993) PRL 71:4401 | Sp(6) symmetry breaking | Foundational paper for algebraic genetic code models |
| Crick (1966) J Mol Biol 19:548-555 | Wobble hypothesis | The biology underlying two-fold filtration |

### 5.2 Context Papers

| Paper | Relevance |
|-------|-----------|
| Sanchez/Grau/Morgado (2005-06) | GF(4) vector space, Lie algebra approach |
| He/Petoukhov/Ricci (2004) Bull Math Biol 66:1405 | Gray code + Hamming + stochastic matrices |
| Rodriguez-Gutierrez et al. (2024) Front Appl Math Stat 10 | Binary metric + toroidal geometry |
| Bashford et al. (1998) PNAS 95:987 | Lie superalgebra sl(6/1) |
| Planat et al. (2020-21) Symmetry 12:1993 | Binary octahedral group model |
| Inouye & Inouye (2020) PNAS 117:28572 | E. coli serine codon usage evidence |
| Bernhardt (2016) Life 6:10 | tRNA^Ser evolution |

### 5.3 Synbio Context

| Paper | Relevance |
|-------|-----------|
| Beattie/Dunkelmann/Chin (2023) Nat Chem 15:903 | Quintuply orthogonal aaRS/tRNA pairs |
| Robertson et al. (2021) Science 372:1057 | Sense codon reassignment |
| Isaacs lab (2025) Nature | Ochre GRO, UAA sole stop |
| McFeely & Hartman (2023) Nat Commun 14:5008 | Breaking degeneracy with ncAAs |

---

## Part VI: Novelty Assessment Summary

### Genuinely Novel (No Prior Publication Found)

1. Two-fold = bit-5 pattern verified across all 25 NCBI tables
2. Four-fold = shared 4-bit prefix formalization
3. Serine disconnection as persistent homology invariant (universal, encoding-invariant)
4. Persistent homology applied to codon Hamming space (no prior TDA application)
5. Disconnection catalogue across variant codes (Thr, Leu, Ala discoveries)
6. Null Model C proving encoding-invariance of Serine disconnection
7. Bit-position bias in reassignment directions (chi-sq p=0.006)
8. Depth calibration negative: epsilon reflects structure not time
9. 83% structural robustness of reassignment landscape
10. Reassignment stability prediction from filtration geometry

### Has Clear Prior Art

1. GF(2)^6 encoding of codons (Nemzer 2017)
2. Purine/pyrimidine + weak/strong bit motivation (Nemzer 2017, Mac Donaill 2003)
3. Serine's two disjoint codon families (Schwartz 2019, Rogozin 2016)
4. Hamming distances between codons (He 2004, Lenstra 2015)
5. Algebraic structure in the genetic code (Hornos & Hornos 1993 and many others)

### Previously Rejected or Problematic

1. PSL(2,7) as genetic code symmetry (Antoneli & Forger 2011)
2. "Fano line" terminology for XOR triples in GF(2)^6 (non-standard)
3. Non-associative multiplication on GF(2)^6 (dimensionally wrong for octonions)
4. KRAS/Fano clinical predictions (failed against data)

---

## Part VII: Strategic Recommendations

### 7.1 Paper Scope

**Include (Paper 1, publishable now)**:
- GF(2)^6 filtration structure (cite Nemzer, differentiate)
- Serine disconnection as persistent homology invariant
- Universal verification across all 25 NCBI tables
- Disconnection catalogue (Thr, Leu, Ala, 3-component Ser)
- Null models A, B, C (encoding-invariance proof)
- Bit-position bias in reassignments
- Depth calibration: epsilon reflects structure not time
- Structural robustness landscape (83%)

**Exclude (Paper 1)**:
- KRAS/Fano clinical prediction (failed)
- Holomorphic embedding (save for math audience)
- PSL(2,7) / Heawood graph / Ihara zeta function
- "Clayworth Algebra" as universal framework
- Quantum error correction, AI, governance applications
- Any reference to Logocentricity Inc. or patent-pending claims in the abstract

**Defer (Paper 2, needs experimental validation)**:
- Synbio feasibility scoring tool
- Meta-analysis of published reassignment experiments vs. filtration predictions
- Experimental collaboration for stability testing

### 7.2 Authorship

Co-first authorship is appropriate and well-justified. Sergey's contributions beyond verification: bit-position bias finding, disconnection catalogue extension, depth calibration insight (epsilon != time), 83% robustness quantification, KRAS negative result (saves Paul from publishing a false positive), synbio application direction.

### 7.3 Target Journals

| Journal | Pros | Cons |
|---------|------|------|
| PLoS Computational Biology | Good fit, open access, accepts independent researchers | Competitive |
| Journal of Theoretical Biology | Mathematical biology, no wet-lab data required | Lower impact |
| BMC Bioinformatics | Computational focus, accessible | Lower prestige |
| Nucleic Acids Research | High impact | Too competitive for first-time unaffiliated authors |

### 7.4 Meeting Preparation

1. Bring the Nemzer (2017) and Antoneli/Forger (2011) citations. Paul may not know these exist.
2. Propose a tight paper scope. The paper is stronger without "Clayworth Algebra" branding.
3. Correct "Fano line" terminology before submission. Use "linearly dependent triple in GF(2)^6."
4. Discuss authorship contribution statement explicitly.
5. Clarify the relationship between the paper and any IP/VC claims.
6. Agree on a target journal.

### 7.5 Risk Factors

1. **Scope creep**: Paul will want to include the broader framework. The paper is stronger without it.
2. **Patent/IP entanglement**: Ensure the paper's scientific claims stand independently of commercial interests.
3. **"Fano line" terminology**: Will cause reviewer pushback if not corrected.
4. **Siddall's objections**: The encoding-dependence critique is addressed by Null Model C, but the "descriptive vs. predictive" question will be raised by any competent reviewer. Address it head-on.
5. **Paul's LinkedIn responses**: His aggressive tone with critics (Siddall) suggests potential difficulty with peer review feedback.

---

## Appendix A: Codebase Verification Summary

- **Total tests**: 280 (all passing)
- **Test coverage**: 96%
- **Core modules verified**: encoding.py, filtration.py, homology.py, embedding.py, fano.py, null_models.py, genetic_codes.py
- **Analysis modules verified**: reassignment_db.py, depth_calibration.py, cosmic_query.py, synbio_feasibility.py
- **Regression tests**: 105 tests reproducing exact values from PRD Appendix 8
- **Technology**: Python 3.11, NumPy, SciPy, pytest + hypothesis for property-based testing

## Appendix B: Key Numerical Results

| Metric | Value | Significance |
|--------|-------|-------------|
| Two-fold filtration pass rate | 100% (all 25 tables) | Universal invariant |
| Four-fold filtration pass rate | 100% (standard code) | 19/25 across all tables |
| Serine disconnection | 25/25 tables, 24/24 encodings | Universal topological invariant |
| Min inter-block Hamming (Ser) | 4 | Reconnects only at epsilon=4 |
| Null Model A p-value (Serine unique) | <0.05 (100k permutations) | Statistically significant |
| KRAS Fano enrichment | p=1.0 (all variants) | No enrichment detected |
| Depth calibration Spearman rho | 0.0 (n=6) | Epsilon != evolutionary time |
| Bit-position bias chi-sq p | 0.006 | Significant non-uniform distribution |
| Structural robustness | 83% of 1,280 reassignments | Large viable design space |

---

*This report was prepared for internal use in advance of co-authorship discussions. All mathematical claims have been independently verified through computational reproduction.*
