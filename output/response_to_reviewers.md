# Response to Reviewers — BIOSYS-D-26-00689

**Manuscript:** Robust error-minimization in the genetic code across physicochemical metrics and variant codes: a graph-theoretic analysis in GF(2)⁶
**Authors:** Paul Clayworth, Sergey Kornilov

---

Dear Dr. Prinz and Reviewers,

We thank the Editor and the reviewers for their detailed and constructive comments. The revision below addresses every point raised, adds the analyses we believe are needed to settle Reviewer 2's principal methodological concerns, and integrates the recent BioSystems literature flagged by Reviewer 1. Section numbers refer to the revised manuscript (`manuscript.typ`) and supplement (`supplement.typ`). All new results are reproducible from the public repository (`codon-topo all --output-dir=./output --seed=135325`) at the commit tagged in the manuscript reproducibility stamp.

A note on Reviewer 3: the decision letter indicated that Reviewer 3's comments were attached; however, no Reviewer 3 report was available in the decision-letter text or in the files accessible to us at the time of the first revision. The report has since been received; our response to it is included below (§ Reviewer 3).

---

## Reviewer 1

> *"The manuscript would also benefit from discussion of recent work by Ekaterina Yurova Axelsson and Andrei Khrennikov on 2-adic and dynamical-system representations of the genetic code published in BioSystems …"*

We thank the reviewer for pointing us to the BioSystems 2-adic literature. We have added a Discussion paragraph (§4.5, "Correspondence with 2-adic codon algebra") comparing our finite graph-theoretic framework with the 2-adic attractor models of Yurova Axelsson and Khrennikov, as well as related p-adic work by Khrennikov and Kozyrev and by Dragovich and Mišić. We now clarify that both approaches treat the genetic code as a structured mathematical object, but they differ in geometry and inferential aim: the 2-adic work models code organisation in an ultrametric dynamical state space in which amino acids appear as attractors, whereas our analysis uses finite Hamming graphs, observed reassignment events, and statistical null models. The two pictures are related but not identical — 2-adic balls and Hamming neighbourhoods are different mathematical objects — and we do not identify them.

To make the relationship concrete, we have added a supplemental Walsh–Hadamard / 2-adic spectral probe as a finite bridge between these perspectives (Supplement §S16; code at `src/codon_topo/analysis/walsh_2adic.py`; tests at `tests/test_walsh_2adic.py`, 17/17 passing). The standard code's total Walsh 2-adic depth across synonymous-block indicator functions is 544, against a block-size matched null mean of 689.3 ± 8.2 (n = 2,000; z = −17.74; 0/2000 nulls beat the standard code). The depth value 544 is mathematically invariant across all 24 base-to-bit encodings. Under a stricter wobble-box-preserving label-permutation null (fixing the 16 first-two-base boxes and their internal partition shapes and randomising only the AA labels), the depth is mathematically constant at 544 — an algebraic invariant of the wobble-box × AA-slot-multiset structure. The result is explicitly interpreted as an algebraic invariant of wobble-box structure, not as evidence for an additional beyond-wobble optimisation signal.

The bibliography now includes Yurova Axelsson & Khrennikov 2024a (*BioSystems* 240:105230), 2024b (*BioSystems* 246:105353), Khrennikov & Kozyrev (2007, *Physica A* 381:265), and Dragovich & Mišić (2019, *BioSystems* 185:104017).

We have added a second Walsh invariant of the same wobble structure (Supplement §S16.3): the wobble-free label-spectrum fraction *S* = 0.7514309076, also encoding-invariant across all 24 base-to-bit bijections. The standard code's amino-acid labeling concentrates Walsh-Hadamard energy on the wobble-free spectral layer to a degree unattainable under a degeneracy-multiset-preserving label permutation (*z* = +30.3, *p* < 10⁻⁴ at *n* = 10,000), but is mathematically invariant under the stricter wobble-box-preserving label permutation — so *S* is a spectral signature of biological wobble-box alignment, not evidence of additional optimisation. Together with spectral_depth = 544, the two encoding-invariant Walsh signatures provide a complete algebraic characterisation of the wobble structure in the discrete shadow of Khrennikov's ultrametric framework. We explicitly disclaim a "beyond wobble" interpretation in §S16.3 and the Discussion.

---

## Reviewer 2

We thank Reviewer 2 for the careful methodological reading. The three concerns are addressed below; in summary, the criticism that the original four metrics overlapped is correct, and we have rewritten the relevant framing throughout. The suggestion to test ProtSub has been implemented in full, and the Di Giulio (2001) caution is now made explicit in the choice of primary measures.

### R2.1 — Overlap among the four physicochemical measures

> *"The list of models taken for the analysis by Clayworth and Kornilov is thus inconsistent; the experimentally obtained polarity variable by Woese overlaps with its statistical correction by Grantham, and 3 variable model by Grantham overlaps with the 2 variable version by Miyata et al. (1979)."*

We agree with the reviewer that the four primary amino-acid measures are not statistically independent. Grantham and Miyata share polarity and volume components, and the full panel reflects overlapping hydrophobicity/polarity structure. We have removed language implying independent replication throughout the manuscript:

- The Introduction phrase "four metrics with disjoint conceptual bases" has been replaced with the explicit acknowledgement that the measures are "established but partially overlapping physicochemical parameterizations" (§1, end of paragraph beginning "Previous work has not exploited this decomposition…").
- The Methods (§2.2) and the Results sentence header (§3.1) now refer to "four established, code-independent physicochemical distance measures with partially overlapping content," not to independent metrics.
- The Supplement claim-hierarchy entry for *hypercube_coloring_optimality* (§S1) has been updated in the same way.

To quantify the degree of overlap, we added metric-correlation analyses at two levels (§3.1; Supplement §S15):

- *Pairwise Spearman across the 190 unordered amino-acid pairs:* range 0.23–0.73; median 0.49. Highest pairs (Grantham–Miyata ρ = 0.71, Miyata–PR ρ = 0.73) confirm the reviewer's point about polarity-component overlap.
- *Pairwise Spearman across 2,000 random codes drawn from the same Freeland–Hurst null:* range 0.71–0.89; median 0.82. Code-level rankings are more correlated than AA-pair distances, as expected when distances are summed over 192 hypercube edges.
- *Partial Spearman correlations of code-level F-scores, controlling for the other measures:* five of ten pairs become statistically independent (|partial ρ| < 0.1) after partialling. Substantial shared variance coexists with non-identical residual structure across the panel.

We now describe the cross-metric concordance as a *convergent sensitivity envelope* rather than independent replication. The claim is that if the standard code lies in the low-cost tail under several established but non-identical physicochemical parameterizations, the result is not confined to a single chosen scale. With five tests at α = 0.05, the Bonferroni-corrected per-test threshold is 0.01; all five p-values pass (Grantham p = 0.006; Miyata p < 0.001; Woese polar requirement p = 0.003; Kyte–Doolittle p = 0.001; ProtSub p = 0.0004).

### R2.2 — ProtSub matrix

> *"It would be better if the authors tested the method on, e.g. ProtSub matrix - derived 'to bring sequence alignments into agreement with structure matches' (Jia K, Jernigan RL. Proteins. 2021; 89:671-682. doi.org/10.1002/prot.26050.)."*

Following the reviewer's recommendation, we added ProtSub as a fifth measure (Methods §2.2; Results §3.1, including the table; Supplement §S15).

ProtSub is provided by Jia & Jernigan as a log-odds substitution matrix in EMBOSS format. We converted it to a positive distance using the standard diagonal-anchored relation *d(i,j) = ½(s(i,i) + s(j,j)) − s(i,j)*, which guarantees *d(i,i) = 0* and is strictly positive off-diagonal for ProtSub (verified). Sanity values are physicochemically sensible: d(I,L) = 2.0, d(K,R) = 4.0, d(F,Y) = 4.0, d(D,E) = 3.5, d(C,W) = 14.5.

Under the same Freeland–Hurst block-preserving null (n = 10,000), the standard code's F-score under ProtSub is 1089.5 versus a null mean of 1180.2 ± 29.7 (z = 3.05; quantile 0.04 %; p = 0.0004). The result satisfies Bonferroni at α = 0.05 across the five-measure panel.

Because ProtSub is alignment-derived, it is subject to the same general class of circularity concerns raised by Di Giulio (2001) — see R2.3 below. We therefore do not use ProtSub as primary evidence for code-origin optimality. Instead, we report it as a structure-aware sensitivity analysis asking whether the standard-code result persists under the matrix suggested by the reviewer. It does.

### R2.3 — Di Giulio (2001)

> *"Considering all, I will add that the general methodological problem was addressed by Massimo Di Giulio (2001) in Article 'The origin of the genetic code cannot be studied using measurements based on the PAM matrix because this matrix reflects the code itself, making any such analyses tautologous'…"*

We agree with Di Giulio's caution as applied to PAM-derived and other alignment-derived measurements; this was the explicit reason for our choice of metric panel, and we have made the rationale explicit in the revised Methods (§2.2). The four primary measures — Woese polar requirement (chromatographic partition coefficient), Grantham 1974 (composition, polarity, volume), Miyata 1979 (polarity and volume only), and Kyte–Doolittle 1982 (hydropathy from transfer free energies) — are all *code-independent* physicochemical scales, derived from properties of amino acids in isolation rather than from observed substitution patterns in coded proteins. We also cite Koonin & Novozhilov (2009; `koonin2009` in the bibliography), who explicitly endorse code-independent metrics for exactly this reason.

The Buschmann et al. (2026, *Bioinformatics* 42(5):btag188, "Much ado about nothing: modeling amino acid replacement with predicted protein structures") study, which we now cite alongside Di Giulio, benchmarked BLOSUM62, an AlphaFold-derived matrix (AFSM), and 16 other substitution matrices and reports that all matrix families perform similarly — their interpretation being that substitution matrices implicitly encode physicochemical reality. We read this as bounding rather than eliminating the Di Giulio circularity concern for alignment-derived matrices: ProtSub still captures genuine physicochemical signal, but its evidentiary status for code-origin questions is distinct from that of measures derived independently of the code.

### R2.4 — Cluster structure of polarity scales (Kosiol et al. 2004)

> *"observed polarity scales consist of basic clusters, specified by the subsets of polar and nonpolar amino acids — see e.g. Fig. 5 & Fig. 6 in Kosiol et al. (J Theor Biol. 2004, 228(1): 97-106 …) for PAM and WAG matrices."*

The Kosiol et al. analysis of PAM and WAG groupings shows a multi-level hierarchical structure — k = 2 (polar/non-polar), k = 3 (aliphatic/aromatic split), k = 4 (charged separated from small-polar), and finer — rather than a single binary partition. The polarity-based measures in our panel (Woese PR, Miyata) reflect this multi-level structure rather than binary polar/non-polar separation, which is why they give different test statistics in our setting (Miyata p < 0.001 vs Grantham p = 0.006, a three-fold range in significance). The Kosiol et al. reference is now included in the bibliography and acknowledged in the §3.1 discussion of metric overlap.

---

## Slavov et al. 2026 — convergent event-level context

We have also added a Discussion paragraph (§4.6) noting the very recent publication of Tsour, Machné, Leduc, Widmer, Koo, Guez, Karczewski & Slavov (2026, *Nature*, doi:10.1038/s41586-026-10678-2, "Alternate RNA decoding results in stable and abundant proteins in mammals"). This work directly observes 8,746 unique amino-acid substitutions across 1,767 genes in human and mouse tissues arising from sense-codon recoding, with codon–anticodon mismatches identified as a primary mechanism.

For each high-confidence event in their Supplementary Data 8 (5,873 events; 5,611 unambiguous-pair events; positional probability > 0.9) we computed the minimum nucleotide-Hamming distance between any source codon and any target codon under the standard genetic code (Supplement §S17). At the *event* level, 65.0 % of observations involve source/target codon pairs at min-NT-1, versus a 39.5 % baseline among all 380 ordered amino-acid pairs (binomial test, one-sided p < 10⁻¹⁰⁰). At the level of unique substitution *types*, the enrichment is modest (41.6 % observed vs 39.5 % baseline; Fisher exact p = 0.26; Mann–Whitney U on NT distances, observed vs unobserved types, p = 0.18); almost every amino-acid substitution type can be observed in their dataset, but the *frequency* of each type is biased toward single-nucleotide-distance codon pairs. We describe this in the Discussion as convergent event-level context for the biological relevance of local codon neighbourhoods, with the unique-pair analyses reported as more conservative.

---

## Reviewer 3

> *"It is quite natural to use the rather simple assignment of 2-bit words to nucleotide bases. However, in some sense this does not mirror the complexity of the genetic code since it implies that all features etc. of codons are implied by properties of their bases. Is this true or a reasonable assumption? For instance, in [1] 'The non-power model of the genetic code: a paradigm for interpreting genomic information' (by Gonzalez, Giannerini, Rosa) a non-power model of the genetic code was developed that is not induced by one of the 24 bijections investigated in the current paper. Could this give even deeper results under the author's methods? I think this should be discussed and also considered as an alternative."*

We thank the reviewer for the recommendation to accept and for the point about non-power codon representations. The observation is well-taken: the 24 base-to-bit bijections we sweep are all *power* representations, in the specific sense that each codon's 6-bit vector is the concatenation of three per-base 2-bit encodings, so codon-level properties are determined by the choice of the base-level bijection. The Gonzalez et al. (2016) non-power model — using a redundant signed-binary numeration with Fibonacci-like positional weights [1,1,2,4,7,8] so that codon-to-integer maps become intentionally many-to-one, mirroring code degeneracy — sits outside this family and is a strictly larger design space.

We have added a paragraph to Discussion §4.5 ("Broader connections") that (i) cites Gonzalez, Giannerini and Rosa (2016) and describes the non-power construction, (ii) makes the point that our primary topology-avoidance result is defined on the encoding-independent Hamming graph *H*(3,4) — which depends only on which nucleotides differ between two codons, not on any codon-level bijection at all, whether power or non-power — so the paper's principal topology claim is not exposed to this concern, and (iii) explicitly flags a non-power extension of the *Q*₆ subgraph analysis (which is already encoding-dependent within the power family, §3.4 and Supplement §S4) as a natural continuation. We have not attempted the full non-power sweep in this revision because it requires reformulating the underlying mutation graph, which changes when codon-to-integer maps become many-to-one, and the paper's principal test is already downstream of any such extension via the *H*(3,4) result. The bibliography now includes Gonzalez, Giannerini and Rosa (2016, *Philos Trans R Soc A* 374(2063):20150062).

---

## Summary of changes to the package

**New analyses:**
1. Pairwise Spearman correlations of the four code-independent measures plus ProtSub across the 190 unordered AA pairs (Supplement §S15).
2. Code-level F-score correlations and partial F-score correlations among the five measures across 2,000 random codes (Supplement §S15).
3. ProtSub as a fifth measure, with full Freeland–Hurst null Monte Carlo (n = 10,000).
4. Walsh–Hadamard / 2-adic spectral probe with two null tests, block-size matched and wobble-box-preserving label permutation (Supplement §S16).
5. Wobble-free label-spectrum invariant S = 0.7514 (Supplement §S16.4).
6. Slavov SAAP cross-analysis (Supplement §S17).

**New code (all in repository):**
- `src/codon_topo/data/protsub.json` (data file derived from Jia & Jernigan EMBOSS matrix)
- `src/codon_topo/analysis/walsh_2adic.py` (~600 LOC; wobble caveat preserved verbatim)
- `tests/test_walsh_2adic.py` (17 tests, all passing)
- Extensions to `src/codon_topo/analysis/coloring_optimality.py` (ProtSub integration; `CODE_INDEPENDENT_METRICS` constant with `multi_metric_sensitivity` now defaulting to the four code-independent measures only).

**New references added to bibliography:**
- Yurova Axelsson & Khrennikov 2024a (*BioSystems* 240:105230)
- Yurova Axelsson & Khrennikov 2024b (*BioSystems* 246:105353)
- Khrennikov & Kozyrev 2007 (*Physica A* 381:265)
- Dragovich & Mišić 2019 (*BioSystems* 185:104017)
- Di Giulio 2001 (*J Theor Biol* 208:141)
- Jia & Jernigan 2021 (*Proteins* 89:671)
- Buschmann et al. 2026 (*Bioinformatics* 42:btag188)
- Kosiol et al. 2004 (*J Theor Biol* 228:97)
- Tsour et al. 2026 (*Nature*, doi:10.1038/s41586-026-10678-2)

**Manuscript edits:**
- Abstract framing reviewed for consistency with the convergent-envelope language.
- Introduction (§1): "disjoint conceptual bases" language removed; convergent-envelope framing added with forward reference to §3.1.
- Methods §2.2: explicit Di Giulio (2001) rationale for the choice of code-independent measures; ProtSub described as a structure-aware sensitivity analysis.
- Results §3.1: new paragraph on metric correlations; ProtSub row added to Table 1; Bonferroni-corrected reporting; "five measures" framing.
- Discussion §4.5: new subsection on correspondence with 2-adic codon algebra and the Walsh–Hadamard bridge.
- Discussion §4.6: new subsection on Tsour/Slavov 2026 alternate-translation evidence.
- Supplement §S1: claim-hierarchy entry updated to reflect partially overlapping content.
- Supplement §S15: ProtSub data conversion + metric correlation tables.
- Supplement §S16: Walsh–Hadamard / 2-adic spectral probe.
- Supplement §S17: Slavov SAAP cross-analysis.

We hope the revision addresses the reviewers' concerns and look forward to any further feedback.

Sincerely,

Paul Clayworth and Sergey Kornilov
