// ============================================================
// Supplementary Material
// Robust error-minimization in the genetic code
// Clayworth & Kornilov
// ============================================================

#set page(paper: "a4", margin: (x: 2.5cm, y: 2.5cm), numbering: "S1")
#set text(font: "Libertinus Serif", size: 10.5pt, lang: "en")
#set par(justify: true, leading: 0.65em)
#set heading(numbering: "S1.1")
#set figure(gap: 0.8em)
#show table.cell.where(y: 0): set text(weight: "bold")
#show table: set text(number-width: "tabular", size: 9.5pt)

#align(center)[
  #text(14pt, weight: "bold")[Supplementary Material]
  #v(0.5em)
  #text(11pt)[
    Robust error-minimization in the genetic code across \
    physicochemical metrics and variant codes: \
    a graph-theoretic analysis in GF(2)#super[6]
  ]
  #v(0.3em)
  #text(10pt)[Paul Clayworth and Sergey Kornilov]
  #v(1.5em)
]

// ============================================================
= Claim hierarchy with full justifications <sec:s-claims>

The following table lists all 15 evaluated claims with their complete justification text, as maintained in the `claim_hierarchy.py` source of truth.

#figure(
  table(
    columns: (auto, auto, 1fr),
    align: (left, center, left),
    inset: 5pt,
    stroke: (x, y) => if y == 0 { (bottom: 0.7pt) } else { none },
    table.header([*Claim ID*], [*Status*], [*Justification*]),
    [hypercube_coloring_optimality], [Supported], [Cross-metric sensitivity: Grantham $p = 0.003$, Miyata $p = 0.001$, polar requirement $p = 0.004$, Kyte--Doolittle $p = 0.001$. Stop penalty sensitivity (0/150/215/300): immaterial.],
    [per_table_optimality_preservation], [Supported], [24/25 tables significant at $p < 0.05$ after BH-FDR. Mean quantile 1.4%. Only Table 3 (yeast mito) at 7.4%.],
    [optimality_rho_robustness], [Supported], [$rho$-sweep: $p < 0.003$ at all values. $z$ increases monotonically from 2.48 ($rho = 0$) to 3.46 ($rho = 1$).],
    [topology_avoidance_depletion], [Supported], [Permutation $p lt.eq 10^(-4)$; hypergeometric $p = 4.8 times 10^(-8)$. Clade exclusion: all $p < 10^(-5)$.],
    [trna_enrichment_reassigned_aa], [Suggestive], [MIS worst-case $p = 0.045$. 18 verified assemblies, 24 pairings, 5 variant codes.],
    [bit_position_bias_weighted], [Exploratory], [Uniform $p = 0.006$ inflated by non-independence. De-duplicated $p = 0.075$.],
    [mechanism_boundary_conditions], [Exploratory], [Three-tier pattern: duplication / stem shortening / modification. Descriptive.],
    [atchley_f3_serine_convergence], [Exploratory], [Serine $F_3 = -4.760$, 2.24 SD below mean. Complementary, not independent.],
    [variant_code_disconnection_catalogue], [Exploratory], [4 variant-code cases: Thr (Table 3), Leu (Table 16), Ala (Table 26), Ser (Table 12).],
    [kras_fano_clinical_prediction], [Falsified], [$p = 1.0$ across all 6 G12 variants. $n = 1{,}670$ MSK-IMPACT mutations.],
    [serine_min_distance_4_invariant], [Rejected], [16/24 encodings give distance 2. Only 8/24 give distance 4.],
    [psl_2_7_symmetry], [Rejected], [No 64-dim irrep. #cite(<antoneli2011>, form: "prose").],
    [holomorphic_embedding], [Rejected], [Domain is finite discrete. Character identity fails: $i^2 = -1 eq.not 1$.],
    [two_fold_bit_5_filtration], [Tautological], [Forced by encoding choice. Holds in 16/24 encodings.],
    [four_fold_prefix_filtration], [Tautological], [Trivial under any bijection from 4 bases to $"GF"(2)^2$.],
  ),
  caption: [Complete claim hierarchy with justifications.],
) <tbl:s-claims>


// ============================================================
= Encoding sensitivity analysis <sec:s-encoding>

The coloring optimality result was tested under all 24 base-to-bit bijections from ${C,U,A,G}$ to $"GF"(2)^2$. All 24 encodings yield significant optimality ($p < 0.05$ under the block-preserving null with $n = 1{,}000$). The mean quantile across encodings is 1.8%, confirming that the result is not an artifact of the default encoding choice.

Properties that are encoding-invariant:
- Serine disconnection at $epsilon = 1$ (holds under all 24 encodings)
- Coloring optimality significance (all 24 significant)
- Four-fold prefix filtration (tautological under any bijection)

Properties that are encoding-dependent:
- Serine inter-family minimum Hamming distance (4 in 8/24, 2 in 16/24)
- Two-fold bit-5 filtration (holds in 16/24)
- Specific score values and rank orderings

Full per-encoding results are available in the `codon-topo` repository via `codon-topo coloring --all-encodings`.


// ============================================================
= Complete tRNA gene count data <sec:s-trna>

== tRNAscan-SE verified organisms

All tRNA gene counts were obtained by running tRNAscan-SE 2.0.12 @chan2019 with Infernal 1.1.4 on NCBI genome assemblies. Eukaryotic organisms were scanned in `-E` mode; _Mycoplasma_ species in `-B` (bacterial) mode.

#figure(
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto),
    align: (left, center, left, right, right, right, left),
    inset: 5pt,
    stroke: (x, y) => if y == 0 { (bottom: 0.7pt) } else { none },
    table.header([*Organism*], [*Tbl*], [*Assembly*], [*Total*], [*Std20*], [*Supp*], [*Reassigned AA*]),
    [_T. thermophila_], [6], [GCF_000189635.1], [718], [672], [37], [Gln: 54 (15+39)],
    [_P. tetraurelia_], [6], [GCF_000165425.1], [216], [202], [11], [Gln: 18 (7+11)],
    [_O. trifallax_], [6], [GCA_000295675.1], [94], [83], [6], [Gln: 8 (2+6)],
    [_P. persalinus_], [6], [GCA_001447515.1], [262], [228], [15], [Gln: 20 (5+15)],
    [_H. grandinella_], [6], [GCA_006369765.1], [130], [121], [3], [Gln: 9 (6+3)],
    [_E. aediculatus_], [10], [GCA_030463445.1], [80], [76], [2], [Cys: 4 (3+1 TCA)],
    [_E. amieti_], [10], [GCA_048569255.1], [120], [103], [6], [Cys: 8 (4+4 TCA)],
    [_E. focardii_], [10], [GCA_001880345.2], [62], [56], [3], [Cys: 3 (1+2 TCA)],
    [_E. parawoodruffi_], [10], [GCA_021440025.1], [149], [128], [5], [Cys: 9 (5+4 TCA)],
    [_E. weissei_], [10], [GCA_021440005.1], [495], [390], [13], [Cys: 20 (17+3 TCA)],
    [_E. woodruffi_], [10], [GCA_027382605.1], [83], [74], [4], [Cys: 4 (2+2 TCA)],
    [_B. stoltei_], [15], [GCA_965603825.1], [169], [165], [0], [Trp: 6],
    [_B. nonstop_ P57], [31], [GCA_028554745.1], [68], [65], [2], [Trp: 2#super[\*]],
    [_M. genitalium_], [4], [GCA_000027325.1], [36], [34], [0], [Trp: 1#super[#sym.dagger]],
    [_M. pneumoniae_], [4], [GCF_910574535.1], [37], [35], [0], [Trp: 1#super[#sym.dagger]],
    [_S. coeruleus_], [1], [GCA_001970955.1], [272], [265], [1], [---],
    [_I. multifiliis_], [1], [GCF_000220395.1], [150], [141], [8], [---],
    [_F. salina_], [1], [GCA_022984795.1], [89], [85], [1], [---],
  ),
  caption: [
    Complete tRNAscan-SE 2.0.12 results. "Std20" = standard 20-AA tRNAs. "Supp" = suppressor tRNAs (CTA/TTA/TCA anticodons). Reassigned AA counts include both standard and suppressor tRNAs counted toward the reassigned amino acid. #super[\*]Anticodon stem shortening. #super[#sym.dagger]Post-transcriptional modification.
  ],
) <tbl:s-trna>


== MIS (maximal independent set) analysis

To address non-independence from shared control organisms, we constructed a conflict graph where edges connect pairings sharing an organism. The Bron--Kerbosch algorithm with pivoting enumerates all maximal independent sets (MIS). Each MIS represents a set of pairings where no two share an organism and no additional pairing can be added without creating a conflict.

With 24 total pairings, 2 MIS of size 6 were identified. Both are significant at $p < 0.05$ under Stouffer's $Z$ combination of per-pairing Fisher exact $p$-values:

- *Best-case MIS* ($p = 0.044$): S. cerevisiae mito/Thr, S. obliquus mito/Leu, P. tannophilus/Ala, P. tetraurelia/Gln, E. aediculatus/Cys, E. weissei/Cys
- *Worst-case MIS* ($p = 0.045$): S. cerevisiae mito/Thr, S. obliquus mito/Leu, C. albicans/Ser, P. tetraurelia/Gln, E. aediculatus/Cys, E. weissei/Cys


// ============================================================
= Complete reassignment database <sec:s-reassignment>

The reassignment database comprises 61 codon reassignment events across the 25 NCBI translation tables, relative to the standard code (Table 1). Each event records the codon, source amino acid, target amino acid, and Hamming distance to the nearest codon already encoding the target amino acid in the standard code. De-duplication to unique (codon, target amino acid) pairs yields 27 events used in the topology avoidance analysis.

The complete database is provided as `output/tables/T10_reassignment_db.csv`.


// ============================================================
= Topology avoidance: clade-exclusion sensitivity <sec:s-clade>

To test whether the topology avoidance result is driven by any single phylogenetic clade, we iteratively excluded:
- All ciliate reassignments
- All metazoan mitochondrial reassignments
- All CUG-clade yeast reassignments
- All chlorophycean reassignments

In every exclusion, the depletion remains highly significant ($p < 10^(-5)$), confirming that the 3.3-fold depletion of topology-breaking changes is a pan-taxonomic pattern, not an artifact of any single lineage.


// ============================================================
= KRAS--Fano clinical prediction: detailed results <sec:s-kras>

The conjecture that XOR ("Fano") relationships in $"GF"(2)^6$ predict enrichment of specific amino acids at KRAS G12 co-mutation sites was tested against 1,670 KRAS mutations from the MSK-IMPACT dataset @zehir2017.

For each of the 6 KRAS G12 variant types (G12D, G12V, G12C, G12A, G12R, G12S), we identified the XOR-predicted amino acid partners and tested for co-mutation enrichment via Fisher's exact test with Bonferroni correction. All 6 tests yielded $p = 1.0$, with odds ratios near 1.0.

This result cleanly separates code-level error-minimization (which operates on the amino acid assignment structure) from mutation-level algebraic predictions (which would require DNA polymerase errors to respect the binary encoding --- a biologically implausible mechanism).


// ============================================================
= Event-level conditional logit model: full details <sec:s-condlogit>

== Model specification and candidate universe

The conditional logit model treats each observed reassignment as a discrete choice from the set of all single-codon reassignments available at the current code state. For a code with 64 codons and $approx 21$ possible amino acid/stop labels, the candidate set $cal(N)(C)$ contains $approx 1{,}280$ moves (64 codons $times$ 20 alternative labels, excluding identity assignments). All 1,280 candidates are evaluated at each step regardless of biological plausibility; the model's purpose is to test whether the observed moves are statistically distinguishable from uniform sampling given the three feature classes. As a conditional logit, the model implies an independence of irrelevant alternatives (IIA) structure within each choice set. Here the goal is not mechanistic completeness but a controlled test of whether observed moves are statistically enriched for low $Delta_"phys"$ and low $Delta_"topo"$ relative to a uniform baseline over the same candidate set; a mixed logit relaxing IIA would be appropriate if the goal were prediction rather than hypothesis testing.

== Feature definitions

*Local physicochemical mismatch change ($Delta_"phys"$)*: For a reassignment of codon $c$ from amino acid $a$ to $a'$, this is the change in the sum of Grantham (1974) distances across all Hamming-1 edges incident to $c$:

$ Delta_"phys" = sum_({c,c'}: d(c,c')=1) [Delta(a', "code"(c')) - Delta(a, "code"(c'))] $

This captures the local impact on error-minimization at the reassigned position. When a neighboring codon is assigned Stop, we set $Delta(a, "Stop")$ equal to the maximum Grantham distance (215), consistent with the stop-penalty convention used in the coloring objective (@sec:objective).

*Topology disruption ($Delta_"topo"$)*: The total increase in connected components (at $epsilon = 1$) summed across all amino acid codon graphs. Stop codons are excluded from the topology sum; $Delta_"topo"$ counts connected-component changes only for amino acid codon families. A move that splits one amino acid's codon family into two components contributes $+1$; a move that fragments two families contributes $+2$; topology-preserving moves contribute $0$.

*tRNA complexity proxy ($Delta_"tRNA"$)*: The Hamming distance from the reassigned codon to the nearest codon already encoding the target amino acid in the current code. This serves as a heuristic for the tRNA repertoire change required to service the reassigned codon, with larger distances implying more novel tRNA machinery needed.

== Fitted coefficients

#figure(
  table(
    columns: (auto, auto, auto, auto),
    align: (left, left, center, center),
    inset: 5pt,
    stroke: (x, y) => if y == 0 { (bottom: 0.7pt) } else { none },
    table.header([Model], [Feature], [$hat(beta)$ (raw)], [$hat(beta)$ (normalized)]),
    [M1], [$Delta_"phys"$], [$-0.0045$], [$-0.65$],
    [M2], [$Delta_"topo"$], [$-3.504$], [$-3.83$],
    [M3], [$Delta_"phys"$], [$-0.0039$], [$-0.56$],
    [M3], [$Delta_"topo"$], [$-3.263$], [$-3.56$],
    [M4], [$Delta_"phys"$], [$-0.0039$], [$-0.56$],
    [M4], [$Delta_"topo"$], [$-2.911$], [$-3.17$],
    [M4], [$Delta_"tRNA"$], [$-0.219$], [$-0.30$],
  ),
  caption: [
    Conditional logit coefficient estimates. Raw coefficients are on the original feature scale; normalized coefficients are on $z$-scored features. All $hat(beta)$ values are negative, indicating that evolution favors moves that reduce physicochemical mismatch, avoid topology disruption, and (weakly) prefer target amino acids already serviced by nearby codons. The tRNA proxy coefficient is small and non-significant (LR $= 0.54$, $p = 0.46$).
  ],
) <tbl:s-condlogit-coefs>

== Likelihood ratio tests

#figure(
  table(
    columns: (auto, auto, auto, auto, auto),
    align: (left, left, center, center, center),
    inset: 5pt,
    stroke: (x, y) => if y == 0 { (bottom: 0.7pt) } else { none },
    table.header([Restricted], [Full], [LR], [df], [$p$]),
    [M1 (phys)], [M3 (phys+topo)], [103.7], [1], [$lt.lt 10^(-10)$],
    [M2 (topo)], [M3 (phys+topo)], [69.1], [1], [$lt.lt 10^(-10)$],
    [M3 (phys+topo)], [M4 (full)], [0.54], [1], [$0.46$],
  ),
  caption: [
    Likelihood-ratio tests for nested conditional logit models. Both topology (added to M1) and physicochemistry (added to M2) provide highly significant improvements. The tRNA-complexity proxy does not improve on the phys+topo model.
  ],
) <tbl:s-condlogit-lrt>

== Confounding diagnostic

Across $approx 78{,}000$ candidate moves pooled from all choice sets, the Spearman correlation between $Delta_"phys"$ and $Delta_"topo"$ is $rho = 0.15$, indicating that the two predictors carry largely independent information. Moves that reduce physicochemical mismatch are only slightly more likely to also preserve topology, and the conditional logit framework accounts for any residual collinearity through simultaneous estimation.

== Observed move percentile ranks

Under the best model (M3), observed natural reassignments rank at the 89.5th percentile on average (mean rank 134 out of $approx 1{,}280$ candidates). Notable individual events: UGA$arrow.r$Trp ranks consistently above the 98th percentile across all tables where it occurs, indicating that this reassignment is among the most "predictable" given the combined phys+topo score. The one clear outlier is the yeast mitochondrial CUU$arrow.r$Thr reassignment (30th percentile), consistent with translation table 3's status as the sole marginal exception in the per-table optimality analysis.

== Order-averaging implementation

For tables with $k > 1$ reassignment events, the temporal ordering is unknown. We marginalized the conditional logit likelihood over all $k!$ orderings, computing $L_"table" = (1\/k!) sum_sigma product_s P(m_(sigma(s))^* | cal(N)(C_(sigma,s)))$ using log-sum-exp for numerical stability. For the largest table (yeast mitochondrial, $k = 6$, $6! = 720$ possible orderings), we used a random subset of 24 orderings. For all tables with $k lt.eq 5$ ($lt.eq 120$ orderings), all orderings were enumerated exactly. We verified that AICc differences and coefficient estimates were stable to increasing the number of sampled orderings (6 vs 24 subsets for the yeast mitochondrial table; results unchanged to two significant digits), indicating that Monte Carlo noise from order-subsampling does not affect the qualitative conclusions.

// ============================================================
= Software and reproducibility <sec:s-software>

All analyses were performed using the `codon-topo` Python package (version 0.3.0), with:
- Python 3.11, NumPy 1.24+, SciPy 1.10+
- R 4.5, ggplot2, ggpubr, viridis, patchwork (for figures)
- tRNAscan-SE 2.0.12 with Infernal 1.1.4
- Random seed: 135325 (all Monte Carlo analyses)

Analyses are fully reproducible via:
```
pip install -e ".[dev]"
codon-topo all --output-dir=./output --seed=135325
python scripts/generate_tables.py
Rscript src/codon_topo/visualization/R/strengthened_figures.R
```

The complete test suite (367 tests) can be run with `python3.11 -m pytest tests/`.

// ============================================================
// REFERENCES (shared with main manuscript)
// ============================================================

#bibliography("references.bib", title: "References", style: "elsevier-harvard")
