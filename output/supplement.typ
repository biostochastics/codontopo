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
    [hypercube_coloring_optimality], [Supported], [Cross-metric sensitivity: Grantham $p = 0.006$, Miyata $p < 0.001$, polar requirement $p = 0.003$, Kyte--Doolittle $p = 0.001$ (block-preserving null, $n = 10,000$). Stop penalty sensitivity (0/150/215/300): immaterial.],
    [per_table_optimality_preservation], [Supported], [Per-table block-preserving null applied to all 27 NCBI tables; significant fraction (BH-FDR $p < 0.05$) and per-table quantiles are reported via pipeline-rendered statistics in the main text and the per-table CSV (`output/tables/T4_per_table_optimality.csv`). Translation table 3 (yeast mito) is the marginal exception.],
    [optimality_rho_robustness], [Supported], [$rho$-sweep at $n = 10,000$: $p lt.eq 0.006$ at all $rho$ values; effect-size $z$ increases monotonically from $rho = 0$ to $rho = 1$.],
    [topology_avoidance_depletion], [Supported], [Permutation $p lt.eq 10^(-4)$; hypergeometric $p = 1.6 times 10^(-8)$ ($Q_6$) and $p = 1.3 times 10^(-6)$ ($K_4^3$). 6 of 28 de-duplicated events break topology vs ~73% of the candidate landscape (3.4-fold depletion). Clade-exclusion sensitivity: all $p < 10^(-5)$.],
    [trna_enrichment_reassigned_aa], [Suggestive], [MIS worst-case $p = 0.045$. 18 tRNAscan-SE-verified assemblies (15 variant + 3 standard controls), 24 pairings across 5 variant codes.],
    [bit_position_bias_weighted], [Exploratory], [Uniform $p = 0.006$ inflated by non-independence. De-duplicated $p = 0.075$.],
    [mechanism_boundary_conditions], [Exploratory], [Three-tier pattern: duplication / stem shortening / modification. Descriptive.],
    [atchley_f3_serine_convergence], [Exploratory], [Serine $F_3 = -4.760$, 2.24 SD below mean. Complementary, not independent.],
    [variant_code_disconnection_catalogue], [Exploratory], [4 variant-code cases at $epsilon = 1$ in $Q_6$ (default encoding): Thr (Table 3, yeast mito), Leu (Table 16, chlorophycean mito), Ala (Table 26, _Pachysolen_), Ser (Table 12, _Candida_). Separately, Table 32 (Balanophoraceae plastid, UAG→Trp) creates a 2-fold Trp pair (UGG, UAG) that remains connected at ε=1 but breaks the bit-5 two-fold filtration — a filtration finding, not a disconnection.],
    [kras_fano_clinical_prediction], [Falsified], [$p = 1.0$ across all 6 G12 variants. $n = 1,670$ MSK-IMPACT mutations.],
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

The coloring optimality result was tested under all 24 base-to-bit bijections from ${C,U,A,G}$ to $"GF"(2)^2$. All 24 encodings yield significant optimality ($p < 0.05$ under the block-preserving null with $n = 1,000$). The mean quantile across encodings is 1.8%, confirming that the result is not an artifact of the default encoding choice.

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

The reassignment database comprises all codon reassignment events across the 27 NCBI translation tables, relative to the standard code (Table 1). Each event records the codon, source amino acid, target amino acid, and Hamming distance to the nearest codon already encoding the target amino acid in the standard code. De-duplication to unique (codon, target amino acid) pairs yields the event set used in the topology-avoidance analysis. Exact counts (which depend on the most recent NCBI gc.prt revision) are reported inline in the main text via pipeline-rendered statistics from `manuscript_stats.json`.

The complete database is provided as `output/tables/T10_reassignment_db.csv`.


// ============================================================
= Topology avoidance: clade-exclusion sensitivity <sec:s-clade>

To test whether the topology avoidance result is driven by any single phylogenetic clade, we iteratively excluded:
- All ciliate reassignments
- All metazoan mitochondrial reassignments
- All CUG-clade yeast reassignments
- All chlorophycean reassignments

In every exclusion, the depletion remains highly significant ($p < 10^(-5)$), confirming that the $approx 3.4$-fold ($Q_6$) and $approx 3.1$-fold ($H(3,4)$) depletion of topology-breaking changes is a pan-taxonomic pattern, not an artifact of any single lineage. Excluding yeast mitochondrial (table 3) actually strengthens the depletion (rate drops from 21.4% to 8.3%; hypergeometric $p$ from $1.6 times 10^(-8)$ to $3.6 times 10^(-11)$), because table 3 contributes 4 of the 6 topology-breaking events. Detailed per-clade counts are available in `output/phylogenetic_sensitivity.json`.


// ============================================================
= Topology-breaking definitions: $2 times 2$ audit <sec:s-topology-defs>

We define two notions of "topology-breaking" candidate move:

+ *New disconnection in a previously connected family* (legacy primary $Q_6$ definition): a candidate move makes some amino acid disconnected at $epsilon = 1$ that was connected in the standard code (i.e., the amino acid was not previously in the disconnection catalogue but is after the move).

+ *Increase in components* ($Delta beta_0 > 0$, the conditional-logit feature): the total number of connected components, summed across amino-acid codon graphs, strictly increases:
  $ Delta_"topo" = sum_a beta_0(G_a^"after") - sum_a beta_0(G_a^"before") > 0. $

Both definitions are reported under both $Q_6$ adjacency (Hamming-1 in the default $"GF"(2)^6$ encoding) and $H(3,4)$ adjacency (full single-nucleotide adjacency, encoding-independent), giving four cells. All four share the same denominators (1{,}280 candidate moves, 28 de-duplicated observed events).

#table(
  columns: (auto, auto, auto, auto, auto, auto),
  align: (left, left, right, right, right, right),
  inset: 6pt,
  stroke: (x, y) => if y == 0 { (bottom: 0.7pt) } else { none },
  table.header(
    [*Adjacency*], [*Definition*], [*$K\/N$ possible*], [*$x\/n$ obs*], [*Hyper. $p$*], [*RR (95% CI)*],
  ),
  [$Q_6$], [new disconnection], [931 / 1280], [6 / 28], [$1.6 times 10^(-8)$], [0.29 (0.14--0.60)],
  [$Q_6$], [$Delta beta_0 > 0$], [963 / 1280], [7 / 28], [$2.2 times 10^(-8)$], [0.33 (0.17--0.63)],
  [$H(3,4)$], [new disconnection], [822 / 1280], [5 / 28], [$5.0 times 10^(-7)$], [0.28 (0.13--0.62)],
  [$H(3,4)$], [$Delta beta_0 > 0$], [846 / 1280], [6 / 28], [$1.3 times 10^(-6)$], [0.32 (0.16--0.66)],
)

The four cells agree qualitatively: depletion is highly significant under all four combinations, with risk ratios in the range 0.28--0.33 and hypergeometric $p < 10^(-5)$ throughout. The main text uses the new-disconnection definition for $Q_6$ (matching the legacy reported counts) and the $Delta beta_0 > 0$ definition for $H(3,4)$ (matching both the conditional-logit feature and the legacy reported counts). The full machine-readable audit is in `output/topology_avoidance.json` under the `definitions_audit` key.


// ============================================================
= Topology avoidance: $Q_6$ encoding-sweep sensitivity <sec:s-encoding-sweep>

The $H(3,4)$ Hamming graph is encoding-independent (every two-bit bijection from ${A, C, G, U}$ to ${0, 1}^2$ produces the same nucleotide-level graph). The $Q_6$ subgraph, however, depends on the encoding because the partition into Hamming-1 (single-bit-change) edges versus Hamming-2 (within-nucleotide diagonal) edges is bijection-specific. To test whether the $Q_6$ topology-avoidance result is robust to the encoding choice, we recomputed the $Q_6$ candidate-landscape rate, observed rate, depletion fold, and hypergeometric $p$-value under all 24 base-to-bit bijections, holding the same 1{,}280 candidate moves and 28 observed events.

Across all 24 encodings:

- candidate-landscape rate: min 36.0%, median 69.9%, max 72.7%
- observed rate: min 21.4%, median 28.6%, max 35.7%
- depletion fold: min 1.01, median 2.45, max 3.39
- hypergeometric $p$: min $1.6 times 10^(-8)$, median $6.1 times 10^(-6)$, max 0.572

The default encoding (C=00, U=01, A=10, G=11) gives the largest depletion and the smallest $p$. Eight of 24 encodings give a candidate-landscape rate of approximately 36%, under which the observed rate of 21--36% does not significantly differ from candidate (p > 0.5). The Q_6 result is therefore not encoding-invariant. We accordingly present the encoding-independent $H(3,4)$ result as the primary topology-avoidance test, with $Q_6$ reported as a representation-specific decomposition for continuity with the broader $"GF"(2)^6$ framework. The full per-encoding sweep is in `output/topology_avoidance.json` under `Q6_encoding_sweep`.


// ============================================================
= Reassignment candidate-universe denominator sensitivity <sec:s-denominator>

The conditional-logit and topology-avoidance analyses both use a candidate-universe of $cal(M)(C) = {(x, y) : x in cal(C), y in cal(A)_20 union {"Stop"}, y eq.not C(x)}$, with $abs(cal(M)(C)) = 64 times 20 = 1280$. Reviewer R1.D requested an explicit comparison of alternative universes:

#table(
  columns: (auto, auto, auto),
  align: (left, right, left),
  inset: 6pt,
  stroke: (x, y) => if y == 0 { (bottom: 0.7pt) } else { none },
  table.header([*Universe*], [*$abs(cal(M))$*], [*Definition*]),
  [U1: 21-label, no identity (primary)], [1{,}280],
    [$y in cal(A)_20 union {"Stop"}$, $y eq.not C(x)$; each codon has exactly 20 alternative labels],
  [U2: AA-only, no identity], [1{,}219],
    [$y in cal(A)_20$, $y eq.not C(x)$; 61 sense codons get 19 AA alternatives, 3 stop codons get 20 AA alternatives],
  [U3: AA, with no-ops], [1{,}280],
    [$y in cal(A)_20$, no identity restriction; 64 candidates are no-ops $y = C(x)$],
  [U4: stop-inclusive with no-ops], [1{,}344],
    [$y in cal(A)_20 union {"Stop"}$, no identity restriction; 64 no-ops included],
)

We adopt U1 as the primary universe because it (i) has a uniform alternative count per codon, (ii) cleanly excludes identity moves which contribute no signal, and (iii) admits stop-codon reassignment which is biologically attested. Topology-avoidance results under U2 and U4 are qualitatively identical (depletion remains $p < 10^(-5)$); the per-cell counts differ only by the constant rescaling $K_2 = K_1 - text("(stop-target candidates)")$ for U2. Detailed numbers are in `output/topology_avoidance.json` under `denominator_sensitivity`.


// ============================================================
= Conditional logit: IIA assumption <sec:s-iia>

Conditional logit assumes Independence of Irrelevant Alternatives (IIA): for any two candidate moves $m_1$ and $m_2$ in a choice set $cal(N)$, the relative probability $P(m_1) \/ P(m_2)$ is the same regardless of which other candidates are in $cal(N)$. In the reassignment context, candidate moves are not exchangeable: a move that targets an amino acid already serviced by many codons is plausibly more substitutable with similar moves than the IIA structure allows.

We adopt IIA here because the goal is _explanatory_ rather than _predictive_: the test asks whether topology adds explanatory value beyond physicochemical cost (LR test, ΔAICc) within the same candidate set, not whether the model accurately predicts which specific reassignment will occur next. The relative-probability ratios that IIA constrains are not the quantities we report; the LR statistics depend only on whether the observed event occupies a high-likelihood position within the candidate set, which is robust to substitution patterns among non-observed candidates. A mixed logit relaxing IIA would be more appropriate if the goal were prediction; future work could pursue this once a larger event set permits identification of mixed-logit covariance parameters.


// ============================================================
= Conditional logit: clade-exclusion sensitivity <sec:s-condlogit-clade>

Reviewer R2.M1 (RIGORA Major #1) noted that the topology-avoidance hypergeometric/permutation tests have a clade-exclusion sensitivity analysis (above, @sec:s-clade) but the conditional logit does not. To address this, we refit M1--M4 under each of the seven clade-exclusion regimes used for the topology-avoidance test, removing the indicated tables and re-running the full conditional-logit pipeline (build candidate sets $arrow.r$ enumerate event orderings $arrow.r$ vectorized fit $arrow.r$ ΔAICc).

The full per-regime ΔAICc(M1$arrow.r$M3), ΔAICc(M2$arrow.r$M3), and ΔAICc(M3$arrow.r$M4) are in `output/condlogit_clade_sensitivity.json`. The interpretation depends on whether the topology coefficient remains decisively favored across all exclusions: if it does, the topology effect is robust to phylogenetic non-independence; if it concentrates in one or two clades, the framing in @sec:res-condlogit and @sec:disc-evidence requires softening.


// ============================================================
= Per-table optimality: standard-code-proximity audit <sec:s-proximity>

The methodological-nuance concern noted by Reviewer R2.6 is that a variant code which differs from the standard code by only a few reassignments may have a per-table block-preserving null distribution dominated by permutations near the standard code. In that case, "table $X$ falls in the bottom 5% of permutations preserving $X$'s block structure" partly tests whether $X$ is close to the standard code (which it is, by construction), not whether $X$ is independently optimal.

To address this, for each NCBI translation table we computed two quantities alongside the unconditional per-table quantile: (i) the Hamming distance from the standard code (number of codons with different AA labels) for both the observed variant code and each block-preserving null draw; and (ii) the variant's null quantile *conditional* on null draws that are within $plus.minus 2$ codons of the variant's distance from standard. If the conditional and unconditional quantiles agree, the per-table optimality is independent of standard-code proximity; if they diverge sharply, the per-table test is largely a proximity test.

Detailed per-table results are in `output/coloring_optimality.json` under `per_table_proximity_audit`. The audit confirms that for variant tables differing from standard by 3 or more codon reassignments, the conditional-on-distance quantile is similar to the unconditional quantile, supporting the interpretation that variant codes preserve error-minimization structure independently of their similarity to the standard code.


// ============================================================
= KRAS--Fano clinical prediction: detailed results <sec:s-kras>

The conjecture that XOR ("Fano") relationships in $"GF"(2)^6$ predict enrichment of specific amino acids at KRAS G12 co-mutation sites was tested against 1,670 KRAS mutations from the MSK-IMPACT dataset @zehir2017.

For each of the 6 KRAS G12 variant types (G12D, G12V, G12C, G12A, G12R, G12S), we identified the XOR-predicted amino acid partners and tested for co-mutation enrichment via Fisher's exact test with Bonferroni correction. All 6 tests yielded $p = 1.0$, with odds ratios near 1.0.

This result cleanly separates code-level error-minimization (which operates on the amino acid assignment structure) from mutation-level algebraic predictions (which would require DNA polymerase errors to respect the binary encoding --- a biologically implausible mechanism).


// ============================================================
= Event-level conditional logit model: full details <sec:s-condlogit>

== Model specification and candidate universe

The conditional logit model treats each observed reassignment as a discrete choice from the set of all single-codon reassignments available at the current code state. For a code with 64 codons and $approx 21$ possible amino acid/stop labels, the candidate set $cal(N)(C)$ contains $approx 1,280$ moves (64 codons $times$ 20 alternative labels, excluding identity assignments). All 1,280 candidates are evaluated at each step regardless of biological plausibility; the model's purpose is to test whether the observed moves are statistically distinguishable from uniform sampling given the three feature classes. As a conditional logit, the model implies an independence of irrelevant alternatives (IIA) structure within each choice set. Here the goal is not mechanistic completeness but a controlled test of whether observed moves are statistically enriched for low $Delta_"phys"$ and low $Delta_"topo"$ relative to a uniform baseline over the same candidate set; a mixed logit relaxing IIA would be appropriate if the goal were prediction rather than hypothesis testing.

== Feature definitions

*Local physicochemical mismatch change ($Delta_"phys"$)*: For a reassignment of codon $c$ from amino acid $a$ to $a'$, this is the change in the sum of Grantham (1974) distances across all Hamming-1 edges incident to $c$:

$ Delta_"phys" = sum_({c,c'}: d(c,c')=1) [Delta(a', "code"(c')) - Delta(a, "code"(c'))] $

This captures the local impact on error-minimization at the reassigned position. When a neighboring codon is assigned Stop, we set $Delta(a, "Stop")$ equal to the maximum Grantham distance (215), consistent with the stop-penalty convention used in the coloring objective (Methods §2.2 of the main manuscript).

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
    [M1], [$Delta_"phys"$], [$-0.0049$], [$-1.33$],
    [M2], [$Delta_"topo"$], [$-3.584$], [$-1.67$],
    [M3], [$Delta_"phys"$], [$-0.0043$], [$-1.17$],
    [M3], [$Delta_"topo"$], [$-3.316$], [$-1.54$],
    [M4], [$Delta_"phys"$], [$-0.0043$], [$-1.17$],
    [M4], [$Delta_"topo"$], [$-2.986$], [$-1.39$],
    [M4], [$Delta_"tRNA"$], [$-0.203$], [$-0.21$],
  ),
  caption: [
    Conditional logit coefficient estimates from the 27-table re-fit. Raw coefficients are on the original feature scale; normalized coefficients are on $z$-scored features. All $hat(beta)$ values are negative, indicating that evolution favors moves that reduce physicochemical mismatch, avoid topology disruption, and (weakly) prefer target amino acids already serviced by nearby codons. The tRNA proxy coefficient is small and non-significant (LR $= 0.12$, $p = 0.73$).
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
    [M1 (phys)], [M3 (phys+topo)], [110.4], [1], [$lt.double 10^(-10)$],
    [M2 (topo)], [M3 (phys+topo)], [91.2], [1], [$lt.double 10^(-10)$],
    [M3 (phys+topo)], [M4 (full)], [0.12], [1], [$0.73$],
  ),
  caption: [
    Likelihood-ratio tests for nested conditional logit models. Both topology (added to M1) and physicochemistry (added to M2) provide highly significant improvements. The tRNA-complexity proxy does not improve on the phys+topo model.
  ],
) <tbl:s-condlogit-lrt>

== Confounding diagnostic

Across $approx 84,000$ candidate moves pooled from all choice sets, the Spearman correlation between $Delta_"phys"$ and $Delta_"topo"$ is $rho = 0.15$, indicating that the two predictors carry largely independent information. Moves that reduce physicochemical mismatch are only slightly more likely to also preserve topology, and the conditional logit framework accounts for any residual collinearity through simultaneous estimation.

== Observed move percentile ranks

Under the best model (M3), observed natural reassignments rank at the 89.5th percentile on average (mean rank 134 out of $approx 1,280$ candidates). Notable individual events: UGA$arrow.r$Trp ranks consistently above the 98th percentile across all tables where it occurs, indicating that this reassignment is among the most "predictable" given the combined phys+topo score. The one clear outlier is the yeast mitochondrial CUU$arrow.r$Thr reassignment (30th percentile), consistent with translation table 3's status as the sole marginal exception in the per-table optimality analysis.

== Order-averaging implementation

For tables with $k > 1$ reassignment events, the temporal ordering is unknown. We marginalized the conditional logit likelihood over all $k!$ orderings, computing $L_"table" = (1\/k!) sum_sigma product_s P(m_(sigma(s))^* | cal(N)(C_(sigma,s)))$ using log-sum-exp for numerical stability. For all tables with $k lt.eq 6$ events (which covers every table in the dataset, including the largest --- yeast mitochondrial, $k = 6$, $6! = 720$ orderings), we enumerate all orderings exactly. The implementation precomputes feature matrices once per (table, ordering, step) into stacked numpy arrays, after which each Nelder--Mead / L-BFGS-B function evaluation reduces to a small batch of matrix multiplies and `scipy.special.logsumexp` calls; the four nested models are fit concurrently on a single shared bundle via `joblib` threads to avoid memory duplication. The conservative random-subset implementation reported in earlier drafts (24 of 720 orderings) is no longer used; AICc and coefficient estimates from the exact-enumeration path agree with the random-subset approximation to within Monte Carlo noise.

// ============================================================
= Software and reproducibility <sec:s-software>

All analyses were performed using the `codon-topo` Python package (version 0.4.0), with:
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
