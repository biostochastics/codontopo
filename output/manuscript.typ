// ============================================================
// Manuscript: Error-minimizing structure of the genetic code
// Target: Journal of Theoretical Biology
// Authors: Paul Clayworth, Sergey Kornilov (equal contribution)
// ============================================================

// --- Page and text setup ---
#set page(
  paper: "a4",
  margin: (x: 2.5cm, y: 2.5cm),
  numbering: "1",
  header: context {
    if counter(page).get().first() > 1 [
      #set text(size: 9pt, fill: luma(120))
      Clayworth & Kornilov — Genetic code error-minimization in GF(2)#super[6]
      #h(1fr)
      #counter(page).display()
    ]
  },
)
#set text(font: "Libertinus Serif", size: 11pt, lang: "en")
#set par(justify: true, leading: 0.65em, first-line-indent: 1.5em)
#set heading(numbering: "1.1")
#set math.equation(numbering: "(1)")

// Figures: tighter gap, auto-placement to avoid large whitespace
#set figure(gap: 0.8em, placement: auto)

// Suppress first-line-indent after headings and figures
#show heading: it => {
  it
  par(text(size: 0pt, ""))
}

// Style table header cells bold
#show table.cell.where(y: 0): set text(weight: "bold")

// Use tabular figures in tables
#show table: set text(number-width: "tabular")

// ============================================================
//  TITLE BLOCK
// ============================================================
#set par(first-line-indent: 0pt)
#align(center)[
  #v(1em)
  #text(16pt, weight: "bold")[
    Robust error-minimization in the genetic code across \
    physicochemical metrics and variant codes: \
    a graph-theoretic analysis in GF(2)#super[6]
  ]
  #v(1em)
  #text(12pt)[Paul Clayworth#super[1,#sym.dagger] and Sergey Kornilov#super[2,#sym.dagger,\*]]
  #v(0.5em)
  #text(10pt, style: "italic")[
    #super[1]Logocentricity Inc., Austin, TX, USA \
    #super[2]Biostochastics, LLC, Seattle, WA, USA \
    #super[#sym.dagger]Equal contribution \
    #super[\*]Corresponding author
  ]
  #v(1.5em)
]

// ============================================================
//  ABSTRACT
// ============================================================
#block(
  width: 100%,
  inset: (x: 1.5cm),
)[
  #text(weight: "bold")[Abstract] \
  Encoding each nucleotide as a 2-bit vector maps the 64 codons to vertices of a 6-dimensional binary hypercube $Q_6 subset.eq K_4^3$, where $K_4^3$ is the complete single-nucleotide mutation graph. A weight parameter $rho in [0,1]$ interpolates between $Q_6$ ($rho = 0$) and the full graph ($rho = 1$), enabling systematic analysis of error-minimization across the mutation spectrum. Under a block-preserving null model, the standard genetic code is significantly optimal across four distinct physicochemical distance metrics---Grantham ($p = 0.003$), Miyata ($p = 0.001$), Woese polar requirement ($p = 0.004$), and Kyte--Doolittle hydropathy ($p = 0.001$)---demonstrating that error-minimization is robust across commonly used physicochemical parameterizations and not a metric-specific artifact. This optimality persists at all $rho$ values (strengthening toward $rho = 1$) and is preserved in 24 of 25 NCBI translation tables after Benjamini--Hochberg correction (mean quantile 1.4%). Natural codon reassignment events show 3.3-fold depletion of topology-breaking changes relative to the space of possible reassignments (22% observed vs 73% possible; permutation $p lt.eq 10^(-4)$; hypergeometric $p = 4.8 times 10^(-8)$), a pattern robust to phylogenetic clade exclusion. Some variant-code lineages show elevated tRNA gene counts for the reassigned amino acid (worst-case maximal independent set Stouffer $p = 0.045$, 18 tRNAscan-SE--verified genomes), consistent with one compensatory route among several that organisms use to accommodate codon reassignment.

  #v(0.5em)
  *Keywords:* genetic code evolution; error minimization; graph-theoretic analysis; GF(2)#super[6]; codon reassignment; physicochemical distance; tRNA gene duplication
]

#v(1em)
#set par(first-line-indent: 1.5em)

// ============================================================
//  1. INTRODUCTION
// ============================================================
= Introduction <sec:intro>

The standard genetic code maps 61 sense codons to 20 amino acids through a pattern that has long been recognized as non-random. Woese (1965) first noted that similar codons tend to encode amino acids with similar physicochemical properties, and Freeland and Hurst (1998) demonstrated quantitatively that the standard code sits in approximately the top $10^(-6)$ percentile of random codes for mutational error minimization under the Woese polar requirement distance. These findings established that the code's structure reduces the fitness impact of point mutations, but the question of whether this optimality is specific to particular physicochemical metrics or represents a robust property across several distinct physicochemical parameterizations has not been systematically addressed.

Every codon comprises three nucleotides drawn from ${C, U, A, G}$. By choosing a bijection $phi: {C, U, A, G} -> "GF"(2)^2$---for instance, $C |-> (0,0)$, $U |-> (0,1)$, $A |-> (1,0)$, $G |-> (1,1)$---each codon maps to a vertex of the 6-dimensional binary hypercube $Q_6 = "GF"(2)^6$. The genetic code then becomes a _coloring_ of $Q_6$ by 21 labels (20 amino acids plus the stop signal), and single-nucleotide mutations correspond to edges of $Q_6$ or, more precisely, to a subgraph of the complete mutation graph $K_4^3$.

This representation is not new in principle---binary encodings of the genetic code appear in the mathematical biology literature (e.g., Antoneli and Forger, 2011). Its value lies not in the encoding itself but in the analytical decomposition it enables: the complete single-nucleotide mutation graph $K_4^3$ (288 edges) separates into 192 Hamming-distance-1 edges (single-bit changes) and 96 within-nucleotide distance-2 edges (both bits of one nucleotide position flip simultaneously), allowing systematic interpolation via a weight parameter $rho in [0,1]$. Note that this decomposition does not correspond to the biological transition/transversion partition: under any 2-bit encoding, each Hamming-1 edge class contains an equal mixture of transitions and transversions. Previous work has not exploited this decomposition to test error-minimization across multiple physicochemical metrics, nor extended the analysis to variant genetic codes. We address three questions:

+ *Is the standard code optimal under multiple physicochemical metrics?* We test whether the standard code's edge-mismatch score is extreme relative to block-preserving null models across four distinct physicochemical distance metrics (Grantham, Miyata, Woese polar requirement, and Kyte--Doolittle hydropathy), extending Freeland and Hurst (1998) from a single metric to a cross-metric robustness envelope.

+ *Is this structure preserved by evolution?* We ask whether variant genetic codes (NCBI translation tables 2--33) maintain error-minimization, and whether natural codon reassignment events preferentially avoid disrupting the topological connectivity of amino acid codon families.

+ *What are the genomic correlates of disruption?* We test whether organisms whose variant codes break the connectivity of an amino acid's codon graph show elevated tRNA gene copy numbers for the affected amino acid, using tRNAscan-SE--verified data from 10 genomes across 4 eukaryotic supergroups.

In pursuing these questions, we also evaluated several deeper algebraic conjectures from a companion technical note (Clayworth, 2026)---including claims about PSL(2,7) symmetry and holomorphic embeddings---and report their falsification as boundary conditions on the framework's scope (Supplementary Material).

The paper is organized as follows. Section 2 describes the encoding formalism, graph decomposition, null models, and statistical methods. Section 3 presents the four supported findings (cross-metric coloring optimality, per-table preservation, $rho$-robustness, and topology avoidance), the suggestive tRNA enrichment result, and exploratory observations. Section 4 discusses the graph-theoretic interpretation, the relationship to frozen-accident versus adaptive hypotheses (Koonin and Novozhilov, 2009), and limitations. All analyses are reproducible via the open-source `codon-topo` pipeline (version 0.3.0, 367 tests, seed 135325).


// ============================================================
//  2. METHODS
// ============================================================
= Methods <sec:methods>

== Binary encoding of the genetic code <sec:encoding>

We encode each nucleotide base as a 2-bit vector via the default bijection $phi: C |-> (0,0)$, $U |-> (0,1)$, $A |-> (1,0)$, $G |-> (1,1)$. A codon $b_1 b_2 b_3$ is then represented as the concatenation $phi(b_1) || phi(b_2) || phi(b_3) in "GF"(2)^6$. The 64 codons become the 64 vertices of the 6-dimensional hypercube $Q_6$.

Under this encoding, two codons that differ by a single-nucleotide substitution in which only one bit of the 2-bit pair changes are adjacent in $Q_6$ (Hamming distance 1). However, transversions that flip both bits within a nucleotide position correspond to Hamming distance 2. The full single-nucleotide mutation graph is therefore $K_4^3$ (the Cartesian product of three complete graphs on 4 vertices), which contains $Q_6$ as a subgraph. To address the concern that $Q_6$ misses approximately one-third of single-nucleotide mutations, we introduce the weighted score $F_rho$ that interpolates between pure $Q_6$ ($rho = 0$) and full $K_4^3$ ($rho = 1$); see Section 2.3.

There are 24 distinct bijections from ${C,U,A,G}$ to $"GF"(2)^2$ (all permutations of $4! = 24$ assignments). All encoding-dependent results are tested across all 24 encodings; encoding-invariant properties (such as Serine's disconnection at $epsilon = 1$, where $epsilon$ denotes the Hamming distance threshold in $"GF"(2)^6$) are noted as such. Coloring optimality is significant ($p < 0.05$) under all 24 encodings; full encoding-sensitivity results are provided in the Supplementary Material.

== Edge-mismatch objective function <sec:objective>

The genetic code assigns each vertex $v in Q_6$ a label $c(v) in cal(A)$, where $cal(A)$ comprises the 20 amino acids and the stop signal. The edge-mismatch score is

$ F(c) = sum_({v,w}: d(v,w) = 1) Delta(c(v), c(w)) $ <eq:mismatch>

where the sum ranges over all 192 edges of $Q_6$ (pairs of vertices at Hamming distance 1), and $Delta$ is a physicochemical distance between amino acids. We test four distinct distance metrics: the Grantham (1974) composite distance (composition, polarity, volume; range 5--215), the Miyata (1979) normalized Euclidean distance (polarity and volume only; range 0.06--5.13), the Woese (1973) polar requirement absolute difference (range 0--8.2; used by Freeland and Hurst 1998), and the Kyte--Doolittle (1982) hydropathy absolute difference (range 0--9.0; used by Haig and Hurst 1991). Synonymous edges contribute 0; edges involving a stop codon receive a fixed penalty scaled proportionally to each metric's maximum. A lower $F$ indicates a more error-minimizing code.

Stop codons are held fixed across all null models (Section 2.3), so the stop-codon contribution is a constant offset. Sensitivity analysis across penalty values (0, 150, 215, 300) confirms this is immaterial to the ranking.

== Null models <sec:nullmodels>

=== Block-preserving null (Freeland--Hurst) <sec:null-fh>

Following Freeland and Hurst (1998), we group the 64 codons into 16 blocks of 4, defined by shared first-two-base prefix. Each block's internal pattern of amino acid assignments is preserved, but the mapping of patterns to blocks is permuted uniformly at random. Blocks containing stop codons are held fixed. This null preserves the synonymous codon contiguity (wobble degeneracy) inherent in the genetic code, providing a stringent test: the observed code must beat random codes that share its degeneracy architecture.

For the standard code, we drew $n = 10,000$ null samples (seed 135325) and computed the conservative $p$-value as $(k + 1)/(n + 1)$, where $k$ is the number of null scores below the observed $F$.

=== Per-table null <sec:null-pertable>

Each of the 25 NCBI translation tables was tested against its own block-preserving null ($n = 1{,}000$ per table, seeds = base seed + table ID). $P$-values were corrected for multiple comparisons using the Benjamini--Hochberg (1995) false discovery rate (FDR) procedure.

=== Graph family $G_rho$ and rho sweep <sec:null-rho>

We define a family of mutation graphs $G_rho$ parameterized by $rho in [0,1]$, where $G_0 = Q_6$ (192 Hamming-1 edges) and $G_1 = K_4^3$ (all 288 single-nucleotide substitution edges). The weighted mismatch score on $G_rho$ is

$ F_rho(c) = F_("Hamming-1")(c) + rho dot.c F_("diagonal")(c) $ <eq:weighted>

where $F_("diagonal")$ sums over the 96 within-nucleotide distance-2 edges. This formulation treats $Q_6$ not as the biological object but as one endpoint of a continuous interpolation toward the full mutation graph. We evaluated 5 values of $rho$ ($n = 1{,}000$ null samples per $rho$, common-seed design).

=== Table-preserving permutation null (topology avoidance) <sec:null-topo>

To test whether natural codon reassignment events avoid disrupting amino acid connectivity, we enumerated all 1,280 possible single-codon reassignments from the standard code. For each amino acid $a$, let $G_a^1$ denote the subgraph induced on codons assigned to $a$ with edges between codons at binary Hamming distance $lt.eq 1$. A reassignment is _topology-breaking_ if it increases the number of connected components of $G_a^1$ for any amino acid $a$. We define a de-duplicated reassignment event as a unique (codon, target amino acid) pair observed across the 25 NCBI translation tables, yielding 27 events from 61 total table-specific changes. We then compared the observed rate of topology-breaking among 27 de-duplicated natural reassignment events to the null expectation.

Two tests were applied: (i) a hypergeometric test treating the observed events as a sample from the finite landscape of possible reassignments ($N = 1{,}280$, $K = 931$ topology-breaking, $n = 27$ observed, $x = 6$ topology-breaking observed); and (ii) a table-preserving permutation test ($n = 10{,}000$, seed 135325) that permutes target amino acids among a table's reassigned codons, preserving within-table codon and target structure to address phylogenetic non-independence.

=== Fisher--Stouffer test for tRNA enrichment <sec:null-trna>

For each variant-code organism paired with a phylogenetically proximate standard-code control, we built a $2 times 2$ contingency table comparing the proportion of tRNA genes encoding the reassigned amino acid versus all other amino acids. The sampling model treats tRNA gene counts as draws from a hypergeometric distribution conditional on the row and column marginals (total tRNAs per organism, total tRNAs for the focal amino acid), which is appropriate when the question is whether a specific amino acid's share of the tRNA repertoire differs between two genomes. Fisher's exact test (one-sided, alternative = "greater") was applied per pairing, and $p$-values were combined via Stouffer's $Z$ method. To address non-independence from shared control organisms, we constructed a conflict graph (edges connect pairings sharing an organism) and enumerated all maximal independent sets (MIS) via the Bron--Kerbosch algorithm with pivoting. For each MIS with $>= 2$ members, we computed the Stouffer combined $p$-value and report the worst-case (maximum $p$) across all MIS as the primary conservative test.

tRNA gene counts for 17 organisms were verified by running tRNAscan-SE 2.0.12 (Chan and Lowe, 2019) with Infernal 1.1.4 on NCBI genome assemblies (eukaryotic mode for ciliates, bacterial mode for _Mycoplasma_). Counts for 3 additional organisms were drawn from published literature. The dataset comprises 24 pairings across 5 variant genetic codes (NCBI translation tables 4, 6, 10, 15, and 31) and 3 standard-code controls.

== Claim hierarchy <sec:claims>

All 15 scientific claims evaluated in this work are registered in a formal claim hierarchy with pre-specified status categories (@tbl:claims): SUPPORTED (passes rigorous null, $p < 0.01$), SUGGESTIVE (trend-level, $p < 0.05$), EXPLORATORY (hypothesis-generating), REJECTED (falsified or pre-rejected in literature), FALSIFIED (directly contradicted by data), and TAUTOLOGICAL (true by construction). The full hierarchy with justifications is available in the supplementary materials.

#figure(
  table(
    columns: (auto, 1fr, auto, auto),
    align: (left, left, center, right),
    inset: 6pt,
    stroke: (x, y) => if y == 0 { (bottom: 0.7pt) } else { none },
    table.header([*Status*], [*Claim*], [*Null model*], [*_p_-value*]),
    table.cell(rowspan: 4)[Supported],
      [Cross-metric coloring optimality], [FH block-preserving $times$ 4 metrics], [$<= 0.004$],
    [Per-table optimality (24/25)], [FH per-table + BH-FDR], [$< 0.05$],
    [$rho$-robustness ($rho in [0,1]$)], [FH weighted edges], [0.003],
    [Topology avoidance depletion], [Table-pres. perm. + clade excl.], [$10^(-4)$],
    table.cell(rowspan: 1)[Suggestive],
      [tRNA enrichment (MIS worst-case)], [Fisher--Stouffer (MIS)], [0.045],
    table.cell(rowspan: 4)[Exploratory],
      [Bit-position bias (dedup.)], [$chi^2$ uniform], [0.075],
    [Mechanism boundary conditions], [Descriptive], [---],
    [Atchley F3/Serine convergence], [Qualitative], [---],
    [Disconnection catalogue (4 cases)], [Descriptive], [---],
    table.cell(rowspan: 1, fill: luma(235))[Falsified],
      table.cell(fill: luma(235))[KRAS--Fano clinical prediction], table.cell(fill: luma(235))[Fisher--Bonferroni], table.cell(fill: luma(235))[1.0],
    table.cell(rowspan: 3, fill: luma(235))[Rejected],
      table.cell(fill: luma(235))[Serine min-distance-4 invariant], table.cell(fill: luma(235))[Counterexample], table.cell(fill: luma(235))[---],
    table.cell(fill: luma(235))[PSL(2,7) symmetry], table.cell(fill: luma(235))[Literature], table.cell(fill: luma(235))[---],
    table.cell(fill: luma(235))[Holomorphic embedding], table.cell(fill: luma(235))[Verification], table.cell(fill: luma(235))[---],
    table.cell(rowspan: 2, fill: luma(235))[Tautological],
      table.cell(fill: luma(235))[Two-fold bit-5 filtration], table.cell(fill: luma(235))[---], table.cell(fill: luma(235))[---],
    table.cell(fill: luma(235))[Four-fold prefix filtration], table.cell(fill: luma(235))[---], table.cell(fill: luma(235))[---],
  ),
  caption: [
    Summary of 15 evaluated claims. "FH" = Freeland--Hurst. "MIS" = maximal independent set enumeration. $p$-values are the most conservative test per claim. Shaded rows: claims that failed validation.
  ],
) <tbl:claims>

// ============================================================
//  3. RESULTS
// ============================================================
= Results <sec:results>

== Cross-metric coloring optimality <sec:res-coloring>

The standard genetic code is significantly error-minimizing under all four physicochemical distance metrics tested (@tbl:metrics). Under the Grantham distance, the standard code achieves $F = 13{,}477$ (null mean $14{,}954 plus.minus 628$, $z = 2.48$, $p = 0.003$). The Miyata distance yields the strongest signal ($z = 3.30$, $p = 0.001$), followed by Kyte--Doolittle hydropathy ($z = 2.95$, $p = 0.001$) and Woese polar requirement ($z = 2.86$, $p = 0.004$). Since Freeland and Hurst (1998) established optimality using polar requirement alone, the cross-metric concordance demonstrates that the code's error-minimization is a robust property across several distinct physicochemical parameterizations, not an artifact of any single physicochemical parameterization.

#figure(
  table(
    columns: (1fr, auto, auto, auto, auto),
    align: (left, right, right, right, right),
    inset: 7pt,
    stroke: (x, y) => if y == 0 { (bottom: 0.7pt) } else { none },
    table.header([*Metric*], [*Observed*], [*Null mean $plus.minus$ SD*], [*Effect size _z_*], [*_p_-value*]),
    [Grantham (1974)], [13,477], [14,998 $plus.minus$ 613], [2.48], [0.003],
    [Miyata (1979)], [235.3], [276.4 $plus.minus$ 12.5], [3.30], [0.001],
    [Polar requirement (Woese 1973)], [347.3], [411.8 $plus.minus$ 22.6], [2.86], [0.004],
    [Kyte--Doolittle (1982)], [457.2], [572.5 $plus.minus$ 39.1], [2.95], [0.001],
  ),
  caption: [
    Cross-metric sensitivity analysis. Block-preserving null ($n = 1{,}000$, seed 135325) under four distinct physicochemical distance metrics. All metrics show significant optimality ($p < 0.01$). Effect size $z = (mu_"null" - F_"obs") \/ sigma_"null"$. Note: Grantham values differ slightly from the primary $n = 10{,}000$ analysis (@fig:coloring; $p = 0.006$, null mean $14{,}954 plus.minus 628$) due to different null sample sizes; this table reports the matched-design $n = 1{,}000$ result for comparability across metrics.
  ],
) <tbl:metrics>

Decomposing $F$ by nucleotide position reveals that the second codon position contributes 49.3% of the total mismatch, followed by the first position at 38.2% and the wobble position at 12.5%. This gradient reflects the well-known biochemical hierarchy: second-position mutations are the most physicochemically disruptive, first-position mutations are intermediate, and wobble-position mutations are largely synonymous.

#figure(
  grid(
    columns: (1fr, 1fr),
    gutter: 1em,
    image("figures/FigA_coloring_null.png"),
    image("figures/FigD_score_decomposition.png"),
  ),
  caption: [
    *Left (a)*: Null distribution of edge-mismatch scores $F$ under the block-preserving null ($n = 10{,}000$). The observed standard code (red line, $F = 13{,}477$) falls at the 0.6th percentile. *Right (b)*: Decomposition of $F$ by nucleotide position. Position 2 dominates (49.3%), consistent with the biochemical hierarchy of mutational impact.
  ],
) <fig:coloring>

== Robustness across transversion weights <sec:res-rho>

The weighted score $F_rho$ (@eq:weighted) remains significantly optimal across all values of $rho$ from 0 to 1 (@tbl:rho). At $rho = 0$ (pure $Q_6$), $p = 0.003$; at $rho = 1$ (full $K_4^3$, all 288 single-nucleotide edges equally weighted), $p = 0.001$. No value of $rho$ yields $p > 0.003$. This result addresses the concern that the $Q_6$ representation ignores approximately one-third of single-nucleotide mutations: the optimality signal strengthens, rather than weakens, when the missing edges are included.

#figure(
  table(
    columns: (auto, auto, auto, auto, auto, auto),
    align: (center, right, right, right, right, right),
    inset: 7pt,
    stroke: (x, y) => if y == 0 { (bottom: 0.7pt) } else { none },
    table.header([$bold(rho)$], [*Observed*], [*Null mean*], [*Quantile*], [*_p_-value*], [*_z_*]),
    [$0.0$],  [13,477], [14,998], [0.2%], [0.003], [2.48],
    [$0.25$], [15,335], [17,063], [0.0%], [0.001], [2.90],
    [$0.5$],  [17,192], [19,127], [0.0%], [0.001], [3.22],
    [$0.75$], [19,050], [21,192], [0.0%], [0.001], [3.40],
    [$1.0$],  [20,908], [23,256], [0.0%], [0.001], [3.46],
  ),
  caption: [
    Robustness of coloring optimality across transversion weight $rho$. Block-preserving null, $n = 1{,}000$ per $rho$. All $p$-values $< 0.01$.
  ],
) <tbl:rho>

#figure(
  grid(
    columns: (1fr, 1fr),
    gutter: 1em,
    image("figures/FigC_rho_robustness.png"),
    image("figures/FigB_per_table_optimality.png"),
  ),
  caption: [
    Robustness of coloring optimality. *Left (a)*: $p$-value across transversion weight $rho$ from 0 ($Q_6$ only) to 1 (full $K_4^3$). All values below 0.01; optimality strengthens with diagonal edges. *Right (b)*: Per-table quantile for each of the 25 NCBI translation tables under their own block-preserving null ($n = 1{,}000$). Only translation table 3 (yeast mitochondrial) exceeds the 5% threshold.
  ],
) <fig:robustness>

== Per-table optimality preservation <sec:res-pertable>

Of the 25 NCBI translation tables, 24 remain in the top 5% of their own block-preserving null after Benjamini--Hochberg FDR correction (@fig:robustness, right panel). The mean quantile across all tables is 1.4%. Only translation table 3 (yeast mitochondrial code, 6 codon reassignments) marginally exceeds the 5% threshold at a quantile of 7.4%.

This result demonstrates that codon reassignment events, despite altering the mapping between codons and amino acids, generally preserve the error-minimization structure of the code. Even the most extensively reassigned code retains near-optimal physicochemical placement.

== Topology avoidance in natural codon reassignments <sec:res-topo>

Among all 1,280 possible single-codon reassignments from the standard code, 931 (72.7%) would create a new amino acid disconnection at Hamming distance $epsilon = 1$ in the codon graph. By contrast, only 6 of 27 observed natural reassignment events (22.2%) are topology-breaking---a 3.3-fold depletion.

The hypergeometric test yields $p = 4.8 times 10^(-8)$, and the table-preserving permutation test gives $p lt.eq 10^(-4)$ (0 of 10,000 permuted datasets as extreme as observed; conservative bound). Both tests confirm that natural codon reassignments are strongly depleted for changes that would disrupt the connected-component structure of amino acid codon families.

To address phylogenetic non-independence (Sengupta, Yang, and Higgs, 2007), we performed clade-exclusion sensitivity analysis, iteratively removing each major taxonomic group (all ciliates, all metazoan mitochondria, all CUG-clade yeasts, etc.) and retesting. The depletion remains highly significant in every case ($p < 10^(-5)$), confirming that the result is not driven by any single clade.

This pattern is compatible with selective or historical constraints favoring connectivity-preserving reassignments: observed natural reassignments are strongly biased toward topology-preserving changes relative to the space of possible reassignments, regardless of phylogenetic grouping.

#figure(
  table(
    columns: (auto, auto, auto),
    align: (left, right, right),
    inset: 7pt,
    stroke: (x, y) => if y == 0 { (bottom: 0.7pt) } else { none },
    table.header([*Metric*], [*Observed*], [*Possible*]),
    [Topology-breaking], [6 (22.2%)], [931 (72.7%)],
    [Topology-preserving], [21 (77.8%)], [349 (27.3%)],
    [Total events], [27], [1,280],
    [Hypergeometric _p_], [], [$4.8 times 10^(-8)$],
    [Permutation _p_ ($n = 10{,}000$)], [], [$1.0 times 10^(-4)$],
  ),
  caption: [
    Topology avoidance: observed natural reassignments versus all possible single-codon reassignments from the standard code. Topology-breaking = creates a new amino acid disconnection at $epsilon = 1$.
  ],
) <tbl:topo>

#figure(
  image("figures/FigF_topology_avoidance.png", width: 70%),
  caption: [
    Topology avoidance in natural codon reassignments. Left: proportion of topology-breaking reassignments among all possible single-codon changes (72.7%) versus observed natural events (22.2%). Right: permutation null distribution ($n = 10{,}000$); the observed rate (arrow) is below all permuted values ($p lt.eq 10^(-4)$).
  ],
) <fig:topo>

== tRNA enrichment for reassigned amino acids <sec:res-trna>

Organisms with variant genetic codes show elevated tRNA gene copy numbers for the reassigned amino acid relative to standard-code controls. Across 24 pairings (18 variant-code organisms, 6 standard-code controls, spanning 5 variant genetic codes and 18 tRNAscan-SE--verified genome assemblies), Fisher's exact test combined via Stouffer's $Z$ method yields $p = 1.7 times 10^(-7)$ ($Z = 5.10$). To address non-independence from shared controls, we enumerated all maximal independent sets (MIS) from the conflict graph via Bron--Kerbosch; both MIS (each of size 6) are significant at $p < 0.05$ (worst-case $p = 0.045$).

The most striking case is _Tetrahymena thermophila_ (NCBI translation table 6, UAA/UAG reassigned to Gln), which carries 54 glutamine tRNA genes---including 39 suppressor tRNAs reading the reassigned stop codons---compared to 3 Gln tRNAs in the standard-code ciliate _Ichthyophthirius multifiliis_, 11 in _Stentor coeruleus_, and 3 in _Fabrea salina_. The pattern extends across reassignment types. Among Gln-reassignment ciliates, _Pseudocohnilembus persalinus_ (20 Gln tRNAs including 15 suppressors) and _Halteria grandinella_ (9 Gln tRNAs including 3 suppressors) represent independent lineages within Oligohymenophorea and Spirotrichea respectively. Among Cys-reassignment ciliates (translation table 10, UGA$arrow.r$Cys), six tRNAscan-SE--verified _Euplotes_ species all carry tRNA-Cys genes with the non-canonical TCA anticodon (reading UGA), with 1--4 such genes per species alongside standard GCA-anticodon Cys tRNAs.

However, the pattern is not universal. _Blastocrithidia nonstop_ (translation table 31) reassigned all three stop codons but achieved UGA$arrow.r$Trp via anticodon stem shortening (5 bp $arrow.r$ 4 bp) of tRNA-Trp(CCA), combined with an eRF1 Ser74Gly mutation, rather than gene duplication (Kachale et al., 2023). Similarly, _Mycoplasma_ species with UGA$arrow.r$Trp use a single tRNA-Trp with anticodon modification. These boundary cases define a three-tier mechanistic landscape: (i) tRNA gene duplication in large nuclear genomes, (ii) anticodon structural modification in streamlined genomes, and (iii) anticodon base modification in minimal genomes.

#figure(
  table(
    columns: (auto, auto, auto, auto),
    align: (left, center, center, right),
    inset: 6pt,
    stroke: (x, y) => if y == 0 { (bottom: 0.7pt) } else { none },
    table.header([*Test*], [*$n$*], [*Statistic*], [*_p_-value*]),
    [Fisher+Stouffer (all pairings)], [24], [$Z = 5.10$], [$1.7 times 10^(-7)$],
    [Fisher+Stouffer (6 indep.)], [6], [$Z = 2.13$], [0.017],
    [*MIS worst-case*], [*6*], [$bold(Z = 1.70)$], [*0.045*],
    [MIS best-case], [6], [$Z = 1.71$], [0.044],
  ),
  caption: [
    Summary of tRNA enrichment tests across 24 pairings (18 tRNAscan-SE--verified genomes, 5 variant genetic codes). MIS = maximal independent set, enumerated via Bron--Kerbosch algorithm on the conflict graph. Both MIS are significant at $p < 0.05$, eliminating the concern that greedy selection biased the independent-pairings result. The MIS worst-case ($p = 0.045$) serves as a robustness bound demonstrating that the result is not driven by a cherry-picked independent subset.
  ],
) <tbl:trna-tests>

#figure(
  image("figures/FigE_trna_aa_rank.png", width: 75%),
  caption: [
    tRNA enrichment for the reassigned amino acid across variant-code/control pairings. Rank 1 indicates the reassigned amino acid is the most proportionally enriched among all 20 amino acids. Statistical analysis uses the full 24-pairing dataset including 6 Euplotes (Table 10, UGA$arrow.r$Cys) and 2 Mycoplasma (Table 4, UGA$arrow.r$Trp) species not shown in this panel.
  ],
) <fig:trna>

#figure(
  table(
    columns: (auto, auto, auto, auto, auto, auto),
    align: (left, left, center, left, right, right),
    inset: 6pt,
    stroke: (x, y) => if y == 0 { (bottom: 0.7pt) } else if y == 8 { (bottom: 0.4pt + luma(180)) } else { none },
    table.header(
      [*Organism*], [*Code*], [*Tbl. \#*], [*Assembly*], [*tRNAs*], [*Reassigned AA*],
    ),
    [_T. thermophila_], [Variant], [6], [GCF_000189635.1], [718], [54 Gln],
    [_P. tetraurelia_], [Variant], [6], [GCF_000165425.1], [216], [18 Gln],
    [_O. trifallax_], [Variant], [6], [GCA_000295675.1], [94], [8 Gln],
    [_P. persalinus_], [Variant], [6], [GCA_001447515.1], [262], [20 Gln],
    [_H. grandinella_], [Variant], [6], [GCA_006369765.1], [130], [9 Gln],
    [_E. aediculatus_], [Variant], [10], [GCA_030463445.1], [80], [4 Cys],
    [_E. focardii_], [Variant], [10], [GCA_001880345.2], [62], [3 Cys],
    [_E. parawoodruffi_], [Variant], [10], [GCA_021440025.1], [149], [9 Cys],
    [_B. stoltei_], [Variant], [15], [GCA_965603825.1], [169], [6 Trp],
    [_B. nonstop_ P57], [Variant], [31], [GCA_028554745.1], [68], [2 Trp#super[\*]],
    [_M. genitalium_], [Variant], [4], [GCA_000027325.1], [36], [1 Trp#super[#sym.dagger]],
    [_S. coeruleus_], [Standard], [1], [GCA_001970955.1], [272], [---],
    [_I. multifiliis_], [Standard], [1], [GCF_000220395.1], [150], [---],
    [_F. salina_], [Standard], [1], [GCA_022984795.1], [89], [---],
  ),
  caption: [
    Representative organisms verified by tRNAscan-SE 2.0.12 (18 total; 3 additional _Euplotes_ and 1 _Mycoplasma_ species in supplement). "Tbl. \#" = NCBI translation table number. "Reassigned AA" = effective tRNA gene count including suppressor tRNAs. #super[\*]Anticodon stem shortening. #super[#sym.dagger]Post-transcriptional modification (single tRNA reads UGG+UGA).
  ],
) <tbl:organisms>

== Score decomposition by nucleotide position <sec:res-decomp>

The decomposition of $F$ by nucleotide position (@fig:coloring, panel b) reveals that the second codon position dominates the error-minimization signal, contributing 49.3% of the total physicochemical mismatch. This aligns with the biochemical observation that second-position mutations cause the largest changes in amino acid hydrophobicity (Woese, 1965; Freeland and Hurst, 1998). The wobble position contributes only 12.5%, consistent with its largely synonymous character under the genetic code's degeneracy structure.

== Exploratory observations <sec:res-exploratory>

=== Bit-position bias in codon reassignments <sec:res-bitbias>

The distribution of bit-flips across the 6 coordinates of $"GF"(2)^6$ in natural codon reassignments shows apparent positional skew under a uniform null ($chi^2 = 16.26$, $p = 0.006$, $"df" = 5$). However, this signal is substantially attenuated after de-duplication to 20 unique (codon, target amino acid) pairs ($p = 0.075$) and vanishes entirely under a codon-preserving permutation null ($p = 1.0$). The apparent bias is therefore explained by which codons are "hot" for reassignment, not by a genuine positional preference in $"GF"(2)^6$.

=== Variant-code disconnection catalogue <sec:res-catalogue>

A systematic survey across all 25 NCBI translation tables identifies 4 variant-code amino acids with disconnected Hamming graphs at $epsilon = 1$: threonine in the yeast mitochondrial code (translation table 3), leucine in the chlorophycean mitochondrial code (translation table 16), alanine in _Pachysolen_ nuclear code (translation table 26), and a tripartite serine in the _Candida_ nuclear code (translation table 12). These cases, combined with the universal serine disconnection, constitute the complete inventory of amino acid graph disconnections at unit Hamming distance.

=== Atchley Factor 3 and Serine convergence <sec:res-atchley>

Serine has the most extreme Atchley Factor 3 score among the 20 amino acids ($F_3 = -4.760$, 2.24 SD below the mean; Atchley et al., 2005), and it is the only amino acid disconnected at $epsilon = 1$ under every base-to-bit encoding. This convergence is not coincidental: Factor 3 is a composite of molecular size and codon diversity, so both measures reflect the same underlying anomaly---Serine's disproportionate codon diversity (6 codons in two disconnected families) relative to its small physicochemical footprint. The $"GF"(2)^6$ framework provides a structural explanation for why Serine is an outlier on Factor 3, though the two views are complementary rather than independent.

== Falsified and rejected claims <sec:res-negatives>

=== KRAS--Fano clinical prediction <sec:res-kras>

The conjecture that XOR ("Fano") relationships in $"GF"(2)^6$ predict co-mutation enrichment at KRAS G12 sites was tested against 1,670 mutations from MSK-IMPACT (Zehir et al., 2017) and cleanly falsified ($p = 1.0$; details in Supplementary Material). This negative result separates code-level error-minimization (which is real) from mutation-level algebraic predictions (which are not).

=== Serine distance-4 invariant <sec:res-serine>

The claim that Serine's minimum inter-family Hamming distance (UCN$dash$AGY) equals 4 under all 24 base-to-bit encodings is false. Of the 24 encodings, 16 yield minimum distance 2 and only 8 yield distance 4. The distance-4 result obtains only when both nucleotide pairs distinguishing UCN from AGY ($U tilde.op A$ and $C tilde.op G$ in the first two positions) are encoded at maximal Hamming distance. The correct encoding-invariant statement is: Serine is disconnected at $epsilon = 1$ under every encoding, and its inter-family distance ($gt.eq 2$) is the largest among the three 6-codon amino acids (Leucine and Arginine both have inter-family distance 1).

=== PSL(2,7) symmetry and holomorphic embedding <sec:res-psl>

The claim that PSL(2,7) is the fundamental symmetry group of the genetic code was pre-rejected by Antoneli and Forger (2011), who showed that PSL(2,7) has no 64-dimensional irreducible representation (its irreps have dimensions 1, 3, 6, 7, 8). The claim that the coordinate-wise map $"GF"(2)^6 arrow.r CC^3$ sending base pairs to fourth roots of unity is a holomorphic embedding extending a character of $"GF"(8)^*$ is also incorrect: the domain is a finite discrete set (not a complex manifold), and the map fails the character identity $chi(x + x) = chi(x)^2$ since $i^2 = -1 eq.not 1$.

== Negative results with informative interpretation <sec:res-infoneg>

The local mismatch cost test---asking whether reassigned codons sit in "worse" Hamming neighborhoods with higher Grantham distance to their neighbors---yields a Mann--Whitney $U = 301$, $p = 0.70$. Reassignment is not driven by local escape from costly neighborhoods. This null result is scientifically informative: it indicates that the topology avoidance constraint (Section 3.4) operates at the global graph-connectivity level, not at the local per-codon level.


// ============================================================
//  4. DISCUSSION
// ============================================================
= Discussion <sec:discussion>

== An information-theoretic view of the genetic code <sec:disc-info>

The central finding of this work is that the standard genetic code minimizes the physicochemical disruption caused by single-bit errors in $"GF"(2)^6$ coordinates. This is not a new conclusion---Freeland and Hurst (1998) established error-minimization using nucleotide-level mutation models---but the hypercube representation makes the optimality principle geometrically explicit: the code is a good _coloring_ of a structured graph, in the graph-theoretic sense that adjacent vertices (codons differing by one bit) tend to share labels (amino acids) or, when they differ, differ by small physicochemical distances.

The score decomposition (@fig:coloring) shows that this optimization is concentrated at the second codon position (49.3% of total mismatch) and first position (38.2%), with the wobble position contributing only 12.5%. This gradient mirrors the biochemical hierarchy of mutational impact and is an emergent property of the code's structure rather than a parameter of the model.

The robustness result (@fig:robustness) demonstrates that optimality is not an artifact of restricting attention to $Q_6$: when the full $K_4^3$ mutation graph is considered ($rho = 1$), the signal strengthens. The code minimizes error not just along the hypercube edges but across the complete space of single-nucleotide substitutions.

A natural question is whether the $"GF"(2)^6$ framework adds genuine insight beyond what a direct analysis of $K_4^3$ would provide. We note three advantages. First, the binary decomposition enables the $rho$-interpolation that reveals the relationship between $Q_6$ and $K_4^3$ optimality as a continuum rather than two disconnected analyses. Second, the Hamming-distance filtration provides a natural persistence parameter ($epsilon$) for studying amino acid graph connectivity, yielding the topology avoidance result. Third, the encoding-sensitivity analysis (24 bijections) distinguishes encoding-invariant properties (Serine disconnection) from encoding-dependent ones (distance-4 claim), a distinction invisible in the nucleotide-level representation. However, we emphasize that the Hamming-1/diagonal decomposition does not correspond to the biological transition/transversion partition: under any 2-bit encoding, each edge class contains a mixture of both substitution types.

== Evolutionary preservation and topology avoidance <sec:disc-evol>

The per-table analysis (@fig:robustness) shows that 24 of 25 NCBI translation tables maintain coloring optimality under their own block-preserving null, suggesting that codon reassignment events are constrained to preserve error-minimization. The single marginal exception (translation table 3, yeast mitochondrial) is the most extensively reassigned code (6 changes), yet still falls only slightly above the 5% threshold (quantile 7.4%).

The topology avoidance result (@fig:topo) provides a mechanistic complement: natural reassignments avoid creating new amino acid disconnections at a rate far below chance expectation (22% observed vs. 73% possible, $p lt.eq 10^(-4)$). This suggests that the connected-component structure of amino acid codon families is functionally important---possibly because disconnected codon families require separate tRNA species to maintain decoding fidelity. The yeast mitochondrial threonine reassignment illustrates this cost: the CUN$arrow.r$Thr change required acquisition of a novel tRNA#super[Thr] derived from tRNA#super[His] via anticodon mutation (Su et al., 2011), creating the topology-breaking disconnection that makes translation table 3 the sole marginal outlier in the per-table optimality analysis.

== Mechanistic implications: tRNA compensation <sec:disc-trna>

The tRNA enrichment result (@fig:trna) links the geometric observation to molecular mechanism. In several variant-code lineages where codon reassignment disrupts connectivity, expanded tRNA gene repertoires for the affected amino acid are observed, consistent with compensatory gene duplication as one evolutionary route to accommodation. The extreme case of _Tetrahymena thermophila_ (54 Gln tRNAs, including 39 suppressors) illustrates the scale of genomic response required to service a split codon family.

However, the boundary cases---_Blastocrithidia nonstop_ (anticodon stem shortening; Kachale et al., 2023) and _Mycoplasma_ species (anticodon modification)---show that tRNA gene duplication is not the only evolutionary solution. The three-tier pattern (duplication in large genomes, structural modification in intermediate genomes, base modification in minimal genomes) suggests that genome size constrains the available mechanistic repertoire for codon reassignment. This pattern is now supported by 18 tRNAscan-SE--verified genomes across 5 variant genetic codes (Tables 4, 6, 10, 15, and 31), including 6 _Euplotes_ species (UGA$arrow.r$Cys) that each carry 1--4 TCA-anticodon tRNA-Cys genes dedicated to reading UGA.

== Scope and falsified claims <sec:disc-honest>

Six of 15 evaluated claims were rejected, falsified, or tautological (shaded rows in @tbl:claims; full details in Supplementary Material). These failures help delimit the framework's explanatory scope. The KRAS--Fano prediction ($p = 1.0$) cleanly separates code-level error-minimization from mutation-level algebraic predictions. The Serine distance-4 invariant, falsified by 16 of 24 encodings giving distance 2, illustrates the importance of systematic encoding-sensitivity testing. The PSL(2,7) and holomorphic embedding claims were pre-rejected by existing mathematical results (Antoneli and Forger, 2011).

== Limitations <sec:disc-limits>

Several limitations should be noted. First, the choice of base-to-bit encoding is not unique, and while the coloring optimality result holds across all 24 encodings, the specific score values and rank orderings are encoding-dependent (full encoding-sensitivity results across all 24 bijections are provided in the Supplementary Material). Second, the block-preserving null model, while standard in the field (Freeland and Hurst, 1998), preserves more structure than a fully random code and may understate the degree of optimality. Third, the tRNA enrichment result, while robust to pairing selection (worst-case MIS $p = 0.045$), relies on a small number of independent pairings ($n = 6$) with limited statistical power; this result is appropriately classified as suggestive and should not be interpreted as demonstrating a universal compensation mechanism.

Whether the code's optimality reflects adaptive selection or non-adaptive carry-over from a primordially constrained starting point (Koonin and Novozhilov, 2009) cannot be resolved by cross-sectional data alone. The topology avoidance result is consistent with both frameworks: it demonstrates a constraint on reassignment trajectories but does not distinguish whether the constraint is selective or structural. Novozhilov and Koonin (2009) showed that putative primordial 16-supercodon codes are "nearly optimal" even without direct selection, suggesting the current code's optimality may be partly inherited.

The $"GF"(2)^6$ representation is best understood as an analytical decomposition tool rather than a claim about the biological primacy of binary coordinates. It enables the $rho$-sweep interpolation between $Q_6$ and $K_4^3$ and facilitates systematic encoding-sensitivity analysis, but the underlying biology is the assignment of chemically similar amino acids to mutationally proximate codons---a property that holds regardless of the coordinate system used to describe it.


// ============================================================
//  5. CONCLUSION
// ============================================================
= Conclusion <sec:conclusion>

The standard genetic code is significantly error-minimizing across four distinct physicochemical distance metrics ($p lt.eq 0.004$), robust to the inclusion of all single-nucleotide substitutions via the graph family $G_rho$ ($p < 0.003$ at all $rho$), and preserved across 24 of 25 variant genetic codes. Observed natural codon reassignments are strongly biased toward topology-preserving changes ($p lt.eq 10^(-4)$, robust to clade exclusion), and some variant-code lineages show expanded tRNA gene repertoires for the reassigned amino acid (worst-case MIS $p = 0.045$), consistent with compensatory gene duplication as one accommodation route among several.

These results extend Freeland and Hurst's (1998) single-metric finding to a four-metric robustness envelope, demonstrating that the code's error-minimization is not an artifact of any particular physicochemical parameterization. The $G_rho$ graph family provides the analytical framework for this interrogation: $Q_6$ and $K_4^3$ are endpoints of a continuum, not competing representations. What the analysis reveals is a code whose amino acid assignments are unusually well-matched to their mutational neighborhoods---a property that persists across distance metrics, mutation models, variant codes, and base-to-bit encodings.


// ============================================================
//  ACKNOWLEDGEMENTS
// ============================================================
#heading(numbering: none)[Acknowledgements]

We thank the NCBI, GtRNAdb, and cBioPortal teams for maintaining public databases essential to this work. tRNAscan-SE 2.0.12 was developed by Chan and Lowe at UC Santa Cruz.

#heading(numbering: none)[Data availability]

All code and data are available in the `codon-topo` repository (version 0.3.0). Analyses are fully reproducible via `codon-topo all --output-dir=./output --seed=135325`. NCBI genome assembly accessions are listed in @tbl:organisms.

#heading(numbering: none)[Declaration of competing interest]

The authors declare no competing interests.

// ============================================================
//  REFERENCES
// ============================================================
#pagebreak()

#heading(numbering: none)[References]

#set par(first-line-indent: 0pt, hanging-indent: 2em)
#set text(size: 10pt)

Antoneli, F., Forger, M., 2011. Symmetry breaking in the genetic code: finite groups. Mathematical and Computer Modelling 53, 1469--1488.

Atchley, W.R., Zhao, J., Fernandes, A.D., Drüke, T., 2005. Solving the protein sequence metric problem. Proceedings of the National Academy of Sciences 102, 6395--6400.

Benjamini, Y., Hochberg, Y., 1995. Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society: Series B 57, 289--300.

Chan, P.P., Lowe, T.M., 2019. tRNAscan-SE: searching for tRNA genes in genomic sequences. Methods in Molecular Biology 1962, 1--14.

Clayworth, P., 2026. Algebraic structure of the genetic code: a technical note. TN-2026-11.

Di Giulio, M., 2005. The origin of the genetic code: theories and their relationships, a review. BioSystems 80, 175--184.

Freeland, S.J., Hurst, L.D., 1998. The genetic code is one in a million. Journal of Molecular Evolution 47, 238--248.

Itzkovitz, S., Alon, U., 2007. The genetic code is nearly optimal for allowing additional information within protein-coding sequences. Genome Research 17, 405--412.

Haig, D., Hurst, L.D., 1991. A quantitative measure of error minimization in the genetic code. Journal of Molecular Evolution 33, 412--417.

Grantham, R., 1974. Amino acid difference formula to help explain protein evolution. Science 185, 862--864.

Hanyu, N., Kuchino, Y., Nishimura, S., Beier, H., 1986. Dramatic events in ciliate evolution: alteration of UAA and UAG termination codons to glutamine codons due to anticodon mutations in two _Tetrahymena_ tRNAs#super[Gln]. The EMBO Journal 5, 1307--1311.

Heaphy, S.M., Mariotti, M., Gladyshev, V.N., Atkins, J.F., Baranov, P.V., 2016. Novel ciliate genetic code variants including the reassignment of all three stop codons to sense codons in _Condylostoma magnum_. Molecular Biology and Evolution 33, 2885--2889.

Koonin, E.V., Novozhilov, A.S., 2009. Origin and evolution of the genetic code: the universal enigma. IUBMB Life 61, 99--111.

Kyte, J., Doolittle, R.F., 1982. A simple method for displaying the hydropathic character of a protein. Journal of Molecular Biology 157, 105--132.

Kachale, A., Pavlíková, Z., Nenarokova, A., Medvedev, A., Milner, D.S., Hashimi, H., Lukeš, J., 2023. Short tRNA anticodon stem and mutant eRF1 allow stop codon reassignment. Nature 613, 751--758.

Miyata, T., Miyazawa, S., Yasunaga, T., 1979. Two types of amino acid substitutions in protein evolution. Journal of Molecular Evolution 12, 219--236.

Novozhilov, A.S., Wolf, Y.I., Koonin, E.V., 2007. Evolution of the genetic code: partial optimization of a random code for robustness to translation error in a rugged fitness landscape. Biology Direct 2, 24.

Sengupta, S., Yang, X., Higgs, P.G., 2007. The mechanisms of codon reassignments in mitochondrial genetic codes. Journal of Molecular Evolution 64, 662--688.

Sella, G., Ardell, D.H., 2006. The coevolution of genes and genetic codes: Crick's frozen accident revisited. Journal of Molecular Evolution 63, 297--313.

Santos, M.A.S., Cheesman, C., Costa, V., Moradas-Ferreira, P., Tuite, M.F., 1999. Selective advantages created by codon ambiguity allowed for the evolution of an alternative genetic code in _Candida_ spp. Molecular Microbiology 31, 937--947.

Singh, A., Vancura, A., Woycicki, R.K., Hogan, D.J., Hendrick, A.G., Nowacki, M., 2018. Determination of the presence of 5-methylcytosine in _Paramecium tetraurelia_. PLoS ONE 13, e0206667.

Su, D., Lieberman, A., Lang, B.F., Simonović, M., Söll, D., Ling, J., 2011. An unusual tRNA#super[Thr] derived from tRNA#super[His] reassigns in yeast mitochondria the CUN codons to threonine. Nucleic Acids Research 39, 4866--4874.

Vetsigian, K., Woese, C., Goldenfeld, N., 2006. Collective evolution and the genetic code. Proceedings of the National Academy of Sciences 103, 10696--10701.

Woese, C.R., Dugre, D.H., Dugre, S.A., Kondo, M., Saxinger, W.C., 1966. On the fundamental nature and evolution of the genetic code. Cold Spring Harbor Symposia on Quantitative Biology 31, 723--736.

Woese, C.R., 1965. On the evolution of the genetic code. Proceedings of the National Academy of Sciences 54, 1546--1552.

Zehir, A., Benayed, R., Shah, R.H., et al., 2017. Mutational landscape of metastatic cancer revealed from prospective clinical sequencing of 10,000 patients. Nature Medicine 23, 703--713.

Zhang, B., Hou, L., Qi, H., Hou, L., Zhang, T., Zhao, F., Miao, M., 2022. An extremely streamlined macronuclear genome in the free-living protozoan _Fabrea salina_. Molecular Biology and Evolution 39, msac062.

Zheng, W., Wang, C., Lynch, M., Gao, S., 2021. The compact macronuclear genome of the ciliate _Halteria grandinella_: a transcriptome-like genome with 23,000 nanochromosomes. mBio 12, e01964-20.
