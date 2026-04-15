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
//  PIPELINE DATA — all inline statistics are read from here
//  Regenerate with: codon-topo all --output-dir=./output
// ============================================================
#let stats = json("manuscript_stats.json")
#let mm = stats.metrics       // per-metric results
#let tq6 = stats.topology_avoidance_q6
#let tk43 = stats.topology_avoidance_k43
#let pt = stats.per_table
#let rho_data = stats.rho_sweep
#let dec = stats.decomposition
#let cl = stats.condlogit
#let trna = stats.trna

// Helper: format a number with commas for thousands
#let fmtk(n) = {
  let v = int(calc.round(n, digits: 0))
  let s = str(v)
  let chars = s.clusters()
  let out = ()
  for (i, c) in chars.rev().enumerate() {
    if i > 0 and calc.rem(i, 3) == 0 { out.push(",") }
    out.push(c)
  }
  out.rev().join()
}

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
  The standard genetic code has long been recognized as reducing the impact of point mutations, but the robustness of this property across physicochemical metrics and across naturally variant codes has not been systematically quantified. We embed the 64 codons in $"GF"(2)^6$, yielding a 6-dimensional hypercube subgraph of the complete single-nucleotide mutation graph $K_4^3$, and introduce a weight parameter $rho in [0,1]$ to interpolate between these mutation neighborhoods. Under a block-preserving null ($n = 10{,}000$ permutations), the standard code is significantly low-cost under four distinct amino acid distance metrics---Grantham ($p = 0.006$), Miyata ($p < 0.001$), polar requirement ($p = 0.003$), and Kyte--Doolittle hydropathy ($p = 0.001$)---and this signal strengthens as all single-nucleotide substitutions are included ($rho arrow.r 1$). Across NCBI translation tables, 24 of 25 variant codes retain comparable optimality after Benjamini--Hochberg correction (mean quantile 1.4%).

  Natural codon reassignment events show $approx$3-fold depletion of topology-breaking moves---changes that fragment an amino acid's codon family---relative to the landscape of possible reassignments (22% observed vs 66--73% possible; RR $= 0.31$--$0.34$; permutation $p lt.eq 10^(-4)$), robust to adjacency definition ($Q_6$ vs $K_4^3$) and phylogenetic clade exclusion. Discrete-choice (conditional logit) modeling shows that topology avoidance is not explained by physicochemical cost alone: each feature class adds substantial independent explanatory power ($Delta"AICc" gt.eq 67$; Spearman $r_s = 0.15$ between predictors). Several variant-code lineages with topology-breaking reassignments show elevated tRNA gene counts for the reassigned amino acid (worst-case Stouffer $p = 0.045$, 18 tRNAscan-SE--verified genomes), consistent with compensatory routes to reassignment accommodation. A retrospective reanalysis of eight genome recoding datasets ($>$248,000 codon-level observations) indicates that codon-family boundary crossing does not predict transcriptomic perturbation, whereas local neighborhood mismatch aligns with established recoding-burden covariates.

  #v(0.5em)
  *Keywords:* genetic code evolution; error minimization; graph-theoretic analysis; GF(2)#super[6]; codon reassignment; physicochemical distance; tRNA gene duplication; conditional logit; genome recoding
]

#v(1em)
#set par(first-line-indent: 1.5em)

// ============================================================
//  1. INTRODUCTION
// ============================================================
= Introduction <sec:intro>

The standard genetic code maps 61 sense codons to 20 amino acids through a pattern that has long been recognized as non-random. #cite(<woese1965>, form: "prose") first noted that similar codons tend to encode amino acids with similar physicochemical properties, and #cite(<freeland1998>, form: "prose") demonstrated quantitatively that the standard code sits in approximately the top $10^(-6)$ percentile of random codes for mutational error minimization under the Woese polar requirement distance. These findings established that the code's structure reduces the fitness impact of point mutations, but the question of whether this optimality is specific to particular physicochemical metrics or represents a robust property across several distinct physicochemical parameterizations has not been systematically addressed.

Every codon comprises three nucleotides drawn from ${C, U, A, G}$. By choosing a bijection $phi: {C, U, A, G} -> "GF"(2)^2$---for instance, $C |-> (0,0)$, $U |-> (0,1)$, $A |-> (1,0)$, $G |-> (1,1)$---each codon maps to a vertex of the 6-dimensional binary hypercube $Q_6 = "GF"(2)^6$. The genetic code then becomes a _coloring_ of $Q_6$ by 21 labels (20 amino acids plus the stop signal), and single-nucleotide mutations correspond to edges of $Q_6$ or, more precisely, to a subgraph of the complete mutation graph $K_4^3$.

This representation is not new in principle---binary encodings of the genetic code appear in the mathematical biology literature (e.g., @antoneli2011). Its value lies not in the encoding itself but in the analytical decomposition it enables: the complete single-nucleotide mutation graph $K_4^3$ (288 edges) separates into 192 Hamming-distance-1 edges (single-bit changes) and 96 within-nucleotide distance-2 edges (both bits of one nucleotide position flip simultaneously), allowing systematic interpolation via a weight parameter $rho in [0,1]$. Note that this decomposition does not correspond to the biological transition/transversion partition: under 16 of the 24 possible 2-bit encodings, the Hamming-1 edges contain an equal mixture of transitions and transversions per nucleotide position (2 of each among 4 Hamming-1 pairs), while the remaining 8 encodings place both transitions on diagonal (Hamming-2) edges. The $rho$ parameter should therefore be interpreted as a diagonal-edge inclusion weight, not a transition/transversion weight. Previous work has not exploited this decomposition to test error-minimization across multiple physicochemical metrics, nor extended the analysis to variant genetic codes. We address three questions:

+ *Is the standard code optimal under multiple physicochemical metrics?* We test whether the standard code's edge-mismatch score is extreme relative to block-preserving null models across four distinct physicochemical distance metrics (Grantham, Miyata, Woese polar requirement, and Kyte--Doolittle hydropathy), extending #cite(<freeland1998>, form: "prose") from a single metric to a cross-metric robustness envelope.

+ *Is this structure preserved by evolution?* We ask whether variant genetic codes (NCBI translation tables 2--33) maintain error-minimization, and whether natural codon reassignment events preferentially avoid disrupting the topological connectivity of amino acid codon families.

+ *What are the genomic correlates of disruption?* We test whether organisms whose variant codes break the connectivity of an amino acid's codon graph show elevated tRNA gene copy numbers for the affected amino acid, using tRNAscan-SE--verified data from 10 genomes across 4 eukaryotic supergroups.

In pursuing these questions, we also evaluated several deeper algebraic conjectures from a companion technical note #cite(<clayworth2026>)---including claims about PSL(2,7) symmetry and holomorphic embeddings---and report their falsification as boundary conditions on the framework's scope (Supplementary Material).

The paper is organized as follows. Section 2 describes the encoding formalism, graph decomposition, null models, and statistical methods, including a new retrospective cross-study reanalysis of eight published genome recoding datasets. Section 3 presents the four supported findings (cross-metric coloring optimality, per-table preservation, $rho$-robustness, and topology avoidance), the suggestive tRNA enrichment result, the cross-study reanalysis of synthetic recoding outcomes, and exploratory observations. Section 4 discusses the graph-theoretic interpretation, the relationship to frozen-accident versus adaptive hypotheses @koonin2009, and a three-layer decomposition of how codon-space structure relates to recoding outcomes. All analyses are reproducible via the open-source `codon-topo` pipeline (version #stats._version, seed 135325).


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

where the sum ranges over all 192 edges of $Q_6$ (pairs of vertices at Hamming distance 1), and $Delta$ is a physicochemical distance between amino acids. We test four distinct distance metrics: the #cite(<grantham1974>, form: "prose") composite distance (composition, polarity, volume; range 5--215), the #cite(<miyata1979>, form: "prose") normalized Euclidean distance (polarity and volume only; range 0.06--5.13), the Woese polar requirement absolute difference (range 0--8.2; @woese1966 @woese1973; used by #cite(<freeland1998>, form: "prose")), and the #cite(<kyte1982>, form: "prose") hydropathy absolute difference (range 0--9.0; used by #cite(<haig1991>, form: "prose")). Synonymous edges contribute 0; edges involving a stop codon receive a fixed penalty scaled proportionally to each metric's maximum. A lower $F$ indicates a more error-minimizing code.

Stop codons are held fixed across all null models (Section 2.3), so the stop-codon contribution is a constant offset. Sensitivity analysis across penalty values (0, 150, 215, 300) confirms this is immaterial to the ranking.

== Null models <sec:nullmodels>

=== Block-preserving null (Freeland--Hurst) <sec:null-fh>

Following #cite(<freeland1998>, form: "prose"), we group the 64 codons into 16 blocks of 4, defined by shared first-two-base prefix. Each block's internal pattern of amino acid assignments is preserved, but the mapping of patterns to blocks is permuted uniformly at random. Blocks containing stop codons are held fixed. This null preserves the synonymous codon contiguity (wobble degeneracy) inherent in the genetic code, providing a stringent test: the observed code must beat random codes that share its degeneracy architecture.

For the standard code, we drew $n = 10,000$ null samples (seed 135325) and computed the conservative $p$-value as $(k + 1)/(n + 1)$, where $k$ is the number of null scores below the observed $F$.

=== Per-table null <sec:null-pertable>

Each of the 25 NCBI translation tables was tested against its own block-preserving null ($n = 1{,}000$ per table, seeds = base seed + table ID). $P$-values were corrected for multiple comparisons using the #cite(<benjamini1995>, form: "prose") false discovery rate (FDR) procedure.

=== Graph family $G_rho$ and rho sweep <sec:null-rho>

We define a family of mutation graphs $G_rho$ parameterized by $rho in [0,1]$, where $G_0 = Q_6$ (192 Hamming-1 edges) and $G_1 = K_4^3$ (all 288 single-nucleotide substitution edges). The weighted mismatch score on $G_rho$ is

$ F_rho(c) = F_("Hamming-1")(c) + rho dot.c F_("diagonal")(c) $ <eq:weighted>

where $F_("diagonal")$ sums over the 96 within-nucleotide distance-2 edges. This formulation treats $Q_6$ not as the biological object but as one endpoint of a continuous interpolation toward the full mutation graph. We evaluated 5 values of $rho$ ($n = 1{,}000$ null samples per $rho$, common-seed design).

=== Table-preserving permutation null (topology avoidance) <sec:null-topo>

To test whether natural codon reassignment events avoid disrupting amino acid connectivity, we enumerated all 1,280 possible single-codon reassignments from the standard code. For each amino acid $a$, let $G_a^1$ denote the subgraph induced on codons assigned to $a$ with edges between codons at binary Hamming distance $lt.eq 1$. A reassignment is _topology-breaking_ if it increases the number of connected components of $G_a^1$ for any amino acid $a$. We define a de-duplicated reassignment event as a unique (codon, target amino acid) pair observed across the 25 NCBI translation tables, yielding 27 events from 61 total table-specific changes. We then compared the observed rate of topology-breaking among 27 de-duplicated natural reassignment events to the null expectation.

Two tests were applied: (i) a hypergeometric test treating the observed events as a sample from the finite landscape of possible reassignments ($N = 1{,}280$, $K = 931$ topology-breaking under $Q_6$, $n = 27$ observed, $x = 6$ topology-breaking observed); and (ii) a table-preserving permutation test ($n = 10{,}000$, seed 135325) that permutes target amino acids among a table's reassigned codons, preserving within-table codon and target structure to address phylogenetic non-independence. Risk ratios with 95% log-normal confidence intervals were computed as $"RR" = (x\/n) \/ (K\/N)$. To test robustness to the adjacency definition, both tests were repeated under $K_4^3$ adjacency ($K = 846$ topology-breaking), where codon-family connectivity is defined by nucleotide-level single-substitution adjacency (encoding-independent).

=== Conditional logit model of reassignment choice <sec:null-condlogit>

To test whether topology avoidance contributes independent explanatory power beyond physicochemical optimization, we fit event-level conditional logit (discrete-choice) models to natural reassignment events. Each observed reassignment is treated as a choice from the set of all $approx 1{,}280$ possible single-codon reassignments available at the current code state. For a table with $k$ reassignment events, the code state evolves sequentially: at each step, the model assigns a probability to each candidate move based on a linear utility score, and the observed move is compared against all alternatives.

For each candidate move $m$ from code state $C$, we computed three features: (i) $Delta_"phys"$, the change in local Grantham mismatch cost summed over Hamming-1 edges incident to the reassigned codon; (ii) $Delta_"topo"$, the total increase in connected components across all amino acid codon graphs at $epsilon = 1$; and (iii) $Delta_"tRNA"$, the Hamming distance from the reassigned codon to the nearest codon already encoding the target amino acid (a heuristic proxy for tRNA repertoire disruption). The conditional logit probability of observing move $m^*$ is

$ P(m^* | cal(N)(C)) = frac(exp(bold(w)^top bold(x)_(m^*)), sum_(m in cal(N)(C)) exp(bold(w)^top bold(x)_m)) $ <eq:condlogit>

where $bold(w)$ is the weight vector and $bold(x)_m$ is the feature vector for move $m$. Features were $z$-scored across all candidates for numerical stability.

Since the temporal ordering of reassignment events within a table is unknown, we marginalized the likelihood over all $k!$ orderings for tables with $k > 1$ changes:

$ L_"table" = frac(1, k!) sum_(sigma in S_k) product_(s=1)^(k) P(m_(sigma(s))^* | cal(N)(C_(sigma, s))) $ <eq:orderavg>

For the largest table (yeast mitochondrial, $k = 6$), we averaged over 24 orderings (a random subset of the 720 possible). The total likelihood is the product across all 23 tables with reassignment events (61 total event-steps).

Four nested models were compared: M1 (physicochemistry only, $w_"topo" = w_"tRNA" = 0$), M2 (topology only, $w_"phys" = w_"tRNA" = 0$), M3 (physicochemistry + topology, $w_"tRNA" = 0$), and M4 (all three features). Weights were estimated by maximum likelihood (scipy.optimize, Nelder--Mead with L-BFGS-B refinement). Model comparison used the corrected Akaike Information Criterion (AICc) and likelihood-ratio tests for nested pairs.

=== Fisher--Stouffer test for tRNA enrichment <sec:null-trna>

For each variant-code organism paired with a phylogenetically proximate standard-code control, we built a $2 times 2$ contingency table comparing the proportion of tRNA genes encoding the reassigned amino acid versus all other amino acids. The sampling model treats tRNA gene counts as draws from a hypergeometric distribution conditional on the row and column marginals (total tRNAs per organism, total tRNAs for the focal amino acid), which is appropriate when the question is whether a specific amino acid's share of the tRNA repertoire differs between two genomes. Fisher's exact test (one-sided, alternative = "greater") was applied per pairing, and $p$-values were combined via Stouffer's $Z$ method. To address non-independence from shared control organisms, we constructed a conflict graph (edges connect pairings sharing an organism) and enumerated all maximal independent sets (MIS) via the Bron--Kerbosch algorithm with pivoting. For each MIS with $>= 2$ members, we computed the Stouffer combined $p$-value and report the worst-case (maximum $p$) across all MIS as the primary conservative test.

tRNA gene counts for 17 organisms were verified by running tRNAscan-SE 2.0.12 @chan2019 with Infernal 1.1.4 on NCBI genome assemblies (eukaryotic mode for ciliates, bacterial mode for _Mycoplasma_). Counts for 3 additional organisms were drawn from published literature. The dataset comprises 24 pairings across 5 variant genetic codes (NCBI translation tables 4, 6, 10, 15, and 31) and 3 standard-code controls.

== Retrospective cross-study reanalysis of synthetic recoding outcomes <sec:methods-codonsafe>

To test whether GF(2)#super[6] topology predicts codon recoding outcomes in synthetic biology experiments, we performed a retrospective cross-study reanalysis across eight published genome recoding datasets. For each codon substitution reported in these studies, we computed three topology features under the default encoding ($C |-> (0,0)$, $U |-> (0,1)$, $A |-> (1,0)$, $G |-> (1,1)$):

+ *Boundary crossing* ($epsilon = 1$): whether the source and target codons lie in different connected components of the amino acid's codon graph at Hamming distance 1. Generalized beyond serine to all amino acids with potential disconnections (including leucine, arginine).

+ *Local mismatch change* ($Delta F_"local"$): the difference in total Grantham distance summed over all Hamming-1 neighbors between the target and source codon positions: $Delta F_"local" = sum_(n in N(t)) Delta(c(t), c(n)) - sum_(n in N(s)) Delta(c(s), c(n))$, where $N(v)$ denotes the 6 Hamming-1 neighbors of vertex $v$ in $Q_6$.

+ *Hamming distance*: the number of bit positions differing between source and target vectors in GF(2)#super[6].

Codon positions were extracted from published GenBank genome files by CDS-level comparison using Biopython @cock2009, with codon frame offsets (`/codon_start`) handled explicitly. CDS features were matched between parent and recoded genomes by locus tag. The primary analysis used the Syn57 dataset @robertson2025syn57, which provides a within-study contrast: 37,146 serine recodings cross the UCN$arrow.l.r$AGY family boundary, while 22,859 alanine recodings remain within the GCN family. Syn57 RNA-seq differential expression data (Data file S10) were used as the outcome measure. Design-to-final genome deviations were identified by comparing the vS33A7 design genome (Data file S2) against the verified Syn57 genome (Data file S8), filtering genes with $>$10 CDS-level differences as structural rearrangements (sensitivity analysis across thresholds 3--50 reported in Supplementary Table SX). The #cite(<ostrov2016>, form: "prose") segment viability data (Table S4) were analyzed as a case-control comparing segments with and without lethal design exceptions; Bonferroni correction was applied for three simultaneous tests.

All cross-study reanalysis code is available in the `codon_topo.analysis.codonsafe` subpackage (version #stats._version) and can be run via `codon-topo codonsafe`. Raw data files and download provenance are documented in `data/codonsafe/DATA_MANIFEST.md`.

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
      [Cross-metric coloring optimality], [FH block-preserving $times$ 4 metrics], [$<= 0.006$],
    [Per-table optimality (24/25)], [FH per-table + BH-FDR], [$< 0.05$],
    [$rho$-robustness ($rho in [0,1]$)], [FH weighted edges], [0.006],
    [Topology avoidance ($Q_6$ + $K_4^3$)], [Table-pres. perm. + clade excl.], [$10^(-4)$],
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

// --- dynamic metric data ---
#let g = mm.at("grantham", default: (:))
#let mi = mm.at("miyata", default: (:))
#let pr = mm.at("polar_requirement", default: (:))
#let kd = mm.at("kyte_doolittle", default: (:))

Under the Freeland--Hurst block-preserving null with $n = #fmtk(stats.coloring.n_samples)$ permutations, the standard genetic code is significantly low-cost across all four physicochemical distance metrics examined (@tbl:metrics). For Grantham distance, the observed score is $F = #fmtk(g.observed)$ versus a null mean of $#fmtk(g.null_mean) plus.minus #str(int(calc.round(g.null_std, digits: 0)))$ ($z = #str(calc.round(g.z, digits: 2))$, quantile #str(calc.round(g.quantile, digits: 2))%, $p = #str(calc.round(g.p, digits: 3))$). The same pattern holds for Miyata distance ($F = #str(calc.round(mi.observed, digits: 1))$ vs $#str(calc.round(mi.null_mean, digits: 1)) plus.minus #str(calc.round(mi.null_std, digits: 1))$, $z = #str(calc.round(mi.z, digits: 2))$, $p < 0.001$), Woese polar requirement ($F = #str(calc.round(pr.observed, digits: 1))$ vs $#str(calc.round(pr.null_mean, digits: 1)) plus.minus #str(calc.round(pr.null_std, digits: 1))$, $z = #str(calc.round(pr.z, digits: 2))$, $p = #str(calc.round(pr.p, digits: 3))$), and Kyte--Doolittle hydropathy ($F = #str(calc.round(kd.observed, digits: 1))$ vs $#str(calc.round(kd.null_mean, digits: 1)) plus.minus #str(calc.round(kd.null_std, digits: 1))$, $z = #str(calc.round(kd.z, digits: 2))$, $p = #str(calc.round(kd.p, digits: 3))$). Relative to null expectations, the observed code lowers the mismatch score by #str(calc.round(g.improvement_pct, digits: 1))%, #str(calc.round(mi.improvement_pct, digits: 1))%, #str(calc.round(pr.improvement_pct, digits: 1))%, and #str(calc.round(kd.improvement_pct, digits: 1))%, respectively. Beta-posterior credible intervals for the empirical $p$-values remain entirely below 0.01 for all four metrics (Supplementary Table S1). Since #cite(<freeland1998>, form: "prose") established optimality using polar requirement alone, the cross-metric concordance demonstrates that error-minimization is robust across multiple distinct physicochemical parameterizations, not an artifact of any single distance metric. Under a weaker degeneracy-only null that preserves only codon family sizes without maintaining block contiguity, the signal strengthens substantially ($z > 9$, $p < 10^(-4)$ for all metrics), confirming that the block-preserving null provides a stringent test.

#figure(
  table(
    columns: (1fr, auto, auto, auto, auto, auto),
    align: (left, right, right, right, right, right),
    inset: 6pt,
    stroke: (x, y) => if y == 0 { (bottom: 0.7pt) } else { none },
    table.header([*Metric*], [*Observed*], [*Null mean $plus.minus$ SD*], [*_z_*], [*_p_*], [*Improv.*]),
    [Grantham~#cite(<grantham1974>)],
      [#fmtk(g.observed)],
      [#fmtk(g.null_mean) $plus.minus$ #str(int(calc.round(g.null_std, digits: 0)))],
      [#str(calc.round(g.z, digits: 2))],
      [#str(calc.round(g.p, digits: 3))],
      [#str(calc.round(g.improvement_pct, digits: 1))%],
    [Miyata~#cite(<miyata1979>)],
      [#str(calc.round(mi.observed, digits: 1))],
      [#str(calc.round(mi.null_mean, digits: 1)) $plus.minus$ #str(calc.round(mi.null_std, digits: 1))],
      [#str(calc.round(mi.z, digits: 2))],
      [#if mi.p < 0.001 [< 0.001] else [#str(calc.round(mi.p, digits: 3))]],
      [#str(calc.round(mi.improvement_pct, digits: 1))%],
    [Polar requirement~#cite(<woese1973>)],
      [#str(calc.round(pr.observed, digits: 1))],
      [#str(calc.round(pr.null_mean, digits: 1)) $plus.minus$ #str(calc.round(pr.null_std, digits: 1))],
      [#str(calc.round(pr.z, digits: 2))],
      [#str(calc.round(pr.p, digits: 3))],
      [#str(calc.round(pr.improvement_pct, digits: 1))%],
    [Kyte--Doolittle~#cite(<kyte1982>)],
      [#str(calc.round(kd.observed, digits: 1))],
      [#str(calc.round(kd.null_mean, digits: 1)) $plus.minus$ #str(calc.round(kd.null_std, digits: 1))],
      [#str(calc.round(kd.z, digits: 2))],
      [#str(calc.round(kd.p, digits: 3))],
      [#str(calc.round(kd.improvement_pct, digits: 1))%],
  ),
  caption: [
    Cross-metric sensitivity analysis. Block-preserving null ($n = #fmtk(stats.coloring.n_samples)$, seed 135325) under four distinct physicochemical distance metrics. All metrics show significant optimality ($p < 0.01$). Effect size $z = (mu_"null" - F_"obs") \/ sigma_"null"$. "Improv." = percent improvement of observed score over null mean. Under a weaker degeneracy-only null (preserving only codon family sizes, not block structure), all four metrics yield $z > 9$ and $p < 10^(-4)$, confirming the result is not an artifact of the Freeland--Hurst @freeland1998 block-preserving design.
  ],
) <tbl:metrics>

// --- dynamic decomposition fractions ---
#let pos2_pct = calc.round(dec.position_fractions.position_2 * 100, digits: 1)
#let pos1_pct = calc.round(dec.position_fractions.position_1 * 100, digits: 1)
#let pos3_pct = calc.round(dec.position_fractions.position_3 * 100, digits: 1)

Decomposing $F$ by nucleotide position reveals that the second codon position contributes #str(pos2_pct)% of the total mismatch, followed by the first position at #str(pos1_pct)% and the wobble position at #str(pos3_pct)%. This gradient reflects the well-known biochemical hierarchy: second-position mutations are the most physicochemically disruptive, first-position mutations are intermediate, and wobble-position mutations are largely synonymous.

// Figure 1 — Core topology (persistent homology + disconnection catalogue)
#figure(
  image("figures/Fig1_core_topology.png"),
  caption: [
    Core topology of the genetic code in $"GF"(2)^6$. *(A)* Persistent homology: connected components ($beta_0$) of amino acid codon graphs as a function of Hamming distance threshold $epsilon$. Serine (bold) is the only amino acid disconnected at $epsilon = 1$ across all 25 NCBI translation tables, reconnecting at $epsilon = 4$. *(B)* Disconnection catalogue across all translation tables. Each tile denotes an amino acid disconnected at $epsilon = 1$; the number is the reconnection $epsilon$. Serine is universally disconnected; Leu, Thr, and Ala appear as variant-code disconnections.
  ],
) <fig:topology>

// Figure 2 — Coloring optimality (null + per-table + rho)
#figure(
  image("figures/Fig2_coloring_optimality.png"),
  caption: [
    Coloring optimality of the genetic code. *(A)* Null distribution of Grantham edge-mismatch scores $F$ under the Freeland--Hurst block-preserving null ($n = #fmtk(stats.coloring.n_samples)$). The observed standard code (red line, $F = #fmtk(g.observed)$) falls at the #str(calc.round(g.quantile, digits: 2))th percentile ($p = #str(calc.round(g.p, digits: 3))$). *(B)* $p$-value across transversion weight $rho$ from 0 ($Q_6$) to 1 (full $K_4^3$); all values below 0.01, with optimality strengthening monotonically. *(C)* Per-table quantile for each of 25 NCBI translation tables under their own block-preserving null ($n = 1{,}000$); only table 3 (yeast mito.) exceeds the 5% threshold.
  ],
) <fig:coloring>

== Robustness across transversion weights <sec:res-rho>

// --- dynamic rho sweep ---
#let rho0 = rho_data.per_rho.at(0)
#let rho4 = rho_data.per_rho.at(rho_data.per_rho.len() - 1)

The weighted score $F_rho$ (@eq:weighted) remains significantly optimal across all values of $rho$ from 0 to 1 (@tbl:rho). At $rho = 0$ (pure $Q_6$), $p = #str(calc.round(rho0.p_value, digits: 3))$ ($z = #str(calc.round(rho0.at("z", default: rho0.at("effect_size_z", default: 0)), digits: 2))$); at $rho = 1$ (full $K_4^3$, all 288 single-nucleotide edges equally weighted), $p < 0.001$ ($z = #str(calc.round(rho4.at("z", default: rho4.at("effect_size_z", default: 0)), digits: 2))$). No value of $rho$ yields $p > #str(calc.round(rho0.p_value, digits: 3))$. The optimality signal strengthens monotonically as diagonal edges are included, addressing the concern that the $Q_6$ representation ignores approximately one-third of single-nucleotide mutations.

#figure(
  table(
    columns: (auto, auto, auto, auto, auto, auto),
    align: (center, right, right, right, right, right),
    inset: 7pt,
    stroke: (x, y) => if y == 0 { (bottom: 0.7pt) } else { none },
    table.header([*ρ*], [*Observed*], [*Null mean*], [*Quantile*], [*_p_-value*], [*_z_*]),
    ..rho_data.per_rho.map(r => (
      [$#str(r.rho)$],
      [#fmtk(r.at("observed", default: r.at("observed_score", default: 0)))],
      [#fmtk(r.at("null_mean", default: 0))],
      [#str(calc.round(r.quantile, digits: 2))%],
      [#if r.p_value < 0.001 [< 0.001] else [#str(calc.round(r.p_value, digits: 3))]],
      [#str(calc.round(r.at("z", default: r.at("effect_size_z", default: 0)), digits: 2))],
    )).flatten(),
  ),
  caption: [
    Robustness of coloring optimality across diagonal-edge weight $rho$. Block-preserving null, $n = #fmtk(stats.coloring.n_samples)$ per $rho$. All $p$-values $< 0.01$. Quantiles reported as empirical rank; values $lt.eq 0.02%$ indicate $lt.eq 2$ null samples below the observed score.
  ],
) <tbl:rho>

// Robustness panels are now in Fig 2 (panels B and C above)

== Per-table optimality preservation <sec:res-pertable>

Of the #pt.n_tables NCBI translation tables, #pt.n_significant_bh remain in the top 5% of their own block-preserving null after Benjamini--Hochberg FDR correction (@fig:coloring, panel C). The mean quantile across all tables is #str(calc.round(pt.mean_quantile, digits: 1))%. Only translation table 3 (yeast mitochondrial code, 6 codon reassignments) marginally exceeds the 5% threshold.

This result demonstrates that codon reassignment events, despite altering the mapping between codons and amino acids, generally preserve the error-minimization structure of the code. Even the most extensively reassigned code retains near-optimal physicochemical placement.

== Topology avoidance in natural codon reassignments <sec:res-topo>

Among all #fmtk(tq6.possible_total) possible single-codon reassignments from the standard code, #fmtk(tq6.possible_breaks) (#str(calc.round(tq6.rate_possible * 100, digits: 1))%) would create a new amino acid disconnection at Hamming distance $epsilon = 1$ in $Q_6$. By contrast, only #tq6.observed_breaks of #tq6.observed_total observed natural reassignment events (#str(calc.round(tq6.rate_observed * 100, digits: 1))%) are topology-breaking---a #str(calc.round(tq6.rate_possible / calc.max(tq6.rate_observed, 0.001), digits: 1))-fold depletion (hypergeometric $p = #str(tq6.hypergeom_p)$; table-preserving permutation $p lt.eq 10^(-4)$).

This depletion is robust to the choice of adjacency graph. Under the fuller $K_4^3$ single-nucleotide mutation graph, where codon-family connectivity is defined by nucleotide-level adjacency (encoding-independent), #fmtk(tk43.possible_breaks) of #fmtk(tk43.possible_total) possible reassignments (#str(calc.round(tk43.rate_possible * 100, digits: 1))%) are topology-breaking---somewhat fewer than under $Q_6$, because some $Q_6$-disconnected pairs become connected when within-nucleotide Hamming-2 edges are admitted. The observed proportion remains #str(calc.round(tk43.rate_observed * 100, digits: 1))% (#tk43.observed_breaks/#tk43.observed_total), yielding a #str(calc.round(tk43.depletion_fold, digits: 1))-fold depletion ($"RR" = #str(calc.round(tk43.risk_ratio, digits: 2))$, 95% CI $[#str(calc.round(tk43.risk_ratio_ci_95.at(0), digits: 2)), #str(calc.round(tk43.risk_ratio_ci_95.at(1), digits: 2))]$; hypergeometric $p = #str(tk43.hypergeom_p)$; permutation $p lt.eq 10^(-4)$). The avoidance of topology-disrupting reassignment trajectories is therefore not specific to the hypercube subgraph.

To address phylogenetic non-independence @sengupta2007, we performed clade-exclusion sensitivity analysis, iteratively removing each major taxonomic group (all ciliates, all metazoan mitochondria, all CUG-clade yeasts, etc.) and retesting. The depletion remains highly significant in every case ($p < 10^(-5)$), confirming that the result is not driven by any single clade.

#figure(
  table(
    columns: (auto, auto, auto, auto, auto),
    align: (left, right, right, right, right),
    inset: 6pt,
    stroke: (x, y) => if y == 0 { (bottom: 0.7pt) } else { none },
    table.header([*Metric*], [*Observed*], [*Possible ($Q_6$)*], [*Possible ($K_4^3$)*], [*Expected (chance)*]),
    [Topology-breaking],
      [#tq6.observed_breaks (#str(calc.round(tq6.rate_observed * 100, digits: 1))%)],
      [#fmtk(tq6.possible_breaks) (#str(calc.round(tq6.rate_possible * 100, digits: 1))%)],
      [#fmtk(tk43.possible_breaks) (#str(calc.round(tk43.rate_possible * 100, digits: 1))%)],
      [---],
    [Topology-preserving],
      [#{ tq6.observed_total - tq6.observed_breaks } (#str(calc.round((1 - tq6.rate_observed) * 100, digits: 1))%)],
      [#{ tq6.possible_total - tq6.possible_breaks } (#str(calc.round((1 - tq6.rate_possible) * 100, digits: 1))%)],
      [#{ tk43.possible_total - tk43.possible_breaks } (#str(calc.round((1 - tk43.rate_possible) * 100, digits: 1))%)],
      [---],
    [Total events], [#tq6.observed_total], [#fmtk(tq6.possible_total)], [#fmtk(tk43.possible_total)], [---],
    [Depletion fold], [---],
      [#str(calc.round(tq6.rate_possible / calc.max(tq6.rate_observed, 0.001), digits: 1))$times$],
      [#str(calc.round(tk43.depletion_fold, digits: 1))$times$],
      [---],
    [RR (95% CI)], [---],
      [---],
      [#str(calc.round(tk43.risk_ratio, digits: 2)) (#str(calc.round(tk43.risk_ratio_ci_95.at(0), digits: 2))--#str(calc.round(tk43.risk_ratio_ci_95.at(1), digits: 2)))],
      [---],
    [Hypergeometric _p_], [---], [$#str(tq6.hypergeom_p)$], [$#str(tk43.hypergeom_p)$], [---],
    [Permutation _p_], [---], [$lt.eq 10^(-4)$], [$lt.eq 10^(-4)$], [---],
  ),
  caption: [
    Topology avoidance under two adjacency definitions. $Q_6$: Hamming-1 edges in $"GF"(2)^6$; $K_4^3$: full single-nucleotide mutation graph (encoding-independent). Topology-breaking = creates a new amino acid graph disconnection. RR = risk ratio with 95% log-normal CI. Both adjacency definitions yield highly significant depletion of topology-breaking reassignments.
  ],
) <tbl:topo>

// Figure 3 — Evolutionary evidence (bit-bias + depth + topology avoidance + tRNA)
#figure(
  image("figures/Fig3_evolutionary_evidence.png"),
  caption: [
    Evolutionary evidence for structural constraints on codon reassignment. *(A)* Bit-position bias: distribution of bit-flips across $"GF"(2)^6$ coordinates in natural reassignment events; dashed line = uniform expectation. *(B)* Evolutionary depth calibration: reconnection $epsilon$ vs estimated divergence age (log scale) for 4 variant-code amino acids. *(C)* Topology avoidance: observed natural reassignments create new disconnections at 22.2% vs 72.7% among all possible single-codon changes (permutation $p lt.eq 10^(-4)$). *(D)* tRNA enrichment: rank of the reassigned amino acid among all 20 AAs by tRNA gene proportion in variant-code vs standard-code organism pairings; rank 1 = most enriched.
  ],
) <fig:topo>

== Explanatory modeling: topology as an independent predictor <sec:res-condlogit>

// --- dynamic conditional logit ---
#let clf = cl.model_fits
#let m3f = clf.at("M3_phys_topo", default: (:))
#let m4f = clf.at("M4_full", default: (:))
#let m2f = clf.at("M2_topo", default: (:))
#let m1f = clf.at("M1_phys", default: (:))
#let best_aicc = m3f.aicc

To test whether topology avoidance is reducible to physicochemical optimization, we fit event-level conditional logit models to the #cl.total_events reassignment event-steps across #cl.n_tables variant-code tables, treating each observed reassignment as a choice among $approx 1{,}280$ candidate single-codon moves (@sec:null-condlogit). The combined physicochemical-plus-topology model (M3) was strongly favored over all alternatives by AICc ($"AICc" = #str(calc.round(m3f.aicc, digits: 1))$), outperforming the physicochemistry-only model (M1; $Delta"AICc" = #str(calc.round(m1f.aicc - m3f.aicc, digits: 1))$) and the topology-only model (M2; $Delta"AICc" = #str(calc.round(m2f.aicc - m3f.aicc, digits: 1))$). Adding a heuristic tRNA-complexity proxy did not improve fit (M3$arrow.r$M4: $Delta"AICc" = #str(calc.round(m4f.aicc - m3f.aicc, digits: 1))$).

// Build model comparison table dynamically, sorted by AICc
#let model_rows = (
  ("M3_phys_topo", "Physicochemistry + topology", m3f),
  ("M4_full", "Physicochemistry + topology + tRNA", m4f),
  ("M2_topo", "Topology only", m2f),
  ("M1_phys", "Physicochemistry only", m1f),
).sorted(key: r => r.at(2).aicc)

#figure(
  table(
    columns: (auto, 1fr, auto, auto, auto, auto),
    align: (center, left, center, center, center, center),
    inset: 5pt,
    stroke: (x, y) => if y == 0 { (bottom: 0.7pt) } else { none },
    table.header([Model], [Description], [$k$], [Log $cal(L)$], [AICc], [$Delta$AICc]),
    ..model_rows.map(r => {
      let f = r.at(2)
      let delta = calc.round(f.aicc - best_aicc, digits: 1)
      let is_best = delta == 0
      (
        [#if is_best [*#r.at(0).split("_").at(0)*] else [#r.at(0).split("_").at(0)]],
        [#r.at(1)],
        [#f.n_params],
        [$#str(calc.round(f.log_likelihood, digits: 1))$],
        [#if is_best [*#str(calc.round(f.aicc, digits: 1))*] else [#str(calc.round(f.aicc, digits: 1))]],
        [#if is_best [*#str(delta)*] else [#str(delta)]],
      )
    }).flatten(),
  ),
  caption: [
    Event-level conditional logit model comparison. Each model assigns a probability to the observed reassignment among $approx 1{,}280$ candidate single-codon moves, with likelihood marginalized over event orderings within each table. $k$ = number of estimated parameters. Bold = best model by AICc.
  ],
) <tbl:condlogit>

// --- dynamic LR test results ---
#let lr_m1m3 = cl.lr_tests.at("M1_vs_M3", default: (:))
#let lr_m2m3 = cl.lr_tests.at("M2_vs_M3", default: (:))

Likelihood-ratio tests confirmed that each feature class adds substantial explanatory value to the other: adding topology to physicochemistry yields LR $= #str(calc.round(lr_m1m3.lr_statistic, digits: 1))$ ($p lt.double 10^(-10)$), and adding physicochemistry to topology yields LR $= #str(calc.round(lr_m2m3.lr_statistic, digits: 1))$ ($p lt.double 10^(-10)$). The two feature classes are only weakly associated across the full candidate landscape, indicating limited confounding.

In the best-fitting model (M3), evolution favored moves that reduce local physicochemical mismatch ($hat(beta)_"phys" = -0.004$ per Grantham unit) and strongly disfavored moves that increase amino acid codon-family disconnection ($hat(beta)_"topo" = -3.26$ per additional connected component). Under M3, observed natural reassignments rank on average at the 89.5th percentile among all candidate moves (@fig:condlogit, panel B), with the recurrent UGA$arrow.r$Trp reassignment consistently ranking above the 98th percentile. The one notable outlier is the yeast mitochondrial CUU$arrow.r$Thr reassignment (30th percentile), consistent with translation table 3 being the sole marginal exception in the per-table optimality analysis (@sec:res-pertable).

// Figure 5 — Conditional logit diagnostics
#figure(
  image("figures/Fig5_condlogit.png"),
  caption: [
    Event-level explanatory modeling of natural codon reassignments. *(A)* AICc comparison of four nested conditional logit models (M1--M4); lower is better. The combined physicochemistry-plus-topology model (M3) is decisively favored ($Delta"AICc" gt.eq 67$). *(B)* Distribution of observed move percentile ranks under M1 (physicochemistry only, light) versus M3 (combined, dark); M3 concentrates ranks toward the top of the candidate set. Dashed line = chance expectation (50th percentile). *(C)* Likelihood-ratio tests for each feature class added to its complement; both topology and physicochemistry contribute highly significant independent information ($p < 0.001$); the heuristic tRNA proxy does not ($p = 0.70$).
  ],
) <fig:condlogit>

== tRNA enrichment for reassigned amino acids <sec:res-trna>

Organisms with variant genetic codes show elevated tRNA gene copy numbers for the reassigned amino acid relative to standard-code controls. Across 24 pairings (18 variant-code organisms, 6 standard-code controls, spanning 5 variant genetic codes and 18 tRNAscan-SE--verified genome assemblies), Fisher's exact test combined via Stouffer's $Z$ method yields $p = 1.7 times 10^(-7)$ ($Z = 5.10$). To address non-independence from shared controls, we enumerated all maximal independent sets (MIS) from the conflict graph via Bron--Kerbosch; both MIS (each of size 6) are significant at $p < 0.05$ (worst-case $p = 0.045$).

The most striking case is _Tetrahymena thermophila_ (NCBI translation table 6, UAA/UAG reassigned to Gln; @hanyu1986), which carries 54 glutamine tRNA genes---including 39 suppressor tRNAs reading the reassigned stop codons---compared to 3 Gln tRNAs in the standard-code ciliate _Ichthyophthirius multifiliis_, 11 in _Stentor coeruleus_ and 3 in _Fabrea salina_ @zhang2022. The pattern extends across reassignment types. Among Gln-reassignment ciliates, _Pseudocohnilembus persalinus_ (20 Gln tRNAs including 15 suppressors) and _Halteria grandinella_ @zheng2021 (9 Gln tRNAs including 3 suppressors) represent independent lineages within Oligohymenophorea and Spirotrichea respectively. Among Cys-reassignment ciliates (translation table 10, UGA$arrow.r$Cys), six tRNAscan-SE--verified _Euplotes_ species all carry tRNA-Cys genes with the non-canonical TCA anticodon (reading UGA), with 1--4 such genes per species alongside standard GCA-anticodon Cys tRNAs.

However, the pattern is not universal. _Blastocrithidia nonstop_ (translation table 31) reassigned all three stop codons (as did _Condylostoma magnum_, where stop-codon function is context-dependent; @heaphy2016) but achieved UGA$arrow.r$Trp via anticodon stem shortening (5 bp $arrow.r$ 4 bp) of tRNA-Trp(CCA), combined with an eRF1 Ser74Gly mutation, rather than gene duplication @kachale2023. Similarly, _Mycoplasma_ species with UGA$arrow.r$Trp use a single tRNA-Trp with anticodon modification. These boundary cases define a three-tier mechanistic landscape: (i) tRNA gene duplication in large nuclear genomes, (ii) anticodon structural modification in streamlined genomes, and (iii) anticodon base modification in minimal genomes.

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

// tRNA enrichment panel is now in Fig 3D above; table-based analysis follows

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

The decomposition of $F$ by nucleotide position (@fig:coloring, panel b) reveals that the second codon position dominates the error-minimization signal, contributing #str(pos2_pct)% of the total physicochemical mismatch. This aligns with the biochemical observation that second-position mutations cause the largest changes in amino acid hydrophobicity @woese1965 #cite(<freeland1998>). The wobble position contributes only #str(pos3_pct)%, consistent with its largely synonymous character under the genetic code's degeneracy structure.

== Exploratory observations <sec:res-exploratory>

=== Bit-position bias in codon reassignments <sec:res-bitbias>

The distribution of bit-flips across the 6 coordinates of $"GF"(2)^6$ in natural codon reassignments shows apparent positional skew under a uniform null ($chi^2 = 16.26$, $p = 0.006$, $"df" = 5$). However, this signal is substantially attenuated after de-duplication to 20 unique (codon, target amino acid) pairs ($p = 0.075$) and vanishes entirely under a codon-preserving permutation null ($p = 1.0$). The apparent bias is therefore explained by which codons are "hot" for reassignment, not by a genuine positional preference in $"GF"(2)^6$.

=== Variant-code disconnection catalogue <sec:res-catalogue>

A systematic survey across all 25 NCBI translation tables identifies 4 variant-code amino acids with disconnected Hamming graphs at $epsilon = 1$: threonine in the yeast mitochondrial code (translation table 3), leucine in the chlorophycean mitochondrial code (translation table 16), alanine in _Pachysolen_ nuclear code (translation table 26), and a tripartite serine in the _Candida_ nuclear code (translation table 12; @santos1999). These cases, combined with the universal serine disconnection, constitute the complete inventory of amino acid graph disconnections at unit Hamming distance.

=== Atchley Factor 3 and Serine convergence <sec:res-atchley>

Serine has the most extreme Atchley Factor 3 score among the 20 amino acids ($F_3 = -4.760$, 2.24 SD below the mean; @atchley2005), and it is the only amino acid disconnected at $epsilon = 1$ under every base-to-bit encoding. This convergence is not coincidental: Factor 3 is a composite of molecular size and codon diversity, so both measures reflect the same underlying anomaly---Serine's disproportionate codon diversity (6 codons in two disconnected families) relative to its small physicochemical footprint. The $"GF"(2)^6$ framework provides a structural explanation for why Serine is an outlier on Factor 3, though the two views are complementary rather than independent.

== Retrospective cross-study reanalysis of synthetic genome recoding outcomes <sec:res-codonsafe>

To assess whether the structural properties identified in Sections 3.1--3.5 translate to measurable phenotypic consequences in synthetic biology, we performed a retrospective cross-study reanalysis of eight published genome recoding datasets comprising over 248,000 codon-level observations (@tbl:codonsafe-datasets). Each codon substitution was classified by its GF(2)#super[6] topology: whether it crosses a connected-component boundary at $epsilon = 1$, the change in local physicochemical mismatch cost ($Delta F_"local"$, Grantham distance), and the Hamming distance between source and target vectors.

#figure(
  table(
    columns: (auto, auto, auto, auto, auto),
    align: (left, left, right, center, left),
    inset: 6pt,
    stroke: (x, y) => if y == 0 { (bottom: 0.7pt) } else { none },
    table.header([*Study*], [*Swap type*], [*$n$*], [*Boundary?*], [*Key result*]),
    [#cite(<robertson2025syn57>, form: "prose") Syn57], [4 Ser + 2 Ala + Stop], [60,240], [Ser=100%, Ala=0%], [Primary contrast],
    [#cite(<fredens2019>, form: "prose") Syn61], [Ser (TCG/TCA$arrow.r$AGC/AGT)], [18,218], [All cross], [Informative null],
    [#cite(<ostrov2016>, form: "prose")], [7 codon types], [62,214], [Mixed], [Segment viability],
    [#cite(<napolitano2016>, form: "prose")], [Arg (AGR$arrow.r$CGN)], [12,888], [None], [SRZ covariates],
    [#cite(<robertson2025ochre>, form: "prose") Ochre], [Stop (TGA$arrow.r$TAA)], [1,195], [N/A], [Growth data],
    [#cite(<nyerges2024>, form: "prose")], [7 codon types], [62,007], [Mixed], [Multi-omics],
    [#cite(<lajoie2013>, form: "prose") C321.$Delta$A], [Stop (UAG$arrow.r$UAA)], [321], [N/A], [Baseline],
    [#cite(<frumkin2018>, form: "prose")], [Arg (CGU/CGC$arrow.r$CGG)], [60], [None], [Translation efficiency],
    [#cite(<ding2024>, form: "prose")], [Ser (TCG, mammalian)], [Variable], [Yes], [Cross-organism],
  ),
  caption: [
    Published genome recoding datasets analyzed. "Boundary?" indicates whether the synonymous codon swaps cross a connected-component boundary at $epsilon = 1$ in GF(2)#super[6]. #super[#sym.dagger]Quantitative extraction pending; included for scope. Only 2 of 8 fully analyzed datasets (Syn57 and Ding) have within-study variation in the boundary-crossing feature.
  ],
) <tbl:codonsafe-datasets>

// Figure 4 — Translational applications (synbio + decomposition + catalogue)
#figure(
  image("figures/Fig4_translational.png"),
  caption: [
    Translational applications and prediction catalogue. *(A)* Feasibility landscape of 1,280 single-codon reassignments from the standard code, colored by degeneracy filtration status. All filtration-preserving variants score $gt.eq 0.8$; all filtration-breaking variants score $lt.eq 0.75$. Dashed line: high-feasibility threshold. *(B)* Grantham mismatch score decomposition by nucleotide position: position 2 dominates (#str(pos2_pct)%), consistent with the biochemical hierarchy of mutational impact. *(C)* Prediction catalogue: 15 evaluated claims distributed across 5 workstreams (WS1--WS6), colored by verification status.
  ],
) <fig:translational>

=== A three-layer decomposition of recoding outcomes <sec:res-threelayer>

The cross-study reanalysis revealed three separable layers of codon-space structure, each with distinct predictive scope.

==== Boundary crossing: validated but not predictive of transcriptomic perturbation

In the #cite(<robertson2025syn57>, form: "prose") Syn57 dataset---the only large-scale dataset with within-study variation in boundary crossing---all 37,146 serine recodings crossed the UCN$arrow.l.r$AGY disconnection boundary at $epsilon = 1$, while all 22,859 alanine recodings remained within the compact GCN family. This provides an exact within-study topology contrast. However, genes containing exclusively serine-type versus exclusively alanine-type recodings did not differ in transcriptomic perturbation (Mann--Whitney $p = 0.40$, $n_"Ser" = 129$, $n_"Ala" = 49$), and the fraction of boundary-crossing recodings per gene was uncorrelated with $|log_2 "FC"|$ ($rho = 0.012$, $p = 0.46$, $n = 3{,}510$ genes). In a separate dataset, all 18,218 Syn61 serine recodings @fredens2019 crossed the boundary, yet 99.96% were viable (7 corrections needed). Boundary crossing is therefore a coarse geometric property of codon-family organization that is validated as a structural distinction but not as a predictor of transcriptomic perturbation.

==== Local mismatch geometry: a fine-scale burden axis with positive signal

The second layer examines whether reassignment changes the physicochemical quality of the immediate mutational neighborhood. In the Syn57 contrast, serine swaps moved codons into better local Grantham neighborhoods ($Delta F_"local" = -37.2 plus.minus 30.5$), whereas alanine swaps moved codons into worse neighborhoods ($Delta F_"local" = +19.0 plus.minus 0.0$); the distributions were completely separated (Mann--Whitney $U = 0$, $p < 10^(-16)$, $n = 60{,}240$). In a topology-fixed setting---arginine synonymous recodings @napolitano2016, where all 12,888 AGR$arrow.r$CGN swaps remain within one connected component at $epsilon = 1$---local mismatch change correlated with established recoding-burden covariates: $Delta F_"local"$ vs RBS deviation ($rho = -0.33$, $p < 10^(-15)$) and vs mRNA folding deviation ($rho = -0.12$, $p < 10^(-15)$). Among the four CGN targets, CGG has the lowest local mismatch cost (246 Grantham units) while CGC has the highest (421), and target-specific RBS deviation spans four orders of magnitude, suggesting that local mismatch and SRZ covariates capture partially independent dimensions of recoding burden.

==== Design deviations: accessibility-dominated, not topology-directed

To test whether engineered systems escape design constraints along topology-defined low-burden directions, we compared the Syn57 design genome (Data file S2; @robertson2025syn57) against the final verified genome (Data file S8) and identified 727 CDS-level codon differences. After filtering 7 genes with $>$10 deviations each (structural rearrangements, predominantly _bioF_ with 311 changes), 326 genuine point deviations remained (267 synonymous, 59 nonsynonymous). Among the 163 deviations that changed local mismatch, 76 moved to a better neighborhood and 87 to a worse one---no directional bias (two-sided exact binomial $p = 0.43$, Clopper--Pearson 95% CI $[0.39, 0.55]$; Wilcoxon signed-rank $p = 0.67$). This null result was stable across all filtering thresholds tested ($p > 0.5$ at cutoffs 3, 5, 10, 15, 20, and 50 deviations per gene). By contrast, deviations strongly favored mutationally proximate moves: 83% were single-bit changes (Hamming distance 1 under the default encoding; note that Hamming distance is encoding-dependent, though nucleotide edit distance is invariant), and the mean Hamming distance was 1.21. Thus, realized deviations from the Syn57 design do not preferentially descend toward lower local mismatch; instead, they follow accessibility-dominated escape routes, favoring nearby codon states regardless of neighborhood quality. The #cite(<ding2024>, form: "prose") mammalian TCG recoding independently confirms the universality of the serine disconnection across standard-code organisms, providing cross-kingdom validation.

==== Segment-level recoding burden in the Ostrov 57-codon design

As an orthogonal test, we examined segment-level outcomes from the #cite(<ostrov2016>, form: "prose") 57-codon genome design. Among 44 segments with fitness data (of 87 total), segments with more recoded codons in essential genes showed a suggestive correlation with worse doubling time ($rho = 0.34$, raw $p = 0.022$, Bonferroni-corrected $p = 0.066$), while total recoding load showed no association ($rho = -0.10$, $p = 0.53$). Problem segments (13 with lethal exceptions) trended toward higher essential-gene recoding load (73.5 vs 39.8 recoded sites, Mann--Whitney $p = 0.086$). Among 44 fully characterized segments, 17 showed spontaneous codon reversions, indicating ongoing design instability at positions where the engineered code was under selection pressure.

=== Summary of cross-study reanalysis scope

The synthetic recoding analyses support a three-layer view: boundary crossing captures a real family-level structural distinction but did not predict transcriptomic perturbation; local mismatch geometry captures a finer-scale burden axis with positive signal; and realized design deviations are dominated by mutational accessibility rather than systematic escape toward lower local mismatch. These layers are partially separable: different structural descriptors appear to matter at different levels of the problem---family topology for the geometry of reassignment space, local mismatch for some aspects of recoding burden, and Hamming accessibility for the short-range routes by which engineered designs deviate or revert.

== Falsified and rejected claims <sec:res-negatives>

=== KRAS--Fano clinical prediction <sec:res-kras>

The conjecture that XOR ("Fano") relationships in $"GF"(2)^6$ predict co-mutation enrichment at KRAS G12 sites was tested against 1,670 mutations from MSK-IMPACT @zehir2017 and cleanly falsified ($p = 1.0$; details in Supplementary Material). This negative result separates code-level error-minimization (which is real) from mutation-level algebraic predictions (which are not).

=== Serine distance-4 invariant <sec:res-serine>

The claim that Serine's minimum inter-family Hamming distance (UCN$dash$AGY) equals 4 under all 24 base-to-bit encodings is false. Of the 24 encodings, 16 yield minimum distance 2 and only 8 yield distance 4. The distance-4 result obtains only when both nucleotide pairs distinguishing UCN from AGY ($U tilde.op A$ and $C tilde.op G$ in the first two positions) are encoded at maximal Hamming distance. The correct encoding-invariant statement is: Serine is disconnected at $epsilon = 1$ under every encoding, and its inter-family distance ($gt.eq 2$) is the largest among the three 6-codon amino acids (Leucine and Arginine both have inter-family distance 1).

=== PSL(2,7) symmetry and holomorphic embedding <sec:res-psl>

The claim that PSL(2,7) is the fundamental symmetry group of the genetic code was pre-rejected by #cite(<antoneli2011>, form: "prose"), who showed that PSL(2,7) has no 64-dimensional irreducible representation (its irreps have dimensions 1, 3, 6, 7, 8). The claim that the coordinate-wise map $"GF"(2)^6 arrow.r CC^3$ sending base pairs to fourth roots of unity is a holomorphic embedding extending a character of $"GF"(8)^*$ is also incorrect: the domain is a finite discrete set (not a complex manifold), and the map fails the character identity $chi(x + x) = chi(x)^2$ since $i^2 = -1 eq.not 1$.

== Negative results with informative interpretation <sec:res-infoneg>

The local mismatch cost test---asking whether reassigned codons sit in "worse" Hamming neighborhoods with higher Grantham distance to their neighbors---yields a Mann--Whitney $U = 301$, $p = 0.70$. Reassignment is not driven by local escape from costly neighborhoods. This null result is scientifically informative: it indicates that the topology avoidance constraint (Section 3.4) operates at the global graph-connectivity level, not at the local per-codon level.


// ============================================================
//  4. DISCUSSION
// ============================================================
= Discussion <sec:discussion>

== An information-theoretic view of the genetic code <sec:disc-info>

The central finding of this work is that the standard genetic code minimizes the physicochemical disruption caused by single-bit errors in $"GF"(2)^6$ coordinates. This is not a new conclusion---#cite(<freeland1998>, form: "prose") established error-minimization using nucleotide-level mutation models---but the hypercube representation makes the optimality principle geometrically explicit: the code is a good _coloring_ of a structured graph, in the graph-theoretic sense that adjacent vertices (codons differing by one bit) tend to share labels (amino acids) or, when they differ, differ by small physicochemical distances.

The score decomposition (@fig:coloring) shows that this optimization is concentrated at the second codon position (49.3% of total mismatch) and first position (38.2%), with the wobble position contributing only 12.5%. This gradient mirrors the biochemical hierarchy of mutational impact and is an emergent property of the code's structure rather than a parameter of the model.

The robustness result (@fig:coloring) demonstrates that optimality is not an artifact of restricting attention to $Q_6$: when the full $K_4^3$ mutation graph is considered ($rho = 1$), the signal strengthens. The code minimizes error not just along the hypercube edges but across the complete space of single-nucleotide substitutions. This error-minimization is complementary to, but distinct from, the finding of #cite(<itzkovitz2007>, form: "prose") that the code is also nearly optimal for carrying parallel regulatory information within protein-coding sequences---a property related to stop-codon identity rather than amino acid physicochemistry.

A natural question is whether the $"GF"(2)^6$ framework adds genuine insight beyond what a direct analysis of $K_4^3$ would provide. We note three advantages. First, the binary decomposition enables the $rho$-interpolation that reveals the relationship between $Q_6$ and $K_4^3$ optimality as a continuum rather than two disconnected analyses. Second, the Hamming-distance filtration provides a natural persistence parameter ($epsilon$) for studying amino acid graph connectivity, yielding the topology avoidance result. Third, the encoding-sensitivity analysis (24 bijections) distinguishes encoding-invariant properties (Serine disconnection) from encoding-dependent ones (distance-4 claim), a distinction invisible in the nucleotide-level representation. However, we emphasize that the Hamming-1/diagonal decomposition does not correspond to the biological transition/transversion partition: 16 of 24 encodings mix transitions and transversions equally among Hamming-1 edges, while 8 encodings place both transitions on diagonal edges (Supplementary Table S2). The $rho$ parameter is best interpreted as a diagonal-edge weight, not a transition/transversion weight.

== Evolutionary preservation and topology avoidance <sec:disc-evol>

The per-table analysis (@fig:coloring) shows that 24 of 25 NCBI translation tables maintain coloring optimality under their own block-preserving null, suggesting that codon reassignment events are constrained to preserve error-minimization. The single marginal exception (translation table 3, yeast mitochondrial) is the most extensively reassigned code (6 changes), yet still falls only slightly above the 5% threshold (quantile 7.4%).

The topology avoidance result (@fig:topo) provides a mechanistic complement: natural reassignments avoid creating new amino acid disconnections at a rate far below chance expectation. Importantly, this depletion is not confined to the $Q_6$ representation: under the full $K_4^3$ single-nucleotide mutation graph, the observed proportion of topology-breaking reassignments remains 22.2% against a possible rate of 66.1% ($"RR" = 0.34$, 95% CI $[0.17, 0.68]$, $p lt.eq 10^(-4)$). The signal is somewhat attenuated relative to the $Q_6$ result ($"RR" = 0.31$) because some $Q_6$-disconnected pairs become connected when all single-nucleotide edges are admitted, but the main qualitative conclusion survives: the connected-component structure of amino acid codon families is functionally important, regardless of whether connectivity is defined on the hypercube subgraph or on the biologically fuller single-substitution graph. The yeast mitochondrial threonine reassignment illustrates this cost: the CUN$arrow.r$Thr change required acquisition of a novel tRNA#super[Thr] derived from tRNA#super[His] via anticodon mutation @su2011, creating the topology-breaking disconnection that makes translation table 3 the sole marginal outlier in the per-table optimality analysis.

Together, these analyses support a four-part picture: the standard code is unusually low-cost, most variant codes preserve that structure, realized natural reassignments avoid topology-disrupting moves, and event-level explanatory modeling shows that this avoidance carries independent information beyond physicochemical optimization. The conditional logit analysis (@tbl:condlogit) makes the last step precise: the topology term retains major explanatory value ($Delta"AICc" = 100$ relative to a physicochemistry-only model) even after accounting for local physicochemical cost, while physicochemical cost also retains value ($Delta"AICc" = 67$) after accounting for topology. The two signals are therefore complementary rather than redundant, consistent with the weak correlation ($r_s = 0.15$) between the feature classes across the candidate landscape. This does not by itself prove direct selection on topology as an abstract graph property; topology may still proxy aspects of decoding architecture---such as tRNA cross-recognition range or ribosomal A-site constraints---not explicitly modeled here. However, within the tested event landscape, topology avoidance cannot be dismissed as a by-product of physicochemical optimization. Topology avoidance is not a geometric restatement of error minimization, but a complementary constraint on how genetic codes can change while remaining evolutionarily accessible.

== Mechanistic implications: tRNA compensation <sec:disc-trna>

The tRNA enrichment result (@fig:topo, panel D) links the geometric observation to molecular mechanism. In several variant-code lineages where codon reassignment disrupts connectivity, expanded tRNA gene repertoires for the affected amino acid are observed, consistent with compensatory gene duplication as one evolutionary route to accommodation. The extreme case of _Tetrahymena thermophila_ (54 Gln tRNAs, including 39 suppressors) illustrates the scale of genomic response required to service a split codon family.

However, the boundary cases---_Blastocrithidia nonstop_ (anticodon stem shortening; @kachale2023) and _Mycoplasma_ species (anticodon modification)---show that tRNA gene duplication is not the only evolutionary solution. The three-tier pattern (duplication in large genomes, structural modification in intermediate genomes, base modification in minimal genomes) suggests that genome size constrains the available mechanistic repertoire for codon reassignment. This pattern is now supported by 18 tRNAscan-SE--verified genomes across 5 variant genetic codes (Tables 4, 6, 10, 15, and 31), including 6 _Euplotes_ species (UGA$arrow.r$Cys) that each carry 1--4 TCA-anticodon tRNA-Cys genes dedicated to reading UGA.

== Synthetic recoding outcomes: a three-layer interpretation <sec:disc-layers>

The retrospective cross-study reanalysis helps delimit which structural quantities in GF(2)#super[6] are functionally relevant and at which biological layer. Boundary crossing is a valid codon-space distinction---serine and alanine provide a clean contrast in which one class always crosses a family boundary and the other never does---yet this distinction did not predict RNA-seq perturbation in Syn57. We therefore do not interpret codon-family boundary crossing as a general transcriptomic burden variable. Its relevance, if any, likely lies at a different biological layer, such as decoding architecture, translational fidelity, or evolutionary accessibility, rather than steady-state expression change.

By contrast, local mismatch geometry showed a clearer positive signal. In both the Syn57 contrast and the topology-fixed arginine setting, codon changes differed systematically in the quality of the physicochemical neighborhood they entered, and these differences aligned with established recoding-burden covariates @napolitano2016. This suggests that local mismatch is best understood as a fine-scale recoding-friction variable: not a universal master predictor, but a meaningful descriptor of neighborhood quality within an already specified codon-family structure.

The design-deviation analysis further clarifies scope. Genuine deviations from the Syn57 design were overwhelmingly short-range (83% single-bit in GF(2)#super[6]), with no directional bias in mismatch change. This indicates that realized escape is governed more by mutational accessibility than by systematic optimization of local neighborhood cost. In other words, once a design is under strain, the available exits appear to be chosen primarily by mutational proximity, not by a global preference for lower mismatch.

These findings argue against treating "topology" as a single explanatory variable. Instead, different structural descriptors appear to matter at different levels: family topology for the geometry of reassignment space (Section 3.4), local mismatch for some aspects of recoding burden (Section 3.7.2), and Hamming accessibility for the short-range routes by which engineered designs deviate. The fact that the #cite(<fredens2019>, form: "prose") Syn61 tolerated 18,218 boundary-crossing serine swaps with 99.96% success, while the NCBI variant codes show 3.3-fold depletion of topology-breaking changes over evolutionary timescales, suggests that boundary crossing may be more relevant to evolutionary trajectory constraints than to acute engineering costs.

== Scope and falsified claims <sec:disc-honest>

Six of 15 evaluated claims were rejected, falsified, or tautological (shaded rows in @tbl:claims; full details in Supplementary Material). These failures help delimit the framework's explanatory scope. The KRAS--Fano prediction ($p = 1.0$) cleanly separates code-level error-minimization from mutation-level algebraic predictions. The Serine distance-4 invariant, falsified by 16 of 24 encodings giving distance 2, illustrates the importance of systematic encoding-sensitivity testing. The PSL(2,7) and holomorphic embedding claims were pre-rejected by existing mathematical results @antoneli2011.

== Limitations <sec:disc-limits>

Several limitations should be noted. First, the choice of base-to-bit encoding is not unique, and while the coloring optimality result holds across all 24 encodings, the specific score values and rank orderings are encoding-dependent (full encoding-sensitivity results across all 24 bijections are provided in the Supplementary Material). Second, the block-preserving null model, while standard in the field @freeland1998, preserves more structure than a fully random code and may understate the degree of optimality. Third, the tRNA enrichment result, while robust to pairing selection (worst-case MIS $p = 0.045$), relies on a small number of independent pairings ($n = 6$) with limited statistical power; this result is appropriately classified as suggestive and should not be interpreted as demonstrating a universal compensation mechanism. Fourth, the conditional logit model (@sec:res-condlogit) is an event-level explanatory model, not direct proof of biological causality; the candidate universe comprises all $approx 1{,}280$ single-codon reassignments regardless of mechanistic plausibility, and the tRNA-complexity proxy (Hamming distance to nearest target-AA codon) is heuristic rather than mechanistically grounded. The finding that this particular proxy did not improve model fit does not exclude the possibility that a richer model of tRNA repertoire change would contribute explanatory value.

Whether the code's optimality reflects adaptive selection or non-adaptive carry-over from a primordially constrained starting point @koonin2009 cannot be resolved by cross-sectional data alone. Multiple theoretical frameworks address this question: #cite(<sella2006>, form: "prose") showed that coevolution of genes and codes generates error-correcting structure resembling the standard code; #cite(<vetsigian2006>, form: "prose") argued that communal evolution via horizontal gene transfer accounts for both universality and optimality; and #cite(<novozhilov2007>, form: "prose") found the standard code sits roughly halfway up a local peak in a rugged fitness landscape (see @digiulio2005 for a review of origin theories). The topology avoidance result is consistent with all these frameworks: it demonstrates a constraint on reassignment trajectories but does not distinguish whether the constraint is selective or structural. #cite(<novozhilov2009>, form: "prose") showed that putative primordial 16-supercodon codes are "nearly optimal" even without direct selection, suggesting the current code's optimality may be partly inherited.

The $"GF"(2)^6$ representation is best understood as an analytical decomposition tool rather than a claim about the biological primacy of binary coordinates. It enables the $rho$-sweep interpolation between $Q_6$ and $K_4^3$ and facilitates systematic encoding-sensitivity analysis. However, we caution that the Hamming-1/diagonal decomposition in $"GF"(2)^6$ does not correspond to the biological transition/transversion partition: 16 of 24 encodings mix transitions and transversions equally among Hamming-1 edges, while 8 encodings place both transitions on diagonal edges. The $rho$ parameter should therefore be interpreted as a diagonal-edge weight, not a proxy for transition/transversion bias. The underlying biology is the assignment of chemically similar amino acids to mutationally proximate codons---a property that holds regardless of the coordinate system, as confirmed by the K4^3 topology avoidance result.


// ============================================================
//  5. CONCLUSION
// ============================================================
= Conclusion <sec:conclusion>

Across four amino acid distance metrics and across a continuum of mutation graphs spanning the hypercube and the complete single-nucleotide substitution graph, the standard genetic code consistently falls in the extreme low-cost tail of block-preserving null models ($p lt.eq 0.006$ under all metrics and all $rho$ values). That near-optimality is largely preserved across naturally variant translation tables (24 of 25 significant after FDR correction), implying that most reassignment histories remain within a constrained region of code space rather than freely exploring it.

Two constraints on reassignment trajectories emerge from the event-level analyses. First, evolution favors physicochemical smoothness: reassigned codons tend to enter neighborhoods with lower local mismatch cost. Second, evolution avoids fragmenting amino acid codon families: topology-breaking moves are $approx$3-fold depleted relative to the landscape of possibilities, and discrete-choice modeling shows that this topology term improves explanatory fit even after accounting for physicochemical cost ($Delta"AICc" gt.eq 67$), and vice versa. This supports the view that topology avoidance is not a geometric restatement of physicochemical optimization, while remaining agnostic about the precise mechanism---whether tRNA cross-recognition range, ribosomal decoding constraints, or other aspects of translational architecture not explicitly modeled here. Several variant-code lineages with topology-breaking reassignments show expanded tRNA gene repertoires for the reassigned amino acid (worst-case MIS Stouffer $p = 0.045$), consistent with compensatory accommodation routes.

In engineered recoding datasets ($>$248,000 codon-level observations across eight studies), these same codon-space descriptors appear to operate at different biological layers: family boundary crossing is a real structural distinction but was not predictive of transcriptomic perturbation, whereas local neighborhood mismatch more often aligned with established recoding-burden covariates and design deviations were dominated by mutational accessibility. Natural code evolution appears constrained along two partly independent axes: physicochemical smoothness and codon-family topological integrity.


// ============================================================
//  ACKNOWLEDGEMENTS
// ============================================================
#heading(numbering: none)[Acknowledgements]

We thank the NCBI, GtRNAdb, and cBioPortal teams for maintaining public databases essential to this work. tRNAscan-SE 2.0.12 was developed by Chan and Lowe at UC Santa Cruz.

#heading(numbering: none)[Data availability]

All code and data are available in the `codon-topo` repository (version #stats._version). Install with `pip install -e ".[all]"` or `uv sync --all-extras`. Analyses are fully reproducible via `codon-topo all --output-dir=./output --seed=135325`, which generates all JSON data files including `manuscript_stats.json` from which every inline statistic in this manuscript is rendered. NCBI genome assembly accessions are listed in @tbl:organisms.

#heading(numbering: none)[Declaration of competing interest]

The authors declare no competing interests.

// ============================================================
//  REFERENCES
// ============================================================
#pagebreak()

#bibliography("references.bib", title: "References", style: "elsevier-harvard")
