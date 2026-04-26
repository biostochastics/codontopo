// ============================================================
// Manuscript: Error-minimizing structure of the genetic code
// Status: pre-submission
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

// Figures: tighter gap. Image figures float (auto) so multi-panel figures
// can take their own page when needed; tables are pinned inline so they
// never land mid-paragraph in a different section.
#set figure(gap: 0.8em, placement: auto)
#show figure.where(kind: table): set figure(placement: none)
// Allow long tables to break across pages so they don't push to a new page
// and leave a large gap at the bottom of the previous page.
#show figure.where(kind: table): set block(breakable: true)

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

// Helper: format a small probability in scientific notation (e.g. 4.78 × 10^-8)
// Uses math mode so the exponent renders as a true superscript without
// the surrounding parentheses that the content-mode form leaves in place.
#let sci(n, sig: 2) = {
  if n == 0 {
    [0]
  } else if calc.abs(n) >= 0.001 and calc.abs(n) < 1 {
    str(calc.round(n, digits: sig + 1))
  } else {
    let exp = int(calc.floor(calc.log(calc.abs(n), base: 10)))
    let mant = calc.round(n / calc.pow(10.0, exp), digits: sig)
    $#mant times 10^(#exp)$
  }
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
    #super[\*]Corresponding author: sergey\@kornilov.bio
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
  The standard genetic code reduces the impact of point mutations, but the robustness of this property across physicochemical metrics, naturally variant codes, and codon-reassignment mechanisms remains incompletely quantified. We embed the 64 codons in $"GF"(2)^6$, representing the hypercube $Q_6$ as a coordinate-dependent subgraph of the encoding-independent single-nucleotide mutation graph $H(3,4)$, which supports continuous $rho$-interpolation between the two and enables joint analysis of physicochemical error minimization and codon-family topology. Under a block-preserving null ($n = $10,000), the standard code is significantly low-cost across four distinct amino-acid distance metrics (Grantham $p = 0.006$; Miyata $p < 0.001$; Woese polar requirement $p = 0.003$; Kyte--Doolittle hydropathy $p = 0.001$), addressing the concern that prior optimality results could be metric-specific; the signal strengthens monotonically as $rho$ moves $Q_6 arrow.r H(3,4)$. Across the 27 NCBI translation tables, near-optimality is broadly preserved: 11 of the 12 informative-distance variants retain top-5% placement after BH--FDR correction (yeast mitochondrial is the sole marginal exception). Natural codon reassignments rarely break codon-family connectivity: under $H(3,4)$, only #tk43.observed_breaks of #tk43.observed_total observed events are topology-breaking versus #str(calc.round(tk43.rate_possible * 100, digits: 0))% of 1,280 candidate moves (RR #str(calc.round(tk43.risk_ratio, digits: 2)), permutation $p lt.eq 10^(-4)$). This depletion is robust to alternative topology definitions, clade exclusions, and base-to-bit encodings, although the small breaker subset (4 of 6 from a single yeast-mitochondrial lineage; denominator effect in clade-exclusion robustness) is underpowered for strong cross-clade inference. Conditional-logit decomposition shows that topology avoidance and local physicochemical cost provide complementary, only weakly correlated signal ($r_s = #str(calc.round(cl.phys_topo_rho, digits: 2))$); a heuristic tRNA-distance proxy does not improve fit, and several variant-code lineages show suggestive tRNA-gene enrichment for reassigned amino acids. Retrospective reanalysis of nine genome-recoding datasets is consistent with---but does not establish---a working hypothesis in which codon-family topology operates at a different biological layer from acute cellular fitness: Syn61 tolerated 18,218 boundary-crossing serine swaps as a class (genome-wide, not per-codon-position viability), yet the same move type is #str(calc.round(tk43.depletion_fold, digits: 1))-fold depleted across natural code evolution. The contribution is the second axis: code evolution is jointly constrained by physicochemical smoothness and codon-family topological integrity, and these two constraints are partly independent.

  #v(0.5em)
  *Keywords:* genetic code evolution; error minimization; GF(2)#super[6] representation; codon reassignment; physicochemical distance; codon-family topology; conditional logit; tRNA gene enrichment; genome recoding
]

#v(0.8em)

#block(
  width: 100%,
  inset: (x: 1.5cm),
)[
  #text(weight: "bold")[Highlights]
  #set list(indent: 0pt, body-indent: 0.6em, marker: [•])
  - Codon-space geometry links genetic-code robustness and reassignment paths
  - Standard and variant codes preserve broad physicochemical error minimization
  - Reassignments are depleted for codon-family topology-breaking moves
  - Conditional-logit models separate topology from physicochemical similarity
  - Synthetic recoding shows boundary conditions for natural-code constraints
]

#v(1em)
#set par(first-line-indent: 1.5em)

// ============================================================
//  1. INTRODUCTION
// ============================================================
= Introduction <sec:intro>

The standard genetic code maps 61 sense codons to 20 amino acids through a pattern that has long been recognized as non-random. #cite(<woese1965>, form: "prose") first noted that similar codons tend to encode amino acids with similar physicochemical properties, and #cite(<freeland1998>, form: "prose") demonstrated quantitatively that the standard code sits among approximately 1 in $10^6$ random codes for mutational error minimization under the Woese polar requirement distance. These findings established that the code's structure reduces the fitness impact of point mutations, but whether this optimality is specific to a single metric or robust across distinct physicochemical parameterizations has not been systematically tested.

Every codon comprises three nucleotides drawn from ${C, U, A, G}$. By choosing a bijection $phi: {C, U, A, G} -> "GF"(2)^2$---for instance, $C |-> (0,0)$, $U |-> (0,1)$, $A |-> (1,0)$, $G |-> (1,1)$---each codon maps to a vertex of the 6-dimensional binary hypercube $Q_6 = "GF"(2)^6$. The genetic code then becomes a _coloring_ of $Q_6$ by 21 labels (20 amino acids plus the stop signal), and single-nucleotide mutations correspond to edges of $Q_6$ or, more precisely, to a subgraph of the complete mutation graph $H(3,4)$.

This representation is not new in principle---binary encodings of the genetic code appear in the mathematical biology literature (e.g., @antoneli2011). Its value lies not in the encoding itself but in the analytical decomposition it enables: the complete single-nucleotide mutation graph $H(3,4)$ (288 edges) separates into 192 Hamming-distance-1 edges (single-bit changes) and 96 within-nucleotide distance-2 edges (both bits of one nucleotide position flip simultaneously), allowing systematic interpolation via a weight parameter $rho in [0,1]$. Note that this decomposition does not correspond to the biological transition/transversion partition: under 16 of the 24 possible 2-bit encodings, the Hamming-1 edges contain an equal mixture of transitions and transversions per nucleotide position (2 of each among 4 Hamming-1 pairs), while the remaining 8 encodings place both transitions on diagonal (Hamming-2) edges. The $rho$ parameter should therefore be interpreted as a diagonal-edge inclusion weight, not a transition/transversion weight. Previous work has not exploited this decomposition to test error-minimization across multiple physicochemical metrics, nor extended the analysis to variant genetic codes. Doing so lets us examine two questions that single-axis analyses cannot. First, whether error-minimization is a property of the code itself rather than of any particular distance function---if four metrics with disjoint conceptual bases (composition--polarity--volume, polarity--volume, chromatographic polar requirement, hydropathy) all place the code in the extreme low-cost tail, the optimization target is amino-acid similarity however operationalized, not a metric-specific quirk. Second, whether codon-family _connectivity_ acts as a second, partly independent axis of constraint on natural reassignment events, distinct from physicochemical cost. If so, code evolution is bounded not only by which codes are low-cost but by which trajectories through code-space leave the decoding substrate intact---a constraint on transitions, not just on states. These two motivating questions are operationalized as three concrete tests:

+ *Is the standard code optimal under multiple physicochemical metrics?* We test whether the standard code's edge-mismatch score is extreme relative to block-preserving null models across four distinct physicochemical distance metrics (Grantham, Miyata, Woese polar requirement, and Kyte--Doolittle hydropathy), extending #cite(<freeland1998>, form: "prose") from a single metric to a cross-metric robustness envelope.

+ *Is this structure preserved by evolution?* We ask whether variant genetic codes (NCBI translation tables 2--33) maintain error-minimization, and whether natural codon reassignment events preferentially avoid disrupting the topological connectivity of amino acid codon families.

+ *What are the genomic correlates of disruption?* We test whether organisms whose variant codes break the connectivity of an amino acid's codon graph show elevated tRNA gene copy numbers for the affected amino acid, using tRNAscan-SE--verified data from 18 genomes spanning 5 variant genetic codes across Alveolata, Opisthokonta, Excavata, and Mollicutes.

To delimit the framework's scope, we also test four conjectural extensions previously associated with this project (one of us, P.C.; @clayworth2026): a Serine distance-4 invariance claim, PSL(2,7) symmetry, a holomorphic embedding extending a character of $"GF"(8)^*$, and a KRAS--Fano clinical prediction. These tests distinguish supported graph-theoretic results from encoding artifacts and unsupported algebraic extensions (@sec:res-negatives; KRAS--Fano detail in Supplement §S13).

The paper is organized as follows. Section 2 describes the encoding formalism, graph decomposition, null models, and statistical methods, including a new retrospective cross-study reanalysis of nine published genome recoding datasets (with quantitative analysis of eight). Section 3 presents the four supported findings (cross-metric coloring optimality, per-table preservation, $rho$-robustness, and topology avoidance), the suggestive tRNA enrichment result, the cross-study reanalysis of synthetic recoding outcomes, and exploratory observations. Section 4 discusses the graph-theoretic interpretation, the relationship to frozen-accident versus adaptive hypotheses @koonin2009, and an exploratory three-layer interpretation of how codon-space structure relates to recoding outcomes. All analyses are reproducible via the open-source `codon-topo` pipeline (version #stats._version, seed 135325).


// ============================================================
//  2. METHODS
// ============================================================
= Methods <sec:methods>

== Binary encoding of the genetic code <sec:encoding>

We encode each nucleotide base as a 2-bit vector via the default bijection $phi: C |-> (0,0)$, $U |-> (0,1)$, $A |-> (1,0)$, $G |-> (1,1)$. A codon $b_1 b_2 b_3$ is then represented as the concatenation $phi(b_1) || phi(b_2) || phi(b_3) in "GF"(2)^6$, and the 64 codons become the 64 vertices of the 6-dimensional hypercube $Q_6$. This default bijection was adopted in companion methodological work (#cite(<clayworth2026>, form: "prose")) because it places the standard code's nine 2-fold-degenerate amino acids on bit-5 differences. The visualization-clarity rationale, however, is mildly circular: the bit-5 two-fold filtration is itself classified as Tautological in the present claim hierarchy (Supplement §S1), so the encoding was chosen to make a tautological property visible. We therefore use the default encoding only as an analytical convenience and make no claim that it is biologically privileged; all encoding-dependent results are reported under all 24 base-to-bit bijections (end of Section 2.1), and the load-bearing topology-avoidance result uses the encoding-independent $H(3,4)$ adjacency (Section 2.3.4).

Under this encoding, two codons that differ by a single-nucleotide substitution in which only one bit of the 2-bit pair changes are adjacent in $Q_6$ (Hamming distance 1). However, transversions that flip both bits within a nucleotide position correspond to Hamming distance 2. The full single-nucleotide mutation graph is therefore the Hamming graph $H(3,4) = K_4 #h(2pt) square.stroked.tiny #h(2pt) K_4 #h(2pt) square.stroked.tiny #h(2pt) K_4$ (the Cartesian product of three complete graphs on the four nucleotide states; equivalently, $K_4^(square.stroked.tiny 3)$), with 64 vertices, regular degree 9, and 288 undirected edges. It contains $Q_6$ as a 192-edge subgraph: $Q_6$ contributes the 192 Hamming-distance-1 edges (single-bit changes within one 2-bit nucleotide), while the remaining 96 within-nucleotide diagonal (Hamming-distance-2) edges complete $H(3,4)$. To address the concern that $Q_6$ misses approximately one-third of single-nucleotide mutations, we introduce the weighted score $F_rho$ that interpolates between pure $Q_6$ ($rho = 0$) and full $H(3,4)$ ($rho = 1$); see Section 2.3. Throughout this paper $H(3,4)$ denotes this nucleotide-level mutation graph, which is encoding-independent (every two-bit bijection from ${A,C,G,U}$ to ${0,1}^2$ yields the same $H(3,4)$); $Q_6$ is encoding-dependent, since the partition into Hamming-1 vs Hamming-2 edges depends on the bijection.

There are $4! = 24$ such bijections. All encoding-dependent results are tested across all 24; encoding-invariant properties (such as Serine's disconnection at $epsilon = 1$, where $epsilon$ denotes the Hamming distance threshold in $"GF"(2)^6$) are noted as such. Coloring optimality is significant ($p < 0.05$) under every encoding (full sweep in the Supplementary Material).

== Edge-mismatch objective function <sec:objective>

The genetic code assigns each vertex $v in Q_6$ a label $c(v) in cal(A)$, where $cal(A)$ comprises the 20 amino acids and the stop signal. The edge-mismatch score is

$ F(c) = sum_({v,w}: d(v,w) = 1) Delta(c(v), c(w)) $ <eq:mismatch>

where the sum ranges over all 192 edges of $Q_6$ (pairs of vertices at Hamming distance 1), and $Delta$ is a physicochemical distance between amino acids. We test four distinct distance metrics: the #cite(<grantham1974>, form: "prose") composite distance (composition, polarity, volume; range 5--215), the #cite(<miyata1979>, form: "prose") normalized Euclidean distance (polarity and volume only; range 0.06--5.13), the Woese polar requirement absolute difference (range 0--8.2; @woese1966 @woese1973; used by #cite(<freeland1998>, form: "prose")), and the #cite(<kyte1982>, form: "prose") hydropathy absolute difference (range 0--9.0; used by #cite(<haig1991>, form: "prose")). Synonymous edges contribute 0; edges involving a stop codon receive a fixed penalty scaled proportionally to each metric's maximum. A lower $F$ indicates a more error-minimizing code.

Stop codons are held fixed across all null models (Section 2.3) because their assignment is constrained by release-factor recognition geometry (eRF1 binding UAR/UGA) rather than tRNA decoding, and permuting stops would conflate two evolutionarily separate optimization problems. The stop-codon contribution is therefore a constant offset; sensitivity analysis across penalty values (0, 150, 215, 300) confirms this is immaterial to the ranking.

== Null models <sec:nullmodels>

=== Block-preserving null (Freeland--Hurst) <sec:null-fh>

Following #cite(<freeland1998>, form: "prose"), we group the 64 codons into 16 blocks of 4, defined by shared first-two-base prefix. Each block's internal pattern of amino acid assignments is preserved, but the mapping of patterns to blocks is permuted uniformly at random. Blocks containing stop codons are held fixed. This null preserves the synonymous codon contiguity (wobble degeneracy) inherent in the genetic code, providing a stringent test: the observed code must beat random codes that share its degeneracy architecture.

For the standard code, we drew $n = $10,000 null samples (seed 135325) and computed the conservative $p$-value as $(k + 1)/(n + 1)$, where $k$ is the number of null scores below the observed $F$.

=== Per-table null <sec:null-pertable>

The NCBI Genetic Codes registry (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi; gc.prt v4.6, retrieved 2026-04-25) currently lists 27 translation tables: codes 1--6, 9--16, and 21--33 (codes 7, 8, and 17--20 are deprecated and absent from the current registry). All 27 are analyzed in this work. Two pairs of tables share identical sense-codon mappings and differ only in their start-codon assignments: tables 1 (Standard) and 11 (Bacterial / Archaeal / Plant Plastid), and tables 27 (Karyorelict Nuclear) and 28 (Condylostoma Nuclear). Both pairs are retained as separate entries to match NCBI numbering, but produce identical results for analyses that depend only on codon$arrow.r$amino-acid mappings; the 27 NCBI tables thus correspond to 25 distinct sense-codon colorings. Each table was tested against its own block-preserving null ($n = $10,000 per table, common-seed design with seeds = base seed + table ID). $P$-values were corrected for multiple comparisons using the #cite(<benjamini1995>, form: "prose") false discovery rate (FDR) procedure. For NCBI tables with dual-function stop/sense codons (e.g., table 27 lists UGA as "Stop or Trp"; table 28 lists UAA, UAG as "Gln or Stop" and UGA as "Trp or Stop"), the primary codon$arrow.r$amino-acid analyses use the amino-acid label given in the NCBI AAs row of the gc.prt definition; stop functionality is handled only in analyses explicitly involving stop labels.

=== Graph family $G_rho$ and rho sweep <sec:null-rho>

We define a family of mutation graphs $G_rho$ parameterized by $rho in [0,1]$, where $G_0 = Q_6$ (192 Hamming-1 edges) and $G_1 = H(3,4)$ (all 288 single-nucleotide substitution edges). The weighted mismatch score on $G_rho$ is

$ F_rho(c) = F_("Hamming-1")(c) + rho dot.c F_("diagonal")(c) $ <eq:weighted>

where $F_("diagonal")$ sums over the 96 within-nucleotide distance-2 edges. This formulation treats $Q_6$ not as the biological object but as one endpoint of a continuous interpolation toward the full mutation graph. We evaluated five values of $rho in {0, 0.25, 0.5, 0.75, 1}$ with $n = $10,000 block-preserving null samples per value (common-seed design). Conservative p-values, computed as $(k+1)/(n+1)$ where $k$ is the number of null samples scoring below the observed code: $p_(rho=0) = 0.0061$ ($k = 60$), $p_(rho=0.25) = 0.0023$ ($k = 22$), $p_(rho=0.5) = 0.0007$ ($k = 6$), $p_(rho=0.75) = 0.0003$ ($k = 2$), $p_(rho=1) = 0.0003$ ($k = 2$). The last two are near the Monte Carlo resolution limit ($1/(n+1) approx 10^(-4)$ at $n = $10,000); the qualitative conclusion---that the observed code lies in the extreme low-cost tail at every $rho$---does not depend on $n$.

=== Table-preserving permutation null (topology avoidance) <sec:null-topo>

To test whether natural codon reassignment events avoid disrupting amino acid connectivity, we define a single primary candidate universe of single-codon reassignments and treat all reasonable alternatives as a sensitivity analysis. From the standard code $C$, the primary candidate set $cal(M)(C)$ pairs each of the 64 codons with each of the 20 alternative labels drawn from $cal(A)_(20) union {"Stop"}$:

$ cal(M)(C) = { (x, y) : x in "codons", y in cal(A)_(20) union {"Stop"}, y eq.not C(x) }, quad abs(cal(M)(C)) = 64 times 20 = 1280. $ <eq:candidate>

This *21-label, identity-excluded* universe (denoted U1 in Supplement §S5) gives every codon exactly 20 alternative labels, admits biologically attested stop-codon reassignment, and excludes identity moves ($y = C(x)$), which contribute no signal. Two strict variants---the *amino-acid-only, identity-excluded* universe with $abs(cal(M)) = 61 times 19 + 3 times 20 = $1,219 (U2) and the *stop-inclusive with no-ops* universe with $abs(cal(M)) = 64 times 21 = $1,344 (U4)---are reported as sensitivity analyses in Supplement §S5; topology-avoidance results are qualitatively identical under all three.

For each amino acid $a$, let $G_a^1$ denote the subgraph induced on codons assigned to $a$ with edges between codons at binary Hamming distance $lt.eq 1$. A reassignment is _topology-breaking_ if it increases the number of connected components of $G_a^1$ for any amino acid $a$ (the *increase-in-components* definition $Delta beta_0 > 0$ used by the conditional logit; the *new-disconnection-in-previously-connected-family* definition is reported in parallel as a sensitivity check; see Supplement §S3). From the 27 NCBI translation tables analyzed in this work (Section 2.3.2), we catalogued the table-specific reassignment events relative to the standard code; de-duplicating by (codon, target amino acid) tuple yields the unique-event set on which the hypergeometric and table-preserving permutation tests operate. The conditional logit model (Section 2.3.5) operates on the full table-specific event-step list, which preserves recurrent reassignments across independent lineages.

Two tests were applied: (i) a hypergeometric test treating the observed events as a sample from the finite candidate landscape ($N = $1,280, $K = 846$ topology-breaking under the primary cell---$H(3,4)$ adjacency with $Delta beta_0 > 0$ definition---against $n = 28$ observed, $x = 6$ topology-breaking observed); and (ii) a table-preserving permutation test ($n = $10,000, seed 135325) that permutes target amino acids among a table's reassigned codons, preserving within-table codon and target structure to address phylogenetic non-independence. Risk ratios with 95% log-normal confidence intervals were computed as $"RR" = (x\/n) \/ (K\/N)$. The full $2 times 2$ definition $times$ adjacency audit ($Q_6$ vs $H(3,4)$ adjacency $times$ new-disconnection vs $Delta beta_0 > 0$ definitions) is given in Supplement §S3; risk ratios fall in 0.28--0.33 across all four cells.

=== Conditional logit model of reassignment choice <sec:null-condlogit>

To test whether topology avoidance contributes independent explanatory power beyond physicochemical optimization, we fit event-level conditional logit (discrete-choice) models to natural reassignment events. Each observed reassignment is treated as a choice from the set of all $approx $1,280 possible single-codon reassignments available at the current code state. For a table with $k$ reassignment events, the code state evolves sequentially: at each step, the model assigns a probability to each candidate move based on a linear utility score, and the observed move is compared against all alternatives.

For each candidate move $m$ from code state $C$, we computed three features: (i) $Delta_"phys"$, the change in local Grantham mismatch cost summed over Hamming-1 edges incident to the reassigned codon; (ii) $Delta_"topo"$, the total increase in connected components across all amino acid codon graphs at $epsilon = 1$; and (iii) $Delta_"tRNA"$, the Hamming distance from the reassigned codon to the nearest codon already encoding the target amino acid (a heuristic proxy for tRNA repertoire disruption). The conditional logit probability of observing move $m^*$ is

$ P(m^* | cal(N)(C)) = frac(exp(bold(w)^top bold(x)_(m^*)), sum_(m in cal(N)(C)) exp(bold(w)^top bold(x)_m)) $ <eq:condlogit>

where $bold(w)$ is the weight vector and $bold(x)_m$ is the feature vector for move $m$. Features were $z$-scored across all candidates for numerical stability.

Since the temporal ordering of reassignment events within a table is unknown, we marginalized the likelihood over all $k!$ orderings for tables with $k > 1$ changes:

$ L_"table" = frac(1, k!) sum_(sigma in S_k) product_(s=1)^(k) P(m_(sigma(s))^* | cal(N)(C_(sigma, s))) $ <eq:orderavg>

For tables with $k lt.eq 6$ events (including the largest, yeast mitochondrial, $k = 6$), we enumerate all $k!$ orderings exactly (up to 720 for $k=6$); for the rare cases of $k > 6$, we sample 720 random orderings with the seeded RNG. The total likelihood is the product across all #cl.n_tables tables with reassignment events (#cl.total_events total event-steps).

Four nested models were compared under $Q_6$ topology: M1 (physicochemistry only, $w_"topo" = w_"tRNA" = 0$), M2 (topology only, $w_"phys" = w_"tRNA" = 0$), M3 (physicochemistry + $Q_6$ topology, $w_"tRNA" = 0$), and M4 (all three features). Two additional verification variants were fit using the encoding-independent $H(3,4)$ topology feature in place of $Q_6$: M2#sub[H(3,4)] (topology only, $H(3,4)$) and M3#sub[H(3,4)] (physicochemistry + $H(3,4)$ topology). The $H(3,4)$ variants test whether M3 dominance is robust to the choice of topology graph; the $Delta_"topo,H(3,4)" = sum_a (beta_0(G_a^"after,H(3,4)") - beta_0(G_a^"before,H(3,4)"))$ feature replaces the $Q_6$ component-count change. Weights were estimated by maximum likelihood (scipy.optimize, Nelder--Mead with L-BFGS-B refinement). Model comparison used the corrected Akaike Information Criterion (AICc) and likelihood-ratio tests for nested pairs.

=== Fisher--Stouffer test for tRNA enrichment <sec:null-trna>

For each variant-code organism paired with a phylogenetically proximate standard-code control, we built a $2 times 2$ contingency table comparing the proportion of tRNA genes encoding the reassigned amino acid versus all other amino acids. The sampling model treats tRNA gene counts as draws from a hypergeometric distribution conditional on the row and column marginals (total tRNAs per organism, total tRNAs for the focal amino acid), which is appropriate when the question is whether a specific amino acid's share of the tRNA repertoire differs between two genomes. Fisher's exact test (one-sided, alternative = "greater") was applied per pairing, and $p$-values were combined via Stouffer's $Z$ method. To address non-independence from shared control organisms, we constructed a conflict graph (edges connect pairings sharing an organism) and enumerated all maximal independent sets (MIS) via the Bron--Kerbosch algorithm with pivoting. For each MIS with $>= 2$ members, we computed the Stouffer combined $p$-value and report the worst-case (maximum $p$) across all MIS as the primary conservative test.

tRNA gene counts for 18 organisms were obtained by running tRNAscan-SE 2.0.12 @chan2019 with Infernal 1.1.4 on NCBI genome assemblies (eukaryotic mode for ciliates and yeast; bacterial mode for _Mycoplasmoides_). The verified set comprises 5 organisms with translation table 6 (UAA/UAG$arrow.r$Gln; ciliates: _Tetrahymena thermophila_, _Paramecium tetraurelia_, _Oxytricha trifallax_, _Pseudocohnilembus persalinus_, _Halteria grandinella_), 6 with table 10 (UGA$arrow.r$Cys; _Euplotes_ species), 1 with table 15 (_Blepharisma stoltei_; the genome tested here reads UGA$arrow.r$Trp via a dedicated suppressor tRNA-Trp(UCA) #cite(<singh2023>),#" "despite the legacy NCBI table-15 UAG$arrow.r$Gln assignment, see Table 8 footnote §), 1 with table 31 (multiple stops reassigned; _Blastocrithidia nonstop_), 2 with table 4 (UGA$arrow.r$Trp; _Mycoplasmoides genitalium_, _M. pneumoniae_), and 3 standard-code controls (_Stentor coeruleus_, _Ichthyophthirius multifiliis_, _Fabrea salina_). One additional standard-code organism (_Saccharomyces cerevisiae_) was sourced from GtRNAdb (https://gtrnadb.ucsc.edu) for context but is not included in the primary 24-pairing analysis. The dataset comprises 24 pairings across 5 variant genetic codes (NCBI translation tables 4, 6, 10, 15, and 31) and 3 tRNAscan-SE--verified standard-code controls.

== Retrospective cross-study reanalysis of synthetic recoding outcomes <sec:methods-codonsafe>

To test whether GF(2)#super[6] topology predicts codon recoding outcomes in synthetic biology experiments, we performed a retrospective cross-study reanalysis across nine published genome recoding datasets, with quantitative analysis of eight (the ninth, @ding2024, is included for cross-kingdom scope without quantitative extraction). For each codon substitution reported in these studies, we computed three topology features under the default encoding ($C |-> (0,0)$, $U |-> (0,1)$, $A |-> (1,0)$, $G |-> (1,1)$):

+ *Boundary crossing* ($epsilon = 1$): whether the source and target codons lie in different connected components of the amino acid's codon graph at Hamming distance 1. Generalized beyond serine to all amino acids with potential disconnections (including leucine, arginine).

+ *Local mismatch change* ($Delta F_"local"$): the difference in total Grantham distance summed over all Hamming-1 neighbors between the target and source codon positions: $Delta F_"local" = sum_(n in N(t)) Delta(c(t), c(n)) - sum_(n in N(s)) Delta(c(s), c(n))$, where $N(v)$ denotes the 6 Hamming-1 neighbors of vertex $v$ in $Q_6$.

+ *Hamming distance*: the number of bit positions differing between source and target vectors in GF(2)#super[6].

Codon positions were extracted from published GenBank genome files by CDS-level comparison using Biopython @cock2009, with codon frame offsets (`/codon_start`) handled explicitly. CDS features were matched between parent and recoded genomes by locus tag. The primary analysis used the Syn57 dataset @robertson2025syn57, which provides a within-study contrast: 37,146 serine recodings cross the UCN$arrow.l.r$AGY family boundary, while 22,859 alanine recodings remain within the GCN family. Syn57 RNA-seq differential expression data (Data file S10) were used as the outcome measure. Design-to-final genome deviations were identified by comparing the vS33A7 design genome (Data file S2) against the verified Syn57 genome (Data file S8), filtering genes with $>$10 CDS-level differences as structural rearrangements (sensitivity analysis across thresholds 3--50 reported in the Supplementary Material). The #cite(<ostrov2016>, form: "prose") segment viability data (Table S4) were analyzed as a case-control comparing segments with and without lethal design exceptions; Bonferroni correction was applied for three simultaneous tests.

All cross-study reanalysis code is released in the `codon_topo.analysis.codonsafe` subpackage of the `codontopo` repository (https://github.com/biostochastics/codontopo, version #stats._version) and can be run via `codon-topo codonsafe`. Raw data files and download provenance are documented in `data/codonsafe/DATA_MANIFEST.md`.

== Synthetic-biology feasibility score (visualization-only) <sec:methods-feasibility>

For Figure 5A we use a heuristic feasibility score $S(m) in [0,1]$ that combines four-fold filtration preservation, a smooth function of local Grantham mismatch change, and Hamming accessibility to the target amino acid. The score is purely a visualization aid for delineating high- versus low-feasibility regions of the candidate single-codon reassignment landscape; it is not used in any inferential test in this paper. The exact functional form, weights, and implementation are given in Supplement §S12 and `src/codon_topo/analysis/synbio_feasibility.py` of the public repository.

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
    [Per-table optimality (#pt.n_significant_bh/#pt.n_tables)], [FH per-table + BH-FDR], [$< 0.05$],
    [$rho$-robustness ($rho in [0,1]$)], [FH weighted edges], [0.006],
    [Topology avoidance ($Q_6$ + $H(3,4)$)], [Table-pres. perm. + clade excl.], [$10^(-4)$],
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

Under the Freeland--Hurst block-preserving null with $n = #fmtk(stats.coloring.n_samples)$ permutations, the standard genetic code is significantly low-cost across all four physicochemical distance metrics examined (@tbl:metrics). For Grantham distance, the observed score is $F = #fmtk(g.observed)$ versus a null mean of $#fmtk(g.null_mean) plus.minus #str(int(calc.round(g.null_std, digits: 0)))$ ($z = #str(calc.round(g.z, digits: 2))$, quantile #str(calc.round(g.quantile, digits: 2))%, $p = #str(calc.round(g.p, digits: 3))$). The same pattern holds for Miyata distance ($F = #str(calc.round(mi.observed, digits: 1))$ vs $#str(calc.round(mi.null_mean, digits: 1)) plus.minus #str(calc.round(mi.null_std, digits: 1))$, $z = #str(calc.round(mi.z, digits: 2))$, $p < 0.001$), Woese polar requirement ($F = #str(calc.round(pr.observed, digits: 1))$ vs $#str(calc.round(pr.null_mean, digits: 1)) plus.minus #str(calc.round(pr.null_std, digits: 1))$, $z = #str(calc.round(pr.z, digits: 2))$, $p = #str(calc.round(pr.p, digits: 3))$), and Kyte--Doolittle hydropathy ($F = #str(calc.round(kd.observed, digits: 1))$ vs $#str(calc.round(kd.null_mean, digits: 1)) plus.minus #str(calc.round(kd.null_std, digits: 1))$, $z = #str(calc.round(kd.z, digits: 2))$, $p = #str(calc.round(kd.p, digits: 3))$). Relative to null expectations, the observed code lowers the mismatch score by #str(calc.round(g.improvement_pct, digits: 1))%, #str(calc.round(mi.improvement_pct, digits: 1))%, #str(calc.round(pr.improvement_pct, digits: 1))%, and #str(calc.round(kd.improvement_pct, digits: 1))%, respectively. Beta-posterior credible intervals for the empirical $p$-values remain entirely below 0.01 for all four metrics (computed from the conservative $(k+1)/(n+1)$ estimator with $"Beta"(k+1, n-k+1)$ posterior; per-metric intervals are listed in `output/coloring_optimality.json` under `beta_credible_intervals`). Since #cite(<freeland1998>, form: "prose") established optimality using polar requirement alone, the cross-metric concordance demonstrates that error-minimization is robust across multiple distinct physicochemical parameterizations, not an artifact of any single distance metric. Under a weaker degeneracy-only null that preserves only codon family sizes without maintaining block contiguity, the signal strengthens substantially ($z > 9$, $p < 10^(-4)$ for all metrics), confirming that the block-preserving null provides a stringent test.

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
#let pos2_pct = calc.round(dec.position_fractions.pos2 * 100, digits: 1)
#let pos1_pct = calc.round(dec.position_fractions.pos1 * 100, digits: 1)
#let pos3_pct = calc.round(dec.position_fractions.pos3_wobble * 100, digits: 1)

Decomposing $F$ by nucleotide position (@fig:translational, panel B) reveals that the second codon position contributes #str(pos2_pct)% of the total mismatch, the first position #str(pos1_pct)%, and the wobble position only #str(pos3_pct)%. This gradient mirrors the biochemical hierarchy of mutational impact: second-position changes are the most physicochemically disruptive, first-position changes are intermediate, and wobble-position changes are largely synonymous @woese1965 @freeland1998.

// Figure 1 — Core topology (persistent homology + disconnection catalogue)
#figure(
  image("figures/Fig1_core_topology.png"),
  caption: [
    Core topology of the genetic code in $"GF"(2)^6$. *(A)* Persistent homology: connected components ($beta_0$) of amino acid codon graphs as a function of Hamming distance threshold $epsilon$. Serine (bold) is the only amino acid disconnected at $epsilon = 1$ across all 27 NCBI translation tables, reconnecting at $epsilon = 4$. *(B)* Disconnection catalogue across all translation tables. Each tile denotes an amino acid disconnected at $epsilon = 1$; the number is the reconnection $epsilon$. Serine is universally disconnected; Leu, Thr, and Ala appear as variant-code disconnections.
  ],
) <fig:topology>

// Figure 2 — Coloring optimality (null + per-table + rho)
#figure(
  image("figures/Fig2_coloring_optimality.png"),
  caption: [
    Coloring optimality of the genetic code. *(A)* Null distribution of Grantham edge-mismatch scores $F$ under the Freeland--Hurst block-preserving null ($n = #fmtk(stats.coloring.n_samples)$). The observed standard code (red line, $F = #fmtk(g.observed)$) falls at the #str(calc.round(g.quantile, digits: 2))th percentile ($p = #str(calc.round(g.p, digits: 3))$). *(B)* $p$-value across diagonal-edge weight $rho$ from 0 ($Q_6$) to 1 (full $H(3,4)$); all values below 0.01, with optimality strengthening monotonically. *(C)* Per-table quantile for each of #pt.n_tables NCBI translation tables under their own block-preserving null ($n = $10,000); only table 3 (yeast mito.) exceeds the 5% threshold.
  ],
) <fig:coloring>

== Robustness across diagonal-edge weight ρ <sec:res-rho>

// --- dynamic rho sweep ---
#let rho0 = rho_data.per_rho.at(0)
#let rho4 = rho_data.per_rho.at(rho_data.per_rho.len() - 1)

The weighted score $F_rho$ (@eq:weighted) remains significantly optimal across all values of $rho$ from 0 to 1 (@tbl:rho). At $rho = 0$ (pure $Q_6$), $p = #str(calc.round(rho0.p_value, digits: 3))$ ($z = #str(calc.round(rho0.at("z", default: rho0.at("effect_size_z", default: 0)), digits: 2))$); at $rho = 1$ (full $H(3,4)$, all 288 single-nucleotide edges equally weighted), $p < 0.001$ ($z = #str(calc.round(rho4.at("z", default: rho4.at("effect_size_z", default: 0)), digits: 2))$). No value of $rho$ yields $p > #str(calc.round(rho0.p_value, digits: 3))$. The optimality signal strengthens monotonically as diagonal edges are included, addressing the concern that the $Q_6$ representation ignores approximately one-third of single-nucleotide mutations.

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

Of the #pt.n_tables NCBI translation tables, #pt.n_significant_bh remain in the top 5% of their own block-preserving null after Benjamini--Hochberg FDR correction (@fig:coloring, panel C). The mean quantile across all tables is #str(calc.round(pt.mean_quantile, digits: 1))%. The single marginal exception is translation table #pt.marginal_table (yeast mitochondrial code, 6 codon reassignments; $p_"BH" = $#str(pt.marginal_p_bh)).

The headline number conflates two distinct evidentiary regimes, and we disaggregate them explicitly. Variant tables that differ from the standard code by $gt.eq 3$ codon reassignments (henceforth "informative-distance" tables) sit far enough from the standard code that the per-table block-preserving null contains many permutations that are not close to the standard code, so the test genuinely probes whether the variant is independently low-cost. Variant tables with $lt.eq 2$ reassignments (henceforth "near-standard" tables) sit so close to the standard code that the per-table null distribution is dominated by permutations near the standard code, and the test cannot reliably distinguish "this variant is independently optimal" from "this variant is close to the standard code, which is optimal." A standard-code-proximity audit (Supplement §S8) makes the structural rationale explicit: every variant in the registry sits at the extreme low-$d_H$ tail of its block-preserving null distribution (every null draw has $d_H gt.eq 30$ from the standard code, while every variant has $d_H lt.eq 6$), so the per-table $p$-value carries an unavoidable proximity-to-standard component for all tables. The "informative-distance" vs "near-standard" disaggregation rests on the absolute $d_H gt.eq 3$ threshold rather than on a conditional null comparison; the substantive claim is that variants whose absolute $d_H$ is large enough to admit a *non-trivial* range of permutations distant from standard remain low-cost on that wider range, while near-standard variants do not generate enough such permutations within the block-preserving null to test independence of proximity.

Of the #pt.informative_total informative-distance tables, #pt.informative_significant retain top-5% placement after BH--FDR correction (the lone marginal case is yeast mitochondrial). Of the #pt.near_standard_total near-standard tables, #pt.near_standard_significant are formally significant under the same correction, but those results should be read as *confirmatory of standard-code geometry* rather than as evidence for independent variant-specific optimization. The substantive evolutionary claim is therefore: of the variant codes that differ enough from the standard code for the test to be informative, all but one (the most extensively reassigned, yeast mitochondrial) preserve near-optimal physicochemical placement. Even the most reassigned code (table 3, 6 codons reassigned) sits at the 7.4th percentile of its block-preserving null, marginally above the 5% threshold but still in the lower tail.

== Topology avoidance in natural codon reassignments <sec:res-topo>

We adopt the encoding-independent $H(3,4)$ Hamming graph as the primary adjacency for the topology-avoidance test, with the $Q_6$ subgraph reported as a coordinate-dependent decomposition (motivated below). Under $H(3,4)$, of the #fmtk(tk43.possible_total) candidate single-codon relabelings (Methods §2.3.4), #fmtk(tk43.possible_breaks) (#str(calc.round(tk43.rate_possible * 100, digits: 1))%) increase the connected-component count of some amino acid's codon graph relative to the standard code ($Delta beta_0 > 0$). Only #tk43.observed_breaks of #tk43.observed_total observed natural reassignment events (#str(calc.round(tk43.rate_observed * 100, digits: 1))%) do so, yielding a #str(calc.round(tk43.depletion_fold, digits: 1))-fold depletion ($"RR" = #str(calc.round(tk43.risk_ratio, digits: 2))$, 95% CI $[#str(calc.round(tk43.risk_ratio_ci_95.at(0), digits: 2)), #str(calc.round(tk43.risk_ratio_ci_95.at(1), digits: 2))]$; hypergeometric $p = #sci(tk43.hypergeom_p)$; table-preserving permutation $p lt.eq 10^(-4)$).

Under $Q_6$ (Hamming-1 adjacency in the default $"GF"(2)^6$ encoding), the corresponding numbers are #fmtk(tq6.possible_breaks) of #fmtk(tq6.possible_total) candidates (#str(calc.round(tq6.rate_possible * 100, digits: 1))%) creating new amino-acid disconnections, against only #tq6.observed_breaks of #tq6.observed_total observed events (#str(calc.round(tq6.rate_observed * 100, digits: 1))%); $Q_6$ depletion is #str(calc.round(tq6.rate_possible / calc.max(tq6.rate_observed, 0.001), digits: 1))-fold (hypergeometric $p = #sci(tq6.hypergeom_p)$; permutation $p lt.eq 10^(-4)$). The Q_6 candidate-rate is somewhat higher than the H(3,4) rate because some $Q_6$-disconnected pairs become connected when within-nucleotide Hamming-2 edges are admitted. Both adjacencies yield highly significant depletion under both topology-breaking definitions (new disconnection in a previously connected family; $Delta beta_0 > 0$ increase in components), with risk ratios in the range 0.28--0.33 across the four cells of the $2 times 2$ definition $times$ adjacency audit (Supplement §S3). The avoidance of topology-disrupting reassignment trajectories is therefore not specific to either the hypercube subgraph or to either definition. We give precedence to the $H(3,4)$ result because $Q_6$ is encoding-dependent: across all 24 base-to-bit bijections, the $Q_6$ candidate-landscape rate ranges from 36% (8 encodings) to 73% (default encoding), with median hypergeometric $p = 6 times 10^(-6)$ but $p > 0.5$ in 8 of 24 encodings (Supplement §S4). $H(3,4)$, by contrast, depends only on nucleotide identity and is robust by construction.

Following #cite(<sengupta2007>, form: "prose"), we performed clade-exclusion sensitivity analysis, iteratively removing each major taxonomic group (all ciliates, all metazoan mitochondria, all CUG-clade yeasts, etc.) and retesting. The depletion is highly significant in every regime ($p < 10^(-5)$), but the structure of the result requires cautious interpretation: the topology-breaking events themselves are not broadly distributed across clades. Four of the 6 breakers under $H(3,4)$ come from a single clade (yeast mitochondrial, NCBI translation table 3; all 4 of its de-duplicated reassignments are breakers); excluding it leaves $2 / 24 = 8.3$% breakers, with hypergeometric $p$ moving from #sci(tq6.hypergeom_p) to $3.6 times 10^(-11)$. The clade-exclusion robustness is thus a denominator effect---removing clades that contributed zero or one breaker each leaves the breakage rate essentially unchanged---rather than evidence of repeated, independently arising topology-breaking events spread across many lineages. The cross-lineage signal is therefore the *avoidance* of topology-breaking moves, supported by 22 of 28 events distributed across many lineages; the cross-clade structure of the breaker set itself is underpowered.

NCBI translation table 3 (yeast mitochondrial) appears as the marginal case in three diagnostics: per-table optimality (§3.3), the lowest-ranking observed move under the M3 conditional logit (§3.5), and the topology-breaker concentration here. These three observations share a $k$-magnitude effect by construction (table 3 has the most reassignments in the registry, and is the only variant code requiring acquisition of a novel tRNA from a different parent---tRNA#super[Thr] derived from tRNA#super[His]; @su2011), so we do not interpret their convergence as independent confirmation of any specific selectionist mechanism. The pattern is consistent with topology-breaking reassignments being feasible but appearing concentrated in unusually perturbed translation systems.

The topology-breaking definition matters less than the choice might suggest: all four cells of the adjacency $times$ definition audit (@tbl:topo-audit) yield hypergeometric $p < 10^(-5)$ with risk ratios in the range 0.28--0.33. The encoding-sensitivity audit is more substantive. The $H(3,4)$ result is constant across all 24 base-to-bit bijections (which all yield the same $H(3,4)$ graph), but the $Q_6$ candidate landscape varies: 8 of the 24 encodings place the $Q_6$ candidate rate near 36% (vs. 73% for the default), and under those the observed rate (21--36%) does not differ from the candidate rate (median $p = 6 times 10^(-6)$ across all 24 encodings, but $p > 0.5$ in 8 of 24). We therefore treat $H(3,4)$ as the primary topology-avoidance result and $Q_6$ as a coordinate-dependent decomposition (Supplement §S4). Denominator sensitivity (alternative candidate universes $abs(M) = 1219$, $abs(M) = 1344$) is reported in Supplement §S5 and does not change the qualitative conclusion.

// Expected counts under the null of random sampling from the candidate
// landscape: E[breaks] = n_obs * (K / N), E[preserving] = n_obs * ((N-K)/N).
#let _exp_breaks = tk43.observed_total * tk43.rate_possible
#let _exp_preserve = tk43.observed_total * (1 - tk43.rate_possible)

#figure(
  table(
    columns: (auto, auto, auto, auto),
    align: (left, right, right, right),
    inset: 6pt,
    stroke: (x, y) => if y == 0 { (bottom: 0.7pt) } else { none },
    table.header([*Quantity*], [*Observed*], [*Candidate landscape*], [*Expected (chance)*]),
    [Topology-breaking ($Delta beta_0 > 0$)],
      [#tk43.observed_breaks (#str(calc.round(tk43.rate_observed * 100, digits: 1))%)],
      [#fmtk(tk43.possible_breaks) (#str(calc.round(tk43.rate_possible * 100, digits: 1))%)],
      [#str(calc.round(_exp_breaks, digits: 1)) (#str(calc.round(tk43.rate_possible * 100, digits: 1))%)],
    [Topology-preserving ($Delta beta_0 = 0$)],
      [#{ tk43.observed_total - tk43.observed_breaks } (#str(calc.round((1 - tk43.rate_observed) * 100, digits: 1))%)],
      [#{ tk43.possible_total - tk43.possible_breaks } (#str(calc.round((1 - tk43.rate_possible) * 100, digits: 1))%)],
      [#str(calc.round(_exp_preserve, digits: 1)) (#str(calc.round((1 - tk43.rate_possible) * 100, digits: 1))%)],
    [Total], [#tk43.observed_total], [#fmtk(tk43.possible_total)], [#tk43.observed_total],
    [Depletion fold (cand. rate / obs. rate)], [---], [#str(calc.round(tk43.depletion_fold, digits: 1))$times$], [1.0$times$],
    [RR (95% CI)], [---],
      [#str(calc.round(tk43.risk_ratio, digits: 2)) (#str(calc.round(tk43.risk_ratio_ci_95.at(0), digits: 2))--#str(calc.round(tk43.risk_ratio_ci_95.at(1), digits: 2)))],
      [1.00 (ref.)],
    [Hypergeometric _p_ (one-sided)], [---], [#sci(tk43.hypergeom_p)], [(null reference)],
    [Permutation _p_ (table-pres.)], [---], [$lt.eq 10^(-4)$], [(null reference)],
  ),
  caption: [
    Topology avoidance in natural codon reassignments under the encoding-independent $H(3,4)$ Hamming graph (primary cell). Topology-breaking = $Delta beta_0 > 0$, an increase in the total number of connected components summed across amino acid codon graphs. The *Expected (chance)* column shows what would be observed under the null hypothesis that the 28 observed events are drawn uniformly at random from the 1,280-move candidate landscape: $E["breaks"] = 28 times K\/N approx 18.5$ versus 6 observed. RR = risk ratio (observed rate / candidate rate) with 95% log-normal CI; depletion fold = (candidate rate)/(observed rate). The remaining three cells of the $2 times 2$ definition $times$ adjacency audit ($Q_6$ vs $H(3,4)$ $times$ new-disconnection vs $Delta beta_0 > 0$) are reported in Supplement §S3 (Table 5 of this manuscript reproduces them inline), with risk ratios in the range 0.28--0.33 and hypergeometric $p < 10^(-5)$ throughout.
  ],
) <tbl:topo>

#figure(
  table(
    columns: (auto, auto, auto, auto, auto),
    align: (left, left, right, right, right),
    inset: 5pt,
    stroke: (x, y) => if y == 0 { (bottom: 0.7pt) } else { none },
    table.header([*Adjacency*], [*Definition*], [*$K\/N$*], [*Hyper. _p_*], [*RR (95% CI)*]),
    [$H(3,4)$ (primary)], [$Delta beta_0 > 0$], [846 / 1,280], [$1.3 times 10^(-6)$], [0.32 (0.16--0.66)],
    [$H(3,4)$], [new disconnection], [822 / 1,280], [$5.0 times 10^(-7)$], [0.28 (0.13--0.62)],
    [$Q_6$], [$Delta beta_0 > 0$], [963 / 1,280], [$2.2 times 10^(-8)$], [0.33 (0.17--0.63)],
    [$Q_6$], [new disconnection], [931 / 1,280], [$1.6 times 10^(-8)$], [0.29 (0.14--0.60)],
  ),
  caption: [
    Sensitivity analysis: the $2 times 2$ topology-breaking definition $times$ adjacency audit. All four cells share the same denominator (1,280 candidates) and the same observed-event total (28 de-duplicated reassignments); the observed topology-breaking count is 5--7 of 28 across cells. All four cells yield depletion in the same direction with comparable risk ratios (0.28--0.33). The $H(3,4)$ Hamming graph is encoding-independent; the $Q_6$ subgraph is encoding-dependent and 8 of 24 base-to-bit bijections give no $Q_6$ depletion at all (Supplement §S4). The main text uses the $H(3,4)$, $Delta beta_0 > 0$ cell as primary.
  ],
) <tbl:topo-audit>

// Figure 3 — Evolutionary evidence (bit-bias + depth + topology avoidance + tRNA)
#figure(
  image("figures/Fig3_evolutionary_evidence.png"),
  caption: [
    Evolutionary evidence for structural constraints on codon reassignment. *(A)* Bit-position bias: distribution of bit-flips across $"GF"(2)^6$ coordinates in natural reassignment events; dashed line = uniform expectation. *(B)* Evolutionary depth calibration: reconnection $epsilon$ vs estimated divergence age (log scale) for 4 variant-code amino acids. *(C)* Topology avoidance shown for the encoding-independent $H(3,4)$ Hamming graph (left panel, primary) and the encoding-dependent $Q_6$ subgraph (right panel, sensitivity). Under $H(3,4)$, observed natural reassignments break topology at #str(calc.round(tk43.rate_observed * 100, digits: 1))% versus #str(calc.round(tk43.rate_possible * 100, digits: 1))% of the candidate landscape (RR #str(calc.round(tk43.risk_ratio, digits: 2)), permutation $p lt.eq 10^(-4)$). *(D)* tRNA enrichment: rank of the reassigned amino acid among all 20 AAs by tRNA gene proportion in variant-code vs standard-code organism pairings; rank 1 = most enriched.
  ],
) <fig:topo>

== Explanatory modeling: topology as an independent predictor <sec:res-condlogit>

// --- dynamic conditional logit ---
#let clf = cl.model_fits
#let m3f  = clf.at("M3_phys_topo", default: (:))
#let m4f  = clf.at("M4_full", default: (:))
#let m2f  = clf.at("M2_topo", default: (:))
#let m1f  = clf.at("M1_phys", default: (:))
#let m2k  = clf.at("M2_topo_k43", default: (:))
#let m3k  = clf.at("M3_phys_topo_k43", default: (:))
#let best_aicc = m3f.aicc

To test whether topology avoidance is reducible to physicochemical optimization, we fit event-level conditional logit models to the #cl.total_events reassignment event-steps across #cl.n_tables variant-code tables, treating each observed reassignment as a choice among $approx $1,280 candidate single-codon moves (@sec:null-condlogit). Six models were fit: four sharing the encoding-dependent $Q_6$ topology feature (M1 phys, M2 topo, M3 phys+topo, M4 phys+topo+tRNA) plus two H(3,4) verification variants that replace the $Delta_"topo,Q_6"$ feature with the encoding-independent $Delta_"topo,H(3,4)"$ (M2#sub[H(3,4)] topo only, M3#sub[H(3,4)] phys+topo). The combined physicochemistry-plus-topology model under $Q_6$ topology (M3) was strongly favored over all alternatives by AICc ($"AICc" = #str(calc.round(m3f.aicc, digits: 1))$), outperforming the physicochemistry-only model (M1; $Delta"AICc" = #str(calc.round(m1f.aicc - m3f.aicc, digits: 1))$) and the topology-only model (M2; $Delta"AICc" = #str(calc.round(m2f.aicc - m3f.aicc, digits: 1))$). Crucially, the encoding-independent verification model M3#sub[H(3,4)] also strongly outperforms the physicochemistry-only baseline ($Delta"AICc"("M1" arrow.r "M3"#sub[H(3,4)]) = #str(calc.round(m1f.aicc - m3k.aicc, digits: 1))$) and the H(3,4)-topology-only baseline ($Delta"AICc"("M2"#sub[H(3,4)] arrow.r "M3"#sub[H(3,4)]) = #str(calc.round(m2k.aicc - m3k.aicc, digits: 1))$), so the conditional-logit thesis does not depend on the choice of topology graph. Adding a heuristic tRNA-complexity proxy did not improve fit (M3$arrow.r$M4: $Delta"AICc" = #str(calc.round(m4f.aicc - m3f.aicc, digits: 1))$).

// Build model comparison table dynamically, sorted by AICc, including
// both Q_6 (legacy primary) and H(3,4) verification variants so the
// reader sees all six models in one place rather than the H(3,4)
// numbers being deferred to the supplement.
#let model_rows = (
  ("M3_phys_topo",     "Phys + topo (Q_6, primary)",          m3f),
  ("M4_full",          "Phys + topo (Q_6) + tRNA proxy",      m4f),
  ("M3_phys_topo_k43", "Phys + topo (H(3,4) verification)",   m3k),
  ("M2_topo",          "Topo only (Q_6)",                     m2f),
  ("M1_phys",          "Phys only",                           m1f),
  ("M2_topo_k43",      "Topo only (H(3,4))",                  m2k),
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
      // Strip the trailing "_k43" so M3_phys_topo_k43 displays as M3 (with
      // the description column carrying the "(H(3,4) verification)" label).
      let short = r.at(0).split("_").at(0)
      (
        [#if is_best [*#short*] else [#short]],
        [#r.at(1)],
        [#f.n_params],
        [$#str(calc.round(f.log_likelihood, digits: 1))$],
        [#if is_best [*#str(calc.round(f.aicc, digits: 1))*] else [#str(calc.round(f.aicc, digits: 1))]],
        [#if is_best [*#str(delta)*] else [#str(delta)]],
      )
    }).flatten(),
  ),
  caption: [
    Event-level conditional logit model comparison. Six models share the same candidate set ($approx $1,280 single-codon moves per event-step) and the same likelihood structure marginalised over event orderings; they differ only in the feature subset used. The first four models use the encoding-dependent $Q_6$ topology feature (the legacy primary, paired with the encoding sweep in §S4); the two M-variants ending in $H(3,4)$ replace $Delta_"topo,Q_6"$ with the encoding-independent $Delta_"topo,H(3,4)"$ feature, providing the encoding-robustness verification. M3 (phys + topo) is favoured under both topology graphs; the gap between M3#sub[Q_6] and M3#sub[H(3,4)] reflects which topology feature carries more information per parameter, not whether topology is itself meaningful. $k$ = number of estimated parameters. Bold = best model by AICc.
  ],
) <tbl:condlogit>

// --- dynamic LR test results ---
#let lr_m1m3 = cl.lr_tests.at("M1_vs_M3", default: (:))
#let lr_m2m3 = cl.lr_tests.at("M2_vs_M3", default: (:))

Likelihood-ratio tests confirmed that each feature class adds substantial explanatory value to the other: adding topology to physicochemistry yields LR $= #str(calc.round(lr_m1m3.lr_statistic, digits: 1))$ ($p lt.double 10^(-10)$), and adding physicochemistry to topology yields LR $= #str(calc.round(lr_m2m3.lr_statistic, digits: 1))$ ($p lt.double 10^(-10)$). The two feature classes are only weakly associated across the full candidate landscape, indicating limited confounding.

In the best-fitting model (M3), observed natural reassignments preferentially populate moves that reduce local physicochemical mismatch ($hat(beta)_"phys" = -0.004$ per Grantham unit) and strongly avoid moves that increase amino acid codon-family disconnection ($hat(beta)_"topo" = -3.26$ per additional connected component, with the conditional-logit feature evaluated under $Q_6$ adjacency using the *increase-in-components* ($Delta beta_0 > 0$) convention). Because $Q_6$ adjacency is encoding-dependent (8 of 24 base-to-bit bijections give no $Q_6$ topology depletion at the landscape level; @sec:res-topo), we additionally fit M3 with the encoding-independent $H(3,4)$ topology feature ($Delta_"topo,H(3,4)"$ rather than $Delta_"topo,Q_6"$). The encoding-robustness comparison gives $Delta"AICc"("M1" arrow.r "M3"#sub[H(3,4)]) = #str(calc.round(cl.at("encoding_robustness", default: (:)).at("delta_aicc_M1_to_M3_k43", default: 0), digits: 1))$ and $Delta"AICc"("M2"#sub[H(3,4)] arrow.r "M3"#sub[H(3,4)]) = #str(calc.round(cl.at("encoding_robustness", default: (:)).at("delta_aicc_M2k43_to_M3k43", default: 0), digits: 1))$, both well above the conventional ΔAICc > 10 threshold treated as strong evidence in the model-comparison literature @burnham2002, and similar in magnitude to the $Q_6$ counterparts ($Delta"AICc"("M1" arrow.r "M3"#sub[Q_6]) = #str(calc.round(cl.at("encoding_robustness", default: (:)).at("delta_aicc_M1_to_M3_q6", default: 0), digits: 1))$). We note that the ΔAICc > 10 threshold was empirically calibrated on linear-regression model-comparison contexts; we use it here as a conventional reference rather than as a formally calibrated cut-off in the conditional-logit setting. The M3 dominance is therefore not an artifact of the $Q_6$ encoding choice. Under M3, observed natural reassignments rank on average at the 89.5th percentile among all candidate moves (@fig:condlogit, panel B), with the recurrent UGA$arrow.r$Trp reassignment consistently ranking above the 98th percentile. The one notable outlier is the yeast mitochondrial CUU$arrow.r$Thr reassignment (30th percentile), consistent with NCBI translation table 3 being the sole marginal exception in the per-table optimality analysis (@sec:res-pertable). The non-significance of the heuristic tRNA-distance proxy (LR = #str(calc.round(cl.lr_tests.at("M3_vs_M4", default: (lr_statistic: 0.0)).lr_statistic, digits: 2)), $p = 0.73$) speaks to the inadequacy of Hamming-distance-to-nearest-target-AA-codon as a proxy for tRNA-mediated mechanistic feasibility; it should not be interpreted as evidence against tRNA effects in general (see @sec:res-trna for direct tRNA-gene-count tests). A posterior-predictive simulation under M3 reproduced the observed topology-breaking rate (observed 0.076 vs simulated mean 0.077; posterior-predictive $p = 0.60$), supporting model calibration rather than merely in-sample AICc improvement.

Like all conditional logit models, M1--M4 assume Independence of Irrelevant Alternatives (IIA): the relative probability between any two candidate moves is unaffected by adding or removing other candidates. We use the model as an explanatory rather than predictive tool---the reported ΔAICc and likelihood-ratio statistics test whether topology adds explanatory value beyond physicochemical cost, not whether the model accurately predicts which specific reassignment will occur next; full IIA discussion appears in Supplement §S6. A separate concern is that the unrestricted candidate set ($approx $1,280 single-codon moves) admits strongly deleterious alternatives (reassigning AUG-Met, multi-codon changes implicit in single-step framing, reassignments to stop in essential codons) that selection has already removed. Models with strongly negative coefficients on $Delta_"topo"$ and $Delta_"phys"$ may therefore partly rediscover that natural reassignments are not catastrophic, inflating apparent ΔAICc gaps. To bound this concern we refit M1--M4 (and the $H(3,4)$ verification variants) on a *restricted candidate set* containing only candidates whose target amino acid is already accessible at Hamming distance $lt.eq d$ from the reassigned codon ($Delta_"tRNA" lt.eq d$, with $d = 2$ as the primary biological-plausibility cut and $d in {1, 3}$ as bracketing thresholds; Supplement §S6.1). Under the recommended $d = 2$ filter the candidate set shrinks from $approx $1,280 to $approx #str(calc.round(cl.at("restricted_candidate", default: (:)).at("by_max_trna", default: (:)).at("2", default: (:)).at("candidate_summary", default: (:)).at("candidates_mean", default: 727), digits: 0))$ candidates per choice set, and M3 still dominates with $Delta"AICc"("M1" arrow.r "M3") = #str(calc.round(cl.at("restricted_candidate", default: (:)).at("by_max_trna", default: (:)).at("2", default: (:)).at("delta_aicc", default: (:)).at("M1_to_M3", default: 0), digits: 0))$ and $Delta"AICc"("M2" arrow.r "M3") = #str(calc.round(cl.at("restricted_candidate", default: (:)).at("by_max_trna", default: (:)).at("2", default: (:)).at("delta_aicc", default: (:)).at("M2_to_M3", default: 0), digits: 0))$, both well above the conventional ΔAICc>10 reference. The qualitative explanatory claim survives even the most stringent $d = 1$ filter ($approx #str(calc.round(cl.at("restricted_candidate", default: (:)).at("by_max_trna", default: (:)).at("1", default: (:)).at("candidate_summary", default: (:)).at("candidates_mean", default: 275), digits: 0))$ candidates per choice set), where $Delta"AICc"("M1" arrow.r "M3") = #str(calc.round(cl.at("restricted_candidate", default: (:)).at("by_max_trna", default: (:)).at("1", default: (:)).at("delta_aicc", default: (:)).at("M1_to_M3", default: 0), digits: 1))$ shrinks but remains above the conventional 10 threshold and $Delta"AICc"("M2" arrow.r "M3") = #str(calc.round(cl.at("restricted_candidate", default: (:)).at("by_max_trna", default: (:)).at("1", default: (:)).at("delta_aicc", default: (:)).at("M2_to_M3", default: 0), digits: 0))$ remains large; under the looser $d = 3$ filter ($approx #str(calc.round(cl.at("restricted_candidate", default: (:)).at("by_max_trna", default: (:)).at("3", default: (:)).at("candidate_summary", default: (:)).at("candidates_mean", default: 1088), digits: 0))$ candidates), the gap recovers to $Delta"AICc"("M1" arrow.r "M3") = #str(calc.round(cl.at("restricted_candidate", default: (:)).at("by_max_trna", default: (:)).at("3", default: (:)).at("delta_aicc", default: (:)).at("M1_to_M3", default: 0), digits: 0))$. The pattern is the expected inflation diagnosis: tighter biological-plausibility cuts shrink ΔAICc gaps as expected, but the qualitative claim "topology adds explanatory value beyond physicochemistry" is robust at every threshold. We therefore read the unrestricted-set ΔAICc magnitudes as upper bounds on the topology--physicochemistry separation, with the $d = 2$ filter giving a more biologically-calibrated effect size. (We do not interpret $Delta"AICc"("M3" arrow.r "M4")$ from the restricted set: because the restriction is itself defined on $Delta_"tRNA"$, the M4 tRNA-feature distribution shifts mechanically with the filter, so the M4 comparison is informative only under the unrestricted set, where it remains uninformative.)

Conditional-logit clade-exclusion sensitivity (refitting M1--M4 with each major clade dropped, matching the regime applied to the topology-avoidance hypergeometric in #cite(<sengupta2007>, form: "prose")) is reported in Supplement §S7. Across all seven exclusion regimes (ciliates, yeast mitochondrial, CUG-clade, metazoan mitochondrial, algal mitochondrial, protist mitochondrial, and hemichordate mitochondrial), M3 is robustly favored over M1 ($Delta"AICc" gt.eq #str(calc.round(cl.at("clade_exclusion", default: (:)).at("delta_M1_M3_min", default: 0), digits: 0))$, median $#str(calc.round(cl.at("clade_exclusion", default: (:)).at("delta_M1_M3_median", default: 0), digits: 0))$, max $#str(calc.round(cl.at("clade_exclusion", default: (:)).at("delta_M1_M3_max", default: 0), digits: 0))$) and over M2 ($Delta"AICc" gt.eq #str(calc.round(cl.at("clade_exclusion", default: (:)).at("delta_M2_M3_min", default: 0), digits: 0))$, median $#str(calc.round(cl.at("clade_exclusion", default: (:)).at("delta_M2_M3_median", default: 0), digits: 0))$). Even the minimum $Delta"AICc"$(M1$arrow.r$M3) of #str(calc.round(cl.at("clade_exclusion", default: (:)).at("delta_M1_M3_min", default: 0), digits: 0)) (from excluding all metazoan mitochondrial codes, leaving 40 events) sits comfortably above the conventional ΔAICc>10 reference threshold; the qualitative conclusion---that topology adds explanatory value beyond physicochemistry---survives every clade-exclusion regime. We rely most heavily on the posterior-predictive calibration check reported above ($p = 0.60$), since it does not require importing a threshold from a different statistical regime.

// Figure 5 — Conditional logit diagnostics
#figure(
  image("figures/Fig5_condlogit.png"),
  caption: [
    Event-level explanatory modeling of natural codon reassignments. *(A)* AICc comparison of all six nested conditional logit models ($Q_6$ + $H(3,4)$ variants); lower is better. The combined physicochemistry-plus-topology model (M3) is robustly favored ($Delta"AICc" gt.eq #str(calc.round(calc.min(cl.model_fits.M1_phys.aicc - cl.model_fits.M3_phys_topo.aicc, cl.model_fits.M2_topo.aicc - cl.model_fits.M3_phys_topo.aicc), digits: 0))$ above the conventional ΔAICc>10 reference threshold). *(B)* Distribution of observed move percentile ranks under M1 (physicochemistry only, light) versus M3 (combined, dark); M3 concentrates ranks toward the top of the candidate set. Dashed line = chance expectation (50th percentile). *(C)* Likelihood-ratio tests for each feature class added to its complement; both topology and physicochemistry contribute highly significant independent information ($p < 0.001$); the heuristic tRNA proxy does not ($p = 0.73$).
  ],
) <fig:condlogit>

== tRNA enrichment for reassigned amino acids <sec:res-trna>

Organisms with variant genetic codes show elevated tRNA gene copy numbers for the reassigned amino acid relative to standard-code controls. Across 24 variant--control pairings derived from 18 tRNAscan-SE--verified genome assemblies (15 variant-code organisms across 5 genetic codes plus 3 standard-code controls), Fisher's exact test combined via Stouffer's $Z$ method yields $p = 1.7 times 10^(-7)$ ($Z = 5.10$). To address non-independence from shared controls, we enumerated all maximal independent sets (MIS) from the conflict graph via Bron--Kerbosch; both MIS (each of size 6) are significant at $p < 0.05$ (worst-case $p = 0.045$).

The most striking case is _Tetrahymena thermophila_ (NCBI translation table 6, UAA/UAG reassigned to Gln; @hanyu1986), which carries 54 glutamine tRNA genes---including 39 suppressor tRNAs reading the reassigned stop codons---compared to 3 Gln tRNAs in the standard-code ciliate _Ichthyophthirius multifiliis_ (assembly GCF_000220395.1), 11 in _Stentor coeruleus_ (assembly GCA_001970955.1), and 3 in _Fabrea salina_ (assembly GCA_022984795.1 from @zhang2022). All four counts were generated in this work by running tRNAscan-SE 2.0.12 on the listed assemblies (full tRNAscan-SE output, including Std20/SeC/Supp/Undet/Pseudo breakdowns, in Supplement §S10). The pattern extends across reassignment types. Among Gln-reassignment ciliates, _Pseudocohnilembus persalinus_ (20 Gln tRNAs including 15 suppressors) and _Halteria grandinella_ (9 Gln tRNAs including 3 suppressors; @zheng2021) represent independent lineages within Oligohymenophorea and Spirotrichea respectively. Among Cys-reassignment ciliates (NCBI translation table 10, UGA$arrow.r$Cys), six tRNAscan-SE--verified _Euplotes_ species all carry tRNA-Cys genes with the non-canonical TCA anticodon (reading UGA), with 1--4 such genes per species alongside standard GCA-anticodon Cys tRNAs.

However, the pattern is not universal. _Blastocrithidia nonstop_ (NCBI translation table 31) reassigned all three stop codons (as did _Condylostoma magnum_, where stop-codon function is context-dependent; @heaphy2016) but achieved UGA$arrow.r$Trp via anticodon stem shortening (5 bp $arrow.r$ 4 bp) of tRNA-Trp(CCA), combined with an eRF1 Ser74Gly mutation, rather than gene duplication @kachale2023. Similarly, _Mycoplasmoides_ species with UGA$arrow.r$Trp use a single tRNA-Trp with anticodon modification. These boundary cases define a three-tier mechanistic landscape: (i) tRNA gene duplication in large nuclear genomes, (ii) anticodon structural modification in streamlined genomes, and (iii) anticodon base modification in minimal genomes.

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
    Summary of tRNA enrichment tests across 24 pairings drawn from 18 tRNAscan-SE--verified genomes (15 variant-code organisms across 5 variant genetic codes; 3 standard-code controls). MIS = maximal independent set, enumerated via Bron--Kerbosch algorithm on the conflict graph. Both MIS are significant at $p < 0.05$, eliminating the concern that greedy selection biased the independent-pairings result. The MIS worst-case ($p = 0.045$) serves as a robustness bound demonstrating that the result is not driven by a cherry-picked independent subset.
  ],
) <tbl:trna-tests>

// tRNA enrichment panel is now in Fig 3D above; table-based analysis follows

#figure(
  table(
    columns: (auto, auto, auto, auto, auto, auto),
    align: (left, left, center, left, right, right),
    inset: 5pt,
    stroke: (x, y) => if y == 0 { (bottom: 0.7pt) } else if y == 16 { (bottom: 0.4pt + luma(180)) } else { none },
    table.header(
      [*Organism*], [*Code*], [*Tbl. \#*], [*Assembly*], [*tRNAs*], [*Reassigned AA*],
    ),
    // Variant-code organisms (table 6, UAA/UAG -> Gln)
    [_T. thermophila_], [Variant], [6], [GCF_000189635.1], [718], [15 Gln (54 incl. supp.)],
    [_P. tetraurelia_], [Variant], [6], [GCF_000165425.1], [216], [7 Gln (18 incl. supp.)],
    [_O. trifallax_], [Variant], [6], [GCA_000295675.1], [94], [2 Gln (8 incl. supp.)],
    [_P. persalinus_], [Variant], [6], [GCA_001447515.1], [262], [5 Gln (20 incl. supp.)],
    [_H. grandinella_], [Variant], [6], [GCA_006369765.1], [130], [6 Gln (9 incl. supp.)],
    // Variant-code organisms (table 10, UGA -> Cys)
    [_E. aediculatus_], [Variant], [10], [GCA_030463445.1], [80], [3 Cys (4 incl. UCA)],
    [_E. amieti_], [Variant], [10], [GCA_048569255.1], [120], [4 Cys (8 incl. UCA)],
    [_E. focardii_], [Variant], [10], [GCA_001880345.2], [62], [1 Cys (3 incl. UCA)],
    [_E. parawoodruffi_], [Variant], [10], [GCA_021440025.1], [149], [5 Cys (9 incl. UCA)],
    [_E. weissei_], [Variant], [10], [GCA_021440005.1], [495], [17 Cys (20 incl. UCA)],
    [_E. woodruffi_], [Variant], [10], [GCA_027382605.1], [83], [2 Cys (4 incl. UCA)],
    // Variant-code organisms (other tables)
    [_B. stoltei_#super[#sym.section]], [Variant], [15], [GCA_965603825.1], [169], [6 Trp],
    [_B. nonstop_ P57], [Variant], [31], [GCA_028554745.1], [68], [2 Trp#super[\*]],
    [_M. genitalium_], [Variant], [4], [GCA_000027325.1], [36], [1 Trp#super[#sym.dagger]],
    [_M. pneumoniae_], [Variant], [4], [GCF_910574535.1], [37], [1 Trp#super[#sym.dagger]],
    // Standard-code controls
    [_S. coeruleus_], [Standard], [1], [GCA_001970955.1], [272], [11 Gln (control)],
    [_I. multifiliis_], [Standard], [1], [GCF_000220395.1], [150], [3 Gln (control)],
    [_F. salina_], [Standard], [1], [GCA_022984795.1], [89], [3 Gln (control)],
  ),
  caption: [
    All 18 organisms verified by tRNAscan-SE 2.0.12 in this work, grouped by genetic code: 5 Gln-reassignment ciliates (table 6, UAA/UAG$arrow.r$Gln), 6 Cys-reassignment _Euplotes_ species (table 10, UGA$arrow.r$Cys), 1 Trp-reassignment ciliate (_Blepharisma stoltei_, see footnote §), 1 stop-codon-reassignment trypanosomatid (_Blastocrithidia nonstop_ P57, table 31), 2 Trp-extension bacteria (_Mycoplasmoides_, table 4, UGA$arrow.r$Trp), and 3 standard-code ciliate controls. Total counts are Infernal 1.1.4-confirmed tRNAs from tRNAscan-SE 2.0.12 (eukaryotic mode for ciliates; bacterial mode for _Mycoplasmoides_); reassigned-AA counts include suppressor tRNAs reading the reassigned codon. #super[#sym.section]_Blepharisma stoltei_ is indexed as NCBI translation table 15 (Blepharisma Nuclear Code, defined as UAG$arrow.r$Gln) by virtue of its genus assignment, but the strain-specific MAC genome (GCA_965603825.1) reads UGA$arrow.r$Trp via a dedicated suppressor tRNA-Trp(UCA) #cite(<singh2023>). The analysis here tests Trp-tRNA enrichment on that empirical reading rather than the legacy table-15 nominal code. #super[\*]Anticodon stem shortening (4-bp stem rather than canonical 5-bp). #super[#sym.dagger]Post-transcriptional modification --- single tRNA-Trp reads both UGG and UGA. _Saccharomyces cerevisiae_ counts (used as a literature-derived control for the yeast-mito Thr disconnection pairing) come from GtRNAdb @chan2016 rather than tRNAscan-SE 2.0.12 and are reported separately in Supplement §S10.3.
  ],
) <tbl:organisms>

== Exploratory observations <sec:res-exploratory>

=== Bit-position bias in codon reassignments <sec:res-bitbias>

The distribution of bit-flips across the 6 coordinates of $"GF"(2)^6$ in natural codon reassignments shows apparent positional skew under a uniform null ($chi^2 = 16.26$, $p = 0.006$, $"df" = 5$). However, this signal is substantially attenuated after de-duplication to 20 unique (codon, target amino acid) pairs ($p = 0.075$) and vanishes entirely under a codon-preserving permutation null ($p = 1.0$). The apparent bias is therefore explained by which codons are recurrently reassigned across lineages, not by a genuine positional preference in $"GF"(2)^6$.

=== Variant-code disconnection catalogue <sec:res-catalogue>

A systematic survey across all 27 NCBI translation tables identifies four lineage-collapsed variant-code amino-acid disconnections at $epsilon = 1$ in $"GF"(2)^6$ under the default encoding: threonine in the yeast mitochondrial code (translation table 3, CUN$arrow.r$Thr); leucine in the chlorophycean mitochondrial codes (translation tables 16 and 22, both with UAG$arrow.r$Leu---table 16 is the chlorophycean mitochondrial code @hayashi1996 and table 22 is the closely related _Scenedesmus obliquus_ mitochondrial code, which additionally reassigns UCA Ser$arrow.r$Stop; both produce equivalent $epsilon = 2$ Leu reconnection profiles, so they collapse to a single algal-mitochondrial event); alanine in _Pachysolen tannophilus_ nuclear code (translation table 26, CUG$arrow.r$Ala); and a tripartite serine in the _Candida_-clade alternative yeast nuclear code (translation table 12, CUG$arrow.r$Ser; @santos1999). These cases, combined with the universal serine disconnection, constitute the complete inventory of amino-acid graph disconnections at unit Hamming distance under the default encoding. A separate and weaker geometric exception---specific to filtration rather than to disconnection---arises in translation table 32 (Balanophoraceae plastid; UAG$arrow.r$Trp). Trp is 1-fold (UGG only) under the standard code; in table 32 it becomes 2-fold (UGG, UAG). The pair differs in the second nucleotide (G$arrow.l.r$A), i.e. at bit position 3 in our 0-based 6-bit indexing (the second bit of the second nucleotide), not at the wobble bit (position 5) where every standard-code 2-fold pair sits. Hamming distance 1, so the pair is connected at $epsilon = 1$ and adds no new entry to the disconnection catalogue. An empirical scan over every 2-fold amino-acid pair in all 27 NCBI tables confirms that this is the unique deviation from the bit-5 two-fold filtration: the standard code's nine 2-fold amino acids and every analogous pair introduced by the other 26 tables differ at bit 5.

=== Atchley Factor 3 and Serine convergence <sec:res-atchley>

Serine has the most extreme Atchley Factor 3 score among the 20 amino acids ($F_3 = -4.760$, 2.24 SD below the mean; @atchley2005), and it is the only amino acid disconnected at $epsilon = 1$ under every base-to-bit encoding. This convergence is unsurprising: Atchley $F_3$ is a composite of molecular size and codon diversity that partly captures codon-family structure, so the $"GF"(2)^6$ topology and $F_3$ are not fully independent views of Serine's anomaly. Both reflect Serine's disproportionate codon diversity (6 codons in two disconnected families) relative to its small physicochemical footprint; the $"GF"(2)^6$ framework provides a complementary structural account rather than independent corroboration.

== Retrospective cross-study reanalysis of synthetic genome recoding outcomes <sec:res-codonsafe>

To assess whether the structural properties identified in Sections 3.1--3.5 translate to measurable phenotypic consequences in synthetic biology, we performed a retrospective cross-study reanalysis of nine published genome recoding datasets, with quantitative analysis of eight ($>$217,000 codon-level observations; @tbl:codonsafe-datasets). The ninth dataset (@ding2024 mammalian Ser TCG recoding) is included for cross-kingdom scope but without quantitative extraction. Each codon substitution was classified by its GF(2)#super[6] topology: whether it crosses a connected-component boundary at $epsilon = 1$, the change in local physicochemical mismatch cost ($Delta F_"local"$, Grantham distance), and the Hamming distance between source and target vectors.

#figure(
  table(
    columns: (auto, auto, auto, auto, auto),
    align: (left, left, right, center, left),
    inset: 6pt,
    stroke: (x, y) => if y == 0 { (bottom: 0.7pt) } else { none },
    table.header([*Study*], [*Swap type*], [*$n$*], [*Boundary?*], [*Key result*]),
    [#cite(<robertson2025syn57>, form: "prose") Syn57], [4 Ser + 2 Ala + Stop], [60,240], [Ser=100%, Ala=0%], [Topology-fixed exposure (gene-level null)],
    [#cite(<fredens2019>, form: "prose") Syn61], [Ser (TCG/TCA$arrow.r$AGC/AGT)], [18,218], [All cross], [Informative null],
    [#cite(<ostrov2016>, form: "prose")], [7 codon types], [62,214], [Mixed], [Segment viability],
    [#cite(<napolitano2016>, form: "prose")], [Arg (AGR$arrow.r$CGN)], [12,888], [None], [SRZ covariates],
    [#cite(<grome2025ochre>, form: "prose") Ochre], [Stop (TGA$arrow.r$TAA)], [1,195], [N/A], [Growth data],
    [#cite(<nyerges2024>, form: "prose")#super[#sym.dagger.double]], [7 codon types], [62,007], [Mixed], [Multi-omics],
    [#cite(<lajoie2013>, form: "prose") C321.$Delta$A], [Stop (UAG$arrow.r$UAA)], [321], [N/A], [Baseline],
    [#cite(<frumkin2018>, form: "prose")], [Arg (CGU/CGC$arrow.r$CGG)], [60], [None], [Translation efficiency],
    [#cite(<ding2024>, form: "prose")], [Ser (TCG, mammalian)], [Variable], [Yes], [Cross-organism],
  ),
  caption: [
    Published genome recoding datasets analyzed. "Boundary?" indicates whether the synonymous codon swaps cross a connected-component boundary at $epsilon = 1$ in GF(2)#super[6]. #super[#sym.dagger]bioRxiv preprint; not peer-reviewed at time of writing. Among quantitatively analyzed datasets, Syn57 provides a *topology-fixed exposure* (every Ser swap crosses the UCN$arrow.l.r$AGY boundary; no Ala swap does) paired with a null gene-level transcriptomic test; the codon-level $Delta F_"local"$ contrast within Syn57 is design-confounded and reported descriptively (see §3.8.1.2). #cite(<ding2024>, form: "prose") is included qualitatively for cross-kingdom scope without quantitative extraction. The Syn57 row total of 60,240 codon-level observations comprises 37,146 serine recodings + 22,859 alanine recodings + 235 stop-codon recodings.
  ],
) <tbl:codonsafe-datasets>

// Figure 4 — Translational applications (synbio + decomposition + catalogue)
#figure(
  image("figures/Fig4_translational.png"),
  caption: [
    Translational applications and prediction catalogue. *(A)* Feasibility landscape of 1,280 single-codon reassignments from the standard code, colored by degeneracy filtration status. All filtration-preserving variants score $gt.eq 0.8$; all filtration-breaking variants score $lt.eq 0.75$. Dashed line: high-feasibility threshold. *(B)* Grantham mismatch score decomposition by nucleotide position: position 2 dominates (#str(pos2_pct)%), consistent with the biochemical hierarchy of mutational impact. *(C)* Prediction catalogue: 15 evaluated claims distributed across 5 workstreams (WS1--WS6), colored by verification status.
  ],
) <fig:translational>

=== A three-layer interpretation of recoding outcomes (exploratory) <sec:res-threelayer>

This subsection is exploratory and hypothesis-generating. We propose a three-layer view of how codon-space structure relates to recoding outcomes, but the partition into evolutionary, recoding-burden, and design-deviation layers was *not* prespecified before the analyses below. The layering is therefore a working synthesis, not a tested model: each contrast below is individually informative, but their integration into a single three-layer schema is interpretive and would require an independent dataset built around a sharper hypothesis to confirm.

==== Boundary crossing: validated but not predictive of transcriptomic perturbation

In the #cite(<robertson2025syn57>, form: "prose") Syn57 dataset---the only large-scale dataset with within-study variation in boundary crossing---all 37,146 serine recodings crossed the UCN$arrow.l.r$AGY disconnection boundary at $epsilon = 1$, while all 22,859 alanine recodings remained within the compact GCN family. This provides an exact within-study topology contrast. However, genes containing exclusively serine-type versus exclusively alanine-type recodings did not differ in transcriptomic perturbation (Mann--Whitney $p = 0.40$, $n_"Ser" = 129$, $n_"Ala" = 49$), and the fraction of boundary-crossing recodings per gene was uncorrelated with $|log_2 "FC"|$ ($rho = 0.012$, $p = 0.46$, $n = $3,510 genes). In a separate dataset, all 18,218 Syn61 serine recodings @fredens2019 crossed the boundary; the seven design-to-final corrections that the Syn61 authors required for cell viability (Supplementary Data 19 of #cite(<fredens2019>, form: "prose")) concentrated at specific functional units---a Ser70 site in _yaaY_, an intergenic separation between _ftsI_ and _murE_, a Ser4 site in _map_, five positions in the _yceQ_ pseudogene needed to preserve viable expression of the downstream essential gene _rne_---rather than at randomly distributed sites. Because all 18,218 forward swaps share the same two AA-identity-preserving swap subtypes (TCG$arrow.r$AGC and TCA$arrow.r$AGT), |Δphys| at the AA-identity level is structurally constrained across all sites, and the 7 corrections do not lie at extreme |Δphys| or unusual Hamming-accessibility positions relative to the 18,211 non-correction sites. We therefore interpret the corrections as protein-level functional pressure at specific residues and operon-context constraints, not a global topology-driven burden, and treat this contrast as descriptive rather than as a hypothesis-test of GF(2)#super[6] features. Together with the gene-level Syn57 null, this indicates that boundary crossing is a coarse geometric property of codon-family organization---validated as a structural distinction over evolutionary timescales but not predictive of acute transcriptomic perturbation in synthetic-biology engineering, where competing ecological pressures are absent and expression constructs are highly optimized.

==== Local mismatch geometry: a fine-scale burden axis with positive signal <sec:local-mismatch-geom>

The second layer examines whether reassignment changes the physicochemical quality of the immediate mutational neighborhood. In the Syn57 contrast, serine swaps moved codons into better local Grantham neighborhoods ($Delta F_"local" = -37.2 plus.minus 30.5$), whereas alanine swaps moved codons into worse neighborhoods at a fixed $Delta F_"local" = +19.0 plus.minus 0.0$ (zero variance). This contrast is design-confounded and we report it descriptively rather than inferentially: every Syn57 alanine swap goes to the same target codon, so all alanine $Delta_"local"$ values are identical, and the Mann--Whitney $U = 0$, $p < 10^(-16)$ comparison reflects design-level non-overlap rather than a within-group inferential test. The methodologically clean inferential test comes from a topology-fixed setting---arginine synonymous recodings @napolitano2016, where all 12,888 AGR$arrow.r$CGN swaps remain within one connected component at $epsilon = 1$---in which local mismatch change correlated with established recoding-burden covariates: $Delta F_"local"$ vs RBS deviation ($rho = -0.33$, $p < 10^(-15)$) and vs mRNA folding deviation ($rho = -0.12$, $p < 10^(-15)$). Among the four CGN targets, CGG has the lowest local mismatch cost (246 Grantham units) while CGC has the highest (421), and target-specific RBS deviation spans four orders of magnitude, suggesting that local mismatch and SRZ covariates capture partially independent dimensions of recoding burden.

==== Design deviations: accessibility-dominated, not topology-directed

To test whether engineered systems escape design constraints along topology-defined low-burden directions, we compared the Syn57 design genome (Data file S2; @robertson2025syn57) against the final verified genome (Data file S8) and identified 727 CDS-level codon differences. After filtering 7 genes with $>$10 deviations each (structural rearrangements, predominantly _bioF_ with 311 changes), 326 genuine point deviations remained (267 synonymous, 59 nonsynonymous). Among the 163 deviations that changed local mismatch, 76 moved to a better neighborhood and 87 to a worse one---no directional bias (two-sided exact binomial $p = 0.43$, Clopper--Pearson 95% CI $[0.39, 0.55]$; Wilcoxon signed-rank $p = 0.67$). This null result was stable across all filtering thresholds tested ($p > 0.5$ at cutoffs 3, 5, 10, 15, 20, and 50 deviations per gene). By contrast, deviations strongly favored mutationally proximate moves: 83% were single-bit changes (Hamming distance 1 under the default encoding; note that Hamming distance is encoding-dependent, though nucleotide edit distance is invariant), and the mean Hamming distance was 1.21. Thus, realized deviations from the Syn57 design do not preferentially descend toward lower local mismatch; instead, they follow accessibility-dominated escape routes, favoring nearby codon states regardless of neighborhood quality. The #cite(<ding2024>, form: "prose") mammalian TCG recoding independently confirms the universality of the serine disconnection across standard-code organisms, providing cross-kingdom validation.

==== Segment-level recoding burden in the Ostrov 57-codon design

As an orthogonal test, we examined segment-level outcomes from the #cite(<ostrov2016>, form: "prose") 57-codon genome design. Among 44 segments with fitness data (of 87 total), segments with more recoded codons in essential genes showed a suggestive correlation with worse doubling time ($rho = 0.34$, raw $p = 0.022$, Bonferroni-corrected $p = 0.066$), while total recoding load showed no association ($rho = -0.10$, $p = 0.53$). Problem segments (13 with lethal exceptions) trended toward higher essential-gene recoding load (73.5 vs 39.8 recoded sites, Mann--Whitney $p = 0.086$). Among 44 fully characterized segments, 17 showed spontaneous codon reversions, indicating ongoing design instability at positions where the engineered code was under selection pressure.

=== Summary of cross-study reanalysis scope

The synthetic recoding analyses, taken together, are consistent with (but do not prove) a three-layer view in which different structural descriptors may operate at different biological layers: family topology for the geometry of reassignment space, local mismatch for some aspects of recoding burden, and Hamming accessibility for the short-range routes by which engineered designs deviate or revert. Specific testable predictions follow: (i) a Syn-style recoding experiment that varies family-boundary crossing while holding local mismatch fixed should fail to predict transcriptomic burden (already partly observed in Syn57); (ii) a recoding scan that varies $Delta F_"local"$ within a single boundary regime should predict RBS-context proxies (consistent with the Napolitano arginine result); (iii) directed evolution starting from a recoded design should produce escape mutations whose Hamming distribution is shaped more by local accessibility than by mismatch optimisation (consistent with the Syn57 design-deviation analysis). Each prediction admits a direct experimental test that would either confirm or falsify the layer structure proposed here.

== Scope of the framework: rejected and null findings <sec:res-negatives>

This section consolidates results that delimit what the framework does and does not do. Four conjectural extensions previously associated with the $"GF"(2)^6$ representation are tested and rejected (§3.9.1--§3.9.3); a separate negative result on source-neighborhood burden clarifies the locus at which the topology-avoidance constraint operates (§3.9.4).

=== KRAS--Fano clinical prediction <sec:res-kras>

The conjecture that XOR ("Fano") relationships in $"GF"(2)^6$ predict co-mutation enrichment at KRAS G12 sites was tested against 1,670 mutations from MSK-IMPACT @zehir2017 and cleanly falsified ($p = 1.0$; per-variant detail in Supplement §S13). This negative result separates code-level error-minimization (which is real) from mutation-level algebraic predictions (which are not).

=== Serine distance-4 invariant <sec:res-serine>

The claim that Serine's minimum inter-family Hamming distance (UCN$dash$AGY) equals 4 under all 24 base-to-bit encodings is false. Of the 24 encodings, 16 yield minimum distance 2 and only 8 yield distance 4. The distance-4 result obtains only when both nucleotide pairs distinguishing UCN from AGY ($U tilde.op A$ and $C tilde.op G$ in the first two positions) are encoded at maximal Hamming distance. The correct encoding-invariant statement is: Serine is disconnected at $epsilon = 1$ under every encoding, and its inter-family distance ($gt.eq 2$) is the largest among the three 6-codon amino acids (Leucine and Arginine both have inter-family distance 1).

=== PSL(2,7) symmetry and holomorphic embedding <sec:res-psl>

The claim that PSL(2,7) is the fundamental symmetry group of the genetic code was pre-rejected by #cite(<antoneli2011>, form: "prose"), who showed that PSL(2,7) has no 64-dimensional irreducible representation (its irreps have dimensions 1, 3, 6, 7, 8). The claim that the coordinate-wise map $"GF"(2)^6 arrow.r CC^3$ sending base pairs to fourth roots of unity is a holomorphic embedding extending a character of $"GF"(8)^*$ is also incorrect: the domain is a finite discrete set (not a complex manifold), and the map fails the character identity $chi(x + x) = chi(x)^2$ since $i^2 = -1 eq.not 1$.

=== Source-neighborhood burden: null result with informative interpretation <sec:res-infoneg>

The source-neighborhood burden test---asking whether reassigned codons sit in worse Hamming neighborhoods with higher Grantham distance to their neighbors---yields Mann--Whitney $U = 301$, $p = 0.70$. Reassignment is therefore not driven by local escape from costly source neighborhoods. This null does not conflict with the conditional-logit finding that natural events favor candidate moves with lower $Delta_"local"$ ($beta_"phys" = -0.004$; @sec:res-condlogit): the two tests probe different quantities. $Delta_"local"$ is the change in mismatch cost induced by a candidate move (a destination-quality measure); the present test asks about the absolute pre-move source-neighborhood burden. Variant codes need not originate from unusually costly source neighborhoods (no absolute source burden), yet still favor candidate moves that lower local mismatch when a reassignment is triggered (opportunistic destination selection). This indicates that the topology-avoidance constraint (@sec:res-topo) operates at the global graph-connectivity level, not at the local per-codon level.


// ============================================================
//  4. DISCUSSION
// ============================================================
= Discussion <sec:discussion>

== An information-theoretic view of the genetic code <sec:disc-info>

The central finding of this work is that the standard genetic code minimizes the physicochemical disruption caused by single-bit errors in $"GF"(2)^6$ coordinates. This is not a new conclusion---#cite(<freeland1998>, form: "prose") established error-minimization using nucleotide-level mutation models---but the hypercube representation makes the optimality principle geometrically explicit: the code is a good _coloring_ of a structured graph, in the graph-theoretic sense that adjacent vertices (codons differing by one bit) tend to share labels (amino acids) or, when they differ, differ by small physicochemical distances.

The score decomposition (@fig:translational, panel B) shows that this optimization is concentrated at the second codon position (49.3% of total mismatch) and first position (38.2%), with the wobble position contributing only 12.5%. This gradient mirrors the biochemical hierarchy of mutational impact and is an emergent property of the code's structure rather than a parameter of the model.

The $rho$-robustness result (@fig:coloring, panel B) demonstrates that optimality is not an artifact of restricting attention to $Q_6$: when the full $H(3,4)$ mutation graph is considered ($rho = 1$), the signal strengthens. The code minimizes error not just along the hypercube edges but across the complete space of single-nucleotide substitutions. This error-minimization is complementary to, but distinct from, the finding of #cite(<itzkovitz2007>, form: "prose") that the code is also nearly optimal for carrying parallel regulatory information within protein-coding sequences---a property related to stop-codon identity rather than amino acid physicochemistry.

Whether the $"GF"(2)^6$ framework adds genuine insight beyond a direct analysis of $H(3,4)$ admits a layered answer: in two specific ways yes, in one circular way no. We list all three explicitly.

*(i) New: the $rho$-sweep.* The decomposition of $H(3,4)$ into the 192 Hamming-distance-1 edges of $Q_6$ plus 96 within-nucleotide diagonal (Hamming-2) edges enables the weighted-mismatch score $F_rho$ (@eq:weighted) to interpolate continuously between $rho = 0$ ($Q_6$ only) and $rho = 1$ (full $H(3,4)$). The monotonic strengthening of optimality as $rho arrow.r 1$ (Table 3, §3.2) is a contribution that does not exist in a pure $H(3,4)$ analysis, where the diagonal/non-diagonal partition is not coordinatized.

*(ii) New: the persistence parameter $epsilon$ for codon-family connectivity.* The Hamming-distance filtration (component count $beta_0(G_a^epsilon)$ as a function of $epsilon$) provides a natural persistence parameter for the disconnection catalogue (@fig:topology, panel A) and for the conditional-logit topology feature ($Delta beta_0$). Equivalent constructions exist on $H(3,4)$ (single-nucleotide-edit distance), but the $epsilon$ axis fits naturally inside the $Q_6$ subgraph framework and connects topology to the cross-metric optimality machinery within a single coordinate system.

*(iii) The encoding-sensitivity audit: a necessary control rather than a contribution.* The 24-encoding sweep distinguishes encoding-invariant properties (Serine disconnection, $H(3,4)$ topology depletion) from encoding-dependent ones (the Serine distance-4 claim, $Q_6$ topology depletion). It identified the $Q_6$ depletion as encoding-dependent (8 of 24 bijections give no signal; @sec:res-topo, Supplement §S4) and prompted promotion of $H(3,4)$ to primary. The audit is corrective rather than additive: it identifies which results are coordinate artifacts, but it does not extend the framework's reach. Pure $H(3,4)$ has no encoding choice and therefore needs no such audit.

The binary representation is therefore best understood as an analytical decomposition tool whose distinct contributions are the $rho$-sweep continuum and the $epsilon$-filtration; the encoding-sensitivity audit is a necessary control for representation-dependent analyses rather than a contribution from the framework. We accordingly present the encoding-independent $H(3,4)$ result as the primary topology-avoidance test and treat $Q_6$ as a coordinate-dependent decomposition useful for the $rho$-sweep continuum but not as a freestanding biological claim ($rho$ being a diagonal-edge weight, not a transition/transversion weight; Methods §2.1, Supplement §S2).

== Evolutionary preservation and topology avoidance <sec:disc-evol>

The per-table analysis (@fig:coloring) shows that #pt.n_significant_bh of #pt.n_tables NCBI translation tables maintain coloring optimality under their own block-preserving null, suggesting that codon reassignment events are constrained to preserve error-minimization. The marginal exception is the most extensively reassigned code: translation table 3 (yeast mitochondrial, 6 codon changes) falls only slightly above the 5% threshold under its own block-preserving null.

The topology avoidance result (@fig:topo) provides a mechanistic complement: natural reassignments avoid creating new amino acid disconnections at a rate far below chance expectation. Importantly, this depletion is not confined to the $Q_6$ representation: under the full $H(3,4)$ single-nucleotide mutation graph, the observed proportion of topology-breaking reassignments remains #str(calc.round(tk43.rate_observed * 100, digits: 1))% against a possible rate of #str(calc.round(tk43.rate_possible * 100, digits: 1))% ($"RR" = #str(calc.round(tk43.risk_ratio, digits: 2))$, 95% CI $[#str(calc.round(tk43.risk_ratio_ci_95.at(0), digits: 2)), #str(calc.round(tk43.risk_ratio_ci_95.at(1), digits: 2))]$, $p lt.eq 10^(-4)$). The signal magnitude is similar to the $Q_6$ result ($"RR" = #str(calc.round(tq6.risk_ratio, digits: 2))$); some $Q_6$-disconnected pairs become connected when all single-nucleotide edges are admitted, slightly lowering the possible-rate denominator, but the main qualitative conclusion survives: the connected-component structure of amino acid codon families is functionally important, regardless of whether connectivity is defined on the hypercube subgraph or on the biologically fuller single-substitution graph. The yeast mitochondrial threonine reassignment illustrates this cost: the CUN$arrow.r$Thr change required acquisition of a novel tRNA#super[Thr] derived from tRNA#super[His] via anticodon mutation @su2011, creating the topology-breaking disconnection that makes translation table 3 the sole marginal outlier in the per-table optimality analysis.

Together these analyses indicate that code evolution is constrained along two partly independent axes, not one. Physicochemical smoothness preserves *protein function under mistranslation*: a reassignment that lands in a physicochemically similar neighborhood is one whose mistranslation errors remain tolerable. Topological integrity preserves *something different*---the connectivity of an amino-acid codon family determines whether existing decoding machinery (a tRNA serving wobble-related codons) can continue to service the family without new molecular hardware. The conditional logit analysis (@tbl:condlogit) makes the separation explicit. The topology term retains major explanatory value ($Delta"AICc" = #str(calc.round(m1f.aicc - m3f.aicc, digits: 0))$ relative to a physicochemistry-only model) after accounting for local physicochemical cost, while physicochemical cost retains comparable value ($Delta"AICc" = #str(calc.round(m2f.aicc - m3f.aicc, digits: 0))$) after accounting for topology; the two feature classes are only weakly correlated across the candidate landscape ($r_s = #str(calc.round(cl.phys_topo_rho, digits: 2))$), so neither term is a redundant restatement of the other.

To rule out a coordinate artifact---given that $Q_6$ topology depletion is itself encoding-dependent---we refit M3 using the encoding-independent $H(3,4)$ topology feature ($Delta_"topo,H(3,4)" = sum_a beta_0(G_a^"after,H(3,4)") - sum_a beta_0(G_a^"before,H(3,4)")$). The topology term retains essentially the same magnitude of explanatory value ($Delta"AICc"("M1" arrow.r "M3"#sub[H(3,4)]) = #str(calc.round(cl.at("encoding_robustness", default: (:)).at("delta_aicc_M1_to_M3_k43", default: 0), digits: 0))$ vs. $Delta"AICc"("M1" arrow.r "M3"#sub[Q_6]) = #str(calc.round(cl.at("encoding_robustness", default: (:)).at("delta_aicc_M1_to_M3_q6", default: 0), digits: 0))$), confirming that the conditional-logit thesis does not depend on the choice of topology graph.

We do not interpret this as direct selection on topology as an abstract graph property; topology may still proxy aspects of decoding architecture---tRNA cross-recognition range, ribosomal A-site constraints---not explicitly modeled here. The substantive claim is structural: within the tested event landscape, topology avoidance cannot be dismissed as a by-product of physicochemical optimization or as a coordinate artifact of the binary representation. The contribution beyond single-axis treatments is the second axis itself---code evolution moves through a low-cost region *and* avoids fragmentation of the decoding substrate, with the two constraints partly orthogonal and neither reducible to the other.

== Mechanistic implications: tRNA compensation <sec:disc-trna>

The tRNA enrichment result (@fig:topo, panel D) links the geometric observation to molecular mechanism. In several variant-code lineages where codon reassignment disrupts connectivity, expanded tRNA gene repertoires for the affected amino acid are observed, consistent with compensatory gene duplication as one evolutionary route to accommodation. The extreme case of _Tetrahymena thermophila_ (54 Gln tRNAs, including 39 suppressors) illustrates the scale of genomic response required to service a split codon family.

However, the boundary cases---_Blastocrithidia nonstop_ (anticodon stem shortening; @kachale2023) and _Mycoplasmoides_ species (anticodon modification)---show that tRNA gene duplication is not the only evolutionary solution. The three-tier pattern (duplication in large genomes, structural modification in intermediate genomes, base modification in minimal genomes) suggests that genome size constrains the available mechanistic repertoire for codon reassignment. This pattern is now supported by 18 tRNAscan-SE--verified genomes across 5 variant genetic codes (Tables 4, 6, 10, 15, and 31), including 6 _Euplotes_ species (UGA$arrow.r$Cys) that each carry 1--4 TCA-anticodon tRNA-Cys genes dedicated to reading UGA.

== Synthetic recoding outcomes: a three-layer interpretation <sec:disc-layers>

The retrospective cross-study reanalysis helps delimit which structural quantities in GF(2)#super[6] are functionally relevant and at which biological layer. Boundary crossing is a valid codon-space distinction---serine and alanine provide a clean contrast in which one class always crosses a family boundary and the other never does---yet this distinction did not predict RNA-seq perturbation in Syn57. We therefore do not interpret codon-family boundary crossing as a general transcriptomic burden variable. Its relevance, if any, likely lies at a different biological layer, such as decoding architecture, translational fidelity, or evolutionary accessibility, rather than steady-state expression change.

By contrast, local mismatch geometry showed a clearer positive signal. In both the Syn57 contrast and the topology-fixed arginine setting, codon changes differed systematically in the quality of the physicochemical neighborhood they entered, and these differences aligned with established recoding-burden covariates @napolitano2016. This suggests that local mismatch is best understood as a fine-scale recoding-friction variable: not a universal master predictor, but a meaningful descriptor of neighborhood quality within an already specified codon-family structure.

The design-deviation analysis further clarifies scope. Genuine deviations from the Syn57 design were overwhelmingly short-range (83% single-bit in GF(2)#super[6]), with no directional bias in mismatch change. This indicates that realized escape is governed more by mutational accessibility than by systematic optimization of local neighborhood cost. In other words, once a design is under strain, the available exits appear to be chosen primarily by mutational proximity, not by a global preference for lower mismatch.

These observations are consistent with---but, critically, do not establish---a working hypothesis that "topology" should not be treated as a single explanatory variable. Under that hypothesis, different structural descriptors would matter at different biological layers: family topology for the geometry of reassignment space (@sec:res-topo), local mismatch for some aspects of recoding burden (@sec:local-mismatch-geom), and Hamming accessibility for the short-range routes by which engineered designs deviate. Each of the three layer claims rests on a single retrospective contrast and the partition into layers was not prespecified, so the schema should be read as hypothesis-generating rather than as a tested model (see §3.8.3 and §3.8.2 for explicit predictions). That the #cite(<fredens2019>, form: "prose") Syn61 design tolerated 18,218 boundary-crossing serine swaps as a class (a genome-wide engineering success rate, not a per-codon-position viability test), while the NCBI variant codes show #str(calc.round(tq6.depletion_fold, digits: 1))-fold depletion of topology-breaking changes over evolutionary timescales, is consistent with this layered reading: boundary crossing may be more relevant to evolutionary trajectory constraints than to acute engineering costs, but a within-study experiment that varies family-boundary crossing while holding local mismatch fixed is required to test that prediction.

== Scope and falsified claims <sec:disc-honest>

Six of 15 evaluated claims were rejected, falsified, or tautological (shaded rows in @tbl:claims; full details in Supplementary Material). These failures help delimit the framework's explanatory scope. The KRAS--Fano prediction ($p = 1.0$) cleanly separates code-level error-minimization from mutation-level algebraic predictions. The Serine distance-4 invariant, falsified by 16 of 24 encodings giving distance 2, illustrates the importance of systematic encoding-sensitivity testing. The PSL(2,7) and holomorphic embedding claims were pre-rejected by existing mathematical results @antoneli2011.

== Limitations <sec:disc-limits>

Several limitations should be noted.

*Encoding choice.* The base-to-bit encoding is not unique. While coloring optimality holds across all 24 encodings, specific score values and rank orderings are encoding-dependent (Supplementary Material).

*Null-model conservatism.* The block-preserving null model, while standard in the field @freeland1998, preserves more structure than a fully random code and may understate the degree of optimality.

*tRNA enrichment fragility.* The tRNA enrichment result is robust to pairing selection (worst-case MIS $p = 0.045$) but rests on a small number of independent pairings ($n = 6$) with limited statistical power; we classify it as suggestive rather than as demonstrating a universal compensation mechanism. The MAC-genome assemblies for the more deeply sampled ciliates vary substantially in fragmentation (predicted-pseudogene fraction in Supplement §S10: 79/495 = 16% in _E. weissei_, 18/262 = 7% in _P. persalinus_, 0--5% elsewhere); removing the two genomes with $> 10$% pseudogene fraction yields 22 pairings, on which the all-pairings Stouffer combination gives $p = 0.041$, so the suggestive result is not driven by the most fragmented assemblies. A small consistency caveat: the _Saccharomyces cerevisiae_ tRNA counts used in the yeast-mitochondrial Thr disconnection pairing come from GtRNAdb @chan2016 rather than from a tRNAscan-SE 2.0.12 run on the assembly; we retain the GtRNAdb values because they are the standard-of-record for yeast and a re-run produces nearly identical numbers.

*Conditional logit scope.* The conditional logit model (@sec:res-condlogit) is event-level explanatory, not direct proof of biological causality. The candidate universe comprises all $approx $1,280 single-codon reassignments regardless of mechanistic plausibility, and the tRNA-complexity proxy (Hamming distance to nearest target-AA codon) is heuristic. The non-significance of this particular proxy does not exclude that a richer model of tRNA repertoire change would contribute explanatory value.

*Survivorship bias.* All event-level analyses operate on reassignments that _persisted_ in extant lineages. Topology-breaking reassignments that produced fitness collapse and lineage extinction are unobservable. The depletion result is therefore consistent with both selection against attempting topology-breaking moves and selection against the lineages that attempted them; cross-sectional NCBI data cannot adjudicate. Both interpretations support the conclusion that codon-family connectivity constrains evolutionary trajectories, with different mechanistic implications for the locus of selection.

*Phylogenetic non-independence ceiling.* The 28 de-duplicated reassignment events used in the topology-avoidance hypergeometric test cluster at the codon level (e.g., the recurrent UGA$arrow.r$Trp event is observed across many distantly related mitochondrial lineages and is collapsed to a single event for the de-duplicated test) but they are *not* independent at the level of evolutionary origin: mitochondrial codes share an ancestral mitochondrial-decoding regime, ciliate nuclear codes share a ciliate-specific eRF1 trajectory, and so on. A conservative lower bound on the number of phylogenetically independent topology-preservation origins is therefore on the order of 4--6 (mitochondrial, ciliate nuclear, _Candida_-clade nuclear, trypanosomatid, _Mycoplasmoides_, and one or two minor lineages), well below the 22 non-breakers the de-duplicated count reports. The hypergeometric $p$ should be read against that ceiling rather than against $n = 22$. The conditional-logit clade-exclusion sensitivity (Supplement §S7; seven phylogenetically motivated regimes) and the within-clade enumerations of breakers and non-breakers in §3.4 jointly bound how strongly the depletion magnitude can be claimed; we have not, in this version, fitted a formal mixed conditional logit with random effects at the clade level, and we flag that as a natural extension.

*Multiple-comparison correction.* Correction was applied within analysis families (metric, $rho$-sweep, per-table, topology-avoidance, synthetic-recoding) rather than across all descriptive, exploratory, and confirmatory quantities. The family structure was the natural organisation for these analyses; we did not file an external pre-registration. The eight test families address conceptually distinct questions (cross-metric robustness; $rho$-interpolation between $Q_6$ and $H(3,4)$; per-NCBI-table preservation; the $2 times 2$ topology-avoidance audit; clade-exclusion sensitivity following #cite(<sengupta2007>, form: "prose"); the M1--M4#sub[H(3,4)] discrete-choice comparison; the 24-pairing tRNA-enrichment Stouffer combination with MIS robustness; and the cross-study recoding reanalysis); their test statistics are non-nested. For transparency: under cross-family Bonferroni at $alpha = 0.05/8 = 6.25 times 10^(-3)$, the headline $H(3,4)$ topology depletion ($p = 1.3 times 10^(-6)$), the cross-metric coloring optimality ($p lt.eq 0.006$ across four metrics, already at threshold), the $rho$-sweep ($p lt.eq 0.006$ across five $rho$), the conditional-logit $Delta$AICc gap of #str(calc.round(m1f.aicc - m3f.aicc, digits: 0)) (M1$arrow.r$M3, well above any threshold), and the per-table BH--FDR result (robust by construction) all survive. The MIS worst-case tRNA enrichment ($p = 0.045$) does *not* survive, consistent with our classification of it as suggestive rather than confirmatory. Within-family corrections also hold: four-metric family $alpha = 0.0125$, $rho$-sweep family $alpha = 0.01$, topology-avoidance $2 times 2$ family $alpha = 0.0125$, Ostrov segment tests $alpha = 0.0167$, and per-table tests use BH--FDR (#pt.n_significant_bh of #pt.n_tables significant). Exploratory and suggestive analyses are labeled as such throughout.

Whether the code's optimality reflects adaptive selection or non-adaptive carry-over from a primordially constrained starting point @koonin2009 cannot be resolved by cross-sectional data alone. Multiple theoretical frameworks address this question: #cite(<sella2006>, form: "prose") showed that coevolution of genes and codes generates error-correcting structure resembling the standard code; #cite(<vetsigian2006>, form: "prose") argued that communal evolution via horizontal gene transfer accounts for both universality and optimality; and #cite(<novozhilov2007>, form: "prose") found the standard code sits roughly halfway up a local peak in a rugged fitness landscape (see @digiulio2005 for a review of origin theories). The topology avoidance result is consistent with all these frameworks: it demonstrates a constraint on reassignment trajectories but does not distinguish whether the constraint is selective or structural. #cite(<novozhilov2009>, form: "prose") showed that putative primordial 16-supercodon codes are "nearly optimal" even without direct selection, suggesting the current code's optimality may be partly inherited.

The $"GF"(2)^6$ representation is best understood as an analytical decomposition tool rather than a claim about the biological primacy of binary coordinates. Its specific contributions---the $rho$-sweep interpolation between $Q_6$ and $H(3,4)$ and the encoding-sensitivity audit---are detailed in §4.1; the underlying biology is the assignment of chemically similar amino acids to mutationally proximate codons, a property that holds regardless of the coordinate system, as confirmed by the encoding-independent $H(3,4)$ topology-avoidance result.


// ============================================================
//  5. CONCLUSION
// ============================================================
= Conclusion <sec:conclusion>

The standard genetic code is robustly low-cost: across four amino acid distance metrics and across a continuum of mutation graphs spanning $Q_6$ and the complete single-nucleotide graph $H(3,4)$, it consistently sits in the extreme tail of block-preserving null models ($p lt.eq 0.006$ under every combination). Across the 27 NCBI translation tables this near-optimality is preserved at the variant level: of the #pt.informative_total codes whose distance from the standard makes the per-table test informative, #pt.informative_significant retain top-5% placement after BH--FDR correction, with only yeast mitochondrial (table 3) marginally above the 5% threshold; the #pt.near_standard_total near-standard tables are formally significant in the same correction but we interpret these as primarily confirmatory of standard-code geometry. Code evolution therefore appears to occur within a constrained low-cost region rather than freely across code space.

Two partly independent constraints govern code evolution within that region. Persistent natural reassignments preferentially enter neighborhoods of lower local Grantham mismatch, and they avoid fragmenting amino-acid codon-family connectivity. Topology-breaking moves are approximately #str(calc.round(tk43.depletion_fold, digits: 1))-fold depleted under the encoding-independent $H(3,4)$ adjacency, and conditional-logit decomposition shows that topology adds explanatory value beyond physicochemistry, and vice versa ($Delta"AICc" gt.eq #str(calc.round(calc.min(cl.model_fits.M1_phys.aicc - cl.model_fits.M3_phys_topo.aicc, cl.model_fits.M2_topo.aicc - cl.model_fits.M3_phys_topo.aicc), digits: 0))$ under $Q_6$ topology; $Delta"AICc" gt.eq #str(calc.round(calc.min(cl.at("encoding_robustness", default: (:)).at("delta_aicc_M1_to_M3_k43", default: 0), cl.at("encoding_robustness", default: (:)).at("delta_aicc_M2k43_to_M3k43", default: 0)), digits: 0))$ under encoding-independent $H(3,4)$ topology). The two axes are only weakly correlated across the candidate landscape ($r_s = #str(calc.round(cl.phys_topo_rho, digits: 2))$); neither is a redundant restatement of the other. The topology coefficient survives all seven phylogenetic clade exclusions ($Delta"AICc"("M1" arrow.r "M3") gt.eq 49$ in the worst regime), and a posterior-predictive simulation under M3 reproduces the observed topology-breaking rate ($p = 0.60$). Several variant-code lineages show suggestive tRNA-gene enrichment for the reassigned amino acid (worst-case MIS Stouffer $p = #str(calc.round(trna.mis_worst_p, digits: 3))$ across 24 pairings from 18 tRNAscan-SE-verified genomes), consistent with compensatory accommodation.

A retrospective cross-study reanalysis of eight published genome-recoding experiments ($>$217,000 codon-level observations) sharpens the scope of the topology constraint. Family-boundary crossing is a real structural distinction over evolutionary timescales but did *not* predict acute transcriptomic perturbation in the Syn57 gene-level contrast, while Syn61 had earlier tolerated 18,218 boundary-crossing serine swaps as a class @fredens2019 (a genome-wide success rate, not a per-codon-position viability measure). This pattern is consistent with topology acting primarily as a constraint on evolutionary trajectories: gradual fragmentation of a codon family could create periods of decoding fragility that are bypassed in synthetic systems where decoding machinery is engineered up-front. Codon-space topology may therefore be a constraint on how genetic codes can change rather than on how they currently function---path-dependent rather than instantaneous. The three-layer integration of these descriptors (evolutionary topology, fine-scale recoding burden, accessibility-driven escape) is exploratory and hypothesis-generating: each layer-claim rests on a single retrospective contrast, the partition was not pre-specified, and confirmation requires the prospective tests outlined in §3.8.3 (@sec:res-threelayer).

The $"GF"(2)^6$ representation is best understood as an analytical decomposition rather than as a claim about biologically privileged binary coordinates. Its main contributions are the $rho$-sweep interpolating between $Q_6$ and $H(3,4)$, the $epsilon$-filtration of codon-family connectivity, and the encoding-sensitivity audit required for representation-dependent analyses. The KRAS--Fano clinical prediction, Serine distance-4 invariant, PSL(2,7) symmetry, and holomorphic-embedding claims were tested and rejected (@tbl:claims): several invariants of the representation turned out to be encoding artifacts. The surviving picture is encoding-independent: the genetic code is a low-cost coloring of $H(3,4)$, and observed reassignment trajectories preserve both physicochemical similarity and codon-family connectivity. Whether the topology axis reflects direct selection or a proxy for a specific feature of decoding architecture remains open. The present results support a working model in which natural genetic-code evolution is constrained by two partly independent axes: physicochemical smoothness and codon-family connectivity.


// ============================================================
//  ACKNOWLEDGEMENTS
// ============================================================
#heading(numbering: none)[Acknowledgements]

We thank the NCBI, GtRNAdb, and cBioPortal teams for maintaining public databases essential to this work. tRNAscan-SE 2.0.12 was developed by Chan and Lowe at UC Santa Cruz.

#heading(numbering: none)[CRediT author contribution statement]

*Paul Clayworth:* Conceptualization, Methodology, Formal analysis, Writing -- original draft, Writing -- review and editing.
*Sergey Kornilov:* Conceptualization, Methodology, Software, Formal analysis, Data curation, Visualization, Validation, Writing -- original draft, Writing -- review and editing.

#heading(numbering: none)[Declaration of generative AI and AI-assisted technologies]

During the preparation of this work the authors used GPT-5.2-Pro, GPT-5.5-Pro, Claude Opus 4.5, Claude Opus 4.7, GLM-5.1, Kimi-K2.5, MiniMax-M2.7, and Gemini 3 Pro to scaffold and review code, review the manuscript outline, and assist in drafting and revising the manuscript text. After using these tools, the authors reviewed and edited the content as needed and take full responsibility for the content of the publication. No generative AI or AI-assisted tools were used to create or alter figures or images in this manuscript.

#heading(numbering: none)[Declaration of competing interest]

The authors declare no competing interests.

#heading(numbering: none)[Ethical statement]

This study did not involve human subjects, animal experiments, or clinical samples. All analyses were performed on publicly available data: NCBI genome assemblies, NCBI gc.prt translation table definitions, GtRNAdb tRNA gene catalogues, the MSK-IMPACT mutation registry as released through cBioPortal, and previously published genome-recoding datasets. No new biological materials were generated.

#heading(numbering: none)[Data and code availability]

All code, raw data, and intermediate analysis outputs are publicly released in the `codontopo` repository at https://github.com/biostochastics/codontopo (version #stats._version, commit #raw("2f1ba6a"), tag #raw("v0.4.0")). Install with `pip install -e ".[all]"` or `uv sync --all-extras`. Analyses are fully reproducible via:

```
git clone https://github.com/biostochastics/codontopo.git
cd codontopo
codon-topo all --output-dir=./output --seed=135325
```

All inline statistics, tables, and figures in this manuscript and supplement are rendered directly from the pipeline outputs by the Typst sources (`manuscript.typ`, `supplement.typ`) included in the same repository. NCBI genome assembly accessions are listed in @tbl:organisms.

// ============================================================
//  REFERENCES
// ============================================================
#pagebreak()

#bibliography("references.bib", title: "References", style: "styles/elsevier-harvard.csl")
