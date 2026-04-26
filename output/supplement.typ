// ============================================================
// Supplementary Material
// Robust error-minimization in the genetic code
// Clayworth & Kornilov
// ============================================================

#set page(paper: "a4", margin: (x: 2.5cm, y: 2.5cm), numbering: "S1")
#set text(font: "Libertinus Serif", size: 10.5pt, lang: "en")
#set par(justify: true, leading: 0.65em)
#set heading(numbering: "S1.1")
// Inline figure placement (no float) keeps figures next to their text.
#set figure(gap: 1.0em, placement: none)
// Add generous breathing room around figures so they don't blend into the
// surrounding prose. Top/bottom margins are intentionally larger than the
// inter-paragraph spacing.
#show figure: it => {
  v(1.2em, weak: false)
  it
  v(1.2em, weak: false)
}
// Numbering: figures and tables get distinct counters with the "S" prefix so
// that the supplement reads "Figure S1", "Table S1", etc. independently of
// the main manuscript's numbering.
#set figure(numbering: n => "S" + str(n))
#show figure.where(kind: table): set figure(supplement: [Table])
#show figure.where(kind: image): set figure(supplement: [Figure])
// Allow long tables to break across pages so they don't push to a new page
// and leave a large gap at the bottom of the previous page.
#show figure.where(kind: table): set block(breakable: true)
#show table.cell.where(y: 0): set text(weight: "bold")
#show table: set text(number-width: "tabular", size: 9.5pt)
// Generous heading spacing so the supplement breathes properly.
#show heading.where(level: 1): it => {
  v(1.2em, weak: false)
  it
  v(0.7em, weak: false)
}
#show heading.where(level: 2): it => {
  v(0.8em, weak: false)
  it
  v(0.4em, weak: false)
}
#show heading.where(level: 3): it => {
  v(0.6em, weak: false)
  it
  v(0.3em, weak: false)
}
// Paragraph indent and spacing for body text
#set par(justify: true, leading: 0.7em, spacing: 0.7em, first-line-indent: 1.2em)

// ============================================================
//  PIPELINE DATA — shared with main manuscript
// ============================================================
#let stats = json("manuscript_stats.json")
#let cl = stats.condlogit
#let trna = stats.trna
#let coloring_enc = stats.coloring.at("encoding_sensitivity", default: (:))

// Helper: format a small probability in scientific notation (matches the
// helper in manuscript.typ so supplement and manuscript render identically).
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

#align(center)[
  #text(14pt, weight: "bold")[Supplementary Material]
  #v(0.3em)
  #text(11pt)[
    Robust error-minimization in the genetic code across \
    physicochemical metrics and variant codes: \
    a graph-theoretic analysis in GF(2)#super[6]
  ]
  #v(0.2em)
  #text(10pt)[Paul Clayworth and Sergey Kornilov]
  #v(0.8em)
]

#v(0.5em)

#block(stroke: (top: 0.5pt, bottom: 0.5pt), inset: (top: 0.6em, bottom: 0.6em, x: 0pt))[
  #set par(first-line-indent: 0pt, leading: 0.65em, spacing: 0.45em)
  #show outline.entry: it => link(it.element.location(), it.indented(it.prefix(), it.inner()))
  #outline(title: text(weight: "bold")[Contents], indent: 1em, depth: 2)
]

#v(0.8em)

*Roadmap.* This supplement is organised around three roles. Sections §S1--§S5 give the *evidentiary scaffolding* underpinning every main-text claim: a registered hierarchy of all 15 claims with status and justification (@sec:s-claims), and the four sensitivity analyses that probe whether each load-bearing result depends on a definitional choice---the 24-encoding sweep for cross-metric coloring optimality (@sec:s-encoding), the $2 times 2$ adjacency $times$ topology-breaking-definition audit (@sec:s-topology-defs), the per-encoding $Q_6$ topology-avoidance sweep (@sec:s-encoding-sweep), and the four candidate-universe denominators (@sec:s-denominator). Sections §S6--§S9 give the *robustness checks for the conditional-logit and topology-avoidance tests*: the IIA discussion and restricted-candidate sensitivity (@sec:s-iia, @sec:s-iia-restricted), conditional-logit clade exclusion (@sec:s-condlogit-clade), the standard-code-proximity audit for the per-table optimality test (@sec:s-proximity), and the hypergeometric/permutation clade-exclusion counterpart (@sec:s-clade). Sections §S10--§S15 give the *complete experimental and computational record*: the tRNAscan-SE 2.0.12 dataset on all 18 verified genomes (@sec:s-trna), the per-table reassignment database (@sec:s-reassignment), the synthetic-biology feasibility score implementation (@sec:s-feasibility), the per-variant KRAS--Fano falsification (@sec:s-kras), the full event-level conditional-logit specification with fitted coefficients and likelihood-ratio tests (@sec:s-condlogit), and software/version/reproducibility metadata (@sec:s-software). All numbers reported here are rendered from the same `manuscript_stats.json` and per-analysis JSON artifacts that drive the main text, so the two documents cannot drift from each other within a single pipeline run.

#v(0.5em)

// ============================================================
= Claim hierarchy with full justifications <sec:s-claims>

The 15 evaluated claims are organised into six evidentiary tiers (Supported / Suggestive / Exploratory / Falsified / Rejected / Tautological) and listed with full justifications in @tbl:s-claims. The hierarchy is registered in code (`src/codon_topo/reports/claim_hierarchy.py` in https://github.com/biostochastics/codontopo) so that every analysis run produces a verifiable status report (`codon-topo claims`) rather than the status table being maintained as prose. The "Falsified" and "Rejected" rows record claims this paper or earlier work has tested and found to be wrong; we list them so that the framework's negative scope is visible alongside its positive results.

#figure(
  table(
    columns: (1fr, auto, 1.1fr),
    align: (left, center, left),
    inset: (x: 6pt, y: 7pt),
    stroke: (x, y) => if y == 0 { (bottom: 0.7pt) } else { none },
    table.header([*Claim ID and statement*], [*Status*], [*Justification*]),
    [*hypercube_coloring_optimality* \ The standard genetic code is significantly error-minimizing under four independent physicochemical distance metrics: Grantham ($p = 0.006$), Miyata ($p < 0.001$), Woese polar requirement ($p = 0.003$), and Kyte--Doolittle hydropathy ($p = 0.001$).], [Supported], [Cross-metric sensitivity: Grantham $p = 0.006$, Miyata $p < 0.001$, polar requirement $p = 0.003$, Kyte--Doolittle $p = 0.001$ (block-preserving null, $n = $10,000). Stop penalty sensitivity (0/150/215/300): immaterial.],
    [*per_table_optimality_preservation* \ 26 of 27 NCBI translation tables remain in the top 5% of their own block-preserving null for Grantham edge-mismatch (BH--FDR corrected); only translation table 3 (yeast mito) exceeds the threshold.], [Supported], [Per-table block-preserving null applied to all 27 NCBI tables; significant fraction (BH--FDR $p < 0.05$) and per-table quantiles are reported in the main text and the per-table CSV `output/tables/T4_per_table_optimality.csv` (released in the `codontopo` repository). Translation table 3 (yeast mito) is the marginal exception.],
    [*optimality_rho_robustness* \ Coloring optimality is robust across all diagonal-edge weights $rho in [0, 1]$, including the full Hamming graph $H(3,4) = K_4 square.stroked K_4 square.stroked K_4$ of single-nucleotide substitutions ($p < 0.01$ at all $rho$ values).], [Supported], [$rho$-sweep at $n = $10,000: $p lt.eq 0.006$ at all $rho$ values; effect-size $z$ increases monotonically from $rho = 0$ to $rho = 1$.],
    [*topology_avoidance_depletion* \ Natural codon reassignments are depleted for topology-breaking changes --- moves that fragment an amino acid's codon family --- at approximately 21% of observed events vs 66--73% of the candidate landscape, robust to adjacency definition ($Q_6$ vs the full Hamming graph $H(3,4)$) and clade exclusion.], [Supported], [Permutation $p lt.eq 10^(-4)$ under both $Q_6$ and $H(3,4)$; hypergeometric $p = 1.6 times 10^(-8)$ ($Q_6$, new disconnection) and $p = 1.3 times 10^(-6)$ ($H(3,4)$, $Delta beta_0 > 0$). 5--7 of 28 de-duplicated events are topology-breaking versus $approx 64$--$75$% of the candidate landscape ($approx 3.0$--$3.4$-fold depletion). $Q_6$ is encoding-dependent (8 of 24 bijections give no depletion); $H(3,4)$ is encoding-independent and is the primary test. Clade-exclusion sensitivity (7 regimes): all $p < 10^(-5)$.],
    [*trna_enrichment_reassigned_aa* \ Organisms with variant genetic codes show elevated tRNA gene copy numbers for the reassigned amino acid, verified by tRNAscan-SE 2.0.12 on 18 NCBI genome assemblies across 5 variant genetic codes (24 pairings, 6 independent).], [Suggestive], [MIS worst-case $p = 0.045$. 18 tRNAscan-SE-verified assemblies (15 variant + 3 standard controls), 24 pairings across 5 variant codes.],
    [*bit_position_bias_weighted* \ Codon reassignment bit-flip distribution shows positional skew in $"GF"(2)^6$ coordinates, but signal is marginal after correcting for non-independence (de-duplicated $p = 0.075$).], [Exploratory], [Uniform $p = 0.006$ inflated by non-independence. De-duplicated $p = 0.075$.],
    [*mechanism_boundary_conditions* \ tRNA gene duplication accompanies codon reassignment in large nuclear genomes (ciliates, yeasts) but not in streamlined genomes (_Blastocrithidia_: anticodon stem shortening; _Mycoplasmoides_: anticodon modification).], [Exploratory], [Three-tier pattern: duplication / stem shortening / modification. Descriptive.],
    [*atchley_f3_serine_convergence* \ Serine's extreme Atchley Factor 3 score ($F_3 = -4.760$, most extreme of 20 amino acids, 2.24 SD below mean) converges with the $"GF"(2)^6$ topological disconnection: $F_3$ captures the mismatch between Serine's small physicochemical footprint and its disproportionate codon diversity, which the geometric framework identifies as maximal inter-family Hamming distance among 6-codon amino acids (4 vs 1 for Leu and Arg).], [Exploratory], [Serine $F_3 = -4.760$, 2.24 SD below mean. Complementary, not independent.],
    [*variant_code_disconnection_catalogue* \ Systematic survey of Hamming-graph connectivity across 27 NCBI translation tables identifies 4 lineage-collapsed variant-code amino-acid disconnections at $epsilon = 1$ under the default encoding: Thr in yeast mitochondrial code (table 3), Leu in chlorophycean mitochondrial codes (tables 16 and 22, collapsed to a single algal-mito event), Ala in _Pachysolen tannophilus_ nuclear code (table 26), and tripartite Ser in _Candida_-clade alternative yeast nuclear code (table 12).], [Exploratory], [4 variant-code cases at $epsilon = 1$ in $Q_6$ (default encoding): Thr (Table 3, yeast mito), Leu (Tables 16/22, chlorophycean and _S. obliquus_ mito; collapsed to one algal-mito event), Ala (Table 26, _Pachysolen_), Ser (Table 12, _Candida_). Separately, Table 32 (Balanophoraceae plastid, UAG→Trp) creates a 2-fold Trp pair (UGG, UAG) that remains connected at ε=1 but breaks the bit-5 two-fold filtration --- a filtration finding, not a disconnection.],
    [*kras_fano_clinical_prediction* \ XOR "Fano" relationships in $"GF"(2)^6$ predict enrichment of specific amino acids at KRAS G12 co-mutation sites.], [Falsified], [$p = 1.0$ across all 6 G12 variants. $n = $1,670 MSK-IMPACT mutations.],
    [*serine_min_distance_4_invariant* \ Serine's minimum UCN--AGY Hamming distance $= 4$ holds invariantly across all 24 base-to-bit encodings.], [Rejected], [16/24 encodings give distance 2. Only 8/24 give distance 4.],
    [*psl_2_7_symmetry* \ $"PSL"(2,7)$ is the fundamental symmetry group of the genetic code.], [Rejected], [No 64-dim irrep. #cite(<antoneli2011>, form: "prose").],
    [*holomorphic_embedding* \ The coordinate-wise map $"GF"(2)^6 arrow.r CC^3$ sending base-pairs to fourth roots of unity is a holomorphic embedding extending a character of $"GF"(8)^*$.], [Rejected], [Domain is finite discrete. Character identity fails: $i^2 = -1 eq.not 1$.],
    [*two_fold_bit_5_filtration* \ All 2-fold synonymous codon pairs in the standard code differ at exactly bit position 5 under the default encoding.], [Tautological], [Forced by encoding choice. Holds in 16/24 encodings.],
    [*four_fold_prefix_filtration* \ All 4-fold synonymous codon groups share a 4-bit prefix with 2-bit suffixes exhausting $"GF"(2)^2$.], [Tautological], [Trivial under any bijection from 4 bases to $"GF"(2)^2$.],
  ),
  caption: [Complete claim hierarchy with full claim statements and justifications.],
) <tbl:s-claims>


// ============================================================
= Encoding sensitivity analysis (coloring optimality) <sec:s-encoding>

There are $4! = 24$ distinct bijections from ${C, U, A, G}$ to $"GF"(2)^2$ (24 ways to assign the four nucleotide letters to the four 2-bit binary patterns 00/01/10/11). The default encoding used in the main text is $C |-> 00$, $U |-> 01$, $A |-> 10$, $G |-> 11$ (rationale stated in main-text §2.1). To check that the coloring-optimality result is not specific to that choice, we re-ran the Grantham block-preserving null analysis under each of the 23 alternative encodings, holding everything else fixed (per-encoding $n = $10,000 null draws, seed 135325). All #coloring_enc.at("n_encodings", default: 24) encodings yield significant optimality ($p < 0.05$ under the block-preserving null), with all per-encoding conservative $p$-values bracketed by $#str(calc.round(coloring_enc.at("p_min", default: 0.0002), digits: 4))$ (best) and $#str(calc.round(coloring_enc.at("p_max", default: 0.0278), digits: 4))$ (worst). The mean per-encoding quantile is #str(calc.round(coloring_enc.at("mean_quantile", default: 1.09), digits: 2))% (range #str(calc.round(coloring_enc.at("quantile_min", default: 0.01), digits: 2))%--#str(calc.round(coloring_enc.at("quantile_max", default: 2.77), digits: 2))%), confirming that the result is not an artifact of the default encoding.

Properties that are encoding-invariant:
- Serine disconnection at $epsilon = 1$ (holds under all 24 encodings)
- Coloring optimality significance (all 24 significant)
- Four-fold prefix filtration (tautological under any bijection)
- $H(3,4)$ adjacency graph (encoding-independent by construction)

Properties that are encoding-dependent:
- Serine inter-family minimum Hamming distance (4 in 8/24 encodings, 2 in 16/24)
- Two-fold bit-5 filtration (holds in 16/24 encodings)
- $Q_6$ (Hamming-1) adjacency partition of $H(3,4)$, and the topology-avoidance signal computed under $Q_6$ adjacency (see §S4)
- Specific score values and rank orderings

Full per-encoding results are written to `output/coloring_optimality.json` (block `encoding_sensitivity`) by the standard `codon-topo all` invocation, and reproducible from the `codontopo` repository.


// ============================================================
= Topology-breaking definitions: 2 × 2 audit <sec:s-topology-defs>

The topology-avoidance result depends on two independent definitional choices: (i) which adjacency graph defines codon-family connectivity ($Q_6$ Hamming-1 in the default $"GF"(2)^6$ encoding, vs the encoding-independent $H(3,4)$ single-nucleotide graph), and (ii) what counts as a "topology-breaking" candidate move (a candidate move creates a *new* disconnection in a previously connected family, vs the candidate move strictly increases $beta_0$ summed across families). To prevent silent dependence on either choice, we report all four combinations.

We define two notions of "topology-breaking" candidate move:

+ *New disconnection in a previously connected family*: a candidate move makes some amino acid disconnected at $epsilon = 1$ that was connected in the standard code (i.e., the amino acid was not previously in the disconnection catalogue but is after the move).

+ *Increase in components* ($Delta beta_0 > 0$, the conditional-logit feature): the total number of connected components, summed across amino-acid codon graphs, strictly increases:
  $ Delta_"topo" = sum_a beta_0(G_a^"after") - sum_a beta_0(G_a^"before") > 0. $

Both definitions are reported under both $Q_6$ adjacency (Hamming-1 in the default $"GF"(2)^6$ encoding) and $H(3,4)$ adjacency (full single-nucleotide adjacency, encoding-independent), giving four cells. All four share the same denominators (1,280 candidate moves, 28 de-duplicated observed events). The full $2 times 2$ result is shown in @tbl:s-topo-audit.

#figure(
  table(
    columns: (auto, 1fr, auto, auto, auto, auto),
    align: (left, left, right, right, right, right),
    inset: (x: 7pt, y: 8pt),
    stroke: (x, y) => if y == 0 { (bottom: 0.7pt) } else { none },
    table.header(
      [*Adjacency*], [*Topology-breaking definition*], [*$K\/N$*], [*$x\/n$*], [*Hyper. $p$*], [*RR (95% CI)*],
    ),
    [$H(3,4)$ (primary)], [$Delta beta_0 > 0$ (increase in components)], [846 / 1,280], [6 / 28], [$1.3 times 10^(-6)$], [0.32 (0.16--0.66)],
    [$H(3,4)$], [new disconnection in connected family], [822 / 1,280], [5 / 28], [$5.0 times 10^(-7)$], [0.28 (0.13--0.62)],
    [$Q_6$], [$Delta beta_0 > 0$ (increase in components)], [963 / 1,280], [7 / 28], [$2.2 times 10^(-8)$], [0.33 (0.17--0.63)],
    [$Q_6$], [new disconnection in connected family], [931 / 1,280], [6 / 28], [$1.6 times 10^(-8)$], [0.29 (0.14--0.60)],
  ),
  caption: [
    *Topology-avoidance $2 times 2$ definition $times$ adjacency audit.* All four cells share the same denominators (1,280 candidate moves, 28 de-duplicated observed events) and differ only in (i) which adjacency graph defines codon-family connectivity ($Q_6$ Hamming-1 vs encoding-independent $H(3,4)$) and (ii) what counts as a topology-breaking move ($Delta beta_0 > 0$ summed over amino acids vs creation of a *new* disconnection in a previously connected family). $K\/N$ = topology-breaking candidates among all candidates; $x\/n$ = topology-breaking observed events among all observed events. All four cells yield depletion in the same direction with comparable risk ratios (0.28--0.33) and hypergeometric $p < 10^(-5)$. The main text uses the $H(3,4)$, $Delta beta_0 > 0$ cell as primary (encoding-independent adjacency, definition matched to the conditional-logit feature).
  ],
) <tbl:s-topo-audit>

The full machine-readable audit is released as part of the `codontopo` repository's `output/` artifacts (the `definitions_audit` block in the topology-avoidance results file).


// ============================================================
= Topology avoidance: $Q_6$ encoding-sweep sensitivity <sec:s-encoding-sweep>

This section documents what motivated the demotion of $Q_6$ from primary to secondary topology adjacency in the manuscript.

The $H(3,4)$ Hamming graph is encoding-independent: every two-bit bijection from ${A, C, G, U}$ to ${0, 1}^2$ produces the same nucleotide-level adjacency, because $H(3,4)$ depends only on which nucleotide letters differ between two codons, not on their binary encoding. The $Q_6$ subgraph, however, depends on the encoding because the partition into Hamming-1 (single-bit-change) edges versus Hamming-2 (within-nucleotide diagonal) edges is bijection-specific---under one encoding two codons that differ at exactly one nucleotide may project to a Hamming-1 edge of $Q_6$, under another encoding they may project to a Hamming-2 edge. To test whether the $Q_6$ topology-avoidance result survives this representation choice, we recomputed the $Q_6$ candidate-landscape rate, observed rate, depletion fold, and hypergeometric $p$-value under all 24 base-to-bit bijections, holding the same 1,280 candidate moves and 28 observed events.

Across all #stats.topology_q6_encoding_sweep.n_encodings encodings:

- candidate-landscape rate: min #str(calc.round(stats.topology_q6_encoding_sweep.rate_poss_min * 100, digits: 1))%, median #str(calc.round(stats.topology_q6_encoding_sweep.rate_poss_median * 100, digits: 1))%, max #str(calc.round(stats.topology_q6_encoding_sweep.rate_poss_max * 100, digits: 1))%
- observed rate: min #str(calc.round(stats.topology_q6_encoding_sweep.rate_obs_min * 100, digits: 1))%, median #str(calc.round(stats.topology_q6_encoding_sweep.rate_obs_median * 100, digits: 1))%, max #str(calc.round(stats.topology_q6_encoding_sweep.rate_obs_max * 100, digits: 1))%
- depletion fold: min #str(calc.round(stats.topology_q6_encoding_sweep.depletion_min, digits: 2)), median #str(calc.round(stats.topology_q6_encoding_sweep.depletion_median, digits: 2)), max #str(calc.round(stats.topology_q6_encoding_sweep.depletion_max, digits: 2))
- hypergeometric $p$: min #sci(stats.topology_q6_encoding_sweep.p_min), median #sci(stats.topology_q6_encoding_sweep.p_median), max #str(calc.round(stats.topology_q6_encoding_sweep.p_max, digits: 3))

The default encoding (C=00, U=01, A=10, G=11) gives the largest depletion and the smallest $p$. Eight of 24 encodings give a candidate-landscape rate of approximately 36%, under which the observed rate of 21--36% does not significantly differ from candidate ($p > 0.5$). The $Q_6$ result is therefore not encoding-invariant. We accordingly present the encoding-independent $H(3,4)$ result as the primary topology-avoidance test, with $Q_6$ reported as a representation-specific decomposition for continuity with the broader $"GF"(2)^6$ framework. The default encoding was adopted in companion methodological work (#cite(<clayworth2026>, form: "prose")) prior to the present encoding sweep, on the grounds of visualization clarity (it places the standard code's nine 2-fold-degenerate amino acids on bit-5 differences). We make no claim that this encoding is biologically privileged. The full per-encoding sweep is released as part of the `codontopo` repository's `output/` artifacts (the `Q6_encoding_sweep` block in the topology-avoidance results file); @fig:s-encoding-sweep visualizes the per-encoding depletion fold and hypergeometric $p$.

#figure(
  image("figures/FigS_encoding_sweep.png", width: 65%),
  placement: auto,
  scope: "parent",
  caption: [
    Per-encoding $Q_6$ topology-avoidance depletion across all 24 base-to-bit bijections. Bars are sorted by depletion fold; the dashed line at $1.0$ marks no-depletion. Bars highlighted in blue are statistically significant ($p < 0.05$); grey bars are not. Eight of 24 encodings (depletion fold $approx 1.0$, $p > 0.5$) place the $Q_6$ candidate-landscape rate at $approx 36$% rather than 73%, eliminating the depletion signal. The default encoding (C=00, U=01, A=10, G=11) yields the largest depletion (3.4-fold). The encoding-independent $H(3,4)$ result is constant across all 24 bijections and is the primary topology-avoidance test in this work.
  ],
) <fig:s-encoding-sweep>


// ============================================================
= Reassignment candidate-universe denominator sensitivity <sec:s-denominator>

The hypergeometric and permutation tests of topology avoidance, and the conditional-logit denominator, both depend on a definition of the candidate-move universe $cal(M)(C)$ from the standard code $C$. Three definitional choices interact: whether stop-codon targets are admitted, whether identity moves ($y = C(x)$) are excluded, and whether the source codon may itself be a stop. The combinations give the variants below; we adopt U1 as primary throughout and report U2 and U4 as sensitivity universes.

Formally, the primary candidate-universe is $cal(M)(C) = {(x, y) : x in cal(C), y in cal(A)_20 union {"Stop"}, y eq.not C(x)}$ with $abs(cal(M)(C)) = 64 times 20 = $1,280. The four definitional variants and their sizes are listed in @tbl:s-denominators.

#figure(
  table(
    columns: (auto, auto, 1fr),
    align: (left, right, left),
    inset: (x: 7pt, y: 7pt),
    stroke: (x, y) => if y == 0 { (bottom: 0.7pt) } else { none },
    table.header([*Universe*], [*$abs(cal(M))$*], [*Definition*]),
    [U1: 21-label, no identity (primary)], [1,280],
      [$y in cal(A)_20 union {"Stop"}$, $y eq.not C(x)$; each codon has exactly 20 alternative labels],
    [U2: AA-only, no identity], [1,219],
      [$y in cal(A)_20$, $y eq.not C(x)$; 61 sense codons get 19 AA alternatives, 3 stop codons get 20 AA alternatives],
    [U3: AA, with no-ops], [1,280],
      [$y in cal(A)_20$, no identity restriction; 64 candidates are no-ops $y = C(x)$],
    [U4: stop-inclusive with no-ops], [1,344],
      [$y in cal(A)_20 union {"Stop"}$, no identity restriction; 64 no-ops included],
  ),
  caption: [
    *Candidate-universe denominator definitions for the topology-avoidance test.* The four variants (U1--U4) interact along three definitional axes: whether stop-codon targets are admitted ($cal(A)_20 union {"Stop"}$ vs $cal(A)_20$), whether identity moves $y = C(x)$ are excluded, and whether the source codon may itself be a stop codon. We adopt U1 as the primary universe because it (i) has a uniform alternative count per codon, (ii) cleanly excludes identity moves which contribute no signal, and (iii) admits stop-codon reassignment, which is biologically attested.
  ],
) <tbl:s-denominators>

Topology-avoidance results under U2 and U4 are qualitatively identical (depletion remains $p < 10^(-5)$); the per-cell counts differ only by the constant rescaling $K_2 = K_1 - text("(stop-target candidates)")$ for U2. Detailed per-universe numbers are released as part of the `codontopo` repository's `output/` artifacts (the `denominator_sensitivity` block in the topology-avoidance results file).


// ============================================================
= Conditional logit: IIA assumption and explanatory framing <sec:s-iia>

The conditional-logit framework assumes Independence of Irrelevant Alternatives (IIA): for any two candidate moves $m_1$ and $m_2$ in a choice set $cal(N)$, the relative probability $P(m_1) \/ P(m_2)$ is the same regardless of which other candidates are in $cal(N)$. In the reassignment context, candidate moves are not exchangeable: a move that targets an amino acid already serviced by many codons is plausibly more substitutable with similar moves than the IIA structure allows.

We adopt IIA here because the goal is _explanatory_ rather than _predictive_: the test asks whether topology adds explanatory value beyond physicochemical cost (LR test, $Delta$AICc) within the same candidate set, not whether the model accurately predicts which specific reassignment will occur next. The relative-probability ratios that IIA constrains are not the quantities we report; the LR statistics depend only on whether the observed event occupies a high-likelihood position within the candidate set, which is robust to substitution patterns among non-observed candidates. A mixed logit relaxing IIA would be more appropriate if the goal were prediction; future work could pursue this once a larger event set permits identification of mixed-logit covariance parameters. The posterior-predictive simulation reported in the main text (§3.5; observed topology-breaking rate 0.076 vs simulated mean 0.077; $p = 0.60$) provides a calibration check that the chosen explanatory model reproduces the marginal feature distribution rather than just maximizing in-sample AICc.

A second, structurally distinct concern about the candidate set is its *composition*---whether the universe of $approx 1{,}280$ moves contains so many biologically implausible alternatives that the fitted topology coefficient is partly identifying their absence rather than a topology-specific avoidance effect. We address this concern in @sec:s-iia-restricted by refitting M1--M4 on candidate sets that have been pruned to biologically plausible alternatives only.

== Candidate-set composition: restricted-candidate sensitivity <sec:s-iia-restricted>

A separate concern about the candidate set is that the universe of $approx $1,280 single-codon moves admits biologically catastrophic alternatives---reassigning AUG-Met, simultaneous multi-codon changes implicit in the single-step framing, reassignments to stop in essential codons---that natural selection has already removed from the option set. Models with strongly negative coefficients on $Delta_"topo"$ and $Delta_"phys"$ may therefore be partly rediscovering that natural reassignments are not biologically catastrophic, inflating ΔAICc magnitudes beyond what the explanatory thesis (topology adds value beyond physicochemistry) strictly requires. The qualitative claim is unaffected by this concern, but the magnitudes need calibration.

We address this by refitting M1--M4 (and the $H(3,4)$ verification variants) on a *restricted candidate set*: at each event-step we retain only candidates whose target amino acid is already serviced by a codon at Hamming distance $lt.eq d$ from the reassigned codon (i.e., $Delta_"tRNA" lt.eq d$), with $d in {1, 2, 3}$. The observed move is always retained regardless of its $Delta_"tRNA"$ so the likelihood remains well-defined. The Hamming-$lt.eq 2$ filter is the primary biological-plausibility criterion (target AA already accessible via at most two nucleotide-position changes); the $lt.eq 1$ filter is the most stringent biological-plausibility cut, and $lt.eq 3$ is included as a looser bound to bracket the choice. @tbl:s-condlogit-restricted reports the resulting $Delta$AICc gaps under all four candidate sets.

#let cl_restr = cl.at("restricted_candidate", default: (:))
#let _restr_block(d) = cl_restr.at("by_max_trna", default: (:)).at(d, default: (:))
#let _restr_field(d, key) = {
  let b = _restr_block(d)
  let da = b.at("delta_aicc", default: (:))
  da.at(key, default: 0)
}
#let _restr_csize(d) = {
  let b = _restr_block(d)
  let cs = b.at("candidate_summary", default: (:))
  cs.at("candidates_mean", default: 0)
}
#let _restr_obs(d) = {
  let b = _restr_block(d)
  let cs = b.at("candidate_summary", default: (:))
  (cs.at("observed_in_filtered_set", default: 0), cs.at("observed_total", default: 0))
}

#if cl_restr.len() > 0 [
  #figure(
    table(
      columns: (auto, auto, auto, auto, auto, auto, auto),
      align: (center, right, right, right, right, right, right),
      inset: (x: 6pt, y: 7pt),
      stroke: (x, y) => if y == 0 { (bottom: 0.7pt) } else { none },
      table.header(
        [*Filter*],
        [*Mean cand.*],
        [*Obs. kept*],
        [*$Delta$AICc \ M1$arrow.r$M3*],
        [*$Delta$AICc \ M2$arrow.r$M3*],
        [*$Delta$AICc \ M3$arrow.r$M4*],
        [*$Delta$AICc \ M1$arrow.r$M3#sub[$H(3,4)$]*],
      ),
      [Full ($d = 7$)],
        [#str(calc.round(cl_restr.at("full_set_summary", default: (:)).at("candidates_mean", default: 0), digits: 0))],
        [---],
        [#str(calc.round(cl.model_fits.M1_phys.aicc - cl.model_fits.M3_phys_topo.aicc, digits: 0))],
        [#str(calc.round(cl.model_fits.M2_topo.aicc - cl.model_fits.M3_phys_topo.aicc, digits: 0))],
        [#str(calc.round(cl.lr_tests.at("M3_vs_M4", default: (lr_statistic: 0)).lr_statistic, digits: 1))],
        [#str(calc.round(cl.at("encoding_robustness", default: (:)).at("delta_aicc_M1_to_M3_k43", default: 0), digits: 0))],
      ..("3", "2", "1").map(d => {
        let (kept, total) = _restr_obs(d)
        (
          [$Delta_"tRNA" lt.eq #d$],
          [#str(calc.round(_restr_csize(d), digits: 0))],
          [#kept / #total],
          [#str(calc.round(_restr_field(d, "M1_to_M3"), digits: 0))],
          [#str(calc.round(_restr_field(d, "M2_to_M3"), digits: 0))],
          [#str(calc.round(_restr_field(d, "M3_to_M4"), digits: 1))],
          [#str(calc.round(_restr_field(d, "M1_to_M3_k43"), digits: 0))],
        )
      }).flatten(),
    ),
    caption: [
      Restricted-candidate sensitivity for the conditional-logit comparison. Each row refits M1, M2, M3, M4 (and the $H(3,4)$ topology variants) on a candidate set filtered to $Delta_"tRNA" lt.eq d$, where $Delta_"tRNA"$ is the Hamming distance from the reassigned codon to the nearest existing codon for the target amino acid. Observed moves are always retained so likelihoods remain comparable. The "Full" row reproduces the unrestricted main-text numbers. ΔAICc gaps shrink as the candidate set is restricted---this is expected, since removing biologically implausible candidates removes contrasts the model otherwise exploits---and the magnitude of the ΔAICc(M1$arrow.r$M3) gap is therefore best read at the primary $d = 2$ threshold (~727 candidates per choice set), where it is ~60. Under the most stringent $d = 1$ filter (~275 candidates), ΔAICc(M1$arrow.r$M3) shrinks to ~14 but stays above the conventional 10-threshold; ΔAICc(M2$arrow.r$M3) stays large at all $d$. The qualitative claim "topology adds explanatory value beyond physicochemistry" is therefore robust to candidate-set composition. The ΔAICc(M3$arrow.r$M4) column under the restricted filters is *not* interpretable, because the filter is defined on $Delta_"tRNA"$ so the M4 tRNA-feature distribution shifts mechanically with the filter; the unrestricted-set ΔAICc(M3$arrow.r$M4) of ~2 is the calibrated reading and shows the heuristic tRNA proxy is uninformative.
    ],
  ) <tbl:s-condlogit-restricted>
] else [
  // Fallback: pipeline did not populate the restricted_candidate block.
  Restricted-candidate sensitivity table will be populated by `codon-topo all` (block `cl.restricted_candidate` in `manuscript_stats.json`).
]

Two structural points are worth flagging beyond the table itself. First, the ΔAICc(M1$arrow.r$M3) gap is exactly the quantity sensitive to candidate-set composition: physicochemistry alone (M1) is a flat baseline that benefits most from the implausible-candidate filler being present, so removing it shrinks the M1$arrow.r$M3 distance more than the M2$arrow.r$M3 distance. Second, ΔAICc(M2$arrow.r$M3) is robust across all three thresholds (≈73--85), reflecting that physicochemistry contributes the same incremental signal regardless of which candidates surround the observed events. The overall reading is therefore: the unrestricted-set ΔAICc(M1$arrow.r$M3) ≈ 110 over-states the *magnitude* of the topology--physicochemistry separation, but does not over-state its existence; the $d = 2$ filter ($approx 60$) gives the biologically-calibrated effect size and the $d = 1$ floor ($approx 14$) shows the qualitative claim survives even under the most stringent biological-plausibility cut available.


// ============================================================
= Conditional logit: clade-exclusion sensitivity <sec:s-condlogit-clade>

The conditional-logit framework can confound a single clade's reassignment events with a global topology-avoidance signal if that clade is unusually extensive (or unusually idiosyncratic). To bound this concern, we refit M1--M4 under each of seven clade-exclusion regimes (@tbl:s-condlogit-clade), removing the indicated NCBI translation tables from the event-step list and re-running the full conditional-logit pipeline end-to-end (build candidate sets $arrow.r$ enumerate event orderings up to $k!$ for $k lt.eq 6$ $arrow.r$ vectorised maximum-likelihood fit $arrow.r$ ΔAICc against M1 and M2). The seven exclusion regimes match those applied to the hypergeometric topology-avoidance test (@sec:s-clade), which themselves follow the phylogenetic-distribution analysis in #cite(<sengupta2007>, form: "prose").

// dynamic table from condlogit_clade_sensitivity.json (read as 'cl.clade_exclusion.rows' if present)
#let clade_rows = cl.at("clade_exclusion", default: (:)).at("rows", default: ())
#let _format_clade_label(s) = {
  let parts = s.replace("_", " ").split(" ")
  parts.slice(1).join(" ")
}

#if clade_rows.len() > 0 [
  #figure(
    table(
      columns: (auto, 1fr, auto, auto, auto, auto),
      align: (left, left, center, right, right, right),
      inset: (x: 6pt, y: 7pt),
      stroke: (x, y) => if y == 0 { (bottom: 0.7pt) } else { none },
      table.header(
        [*Excluded clade*],
        [*Tables out*],
        [*$n$*],
        [*$Delta$AICc \ M1$arrow.r$M3*],
        [*$Delta$AICc \ M2$arrow.r$M3*],
        [*$Delta$AICc \ M3$arrow.r$M4*],
      ),
      ..clade_rows.map(r => (
        [#_format_clade_label(r.excluded_clade)],
        [#r.excluded_tables.map(str).join(", ")],
        [#r.n_events_remaining],
        [#str(calc.round(r.delta_aicc_M1_to_M3, digits: 1))],
        [#str(calc.round(r.delta_aicc_M2_to_M3, digits: 1))],
        [#str(calc.round(r.delta_aicc_M3_to_M4, digits: 1))],
      )).flatten(),
    ),
    caption: [
      *Per-regime conditional-logit clade-exclusion sensitivity (7 regimes).* Each row refits M1--M4 with the indicated NCBI translation tables removed; $n$ = remaining event-steps. M3 (physicochemistry + $Q_6$ topology) remains robustly favored over both M1 (physicochemistry only) and M2 (topology only) in every regime. The minimum $Delta$AICc(M1$arrow.r$M3) across all seven regimes is #str(calc.round(cl.at("clade_exclusion", default: (:)).at("delta_M1_M3_min", default: 0), digits: 0)); the maximum is #str(calc.round(cl.at("clade_exclusion", default: (:)).at("delta_M1_M3_max", default: 0), digits: 0)). The conventional Burnham--Anderson reference of $Delta$AICc$>10$ was calibrated on linear-regression contexts, so we use it as a reference threshold rather than a formally calibrated cut-off in this conditional-logit setting; the posterior-predictive simulation reported in main-text §3.5 (observed topology-breaking rate 0.076 vs simulated 0.077, $p = 0.60$) is the more directly interpretable calibration check. $Delta$AICc(M3$arrow.r$M4) values near 2 indicate that adding the heuristic tRNA-distance proxy provides no incremental fit.
    ],
  ) <tbl:s-condlogit-clade>
] else [
  // Fallback static table (in case manuscript_stats.json was generated by an
  // older pipeline that omitted the per-regime rows).
  #figure(
    table(
      columns: (auto, 1fr, auto, auto, auto, auto),
      align: (left, left, center, right, right, right),
      inset: (x: 6pt, y: 7pt),
      stroke: (x, y) => if y == 0 { (bottom: 0.7pt) } else { none },
      table.header(
        [*Excluded clade*],
        [*Tables out*],
        [*$n$*],
        [*$Delta$AICc \ M1$arrow.r$M3*],
        [*$Delta$AICc \ M2$arrow.r$M3*],
        [*$Delta$AICc \ M3$arrow.r$M4*],
      ),
      [ciliates], [6, 10, 15, 27, 28, 29, 30], [52], [83.7], [31.6], [1.8],
      [yeast mito], [3], [60], [103.1], [87.7], [2.2],
      [CUG-clade], [12, 26], [64], [116.7], [94.9], [2.2],
      [metazoan mito], [2, 5, 9, 13, 14, 21], [40], [48.8], [88.6], [0.4],
      [algal mito], [16, 22], [63], [115.5], [85.7], [2.0],
      [protist mito], [4, 23], [63], [100.7], [92.4], [1.9],
      [hemichordate mito], [24, 33], [59], [92.3], [83.5], [2.1],
    ),
    caption: [
      *Per-regime conditional-logit clade-exclusion sensitivity (7 regimes).* M3 remains robustly favored across all regimes. Static fallback table (a current pipeline run will populate the dynamic version above).
    ],
  ) <tbl:s-condlogit-clade>
]


// ============================================================
= Per-table optimality: standard-code-proximity audit <sec:s-proximity>

This audit grounds the disaggregation reported in main-text §3.3, where we separated the per-table BH-FDR result into "informative-distance" tables ($gt.eq 3$ reassignments from standard) and "near-standard" tables ($lt.eq 2$ reassignments).

The methodological concern: a variant code that differs from the standard by only a few reassignments has a per-table block-preserving null distribution dominated by permutations very close to the standard code, simply because few permutations of the variant's block structure yield codes that are far from standard. In that limit, "table $X$ falls in the bottom 5% of permutations preserving $X$'s block structure" partly tests whether $X$ is close to the standard code (which it is, by construction), not whether $X$ is independently optimal.

For each NCBI translation table we computed three quantities alongside the unconditional per-table quantile: (i) the Hamming distance $d_H$ from the standard code (number of codons with a different AA label) for both the observed variant code and each block-preserving null draw, (ii) the variant's null quantile *conditional* on null draws within $plus.minus 2$ codons of the variant's $d_H$, and (iii) the fraction of null draws whose $d_H$ is at most the variant's $d_H$ (a complementary "proximity rank" of the variant against the null). The motivation for the conditional quantile was to ask whether variants with low unconditional $p$-values would still appear unusually low-cost when restricted to null draws of equivalent $d_H$. As @tbl:s-proximity shows, the conditional bucket turns out to be empty for every variant: under block-preserving permutation, every null draw has $d_H gt.eq 30$ from the standard code (the lowest observed $d_H$ in the null distribution across all 27 tables is 30, while the highest variant $d_H$ in the registry is 6), so no null draws fall within $plus.minus 2$ codons of any variant's $d_H$. The conditional quantile is therefore *not informative* about whether a variant's optimality is intrinsic versus standard-code-proximity-driven.

The proximity-rank diagnostic in the rightmost column is informative: every variant in the registry has $d_H$ smaller than every null draw, so all variants sit at the extreme low-$d_H$ tail of the block-preserving null distribution. This means each variant's per-table $p$-value is, structurally, partly a proximity-to-standard-code measurement, and the "informative-distance" vs "near-standard" disaggregation in main-text §3.3 is therefore based on absolute $d_H$ thresholds (3 reassignments) rather than on this audit's conditional quantile.

#let prox_audit = stats.coloring.at("per_table_proximity_audit", default: (:))
#let prox_rows = prox_audit.at("per_table", default: ())
// Drop the standard code itself (table 1, dH=0) from the audit table
#let _prox_filtered = prox_rows.filter(r => r.table_id != 1)
#let _prox_informative = _prox_filtered.filter(r => r.n_reassignments_from_standard >= 3).sorted(key: r => r.n_reassignments_from_standard).rev()
#let _prox_near = _prox_filtered.filter(r => r.n_reassignments_from_standard <= 2).sorted(key: r => r.n_reassignments_from_standard).rev()
#let _fmt_q(q) = if q == none [---] else [#str(calc.round(q, digits: 1))%]

#if _prox_filtered.len() > 0 [
  #figure(
    table(
      columns: (auto, auto, auto, auto, auto, auto),
      align: (center, center, right, right, right, right),
      inset: (x: 6pt, y: 6pt),
      stroke: (x, y) => if y == 0 { (bottom: 0.7pt) } else { none },
      table.header(
        [*Table*],
        [*$d_H$ (variant)*],
        [*Null $d_H$ range*],
        [*Quantile (uncond.)*],
        [*$n$ in $d_H plus.minus 2$ bucket*],
        [*Frac. null $d_H lt.eq d_H^"obs"$*],
      ),
      table.cell(colspan: 6, fill: luma(240), align: left)[
        *Informative-distance tables ($d_H gt.eq 3$ codon reassignments from standard).*
      ],
      ..(_prox_informative.map(r => (
        [#r.table_id],
        [#r.n_reassignments_from_standard],
        [#{r.null_dH_min}--#{r.null_dH_max}],
        [#_fmt_q(r.quantile_unconditional)],
        [#r.n_null_in_dH_bucket],
        [#str(calc.round(r.frac_null_dH_le_obs * 100, digits: 1))%],
      )).flatten()),
      table.cell(colspan: 6, fill: luma(240), align: left)[
        *Near-standard tables ($d_H lt.eq 2$ codon reassignments from standard).*
      ],
      ..(_prox_near.map(r => (
        [#r.table_id],
        [#r.n_reassignments_from_standard],
        [#{r.null_dH_min}--#{r.null_dH_max}],
        [#_fmt_q(r.quantile_unconditional)],
        [#r.n_null_in_dH_bucket],
        [#str(calc.round(r.frac_null_dH_le_obs * 100, digits: 1))%],
      )).flatten()),
    ),
    caption: [
      *Per-table standard-code-proximity audit.* For each NCBI translation table other than the standard code (table 1), $d_H$ is the number of codons with a different AA label from the standard code. *Null $d_H$ range* is the (min, max) of $d_H$ across the 1,000 block-preserving null draws for that table. *Quantile (uncond.)* is the variant's per-table block-preserving-null quantile from main-text §3.3 (lower = more error-minimising). *$n$ in $d_H plus.minus 2$ bucket* counts null draws with $d_H$ within $plus.minus 2$ of the variant's $d_H$ (would have been the conditional-quantile denominator); this is zero for every table because every block-preserving null draw has $d_H gt.eq 30$ while every variant has $d_H lt.eq 6$. *Frac. null $d_H lt.eq d_H^"obs"$* is the fraction of null draws closer to (or as close as) the variant in $d_H$; this is uniformly $0$%, confirming that every variant sits at the extreme low-$d_H$ tail of its block-preserving null. The conditional-quantile diagnostic is therefore not informative on this data, and the disaggregation in main-text §3.3 rests on the absolute $d_H gt.eq 3$ threshold rather than on this conditional analysis.
    ],
  ) <tbl:s-proximity>
] else [
  Per-table proximity-audit table will be populated by `codon-topo all` (block `coloring.per_table_proximity_audit` in `manuscript_stats.json`).
]

The full per-table audit, including $d_H$ summary statistics for the null distribution and the per-table observed score, is released in `output/coloring_optimality.json` under the `per_table_proximity_audit` key.


// ============================================================
= Topology avoidance: clade-exclusion sensitivity (hypergeometric/permutation) <sec:s-clade>

This section provides the hypergeometric/permutation counterpart to the conditional-logit clade-exclusion analysis in @sec:s-condlogit-clade. The two tests probe related but distinct quantities: @sec:s-condlogit-clade tests whether the $Delta_"topo"$ regression coefficient survives clade exclusion; this section tests whether the global landscape-vs-observed depletion ratio survives clade exclusion. Following the phylogenetic-distribution analysis of mitochondrial reassignment events by #cite(<sengupta2007>, form: "prose"), we iteratively excluded each of seven major taxonomic groups (@tbl:s-phylo-clade) and re-ran the $Q_6$ topology-avoidance hypergeometric on the reduced event set, holding the 1,280-move candidate landscape fixed.

#let phylo_rows = stats.phylo.at("clade_exclusion", default: ())
#let _phylo_label(s) = {
  let parts = s.replace("_", " ").split(" ")
  parts.slice(1).join(" ")
}

#if phylo_rows.len() > 0 [
  #figure(
    table(
      columns: (auto, 1fr, auto, auto, auto, auto),
      align: (left, left, center, right, right, right),
      inset: (x: 6pt, y: 7pt),
      stroke: (x, y) => if y == 0 { (bottom: 0.7pt) } else { none },
      table.header(
        [*Excluded clade*],
        [*Tables out*],
        [*$n$*],
        [*Breakers*],
        [*Rate (%)*],
        [*Hyper. $p$*],
      ),
      ..phylo_rows.map(r => (
        [#_phylo_label(r.excluded_clade)],
        [#r.excluded_tables.map(str).join(", ")],
        [#r.n_events_remaining],
        [#r.creates_disc],
        [#str(calc.round(r.rate_observed * 100, digits: 1))],
        [#sci(r.hypergeom_p)],
      )).flatten(),
    ),
    caption: [
      *Hypergeometric clade-exclusion sensitivity for the $Q_6$ topology-avoidance test.* Each row removes the indicated NCBI translation tables (matching the clade definitions of #cite(<sengupta2007>, form: "prose")) from the de-duplicated event list and re-runs the hypergeometric depletion test against the same 1,280-move candidate landscape. *Breakers* = topology-breaking events (new disconnection, $Q_6$ adjacency); *Rate* = observed breaker rate (vs candidate-landscape rate of 72.7% under the default encoding). All seven exclusions yield $p < 10^(-5)$. Excluding yeast mitochondrial (NCBI translation table 3) *strengthens* the depletion to $p approx 3.6 times 10^(-11)$ because that table contributes 4 of the 6 lineage-collapsed topology-breakers; this is the denominator-effect noted in main-text §3.4.
    ],
  ) <tbl:s-phylo-clade>
] else [
  // Fallback static table
  #figure(
    table(
      columns: (auto, 1fr, auto, auto, auto, auto),
      align: (left, left, center, right, right, right),
      inset: (x: 6pt, y: 7pt),
      stroke: (x, y) => if y == 0 { (bottom: 0.7pt) } else { none },
      table.header(
        [*Excluded clade*],
        [*Tables out*],
        [*$n$*],
        [*Breakers*],
        [*Rate (%)*],
        [*Hyper. $p$*],
      ),
      [ciliates], [6, 10, 15, 27, 28, 29, 30], [24], [6], [25.0], [$1.2 times 10^(-6)$],
      [yeast mito], [3], [24], [2], [8.3], [$3.6 times 10^(-11)$],
      [CUG-clade], [12, 26], [26], [5], [19.2], [$1.4 times 10^(-8)$],
      [metazoan mito], [2, 5, 9, 13, 14, 21], [22], [6], [27.3], [$9.8 times 10^(-6)$],
      [algal mito], [16, 22], [26], [5], [19.2], [$1.4 times 10^(-8)$],
      [protist mito], [4, 23], [27], [6], [22.2], [$4.8 times 10^(-8)$],
      [hemichordate mito], [24, 33], [27], [6], [22.2], [$4.8 times 10^(-8)$],
    ),
    caption: [
      *Hypergeometric clade-exclusion sensitivity for the $Q_6$ topology-avoidance test.* All seven exclusions yield $p < 10^(-5)$. Excluding yeast mitochondrial strengthens the depletion because NCBI translation table 3 contributes 4 of the 6 lineage-collapsed topology-breakers. Static fallback table; a current pipeline run will populate the dynamic version above.
    ],
  ) <tbl:s-phylo-clade>
]

In every exclusion, the depletion remains highly significant ($p < 10^(-5)$), confirming that the $approx 3.4$-fold ($Q_6$) and $approx 3.1$-fold ($H(3,4)$) depletion of topology-breaking changes is a pan-taxonomic pattern, not an artifact of any single lineage. The clade-exclusion robustness is a denominator effect (removing clades that contributed zero or one topology-breaking event each leaves the breakage rate essentially unchanged) rather than a numerator effect; the *avoidance* of topology-breaking moves is what propagates across many lineages, while the *attempts* to break topology are concentrated in yeast mitochondrial. Detailed per-clade counts are released as part of the `codontopo` repository's `output/` artifacts (the `phylogenetic_sensitivity.json` file).


// ============================================================
= Complete tRNA gene count data <sec:s-trna>

This section documents the empirical foundation for the tRNA enrichment analysis in main-text §3.6. We re-ran tRNAscan-SE 2.0.12 on the 18 NCBI assemblies listed below to produce a uniform tRNA-prediction protocol across the 24 pairings; main-text Table 8 reports the Total column from this analysis, while §S10.1 below resolves Total into its five mutually-exclusive components (Std20 + SeC + Supp + Undet + Pseudo) so the reader can see the per-genome composition.

== tRNAscan-SE verified organisms

All tRNA gene counts were obtained by running tRNAscan-SE 2.0.12 @chan2019 with Infernal 1.1.4 on NCBI genome assemblies. Eukaryotic organisms were scanned in `-E` mode; _Mycoplasmoides_ species in `-B` (bacterial) mode. Total counts are Infernal-confirmed (the more conservative count after the second-pass filter) and match the Total column of main-text Table 8. The five-column breakdown in @tbl:s-trna sums exactly to Total: Std20 (decoding the standard 20 amino acids) + SeC (selenocysteine, anticodon UCA, scanned via the dedicated `TRNAinf-euk-SeC.cm` model) + Supp (possible suppressor tRNAs with CTA/TTA/UCA anticodons) + Undet (predicted tRNAs whose isotype could not be determined) + Pseudo (predicted pseudogenes filtered by the Infernal second-pass isotype validation).

#figure(
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto, auto, auto, 1fr),
    align: (left, center, left, right, right, right, right, right, right, left),
    inset: (x: 5pt, y: 6pt),
    stroke: (x, y) => if y == 0 { (bottom: 0.7pt) } else { none },
    table.header([*Organism*], [*Tbl*], [*Assembly*], [*Total*], [*Std20*], [*SeC*], [*Supp*], [*Undet*], [*Pseudo*], [*Reassigned AA*]),
    [_T. thermophila_], [6], [GCF_000189635.1], [718], [672], [1], [37], [4], [4], [Gln: 54 (15+39)],
    [_P. tetraurelia_], [6], [GCF_000165425.1], [216], [202], [1], [11], [0], [2], [Gln: 18 (7+11)],
    [_O. trifallax_], [6], [GCA_000295675.1], [94], [83], [2], [6], [2], [1], [Gln: 8 (2+6)],
    [_P. persalinus_], [6], [GCA_001447515.1], [262], [228], [1], [15], [0], [18], [Gln: 20 (5+15)],
    [_H. grandinella_], [6], [GCA_006369765.1], [130], [121], [1], [3], [0], [5], [Gln: 9 (6+3)],
    [_E. aediculatus_], [10], [GCA_030463445.1], [80], [76], [1], [2], [1], [0], [Cys: 4 (3+1 UCA)],
    [_E. amieti_], [10], [GCA_048569255.1], [120], [103], [1], [6], [1], [9], [Cys: 8 (4+4 UCA)],
    [_E. focardii_], [10], [GCA_001880345.2], [62], [56], [1], [3], [0], [2], [Cys: 3 (1+2 UCA)],
    [_E. parawoodruffi_], [10], [GCA_021440025.1], [149], [128], [0], [3], [3], [15], [Cys: 9 (5+4 UCA)],
    [_E. weissei_], [10], [GCA_021440005.1], [495], [390], [0], [12], [14], [79], [Cys: 20 (17+3 UCA)],
    [_E. woodruffi_], [10], [GCA_027382605.1], [83], [74], [1], [3], [1], [4], [Cys: 4 (2+2 UCA)],
    [_B. stoltei_#super[#sym.section]], [15], [GCA_965603825.1], [169], [165], [1], [0], [1], [2], [Trp: 6],
    [_B. nonstop_ P57], [31], [GCA_028554745.1], [68], [65], [1], [2], [0], [0], [Trp: 2#super[\*]],
    [_M. genitalium_], [4], [GCA_000027325.1], [36], [35], [0], [1], [0], [0], [Trp: 1#super[#sym.dagger]],
    [_M. pneumoniae_], [4], [GCF_910574535.1], [37], [36], [0], [1], [0], [0], [Trp: 1#super[#sym.dagger]],
    [_S. coeruleus_], [1], [GCA_001970955.1], [272], [265], [1], [1], [3], [2], [---],
    [_I. multifiliis_], [1], [GCF_000220395.1], [150], [141], [1], [8], [0], [0], [---],
    [_F. salina_], [1], [GCA_022984795.1], [89], [85], [1], [1], [1], [1], [---],
  ),
  caption: [
    Complete tRNAscan-SE 2.0.12 results (Infernal-confirmed counts). The five tRNA categories sum exactly to Total: *Std20* = standard 20-AA tRNAs; *SeC* = selenocysteine (anticodon UCA); *Supp* = possible suppressor tRNAs (CTA/TTA/UCA anticodons that could read stop codons); *Undet* = predicted tRNAs whose isotype could not be determined; *Pseudo* = predicted pseudogenes filtered by the Infernal second-pass isotype validation. Reassigned-AA counts include standard and suppressor tRNAs counted toward the reassigned amino acid. #super[#sym.section]_Blepharisma stoltei_ is indexed as NCBI translation table 15 (Blepharisma Nuclear Code, defined as UAG$arrow.r$Gln) by virtue of its genus assignment, but the strain-specific MAC genome reads UGA$arrow.r$Trp via a dedicated suppressor tRNA-Trp(UCA) #cite(<singh2023>). The analysis here tests Trp-tRNA enrichment on that empirical reading rather than the legacy table-15 nominal code. #super[\*]Anticodon stem shortening (4-bp stem rather than canonical 5-bp). #super[#sym.dagger]Post-transcriptional modification (a single tRNA-Trp(CCA) reads both UGG and UGA after base modification).
  ],
) <tbl:s-trna>


== MIS (maximal independent set) analysis

To address non-independence from shared control organisms, we constructed a conflict graph where edges connect pairings sharing an organism. The Bron--Kerbosch algorithm with pivoting enumerates all maximal independent sets (MIS). Each MIS represents a set of pairings where no two share an organism and no additional pairing can be added without creating a conflict.

With 24 total pairings, 2 MIS of size 6 were identified. Both are significant at $p < 0.05$ under Stouffer's $Z$ combination of per-pairing Fisher exact $p$-values:

- *Best-case MIS* ($p = 0.044$): _S. cerevisiae_ mito/Thr, _S. obliquus_ mito/Leu, _P. tannophilus_/Ala, _P. tetraurelia_/Gln, _E. aediculatus_/Cys, _E. weissei_/Cys
- *Worst-case MIS* ($p = 0.045$): _S. cerevisiae_ mito/Thr, _S. obliquus_ mito/Leu, _C. albicans_/Ser, _P. tetraurelia_/Gln, _E. aediculatus_/Cys, _E. weissei_/Cys

== _Saccharomyces cerevisiae_ literature-derived control

The _Saccharomyces cerevisiae_ tRNA gene counts used in the yeast-mito Thr disconnection pairing come from GtRNAdb @chan2016 rather than tRNAscan-SE 2.0.12 in this work. We retain them to keep the well-characterized yeast-mito Thr case in the pairing set and flag the difference in source explicitly.


// ============================================================
= Complete reassignment database <sec:s-reassignment>

The reassignment database underpins both the topology-avoidance test (main-text §3.4) and the conditional-logit fits (main-text §3.5). It comprises every codon reassignment event across the 27 NCBI translation tables, relative to the standard code (table 1). Each event records the codon, source amino acid, target amino acid, and Hamming distance to the nearest codon already encoding the target amino acid in the standard code. De-duplication to unique (codon, target amino acid) pairs yields the event set used in the topology-avoidance hypergeometric and table-preserving permutation tests. The conditional-logit model operates on the full per-table event-step list, which preserves recurrent reassignments across independent lineages. Exact counts depend on the most recent NCBI gc.prt revision (currently v4.6 retrieved 2026-04-25); a current pipeline run yields #stats.reassignment.n_events raw events de-duplicating to #stats.reassignment.n_dedup unique (codon, target-AA) pairs.

#let reass_rows = stats.reassignment.at("per_table_counts", default: ())

#if reass_rows.len() > 0 [
  #figure(
    table(
      columns: (auto, 1fr, auto),
      align: (center, left, center),
      inset: (x: 7pt, y: 6pt),
      stroke: (x, y) => if y == 0 { (bottom: 0.7pt) } else { none },
      table.header([*Table*], [*NCBI translation table*], [*Events*]),
      ..reass_rows.map(r => (
        [#r.table_id],
        [#r.table_name],
        [#r.n_events],
      )).flatten(),
    ),
    caption: [
      *Per-table reassignment-event counts (raw, before de-duplication).* The 25 NCBI translation tables that differ from the standard code account for #stats.reassignment.n_events raw codon-reassignment events. Tables 1 (Standard) and 11 (Bacterial / Archaeal / Plant Plastid) share identical sense-codon mappings and contribute zero reassignment events. Tables 27 (Karyorelict Nuclear) and 28 (Condylostoma Nuclear) likewise share identical sense-codon mappings and are listed separately to match NCBI numbering. The 25-table counts sum to the all-events total, which de-duplicates to #stats.reassignment.n_dedup unique (codon, target-AA) pairs across the registry. The complete event-level database (with codon, source-AA, target-AA, Hamming-to-nearest-target columns) is released as `output/tables/T10_reassignment_db.csv` in the `codontopo` repository.
    ],
  ) <tbl:s-reassignment>
] else [
  // Fallback when manuscript_stats.json predates the per-table-counts field.
  Per-table reassignment-event summary will be populated by `codon-topo all` (block `reassignment.per_table_counts` in `manuscript_stats.json`). The complete event-level database is released as `output/tables/T10_reassignment_db.csv` in the `codontopo` repository.
]


// ============================================================
= Synthetic-biology feasibility score (visualization-only) <sec:s-feasibility>

For Figure 5A of the main text we use a heuristic feasibility score $S(m) in [0, 1]$ for each candidate single-codon reassignment $m$ from the standard code. The score is *not* used in any inferential test in the manuscript; it is a visualization aid for delineating high- versus low-feasibility regions of the 1,280-move candidate landscape. We report it in detail here for completeness.

The score combines three structural factors:

+ *Filtration preservation* $I_F(m) in {0, 1}$: indicator of whether the post-reassignment code retains a uniform 4-bit prefix for any 4-fold-degenerate amino acid (i.e., whether the reassignment preserves the four-fold prefix filtration described in main text §2.2).

+ *Local mismatch smoothness* $f(Delta_"phys"(m))$: a smooth decay of the move's local Grantham mismatch change, $f(Delta_"phys") = exp(- max(0, Delta_"phys") / d_0)$ with $d_0 = 100$ Grantham units. Negative $Delta_"phys"$ (improving local mismatch) saturates to 1; positive $Delta_"phys"$ decays toward 0.

+ *Accessibility* $g(Delta_"acc"(m))$: $g(Delta_"acc") = (7 - Delta_"acc") / 6$, where $Delta_"acc"(m)$ is the minimum Hamming distance from the reassigned codon to any codon currently encoding the target amino acid. Range $[0, 1]$ with full credit at distance 1, zero credit at distance 7.

The composite score is the equal-weighted average:

$ S(m) = (1/3) I_F(m) + (1/3) f(Delta_"phys"(m)) + (1/3) g(Delta_"acc"(m)). $

Filtration-preserving variants ($I_F = 1$) score $gt.eq 0.8$ under this composite; filtration-breaking variants ($I_F = 0$) score $lt.eq 0.75$. The visual gap between the two groups in Figure 5A is therefore mostly driven by $I_F$. Exact implementation: `src/codon_topo/analysis/synbio_feasibility.py` in the `codontopo` repository.


// ============================================================
= KRAS--Fano clinical prediction: detailed results <sec:s-kras>

The KRAS--Fano clinical prediction is the most concrete biomedical extrapolation considered for the $"GF"(2)^6$ framework, and we tested it directly to delimit scope. The conjecture: XOR ("Fano") relationships in $"GF"(2)^6$ predict enrichment of specific amino acids at KRAS G12 co-mutation sites. For each of the six common somatic KRAS G12 variants, we computed the predicted Fano partner codon as $G G U xor m_i = "partner"_i$ where $m_i$ is the mutant codon for that variant and $xor$ is bitwise XOR in $"GF"(2)^6$ (@tbl:s-kras). We then tested for co-mutation enrichment of the predicted partner amino acid against 1,670 KRAS-mutant samples from the MSK-IMPACT pan-cancer dataset @zehir2017 using Fisher's exact test with Bonferroni correction across the six variants.

#figure(
  table(
    columns: (auto, auto, auto, auto, auto, auto),
    align: (center, center, center, center, center, right),
    inset: (x: 7pt, y: 7pt),
    stroke: (x, y) => if y == 0 { (bottom: 0.7pt) } else { none },
    table.header(
      [*Variant*],
      [*WT codon*],
      [*Mutant codon*],
      [*Fano partner codon*],
      [*Predicted partner AA*],
      [*Fisher $p$ (Bonf.)*],
    ),
    [G12V], [GGU], [GUU], [CAC], [His (H)], [1.0],
    [G12D], [GGU], [GAU], [CUC], [Leu (L)], [1.0],
    [G12A], [GGU], [GCU], [CGC], [Arg (R)], [1.0],
    [G12R], [GGU], [CGU], [GCC], [Ala (A)], [1.0],
    [G12C], [GGU], [UGU], [ACC], [Thr (T)], [1.0],
    [G12S], [GGU], [AGU], [UCC], [Ser (S)], [1.0],
  ),
  caption: [
    *KRAS--Fano clinical prediction tested against MSK-IMPACT.* For each KRAS G12 variant, the Fano partner codon is computed as the bitwise XOR (in $"GF"(2)^6$) of the wild-type and mutant codons against a fixed reference; the predicted partner amino acid is the standard-code translation of that codon. Fisher's exact test for co-mutation enrichment of the predicted partner amino acid in MSK-IMPACT samples bearing each KRAS G12 variant ($n = $1,670 KRAS-mutant tumours; @zehir2017) yields $p = 1.0$ for all six variants after Bonferroni correction across the six tests. Odds ratios are near 1.0 throughout. The conjecture is cleanly falsified.
  ],
) <tbl:s-kras>

This result cleanly separates code-level error-minimization, which operates on the amino acid assignment structure, from mutation-level algebraic predictions, which would require DNA polymerase errors to respect the binary encoding---a biologically implausible mechanism. The Fano-line construction was a meaningful test of whether the binary representation has predictive power at the *mutation* level (in addition to the *assignment* level where it does carry signal); the answer is no, and we record this falsification explicitly in @tbl:s-claims.


// ============================================================
= Event-level conditional logit model: full details <sec:s-condlogit>

This section gives the implementation-level detail for the discrete-choice analysis reported in main-text §3.5: the candidate universe, the three feature definitions, the optimisation routine, the order-averaging procedure for unknown event sequences, the fitted raw and normalised coefficients, the likelihood-ratio tests, the predictor-confounding diagnostic, and the per-event observed-move ranks. The complete fitted JSON output is also released as `output/evolutionary_simulation.json` in the public repository.

== Model specification and candidate universe

The conditional logit model treats each observed reassignment as a discrete choice from the set of all single-codon reassignments available at the current code state. For a code with 64 codons and 21 possible amino acid/stop labels, the candidate set $cal(N)(C)$ contains exactly 1,280 moves (universe U1 in §S5: 64 codons $times$ 20 alternative labels, excluding identity assignments). All 1,280 candidates are evaluated at each step regardless of biological plausibility; the model's purpose is to test whether the observed moves are statistically distinguishable from uniform sampling given the three feature classes. As a conditional logit, the model implies an independence of irrelevant alternatives (IIA) structure within each choice set; the explanatory--rather--than--predictive framing (§S6) makes IIA tolerable for our purposes.

== Feature definitions

*Local physicochemical mismatch change ($Delta_"phys"$)*: For a reassignment of codon $c$ from amino acid $a$ to $a'$, this is the change in the sum of #cite(<grantham1974>, form: "prose") distances across all Hamming-1 edges incident to $c$:

$ Delta_"phys" = sum_({c,c'}: d(c,c')=1) [Delta(a', "code"(c')) - Delta(a, "code"(c'))] $

This captures the local impact on error-minimization at the reassigned position. When a neighboring codon is assigned Stop, we set $Delta(a, "Stop")$ equal to the maximum Grantham distance (215), consistent with the stop-penalty convention used in the coloring objective (Methods §2.2 of the main manuscript).

*Topology disruption ($Delta_"topo,Q_6"$)*: The total increase in connected components (at $epsilon = 1$) summed across all amino acid codon graphs under $Q_6$ adjacency. Stop codons are excluded from the topology sum; $Delta_"topo"$ counts connected-component changes only for amino acid codon families. A move that splits one amino acid's codon family into two components contributes $+1$; a move that fragments two families contributes $+2$; topology-preserving moves contribute $0$. The default $Delta_"topo,Q_6"$ uses $Q_6$ adjacency (Hamming-1 in the default $"GF"(2)^6$ encoding $C arrow.r 00$, $U arrow.r 01$, $A arrow.r 10$, $G arrow.r 11$). Because $Q_6$ adjacency is encoding-dependent (8 of 24 base-to-bit bijections give no $Q_6$ topology depletion at the candidate-landscape level; @sec:s-encoding-sweep), we additionally compute $Delta_"topo,H(3,4)"$ using the encoding-independent $H(3,4)$ adjacency (two codons are neighbors iff they differ at exactly one nucleotide position) and refit two model variants: M2#sub[$H(3,4)$] (topology-only with $Delta_"topo,H(3,4)"$) and M3#sub[$H(3,4)$] (physicochemistry + $Delta_"topo,H(3,4)"$). The encoding-robustness comparison $Delta$AICc(M1 $arrow.r$ M3#sub[$H(3,4)$]) vs $Delta$AICc(M1 $arrow.r$ M3) is reported in the next subsection.

*tRNA complexity proxy ($Delta_"tRNA"$)*: The Hamming distance from the reassigned codon to the nearest codon already encoding the target amino acid in the current code. This serves as a heuristic for the tRNA repertoire change required to service the reassigned codon, with larger distances implying more novel tRNA machinery needed.

== Fitted coefficients

@tbl:s-condlogit-coefs reports the maximum-likelihood coefficient estimates for all six model variants (M1, M2, M3, M4 under $Q_6$ topology, plus M2#sub[$H(3,4)$] and M3#sub[$H(3,4)$] under encoding-independent $H(3,4)$ topology).

#let mf = cl.model_fits
#let _raw = (name, idx) => {
  let f = mf.at(name, default: (:))
  let raws = f.at("weights_raw", default: ())
  if raws.len() > idx [
    #str(calc.round(raws.at(idx), digits: 4))
  ] else [
    ---
  ]
}
#let _norm = (name, idx) => {
  let f = mf.at(name, default: (:))
  let norms = f.at("weights_normalized", default: ())
  if norms.len() > idx [
    #str(calc.round(norms.at(idx), digits: 2))
  ] else [
    ---
  ]
}

#figure(
  table(
    columns: (auto, 1fr, auto, auto),
    align: (left, left, center, center),
    inset: (x: 6pt, y: 6pt),
    stroke: (x, y) => if y == 0 { (bottom: 0.7pt) } else { none },
    table.header([*Model*], [*Feature*], [*$hat(beta)$ (raw)*], [*$hat(beta)$ (normalized)*]),
    // Q_6 topology variants (legacy primary)
    [M1], [$Delta_"phys"$], _raw("M1_phys", 0), _norm("M1_phys", 0),
    [M2 ($Q_6$)], [$Delta_"topo,Q_6"$], _raw("M2_topo", 0), _norm("M2_topo", 0),
    [M3 ($Q_6$)], [$Delta_"phys"$], _raw("M3_phys_topo", 0), _norm("M3_phys_topo", 0),
    [M3 ($Q_6$)], [$Delta_"topo,Q_6"$], _raw("M3_phys_topo", 1), _norm("M3_phys_topo", 1),
    [M4], [$Delta_"phys"$], _raw("M4_full", 0), _norm("M4_full", 0),
    [M4], [$Delta_"topo,Q_6"$], _raw("M4_full", 1), _norm("M4_full", 1),
    [M4], [$Delta_"tRNA"$], _raw("M4_full", 2), _norm("M4_full", 2),
    // K_4^3 topology verification variants (encoding-independent)
    [M2 ($H(3,4)$)], [$Delta_"topo,H(3,4)"$], _raw("M2_topo_k43", 0), _norm("M2_topo_k43", 0),
    [M3 ($H(3,4)$)], [$Delta_"phys"$], _raw("M3_phys_topo_k43", 0), _norm("M3_phys_topo_k43", 0),
    [M3 ($H(3,4)$)], [$Delta_"topo,H(3,4)"$], _raw("M3_phys_topo_k43", 1), _norm("M3_phys_topo_k43", 1),
  ),
  caption: [
    Conditional logit coefficient estimates from the 27-table re-fit. Raw coefficients are on the original feature scale; normalized coefficients are on $z$-scored features (multiplying raw $hat(beta)$ by the global feature standard deviation). All $hat(beta)$ values are negative, indicating that observed reassignment histories preferentially populate moves that reduce physicochemical mismatch, avoid topology disruption, and (weakly) prefer target amino acids already serviced by nearby codons. The tRNA proxy coefficient is small and non-significant (LR $= 0.12$, $p = 0.73$). The $H(3,4)$ verification variants (M2#sub[H(3,4)], M3#sub[H(3,4)]) replace the encoding-dependent $Delta_"topo,Q_6"$ feature with the encoding-independent $Delta_"topo,H(3,4)"$. Their fitted ΔAICc values match the main-text encoding-robustness numbers.
  ],
) <tbl:s-condlogit-coefs>

== Likelihood ratio tests

@tbl:s-condlogit-lrt reports nested-model likelihood-ratio statistics under both $Q_6$ and $H(3,4)$ topology variants.

#let lr = cl.at("lr_tests", default: (:))
#let _lr_row = (key) => {
  let row = lr.at(key, default: (:))
  let s = str(calc.round(row.at("lr_statistic", default: 0), digits: 1))
  s
}

#figure(
  table(
    columns: (1fr, 1fr, auto, auto, auto),
    align: (left, left, center, center, center),
    inset: (x: 6pt, y: 6pt),
    stroke: (x, y) => if y == 0 { (bottom: 0.7pt) } else { none },
    table.header([*Restricted*], [*Full*], [*LR*], [*df*], [*$p$*]),
    [M1 (phys)], [M3 (phys+topo, $Q_6$)], [#_lr_row("M1_vs_M3")], [1], [$lt.double 10^(-10)$],
    [M2 (topo, $Q_6$)], [M3 (phys+topo, $Q_6$)], [#_lr_row("M2_vs_M3")], [1], [$lt.double 10^(-10)$],
    [M3 (phys+topo, $Q_6$)], [M4 (full)], [#_lr_row("M3_vs_M4")], [1], [$0.73$],
    [M1 (phys)], [M3#sub[$H(3,4)$] (phys+topo, $H(3,4)$)], [#_lr_row("M1_vs_M3_k43")], [1], [$lt.double 10^(-10)$],
    [M2#sub[$H(3,4)$] (topo, $H(3,4)$)], [M3#sub[$H(3,4)$] (phys+topo, $H(3,4)$)], [#_lr_row("M2_k43_vs_M3_k43")], [1], [$lt.double 10^(-10)$],
  ),
  caption: [
    Likelihood-ratio tests for nested conditional logit models, including both $Q_6$ topology variants (legacy primary) and $H(3,4)$ topology verification variants. Both topology features (Q_6 added to M1, H(3,4) added to M1) and physicochemistry (added to M2 / M2#sub[H(3,4)]) provide highly significant improvements. The tRNA-complexity proxy does not improve on the phys+topo model.
  ],
) <tbl:s-condlogit-lrt>

== Confounding diagnostic

Across $approx $84,000 candidate moves pooled from all choice sets, the Spearman correlation between $Delta_"phys"$ and $Delta_"topo,Q_6"$ is $rho = #str(calc.round(cl.at("phys_topo_rho", default: 0.15), digits: 2))$, indicating that the two predictors carry largely independent information. Moves that reduce physicochemical mismatch are only slightly more likely to also preserve topology, and the conditional logit framework accounts for any residual collinearity through simultaneous estimation.

== Observed move percentile ranks

Under the best model (M3), observed natural reassignments rank at the 89.5th percentile on average (mean rank 134 out of 1,280 candidates). Notable individual events: UGA$arrow.r$Trp ranks consistently above the 98th percentile across all tables where it occurs, indicating that this reassignment is the highest-ranked under the combined phys+topo score. The one clear outlier is the yeast mitochondrial CUU$arrow.r$Thr reassignment (30th percentile), consistent with NCBI translation table 3's status as the sole marginal exception in the per-table optimality analysis.

== Order-averaging implementation

For tables with $k > 1$ reassignment events, the temporal ordering is unknown. We marginalized the conditional logit likelihood over all $k!$ orderings, computing $L_"table" = (1\/k!) sum_sigma product_s P(m_(sigma(s))^* | cal(N)(C_(sigma,s)))$ using log-sum-exp for numerical stability. For all tables with $k lt.eq 6$ events (which covers every table in the dataset, including the largest --- yeast mitochondrial, $k = 6$, $6! = 720$ orderings), we enumerate all orderings exactly. The implementation precomputes feature matrices once per (table, ordering, step) into stacked numpy arrays, after which each Nelder--Mead / L-BFGS-B function evaluation reduces to a small batch of matrix multiplies and `scipy.special.logsumexp` calls; the four nested models are fit concurrently on a single shared bundle via `joblib` threads to avoid memory duplication.


// ============================================================
= Software and reproducibility <sec:s-software>

This section provides the metadata needed to reproduce every number in the manuscript and supplement bit-for-bit from the public repository. Every figure, table, and inline statistic is rendered by the Typst sources `manuscript.typ` and `supplement.typ` (also in the repository) from the JSON outputs of a single `codon-topo all` invocation, so the manuscript and supplement cannot drift from each other within a single pipeline run.

All analyses were performed using the `codon-topo` Python package (version #stats._version, commit #raw("bd92dfe"), tag #raw("v0.4.0")). The code is publicly released at https://github.com/biostochastics/codontopo. Dependencies and runtime requirements:

- Python 3.11, NumPy 1.24+, SciPy 1.10+
- R 4.5, ggplot2, ggpubr, viridis, patchwork (for figures)
- tRNAscan-SE 2.0.12 with Infernal 1.1.4
- Random seed: 135325 (all Monte Carlo analyses)

Analyses are fully reproducible via:
```
git clone https://github.com/biostochastics/codontopo.git
cd codontopo
pip install -e ".[dev]"
codon-topo all --output-dir=./output --seed=135325
python scripts/finalize_manuscript_stats.py
python scripts/generate_tables.py
Rscript src/codon_topo/visualization/R/strengthened_figures.R
```

The complete test suite (432 tests) can be run with `python3.11 -m pytest tests/ -m "not slow"`.

// ============================================================
// REFERENCES (shared with main manuscript)
// ============================================================

#bibliography("references.bib", title: "References", style: "elsevier-harvard")
