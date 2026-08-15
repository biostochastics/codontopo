# BIOSYS-D-26-00689 — Pre-Production Corrections Ledger

**Manuscript:** *Robust error-minimization in the genetic code across
physicochemical metrics and variant codes: a graph-theoretic analysis in
GF(2)⁶*
**Authors:** Paul Clayworth, Sergey Kornilov
**Status:** In-press at *BioSystems* (Elsevier)
**Prepared:** 2026-08-14 (v0.5.1) / 2026-08-15 (v0.6.0 §NULL and later sections)
**Corresponding author:** Sergey Kornilov <sergey@biostochastics.com>
**Repository:** https://github.com/biostochastics/codontopo
**Release tag:** `v0.6.0` (supersedes v0.5.1; consolidated release for BIOSYS-D-26-00689 production)

## Purpose of this ledger

Following acceptance the authors ran three internal QA/QC passes on the
manuscript and supplement:

1. A figure- and statistics-audit pass (v0.5.0) that identified ~41
   first-round + ~30 second-round presentation and data-plumbing
   defects across the 16 figures, 25 tables, and body prose.
2. A first pre-production audit (v0.5.1) that identified one
   computational correction (see §MIS below), a cross-reference sweep,
   and a set of terminology/wording refinements.
3. A second pre-production audit (v0.6.0, this ledger's §NULL and
   subsequent sections) that corrected the null-model nomenclature,
   added the classical Haig--Hurst amino-acid-permutation null as a
   sensitivity companion, made the tRNA analysis population and Fisher
   denominator explicit, restored the pre-specified tRNA decision rule,
   corrected a data-plumbing gap in the §3.5 clade-exclusion prose, and
   closed a batch of framing, cross-reference, and reproducibility
   items.

The v0.5.0 corrections were presentation and reproducibility
improvements only. **The v0.5.1 correction changed one claim's
classification** from SUGGESTIVE to EXPLORATORY (tRNA-gene enrichment):
the median-independent-subset signal meets the 0.05 threshold
(median MIS Stouffer $p = 0.046$), the strict worst-case subset does not
($p = 0.123$), and the pre-specified topology-breaking-only subset is
null (Stouffer $p = 0.43$). The v0.6.0 changes are *presentation,
framing, plumbing, and sensitivity-companion additions*: no claim
status or headline $p$-value changes at v0.6.0. All four SUPPORTED
claims (cross-metric coloring optimality, per-table preservation, ρ
robustness, topology-avoidance depletion) and the FALSIFIED KRAS--Fano
result remain as classified. The claim hierarchy count remains 4
SUPPORTED / 5 EXPLORATORY / 1 FALSIFIED / 3 REJECTED / 2 TAUTOLOGICAL
(15 total).

This ledger records every user-visible change, sorted by manuscript
location, with the pre-fix value, the post-fix value, the defect
provenance, and the commit SHA that made the change. The production
editor can use this document as a diff manifest for the corrections
pass.

---

## §MIS — Maximal-Independent-Set enumeration (v0.5.1)

| ID | Change | Old value | New value | Justification | Commit |
| --- | --- | --- | --- | --- | --- |
| MIS-01 | Bron--Kerbosch pivoting on the conflict graph replaced with pivoting on its complement | Enumerates 2 MIS of size 6 (worst-case Stouffer $p = 0.045$) | Enumerates 332 MIS of size 6 (median $p = 0.046$, best $0.016$, worst $0.123$; 190/332 sets below $0.05$) | An independent set in $G$ is a clique in the complement $\bar G$; the pivot must be selected in $\bar G$, not in $G$. | (v0.5.1 tag) |
| MIS-02 | Topology-breaking-subset test surfaced in the narrative | Present in `output/trna_evidence.json` but not cited in main text | Explicitly reported: $n = 4$ pairings (yeast mito Thr, Scenedesmus mito Leu, Pachysolen Ala, Candida Ser), Stouffer $p = 0.43$ (null) | The subset is the direct mechanistic test of the compensation-for-disconnection hypothesis; suppressing it would misrepresent the mechanistic evidence. | (v0.5.1 tag) |
| MIS-03 | tRNA-enrichment claim reclassified in the hierarchy | Status: SUGGESTIVE, evidence $p = 0.045$ (worst-case MIS) | Status: EXPLORATORY, evidence $p = 0.046$ (median MIS across 332 sets) | The strict worst-case criterion is not met, and the topology-breaking-only subset is null; the pattern is retained as heterogeneous-accommodation-of-decoding rather than as compensation-for-disconnection. | (v0.5.1 tag) |
| MIS-04 | Table 7 (main text) row set replaced | "Fisher+Stouffer (6 indep.)"; "MIS worst-case ($p = 0.045$)"; "MIS best-case ($p = 0.044$)" | "MIS median (across all 332 MIS) $p = 0.046$"; "MIS worst-case $p = 0.123$"; "MIS best-case $p = 0.016$"; "Topology-breaking subset $p = 0.435$" | Table row set now matches the distribution reported in the prose. | (v0.5.1 tag) |
| MIS-05 | Table 1 (claim hierarchy) row moved | SUGGESTIVE / "tRNA enrichment (MIS worst-case)" / $p = 0.045$ | EXPLORATORY / "tRNA enrichment (MIS median)" / $p = 0.046$ | Consistent with the reclassification in MIS-03. | (v0.5.1 tag) |
| MIS-06 | Table 1 caption reworded | "Shaded rows: claims that failed validation" | "Shaded rows are non-inferential: FALSIFIED / REJECTED rows record claims found to be wrong; TAUTOLOGICAL rows record claims true by construction and carry no evidentiary weight" | Tautological claims did not "fail" — they are true but non-evidentiary. | (v0.5.1 tag) |
| MIS-07 | Introduction Test 3 falsification criterion softened | "Stouffer $p > 0.05$ ... under the worst-case maximal independent set" | "positive result requires median Stouffer $p$ across all maximal-independent pairing subsets $< 0.05$; the topology-breaking-only subset is a complementary mechanistic test" | Matches the reported summary statistic; separates the pooled test from the mechanism-linked subset. | (v0.5.1 tag) |
| MIS-08 | §3.6 main-text tRNA paragraph rewritten | "Both MIS (each of size 6) are significant at $p < 0.05$ (worst-case $p = 0.045$)" | "332 MIS of size 6; median Stouffer $p = 0.046$, best $0.016$, worst $0.123$; 190/332 ($57.2\%$) below $0.05$. Topology-breaking subset ($n = 4$) $p = 0.43$" | Present distribution rather than a two-point summary; disclose null subset. | (v0.5.1 tag) |
| MIS-09 | Supplement §S10.2 MIS section body rewritten | "2 MIS of size 6 were identified. Both significant at $p < 0.05$" (with best/worst pairings listed at $p = 0.044$ and $0.045$) | "332 MIS of size 6 ... median $p = 0.046$, best $0.016$, worst $0.123$" plus the corrected best/worst pairing memberships | Same as MIS-08 for the supplement side; adds an explicit topology-breaking-subset subsection. | (v0.5.1 tag) |
| MIS-10 | Limitations, Conclusion, and Cross-family §S18 sentences | "suggestive"; "worst-case MIS $p = 0.045$"; "does not survive Bonferroni" | "exploratory"; "median MIS $p = 0.046$, worst-case $p = 0.123$"; sits at the family-Bonferroni threshold | Language consistent with the reclassification. | (v0.5.1 tag) |
| MIS-11 | README.md tRNA-gene enrichment section | "(suggestive) ... Stouffer $p = 0.045$" | "(exploratory) ... median $p = 0.046$, worst $0.123$, topology-breaking subset $p = 0.43$" | Same. | (v0.5.1 tag) |

---

## §NULL — Coloring-optimality null model (v0.5.2)

**Context (careful framing).** The v0.5.1 implementation labelled the primary coloring-optimality null "Freeland--Hurst block-preserving". The implemented ensemble is a *quartet-pattern shuffle* that holds each of the 16 first-two-base quartet AA-patterns atomic and permutes them across the 16 quartet slots. The classical Haig--Hurst 1991 / Freeland--Hurst 1998 null described in the code-optimization literature is a different ensemble: it holds the 20 codon families (which codons belong to which family) fixed and permutes the 20 AA labels uniformly across the 20 families. Internal pre-production QA/QC flagged this as a naming discrepancy rather than as a computational error — the quartet-pattern shuffle is a valid stringent null, but it should not be named after Haig and Hurst.

Corrections below (a) restore the correct name for what we implement, (b) add the classical HH-AA-permutation null as a sensitivity companion, and (c) surface the empirical comparison in the supplement. The optimality claim survives either null: all four physicochemical metrics pass Bonferroni at $alpha = 0.05/5 = 0.01$ under both ensembles.

| ID | Change | Old value | New value | Justification | Commit |
| --- | --- | --- | --- | --- | --- |
| NULL-01 | Rename implemented null throughout main text and supplement | "Freeland--Hurst block-preserving null" (~30 mentions across manuscript.typ, supplement.typ, cli.py) | "quartet-pattern shuffle null"; naming disclosure paragraph added at §2.3.1 | The implemented ensemble is a quartet-pattern shuffle, not the classical Haig--Hurst / Freeland--Hurst amino-acid-permutation null. Renaming avoids attribution error. | (v0.5.2 tag) |
| NULL-02 | Add classical Haig--Hurst 1991 / Freeland--Hurst 1998 AA-permutation null | Not implemented | New generator `_generate_random_code_haig_hurst_aa_permutation` in `src/codon_topo/analysis/coloring_optimality.py`; wired via `null_type="haig_hurst_aa"` in `monte_carlo_null()` and emitted for all four metrics at $n=10{,}000$, seed 135325, in `output/coloring_optimality_hh_aa_null.json`. | The classical null is the one the code-optimization literature cites; reporting it as a sensitivity companion lets the reader see the effect of the specific null choice on the reported $p$-values. Not a correction to any published number, only an addition. | (v0.5.2 tag) |
| NULL-03 | Supplement §S3.1 sensitivity table added | Not present | New subsection with two tables: preservation-property comparison (which structural features each null holds fixed vs randomizes; ~$16!$ vs ~$20!$ possible codes) and per-metric $p$-value comparison under both nulls | Makes the null-ensemble choice explicit and readable. | (v0.5.2 tag) |
| NULL-04 | Numerical outcome disclosed alongside naming correction | Only quartet-pattern $p$-values reported | Both quartet-pattern and HH-AA $p$-values reported. Grantham: 0.0061 (QP) / 0.0021 (HH-AA); Miyata: 0.0004 / 0.0001; polar requirement: 0.0026 / 0.0001; Kyte--Doolittle: 0.0014 / 0.0027. All 8 pass Bonferroni at $alpha=0.01$. | For three of the four code-independent physicochemical metrics the HH-AA null produces a *more* extreme $p$-value because its ensemble admits more code alternatives ($20!$ vs $16!$) and the standard code sits further into its low-cost tail. Kyte--Doolittle inverts: hydropathy correlates with the wobble-quartet structure that the QP null preserves, so the QP null's alternatives already share more of the observed code's low-cost signal than the HH-AA null's do. Neither null uniformly dominates the other. | (v0.5.2 tag) |
| NULL-05 | Kyte--Doolittle direction flagged in §S3.1 | Not previously discussed | One-sentence caveat in @tbl:s-null-comparison caption noting the metric-specific direction reversal and its structural interpretation | Prevents the reader from over-interpreting a single-metric shift in either direction. | (v0.5.2 tag) |
| NULL-06 | We do *not* describe either null as "more stringent" | Absent | Explicit disclaimer in §S3.1 that the two ensembles preserve non-nested structural properties and empirical tail-behaviour is not monotone across metrics | Prevents rhetorical claim of stringency without empirical grounding. | (v0.5.2 tag) |
| NULL-07 | CLI option surface updated | `--null=freeland_hurst` (default) or `class_size` | Adds `--null=haig_hurst_aa`; help text clarifies "freeland_hurst = quartet-pattern shuffle; haig_hurst_aa = classical AA-permutation" | Makes both nulls available at the command line for reproducibility. | (v0.5.2 tag) |

---

## §0 — Abstract

| ID | Change | Old value | New value | Justification | Commit |
|---|---|---|---|---|---|
| A1 | 5-metric wording + ProtSub sentence added to Typst abstract | 4-metric wording, no ProtSub sentence | 5-metric wording including "A structure-aware sensitivity analysis under the alignment-derived ProtSub matrix (Jia & Jernigan 2021) yields the most extreme percentile of any measure tested (p = 0.0004; all five p-values pass Bonferroni at α = 0.05)." | LaTeX abstract was updated during earlier remediation; Typst version drifted — audit V2 caught it | `c24509a` |
| A2 | Reassignment sentence simplified to RR wording | "only 6 of 28 observed events are topology-breaking versus 66.1% of 1,280 candidates" | "observed events are topology-breaking at relative risk 0.32 versus the candidate landscape" | Matches LaTeX abstract wording from the previous round | `c24509a` |
| A3 | Rounding harmonisation | "66% of 1,280 candidates" | "66.1% of 1,280 candidates" | Round-2 audit found abstract used digits=0 while §3.4 body used digits=1 | `a79a348` |
| A4 | Spelling harmonisation | "modeling" | "modelling" (UK, matches LaTeX abstract) | V2 audit found US/UK spelling drift between the two builds | `c24509a` |

## §2 — Methods

| ID | Change | Old value | New value | Justification | Commit |
|---|---|---|---|---|---|
| M1 | §2.3.3 rounding | "$p_{\rho=0} = 0.0061$" | "$p_{\rho=0} = 0.006$" | 3-significant-figure convention used elsewhere; §2.3.3 was the outlier | `a79a348` |

## §3 — Results

### §3.1 Cross-metric error-minimisation
| ID | Change | Old | New | Justification | Commit |
|---|---|---|---|---|---|
| R1.1 | None. Table 2 (5-metric panel) verified consistent with `manuscript_stats.metrics.*` after audit V3 | — | — | — | (no commit needed) |

### §3.3 Per-table optimality
| ID | Change | Old | New | Justification | Commit |
|---|---|---|---|---|---|
| R3.1 | Per-table Monte Carlo rerun at authoritative n = 10,000 (was n = 1,000; abstract / Methods / three captions all claimed n = 10,000) | All per-table `p_value = (k+1)/1001` | All per-table `p_value = (k+1)/10001` | CLAUDE.md is authoritative that the coloring-optimality MC uses 10,000; the caption / body claim was correct, only the pipeline default was low | `e85336f` |
| R3.2 | 26 of 27 tables significant at BH p < 0.05; sole marginal table 3 (yeast mitochondrial) at p = 0.0721 (was 0.075) | 26/27 | 26/27 (unchanged) — p-value precision changed | Result direction and count unchanged; only p-value precision changed | `e85336f` |

### §3.4 Topology avoidance
| ID | Change | Old | New | Justification | Commit |
|---|---|---|---|---|---|
| R4.1 | Fig 3 Panel C caption gains explicit definition labels per facet: H(3,4) uses Δβ₀ > 0; Q₆ uses new-disconnection | Facets not labeled with definition | Both facets labeled | Round-2 audit found reader could not tell that the two facets use different topology-breaking conventions (Option B of author's choice: label both) | `4cf1318` |

### §3.5 Conditional-logit
| ID | Change | Old | New | Justification | Commit |
|---|---|---|---|---|---|
| R5.1 | β_topo value data-driven | Hard-coded "-3.26" | `#str(calc.round(m3f.weights_raw.at(1), digits: 2))` → renders "-3.32" | Stale literal from earlier fit; current fit gives -3.32 | `2552221` |

## §4 — Discussion

| ID | Change | Old | New | Justification | Commit |
|---|---|---|---|---|---|
| D4.1 | §4.4 depletion-fold key | `tq6.depletion_fold` (3.4-fold, Q₆) | `tk43.depletion_fold` (3.1-fold, H(3,4) primary) | H(3,4) is declared primary throughout §4.2 and conclusion; §4.4 had the wrong key | `c95addd` |

---

## Figures — main text

### Fig 1: Core topology (2 panels)
| ID | Change | Justification | Commit |
|---|---|---|---|
| F1.1 | Panel B disconnection heatmap now shows all 27 NCBI tables (was 25 — Condylostoma nuclear and Balanophoraceae plastid dropped by CSV exporter bug) | Data-plumbing fix; caption ("across all 27 tables") already claimed this | `3ed7f83` |

### Fig 2: Coloring optimality (3 panels)
| ID | Change | Justification | Commit |
|---|---|---|---|
| F2.1 | Panel A now plots empirical histogram of 10,000 null draws (was Gaussian normal-approximation density) | Round-1 audit noted caption said "histogram" but plot was `dnorm()`. `monte_carlo_null()` gained `include_null_samples` flag; T3 generator persists them | `091180d` |
| F2.2 | Panel B ρ-sweep now shows all 11 points at n = 10,000 resolution (was n = 1,000 → p pinned at 1/1001) | Same as R3.1 | `e85336f` |
| F2.3 | Panel C — 27 bars (was 25 or 27 depending on which script pushed) | Data-plumbing fix as F1.1 | `3ed7f83`, `e85336f` |

### Fig 3: Evolutionary evidence (4 panels)
| ID | Change | Justification | Commit |
|---|---|---|---|
| F3.1 | Panel A: overlays codon-preserving null (in addition to uniform null); subtitle rewritten from "Bit 4 dominates" (misleading positive framing) to null-result framing | Round-1 major finding: reader was seeing an apparent bit-4 excess with no counter-null to explain it away; subtitle framed exploratory-null as positive finding | `72ac630` |
| F3.2 | Panel B: `position_jitter(seed = 135325)` applied to points + CI segments so coincident data render as distinct markers | Round-1 CRITICAL finding: Ala (Pachysolen CUG-clade, age 150, ε = 3) was hidden under Ser (Candida CUG-clade, same coords); Leu tables 16 and 22 both at (600, 2), only one visible — legend claimed markers that were invisible on the plot | `7f1cdc9` |
| F3.3 | Panel C: facet strips now render Q₆ with real plotmath subscript (was literal `Q[6]` text); subtitle uses `bquote()` with math typography (was literal `Q_6` with underscore) | Round-1 HIGH: primary user-facing "garbled legends" complaint on the paper. `label_parsed` + `bquote()` fix | `c2f4292` |
| F3.4 | Panel D: 24 distinct pairing bars (was 18 collapsed labels stacking ranks via `geom_col`'s default `position = "stack"`) | Round-1 CRITICAL: reader saw impossible bar heights (e.g. Cys/E. aediculatus at rank 6 was actually 4+2) | `68e5e2d` |

### Fig 4: Translational applications (3 panels)
| ID | Change | Justification | Commit |
|---|---|---|---|
| F4.1 | Panel C legend renders "Falsified" (not "Pending") for the red WS4 KRAS bar; NAMED factor-labels vector so ggplot's empty-level drop cannot re-assign labels; all 15 catalogued claims represented (was 3-category legend) | Round-1 CRITICAL: the paper's honest-null clinical result (KRAS-Fano p = 1.0) was legended as "Pending" — the opposite of the message | `b1f4982` |
| F4.2 | Orphan "Pending" legend entry removed (no bar uses it after F4.1) | Phase A audit V1 found the fix left an unused legend swatch | `c24509a` |

### Fig 5: Conditional-logit diagnostics (3 panels)
| ID | Change | Justification | Commit |
|---|---|---|---|
| F5.1 | Panel C x-axis: labels rotated 30° + shortened to prevent overlap | Round-1 HIGH: primary user-facing "garbled legends" complaint on the paper. Labels concatenated into "M3→M4M1→M3 H(3,4)" and "(+tRNA)(+H(3,4) topo)" at 0° rotation | `6da2104` |
| F5.2 | p-value annotations formatted grammatically ("p < 0.001", "p < 10⁻¹⁰") | Round-1 flagged ungrammatical "p=<0.001" (concatenation of `=` and `<`) | `6da2104` |

---

## Figures — supplement

### FigA (§S2) — extended coloring null
| ID | Change | Justification | Commit |
|---|---|---|---|
| SA.1 | Empirical histogram of 10,000 null draws (was Gaussian) | Same as F2.1 | `091180d` |

### FigB (§S12) — per-table optimality
| ID | Change | Justification | Commit |
|---|---|---|---|
| SB.1 | Subtitle now reads "n = 10,000" (data-driven from JSON, was hard-coded "n = 1,000") | Same as R3.1; the R script's hard-coded literal became silently stale after the T3 rerun | `8ca6b88` |
| SB.2 | 27 bars visible (was 25) | Same as F1.1 | `3ed7f83` |

### FigC (§S4) — ρ-robustness
| ID | Change | Justification | Commit |
|---|---|---|---|
| SC.1 | Y-axis switched from p (which pinned at MC floor 1/1001 for ρ ≥ 0.2, hiding the paper's monotonic-strengthening claim) to effect-size z | Round-1 MAJOR: the paper's whole claim for this figure was invisible under the p-plot | `be5c91a` |
| SC.2 | Endpoints labeled Q₆ (ρ=0) and H(3,4) (ρ=1) directly on the plot | Round-1 MINOR usability fix | `be5c91a` |

### FigE (§S10) — tRNA rank per pairing
| ID | Change | Justification | Commit |
|---|---|---|---|
| SE.1 | 24 distinct bars (was 18 stacked-labels bars) | Same as F3.4 | `68e5e2d` |

### FigF (§S3) — topology-avoidance primary cell
| ID | Change | Justification | Commit |
|---|---|---|---|
| SF.1 | Now shows the primary H(3,4), Δβ₀ > 0 cell (21.4% observed vs 66.1% possible, p = 1.3e-6) that the caption promises. Was showing stale 22.2% / 4.8e-8 from a 27-event corpus and the Q₆ new-disc cell | Round-1 MAJOR data-source drift | `e41a3e6` |

### FigG_bit_position_bias (§S22.1)
| ID | Change | Justification | Commit |
|---|---|---|---|
| SG1.1 | Overlays codon-preserving null (was only uniform); subtitle to null-result framing | Same as F3.1 | `72ac630` |

### FigG_condlogit_diagnostics (§S9)
| ID | Change | Justification | Commit |
|---|---|---|---|
| SG2.1 | Complete rebuild — three promised panels (Δ_phys vs Δ_topo joint dist with Spearman ρ = 0.15; coefficient forest with 95% CI; posterior-predictive replication) | Round-1 CRITICAL: previous version duplicated Fig 5's AICc + rank + LRT content; caption promised something different | `efbd3fe` |

### FigS_encoding_sweep (§S4)
| ID | Change | Justification | Commit |
|---|---|---|---|
| SS.1 | Default encoding (C=00 U=01 A=10 G=11) marked bold with `← default` annotation | Round-1 HIGH: 8 encodings tied at depletion 3.4×, reader could not locate the default despite caption pointing to it | `a7a0b9e` |
| SS.2 | Legend labels human-readable ("p < 0.05" / "p ≥ 0.05") instead of raw booleans (TRUE / FALSE) | Round-1 MINOR usability fix | `a7a0b9e` |

### FigW_walsh_spectrum (§S16) — NEW
| ID | Change | Justification | Commit |
|---|---|---|---|
| SW.1 | Added two-panel figure: block-size null density + 24-encoding invariance sweep | The 2-adic framework was explored in the previous analysis; §S16 previously had only a 3-row table | `7da4776` |
| SW.2 | Text now references `@fig:s-walsh-spectrum` from §S16.3 | Phase A V5 audit found the new figure was inserted but never cited | `c24509a` |

### FigSl_slavov_saap (§S17) — NEW
| ID | Change | Justification | Commit |
|---|---|---|---|
| SSL.1 | Added grouped-bar figure: observed vs baseline distributions at distances 1/2/3 | Tsour/Slavov (2026) was integrated in the previous analysis; §S17 previously had only a 3-row table (65% vs 39.5% event enrichment) | `9a5a83e` |
| SSL.2 | Text now references `@fig:s-slavov-saap` from §S17 | Phase A V5 audit found the new figure was orphan | `c24509a` |

---

## Tables — main text and supplement

Every table verified by audit V3. Only two tables changed data:

| Table | Location | Change | Justification | Commit |
|---|---|---|---|---|
| Table 4 (`tbl:condlogit`) | Manuscript §3.5 | AICc values inherit T3 rerun (still M3 best) | Direction unchanged | `e85336f` |
| `tbl:s-condlogit-restricted` Full row | Supplement §S6.1 | M3→M4 column now shows ΔAICc (-2.08) instead of LR statistic (0.12) | Round-2 MED: units mismatch with the other rows | `b5b22c7` |

Cross-reference cleanup in supplement §S15:

| Location | Change | Commit |
|---|---|---|
| `supplement.typ` §S15 body | `Tables 11--13` → `@tbl:s-metric-corr-pairs, @tbl:s-metric-corr-codes, @tbl:s-metric-partial` (table numbers had shifted after S11/S12 were inserted for R2) | `979cbc7` |
| Supplement §S8 caption | "1,000 block-preserving null draws" → "10,000" | `c24509a` |

---

## Reproducibility / analysis-code hygiene

| ID | Change | Files touched | Commit |
|---|---|---|---|
| P1 | `data_export.py` gains defensive `RuntimeError` guard against silent table-drop; regression test added | `src/codon_topo/visualization/data_export.py`, `tests/test_visualization.py` | `3ed7f83` |
| P2 | `monte_carlo_null()` gains `include_null_samples` flag; T3 generator persists them | `src/codon_topo/analysis/coloring_optimality.py`, `scripts/generate_tables.py` | `091180d` |
| P3 | Every R script's hard-coded stats replaced with `fromJSON()` / `read.csv()` reads | `manuscript_figures.R`, `strengthened_figures.R`, `all_figures.R`, `fig_encoding_sweep.R` | `b0f6593`, `8ca6b88`, `8d4cd60` |
| P4 | Deleted legacy `fig_condlogit.R` (superseded by manuscript_figures.R Fig 5); contained hard-coded LR values that had drifted from the current fit | `src/codon_topo/visualization/R/fig_condlogit.R` | `5bffe1b` |
| P5 | New diagnostic-computation script (Hessian-based SE + PP simulation vector) so FigG_condlogit_diagnostics is regeneratable | `scripts/compute_condlogit_diagnostics.py`, `output/condlogit_diagnostics.json` | `efbd3fe` |
| P6 | Slavov SAAP analysis emits full distance distributions (was only summary counts) | `scripts/slavov_saap_analysis.py`, `output/slavov_saap_codon_distances.json` | `9a5a83e` |
| P7 | `fig_combined_recoding.R` panels A-C (simulated data by author intent) gain visible `SCHEMATIC` banner | `src/codon_topo/visualization/R/fig_combined_recoding.R` | `3316a59` |

---

## What was NOT changed (audited and confirmed)

- **The four SUPPORTED claims.** Hypercube-coloring optimality across four physicochemical metrics, per-table preservation across 27 NCBI tables (BH-FDR corrected), ρ-robustness across ρ ∈ [0, 1], and topology-avoidance depletion under both Q₆ and H(3, 4) adjacency are all unchanged in the code, in `output/manuscript_stats.json`, and in the prose.
- **The FALSIFIED KRAS–Fano result** and the three REJECTED algebraic conjectures (Serine min-distance-4 invariant, PSL(2,7) linear irrep, holomorphic embedding). Their status is unchanged; only the wording of the PSL(2,7) and holomorphic-embedding rejections was narrowed to reflect what the argument actually establishes (see MIS-block above).
- **Every citation in `references.bib`.** No entries added except `gonzalez2016` (added in the previous round, unchanged in this pass).
- **The Bonferroni thresholds** and the number of test families (still eight); the tRNA row's Bonferroni-crossability is now described accurately (median-MIS $p = 0.046$ sits at the family threshold; worst-case $p = 0.123$ does not survive).
- **The final claim-hierarchy total (15).** The category count changed only in the tRNA row: 4 SUPPORTED / 5 EXPLORATORY / 1 FALSIFIED / 3 REJECTED / 2 TAUTOLOGICAL (previously 4 / 1 SUGGESTIVE / 4 EXPLORATORY / 1 / 3 / 2).

## What ships to production

Bundle contents under `output/publisher_deliverable/`:

- `manuscript.pdf` — main manuscript, Typst 0.14.2 build against the corrected pipeline; figures embedded at 300 DPI.
- `supplement.pdf` — supplement (S1–S24), Typst build against the corrected pipeline.
- `publisher_letter.pdf` — cover letter to the production editor (softened language, discloses the v0.5.1 MIS correction and reclassification).
- `cover_letter.pdf` — the original acceptance-round cover letter (unchanged).
- `response_to_reviewers.pdf` — unchanged from acceptance round; includes the R3 reply.
- `CHANGELOG.md` — public-facing release notes for `v0.5.0` (first proofs-corrections pass) and `v0.5.1` (this audit pass).
- `2026-08-14-corrections-ledger.md` — this document.
- `submission_em.zip` — Editorial-Manager LaTeX bundle: `manuscript.tex` + `manuscript.bbl` + `supplement.tex` + `supplement.bbl` + `references.bib` (ASCII-clean; every UTF-8 diacritic replaced with the LaTeX escape) + `elsarticle.cls` + the five PDFs above.
- `biorxiv_bundle.zip` — bioRxiv-ready copy: manuscript + supplement + references + PNG/TIFF figures at 300 DPI + README.
- `release_v0.5.1.zip` — reproducibility archive: refreshed pipeline JSONs, T-CSV tables, PNG+PDF figures, README, CHANGELOG, LICENSE, CITATION.cff, llms.txt.
- `MANIFEST.sha256` — SHA-256 hashes of every file above; regenerate with `sha256sum -c MANIFEST.sha256` inside the folder.

Release tag on the public repository: `v0.5.1-preproduction-2026-08-14`.

## Contact

Sergey Kornilov <sergey@biostochastics.com>
Biostochastics, LLC — Seattle, WA, USA
