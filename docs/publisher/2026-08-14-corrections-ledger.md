# BIOSYS-D-26-00689 (BIO 105919) — Pre-publication QA/QC Ledger

**Manuscript:** *Robust error-minimization in the genetic code across
physicochemical metrics and variant codes: a graph-theoretic analysis in
GF(2)⁶*
**Authors:** Paul Clayworth, Sergey Kornilov
**Status:** In press, *BioSystems* (Elsevier). Article reference
BIO_BIOSYS-D-26-00689; production reference BIO 105919.
**Prepared:** 15 August 2026
**Corresponding author:** Sergey Kornilov <sergey@biostochastics.com>
**Repository:** https://github.com/biostochastics/codontopo
**Release tag:** `v0.6.1` — the single consolidated release accompanying
this ledger.

## Purpose of this ledger

Between acceptance and typesetting the authors carried out a
pre-publication QA/QC audit of every figure, table, and inferential
quantity across the manuscript and the supplement. This ledger is the
diff manifest for that audit. It records every user-visible change,
sorted by manuscript location, with the pre-fix value, the post-fix
value, and the justification.

**All headline scientific conclusions of the accepted manuscript are
unchanged.** Four SUPPORTED claims — cross-metric coloring optimality,
per-table preservation, ρ-robustness, and topology-avoidance depletion —
the EXPLORATORY tRNA-enrichment classification, the FALSIFIED KRAS–Fano
prediction, and the three REJECTED / two TAUTOLOGICAL entries are as
reported in the accepted version. The claim hierarchy count remains 15
(4 SUPPORTED / 5 EXPLORATORY / 1 FALSIFIED / 3 REJECTED / 2 TAUTOLOGICAL).

Two substantive computational corrections were made during the audit and
their effects on individual displayed values are documented explicitly
below (§MIS and §NULL). Several displayed p-values shifted at the third
and fourth decimal as a result; none crossed a decision threshold. In
addition, a number of clarifications, framing corrections, and editorial
fixes were applied across the manuscript, supplement, response letter,
CHANGELOG, and public README, and the delivery inventory was regenerated
from a single source-of-truth.

---

## §MIS — Maximal-Independent-Set enumeration for tRNA enrichment

**Root cause.** The initial Bron–Kerbosch pivoting was applied to the
conflict graph *G* rather than to its complement *Ḡ*. An independent set
in *G* is a clique in *Ḡ*, so the pivot must be chosen in *Ḡ*. Fixing
the pivot produced the complete MIS enumeration below.

| ID | Change | Old value | New value | Justification |
| --- | --- | --- | --- | --- |
| MIS-01 | Bron–Kerbosch pivot moved to complement graph | Enumerates 2 MIS of size 6 (worst-case Stouffer *p* = 0.045) | Enumerates 332 MIS of size 6; median *p* = 0.037, best 0.012, worst 0.104; 264/332 (79.5%) sets below 0.05 | An independent set in *G* is a clique in *Ḡ*; pivot must be selected in *Ḡ*, not in *G*. Present numbers reflect the post-audit tRNA vector reconciliation (§MIS-11) applied to the complete MIS enumeration. |
| MIS-02 | Topology-breaking-subset test surfaced | Present in `output/trna_evidence.json` but not cited in main text | Explicitly reported in §3.6 and §S10.3: *n* = 4 pairings (yeast mito Thr, *Scenedesmus* mito Leu, *Pachysolen* Ala, *Candida* Ser), Stouffer *p* = 0.387 (null) | The subset is the direct mechanistic test of the compensation-for-disconnection hypothesis; suppressing it would misrepresent the mechanistic evidence. |
| MIS-03 | tRNA-enrichment claim reclassified | Status: SUGGESTIVE, evidence *p* = 0.045 (worst-case MIS) | Status: EXPLORATORY, evidence median MIS *p* = 0.037, worst 0.104; topology-breaking subset *p* = 0.387 (null) | The strict worst-case criterion is not met, and the mechanism-linked subset is null; the signal is retained as heterogeneous accommodation-of-decoding rather than compensation-for-disconnection. |
| MIS-04 | Main-text Table 7 templated from stats | All rows literal or partly templated | All five rows (Fisher+Stouffer, median MIS, worst MIS, best MIS, topology-breaking subset) render from `stats.trna.*` fields | Prevents future drift between prose and table; the same fields populate the corresponding supplement rows and Table S1 claim justification. |
| MIS-05 | Main-text §3.6 tRNA paragraph rewritten | Two-point summary and stale numbers | Distribution across 332 MIS of size 6 with median/best/worst and 264/332 fraction; topology-breaking-subset null called out | Presents the distribution rather than a two-point summary; discloses the null subset explicitly. |
| MIS-06 | Supplement §S10.2 (MIS) rewritten | Two MIS reported | 332 MIS reported; best/worst pairings identified; topology-breaking-subset subsection added | Same as MIS-05 for the supplement. |
| MIS-07 | tRNA legacy repertoires reconciled against primary literature (per pre-publication audit) | *S. cerevisiae* mitochondrial and *Yarrowia lipolytica* mitochondrial counts inconsistent with primary sources | *S. cerevisiae* mitochondrial: 24 tRNAs (2 Arg tRNAs decoding AGR + CGN) after Bonitz *et al.* 1980 (PNAS 77:3167). *Y. lipolytica* mitochondrial: 27 tRNAs (no CGN-Arg) after Kerscher *et al.* 2001 (Comp Funct Genomics 2:80, mould-mitochondrial table 4). | Primary-source anchoring; both references added to the manuscript and supplement bibliographies (§Bibliography). |
| MIS-08 | Main-text §2.3.6 tRNA analysis population paragraph rewritten | Named "5 variant codes (tables 4, 6, 10, 15, 31)" and said *S. cerevisiae* was "not included in the primary 24-pairing analysis" | Correctly names primary tables 3 (yeast mito CUN→Thr), 6, 10, 12, 15, 22 (chlorophycean mito equivalent of table 16), and 26; names *Blastocrithidia nonstop* + 2 *Mycoplasmoides* species as mechanistic boundary cases excluded from the 24-pairing analysis; states that *S. cerevisiae* mito Thr appears in every MIS | The prior text was inconsistent with the 24-row provenance table and with the MIS enumeration; the corrected paragraph is generated from the same 24-row input. |
| MIS-09 | Main-text §4.3 signal-driver list corrected | Said signal driven by ciliates, *Euplotes*, "*Blepharisma* and *Mycoplasmoides*" | Signal driven by UAR→Gln ciliates, UGA→Cys *Euplotes*, and *Blepharisma stoltei* UGA→Trp; *Blastocrithidia* and *Mycoplasmoides* are excluded boundary cases (single-tRNA anticodon-stem or base-modification routes without gene-copy expansion) | *Mycoplasmoides* is outside the 24-pairing set precisely because its mechanism does not produce a copy-number signal; naming it as a driver reversed the evidentiary role. |
| MIS-10 | Supplement Table S11 gains a companion row-level provenance table | Only analysis-input columns visible | New Table S11b (Row-level provenance) with compartment, NCBI table id, source class, and assembly accession or primary citation for each of the 48 organism-slots; machine-readable CSV `output/tables/T-trna-24row-provenance.csv` extended with `variant_source_full` and `control_source_full` columns | Aggregate source-class counts alone were insufficient to audit each 2×2 table; the audit required row-level provenance. |
| MIS-11 | Limitations, Conclusion, and README | "(suggestive)"; "*p* = 0.045" | "(exploratory)"; "median MIS *p* = 0.037, worst 0.104, topology-breaking subset *p* = 0.387" | Same reclassification as MIS-03, propagated to every venue that quotes the tRNA numbers. |

---

## §NULL — Coloring-optimality null model

**Context.** The initial implementation labelled the primary
coloring-optimality null "Freeland–Hurst block-preserving". The
implemented ensemble is a *quartet-pattern shuffle* that holds each
mobile first-two-base quartet AA-pattern atomic and permutes them across
the 14 mobile quartet slots. The classical Haig–Hurst 1991 /
Freeland–Hurst 1998 null described in the code-optimization literature
is a different ensemble: it holds the 20 sense-codon families fixed and
permutes the 20 amino-acid labels uniformly across the 20 families. The
audit identified this as a naming discrepancy rather than a computational
error — the quartet-pattern shuffle is a valid stringent null, but it
should not be named after Haig and Hurst — and separately corrected the
Monte Carlo tail convention.

**Preservation properties.** The quartet-pattern shuffle preserves each
quartet's internal four-position amino-acid split shape (whether it is a
(4), (2, 2), (3, 1), or singleton block) and each named amino acid's
total codon count; it permutes which codons carry those labels. It does
NOT fix the codon-family partition. The classical Haig–Hurst
amino-acid permutation DOES fix the codon-family partition; it does NOT
fix per-named-AA codon count. The state-space count for the
quartet-pattern shuffle is ≈ 14! ≈ 8.7 × 10¹⁰ (14 non-stop quartets are
mobile; the 2 stop-containing quartets are held fixed), not ≈ 16! as
previously stated.

**Monte Carlo tail convention.** `monte_carlo_null` and its three
sibling functions previously counted null draws with *f* < *F*_obs
rather than the standard lower-tail permutation-test convention
*f* ≤ *F*_obs (i.e., they excluded ties at the observed value). The
Grantham null contains exactly one draw equal to *F*_obs = 13 477;
conservative *p*_cons moves from 61/10 001 = 0.006099 (strict) to
62/10 001 = 0.006199 (ties included). The other four metrics have no
ties and are unchanged. All eight quartet-pattern-shuffle + Haig–Hurst
*p*-values pass Bonferroni at α = 0.05/8 = 0.00625.

| ID | Change | Old value | New value | Justification |
| --- | --- | --- | --- | --- |
| NULL-01 | Renamed implemented null throughout main text and supplement | "Freeland–Hurst block-preserving null" | "quartet-pattern shuffle null" with an explicit naming/method paragraph in §2.3.1 and full comparison in §S2.1 | Attribution accuracy: the implemented ensemble is not the classical AA-permutation null. |
| NULL-02 | Added classical Haig–Hurst 1991 / Freeland–Hurst 1998 AA-permutation null | Not implemented | New generator `_generate_random_code_haig_hurst_aa_permutation` in `src/codon_topo/analysis/coloring_optimality.py`; wired via `null_type="haig_hurst_aa"` in `monte_carlo_null()`; emitted for all four metrics at *n* = 10 000, seed 135325, in `output/coloring_optimality_hh_aa_null.json` | Companion sensitivity ensemble; lets readers see the effect of the null-ensemble choice on displayed *p*-values. |
| NULL-03 | Supplement §S2.1 preservation-property table added | Not present | New comparison of what each null holds fixed vs. randomises, with per-metric *p*-value comparison and state-space counts (14! vs 20!) | Makes the null-ensemble choice explicit and readable. |
| NULL-04 | Monte Carlo tail convention corrected (audit surface) | `k` = number of null draws with `f < F_obs` in four `monte_carlo_null()` functions | `k` = number of null draws with `f ≤ F_obs` (ties at the observed value included); `p_conservative = (k+1)/(n+1)`; result dict gains `n_at_or_below_observed`, `n_strictly_below_observed`, `n_ties_at_observed` fields | Standard lower-tail permutation-test convention. Only Grantham has a tie at *F*_obs; other metrics unchanged. |
| NULL-05 | Main-text §2.3.1 methods rewritten for tail convention | "the number of null scores *below* the observed *F*" | "the number of null scores satisfying *F*_null ≤ *F*_obs; ties at the observed score are included in *k*" | Same convention as the corrected implementation. |
| NULL-06 | Main-text §2.3.3 ρ-sweep literal updated | *p*_ρ=0 = 0.006, *k* = 60 | *p*_ρ=0 = 0.0062, *k* = 61 (60 strict + 1 tie); other ρ values unchanged (no ties) | Same convention. |
| NULL-07 | Main-text Fig 2A caption resolution increased | Caption rendered *p* at 3 significant digits (0.006) | Caption renders at 4 significant digits (0.0062), matching the figure's baked-in label | Prevents caption/figure mismatch given the Bonferroni threshold 0.00625 is close to the observed value. |
| NULL-08 | Abstract and §3.1 Grantham *p* literals updated | Grantham *p* = 0.006 | Grantham *p* = 0.0062 (four-digit display; important because the Bonferroni-corrected threshold at α = 0.05/8 = 0.00625 is close) | Editorial finding E13. |
| NULL-09 | CLI option surface updated | `--null=freeland_hurst` (default) or `class_size` | Adds `--null=haig_hurst_aa`; help text clarifies "freeland_hurst = quartet-pattern shuffle; haig_hurst_aa = classical AA-permutation" | Makes both nulls reproducible from the command line. |
| NULL-10 | §S18 cross-family multiplicity narrative rewritten | Placed raw *p*-values, within-family BH-FDR results, permutation *p*-values, and AICc gaps in one cross-family Bonferroni list | Recast as a descriptive conservative-threshold sensitivity: primary/smallest *p*-values reported against a common conservative reference threshold α* = 0.05/8 = 6.25 × 10⁻³; conditional-logit contribution reported as the M1→M3 likelihood-ratio test (χ² = 110.4 on 1 df, *p* ≪ 10⁻²⁰) with the accompanying ΔAICc reported as a descriptive model-selection statistic; within-family corrections retained | AICc is not a *p*-value and cannot "survive" an α threshold; the reported quantities are non-exchangeable and do not form a single family-level list. |

---

## §H34 — H(3,4) clade-exclusion overclaim (main-text §3.4)

| ID | Change | Old value | New value | Justification |
| --- | --- | --- | --- | --- |
| H34-01 | Main-text §3.4 clade-exclusion global threshold | "The depletion is highly significant in every regime (*p* < 10⁻⁵)" | "Under the primary *H*(3,4)/Δβ₀ > 0 definition all seven exclusions satisfy *p* < 10⁻³ (largest *p* ≈ 2.03 × 10⁻⁴, from excluding all metazoan mitochondria); under the sensitivity *Q*₆/new-disconnection definition all seven satisfy *p* < 10⁻⁵" | Under the primary *H*(3,4) adjacency, 2 of 7 regimes exceed 10⁻⁵ (ciliates 3.96 × 10⁻⁵, metazoan mitochondria 2.03 × 10⁻⁴); the *p* < 10⁻⁵ statement is true only for the *Q*₆/new-disconnection sensitivity. |
| H34-02 | Table S1 (claim hierarchy) topology-avoidance row aligned | "Clade-exclusion sensitivity (7 regimes): all *p* < 10⁻⁵" | "Clade-exclusion sensitivity (7 regimes): all *p* < 10⁻³ under *H*(3,4)/Δβ₀ > 0 (largest *p* ≈ 2.03 × 10⁻⁴); all *p* < 10⁻⁵ under *Q*₆/new-disconnection" | Same fix in the registered claim hierarchy. |
| H34-03 | Supplement §S9 (H(3,4) clade-exclusion table) retained | 7-row table already present | Referenced from the corrected main-text §3.4 statement | Provides the underlying values that support the corrected statement. |

---

## §PARTIAL — Metric partial-correlation count

| ID | Change | Old value | New value | Justification |
| --- | --- | --- | --- | --- |
| PC-01 | Main-text §3.1 count | "five of ten pairs become approximately uncorrelated after adjustment (\|ρ_partial\| < 0.1)" naming Grantham~PR, PR~KD, PR~ProtSub | "four of ten pairs" naming Grantham~PR (+0.03), Miyata~ProtSub (−0.06), PR~KD (+0.01), PR~ProtSub (+0.08); the descriptor "a descriptive point-estimate statement, not a test of statistical independence" added | Table S19 lists exactly four pairs with \|ρ\| < 0.1; main-text and CHANGELOG count corrected. |
| PC-02 | CHANGELOG partial-correlation line | "five of ten" | "four of ten" | Same. |
| PC-03 | README §Coloring optimality | "five of ten pairs are statistically independent" | "four of ten pairs are near-zero after adjustment … (a descriptive point-estimate statement)" | Same. |

---

## §NOTATION — Topology-graph notation and quartet-null description

| ID | Change | Old value | New value | Justification |
| --- | --- | --- | --- | --- |
| N-01 | Main-text §2.3.4 graph notation | "let *G*_a^1 denote the subgraph induced on codons assigned to *a* with edges between codons at binary Hamming distance ≤ 1" (Q₆-only) | "for each amino acid *a* and each adjacency *A* ∈ {Q₆, H(3,4)}, let *G*_A^a denote the subgraph induced on codons assigned to *a* with edges given by adjacency *A*; the topology-breaking indicator Δβ₀(*G*_A^a) > 0 is defined separately for each adjacency; the primary test uses *A* = H(3,4)" | Prior notation defined only the Q₆ variant while declaring H(3,4) primary; the corrected definition indexes the graph by adjacency. |
| N-02 | Main-text §2.3.1 quartet-null description | "This null preserves both the synonymous codon contiguity (wobble degeneracy) and each quartet's internal split structure" | "The null preserves each quartet's internal four-position amino-acid split shape … and each named amino acid's total codon count, while permuting *which* codons carry those labels; no quartet is broken up, but amino-acid labels move between quartet slots" | The prior phrasing suggested named-amino-acid to codon-set membership was fixed; the exact preservation properties are as stated. |
| N-03 | Main-text §2.1 encoding-sweep sentence | "All encoding-dependent results are tested across all 24" | "The two encoding-dependent analyses in this work (cross-metric coloring optimality and the *Q*₆ topology-avoidance depletion) are evaluated under all 24 bijections; the M1–M4 conditional-logit models and the primary *H*(3,4) topology-avoidance test are encoding-invariant by construction and are not swept" | Prior sentence overstated the scope of the encoding sweep; only two encoding-dependent analyses are swept. |
| N-04 | Abstract encoding-robustness wording | "robust to alternative topology definitions and base-to-bit encodings" | "The *H*(3,4) result is stable by construction; the *Q*₆ decomposition is representation-specific and fails to show depletion under 8 of 24 base-to-bit encodings, so we report *H*(3,4) as the primary test and *Q*₆ as a sensitivity" | Distinguishes the encoding-independent *H*(3,4) result from the encoding-dependent *Q*₆ decomposition. |

---

## §FRAMING — Discussion framing

| ID | Change | Old value | New value | Justification |
| --- | --- | --- | --- | --- |
| F-01 | Main-text §4.2 per-table framing | "26 of 27 tables … suggesting that codon reassignment events are constrained to preserve error-minimization" | "26 of 27 tables remain low-cost relative to their own quartet-shuffle ensembles. Because the per-table null lacks distance-matched draws (§S8), this result is *compatible* with preservation of error-minimization during reassignment but does not identify independent variant-specific optimization" | The per-table null lacks distance-matched draws; the result is compatible with preservation but does not evidence independent optimization. |
| F-02 | Supplement §S6 (IIA) framing | "the explanatory–rather–than–predictive framing makes IIA tolerable for our purposes" | "We use the model for explanatory comparison, not as an unconditional endorsement of IIA: IIA remains an untested structural assumption, the restricted-candidate refits are sensitivity analyses rather than evidence that IIA is harmless, and coefficients and likelihoods could shift under a mixed-logit relaxation" | Aligns §S19.1 with the more careful §S6 caveat that IIA is untested. |
| F-03 | Supplement §S22 Slavov event-level framing | Event-weighted binomial *p* < 10⁻¹⁰⁰ presented as an "empirical signature" | Frequency contrast (65.0% vs 39.5%) framed as the primary descriptor; the event-weighted binomial *p*-value is labelled a descriptive naive summary that ignores clustering by substitution type, gene, sample, and measurement opportunity; the substitution-type analysis (OR 1.17, *p* = 0.26) is placed immediately beside it | Harmonises the supplement with the main text's more careful phrasing. |
| F-04 | Main-text Limitations 22-pair sentence | "Assembly-fragmentation robustness (removing two genomes with >10% pseudogene fraction yields Stouffer *p* = 0.041 on the remaining 22 pairings) is reported in Supplement §S10" | "Per-genome tRNAscan-SE pseudogene fractions vary across the 18 verified assemblies and are reported in Supplement §S10 (Table S12, *Pseudo* column); we do not report a formal fragmentation-based sensitivity test in this manuscript" | The prior sentence pointed to §S10 for a fragmentation-based sensitivity that the supplement did not document; the corrected sentence removes an unsupported claim. |

---

## §CLAIM — Claim-hierarchy definitions

| ID | Change | Old value | New value | Justification |
| --- | --- | --- | --- | --- |
| C-01 | Main-text §2.6 SUPPORTED definition | "SUPPORTED (passes rigorous null, *p* < 0.01)" | "SUPPORTED (meets its prespecified claim-specific decision rule under the stated null or model and multiplicity control); EXPLORATORY (descriptive or hypothesis-generating; does not meet the relevant confirmatory gate); REJECTED; FALSIFIED; TAUTOLOGICAL — Each claim's decision rule is stated in the Justification column of the full hierarchy in the supplementary materials" | The prior blanket *p* < 0.01 rule did not fit the BH-FDR-based per-table claim; a per-claim decision rule is used. |
| C-02 | `claim_hierarchy.py` SUPPORTED docstring | "Passes rigorous null; cite as finding" | "Meets its prespecified claim-specific decision rule under the stated null or model and multiplicity control" | Same. |
| C-03 | Fig 5C caption + status labels | Labels "Verified / Tested / Pending / Falsified / Tautological / Exploratory" without cross-reference to Table S1 | Labels renamed to registered categories where possible ("Supported / Rejected / Pending / Falsified / Tautological / Exploratory"); caption explains the per-prediction vs per-claim distinction | Aligns Figure 5C with the registered claim hierarchy. |
| C-04 | Fig 5A axis label | "Feasibility score" | "Structural-preservation index *S* (visualization-only)" | Aligns Fig 5A with the Methods §2.5 rename to "Structural-preservation index (visualization-only)". |

---

## §BIB — Bibliography additions and citation cleanups

| ID | Change | Justification |
| --- | --- | --- |
| B-01 | Added `@bonitz1980` (Bonitz *et al.* 1980, PNAS 77:3167, doi:10.1073/pnas.77.6.3167) to `output/references.bib`; cited in §2.3.6 as the primary source for the *S. cerevisiae* mitochondrial tRNA vector. | Primary source for the corrected yeast-mito tRNA counts. |
| B-02 | Added `@kerscher2001` (Kerscher *et al.* 2001, Comp Funct Genomics 2:80, doi:10.1002/cfg.72) to `output/references.bib`; cited in §2.3.6 as the primary source for the *Y. lipolytica* mitochondrial tRNA vector. | Primary source for the corrected *Y. lipolytica* mito tRNA counts. |
| B-03 | Added explicit `@atchley2005` in-text citation at the Serine/Factor 3 exploratory-observation paragraph (main-text §3.7). | Prior text used *F*_3 = -4.760 without citing the Atchley source. |
| B-04 | Retained existing `@chan2016` GtRNAdb citation and referenced it at the GtRNAdb-sourced organism-slots in §2.3.6. | Replaces bare URL. |

---

## §AI — Generative-AI-assisted code disclosure

| ID | Change | Justification |
| --- | --- | --- |
| AI-01 | New Methods subsection §2.6 "Software, validation, and reproducibility" added, describing Python 3.11 / R 4.4 stack, the use of generative-AI systems to scaffold portions of analysis and figure-generation code, human review, version control, unit / regression / property-based test coverage (>400 pytest tests), independent numerical anchor-value checks, and confirmation that no AI system generated or altered research data, inferential outputs, or figures. | Elsevier's current journal policy calls for Methods-level disclosure of AI-assisted code development in addition to the end-of-manuscript declaration. |

---

## §FIG — Figure fixes carried into this release

| ID | Change | Justification |
| --- | --- | --- |
| FIG-01 | Figure 3B: added deterministic pre-jitter so that coincident (age, ε) markers render as visually distinct points in both PNG and PDF renders | Prior versions had visible jitter in the PDF render but not the PNG; the fix makes both consistent. |
| FIG-02 | `save_figure()` in `theme_codon.R` now emits PNG + PDF + TIFF (LZW-compressed) for every figure. | Elsevier production requires TIFFs for production upload. |

---

## §DELIVERY — Delivery inventory (unified)

The single authoritative delivery inventory for this release, sorted by
recipient use. Every entry is hashed in `MANIFEST.sha256`.

| File | Purpose |
| --- | --- |
| `manuscript.pdf` | Main manuscript, Typst build; figures embedded at 300 DPI; PDF metadata (title, authors, keywords, date) written to the Info dictionary. |
| `supplement.pdf` | Supplement, Typst build; PDF metadata written to the Info dictionary. |
| `elsevier_response.pdf` | Response letter accompanying this delivery. |
| `corrections_ledger.pdf` | This document, rendered as a PDF. |
| `submission_bundle.zip` | Editorial-Manager-ready bundle: `manuscript.tex` + `supplement.tex` (bibliography inlined via `\begin{thebibliography}` for self-contained pdflatex compilation) + `elsarticle.cls` + 300-DPI PNGs and TIFFs for every figure + `references.bib` + the compiled `manuscript.pdf` / `supplement.pdf`. |
| `submission_em_v0.6.1/` | Unpacked source of `submission_bundle.zip` for reviewers preferring an unarchived layout. |
| `CHANGELOG.md` | Release notes: consolidated v0.6.1 entry for the pre-publication QA/QC pass. |
| `2026-08-14-corrections-ledger.md` | Markdown source of this ledger. |
| `MANIFEST.sha256` | SHA-256 hashes of every file above. Verify with `shasum -a 256 -c MANIFEST.sha256` inside the folder. |

Release tag on the public repository: `v0.6.1` (annotated).

---

## §BUILD — Canonical build command

One canonical command reproduces every displayed value from a clean
clone of the public repository at tag `v0.6.1`:

```bash
git clone https://github.com/biostochastics/codontopo && cd codontopo
git checkout v0.6.1
pip install -e ".[dev]"
bash scripts/build_publisher_release.sh
```

`scripts/build_publisher_release.sh` runs the analysis pipeline
(`codon-topo all --seed=135325 --n=10000`), the provenance-emit and
table-generation scripts, both R figure scripts, and the two Typst PDF
compilations, in that order. The same four-line command block is
referenced verbatim in the Elsevier response letter, `CHANGELOG.md`,
main-text Data-and-code availability, supplement §S24, and `README.md`.
Python dependency versions are pinned in `requirements.lock`; a
reference `Dockerfile` provides a container floor.

---

## §UNCHANGED — What was audited and NOT changed

- **The four SUPPORTED claims.** Cross-metric coloring optimality (all four code-independent metrics below α*, Grantham *p* = 0.0062 after tie correction), per-table preservation (11 of 12 informative-distance variants below BH-FDR α = 0.05), ρ-robustness (*p* ≤ 0.006 at every ρ), and *H*(3,4) topology-avoidance depletion (hypergeometric *p* = 1.28 × 10⁻⁶; RR ≈ 0.32) all retain SUPPORTED status.
- **The FALSIFIED KRAS–Fano result** and the three REJECTED algebraic conjectures (Serine min-distance-4 invariant, PSL(2, 7) linear irrep, holomorphic embedding). Status unchanged; wording narrowed to reflect what the argument establishes.
- **The final claim-hierarchy total (15).** Category count: 4 SUPPORTED / 5 EXPLORATORY / 1 FALSIFIED / 3 REJECTED / 2 TAUTOLOGICAL, as in the accepted version.
- **All figure captions and legends** except those explicitly listed in §FIG.
- **All bibliography entries** except those in §BIB.

---

## Contact

Sergey Kornilov <sergey@biostochastics.com>
Biostochastics, LLC — Seattle, WA, USA
