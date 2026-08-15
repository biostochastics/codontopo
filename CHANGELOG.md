# CHANGELOG — CODON-TOPO

Release notes for the `codon-topo` pipeline and the accompanying
manuscript **BIOSYS-D-26-00689** (Clayworth & Kornilov, *BioSystems*, in
press). Semantic versioning; each release corresponds to a manuscript
revision or a proofs-corrections pass.

---

## [v0.6.0] — 2026-08-15 — Second pre-production QA/QC pass

Consolidated close-out of a second internal QA/QC read. All v0.5.1 corrections preserved. The four SUPPORTED claims (cross-metric coloring optimality, per-table preservation, ρ-robustness, topology-avoidance depletion) are unchanged; tRNA-enrichment remains EXPLORATORY. Substantive changes:

- **Null-model naming.** Renamed implemented null everywhere as *quartet-pattern shuffle* (16 quartet AA-patterns permuted across quartet slots); it is not the classical Haig-Hurst 1991 / Freeland-Hurst 1998 AA-permutation null. Added the classical HH-AA null as sensitivity companion (n=10 000, all four metrics; Supplement §S3.1). All eight *p*-values (4 metrics × 2 nulls) pass Bonferroni at α=0.01.
- **tRNA population.** Complete 24-row pairing/provenance table with exact 2×2 counts added as Supplement §S10.1; Fisher denominator stated as by-AA-sum over Std20. Wording corrected to "15 of 18 tRNAscan-SE-verified genomes populate the pairings; 3 (Blastocrithidia, Mycoplasmoides ×2) are mechanistic boundary cases discussed in §4.3".
- **Decision rules and framing.** §2.6 tRNA gate restored to worst-case-MIS (failed at *p*=0.123); §4.3 retitled "Exploratory tRNA accommodation"; §2.5 replaced with the discrete 4-indicator structural-preservation index that Fig 5A actually uses; per-table optimality reframed as "relative to quartet-shuffle null"; encoding-sweep claim narrowed to the two analyses actually swept; IIA acknowledged as untested; §S18 Bonferroni statement corrected; ρ-sweep 5 vs 11 grid reconciled; parametric-predictive rate stated with exact numerator/denominator/adjacency (5/66 event-step, Q₆, Δβ₀>0).
- **Data plumbing.** `manuscript_stats.condlogit.clade_exclusion` now populated (previously silently defaulted to 0, rendering "ΔAICc ≥ 0, median 0, max 0" in §3.5); "Static fallback table" branches deleted and replaced with `#assert` so future gaps error visibly at build.
- **Bibliography and reproducibility.** Kosiol 2004 and second Yurova Axelsson-Khrennikov 2024 paper now cited; reproducibility statement downgraded to "within numerical/rendering tolerance"; `requirements.lock` and `Dockerfile` added.
- **Rendering.** Typst source-escape defects fixed (`10^{-5}`, `1{,}280`, two `10{,}000`, `{≈}5×`, `{≈}100%`, `e.g._T. thermophila_`, `Table @tbl:s-claims` doubling, Chan & Lowe DOI); cross-reference sweep.
- **Release engineering.** Version 0.5.1 → 0.6.0 propagated across `pyproject.toml`, `CITATION.cff`, `requirements.lock`, `Dockerfile`, `llms.txt`, manuscript and supplement reproducibility stamps. Ruff format + check clean. Full pipeline rerun at seed 135325. Tag `v0.6.0`; SHA-256 manifest.

Reproduction identical to v0.5.1 (same three-command sequence), with `git checkout v0.6.0` in place of the v0.5.1 tag.

---

## [v0.5.1] — 2026-08-14 — First pre-production QA/QC pass

Second audit ahead of publisher proofs. One computational correction,
plus a broad presentation, cross-reference, and reproducibility sweep.

### Computational correction (disclosed)

The maximal-independent-set (MIS) enumeration used for the tRNA-
enrichment robustness bound in Supplement §S10.2 was applying
Bron–Kerbosch clique-branching pivoting to a graph whose recursion was
set up for independent-set search. Correcting the algorithm (standard
Bron–Kerbosch on the complement of the conflict graph) changes the
enumerated MIS count from 2 to 332 (all of size 6) and expands the
reported distribution: median Stouffer p = 0.046, best 0.016, worst
0.123; 190 of 332 sets fall below 0.05, 0 below 0.01. A pre-specified
topology-breaking-only subset (n = 4 pairings) that was already computed
in the repository is null (Stouffer p = 0.43).

**Effect on scientific claims.** The all-pairings Fisher–Stouffer
result (p = 1.7 × 10⁻⁷) and the four SUPPORTED claims (hypercube
coloring optimality, per-table preservation, ρ robustness, topology-
avoidance depletion) are unchanged. The tRNA-enrichment claim is
reclassified from SUGGESTIVE to EXPLORATORY: the median-independent-
subset signal is present but the strict worst-case bound is not met and
the topology-breaking-only subset is null. The claim hierarchy now
enumerates 4 SUPPORTED / 5 EXPLORATORY / 1 FALSIFIED / 3 REJECTED / 2
TAUTOLOGICAL claims (unchanged total of 15).

### Manuscript, supplement, and README

- Abstract, Introduction Test 3 criterion, §3.6, Table 7 & caption,
  Discussion, Limitations, and Conclusion — tRNA text rewritten with
  the corrected MIS distribution and the topology-breaking-subset
  disclosure.
- Table 1 (claim hierarchy) — tRNA row moved from SUGGESTIVE to
  EXPLORATORY, p updated (MIS median 0.046). Shaded-rows caption
  reworded so FALSIFIED / REJECTED / TAUTOLOGICAL rows are described as
  non-inferential rather than "failed validation".
- Cross-reference sweep — main-manuscript pointers to §S13 / §S15 /
  §S16 / §S17 / §S18 updated to the current supplement numbering
  (§S13 / §S14 / §S15 / §S16 / §S17 / §S18 / §S20 / §S21 / §S22 / §S23).
- §5.10 H(3, 4) clade-exclusion paragraph — p-value now quoted from
  the matching primary-cell calculation (H(3, 4), Δβ₀ > 0, N = 1280,
  K = 846, n = 24, x = 2 → p ≈ 4.2 × 10⁻⁹), with the Q₆ /
  new-disconnection sensitivity variant (≈ 3.6 × 10⁻¹¹) called out
  separately.
- Figure 4 caption — six *candidate* conditional-logit models (LR
  tests applied only to nested pairs); AICc gaps itemised (M1 / M2
  baselines 89.1–112.1 vs M3; H(3, 4) variant 16.9; M4 2.1).
- §S12 — feasibility-score body corrected to describe the discrete
  four-indicator implementation (S ∈ {0.55, 0.75, 0.80, 1.00}), and
  the section renamed to "Structural-preservation index".
- Table S8 caption — discloses the tRNAscan-SE two-pass distinction
  (Infernal-filtered `Supp` column vs first-pass Isotype/Anticodon
  totals in the Reassigned-AA parenthetical).
- "posterior-predictive" renamed to "parametric predictive" throughout
  (the estimation is MLE + AICc, not Bayesian).
- PSL(2,7) rejection narrowed to "64-dimensional irreducible linear
  representation"; holomorphic-embedding rejection grounded in the
  character-identity failure rather than in an incorrect appeal to
  domain type.
- README — tRNA-enrichment section rewritten with the corrected
  distribution.

### Code

- `src/codon_topo/analysis/trna_evidence.py::_bron_kerbosch` — pivoting
  applied to the complement of the conflict graph via `non_conflicts`
  adjacency; enumerates all 332 MIS of size 6.
- `src/codon_topo/reports/claim_hierarchy.py` — tRNA-enrichment entry
  restated with the new distribution and reclassified EXPLORATORY.
- `src/codon_topo/cli.py`
  - `manuscript_stats.json.trna` gains `mis_median_p`, `n_mis_total`,
    `mis_frac_significant_p05`, `topology_breaking_p`,
    `topology_breaking_z`, `topology_breaking_n`.
  - `manuscript_stats.json.topology_avoidance_q6` gains
    `depletion_fold`, `risk_ratio`, `risk_ratio_ci_95` (matching the
    tk43 shape).
  - `_write_json` now sanitises NaN / Inf to `null` before serialisation
    and uses `allow_nan=False`, so every JSON parses under strict
    readers (R jsonlite, Typst `json()`).
- `scripts/generate_tables.py::T7` — summary CSV emits MIS-median,
  MIS-worst-case, and topology-breaking-subset rows alongside the
  existing Fisher and AA-label rows.

### Delivery

- Rebuilt manuscript.pdf, supplement.pdf, publisher_letter.pdf, and
  the SHA-256 manifest under `output/publisher_deliverable/`.
- Release tag: `v0.5.1-preproduction-2026-08-14`.

### Reproduction

```
git checkout v0.5.1-preproduction-2026-08-14
pip install -e ".[dev]"
codon-topo all --output-dir=./output --seed=135325 --n=10000
python3.11 scripts/generate_tables.py
Rscript src/codon_topo/visualization/R/manuscript_figures.R
Rscript src/codon_topo/visualization/R/strengthened_figures.R
Rscript src/codon_topo/visualization/R/fig_encoding_sweep.R
Rscript src/codon_topo/visualization/R/fig_condlogit_diagnostics.R
Rscript src/codon_topo/visualization/R/fig_walsh_spectrum.R
Rscript src/codon_topo/visualization/R/fig_slavov_saap.R
typst compile output/manuscript.typ output/manuscript.pdf
typst compile output/supplement.typ output/supplement.pdf
```

Seed = 135325 across every Monte Carlo path.

---

## [v0.5.0] — 2026-08-14 — Proofs-corrections pass

First internal figure- and statistics-audit pass. Resolved ~41 first-
round + ~30 second-round presentation and data-plumbing defects across
the 16 figures, 25 tables, and body prose of the accepted manuscript
and supplement. All corrections were presentation, plumbing, or
reproducibility improvements (no scientific claim changed at v0.5.0;
the tRNA-classification change was identified later and shipped as
v0.5.1 above).

Highlights:

- CSV exporter for the disconnection catalogue no longer silently drops
  NCBI translation tables 28 (Condylostoma nuclear) and 32
  (Balanophoraceae plastid); Fig 1B and supplement FigB gain the two
  missing rows.
- Coloring-optimality per-table Monte Carlo and ρ-sweep reran at
  n = 10,000 (the abstract / Methods / three figure captions had
  claimed that count; the code had been at 1,000). The 26/27
  significance count and the sole marginal (table 3 yeast mito at
  p = 0.072) are unchanged.
- `output/tables/T9_topology_avoidance.csv` now emits the primary
  H(3, 4), Δβ₀ > 0 cell rather than the Q₆/new-disconnection cell.
- Fig 3 depth-calibration panel — reproducible-seed jitter so
  coincident CUG-clade markers are all visible.
- Fig 3 tRNA-rank panel and supplement FigE — labels now include the
  control organism so all 24 pairings render as 24 distinct bars.
- Fig 4 catalogue panel — named-vector label mapping; KRAS bar now
  correctly labelled "Falsified" (was showing "Pending").
- Two new supplement figures: Walsh–Hadamard / 2-adic spectral panel
  (§S21) and Slavov / Tsour SAAP grouped-bar figure (§S22).
- Reproducibility: every R visualization script now reads statistics
  from `manuscript_stats.json` / per-analysis JSONs / T-CSVs rather
  than embedding hard-coded literals.

---

## [v0.4.x] — earlier — Additional bridging analyses

Additional analyses: ProtSub sensitivity, Walsh–Hadamard §S21, non-power-
representation discussion in §4.5, Slavov / Tsour SAAP §S22, metric
correlation analysis §S20. See git log for per-commit detail.
