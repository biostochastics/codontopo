# codon-topo v0.6.0 — 2026-08-15

**Consolidated pre-production release for BIOSYS-D-26-00689** (Clayworth & Kornilov, *Robust error-minimization in the genetic code across physicochemical metrics and variant codes: a graph-theoretic analysis in GF(2)⁶*, *BioSystems*, in press).

This release supersedes v0.5.1. It closes out a second internal QA/QC read of the manuscript, supplement, and pipeline. **The four SUPPORTED scientific claims and the falsified KRAS–Fano result are unchanged in status and p-value; the tRNA-enrichment claim remains EXPLORATORY (as reclassified in v0.5.1). The changes in this release are additive, framing, plumbing, and reproducibility improvements.**

## What changed

### Null-model naming corrected + sensitivity companion added
- The coloring-optimality null implemented throughout the paper is a **quartet-pattern shuffle** (16 first-two-base codon quartets, each quartet's internal AA-pattern held atomic, permuted across quartet slots). Prior text called this "Freeland–Hurst block-preserving", which conflated it with the classical Haig & Hurst 1991 / Freeland & Hurst 1998 amino-acid-permutation null. Renamed consistently across manuscript, supplement, and code.
- Added the classical **Haig–Hurst amino-acid-permutation null** (`null_type="haig_hurst_aa"`) as a sensitivity companion (n = 10 000 per metric). Both nulls are reported side-by-side in Supplement §S3.1 with a preservation-property comparison and per-metric p-value table. All eight p-values (four metrics × two nulls) pass Bonferroni at α = 0.01; three of four code-independent metrics give a *more* extreme p under the classical null. The optimality claim is robust to the null-ensemble choice.

### tRNA analysis population made explicit
- New **24-row provenance table** in Supplement §S10.1: for each pairing, the variant + control organism, compartment, assembly accession, reassignment, Q₆ + H(3,4) topology status, exact 2×2 Fisher counts, per-organism source class (tRNAscan-SE / GtRNAdb / literature / annotation), and per-pair Fisher p-value.
- **Fisher denominator convention** now stated explicitly: by-amino-acid sum over the Std20 column of Table S8 (excludes SeC, undetermined, suppressor, and pseudogene rows).
- Wording corrected: *15 of 18 tRNAscan-SE-verified genomes populate the 24-pairing enrichment analysis; the remaining 3 (*Blastocrithidia nonstop*, *Mycoplasmoides genitalium*, *Mycoplasmoides pneumoniae*) are boundary-case mechanisms discussed only in §4.3*. Their reassignment routes (anticodon stem shortening; single-tRNA anticodon modification) act without changing gene counts, so a Fisher-count enrichment test would misinterpret their mechanism as absence-of-response.

### Data-plumbing fixes (silent-zero rendering bug)
- `manuscript_stats.condlogit.clade_exclusion` aggregate keys (min / median / max ΔAICc M1→M3) were absent from the emitted stats blob, causing the §3.5 template's `.at(..., default: 0)` fallback to silently render zeros ("ΔAICc ≥ 0, median 0, max 0"). Corrected. §3.5 now renders the true values (min = 48.8, median = 100.7, max = 116.7).
- "Static fallback table" if/else branches in Supplement §§S7 and S8 replaced with `#assert` guards so a missing data lift errors visibly at build time rather than silently emitting stale content.

### Framing and methods corrections
- §2.5 replaced with the discrete four-indicator structural-preservation index that Figure 5A actually uses (renamed *Structural-preservation index (visualization-only)*); the "synthetic-biology feasibility" phrasing dropped since no independently validated feasibility endpoint exists.
- §2.6 tRNA decision-gate restored to the originally specified worst-case-MIS criterion (not met at p = 0.123); median MIS p reported descriptively as a summary statistic, not as a re-specified criterion.
- §4.3 retitled *Exploratory tRNA accommodation across variant-code systems*; compensation-for-disconnection claims removed; genome-size claim explicitly framed as a hypothesis motivated by three case examples rather than a tested association.
- Per-table optimality reframed as "significantly low-cost *relative to its own quartet-shuffle null*"; the per-table quantile does not identify a boundary between proximity-driven and independently-optimized codes (every variant has d_H ≤ 6 while every null draw has d_H ≥ 30). Both registry-entry (11/12) and distinct-coloring (10/11) counts now reported.
- Restricted-candidate sensitivity: no *d* value designated the "primary" biological cut; three values (d ∈ {1, 2, 3}) reported as a bracket with d = 1 as the most defensible biological floor.
- IIA discussion rewritten: IIA remains an untested structural assumption even for the explanatory reading of the model.
- ρ-sweep clarified: 5 pre-specified inferential grid points (multiplicity-corrected); 6 additional plotted points are descriptive interpolation.
- Cross-family Bonferroni statement in §S18 corrected: the MIS-median p = 0.046 does **not** survive the α = 0.00625 threshold (it exceeds it by > 7×); consistent with the EXPLORATORY classification.
- Candidate-universe Table S3: U3 no-op count corrected to 61 (only sense codons carry an AA identity); source-stop axis removed from the introductory framing.
- Encoding-sweep claim narrowed: "all encoding-dependent results tested across 24 encodings" replaced with the exact list (coloring optimality + Q₆ landscape depletion). M3-H(3,4) renamed a *topology-verification* variant (its local physicochemical feature is still Q₆-based).
- Parametric-predictive rate reported with explicit numerator / denominator / adjacency / definition (5/66 event-steps under Q₆ / Δβ₀ > 0), distinguished from the lineage-collapsed 6/28 rate under H(3,4) / new-disconnection quoted in §3.4.
- Bibliography additions completed: Kosiol 2004; second Yurova Axelsson-Khrennikov 2024 paper.

### Reproducibility
- Version bumped 0.5.1 → 0.6.0 across `pyproject.toml`, `CITATION.cff`, `requirements.lock`, `Dockerfile`, `llms.txt`, manuscript and supplement reproducibility stamps.
- New `requirements.lock` (pinned direct + dev dependencies for CPython 3.11.14 reference build).
- New `Dockerfile` (Python 3.11.14-slim + Typst 0.14.2 + pinned pip install; regenerates every JSON artifact deterministically at seed 135325).
- Reproducibility claim downgraded from "bit-for-bit" to "within numerical / rendering tolerance" (§S24) since transitive dependencies (BLAS routing, libgfortran ABI, Typst layout across versions) are not fully pinned.
- CLI: `--null=haig_hurst_aa` option added to `codon-topo coloring`; `codon-topo all` now emits `output/coloring_optimality_hh_aa_null.json` and `output/tables/T-trna-24row-provenance.csv`.

### Rendering
- Fixed Typst source-escape defects that had rendered as literal glyphs in the previous PDF (`10^{-5}`, `1{,}280`, two `10{,}000`, `{≈}5×`, `{≈}100%`, `e.g._T. thermophila_` spacing, doubled "Table Table S1", Chan & Lowe DOI de-escaped for hayagriva).
- Cross-reference sweep: §S10.1 "main-text Table 8" → correct tRNA reference; Fig S7 label "disconnection × organism" → "variant-code × control"; Fig 1 caption gains "under the default encoding" for the ε = 4 Serine reconnection.
- Table S15 "statistically independent" → "approximately uncorrelated after adjustment".
- Walsh §S9: z-magnitude-ratio removed; conservative p = 1/(n+1) for zero exceedances.

### Code hygiene
- Ruff format + check clean across `src/` and `scripts/`.
- Pyright clean across `src/codon_topo` (0 errors, 0 warnings, 0 info).
- Full pipeline rerun at seed 135325 with the corrected code base.

## Scientific status (unchanged)

| Claim | Status | p-value |
|:--|:--|--:|
| Cross-metric coloring optimality (four metrics) | Supported | ≤ 0.006 (all metrics; quartet-pattern shuffle) |
| Per-table preservation across 27 NCBI tables | Supported | BH–FDR; 11/12 informative + 10/11 distinct colorings retain top-5% |
| ρ-robustness | Supported | 0.006 (Q₆) → 0.0003 (H(3,4)) |
| Topology avoidance under H(3,4) | Supported | RR 0.32; permutation p ≤ 10⁻⁴; hypergeometric p = 1.3 × 10⁻⁶ |
| tRNA enrichment | **Exploratory** | worst-case MIS 0.123 (gate not met); median 0.046 (descriptive); topology-breaking subset 0.43 (null) |
| KRAS–Fano clinical prediction | Falsified | 1.0 (MSK-IMPACT n = 1670) |

## Reproduction

```bash
git clone https://github.com/biostochastics/codontopo.git
cd codontopo
git checkout v0.6.0
pip install -r requirements.lock
pip install -e ".[dev]"
codon-topo all --output-dir=./output --seed=135325
python3.11 scripts/patch_pt_keys.py
python3.11 scripts/emit_trna_provenance_table.py
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

Or via container:

```bash
docker build -t codontopo:v0.6.0 .
docker run --rm -v "$PWD/output:/repo/output" codontopo:v0.6.0 \
    codon-topo all --output-dir=/repo/output --seed=135325
```

SHA-256 checksums of the reference-build artifacts are in `output/publisher_deliverable/checksums.sha256`.

## Downloads

- `manuscript.pdf` — main text
- `supplement.pdf` — supplementary material
- `publisher_letter.pdf` — cover letter to the BioSystems production editor
- `release_v0.6.0.zip` — reference-build artifacts (JSON + CSV + figure PNGs + PDFs)
- `biorxiv_bundle.zip` — bioRxiv-format submission bundle
- `submission_em.zip` — Editorial Manager LaTeX submission bundle
- `requirements.lock`, `Dockerfile` — pinned reference-build environment

## Correspondence & documentation

- `docs/publisher/2026-08-14-corrections-ledger.md` — full per-item ledger with commit provenance for v0.5.1 and v0.6.0 corrections
- `CHANGELOG.md` — release notes for v0.4.x / v0.5.0 / v0.5.1 / v0.6.0
- `README.md` — project overview
- `llms.txt` — LLM-oriented repository map
- `CITATION.cff` — machine-readable citation metadata

---

*Corresponding author:* Sergey Kornilov <sergey@biostochastics.com>
*License:* CC-BY-NC-4.0
