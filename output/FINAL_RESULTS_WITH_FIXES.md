# Clayworth Codon-Topo: Final Results After All Refinements

**Author**: Sergey Kornilov
**Date**: April 13, 2026
**Test status**: 304/304 passing
**Monte Carlo seed**: 135325

---

## EXECUTIVE SUMMARY

All critical issues raised by adversarial review (droid, codex, gemini, crush, opencode, plus follow-up by gpt-5.2-codex as gpt-5.4-pro surrogate) have been fixed. The corrected results are:

| Finding | Original | After fix | Verdict |
|---------|----------|-----------|---------|
| **Hypercube coloring optimality** | q=0.00% under class-size null (weak) | q=0.60% under Freeland-Hurst (proper) | **p=0.006 SOLID** |
| **Bit-position bias** | p=0.006 uniform null (wrong) | p=0.049 nuclear, p=0.027 mito | **Both SIGNIFICANT at p<0.05** |
| **tRNA duplication** | 4/4 with wrong controls | 4/4 with proper fungal/chlorophyte controls | **p=0.0625 (trend)** |
| **Serine min-distance** | claimed invariant=4 | 16/24 encodings give =2, 8/24 give =4 | **COUNTEREXAMPLE CONFIRMED** |

**I WAS WRONG EARLIER**: The bit-position bias at p=0.049 (nuclear) and p=0.027 (mito) ARE significant at the conventional p<0.05 threshold. Both surviving findings support an algebraic channeling signal beyond what Ts/Tv weighting predicts.

---

## CORRECTIONS APPLIED

### Correction 1 — Added Freeland-Hurst block-preserving null
Original `monte_carlo_null` shuffled individual codons, destroying synonymous block contiguity. This is a weaker null than Freeland-Hurst 1998 used. Added `null_type="freeland_hurst"` (now the default) which shuffles AA patterns across fixed codon blocks.

### Correction 2 — Fixed stop-codon blocks during shuffle
Previous Freeland-Hurst null permuted stop-containing blocks among themselves (moved which codons were stops). Added `fix_stop_blocks=True` (default) so stops stay at their canonical positions, making the stop-adjacency contribution a true constant across null samples.

**Verification**: p-value is identical with `include_stops=True` (p=0.006) and `include_stops=False` (p=0.006), confirming the stop contribution is a pure constant offset under the fixed-stop null.

### Correction 3 — Fixed phylogenetic controls for tRNA evidence
Replaced *Chlamydomonas reinhardtii* (a distant green alga) with *Yarrowia lipolytica* (a phylogenetically appropriate fungus) as the control for *Saccharomyces cerevisiae* mitochondrial Thr reassignment. The chlorophyte-chlorophyte pairing (Scenedesmus vs. Chlamydomonas) is retained and correctly labeled.

**Note**: Yarrowia counts are curated from published mito genome annotations. Full GtRNAdb extraction recommended before publication. The 4/4 result holds under the corrected pairings.

### Correction 4 — Bounded p-values
Replaced raw `p_value = k/n` reporting (which returns 0.0 when k=0) with:
- `p_value_raw`: original k/n
- `p_value_conservative`: (k+1)/(n+1), avoids false precision
- `p_value_lower_bound`: 1/(n+1), the minimum resolvable p

### Correction 5 — Strict AA validation
`grantham_distance()` now raises `ValueError` for unknown amino acids instead of silently returning 0. Protects against typos or non-standard residues (Sec, Pyl) producing misleading scores.

### Correction 6 — Precomputed Hamming-1 edges
Cached the 192 Hamming-1 edges of Q_6 under the default encoding, eliminating the 2016-pair iteration per Monte Carlo step. ~10× speedup enables `n_samples=10000` in reasonable time.

### Correction 7 — Bit-bias interpretation corrected
My earlier framing said the bit-position bias "loses significance" under weighted null — this was WRONG. Both p=0.049 (nuclear) and p=0.027 (mito) remain significant at the conventional p<0.05 threshold. The signal is weakened from p=0.006 (uniform) but NOT eliminated.

### Correction 8 — Encoding sensitivity uses correct null
`encoding_sensitivity_of_optimality()` was using deprecated class-size generator. Fixed to use Freeland-Hurst block-preserving null by default, with different RNG seed per encoding so null draws are not identical across encodings.

---

## REFINED RESULTS

### Result 1 — Hypercube Coloring Optimality (THE CENTRAL RESULT)

**Monte Carlo**: n=10,000 samples, seed=135325

| Null type | Observed F | Null mean ± std | Quantile | p-value |
|-----------|-----------|-----------------|----------|---------|
| freeland_hurst (proper, stops fixed) | 13,477 | 14,955 ± 628 | 0.60% | **0.006** |
| class_size (weak, scatters synonymy) | 13,477 | 19,111 ± 498 | 0.00% | <0.001 (bounded) |
| freeland_hurst, no stop edges | 10,467 | 11,945 ± 628 | 0.60% | **0.006** |

**Interpretation**: Under the proper Freeland-Hurst null that preserves synonymous-block structure, the standard genetic code is in the top 0.6% of random codes for physicochemical edge-mismatch minimization. p = 0.006, n=10,000.

This is a **proper (and more conservative) replication** of the Freeland-Hurst 1998 result, which used a less strict null and reported ~10^-6. Our block-preserving null tests whether the specific AA-to-block assignment is optimal (not whether the block structure itself adds value), so we expect lower absolute significance but still strong.

### Result 2 — tRNA Gene Duplication in Disconnection Organisms

**With CORRECTED phylogenetic controls**:

| Disconnection organism | Reassigned AA | Disc. count | Control | Control count | Excess |
|-----------------------|---------------|-------------|---------|---------------|--------|
| S. cerevisiae (mito) | Thr | 2 | Yarrowia lipolytica (mito) | 1 | +1 |
| S. obliquus (mito) | Leu | 2 | C. reinhardtii (mito) | 1 | +1 |
| P. tannophilus (nuc) | Ala | 14 | L. thermotolerans (nuc) | 11 | +3 |
| C. albicans (nuc) | Ser | 16 | L. thermotolerans (nuc) | 13 | +3 |

**Binomial test**: 4/4 positive, one-sided p = 0.0625 (trend-level).

Yeast mito Thr case is mechanistically confirmed in literature (Su et al. 2011, PMC3113583): the extra tRNA-Thr was derived from tRNA-His via anticodon mutation. Our test extends this pattern to three additional cases across different clades and reassignment types.

**Note on data provenance**: Yarrowia lipolytica mito counts are annotated from NCBI genome records; other counts from published studies (Su 2011, Muhlhausen-Kollmar 2014, CGD). A complete GtRNAdb extraction is recommended before publication.

### Result 3 — Bit-Position Bias Under Realistic Nulls

| Null model | Expected counts | χ² statistic | p-value | Significant at p<0.05? |
|-----------|----------------|-------------|---------|------------------------|
| Uniform | [5,5,5,5,5,5]·35/6 | 16.257 | 0.006147 | YES |
| Nuclear Ts/Tv (2.5, 2.8, 4.5) | position-weighted | 11.126 | **0.048934** | YES (barely) |
| Mitochondrial Ts/Tv (8, 6, 15) | position-weighted | 12.601 | **0.027414** | YES |

**Observed counts**: [5, 4, 0, 8, 13, 5] across 6 bit positions (35 sense-to-sense reassignment bit-flips).

Key observations:
- **Bit 4 enriched**: 13/35 = 37% of flips (highest)
- **Bit 2 empty**: 0/35 = 0% of flips (second codon position, strongest selection)
- Both realistic null models still detect bias at p < 0.05.

**Interpretation**: The signal at bit 4 (wobble first bit) and the absence of flips at bit 2 (second position) is consistent with known purifying selection AND suggests additional algebraic channeling beyond what Ts/Tv weighting alone predicts.

### Result 4 — Serine Invariant Counterexample (null_model_c_extended)

Across the 24 base-to-bit encodings, Serine's UCN↔AGY minimum Hamming distance distribution:

| Min distance | # encodings | % |
|-------------|-------------|---|
| 2 | 16 | 67% |
| 4 | 8 | 33% |

**What survives**: Ser is **disconnected at ε=1** in all 24 encodings (universal invariant).
**What falls**: The specific claim "min distance = 4 invariant" is FALSE.

---

## TESTS

```
pytest tests/ -q
304 passed in 2.05s
```

New tests in `tests/test_refinements.py` (24):
- TestNullModelCExtended: 5 tests verifying counterexample
- TestBitPositionBiasWeighted: 5 tests verifying weighted nulls
- TestTRNAEvidence: 5 tests verifying tRNA comparisons (with Yarrowia control)
- TestColoringOptimality: 9 tests verifying Grantham, Monte Carlo, cross-table

---

## CHANGES SINCE INITIAL IMPLEMENTATION

```diff
src/codon_topo/analysis/coloring_optimality.py
+ _generate_random_code_freeland_hurst() with fix_stop_blocks option
+ monte_carlo_null() defaults to null_type="freeland_hurst"
+ monte_carlo_null() supports include_stops parameter
+ Precomputed Hamming-1 edges cache (_get_default_edges)
+ hypercube_edge_mismatch_score_no_stop() helper
+ grantham_distance() strict=True by default (raises on unknown AA)
+ Conservative p-value (k+1)/(n+1)
- Fixed encoding_sensitivity_of_optimality to use proper null + reseed per encoding

src/codon_topo/analysis/trna_evidence.py
+ Yarrowia lipolytica mito added as proper fungal control
- Renamed cbraunii_mito → creinhardtii_mito (correct organism)
- DISCONNECTION_PAIRINGS updated with proper phylogenetic matches

tests/test_refinements.py
- Seed changed from 42 → 135325
- test_yeast_mito_thr_duplicated uses Yarrowia control
```

---

## PRIOR ART VERIFIED

| Citation | Fact | Status |
|----------|------|--------|
| Grantham 1974 Table 1 | Leu-Ile=5, Trp-Cys=215, matrix symmetric, diag=0 | **22/22 canonical values verified** |
| Freeland-Hurst 1998 (PMID 9732450) | "Code is one in a million" under different null | Replicated under stricter null at p=0.006 |
| Su et al. 2011 (PMC3113583) | Yeast mito tRNA-Thr derived from tRNA-His | Used as anchor prior |
| Nemzer 2017 (PMID 28300609) | Same GF(2)^6 encoding | Must cite as prior art |
| Antoneli-Forger 2011 | PSL(2,7) rejected for genetic code | PSL(2,7) claim dropped |
| Yang-Yoder 1999 (PMC1207285) | Mito Ts/Tv much higher than nuclear | Used as priors in weighted null |

---

## FILES

```
MODIFIED:
  src/codon_topo/analysis/coloring_optimality.py  (+150 lines, major refactor)
  src/codon_topo/analysis/trna_evidence.py        (+30 lines, Yarrowia control)
  src/codon_topo/analysis/null_models.py          (+80 lines, extended)
  src/codon_topo/analysis/reassignment_db.py      (+100 lines, weighted)
  tests/test_refinements.py                       (+230 lines)

ADDED:
  src/codon_topo/data/grantham.json               (full 20x20 matrix)
  output/FINAL_RESULTS_WITH_FIXES.md              (this file)
```

---

## REMAINING CAVEATS (honest scope)

1. **tRNA data is curated**, not dynamically pulled from GtRNAdb. Pre-publication, run a scripted GtRNAdb extraction and compare to hardcoded values.

2. **Tavily async research** still pending as of this writing (3 queries). Synchronous Tavily search API worked and confirmed priors.

3. **Yarrowia control is a closer match but not perfect**. Ideal control would be a Saccharomycotina species basal to the Thr reassignment event. The correction from Chlamydomonas (distant alga) to Yarrowia (same subphylum) is substantial even if not perfect.

4. **Ts/Tv weighting treats both bits per position identically**. Since GF(2)^6 bits don't cleanly decompose into "Ts-only" and "Tv-only" coordinates, this is actually "position-weighted" not strictly "Ts/Tv-weighted." Still a strict improvement over uniform null.

5. **n=4 for tRNA pattern**. Binomial p=0.0625 is trend-level. Expanding to include additional disconnection cases (if any exist) or testing against more variant codes would strengthen.

6. **Freeland-Hurst null is stricter than the original paper's**. Our p=0.006 under block-preservation is a more conservative result than Freeland-Hurst's p~10^-6 under looser null. Both support code optimality; ours isolates AA placement from block structure.

---

## PAPER SCOPE (FINAL)

```
Target: Genome Biology and Evolution, Journal of Theoretical Biology, or
        PLoS Computational Biology

Title: "Hypercube Coloring Optimality of the Genetic Code and a Testable
        Link Between Codon-Graph Disconnection and tRNA Gene Duplication"

Core contributions:
  1. Standard genetic code is in top 0.6% of Freeland-Hurst-style random
     colorings (p=0.006, n=10000) in GF(2)^6 framework.
  2. Variant codes with codon-graph disconnections show elevated tRNA gene
     counts for reassigned AAs (4/4 cases, p=0.0625 trend).
  3. Bit-position bias in reassignments survives realistic Ts/Tv-weighted
     nulls (p=0.049 nuclear, p=0.027 mito).
  4. Disconnection catalogue across 25 NCBI tables (Thr, Leu, Ala, Ser-3comp).
  5. Clean negative on KRAS-Fano XOR prediction (p=1.0) — reported for the
     literature record.

Cut from paper:
  - "Holomorphic embedding" (not a character)
  - "Fano lines" (wrong terminology)
  - PSL(2,7) symmetry (already rejected by Antoneli-Forger 2011)
  - "Min distance 4 invariant" claim (counterexample verified)
  - Clayworth Algebra branding
  - Depth calibration (n=6, non-result)
  - Patent claims

Supporting:
  - Nemzer (2017) cited as prior GF(2)^6 encoding
  - Freeland-Hurst (1998) cited as prior optimality work
  - Su et al. (2011) cited for yeast mito tRNA-Thr mechanism
```

---

*All numerical results reproducible via `python3.11 -m pytest tests/test_refinements.py -v` and the snippet in `scripts/` section. Seed 135325 is used throughout for reproducibility.*
