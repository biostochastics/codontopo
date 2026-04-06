# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**CODON-TOPO** — Codon Geometry Validation & Prediction Engine. A computational pipeline validating the algebraic structure of genetic codes when encoded as 6-bit binary vectors in GF(2)^6. The project verifies and extends claims from Clayworth TN-2026-11 across all 25 NCBI translation tables.

The PRD is in `CODON_TOPO_PRD_v1.docx` at the repository root.

## Core Domain Concepts

- **GF(2)^6 encoding**: Each codon is a 6-bit binary vector (default: C=00, U=01, A=10, G=11, two bits per nucleotide position)
- **Filtration**: Two-fold degeneracy = bit-5 difference; four-fold degeneracy = shared 4-bit prefix with last 2 bits exhausting GF(2)^2
- **Persistent homology**: Connected components of amino acid codon graphs at increasing Hamming distance threshold epsilon
- **Serine invariant**: Serine is topologically disconnected (epsilon=4) in every known genetic code — a universal invariant
- **Fano lines**: XOR triples in GF(2)^6 (e.g., GGU XOR GUU XOR CAC = 0)
- **Holomorphic embedding**: phi: GF(2)^6 -> C^3, maps codons to complex 3-space preserving degeneracy structure

## Preliminary Code

`codon_topo_preliminary/` contains 6 verified scripts (stdlib-only Python, no external deps) that form the foundation for the `codon_topo` package. All produce correct output — these are the regression reference.

| Script | What it does | Refactors into |
|---|---|---|
| `verify_filtration.py` | Claims 1-3: bit-5 two-fold, prefix four-fold, Serine persistent homology | `core/filtration.py`, `core/homology.py` |
| `verify_embedding.py` | Holomorphic embedding phi: GF(2)^6 -> C^3, base->fourth-root-of-unity | `core/embedding.py` |
| `verify_kras.py` | KRAS G12V Fano line (GGU XOR GUU XOR CAC = 0), all single-bit mutations | `core/fano.py`, `analysis/cosmic_query.py` |
| `verify_mitochondrial.py` | 6-code cross-validation, Serine persistence, Thr disconnection discovery | `core/genetic_codes.py` + analysis pipeline |
| `all_codes.py` | All 24 NCBI tables, disconnection catalogue, filtration breakage | Integration tests, `analysis/reassignment_db.py` |
| `summary_table.py` | Disconnection-depth catalogue with evolutionary dating context | Reports |

Shared utilities duplicated across scripts (must be consolidated):
- `BASE_TO_BITS = {'C': (0,0), 'U': (0,1), 'A': (1,0), 'G': (1,1)}`
- `codon_to_bits(codon)` — codon string to 6-bit tuple
- `hamming_distance(a, b)` — bitwise Hamming distance
- `connected_components(bits_list, epsilon)` — union-find at threshold
- `STANDARD` / `GENETIC_CODE` — standard code dict (duplicated 4 times)
- `all_codes.py` has all 24 NCBI tables via `make_code(changes)` pattern

## Package Structure (Implemented)

```
src/codon_topo/
  __init__.py            # Public API re-exports
  core/
    encoding.py          # Base->bit mappings, codon->vector conversion, all 24 encodings
    genetic_codes.py     # All 25 NCBI translation tables as dictionaries
    filtration.py        # Two-fold (bit-5) and four-fold (prefix) degeneracy checks
    homology.py          # Connected components, persistent homology, disconnection catalogue
    embedding.py         # Holomorphic embedding phi: GF(2)^6 -> C^3
    fano.py              # Fano-line computation, XOR triples
  analysis/
    null_models.py       # Null Model A (random), B (block shuffle), C (24 encodings)
    reassignment_db.py   # Reassignment database, Hamming path analysis, directionality stats
    cosmic_query.py      # cBioPortal API client, KRAS Fano-line co-occurrence, WS4 gate
    depth_calibration.py # Evolutionary depth calibration, epsilon-age Spearman correlation
    synbio_feasibility.py # Synthetic biology feasibility scoring, reassignment landscape
  visualization/
    data_export.py       # CSV export for R visualization (all workstreams)
  reports/
    catalogue.py         # Prediction catalogue with evidence grading (WS5 synthesis)
tests/
  test_encoding.py       # Encoding primitives + hypothesis property tests
  test_genetic_codes.py  # All 25 NCBI tables
  test_filtration.py     # Two-fold/four-fold checks across all tables
  test_homology.py       # Serine invariant, novel disconnections
  test_embedding.py      # C^3 embedding properties
  test_fano.py           # KRAS Fano line, XOR triples
  test_null_models.py    # Null model smoke tests
  test_regression.py     # Full PRD Appendix 8 regression (105 tests)
  test_reassignment_db.py       # WS2 reassignment database
  test_reassignment_stats.py    # WS2 directionality statistics
  test_depth_calibration.py     # WS3 evolutionary depth calibration
  test_cosmic_query.py          # WS4 cBioPortal client
  test_kras_enrichment.py       # WS4 Fano enrichment test
  test_synbio_feasibility.py    # WS6 feasibility scoring
  test_catalogue.py             # WS5 prediction catalogue
  test_ws_exports.py            # WS2-WS6 data exports
  test_integration_ws2_ws6.py   # Cross-workstream integration tests
```

## Technology Stack

- **Python 3.11+**, NumPy, SciPy for core computation
- **GUDHI or Ripser.py** for production persistent homology (hand-rolled code is verification-only)
- **ggplot2 + ggpubr** (R) for all figures — publication-quality, statistical annotations via ggpubr
- **requests + pandas** for COSMIC REST API, cBioPortal API
- **SciPy.stats, statsmodels** for permutation tests, bootstrap CIs, multiple testing correction
- **pytest + hypothesis** for property-based testing of mathematical invariants
- **Sphinx + NumPy-style docstrings** for documentation

## Build & Test Commands

```bash
pip install -e ".[dev]"                              # Install in development mode
python3.11 -m pytest                                  # Run all tests (280 tests)
python3.11 -m pytest tests/test_encoding.py -v        # Run a single test file
python3.11 -m pytest tests/test_regression.py -v      # Run regression suite (105 tests)
python3.11 -m pytest --cov=codon_topo --cov-report=term-missing  # Coverage (96%)
```

Note: Use `python3.11 -m pytest` (not bare `pytest`) because the system default Python is 3.14 but dev dependencies are installed under 3.11.

## Workstreams (WS1-WS6)

Dependency order: WS1 (foundation) -> WS4 (fast-fail KRAS test) -> WS2/WS3 (independent) -> WS6 -> WS5 (synthesis).

- **WS1**: Core replication, null models (A: random assignment 100k permutations, B: block-preserving shuffle, C: all 24 base encodings)
- **WS2**: Reassignment directionality — Hamming path analysis of ~40-60 known codon reassignment events
- **WS3**: Evolutionary depth calibration — epsilon vs. divergence age correlation
- **WS4**: KRAS/COSMIC Fano-line co-occurrence test (critical decision gate: null result kills clinical track only)
- **WS5**: Prediction catalogue synthesis
- **WS6**: Synthetic biology feasibility + commercial assessment

## Regression Test Values (from Appendix 8)

These must be reproduced exactly:
- Two-fold filtration: 100% bit-5 pattern across all 24 NCBI tables, zero exceptions
- Four-fold filtration: 100% prefix uniformity, breaks only on stop->amino acid reassignments
- Serine: disconnected at epsilon=1 in all 24 tables, min inter-block Hamming = 4, reconnects at epsilon=4
- Non-Serine disconnections: Thr (yeast mito, epsilon=2), Leu (chlorophycean mito, epsilon=2), Ala (Pachysolen nuclear, epsilon=3), three-component Ser (Candida nuclear, epsilon=3)
- KRAS Fano line: GGU XOR GUU XOR CAC = 0 in GF(2)^6

## Key Constraints

- All figures: ggplot2 + ggpubr (R), 300 DPI minimum, vector where possible, colorblind-friendly palette
- Null Model A requires 100,000 permutations
- Success threshold: p < 0.01 for null model rejection
- 100% unit test coverage on core computation functions
