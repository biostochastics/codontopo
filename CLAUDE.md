# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**CODON-TOPO** — Codon Geometry Validation & Prediction Engine. A computational pipeline validating the algebraic structure of genetic codes when encoded as 6-bit binary vectors in GF(2)^6. The project verifies and extends claims from Clayworth TN-2026-11 across all 23 NCBI translation tables.

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
| `all_codes.py` | All 23 NCBI tables, disconnection catalogue, filtration breakage | Integration tests, `analysis/reassignment_db.py` |
| `summary_table.py` | Disconnection-depth catalogue with evolutionary dating context | Reports |

Shared utilities duplicated across scripts (must be consolidated):
- `BASE_TO_BITS = {'C': (0,0), 'U': (0,1), 'A': (1,0), 'G': (1,1)}`
- `codon_to_bits(codon)` — codon string to 6-bit tuple
- `hamming_distance(a, b)` — bitwise Hamming distance
- `connected_components(bits_list, epsilon)` — union-find at threshold
- `STANDARD` / `GENETIC_CODE` — standard code dict (duplicated 4 times)
- `all_codes.py` has all 23 NCBI tables via `make_code(changes)` pattern

## Target Package Structure

```
codon_topo/
  core/
    encoding.py        # Base->bit mappings, codon->vector conversion, all 24 encodings
    genetic_codes.py   # All NCBI translation tables as dictionaries
    filtration.py      # Projection tower, fibre computation, filtration checks
    homology.py        # Connected components, persistent homology (Betti numbers)
    embedding.py       # Holomorphic embedding phi: GF(2)^6 -> C^3
    fano.py            # Fano-line computation, XOR triples
  analysis/
    null_models.py     # Random assignment, block shuffle, encoding permutation
    reassignment_db.py # Reassignment database construction and querying
    cosmic_query.py    # COSMIC/cBioPortal API integration
    synbio_meta.py     # Synthetic biology literature meta-analysis
  visualization/
    hypercube.py       # 6D->2D/3D projections with amino acid coloring
    barcodes.py        # Persistent homology barcode plots
    embedding_viz.py   # C^3 scatter plots with degeneracy coloring
    reassignment_flow.py # Hamming-path flow diagrams
  reports/
    templates/         # Report generation templates
  tests/
    test_core.py       # Unit tests for all core computations
    test_null.py       # Null model validation tests
    test_regression.py # Regression tests against known results from Appendix 8
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
pip install -e ".[dev]"          # Install in development mode
pytest                           # Run all tests
pytest tests/test_core.py        # Run core computation tests only
pytest tests/test_core.py::test_name -v  # Run a single test
pytest --cov=codon_topo          # Run with coverage
```

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
- Two-fold filtration: 100% bit-5 pattern across all 23 NCBI tables, zero exceptions
- Four-fold filtration: 100% prefix uniformity, breaks only on stop->amino acid reassignments
- Serine: disconnected at epsilon=1 in all 23 tables, min inter-block Hamming = 4, reconnects at epsilon=4
- Non-Serine disconnections: Thr (yeast mito, epsilon=2), Leu (chlorophycean mito, epsilon=2), Ala (Pachysolen nuclear, epsilon=3), three-component Ser (Candida nuclear, epsilon=3)
- KRAS Fano line: GGU XOR GUU XOR CAC = 0 in GF(2)^6

## Key Constraints

- All figures: ggplot2 + ggpubr (R), 300 DPI minimum, vector where possible, colorblind-friendly palette
- Null Model A requires 100,000 permutations
- Success threshold: p < 0.01 for null model rejection
- 100% unit test coverage on core computation functions
