# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**CODON-TOPO** — Codon Geometry Validation & Prediction Engine. A computational pipeline validating the algebraic structure of genetic codes when encoded as 6-bit binary vectors in GF(2)^6. The project verifies and extends claims from Clayworth TN-2026-11 across all 27 NCBI translation tables (codes 1–6, 9–16, 21–33; codes 7, 8, 17–20 deprecated).

Authoritative scientific scope: see `src/codon_topo/reports/claim_hierarchy.py` (15 claims, status + p-values + justifications) and the manuscript at `output/manuscript.typ`. The original product-requirements document was internal-only and has been removed from the public release.

## Core Domain Concepts

- **GF(2)^6 encoding**: Each codon is a 6-bit binary vector (default: C=00, U=01, A=10, G=11, two bits per nucleotide position)
- **Filtration**: Two-fold degeneracy = bit-5 difference; four-fold degeneracy = shared 4-bit prefix with last 2 bits exhausting GF(2)^2
- **Persistent homology**: Connected components of amino acid codon graphs at increasing Hamming distance threshold epsilon
- **Serine invariant**: Serine is topologically disconnected at epsilon=1 in every known genetic code and every base-to-bit encoding — the only encoding-invariant disconnection. Under the default encoding, min inter-family Hamming distance = 4 (UCN vs AGY), uniquely extreme among 6-codon AAs (Leu and Arg = 1). Distance-4 holds for 8/24 encodings; 16/24 give distance 2. Converges with Serine's extreme Atchley Factor 3 score (Atchley et al. 2005)
- **Fano lines**: XOR triples in GF(2)^6 (e.g., GGU XOR GUU XOR CAC = 0)
- **Root-of-unity map**: phi: GF(2)^6 -> C^3, coordinate-wise bijection to fourth roots of unity (not holomorphic — domain is finite discrete)

## Claim Hierarchy

The single source of truth for what the paper claims is `src/codon_topo/reports/claim_hierarchy.py`. Run `codon-topo claims` to view it. See also `ARCHITECTURE.md` for the full module dependency graph.

Current status (15 claims):
- 4 SUPPORTED: hypercube coloring optimality (p=0.006), per-table preservation across the 27 NCBI tables (BH-FDR; see manuscript_stats.json for current count), rho robustness (p=0.003), topology avoidance depletion (permutation p=0.0001)
- 1 SUGGESTIVE: tRNA enrichment (worst-case MIS Stouffer p=0.045; 18 organisms, 17 tRNAscan-SE verified)
- 4 EXPLORATORY: bit-position bias, mechanism boundary conditions, Atchley F3/Serine convergence, variant-code disconnection catalogue
- 3 REJECTED: Serine min-distance-4 invariant, PSL(2,7), holomorphic embedding
- 1 FALSIFIED: KRAS-Fano clinical prediction (p=1.0)
- 2 TAUTOLOGICAL: two-fold bit-5 filtration, four-fold prefix filtration

## Package Structure (Implemented)

```
src/codon_topo/
  __init__.py            # Public API re-exports, DEFAULT_SEED=135325
  cli.py                 # Click-based CLI: codon-topo {filtration,disconnections,coloring,...}
  core/
    encoding.py          # Base->bit mappings, codon->vector conversion, all 24 encodings
    genetic_codes.py     # All 27 NCBI translation tables as dictionaries (codes 1–6, 9–16, 21–33)
    filtration.py        # Two-fold (bit-5) and four-fold (prefix) degeneracy checks
    homology.py          # Connected components, persistent homology, disconnection catalogue
    embedding.py         # Coordinate-wise root-of-unity map GF(2)^6 -> C^3
    fano.py              # XOR triple computation
  analysis/
    null_models.py       # Null Model A (random), B (block shuffle), C (24 encodings), C_extended
    reassignment_db.py   # Reassignment database, Hamming paths, bit-position bias (uniform + weighted)
    cosmic_query.py      # cBioPortal API client, KRAS Fano co-occurrence, WS4 gate
    depth_calibration.py # Evolutionary depth calibration, epsilon-age Spearman correlation
    synbio_feasibility.py # Synthetic biology feasibility scoring, reassignment landscape
    coloring_optimality.py # Hypercube coloring Monte Carlo (primary publishable result)
    trna_evidence.py     # tRNA enrichment test (18 organisms, 17 tRNAscan-SE verified + 1 from literature)
  data/
    grantham.json        # Grantham 1974 physicochemical distance matrix
    assembly_accessions.tsv  # NCBI genome assemblies for tRNAscan-SE verification
    trnascan_results/    # Raw tRNAscan-SE .out/.stats files
  visualization/
    data_export.py       # CSV export for R visualization (all workstreams)
  reports/
    catalogue.py         # Prediction catalogue with evidence grading (WS5 synthesis)
    claim_hierarchy.py   # Single source of truth for claim status (15 claims)
tests/
  test_encoding.py       # Encoding primitives + hypothesis property tests
  test_genetic_codes.py  # All 27 NCBI tables
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
  test_claim_hierarchy.py       # Claim hierarchy tests
  test_refinements.py           # Refinement-round tests
  test_cli.py                   # CLI subcommand tests
  test_ws_exports.py            # WS2-WS6 data exports
  test_integration_ws2_ws6.py   # Cross-workstream integration tests
```

## Technology Stack

- **Python 3.11+**, NumPy, SciPy for core computation
- **GUDHI or Ripser.py** for production persistent homology (hand-rolled code is verification-only)
- **ggplot2 + ggpubr** (R) for all figures — publication-quality, statistical annotations via ggpubr
- **click + rich** for CLI (`codon-topo` command)
- **requests** for cBioPortal API
- **SciPy.stats, statsmodels** for permutation tests, bootstrap CIs, multiple testing correction
- **pytest + hypothesis** for property-based testing of mathematical invariants
- **Sphinx + NumPy-style docstrings** for documentation

## Build & Test Commands

```bash
pip install -e ".[dev]"                              # Install in development mode
python3.11 -m pytest                                  # Run all tests (432 tests)
python3.11 -m pytest tests/test_encoding.py -v        # Run a single test file
python3.11 -m pytest tests/test_regression.py -v      # Run regression suite (105 tests)
python3.11 -m pytest --cov=codon_topo --cov-report=term-missing  # Coverage (≥96%)
codon-topo --help                                     # CLI usage
codon-topo claims                                     # View claim hierarchy
codon-topo all --output-dir=./output                  # Run all analyses
```

Note: Use `python3.11 -m pytest` (not bare `pytest`) because the system default Python is 3.14 but dev dependencies are installed under 3.11.

## Workstreams (WS1-WS6)

Dependency order: WS1 (foundation) -> WS4 (fast-fail KRAS test) -> WS2/WS3 (independent) -> WS6 -> WS5 (synthesis).

- **WS1**: Core replication, null models (A: random assignment 100k permutations, B: block-preserving shuffle, C: all 24 base encodings)
- **WS2**: Reassignment directionality — Hamming path analysis of ~40-60 known codon reassignment events
- **WS3**: Evolutionary depth calibration — epsilon vs. divergence age correlation
- **WS4**: KRAS/COSMIC Fano-line co-occurrence test (critical decision gate: null result kills clinical track only)
- **WS5**: Prediction catalogue synthesis
- **WS6**: Topology avoidance + conditional logit modeling + cross-study CodonSafe reanalysis (`topology-avoidance`, `topology-avoidance-k43`, `condlogit`, `phylo-sensitivity`, `codonsafe`)

## Regression Test Values (from Appendix 8)

These must be reproduced exactly:
- Two-fold filtration: bit-5 pattern across all 27 NCBI tables for the standard code's nine 2-fold AAs (zero exceptions among them); Table 32 (Balanophoraceae plastid, UAG→Trp) creates a new 2-fold Trp pair (UGG, UAG) at Hamming distance 2 (bits 3-4 differ rather than bit 5) — documented exception in `tests/test_filtration.py::TWOFOLD_FILTRATION_EXCEPTIONS`
- Four-fold filtration: 100% prefix uniformity, breaks only on stop->amino acid reassignments
- Serine: disconnected at epsilon=1 in all 27 tables, min inter-block Hamming = 4 under the default encoding, reconnects at epsilon=4
- Non-Serine disconnections: Thr (yeast mito, epsilon=2), Leu (chlorophycean mito, epsilon=2), Ala (Pachysolen nuclear, epsilon=3), three-component Ser (Candida nuclear, epsilon=3), Trp (Balanophoraceae plastid table 32, epsilon=2 — new in 27-table analysis)
- KRAS Fano line: GGU XOR GUU XOR CAC = 0 in GF(2)^6

## External Dependencies

- **tRNAscan-SE 2.0.12** at `~/.local/bin/tRNAscan-SE` with Infernal 1.1.4 at `~/.local/bin/`
- Run `bash scripts/run_trnascan.sh` to reproduce all tRNA scans from NCBI assemblies
- Assembly accessions in `data/assembly_accessions.tsv`; raw results in `data/trnascan_results/`

## Key Constraints

- All figures: ggplot2 + ggpubr (R), 300 DPI minimum, vector where possible, colorblind-friendly palette
- Null Model A (filtration/homology) requires 100,000 permutations; the coloring-optimality Monte Carlo (the headline result) uses 10,000 — these are distinct sample sizes for different tests, do not confuse them
- Success threshold: p < 0.01 for null model rejection
- 100% unit test coverage on core computation functions
