# CODON-TOPO

Codon Geometry Validation & Prediction Engine. Analyzes the algebraic structure of genetic codes encoded as 6-bit binary vectors in GF(2)^6.

**Version**: 0.2.0 | **Tests**: 344 passing | **Coverage**: ≥96% | **Target journal**: Journal of Theoretical Biology

## Project Status

Following multi-model adversarial review (10 LLMs, 12 passes), the claim hierarchy is:

| Status | Count | Key finding |
|--------|-------|-------------|
| **Supported** | 1 | Hypercube coloring optimality (p=0.006, Freeland-Hurst null) |
| **Suggestive** | 1 | tRNA duplication correlation (4/4 cases, p=0.0625) |
| **Exploratory** | 2 | Bit-position bias; variant-code disconnection catalogue |
| **Rejected** | 3 | Serine min-distance-4 invariant; PSL(2,7) symmetry; holomorphic embedding |
| **Falsified** | 1 | KRAS-Fano clinical prediction (p=1.0) |
| **Tautological** | 2 | Two-fold bit-5 filtration; four-fold prefix filtration |

See `codon-topo claims` or `src/codon_topo/reports/claim_hierarchy.py` for details.

### Scientific caveats

- The coloring optimality result replicates the *phenomenon* of Freeland-Hurst (1998) optimality under a stricter null, but not the original 10^-6 magnitude.
- The tRNA pattern is suggestive (n=4, binomial p=0.0625 trend). Su et al. 2011 mechanistically confirmed the yeast mito case.
- PSL(2,7) symmetry claim was dropped (pre-rejected by Antoneli & Forger 2011).
- KRAS-Fano enrichment test returned a clean negative (p=1.0 across all G12 variants).
- The coordinate-wise root-of-unity map is a bijection, not a holomorphic embedding (domain is finite discrete).
- Two-fold and four-fold filtration properties are tautological under any base→GF(2)^2 bijection.

## Quick Start

```bash
pip install -e ".[dev]"

python3.11 -m pytest                 # Run all 344 tests
python3.11 -m pytest --cov=codon_topo --cov-report=term-missing  # Coverage

codon-topo --help                    # CLI usage
codon-topo claims                    # View claim hierarchy
codon-topo coloring --n=10000        # Run optimality Monte Carlo
codon-topo all --output-dir=./output # Run everything, write JSON reports
```

> **Note**: Use `python3.11 -m pytest` if your system default Python differs from where dev dependencies are installed.

## CLI

```
codon-topo filtration [--table=1] [--all-tables] [--json]
codon-topo disconnections [--table=1] [--all-tables] [--extended] [--json]
codon-topo coloring [--null=freeland_hurst|class_size] [--n=10000] [--seed=135325] [--no-stops] [--json]
codon-topo bit-bias [--compartment=uniform|nuclear|mitochondrial] [--json]
codon-topo trna [--json]
codon-topo kras [--offline] [--json]
codon-topo claims [--json]
codon-topo all [--output-dir=./output] [--seed=135325] [--n=10000]
```

All subcommands support `--json` for machine-readable output. Interactive mode uses rich tables.

## Workstreams

| WS | Name | Status | Description |
|----|------|--------|-------------|
| **WS1** | Core Replication | Complete | Filtration, homology, embedding, Fano lines, null models A/B/C/C_extended |
| **WS2** | Reassignment Directionality | Complete | Database of codon reassignment events, Hamming paths, bit-position bias (uniform + Ts/Tv weighted) |
| **WS3** | Evolutionary Depth Calibration | Complete | Epsilon-age Spearman correlation across 6 calibration points (non-result: rho=0.0) |
| **WS4** | KRAS/COSMIC Fast-Fail Gate | Complete | cBioPortal API client, Fano enrichment test — **clean negative** (p=1.0) |
| **WS5** | Prediction Catalogue | Complete | 15 predictions with evidence grading and claim hierarchy |
| **WS6** | Synthetic Biology Feasibility | Complete | Feasibility scoring for 1280 single-codon reassignment variants |

## Package Structure

```
src/codon_topo/
  cli.py                 # Click-based CLI (codon-topo command)
  core/
    encoding.py          # GF(2)^6 binary encoding, Hamming distance, all 24 encodings
    genetic_codes.py     # All 25 NCBI translation tables
    filtration.py        # Two-fold (bit-5) and four-fold (prefix) degeneracy checks
    homology.py          # Persistent homology via union-find, disconnection catalogue
    embedding.py         # Coordinate-wise root-of-unity map GF(2)^6 -> C^3
    fano.py              # XOR triple computation
  analysis/
    null_models.py       # Null models A (random), B (block shuffle), C (24 encodings), C_extended
    reassignment_db.py   # Reassignment database, Hamming paths, bit-position bias (uniform + weighted)
    cosmic_query.py      # cBioPortal API client, KRAS Fano co-occurrence, WS4 gate
    depth_calibration.py # Evolutionary depth calibration, epsilon-age correlation
    synbio_feasibility.py # Feasibility scoring for alternative genetic codes
    coloring_optimality.py # Hypercube coloring Monte Carlo (primary publishable result)
    trna_evidence.py     # tRNA duplication correlation test (4 organisms + controls)
  data/
    grantham.json        # Grantham 1974 physicochemical distance matrix
  visualization/
    data_export.py       # CSV export for R visualization scripts
  reports/
    catalogue.py         # Prediction catalogue with evidence grading
    claim_hierarchy.py   # Single source of truth for claim status (10 claims)
```

## Usage Examples

```python
from codon_topo import (
    codon_to_vector, hamming_distance, STANDARD, get_code,
    analyze_filtration, disconnection_catalogue, embed_codon,
    is_fano_line, fano_partner, monte_carlo_null,
    CLAIM_HIERARCHY, supported_claims,
)

# Encode a codon as a 6-bit vector
codon_to_vector('GGU')  # (1, 1, 1, 1, 0, 1)

# Check the KRAS Fano line (XOR = 0)
is_fano_line('GGU', 'GUU', 'CAC')  # True
fano_partner('GGU', 'GUU')  # 'CAC'

# Analyze filtration for any genetic code
analyze_filtration(STANDARD)
# {'twofold_pass': 9, 'twofold_fail': 0, 'fourfold_pass': 5, ...}

# Find all disconnected amino acids
disconnection_catalogue(STANDARD)
# [{'aa': 'Ser', 'n_components': 2, 'reconnect_eps': 4, ...}]

# Run the coloring optimality Monte Carlo
result = monte_carlo_null(n_samples=10000, seed=135325)
# {'quantile_of_observed': 0.6, 'p_value_conservative': 0.006, ...}

# Query the claim hierarchy
for claim in supported_claims():
    print(claim.id, claim.evidence_p_value)
```

## Null Models

| Model | What it tests | Method |
|-------|---------------|--------|
| **A** | Is the codon→AA assignment special? | Fix degeneracy structure, randomly reassign codons (100k permutations) |
| **B** | Is the block structure special? | Preserve 4-codon blocks, shuffle block→AA assignments |
| **C** | Is the binary encoding special? | Test all 24 possible base→bit-pair mappings |
| **C_extended** | Per-encoding min distances | Emit per-encoding min inter-block Hamming for each disconnected AA |
| **Freeland-Hurst** | Is the coloring optimal? | Block-preserving shuffle with Grantham distance scoring |
| **Class-size** | Weaker coloring null | Class-size-preserving shuffle (no block contiguity) |

## Technology

- **Python 3.11+**, NumPy, SciPy, requests, click, rich
- **pytest + hypothesis** for property-based testing
- **ggplot2 + ggpubr** (R) for publication figures
- **SciPy.stats** for Spearman correlation, Fisher's exact test, chi-squared, bootstrap CIs
