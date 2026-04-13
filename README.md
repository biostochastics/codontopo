# CODON-TOPO

Codon Geometry Validation & Prediction Engine. Analyzes the algebraic structure of genetic codes encoded as 6-bit binary vectors in GF(2)^6.

**Version**: 0.3.0 | **Tests**: 367 passing | **Coverage**: >=96% | **Target journal**: Journal of Theoretical Biology

## Project Status

15 claims evaluated. Single source of truth: `src/codon_topo/reports/claim_hierarchy.py`.

| Status | Count | Key findings |
|--------|-------|--------------|
| **Supported** | 4 | Cross-metric coloring optimality (p<0.004 across Grantham, Miyata, polar requirement, Kyte-Doolittle); per-table preservation (24/25); rho robustness (p=0.003); topology avoidance depletion (perm p=0.0001, robust to phylogenetic clade exclusion) |
| **Suggestive** | 1 | tRNA enrichment (worst-case MIS Stouffer p=0.042 across 17 pairings, 10 tRNAscan-SE verified genomes) |
| **Exploratory** | 4 | Bit-position bias; mechanism boundary conditions; Atchley F3/Serine convergence; variant-code disconnection catalogue |
| **Rejected** | 3 | Serine min-distance-4 invariant; PSL(2,7) symmetry; holomorphic embedding |
| **Falsified** | 1 | KRAS-Fano clinical prediction (p=1.0) |
| **Tautological** | 2 | Two-fold bit-5 filtration; four-fold prefix filtration |

Run `codon-topo claims` for the full hierarchy with p-values and justifications.

### Scientific summary

The standard genetic code is significantly error-minimizing under four independent physicochemical distance metrics (Grantham 1974, Miyata 1979, Woese polar requirement, Kyte-Doolittle hydropathy; all p < 0.004), establishing this as a general structural property rather than a metric-specific artifact. Freeland and Hurst (1998) demonstrated optimality under polar requirement alone; the GF(2)^6 graph decomposition enables systematic interpolation across the full codon mutation graph K_4^3 via a transversion weight parameter rho, with optimality strengthening at rho=1 (all p < 0.003). This structure is preserved in 24 of 25 NCBI variant codes. Natural codon reassignments are 3.3-fold depleted for topology-breaking changes (22% vs 73%; permutation p=0.0001), a finding robust to phylogenetic clade exclusion (all p < 10^-5). Organisms with variant codes show elevated tRNA gene copy numbers for the reassigned amino acid (worst-case maximal independent set Stouffer p=0.042, 10 tRNAscan-SE verified genomes). Several algebraic conjectures (PSL(2,7), holomorphic embedding, KRAS-Fano) were cleanly falsified.

## Quick Start

```bash
pip install -e ".[dev]"

python3.11 -m pytest                 # Run all 367 tests
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
codon-topo metric-sensitivity [--n=1000] [--seed=135325] [--json]
codon-topo rho-sweep [--n=1000] [--seed=135325] [--json]
codon-topo per-table [--n=1000] [--seed=135325] [--json]
codon-topo topology-avoidance [--json]
codon-topo phylo-sensitivity [--json]
codon-topo trna [--json]
codon-topo mis-analysis [--json]
codon-topo bit-bias [--compartment=uniform|nuclear|mitochondrial] [--json]
codon-topo kras [--offline] [--json]
codon-topo claims [--json]
codon-topo decompose [--json]
codon-topo all [--output-dir=./output] [--seed=135325] [--n=10000]
```

All subcommands support `--json` for machine-readable output. Interactive mode uses rich tables.

## Workstreams

| WS | Name | Status | Description |
|----|------|--------|-------------|
| **WS1** | Core Replication | Complete | Filtration, homology, embedding, Fano lines, null models A/B/C/C_extended |
| **WS2** | Reassignment Directionality | Complete | Database of codon reassignment events, Hamming paths, bit-position bias (uniform + Ts/Tv weighted) |
| **WS3** | Evolutionary Depth Calibration | Complete | Epsilon-age Spearman correlation across 6 calibration points (non-result: rho=0.0) |
| **WS4** | KRAS/COSMIC Fast-Fail Gate | Complete | cBioPortal API client, Fano enrichment test — clean negative (p=1.0) |
| **WS5** | Prediction Catalogue | Complete | 15 predictions with evidence grading and claim hierarchy |
| **WS6** | Synthetic Biology Feasibility | Complete | Feasibility scoring, topology avoidance test (perm p=0.0001) |

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
    reassignment_db.py   # Reassignment database, Hamming paths, bit-position bias
    cosmic_query.py      # cBioPortal API client, KRAS Fano co-occurrence
    depth_calibration.py # Evolutionary depth calibration, epsilon-age correlation
    synbio_feasibility.py # Feasibility scoring, topology avoidance permutation test
    coloring_optimality.py # Hypercube coloring Monte Carlo, rho robustness, per-table
    trna_evidence.py     # tRNA enrichment test (19 organisms, 7 tRNAscan-SE verified)
  data/
    grantham.json              # Grantham 1974 physicochemical distance matrix
    assembly_accessions.tsv    # NCBI genome assemblies for tRNAscan-SE verification
    trnascan_results/          # Raw tRNAscan-SE .out/.stats files (7 organisms)
  visualization/
    data_export.py       # CSV export for R visualization scripts
    R/                   # ggplot2 + ggpubr publication figures (7 panels + combined)
  reports/
    catalogue.py         # Prediction catalogue with evidence grading
    claim_hierarchy.py   # Single source of truth for claim status (15 claims)
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
| **A** | Is the codon->AA assignment special? | Fix degeneracy structure, randomly reassign codons (100k permutations) |
| **B** | Is the block structure special? | Preserve 4-codon blocks, shuffle block->AA assignments |
| **C** | Is the binary encoding special? | Test all 24 possible base->bit-pair mappings |
| **C_extended** | Per-encoding min distances | Emit per-encoding min inter-block Hamming for each disconnected AA |
| **Freeland-Hurst** | Is the coloring optimal? | Block-preserving shuffle with Grantham distance scoring |
| **Class-size** | Weaker coloring null | Class-size-preserving shuffle (no block contiguity) |
| **Table-preserving permutation** | Does evolution avoid topology disruption? | Shuffle reassignment targets within each NCBI table (n=10,000) |

## Technology

- **Python 3.11+**, NumPy, SciPy, requests, click, rich
- **pytest + hypothesis** for property-based testing
- **ggplot2 + ggpubr** (R) for publication figures (300 DPI, colorblind-friendly viridis palette)
- **tRNAscan-SE 2.0.12** + Infernal 1.1.4 for tRNA gene verification
- **SciPy.stats** for permutation tests, Fisher's exact, Stouffer's Z, bootstrap CIs, BH-FDR
