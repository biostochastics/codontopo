<div align="center">

  # CODON-TOPO

  **Codon Geometry Validation & Prediction Engine**

  [![Version](https://img.shields.io/badge/version-0.4.0-blue)]()
  [![Tests](https://img.shields.io/badge/tests-416%20passing-success)]()
  [![Coverage](https://img.shields.io/badge/coverage-%E2%89%A596%25-brightgreen)]()
  [![Python](https://img.shields.io/badge/python-3.11%2B-yellow)]()
  [![License: CC BY-NC 4.0](https://img.shields.io/badge/license-CC%20BY--NC%204.0-lightgrey)](LICENSE)

</div>

---

## What is CODON-TOPO?

CODON-TOPO validates the algebraic structure of genetic codes when encoded as 6-bit binary vectors in GF(2)^6. It provides a complete, reproducible pipeline for the analyses described in:

> **Robust error-minimization in the genetic code across physicochemical metrics and variant codes: a graph-theoretic analysis in GF(2)^6**
> Paul Clayworth & Sergey Kornilov (2026). Pre-submission manuscript; PDF compiles from `output/manuscript.typ`.

### Key Findings

| Status | Count | Highlights |
|--------|-------|------------|
| **Supported** | 4 | Cross-metric coloring optimality (4 metrics, p ≤ 0.006); per-table preservation (**26 of 27** NCBI tables, mean quantile 1.4%; standard-code-proximity audit confirms variant tables are independently optimal); ρ-robustness across the full Hamming graph H(3,4) = K₄ □ K₄ □ K₄; topology-avoidance depletion under both Q₆ (encoding-dependent) and **encoding-independent H(3,4)** adjacency (RR 0.28–0.33, permutation p ≤ 10⁻⁴, robust to clade exclusion and to both new-disconnection and Δβ₀>0 definitions) |
| **Suggestive** | 1 | tRNA enrichment for reassigned amino acid (worst-case MIS Stouffer p = 0.045 across 24 pairings; 18 tRNAscan-SE–verified genomes); the 4-pairing topology-breaking-restricted subset alone is underpowered (Stouffer p = 0.43) |
| **Exploratory** | 4 | Bit-position bias (deduplicated p = 0.075); mechanism boundary conditions (3-tier: gene duplication / stem shortening / anticodon modification); Atchley F3/Serine convergence; disconnection catalogue (Thr / Leu / Ala / Ser; Trp Table 32 = filtration-only exception) |
| **Rejected** | 3 | Serine min-distance-4 invariant (encoding-dependent); PSL(2,7); holomorphic embedding |
| **Falsified** | 1 | KRAS-Fano clinical prediction (p = 1.0 on n = 1,670 MSK-IMPACT mutations) |
| **Tautological** | 2 | Two-fold bit-5 filtration (encoding-dependent); four-fold prefix filtration |

**Notation:** The full single-nucleotide mutation graph is consistently written as **H(3,4) = K₄ □ K₄ □ K₄** (the Hamming graph; 64 vertices, regular degree 9, 288 undirected edges) rather than the ambiguous K₄³. Q₆ is a 192-edge subgraph of H(3,4); the remaining 96 within-nucleotide diagonal edges complete H(3,4). CLI flags retain the legacy `k43` spelling (e.g. `topology-avoidance-k43`) for backward compatibility.

**Encoding sensitivity (24 base-to-bit bijection sweep):** The Q₆ topology-avoidance result is encoding-dependent — 8 of 24 bijections give a Q₆ candidate-landscape rate near 36% (rather than 73% under the default encoding) and no statistically significant depletion. The H(3,4) result is encoding-independent and **is reported as the primary topology-avoidance test**; Q₆ is now framed as a coordinate-dependent decomposition. Q₆ remains useful for the ρ-sweep (continuous interpolation between Q₆ and H(3,4)).

**Conditional logit (M3 phys+topo) under both topology encodings:** Decisively favored over single-feature models. Under encoding-dependent Q₆ topology: ΔAICc(M1→M3) = 108.2, ΔAICc(M2→M3) = 89.1. Under encoding-independent H(3,4) topology (verifying that the result is not an artifact of the Q₆ encoding): ΔAICc(M1→M3_H(3,4)) = **91.3**, ΔAICc(M2_H(3,4)→M3_H(3,4)) = **95.1** — both decisive (>10) and similar in magnitude to the Q₆ counterparts. Adding the tRNA-distance proxy (M4) does not improve fit (LR = 0.12, p = 0.73). Spearman ρ between Δ_phys and Δ_topo across the 1,280-move candidate landscape = 0.15 (largely independent predictors). Conditional-logit clade-exclusion sensitivity (per Sengupta et al. 2007, refitting M1-M4 with each major clade dropped) and posterior-predictive validation (observed 0.076 vs simulated 0.077; pp p = 0.60) confirm robustness.

**Restricted-candidate sensitivity:** Refitting M1-M4 on candidate sets restricted to biologically plausible moves (target AA already accessible at Hamming distance ≤ d) shows the qualitative claim "topology adds value beyond physicochemistry" survives at every threshold tested. Under the primary d=2 filter (≈727 candidates per choice set), ΔAICc(M1→M3) = 60 and ΔAICc(M2→M3) = 77, both well above the conventional ΔAICc>10 reference. Under the most stringent d=1 filter (≈275 candidates), ΔAICc(M1→M3) shrinks to 14 but stays above 10; ΔAICc(M2→M3) stays at 73. The unrestricted ΔAICc magnitudes are upper bounds; the d=2 filter gives a more biologically-calibrated effect size.

**Methodological caveats explicitly disclosed in Limitations:**
- Survivorship bias: cross-sectional NCBI data cannot distinguish "selection against attempting topology-breaking moves" from "selection against the lineages that attempted them"
- Independence-of-irrelevant-alternatives (IIA) assumption in conditional logit (used as explanatory rather than predictive tool)
- Family-wise multiple-comparison correction within prespecified analysis families (no spurious global-Bonferroni claim)
- Tables 1/11 and 27/28 share identical sense-codon mappings (27 NCBI tables = 25 distinct sense-codon colorings)
- Per-table block-preserving null is partly dominated by near-standard permutations for variants with few reassignments — addressed by standard-code-proximity audit (Supplement)

Run `codon-topo claims` for the full hierarchy with p-values and justifications.

---

## Quick Start

### Prerequisites

- **Python**: 3.11+
- **Package manager**: [uv](https://docs.astral.sh/uv/) (recommended) or pip
- **Optional**: R 4.5+ with `ggplot2`, `ggpubr`, `viridis`, `patchwork` for publication figures
- **Optional**: tRNAscan-SE 2.0.12 for tRNA gene verification

### Installation

```bash
git clone https://github.com/biostochastics/codontopo.git
cd codontopo

# With uv (recommended)
uv sync --all-extras
uv run codon-topo --help

# With pip
pip install -e ".[dev]"
codon-topo --help
```

### Run the Full Pipeline

```bash
# Run everything and generate manuscript_stats.json
codon-topo all --output-dir=./output --seed=135325

# Individual analyses
codon-topo coloring --n=10000          # Coloring optimality Monte Carlo
codon-topo metric-sensitivity          # Cross-metric (Grantham, Miyata, PR, KD)
codon-topo rho-sweep                   # Rho robustness (Q6 -> K4^3)
codon-topo per-table                   # All 27 NCBI translation tables
codon-topo topology-avoidance          # Topology avoidance (Q6)
codon-topo topology-avoidance-k43      # Topology avoidance (K4^3, encoding-independent)
codon-topo condlogit                   # Conditional logit models (M1-M4)
codon-topo condlogit-restricted        # Restricted-candidate sensitivity (delta_trna<=1,2,3)
codon-topo trna                        # tRNA enrichment test
codon-topo mis-analysis                # Maximal independent set analysis
codon-topo phylo-sensitivity           # Clade-exclusion robustness
codon-topo claims                      # View claim hierarchy
```

### Run the CodonSafe Cross-Study Reanalysis

```bash
# Requires raw data in data/codonsafe/ (see DATA_MANIFEST.md)
pip install -e ".[codonsafe]"
codon-topo codonsafe
```

### Run the Test Suite

```bash
python3.11 -m pytest tests/ -q                    # all tests
python3.11 -m pytest tests/ --cov=codon_topo      # with coverage
python3.11 -m pytest tests/test_regression.py -v   # regression suite (105 tests)
```

> **Note**: Use `python3.11 -m pytest` if your system default Python differs from where dev dependencies are installed.

### Generate Publication Figures

```bash
Rscript src/codon_topo/visualization/R/strengthened_figures.R
```

---

## Reproducibility

The core design principle: **a user who clones this repo should be able to regenerate every number in the manuscript.**

```bash
# Full reproducibility from scratch
git clone https://github.com/biostochastics/codontopo.git
cd codontopo
uv sync --all-extras
uv run codon-topo all --output-dir=./output --seed=135325
# -> generates output/manuscript_stats.json
# -> manuscript.typ reads all inline statistics from this JSON
```

The `manuscript_stats.json` file contains every statistic cited in the paper. The Typst manuscript (`output/manuscript.typ`) reads this file and renders all inline numbers dynamically:

```typst
#let stats = json("manuscript_stats.json")
// All tables and inline stats reference stats.* fields
```

Random seed: **135325** (all Monte Carlo analyses).

---

## Architecture

```
codon-topo all
    |
    +-- Filtration (WS1) .................. Two-fold/four-fold degeneracy checks
    +-- Disconnections (WS1) .............. Persistent homology catalogue
    +-- Coloring Optimality (WS1) ......... Block-preserving Monte Carlo
    |     +-- Multi-metric sensitivity .... Grantham, Miyata, PR, KD
    |     +-- Rho robustness sweep ........ Q6 -> H(3,4) interpolation
    |     +-- Per-table optimality ........ 27 NCBI tables + BH-FDR
    |     +-- Per-table proximity audit ... dH-conditional vs unconditional quantile
    |     +-- Score decomposition ......... By nucleotide position
    +-- Reassignment Analysis (WS2) ....... Database, Hamming paths, bit bias
    +-- Topology Avoidance (WS6) .......... Q6 + H(3,4), 2x2 definitions audit,
    |                                       24-encoding sweep, denominator sensitivity
    +-- tRNA Evidence (WS1) ............... Fisher-Stouffer + MIS enumeration
    |                                       + topology-breaking subset (n=4)
    +-- Phylogenetic Sensitivity (WS6) .... Clade-exclusion robustness
    +-- Conditional Logit (WS6) ........... M1-M4 (Q6) + M2k43, M3k43 (H(3,4))
    |     +-- Encoding robustness ......... Q6 vs H(3,4) ΔAICc comparison
    |     +-- Clade-exclusion sensitivity . 7 clade regimes per Sengupta et al. 2007
    |     +-- Restricted-candidate sens. .. delta_trna<=1,2,3 biological-plausibility filter
    |     +-- Posterior predictive ........ Observed vs simulated topology rate
    +-- Depth Calibration (WS3) ........... Epsilon-age correlation
    +-- KRAS-Fano (WS4) .................. cBioPortal enrichment (negative)
    +-- Claims + Catalogue (WS5) .......... 15 claims, evidence grading
    |
    +-> output/manuscript_stats.json ...... Consolidated stats for Typst
    +-> output/*.json ..................... Per-analysis detailed results
```

### Package Structure

| Component | Path | Role |
|-----------|------|------|
| CLI | `src/codon_topo/cli.py` | Click-based CLI with 18 subcommands |
| Encoding | `src/codon_topo/core/encoding.py` | GF(2)^6, Hamming distance, all 24 encodings |
| Genetic codes | `src/codon_topo/core/genetic_codes.py` | All 27 NCBI translation tables (codes 1-6, 9-16, 21-33) |
| Filtration | `src/codon_topo/core/filtration.py` | Two-fold (bit-5) and four-fold (prefix) checks |
| Homology | `src/codon_topo/core/homology.py` | Connected components, disconnection catalogue |
| Embedding | `src/codon_topo/core/embedding.py` | Root-of-unity map GF(2)^6 -> C^3 |
| Fano | `src/codon_topo/core/fano.py` | XOR triple computation |
| Coloring optimality | `src/codon_topo/analysis/coloring_optimality.py` | Monte Carlo, rho sweep, per-table, multi-metric |
| Null models | `src/codon_topo/analysis/null_models.py` | Models A/B/C/C_extended |
| Reassignment DB | `src/codon_topo/analysis/reassignment_db.py` | Database, Hamming paths, bit-position bias |
| Topology avoidance | `src/codon_topo/analysis/synbio_feasibility.py` | Q6 + K4^3 tests, phylogenetic sensitivity |
| Evolutionary simulation | `src/codon_topo/analysis/evolutionary_simulation.py` | Conditional logit M1-M4, order-averaging |
| tRNA evidence | `src/codon_topo/analysis/trna_evidence.py` | Fisher-Stouffer, MIS via Bron-Kerbosch |
| CodonSafe | `src/codon_topo/analysis/codonsafe/` | Cross-study reanalysis of 8 recoding datasets |
| Statistical utils | `src/codon_topo/analysis/statistical_utils.py` | Beta CIs, risk ratios, quantile CIs |
| Visualization | `src/codon_topo/visualization/` | CSV export + R ggplot2 scripts |
| Claims | `src/codon_topo/reports/claim_hierarchy.py` | Single source of truth for 15 claims |
| Catalogue | `src/codon_topo/reports/catalogue.py` | Evidence grading across workstreams |

---

## CLI Reference

| Command | Description |
|---------|-------------|
| `codon-topo all` | Run everything, generate `manuscript_stats.json` |
| `codon-topo filtration` | Two-fold/four-fold filtration checks |
| `codon-topo disconnections` | Disconnection catalogue (persistent homology) |
| `codon-topo coloring` | Hypercube coloring Monte Carlo |
| `codon-topo metric-sensitivity` | Cross-metric sensitivity (4 metrics) |
| `codon-topo rho-sweep` | Rho robustness (Q6 -> K4^3) |
| `codon-topo per-table` | Per-table optimality (27 NCBI tables) |
| `codon-topo decompose` | Score decomposition by nucleotide position |
| `codon-topo topology-avoidance` | Topology avoidance test (Q6) |
| `codon-topo topology-avoidance-k43` | Topology avoidance test (K4^3) |
| `codon-topo condlogit` | Conditional logit model comparison (M1-M4) |
| `codon-topo condlogit-restricted` | Restricted-candidate-set sensitivity (delta_trna ≤ d) |
| `codon-topo phylo-sensitivity` | Clade-exclusion sensitivity analysis |
| `codon-topo trna` | tRNA enrichment test |
| `codon-topo mis-analysis` | Maximal independent set enumeration |
| `codon-topo bit-bias` | Bit-position bias test |
| `codon-topo kras` | KRAS-Fano enrichment test (negative) |
| `codon-topo codonsafe` | CodonSafe cross-study reanalysis |
| `codon-topo claims` | View claim hierarchy |

All subcommands support `--json` for machine-readable output. Interactive mode uses rich tables.

---

## Workstreams

| WS | Name | CLI commands | Status |
|----|------|-------------|--------|
| **WS1** | Core Replication | `filtration`, `disconnections`, `coloring`, `metric-sensitivity`, `rho-sweep`, `per-table`, `decompose` | Complete |
| **WS2** | Reassignment Directionality | `bit-bias` | Complete |
| **WS3** | Evolutionary Depth | _(in `all`)_ | Complete |
| **WS4** | KRAS/COSMIC | `kras` | Complete (negative) |
| **WS5** | Prediction Catalogue | `claims` | Complete |
| **WS6** | Topology & Synbio | `topology-avoidance`, `topology-avoidance-k43`, `condlogit`, `phylo-sensitivity`, `codonsafe` | Complete |

---

## Null Models

| Model | What it tests | CLI |
|-------|---------------|-----|
| **Freeland-Hurst** | Is the coloring optimal? Block-preserving shuffle | `coloring` |
| **Class-size** | Weaker null (degeneracy-only, no block contiguity) | `coloring --null=class_size` |
| **Model C** | Is the encoding special? All 24 base-to-bit mappings | `disconnections --extended` |
| **Table-preserving permutation** | Does evolution avoid topology disruption? | `topology-avoidance` |
| **Conditional logit** | Is topology an independent predictor? | `condlogit` |

---

## Technology Stack

- **Python 3.11+**, NumPy, SciPy for core computation
- **click + rich** for CLI
- **pytest + hypothesis** for property-based testing (416 tests, >=96% coverage)
- **ggplot2 + ggpubr** (R) for publication figures (300 DPI, colorblind-friendly viridis)
- **Typst** for manuscript typesetting (reads `manuscript_stats.json` for dynamic stats)
- **tRNAscan-SE 2.0.12** + Infernal 1.1.4 for tRNA verification (18 genomes across 5 variant codes + 3 standard-code controls)
- **Biopython** for GenBank parsing (CodonSafe reanalysis)

---

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

# Run the coloring optimality Monte Carlo
result = monte_carlo_null(n_samples=10000, seed=135325)
# {'quantile_of_observed': 0.6, 'p_value_conservative': 0.006, ...}

# Query the claim hierarchy
for claim in supported_claims():
    print(claim.id, claim.evidence_p_value)
```

---

## Documentation

| Document | Purpose |
|----------|---------|
| [`CLAUDE.md`](CLAUDE.md) | AI/contributor guidance |
| [`ARCHITECTURE.md`](ARCHITECTURE.md) | Module dependency graph |
| [`data/codonsafe/DATA_MANIFEST.md`](data/codonsafe/DATA_MANIFEST.md) | Raw data provenance for cross-study reanalysis |

---

## License

Released under the [Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)](LICENSE) license. You may share and adapt this work with attribution, but commercial use requires a separate license — contact the authors.

To cite, see [`CITATION.cff`](CITATION.cff) or the bibliography entry generated by GitHub's "Cite this repository" button.

---

<div align="center">

**[Quick Start](#quick-start)** &bull;
**[CLI Reference](#cli-reference)** &bull;
**[Reproducibility](#reproducibility)** &bull;
**[Architecture](#architecture)**

</div>
