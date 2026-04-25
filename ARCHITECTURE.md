# Architecture

## Module Dependency Graph

```
cli.py
  -> core.filtration
  -> core.homology
  -> core.genetic_codes
  -> analysis.coloring_optimality
  -> analysis.null_models
  -> analysis.reassignment_db
  -> analysis.trna_evidence
  -> analysis.depth_calibration
  -> analysis.cosmic_query
  -> analysis.synbio_feasibility
  -> reports.claim_hierarchy
  -> reports.catalogue

core.filtration       -> core.encoding
core.homology         -> core.encoding
core.embedding        -> core.encoding
core.fano             -> core.encoding
core.genetic_codes    (standalone, no internal deps)
core.encoding         (standalone, no internal deps)

analysis.null_models           -> core.encoding, core.genetic_codes, core.homology
analysis.coloring_optimality   -> core.encoding, core.genetic_codes
analysis.reassignment_db       -> core.encoding, core.genetic_codes
analysis.trna_evidence         (standalone, curated data only)
analysis.depth_calibration     (standalone, curated data only)
analysis.cosmic_query          -> core.genetic_codes, core.fano
analysis.synbio_feasibility    -> core.filtration, core.homology, core.genetic_codes,
                                  analysis.reassignment_db

reports.claim_hierarchy  (standalone, no internal deps)
reports.catalogue        (standalone, no internal deps)

visualization.data_export -> analysis.*, reports.*
visualization.R/          -> output/tables/*.csv, output/tables/*.json
```

## Data Flow

```
NCBI Translation Tables (27 tables: 1–6, 9–16, 21–33)
  |
  v
core.genetic_codes.get_code(table_id) -> dict[str, str]
  |
  +----> core.encoding.codon_to_vector() -> tuple[int, ...]
  |        |
  |        +----> core.filtration: two-fold / four-fold checks
  |        +----> core.homology: connected components, disconnection catalogue
  |        +----> core.embedding: coordinate-wise root-of-unity map
  |        +----> core.fano: XOR triples
  |
  +----> analysis.null_models: A (random), B (block shuffle), C (encodings)
  +----> analysis.coloring_optimality: Grantham-weighted edge mismatch, rho sweep, per-table
  +----> analysis.reassignment_db: cross-table diff -> Hamming paths, bit-position bias
  +----> analysis.synbio_feasibility: topology avoidance (permutation null)
  +----> analysis.trna_evidence: tRNA enrichment (Fisher+Stouffer, 18 organisms; 17 tRNAscan-SE verified + 1 from literature)
  |
  v
reports.claim_hierarchy: single source of truth for 15 claims
reports.catalogue: structured prediction list with evidence grading
  |
  v
cli.py: subcommands -> JSON + rich tables
scripts/generate_tables.py: CSV/JSON -> output/tables/
visualization.R/: ggplot2 + ggpubr -> output/figures/
```

## Null Model Taxonomy

| ID | Name | Preserves | Scores | Usage |
|----|------|-----------|--------|-------|
| A | Random assignment | Degeneracy sizes | Filtration + homology | `null_models.null_model_a()` |
| B | Block shuffle | 4-codon block structure | Serine uniqueness | `null_models.null_model_b()` |
| C | Encoding sweep | Everything except encoding | Filtration invariance | `null_models.null_model_c()` |
| C_ext | Encoding + distances | Everything except encoding | Per-AA min inter-block distance | `null_models.null_model_c_extended()` |
| FH | Freeland-Hurst block-preserving | Block contiguity, stop positions | Grantham edge mismatch | `coloring_optimality.monte_carlo_null(null_type="freeland_hurst")` |
| CS | Class-size-preserving | Degeneracy sizes (no block contiguity) | Grantham edge mismatch | `coloring_optimality.monte_carlo_null(null_type="class_size")` |
| TP | Table-preserving permutation | Within-table codon set | Novel disconnection rate | `synbio_feasibility.topology_avoidance_test()` |
| CL | Conditional logit (M1-M4) | Event-level candidate features (Δ_local, topology, Hamming, tRNA-proxy) | Independent explanatory power via AICc / LRT | `evolutionary_simulation.run_evolutionary_simulation_analysis()` |

## Claim Hierarchy

`reports/claim_hierarchy.py` is the single source of truth. Each claim has:
- `ClaimStatus`: SUPPORTED, SUGGESTIVE, EXPLORATORY, REJECTED, FALSIFIED, TAUTOLOGICAL
- Evidence: p-value, null model, sample size
- Justification: why this status was assigned

Current state (15 claims): 4 SUPPORTED, 1 SUGGESTIVE, 4 EXPLORATORY, 3 REJECTED, 1 FALSIFIED, 2 TAUTOLOGICAL.

Run `codon-topo claims` to view.

## Data Provenance

tRNA gene counts verified by tRNAscan-SE 2.0.12 on NCBI genome assemblies (representative sample shown; complete 18-organism dataset is in `data/assembly_accessions.tsv` and `data/trnascan_results/`):

| Organism | Assembly | tRNAscan Total | Role |
|----------|----------|----------------|------|
| T. thermophila | GCF_000189635.1 | 718 | Variant (Table 6, Gln) |
| P. tetraurelia | GCF_000165425.1 | 216 | Variant (Table 6, Gln) |
| O. trifallax | GCA_000295675.1 | 94 | Variant (Table 6, Gln) |
| S. coeruleus | GCA_001970955.1 | 272 | Control (Standard) |
| I. multifiliis | GCF_000220395.1 | 150 | Control (Standard) |
| B. nonstop P57 | GCA_028554745.1 | 68 | Boundary (Table 31, Trp) |

Raw results in `data/trnascan_results/`. Reproduce with `bash scripts/run_trnascan.sh`.

## Output Artifacts

| Path | Description |
|------|-------------|
| `output/tables/T1_claim_hierarchy.csv` | All 15 claims with status and p-values |
| `output/tables/T3_coloring_optimality.json` | Monte Carlo null distribution (n=10,000) |
| `output/tables/T4_per_table_optimality.csv` | Per-table quantiles (27 tables) |
| `output/tables/T5_rho_robustness.csv` | Rho sweep (11 values, 0.0 to 1.0) |
| `output/tables/T7_trna_per_pairing.csv` | Fisher exact per pairing (11 rows) |
| `output/tables/T9_topology_avoidance.csv` | Permutation + hypergeometric p-values |
| `output/figures/FigA-G_*.png` | 7 publication panels (ggplot2, 300 DPI) |
| `output/figures/panel_strengthened.png` | Combined 6-panel figure |
| `output/hypercube_coloring.html` | Interactive D3 visualization |
