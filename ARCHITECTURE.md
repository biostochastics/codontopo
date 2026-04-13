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
analysis.synbio_feasibility    -> core.filtration, core.homology, core.genetic_codes

reports.claim_hierarchy  (standalone, no internal deps)
reports.catalogue        (standalone, no internal deps)

visualization.data_export -> analysis.*, reports.*
```

## Data Flow

```
NCBI Translation Tables (25 tables)
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
  +----> analysis.coloring_optimality: Grantham-weighted edge mismatch
  +----> analysis.reassignment_db: cross-table diff -> Hamming paths
  |
  v
reports.claim_hierarchy: single source of truth for claim status
reports.catalogue: structured prediction list with evidence grading
  |
  v
cli.py: subcommands -> JSON + rich tables
visualization.data_export: CSV -> R/ggplot2
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

## Claim Hierarchy

`reports/claim_hierarchy.py` is the single source of truth for what the paper claims. Each claim has:
- `ClaimStatus`: SUPPORTED, SUGGESTIVE, EXPLORATORY, REJECTED, FALSIFIED, TAUTOLOGICAL
- Evidence: p-value, null model, sample size
- Justification: why this status was assigned

The hierarchy was established through multi-model adversarial review (April 2026). Run `codon-topo claims` to view it.
