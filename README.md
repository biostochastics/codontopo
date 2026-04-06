# CODON-TOPO

Codon Geometry Validation & Prediction Engine. Validates the algebraic structure of genetic codes when encoded as 6-bit binary vectors in GF(2)^6.

## What This Does

The genetic code maps 64 codons to 20 amino acids + Stop. When you encode each nucleotide as a 2-bit pair (C=00, U=01, A=10, G=11), each codon becomes a 6-bit vector in GF(2)^6. This package proves that the resulting algebraic structure has remarkable properties:

- **Two-fold filtration**: All 9 two-fold degenerate amino acids differ at exactly bit 5 (the wobble position parity bit) — across ALL 25 known genetic codes
- **Four-fold filtration**: All 5 four-fold degenerate amino acids share an identical 4-bit prefix, with the last 2 bits exhausting GF(2)^2
- **Serine invariant**: Serine is the *only* amino acid whose codon graph is topologically disconnected at Hamming distance 1 — a universal invariant across every known genetic code
- **KRAS Fano line**: The cancer-relevant KRAS G12V mutation (GGU->GUU) forms a Fano line with CAC (Histidine): GGU XOR GUU XOR CAC = 0 in GF(2)^6
- **Holomorphic embedding**: Mapping bases to fourth roots of unity (C->1, U->i, A->-1, G->-i) sends codons to C^3, preserving degeneracy structure

## Quick Start

```bash
# Install
pip install -e ".[dev]"

# Run all 280 tests
python3.11 -m pytest

# Run with coverage (96%)
python3.11 -m pytest --cov=codon_topo --cov-report=term-missing

# Run just the regression suite (105 tests against PRD Appendix 8)
python3.11 -m pytest tests/test_regression.py -v
```

> **Note**: Use `python3.11 -m pytest` if your system default Python differs from where dev dependencies are installed.

## Workstreams

| WS | Name | Status | Description |
|----|------|--------|-------------|
| **WS1** | Core Replication | Complete | Filtration, homology, embedding, Fano lines, null models A/B/C |
| **WS2** | Reassignment Directionality | Complete | Database of ~55 codon reassignment events, Hamming path analysis, bit-position bias test |
| **WS3** | Evolutionary Depth Calibration | Complete | Epsilon-age Spearman correlation across 6 calibration points (Serine 3.5 Gya to CUG clade 150 Mya) |
| **WS4** | KRAS/COSMIC Fast-Fail Gate | Complete | cBioPortal API client, Fano-line co-occurrence prediction, Fisher's exact test with Bonferroni correction |
| **WS5** | Prediction Catalogue | Complete | 15 predictions across all workstreams with evidence grading (verified/tested/pending) |
| **WS6** | Synthetic Biology Feasibility | Complete | Feasibility scoring for ~1280 single-codon reassignment variants |

## Package Structure

```
src/codon_topo/
  core/
    encoding.py          # GF(2)^6 binary encoding, Hamming distance, all 24 encodings
    genetic_codes.py     # All 25 NCBI translation tables
    filtration.py        # Two-fold (bit-5) and four-fold (prefix) degeneracy checks
    homology.py          # Persistent homology via union-find, disconnection catalogue
    embedding.py         # Holomorphic embedding phi: GF(2)^6 -> C^3
    fano.py              # Fano-line (XOR triple) computation
  analysis/
    null_models.py       # Statistical null models A (random), B (block shuffle), C (24 encodings)
    reassignment_db.py   # Reassignment database, Hamming path analysis, directionality stats
    cosmic_query.py      # cBioPortal API client, KRAS Fano co-occurrence, WS4 gate decision
    depth_calibration.py # Evolutionary depth calibration, epsilon-age correlation
    synbio_feasibility.py # Feasibility scoring for alternative genetic codes
  visualization/
    data_export.py       # CSV export for R visualization scripts (all workstreams)
  reports/
    catalogue.py         # Prediction catalogue with evidence grading
```

## Usage Examples

```python
from codon_topo import (
    codon_to_vector, hamming_distance, STANDARD, get_code,
    analyze_filtration, disconnection_catalogue, embed_codon,
    is_fano_line, fano_partner,
)

# Encode a codon as a 6-bit vector
codon_to_vector('GGU')  # (1, 1, 1, 1, 0, 1)

# Check the KRAS Fano line
is_fano_line('GGU', 'GUU', 'CAC')  # True (XOR = 0)
fano_partner('GGU', 'GUU')  # 'CAC'

# Analyze filtration for any genetic code
report = analyze_filtration(STANDARD)
# {'twofold_pass': 9, 'twofold_fail': 0, 'fourfold_pass': 5, ...}

# Find all disconnected amino acids
catalogue = disconnection_catalogue(STANDARD)
# [{'aa': 'Ser', 'n_components': 2, 'reconnect_eps': 4, 'min_inter_distance': 4, ...}]

# Embed a codon in C^3
embed_codon('GGU')  # (-i, -i, i)
```

### WS2-WS6 Examples

```python
# WS2: Reassignment database
from codon_topo import build_reassignment_db
db = build_reassignment_db()
# 55 events: ReassignmentEvent(table_id=2, codon='AGA', source_aa='Arg', target_aa='Stop', ...)

# WS3: Evolutionary depth calibration
from codon_topo import compute_correlation
result = compute_correlation()
# {'spearman_rho': ..., 'spearman_p': ..., 'n_points': 6}

# WS4: KRAS Fano predictions
from codon_topo import fano_predictions_for_kras
preds = fano_predictions_for_kras()
# {'G12V': {'fano_partner_codon': 'CAC', 'fano_partner_aa': 'His', ...}, ...}

# WS4: Run the decision gate (requires cBioPortal data or mock)
from codon_topo import ws4_gate_decision
result = ws4_gate_decision(mutation_data, p_threshold=0.01)
# {'pass': True/False, 'corrected_threshold': 0.00167, 'significant_variants': [...]}

# WS6: Score a variant genetic code
from codon_topo import score_variant_code, get_code
score = score_variant_code(get_code(3))  # Yeast mitochondrial
# {'feasibility_score': 0.55, 'serine_disconnected': True, 'n_disconnected_aas': 2, ...}

# WS5: Full prediction catalogue
from codon_topo import build_catalogue
cat = build_catalogue()  # 15 predictions with evidence grading

# Export any workstream data for R visualization
from codon_topo.visualization.data_export import (
    export_reassignment_db, export_depth_calibration,
    export_fano_predictions, export_catalogue,
)
export_reassignment_db('output/reassignments.csv')
export_depth_calibration('output/depth_calibration.csv')
export_fano_predictions('output/fano_predictions.csv')
export_catalogue('output/catalogue.csv')
```

## Null Models

| Model | What it tests | Method |
|-------|---------------|--------|
| **A** | Is the codon->AA assignment special? | Fix degeneracy structure, randomly reassign codons (100k permutations) |
| **B** | Is the block structure special? | Preserve 4-codon blocks, shuffle block->AA assignments |
| **C** | Is the binary encoding special? | Test all 24 possible base->bit-pair mappings |

## Technology

- **Python 3.11+**, NumPy, SciPy, requests
- **pytest + hypothesis** for property-based testing
- **ggplot2 + ggpubr** (R) for publication figures
- **ruff** for linting/formatting, **mypy** for type checking
- **SciPy.stats** for Spearman correlation, Fisher's exact test, chi-squared, bootstrap CIs
