# CODON-TOPO

Codon Geometry Validation & Prediction Engine. Validates the algebraic structure of genetic codes when encoded as 6-bit binary vectors in GF(2)^6.

## What This Does

The genetic code maps 64 codons to 20 amino acids + Stop. When you encode each nucleotide as a 2-bit pair (C=00, U=01, A=10, G=11), each codon becomes a 6-bit vector in GF(2)^6. This package proves that the resulting algebraic structure has remarkable properties:

- **Two-fold filtration**: All 9 two-fold degenerate amino acids differ at exactly bit 5 (the wobble position parity bit) — across ALL 25 known genetic codes
- **Four-fold filtration**: All 5 four-fold degenerate amino acids share an identical 4-bit prefix, with the last 2 bits exhausting GF(2)^2
- **Serine invariant**: Serine is the *only* amino acid whose codon graph is topologically disconnected at Hamming distance 1 — a universal invariant across every known genetic code
- **KRAS Fano line**: The cancer-relevant KRAS G12V mutation (GGU→GUU) forms a Fano line with CAC (Histidine): GGU ⊕ GUU ⊕ CAC = 0 in GF(2)^6
- **Holomorphic embedding**: Mapping bases to fourth roots of unity (C→1, U→i, A→-1, G→-i) sends codons to C^3, preserving degeneracy structure

## Quick Start

```bash
# Install
pip install -e ".[dev]"

# Run all 221 tests
python3.11 -m pytest

# Run with coverage (99%)
python3.11 -m pytest --cov=codon_topo --cov-report=term-missing

# Run just the regression suite (105 tests against PRD Appendix 8)
python3.11 -m pytest tests/test_regression.py -v
```

> **Note**: Use `python3.11 -m pytest` if your system default Python differs from where dev dependencies are installed.

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
  visualization/
    data_export.py       # CSV export for R visualization scripts
    R/                   # ggplot2 + ggpubr scripts for publication figures
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

# Run null model A (statistical validation)
from codon_topo.analysis.null_models import null_model_a
result = null_model_a(n_permutations=10_000, seed=42)
print(f"p(Serine unique) = {result['p_value_serine_unique']}")
print(f"p(bit-5 uniform) = {result['p_value_bit5_uniform']}")

# Export data for R visualization
from codon_topo.visualization.data_export import export_persistent_homology
export_persistent_homology('output/persistent_homology.csv')
# Then: Rscript src/codon_topo/visualization/R/barcode_plot.R output/persistent_homology.csv output/barcode.pdf
```

## Null Models

| Model | What it tests | Method |
|-------|---------------|--------|
| **A** | Is the codon→AA assignment special? | Fix degeneracy structure, randomly reassign codons (100k permutations) |
| **B** | Is the block structure special? | Preserve 4-codon blocks, shuffle block→AA assignments |
| **C** | Is the binary encoding special? | Test all 24 possible base→bit-pair mappings |

## Technology

- **Python 3.11+**, NumPy, SciPy
- **pytest + hypothesis** for property-based testing
- **ggplot2 + ggpubr** (R) for publication figures
- **ruff** for linting/formatting, **mypy** for type checking
