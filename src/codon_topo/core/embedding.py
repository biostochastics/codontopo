"""Holomorphic embedding phi: GF(2)^6 -> C^3.

Maps each nucleotide base to a fourth root of unity (i^k),
then a codon (3 bases) maps to a point in C^3.

The mapping derives from the binary encoding: a base with bit pair
(b0, b1) maps to i^(2*b0 + b1), where i is the imaginary unit.
"""
from codon_topo.core.encoding import DEFAULT_ENCODING

_I = complex(0, 1)
_ROOTS = [_I ** k for k in range(4)]  # [1, i, -1, -i]

# Default base->complex for the standard encoding C=(0,0)->1, U=(0,1)->i, etc.
BASE_TO_COMPLEX: dict[str, complex] = {
    'C': complex(1, 0),    # i^0 = 1
    'U': complex(0, 1),    # i^1 = i
    'A': complex(-1, 0),   # i^2 = -1
    'G': complex(0, -1),   # i^3 = -i
}


def _base_to_complex(
    encoding: dict[str, tuple[int, int]],
) -> dict[str, complex]:
    """Derive base->complex mapping from a binary encoding."""
    return {base: _ROOTS[2 * b0 + b1] for base, (b0, b1) in encoding.items()}


def embed_codon(
    codon: str,
    encoding: dict[str, tuple[int, int]] | None = None,
) -> tuple[complex, complex, complex]:
    """Map a 3-letter codon to C^3 via fourth roots of unity."""
    b2c = _base_to_complex(encoding) if encoding else BASE_TO_COMPLEX
    return (b2c[codon[0]], b2c[codon[1]], b2c[codon[2]])
