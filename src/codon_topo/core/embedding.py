"""Holomorphic embedding phi: GF(2)^6 -> C^3.

Maps each nucleotide base to a fourth root of unity (i^k),
then a codon (3 bases) maps to a point in C^3.
"""

BASE_TO_COMPLEX: dict[str, complex] = {
    'C': complex(1, 0),    # i^0 = 1
    'U': complex(0, 1),    # i^1 = i
    'A': complex(-1, 0),   # i^2 = -1
    'G': complex(0, -1),   # i^3 = -i
}


def embed_codon(codon: str) -> tuple[complex, complex, complex]:
    """Map a 3-letter codon to C^3 via fourth roots of unity."""
    return (
        BASE_TO_COMPLEX[codon[0]],
        BASE_TO_COMPLEX[codon[1]],
        BASE_TO_COMPLEX[codon[2]],
    )
