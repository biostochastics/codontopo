"""Fano-line computation in GF(2)^6.

A Fano line is a triple of vectors (a, b, c) such that a XOR b XOR c = 0.
Equivalently, c = a XOR b.
"""
from codon_topo.core.encoding import (
    codon_to_vector, ALL_CODONS, DEFAULT_ENCODING,
)

# Reverse lookup: vector -> codon (default encoding)
_VECTOR_TO_CODON: dict[tuple[int, ...], str] = {
    codon_to_vector(c): c for c in ALL_CODONS
}


def _xor(a: tuple[int, ...], b: tuple[int, ...]) -> tuple[int, ...]:
    return tuple((x + y) % 2 for x, y in zip(a, b))


def is_fano_line(c1: str, c2: str, c3: str) -> bool:
    """Check if three codons form a Fano line (XOR = 0)."""
    v1, v2, v3 = codon_to_vector(c1), codon_to_vector(c2), codon_to_vector(c3)
    return all(x == 0 for x in _xor(_xor(v1, v2), v3))


def fano_partner(c1: str, c2: str) -> str:
    """Return the third codon completing the Fano line through c1 and c2."""
    v1, v2 = codon_to_vector(c1), codon_to_vector(c2)
    v3 = _xor(v1, v2)
    return _VECTOR_TO_CODON[v3]


def all_single_bit_fano_lines(codon: str) -> list[dict]:
    """For each single-bit mutation from codon, compute the Fano partner.

    Returns list of dicts: {bit_pos, mutant_codon, fano_partner}.
    """
    v = codon_to_vector(codon)
    results = []
    for bit_pos in range(6):
        mutant = list(v)
        mutant[bit_pos] = 1 - mutant[bit_pos]
        mutant = tuple(mutant)
        partner = _xor(v, mutant)
        results.append({
            'bit_pos': bit_pos,
            'mutant_codon': _VECTOR_TO_CODON[mutant],
            'fano_partner': _VECTOR_TO_CODON[partner],
        })
    return results
