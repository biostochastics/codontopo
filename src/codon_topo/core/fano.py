"""Fano-line computation in GF(2)^6.

A Fano line is a triple of vectors (a, b, c) such that a XOR b XOR c = 0.
Equivalently, c = a XOR b.
"""
from codon_topo.core.encoding import (
    codon_to_vector, ALL_CODONS, DEFAULT_ENCODING,
)

# Cache for reverse lookups per encoding (encoding values are unhashable tuples,
# so we key on the frozenset of items).
_REVERSE_CACHE: dict[frozenset, dict[tuple[int, ...], str]] = {}


def _reverse_lookup(
    encoding: dict[str, tuple[int, int]] | None = None,
) -> dict[tuple[int, ...], str]:
    """Return vector->codon mapping for a given encoding, with caching."""
    enc = encoding or DEFAULT_ENCODING
    key = frozenset(enc.items())
    if key not in _REVERSE_CACHE:
        _REVERSE_CACHE[key] = {
            codon_to_vector(c, enc): c for c in ALL_CODONS
        }
    return _REVERSE_CACHE[key]


def _xor(a: tuple[int, ...], b: tuple[int, ...]) -> tuple[int, ...]:
    return tuple((x + y) % 2 for x, y in zip(a, b))


def is_fano_line(
    c1: str, c2: str, c3: str,
    encoding: dict[str, tuple[int, int]] | None = None,
) -> bool:
    """Check if three codons form a Fano line (XOR = 0)."""
    enc = encoding or DEFAULT_ENCODING
    v1, v2, v3 = codon_to_vector(c1, enc), codon_to_vector(c2, enc), codon_to_vector(c3, enc)
    return all(x == 0 for x in _xor(_xor(v1, v2), v3))


def fano_partner(
    c1: str, c2: str,
    encoding: dict[str, tuple[int, int]] | None = None,
) -> str:
    """Return the third codon completing the Fano line through c1 and c2."""
    enc = encoding or DEFAULT_ENCODING
    v1, v2 = codon_to_vector(c1, enc), codon_to_vector(c2, enc)
    v3 = _xor(v1, v2)
    return _reverse_lookup(enc)[v3]


def all_single_bit_fano_lines(
    codon: str,
    encoding: dict[str, tuple[int, int]] | None = None,
) -> list[dict]:
    """For each single-bit mutation from codon, compute the Fano partner.

    Returns list of dicts: {bit_pos, mutant_codon, fano_partner}.
    """
    enc = encoding or DEFAULT_ENCODING
    rev = _reverse_lookup(enc)
    v = codon_to_vector(codon, enc)
    results = []
    for bit_pos in range(6):
        mutant = list(v)
        mutant[bit_pos] = 1 - mutant[bit_pos]
        mutant = tuple(mutant)
        partner = _xor(v, mutant)
        results.append({
            'bit_pos': bit_pos,
            'mutant_codon': rev[mutant],
            'fano_partner': rev[partner],
        })
    return results
