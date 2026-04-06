"""Binary encoding of codons as vectors in GF(2)^6.

The default encoding maps each nucleotide to a 2-bit pair:
C=(0,0), U=(0,1), A=(1,0), G=(1,1). A codon (3 bases) becomes
a 6-bit tuple. All 24 permutations of the 4 possible 2-bit values
across the 4 bases are available for Null Model C analysis.
"""

from itertools import permutations

BASES = ("C", "U", "A", "G")
_BIT_PAIRS = ((0, 0), (0, 1), (1, 0), (1, 1))

DEFAULT_ENCODING: dict[str, tuple[int, int]] = dict(zip(BASES, _BIT_PAIRS))

ALL_CODONS: list[str] = [b1 + b2 + b3 for b1 in BASES for b2 in BASES for b3 in BASES]


def codon_to_vector(
    codon: str,
    encoding: dict[str, tuple[int, int]] | None = None,
) -> tuple[int, ...]:
    """Convert a 3-letter codon to a 6-bit tuple in GF(2)^6."""
    enc = encoding or DEFAULT_ENCODING
    bits: list[int] = []
    for base in codon:
        bits.extend(enc[base])
    return tuple(bits)


def hamming_distance(a: tuple[int, ...], b: tuple[int, ...]) -> int:
    """Hamming distance between two bit-tuples."""
    return sum(x != y for x, y in zip(a, b))


def all_encodings() -> list[dict[str, tuple[int, int]]]:
    """Return all 24 bijections from {C,U,A,G} to {00,01,10,11}."""
    return [dict(zip(BASES, perm)) for perm in permutations(_BIT_PAIRS)]
