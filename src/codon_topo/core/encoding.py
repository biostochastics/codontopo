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


def nucleotide_distance(codon1: str, codon2: str) -> int:
    """Number of nucleotide positions at which two codons differ.

    This is the distance in K4^3 (the complete single-nucleotide mutation
    graph), which is encoding-independent. Two codons are K4^3-adjacent
    iff nucleotide_distance == 1.
    """
    return sum(a != b for a, b in zip(codon1, codon2))


def all_encodings() -> list[dict[str, tuple[int, int]]]:
    """Return all 24 bijections from {C,U,A,G} to {00,01,10,11}."""
    return [dict(zip(BASES, perm)) for perm in permutations(_BIT_PAIRS)]


# Transitions: purine<->purine or pyrimidine<->pyrimidine
_TRANSITIONS = {frozenset(("A", "G")), frozenset(("C", "U"))}


def ts_tv_classification_per_encoding() -> list[dict]:
    """For each of the 24 encodings, classify Hamming-1 edges as Ts or Tv.

    At a single nucleotide position, two bases are Hamming-1 neighbors
    (in the 2-bit encoding) iff they differ at exactly one bit. Each
    encoding produces 4 such Hamming-1 pairs out of the 6 possible base
    pairs. The remaining 2 base pairs are Hamming-2 ("diagonal") edges.

    Returns per-encoding:
      - n_ts_hamming1: how many transition pairs are Hamming-1
      - n_tv_hamming1: how many transversion pairs are Hamming-1
      - n_ts_diagonal: how many transition pairs are Hamming-2
      - n_tv_diagonal: how many transversion pairs are Hamming-2
      - hamming1_pairs: which base pairs are Hamming-1
    """
    results = []
    for enc in all_encodings():
        n_ts_h1 = 0
        n_tv_h1 = 0
        n_ts_diag = 0
        n_tv_diag = 0
        h1_pairs = []
        diag_pairs = []

        for i, b1 in enumerate(BASES):
            for b2 in BASES[i + 1 :]:
                v1 = enc[b1]
                v2 = enc[b2]
                hdist = sum(x != y for x, y in zip(v1, v2))
                is_ts = frozenset((b1, b2)) in _TRANSITIONS

                if hdist == 1:
                    h1_pairs.append((b1, b2))
                    if is_ts:
                        n_ts_h1 += 1
                    else:
                        n_tv_h1 += 1
                else:  # hdist == 2
                    diag_pairs.append((b1, b2))
                    if is_ts:
                        n_ts_diag += 1
                    else:
                        n_tv_diag += 1

        results.append(
            {
                "encoding": dict(enc),
                "n_ts_hamming1": n_ts_h1,
                "n_tv_hamming1": n_tv_h1,
                "n_ts_diagonal": n_ts_diag,
                "n_tv_diagonal": n_tv_diag,
                "hamming1_pairs": h1_pairs,
                "diagonal_pairs": diag_pairs,
            }
        )
    return results
