from hypothesis import given, strategies as st
from codon_topo.core.encoding import (
    DEFAULT_ENCODING,
    codon_to_vector,
    hamming_distance,
    ALL_CODONS,
    BASES,
    all_encodings,
)


def test_default_encoding():
    assert DEFAULT_ENCODING == {"C": (0, 0), "U": (0, 1), "A": (1, 0), "G": (1, 1)}


def test_codon_to_vector_ggu():
    """GGU = (1,1,1,1,0,1) per verify_kras.py line 30-31."""
    assert codon_to_vector("GGU") == (1, 1, 1, 1, 0, 1)


def test_codon_to_vector_cac():
    """CAC = (0,0,1,0,0,0) per verify_kras.py."""
    assert codon_to_vector("CAC") == (0, 0, 1, 0, 0, 0)


def test_all_codons_count():
    assert len(ALL_CODONS) == 64


def test_all_codons_unique():
    assert len(set(ALL_CODONS)) == 64


def test_hamming_self():
    v = codon_to_vector("GGU")
    assert hamming_distance(v, v) == 0


def test_hamming_ggu_guu():
    """d(GGU, GUU) = 1 per verify_kras.py line 34."""
    ggu = codon_to_vector("GGU")
    guu = codon_to_vector("GUU")
    assert hamming_distance(ggu, guu) == 1


@given(st.sampled_from(["UUU", "UUC", "UUA", "UUG", "GGG", "AAA", "CCC"]))
def test_vector_length_is_6(codon):
    v = codon_to_vector(codon)
    assert len(v) == 6
    assert all(b in (0, 1) for b in v)


@given(
    st.sampled_from(["UUU", "GGG", "AAA", "CCC", "AUG", "UAG"]),
    st.sampled_from(["UUU", "GGG", "AAA", "CCC", "AUG", "UAG"]),
)
def test_hamming_symmetry(c1, c2):
    v1, v2 = codon_to_vector(c1), codon_to_vector(c2)
    assert hamming_distance(v1, v2) == hamming_distance(v2, v1)


@given(
    st.sampled_from(["UUU", "GGG", "AAA", "CCC", "AUG"]),
    st.sampled_from(["UUU", "GGG", "AAA", "CCC", "AUG"]),
)
def test_hamming_range(c1, c2):
    d = hamming_distance(codon_to_vector(c1), codon_to_vector(c2))
    assert 0 <= d <= 6


def test_all_encodings_count():
    """4! = 24 possible base-to-bit encodings."""
    encs = all_encodings()
    assert len(encs) == 24


def test_all_encodings_are_bijections():
    for enc in all_encodings():
        assert len(enc) == 4
        assert set(enc.keys()) == set(BASES)
        assert len(set(enc.values())) == 4  # all values distinct
