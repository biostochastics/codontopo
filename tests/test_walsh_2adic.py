"""Tests for the Walsh-Hadamard / 2-adic spectral probe."""

import numpy as np

from codon_topo.analysis.walsh_2adic import (
    WOBBLE_CAVEAT,
    codon_to_index,
    encoding_invariance_test,
    label_spectrum_null_test,
    label_spectrum_S,
    spectral_depth,
    spectral_depth_null_test,
    synonymous_blocks,
    v2,
    walsh_hadamard,
    wobble_free_walsh_mask,
)


# ---------------------------------------------------------------------------
# codon_to_index: respect project encoding
# ---------------------------------------------------------------------------
def test_codon_to_index_ccc_is_zero():
    """C = (0,0): CCC -> 000000 -> 0."""
    assert codon_to_index("CCC") == 0


def test_codon_to_index_ggg_is_63():
    """G = (1,1): GGG -> 111111 -> 63."""
    assert codon_to_index("GGG") == 63


def test_codon_to_index_aug_is_39():
    """A=10, U=01, G=11: AUG -> 100111 -> 39."""
    assert codon_to_index("AUG") == 39


def test_codon_to_index_all_64_unique():
    from codon_topo.core.encoding import ALL_CODONS

    idxs = [codon_to_index(c) for c in ALL_CODONS]
    assert sorted(idxs) == list(range(64))


# ---------------------------------------------------------------------------
# 2-adic valuation
# ---------------------------------------------------------------------------
def test_v2_powers_of_two():
    for k in range(10):
        assert v2(1 << k) == k


def test_v2_odd_is_zero():
    for n in [1, 3, 5, 7, 9, 11, 13, 15, 17]:
        assert v2(n) == 0


def test_v2_zero_is_none():
    assert v2(0) is None


def test_v2_negative_matches_positive():
    for n in [12, 24, 48, 96]:
        assert v2(n) == v2(-n)


# ---------------------------------------------------------------------------
# Walsh-Hadamard transform: structural properties
# ---------------------------------------------------------------------------
def test_walsh_hadamard_inverse():
    """W is self-inverse up to scale 64."""
    rng = np.random.default_rng(0)
    f = rng.standard_normal(64)
    g = walsh_hadamard(walsh_hadamard(f))
    np.testing.assert_allclose(g, 64.0 * f, atol=1e-9)


def test_walsh_hadamard_delta_at_origin():
    """W applied to delta_0 gives the constant-1 vector."""
    f = np.zeros(64)
    f[0] = 1.0
    out = walsh_hadamard(f)
    np.testing.assert_array_equal(out, np.ones(64))


def test_walsh_hadamard_constant_input():
    """W applied to all-ones gives 64 * delta_0."""
    f = np.ones(64)
    out = walsh_hadamard(f)
    expected = np.zeros(64)
    expected[0] = 64.0
    np.testing.assert_array_equal(out, expected)


# ---------------------------------------------------------------------------
# Synonymous blocks: respect the standard code
# ---------------------------------------------------------------------------
def test_synonymous_blocks_partition_64():
    blocks = synonymous_blocks()
    flat = [i for b in blocks for i in b]
    assert sorted(flat) == list(range(64))


def test_synonymous_blocks_methionine_singleton():
    """Met has exactly one codon under the standard code."""
    from codon_topo.core.genetic_codes import STANDARD
    from codon_topo.core.encoding import ALL_CODONS

    blocks = synonymous_blocks()
    sizes = sorted(len(b) for b in blocks)
    expected = sorted(
        [
            v
            for v in __import__("collections")
            .Counter(STANDARD.get(c, "Stop") for c in ALL_CODONS)
            .values()
        ]
    )
    assert sizes == expected


# ---------------------------------------------------------------------------
# spectral_depth: encoding-invariant
# ---------------------------------------------------------------------------
def test_spectral_depth_encoding_invariant():
    """Under the standard code, spectral depth must be identical across
    all 24 base-to-bit encodings."""
    from codon_topo.core.encoding import all_encodings

    depths = {spectral_depth(synonymous_blocks(encoding=e)) for e in all_encodings()}
    assert len(depths) == 1
    # Numerical value with our encoding (C=00,U=01,A=10,G=11): 544
    assert depths == {544}


def test_spectral_depth_null_test_returns_significant_z():
    """Standard code should be many SDs shallower than the block-size null."""
    res = spectral_depth_null_test(n_null=300, seed=42)
    assert res["observed_depth"] == 544
    assert res["z_score"] < -10  # Paul reports z ~ -18; n=300 still strongly negative


def test_encoding_invariance_test_marks_invariant_and_all_negative():
    res = encoding_invariance_test(n_null=200, seed=42)
    assert res["invariant"] is True
    assert res["observed_depths_unique"] == [544]
    assert res["all_z_negative"] is True


# ---------------------------------------------------------------------------
# Label spectrum S: observed value, encoding invariance, null behaviour
# ---------------------------------------------------------------------------
def test_wobble_free_walsh_mask_count():
    """16 of 64 Walsh frequencies are wobble-free (no support on bits 4,5)."""
    wf = wobble_free_walsh_mask(64)
    assert wf.sum() == 16
    assert wf[0]  # y=0 (DC) is always wobble-free
    # All wobble-free frequencies satisfy (y & 3) == 0
    assert all((y & 3) == 0 for y in range(64) if wf[y])


def test_label_spectrum_S_observed_value():
    """Standard code S = 0.7514309 under the project encoding."""
    s = label_spectrum_S()
    assert abs(s - 0.7514309076042518) < 1e-9


def test_label_spectrum_S_encoding_invariant():
    """S is mathematically invariant across all 24 base-to-bit bijections."""
    from codon_topo.core.encoding import all_encodings

    vals = {round(label_spectrum_S(encoding=enc), 10) for enc in all_encodings()}
    assert len(vals) == 1


def test_label_spectrum_null_extreme_z():
    """Under the label-permutation null, S is highly anomalous (z >> 0)."""
    res = label_spectrum_null_test(n_null=500, seed=42)
    assert res["observed_S"] > 0.75
    assert res["z_score"] > 20
    assert res["fraction_null_at_or_above_observed"] == 0.0


def test_label_spectrum_null_is_distinct_from_block_geometry():
    """The label-permutation null produces genuine variance in S, in contrast
    to the wobble-box-preserving null which leaves spectral_depth constant.
    This documents that the two tests probe orthogonal aspects of the code."""
    res = label_spectrum_null_test(n_null=500, seed=42)
    assert res["null_std"] > 0  # genuine variance, unlike spectral_depth's null
    # ~15 distinct values is the known coset-collapse property of this statistic
    assert res["n_distinct_null_S_values"] <= 30


# ---------------------------------------------------------------------------
# Wobble caveat is preserved verbatim
# ---------------------------------------------------------------------------
def test_wobble_caveat_present():
    assert "wobble" in WOBBLE_CAVEAT.lower()
    assert "label-permutation null" in WOBBLE_CAVEAT
    # Caveat is also embedded in the null-test output
    res = spectral_depth_null_test(n_null=100, seed=0)
    assert res["wobble_caveat"] == WOBBLE_CAVEAT
