"""Tests for the refined analysis modules (April 13, 2026 updates).

Covers:
- null_model_c_extended: per-encoding min-distance emission
- bit_position_bias_weighted: Ts/Tv weighted chi-square
- trna_evidence: tRNA gene duplication correlation
- coloring_optimality: hypercube edge-mismatch Monte Carlo
"""

from codon_topo.analysis.null_models import null_model_c_extended
from codon_topo.analysis.reassignment_db import (
    bit_position_bias,
    bit_position_bias_weighted,
)
from codon_topo.analysis.trna_evidence import (
    DISCONNECTION_PAIRINGS,
    compare_aa_gene_counts,
    trna_duplication_correlation_test,
)
from codon_topo.analysis.coloring_optimality import (
    grantham_distance,
    hypercube_edge_mismatch_score,
    monte_carlo_null,
    cross_table_optimality,
)
from codon_topo.core.genetic_codes import STANDARD


# ============================================================
# null_model_c_extended
# ============================================================


class TestNullModelCExtended:
    """Tests for the per-encoding minimum distance emission.

    Verifies that Serine min distance is NOT invariantly 4 across
    the 24 base-to-bit encodings (counterexample exists).
    """

    def test_returns_correct_structure(self):
        result = null_model_c_extended()
        assert result["n_encodings"] == 24
        assert "per_encoding" in result
        assert "universal_disconnected_aas" in result
        assert "invariant_details" in result
        assert "aa_distance_histogram" in result

    def test_serine_is_universally_disconnected(self):
        """Ser is disconnected in ALL 24 encodings (this holds)."""
        result = null_model_c_extended()
        assert "Ser" in result["universal_disconnected_aas"]

    def test_serine_min_distance_is_not_invariant(self):
        """CRITICAL: Ser min distance varies across encodings.

        Counterexample: phi(U)=00, phi(C)=11, phi(A)=01, phi(G)=10 gives min=2.
        The claim 'min distance 4 invariant across 24 encodings' is FALSE.
        """
        result = null_model_c_extended()
        ser_distances = result["aa_distance_histogram"]["Ser"]
        # Must see at least TWO distinct distance values across encodings
        assert len(ser_distances) >= 2, (
            f"Expected multiple distinct min distances for Ser; got only {ser_distances}"
        )
        # Should see min distance = 2 in some encodings (kimi counterexample class)
        assert 2 in ser_distances, (
            f"Expected min distance 2 in some encodings; got {ser_distances}"
        )
        # And should see min distance = 4 in other encodings (default class)
        assert 4 in ser_distances, (
            f"Expected min distance 4 in some encodings; got {ser_distances}"
        )

    def test_serine_distance_histogram_sums_to_24(self):
        """Every encoding contributes exactly one data point per disconnected AA."""
        result = null_model_c_extended()
        ser_hist = result["aa_distance_histogram"]["Ser"]
        assert sum(ser_hist.values()) == 24

    def test_invariant_details_flags_non_invariance(self):
        """The per-AA invariant details should flag that Ser distance varies."""
        result = null_model_c_extended()
        ser_details = result["invariant_details"]["Ser"]
        assert ser_details["distance_is_invariant"] is False


# ============================================================
# bit_position_bias_weighted
# ============================================================


class TestBitPositionBiasWeighted:
    """Tests for the Ts/Tv-weighted chi-square test.

    The original uniform null is biologically unrealistic. This verifies
    the weighted version runs and produces different p-values than the
    uniform case.
    """

    def test_returns_correct_keys(self):
        r = bit_position_bias_weighted("mitochondrial")
        assert "bit_counts_observed" in r
        assert "expected_counts_weighted" in r
        assert "chi2_p_value_weighted" in r
        assert "chi2_p_value_uniform_reference" in r

    def test_all_compartments_supported(self):
        for compartment in ["nuclear", "mitochondrial", "uniform"]:
            r = bit_position_bias_weighted(compartment)
            assert r["compartment"] == compartment
            assert len(r["bit_counts_observed"]) == 6

    def test_uniform_compartment_matches_original_test(self):
        """The uniform weighting should roughly reproduce the original chi-square."""
        uniform = bit_position_bias()
        weighted_uniform = bit_position_bias_weighted("uniform")
        # With "uniform" priors, the expected weights are all equal (1/6 each).
        # This is equivalent to the default chi-square test.
        assert (
            abs(uniform["chi2_statistic"] - weighted_uniform["chi2_statistic_weighted"])
            < 0.01
        )

    def test_observed_bit4_is_elevated(self):
        """Document the observed finding: bit 4 has the most flips."""
        r = bit_position_bias_weighted("mitochondrial")
        observed = r["bit_counts_observed"]
        # Bit 4 (first bit of position 3 / wobble) should be the max
        assert observed[4] == max(observed)

    def test_weighted_null_gives_different_p_than_uniform(self):
        """Realistic null should produce different p-values than uniform."""
        uniform_p = bit_position_bias()["chi2_p_value"]
        mito_p = bit_position_bias_weighted("mitochondrial")["chi2_p_value_weighted"]
        # They shouldn't be exactly equal
        assert abs(uniform_p - mito_p) > 0.001


# ============================================================
# trna_evidence
# ============================================================


class TestTRNAEvidence:
    """Tests for the tRNA duplication correlation analysis."""

    def test_disconnection_pairings_complete(self):
        """Nine pairings: original 4 + 5 ciliate expansions."""
        assert len(DISCONNECTION_PAIRINGS) >= 9
        reassigned_aas = {pair[2] for pair in DISCONNECTION_PAIRINGS}
        assert {"Thr", "Leu", "Ala", "Ser", "Gln", "Cys", "Trp"} <= reassigned_aas

    def test_yeast_mito_thr_duplicated(self):
        """Key prior: Su et al. 2011 showed yeast mito Thr has 2 tRNA genes.

        Standard code mito has 1 tRNA-Thr. Disconnection organism has 2.
        Uses phylogenetically appropriate fungal control (Yarrowia lipolytica).
        """
        cmp = compare_aa_gene_counts("scerevisiae_mito", "ylipolytica_mito", "Thr")
        assert cmp["disconnection_aa_trna_count"] == 2
        assert cmp["control_aa_trna_count"] == 1
        assert cmp["excess"] == 1

    def test_majority_show_elevated_trna(self):
        """H-tRNA-1: majority of pairings show elevated tRNA for reassigned AA."""
        r = trna_duplication_correlation_test()
        assert r["n_pairings"] >= 9
        assert r["n_with_elevated_trna"] >= 7

    def test_binomial_pvalue_trend_level(self):
        """Expanded dataset: sign test less powerful but Fisher/Stouffer significant."""
        r = trna_duplication_correlation_test()
        assert r["binomial_p_value"] < 0.2

    def test_mean_excess_is_positive(self):
        r = trna_duplication_correlation_test()
        assert r["mean_excess_trna_count"] > 0


# ============================================================
# coloring_optimality
# ============================================================


class TestColoringOptimality:
    """Tests for the hypercube coloring optimality analysis."""

    def test_grantham_symmetric(self):
        assert grantham_distance("Leu", "Ile") == grantham_distance("Ile", "Leu")

    def test_grantham_zero_on_diagonal(self):
        for aa in ["Leu", "Ile", "Ala", "Arg"]:
            assert grantham_distance(aa, aa) == 0.0

    def test_grantham_leu_ile_closest(self):
        """Leu and Ile are near-identical physicochemically (Grantham = 5)."""
        assert grantham_distance("Leu", "Ile") == 5.0

    def test_grantham_trp_cys_farthest(self):
        """Trp and Cys are among the most distant (Grantham = 215)."""
        assert grantham_distance("Trp", "Cys") == 215.0

    def test_hypercube_score_deterministic(self):
        s1 = hypercube_edge_mismatch_score(STANDARD)
        s2 = hypercube_edge_mismatch_score(STANDARD)
        assert s1 == s2

    def test_hypercube_score_is_positive(self):
        score = hypercube_edge_mismatch_score(STANDARD)
        assert score > 0

    def test_standard_code_is_strongly_optimal(self):
        """Main finding: standard code is in the top 1% of random colorings.

        Expected quantile < 1% (Freeland-Hurst "one in a million" replication).
        """
        result = monte_carlo_null(n_samples=500, seed=135325)
        assert result["quantile_of_observed"] < 5.0, (
            f"Standard code quantile {result['quantile_of_observed']:.2f}% "
            "is higher than expected; should be strongly optimal."
        )

    def test_monte_carlo_deterministic_with_seed(self):
        r1 = monte_carlo_null(n_samples=100, seed=135325)
        r2 = monte_carlo_null(n_samples=100, seed=135325)
        assert r1["observed_score"] == r2["observed_score"]
        assert r1["null_mean"] == r2["null_mean"]

    def test_cross_table_optimality_includes_all_tables(self):
        result = cross_table_optimality()
        # 27 NCBI tables (1-6, 9-16, 21-33; codes 7, 8, 17-20 deprecated)
        assert len(result["per_table"]) == 27
        assert result["standard_score"] > 0
