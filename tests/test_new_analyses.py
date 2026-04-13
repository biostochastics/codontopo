"""Tests for the new analyses added in the April 2026 strengthening round."""

from codon_topo.analysis.trna_evidence import (
    fisher_exact_per_pairing,
    aa_label_permutation_test,
    NEGATIVE_CONTROL_PAIRINGS,
    validate_repertoires,
)
from codon_topo.analysis.reassignment_db import (
    bit_bias_permutation_null,
    bit_bias_deduplicated,
    nucleotide_position_bias,
)
from codon_topo.analysis.coloring_optimality import (
    local_mismatch_by_codon,
    reassignment_local_cost_test,
    per_table_optimality,
    weighted_mismatch_score,
    rho_robustness_sweep,
    score_decomposition_by_position,
)
from codon_topo.analysis.synbio_feasibility import topology_avoidance_test
from codon_topo.core.genetic_codes import STANDARD


class TestTRNAExpanded:
    def test_fisher_stouffer_returns_combined_p(self):
        r = fisher_exact_per_pairing()
        assert "stouffer_z" in r
        assert "stouffer_p" in r
        assert r["n_pairings"] >= 9

    def test_fisher_stouffer_significant(self):
        r = fisher_exact_per_pairing()
        assert r["stouffer_p"] < 0.05

    def test_permutation_test_runs(self):
        r = aa_label_permutation_test(n_permutations=100, seed=135325)
        assert r["n_pairings"] >= 9
        assert "stouffer_p" in r

    def test_negative_control_exists(self):
        assert len(NEGATIVE_CONTROL_PAIRINGS) >= 1

    def test_repertoire_validation_passes(self):
        problems = validate_repertoires()
        # Only the intentional negative control mismatch
        assert len(problems) <= 1


class TestBitBiasExpanded:
    def test_permutation_null_table_preserving(self):
        r = bit_bias_permutation_null(
            n_permutations=100, seed=42, mode="table_preserving"
        )
        assert "p_value" in r
        assert 0 <= r["p_value"] <= 1

    def test_permutation_null_codon_preserving(self):
        r = bit_bias_permutation_null(
            n_permutations=100, seed=42, mode="codon_preserving"
        )
        assert "p_value" in r

    def test_deduplicated(self):
        r = bit_bias_deduplicated()
        assert r["n_unique"] <= r["n_all"]
        assert r["n_unique"] > 0

    def test_nucleotide_position(self):
        r = nucleotide_position_bias()
        assert len(r["position_counts"]) == 3
        assert sum(r["position_counts"]) > 0


class TestLocalMismatch:
    def test_local_costs_computed(self):
        costs = local_mismatch_by_codon()
        assert len(costs) == 64
        assert all(v >= 0 for v in costs.values())

    def test_reassignment_test_runs(self):
        r = reassignment_local_cost_test()
        assert r["n_reassigned"] > 0
        assert "mann_whitney_p" in r


class TestPerTableOptimality:
    def test_runs_with_small_n(self):
        r = per_table_optimality(n_samples=50, seed=135325)
        assert r["n_tables"] >= 20
        assert "n_significant_p05_raw" in r
        assert "n_significant_p05_bh" in r


class TestWeightedEdges:
    def test_rho_zero_matches_original(self):
        from codon_topo.analysis.coloring_optimality import (
            hypercube_edge_mismatch_score,
        )

        s_original = hypercube_edge_mismatch_score(STANDARD)
        s_rho0 = weighted_mismatch_score(STANDARD, rho=0.0)
        assert s_original == s_rho0

    def test_rho_one_higher_than_zero(self):
        s0 = weighted_mismatch_score(STANDARD, rho=0.0)
        s1 = weighted_mismatch_score(STANDARD, rho=1.0)
        assert s1 > s0

    def test_sweep_runs(self):
        r = rho_robustness_sweep(rho_values=[0.0, 0.5, 1.0], n_samples=50, seed=135325)
        assert len(r["per_rho"]) == 3


class TestScoreDecomposition:
    def test_positions_sum_to_total(self):
        r = score_decomposition_by_position()
        pos_sum = sum(r["by_nucleotide_position"].values())
        assert abs(pos_sum - r["total_score"]) < 0.01

    def test_top_pairs_present(self):
        r = score_decomposition_by_position()
        assert len(r["top_aa_pairs"]) > 0


class TestTopologyAvoidance:
    def test_observed_less_than_possible(self):
        r = topology_avoidance_test()
        assert r["rate_observed"] < r["rate_possible"]

    def test_significant(self):
        r = topology_avoidance_test()
        assert r["fisher_p"] < 0.001
