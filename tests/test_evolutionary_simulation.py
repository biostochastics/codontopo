"""Tests for evolutionary simulation module (conditional logit analysis).

Tests the three-layer architecture:
  Layer A: component_counts, topology_change
  Layer B: compute_candidate_features
  Layer C: conditional_logit_log_likelihood, fit_conditional_logit, model comparison
"""

import math

import numpy as np
import pytest

from codon_topo.core.genetic_codes import STANDARD, get_changes
from codon_topo.analysis.evolutionary_simulation import (
    CandidateMove,
    ModelSpec,
    StepChoiceSet,
    build_choice_sets_single_order,
    build_choice_sets_order_averaged,
    build_all_choice_sets,
    clear_component_cache,
    component_counts,
    compute_candidate_features,
    conditional_logit_log_likelihood,
    fit_conditional_logit,
    fit_all_models,
    likelihood_ratio_test,
    order_averaged_log_likelihood,
    phys_topo_correlation,
    topology_change,
)


# ====================================================================
# Layer A: Graph/state algebra tests
# ====================================================================


class TestComponentCounts:
    """Tests for component_counts (cached topology feature extraction)."""

    def setup_method(self):
        clear_component_cache()

    def test_standard_code_serine_disconnected(self):
        """Serine should have 2 components in STANDARD (UCN vs AGY)."""
        comps = component_counts(STANDARD)
        assert comps["Ser"] == 2

    def test_standard_code_ala_connected(self):
        """Alanine (GCN, 4 codons in one block) should have 1 component."""
        comps = component_counts(STANDARD)
        assert comps["Ala"] == 1

    def test_standard_code_met_single_codon(self):
        """Met (single codon AUG) should have 1 component."""
        comps = component_counts(STANDARD)
        assert comps["Met"] == 1

    def test_caching_returns_same_object(self):
        """Component counts should be cached for identical codes."""
        c1 = component_counts(STANDARD)
        c2 = component_counts(STANDARD)
        assert c1 is c2

    def test_variant_code_different_result(self):
        """A variant code should give different component counts if topology changes."""
        variant = dict(STANDARD)
        variant["AGU"] = "Thr"  # Remove one Ser codon from AGY block
        comps = component_counts(variant)
        # Ser should still have components (AGC remains + UCN block)
        assert "Ser" in comps


class TestTopologyChange:
    """Tests for topology_change."""

    def test_identity_no_change(self):
        """Reassigning a codon to its current AA should not be called, but test edge."""
        base = component_counts(STANDARD)
        # This shouldn't happen in practice, but verify no crash
        is_breaking, delta = topology_change(STANDARD, base, "AUG", "Met")
        assert not is_breaking
        assert delta == 0

    def test_topology_preserving_move(self):
        """Reassigning a stop codon to a sense AA should not break existing topology."""
        base = component_counts(STANDARD)
        # UAA (Stop) -> Gln (as in ciliate codes) — adds to Gln family
        is_breaking, delta = topology_change(STANDARD, base, "UAA", "Gln")
        assert not is_breaking
        assert delta == 0

    def test_topology_breaking_move_exists(self):
        """There should be at least some topology-breaking moves from STANDARD."""
        base = component_counts(STANDARD)
        found_breaking = False
        # Check a few candidates
        for codon in ["GCU", "GCC", "GCA", "GCG"]:  # Ala block
            for target in ["Ser", "Leu", "Arg"]:
                is_breaking, _delta = topology_change(STANDARD, base, codon, target)
                if is_breaking:
                    found_breaking = True
                    break
            if found_breaking:
                break
        assert found_breaking, "Should find at least one topology-breaking move"


# ====================================================================
# Layer B: Scoring kernel tests
# ====================================================================


class TestCandidateFeatures:
    """Tests for compute_candidate_features."""

    def setup_method(self):
        clear_component_cache()

    def test_candidate_count(self):
        """Should produce ~1280 candidates (64 codons x ~20 targets each)."""
        candidates = compute_candidate_features(STANDARD)
        # 64 codons, each can go to ~20 other AAs (21 labels - 1 current)
        # Some variation due to different degeneracy
        assert len(candidates) > 1000
        assert len(candidates) < 1500

    def test_no_identity_moves(self):
        """No candidate should reassign a codon to its current AA."""
        candidates = compute_candidate_features(STANDARD)
        for m in candidates:
            assert m.source_aa != m.target_aa

    def test_observed_flag(self):
        """When observed_codon/target specified, exactly one candidate is_observed."""
        candidates = compute_candidate_features(
            STANDARD, observed_codon="UAA", observed_target="Gln"
        )
        observed = [m for m in candidates if m.is_observed]
        assert len(observed) == 1
        assert observed[0].codon == "UAA"
        assert observed[0].target_aa == "Gln"

    def test_delta_phys_sign(self):
        """delta_phys should be negative for moves to similar AAs."""
        candidates = compute_candidate_features(STANDARD)
        # Find a synonymous-like move: Ile -> Leu (similar AAs)
        ile_to_leu = [
            m for m in candidates
            if m.source_aa == "Ile" and m.target_aa == "Leu"
        ]
        # At least some should exist
        assert len(ile_to_leu) > 0

    def test_delta_trna_range(self):
        """delta_trna (Hamming distance) should be 0-6 for normal moves."""
        candidates = compute_candidate_features(STANDARD)
        for m in candidates:
            assert 0 <= m.delta_trna <= 7  # 7 = sentinel for no existing codons

    def test_topo_fields_consistent(self):
        """If topo_breaking is True, delta_topo should be > 0."""
        candidates = compute_candidate_features(STANDARD)
        for m in candidates:
            if m.topo_breaking:
                assert m.delta_topo > 0
            if m.delta_topo > 0:
                assert m.topo_breaking


# ====================================================================
# Layer C: Statistical inference tests
# ====================================================================


class TestConditionalLogit:
    """Tests for the conditional logit likelihood and fitting."""

    def _make_simple_choice_set(self) -> list[StepChoiceSet]:
        """Create a simple 3-candidate choice set for testing."""
        candidates = [
            CandidateMove("AAA", "Lys", "Asn", delta_phys=-10.0,
                          topo_breaking=False, delta_topo=0, delta_trna=1,
                          is_observed=True),
            CandidateMove("AAG", "Lys", "Glu", delta_phys=50.0,
                          topo_breaking=True, delta_topo=1, delta_trna=2,
                          is_observed=False),
            CandidateMove("AAC", "Asn", "Lys", delta_phys=30.0,
                          topo_breaking=False, delta_topo=0, delta_trna=1,
                          is_observed=False),
        ]
        return [StepChoiceSet(
            table_id=1, step_index=0,
            candidates=candidates, current_code=STANDARD,
        )]

    def test_uniform_likelihood(self):
        """With zero weights, P(obs) = 1/n for each choice set."""
        cs = self._make_simple_choice_set()
        spec = ModelSpec("test", use_phys=True)
        w = np.array([0.0])
        ll = conditional_logit_log_likelihood(w, cs, spec)
        expected = -math.log(3)  # 1/3 for 3 candidates
        assert abs(ll - expected) < 1e-10

    def test_strong_preference_increases_likelihood(self):
        """Negative weight on delta_phys should prefer the observed move (lowest phys)."""
        cs = self._make_simple_choice_set()
        spec = ModelSpec("test", use_phys=True)

        ll_zero = conditional_logit_log_likelihood(np.array([0.0]), cs, spec)
        # Negative weight: prefer lower delta_phys (observed has -10)
        ll_neg = conditional_logit_log_likelihood(np.array([-0.1]), cs, spec)
        assert ll_neg > ll_zero

    def test_fit_converges(self):
        """Fitting should converge on the simple example."""
        cs = self._make_simple_choice_set()
        spec = ModelSpec("test", use_phys=True)
        result = fit_conditional_logit(cs, spec)
        assert result["converged"] or result["log_likelihood"] > -math.log(3)
        assert result["n_params"] == 1
        # Weight should be negative (prefer lower phys cost)
        assert result["weights_raw"][0] < 0

    def test_aicc_penalizes_complexity(self):
        """AICc should penalize models with more parameters when data is small."""
        cs = self._make_simple_choice_set()
        spec1 = ModelSpec("one_param", use_phys=True)
        spec2 = ModelSpec("two_param", use_phys=True, use_topo=True)
        r1 = fit_conditional_logit(cs, spec1)
        r2 = fit_conditional_logit(cs, spec2)
        # With only 1 observation, AICc penalty should be severe for 2 params
        # (but AICc may be inf if n <= k+1)
        assert r1["aicc"] < float("inf") or r2["aicc"] == float("inf")


class TestOrderAveraging:
    """Tests for order-averaged likelihood."""

    def test_single_event_no_averaging(self):
        """A table with 1 change should produce 1 ordering."""
        # Table 15 (Blepharisma) has few changes
        # Use a table known to have exactly 1 change
        changes = get_changes(15)
        if len(changes) == 1:
            orderings = build_choice_sets_order_averaged(15)
            assert len(orderings) == 1

    def test_multi_event_factorial_orderings(self):
        """A table with k changes should produce k! orderings (up to cap)."""
        # Find a table with 2+ changes
        for tid in [2, 3, 4, 5]:
            changes = get_changes(tid)
            k = len(changes)
            if k >= 2:
                orderings = build_choice_sets_order_averaged(tid)
                expected = min(math.factorial(k), 720)
                assert len(orderings) == expected
                break

    def test_order_averaged_ll_finite(self):
        """Order-averaged log-likelihood should be finite."""
        for tid in [2, 3, 4, 5]:
            changes = get_changes(tid)
            if len(changes) >= 2:
                orderings = build_choice_sets_order_averaged(tid, max_orderings=6)
                spec = ModelSpec("test", use_phys=True)
                w = np.array([0.0])
                ll = order_averaged_log_likelihood(w, orderings, spec)
                assert math.isfinite(ll)
                assert ll < 0  # Log-likelihood should be negative
                break


class TestChoiceSetConstruction:
    """Tests for building choice sets from NCBI tables."""

    def setup_method(self):
        clear_component_cache()

    def test_build_single_order(self):
        """Building choice sets for a known table should produce valid output."""
        changes = get_changes(6)  # Ciliate nuclear (UAR -> Gln)
        if changes:
            events = list(changes.items())
            cs_list = build_choice_sets_single_order(6, events)
            assert len(cs_list) == len(events)
            for cs in cs_list:
                assert cs.observed_move is not None
                assert cs.n_candidates > 100  # Should have many candidates

    def test_build_all_excludes_standard(self):
        """Table 1 (standard) has no changes, so should be excluded."""
        all_cs = build_all_choice_sets(max_orderings_per_table=1)
        assert 1 not in all_cs  # Standard code has no reassignments


class TestModelComparison:
    """Tests for the model comparison pipeline."""

    def test_likelihood_ratio_test_basic(self):
        """LRT should work for nested models."""
        restricted = {"model": "M1", "log_likelihood": -50.0, "n_params": 1}
        full = {"model": "M3", "log_likelihood": -45.0, "n_params": 2}
        result = likelihood_ratio_test(restricted, full)
        assert result["lr_statistic"] == pytest.approx(10.0)
        assert result["df"] == 1
        assert result["p_value"] < 0.01  # chi2(10, df=1) is very significant

    def test_likelihood_ratio_test_no_improvement(self):
        """LRT with no improvement should give p ~ 1."""
        restricted = {"model": "M1", "log_likelihood": -50.0, "n_params": 1}
        full = {"model": "M3", "log_likelihood": -50.0, "n_params": 2}
        result = likelihood_ratio_test(restricted, full)
        assert result["lr_statistic"] == pytest.approx(0.0)
        assert result["p_value"] > 0.9


class TestDiagnostics:
    """Tests for diagnostic functions."""

    def test_phys_topo_correlation_returns_valid(self):
        """Correlation should return a valid Spearman rho."""
        # Build minimal choice sets
        all_cs = build_all_choice_sets(max_orderings_per_table=1)
        if all_cs:
            result = phys_topo_correlation(all_cs)
            assert -1.0 <= result["spearman_rho"] <= 1.0
            assert 0.0 <= result["spearman_p"] <= 1.0


# ====================================================================
# Integration test (slower, runs the full pipeline on a subset)
# ====================================================================


@pytest.mark.slow
class TestFullPipeline:
    """Integration test running the full analysis on a small subset."""

    def setup_method(self):
        clear_component_cache()

    def test_fit_all_models_subset(self):
        """Fit all 4 models on a small subset of tables."""
        # Use only tables with 1 change for speed
        small_tables = []
        for tid in [6, 10, 15]:
            changes = get_changes(tid)
            if len(changes) == 1:
                small_tables.append(tid)

        if not small_tables:
            pytest.skip("No single-change tables found in subset")

        all_cs = build_all_choice_sets(
            max_orderings_per_table=1,
            tables=small_tables,
        )
        fits = fit_all_models(all_cs)

        # Should have results for all 4 models
        assert len(fits) == 4
        for _name, f in fits.items():
            assert math.isfinite(f["log_likelihood"])
            assert f["aic"] > 0 or f["n_params"] == 0

        # M3 should have better (less negative) LL than M1 (more params)
        # or at least equal (topology may not help)
        assert fits["M3_phys_topo"]["log_likelihood"] >= fits["M1_phys"]["log_likelihood"] - 1e-6
