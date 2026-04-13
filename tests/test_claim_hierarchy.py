"""Tests for the claim hierarchy module.

The paper has an explicit claim hierarchy encoded as structured data.
Tests enforce consistency (no orphaned statuses, no missing justifications).
"""

from codon_topo.reports.claim_hierarchy import (
    CLAIM_HIERARCHY,
    ClaimStatus,
    claims_by_status,
    supported_claims,
    suggestive_claims,
    exploratory_claims,
    rejected_claims,
    publishable_claims,
    hierarchy_summary_table,
    abstract_ready_paragraph,
)


class TestClaimHierarchy:
    def test_all_claims_have_unique_ids(self):
        ids = [c.id for c in CLAIM_HIERARCHY]
        assert len(ids) == len(set(ids)), f"Duplicate claim IDs: {ids}"

    def test_all_claims_have_justification(self):
        for c in CLAIM_HIERARCHY:
            assert c.justification, f"{c.id}: justification is empty"
            assert len(c.justification) > 20, f"{c.id}: justification too short"

    def test_all_claims_have_statement(self):
        for c in CLAIM_HIERARCHY:
            assert c.statement, f"{c.id}: statement is empty"

    def test_status_values_are_valid(self):
        valid_statuses = set(ClaimStatus)
        for c in CLAIM_HIERARCHY:
            assert c.status in valid_statuses

    def test_supported_claims_not_empty(self):
        """We should have at least one supported finding for a paper."""
        assert len(supported_claims()) >= 1, (
            "No supported claims — cannot write a paper with nothing solid"
        )

    def test_rejected_claims_documented(self):
        """All rejected/falsified/tautological claims should have citations or data."""
        for c in rejected_claims():
            # These should have an explanation of WHY they're rejected
            assert (
                c.citation_key is not None
                or c.null_model is not None
                or "counterexample" in c.justification.lower()
                or "rejected" in c.justification.lower()
                or "tautological" in c.justification.lower()
                or "trivial" in c.justification.lower()
                or "forced" in c.justification.lower()
                or "not holomorphic" in c.justification.lower()
            ), f"{c.id}: rejection not clearly justified"

    def test_publishable_subset_correct(self):
        """Only supported, suggestive, and exploratory go into the paper."""
        pub = publishable_claims()
        for c in pub:
            assert c.status in (
                ClaimStatus.SUPPORTED,
                ClaimStatus.SUGGESTIVE,
                ClaimStatus.EXPLORATORY,
            )
        # Everything else must be filtered out
        for c in CLAIM_HIERARCHY:
            if c.status in (
                ClaimStatus.REJECTED,
                ClaimStatus.FALSIFIED,
                ClaimStatus.TAUTOLOGICAL,
            ):
                assert c not in pub

    def test_hypercube_coloring_is_supported(self):
        """The central result must be in the supported tier."""
        sup = supported_claims()
        assert any(c.id == "hypercube_coloring_optimality" for c in sup)

    def test_kras_is_falsified(self):
        """Clinical KRAS prediction must be flagged as falsified."""
        fal = claims_by_status(ClaimStatus.FALSIFIED)
        assert any(c.id == "kras_fano_clinical_prediction" for c in fal)

    def test_holomorphic_is_rejected(self):
        """Holomorphic embedding must be rejected."""
        rej = claims_by_status(ClaimStatus.REJECTED)
        assert any(c.id == "holomorphic_embedding" for c in rej)

    def test_serine_min_4_is_rejected(self):
        """The 'invariant min distance 4' claim must be rejected."""
        rej = claims_by_status(ClaimStatus.REJECTED)
        assert any(c.id == "serine_min_distance_4_invariant" for c in rej)

    def test_psl_2_7_is_rejected(self):
        rej = claims_by_status(ClaimStatus.REJECTED)
        assert any(c.id == "psl_2_7_symmetry" for c in rej)

    def test_tautological_claims_flagged(self):
        """Two-fold and four-fold filtrations are tautological by construction."""
        taut = claims_by_status(ClaimStatus.TAUTOLOGICAL)
        assert any(c.id == "two_fold_bit_5_filtration" for c in taut)
        assert any(c.id == "four_fold_prefix_filtration" for c in taut)

    def test_trna_is_suggestive_not_supported(self):
        """tRNA enrichment with p=0.033 (independent) is suggestive, not primary."""
        sug = suggestive_claims()
        assert any(c.id == "trna_enrichment_reassigned_aa" for c in sug)
        sup = supported_claims()
        assert not any(c.id == "trna_enrichment_reassigned_aa" for c in sup)

    def test_bit_position_bias_is_exploratory(self):
        """p=0.075 after de-duplication is exploratory, not supported."""
        exp = exploratory_claims()
        assert any(c.id == "bit_position_bias_weighted" for c in exp)

    def test_hierarchy_summary_renders(self):
        summary = hierarchy_summary_table()
        assert "SUPPORTED" in summary
        assert "REJECTED" in summary

    def test_abstract_paragraph_renders(self):
        p = abstract_ready_paragraph()
        assert "rigorously evaluates" in p
        assert "supported" in p
        assert "suggestive" in p or "suggestive" in p.lower()


class TestClaimEvidence:
    def test_supported_claims_have_evidence(self):
        for c in supported_claims():
            assert c.null_model is not None, f"{c.id}: no null model"
            assert c.sample_size is not None, f"{c.id}: no sample size"
            # p_value can be None for aggregate claims (e.g., 24/25 FDR)
            # but justification must explain the evidence
            assert len(c.justification) > 50, f"{c.id}: justification too short"

    def test_supported_p_values_are_significant(self):
        """Supported claims with explicit p-values must have p < 0.05."""
        for c in supported_claims():
            if c.evidence_p_value is not None:
                assert c.evidence_p_value < 0.05, (
                    f"{c.id}: p={c.evidence_p_value} is not significant"
                )

    def test_falsified_claims_have_p_value_1(self):
        """Falsified claims should document the failure (p=1.0 for KRAS)."""
        for c in claims_by_status(ClaimStatus.FALSIFIED):
            # Not a strict requirement but documents the failure
            if c.evidence_p_value is not None:
                assert c.evidence_p_value >= 0.5, (
                    f"{c.id}: p={c.evidence_p_value} not reflective of failure"
                )
