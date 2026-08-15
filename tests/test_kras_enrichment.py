"""Tests for KRAS Fano-line enrichment analysis."""

from codon_topo.analysis.cosmic_query import (
    fano_enrichment_test,
    fano_predictions_for_kras,
    ws4_gate_decision,
)


# Synthetic data where His IS enriched in G12V samples
ENRICHED_DATA = [
    # G12V samples — His mutations present
    *[
        {
            "sampleId": f"V{i}",
            "proteinChange": "G12V",
            "mutationType": "Missense_Mutation",
        }
        for i in range(50)
    ],
    *[
        {
            "sampleId": f"V{i}",
            "proteinChange": f"H{100 + i}R",
            "mutationType": "Missense_Mutation",
        }
        for i in range(30)
    ],  # 30 of 50 G12V samples have His mutations
    # G12V samples also have some non-His co-mutations
    *[
        {
            "sampleId": f"V{i}",
            "proteinChange": f"T{300 + i}A",
            "mutationType": "Missense_Mutation",
        }
        for i in range(20)
    ],  # 20 non-His co-mutations in G12V samples
    # G12D samples — fewer His mutations
    *[
        {
            "sampleId": f"D{i}",
            "proteinChange": "G12D",
            "mutationType": "Missense_Mutation",
        }
        for i in range(50)
    ],
    *[
        {
            "sampleId": f"D{i}",
            "proteinChange": f"H{200 + i}L",
            "mutationType": "Missense_Mutation",
        }
        for i in range(5)
    ],  # Only 5 of 50 G12D samples have His mutations
    # G12D samples also have non-His co-mutations
    *[
        {
            "sampleId": f"D{i}",
            "proteinChange": f"P{400 + i}S",
            "mutationType": "Missense_Mutation",
        }
        for i in range(30)
    ],  # 30 non-His co-mutations in G12D samples
]


def test_fano_enrichment_returns_result():
    result = fano_enrichment_test(ENRICHED_DATA, variant="G12V")
    assert "variant" in result
    assert "fano_predicted_aa" in result
    assert "observed_count" in result
    assert "total_co_mutations" in result
    assert "fishers_p" in result
    assert "odds_ratio" in result


def test_fano_enrichment_detects_enrichment():
    result = fano_enrichment_test(ENRICHED_DATA, variant="G12V")
    assert result["fano_predicted_aa"] == "H"
    assert result["observed_count"] == 30
    assert result["fishers_p"] < 0.01  # Should be significant


def test_fano_enrichment_g12d():
    preds = fano_predictions_for_kras()
    g12d_aa = preds["G12D"]["fano_partner_aa_1letter"]
    result = fano_enrichment_test(ENRICHED_DATA, variant="G12D")
    assert result["variant"] == "G12D"
    assert result["fano_predicted_aa"] == g12d_aa


def test_ws4_gate_decision_structure():
    result = ws4_gate_decision(ENRICHED_DATA)
    assert "pass" in result
    assert "variants_tested" in result
    assert "significant_variants" in result
    assert isinstance(result["pass"], bool)
