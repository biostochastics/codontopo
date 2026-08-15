"""Tests for cBioPortal API client and KRAS co-occurrence analysis.

All tests use mock data — no live API calls in the test suite.
"""

from codon_topo.analysis.cosmic_query import (
    KRASMutation,
    KRAS_G12_VARIANTS,
    fano_predictions_for_kras,
    parse_mutation_data,
    compute_cooccurrence,
    fetch_kras_mutations,
    CBioPortalClient,
)


# ---- Mock data fixtures ----

MOCK_MUTATIONS = [
    {"sampleId": "S1", "proteinChange": "G12V", "mutationType": "Missense_Mutation"},
    {"sampleId": "S1", "proteinChange": "H47R", "mutationType": "Missense_Mutation"},
    {"sampleId": "S2", "proteinChange": "G12V", "mutationType": "Missense_Mutation"},
    {"sampleId": "S2", "proteinChange": "T58A", "mutationType": "Missense_Mutation"},
    {"sampleId": "S3", "proteinChange": "G12D", "mutationType": "Missense_Mutation"},
    {"sampleId": "S3", "proteinChange": "H100L", "mutationType": "Missense_Mutation"},
    {"sampleId": "S4", "proteinChange": "G12D", "mutationType": "Missense_Mutation"},
    {"sampleId": "S5", "proteinChange": "G12V", "mutationType": "Missense_Mutation"},
    {"sampleId": "S5", "proteinChange": "P150S", "mutationType": "Missense_Mutation"},
    {"sampleId": "S5", "proteinChange": "H200Y", "mutationType": "Missense_Mutation"},
]


def test_kras_g12_variants_defined():
    assert "G12V" in KRAS_G12_VARIANTS
    assert "G12D" in KRAS_G12_VARIANTS
    assert "G12A" in KRAS_G12_VARIANTS
    assert len(KRAS_G12_VARIANTS) >= 4


def test_fano_predictions_for_kras():
    preds = fano_predictions_for_kras()
    assert "G12V" in preds
    # G12V: GGU->GUU, Fano partner = CAC -> His
    assert preds["G12V"]["fano_partner_codon"] == "CAC"
    assert preds["G12V"]["fano_partner_aa"] == "His"


def test_parse_mutation_data():
    mutations = parse_mutation_data(MOCK_MUTATIONS)
    assert len(mutations) == 10
    assert all(isinstance(m, KRASMutation) for m in mutations)


def test_compute_cooccurrence_mock():
    result = compute_cooccurrence(MOCK_MUTATIONS)
    assert isinstance(result, dict)
    assert "G12V" in result
    g12v = result["G12V"]
    assert "n_samples" in g12v
    assert g12v["n_samples"] == 3  # S1, S2, S5
    assert "co_mutations" in g12v


def test_cooccurrence_his_enrichment_in_mock():
    """In mock data, G12V samples S1 and S5 have His mutations."""
    result = compute_cooccurrence(MOCK_MUTATIONS)
    g12v = result["G12V"]
    his_count = sum(1 for m in g12v["co_mutations"] if m["ref_aa"] == "H")
    assert his_count == 2  # H47R (S1) and H200Y (S5)


def test_cbio_client_offline_mode():
    """Client in offline mode returns empty results, no network calls."""
    client = CBioPortalClient(offline=True)
    result = client.get_kras_mutations(study_id="test")
    assert result == []


def test_fetch_kras_mutations_with_mock():
    """fetch_kras_mutations accepts pre-loaded data."""
    result = fetch_kras_mutations(mock_data=MOCK_MUTATIONS)
    assert len(result) > 0
