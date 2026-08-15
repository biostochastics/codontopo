import pytest
from codon_topo.core.homology import (
    connected_components,
    persistent_homology,
    disconnection_catalogue,
)
from codon_topo.core.encoding import codon_to_vector
from codon_topo.core.genetic_codes import STANDARD, get_code, all_table_ids


def test_single_codon_connected():
    v = [codon_to_vector("AUG")]
    assert connected_components(v, 1) == 1


def test_serine_standard_disconnected_eps1():
    """Claim 3: Serine uniquely disconnected at eps=1."""
    ser_codons = [c for c, aa in STANDARD.items() if aa == "Ser"]
    vectors = [codon_to_vector(c) for c in ser_codons]
    assert connected_components(vectors, 1) == 2


def test_serine_reconnects_at_eps4():
    ser_codons = [c for c, aa in STANDARD.items() if aa == "Ser"]
    vectors = [codon_to_vector(c) for c in ser_codons]
    assert connected_components(vectors, 3) == 2
    assert connected_components(vectors, 4) == 1


def test_glycine_connected():
    gly_codons = [c for c, aa in STANDARD.items() if aa == "Gly"]
    vectors = [codon_to_vector(c) for c in gly_codons]
    assert connected_components(vectors, 1) == 1


def test_persistent_homology_serine():
    ser_codons = [c for c, aa in STANDARD.items() if aa == "Ser"]
    vectors = [codon_to_vector(c) for c in ser_codons]
    ph = persistent_homology(vectors, max_eps=6)
    assert ph[1] == 2  # disconnected
    assert ph[4] == 1  # reconnected


@pytest.mark.parametrize("table_id", all_table_ids())
def test_serine_always_disconnected(table_id):
    """Universal Serine invariant: disconnected at eps<=4 in all codes."""
    code = get_code(table_id)
    ser_codons = [c for c, aa in code.items() if aa == "Ser"]
    if len(ser_codons) < 2:
        pytest.skip("Serine has <2 codons")
    vectors = [codon_to_vector(c) for c in ser_codons]
    assert connected_components(vectors, 1) > 1, (
        f"Table {table_id}: Serine connected at eps=1"
    )


def test_thr_yeast_mito_disconnected():
    """Novel finding: Thr disconnected in yeast mito, eps=2."""
    code = get_code(3)
    thr_codons = [c for c, aa in code.items() if aa == "Thr"]
    vectors = [codon_to_vector(c) for c in thr_codons]
    assert connected_components(vectors, 1) == 2
    assert connected_components(vectors, 2) == 1


def test_leu_chlorophycean_disconnected():
    """Novel finding: Leu disconnected in chlorophycean mito, eps=2."""
    code = get_code(16)
    leu_codons = [c for c, aa in code.items() if aa == "Leu"]
    vectors = [codon_to_vector(c) for c in leu_codons]
    assert connected_components(vectors, 1) == 2
    assert connected_components(vectors, 2) == 1


def test_ala_pachysolen_disconnected():
    """Novel finding: Ala disconnected in Pachysolen nuclear, eps=3."""
    code = get_code(26)
    ala_codons = [c for c, aa in code.items() if aa == "Ala"]
    vectors = [codon_to_vector(c) for c in ala_codons]
    assert connected_components(vectors, 1) == 2
    assert connected_components(vectors, 3) == 1


def test_ser_candida_three_components():
    """Novel finding: Serine has 3 components in Candida (table 12)."""
    code = get_code(12)
    ser_codons = [c for c, aa in code.items() if aa == "Ser"]
    vectors = [codon_to_vector(c) for c in ser_codons]
    assert connected_components(vectors, 1) == 3


def test_disconnection_catalogue_standard():
    cat = disconnection_catalogue(STANDARD)
    aa_names = [entry["aa"] for entry in cat]
    assert "Ser" in aa_names
    ser = [e for e in cat if e["aa"] == "Ser"][0]
    assert ser["reconnect_eps"] == 4
    assert ser["min_inter_distance"] == 4
    assert ser["n_components"] == 2
