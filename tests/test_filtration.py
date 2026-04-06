import pytest
from codon_topo.core.filtration import (
    classify_degeneracy, check_twofold, check_fourfold, analyze_filtration,
)
from codon_topo.core.genetic_codes import STANDARD, get_code, all_table_ids

def test_classify_degeneracy_standard():
    deg = classify_degeneracy(STANDARD)
    assert sorted(deg[2]) == sorted(['Phe','Tyr','His','Gln','Asn','Lys','Asp','Glu','Cys'])
    assert sorted(deg[4]) == sorted(['Val','Pro','Thr','Ala','Gly'])
    assert sorted(deg[6]) == sorted(['Leu','Arg','Ser'])
    assert deg[1] == ['Met', 'Trp'] or set(deg[1]) == {'Met','Trp'}

def test_twofold_standard_all_pass():
    """Claim 1: all 9 two-fold AAs differ at exactly bit 5."""
    results = check_twofold(STANDARD)
    for aa, passed, differing in results:
        assert passed, f"{aa} failed: differs at {differing}"

def test_fourfold_standard_all_pass():
    """Claim 2: all 5 four-fold AAs have uniform prefix."""
    results = check_fourfold(STANDARD)
    for aa, passed in results:
        assert passed, f"{aa} failed four-fold check"

@pytest.mark.parametrize("table_id", all_table_ids())
def test_twofold_universal(table_id):
    """Two-fold filtration is 100% invariant across all 24 codes."""
    code = get_code(table_id)
    results = check_twofold(code)
    for aa, passed, differing in results:
        assert passed, f"Table {table_id}, {aa}: differs at {differing}"

def test_analyze_filtration_returns_summary():
    report = analyze_filtration(STANDARD)
    assert report['twofold_pass'] == 9
    assert report['twofold_fail'] == 0
    assert report['fourfold_pass'] == 5
    assert report['fourfold_fail'] == 0
