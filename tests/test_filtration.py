import pytest
from codon_topo.core.filtration import (
    classify_degeneracy,
    check_twofold,
    check_fourfold,
    analyze_filtration,
)
from codon_topo.core.genetic_codes import STANDARD, get_code, all_table_ids


def test_classify_degeneracy_standard():
    deg = classify_degeneracy(STANDARD)
    assert sorted(deg[2]) == sorted(
        ["Phe", "Tyr", "His", "Gln", "Asn", "Lys", "Asp", "Glu", "Cys"]
    )
    assert sorted(deg[4]) == sorted(["Val", "Pro", "Thr", "Ala", "Gly"])
    assert sorted(deg[6]) == sorted(["Leu", "Arg", "Ser"])
    assert deg[1] == ["Met", "Trp"] or set(deg[1]) == {"Met", "Trp"}


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


# Variant codes that create new 2-fold pairs that break bit-5 filtration.
# This is a substantive biological finding: when a stop codon is reassigned
# to an amino acid that was previously 1-fold (Trp) or that creates a new
# pair partner not differing at bit 5, the tautological "bit-5 filtration"
# property of the standard code is violated.
#
# Documented exceptions (NCBI gc.prt v4.6, default encoding C=00, U=01, A=10, G=11):
#   Table 32 (Balanophoraceae Plastid): UAG -> Trp creates pair (UGG, UAG)
#     differing at bit position 3, not bit 5.
TWOFOLD_FILTRATION_EXCEPTIONS = {
    32: {"Trp"},
}


@pytest.mark.parametrize("table_id", all_table_ids())
def test_twofold_universal(table_id):
    """Two-fold filtration holds across all NCBI tables, with documented exceptions.

    The bit-5 filtration is tautological for the standard code's nine 2-fold
    amino acids. Variant codes that reassign a stop codon to a previously
    1-fold amino acid (e.g. Trp) can create new 2-fold pairs whose Hamming
    separation differs from bit 5. These exceptions are recorded explicitly.
    """
    code = get_code(table_id)
    results = check_twofold(code)
    expected_failures = TWOFOLD_FILTRATION_EXCEPTIONS.get(table_id, set())
    for aa, passed, differing in results:
        if aa in expected_failures:
            assert not passed, (
                f"Table {table_id}, {aa}: expected filtration-breaking pair, "
                f"but it passed. Exception list may be stale."
            )
        else:
            assert passed, f"Table {table_id}, {aa}: differs at {differing}"


def test_analyze_filtration_returns_summary():
    report = analyze_filtration(STANDARD)
    assert report["twofold_pass"] == 9
    assert report["twofold_fail"] == 0
    assert report["fourfold_pass"] == 5
    assert report["fourfold_fail"] == 0
