from codon_topo.core.genetic_codes import STANDARD, get_code, all_table_ids
from codon_topo.core.encoding import ALL_CODONS

def test_standard_has_64_entries():
    assert len(STANDARD) == 64

def test_standard_codons_match():
    assert set(STANDARD.keys()) == set(ALL_CODONS)

def test_standard_ggu_is_gly():
    assert STANDARD['GGU'] == 'Gly'

def test_standard_uaa_is_stop():
    assert STANDARD['UAA'] == 'Stop'

def test_standard_aug_is_met():
    assert STANDARD['AUG'] == 'Met'

def test_all_table_ids():
    ids = all_table_ids()
    assert 1 in ids
    assert 2 in ids
    assert len(ids) == 24

def test_get_code_table1_is_standard():
    assert get_code(1) == STANDARD

def test_get_code_vert_mito_uga_trp():
    code = get_code(2)
    assert code['UGA'] == 'Trp'
    assert code['AGA'] == 'Stop'

def test_get_code_yeast_mito_cun_thr():
    code = get_code(3)
    for codon in ('CUU', 'CUC', 'CUA', 'CUG'):
        assert code[codon] == 'Thr'

def test_every_code_has_64_entries():
    for tid in all_table_ids():
        assert len(get_code(tid)) == 64, f"Table {tid} missing codons"

def test_invalid_table_raises():
    import pytest
    with pytest.raises(KeyError):
        get_code(999)
