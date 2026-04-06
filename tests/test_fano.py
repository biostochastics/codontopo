from codon_topo.core.fano import is_fano_line, fano_partner, all_single_bit_fano_lines
from codon_topo.core.encoding import codon_to_vector

def test_ggu_guu_cac_is_fano_line():
    """KRAS: GGU XOR GUU XOR CAC = 0."""
    assert is_fano_line('GGU', 'GUU', 'CAC')

def test_fano_partner_ggu_guu():
    assert fano_partner('GGU', 'GUU') == 'CAC'

def test_non_fano_line():
    assert not is_fano_line('GGU', 'GUU', 'AAA')

def test_fano_line_symmetric():
    assert is_fano_line('GGU', 'GUU', 'CAC')
    assert is_fano_line('GUU', 'CAC', 'GGU')
    assert is_fano_line('CAC', 'GGU', 'GUU')

def test_all_single_bit_fano_from_ggu():
    """Every single-bit mutation from GGU produces a Fano line."""
    lines = all_single_bit_fano_lines('GGU')
    assert len(lines) == 6  # 6 bit positions
    # One of them should be GGU->GUU->CAC
    codons_and_partners = [(l['mutant_codon'], l['fano_partner']) for l in lines]
    assert ('GUU', 'CAC') in codons_and_partners
