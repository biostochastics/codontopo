import pytest
from codon_topo.analysis.null_models import (
    null_model_a, null_model_b, null_model_c,
)

def test_null_model_a_returns_results():
    """Smoke test with small n_permutations."""
    result = null_model_a(n_permutations=100, seed=42)
    assert 'n_permutations' in result
    assert 'p_value_serine_unique' in result
    assert 'p_value_bit5_uniform' in result
    assert 'p_value_fourfold_uniform' in result
    assert 0.0 <= result['p_value_serine_unique'] <= 1.0
    assert 0.0 <= result['p_value_bit5_uniform'] <= 1.0
    assert 0.0 <= result['p_value_fourfold_uniform'] <= 1.0

def test_null_model_a_deterministic():
    r1 = null_model_a(n_permutations=50, seed=123)
    r2 = null_model_a(n_permutations=50, seed=123)
    assert r1 == r2

def test_null_model_a_observed_flags():
    result = null_model_a(n_permutations=10, seed=1)
    assert result['observed_exactly_one_disconnected'] is True
    assert result['observed_twofold_all_pass'] is True
    assert result['observed_fourfold_all_pass'] is True

def test_null_model_a_accepts_code():
    from codon_topo.core.genetic_codes import get_code
    result = null_model_a(n_permutations=10, code=get_code(2), seed=42)
    assert result['n_permutations'] == 10

def test_null_model_b_returns_results():
    result = null_model_b(n_permutations=100, seed=42)
    assert 'p_value_serine_unique' in result
    assert result['exclude_stops'] is False

def test_null_model_b_exclude_stops():
    result = null_model_b(n_permutations=100, seed=42, exclude_stops=True)
    assert result['exclude_stops'] is True
    assert 'p_value_serine_unique' in result

def test_null_model_b_accepts_code():
    from codon_topo.core.genetic_codes import get_code
    result = null_model_b(n_permutations=10, code=get_code(2), seed=42)
    assert result['n_permutations'] == 10

def test_null_model_c_returns_results():
    """Tests all 24 encodings. No permutations needed."""
    result = null_model_c()
    assert 'n_encodings' in result
    assert result['n_encodings'] == 24
    assert 'twofold_invariant_count' in result
    assert 'encoding_results' in result

def test_null_model_c_has_disconnected_aas():
    result = null_model_c()
    for entry in result['encoding_results']:
        assert 'disconnected_aas' in entry
        assert 'fourfold_all_pass' in entry

def test_null_model_c_accepts_code():
    from codon_topo.core.genetic_codes import get_code
    result = null_model_c(code=get_code(2))
    assert result['n_encodings'] == 24
