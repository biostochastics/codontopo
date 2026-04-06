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
    assert 0.0 <= result['p_value_serine_unique'] <= 1.0
    assert 0.0 <= result['p_value_bit5_uniform'] <= 1.0

def test_null_model_a_deterministic():
    r1 = null_model_a(n_permutations=50, seed=123)
    r2 = null_model_a(n_permutations=50, seed=123)
    assert r1 == r2

def test_null_model_b_returns_results():
    result = null_model_b(n_permutations=100, seed=42)
    assert 'p_value_serine_unique' in result

def test_null_model_c_returns_results():
    """Tests all 24 encodings. No permutations needed."""
    result = null_model_c()
    assert 'n_encodings' in result
    assert result['n_encodings'] == 24
    assert 'twofold_invariant_count' in result
    assert 'encoding_results' in result
