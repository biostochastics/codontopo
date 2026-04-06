"""Tests for reassignment directionality statistics."""

from codon_topo.analysis.reassignment_db import (
    bit_position_bias,
    directionality_summary,
)


def test_bit_position_bias_returns_6_positions():
    result = bit_position_bias()
    assert len(result["bit_counts"]) == 6
    assert all(isinstance(c, int) for c in result["bit_counts"])
    assert sum(result["bit_counts"]) > 0


def test_bit_position_bias_has_chi2():
    result = bit_position_bias()
    assert "chi2_statistic" in result
    assert "chi2_p_value" in result
    assert 0.0 <= result["chi2_p_value"] <= 1.0


def test_directionality_summary_structure():
    result = directionality_summary()
    assert "total_events" in result
    assert "sense_to_sense" in result
    assert "sense_to_stop" in result
    assert "stop_to_sense" in result
    assert "unique_codons" in result
    assert "unique_source_aas" in result
    assert "unique_target_aas" in result
    assert result["total_events"] > 0


def test_directionality_most_reassigned_codon():
    """UGA should be among the most frequently reassigned codons."""
    result = directionality_summary()
    assert "UGA" in result["top_codons"]
