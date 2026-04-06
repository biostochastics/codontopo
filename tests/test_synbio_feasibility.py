"""Tests for synthetic biology feasibility scoring."""

from codon_topo.analysis.synbio_feasibility import (
    score_variant_code,
    single_reassignment_landscape,
    feasibility_summary,
)
from codon_topo.core.genetic_codes import STANDARD, get_code


def test_score_standard_code():
    """Standard code should have perfect scores."""
    score = score_variant_code(STANDARD)
    assert score["twofold_intact"] is True
    assert score["fourfold_intact"] is True
    assert score["serine_disconnected"] is True
    assert score["n_disconnected_aas"] == 1  # Only Ser
    assert score["feasibility_score"] == 1.0


def test_score_yeast_mito():
    """Yeast mito (Table 3) gains Thr disconnection."""
    score = score_variant_code(get_code(3))
    assert score["serine_disconnected"] is True
    assert score["n_disconnected_aas"] >= 2  # Ser + Thr


def test_single_reassignment_landscape_returns_list():
    landscape = single_reassignment_landscape()
    assert isinstance(landscape, list)
    assert len(landscape) > 0
    entry = landscape[0]
    assert "codon" in entry
    assert "original_aa" in entry
    assert "new_aa" in entry
    assert "feasibility_score" in entry


def test_single_reassignment_landscape_excludes_noop():
    """Should not include identity reassignments (AA -> same AA)."""
    landscape = single_reassignment_landscape()
    for entry in landscape:
        assert entry["original_aa"] != entry["new_aa"]


def test_feasibility_summary_structure():
    summary = feasibility_summary()
    assert "total_variants" in summary
    assert "high_feasibility" in summary
    assert "medium_feasibility" in summary
    assert "low_feasibility" in summary
    assert "best_variants" in summary
    assert summary["total_variants"] > 0


def test_stop_reassignment_counted():
    """Reassigning a stop codon should produce a scored variant."""
    landscape = single_reassignment_landscape()
    stop_entries = [e for e in landscape if e["original_aa"] == "Stop"]
    assert len(stop_entries) > 0
