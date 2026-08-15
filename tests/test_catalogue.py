"""Tests for prediction catalogue assembly."""

from codon_topo.reports.catalogue import (
    Prediction,
    build_catalogue,
    catalogue_summary,
)


def test_prediction_fields():
    cat = build_catalogue()
    assert len(cat) > 0
    p = cat[0]
    assert isinstance(p, Prediction)
    assert hasattr(p, "id")
    assert hasattr(p, "claim")
    assert hasattr(p, "workstream")
    assert hasattr(p, "status")
    assert hasattr(p, "evidence_strength")
    assert hasattr(p, "implications")


def test_catalogue_has_ws1_predictions():
    cat = build_catalogue()
    ws1 = [p for p in cat if p.workstream == "WS1"]
    assert len(ws1) >= 3  # Claims 1, 2, 3 at minimum


def test_catalogue_has_ws4_prediction():
    cat = build_catalogue()
    ws4 = [p for p in cat if p.workstream == "WS4"]
    assert len(ws4) >= 1


def test_catalogue_has_serine_invariant():
    cat = build_catalogue()
    ser = [p for p in cat if "Serine" in p.claim or "serine" in p.claim]
    assert len(ser) >= 1
    assert ser[0].evidence_strength in ("strong", "very_strong")


def test_catalogue_summary_structure():
    s = catalogue_summary()
    assert "total_predictions" in s
    assert "by_workstream" in s
    assert "by_status" in s
    assert "by_strength" in s
    assert s["total_predictions"] > 0


def test_prediction_ids_unique():
    cat = build_catalogue()
    ids = [p.id for p in cat]
    assert len(ids) == len(set(ids))
