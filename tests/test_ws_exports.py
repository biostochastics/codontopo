"""Tests for WS2-WS6 data exports."""

import csv

from codon_topo.visualization.data_export import (
    export_reassignment_db,
    export_depth_calibration,
    export_fano_predictions,
    export_synbio_landscape,
    export_catalogue,
)


def test_export_reassignment_db(tmp_path):
    out = export_reassignment_db(tmp_path / "reassignments.csv")
    assert out.exists()
    with open(out) as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    assert len(rows) > 30
    assert "table_id" in rows[0]
    assert "codon" in rows[0]
    assert "source_aa" in rows[0]
    assert "target_aa" in rows[0]
    assert "hamming_to_nearest_target" in rows[0]


def test_export_depth_calibration(tmp_path):
    out = export_depth_calibration(tmp_path / "depth.csv")
    assert out.exists()
    with open(out) as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    assert len(rows) >= 5
    assert "reconnect_eps" in rows[0]
    assert "age_midpoint_mya" in rows[0]


def test_export_fano_predictions(tmp_path):
    out = export_fano_predictions(tmp_path / "fano.csv")
    assert out.exists()
    with open(out) as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    assert len(rows) >= 4
    assert "variant" in rows[0]
    assert "fano_partner_aa" in rows[0]


def test_export_synbio_landscape_small(tmp_path):
    """Export with a small subset to keep test fast."""
    out = export_synbio_landscape(tmp_path / "synbio.csv", max_variants=20)
    assert out.exists()
    with open(out) as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    assert len(rows) <= 20
    assert "codon" in rows[0]
    assert "feasibility_score" in rows[0]


def test_export_catalogue(tmp_path):
    out = export_catalogue(tmp_path / "catalogue.csv")
    assert out.exists()
    with open(out) as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    assert len(rows) > 0
    assert "id" in rows[0]
    assert "claim" in rows[0]
    assert "evidence_strength" in rows[0]
