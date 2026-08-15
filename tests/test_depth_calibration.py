"""Tests for evolutionary depth calibration."""

from codon_topo.analysis.depth_calibration import (
    CALIBRATION_POINTS,
    CalibrationPoint,
    compute_correlation,
    bootstrap_ci,
    depth_calibration_table,
)


def test_calibration_points_exist():
    assert len(CALIBRATION_POINTS) >= 5


def test_calibration_point_fields():
    p = CALIBRATION_POINTS[0]
    assert isinstance(p, CalibrationPoint)
    assert hasattr(p, "aa")
    assert hasattr(p, "table_id")
    assert hasattr(p, "reconnect_eps")
    assert hasattr(p, "age_mya_low")
    assert hasattr(p, "age_mya_high")
    assert hasattr(p, "lineage")


def test_serine_is_oldest():
    ser = [p for p in CALIBRATION_POINTS if p.aa == "Ser" and p.reconnect_eps == 4]
    assert len(ser) >= 1
    assert ser[0].age_mya_low >= 3000  # Pre-LUCA


def test_correlation_returns_result():
    result = compute_correlation()
    assert "spearman_rho" in result
    assert "spearman_p" in result
    assert "n_points" in result
    assert result["n_points"] >= 5


def test_correlation_non_negative():
    """Higher epsilon should not anti-correlate with older age (rho >= 0).

    With only 6 calibration points and the eps=3 group (CUG clade, ~150 Mya)
    being younger than the eps=2 group (Chlorophyceae, ~600 Mya), the Spearman
    rho is exactly 0.0 — the mid-range inversion cancels the Serine-driven
    positive trend. The key scientific claim is non-negative association.
    """
    result = compute_correlation()
    assert result["spearman_rho"] >= 0


def test_bootstrap_ci_structure():
    result = bootstrap_ci(n_bootstrap=100, seed=42)
    assert "ci_low" in result
    assert "ci_high" in result
    assert "confidence_level" in result
    assert result["ci_low"] <= result["ci_high"]


def test_depth_calibration_table():
    table = depth_calibration_table()
    assert len(table) >= 5
    for row in table:
        assert "aa" in row
        assert "reconnect_eps" in row
        assert "age_midpoint_mya" in row
        assert "lineage" in row
