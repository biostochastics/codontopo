"""Evolutionary depth calibration: epsilon vs divergence age correlation.

Maps disconnection reconnection epsilon values to approximate evolutionary
divergence dates from the literature. Tests the hypothesis that higher
epsilon (deeper disconnection) correlates with older evolutionary origin.

Calibration data sources:
- Serine (universal): Pre-LUCA, >3.5 Gya (Koonin & Novozhilov 2009)
- Thr (yeast mito, Table 3): Saccharomycetales, ~200-300 Mya
- Leu (chlorophycean mito, Table 16): Chlorophyceae, ~500-700 Mya
- Ala (Pachysolen, Table 26): CUG clade, ~100-200 Mya
- Ser 3-component (Candida, Table 12): CUG clade, ~100-200 Mya
"""

import random
from dataclasses import dataclass
from statistics import mean

from scipy.stats import spearmanr


@dataclass(frozen=True)
class CalibrationPoint:
    """A disconnection event with approximate evolutionary dating."""

    aa: str
    table_id: int
    table_name: str
    reconnect_eps: int
    age_mya_low: float
    age_mya_high: float
    lineage: str
    notes: str = ""


CALIBRATION_POINTS: list[CalibrationPoint] = [
    CalibrationPoint(
        aa="Ser",
        table_id=1,
        table_name="Universal (all tables)",
        reconnect_eps=4,
        age_mya_low=3500.0,
        age_mya_high=4000.0,
        lineage="Pre-LUCA, all domains of life",
        notes="UCN-AGY split, min inter-block Hamming=4, present in all 27 NCBI tables",
    ),
    CalibrationPoint(
        aa="Thr",
        table_id=3,
        table_name="Yeast Mitochondrial",
        reconnect_eps=2,
        age_mya_low=200.0,
        age_mya_high=300.0,
        lineage="Saccharomycetales (budding yeasts)",
        notes="CUN block reassigned Leu->Thr, creates ACN/CUN split",
    ),
    CalibrationPoint(
        aa="Leu",
        table_id=16,
        table_name="Chlorophycean Mito",
        reconnect_eps=2,
        age_mya_low=500.0,
        age_mya_high=700.0,
        lineage="Chlorophyceae (green algae)",
        notes="UAG reassigned Stop->Leu",
    ),
    CalibrationPoint(
        aa="Leu",
        table_id=22,
        table_name="Scenedesmus obliquus Mito",
        reconnect_eps=2,
        age_mya_low=500.0,
        age_mya_high=700.0,
        lineage="Chlorophyceae (Scenedesmus)",
        notes="UAG reassigned Stop->Leu, UCA reassigned Ser->Stop",
    ),
    CalibrationPoint(
        aa="Ala",
        table_id=26,
        table_name="Pachysolen tannophilus Nuclear",
        reconnect_eps=3,
        age_mya_low=100.0,
        age_mya_high=200.0,
        lineage="CUG clade (Pachysolen)",
        notes="CUG reassigned Leu->Ala",
    ),
    CalibrationPoint(
        aa="Ser",
        table_id=12,
        table_name="Alternative Yeast Nuclear",
        reconnect_eps=3,
        age_mya_low=100.0,
        age_mya_high=200.0,
        lineage="CUG clade (Candida)",
        notes="CUG reassigned Leu->Ser, creates 3-component Ser (unique)",
    ),
]


def _age_midpoint(p: CalibrationPoint) -> float:
    return (p.age_mya_low + p.age_mya_high) / 2.0


def compute_correlation() -> dict:
    """Compute Spearman rank correlation between reconnect_eps and age midpoint."""
    eps_values = [p.reconnect_eps for p in CALIBRATION_POINTS]
    age_values = [_age_midpoint(p) for p in CALIBRATION_POINTS]

    sr = spearmanr(eps_values, age_values)
    return {
        "spearman_rho": float(sr.statistic),  # type: ignore[attr-defined]
        "spearman_p": float(sr.pvalue),  # type: ignore[attr-defined]
        "n_points": len(CALIBRATION_POINTS),
    }


def bootstrap_ci(
    n_bootstrap: int = 10_000,
    confidence_level: float = 0.95,
    seed: int | None = None,
) -> dict:
    """Bootstrap confidence interval for the Spearman correlation."""
    rng = random.Random(seed)
    n = len(CALIBRATION_POINTS)
    rhos: list[float] = []

    for _ in range(n_bootstrap):
        sample = [rng.choice(CALIBRATION_POINTS) for _ in range(n)]
        eps_vals = [p.reconnect_eps for p in sample]
        age_vals = [_age_midpoint(p) for p in sample]
        if len(set(eps_vals)) < 2 or len(set(age_vals)) < 2:
            continue
        sr = spearmanr(eps_vals, age_vals)
        rhos.append(float(sr.statistic))  # type: ignore[attr-defined]

    if not rhos:
        return {
            "ci_low": 0.0,
            "ci_high": 0.0,
            "confidence_level": confidence_level,
            "n_bootstrap": n_bootstrap,
            "n_valid": 0,
        }

    rhos.sort()
    alpha = 1 - confidence_level
    lo_idx = int(len(rhos) * alpha / 2)
    hi_idx = int(len(rhos) * (1 - alpha / 2)) - 1

    return {
        "ci_low": rhos[max(0, lo_idx)],
        "ci_high": rhos[min(len(rhos) - 1, hi_idx)],
        "confidence_level": confidence_level,
        "n_bootstrap": n_bootstrap,
        "n_valid": len(rhos),
        "mean_rho": mean(rhos),
    }


def depth_calibration_table() -> list[dict]:
    """Return calibration points as a list of flat dicts (for CSV export)."""
    return [
        {
            "aa": p.aa,
            "table_id": p.table_id,
            "table_name": p.table_name,
            "reconnect_eps": p.reconnect_eps,
            "age_mya_low": p.age_mya_low,
            "age_mya_high": p.age_mya_high,
            "age_midpoint_mya": _age_midpoint(p),
            "lineage": p.lineage,
            "notes": p.notes,
        }
        for p in CALIBRATION_POINTS
    ]
