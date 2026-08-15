#!/usr/bin/env python3
"""Generate manuscript_stats.json — all inline statistics for the Typst manuscript.

Every number in the manuscript is rendered from this JSON. No hardcoded values.

Usage:
    python3.11 scripts/generate_manuscript_stats.py [--output output/manuscript_stats.json]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


class _Encoder(json.JSONEncoder):
    """Handle numpy types transparently."""

    def default(self, o):  # noqa: ANN001
        import numpy as np

        if isinstance(o, np.integer):
            return int(o)
        if isinstance(o, np.floating):
            return float(o)
        if isinstance(o, np.bool_):
            return bool(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        return super().default(o)


def _f(x) -> float:
    """Coerce to native float."""
    return float(x)


def _i(x) -> int:
    """Coerce to native int."""
    return int(x)


def compute_metrics(n_samples: int = 10_000, seed: int = 135325) -> dict:
    """Cross-metric coloring optimality (Section 3.1)."""
    from codon_topo.analysis.coloring_optimality import multi_metric_sensitivity

    raw = multi_metric_sensitivity(n_samples=n_samples, seed=seed)
    metrics = {}
    for r in raw["per_metric"]:
        m = r["metric"]
        metrics[m] = {
            "observed": _f(r["observed_score"]),
            "null_mean": _f(r["null_mean"]),
            "null_std": _f(r["null_std"]),
            "z": _f(r["effect_size_z"]),
            "p": _f(r["p_value_conservative"]),
            "quantile": _f(r["quantile"]),
            "improvement_pct": _f(
                100.0 * (r["null_mean"] - r["observed_score"]) / r["null_mean"]
            ),
        }
    return metrics


def compute_rho_sweep(n_samples: int = 1_000, seed: int = 135325) -> dict:
    """Rho robustness sweep (Section 3.2)."""
    from codon_topo.analysis.coloring_optimality import (
        rho_robustness_sweep,
    )

    raw = rho_robustness_sweep(seed=seed, n_samples=n_samples)
    return {
        "per_rho": [
            {
                "rho": _f(r["rho"]),
                "p_value": _f(r["p_value"]),
                "effect_size_z": _f(r["effect_size_z"]),
                "quantile": _f(r["quantile"]),
                "observed_score": _f(r["observed_score"]),
                "null_mean": _f(r["null_mean"]),
            }
            for r in raw["per_rho"]
        ]
    }


def compute_per_table(n_samples: int = 1_000, seed: int = 135325) -> dict:
    """Per-table optimality preservation (Section 3.3)."""
    from codon_topo.analysis.coloring_optimality import per_table_optimality

    raw = per_table_optimality(n_samples=n_samples, seed=seed)
    return {
        "n_significant_bh": _i(raw["n_significant_p05_bh"]),
        "n_tables": _i(raw["n_tables"]),
        "mean_quantile": _f(raw["mean_quantile"]),
        "max_quantile": _f(raw["max_quantile"]),
    }


def compute_decomposition() -> dict:
    """Score decomposition by nucleotide position."""
    from codon_topo.analysis.coloring_optimality import score_decomposition_by_position

    raw = score_decomposition_by_position()
    return {
        "total_score": _f(raw["total_score"]),
        "position_fractions": {
            "position_1": _f(raw["position_fractions"]["pos1"]),
            "position_2": _f(raw["position_fractions"]["pos2"]),
            "position_3": _f(raw["position_fractions"]["pos3_wobble"]),
        },
    }


def compute_topology_avoidance_q6() -> dict:
    """Topology avoidance under Q6 adjacency (Section 3.4)."""
    from codon_topo.analysis.synbio_feasibility import topology_avoidance_test

    raw = topology_avoidance_test()
    rate_obs = _f(raw["rate_observed"])
    rate_poss = _f(raw["rate_possible"])
    return {
        "observed_breaks": _i(raw["observed_creates_disc"]),
        "observed_total": _i(raw["observed_total"]),
        "rate_observed": rate_obs,
        "possible_breaks": _i(raw["possible_creates_disc"]),
        "possible_total": _i(raw["possible_total"]),
        "rate_possible": rate_poss,
        "depletion_fold": _f(rate_poss / max(rate_obs, 0.001)),
        "hypergeom_p": _f(raw["hypergeom_p"]),
        "permutation_p": _f(raw["permutation_p"]),
    }


def compute_topology_avoidance_k43(seed: int = 135325) -> dict:
    """Topology avoidance under K4^3 adjacency (Section 3.4)."""
    from codon_topo.analysis.synbio_feasibility import topology_avoidance_k43

    raw = topology_avoidance_k43(seed=seed)
    rate_obs = _f(raw["rate_observed"])
    rate_poss = _f(raw["rate_possible"])
    return {
        "observed_breaks": _i(raw["observed_breaks"]),
        "observed_total": _i(raw["observed_total"]),
        "rate_observed": rate_obs,
        "possible_breaks": _i(raw["possible_breaks"]),
        "possible_total": _i(raw["possible_total"]),
        "rate_possible": rate_poss,
        "depletion_fold": _f(rate_poss / max(rate_obs, 0.001)),
        "risk_ratio": _f(raw["risk_ratio"]),
        "risk_ratio_ci_95": [_f(x) for x in raw["risk_ratio_ci_95"]],
        "hypergeom_p": _f(raw["hypergeom_p"]),
        "permutation_p": _f(raw["permutation_p"]),
    }


def compute_condlogit(max_orderings: int = 6) -> dict:
    """Conditional logit event-level model (Section 3.5)."""
    from codon_topo.analysis.evolutionary_simulation import (
        MODELS,
        build_all_choice_sets,
        clear_component_cache,
        fit_all_models,
        likelihood_ratio_test,
        observed_move_ranks,
        phys_topo_correlation,
    )

    clear_component_cache()
    all_cs = build_all_choice_sets(max_orderings_per_table=max_orderings)
    fits = fit_all_models(all_cs)
    corr = phys_topo_correlation(all_cs)

    model_fits: dict = {}
    for name, f in fits.items():
        wr = f.get("weights_raw", [])
        if hasattr(wr, "tolist"):
            wr = wr.tolist()
        model_fits[name] = {
            "ll": _f(f["log_likelihood"]),
            "log_likelihood": _f(f["log_likelihood"]),
            "k": _i(f["n_params"]),
            "n_params": _i(f["n_params"]),
            "n": _i(f["n_obs"]),
            "n_obs": _i(f["n_obs"]),
            "aic": _f(f["aic"]),
            "aicc": _f(f["aicc"]),
            "weights_raw": [_f(x) for x in wr],
            "weight_labels": f.get("weight_labels", []),
        }

    lr_tests: dict = {}
    for rmodel, fmodel, label in [
        ("M1_phys", "M3_phys_topo", "M1_vs_M3"),
        ("M2_topo", "M3_phys_topo", "M2_vs_M3"),
        ("M3_phys_topo", "M4_full", "M3_vs_M4"),
    ]:
        lrt = likelihood_ratio_test(fits[rmodel], fits[fmodel])
        lr_tests[label] = {
            "lr_statistic": _f(lrt["lr_statistic"]),
            "df": _i(lrt["df"]),
            "p_value": _f(lrt["p_value"]),
            "significant": bool(lrt["significant_p05"]),
        }

    ranks_best = observed_move_ranks(
        all_cs, fits["M3_phys_topo"], MODELS["M3_phys_topo"]
    )
    avg_pct = _f(sum(r["percentile"] for r in ranks_best) / len(ranks_best))

    return {
        "model_fits": model_fits,
        "models": model_fits,  # alias for backward compat
        "lr_tests": lr_tests,
        "phys_topo_rho": _f(corr["spearman_rho"]),
        "phys_topo_p": _f(corr["spearman_p"]),
        "avg_percentile": avg_pct,
        "n_tables": _i(len(all_cs)),
        "total_events": _i(sum(len(o[0]) for o in all_cs.values())),
    }


def compute_trna() -> dict:
    """tRNA enrichment statistics (Section 3.6)."""
    from codon_topo.analysis.trna_evidence import (
        fisher_exact_per_pairing,
        maximal_independent_set_analysis,
        trna_evidence_summary,
    )

    fisher = fisher_exact_per_pairing()
    mis = maximal_independent_set_analysis()
    summary = trna_evidence_summary()

    return {
        "all_pairings_p": _f(fisher["stouffer_p"]),
        "all_pairings_z": _f(fisher["stouffer_z"]),
        "all_pairings_n": _i(fisher["n_pairings"]),
        "independent_p": _f(fisher["stouffer_p_independent"]),
        "independent_z": _f(fisher["stouffer_z_independent"]),
        "independent_n": _i(fisher["n_independent_pairings"]),
        "mis_worst_p": _f(mis["worst_case_stouffer_p"]),
        "mis_best_p": _f(mis["best_case_stouffer_p"]),
        "mis_count": _i(mis["n_mis_size_ge2"]),
        "n_organisms": _i(
            summary["n_disconnection_organisms"] + summary["n_control_organisms"]
        ),
        "n_pairings": _i(fisher["n_pairings"]),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="output/manuscript_stats.json",
        help="Output JSON path",
    )
    parser.add_argument("--seed", type=int, default=135325, help="Random seed")
    parser.add_argument(
        "--n-metric-samples",
        type=int,
        default=10_000,
        help="Monte Carlo samples for metric sensitivity",
    )
    args = parser.parse_args()

    print("Computing metrics...", flush=True)
    metrics = compute_metrics(n_samples=args.n_metric_samples, seed=args.seed)

    print("Computing rho sweep...", flush=True)
    rho_sweep = compute_rho_sweep(seed=args.seed)

    print("Computing per-table optimality...", flush=True)
    per_table = compute_per_table(seed=args.seed)

    print("Computing score decomposition...", flush=True)
    decomposition = compute_decomposition()

    print("Computing topology avoidance (Q6)...", flush=True)
    tq6 = compute_topology_avoidance_q6()

    print("Computing topology avoidance (K4^3)...", flush=True)
    tk43 = compute_topology_avoidance_k43(seed=args.seed)

    print("Computing conditional logit...", flush=True)
    condlogit = compute_condlogit()

    print("Computing tRNA enrichment...", flush=True)
    try:
        trna = compute_trna()
    except Exception as e:
        print(f"  tRNA computation failed ({e}), using cached values", flush=True)
        trna = {
            "all_pairings_p": None,
            "all_pairings_z": None,
            "mis_worst_p": None,
            "n_organisms": None,
            "n_pairings": None,
        }

    stats = {
        "_version": "0.3.1",
        "coloring": {
            "n_samples": args.n_metric_samples,
            "seed": args.seed,
        },
        "metrics": metrics,
        "rho_sweep": rho_sweep,
        "per_table": per_table,
        "decomposition": decomposition,
        "topology_avoidance_q6": tq6,
        "topology_avoidance_k43": tk43,
        "condlogit": condlogit,
        "trna": trna,
    }

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        json.dump(stats, f, indent=2, cls=_Encoder)

    print(f"Wrote {out} ({out.stat().st_size:,} bytes)", flush=True)


if __name__ == "__main__":
    main()
