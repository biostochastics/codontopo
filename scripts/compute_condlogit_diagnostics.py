"""Compute supplementary diagnostics for FigG_condlogit_diagnostics.

Produces `output/condlogit_diagnostics.json` with:
  1. `phys_topo_scatter_sample` — subsample of (delta_phys, delta_topo)
     pairs pooled across the candidate landscape, for the joint-distribution
     panel.
  2. `model_ses` — per-model normalized-scale coefficient standard errors
     computed from the numerical Hessian at the fitted weights, plus the
     corresponding raw-scale SEs (dividing by feature stds), for the
     forest-plot panel with 95% Wald CIs.
  3. `posterior_predictive_simulated_rates` — the full array of simulated
     topology-breaking rates from the M3 posterior-predictive check, for
     the histogram panel.

Reuses fitted weights already persisted in `output/manuscript_stats.json`
so refitting is skipped; only the choice sets and Hessian/simulation
computations are done here.
"""

from __future__ import annotations

import json
import random
from pathlib import Path

import numpy as np

from codon_topo.analysis.evolutionary_simulation import (
    MODELS,
    _feature_vector,
    _global_normalization,
    _precompute_feature_bundle,
    _spec_columns,
    _total_log_likelihood_vec,
    build_all_choice_sets,
)


REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = REPO_ROOT / "output"
STATS_PATH = OUT_DIR / "manuscript_stats.json"
OUT_PATH = OUT_DIR / "condlogit_diagnostics.json"

SEED = 135325
MAX_ORDERINGS = 720
N_PP_SIMULATIONS = 10_000
SCATTER_SAMPLE_SIZE = 5_000  # subsample pooled candidate deltas for scatter


def numerical_hessian(f, x: np.ndarray, eps: float = 1e-4) -> np.ndarray:
    """Central-difference approximation of the Hessian of scalar f at x.

    H[i,j] = (f(x + e_i + e_j) - f(x + e_i - e_j)
            - f(x - e_i + e_j) + f(x - e_i - e_j)) / (4 eps^2)
    """
    n = len(x)
    H = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        for j in range(i, n):
            xpp = x.copy()
            xpp[i] += eps
            xpp[j] += eps
            xpm = x.copy()
            xpm[i] += eps
            xpm[j] -= eps
            xmp = x.copy()
            xmp[i] -= eps
            xmp[j] += eps
            xmm = x.copy()
            xmm[i] -= eps
            xmm[j] -= eps
            H[i, j] = (f(xpp) - f(xpm) - f(xmp) + f(xmm)) / (4.0 * eps * eps)
            H[j, i] = H[i, j]
    return H


def compute_model_ses(
    bundle, feat_means_full, feat_stds_full, model_fits: dict
) -> dict:
    """Compute SEs for each fitted model via numerical Hessian at optimum.

    Uses the same `_total_log_likelihood_vec` used at fit time. The
    covariance matrix on the normalized scale is Cov = (-Hessian(logL))^-1;
    SE on the raw scale is SE_normalized / feat_stds_active.
    """
    ses: dict[str, dict] = {}
    for name, fit in model_fits.items():
        if not fit.get("weight_labels"):
            continue
        spec = MODELS.get(name)
        if spec is None:
            continue
        w_norm = np.asarray(fit["weights_normalized"], dtype=np.float64)
        if w_norm.size == 0:
            continue
        cols = _spec_columns(spec)
        feat_stds_active = feat_stds_full[cols]

        def logL(w: np.ndarray) -> float:
            return _total_log_likelihood_vec(
                w, bundle, cols, feat_means_full, feat_stds_full
            )

        # -Hessian of logL = observed information matrix
        H = numerical_hessian(logL, w_norm, eps=1e-3)
        info = -H
        try:
            cov = np.linalg.inv(info)
            se_norm = np.sqrt(np.clip(np.diag(cov), 0.0, None))
        except np.linalg.LinAlgError:
            se_norm = np.full(w_norm.shape[0], float("nan"))
        se_raw = se_norm / feat_stds_active

        ses[name] = {
            "weight_labels": fit["weight_labels"],
            "weights_normalized": w_norm.tolist(),
            "weights_raw": (w_norm / feat_stds_active).tolist(),
            "weights_se_normalized": se_norm.tolist(),
            "weights_se_raw": se_raw.tolist(),
        }
    return ses


def sample_phys_topo_pairs(all_cs, n_sample: int, seed: int) -> dict:
    """Pool (delta_phys, delta_topo) across all candidate moves; subsample."""
    phys: list[float] = []
    topo: list[float] = []
    for orderings in all_cs.values():
        for cs in orderings[0]:
            for m in cs.candidates:
                phys.append(m.delta_phys)
                topo.append(float(m.delta_topo))
    n_total = len(phys)
    rng = random.Random(seed)
    if n_total > n_sample:
        idx = rng.sample(range(n_total), n_sample)
        phys_s = [phys[i] for i in idx]
        topo_s = [topo[i] for i in idx]
    else:
        phys_s = phys
        topo_s = topo
    return {
        "n_pool": n_total,
        "n_sample": len(phys_s),
        "delta_phys": phys_s,
        "delta_topo": topo_s,
    }


def posterior_predictive_with_rates(all_cs, fit_result: dict, spec) -> dict:
    """Vectorized posterior-predictive simulation matching evolutionary_simulation.py.

    Precomputes per-choice-set utilities and softmax probabilities once, then
    draws all N_PP_SIMULATIONS samples per choice set via a single vectorized
    numpy call. Result is bitwise identical (modulo RNG stream) to the
    per-move loop in `posterior_predictive_topo_rate`.
    """
    weights = np.asarray(fit_result["weights_normalized"], dtype=np.float64)
    feat_stds = np.asarray(fit_result["feat_stds"], dtype=np.float64)
    # Softmax is shift-invariant, so subtracting any constant (e.g. the
    # global mean) from every candidate's utility in a fixed choice set
    # leaves the probabilities unchanged. We therefore skip feat_means.

    cs_list = []
    for orderings in all_cs.values():
        cs_list.extend(orderings[0])
    cs_list = [cs for cs in cs_list if cs.observed_move is not None]
    n_cs = len(cs_list)

    # Observed breaking rate
    n_obs_breaking = sum(
        1 for cs in cs_list if cs.observed_move and cs.observed_move.topo_breaking
    )
    obs_rate = n_obs_breaking / max(n_cs, 1)

    # Precompute per-choice-set utilities and topo_breaking flags
    # Utilities: dot(w, x / feat_stds) for each candidate
    per_cs_probs: list[np.ndarray] = []
    per_cs_break: list[np.ndarray] = []
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        for cs in cs_list:
            # Feature matrix (n_cand, k) via _feature_vector applied per candidate
            feat_mat = np.stack([_feature_vector(m, spec) for m in cs.candidates])
            # z-score by feat_stds
            feat_mat = feat_mat / feat_stds  # broadcasting over rows
            util = feat_mat @ weights
            util -= util.max()  # numerical stability
            expu = np.exp(util)
            probs = expu / expu.sum()
            per_cs_probs.append(probs)
            per_cs_break.append(
                np.array([bool(m.topo_breaking) for m in cs.candidates], dtype=bool)
            )

    # Draw N × n_cs samples: for each choice set, N multinomial draws returning
    # the chosen candidate index. Use numpy Generator with SEED.
    rng = np.random.default_rng(SEED)
    # Cumulative probs per choice set → vectorized inverse-CDF sampling
    n_breaking_per_sim = np.zeros(N_PP_SIMULATIONS, dtype=np.int64)
    for probs, break_flags in zip(per_cs_probs, per_cs_break):
        cdf = np.cumsum(probs)
        # Draw N uniforms and searchsorted into the CDF → chosen indices
        u = rng.random(N_PP_SIMULATIONS)
        idx = np.searchsorted(cdf, u, side="left")
        # Clip to guard against numerical edge (searchsorted can return len(cdf))
        np.clip(idx, 0, len(probs) - 1, out=idx)
        n_breaking_per_sim += break_flags[idx].astype(np.int64)

    sim_rates_arr = n_breaking_per_sim / n_cs
    sim_rates = sim_rates_arr.tolist()

    sim_mean = float(sim_rates_arr.mean())
    sim_std = float(sim_rates_arr.std())
    pp_p = float((sim_rates_arr <= obs_rate).sum()) / N_PP_SIMULATIONS

    return {
        "model": fit_result["model"],
        "observed_topo_breaking_rate": obs_rate,
        "simulated_mean_rate": sim_mean,
        "simulated_std_rate": sim_std,
        "posterior_predictive_p": pp_p,
        "n_simulations": N_PP_SIMULATIONS,
        "n_choice_sets": n_cs,
        "simulated_rates": sim_rates,
    }


def main() -> None:
    print("Loading manuscript_stats.json ...")
    stats = json.loads(STATS_PATH.read_text())
    cl = stats["condlogit"]

    print(f"Building choice sets (max_orderings={MAX_ORDERINGS}) ...")
    all_cs = build_all_choice_sets(MAX_ORDERINGS)
    print(f"  Tables: {len(all_cs)}, events: {sum(len(o[0]) for o in all_cs.values())}")

    print("Precomputing feature bundle for Hessian computation ...")
    bundle = _precompute_feature_bundle(all_cs)
    feat_means_full, feat_stds_full = _global_normalization(bundle)

    print("Computing standard errors via numerical Hessian ...")
    ses = compute_model_ses(bundle, feat_means_full, feat_stds_full, cl["model_fits"])
    for name, block in ses.items():
        print(
            f"  {name}: SE(normalized) = "
            + ", ".join(f"{s:.4f}" for s in block["weights_se_normalized"])
        )

    print("Sampling pooled (delta_phys, delta_topo) pairs ...")
    scatter = sample_phys_topo_pairs(all_cs, SCATTER_SAMPLE_SIZE, SEED)
    print(f"  n_pool = {scatter['n_pool']}, n_sample = {scatter['n_sample']}")

    print("Running posterior-predictive simulation (persisting rates array) ...")
    best_name = cl["aicc_ranking"][0][0]
    best_fit = cl["model_fits"][best_name]
    pp = posterior_predictive_with_rates(all_cs, best_fit, MODELS[best_name])
    print(
        f"  model={pp['model']}  observed={pp['observed_topo_breaking_rate']:.4f}  "
        f"sim_mean={pp['simulated_mean_rate']:.4f}  pp_p={pp['posterior_predictive_p']:.3f}"
    )

    out = {
        "_doc": (
            "Supplementary diagnostics for FigG_condlogit_diagnostics: "
            "pooled candidate (delta_phys, delta_topo) sample, per-model "
            "standard errors from numerical Hessian, and posterior-predictive "
            "simulated topology-breaking rates (10 000 draws)."
        ),
        "seed": SEED,
        "phys_topo_rho": cl.get("phys_topo_rho"),
        "phys_topo_correlation": cl.get("phys_topo_correlation"),
        "phys_topo_scatter_sample": scatter,
        "model_ses": ses,
        "posterior_predictive": pp,
    }
    OUT_PATH.write_text(json.dumps(out))
    print(f"Wrote {OUT_PATH} ({OUT_PATH.stat().st_size / 1024:.1f} KB)")


if __name__ == "__main__":
    main()
