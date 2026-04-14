"""Evolutionary simulation: conditional logit model of codon reassignment.

Tests whether topology avoidance is an independent evolutionary constraint
beyond physicochemical error-minimization, using a discrete-choice framework
(conditional logit) on the 27 observed natural reassignment events.

Architecture (three-layer design, cf. hybrid stack recommendation):
  Layer A — Graph/state algebra: component-count features via homology.py
  Layer B — Scoring kernels: per-move feature vectors (delta_phys, delta_topo, delta_trna)
  Layer C — Statistical inference: conditional logit MLE, AICc model comparison

The conditional logit framing treats each observed reassignment as a choice
among ~1,280 candidate single-codon reassignments. This converts 27 data
points into 27 choices with thousands of contrasted alternatives, dramatically
improving statistical power over binary (break/no-break) analyses.

For tables with k>1 changes, the likelihood is marginalized over all k!
orderings of events (order-averaging), since temporal ordering is unknown.

Reference: GPT-5.2-pro consultation (2026-04-13) recommended SSWM origin-
fixation dynamics with conditional logit as the retrospective statistical test.

Models compared:
  M1: physicochemistry only (w_topo=0, w_trna=0)
  M2: topology only (w_phys=0, w_trna=0)
  M3: physicochemistry + topology (w_trna=0)
  M4: physicochemistry + topology + tRNA complexity (all free)
"""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass, field
from itertools import permutations
from typing import Sequence

import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2

from codon_topo.core.encoding import ALL_CODONS, codon_to_vector, hamming_distance
from codon_topo.core.genetic_codes import STANDARD, all_table_ids, get_changes
from codon_topo.core.homology import connected_components

# ====================================================================
# Layer A: Graph/state algebra — fast topology features with caching
# ====================================================================

# Cache: frozenset of (codon, aa) items -> {aa: n_components}
_COMPONENT_CACHE: dict[frozenset, dict[str, int]] = {}


def _code_to_key(code: dict[str, str]) -> frozenset:
    """Convert code dict to hashable key for caching."""
    return frozenset(code.items())


def component_counts(code: dict[str, str]) -> dict[str, int]:
    """Count connected components per AA at epsilon=1 (Hamming distance).

    Uses aggressive caching since many intermediate codes are revisited
    during order-averaging over k! permutations.

    Returns dict mapping AA -> number of connected components.
    """
    key = _code_to_key(code)
    if key in _COMPONENT_CACHE:
        return _COMPONENT_CACHE[key]

    aa_codons: dict[str, list[str]] = defaultdict(list)
    for codon, aa in code.items():
        if aa != "Stop":
            aa_codons[aa].append(codon)

    result: dict[str, int] = {}
    for aa, codons in aa_codons.items():
        if len(codons) < 2:
            result[aa] = 1 if codons else 0
            continue
        vectors = [codon_to_vector(c) for c in codons]
        result[aa] = connected_components(vectors, 1)

    _COMPONENT_CACHE[key] = result
    return result


def clear_component_cache() -> None:
    """Clear the component count cache (for testing/memory management)."""
    _COMPONENT_CACHE.clear()


def topology_change(
    base_code: dict[str, str],
    base_components: dict[str, int],
    codon: str,
    new_aa: str,
) -> tuple[bool, int]:
    """Compute topology change from a single-codon reassignment.

    Returns:
        (is_breaking, delta_components):
            is_breaking: True if any AA gains a new connected component
            delta_components: total increase in components across all AAs
                (graded measure; 0 = topology-preserving, >0 = breaking)
    """
    variant = dict(base_code)
    variant[codon] = new_aa
    var_components = component_counts(variant)

    delta = 0
    is_breaking = False
    # Check all AAs that exist in either code
    all_aas = set(base_components.keys()) | set(var_components.keys())
    for aa in all_aas:
        diff = var_components.get(aa, 0) - base_components.get(aa, 0)
        if diff > 0:
            is_breaking = True
            delta += diff

    return is_breaking, delta


# ====================================================================
# Layer B: Scoring kernels — per-move feature computation
# ====================================================================


@dataclass(frozen=True)
class CandidateMove:
    """A single candidate reassignment move from a given code state."""

    codon: str
    source_aa: str
    target_aa: str
    delta_phys: float       # Change in local Grantham mismatch
    topo_breaking: bool     # Binary: creates new disconnection?
    delta_topo: int         # Graded: total component count increase
    delta_trna: int         # Hamming distance to nearest target-AA codon
    is_observed: bool       # Was this the actually observed move?


@dataclass
class StepChoiceSet:
    """The choice set at one step of a reassignment walk.

    Contains all candidate moves from the current code state,
    with one marked as the observed choice.
    """

    table_id: int
    step_index: int
    candidates: list[CandidateMove]
    current_code: dict[str, str] = field(repr=False)

    @property
    def observed_move(self) -> CandidateMove | None:
        for m in self.candidates:
            if m.is_observed:
                return m
        return None

    @property
    def n_candidates(self) -> int:
        return len(self.candidates)


def _build_hamming1_edges() -> list[tuple[str, str]]:
    """Build list of all Hamming-1 codon pairs (192 edges of Q6)."""
    from itertools import combinations

    codon_vecs = [(c, codon_to_vector(c)) for c in ALL_CODONS]
    return [
        (c1, c2)
        for (c1, v1), (c2, v2) in combinations(codon_vecs, 2)
        if hamming_distance(v1, v2) == 1
    ]


# Precomputed edges (module-level lazy init)
_EDGES: list[tuple[str, str]] | None = None


def _get_edges() -> list[tuple[str, str]]:
    global _EDGES
    if _EDGES is None:
        _EDGES = _build_hamming1_edges()
    return _EDGES


# Precomputed per-codon edge index for fast local cost
_CODON_EDGES: dict[str, list[tuple[str, str]]] | None = None


def _get_codon_edges() -> dict[str, list[tuple[str, str]]]:
    """Index: codon -> list of edges incident to that codon."""
    global _CODON_EDGES
    if _CODON_EDGES is None:
        edges = _get_edges()
        idx: dict[str, list[tuple[str, str]]] = defaultdict(list)
        for c1, c2 in edges:
            idx[c1].append((c1, c2))
            idx[c2].append((c1, c2))
        _CODON_EDGES = dict(idx)
    return _CODON_EDGES


def _nearest_hamming_to_aa(codon: str, target_aa: str, code: dict[str, str]) -> int:
    """Hamming distance from codon to nearest codon already encoding target_aa.

    Uses the given code (not necessarily STANDARD) as reference.
    Returns 0 if the codon itself encodes target_aa, 7 if no codons encode it.
    """
    target_codons = [c for c, aa in code.items() if aa == target_aa and c != codon]
    if not target_codons:
        return 7  # Sentinel: no existing codons for this AA
    v = codon_to_vector(codon)
    return min(hamming_distance(v, codon_to_vector(tc)) for tc in target_codons)


def compute_candidate_features(
    code: dict[str, str],
    observed_codon: str | None = None,
    observed_target: str | None = None,
) -> list[CandidateMove]:
    """Compute features for all single-codon reassignment candidates.

    For each of the 64 codons x ~20 possible target AAs (excluding identity),
    compute delta_phys, delta_topo, and delta_trna.

    Args:
        code: Current genetic code state.
        observed_codon: The codon that was actually reassigned (if known).
        observed_target: The AA it was reassigned to (if known).

    Returns:
        List of CandidateMove with features populated.
    """
    from codon_topo.analysis.coloring_optimality import grantham_distance

    codon_edges = _get_codon_edges()
    base_components = component_counts(code)
    all_aas = sorted(set(code.values()))

    candidates: list[CandidateMove] = []
    for codon in ALL_CODONS:
        original_aa = code[codon]
        incident_edges = codon_edges.get(codon, [])

        # Pre-compute base local cost for this codon
        base_local = 0.0
        for c1, c2 in incident_edges:
            other = c2 if c1 == codon else c1
            other_aa = code[other]
            if original_aa != other_aa:
                base_local += grantham_distance(original_aa, other_aa)

        for new_aa in all_aas:
            if new_aa == original_aa:
                continue

            # delta_phys: change in local Grantham cost at this codon
            new_local = 0.0
            for c1, c2 in incident_edges:
                other = c2 if c1 == codon else c1
                other_aa = code[other]
                if new_aa != other_aa:
                    new_local += grantham_distance(new_aa, other_aa)
            delta_phys = new_local - base_local

            # delta_topo: topology change
            is_breaking, delta_topo = topology_change(
                code, base_components, codon, new_aa
            )

            # delta_trna: Hamming distance proxy for tRNA complexity
            delta_trna = _nearest_hamming_to_aa(codon, new_aa, code)

            is_observed = (
                codon == observed_codon and new_aa == observed_target
            )

            candidates.append(
                CandidateMove(
                    codon=codon,
                    source_aa=original_aa,
                    target_aa=new_aa,
                    delta_phys=delta_phys,
                    topo_breaking=is_breaking,
                    delta_topo=delta_topo,
                    delta_trna=delta_trna,
                    is_observed=is_observed,
                )
            )

    return candidates


# ====================================================================
# Layer C: Statistical inference — conditional logit + model comparison
# ====================================================================


@dataclass
class ModelSpec:
    """Specification of which features to include in the conditional logit."""

    name: str
    use_phys: bool = True
    use_topo: bool = False
    use_trna: bool = False

    @property
    def n_params(self) -> int:
        return sum([self.use_phys, self.use_topo, self.use_trna])


# The four nested models
MODELS = {
    "M1_phys": ModelSpec("M1_phys", use_phys=True, use_topo=False, use_trna=False),
    "M2_topo": ModelSpec("M2_topo", use_phys=False, use_topo=True, use_trna=False),
    "M3_phys_topo": ModelSpec("M3_phys_topo", use_phys=True, use_topo=True, use_trna=False),
    "M4_full": ModelSpec("M4_full", use_phys=True, use_topo=True, use_trna=True),
}


def _feature_vector(move: CandidateMove, spec: ModelSpec) -> np.ndarray:
    """Extract the feature vector for a move under a model specification."""
    feats = []
    if spec.use_phys:
        feats.append(move.delta_phys)
    if spec.use_topo:
        feats.append(float(move.delta_topo))
    if spec.use_trna:
        feats.append(float(move.delta_trna))
    return np.array(feats, dtype=np.float64)


def _normalize_features(choice_sets: list[StepChoiceSet], spec: ModelSpec) -> tuple[np.ndarray, np.ndarray]:
    """Compute means and stds for feature normalization across all candidates.

    Returns (means, stds) arrays of shape (n_features,).
    Normalization improves optimizer convergence for conditional logit.
    """
    all_feats = []
    for cs in choice_sets:
        for m in cs.candidates:
            all_feats.append(_feature_vector(m, spec))
    if not all_feats:
        k = spec.n_params
        return np.zeros(k), np.ones(k)
    mat = np.array(all_feats)
    means = mat.mean(axis=0)
    stds = mat.std(axis=0)
    stds[stds == 0] = 1.0  # Avoid division by zero
    return means, stds


def conditional_logit_log_likelihood(
    weights: np.ndarray,
    choice_sets: list[StepChoiceSet],
    spec: ModelSpec,
    feat_means: np.ndarray | None = None,
    feat_stds: np.ndarray | None = None,
) -> float:
    """Compute log-likelihood of observed choices under conditional logit.

    For each choice set, P(observed | alternatives) = exp(w'x_obs) / sum_j exp(w'x_j).
    Total log-likelihood = sum over choice sets of log P(observed).

    Features are optionally z-scored using provided means/stds.
    """
    ll = 0.0
    for cs in choice_sets:
        obs = cs.observed_move
        if obs is None:
            continue

        # Compute utilities for all candidates
        utilities = []
        obs_utility = None
        for m in cs.candidates:
            x = _feature_vector(m, spec)
            if feat_means is not None and feat_stds is not None:
                x = (x - feat_means) / feat_stds
            u = float(np.dot(weights, x))
            utilities.append(u)
            if m.is_observed:
                obs_utility = u

        if obs_utility is None:
            continue

        # Log-sum-exp for numerical stability
        max_u = max(utilities)
        log_denom = max_u + math.log(sum(math.exp(u - max_u) for u in utilities))
        ll += obs_utility - log_denom

    return ll


def fit_conditional_logit(
    choice_sets: list[StepChoiceSet],
    spec: ModelSpec,
) -> dict:
    """Fit conditional logit model via MLE (scipy.optimize).

    Returns dict with:
        weights: fitted coefficients (on normalized scale)
        weights_raw: coefficients on original feature scale
        log_likelihood: maximized log-likelihood
        n_params: number of parameters
        n_obs: number of choice sets (observed events)
        aic: Akaike Information Criterion
        aicc: corrected AIC (small-sample)
        converged: whether optimizer converged
    """
    k = spec.n_params
    if k == 0:
        # Null model: uniform choice
        ll_null = 0.0
        for cs in choice_sets:
            if cs.observed_move is not None:
                ll_null += -math.log(cs.n_candidates)
        return {
            "model": spec.name,
            "weights": np.array([]),
            "weights_raw": np.array([]),
            "log_likelihood": ll_null,
            "n_params": 0,
            "n_obs": len(choice_sets),
            "aic": -2 * ll_null,
            "aicc": -2 * ll_null,
            "converged": True,
        }

    feat_means, feat_stds = _normalize_features(choice_sets, spec)

    def neg_ll(w: np.ndarray) -> float:
        return -conditional_logit_log_likelihood(
            w, choice_sets, spec, feat_means, feat_stds
        )

    # Start from zeros (uniform prior)
    w0 = np.zeros(k)
    result = minimize(neg_ll, w0, method="Nelder-Mead", options={"maxiter": 5000, "xatol": 1e-8, "fatol": 1e-10})

    # Also try L-BFGS-B if Nelder-Mead finds suboptimal
    try:
        result2 = minimize(neg_ll, result.x, method="L-BFGS-B", options={"maxiter": 5000})
        if result2.fun < result.fun:
            result = result2
    except (ValueError, RuntimeError):
        pass  # Fallback to Nelder-Mead result

    w_norm = result.x
    ll = -result.fun
    n = len([cs for cs in choice_sets if cs.observed_move is not None])

    # Convert weights back to original scale
    w_raw = w_norm / feat_stds

    # AIC and AICc
    aic = -2 * ll + 2 * k
    if n > k + 1:
        aicc = aic + 2 * k * (k + 1) / (n - k - 1)
    else:
        aicc = float("inf")

    return {
        "model": spec.name,
        "weights": w_norm,
        "weights_raw": w_raw,
        "feat_means": feat_means,
        "feat_stds": feat_stds,
        "log_likelihood": ll,
        "n_params": k,
        "n_obs": n,
        "aic": aic,
        "aicc": aicc,
        "converged": result.success,
    }


# ====================================================================
# Choice set construction with order-averaging
# ====================================================================


def _table_events(table_id: int) -> list[tuple[str, str]]:
    """Get (codon, target_aa) pairs for a table's changes from STANDARD."""
    changes = get_changes(table_id)
    return [(codon, aa) for codon, aa in changes.items()]


def build_choice_sets_single_order(
    table_id: int,
    event_order: Sequence[tuple[str, str]],
) -> list[StepChoiceSet]:
    """Build choice sets for a specific ordering of events within a table.

    Starting from STANDARD, apply events in the given order. At each step,
    compute features for all candidate moves from the current code state,
    marking the actual next event as the observed choice.
    """
    code = dict(STANDARD)
    choice_sets: list[StepChoiceSet] = []

    for step_idx, (obs_codon, obs_target) in enumerate(event_order):
        candidates = compute_candidate_features(
            code,
            observed_codon=obs_codon,
            observed_target=obs_target,
        )

        choice_sets.append(
            StepChoiceSet(
                table_id=table_id,
                step_index=step_idx,
                candidates=candidates,
                current_code=dict(code),
            )
        )

        # Apply the observed move
        code[obs_codon] = obs_target

    return choice_sets


def build_choice_sets_order_averaged(
    table_id: int,
    max_orderings: int = 720,
) -> list[list[StepChoiceSet]]:
    """Build choice sets for all orderings of a table's events.

    For a table with k events, generates up to k! orderings (capped at
    max_orderings). Each ordering produces a sequence of k StepChoiceSets.

    Returns list of orderings, each being a list of StepChoiceSets.
    """
    events = _table_events(table_id)
    k = len(events)
    if k == 0:
        return []
    if k == 1:
        return [build_choice_sets_single_order(table_id, events)]

    all_orderings: list[list[StepChoiceSet]] = []
    for i, perm in enumerate(permutations(events)):
        if i >= max_orderings:
            break
        cs = build_choice_sets_single_order(table_id, list(perm))
        all_orderings.append(cs)

    return all_orderings


def order_averaged_log_likelihood(
    weights: np.ndarray,
    orderings: list[list[StepChoiceSet]],
    spec: ModelSpec,
    feat_means: np.ndarray | None = None,
    feat_stds: np.ndarray | None = None,
) -> float:
    """Log-likelihood marginalized over orderings for one table.

    L_table = (1/|orderings|) * sum_o prod_s P(obs_s | candidates_s; w)

    In log space: log L_table = log_sum_exp(ll_o) - log(|orderings|)
    where ll_o is the log-likelihood under ordering o.
    """
    if not orderings:
        return 0.0

    log_liks = []
    for ordering_cs in orderings:
        ll = conditional_logit_log_likelihood(
            weights, ordering_cs, spec, feat_means, feat_stds
        )
        log_liks.append(ll)

    # Log-sum-exp
    max_ll = max(log_liks)
    lse = max_ll + math.log(sum(math.exp(ll - max_ll) for ll in log_liks))
    return lse - math.log(len(orderings))


# ====================================================================
# Full pipeline: build all choice sets, fit models, compare
# ====================================================================


def build_all_choice_sets(
    max_orderings_per_table: int = 720,
    tables: list[int] | None = None,
) -> dict[int, list[list[StepChoiceSet]]]:
    """Build choice sets for all tables with reassignment events.

    Returns dict: table_id -> list of orderings, each an ordered list
    of StepChoiceSets.

    Tables with 0 changes are excluded. Tables with 1 change have a
    single ordering. Tables with k>1 changes have up to k! orderings.
    """
    result: dict[int, list[list[StepChoiceSet]]] = {}
    table_ids = tables or all_table_ids()

    for tid in table_ids:
        events = _table_events(tid)
        if not events:
            continue
        orderings = build_choice_sets_order_averaged(tid, max_orderings_per_table)
        if orderings:
            result[tid] = orderings

    return result


def fit_all_models(
    all_choice_sets: dict[int, list[list[StepChoiceSet]]],
    models: dict[str, ModelSpec] | None = None,
) -> dict[str, dict]:
    """Fit all models and return comparison results.

    For tables with single ordering, uses standard conditional logit.
    For tables with multiple orderings, uses order-averaged likelihood.

    The likelihood for the full dataset is the product across tables.
    """
    if models is None:
        models = MODELS

    results: dict[str, dict] = {}

    for model_name, spec in models.items():
        k = spec.n_params

        # Collect all single-ordering choice sets (for normalization)
        all_cs_flat: list[StepChoiceSet] = []
        for orderings in all_choice_sets.values():
            # Use first ordering for normalization stats
            all_cs_flat.extend(orderings[0])

        if k == 0:
            # Null model
            total_ll = 0.0
            for cs in all_cs_flat:
                if cs.observed_move is not None:
                    total_ll += -math.log(cs.n_candidates)
            n_obs = len([cs for cs in all_cs_flat if cs.observed_move is not None])
            results[model_name] = {
                "model": model_name,
                "weights": np.array([]),
                "log_likelihood": total_ll,
                "n_params": 0,
                "n_obs": n_obs,
                "aic": -2 * total_ll,
                "aicc": -2 * total_ll,
                "converged": True,
            }
            continue

        feat_means, feat_stds = _normalize_features(all_cs_flat, spec)

        def neg_ll_total(w: np.ndarray) -> float:
            total = 0.0
            for orderings in all_choice_sets.values():
                if len(orderings) == 1:
                    # Single ordering: standard conditional logit
                    total += conditional_logit_log_likelihood(
                        w, orderings[0], spec, feat_means, feat_stds
                    )
                else:
                    # Multiple orderings: order-averaged
                    total += order_averaged_log_likelihood(
                        w, orderings, spec, feat_means, feat_stds
                    )
            return -total

        w0 = np.zeros(k)
        opt = minimize(neg_ll_total, w0, method="Nelder-Mead",
                       options={"maxiter": 10000, "xatol": 1e-8, "fatol": 1e-10})
        try:
            opt2 = minimize(neg_ll_total, opt.x, method="L-BFGS-B",
                            options={"maxiter": 10000})
            if opt2.fun < opt.fun:
                opt = opt2
        except (ValueError, RuntimeError):
            pass

        ll = -opt.fun
        w_norm = opt.x
        w_raw = w_norm / feat_stds

        n_obs = len([cs for cs in all_cs_flat if cs.observed_move is not None])
        aic = -2 * ll + 2 * k
        aicc = aic + 2 * k * (k + 1) / (n_obs - k - 1) if n_obs > k + 1 else float("inf")

        # Feature labels
        labels = []
        if spec.use_phys:
            labels.append("delta_phys")
        if spec.use_topo:
            labels.append("delta_topo")
        if spec.use_trna:
            labels.append("delta_trna")

        results[model_name] = {
            "model": model_name,
            "weights": w_norm,
            "weights_raw": w_raw,
            "weight_labels": labels,
            "feat_means": feat_means,
            "feat_stds": feat_stds,
            "log_likelihood": ll,
            "n_params": k,
            "n_obs": n_obs,
            "aic": aic,
            "aicc": aicc,
            "converged": opt.success,
        }

    return results


# ====================================================================
# Diagnostics
# ====================================================================


def observed_move_ranks(
    all_choice_sets: dict[int, list[list[StepChoiceSet]]],
    fit_result: dict,
    spec: ModelSpec,
) -> list[dict]:
    """Compute rank of observed move among candidates under fitted model.

    Uses the first ordering for each table. Rank 1 = highest utility.
    Lower average rank = better model fit.
    """
    weights = fit_result["weights"]
    feat_means = fit_result.get("feat_means")
    feat_stds = fit_result.get("feat_stds")

    ranks: list[dict] = []
    for tid, orderings in all_choice_sets.items():
        for cs in orderings[0]:  # First ordering
            obs = cs.observed_move
            if obs is None:
                continue

            # Score all candidates
            scores = []
            for m in cs.candidates:
                x = _feature_vector(m, spec)
                if feat_means is not None and feat_stds is not None:
                    x = (x - feat_means) / feat_stds
                u = float(np.dot(weights, x))
                scores.append((u, m.is_observed, m.codon, m.target_aa))

            # Rank by utility (descending)
            scores.sort(key=lambda s: -s[0])
            obs_rank = next(
                i + 1 for i, (_, is_obs, _, _) in enumerate(scores) if is_obs
            )
            ranks.append({
                "table_id": tid,
                "step": cs.step_index,
                "codon": obs.codon,
                "target_aa": obs.target_aa,
                "rank": obs_rank,
                "n_candidates": len(scores),
                "percentile": 100.0 * (1 - obs_rank / len(scores)),
            })

    return ranks


def phys_topo_correlation(
    all_choice_sets: dict[int, list[list[StepChoiceSet]]],
) -> dict:
    """Compute correlation between delta_phys and delta_topo across candidates.

    Quantifies confounding: if r is high, topology may be epiphenomenal
    (correlated with physicochemistry). If r is low/moderate, they carry
    independent information.
    """
    from scipy.stats import spearmanr

    phys_vals: list[float] = []
    topo_vals: list[float] = []

    for orderings in all_choice_sets.values():
        for cs in orderings[0]:
            for m in cs.candidates:
                phys_vals.append(m.delta_phys)
                topo_vals.append(float(m.delta_topo))

    if len(phys_vals) < 3:
        return {"spearman_rho": 0.0, "spearman_p": 1.0, "n": len(phys_vals)}

    sr = spearmanr(phys_vals, topo_vals)
    rho_val = float(sr.statistic)  # type: ignore[union-attr]
    p_val = float(sr.pvalue)  # type: ignore[union-attr]
    return {
        "spearman_rho": rho_val,
        "spearman_p": p_val,
        "n": len(phys_vals),
        "interpretation": (
            "High |rho| indicates confounding: topology and physicochemistry "
            "are correlated across the candidate landscape. Low |rho| suggests "
            "they carry independent information."
        ),
    }


def likelihood_ratio_test(
    restricted: dict,
    full: dict,
) -> dict:
    """Likelihood ratio test for nested models.

    H0: restricted model is adequate (extra parameters = 0).
    Test statistic: -2(ll_restricted - ll_full) ~ chi2(df).
    """
    ll_r = restricted["log_likelihood"]
    ll_f = full["log_likelihood"]
    k_r = restricted["n_params"]
    k_f = full["n_params"]
    df = k_f - k_r

    if df <= 0:
        return {"error": "Full model must have more parameters than restricted"}

    lr_stat = -2 * (ll_r - ll_f)
    p_value = 1.0 - chi2.cdf(lr_stat, df)

    return {
        "restricted_model": restricted["model"],
        "full_model": full["model"],
        "lr_statistic": float(lr_stat),
        "df": df,
        "p_value": float(p_value),
        "significant_p05": p_value < 0.05,
    }


def posterior_predictive_topo_rate(
    all_choice_sets: dict[int, list[list[StepChoiceSet]]],
    fit_result: dict,
    spec: ModelSpec,
    n_simulations: int = 10_000,
    seed: int = 135325,
) -> dict:
    """Simulate from fitted model and check topology-breaking rate.

    For each choice set, sample a move proportional to exp(w'x),
    record whether it's topology-breaking. Repeat n_simulations times.
    Compare simulated rate to observed 22%.

    This is the key posterior predictive check: does the model
    reproduce the observed depletion of topology-breaking events?
    """
    import random as rng_mod

    rng = rng_mod.Random(seed)
    weights = fit_result["weights"]
    feat_means = fit_result.get("feat_means")
    feat_stds = fit_result.get("feat_stds")

    # Collect choice sets (first ordering)
    cs_list: list[StepChoiceSet] = []
    for orderings in all_choice_sets.values():
        cs_list.extend(orderings[0])
    cs_list = [cs for cs in cs_list if cs.observed_move is not None]

    # Observed topology-breaking rate
    n_obs_breaking = sum(
        1 for cs in cs_list if cs.observed_move and cs.observed_move.topo_breaking
    )
    obs_rate = n_obs_breaking / max(len(cs_list), 1)

    # Simulate
    sim_rates: list[float] = []
    for _ in range(n_simulations):
        n_breaking = 0
        for cs in cs_list:
            # Compute utilities
            utilities: list[float] = []
            for m in cs.candidates:
                x = _feature_vector(m, spec)
                if feat_means is not None and feat_stds is not None:
                    x = (x - feat_means) / feat_stds
                u = float(np.dot(weights, x))
                utilities.append(u)

            # Softmax -> probabilities
            max_u = max(utilities)
            exp_u = [math.exp(u - max_u) for u in utilities]
            total_exp = sum(exp_u)
            probs = [e / total_exp for e in exp_u]

            # Sample
            r = rng.random()
            cum = 0.0
            chosen_idx = len(probs) - 1
            for j, p in enumerate(probs):
                cum += p
                if r <= cum:
                    chosen_idx = j
                    break

            if cs.candidates[chosen_idx].topo_breaking:
                n_breaking += 1

        sim_rates.append(n_breaking / len(cs_list))

    sim_mean = float(np.mean(sim_rates))
    sim_std = float(np.std(sim_rates))

    # Posterior predictive p-value: P(sim_rate <= obs_rate)
    pp_p = sum(1 for r in sim_rates if r <= obs_rate) / n_simulations

    return {
        "model": fit_result["model"],
        "observed_topo_breaking_rate": obs_rate,
        "simulated_mean_rate": sim_mean,
        "simulated_std_rate": sim_std,
        "posterior_predictive_p": pp_p,
        "n_simulations": n_simulations,
        "n_choice_sets": len(cs_list),
        "interpretation": (
            f"Observed rate: {obs_rate:.3f}. Model predicts: {sim_mean:.3f} "
            f"(+/- {sim_std:.3f}). "
            + (
                "Model reproduces observed depletion (pp p > 0.05)."
                if pp_p > 0.05
                else "Model fails to reproduce observed depletion (pp p < 0.05)."
            )
        ),
    }


# ====================================================================
# Top-level analysis function
# ====================================================================


def run_evolutionary_simulation_analysis(
    max_orderings_per_table: int = 720,
    n_pp_simulations: int = 10_000,
    seed: int = 135325,
) -> dict:
    """Run the complete evolutionary simulation analysis.

    1. Build choice sets for all tables with reassignment events
    2. Fit 4 nested models (M1-M4) via conditional logit MLE
    3. Compare models via AICc and likelihood ratio tests
    4. Compute diagnostics: ranks, confounding, posterior predictive

    Returns comprehensive results dict.
    """
    # Step 1: Build choice sets
    all_cs = build_all_choice_sets(max_orderings_per_table)

    # Step 2: Fit models
    fits = fit_all_models(all_cs)

    # Step 3: Model comparison
    comparisons = {}
    # M1 vs M3 (does topology help beyond phys?)
    if "M1_phys" in fits and "M3_phys_topo" in fits:
        comparisons["M1_vs_M3"] = likelihood_ratio_test(
            fits["M1_phys"], fits["M3_phys_topo"]
        )
    # M2 vs M3 (does phys help beyond topology?)
    if "M2_topo" in fits and "M3_phys_topo" in fits:
        comparisons["M2_vs_M3"] = likelihood_ratio_test(
            fits["M2_topo"], fits["M3_phys_topo"]
        )
    # M3 vs M4 (does tRNA help beyond phys+topo?)
    if "M3_phys_topo" in fits and "M4_full" in fits:
        comparisons["M3_vs_M4"] = likelihood_ratio_test(
            fits["M3_phys_topo"], fits["M4_full"]
        )

    # Step 4: AICc ranking
    aicc_ranking = sorted(
        [(name, f["aicc"]) for name, f in fits.items()],
        key=lambda x: x[1],
    )

    # Step 5: Diagnostics
    diagnostics = {}
    diagnostics["phys_topo_correlation"] = phys_topo_correlation(all_cs)

    for model_name, spec in MODELS.items():
        if model_name in fits:
            diagnostics[f"ranks_{model_name}"] = observed_move_ranks(
                all_cs, fits[model_name], spec
            )

    # Step 6: Posterior predictive for best model
    best_model_name = aicc_ranking[0][0] if aicc_ranking else None
    if best_model_name and best_model_name in fits:
        diagnostics["posterior_predictive"] = posterior_predictive_topo_rate(
            all_cs,
            fits[best_model_name],
            MODELS[best_model_name],
            n_simulations=n_pp_simulations,
            seed=seed,
        )

    # Compile summary
    tables_used = sorted(all_cs.keys())
    total_events = sum(len(orderings[0]) for orderings in all_cs.values())

    return {
        "tables_used": tables_used,
        "n_tables": len(tables_used),
        "total_events": total_events,
        "model_fits": {
            name: {
                "model": f["model"],
                "log_likelihood": f["log_likelihood"],
                "n_params": f["n_params"],
                "aic": f["aic"],
                "aicc": f["aicc"],
                "converged": f["converged"],
                "weights_raw": f.get("weights_raw", np.array([])).tolist(),
                "weight_labels": f.get("weight_labels", []),
            }
            for name, f in fits.items()
        },
        "aicc_ranking": aicc_ranking,
        "likelihood_ratio_tests": comparisons,
        "diagnostics": diagnostics,
    }
