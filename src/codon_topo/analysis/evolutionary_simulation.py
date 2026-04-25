"""Evolutionary simulation: conditional logit model of codon reassignment.

Tests whether topology avoidance is an independent evolutionary constraint
beyond physicochemical error-minimization, using a discrete-choice framework
(conditional logit) on the 27 observed natural reassignment events.

Architecture (three-layer design):
  Layer A — Graph/state algebra: component-count features via homology.py
  Layer B — Scoring kernels: per-move feature vectors (delta_phys, delta_topo, delta_trna)
  Layer C — Statistical inference: conditional logit MLE, AICc model comparison

The conditional logit framing treats each observed reassignment as a choice
among ~1,280 candidate single-codon reassignments. This converts 27 data
points into 27 choices with thousands of contrasted alternatives, dramatically
improving statistical power over binary (break/no-break) analyses. The
underlying evolutionary assumption is strong-selection-weak-mutation (SSWM)
origin-fixation dynamics, with the conditional logit as the retrospective
statistical test.

For tables with k>1 changes, the likelihood is marginalized over all k!
orderings of events (order-averaging), since temporal ordering is unknown.

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
from typing import Any, Sequence

import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp as _logsumexp
from scipy.stats import chi2

from codon_topo.core.encoding import ALL_CODONS, codon_to_vector, hamming_distance
from codon_topo.core.genetic_codes import STANDARD, all_table_ids, get_changes
from codon_topo.core.homology import connected_components

try:
    from joblib import Parallel, delayed

    _HAS_JOBLIB = True
except ImportError:  # pragma: no cover
    Parallel = None  # type: ignore[assignment]
    delayed = None  # type: ignore[assignment]
    _HAS_JOBLIB = False

# ====================================================================
# Layer A: Graph/state algebra — fast topology features with caching
# ====================================================================

# Cache: frozenset of (codon, aa) items -> {aa: n_components}
_COMPONENT_CACHE: dict[frozenset, dict[str, int]] = {}


def _code_to_key(code: dict[str, str]) -> frozenset:
    """Convert code dict to hashable key for caching."""
    return frozenset(code.items())


def component_counts(code: dict[str, str]) -> dict[str, int]:
    """Count connected components per AA at epsilon=1 in Q_6 (encoding-dependent).

    Uses aggressive caching since many intermediate codes are revisited
    during order-averaging over k! permutations.

    This is the *Q_6* component-count: Hamming-1 adjacency in the default
    GF(2)^6 encoding (C=00, U=01, A=10, G=11). For the encoding-independent
    K_4^3 / H(3,4) component-count, see component_counts_k43.

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


# Cache for K_4^3 (nucleotide-level, encoding-independent) component counts.
_COMPONENT_CACHE_K43: dict[frozenset, dict[str, int]] = {}


def component_counts_k43(code: dict[str, str]) -> dict[str, int]:
    """Count connected components per AA under K_4^3 (encoding-independent).

    Uses nucleotide-level adjacency: two codons are neighbors iff they
    differ at exactly one nucleotide position. This is the encoding-
    independent counterpart of component_counts (which uses Q_6 / Hamming-1
    in the GF(2)^6 encoding). Reviewer R2.M1 / glm-5.1 verification: the
    conditional-logit topology feature delta_topo uses Q_6 by default,
    which is encoding-dependent (8/24 encodings give no Q_6 depletion at
    the landscape level). The K_4^3 variant verifies that M3's dominance
    is not an artifact of the chosen Q_6 encoding.

    Returns dict mapping AA -> number of connected components.
    """
    from codon_topo.core.encoding import nucleotide_distance

    key = _code_to_key(code)
    if key in _COMPONENT_CACHE_K43:
        return _COMPONENT_CACHE_K43[key]

    aa_codons: dict[str, list[str]] = defaultdict(list)
    for codon, aa in code.items():
        if aa != "Stop":
            aa_codons[aa].append(codon)

    result: dict[str, int] = {}
    for aa, codons in aa_codons.items():
        if not codons:
            result[aa] = 0
            continue
        if len(codons) == 1:
            result[aa] = 1
            continue
        parent = list(range(len(codons)))

        def find(x: int) -> int:
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for i in range(len(codons)):
            for j in range(i + 1, len(codons)):
                if nucleotide_distance(codons[i], codons[j]) == 1:
                    pi, pj = find(i), find(j)
                    if pi != pj:
                        parent[pi] = pj
        result[aa] = len({find(i) for i in range(len(codons))})

    _COMPONENT_CACHE_K43[key] = result
    return result


def clear_component_cache() -> None:
    """Clear all component count caches (for testing/memory management)."""
    _COMPONENT_CACHE.clear()
    _COMPONENT_CACHE_K43.clear()


def topology_change(
    base_code: dict[str, str],
    base_components: dict[str, int],
    codon: str,
    new_aa: str,
) -> tuple[bool, int]:
    """Compute topology change from a single-codon reassignment under Q_6.

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


def topology_change_k43(
    base_code: dict[str, str],
    base_components_k43: dict[str, int],
    codon: str,
    new_aa: str,
) -> tuple[bool, int]:
    """Compute K_4^3 topology change from a single-codon reassignment.

    Encoding-independent counterpart of topology_change. Used to verify
    that the conditional-logit M3 dominance is not an artifact of the
    Q_6 encoding choice (per glm-5.1 review).

    Returns:
        (is_breaking_k43, delta_components_k43):
            is_breaking_k43: True if any AA gains a new connected component
                under K_4^3 nucleotide-level adjacency
            delta_components_k43: total increase in components across all AAs
    """
    variant = dict(base_code)
    variant[codon] = new_aa
    var_components = component_counts_k43(variant)

    delta = 0
    is_breaking = False
    all_aas = set(base_components_k43.keys()) | set(var_components.keys())
    for aa in all_aas:
        diff = var_components.get(aa, 0) - base_components_k43.get(aa, 0)
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
    delta_phys: float  # Change in local Grantham mismatch
    topo_breaking: bool  # Binary: creates new component (Q_6)?
    delta_topo: int  # Graded Q_6: total component count increase
    delta_trna: int  # Hamming distance to nearest target-AA codon
    is_observed: bool  # Was this the actually observed move?
    # Encoding-independent (K_4^3) topology — for sensitivity vs Q_6
    topo_breaking_k43: bool = False
    delta_topo_k43: int = 0


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
    base_components_k43 = component_counts_k43(code)
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

            # delta_topo: topology change under Q_6 (encoding-dependent)
            is_breaking, delta_topo = topology_change(
                code, base_components, codon, new_aa
            )

            # delta_topo_k43: topology change under K_4^3 (encoding-independent)
            is_breaking_k43, delta_topo_k43 = topology_change_k43(
                code, base_components_k43, codon, new_aa
            )

            # delta_trna: Hamming distance proxy for tRNA complexity
            delta_trna = _nearest_hamming_to_aa(codon, new_aa, code)

            is_observed = codon == observed_codon and new_aa == observed_target

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
                    topo_breaking_k43=is_breaking_k43,
                    delta_topo_k43=delta_topo_k43,
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
    use_topo: bool = False  # Q_6 topology (delta_topo)
    use_trna: bool = False
    use_topo_k43: bool = False  # K_4^3 topology (delta_topo_k43, encoding-independent)

    @property
    def n_params(self) -> int:
        return sum([self.use_phys, self.use_topo, self.use_trna, self.use_topo_k43])


# The original four nested models (Q_6 topology) plus K_4^3 verification variants.
# Per glm-5.1 review: M3 dominance must be confirmed under encoding-independent
# topology before treating "topology adds independent power" as a primary claim.
MODELS = {
    "M1_phys": ModelSpec("M1_phys", use_phys=True, use_topo=False, use_trna=False),
    "M2_topo": ModelSpec("M2_topo", use_phys=False, use_topo=True, use_trna=False),
    "M3_phys_topo": ModelSpec(
        "M3_phys_topo", use_phys=True, use_topo=True, use_trna=False
    ),
    "M4_full": ModelSpec("M4_full", use_phys=True, use_topo=True, use_trna=True),
    # K_4^3 (encoding-independent) variants:
    "M2_topo_k43": ModelSpec(
        "M2_topo_k43",
        use_phys=False,
        use_topo=False,
        use_trna=False,
        use_topo_k43=True,
    ),
    "M3_phys_topo_k43": ModelSpec(
        "M3_phys_topo_k43",
        use_phys=True,
        use_topo=False,
        use_trna=False,
        use_topo_k43=True,
    ),
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
    if spec.use_topo_k43:
        feats.append(float(move.delta_topo_k43))
    return np.array(feats, dtype=np.float64)


def _normalize_features(
    choice_sets: list[StepChoiceSet], spec: ModelSpec
) -> tuple[np.ndarray, np.ndarray]:
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
    result = minimize(
        neg_ll,
        w0,
        method="Nelder-Mead",
        options={"maxiter": 5000, "xatol": 1e-8, "fatol": 1e-10},
    )

    # Also try L-BFGS-B if Nelder-Mead finds suboptimal
    try:
        result2 = minimize(
            neg_ll, result.x, method="L-BFGS-B", options={"maxiter": 5000}
        )
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
# Vectorized fast path: precompute feature arrays once, then use
# pure-numpy log-likelihood inside the optimizer's neg-LL evaluation.
#
# This replaces the O(N_iter * N_orderings * N_steps * N_candidates)
# Python loop in `neg_ll_total` with batched numpy operations, yielding
# 100-1000x speedup on the 27-table fit. The math is identical to the
# loop-based functions above; those are retained for testing/clarity.
# ====================================================================


@dataclass
class _PrecomputedTable:
    """Stacked feature arrays for one table, all orderings, all steps.

    For a table with `n_ord` orderings (each a sequence of `k` steps), we
    flatten to a single matrix by concatenating all candidate-feature
    vectors. `step_offsets` records where each step begins in the flat
    matrix, and `obs_flat_idx` records the row index of the observed move
    within each step. `ord_step_starts` records, for each ordering, the
    flat index of its first step.

    Shape conventions:
        feat_matrix: (total_candidates, 3)  — columns are [phys, topo, trna]
        step_offsets: (total_steps + 1,)    — half-open [start, end) per step
        obs_flat_idx: (total_steps,)        — observed-row index in feat_matrix
                                             (-1 if no observed move)
        ord_step_starts: (n_ord + 1,)       — half-open [start, end) per ordering
                                             (in units of steps, not rows)
    """

    feat_matrix: np.ndarray
    step_offsets: np.ndarray
    obs_flat_idx: np.ndarray
    ord_step_starts: np.ndarray
    n_orderings: int


def _build_table_bundle(
    orderings: list[list[StepChoiceSet]],
) -> _PrecomputedTable:
    """Pack all (codon, target_aa) candidate features for one table into
    contiguous numpy arrays, ready for vectorized scoring.
    """
    feat_rows: list[np.ndarray] = []
    step_offsets: list[int] = [0]
    obs_flat_idx: list[int] = []
    ord_step_starts: list[int] = [0]

    cumulative_rows = 0
    cumulative_steps = 0
    for ordering_cs in orderings:
        for cs in ordering_cs:
            obs_local: int = -1
            for j, m in enumerate(cs.candidates):
                feat_rows.append(
                    np.array(
                        [
                            m.delta_phys,
                            float(m.delta_topo),
                            float(m.delta_trna),
                            float(m.delta_topo_k43),
                        ],
                        dtype=np.float64,
                    )
                )
                if m.is_observed:
                    obs_local = j
            n_cands = len(cs.candidates)
            obs_flat_idx.append(obs_local + cumulative_rows if obs_local >= 0 else -1)
            cumulative_rows += n_cands
            step_offsets.append(cumulative_rows)
            cumulative_steps += 1
        ord_step_starts.append(cumulative_steps)

    if feat_rows:
        feat_matrix = np.vstack(feat_rows)
    else:
        feat_matrix = np.zeros((0, 4), dtype=np.float64)

    return _PrecomputedTable(
        feat_matrix=feat_matrix,
        step_offsets=np.asarray(step_offsets, dtype=np.int64),
        obs_flat_idx=np.asarray(obs_flat_idx, dtype=np.int64),
        ord_step_starts=np.asarray(ord_step_starts, dtype=np.int64),
        n_orderings=len(orderings),
    )


def _precompute_feature_bundle(
    all_choice_sets: dict[int, list[list[StepChoiceSet]]],
) -> dict[int, _PrecomputedTable]:
    """Build a vectorized bundle for every table at once."""
    return {tid: _build_table_bundle(ords) for tid, ords in all_choice_sets.items()}


def _global_normalization(
    bundle: dict[int, _PrecomputedTable],
) -> tuple[np.ndarray, np.ndarray]:
    """Compute (means, stds) across the union of all candidates in the bundle.

    Returns means/stds of shape (4,) covering [phys, topo, trna, topo_k43]
    columns. Mirrors `_normalize_features` but on the flat numpy matrix.
    """
    rows = [pt.feat_matrix for pt in bundle.values() if pt.feat_matrix.size]
    if not rows:
        return np.zeros(4), np.ones(4)
    mat = np.vstack(rows)
    means = mat.mean(axis=0)
    stds = mat.std(axis=0)
    stds[stds == 0] = 1.0
    return means, stds


# Indices into the (phys, topo, trna, topo_k43) feature columns.
_FEATURE_INDEX = {"phys": 0, "topo": 1, "trna": 2, "topo_k43": 3}


def _spec_columns(spec: ModelSpec) -> list[int]:
    """Map a ModelSpec to the feature-matrix column indices it uses."""
    cols: list[int] = []
    if spec.use_phys:
        cols.append(_FEATURE_INDEX["phys"])
    if spec.use_topo:
        cols.append(_FEATURE_INDEX["topo"])
    if spec.use_trna:
        cols.append(_FEATURE_INDEX["trna"])
    if spec.use_topo_k43:
        cols.append(_FEATURE_INDEX["topo_k43"])
    return cols


def _table_log_likelihood_vec(
    w: np.ndarray,
    pt: _PrecomputedTable,
    cols: list[int],
    feat_means: np.ndarray,
    feat_stds: np.ndarray,
) -> float:
    """Vectorized order-averaged log-likelihood for one precomputed table.

    Replaces the per-step Python loop in
    `conditional_logit_log_likelihood` + `order_averaged_log_likelihood`
    with batched numpy ops:

      1. Slice the feature matrix to the active columns and z-score in
         a single broadcast op.
      2. Compute utilities = X @ w in one matmul over the entire flat
         matrix (not per step).
      3. For each step, the per-step log-likelihood is
            u[obs] - logsumexp(u[step_slice])
         which we evaluate via a numpy-friendly group-logsumexp over
         contiguous index ranges.
      4. Order-average via logsumexp over the per-ordering totals minus
         log(n_orderings).
    """
    if pt.feat_matrix.size == 0:
        return 0.0

    means_active = feat_means[cols]
    stds_active = feat_stds[cols]
    with np.errstate(invalid="ignore", over="ignore", divide="ignore"):
        X = (pt.feat_matrix[:, cols] - means_active) / stds_active  # (N, k)
        u = X @ w  # (N,)

    n_steps = pt.obs_flat_idx.shape[0]
    obs_u = np.where(
        pt.obs_flat_idx >= 0,
        u[np.clip(pt.obs_flat_idx, 0, u.shape[0] - 1)],
        0.0,
    )

    # logsumexp over each step's contiguous slice
    log_denoms = np.empty(n_steps, dtype=np.float64)
    starts = pt.step_offsets[:-1]
    ends = pt.step_offsets[1:]
    for s in range(n_steps):
        log_denoms[s] = _logsumexp(u[starts[s] : ends[s]])
    step_lls = np.where(pt.obs_flat_idx >= 0, obs_u - log_denoms, 0.0)

    # Sum step log-likelihoods within each ordering
    ord_starts = pt.ord_step_starts[:-1]
    ord_ends = pt.ord_step_starts[1:]
    n_ord = pt.n_orderings
    ord_lls = np.empty(n_ord, dtype=np.float64)
    for i in range(n_ord):
        ord_lls[i] = step_lls[ord_starts[i] : ord_ends[i]].sum()

    if n_ord == 1:
        return float(ord_lls[0])
    return float(_logsumexp(ord_lls)) - math.log(n_ord)  # type: ignore[arg-type]


def _total_log_likelihood_vec(
    w_active: np.ndarray,
    bundle: dict[int, _PrecomputedTable],
    cols: list[int],
    feat_means: np.ndarray,
    feat_stds: np.ndarray,
) -> float:
    """Sum vectorized log-likelihood across all tables in the bundle."""
    return sum(
        _table_log_likelihood_vec(w_active, pt, cols, feat_means, feat_stds)
        for pt in bundle.values()
    )


def _fit_one_model_vec(
    model_name: str,
    spec: ModelSpec,
    bundle: dict[int, _PrecomputedTable],
    feat_means_full: np.ndarray,
    feat_stds_full: np.ndarray,
    n_obs: int,
) -> dict:
    """Fit a single conditional-logit model on the precomputed bundle.

    Mirrors the per-model branch of the original `fit_all_models` but
    uses `_total_log_likelihood_vec` so each Nelder-Mead / L-BFGS-B
    function evaluation is a few matrix multiplies rather than 10^5
    Python-level dot products.
    """
    k = spec.n_params
    if k == 0:
        # Null model: uniform choice per step.
        total_ll = 0.0
        for pt in bundle.values():
            for s in range(pt.obs_flat_idx.shape[0] // max(pt.n_orderings, 1)):
                # In the null model order-averaging is irrelevant: every
                # ordering gives the same step structure for the same table,
                # so use the first ordering's step sizes.
                start = pt.step_offsets[s]
                end = pt.step_offsets[s + 1]
                if pt.obs_flat_idx[s] >= 0:
                    total_ll += -math.log(end - start)
        return {
            "model": model_name,
            "weights": np.array([]),
            "weights_raw": np.array([]),
            "log_likelihood": total_ll,
            "n_params": 0,
            "n_obs": n_obs,
            "aic": -2 * total_ll,
            "aicc": -2 * total_ll,
            "converged": True,
        }

    cols = _spec_columns(spec)

    def neg_ll(w: np.ndarray) -> float:
        return -_total_log_likelihood_vec(
            w, bundle, cols, feat_means_full, feat_stds_full
        )

    w0 = np.zeros(k)
    opt = minimize(
        neg_ll,
        w0,
        method="Nelder-Mead",
        options={"maxiter": 10_000, "xatol": 1e-8, "fatol": 1e-10},
    )
    try:
        opt2 = minimize(neg_ll, opt.x, method="L-BFGS-B", options={"maxiter": 10_000})
        if opt2.fun < opt.fun:
            opt = opt2
    except (ValueError, RuntimeError):
        pass

    ll = -opt.fun
    w_norm = opt.x
    feat_stds_active = feat_stds_full[cols]
    feat_means_active = feat_means_full[cols]
    w_raw = w_norm / feat_stds_active

    aic = -2 * ll + 2 * k
    aicc = aic + 2 * k * (k + 1) / (n_obs - k - 1) if n_obs > k + 1 else float("inf")

    labels: list[str] = []
    if spec.use_phys:
        labels.append("delta_phys")
    if spec.use_topo:
        labels.append("delta_topo")
    if spec.use_trna:
        labels.append("delta_trna")
    if spec.use_topo_k43:
        labels.append("delta_topo_k43")

    return {
        "model": model_name,
        "weights": w_norm,
        "weights_raw": w_raw,
        "weight_labels": labels,
        "feat_means": feat_means_active,
        "feat_stds": feat_stds_active,
        "log_likelihood": ll,
        "n_params": k,
        "n_obs": n_obs,
        "aic": aic,
        "aicc": aicc,
        "converged": opt.success,
    }


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
    n_jobs: int = 2,
) -> dict[str, dict]:
    """Fit all models and return comparison results.

    Uses the vectorized fast path: features are precomputed once into
    numpy matrices, then each Nelder-Mead / L-BFGS-B function evaluation
    becomes a small batch of matmul + logsumexp calls.

    Models are fit concurrently using a SHARED-MEMORY threading backend
    (`prefer="threads"`) rather than `processes`, because the precomputed
    feature bundle is several hundred MB and process-level parallelism
    duplicates it per worker (memory-blowing). scipy.optimize releases
    the GIL during BLAS calls, so threads still give meaningful overlap.

    Default `n_jobs=2` keeps RAM bounded; pass `n_jobs=1` for fully
    serial / debug execution. The math is identical to the loop-based
    implementation: for tables with a single ordering, standard
    conditional logit; for multiple orderings, order-averaged
    likelihood. Total likelihood is the product across tables.
    """
    if models is None:
        models = MODELS

    # Precompute feature bundles once and global normalization stats.
    bundle = _precompute_feature_bundle(all_choice_sets)
    feat_means_full, feat_stds_full = _global_normalization(bundle)

    # Count observed events on the first ordering of each table.
    n_obs = 0
    for orderings in all_choice_sets.values():
        for cs in orderings[0]:
            if cs.observed_move is not None:
                n_obs += 1

    items = list(models.items())

    use_parallel = (
        _HAS_JOBLIB
        and n_jobs != 1
        and len(items) > 1
        and Parallel is not None
        and delayed is not None
    )
    if use_parallel:
        assert Parallel is not None and delayed is not None  # narrowed via use_parallel
        # threads (not processes): bundle is shared, no duplication
        outputs = Parallel(n_jobs=n_jobs, prefer="threads", verbose=0)(
            delayed(_fit_one_model_vec)(
                name, spec, bundle, feat_means_full, feat_stds_full, n_obs
            )
            for name, spec in items
        )
    else:  # pragma: no cover  - serial fallback
        outputs = [
            _fit_one_model_vec(
                name, spec, bundle, feat_means_full, feat_stds_full, n_obs
            )
            for name, spec in items
        ]

    out: dict[str, dict] = {}
    for (name, _spec), result in zip(items, outputs):
        if result is None:
            continue
        out[name] = result
    return out


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
            ranks.append(
                {
                    "table_id": tid,
                    "step": cs.step_index,
                    "codon": obs.codon,
                    "target_aa": obs.target_aa,
                    "rank": obs_rank,
                    "n_candidates": len(scores),
                    "percentile": 100.0 * (1 - obs_rank / len(scores)),
                }
            )

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


def run_clade_exclusion_sensitivity(
    max_orderings_per_table: int = 720,
) -> dict:
    """Refit M1-M4 under each of 7 clade-exclusion regimes.

    Reviewer R1.C / R2.M1 (RIGORA Major #1) noted that the topology-avoidance
    hypergeometric/permutation tests have a clade-exclusion sensitivity
    analysis (Sengupta et al. 2007), but the conditional logit does not.
    This routine fits M1-M4 with each major clade dropped and reports
    ΔAICc(M1->M3) and ΔAICc(M2->M3) for each. If both deltas remain large
    across all exclusions, the topology coefficient is robust to phylogenetic
    non-independence; if the effect concentrates in one or two clades, the
    Section 3.5 / 4.2 framing needs softening.

    Clade definitions match `synbio_feasibility.CLADE_GROUPS` so the
    sensitivity is directly comparable to `phylogenetic_sensitivity.json`.
    """
    from codon_topo.analysis.synbio_feasibility import CLADE_GROUPS

    base_table_ids = list(all_table_ids())
    rows: list[dict] = []
    for clade_name, excluded in CLADE_GROUPS.items():
        keep = [t for t in base_table_ids if t not in set(excluded)]
        all_cs = build_all_choice_sets(max_orderings_per_table, tables=keep)
        if not all_cs:
            rows.append(
                {
                    "excluded_clade": clade_name,
                    "excluded_tables": sorted(excluded),
                    "n_tables_remaining": 0,
                    "n_events_remaining": 0,
                    "skipped": "no remaining choice sets",
                }
            )
            continue
        fits = fit_all_models(all_cs)
        m1 = fits.get("M1_phys")
        m2 = fits.get("M2_topo")
        m3 = fits.get("M3_phys_topo")
        m4 = fits.get("M4_full")
        n_events = sum(len(o[0]) for o in all_cs.values())
        if m1 is None or m2 is None or m3 is None:
            rows.append(
                {
                    "excluded_clade": clade_name,
                    "excluded_tables": sorted(excluded),
                    "n_tables_remaining": len(all_cs),
                    "n_events_remaining": n_events,
                    "skipped": "model fit returned None",
                }
            )
            continue
        delta_m1_m3 = m1["aicc"] - m3["aicc"]
        delta_m2_m3 = m2["aicc"] - m3["aicc"]
        delta_m3_m4 = (m4["aicc"] - m3["aicc"]) if m4 is not None else None
        rows.append(
            {
                "excluded_clade": clade_name,
                "excluded_tables": sorted(excluded),
                "n_tables_remaining": len(all_cs),
                "n_events_remaining": n_events,
                "M1_aicc": m1["aicc"],
                "M2_aicc": m2["aicc"],
                "M3_aicc": m3["aicc"],
                "M4_aicc": m4["aicc"] if m4 is not None else None,
                "delta_aicc_M1_to_M3": delta_m1_m3,
                "delta_aicc_M2_to_M3": delta_m2_m3,
                "delta_aicc_M3_to_M4": delta_m3_m4,
                "topo_coef_raw_M3": (
                    m3.get("weights_raw", np.array([])).tolist()
                    if isinstance(m3.get("weights_raw"), np.ndarray)
                    else m3.get("weights_raw")
                ),
                "topo_coef_labels_M3": m3.get("weight_labels", []),
                "M3_favored_over_M1": delta_m1_m3 > 2.0,
                "M3_favored_over_M2": delta_m2_m3 > 2.0,
            }
        )

    # Summary
    valid = [r for r in rows if "delta_aicc_M1_to_M3" in r]
    summary: dict[str, Any] = {
        "method": (
            "Conditional-logit clade-exclusion sensitivity. For each of 7 "
            "phylogenetic-clade exclusion regimes (matching "
            "phylogenetic_sensitivity.json), refit M1-M4 with the "
            "indicated tables removed and report ΔAICc(M1->M3) and "
            "ΔAICc(M2->M3). Reviewer R1.C / R2.M1: tests whether the "
            "topology coefficient is robust to phylogenetic "
            "non-independence."
        ),
        "rows": rows,
        "all_M3_favored_over_M1": (
            all(r["M3_favored_over_M1"] for r in valid) if valid else False
        ),
        "all_M3_favored_over_M2": (
            all(r["M3_favored_over_M2"] for r in valid) if valid else False
        ),
    }
    if valid:
        deltas_m1 = [r["delta_aicc_M1_to_M3"] for r in valid]
        deltas_m2 = [r["delta_aicc_M2_to_M3"] for r in valid]
        summary.update(
            {
                "delta_M1_M3_min": min(deltas_m1),
                "delta_M1_M3_max": max(deltas_m1),
                "delta_M1_M3_median": sorted(deltas_m1)[len(deltas_m1) // 2],
                "delta_M2_M3_min": min(deltas_m2),
                "delta_M2_M3_max": max(deltas_m2),
                "delta_M2_M3_median": sorted(deltas_m2)[len(deltas_m2) // 2],
            }
        )
    return summary


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
    # --- Q_6 (encoding-dependent) topology comparisons (legacy primary) ---
    # M1 vs M3 (does Q_6 topology help beyond phys?)
    if "M1_phys" in fits and "M3_phys_topo" in fits:
        comparisons["M1_vs_M3"] = likelihood_ratio_test(
            fits["M1_phys"], fits["M3_phys_topo"]
        )
    # M2 vs M3 (does phys help beyond Q_6 topology?)
    if "M2_topo" in fits and "M3_phys_topo" in fits:
        comparisons["M2_vs_M3"] = likelihood_ratio_test(
            fits["M2_topo"], fits["M3_phys_topo"]
        )
    # M3 vs M4 (does tRNA help beyond phys+Q_6-topo?)
    if "M3_phys_topo" in fits and "M4_full" in fits:
        comparisons["M3_vs_M4"] = likelihood_ratio_test(
            fits["M3_phys_topo"], fits["M4_full"]
        )
    # --- K_4^3 (encoding-independent) topology comparisons ---
    # Per glm-5.1 review: M3 dominance must be confirmed under
    # encoding-independent topology before treating "topology adds
    # independent power" as a primary claim.
    if "M1_phys" in fits and "M3_phys_topo_k43" in fits:
        comparisons["M1_vs_M3_k43"] = likelihood_ratio_test(
            fits["M1_phys"], fits["M3_phys_topo_k43"]
        )
    if "M2_topo_k43" in fits and "M3_phys_topo_k43" in fits:
        comparisons["M2_k43_vs_M3_k43"] = likelihood_ratio_test(
            fits["M2_topo_k43"], fits["M3_phys_topo_k43"]
        )

    # Step 4: AICc ranking
    aicc_ranking = sorted(
        [(name, f["aicc"]) for name, f in fits.items()],
        key=lambda x: x[1],
    )

    # Step 5: Diagnostics
    diagnostics: dict[str, Any] = {}
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

    # Q_6 vs K_4^3 topology comparison block (manuscript-readable)
    encoding_robustness: dict[str, Any] = {}
    m1 = fits.get("M1_phys")
    m3_q6 = fits.get("M3_phys_topo")
    m3_k43 = fits.get("M3_phys_topo_k43")
    m2_q6 = fits.get("M2_topo")
    m2_k43 = fits.get("M2_topo_k43")
    if m1 and m3_q6 and m3_k43:
        encoding_robustness["delta_aicc_M1_to_M3_q6"] = m1["aicc"] - m3_q6["aicc"]
        encoding_robustness["delta_aicc_M1_to_M3_k43"] = m1["aicc"] - m3_k43["aicc"]
        encoding_robustness["m3_q6_aicc"] = m3_q6["aicc"]
        encoding_robustness["m3_k43_aicc"] = m3_k43["aicc"]
        encoding_robustness["m3_q6_log_likelihood"] = m3_q6["log_likelihood"]
        encoding_robustness["m3_k43_log_likelihood"] = m3_k43["log_likelihood"]
        encoding_robustness["m3_q6_minus_m3_k43_aicc"] = m3_q6["aicc"] - m3_k43["aicc"]
        encoding_robustness["m3_k43_favored_over_m1"] = (
            m1["aicc"] - m3_k43["aicc"]
        ) > 2.0
        encoding_robustness["m3_q6_favored_over_m1"] = (
            m1["aicc"] - m3_q6["aicc"]
        ) > 2.0
    if m2_q6 and m2_k43 and m3_q6 and m3_k43:
        encoding_robustness["delta_aicc_M2q6_to_M3q6"] = m2_q6["aicc"] - m3_q6["aicc"]
        encoding_robustness["delta_aicc_M2k43_to_M3k43"] = (
            m2_k43["aicc"] - m3_k43["aicc"]
        )
    encoding_robustness["interpretation"] = (
        "If delta_aicc_M1_to_M3_k43 > 10 (K_4^3 phys+topo decisively beats "
        "phys alone), M3 dominance is robust to the choice of topology "
        "encoding (Q_6 vs encoding-independent K_4^3 / H(3,4)). If "
        "delta_aicc_M1_to_M3_k43 << delta_aicc_M1_to_M3_q6, the headline "
        "topology effect is partly an artifact of the Q_6 encoding."
    )

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
        "encoding_robustness": encoding_robustness,
    }
