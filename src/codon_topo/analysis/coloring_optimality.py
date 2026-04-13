"""Hypercube coloring optimality: the primary publishable direction.

Frames the genetic code as a coloring of the 6-dimensional hypercube Q_6 by
22 classes (21 amino acids + stop) and tests whether this coloring is
locally optimal for minimizing physicochemical mismatch across 1-Hamming
edges.

Key prior: Freeland & Hurst (1998) "The genetic code is one in a million"
(PMID 9732450) showed the standard code is in the top ~0.0001% of random
codes for mutational error minimization. This module replicates that result
in the explicit GF(2)^6 coloring framework and extends it across all 25
NCBI translation tables.

Proposed theorem (gemini3 / glm-5 consensus, April 13 2026):
  Among all colorings of Q_6 with class-size distribution matching the
  standard genetic code, the canonical coloring achieves optimality in the
  top X% for the edge-mismatch objective
      F(C) = sum over {v,w} with d(v,w)=1: Delta(color(v), color(w))
  where Delta is the Grantham (1974) physicochemical distance.
"""

import json
import random
from collections import defaultdict
from itertools import combinations
from pathlib import Path
from statistics import mean, stdev

from codon_topo.core.encoding import (
    ALL_CODONS,
    DEFAULT_ENCODING,
    codon_to_vector,
    hamming_distance,
)
from codon_topo.core.genetic_codes import STANDARD, all_table_ids, get_code

# Three-letter to one-letter AA mapping for Grantham lookups.
AA3_TO_1 = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
}

# Large penalty for stop-codon adjacency - representing lethal mistranslation.
STOP_ADJACENCY_PENALTY = 215.0  # Matches Grantham's maximum distance

# Cache the Grantham matrix at module level to avoid repeated file reads.
_GRANTHAM_CACHE: dict | None = None


def _load_grantham() -> dict:
    """Load the Grantham distance matrix from data/grantham.json."""
    global _GRANTHAM_CACHE
    if _GRANTHAM_CACHE is not None:
        return _GRANTHAM_CACHE

    data_path = Path(__file__).parent.parent / "data" / "grantham.json"
    with open(data_path) as f:
        raw = json.load(f)
    # Strip metadata
    clean = {k: v for k, v in raw.items() if not k.startswith("_")}
    _GRANTHAM_CACHE = clean
    return clean


def grantham_distance(aa1: str, aa2: str, strict: bool = True) -> float:
    """Return Grantham physicochemical distance between two amino acids.

    Accepts either 3-letter or 1-letter codes. Stop codons return a large
    penalty to represent lethal mistranslation cost.

    Caveat: Because stop codons are fixed at their canonical positions in the
    null model, the stop penalty acts as a CONSTANT offset shared across all
    codes. It does NOT affect rank-based quantile comparisons. Use
    `hypercube_edge_mismatch_score_no_stop()` to obtain the stop-adjusted score.

    Args:
        strict: if True (default), raise ValueError for unknown amino acid
            codes. Set False to return 0.0 (legacy behavior, not recommended).
    """
    if aa1 == "Stop" or aa2 == "Stop":
        return STOP_ADJACENCY_PENALTY if aa1 != aa2 else 0.0

    matrix = _load_grantham()
    a1 = AA3_TO_1.get(aa1, aa1)
    a2 = AA3_TO_1.get(aa2, aa2)

    if a1 not in matrix or a2 not in matrix.get(a1, {}):
        if strict:
            raise ValueError(
                f"Unknown amino acid pair: {aa1!r} -> {aa2!r} "
                f"(canonicalized: {a1!r}, {a2!r})"
            )
        return 0.0
    return float(matrix[a1][a2])


# Precomputed Hamming-1 codon pairs under the default encoding — avoids
# recomputing 2016 pair comparisons per Monte Carlo iteration.
_HAMMING_1_EDGES_DEFAULT: list[tuple[str, str]] | None = None


def _build_edges(encoding: dict | None = None) -> list[tuple[str, str]]:
    """Precompute the 192 Hamming-1 edges of Q_6 under a given encoding."""
    enc = encoding or DEFAULT_ENCODING
    codon_vecs = [(c, codon_to_vector(c, enc)) for c in ALL_CODONS]
    edges = [
        (c1, c2)
        for (c1, v1), (c2, v2) in combinations(codon_vecs, 2)
        if hamming_distance(v1, v2) == 1
    ]
    return edges


def _get_default_edges() -> list[tuple[str, str]]:
    global _HAMMING_1_EDGES_DEFAULT
    if _HAMMING_1_EDGES_DEFAULT is None:
        _HAMMING_1_EDGES_DEFAULT = _build_edges(DEFAULT_ENCODING)
    return _HAMMING_1_EDGES_DEFAULT


def hypercube_edge_mismatch_score(
    code: dict[str, str],
    encoding: dict | None = None,
    include_stops: bool = True,
) -> float:
    """Compute total physicochemical mismatch across all Hamming-1 edges.

    F(code) = sum over {v,w} with d(v,w)=1: delta(color(v), color(w))

    Synonymous adjacencies (same AA) contribute 0. Non-synonymous adjacencies
    contribute the Grantham distance between the two amino acids.

    A LOWER score means a more error-correcting coloring (smaller penalty
    when a single-nucleotide mutation happens).

    Args:
        include_stops: if False, skip edges involving stop codons entirely.
            Since stops are typically fixed across null models, their
            contribution is a CONSTANT offset that does not affect relative
            ranking — setting False gives a cleaner score isolating the
            sense-codon placement.
    """
    enc = encoding or DEFAULT_ENCODING
    # Fast path: precomputed edges for the default encoding
    if enc == DEFAULT_ENCODING:
        edges = _get_default_edges()
    else:
        edges = _build_edges(enc)

    total = 0.0
    for c1, c2 in edges:
        aa1, aa2 = code[c1], code[c2]
        if aa1 == aa2:
            continue
        if not include_stops and (aa1 == "Stop" or aa2 == "Stop"):
            continue
        total += grantham_distance(aa1, aa2)
    return total


def hypercube_edge_mismatch_score_no_stop(
    code: dict[str, str],
    encoding: dict | None = None,
) -> float:
    """Score excluding stop-codon adjacencies. Useful for comparing across
    codes that differ in stop placement (e.g., UAG->Leu reassignments).
    """
    return hypercube_edge_mismatch_score(code, encoding, include_stops=False)


def _build_diagonal_edges(encoding: dict | None = None) -> list[tuple[str, str]]:
    """Compute the 96 within-nucleotide distance-2 ('diagonal') edges.

    These are codon pairs that differ at exactly 2 bits, where both
    differing bits fall within the same nucleotide position (bits 0-1,
    2-3, or 4-5). In the codon mutation graph K₄³, these correspond to
    the single-nucleotide substitutions that the Q₆ Hamming-1 model misses.
    """
    enc = encoding or DEFAULT_ENCODING
    codon_vecs = [(c, codon_to_vector(c, enc)) for c in ALL_CODONS]
    edges = []
    for (c1, v1), (c2, v2) in combinations(codon_vecs, 2):
        if hamming_distance(v1, v2) != 2:
            continue
        diff_bits = [i for i in range(6) if v1[i] != v2[i]]
        if len(diff_bits) == 2:
            b1, b2 = diff_bits
            if b2 - b1 == 1 and b1 % 2 == 0:
                edges.append((c1, c2))
    return edges


def weighted_mismatch_score(
    code: dict[str, str],
    rho: float = 1.0,
    encoding: dict | None = None,
    include_stops: bool = True,
) -> float:
    """Mismatch score including diagonal (within-nucleotide d=2) edges.

    F_ρ(code) = F_hamming1(code) + ρ · F_diagonal(code)

    ρ=0 recovers the pure Q₆ model (192 Hamming-1 edges only).
    ρ=1 weights all single-nucleotide substitutions equally (K₄³ model).
    Biologically, ρ < 1 models the typical transition bias (Hamming-1
    edges tend to represent transitions; diagonal edges tend to represent
    the "opposite" transversion).
    """
    h1_score = hypercube_edge_mismatch_score(code, encoding, include_stops)
    if rho == 0.0:
        return h1_score

    enc = encoding or DEFAULT_ENCODING
    diag_edges = _build_diagonal_edges(enc)
    diag_score = 0.0
    for c1, c2 in diag_edges:
        aa1, aa2 = code[c1], code[c2]
        if aa1 == aa2:
            continue
        if not include_stops and (aa1 == "Stop" or aa2 == "Stop"):
            continue
        diag_score += grantham_distance(aa1, aa2)

    return h1_score + rho * diag_score


def rho_robustness_sweep(
    rho_values: list[float] | None = None,
    n_samples: int = 1_000,
    seed: int | None = None,
) -> dict:
    """Sweep over ρ values to test optimality robustness.

    For each ρ, runs a Monte Carlo test under the Freeland-Hurst null
    using the weighted score F_ρ. Reports the quantile and p-value at
    each ρ.

    Addresses reviewer concern: "you ignored 1/3 of mutations."
    """
    if rho_values is None:
        rho_values = [0.0, 0.25, 0.5, 0.75, 1.0]

    results = []
    for rho in rho_values:
        rng = random.Random(seed)
        ref = STANDARD

        observed_f = weighted_mismatch_score(ref, rho=rho)

        null_scores = []
        n_below = 0
        for _ in range(n_samples):
            rc = _generate_random_code_freeland_hurst(ref, rng)
            f = weighted_mismatch_score(rc, rho=rho)
            null_scores.append(f)
            if f < observed_f:
                n_below += 1

        quantile = 100.0 * n_below / n_samples
        p_cons = (n_below + 1) / (n_samples + 1)

        results.append(
            {
                "rho": rho,
                "observed_score": observed_f,
                "null_mean": mean(null_scores),
                "null_std": stdev(null_scores) if len(null_scores) > 1 else 0.0,
                "quantile": quantile,
                "p_value": p_cons,
            }
        )

    return {
        "rho_values": rho_values,
        "per_rho": results,
        "all_significant_p05": all(r["p_value"] < 0.05 for r in results),
        "all_significant_p01": all(r["p_value"] < 0.01 for r in results),
    }


def _degeneracy_profile(code: dict[str, str]) -> list[int]:
    """Extract sorted class-size distribution (excluding stops)."""
    counts: dict[str, int] = defaultdict(int)
    for aa in code.values():
        if aa != "Stop":
            counts[aa] += 1
    return sorted(counts.values(), reverse=True)


def _generate_random_code_class_size_preserving(
    reference_code: dict[str, str],
    rng: random.Random,
) -> dict[str, str]:
    """Generate a random code with matching class-size distribution.

    Shuffles individual codons to amino acid classes with the same
    class-size distribution as the reference. DOES NOT preserve block
    (prefix) structure — synonymous codons are scattered across the
    hypercube rather than kept contiguous. This is a WEAKER null than
    Freeland-Hurst 1998 which preserved block structure.

    Stop codons remain at their fixed positions (same codons remain stops)
    to isolate the AA placement effect.
    """
    sense_codons = [c for c in ALL_CODONS if reference_code[c] != "Stop"]
    stop_codons = [c for c in ALL_CODONS if reference_code[c] == "Stop"]

    # Extract AA->size, sorted for determinism
    aa_to_size: dict[str, int] = defaultdict(int)
    for c in sense_codons:
        aa_to_size[reference_code[c]] += 1
    sorted_items = sorted(aa_to_size.items(), key=lambda x: (-x[1], x[0]))

    shuffled = list(sense_codons)
    rng.shuffle(shuffled)

    perm_code = {c: "Stop" for c in stop_codons}
    idx = 0
    for aa, size in sorted_items:
        for c in shuffled[idx : idx + size]:
            perm_code[c] = aa
        idx += size
    return perm_code


def _generate_random_code_freeland_hurst(
    reference_code: dict[str, str],
    rng: random.Random,
    fix_stop_blocks: bool = True,
) -> dict[str, str]:
    """Generate a random code preserving the BLOCK structure.

    Groups codons by first-2-base prefix (16 blocks of 4 codons each),
    extracts the per-block AA pattern from the reference, then shuffles
    which block gets which pattern. This is the Freeland-Hurst 1998 null
    that preserves synonymous-codon contiguity.

    Args:
        fix_stop_blocks: if True (default), blocks containing stop codons
            are held at their fixed positions (not permuted). This ensures
            stop-codon adjacency contribution is a pure constant across
            null samples (addressing Issue 2 from codex review).
            If False, stop-containing blocks permute among themselves
            (original behavior, moves which codons are stops).
    """
    # Group codons by first-2-base prefix
    blocks: dict[str, list[str]] = defaultdict(list)
    for codon in ALL_CODONS:
        blocks[codon[:2]].append(codon)
    block_keys = sorted(blocks.keys())

    # Extract AA pattern per block (tuple of 4 AAs, sorted by suffix)
    block_patterns: list[tuple[str, ...]] = []
    for prefix in block_keys:
        codons_in_block = sorted(blocks[prefix])
        pattern = tuple(reference_code[c] for c in codons_in_block)
        block_patterns.append(pattern)

    # Separate stop-containing blocks from pure-sense blocks
    stop_blocks = [i for i, p in enumerate(block_patterns) if "Stop" in p]
    sense_blocks = [i for i, p in enumerate(block_patterns) if "Stop" not in p]

    sense_perms = [block_patterns[i] for i in sense_blocks]
    rng.shuffle(sense_perms)

    new_patterns = list(block_patterns)
    if fix_stop_blocks:
        # Leave stop-containing blocks exactly in place — their contribution
        # to edge mismatch is a PROPER constant across null samples
        pass
    else:
        stop_perms = [block_patterns[i] for i in stop_blocks]
        rng.shuffle(stop_perms)
        for new_pat, orig_idx in zip(stop_perms, stop_blocks):
            new_patterns[orig_idx] = new_pat

    for new_pat, orig_idx in zip(sense_perms, sense_blocks):
        new_patterns[orig_idx] = new_pat

    # Map back to codons
    perm_code = {}
    for prefix, pattern in zip(block_keys, new_patterns):
        codons_in_block = sorted(blocks[prefix])
        for c, aa in zip(codons_in_block, pattern):
            perm_code[c] = aa
    return perm_code


def monte_carlo_null(
    n_samples: int = 10_000,
    reference_code: dict[str, str] | None = None,
    seed: int | None = None,
    null_type: str = "freeland_hurst",  # "class_size" or "freeland_hurst" (default)
    include_stops: bool = True,
) -> dict:
    """Sample random colorings and compute F for each.

    null_type="class_size": shuffle individual codons preserving class sizes
        (does NOT preserve synonymous block contiguity — weaker null)
    null_type="freeland_hurst": shuffle AA patterns across fixed blocks
        (preserves block contiguity — matches Freeland & Hurst 1998)

    Returns the observed F(reference_code), null distribution stats, and the
    quantile of the observed score. Lower score = better (less edge mismatch).

    P-value is bounded below by 1/(n_samples+1) to avoid false precision
    when no null samples beat the observed score (conservative estimate).
    """
    ref = reference_code if reference_code is not None else STANDARD
    rng = random.Random(seed)

    if null_type not in ("freeland_hurst", "class_size"):
        raise ValueError(f"Unknown null_type: {null_type!r}")

    observed_f = hypercube_edge_mismatch_score(ref, include_stops=include_stops)

    null_scores: list[float] = []
    n_below_observed = 0
    for _ in range(n_samples):
        if null_type == "freeland_hurst":
            random_code = _generate_random_code_freeland_hurst(ref, rng)
        else:
            random_code = _generate_random_code_class_size_preserving(ref, rng)
        f = hypercube_edge_mismatch_score(random_code, include_stops=include_stops)
        null_scores.append(f)
        if f < observed_f:
            n_below_observed += 1

    null_mean = mean(null_scores)
    null_std = stdev(null_scores) if len(null_scores) > 1 else 0.0
    quantile = 100.0 * n_below_observed / n_samples

    # Conservative p-value: (k + 1) / (n + 1). Avoids reporting p = 0 exactly.
    p_value_conservative = (n_below_observed + 1) / (n_samples + 1)
    p_value_raw = n_below_observed / n_samples
    p_bound = 1.0 / (n_samples + 1)  # minimum resolvable p-value

    if p_value_raw < 0.001:
        interp = (
            f"Strongly optimal (p < {max(p_bound, 0.001):.4f}, bounded by sample size)"
        )
    elif p_value_raw < 0.01:
        interp = f"Strongly optimal (p = {p_value_raw:.4f})"
    elif p_value_raw < 0.05:
        interp = f"Significantly optimal (p = {p_value_raw:.4f})"
    else:
        interp = f"Not distinguishable from null (p = {p_value_raw:.4f})"

    return {
        "null_type": null_type,
        "include_stops": include_stops,
        "observed_score": observed_f,
        "null_mean": null_mean,
        "null_std": null_std,
        "null_min": min(null_scores),
        "null_max": max(null_scores),
        "quantile_of_observed": quantile,
        "n_beaten_observed": n_below_observed,
        "n_samples": n_samples,
        "p_value_raw": p_value_raw,
        "p_value_conservative": p_value_conservative,  # (k+1)/(n+1)
        "p_value_lower_bound": p_bound,
        "interpretation": interp,
    }


def cross_table_optimality() -> dict:
    """Compute F score for every NCBI translation table.

    Tests whether variant codes maintain similar optimality, or whether
    reassignment trades off error-correction for decoding flexibility.
    """
    table_ids = all_table_ids()
    scores: list[float] = []
    profiles: list[list[int]] = []
    for tid in table_ids:
        code = get_code(tid)
        scores.append(hypercube_edge_mismatch_score(code))
        profiles.append(_degeneracy_profile(code))

    standard_score = scores[0]
    per_table = [
        {"table_id": tid, "score": s, "degeneracy_profile": p}
        for tid, s, p in zip(table_ids, scores, profiles)
    ]
    return {
        "per_table": per_table,
        "standard_score": standard_score,
        "min_score": min(scores),
        "max_score": max(scores),
        "relative_to_standard": [
            {
                "table_id": tid,
                "score": s,
                "pct_of_standard": 100.0 * s / standard_score,
            }
            for tid, s in zip(table_ids, scores)
        ],
    }


def encoding_sensitivity_of_optimality(
    n_samples: int = 1_000,
    seed: int | None = None,
    null_type: str = "freeland_hurst",
) -> dict:
    """Check if the optimality result depends on the choice of encoding.

    Runs the Monte Carlo test under each of the 24 base-to-bit encodings
    and reports the quantile range. If the optimality is robust, all 24
    encodings should give similar quantiles.

    Args:
        null_type: "freeland_hurst" (block-preserving, default) or
            "class_size" (class-size-preserving only).
    """
    from codon_topo.core.encoding import all_encodings

    if null_type not in ("freeland_hurst", "class_size"):
        raise ValueError(f"Unknown null_type: {null_type!r}")

    per_encoding_results = []
    quantiles: list[float] = []
    p_values: list[float] = []
    for enc_idx, enc in enumerate(all_encodings()):
        observed_f = hypercube_edge_mismatch_score(STANDARD, encoding=enc)
        # Re-seed per encoding so different encodings get different null draws
        # (rather than identical draws that only differ by the scoring map).
        rng = random.Random((seed or 0) + enc_idx)

        n_below = 0
        enc_scores: list[float] = []
        for _ in range(n_samples):
            if null_type == "freeland_hurst":
                rc = _generate_random_code_freeland_hurst(STANDARD, rng)
            else:
                rc = _generate_random_code_class_size_preserving(STANDARD, rng)
            f = hypercube_edge_mismatch_score(rc, encoding=enc)
            enc_scores.append(f)
            if f < observed_f:
                n_below += 1

        p_raw = n_below / n_samples
        p_cons = (n_below + 1) / (n_samples + 1)
        q = 100.0 * n_below / n_samples
        quantiles.append(q)
        p_values.append(p_raw)
        per_encoding_results.append(
            {
                "encoding_index": enc_idx,
                "encoding": dict(enc),
                "observed_score": observed_f,
                "null_mean": mean(enc_scores),
                "quantile": q,
                "p_value_raw": p_raw,
                "p_value_conservative": p_cons,
            }
        )

    return {
        "n_encodings": len(per_encoding_results),
        "null_type": null_type,
        "per_encoding": per_encoding_results,
        "quantile_range": (min(quantiles), max(quantiles)),
        "mean_quantile": mean(quantiles),
        "all_encodings_significant_p05": all(p < 0.05 for p in p_values),
    }


def local_mismatch_by_codon(
    code: dict[str, str] | None = None,
    encoding: dict | None = None,
) -> dict[str, float]:
    """Compute per-codon local Grantham mismatch cost.

    For each codon, sums the Grantham distance to every Hamming-1 neighbor's
    amino acid. A high local cost means that codon sits in a "bad neighborhood"
    where single-bit mutations cause large physicochemical jumps.
    """
    ref = code if code is not None else STANDARD
    enc = encoding or DEFAULT_ENCODING
    if enc == DEFAULT_ENCODING:
        edges = _get_default_edges()
    else:
        edges = _build_edges(enc)

    local_cost: dict[str, float] = {c: 0.0 for c in ALL_CODONS}
    for c1, c2 in edges:
        aa1, aa2 = ref[c1], ref[c2]
        if aa1 == aa2:
            continue
        d = grantham_distance(aa1, aa2)
        local_cost[c1] += d
        local_cost[c2] += d
    return local_cost


def reassignment_local_cost_test() -> dict:
    """Test whether reassigned codons have higher local Grantham cost.

    Prediction: codons reassigned across NCBI tables sit in "worse"
    neighborhoods (higher local mismatch) than non-reassigned codons.
    Uses Mann-Whitney U test (one-sided, greater).
    """
    from scipy.stats import mannwhitneyu

    from codon_topo.analysis.reassignment_db import build_reassignment_db

    local_costs = local_mismatch_by_codon()
    db = build_reassignment_db()
    reassigned_codons = {e.codon for e in db}

    costs_reassigned = [local_costs[c] for c in reassigned_codons]
    costs_non_reassigned = [
        local_costs[c] for c in ALL_CODONS if c not in reassigned_codons
    ]

    stat, p = mannwhitneyu(
        costs_reassigned, costs_non_reassigned, alternative="greater"
    )

    return {
        "n_reassigned": len(costs_reassigned),
        "n_non_reassigned": len(costs_non_reassigned),
        "mean_cost_reassigned": sum(costs_reassigned) / len(costs_reassigned),
        "mean_cost_non_reassigned": (
            sum(costs_non_reassigned) / len(costs_non_reassigned)
        ),
        "median_cost_reassigned": sorted(costs_reassigned)[len(costs_reassigned) // 2],
        "median_cost_non_reassigned": sorted(costs_non_reassigned)[
            len(costs_non_reassigned) // 2
        ],
        "mann_whitney_U": float(stat),
        "mann_whitney_p": float(p),
        "hypothesis": (
            "Reassigned codons have higher local Grantham mismatch cost "
            "than non-reassigned codons (one-sided, greater)"
        ),
    }


def per_table_optimality(
    n_samples: int = 1_000,
    seed: int | None = None,
) -> dict:
    """Test whether variant codes preserve coloring optimality.

    For each NCBI table, runs a Monte Carlo test using that table's own
    code as reference. Prediction: most variant codes remain in the top
    X% of their own block-preserving null, i.e., reassignments don't
    destroy error minimization.
    """
    results = []
    for tid in all_table_ids():
        code = get_code(tid)
        r = monte_carlo_null(
            n_samples=n_samples,
            reference_code=code,
            seed=(seed or 0) + tid,
            null_type="freeland_hurst",
        )
        results.append(
            {
                "table_id": tid,
                "observed_score": r["observed_score"],
                "quantile": r["quantile_of_observed"],
                "p_value": r["p_value_conservative"],
            }
        )

    significant_raw = sum(1 for r in results if r["p_value"] < 0.05)
    quantiles = [r["quantile"] for r in results]

    # Benjamini-Hochberg FDR correction
    n_tables = len(results)
    sorted_idx = sorted(range(n_tables), key=lambda i: results[i]["p_value"])
    bh_adjusted = [0.0] * n_tables
    for rank_pos, idx in enumerate(sorted_idx):
        bh_adjusted[idx] = results[idx]["p_value"] * n_tables / (rank_pos + 1)
    # Enforce monotonicity (step-up)
    for i in range(len(sorted_idx) - 2, -1, -1):
        idx_curr = sorted_idx[i]
        idx_next = sorted_idx[i + 1]
        bh_adjusted[idx_curr] = min(bh_adjusted[idx_curr], bh_adjusted[idx_next])
    # Cap at 1.0
    bh_adjusted = [min(p, 1.0) for p in bh_adjusted]

    for i, r in enumerate(results):
        r["p_value_bh"] = bh_adjusted[i]

    significant_bh = sum(1 for r in results if r["p_value_bh"] < 0.05)

    return {
        "per_table": results,
        "n_tables": n_tables,
        "n_significant_p05_raw": significant_raw,
        "n_significant_p05_bh": significant_bh,
        "mean_quantile": mean(quantiles),
        "max_quantile": max(quantiles),
        "all_below_5pct": all(q < 5.0 for q in quantiles),
        "hypothesis": (
            "Most variant genetic codes remain in the top 5% of "
            "their own block-preserving null for Grantham edge-mismatch"
        ),
    }


def score_decomposition_by_position(
    code: dict[str, str] | None = None,
    encoding: dict | None = None,
) -> dict:
    """Decompose the mismatch score by nucleotide position and top AA pairs.

    Breaks F(code) into contributions from edges flipping bits in each
    nucleotide position. Also identifies which AA pairs contribute the most
    to the total score. Helps reviewers understand what biology drives the
    statistic.
    """
    ref = code if code is not None else STANDARD
    enc = encoding or DEFAULT_ENCODING
    if enc == DEFAULT_ENCODING:
        edges = _get_default_edges()
    else:
        edges = _build_edges(enc)

    by_position = [0.0, 0.0, 0.0]  # nucleotide positions 1, 2, 3
    by_aa_pair: dict[tuple[str, str], float] = defaultdict(float)

    for c1, c2 in edges:
        aa1, aa2 = ref[c1], ref[c2]
        if aa1 == aa2:
            continue
        d = grantham_distance(aa1, aa2)

        # Which bit differs?
        v1, v2 = codon_to_vector(c1, enc), codon_to_vector(c2, enc)
        for i in range(6):
            if v1[i] != v2[i]:
                by_position[i // 2] += d
                break

        pair_key: tuple[str, str] = (min(aa1, aa2), max(aa1, aa2))
        by_aa_pair[pair_key] += d

    total = sum(by_position)
    top_pairs = sorted(by_aa_pair.items(), key=lambda x: -x[1])[:10]

    return {
        "total_score": total,
        "by_nucleotide_position": {
            "pos1": by_position[0],
            "pos2": by_position[1],
            "pos3_wobble": by_position[2],
        },
        "position_fractions": {
            "pos1": by_position[0] / max(total, 1),
            "pos2": by_position[1] / max(total, 1),
            "pos3_wobble": by_position[2] / max(total, 1),
        },
        "top_aa_pairs": [
            {"pair": list(pair), "score": score, "fraction": score / max(total, 1)}
            for pair, score in top_pairs
        ],
    }


def bootstrap_quantile_ci(
    n_samples: int = 10_000,
    n_bootstrap: int = 100,
    seed: int | None = None,
    confidence_level: float = 0.95,
) -> dict:
    """Bootstrap confidence interval for the observed quantile.

    Repeats the Monte Carlo test with different seeds and reports the
    distribution of quantile estimates. This gives uncertainty on the
    p-value / quantile that a single run cannot.
    """
    quantiles = []
    for i in range(n_bootstrap):
        boot_seed = (seed or 0) + i * 7919  # distinct seeds
        r = monte_carlo_null(
            n_samples=n_samples,
            seed=boot_seed,
            null_type="freeland_hurst",
        )
        quantiles.append(r["quantile_of_observed"])

    quantiles.sort()
    alpha = 1 - confidence_level
    lo_idx = int(len(quantiles) * alpha / 2)
    hi_idx = int(len(quantiles) * (1 - alpha / 2)) - 1

    return {
        "n_bootstrap": n_bootstrap,
        "n_samples_per_boot": n_samples,
        "mean_quantile": mean(quantiles),
        "median_quantile": quantiles[len(quantiles) // 2],
        "ci_low": quantiles[max(0, lo_idx)],
        "ci_high": quantiles[min(len(quantiles) - 1, hi_idx)],
        "confidence_level": confidence_level,
        "z_score": (mean(quantiles) - 50.0) / max(stdev(quantiles), 0.001),
    }


def affine_subspace_constrained_score(code: dict[str, str]) -> dict:
    """Score contribution decomposition: synonymous-block-internal vs inter-block.

    Tests the gemini3/glm-5 proposed refinement: how much of the optimality
    comes from synonymous blocks being affine subspaces (required by four-fold
    filtration) vs. from the specific AA placement on those subspaces?
    """
    total = hypercube_edge_mismatch_score(code)

    # Count synonymous edges vs non-synonymous edges
    n_syn_edges = 0
    n_nonsyn_edges = 0
    for c1, c2 in combinations(ALL_CODONS, 2):
        if hamming_distance(codon_to_vector(c1), codon_to_vector(c2)) == 1:
            if code[c1] == code[c2]:
                n_syn_edges += 1
            else:
                n_nonsyn_edges += 1

    return {
        "total_score": total,
        "n_synonymous_edges": n_syn_edges,
        "n_nonsynonymous_edges": n_nonsyn_edges,
        "total_hamming_1_edges": n_syn_edges + n_nonsyn_edges,
        "synonymous_fraction": n_syn_edges / (n_syn_edges + n_nonsyn_edges),
        "mean_mismatch_per_nonsyn_edge": total / max(n_nonsyn_edges, 1),
    }
