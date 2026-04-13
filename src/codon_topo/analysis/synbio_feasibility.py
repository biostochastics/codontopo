"""Synthetic biology feasibility scoring for alternative genetic codes.

Scores each possible single-codon reassignment from the standard code
by how well it preserves the algebraic structure discovered in WS1:
two-fold filtration, four-fold filtration, Serine disconnection, and
overall disconnection profile.

The feasibility score is a weighted composite:
- Two-fold filtration intact: 0.25
- Four-fold filtration intact: 0.25
- Serine disconnected: 0.30
- Single disconnected AA (like standard): 0.20
"""

from codon_topo.core.encoding import ALL_CODONS
from codon_topo.core.genetic_codes import STANDARD
from codon_topo.core.filtration import check_twofold, check_fourfold
from codon_topo.core.homology import disconnection_catalogue


def score_variant_code(code: dict[str, str]) -> dict:
    """Score a genetic code variant for algebraic structure preservation.

    Returns dict with boolean flags and a composite feasibility score [0, 1].
    """
    # Two-fold check
    tw = check_twofold(code)
    tw_intact = all(ok for _, ok, _ in tw) if tw else True

    # Four-fold check
    ff = check_fourfold(code)
    ff_intact = all(ok for _, ok in ff) if ff else True

    # Disconnection analysis
    cat = disconnection_catalogue(code)
    disconnected_aas = [e["aa"] for e in cat]
    ser_disc = "Ser" in disconnected_aas
    n_disc = len(disconnected_aas)

    # Composite score
    score = 0.0
    if tw_intact:
        score += 0.25
    if ff_intact:
        score += 0.25
    if ser_disc:
        score += 0.30
    if n_disc == 1:  # Exactly one disconnected (like standard)
        score += 0.20
    elif n_disc == 0:
        score += 0.10  # Partial credit

    return {
        "twofold_intact": tw_intact,
        "fourfold_intact": ff_intact,
        "serine_disconnected": ser_disc,
        "n_disconnected_aas": n_disc,
        "disconnected_aas": disconnected_aas,
        "feasibility_score": round(score, 4),
    }


def single_reassignment_landscape(
    base_code: dict[str, str] | None = None,
) -> list[dict]:
    """Score every possible single-codon reassignment from the base code.

    For each of the 64 codons, try reassigning it to every other AA
    (including Stop) and score the result.

    Excludes identity reassignments (same AA).
    """
    ref = base_code or STANDARD
    all_aas = sorted(set(ref.values()))  # Includes 'Stop'

    landscape: list[dict] = []
    for codon in ALL_CODONS:
        original_aa = ref[codon]
        for new_aa in all_aas:
            if new_aa == original_aa:
                continue
            variant = dict(ref)
            variant[codon] = new_aa
            sc = score_variant_code(variant)
            landscape.append(
                {
                    "codon": codon,
                    "original_aa": original_aa,
                    "new_aa": new_aa,
                    "feasibility_score": sc["feasibility_score"],
                    "twofold_intact": sc["twofold_intact"],
                    "fourfold_intact": sc["fourfold_intact"],
                    "serine_disconnected": sc["serine_disconnected"],
                    "n_disconnected_aas": sc["n_disconnected_aas"],
                }
            )
    return landscape


def feasibility_summary(
    base_code: dict[str, str] | None = None,
) -> dict:
    """Summary statistics for the single-reassignment landscape.

    Classifies variants into high (>=0.8), medium (0.5-0.8), low (<0.5)
    feasibility tiers.
    """
    landscape = single_reassignment_landscape(base_code)
    high = [e for e in landscape if e["feasibility_score"] >= 0.8]
    medium = [e for e in landscape if 0.5 <= e["feasibility_score"] < 0.8]
    low = [e for e in landscape if e["feasibility_score"] < 0.5]

    best = sorted(landscape, key=lambda x: x["feasibility_score"], reverse=True)[:10]

    return {
        "total_variants": len(landscape),
        "high_feasibility": len(high),
        "medium_feasibility": len(medium),
        "low_feasibility": len(low),
        "best_variants": best,
    }


def topology_avoidance_test() -> dict:
    """Test whether natural reassignments avoid creating new disconnections.

    Compares the rate of disconnection-creating changes among OBSERVED
    natural reassignments vs ALL POSSIBLE single-codon reassignments.

    Prediction: observed events are depleted for topology-breaking changes
    (i.e., observed rate < possible rate).
    """
    from scipy.stats import fisher_exact, hypergeom

    from codon_topo.analysis.reassignment_db import build_reassignment_db

    # 1. What fraction of POSSIBLE reassignments create new disconnections?
    # (Use the standard code as base, matching the original landscape)
    landscape = single_reassignment_landscape(STANDARD)
    standard_cat = disconnection_catalogue(STANDARD)
    standard_disc_aas = {e["aa"] for e in standard_cat}

    possible_creates_disc = 0
    possible_no_disc = 0
    for entry in landscape:
        # Does this variant have MORE disconnected AAs than standard?
        new_disc_aas = set()
        variant = dict(STANDARD)
        variant[entry["codon"]] = entry["new_aa"]
        new_cat = disconnection_catalogue(variant)
        new_disc_aas = {e["aa"] for e in new_cat}
        novel_discs = new_disc_aas - standard_disc_aas
        if novel_discs:
            possible_creates_disc += 1
        else:
            possible_no_disc += 1

    # 2. What fraction of OBSERVED natural reassignments create disconnections?
    db = build_reassignment_db()
    observed_creates_disc = 0
    observed_no_disc = 0
    seen: set[tuple[str, str]] = set()  # de-duplicate
    for e in db:
        key = (e.codon, e.target_aa)
        if key in seen:
            continue
        seen.add(key)
        variant = dict(STANDARD)
        variant[e.codon] = e.target_aa
        new_cat = disconnection_catalogue(variant)
        new_disc_aas = {entry["aa"] for entry in new_cat}
        novel_discs = new_disc_aas - standard_disc_aas
        if novel_discs:
            observed_creates_disc += 1
        else:
            observed_no_disc += 1

    a = observed_creates_disc
    n_obs = observed_creates_disc + observed_no_disc
    N = possible_creates_disc + possible_no_disc
    K = possible_creates_disc

    rate_observed = a / max(n_obs, 1)
    rate_possible = K / max(N, 1)

    # Hypergeometric test: n_obs draws from landscape of N items,
    # K of which are "disc-creating". P(X <= a) tests depletion.
    hypergeom_p = float(hypergeom.cdf(a, N, K, n_obs))

    # Also report Fisher's exact for comparison.
    # Observed is a subset of possible, so the second row is the complement.
    b = observed_no_disc
    c_complement = K - a  # possible disc-creating NOT in observed
    d_complement = (N - K) - b  # possible non-disc NOT in observed
    odds_ratio, fisher_p = fisher_exact(
        [[a, b], [max(c_complement, 0), max(d_complement, 0)]],
        alternative="less",
    )

    # Table-preserving permutation null: within each table, keep which
    # codons change but permute their target AAs among that table's
    # observed targets. Recompute "creates novel disconnections."
    import random as _rng_mod

    from codon_topo.analysis.reassignment_db import build_reassignment_db as _build_db

    by_table: dict[int, list[tuple[str, str]]] = {}
    for e in _build_db():
        by_table.setdefault(e.table_id, []).append((e.codon, e.target_aa))

    rng = _rng_mod.Random(135325)
    n_perm = 10_000
    n_perm_extreme = 0
    for _ in range(n_perm):
        perm_creates = 0
        perm_seen: set[tuple[str, str]] = set()
        for _tid, events in by_table.items():
            codons = [c for c, _t in events]
            targets = [t for _c, t in events]
            rng.shuffle(targets)
            for c, t in zip(codons, targets):
                key = (c, t)
                if key in perm_seen:
                    continue
                perm_seen.add(key)
                variant = dict(STANDARD)
                variant[c] = t
                new_cat = disconnection_catalogue(variant)
                new_disc = {entry["aa"] for entry in new_cat}
                if new_disc - standard_disc_aas:
                    perm_creates += 1
        perm_rate = perm_creates / max(len(perm_seen), 1)
        if perm_rate <= rate_observed:
            n_perm_extreme += 1

    perm_p = (n_perm_extreme + 1) / (n_perm + 1)

    return {
        "observed_creates_disc": a,
        "observed_no_disc": b,
        "observed_total": n_obs,
        "possible_creates_disc": K,
        "possible_no_disc": N - K,
        "possible_total": N,
        "rate_observed": rate_observed,
        "rate_possible": rate_possible,
        "hypergeom_p": hypergeom_p,
        "permutation_p": perm_p,
        "n_permutations": n_perm,
        "odds_ratio": float(odds_ratio),
        "fisher_p": float(fisher_p),
        "hypothesis": (
            "Natural reassignments are depleted for topology-breaking "
            "changes relative to the space of possible reassignments"
        ),
        "caveat": (
            "Hypergeometric p assumes iid draws; permutation p preserves "
            "table structure. Both reported for transparency."
        ),
    }
