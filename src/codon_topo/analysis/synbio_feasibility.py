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
