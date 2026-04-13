"""Reassignment database: structured catalogue of codon reassignment events.

Extracts all codon reassignments across NCBI translation tables relative
to the standard genetic code (Table 1). Each event records the codon,
source AA, target AA, and Hamming distance to the nearest codon already
encoding the target AA in the standard code.
"""

from collections import defaultdict
from dataclasses import dataclass
from statistics import mean, median

from scipy.stats import chisquare

from codon_topo.core.encoding import codon_to_vector, hamming_distance
from codon_topo.core.genetic_codes import (
    STANDARD,
    all_table_ids,
    get_changes,
    get_code_name,
)


@dataclass(frozen=True)
class ReassignmentEvent:
    """A single codon reassignment event."""

    table_id: int
    table_name: str
    codon: str
    source_aa: str
    target_aa: str
    hamming_to_nearest_target: int


def _nearest_hamming_to_aa(codon: str, target_aa: str) -> int | None:
    """Hamming distance from codon to nearest codon encoding target_aa in STANDARD.

    Returns None if target_aa has no codons in STANDARD (e.g., a synthetic AA
    label such as Sec / Pyl / O / U used in variant tables). Callers should
    filter None values out before aggregating, otherwise distance statistics
    get biased downward.

    Stop codons are present in STANDARD, so Stop targets get real distances.
    """
    target_codons = [c for c, aa in STANDARD.items() if aa == target_aa]
    if not target_codons:
        return None
    v = codon_to_vector(codon)
    return min(hamming_distance(v, codon_to_vector(tc)) for tc in target_codons)


def build_reassignment_db() -> list[ReassignmentEvent]:
    """Build complete reassignment database from all NCBI translation tables.

    Events where target_aa has no STANDARD codons (non-standard residues
    like Sec/Pyl) are stored with hamming_to_nearest_target=-1 as a sentinel
    so downstream aggregators can filter them. Use hamming_path_lengths()
    for statistics that automatically exclude these.
    """
    events: list[ReassignmentEvent] = []
    for tid in all_table_ids():
        name = get_code_name(tid)
        changes = get_changes(tid)
        for codon, target_aa in changes.items():
            source_aa = STANDARD[codon]
            h_raw = _nearest_hamming_to_aa(codon, target_aa)
            # Sentinel -1 when target has no standard codons (Sec, Pyl, etc.)
            h: int = -1 if h_raw is None else h_raw
            events.append(
                ReassignmentEvent(
                    table_id=tid,
                    table_name=name,
                    codon=codon,
                    source_aa=source_aa,
                    target_aa=target_aa,
                    hamming_to_nearest_target=h,
                )
            )
    return events


def reassignments_by_table() -> dict[int, list[ReassignmentEvent]]:
    """Group reassignment events by NCBI table ID."""
    db = build_reassignment_db()
    result: dict[int, list[ReassignmentEvent]] = defaultdict(list)
    for e in db:
        result[e.table_id].append(e)
    return dict(result)


def reassignments_by_codon() -> dict[str, list[ReassignmentEvent]]:
    """Group reassignment events by codon."""
    db = build_reassignment_db()
    result: dict[str, list[ReassignmentEvent]] = defaultdict(list)
    for e in db:
        result[e.codon].append(e)
    return dict(result)


def hamming_path_lengths() -> dict:
    """Compute Hamming distance statistics for all reassignment paths.

    For each sense-to-sense reassignment, compute the Hamming distance from
    the reassigned codon to the nearest codon already encoding the target
    AA in the standard code. Events with non-standard targets (stored with
    sentinel -1 from _nearest_hamming_to_aa) are excluded to avoid biasing
    statistics downward.
    """
    db = build_reassignment_db()
    sense_events = [
        e
        for e in db
        if e.source_aa != "Stop"
        and e.target_aa != "Stop"
        and e.hamming_to_nearest_target >= 0  # Exclude non-standard-target sentinels
    ]
    excluded_non_standard = sum(
        1
        for e in db
        if e.source_aa != "Stop"
        and e.target_aa != "Stop"
        and e.hamming_to_nearest_target < 0
    )
    event_dicts = [
        {
            "table_id": e.table_id,
            "codon": e.codon,
            "source_aa": e.source_aa,
            "target_aa": e.target_aa,
            "hamming_to_nearest_target": e.hamming_to_nearest_target,
        }
        for e in sense_events
    ]
    distances = [e.hamming_to_nearest_target for e in sense_events]
    return {
        "events": event_dicts,
        "mean_hamming": mean(distances) if distances else 0.0,
        "median_hamming": median(distances) if distances else 0.0,
        "n_events_included": len(sense_events),
        "n_events_excluded_non_standard_target": excluded_non_standard,
    }


def bit_position_bias() -> dict:
    """Test whether reassignments show bias toward specific bit positions.

    For each sense-to-sense reassignment, compute which bits differ between
    the reassigned codon and its nearest target codon. Under null, all 6
    bit positions are equally likely.
    """
    db = build_reassignment_db()
    bit_counts = [0] * 6

    for e in db:
        if e.target_aa == "Stop" or e.source_aa == "Stop":
            continue
        v_source = codon_to_vector(e.codon)
        target_codons = [c for c, aa in STANDARD.items() if aa == e.target_aa]
        if not target_codons:
            continue
        nearest = min(
            target_codons,
            key=lambda tc: hamming_distance(v_source, codon_to_vector(tc)),
        )
        v_target = codon_to_vector(nearest)
        for i in range(6):
            if v_source[i] != v_target[i]:
                bit_counts[i] += 1

    total = sum(bit_counts)
    if total == 0:
        return {"bit_counts": bit_counts, "chi2_statistic": 0.0, "chi2_p_value": 1.0}

    stat, p = chisquare(bit_counts)
    return {
        "bit_counts": bit_counts,
        "chi2_statistic": float(stat),
        "chi2_p_value": float(p),
    }


# Position-level mutation-rate priors (NOT pure Ts/Tv decomposition).
#
# What this model actually does:
#   Each codon position gets a single rate weight reflecting its RELATIVE
#   mutation rate across Ts+Tv. Both GF(2)^6 bits within one nucleotide
#   position receive the SAME weight because the default encoding's per-bit
#   semantics (pyr/pur, amino/keto) do not cleanly partition substitutions
#   into "Ts-only bit" vs "Tv-only bit" — a C↔U transition flips bit 2 only,
#   an A↔U transversion flips both bits, etc. A strict Ts/Tv decomposition
#   would require a different encoding mapping.
#
# The numeric values below are DERIVED from published Ts/Tv ratios by
# combining Ts and Tv rates at each position, but the resulting null is
# best described as "position-weighted" not "Ts/Tv-weighted."
#
# Source: Yang & Yoder 1999 (PMC1207285) and subsequent meta-analyses.
TS_TV_PRIORS = {
    "nuclear": {
        "pos1_ts_tv": 2.5,
        "pos2_ts_tv": 2.8,
        "pos3_ts_tv": 4.5,
    },
    "mitochondrial": {
        "pos1_ts_tv": 8.0,
        "pos2_ts_tv": 6.0,
        "pos3_ts_tv": 15.0,
    },
    "uniform": {
        # Control: equal probability of all bit flips (original test's null)
        "pos1_ts_tv": 1.0,
        "pos2_ts_tv": 1.0,
        "pos3_ts_tv": 1.0,
    },
}


def _expected_bit_weights(compartment: str) -> list[float]:
    """Expected per-bit frequency under the position-weighted neutral null.

    Assigns the same per-position weight to both GF(2)^6 bits within one
    nucleotide position. See the TS_TV_PRIORS comment above for why this
    is not a strict Ts/Tv decomposition.

    Position→bit mapping:
      - bits 0,1 → nucleotide position 1
      - bits 2,3 → nucleotide position 2
      - bits 4,5 → nucleotide position 3 (wobble)
    """
    priors = TS_TV_PRIORS[compartment]
    pos_weights = [
        priors["pos1_ts_tv"],
        priors["pos1_ts_tv"],
        priors["pos2_ts_tv"],
        priors["pos2_ts_tv"],
        priors["pos3_ts_tv"],
        priors["pos3_ts_tv"],
    ]
    total = sum(pos_weights)
    return [w / total for w in pos_weights]


def bit_position_bias_weighted(compartment: str = "mitochondrial") -> dict:
    """Chi-square test for bit-position bias under Ts/Tv-weighted null.

    The original bit_position_bias() test uses a uniform null that is
    biologically unrealistic. This version weights expected bit-flip
    frequencies by position-specific transition/transversion ratios.

    Key question: does the observed p = 0.006 (uniform null) survive a
    realistic null? If yes, the signal is algebraic channeling beyond
    known purifying selection. If no, the signal is mutational spectrum
    rediscovered.

    Args:
        compartment: "nuclear", "mitochondrial", or "uniform" (control).
    """
    db = build_reassignment_db()
    observed = [0] * 6

    for e in db:
        if e.target_aa == "Stop" or e.source_aa == "Stop":
            continue
        v_source = codon_to_vector(e.codon)
        target_codons = [c for c, aa in STANDARD.items() if aa == e.target_aa]
        if not target_codons:
            continue
        nearest = min(
            target_codons,
            key=lambda tc: hamming_distance(v_source, codon_to_vector(tc)),
        )
        v_target = codon_to_vector(nearest)
        for i in range(6):
            if v_source[i] != v_target[i]:
                observed[i] += 1

    total = sum(observed)
    if total == 0:
        return {
            "compartment": compartment,
            "bit_counts_observed": observed,
            "expected_weights": _expected_bit_weights(compartment),
            "chi2_statistic_weighted": 0.0,
            "chi2_p_value_weighted": 1.0,
            "chi2_p_value_uniform_reference": 1.0,
        }

    weights = _expected_bit_weights(compartment)
    expected = [w * total for w in weights]

    # Manual chi-square to avoid scipy issues with zero expected values
    # (which shouldn't happen with Ts/Tv priors but be safe)
    stat_weighted = sum(
        (o - e) ** 2 / e if e > 0 else 0.0 for o, e in zip(observed, expected)
    )
    # Degrees of freedom = k - 1 = 5
    from scipy.stats import chi2

    p_weighted = 1 - chi2.cdf(stat_weighted, df=5)

    # Also report the uniform test for reference
    uniform_result = bit_position_bias()

    return {
        "compartment": compartment,
        "bit_counts_observed": observed,
        "expected_counts_weighted": expected,
        "expected_weights": weights,
        "chi2_statistic_weighted": float(stat_weighted),
        "chi2_p_value_weighted": float(p_weighted),
        "chi2_statistic_uniform_reference": uniform_result["chi2_statistic"],
        "chi2_p_value_uniform_reference": uniform_result["chi2_p_value"],
        "n_events": total,
    }


def directionality_summary() -> dict:
    """Summary statistics for reassignment directionality."""
    db = build_reassignment_db()
    sense_to_sense = [e for e in db if e.source_aa != "Stop" and e.target_aa != "Stop"]
    sense_to_stop = [e for e in db if e.source_aa != "Stop" and e.target_aa == "Stop"]
    stop_to_sense = [e for e in db if e.source_aa == "Stop" and e.target_aa != "Stop"]

    codon_counts: dict[str, int] = defaultdict(int)
    for e in db:
        codon_counts[e.codon] += 1
    top = sorted(codon_counts, key=lambda c: codon_counts[c], reverse=True)[:5]

    return {
        "total_events": len(db),
        "sense_to_sense": len(sense_to_sense),
        "sense_to_stop": len(sense_to_stop),
        "stop_to_sense": len(stop_to_sense),
        "unique_codons": len(set(e.codon for e in db)),
        "unique_source_aas": len(set(e.source_aa for e in db)),
        "unique_target_aas": len(set(e.target_aa for e in db)),
        "top_codons": top,
    }


def _compute_bit_histogram(events: list[ReassignmentEvent]) -> list[int]:
    """Shared helper: compute 6-bin bit-flip histogram for a set of events."""
    counts = [0] * 6
    for e in events:
        if e.target_aa == "Stop" or e.source_aa == "Stop":
            continue
        v_source = codon_to_vector(e.codon)
        target_codons = [c for c, aa in STANDARD.items() if aa == e.target_aa]
        if not target_codons:
            continue
        nearest = min(
            target_codons,
            key=lambda tc: hamming_distance(v_source, codon_to_vector(tc)),
        )
        v_target = codon_to_vector(nearest)
        for i in range(6):
            if v_source[i] != v_target[i]:
                counts[i] += 1
    return counts


def bit_bias_permutation_null(
    n_permutations: int = 10_000,
    seed: int | None = None,
    mode: str = "table_preserving",
) -> dict:
    """Empirical permutation null for bit-position bias.

    Addresses non-iid structure: the same codon can be reassigned across
    multiple tables, and lineage-specific patterns create dependencies.

    Two modes:
      - "table_preserving": within each table, permute target AAs among that
        table's reassigned codons. Preserves which codons moved per table.
      - "codon_preserving": for each codon, permute its targets across tables.
        Preserves codon "hotness".

    Computes the chi-square statistic on the 6-bin bit histogram under each
    permuted dataset and returns an empirical p-value.
    """
    import random

    rng = random.Random(seed)
    db = build_reassignment_db()
    sense_events = [
        e
        for e in db
        if e.source_aa != "Stop"
        and e.target_aa != "Stop"
        and e.hamming_to_nearest_target >= 0
    ]

    observed_hist = _compute_bit_histogram(sense_events)
    total = sum(observed_hist)
    if total == 0:
        return {"observed": observed_hist, "p_value": 1.0, "mode": mode}

    expected_uniform = total / 6
    observed_chi2 = sum(
        (o - expected_uniform) ** 2 / expected_uniform for o in observed_hist
    )

    # Group events for permutation
    if mode == "table_preserving":
        by_table: dict[int, list[ReassignmentEvent]] = defaultdict(list)
        for e in sense_events:
            by_table[e.table_id].append(e)
    elif mode == "codon_preserving":
        by_codon: dict[str, list[ReassignmentEvent]] = defaultdict(list)
        for e in sense_events:
            by_codon[e.codon].append(e)
    else:
        raise ValueError(f"Unknown mode: {mode!r}")

    n_extreme = 0
    for _ in range(n_permutations):
        perm_events = []
        if mode == "table_preserving":
            for _tid, events in by_table.items():
                targets = [e.target_aa for e in events]
                rng.shuffle(targets)
                for e, new_target in zip(events, targets):
                    perm_events.append(
                        ReassignmentEvent(
                            table_id=e.table_id,
                            table_name=e.table_name,
                            codon=e.codon,
                            source_aa=e.source_aa,
                            target_aa=new_target,
                            hamming_to_nearest_target=e.hamming_to_nearest_target,
                        )
                    )
        else:  # codon_preserving
            for _codon, events in by_codon.items():
                targets = [e.target_aa for e in events]
                rng.shuffle(targets)
                for e, new_target in zip(events, targets):
                    perm_events.append(
                        ReassignmentEvent(
                            table_id=e.table_id,
                            table_name=e.table_name,
                            codon=e.codon,
                            source_aa=e.source_aa,
                            target_aa=new_target,
                            hamming_to_nearest_target=e.hamming_to_nearest_target,
                        )
                    )

        perm_hist = _compute_bit_histogram(perm_events)
        perm_total = sum(perm_hist)
        if perm_total == 0:
            continue
        perm_exp = perm_total / 6
        perm_chi2 = sum((o - perm_exp) ** 2 / perm_exp for o in perm_hist)
        if perm_chi2 >= observed_chi2:
            n_extreme += 1

    p_value = (n_extreme + 1) / (n_permutations + 1)
    return {
        "observed_histogram": observed_hist,
        "observed_chi2": observed_chi2,
        "n_permutations": n_permutations,
        "n_extreme": n_extreme,
        "p_value": p_value,
        "mode": mode,
        "n_events": len(sense_events),
    }


def bit_bias_deduplicated() -> dict:
    """Bit-position bias on unique (codon, target_aa) pairs only.

    The raw analysis counts the same codon→target multiple times if it
    appears in multiple NCBI tables (e.g., AGA→Stop in 5+ tables).
    De-duplication tests whether the signal survives when each event
    type is counted only once.
    """
    db = build_reassignment_db()
    seen: set[tuple[str, str]] = set()
    unique_events = []
    for e in db:
        if e.source_aa == "Stop" or e.target_aa == "Stop":
            continue
        if e.hamming_to_nearest_target < 0:
            continue
        key = (e.codon, e.target_aa)
        if key not in seen:
            seen.add(key)
            unique_events.append(e)

    hist_all = _compute_bit_histogram(list(db))
    hist_unique = _compute_bit_histogram(unique_events)
    total_unique = sum(hist_unique)

    if total_unique == 0:
        return {
            "all_events_histogram": hist_all,
            "unique_events_histogram": hist_unique,
            "n_unique": 0,
            "chi2_unique": 0.0,
            "p_value_unique": 1.0,
        }

    stat, p = chisquare(hist_unique)
    return {
        "all_events_histogram": hist_all,
        "n_all": sum(hist_all),
        "unique_events_histogram": hist_unique,
        "n_unique": total_unique,
        "chi2_unique": float(stat),
        "p_value_unique": float(p),
    }


def nucleotide_position_bias() -> dict:
    """Collapse 6 GF(2)^6 bits to 3 nucleotide positions and test bias.

    Bits (0,1) → nucleotide position 1, (2,3) → position 2, (4,5) → position 3.
    A 3-bin test that's biologically legible and encoding-agnostic.
    """
    db = build_reassignment_db()
    bit_hist = _compute_bit_histogram(list(db))
    pos_counts = [
        bit_hist[0] + bit_hist[1],  # position 1
        bit_hist[2] + bit_hist[3],  # position 2
        bit_hist[4] + bit_hist[5],  # position 3 (wobble)
    ]
    total = sum(pos_counts)
    if total == 0:
        return {
            "position_counts": pos_counts,
            "chi2_statistic": 0.0,
            "chi2_p_value": 1.0,
        }

    stat, p = chisquare(pos_counts)
    return {
        "position_counts": pos_counts,
        "position_labels": [
            "pos1 (bits 0-1)",
            "pos2 (bits 2-3)",
            "pos3/wobble (bits 4-5)",
        ],
        "n_events": total,
        "chi2_statistic": float(stat),
        "chi2_p_value": float(p),
        "interpretation": (
            f"Position 2 has {pos_counts[1]} flips, wobble has {pos_counts[2]}. "
            f"Under uniform null: chi2={stat:.2f}, p={p:.4f}."
        ),
    }
