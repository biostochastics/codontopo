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


def _nearest_hamming_to_aa(codon: str, target_aa: str) -> int:
    """Hamming distance from codon to nearest codon encoding target_aa in STANDARD.

    Returns 0 if target_aa has no codons in STANDARD (e.g., a synthetic AA label).
    Stop codons are present in STANDARD, so Stop targets get real distances.
    """
    target_codons = [c for c, aa in STANDARD.items() if aa == target_aa]
    if not target_codons:
        return 0
    v = codon_to_vector(codon)
    return min(hamming_distance(v, codon_to_vector(tc)) for tc in target_codons)


def build_reassignment_db() -> list[ReassignmentEvent]:
    """Build complete reassignment database from all NCBI translation tables."""
    events: list[ReassignmentEvent] = []
    for tid in all_table_ids():
        name = get_code_name(tid)
        changes = get_changes(tid)
        for codon, target_aa in changes.items():
            source_aa = STANDARD[codon]
            h = _nearest_hamming_to_aa(codon, target_aa)
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

    For each reassignment event where both source and target are amino acids
    (not Stop), compute the Hamming distance from the reassigned codon to
    the nearest codon already encoding the target AA in the standard code.
    """
    db = build_reassignment_db()
    sense_events = [e for e in db if e.source_aa != "Stop" and e.target_aa != "Stop"]
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
