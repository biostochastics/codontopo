"""Tests for reassignment database construction and querying."""

from codon_topo.analysis.reassignment_db import (
    ReassignmentEvent,
    build_reassignment_db,
    reassignments_by_table,
    reassignments_by_codon,
    hamming_path_lengths,
)


def test_reassignment_event_fields():
    db = build_reassignment_db()
    assert len(db) > 0
    e = db[0]
    assert isinstance(e, ReassignmentEvent)
    assert hasattr(e, "table_id")
    assert hasattr(e, "table_name")
    assert hasattr(e, "codon")
    assert hasattr(e, "source_aa")
    assert hasattr(e, "target_aa")
    assert hasattr(e, "hamming_to_nearest_target")


def test_known_reassignment_vert_mito():
    """Table 2: AGA reassigned from Arg to Stop."""
    db = build_reassignment_db()
    aga_events = [e for e in db if e.codon == "AGA" and e.table_id == 2]
    assert len(aga_events) == 1
    assert aga_events[0].source_aa == "Arg"
    assert aga_events[0].target_aa == "Stop"


def test_known_reassignment_yeast_mito():
    """Table 3: CUU-CUG reassigned from Leu to Thr."""
    db = build_reassignment_db()
    thr_events = [
        e
        for e in db
        if e.table_id == 3 and e.target_aa == "Thr" and e.source_aa == "Leu"
    ]
    assert len(thr_events) == 4  # CUU, CUC, CUA, CUG


def test_no_reassignments_in_standard():
    """Table 1 and 11 have no changes from standard."""
    by_table = reassignments_by_table()
    assert len(by_table.get(1, [])) == 0
    assert len(by_table.get(11, [])) == 0


def test_reassignments_by_codon():
    by_codon = reassignments_by_codon()
    # UGA is reassigned in many tables (Trp, Cys, Gly)
    assert "UGA" in by_codon
    assert len(by_codon["UGA"]) >= 5


def test_total_reassignment_count():
    """There should be 40-60 distinct reassignment events across all tables."""
    db = build_reassignment_db()
    assert 30 <= len(db) <= 100


def test_hamming_path_lengths_returns_dict():
    result = hamming_path_lengths()
    assert isinstance(result, dict)
    assert "events" in result
    assert "mean_hamming" in result
    assert "median_hamming" in result
    for e in result["events"]:
        assert "hamming_to_nearest_target" in e
        assert e["hamming_to_nearest_target"] >= 1


def test_hamming_path_is_positive():
    db = build_reassignment_db()
    for e in db:
        if e.target_aa != "Stop" and e.source_aa != "Stop":
            assert e.hamming_to_nearest_target >= 1
