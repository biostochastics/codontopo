"""Integration tests for WS2-WS6 pipeline."""

import pytest
from codon_topo.core.genetic_codes import STANDARD, get_code, get_changes
from codon_topo.core.encoding import codon_to_vector, hamming_distance
from codon_topo.core.homology import disconnection_catalogue
from codon_topo.core.fano import is_fano_line
from codon_topo.analysis.reassignment_db import (
    build_reassignment_db,
)
from codon_topo.analysis.depth_calibration import (
    CALIBRATION_POINTS,
)
from codon_topo.analysis.cosmic_query import (
    fano_predictions_for_kras,
    ws4_gate_decision,
)
from codon_topo.analysis.synbio_feasibility import score_variant_code
from codon_topo.reports.catalogue import build_catalogue, catalogue_summary


class TestWS2Integration:
    """Reassignment database cross-validates with genetic_codes."""

    def test_all_reassignments_match_genetic_codes(self):
        """Every reassignment event should correspond to a get_changes diff."""
        db = build_reassignment_db()
        for e in db:
            changes = get_changes(e.table_id)
            assert e.codon in changes, f"Table {e.table_id}: {e.codon} not in changes"
            assert changes[e.codon] == e.target_aa
            assert STANDARD[e.codon] == e.source_aa

    def test_hamming_paths_consistent_with_encoding(self):
        """Hamming distances in reassignment DB should match encoding module."""
        db = build_reassignment_db()
        for e in db:
            if e.target_aa == "Stop" or e.source_aa == "Stop":
                continue
            v = codon_to_vector(e.codon)
            target_codons = [c for c, aa in STANDARD.items() if aa == e.target_aa]
            min_h = min(
                hamming_distance(v, codon_to_vector(tc)) for tc in target_codons
            )
            assert e.hamming_to_nearest_target == min_h


class TestWS3Integration:
    """Depth calibration cross-validates with disconnection catalogue."""

    def test_calibration_points_match_catalogue(self):
        """Each calibration point's epsilon should match disconnection_catalogue."""
        for p in CALIBRATION_POINTS:
            cat = disconnection_catalogue(get_code(p.table_id))
            matches = [e for e in cat if e["aa"] == p.aa]
            assert len(matches) >= 1, f"Expected {p.aa} in table {p.table_id}"
            assert matches[0]["reconnect_eps"] == p.reconnect_eps


class TestWS4Integration:
    """Fano predictions cross-validate with fano module."""

    def test_fano_predictions_are_valid_fano_lines(self):
        preds = fano_predictions_for_kras()
        for variant, data in preds.items():
            assert is_fano_line(
                data["wt_codon"], data["mutant_codon"], data["fano_partner_codon"]
            ), f"{variant}: not a valid Fano line"

    def test_ws4_gate_with_empty_data(self):
        """Gate should not pass with empty data."""
        result = ws4_gate_decision([])
        assert result["pass"] is False


class TestWS5Integration:
    """Catalogue cross-validates with all workstreams."""

    def test_catalogue_count_stable(self):
        cat = build_catalogue()
        s = catalogue_summary()
        assert s["total_predictions"] == len(cat)

    def test_all_workstreams_represented(self):
        s = catalogue_summary()
        for ws in ["WS1", "WS2", "WS3", "WS4", "WS6"]:
            assert ws in s["by_workstream"], f"{ws} missing from catalogue"


class TestWS6Integration:
    """Synbio scoring cross-validates with WS1 filtration."""

    def test_standard_code_is_optimal(self):
        score = score_variant_code(STANDARD)
        assert score["feasibility_score"] == 1.0

    @pytest.mark.parametrize("table_id", [2, 3, 5, 6])
    def test_variant_codes_have_valid_scores(self, table_id):
        code = get_code(table_id)
        score = score_variant_code(code)
        assert 0.0 <= score["feasibility_score"] <= 1.0
        assert score["serine_disconnected"] is True
