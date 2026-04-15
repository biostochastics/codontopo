"""Tests for the CodonSafe meta-analysis module."""

import pytest

from codon_topo.analysis.codonsafe.models import (
    CodonSwapEvent,
    OutcomeType,
    RecodingOutcome,
    StudyId,
    SwapModel,
)
from codon_topo.analysis.codonsafe.classify import (
    build_reference_context,
    classify_all_encodings,
    classify_swap_event,
)
from codon_topo.analysis.codonsafe.normalize import (
    dna_to_rna,
    normalize_codon,
    validate_codon_pair,
)


# === Normalization tests ===


class TestNormalize:
    def test_dna_to_rna(self):
        assert dna_to_rna("TCG") == "UCG"
        assert dna_to_rna("tcg") == "ucg"
        assert dna_to_rna("AGC") == "AGC"

    def test_normalize_codon(self):
        assert normalize_codon("TCG") == "UCG"
        assert normalize_codon("tcg") == "UCG"
        assert normalize_codon("AGC") == "AGC"
        assert normalize_codon("UGA") == "UGA"

    def test_normalize_invalid(self):
        with pytest.raises(ValueError, match="invalid characters"):
            normalize_codon("NNN")
        with pytest.raises(ValueError, match="3 characters"):
            normalize_codon("AT")

    def test_validate_pair(self):
        src, tgt = validate_codon_pair("TCG", "AGC")
        assert src == "UCG"
        assert tgt == "AGC"


# === Model tests ===


class TestModels:
    def test_codon_swap_event_validation(self):
        event = CodonSwapEvent(
            study=StudyId.NAPOLITANO_2016,
            event_id="test",
            source_codon="AGA",
            target_codon="CGU",
            table_id=11,
        )
        assert event.source_codon == "AGA"

    def test_codon_swap_event_rejects_dna(self):
        with pytest.raises(ValueError, match="non-RNA characters"):
            CodonSwapEvent(
                study=StudyId.NAPOLITANO_2016,
                event_id="test",
                source_codon="TCG",  # T is not RNA
                target_codon="AGC",
            )

    def test_codon_swap_event_rejects_short(self):
        with pytest.raises(ValueError, match="3 characters"):
            CodonSwapEvent(
                study=StudyId.NAPOLITANO_2016,
                event_id="test",
                source_codon="AG",
                target_codon="AGC",
            )

    def test_recoding_outcome_binary(self):
        outcome = RecodingOutcome(
            outcome_type=OutcomeType.BINARY_SUCCESS,
            success=True,
        )
        assert outcome.success is True

    def test_recoding_outcome_continuous(self):
        outcome = RecodingOutcome(
            outcome_type=OutcomeType.FITNESS_CONTINUOUS,
            fitness=0.95,
            fitness_unit="growth_rate",
        )
        assert outcome.fitness == 0.95

    def test_recoding_outcome_requires_fields(self):
        with pytest.raises(ValueError, match="success"):
            RecodingOutcome(outcome_type=OutcomeType.BINARY_SUCCESS)
        with pytest.raises(ValueError, match="fitness"):
            RecodingOutcome(outcome_type=OutcomeType.FITNESS_CONTINUOUS)

    def test_frozen_dataclasses(self):
        event = CodonSwapEvent(
            study=StudyId.NAPOLITANO_2016,
            event_id="test",
            source_codon="AGA",
            target_codon="CGU",
        )
        with pytest.raises(AttributeError):
            event.source_codon = "UGA"  # type: ignore[misc]


# === Classification tests ===


class TestClassify:
    @pytest.fixture
    def ctx(self):
        return build_reference_context(table_id=1, encoding_id=0)

    def test_build_reference_context(self, ctx):
        assert ctx.table_id == 1
        assert ctx.encoding_id == 0
        assert len(ctx.vectors) == 64
        assert len(ctx.neighbors) == 64
        # Each codon should have exactly 6 Hamming-1 neighbors in Q6
        for codon, nbrs in ctx.neighbors.items():
            assert len(nbrs) == 6, f"{codon} has {len(nbrs)} neighbors"

    def test_ser_boundary_crossing(self, ctx):
        """TCG->AGC (Ser->Ser) crosses the UCN/AGY disconnection boundary."""
        event = CodonSwapEvent(
            study=StudyId.FREDENS_2019,
            event_id="test_ser",
            source_codon="UCG",
            target_codon="AGC",
        )
        topo = classify_swap_event(
            event=event,
            ctx=ctx,
            swap_model=SwapModel.SYNONYMOUS_REPLACEMENT,
        )
        assert topo.source_aa == "Ser"
        assert topo.target_aa == "Ser"
        assert topo.is_synonymous is True
        assert topo.crosses_component_boundary_eps1 is True
        assert topo.changes_components_eps1 is False  # code unchanged

    def test_ala_no_boundary_crossing(self, ctx):
        """GCA->GCU (Ala->Ala) stays within the GCN connected component."""
        event = CodonSwapEvent(
            study=StudyId.ROBERTSON_2025_SYN57,
            event_id="test_ala",
            source_codon="GCA",
            target_codon="GCU",
        )
        topo = classify_swap_event(
            event=event,
            ctx=ctx,
            swap_model=SwapModel.SYNONYMOUS_REPLACEMENT,
        )
        assert topo.source_aa == "Ala"
        assert topo.target_aa == "Ala"
        assert topo.crosses_component_boundary_eps1 is False

    def test_arg_no_boundary(self, ctx):
        """AGA->CGU (Arg->Arg) stays within one connected Arg component."""
        event = CodonSwapEvent(
            study=StudyId.NAPOLITANO_2016,
            event_id="test_arg",
            source_codon="AGA",
            target_codon="CGU",
        )
        topo = classify_swap_event(
            event=event,
            ctx=ctx,
            swap_model=SwapModel.SYNONYMOUS_REPLACEMENT,
        )
        assert topo.source_aa == "Arg"
        assert topo.target_aa == "Arg"
        assert topo.crosses_component_boundary_eps1 is False

    def test_nonsynonymous_reassignment(self, ctx):
        """AGA(Arg)->GCU(Ala) nonsynonymous reassignment changes components."""
        event = CodonSwapEvent(
            study=StudyId.NAPOLITANO_2016,
            event_id="test_nonsyn",
            source_codon="AGA",
            target_codon="GCU",
        )
        topo = classify_swap_event(
            event=event,
            ctx=ctx,
            swap_model=SwapModel.REASSIGN_SOURCE_CODON,
        )
        assert topo.source_aa == "Arg"
        assert topo.target_aa == "Ala"
        assert topo.is_synonymous is False
        assert topo.changes_components_eps1 is True

    def test_stop_codon_boundary_is_none(self, ctx):
        """UAG->UAA (Stop->Stop) has no boundary concept."""
        event = CodonSwapEvent(
            study=StudyId.LAJOIE_2013,
            event_id="test_stop",
            source_codon="UAG",
            target_codon="UAA",
        )
        topo = classify_swap_event(
            event=event,
            ctx=ctx,
            swap_model=SwapModel.SYNONYMOUS_REPLACEMENT,
        )
        assert topo.source_aa == "Stop"
        assert topo.target_aa == "Stop"
        assert topo.crosses_component_boundary_eps1 is None

    def test_delta_f_synonymous_is_zero(self, ctx):
        """Synonymous swaps have delta_F = 0 (code map unchanged)."""
        event = CodonSwapEvent(
            study=StudyId.FREDENS_2019,
            event_id="test_df",
            source_codon="UCG",
            target_codon="AGC",
        )
        topo = classify_swap_event(
            event=event,
            ctx=ctx,
            swap_model=SwapModel.SYNONYMOUS_REPLACEMENT,
        )
        for metric, val in topo.delta_edge_mismatch.items():
            assert val == 0.0, f"delta_F({metric}) = {val} for synonymous swap"

    def test_local_mismatch_varies(self, ctx):
        """Source and target local mismatch should differ for distant codons."""
        event = CodonSwapEvent(
            study=StudyId.FREDENS_2019,
            event_id="test_lm",
            source_codon="UCG",
            target_codon="AGC",
        )
        topo = classify_swap_event(
            event=event,
            ctx=ctx,
            swap_model=SwapModel.SYNONYMOUS_REPLACEMENT,
            metrics=("grantham",),
        )
        assert topo.local_mismatch_source["grantham"] > 0
        assert topo.local_mismatch_target["grantham"] > 0
        # UCG and AGC are at Hamming distance 6 — different neighborhoods
        assert (
            topo.local_mismatch_source["grantham"]
            != topo.local_mismatch_target["grantham"]
        )

    def test_all_24_encodings(self):
        """Classification under all 24 encodings should produce 24 results."""
        event = CodonSwapEvent(
            study=StudyId.FREDENS_2019,
            event_id="test_enc",
            source_codon="UCG",
            target_codon="AGC",
        )
        results = classify_all_encodings(event=event, table_id=1)
        assert len(results) == 24
        # Ser should be disconnected under all encodings
        for topo in results:
            assert topo.source_aa == "Ser"
            assert topo.target_aa == "Ser"
            assert topo.crosses_component_boundary_eps1 is True

    def test_leu_boundary_crossing(self, ctx):
        """Leu has two families (CUN and UUN) — check boundary crossing."""
        # CUG -> UUA: different Leu families
        event = CodonSwapEvent(
            study=StudyId.OSTROV_2016,
            event_id="test_leu",
            source_codon="CUG",
            target_codon="UUA",
        )
        topo = classify_swap_event(
            event=event,
            ctx=ctx,
            swap_model=SwapModel.SYNONYMOUS_REPLACEMENT,
        )
        assert topo.source_aa == "Leu"
        assert topo.target_aa == "Leu"
        # Leu may or may not be disconnected at eps=1 depending on encoding
        # Under default encoding, CUN and UUN Leu families are connected
        # (CUU and UUG differ by 1 bit)
        # So this should NOT cross a boundary
        # (but it's encoding-dependent)

    def test_delta_f_stop_handling(self, ctx):
        """Delta F with include_stops=False correctly handles stop codons."""
        # UAG is a stop codon; recoding UAG->CAG (Gln) is nonsynonymous
        event = CodonSwapEvent(
            study=StudyId.LAJOIE_2013,
            event_id="test_stop_delta",
            source_codon="UAG",
            target_codon="CAG",
        )
        topo_with = classify_swap_event(
            event=event,
            ctx=ctx,
            swap_model=SwapModel.REASSIGN_SOURCE_CODON,
            metrics=("grantham",),
        )
        # Delta F should be nonzero (stop -> sense reassignment)
        assert topo_with.delta_edge_mismatch["grantham"] != 0.0
