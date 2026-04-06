"""Regression tests: reproduce all preliminary results exactly.

These values come from the verified output of the preliminary scripts
and PRD Appendix 8. Any change that breaks these tests means the
refactored code disagrees with the original verification.
"""
import pytest
from codon_topo.core.encoding import codon_to_vector
from codon_topo.core.genetic_codes import STANDARD, get_code, all_table_ids
from codon_topo.core.filtration import check_twofold, check_fourfold, analyze_filtration
from codon_topo.core.homology import (
    connected_components, disconnection_catalogue, persistent_homology,
)
from codon_topo.core.embedding import embed_codon
from codon_topo.core.fano import is_fano_line


class TestFiltrationInvariants:
    """PRD Appendix 8, Claims 1-2."""

    @pytest.mark.parametrize("table_id", all_table_ids())
    def test_twofold_100pct(self, table_id):
        code = get_code(table_id)
        results = check_twofold(code)
        failures = [(aa, diff) for aa, ok, diff in results if not ok]
        assert failures == [], f"Table {table_id}: {failures}"

    def test_fourfold_standard(self):
        results = check_fourfold(STANDARD)
        assert all(ok for _, ok in results)

    @pytest.mark.parametrize("table_id", all_table_ids())
    def test_fourfold_breaks_only_on_stop_reassignment(self, table_id):
        """Four-fold breaks only when stops are reassigned to AAs."""
        code = get_code(table_id)
        results = check_fourfold(code)
        failures = [aa for aa, ok in results if not ok]
        if failures:
            # Verify these failures involve stop codon reassignments
            stop_reassigned = {
                c for c in ('UAA','UAG','UGA')
                if code[c] != 'Stop' and STANDARD[c] == 'Stop'
            }
            assert len(stop_reassigned) > 0, (
                f"Table {table_id}: four-fold failures {failures} "
                f"without stop reassignment"
            )


class TestSerineInvariant:
    """PRD Appendix 8, Claim 3."""

    @pytest.mark.parametrize("table_id", all_table_ids())
    def test_serine_disconnected(self, table_id):
        code = get_code(table_id)
        ser = [c for c, aa in code.items() if aa == 'Ser']
        if len(ser) < 2:
            pytest.skip("Serine has <2 codons")
        vectors = [codon_to_vector(c) for c in ser]
        assert connected_components(vectors, 1) > 1

    @pytest.mark.parametrize("table_id", all_table_ids())
    def test_serine_reconnects_at_4(self, table_id):
        code = get_code(table_id)
        ser = [c for c, aa in code.items() if aa == 'Ser']
        if len(ser) < 2:
            pytest.skip("Serine has <2 codons")
        vectors = [codon_to_vector(c) for c in ser]
        assert connected_components(vectors, 4) == 1

    def test_serine_min_inter_block_distance_4(self):
        cat = disconnection_catalogue(STANDARD)
        ser = [e for e in cat if e['aa'] == 'Ser'][0]
        assert ser['min_inter_distance'] == 4


class TestNovelDisconnections:
    """PRD Appendix 8.2."""

    def test_thr_yeast_mito(self):
        cat = disconnection_catalogue(get_code(3))
        thr = [e for e in cat if e['aa'] == 'Thr']
        assert len(thr) == 1
        assert thr[0]['reconnect_eps'] == 2

    def test_leu_chlorophycean(self):
        cat = disconnection_catalogue(get_code(16))
        leu = [e for e in cat if e['aa'] == 'Leu']
        assert len(leu) == 1
        assert leu[0]['reconnect_eps'] == 2

    def test_ala_pachysolen(self):
        cat = disconnection_catalogue(get_code(26))
        ala = [e for e in cat if e['aa'] == 'Ala']
        assert len(ala) == 1
        assert ala[0]['reconnect_eps'] == 3

    def test_ser_candida_three_components(self):
        cat = disconnection_catalogue(get_code(12))
        ser = [e for e in cat if e['aa'] == 'Ser']
        assert len(ser) == 1
        assert ser[0]['n_components'] == 3
        assert ser[0]['reconnect_eps'] == 3


class TestEmbedding:
    """PRD Appendix 8.1, embedding properties."""

    def test_twofold_coords_12_shared(self):
        from collections import defaultdict
        aa_codons = defaultdict(list)
        for c, aa in STANDARD.items():
            if aa != 'Stop':
                aa_codons[aa].append(c)
        for aa, codons in aa_codons.items():
            if len(codons) != 2:
                continue
            z1, z2 = embed_codon(codons[0]), embed_codon(codons[1])
            assert z1[0] == z2[0] and z1[1] == z2[1], f"{aa}"

    def test_fourfold_third_coord_exhausts_roots(self):
        from collections import defaultdict
        aa_codons = defaultdict(list)
        for c, aa in STANDARD.items():
            if aa != 'Stop':
                aa_codons[aa].append(c)
        roots = {complex(1,0), complex(0,1), complex(-1,0), complex(0,-1)}
        for aa, codons in aa_codons.items():
            if len(codons) != 4:
                continue
            thirds = set(embed_codon(c)[2] for c in codons)
            assert thirds == roots, f"{aa}"


class TestKRAS:
    """PRD Appendix 8.1, KRAS Fano line."""

    def test_ggu_guu_cac(self):
        assert is_fano_line('GGU', 'GUU', 'CAC')


class TestSerineEncodingInvariance:
    """Serine disconnection should hold across all 24 base-to-bit encodings."""

    def test_serine_disconnected_all_encodings(self):
        from codon_topo.core.encoding import all_encodings
        for enc in all_encodings():
            ser_codons = [c for c, aa in STANDARD.items() if aa == 'Ser']
            vectors = [codon_to_vector(c, enc) for c in ser_codons]
            n_comp = connected_components(vectors, 1)
            assert n_comp > 1, f"Serine connected at eps=1 with encoding {enc}"


class TestNullModelValidation:
    """Validate that null models produce low p-values for observed invariants."""

    @pytest.mark.slow
    def test_model_a_serine_rare(self):
        from codon_topo.analysis.null_models import null_model_a
        result = null_model_a(n_permutations=1000, seed=42)
        assert result['p_value_serine_unique'] < 0.05, (
            f"Serine uniqueness p={result['p_value_serine_unique']} >= 0.05"
        )
