import csv

import pytest

from codon_topo.core.genetic_codes import all_table_ids
from codon_topo.visualization.data_export import (
    export_persistent_homology,
    export_embedding_coords,
    export_disconnection_catalogue,
    export_hamming_matrix,
)


def test_export_persistent_homology(tmp_path):
    out = export_persistent_homology(tmp_path / "ph.csv")
    assert out.exists()
    with open(out) as f:
        rows = list(csv.DictReader(f))
    assert len(rows) > 0
    assert "aa" in rows[0]
    assert "epsilon" in rows[0]
    assert "beta_0" in rows[0]
    # Serine should be present
    ser_rows = [r for r in rows if r["aa"] == "Ser"]
    assert len(ser_rows) == 6  # eps 1..6


def test_export_embedding_coords(tmp_path):
    out = export_embedding_coords(tmp_path / "embed.csv")
    assert out.exists()
    with open(out) as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 64


def test_export_disconnection_catalogue(tmp_path):
    out = export_disconnection_catalogue(tmp_path / "disc.csv")
    assert out.exists()
    with open(out) as f:
        rows = list(csv.DictReader(f))
    assert len(rows) > 0
    # Serine should appear in multiple tables
    ser_rows = [r for r in rows if r["aa"] == "Ser"]
    assert len(ser_rows) >= 20


def test_export_disconnection_catalogue_covers_all_27_ncbi_tables(tmp_path):
    """Fig 1B and FigB require every NCBI table id (including 28 and 32) present.

    Regression guard: an earlier disk artifact only listed 25 of 27 tables
    (28 Condylostoma Nuclear and 32 Balanophoraceae Plastid were missing),
    silently contradicting the "all 27 NCBI tables" claim in the manuscript.
    Every table in ``all_table_ids()`` must contribute at least one row to
    the exported catalogue.
    """
    out = export_disconnection_catalogue(tmp_path / "disc.csv")
    with open(out) as f:
        rows = list(csv.DictReader(f))

    expected = set(all_table_ids())
    # sanity: the codebase encodes all 27 NCBI translation tables
    # (codes 1-6, 9-16, 21-33; codes 7, 8, 17-20 are deprecated).
    assert len(expected) == 27
    assert {28, 32}.issubset(expected)

    seen = {int(r["table_id"]) for r in rows}
    missing = expected - seen
    assert not missing, (
        f"disconnection catalogue CSV is missing NCBI tables {sorted(missing)}; "
        "Fig 1B (main text) and FigB (supplement) will render incomplete axes."
    )


def test_export_disconnection_catalogue_raises_when_a_table_would_be_dropped(
    tmp_path, monkeypatch
):
    """The exporter must fail loudly if any table would be silently dropped.

    Simulates a regression (e.g. a filter, or an encoding change that makes
    Ser connected in some table) by monkey-patching the catalogue function
    to return an empty list. The exporter should raise, not silently write
    a truncated CSV.
    """
    from codon_topo.visualization import data_export

    def _empty_catalogue(*args, **kwargs):
        _ = (args, kwargs)  # accept any signature; silence unused-param linters
        return []

    monkeypatch.setattr(data_export, "disconnection_catalogue", _empty_catalogue)
    with pytest.raises(RuntimeError, match="dropped NCBI translation tables"):
        export_disconnection_catalogue(tmp_path / "disc.csv")


def test_export_hamming_matrix(tmp_path):
    out = export_hamming_matrix(tmp_path / "ser_hamming.csv", aa="Ser")
    assert out.exists()
    with open(out) as f:
        reader = csv.reader(f)
        header = next(reader)
    assert len(header) == 7  # '' + 6 serine codons
