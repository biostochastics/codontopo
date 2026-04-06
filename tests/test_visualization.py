import csv

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


def test_export_hamming_matrix(tmp_path):
    out = export_hamming_matrix(tmp_path / "ser_hamming.csv", aa="Ser")
    assert out.exists()
    with open(out) as f:
        reader = csv.reader(f)
        header = next(reader)
    assert len(header) == 7  # '' + 6 serine codons
