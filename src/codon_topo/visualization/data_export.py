"""Export codon_topo results as CSV/JSON for R visualization scripts.

All figures are rendered in R (ggplot2 + ggpubr) per PRD requirements.
This module provides data export functions that produce tidy data frames
as CSV files ready for consumption by the R scripts in visualization/R/.
"""
import csv
import json
from pathlib import Path
from collections import defaultdict

from codon_topo.core.encoding import (
    codon_to_vector, hamming_distance, ALL_CODONS, DEFAULT_ENCODING,
)
from codon_topo.core.genetic_codes import STANDARD, get_code, get_code_name, all_table_ids
from codon_topo.core.homology import persistent_homology, disconnection_catalogue
from codon_topo.core.embedding import embed_codon


def export_persistent_homology(
    output_path: str | Path,
    code: dict[str, str] | None = None,
    encoding: dict | None = None,
) -> Path:
    """Export persistent homology data for barcode plots.

    Writes CSV with columns: aa, degeneracy, epsilon, beta_0.
    """
    ref = code or STANDARD
    enc = encoding or DEFAULT_ENCODING
    output = Path(output_path)

    aa_codons: dict[str, list[str]] = defaultdict(list)
    for codon, aa in ref.items():
        if aa != 'Stop':
            aa_codons[aa].append(codon)

    rows = []
    for aa in sorted(aa_codons):
        codons = aa_codons[aa]
        if len(codons) < 2:
            continue
        vectors = [codon_to_vector(c, enc) for c in codons]
        ph = persistent_homology(vectors, max_eps=6)
        for eps, beta_0 in ph.items():
            rows.append({
                'aa': aa,
                'degeneracy': len(codons),
                'epsilon': eps,
                'beta_0': beta_0,
            })

    with open(output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['aa', 'degeneracy', 'epsilon', 'beta_0'])
        writer.writeheader()
        writer.writerows(rows)
    return output


def export_embedding_coords(
    output_path: str | Path,
    code: dict[str, str] | None = None,
) -> Path:
    """Export C^3 embedding coordinates for scatter plots.

    Writes CSV with columns: codon, aa, z1_re, z1_im, z2_re, z2_im, z3_re, z3_im.
    """
    ref = code or STANDARD
    output = Path(output_path)

    rows = []
    for codon in ALL_CODONS:
        aa = ref[codon]
        z = embed_codon(codon)
        rows.append({
            'codon': codon,
            'aa': aa,
            'z1_re': z[0].real, 'z1_im': z[0].imag,
            'z2_re': z[1].real, 'z2_im': z[1].imag,
            'z3_re': z[2].real, 'z3_im': z[2].imag,
        })

    with open(output, 'w', newline='') as f:
        writer = csv.DictWriter(
            f, fieldnames=['codon', 'aa', 'z1_re', 'z1_im', 'z2_re', 'z2_im', 'z3_re', 'z3_im'],
        )
        writer.writeheader()
        writer.writerows(rows)
    return output


def export_disconnection_catalogue(
    output_path: str | Path,
) -> Path:
    """Export disconnection catalogue across all NCBI tables.

    Writes CSV with columns: table_id, table_name, aa, n_components,
    reconnect_eps, min_inter_distance.
    """
    output = Path(output_path)
    rows = []

    for tid in all_table_ids():
        code = get_code(tid)
        name = get_code_name(tid)
        cat = disconnection_catalogue(code)
        for entry in cat:
            rows.append({
                'table_id': tid,
                'table_name': name,
                'aa': entry['aa'],
                'n_components': entry['n_components'],
                'reconnect_eps': entry['reconnect_eps'],
                'min_inter_distance': entry['min_inter_distance'],
            })

    with open(output, 'w', newline='') as f:
        writer = csv.DictWriter(
            f, fieldnames=['table_id', 'table_name', 'aa', 'n_components',
                           'reconnect_eps', 'min_inter_distance'],
        )
        writer.writeheader()
        writer.writerows(rows)
    return output


def export_hamming_matrix(
    output_path: str | Path,
    aa: str,
    code: dict[str, str] | None = None,
) -> Path:
    """Export pairwise Hamming distance matrix for an amino acid's codons.

    Writes CSV suitable for heatmap visualization.
    """
    ref = code or STANDARD
    output = Path(output_path)
    codons = sorted(c for c, a in ref.items() if a == aa)
    vectors = [codon_to_vector(c) for c in codons]

    with open(output, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([''] + codons)
        for i, c1 in enumerate(codons):
            row = [c1] + [hamming_distance(vectors[i], vectors[j]) for j in range(len(codons))]
            writer.writerow(row)
    return output
