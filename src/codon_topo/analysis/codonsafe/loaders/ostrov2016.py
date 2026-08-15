"""Loader for Ostrov et al. 2016 (Science) 57-codon genome data.

Reference: Ostrov N et al. "Design, synthesis, and testing toward
a 57-codon genome." Science 353:819-822.
doi:10.1126/science.aaf3639

Data sources:
  - Data File S1 (aaf3639_ostrov_re.coli-57.design.gb.zip): GenBank
    file of the designed recoded genome.
  - Table S3 (table_s3.xlsx): Segment-level viability data.
  - Table S4 (table_s4.xlsx): 13 lethal design exceptions.

Key characteristics:
  - 62,214 codons across 7 types recoded:
    AGA->CGU, AGG->CGU (Arg), AGC->UCU, AGU->UCU (Ser),
    UUA->CUG, UUG->CUG (Leu), UAG->UAA (Stop).
  - Ser swaps cross the disconnection boundary (AGY->UCN).
  - Arg swaps stay within one connected component.
  - 55 segments of ~50kb tested; 91% of essential genes retained function.
  - 13 lethal design exceptions identified.

Note: Download supplementary data and place in:
  data/codonsafe/ostrov2016/
"""

from __future__ import annotations

from pathlib import Path

from codon_topo.analysis.codonsafe.models import (
    CodonSwapEvent,
    OutcomeType,
    RecodingOutcome,
    StudyId,
)
from codon_topo.analysis.codonsafe.normalize import normalize_codon

DATA_DIR = (
    Path(__file__).parent.parent.parent.parent / "data" / "codonsafe" / "ostrov2016"
)

# Recoding rules from the paper (DNA alphabet)
OSTROV_RECODE_MAP: dict[str, str] = {
    "AGA": "CGT",  # Arg → Arg (synonymous, no boundary crossing)
    "AGG": "CGT",  # Arg → Arg (synonymous, no boundary crossing)
    "AGC": "TCT",  # Ser → Ser (synonymous, AGY→UCN = boundary crossing)
    "AGT": "TCT",  # Ser → Ser (synonymous, AGY→UCN = boundary crossing)
    "TTA": "CTG",  # Leu → Leu (synonymous)
    "TTG": "CTG",  # Leu → Leu (synonymous)
    "TAG": "TAA",  # Stop → Stop
}


def load_segment_viability() -> list[tuple[CodonSwapEvent, RecodingOutcome]]:
    """Load segment-level viability data from Table S3.

    Each row represents a ~50kb segment with its viability status.
    Individual codon positions within each segment are grouped by unit_id.
    """
    ts3_path = DATA_DIR / "table_s3.xlsx"
    if not ts3_path.exists():
        raise FileNotFoundError(
            f"Ostrov 2016 Table S3 not found at {ts3_path}. "
            f"Download from https://www.science.org/doi/10.1126/science.aaf3639 "
            f"and place in {DATA_DIR}/"
        )

    import pandas as pd

    df = pd.read_excel(ts3_path)
    results: list[tuple[CodonSwapEvent, RecodingOutcome]] = []

    for _, row in df.iterrows():
        segment = str(row.get("Segment", row.get("segment", "")))
        viable = row.get("Viable", row.get("viable", True))
        n_codons = int(row.get("N_codons", row.get("n_codons", 0)) or 0)

        # Use AGA->CGU as representative swap for segment-level events.
        # AGA is the most common forbidden codon in the Ostrov design.
        # Topology classification for segments should use aggregate.py
        # to combine per-codon classifications within each segment.
        event = CodonSwapEvent(
            study=StudyId.OSTROV_2016,
            event_id=f"ostrov_segment_{segment}",
            source_codon=normalize_codon("AGA"),
            target_codon=normalize_codon("CGT"),
            table_id=11,
            organism="E_coli",
            strain="rE.coli-57",
            unit_id=f"segment_{segment}",
            covariates={"n_codons_in_segment": n_codons},
        )

        outcome = RecodingOutcome(
            outcome_type=OutcomeType.BINARY_SUCCESS,
            success=bool(viable),
        )

        results.append((event, outcome))

    return results


def load_lethal_exceptions() -> list[tuple[CodonSwapEvent, RecodingOutcome]]:
    """Load the 13 lethal design exceptions from Table S4.

    These are individual recoded genes that failed to support viability.
    """
    ts4_path = DATA_DIR / "table_s4.xlsx"
    if not ts4_path.exists():
        raise FileNotFoundError(
            f"Ostrov 2016 Table S4 not found at {ts4_path}. "
            f"Download from https://www.science.org/doi/10.1126/science.aaf3639 "
            f"and place in {DATA_DIR}/"
        )

    import pandas as pd

    df = pd.read_excel(ts4_path)
    results: list[tuple[CodonSwapEvent, RecodingOutcome]] = []

    for _, row in df.iterrows():
        gene = str(row.get("Gene", row.get("gene", "")))
        source_dna = str(row.get("WT_Codon", row.get("codon", "NNN")))
        pos = int(row.get("Position", row.get("position", 0)) or 0)

        target_dna = OSTROV_RECODE_MAP.get(source_dna, "NNN")

        try:
            source_rna = normalize_codon(source_dna)
            target_rna = normalize_codon(target_dna)
        except ValueError:
            continue

        event = CodonSwapEvent(
            study=StudyId.OSTROV_2016,
            event_id=f"ostrov_lethal_{gene}_{pos}",
            source_codon=source_rna,
            target_codon=target_rna,
            table_id=11,
            organism="E_coli",
            strain="rE.coli-57",
            gene=gene,
            codon_index_in_cds=pos,
            is_essential_gene=True,
        )

        outcome = RecodingOutcome(
            outcome_type=OutcomeType.BINARY_SUCCESS,
            success=False,
        )

        results.append((event, outcome))

    return results
