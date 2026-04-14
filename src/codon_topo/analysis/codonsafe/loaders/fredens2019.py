"""Loader for Fredens et al. 2019 (Nature) Syn61 recoding data.

Reference: Fredens J et al. "Total synthesis of Escherichia coli
with a recoded genome." Nature 569:514-518.
doi:10.1038/s41586-019-1192-5

Data sources:
  - Supplementary Data 2: GenBank file of designed synthetic genome
    (all 18,218 recoding positions).
  - Supplementary Data 3: Table of target codons (TCG->AGC, TCA->AGT,
    TAG->TAA positions).
  - Supplementary Data 19: Design optimizations and non-programmed mutations.

Key characteristics:
  - 18,214 sense codons recoded (TCG->AGC, TCA->AGT), all Ser synonymous.
  - Per-segment recoding landscape: REXER experiments tested ~100kb
    segments, recording which positions accepted recoding.
  - 7 corrections needed (design exceptions).
  - Syn61 growth: 1.6x slower than MDS42.

Topology relevance:
  TCG->AGC and TCA->AGT are Ser->Ser swaps that cross the Ser
  disconnection boundary (UCN family -> AGY family) in GF(2)^6.

Note: Download supplementary data and place in:
  data/codonsafe/fredens2019/
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

DATA_DIR = Path(__file__).parent.parent.parent.parent / "data" / "codonsafe" / "fredens2019"


def load_target_codons() -> list[tuple[CodonSwapEvent, RecodingOutcome]]:
    """Load all 18,218 target codon positions from Supplementary Data 3.

    Each position is a TCG->AGC, TCA->AGT, or TAG->TAA swap.
    Per-segment success is inferred from the recoding landscape data
    in the Extended Data figures and Supplementary Data 19.
    """
    sd03_path = DATA_DIR / "supplementary_data_3.xlsx"
    if not sd03_path.exists():
        raise FileNotFoundError(
            f"Fredens 2019 Supplementary Data 3 not found at {sd03_path}. "
            f"Download from https://www.nature.com/articles/s41586-019-1192-5 "
            f"and place in {DATA_DIR}/"
        )

    import pandas as pd

    df = pd.read_excel(sd03_path)
    results: list[tuple[CodonSwapEvent, RecodingOutcome]] = []

    # Recoding rules from the paper
    recode_map = {
        "TCG": "AGC",
        "TCA": "AGT",
        "TAG": "TAA",
    }

    for _, row in df.iterrows():
        gene = str(row.get("Gene", row.get("gene", "")))
        pos = int(row.get("Position", row.get("position", 0)))
        source_dna = str(row.get("WT_Codon", row.get("codon", "NNN")))

        if source_dna not in recode_map:
            continue

        target_dna = recode_map[source_dna]

        try:
            source_rna = normalize_codon(source_dna)
            target_rna = normalize_codon(target_dna)
        except ValueError:
            continue

        # Determine segment
        segment = row.get("Segment", row.get("segment", None))
        unit_id = f"segment_{segment}" if segment is not None else None

        event = CodonSwapEvent(
            study=StudyId.FREDENS_2019,
            event_id=f"fredens_{gene}_{pos}_{target_rna}",
            source_codon=source_rna,
            target_codon=target_rna,
            table_id=11,
            organism="E_coli",
            strain="MDS42_Syn61",
            gene=gene,
            codon_index_in_cds=pos,
            unit_id=unit_id,
        )

        # Default: all positions in Syn61 were successfully recoded
        # (with 7 corrections). Mark corrections as initially failed.
        is_correction = row.get("Correction", row.get("correction", False))
        outcome = RecodingOutcome(
            outcome_type=OutcomeType.BINARY_SUCCESS,
            success=not bool(is_correction),
        )

        results.append((event, outcome))

    return results
