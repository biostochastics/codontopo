"""Loader for Robertson et al. 2025 (Science) Syn57 recoding data.

Reference: Robertson WE et al. "Escherichia coli with a 57-codon
genetic code." Science 390:eady4368 (2025).
bioRxiv: 10.1101/2025.05.02.651837

This is the most important new dataset for CodonSafe:
  - ~101,302 codons recoded (4 Ser + 2 Ala + 1 Stop = 7 codons removed)
  - 38 fragments of ~100kb each with per-fragment fitness landscape
  - Recoding-fitness linkage map (fig. S20) maps each position to fitness
  - 4 of 6 Ser codons eliminated — massive perturbation of Ser's
    already-disconnected codon graph in GF(2)^6

Recoding scheme (vS33A7):
  TCG -> AGC (Ser)   TCT -> AGC (Ser)
  TCC -> AGC (Ser)   TCA -> AGT (Ser)
  GCA -> GCT (Ala)   GCG -> GCC (Ala)
  TAG -> TAA (Stop)

Topology significance:
  - Ser: 4 codons removed from both UCN and AGY families. Remaining
    Ser codons are only UCC and UCU (both in the UCN family).
    This collapses the disconnected Ser graph to a single connected
    component — a natural experiment in topology disruption.
  - Ala: 2 of 4 codons removed (GCA, GCG). Remaining: GCU, GCC.
    The Ala codon graph remains connected (all in GCN box).

Note: Download supplementary data and place in:
  data/codonsafe/robertson2025_syn57/
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
    Path(__file__).parent.parent.parent.parent
    / "data"
    / "codonsafe"
    / "robertson2025_syn57"
)

# Recoding rules for Syn57 (vS33A7 scheme), DNA alphabet
SYN57_RECODE_MAP: dict[str, str] = {
    "TCG": "AGC",
    "TCT": "AGC",
    "TCC": "AGC",
    "TCA": "AGT",
    "GCA": "GCT",
    "GCG": "GCC",
    "TAG": "TAA",
}

# The 38 fragments (~100kb each)
N_FRAGMENTS = 38


def load_recoding_landscape() -> list[tuple[CodonSwapEvent, RecodingOutcome]]:
    """Load per-fragment recoding landscape data from Supplementary Data.

    Each fragment has a fitness score derived from the recoding-fitness
    linkage map (fig. S20, data S8). Fragments that required troubleshooting
    are marked as initially failed (success=False).

    The actual data format will be finalized once the supplementary
    Excel files are downloaded.
    """
    # Try multiple possible filenames
    candidates = [
        DATA_DIR / "supplementary_data_1.xlsx",
        DATA_DIR / "data_s8.xlsx",
        DATA_DIR / "syn57_recoding_landscape.xlsx",
    ]

    data_path = None
    for c in candidates:
        if c.exists():
            data_path = c
            break

    if data_path is None:
        raise FileNotFoundError(
            f"Syn57 supplementary data not found in {DATA_DIR}/. "
            f"Download from https://www.biorxiv.org/content/10.1101/2025.05.02.651837v1 "
            f"or https://www.science.org/doi/10.1126/science.ady4368 "
            f"and place in {DATA_DIR}/"
        )

    import pandas as pd

    df = pd.read_excel(data_path)
    results: list[tuple[CodonSwapEvent, RecodingOutcome]] = []

    for _, row in df.iterrows():
        gene = str(row.get("Gene", row.get("gene", "")))
        pos = int(row.get("Position", row.get("position", 0)))
        source_dna = str(row.get("WT_Codon", row.get("codon", "")))

        if source_dna.upper() not in SYN57_RECODE_MAP:
            continue

        target_dna = SYN57_RECODE_MAP[source_dna.upper()]

        try:
            source_rna = normalize_codon(source_dna)
            target_rna = normalize_codon(target_dna)
        except ValueError:
            continue

        fragment = row.get("Fragment", row.get("fragment", None))
        unit_id = f"fragment_{fragment}" if fragment is not None else None

        # Fitness from linkage map (if available)
        fitness = row.get("Fitness", row.get("fitness", None))
        recoding_maintained = row.get(
            "Recoding_Maintained", row.get("recoding_maintained", None)
        )

        event = CodonSwapEvent(
            study=StudyId.ROBERTSON_2025_SYN57,
            event_id=f"syn57_{gene}_{pos}_{target_rna}",
            source_codon=source_rna,
            target_codon=target_rna,
            table_id=11,
            organism="E_coli",
            strain="Syn57",
            gene=gene,
            codon_index_in_cds=pos,
            unit_id=unit_id,
        )

        if fitness is not None:
            outcome = RecodingOutcome(
                outcome_type=OutcomeType.FITNESS_CONTINUOUS,
                fitness=float(fitness),
                fitness_unit="linkage_fitness_score",
            )
        elif recoding_maintained is not None:
            outcome = RecodingOutcome(
                outcome_type=OutcomeType.BINARY_SUCCESS,
                success=bool(recoding_maintained),
            )
        else:
            continue  # skip events with no outcome data

        results.append((event, outcome))

    return results


def get_ser_topology_summary() -> dict:
    """Summarize the Ser topology changes in Syn57 vs standard code.

    In the standard code, Ser has 6 codons in two disconnected families:
      UCN: UCU, UCC, UCA, UCG (4 codons)
      AGY: AGU, AGC (2 codons)

    Syn57 removes TCG, TCT, TCC (from UCN) and TCA (also UCN).
    Remaining Ser codons: UCC and UCU only (both in UCN family).
    The AGY family codons (AGU, AGC) are now TARGET codons, not Ser sources.

    This means Syn57 collapses Ser to a single connected component,
    but only has 2 Ser codons left (both in the UCN box).
    """
    from codon_topo.core.encoding import codon_to_vector, hamming_distance
    from codon_topo.core.genetic_codes import STANDARD
    from codon_topo.core.homology import connected_components

    # Standard Ser codons
    std_ser = [c for c in STANDARD if STANDARD[c] == "Ser"]
    std_vecs = [codon_to_vector(c) for c in sorted(std_ser)]
    std_components = connected_components(std_vecs, 1)

    # Syn57 removes all UCN Ser codons (UCG, UCU, UCC, UCA).
    # Remaining Ser codons are AGC and AGU (the AGY family).
    # These are the TARGET codons of the recoding (TCG/TCT/TCC→AGC, TCA→AGT).
    syn57_ser_remaining = ["AGC", "AGU"]
    syn57_vecs = [codon_to_vector(c) for c in syn57_ser_remaining]
    syn57_components = connected_components(syn57_vecs, 1)

    # Hamming distance between remaining Ser codons
    hd = hamming_distance(syn57_vecs[0], syn57_vecs[1])

    return {
        "standard_ser_codons": sorted(std_ser),
        "standard_n_components": std_components,
        "syn57_ser_remaining": syn57_ser_remaining,
        "syn57_n_components": syn57_components,
        "syn57_ser_hamming": hd,
        "ser_codons_removed": ["UCG", "UCU", "UCC", "UCA"],
        "interpretation": (
            f"Standard code: {std_components} Ser components (UCN vs AGY). "
            f"Syn57: {syn57_components} component(s) with only "
            f"{len(syn57_ser_remaining)} codons remaining "
            f"(Hamming distance = {hd}). "
            f"Syn57 collapses Ser disconnection by removing 4 of 6 codons."
        ),
    }
