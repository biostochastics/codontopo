"""Loader for Napolitano et al. 2016 (PNAS) CRAM fitness data.

Reference: Napolitano MG et al. "Emergent rules for codon choice
elucidated by editing rare arginine codons in Escherichia coli."
PNAS 113(38):E5588-E5597. doi:10.1073/pnas.1605856113

Data sources:
  - Dataset S6 (pnas.1605856113.sd06.xlsx): CRAM fitness tracking
    for all 64 codon alternatives at 14 AGR positions in essential genes.
  - Dataset S1 (pnas.1605856113.sd01.xlsx): 123 AGR codon positions,
    recalcitrance classification, SRZ metrics.

Key data structure:
  - 14 positions x 64 alternatives = 896 datapoints
  - Each row: gene, codon position, alternative codon, CRAM abundance
    at timepoints T0, T1, T2, T3
  - Success/failure determined by depletion: codons outside SRZ
    are rapidly depleted from the growing population.

Note: Download supplementary data from PNAS and place in:
  data/codonsafe/napolitano2016/
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from codon_topo.analysis.codonsafe.models import (
    CodonSwapEvent,
    OutcomeType,
    RecodingOutcome,
    StudyId,
)
from codon_topo.analysis.codonsafe.normalize import normalize_codon

if TYPE_CHECKING:
    pass

# The 14 AGR positions tested by CRAM (from Napolitano Fig. 6 / Dataset S6)
# Format: (gene, codon_index_in_cds, wild_type_codon_dna)
CRAM_POSITIONS: list[tuple[str, int, str]] = [
    # These will be populated from Dataset S6 upon data download.
    # Placeholder structure based on paper text.
]

DATA_DIR = (
    Path(__file__).parent.parent.parent.parent / "data" / "codonsafe" / "napolitano2016"
)


def load_cram_data() -> list[tuple[CodonSwapEvent, RecodingOutcome]]:
    """Load CRAM fitness data from Dataset S6.

    Returns list of (event, outcome) pairs for all 14x64 codon alternatives.

    Raises FileNotFoundError if the supplementary data has not been downloaded.
    """
    sd06_path = DATA_DIR / "pnas.1605856113.sd06.xlsx"
    if not sd06_path.exists():
        raise FileNotFoundError(
            f"Napolitano 2016 Dataset S6 not found at {sd06_path}. "
            f"Download from https://www.pnas.org/doi/10.1073/pnas.1605856113 "
            f"and place in {DATA_DIR}/"
        )

    import pandas as pd

    df = pd.read_excel(sd06_path)

    results: list[tuple[CodonSwapEvent, RecodingOutcome]] = []

    # Iterate over each position x alternative
    # Column mapping will be finalized after inspecting the actual Excel structure.
    # Expected columns based on paper: Gene, Position, Codon, T0, T1, T2, T3
    for _, row in df.iterrows():
        gene = str(row.get("Gene", row.get("gene", "")))
        pos = int(row.get("Position", row.get("position", 0)))

        source_dna = str(row.get("WT_Codon", row.get("wt_codon", "NNN")))
        target_dna = str(row.get("Codon", row.get("codon", "NNN")))

        try:
            source_rna = normalize_codon(source_dna)
            target_rna = normalize_codon(target_dna)
        except ValueError:
            continue

        event_id = f"napolitano_{gene}_{pos}_{target_rna}"

        # Extract fitness (CRAM abundance at intermediate timepoint)
        fitness_val = row.get("T2", row.get("t2", None))
        if fitness_val is not None:
            fitness_val = float(fitness_val)

        # Determine success: codons within SRZ are viable
        srz_ok = row.get("SRZ", row.get("in_srz", None))

        covariates: dict = {}
        for cov_key in [
            "mRNA_structure_deviation",
            "rbs_strength_deviation",
            "mRNA_fold_energy",
            "RBS_strength",
        ]:
            if cov_key in row and pd.notna(row[cov_key]):
                covariates[cov_key] = float(row[cov_key])

        event = CodonSwapEvent(
            study=StudyId.NAPOLITANO_2016,
            event_id=event_id,
            source_codon=source_rna,
            target_codon=target_rna,
            table_id=11,
            organism="E_coli",
            strain="MDS42",
            gene=gene,
            codon_index_in_cds=pos,
            is_essential_gene=True,
            covariates=covariates,
        )

        if fitness_val is not None:
            outcome = RecodingOutcome(
                outcome_type=OutcomeType.FITNESS_CONTINUOUS,
                fitness=fitness_val,
                fitness_unit="cram_abundance_t2",
            )
        elif srz_ok is not None:
            outcome = RecodingOutcome(
                outcome_type=OutcomeType.BINARY_SUCCESS,
                success=bool(srz_ok),
            )
        else:
            # Skip events with no outcome data to avoid positive bias
            continue

        results.append((event, outcome))

    return results


def load_recalcitrant_codons() -> list[tuple[CodonSwapEvent, RecodingOutcome]]:
    """Load the 13 recalcitrant AGR codons from Dataset S1.

    These are the AGR->CGU swaps that failed (lethal). Each gets
    outcome success=False.
    """
    sd01_path = DATA_DIR / "pnas.1605856113.sd01.xlsx"
    if not sd01_path.exists():
        raise FileNotFoundError(
            f"Napolitano 2016 Dataset S1 not found at {sd01_path}. "
            f"Download from https://www.pnas.org/doi/10.1073/pnas.1605856113"
        )

    import pandas as pd

    df = pd.read_excel(sd01_path)
    results: list[tuple[CodonSwapEvent, RecodingOutcome]] = []

    for _, row in df.iterrows():
        gene = str(row.get("Gene", row.get("gene", "")))
        pos = int(row.get("Position", row.get("position", 0)))
        source_dna = str(row.get("Codon", row.get("codon", "NNN")))
        recalcitrant = row.get("Recalcitrant", row.get("recalcitrant", False))

        try:
            source_rna = normalize_codon(source_dna)
        except ValueError:
            continue

        # Default replacement was CGU
        target_rna = normalize_codon("CGT")

        event = CodonSwapEvent(
            study=StudyId.NAPOLITANO_2016,
            event_id=f"napolitano_recalc_{gene}_{pos}",
            source_codon=source_rna,
            target_codon=target_rna,
            table_id=11,
            organism="E_coli",
            strain="MDS42",
            gene=gene,
            codon_index_in_cds=pos,
            is_essential_gene=True,
        )

        outcome = RecodingOutcome(
            outcome_type=OutcomeType.BINARY_SUCCESS,
            success=not bool(recalcitrant),
        )

        results.append((event, outcome))

    return results
