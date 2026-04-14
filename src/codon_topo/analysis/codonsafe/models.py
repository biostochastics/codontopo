"""Data models for CodonSafe meta-analysis.

Defines the core dataclasses representing codon swap events, their
experimental outcomes, and their topology classifications in GF(2)^6.

Design principles:
  - Frozen dataclasses for immutability (mathematical results should not mutate).
  - Study-faithful: raw data stored as-is; classification computed separately.
  - Flexible covariates dict for study-specific confounders (mRNA structure,
    RBS disruption, codon position) without schema lock-in.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any


class StudyId(str, Enum):
    """Published genome recoding studies with per-codon fitness data."""

    NAPOLITANO_2016 = "napolitano_2016"
    FREDENS_2019 = "fredens_2019"
    OSTROV_2016 = "ostrov_2016"
    NYERGES_2024 = "nyerges_2024"
    FRUMKIN_2018 = "frumkin_2018"
    ROBERTSON_2025_OCHRE = "robertson_2025_ochre"
    ROBERTSON_2025_SYN57 = "robertson_2025_syn57"
    LAJOIE_2013 = "lajoie_2013"
    WANNIER_2018 = "wannier_2018"
    DING_2024 = "ding_2024"


class OutcomeType(str, Enum):
    """Type of fitness measurement for a recoding event."""

    BINARY_SUCCESS = "binary_success"
    FITNESS_CONTINUOUS = "fitness_continuous"


class SwapModel(str, Enum):
    """How a source->target codon swap is interpreted on the genetic code.

    REASSIGN_SOURCE_CODON:
      Construct C' identical to C except C'[source_codon] = aa(target_codon).
      One-vertex recoloring. Appropriate for nonsynonymous Napolitano CRAM data
      where all 64 alternatives are tested at a single position.

    SYNONYMOUS_REPLACEMENT:
      Source and target encode the same amino acid. The code map is unchanged,
      but local Hamming neighborhood geometry is affected. Appropriate for
      Fredens Syn61 (TCG->AGC, both Ser) and Ostrov 57-codon swaps.
    """

    REASSIGN_SOURCE_CODON = "reassign_source_codon"
    SYNONYMOUS_REPLACEMENT = "synonymous_replacement"


@dataclass(frozen=True, slots=True)
class CodonSwapEvent:
    """A single source->target codon edit attempted in a recoding study.

    Codons use RNA alphabet {A,C,G,U}. DNA-alphabet inputs (T) should be
    normalized to U in the loader before constructing this object.

    Parameters
    ----------
    study : StudyId
        Which published study this event comes from.
    event_id : str
        Unique identifier within the study (e.g., "rplN_pos5_GCU").
    source_codon : str
        Wild-type codon (RNA alphabet, 3 characters).
    target_codon : str
        Replacement codon (RNA alphabet, 3 characters).
    table_id : int
        NCBI translation table for the organism. Default 11 (bacterial).
    gene : str or None
        Gene name where the swap occurs.
    locus_tag : str or None
        Locus tag (e.g., "b0169").
    genome_pos_nt : int or None
        Genomic position of the first nucleotide (0-based).
    codon_index_in_cds : int or None
        Position of the codon within the CDS (0-based).
    is_essential_gene : bool or None
        Whether the gene is essential for viability.
    unit_id : str or None
        Grouping key for segment-level outcomes (e.g., "segment_12").
    covariates : dict
        Study-specific covariates. Common keys:
          - "mRNA_structure_deviation": float (Napolitano SRZ metric)
          - "rbs_strength_deviation": float (Napolitano SRZ metric)
          - "codon_usage_freq": float (codons per 1000)
          - "trna_gene_copies": int
    """

    study: StudyId
    event_id: str
    source_codon: str
    target_codon: str

    table_id: int = 11
    organism: str | None = None
    strain: str | None = None

    gene: str | None = None
    locus_tag: str | None = None
    genome_pos_nt: int | None = None
    codon_index_in_cds: int | None = None
    is_essential_gene: bool | None = None

    unit_id: str | None = None
    covariates: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if len(self.source_codon) != 3 or len(self.target_codon) != 3:
            raise ValueError(
                f"Codons must be 3 characters: "
                f"source={self.source_codon!r}, target={self.target_codon!r}"
            )
        valid = set("ACGU")
        for label, codon in [
            ("source", self.source_codon),
            ("target", self.target_codon),
        ]:
            if not set(codon).issubset(valid):
                raise ValueError(
                    f"{label}_codon {codon!r} contains non-RNA characters. "
                    f"Normalize T->U before constructing CodonSwapEvent."
                )


@dataclass(frozen=True, slots=True)
class RecodingOutcome:
    """Experimental fitness outcome for a codon swap event.

    Parameters
    ----------
    outcome_type : OutcomeType
        Whether this is a binary success/failure or continuous fitness.
    success : bool or None
        For BINARY_SUCCESS: True if the swap was viable.
    fitness : float or None
        For FITNESS_CONTINUOUS: quantitative fitness measure.
    fitness_unit : str or None
        Unit of the fitness measure (e.g., "doublings_per_hr", "log2FC",
        "cram_abundance_t3").
    stderr : float or None
        Standard error of the fitness estimate.
    replicate_n : int or None
        Number of biological replicates.
    timepoint : str or None
        For time-series data (CRAM), which timepoint this represents.
    """

    outcome_type: OutcomeType
    success: bool | None = None
    fitness: float | None = None
    fitness_unit: str | None = None
    stderr: float | None = None
    replicate_n: int | None = None
    timepoint: str | None = None

    def __post_init__(self) -> None:
        if self.outcome_type == OutcomeType.BINARY_SUCCESS and self.success is None:
            raise ValueError("BINARY_SUCCESS outcome requires success field")
        if (
            self.outcome_type == OutcomeType.FITNESS_CONTINUOUS
            and self.fitness is None
        ):
            raise ValueError("FITNESS_CONTINUOUS outcome requires fitness field")


@dataclass(frozen=True, slots=True)
class TopologyClassification:
    """GF(2)^6 topology and mismatch classification for a codon swap.

    This is the mathematical core of the CodonSafe analysis. Each field
    is deterministically computed from (event, encoding, translation_table).

    Parameters
    ----------
    encoding_id : int
        Index 0..23 in all_encodings().
    encoding : dict
        The base-to-bit mapping used.
    source_vec : tuple of int
        6-bit vector for the source codon.
    target_vec : tuple of int
        6-bit vector for the target codon.
    hamming : int
        Hamming distance between source and target vectors.
    source_aa : str
        Amino acid encoded by source codon under the reference table.
    target_aa : str
        Amino acid encoded by target codon under the reference table.
    is_synonymous : bool
        Whether source_aa == target_aa.
    changes_components_eps1 : bool
        Whether the swap changes the number of connected components at
        epsilon=1 for any affected amino acid.
    affected_aas : tuple of str
        Amino acids whose component structure is checked.
    delta_components_by_aa : dict
        {aa: beta0_after - beta0_before} at epsilon=1.
    crosses_component_boundary_eps1 : bool or None
        For Ser->Ser swaps: whether source and target lie in different
        connected components of the Ser codon graph at epsilon=1.
        None if the swap does not involve Ser.
    delta_edge_mismatch : dict
        {metric_name: delta_F} where delta_F = F(C') - F(C).
        Metrics: "grantham", "miyata", "polar_requirement", "kyte_doolittle".
    local_mismatch_source : dict
        {metric_name: local_cost} for the source codon position.
    local_mismatch_target : dict
        {metric_name: local_cost} for the target codon position.
    """

    encoding_id: int
    encoding: dict[str, tuple[int, int]]

    source_vec: tuple[int, ...]
    target_vec: tuple[int, ...]
    hamming: int

    source_aa: str
    target_aa: str
    is_synonymous: bool

    changes_components_eps1: bool
    affected_aas: tuple[str, ...]
    delta_components_by_aa: dict[str, int]

    crosses_component_boundary_eps1: bool | None
    delta_edge_mismatch: dict[str, float]
    local_mismatch_source: dict[str, float]
    local_mismatch_target: dict[str, float]


@dataclass(frozen=True, slots=True)
class AnnotatedSwap:
    """A fully annotated codon swap: event + outcome + topology classification.

    This is the unit of analysis for the meta-analysis. One CodonSwapEvent
    may produce multiple AnnotatedSwaps (one per encoding) for encoding
    sensitivity analysis.
    """

    event: CodonSwapEvent
    outcome: RecodingOutcome | None
    topo: TopologyClassification
