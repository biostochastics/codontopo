"""CodonSafe: Retrospective meta-analysis of genome recoding experiments.

Tests whether codons classified as topology-breaking in GF(2)^6 show
higher rates of recoding failure or fitness reduction than topology-
preserving codons, across published genome recoding datasets.

Tier 1 datasets (per-codon fitness resolution):
  - Napolitano et al. 2016 (PNAS): CRAM fitness, 14 loci x 64 alternatives
  - Fredens et al. 2019 (Nature): Syn61, 18,214 TCG/TCA->AGC swaps
  - Ostrov et al. 2016 (Science): 57-codon genome, 62,214 codons
  - Nyerges et al. 2024 (Nature/bioRxiv): Ec_Syn57 multi-omics
  - Robertson et al. 2025 (Nature): Ochre, 1,195 TGA->TAA
  - Frumkin et al. 2018 (PNAS): CGG recoding, ribosome profiling

Reference: Clayworth & Kornilov, CodonSafe validation plan v0.2 (2026).
"""

from codon_topo.analysis.codonsafe.models import (
    AnnotatedSwap as AnnotatedSwap,
    CodonSwapEvent as CodonSwapEvent,
    OutcomeType as OutcomeType,
    RecodingOutcome as RecodingOutcome,
    StudyId as StudyId,
    SwapModel as SwapModel,
    TopologyClassification as TopologyClassification,
)
from codon_topo.analysis.codonsafe.classify import (
    ReferenceContext as ReferenceContext,
    build_reference_context as build_reference_context,
    classify_swap_event as classify_swap_event,
)

__all__ = [
    "AnnotatedSwap",
    "CodonSwapEvent",
    "OutcomeType",
    "RecodingOutcome",
    "StudyId",
    "SwapModel",
    "TopologyClassification",
    "ReferenceContext",
    "build_reference_context",
    "classify_swap_event",
]
