"""Prediction catalogue: structured synthesis of WS1-WS6 findings.

Assembles all validated claims, statistical results, and predictions
into a structured catalogue with evidence grading and implications.
"""

from collections import Counter
from dataclasses import dataclass


@dataclass(frozen=True)
class Prediction:
    """A single prediction or validated claim."""

    id: str
    claim: str
    workstream: str
    status: str  # 'verified', 'tested', 'pending', 'null'
    evidence_strength: str  # 'very_strong', 'strong', 'moderate', 'weak', 'untested'
    p_value: float | None = None
    implications: str = ""
    notes: str = ""


def build_catalogue() -> list[Prediction]:
    """Build the complete prediction catalogue from WS1-WS6 results."""
    return [
        Prediction(
            id="WS1-C1",
            claim="Two-fold degeneracy is uniformly explained by bit-5 difference in GF(2)^6",
            workstream="WS1",
            status="verified",
            evidence_strength="very_strong",
            implications="All 9 two-fold AAs show identical algebraic structure across all 25 NCBI tables",
            notes="100% pass rate, zero exceptions, verified by test_regression.py",
        ),
        Prediction(
            id="WS1-C2",
            claim="Four-fold degeneracy corresponds to shared 4-bit prefix exhausting GF(2)^2",
            workstream="WS1",
            status="verified",
            evidence_strength="very_strong",
            implications="Breaks only on stop->AA reassignments (structural necessity)",
        ),
        Prediction(
            id="WS1-C3",
            claim="Serine is uniquely disconnected at epsilon=1 in all known genetic codes",
            workstream="WS1",
            status="verified",
            evidence_strength="very_strong",
            implications="Universal invariant: UCN-AGY split with min Hamming distance 4, reconnects at epsilon=4",
        ),
        Prediction(
            id="WS1-NULL-A",
            claim="Serine uniqueness is statistically unlikely under random codon assignment",
            workstream="WS1",
            status="verified",
            evidence_strength="strong",
            p_value=0.05,
            implications="Null Model A rejects random assignment hypothesis",
        ),
        Prediction(
            id="WS1-NULL-C",
            claim="Serine disconnection is encoding-invariant (all 24 base-to-bit mappings)",
            workstream="WS1",
            status="verified",
            evidence_strength="very_strong",
            implications="Algebraic structure is not an artifact of encoding choice",
        ),
        Prediction(
            id="WS1-NOVEL-1",
            claim="Threonine is disconnected in yeast mitochondrial code (Table 3, epsilon=2)",
            workstream="WS1",
            status="verified",
            evidence_strength="strong",
            implications="CUN block reassignment Leu->Thr creates lineage-specific topological signature",
        ),
        Prediction(
            id="WS1-NOVEL-2",
            claim="Leucine is disconnected in chlorophycean mito (Table 16, epsilon=2)",
            workstream="WS1",
            status="verified",
            evidence_strength="strong",
            implications="UAG reassignment Stop->Leu creates algal-specific signature",
        ),
        Prediction(
            id="WS1-NOVEL-3",
            claim="Alanine is disconnected in Pachysolen nuclear code (Table 26, epsilon=3)",
            workstream="WS1",
            status="verified",
            evidence_strength="strong",
            implications="CUG reassignment Leu->Ala in CUG clade",
        ),
        Prediction(
            id="WS1-NOVEL-4",
            claim="Serine has three components in alternative yeast nuclear code (Table 12)",
            workstream="WS1",
            status="verified",
            evidence_strength="strong",
            implications="Unique 3-component topology from CUG->Ser reassignment",
        ),
        Prediction(
            id="WS2-DIR",
            claim="Codon reassignment events show directional bias in GF(2)^6 bit positions",
            workstream="WS2",
            status="tested",
            evidence_strength="moderate",
            implications="Certain bit positions are more susceptible to reassignment than others",
        ),
        Prediction(
            id="WS2-PATH",
            claim="Reassignment Hamming paths have characteristic length distribution",
            workstream="WS2",
            status="tested",
            evidence_strength="moderate",
            implications="Mean Hamming distance constrains the 'reach' of evolutionary reassignment",
        ),
        Prediction(
            id="WS3-CORR",
            claim="Disconnection epsilon correlates positively with evolutionary age",
            workstream="WS3",
            status="tested",
            evidence_strength="moderate",
            implications="Higher epsilon = deeper evolutionary origin; ordinal dating framework",
        ),
        Prediction(
            id="WS4-FANO",
            claim="KRAS G12V Fano partner (His/CAC) is enriched in co-occurring mutations",
            workstream="WS4",
            status="pending",
            evidence_strength="untested",
            implications="If confirmed: Fano geometry predicts somatic co-mutation patterns. If null: clinical track closed, math stands.",
        ),
        Prediction(
            id="WS6-LANDSCAPE",
            claim="Most single-codon reassignments preserve algebraic structure (feasibility > 0.5)",
            workstream="WS6",
            status="tested",
            evidence_strength="moderate",
            implications="Synthetic code engineering has a large viable design space",
        ),
        Prediction(
            id="WS6-SERINE",
            claim="Serine disconnection is robust to most single-codon perturbations",
            workstream="WS6",
            status="tested",
            evidence_strength="strong",
            implications="UCN-AGY gap is structurally stable — difficult to bridge by single reassignment",
        ),
    ]


def catalogue_summary() -> dict:
    """Summary statistics for the prediction catalogue."""
    cat = build_catalogue()
    by_ws = Counter(p.workstream for p in cat)
    by_status = Counter(p.status for p in cat)
    by_strength = Counter(p.evidence_strength for p in cat)

    return {
        "total_predictions": len(cat),
        "by_workstream": dict(by_ws),
        "by_status": dict(by_status),
        "by_strength": dict(by_strength),
    }
