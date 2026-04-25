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
            status="tautological",
            evidence_strength="encoding_dependent",
            implications=(
                "Tautological under default encoding: all 2-fold synonymous "
                "pairs in the standard code end in {U,C} or {A,G}, so any "
                "bijection placing those pairs at identical first-bit positions "
                "yields the bit-5 pattern by construction. Holds in 16 of 24 "
                "encodings; not a discovery about the code itself."
            ),
            notes=(
                "Reclassified TAUTOLOGICAL per claim_hierarchy.py "
                "two_fold_bit_5_filtration. Manuscript Table 1 lists this as "
                "TAUTOLOGICAL; catalogue.json now agrees."
            ),
        ),
        Prediction(
            id="WS1-C2",
            claim="Four-fold degeneracy corresponds to shared 4-bit prefix exhausting GF(2)^2",
            workstream="WS1",
            status="tautological",
            evidence_strength="encoding_dependent",
            implications=(
                "Trivial under any bijection from {A,C,G,U} to GF(2)^2: "
                "four-fold degenerate codons {XYN} have third-base bits that "
                "exhaust GF(2)^2 by definition. Not a claim about the genetic "
                "code, just the definition of four-fold degeneracy."
            ),
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
            claim=(
                "Leucine is disconnected in chlorophycean mitochondrial codes "
                "(translation tables 16 and 22, epsilon=2)"
            ),
            workstream="WS1",
            status="verified",
            evidence_strength="strong",
            implications=(
                "UAG reassignment Stop->Leu creates an algal-mitochondrial "
                "signature shared by NCBI table 16 (chlorophycean mito) and "
                "table 22 (Scenedesmus obliquus mito; table 22 additionally "
                "reassigns UCA Ser->Stop). When de-duplicated to lineage-level "
                "events, the two tables collapse to a single algal-mito "
                "Leu-disconnection event."
            ),
            notes=(
                "depth_calibration.py CALIBRATION_POINTS includes both "
                "table 16 and table 22 with the same age range (500-700 Mya). "
                "Manuscript Section 3.2 Disconnection catalogue lists both."
            ),
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
            claim=(
                "Codon reassignment events show directional bias in GF(2)^6 "
                "bit positions, but signal is marginal after correcting for "
                "non-independence (de-duplicated p=0.075)"
            ),
            workstream="WS2",
            status="exploratory",
            evidence_strength="weak",
            p_value=0.075,
            implications=(
                "Original 6-bin chi-squared p=0.006 was inflated by "
                "non-independence (same codons reassigned across multiple "
                "tables). After de-duplication to 20 unique (codon, target_aa) "
                "events the marginal-position signal degrades to p=0.075; "
                "table-preserving permutation p=0.27; codon-preserving "
                "permutation p=1.0. The apparent positional bias is entirely "
                "explained by which codons are 'hot' for reassignment, not "
                "by an intrinsic preference for particular GF(2)^6 bit "
                "positions."
            ),
            notes=(
                "Reclassified EXPLORATORY/weak per claim_hierarchy.py "
                "bit_position_bias_weighted. Manuscript Table 1 lists this "
                "as EXPLORATORY with deduplicated p=0.075."
            ),
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
            evidence_strength="weak",
            p_value=1.0,
            implications="Relationship is ordinal, not linear: eps=3 events (CUG clade, 150 Mya) are younger than eps=2 events (Chlorophyceae, 600 Mya). Structural position in GF(2)^6 matters more than raw time.",
            notes="Spearman rho=0.0 on 6 calibration points. Bootstrap 95% CI spans [-1, 1]. The epsilon-age mapping is a partial order, not a monotonic function.",
        ),
        Prediction(
            id="WS4-FANO",
            claim="KRAS G12V Fano partner (His/CAC) is enriched in co-occurring mutations",
            workstream="WS4",
            status="null",
            evidence_strength="moderate",
            p_value=1.0,
            implications="Null result: Fano geometry does not predict somatic co-mutation patterns in MSK-IMPACT data. Clinical prediction track closed. Algebraic framework remains valid.",
            notes="Tested on 1670 KRAS mutations (371 G12V) from MSK-IMPACT 2017 via cBioPortal. Zero Fano-predicted AAs observed at co-mutation sites across all 6 G12 variants. Bonferroni-corrected threshold 0.00167.",
        ),
        Prediction(
            id="WS6-LANDSCAPE",
            claim="Most single-codon reassignments preserve algebraic structure (feasibility > 0.5)",
            workstream="WS6",
            status="tested",
            evidence_strength="strong",
            implications="83% of 1280 single-codon variants score >=0.8 feasibility. Zero score below 0.5. The algebraic structure is remarkably robust to perturbation.",
            notes="Scored all 64 codons x 20 possible reassignments. Feasibility = weighted composite of filtration integrity, Serine disconnection, and disconnection profile.",
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
