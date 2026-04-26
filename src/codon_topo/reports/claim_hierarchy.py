"""Claim hierarchy for the Clayworth codon-topo framework.

Single source of truth for every scientific claim's status: supported,
suggestive, exploratory, rejected, falsified, or tautological. Each claim
carries its p-value, null model, sample size, and justification.

Usage:
    from codon_topo.reports.claim_hierarchy import CLAIM_HIERARCHY, supported_claims
    for claim in supported_claims():
        print(claim.id, claim.statement)
"""

from dataclasses import dataclass
from enum import Enum


class ClaimStatus(Enum):
    """Status of a scientific claim."""

    SUPPORTED = "supported"  # Passes rigorous null; cite as finding
    SUGGESTIVE = "suggestive"  # Trend-level; cite as supporting evidence
    EXPLORATORY = "exploratory"  # Borderline; hypothesis-generating
    REJECTED = "rejected"  # Falsified or pre-rejected in literature
    FALSIFIED = "falsified"  # Directly contradicted by data
    TAUTOLOGICAL = "tautological"  # True by construction; not a discovery


@dataclass(frozen=True)
class Claim:
    """A scientific claim with its current status."""

    id: str
    statement: str
    status: ClaimStatus
    evidence_p_value: float | None
    null_model: str | None
    sample_size: int | None
    justification: str
    citation_key: str | None = None

    @property
    def is_publishable(self) -> bool:
        """Whether the claim should appear in the paper (with appropriate framing)."""
        return self.status in (
            ClaimStatus.SUPPORTED,
            ClaimStatus.SUGGESTIVE,
            ClaimStatus.EXPLORATORY,
        )


CLAIM_HIERARCHY: list[Claim] = [
    # ============================================================
    # SUPPORTED — cite as primary findings
    # ============================================================
    Claim(
        id="hypercube_coloring_optimality",
        statement=(
            "The standard genetic code is significantly error-minimizing "
            "under four independent physicochemical distance metrics: "
            "Grantham (p=0.006), Miyata (p<0.001), Woese polar requirement "
            "(p=0.003), and Kyte-Doolittle hydropathy (p=0.001)"
        ),
        status=ClaimStatus.SUPPORTED,
        evidence_p_value=0.006,  # max p across all 4 metrics
        null_model="freeland_hurst_block_preserving_x4_metrics",
        sample_size=10_000,
        justification=(
            "Cross-metric sensitivity analysis under block-preserving null. "
            "Freeland & Hurst (1998) used polar requirement only; we confirm "
            "optimality holds across 4 independent metrics (Grantham, Miyata, "
            "polar requirement, Kyte-Doolittle). Stop penalty sensitivity "
            "tested (0/150/215/300): immaterial. Establishes code error-"
            "minimization as a general structural property."
        ),
        citation_key="FreelandHurst1998",
    ),
    Claim(
        id="per_table_optimality_preservation",
        statement=(
            "26 of 27 NCBI translation tables remain in the top 5% of "
            "their own block-preserving null for Grantham edge-mismatch "
            "(BH-FDR corrected); only translation table 3 (yeast mito) "
            "exceeds the threshold"
        ),
        status=ClaimStatus.SUPPORTED,
        evidence_p_value=None,  # Aggregate: count from manuscript_stats.json per_table.n_significant_bh
        null_model="freeland_hurst_block_preserving_per_table",
        sample_size=10_000,
        justification=(
            "Each variant code tested against its own block-preserving null "
            "(n=10,000 per table, common-seed design, BH-FDR corrected). "
            "Across 27 NCBI tables (codes 1-6, 9-16, 21-33; codes 7, 8, "
            "17-20 deprecated): 26 significant at adjusted p<0.05. Mean "
            "quantile ~1.4%. Only table 3 (yeast mitochondrial, 6 codon "
            "reassignments including the topology-breaking CUN-Thr) "
            "exceeds the threshold. Tables 1/11 and 27/28 share identical "
            "sense-codon mappings; the per-table null seeds differ so "
            "rendered quantiles differ slightly even when the codes are "
            "sense-codon-identical. A standard-code-proximity audit "
            "(Supplement) confirms the variant-code quantile is not driven "
            "purely by closeness to the standard code."
        ),
    ),
    Claim(
        id="optimality_rho_robustness",
        statement=(
            "Coloring optimality is robust across all diagonal-edge weights "
            "rho in [0,1], including the full Hamming graph H(3,4) = "
            "K_4 □ K_4 □ K_4 of single-nucleotide substitutions (p<0.01 "
            "at all rho values)"
        ),
        status=ClaimStatus.SUPPORTED,
        evidence_p_value=0.0061,  # max p across all rho values (rho=0)
        null_model="freeland_hurst_weighted_edges",
        sample_size=10_000,
        justification=(
            "Adding within-nucleotide distance-2 'diagonal' edges (completing "
            "H(3,4) from Q_6) preserves and strengthens optimality. At rho=1 "
            "(all 288 single-nucleotide edges equally weighted), p=3e-4. "
            "Common-seed design with n=10,000 block-preserving null samples "
            "per rho value. Exact p-values: rho=0 p=0.0061, rho=0.25 "
            "p=0.0023, rho=0.5 p=0.0007, rho=0.75 p=0.0003, rho=1 p=0.0003. "
            "Confirms the optimality is not an artifact of the Q_6 subgraph "
            "omitting the 96 within-nucleotide diagonal edges (~1/3 of "
            "single-nucleotide mutations). The rho parameter is a "
            "diagonal-edge weight, NOT a transition/transversion weight: "
            "16 of 24 encodings mix transitions and transversions equally "
            "among Hamming-1 edges."
        ),
    ),
    Claim(
        id="topology_avoidance_depletion",
        statement=(
            "Natural codon reassignments are depleted for topology-breaking "
            "changes — moves that fragment an amino acid's codon family — "
            "at approximately 21% of observed events vs 66-73% of the "
            "candidate landscape, robust to adjacency definition (Q_6 vs "
            "the full Hamming graph H(3,4)) and clade exclusion"
        ),
        status=ClaimStatus.SUPPORTED,
        evidence_p_value=1e-4,  # Permutation p across both adjacency definitions
        null_model="hypergeometric + table_preserving_permutation_10k",
        sample_size=28,
        justification=(
            "Candidate universe N=1,280 single-codon relabelings "
            "((x,y): x in 64 codons, y in 20 amino acids + Stop, y != C(x)). "
            "Q_6 adjacency: 931 of 1,280 (72.7%) candidates create a new "
            "disconnection in a previously connected family; observed 6 of "
            "28 events (21.4%); hypergeometric p=1.6e-8, permutation "
            "p<=1e-4, RR=0.29 (95% CI [0.16, 0.66]). H(3,4) adjacency "
            "(encoding-independent): 846 of 1,280 (66.1%); same 6 of 28 "
            "events; hypergeometric p=1.3e-6, permutation p<=1e-4, RR=0.32 "
            "(95% CI [0.16, 0.66]). Phylogenetic clade-exclusion sensitivity "
            "(Sengupta et al. 2007): all 7 exclusion regimes remain "
            "significant; excluding yeast mito (Table 3) actually "
            "strengthens depletion (rate drops 21.4%->8.3%, p drops "
            "1.6e-8->3.6e-11). Yeast mito accounts for 4 of 6 "
            "topology-breaking events and is the marginal outlier across "
            "three independent analyses (per-table optimality, "
            "topology-avoidance counts, conditional-logit case ranks)."
        ),
    ),
    # ============================================================
    # SUGGESTIVE — cite as supporting evidence, not primary
    # ============================================================
    Claim(
        id="trna_enrichment_reassigned_aa",
        statement=(
            "Organisms with variant genetic codes show elevated tRNA gene "
            "copy numbers for the reassigned amino acid, verified by "
            "tRNAscan-SE 2.0.12 on 18 NCBI genome assemblies across "
            "5 variant genetic codes (24 pairings, 6 independent)"
        ),
        status=ClaimStatus.SUGGESTIVE,
        evidence_p_value=0.045,  # MIS worst-case Stouffer (conservative primary)
        null_model="fisher_exact_stouffer_MIS_worst_case",
        sample_size=6,
        justification=(
            "Conservative primary: MIS (maximal independent set) worst-case "
            "Stouffer p=0.045. Bron-Kerbosch enumeration of all 2 MIS from "
            "conflict graph; both significant at p<0.05. Eliminates greedy-"
            "selection bias. All-pairings (n=24): p=1.7e-7. Independent "
            "greedy (n=6): p=0.017. 18 tRNAscan-SE verified assemblies: "
            "5 Table 6 (Gln), 6 Table 10 (Cys), 1 Table 15 (Trp), "
            "1 Table 31 (Trp/Glu), 2 Table 4 (Trp/bacterial), plus "
            "3 standard-code controls. Covers 5 variant genetic codes "
            "across Alveolata, Opisthokonta, Excavata, and Mollicutes."
        ),
        citation_key="Su2011PMC3113583",
    ),
    # ============================================================
    # EXPLORATORY — include with caveats, frame as hypothesis
    # ============================================================
    Claim(
        id="bit_position_bias_weighted",
        statement=(
            "Codon reassignment bit-flip distribution shows positional skew "
            "in GF(2)^6 coordinates, but signal is marginal after correcting "
            "for non-independence (de-duplicated p=0.075)"
        ),
        status=ClaimStatus.EXPLORATORY,
        evidence_p_value=0.075,  # De-duplicated unique events
        null_model="chi_square_deduplicated",
        sample_size=20,
        justification=(
            "Original 6-bin p=0.006 inflated by non-independence (same codons "
            "across multiple tables). After de-duplication to 20 unique "
            "(codon, target_aa) events: p=0.075. Table-preserving permutation "
            "p=0.27; codon-preserving p=1.0. Signal is entirely explained by "
            "which codons are hot for reassignment, not by positional preference."
        ),
        citation_key=None,
    ),
    Claim(
        id="mechanism_boundary_conditions",
        statement=(
            "tRNA gene duplication accompanies codon reassignment in large "
            "nuclear genomes (ciliates, yeasts) but not in streamlined "
            "genomes (Blastocrithidia: anticodon stem shortening; "
            "Mycoplasma: anticodon modification)"
        ),
        status=ClaimStatus.EXPLORATORY,
        evidence_p_value=None,
        null_model=None,
        sample_size=7,
        justification=(
            "Three-tier pattern: (1) Ciliates (T. thermophila 54 Gln tRNAs, "
            "39 suppressor) — massive gene duplication; (2) Blastocrithidia "
            "nonstop P57 (GCA_028554745.1, 24.7 Mb genome, 68 tRNAs: "
            "UGA→Trp via tRNA-Trp(CCA) anticodon stem shortening 5bp→4bp "
            "+ eRF1 Ser74Gly mutation; UAA/UAG→Glu via 2 new suppressor "
            "tRNAs; Kachale 2023 Nature 613:751-758); (3) Mycoplasma "
            "(1 Trp tRNA reads both UGA+UGG — anticodon modification). "
            "Consistent with multiple evolutionary routes to reassignment. "
            "Descriptive, not inferential."
        ),
    ),
    Claim(
        id="atchley_f3_serine_convergence",
        statement=(
            "Serine's extreme Atchley Factor 3 score (F3 = −4.760, most "
            "extreme of 20 amino acids, 2.24 SD below mean) converges with "
            "the GF(2)^6 topological disconnection: F3 captures the mismatch "
            "between Serine's small physicochemical footprint and its "
            "disproportionate codon diversity, which the geometric framework "
            "identifies as maximal inter-family Hamming distance among "
            "6-codon amino acids (4 vs 1 for Leu and Arg)"
        ),
        status=ClaimStatus.EXPLORATORY,
        evidence_p_value=None,
        null_model=None,
        sample_size=None,
        justification=(
            "Atchley et al. (PNAS 2005) F3 is a factor-analytic composite "
            "of molecular size and codon diversity derived from 54 amino acid "
            "attributes. The convergence is real but not 'independent "
            "confirmation' in the strict sense: F3 partially reflects codon "
            "usage, so both measurements share a proximate cause (Serine's "
            "codon distribution). Correctly framed as: the geometric "
            "framework provides a structural explanation for why Serine is "
            "an F3 outlier — two complementary views of the same anomaly, "
            "not two independent detections. Encoding-dependence caveat (distance-4 holds for "
            "8/24 encodings; disconnection at ε=1 is universal) applies."
        ),
        citation_key="Atchley2005",
    ),
    Claim(
        id="variant_code_disconnection_catalogue",
        statement=(
            "Systematic survey of Hamming-graph connectivity across 27 NCBI "
            "translation tables identifies 4 lineage-collapsed variant-code "
            "amino-acid disconnections at epsilon=1 under the default "
            "encoding: Thr in yeast mitochondrial code (table 3), Leu in "
            "chlorophycean mitochondrial codes (tables 16 and 22, collapsed "
            "to a single algal-mito event), Ala in Pachysolen tannophilus "
            "nuclear code (table 26), and tripartite Ser in Candida-clade "
            "alternative yeast nuclear code (table 12)"
        ),
        status=ClaimStatus.EXPLORATORY,
        evidence_p_value=None,  # Descriptive catalogue, not inferential
        null_model=None,
        sample_size=27,
        justification=(
            "Descriptive survey across all 27 current NCBI translation tables "
            "(codes 1-6, 9-16, 21-33; codes 7, 8, 17-20 deprecated). "
            "Tables 16 and 22 (chlorophycean and Scenedesmus obliquus mito "
            "respectively) both reassign UAG-Leu and produce equivalent "
            "epsilon=2 reconnection profiles for Leu; lineage-collapsed they "
            "are a single event. Table 22 additionally reassigns UCA Ser->Stop. "
            "Links to literature: yeast mito CUN-Thr (Miranda 2006), "
            "chlorophycean UAG-Leu (Knaap 2002), Candida CUG-Ser "
            "(Santos 1996), Pachysolen CUG-Ala (Muhlhausen 2018). A "
            "separate filtration-only deviation (translation table 32, "
            "Balanophoraceae plastid, UAG-Trp at bit 3) is documented but "
            "does not introduce a new disconnection."
        ),
    ),
    # ============================================================
    # FALSIFIED — tested against data, failed
    # ============================================================
    Claim(
        id="kras_fano_clinical_prediction",
        statement=(
            "XOR 'Fano' relationships in GF(2)^6 predict enrichment of "
            "specific amino acids at KRAS G12 co-mutation sites"
        ),
        status=ClaimStatus.FALSIFIED,
        evidence_p_value=1.0,  # Fisher's exact with Bonferroni
        null_model="fisher_exact_bonferroni_corrected",
        sample_size=1670,
        justification=(
            "Tested against MSK-IMPACT 2017 (cBioPortal, 1670 KRAS "
            "mutations). Zero enrichment across all 6 G12 variants. "
            "No biological mechanism would make DNA polymerase respect "
            "XOR structure. Report as clean negative for literature record."
        ),
    ),
    # ============================================================
    # REJECTED — incorrect math or pre-rejected in literature
    # ============================================================
    Claim(
        id="serine_min_distance_4_invariant",
        statement=(
            "Serine's minimum UCN-AGY Hamming distance = 4 holds "
            "invariantly across all 24 base-to-bit encodings"
        ),
        status=ClaimStatus.REJECTED,
        evidence_p_value=None,
        null_model="counterexample_verification",
        sample_size=24,
        justification=(
            "Counterexample: phi(U)=00, phi(C)=11, phi(A)=01, "
            "phi(G)=10 gives min distance 2. Verified in repository code. "
            "Across 24 encodings: 16 give min=2, 8 give min=4. Distance-4 "
            "obtains exactly when U⊕A = (1,1) and C⊕G = (1,1), i.e. both "
            "nucleotide pairs distinguishing the first two positions of UCN "
            "vs AGY are encoded at maximal Hamming distance. The only "
            "universal (encoding-invariant) Serine result is disconnection "
            "at ε=1. Even so, Serine's inter-family distance (4 under the "
            "default encoding, ≥2 under all encodings) is uniquely extreme "
            "among 6-codon amino acids: Leu (CUN-UUA/UUG) and Arg "
            "(CGN-AGA/AGG) both have min inter-family distance = 1."
        ),
    ),
    Claim(
        id="psl_2_7_symmetry",
        statement="PSL(2,7) is the fundamental symmetry group of the genetic code",
        status=ClaimStatus.REJECTED,
        evidence_p_value=None,
        null_model=None,
        sample_size=None,
        justification=(
            "Pre-rejected by Antoneli & Forger (2011, Math Comp Model 53). "
            "PSL(2,7) has no 64-dimensional irreducible representation "
            "(its irreps are dim 1, 3, 6, 7, 8). Cannot accommodate the "
            "codon multiplet structure."
        ),
        citation_key="AntoneliForger2011",
    ),
    Claim(
        id="holomorphic_embedding",
        statement=(
            "The coordinate-wise map GF(2)^6 -> C^3 sending base-pairs to "
            "fourth roots of unity is a holomorphic embedding extending a "
            "character of GF(8)*"
        ),
        status=ClaimStatus.REJECTED,
        evidence_p_value=None,
        null_model=None,
        sample_size=None,
        justification=(
            "Not holomorphic: domain is a finite discrete set. "
            "Not a character: GF(2)^2 has exponent 2, so characters to C* "
            "must take values in {+1, -1}; but (0,1) -> i has order 4 and "
            "fails chi(x+x) = chi(x)^2 (i^2 = -1 != 1). The map is simply "
            "a coordinate-wise bijection to the fourth roots of unity."
        ),
    ),
    # ============================================================
    # TAUTOLOGICAL — true by construction, cite as 'expected'
    # ============================================================
    Claim(
        id="two_fold_bit_5_filtration",
        statement=(
            "All 2-fold synonymous codon pairs in the standard code "
            "differ at exactly bit position 5 under the default encoding"
        ),
        status=ClaimStatus.TAUTOLOGICAL,
        evidence_p_value=None,
        null_model=None,
        sample_size=None,
        justification=(
            "Forced by the Gray-like encoding. Two-fold degenerate codons "
            "in the standard code always end in {U,C} or {A,G}; the default "
            "encoding places these pairs at identical first-bit positions. "
            "Not a discovery about the code; a property of the encoding "
            "choice. Holds in 16/24 encodings."
        ),
    ),
    Claim(
        id="four_fold_prefix_filtration",
        statement=(
            "All 4-fold synonymous codon groups share a 4-bit prefix "
            "with 2-bit suffixes exhausting GF(2)^2"
        ),
        status=ClaimStatus.TAUTOLOGICAL,
        evidence_p_value=None,
        null_model=None,
        sample_size=None,
        justification=(
            "Trivial under any bijection from {A,C,G,U} to GF(2)^2. Four-"
            "fold degenerate amino acids have codons {XYN}; under any "
            "bijection the third-base bits exhaust GF(2)^2 by definition. "
            "Not a claim about the genetic code, just the definition of "
            "four-fold degeneracy."
        ),
    ),
]


def claims_by_status(status: ClaimStatus) -> list[Claim]:
    """Return all claims with a given status."""
    return [c for c in CLAIM_HIERARCHY if c.status == status]


def supported_claims() -> list[Claim]:
    return claims_by_status(ClaimStatus.SUPPORTED)


def suggestive_claims() -> list[Claim]:
    return claims_by_status(ClaimStatus.SUGGESTIVE)


def exploratory_claims() -> list[Claim]:
    return claims_by_status(ClaimStatus.EXPLORATORY)


def rejected_claims() -> list[Claim]:
    """Claims that failed (falsified, rejected, or tautological)."""
    return [
        c
        for c in CLAIM_HIERARCHY
        if c.status
        in (ClaimStatus.REJECTED, ClaimStatus.FALSIFIED, ClaimStatus.TAUTOLOGICAL)
    ]


def publishable_claims() -> list[Claim]:
    return [c for c in CLAIM_HIERARCHY if c.is_publishable]


def hierarchy_summary_table() -> str:
    """Render a human-readable summary of all claims grouped by status."""
    lines = []
    for status in ClaimStatus:
        claims = claims_by_status(status)
        if not claims:
            continue
        lines.append(f"\n## {status.value.upper()} ({len(claims)} claims)\n")
        for c in claims:
            p_str = f" p={c.evidence_p_value}" if c.evidence_p_value is not None else ""
            n_str = f" n={c.sample_size}" if c.sample_size is not None else ""
            lines.append(f"- **{c.id}**:{p_str}{n_str}")
            lines.append(f"  {c.statement}")
    return "\n".join(lines)


def abstract_ready_paragraph() -> str:
    """Generate an abstract-ready paragraph summarizing the hierarchy."""
    sup = supported_claims()
    sug = suggestive_claims()
    exp = exploratory_claims()
    rej = rejected_claims()
    return (
        f"This work rigorously evaluates {len(CLAIM_HIERARCHY)} claims "
        f"about the algebraic structure of the genetic code. "
        f"{len(sup)} claim(s) are supported under strict null models, "
        f"{len(sug)} claim(s) are suggestive (trend-level), "
        f"{len(exp)} claim(s) remain exploratory, and "
        f"{len(rej)} claim(s) were rejected as tautological, falsified, "
        f"or pre-rejected in prior literature."
    )
