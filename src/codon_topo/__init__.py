"""Codon Geometry Validation & Prediction Engine."""

from codon_topo.core.encoding import (
    ALL_CODONS as ALL_CODONS,
    codon_to_vector as codon_to_vector,
    hamming_distance as hamming_distance,
)
from codon_topo.core.embedding import embed_codon as embed_codon
from codon_topo.core.fano import (
    fano_partner as fano_partner,
    is_fano_line as is_fano_line,
)
from codon_topo.core.filtration import analyze_filtration as analyze_filtration
from codon_topo.core.genetic_codes import (
    STANDARD as STANDARD,
    all_table_ids as all_table_ids,
    get_code as get_code,
)
from codon_topo.core.homology import (
    disconnection_catalogue as disconnection_catalogue,
    persistent_homology as persistent_homology,
)
from codon_topo.analysis.reassignment_db import (
    build_reassignment_db as build_reassignment_db,
    ReassignmentEvent as ReassignmentEvent,
)
from codon_topo.analysis.depth_calibration import (
    compute_correlation as compute_correlation,
    CALIBRATION_POINTS as CALIBRATION_POINTS,
)
from codon_topo.analysis.cosmic_query import (
    fano_predictions_for_kras as fano_predictions_for_kras,
    CBioPortalClient as CBioPortalClient,
    ws4_gate_decision as ws4_gate_decision,
)
from codon_topo.analysis.synbio_feasibility import (
    score_variant_code as score_variant_code,
)
from codon_topo.analysis.coloring_optimality import (
    monte_carlo_null as monte_carlo_null,
    hypercube_edge_mismatch_score as hypercube_edge_mismatch_score,
    cross_table_optimality as cross_table_optimality,
)
from codon_topo.analysis.trna_evidence import (
    trna_duplication_correlation_test as trna_duplication_correlation_test,
)
from codon_topo.reports.catalogue import (
    build_catalogue as build_catalogue,
    Prediction as Prediction,
)
from codon_topo.reports.claim_hierarchy import (
    CLAIM_HIERARCHY as CLAIM_HIERARCHY,
    ClaimStatus as ClaimStatus,
    supported_claims as supported_claims,
    hierarchy_summary_table as hierarchy_summary_table,
)

DEFAULT_SEED: int = 135325
