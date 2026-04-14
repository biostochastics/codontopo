"""Topology and mismatch classification engine for CodonSafe.

Classifies each codon swap event by its effect on GF(2)^6 topology:
  1. Does the swap change connected-component structure at epsilon=1?
  2. For Ser swaps, does it cross the disconnection boundary?
  3. What is the change in edge-mismatch score (delta_F) under 4 metrics?
  4. What is the local mismatch cost at source and target positions?

Performance note: For the CRAM dataset (~896 events x 24 encodings x 4 metrics
= ~86k classifications), the fast delta_F computation avoids recomputing the
full 192-edge sum. Instead, only the 6 edges incident to the recolored codon
are recalculated.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from itertools import combinations

from codon_topo.analysis.codonsafe.models import (
    CodonSwapEvent,
    SwapModel,
    TopologyClassification,
)
from codon_topo.analysis.coloring_optimality import (
    AVAILABLE_METRICS,
    get_distance_func,
)
from codon_topo.core.encoding import (
    ALL_CODONS,
    all_encodings,
    codon_to_vector,
    hamming_distance,
)
from codon_topo.core.genetic_codes import get_code
from codon_topo.core.homology import connected_components


@dataclass(frozen=True)
class ReferenceContext:
    """Precomputed context for a (table_id, encoding) pair.

    Caching this avoids recomputing codon vectors and neighbor lists
    for every event in a study.

    Parameters
    ----------
    table_id : int
        NCBI translation table.
    encoding_id : int
        Index 0..23 in all_encodings().
    encoding : dict
        Base-to-bit mapping.
    code : dict
        Codon->amino acid mapping for the translation table.
    vectors : dict
        Codon->6-bit vector mapping under this encoding.
    neighbors : dict
        Codon->list of Hamming-1 neighbor codons under this encoding.
    """

    table_id: int
    encoding_id: int
    encoding: dict[str, tuple[int, int]]
    code: dict[str, str]
    vectors: dict[str, tuple[int, ...]]
    neighbors: dict[str, list[str]]


def build_reference_context(
    *,
    table_id: int = 1,
    encoding_id: int = 0,
    encoding: dict[str, tuple[int, int]] | None = None,
) -> ReferenceContext:
    """Build a ReferenceContext for a (table_id, encoding) pair.

    Parameters
    ----------
    table_id : int
        NCBI translation table (default: 1, standard code).
    encoding_id : int
        Index in all_encodings(). Ignored if encoding is provided directly.
    encoding : dict or None
        Base-to-bit mapping. If None, uses all_encodings()[encoding_id].
    """
    if encoding is None:
        enc = all_encodings()[encoding_id]
    else:
        enc = encoding

    code = get_code(table_id)
    vectors = {c: codon_to_vector(c, enc) for c in ALL_CODONS}

    # Build Hamming-1 neighbor map
    neighbors: dict[str, list[str]] = defaultdict(list)
    for c1, c2 in combinations(ALL_CODONS, 2):
        if hamming_distance(vectors[c1], vectors[c2]) == 1:
            neighbors[c1].append(c2)
            neighbors[c2].append(c1)

    return ReferenceContext(
        table_id=table_id,
        encoding_id=encoding_id,
        encoding=enc,
        code=code,
        vectors=dict(vectors),
        neighbors=dict(neighbors),
    )


def _beta0_eps1_for_aa(
    code: dict[str, str],
    aa: str,
    vectors: dict[str, tuple[int, ...]],
) -> int:
    """Count connected components at epsilon=1 for codons encoding `aa`."""
    aa_codons = [c for c in ALL_CODONS if code.get(c) == aa]
    if len(aa_codons) <= 1:
        return len(aa_codons)
    aa_vecs = [vectors[c] for c in aa_codons]
    return connected_components(aa_vecs, epsilon=1)


def _crosses_component_boundary(
    code: dict[str, str],
    vectors: dict[str, tuple[int, ...]],
    source_codon: str,
    target_codon: str,
) -> bool | None:
    """Check if a synonymous swap crosses a disconnection boundary at eps=1.

    Works for any amino acid (Ser, Leu, Arg, etc.), not just Serine.
    Returns None if source and target encode different amino acids.
    Returns True if source and target lie in different connected components
    of the amino acid's codon graph at epsilon=1.
    """
    source_aa = code.get(source_codon)
    target_aa = code.get(target_codon)
    if source_aa != target_aa or source_aa is None or source_aa == "Stop":
        return None

    aa = source_aa
    aa_codons = sorted(c for c in ALL_CODONS if code.get(c) == aa)
    if len(aa_codons) < 2:
        return None
    aa_vecs = [vectors[c] for c in aa_codons]

    # Find which component each codon belongs to using union-find
    n = len(aa_codons)
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x: int, y: int) -> None:
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py

    for i in range(n):
        for j in range(i + 1, n):
            if hamming_distance(aa_vecs[i], aa_vecs[j]) <= 1:
                union(i, j)

    try:
        src_idx = aa_codons.index(source_codon)
        tgt_idx = aa_codons.index(target_codon)
    except ValueError:
        return None
    return find(src_idx) != find(tgt_idx)


def _local_mismatch_at_codon(
    codon: str,
    code: dict[str, str],
    ctx: ReferenceContext,
    metric: str,
) -> float:
    """Compute local mismatch cost at a single codon position.

    Sum of physicochemical distances to all Hamming-1 neighbors.
    """
    dist_fn = get_distance_func(metric)
    aa = code[codon]
    total = 0.0
    for nbr in ctx.neighbors[codon]:
        nbr_aa = code[nbr]
        if aa != nbr_aa:
            total += dist_fn(aa, nbr_aa)
    return total


def _delta_f_reassign(
    ctx: ReferenceContext,
    source_codon: str,
    new_aa: str,
    metric: str,
    include_stops: bool = True,
) -> float:
    """Fast delta_F for one-vertex recoloring: C'[source_codon] = new_aa.

    Only edges incident to source_codon change. For each neighbor n:
      delta += dist(new_aa, C[n]) - dist(old_aa, C[n])

    This is O(6) per codon instead of O(192) for the full score.
    """
    dist_fn = get_distance_func(metric)
    old_aa = ctx.code[source_codon]

    if old_aa == new_aa:
        return 0.0

    delta = 0.0
    for nbr in ctx.neighbors[source_codon]:
        nbr_aa = ctx.code[nbr]

        # Old contribution: edge (old_aa, nbr_aa)
        old_contrib = 0.0
        if old_aa != nbr_aa:
            if include_stops or (old_aa != "Stop" and nbr_aa != "Stop"):
                old_contrib = dist_fn(old_aa, nbr_aa)

        # New contribution: edge (new_aa, nbr_aa)
        new_contrib = 0.0
        if new_aa != nbr_aa:
            if include_stops or (new_aa != "Stop" and nbr_aa != "Stop"):
                new_contrib = dist_fn(new_aa, nbr_aa)

        delta += new_contrib - old_contrib

    return delta


def classify_swap_event(
    *,
    event: CodonSwapEvent,
    ctx: ReferenceContext,
    swap_model: SwapModel = SwapModel.REASSIGN_SOURCE_CODON,
    metrics: tuple[str, ...] | None = None,
    include_stops: bool = True,
) -> TopologyClassification:
    """Classify a single codon swap event under a given encoding.

    Parameters
    ----------
    event : CodonSwapEvent
        The swap to classify.
    ctx : ReferenceContext
        Precomputed context for (table_id, encoding).
    swap_model : SwapModel
        How the swap modifies the code map.
    metrics : tuple of str or None
        Distance metrics to use. Default: all four available metrics.
    include_stops : bool
        Whether to include stop-codon edges in delta_F.

    Returns
    -------
    TopologyClassification
        Complete topology and mismatch annotation.
    """
    if metrics is None:
        metrics = AVAILABLE_METRICS

    source_vec = ctx.vectors[event.source_codon]
    target_vec = ctx.vectors[event.target_codon]
    hamming = hamming_distance(source_vec, target_vec)

    source_aa = ctx.code[event.source_codon]
    target_aa = ctx.code[event.target_codon]
    is_synonymous = source_aa == target_aa

    # --- Connected component analysis ---
    if swap_model == SwapModel.SYNONYMOUS_REPLACEMENT:
        # For synonymous swaps, the code map doesn't change, so
        # beta_0 is unchanged by definition.
        changes = False
        delta_by_aa: dict[str, int] = {}
        affected = ()
    else:
        # REASSIGN_SOURCE_CODON: C'[source] = target_aa
        # Check which AAs are affected
        affected_set = {source_aa, target_aa} - {"Stop"}
        affected = tuple(sorted(affected_set))

        # Compute beta_0 before and after
        code_before = ctx.code
        code_after = dict(code_before)
        code_after[event.source_codon] = target_aa

        delta_by_aa = {}
        changes = False
        for aa in affected:
            b0_before = _beta0_eps1_for_aa(code_before, aa, ctx.vectors)
            b0_after = _beta0_eps1_for_aa(code_after, aa, ctx.vectors)
            delta_by_aa[aa] = b0_after - b0_before
            if b0_after != b0_before:
                changes = True

    # --- Ser boundary crossing ---
    crosses_ser = _crosses_component_boundary(
        ctx.code, ctx.vectors, event.source_codon, event.target_codon
    )

    # --- Delta F under each metric ---
    delta_f: dict[str, float] = {}
    local_src: dict[str, float] = {}
    local_tgt: dict[str, float] = {}

    for m in metrics:
        if swap_model == SwapModel.SYNONYMOUS_REPLACEMENT:
            # No code change => delta_F = 0
            delta_f[m] = 0.0
        else:
            delta_f[m] = _delta_f_reassign(
                ctx,
                event.source_codon,
                target_aa,
                m,
                include_stops=include_stops,
            )

        local_src[m] = _local_mismatch_at_codon(
            event.source_codon, ctx.code, ctx, m
        )
        local_tgt[m] = _local_mismatch_at_codon(
            event.target_codon, ctx.code, ctx, m
        )

    return TopologyClassification(
        encoding_id=ctx.encoding_id,
        encoding=ctx.encoding,
        source_vec=source_vec,
        target_vec=target_vec,
        hamming=hamming,
        source_aa=source_aa,
        target_aa=target_aa,
        is_synonymous=is_synonymous,
        changes_components_eps1=changes,
        affected_aas=affected,
        delta_components_by_aa=delta_by_aa,
        crosses_component_boundary_eps1=crosses_ser,
        delta_edge_mismatch=delta_f,
        local_mismatch_source=local_src,
        local_mismatch_target=local_tgt,
    )


def classify_all_encodings(
    *,
    event: CodonSwapEvent,
    table_id: int = 1,
    swap_model: SwapModel = SwapModel.REASSIGN_SOURCE_CODON,
    metrics: tuple[str, ...] | None = None,
) -> list[TopologyClassification]:
    """Classify a swap under all 24 base-to-bit encodings.

    Returns a list of 24 TopologyClassification objects, one per encoding.
    Used for encoding sensitivity analysis.
    """
    results = []
    for enc_id, enc in enumerate(all_encodings()):
        ctx = build_reference_context(
            table_id=table_id,
            encoding_id=enc_id,
            encoding=enc,
        )
        topo = classify_swap_event(
            event=event,
            ctx=ctx,
            swap_model=swap_model,
            metrics=metrics,
        )
        results.append(topo)
    return results
