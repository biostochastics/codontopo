"""Persistent homology: connected components at Hamming distance thresholds."""
from collections import defaultdict
from codon_topo.core.encoding import codon_to_vector, hamming_distance, DEFAULT_ENCODING


def _partition(vectors: list[tuple[int, ...]], epsilon: int) -> list[list[int]]:
    """Partition vector indices into connected components at Hamming threshold.

    Returns list of clusters, where each cluster is a list of indices.
    """
    n = len(vectors)
    if n == 0:
        return []
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
            if hamming_distance(vectors[i], vectors[j]) <= epsilon:
                union(i, j)

    clusters: dict[int, list[int]] = defaultdict(list)
    for i in range(n):
        clusters[find(i)].append(i)
    return list(clusters.values())


def connected_components(vectors: list[tuple[int, ...]], epsilon: int) -> int:
    """Count connected components at Hamming distance threshold epsilon.

    Two vectors are in the same component if there is a path of vectors
    where each consecutive pair has Hamming distance <= epsilon.
    """
    return len(_partition(vectors, epsilon))


def persistent_homology(
    vectors: list[tuple[int, ...]],
    max_eps: int = 6,
) -> dict[int, int]:
    """Compute beta_0 (number of connected components) at each epsilon.

    Returns {epsilon: beta_0} for epsilon in 1..max_eps.
    """
    return {eps: connected_components(vectors, eps) for eps in range(1, max_eps + 1)}


def disconnection_catalogue(
    code: dict[str, str],
    encoding: dict | None = None,
) -> list[dict]:
    """Find all amino acids with disconnected codon graphs at epsilon=1.

    Returns list of dicts with keys: aa, n_components, reconnect_eps,
    min_inter_distance, blocks.
    """
    enc = encoding or DEFAULT_ENCODING
    aa_codons: dict[str, list[str]] = defaultdict(list)
    for codon, aa in code.items():
        if aa != 'Stop':
            aa_codons[aa].append(codon)

    catalogue = []
    for aa in sorted(aa_codons):
        codons = sorted(aa_codons[aa])
        if len(codons) < 2:
            continue
        vectors = [codon_to_vector(c, enc) for c in codons]
        clusters = _partition(vectors, 1)
        if len(clusters) <= 1:
            continue

        blocks = [[codons[i] for i in cluster] for cluster in clusters]

        # Build index-to-cluster map for inter-block distance
        idx_to_cluster = {}
        for ci, cluster in enumerate(clusters):
            for i in cluster:
                idx_to_cluster[i] = ci

        n = len(vectors)
        min_inter = min(
            hamming_distance(vectors[i], vectors[j])
            for i in range(n) for j in range(i + 1, n)
            if idx_to_cluster[i] != idx_to_cluster[j]
        )

        # Reconnection epsilon
        reconnect_eps = None
        for eps in range(2, 7):
            if connected_components(vectors, eps) == 1:
                reconnect_eps = eps
                break

        catalogue.append({
            'aa': aa,
            'n_components': len(clusters),
            'reconnect_eps': reconnect_eps,
            'min_inter_distance': min_inter,
            'blocks': blocks,
        })

    return catalogue
