"""Persistent homology: connected components at Hamming distance thresholds."""
from collections import defaultdict
from codon_topo.core.encoding import codon_to_vector, hamming_distance, DEFAULT_ENCODING


def connected_components(vectors: list[tuple[int, ...]], epsilon: int) -> int:
    """Count connected components at Hamming distance threshold epsilon.

    Two vectors are in the same component if there is a path of vectors
    where each consecutive pair has Hamming distance <= epsilon.
    """
    n = len(vectors)
    if n == 0:
        return 0
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

    return len(set(find(i) for i in range(n)))


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
        n_comp = connected_components(vectors, 1)
        if n_comp <= 1:
            continue

        # Find blocks via union-find at eps=1
        n = len(vectors)
        parent = list(range(n))
        def find(x, p=parent):
            while p[x] != x:
                p[x] = p[p[x]]
                x = p[x]
            return x
        def union(x, y, p=parent):
            px, py = find(x), find(y)
            if px != py:
                p[px] = py
        for i in range(n):
            for j in range(i + 1, n):
                if hamming_distance(vectors[i], vectors[j]) <= 1:
                    union(i, j)

        clusters: dict[int, list[str]] = defaultdict(list)
        for i in range(n):
            clusters[find(i)].append(codons[i])
        blocks = list(clusters.values())

        # Min inter-block Hamming distance
        min_inter = min(
            hamming_distance(vectors[i], vectors[j])
            for i in range(n) for j in range(i + 1, n)
            if find(i) != find(j)
        )

        # Reconnection epsilon
        reconnect_eps = None
        for eps in range(2, 7):
            if connected_components(vectors, eps) == 1:
                reconnect_eps = eps
                break

        catalogue.append({
            'aa': aa,
            'n_components': n_comp,
            'reconnect_eps': reconnect_eps,
            'min_inter_distance': min_inter,
            'blocks': blocks,
        })

    return catalogue
