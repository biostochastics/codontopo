"""Null models for statistical validation of filtration invariants.

Model A: Fix degeneracy structure, randomly assign codons to AAs.
Model B: Preserve block structure, shuffle block-to-AA assignments.
Model C: Test all 24 base-to-bit encodings.
"""

import random
from collections import defaultdict

from codon_topo.core.encoding import (
    ALL_CODONS,
    codon_to_vector,
    all_encodings,
    DEFAULT_ENCODING,
)
from codon_topo.core.genetic_codes import STANDARD
from codon_topo.core.homology import connected_components


def _degeneracy_structure(code: dict[str, str]) -> list[int]:
    """Extract degeneracy counts: how many codons each AA has."""
    aa_counts: dict[str, int] = defaultdict(int)
    for codon, aa in code.items():
        if aa != "Stop":
            aa_counts[aa] += 1
    return sorted(aa_counts.values())


_VECTOR_CACHE: dict[frozenset, dict[str, tuple[int, ...]]] = {}


def _get_vectors(
    encoding: dict[str, tuple[int, int]] | None = None,
) -> dict[str, tuple[int, ...]]:
    """Return cached codon->vector mapping for a given encoding."""
    enc = encoding or DEFAULT_ENCODING
    key = frozenset(enc.items())
    if key not in _VECTOR_CACHE:
        _VECTOR_CACHE[key] = {c: codon_to_vector(c, enc) for c in ALL_CODONS}
    return _VECTOR_CACHE[key]


def _score_code(
    code: dict[str, str],
    encoding: dict[str, tuple[int, int]] | None = None,
) -> dict:
    """Compute filtration and topology metrics for a code.

    Uses cached vector lookups to avoid redundant codon_to_vector calls
    inside permutation loops.
    """
    vec_map = _get_vectors(encoding)
    expected_suffixes = {(0, 0), (0, 1), (1, 0), (1, 1)}

    # Group codon vectors by AA
    groups: dict[str, list[tuple[int, ...]]] = defaultdict(list)
    for codon, aa in code.items():
        if aa != "Stop":
            groups[aa].append(vec_map[codon])

    tw_all_pass = True
    ff_all_pass = True
    disconnected = []

    for aa, vectors in groups.items():
        n_codons = len(vectors)

        if n_codons == 2:
            v1, v2 = vectors
            differing = [i for i in range(6) if v1[i] != v2[i]]
            if differing != [5]:
                tw_all_pass = False

        elif n_codons == 4:
            prefixes = set(v[:4] for v in vectors)
            suffixes = set(v[4:] for v in vectors)
            if len(prefixes) != 1 or suffixes != expected_suffixes:
                ff_all_pass = False

        if n_codons >= 2:
            if connected_components(vectors, 1) > 1:
                disconnected.append(aa)

    return {
        "twofold_all_pass": tw_all_pass,
        "fourfold_all_pass": ff_all_pass,
        "n_disconnected": len(disconnected),
        "exactly_one_disconnected": len(disconnected) == 1,
        "disconnected_aas": sorted(disconnected),
    }


def null_model_a(
    n_permutations: int = 100_000,
    code: dict[str, str] | None = None,
    seed: int | None = None,
) -> dict:
    """Null Model A: fix degeneracy structure, randomly assign codons.

    Keeps the same number of AAs with 1, 2, 3, 4, 6 codons but
    randomly assigns which codons go to which AA.

    Args:
        code: Genetic code to use as reference. Defaults to STANDARD.
    """
    rng = random.Random(seed)
    ref = code or STANDARD

    aa_sizes = _degeneracy_structure(ref)
    sense_codons = [c for c in ALL_CODONS if ref[c] != "Stop"]
    stop_codons = [c for c in ALL_CODONS if ref[c] == "Stop"]

    obs = _score_code(ref)

    count_serine_unique = 0
    count_bit5_uniform = 0
    count_fourfold_uniform = 0

    for _ in range(n_permutations):
        shuffled = list(sense_codons)
        rng.shuffle(shuffled)

        perm_code = {c: "Stop" for c in stop_codons}
        idx = 0
        for i, size in enumerate(aa_sizes):
            aa_name = f"AA{i}"
            for c in shuffled[idx : idx + size]:
                perm_code[c] = aa_name
            idx += size

        score = _score_code(perm_code)
        if score["exactly_one_disconnected"]:
            count_serine_unique += 1
        if score["twofold_all_pass"]:
            count_bit5_uniform += 1
        if score["fourfold_all_pass"]:
            count_fourfold_uniform += 1

    return {
        "n_permutations": n_permutations,
        "p_value_serine_unique": count_serine_unique / n_permutations,
        "p_value_bit5_uniform": count_bit5_uniform / n_permutations,
        "p_value_fourfold_uniform": count_fourfold_uniform / n_permutations,
        "observed_exactly_one_disconnected": obs["exactly_one_disconnected"],
        "observed_twofold_all_pass": obs["twofold_all_pass"],
        "observed_fourfold_all_pass": obs["fourfold_all_pass"],
    }


def null_model_b(
    n_permutations: int = 100_000,
    code: dict[str, str] | None = None,
    include_stops: bool = True,
    seed: int | None = None,
) -> dict:
    """Null Model B: preserve block structure, shuffle block-to-AA assignment.

    A 'block' is a group of codons sharing the first two nucleotide positions
    (16 blocks of 4 codons each). This shuffles which block gets which
    AA assignment pattern.

    Note on stop codons: by default, blocks containing Stop codons (UA, UG)
    participate in the shuffle, testing joint randomness of AA AND stop
    placement. Set include_stops=False to fix stop-containing blocks in
    place and test only sense-codon block assignments.

    Args:
        code: Genetic code to use as reference. Defaults to STANDARD.
        include_stops: If True (default), stop-containing blocks participate
            in the shuffle. If False, they are held fixed.
    """
    rng = random.Random(seed)
    ref = code or STANDARD

    # Build blocks: groups by first-2-base prefix
    blocks: dict[str, list[str]] = defaultdict(list)
    for codon in ALL_CODONS:
        blocks[codon[:2]].append(codon)
    block_list = list(blocks.values())  # 16 blocks

    # Get AA pattern per block
    std_block_aas = []
    for block_codons in block_list:
        aas = [ref[c] for c in block_codons]
        std_block_aas.append(tuple(aas))

    # Identify which blocks contain stops
    has_stop = [any(aa == "Stop" for aa in pattern) for pattern in std_block_aas]

    count_serine_unique = 0

    for _ in range(n_permutations):
        if not include_stops:
            # Only shuffle non-stop blocks
            sense_indices = [i for i, s in enumerate(has_stop) if not s]
            sense_patterns = [std_block_aas[i] for i in sense_indices]
            rng.shuffle(sense_patterns)
            shuffled_patterns = list(std_block_aas)
            for idx, pat in zip(sense_indices, sense_patterns):
                shuffled_patterns[idx] = pat
        else:
            shuffled_patterns = list(std_block_aas)
            rng.shuffle(shuffled_patterns)

        perm_code = {}
        for block_codons, pattern in zip(block_list, shuffled_patterns):
            for codon, aa in zip(block_codons, pattern):
                perm_code[codon] = aa

        score = _score_code(perm_code)
        if score["exactly_one_disconnected"]:
            count_serine_unique += 1

    return {
        "n_permutations": n_permutations,
        "p_value_serine_unique": count_serine_unique / n_permutations,
        "include_stops": include_stops,
    }


def null_model_c(
    code: dict[str, str] | None = None,
) -> dict:
    """Null Model C: test all 24 base-to-bit encodings.

    Determines whether filtration properties are encoding-dependent
    or encoding-invariant.

    Args:
        code: Genetic code to test. Defaults to STANDARD.
    """
    ref = code or STANDARD
    encodings = all_encodings()
    results = []
    twofold_invariant = 0

    for enc in encodings:
        score = _score_code(ref, encoding=enc)
        results.append(
            {
                "encoding": enc,
                "twofold_all_pass": score["twofold_all_pass"],
                "fourfold_all_pass": score["fourfold_all_pass"],
                "exactly_one_disconnected": score["exactly_one_disconnected"],
                "n_disconnected": score["n_disconnected"],
                "disconnected_aas": score["disconnected_aas"],
            }
        )
        if score["twofold_all_pass"]:
            twofold_invariant += 1

    return {
        "n_encodings": len(encodings),
        "twofold_invariant_count": twofold_invariant,
        "encoding_results": results,
    }


def null_model_c_extended(
    code: dict[str, str] | None = None,
) -> dict:
    """Null Model C EXTENDED: emit per-encoding MIN INTER-BLOCK HAMMING DISTANCE
    for every disconnected amino acid.

    The claim "Serine min distance 4 across all 24 encodings" is false (16/24
    give min=2). What is universal: disconnected at eps=1 in all encodings.

    For each of the 24 encodings and each amino acid whose codon graph is
    disconnected at eps=1, this emits:
        - n_components
        - min_inter_distance (can vary across encodings)
        - reconnect_eps (smallest eps at which the graph becomes connected)

    Returns a structure that lets you aggregate: for how many encodings does
    Serine have min distance 4 vs 2? vs 3?
    """
    from codon_topo.core.homology import disconnection_catalogue

    ref = code or STANDARD
    encodings = all_encodings()

    per_encoding: list[dict] = []
    # Aggregate: {aa: {min_dist_value: count_across_encodings}}
    aa_distance_histogram: dict[str, dict[int, int]] = defaultdict(
        lambda: defaultdict(int)
    )
    # Aggregate: {aa: {reconnect_eps: count_across_encodings}}
    aa_reconnect_histogram: dict[str, dict[int, int]] = defaultdict(
        lambda: defaultdict(int)
    )

    for enc_idx, enc in enumerate(encodings):
        cat = disconnection_catalogue(ref, encoding=enc)
        encoding_summary = {
            "encoding_index": enc_idx,
            "encoding_map": dict(enc),
            "disconnected": [
                {
                    "aa": entry["aa"],
                    "n_components": entry["n_components"],
                    "min_inter_distance": entry["min_inter_distance"],
                    "reconnect_eps": entry["reconnect_eps"],
                }
                for entry in cat
            ],
        }
        per_encoding.append(encoding_summary)

        for entry in cat:
            aa = entry["aa"]
            aa_distance_histogram[aa][entry["min_inter_distance"]] += 1
            if entry["reconnect_eps"] is not None:
                aa_reconnect_histogram[aa][entry["reconnect_eps"]] += 1

    # Universal invariants: AA that is disconnected in ALL 24 encodings
    universal_disconnected_aas = [
        aa
        for aa in aa_distance_histogram
        if sum(aa_distance_histogram[aa].values()) == len(encodings)
    ]

    # For each universal AA, enumerate the distance values it takes
    invariant_details: dict[str, dict] = {}
    for aa in universal_disconnected_aas:
        invariant_details[aa] = {
            "min_distance_values": dict(aa_distance_histogram[aa]),
            "reconnect_eps_values": dict(aa_reconnect_histogram[aa]),
            "distance_is_invariant": len(aa_distance_histogram[aa]) == 1,
            "reconnect_is_invariant": len(aa_reconnect_histogram[aa]) == 1,
        }

    return {
        "n_encodings": len(encodings),
        "per_encoding": per_encoding,
        "universal_disconnected_aas": sorted(universal_disconnected_aas),
        "invariant_details": invariant_details,
        "aa_distance_histogram": {
            aa: dict(counts) for aa, counts in aa_distance_histogram.items()
        },
        "aa_reconnect_histogram": {
            aa: dict(counts) for aa, counts in aa_reconnect_histogram.items()
        },
    }
