"""Null models for statistical validation of filtration invariants.

Model A: Fix degeneracy structure, randomly assign codons to AAs.
Model B: Preserve block structure, shuffle block-to-AA assignments.
Model C: Test all 24 base-to-bit encodings.
"""
import random
from collections import defaultdict

from codon_topo.core.encoding import (
    ALL_CODONS, codon_to_vector, all_encodings, DEFAULT_ENCODING,
)
from codon_topo.core.genetic_codes import STANDARD
from codon_topo.core.filtration import check_twofold, check_fourfold
from codon_topo.core.homology import connected_components


def _degeneracy_structure(code: dict[str, str]) -> list[int]:
    """Extract degeneracy counts: how many codons each AA has."""
    aa_counts: dict[str, int] = defaultdict(int)
    for codon, aa in code.items():
        if aa != 'Stop':
            aa_counts[aa] += 1
    return sorted(aa_counts.values())


def _score_code(code: dict[str, str], encoding=None) -> dict:
    """Compute filtration and topology metrics for a code."""
    enc = encoding or DEFAULT_ENCODING
    # Group codons by AA
    aa_codons: dict[str, list[str]] = defaultdict(list)
    for codon, aa in code.items():
        if aa != 'Stop':
            aa_codons[aa].append(codon)

    # Two-fold bit-5 check
    tw_results = check_twofold(code, enc)
    tw_all_pass = all(ok for _, ok, _ in tw_results) if tw_results else True

    # Four-fold prefix check
    ff_results = check_fourfold(code, enc)
    ff_all_pass = all(ok for _, ok in ff_results) if ff_results else True

    # Count disconnected AAs at eps=1
    disconnected = []
    for aa, codons in aa_codons.items():
        if len(codons) < 2:
            continue
        vectors = [codon_to_vector(c, enc) for c in codons]
        n_comp = connected_components(vectors, 1)
        if n_comp > 1:
            disconnected.append(aa)

    return {
        'twofold_all_pass': tw_all_pass,
        'fourfold_all_pass': ff_all_pass,
        'n_disconnected': len(disconnected),
        'exactly_one_disconnected': len(disconnected) == 1,
        'disconnected_aas': sorted(disconnected),
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
    sense_codons = [c for c in ALL_CODONS if ref[c] != 'Stop']
    stop_codons = [c for c in ALL_CODONS if ref[c] == 'Stop']

    obs = _score_code(ref)

    count_serine_unique = 0
    count_bit5_uniform = 0
    count_fourfold_uniform = 0

    for _ in range(n_permutations):
        shuffled = list(sense_codons)
        rng.shuffle(shuffled)

        perm_code = {c: 'Stop' for c in stop_codons}
        idx = 0
        for i, size in enumerate(aa_sizes):
            aa_name = f'AA{i}'
            for c in shuffled[idx:idx + size]:
                perm_code[c] = aa_name
            idx += size

        score = _score_code(perm_code)
        if score['exactly_one_disconnected']:
            count_serine_unique += 1
        if score['twofold_all_pass']:
            count_bit5_uniform += 1
        if score['fourfold_all_pass']:
            count_fourfold_uniform += 1

    return {
        'n_permutations': n_permutations,
        'p_value_serine_unique': count_serine_unique / n_permutations,
        'p_value_bit5_uniform': count_bit5_uniform / n_permutations,
        'p_value_fourfold_uniform': count_fourfold_uniform / n_permutations,
        'observed_exactly_one_disconnected': obs['exactly_one_disconnected'],
        'observed_twofold_all_pass': obs['twofold_all_pass'],
        'observed_fourfold_all_pass': obs['fourfold_all_pass'],
    }


def null_model_b(
    n_permutations: int = 100_000,
    code: dict[str, str] | None = None,
    exclude_stops: bool = False,
    seed: int | None = None,
) -> dict:
    """Null Model B: preserve block structure, shuffle block-to-AA assignment.

    A 'block' is a group of codons sharing the first two nucleotide positions
    (16 blocks of 4 codons each). This shuffles which block gets which
    AA assignment pattern.

    Note on stop codons: by default, blocks containing Stop codons (UA, UG)
    participate in the shuffle, testing joint randomness of AA AND stop
    placement. Set exclude_stops=True to fix stop-containing blocks in
    place and test only sense-codon block assignments.

    Args:
        code: Genetic code to use as reference. Defaults to STANDARD.
        exclude_stops: If True, blocks containing Stop codons are not shuffled.
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
    has_stop = [any(aa == 'Stop' for aa in pattern) for pattern in std_block_aas]

    count_serine_unique = 0

    for _ in range(n_permutations):
        if exclude_stops:
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
        if score['exactly_one_disconnected']:
            count_serine_unique += 1

    return {
        'n_permutations': n_permutations,
        'p_value_serine_unique': count_serine_unique / n_permutations,
        'exclude_stops': exclude_stops,
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
        results.append({
            'encoding': enc,
            'twofold_all_pass': score['twofold_all_pass'],
            'fourfold_all_pass': score['fourfold_all_pass'],
            'exactly_one_disconnected': score['exactly_one_disconnected'],
            'n_disconnected': score['n_disconnected'],
            'disconnected_aas': score['disconnected_aas'],
        })
        if score['twofold_all_pass']:
            twofold_invariant += 1

    return {
        'n_encodings': len(encodings),
        'twofold_invariant_count': twofold_invariant,
        'encoding_results': results,
    }
