"""Filtration analysis: two-fold (bit-5) and four-fold (prefix uniformity) checks."""
from collections import defaultdict
from codon_topo.core.encoding import codon_to_vector, DEFAULT_ENCODING


def classify_degeneracy(
    code: dict[str, str],
) -> dict[int, list[str]]:
    """Group amino acids by degeneracy (number of codons)."""
    aa_codons: dict[str, list[str]] = defaultdict(list)
    for codon, aa in code.items():
        if aa != 'Stop':
            aa_codons[aa].append(codon)
    result: dict[int, list[str]] = defaultdict(list)
    for aa, codons in sorted(aa_codons.items()):
        result[len(codons)].append(aa)
    return dict(result)


def _group_by_aa(
    code: dict[str, str],
    encoding: dict | None = None,
) -> dict[str, list[tuple[str, tuple[int, ...]]]]:
    """Return {aa: [(codon, vector), ...]} excluding stops."""
    enc = encoding or DEFAULT_ENCODING
    groups: dict[str, list[tuple[str, tuple[int, ...]]]] = defaultdict(list)
    for codon, aa in code.items():
        if aa != 'Stop':
            groups[aa].append((codon, codon_to_vector(codon, enc)))
    return dict(groups)


def check_twofold(
    code: dict[str, str],
    encoding: dict | None = None,
) -> list[tuple[str, bool, list[int]]]:
    """Check two-fold degeneracy: synonymous codons differ at exactly bit 5.

    Returns list of (aa, passed, differing_positions).
    """
    groups = _group_by_aa(code, encoding)
    results = []
    for aa in sorted(groups):
        members = groups[aa]
        if len(members) != 2:
            continue
        (_, v1), (_, v2) = members
        differing = [i for i in range(6) if v1[i] != v2[i]]
        results.append((aa, differing == [5], differing))
    return results


def check_fourfold(
    code: dict[str, str],
    encoding: dict | None = None,
) -> list[tuple[str, bool]]:
    """Check four-fold degeneracy: identical 4-bit prefix, last 2 exhaust GF(2)^2.

    Returns list of (aa, passed).
    """
    groups = _group_by_aa(code, encoding)
    expected_suffixes = {(0,0),(0,1),(1,0),(1,1)}
    results = []
    for aa in sorted(groups):
        members = groups[aa]
        if len(members) != 4:
            continue
        vectors = [v for _, v in members]
        prefixes = set(v[:4] for v in vectors)
        suffixes = set(v[4:] for v in vectors)
        passed = len(prefixes) == 1 and suffixes == expected_suffixes
        results.append((aa, passed))
    return results


def analyze_filtration(
    code: dict[str, str],
    encoding: dict | None = None,
) -> dict:
    """Full filtration report for a genetic code."""
    tw = check_twofold(code, encoding)
    ff = check_fourfold(code, encoding)
    return {
        'twofold_pass': sum(1 for _, p, _ in tw if p),
        'twofold_fail': sum(1 for _, p, _ in tw if not p),
        'fourfold_pass': sum(1 for _, p in ff if p),
        'fourfold_fail': sum(1 for _, p in ff if not p),
        'twofold_detail': tw,
        'fourfold_detail': ff,
    }
