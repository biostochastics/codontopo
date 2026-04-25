# WS1: Core Package Implementation Plan

**Goal:** Refactor the 6 preliminary verification scripts into a tested, pip-installable Python package `codon_topo` that reproduces all verified results as regression tests, then implement the three null models.

**Architecture:** Bottom-up build. First consolidate shared primitives (encoding, Hamming distance) into `core/`, then layer filtration, homology, embedding, and Fano modules on top. Each module gets property-based tests via hypothesis. Null models go in `analysis/`. Visualization is R/ggplot2 scripts in `visualization/R/`. The preliminary scripts stay untouched as the reference oracle.

**Tech Stack:** Python 3.11+ (stdlib + NumPy + SciPy), pytest + hypothesis for testing, ggplot2 + ggpubr (R) for figures.

---

### Task 0: Project scaffolding

**Files:**
- Create: `pyproject.toml`
- Create: `src/codon_topo/__init__.py`
- Create: `src/codon_topo/core/__init__.py`
- Create: `src/codon_topo/analysis/__init__.py`
- Create: `tests/__init__.py`
- Create: `tests/conftest.py`

**Step 1: Create pyproject.toml**

```toml
[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[project]
name = "codon-topo"
version = "0.1.0"
description = "Codon Geometry Validation & Prediction Engine"
requires-python = ">=3.11"
dependencies = ["numpy>=1.24", "scipy>=1.10"]

[project.optional-dependencies]
dev = ["pytest>=7.0", "hypothesis>=6.0", "pytest-cov"]

[tool.pytest.ini_options]
testpaths = ["tests"]
```

**Step 2: Create empty package init files**

`src/codon_topo/__init__.py`:
```python
"""Codon Geometry Validation & Prediction Engine."""
```

`src/codon_topo/core/__init__.py` and `src/codon_topo/analysis/__init__.py`: empty files.

`tests/__init__.py`: empty file.

`tests/conftest.py`:
```python
"""Shared test fixtures for codon_topo."""
```

**Step 3: Install in dev mode and verify**

Run: `pip install -e ".[dev]"`
Run: `pytest --co -q` (should collect 0 tests, no errors)

**Step 4: Commit**

```bash
git add pyproject.toml src/ tests/
git commit -m "scaffold: codon_topo package structure with pyproject.toml"
```

---

### Task 1: Core encoding module

**Files:**
- Create: `src/codon_topo/core/encoding.py`
- Create: `tests/test_encoding.py`

This consolidates `BASE_TO_BITS`, `codon_to_bits`, `hamming_distance`, and `ALL_CODONS` from the preliminary scripts. Also adds support for all 24 base-to-bit encodings (WS1 Null Model C requirement).

**Step 1: Write failing tests**

`tests/test_encoding.py` — test the core primitives:
```python
import pytest
from hypothesis import given, strategies as st
from codon_topo.core.encoding import (
    DEFAULT_ENCODING, codon_to_vector, hamming_distance,
    ALL_CODONS, BASES, all_encodings,
)

def test_default_encoding():
    assert DEFAULT_ENCODING == {'C': (0,0), 'U': (0,1), 'A': (1,0), 'G': (1,1)}

def test_codon_to_vector_ggu():
    """GGU = (1,1,1,1,0,1) per verify_kras.py line 30-31."""
    assert codon_to_vector('GGU') == (1,1,1,1,0,1)

def test_codon_to_vector_cac():
    """CAC = (0,0,1,0,0,0) per verify_kras.py."""
    assert codon_to_vector('CAC') == (0,0,1,0,0,0)

def test_all_codons_count():
    assert len(ALL_CODONS) == 64

def test_all_codons_unique():
    assert len(set(ALL_CODONS)) == 64

def test_hamming_self():
    v = codon_to_vector('GGU')
    assert hamming_distance(v, v) == 0

def test_hamming_ggu_guu():
    """d(GGU, GUU) = 1 per verify_kras.py line 34."""
    ggu = codon_to_vector('GGU')
    guu = codon_to_vector('GUU')
    assert hamming_distance(ggu, guu) == 1

@given(st.sampled_from(['UUU','UUC','UUA','UUG','GGG','AAA','CCC']))
def test_vector_length_is_6(codon):
    v = codon_to_vector(codon)
    assert len(v) == 6
    assert all(b in (0,1) for b in v)

@given(
    st.sampled_from(['UUU','GGG','AAA','CCC','AUG','UAG']),
    st.sampled_from(['UUU','GGG','AAA','CCC','AUG','UAG']),
)
def test_hamming_symmetry(c1, c2):
    v1, v2 = codon_to_vector(c1), codon_to_vector(c2)
    assert hamming_distance(v1, v2) == hamming_distance(v2, v1)

@given(
    st.sampled_from(['UUU','GGG','AAA','CCC','AUG']),
    st.sampled_from(['UUU','GGG','AAA','CCC','AUG']),
)
def test_hamming_range(c1, c2):
    d = hamming_distance(codon_to_vector(c1), codon_to_vector(c2))
    assert 0 <= d <= 6

def test_all_encodings_count():
    """4! = 24 possible base-to-bit encodings."""
    encs = all_encodings()
    assert len(encs) == 24

def test_all_encodings_are_bijections():
    for enc in all_encodings():
        assert len(enc) == 4
        assert set(enc.keys()) == set(BASES)
        assert len(set(enc.values())) == 4  # all values distinct
```

**Step 2: Run tests to verify failure**

Run: `pytest tests/test_encoding.py -v`
Expected: ImportError — module doesn't exist yet.

**Step 3: Implement encoding.py**

`src/codon_topo/core/encoding.py`:
```python
"""Binary encoding of codons as vectors in GF(2)^6.

The default encoding maps each nucleotide to a 2-bit pair:
C=(0,0), U=(0,1), A=(1,0), G=(1,1). A codon (3 bases) becomes
a 6-bit tuple. All 24 permutations of the 4 possible 2-bit values
across the 4 bases are available for Null Model C analysis.
"""
from itertools import permutations

BASES = ('C', 'U', 'A', 'G')
_BIT_PAIRS = ((0,0), (0,1), (1,0), (1,1))

DEFAULT_ENCODING: dict[str, tuple[int,int]] = dict(zip(BASES, _BIT_PAIRS))

ALL_CODONS: list[str] = [
    b1 + b2 + b3
    for b1 in BASES for b2 in BASES for b3 in BASES
]


def codon_to_vector(
    codon: str,
    encoding: dict[str, tuple[int,int]] | None = None,
) -> tuple[int, ...]:
    """Convert a 3-letter codon to a 6-bit tuple in GF(2)^6."""
    enc = encoding or DEFAULT_ENCODING
    bits: list[int] = []
    for base in codon:
        bits.extend(enc[base])
    return tuple(bits)


def hamming_distance(a: tuple[int, ...], b: tuple[int, ...]) -> int:
    """Hamming distance between two bit-tuples."""
    return sum(x != y for x, y in zip(a, b))


def all_encodings() -> list[dict[str, tuple[int,int]]]:
    """Return all 24 bijections from {C,U,A,G} to {00,01,10,11}."""
    return [dict(zip(BASES, perm)) for perm in permutations(_BIT_PAIRS)]
```

**Step 4: Run tests to verify all pass**

Run: `pytest tests/test_encoding.py -v`
Expected: All pass.

**Step 5: Commit**

```bash
git add src/codon_topo/core/encoding.py tests/test_encoding.py
git commit -m "feat: core encoding module — codon vectors, Hamming distance, 24 encodings"
```

---

### Task 2: Genetic codes module

**Files:**
- Create: `src/codon_topo/core/genetic_codes.py`
- Create: `tests/test_genetic_codes.py`

Consolidates the STANDARD dict and all 23 NCBI tables from `all_codes.py`.

**Step 1: Write failing tests**

`tests/test_genetic_codes.py`:
```python
from codon_topo.core.genetic_codes import STANDARD, get_code, all_table_ids
from codon_topo.core.encoding import ALL_CODONS

def test_standard_has_64_entries():
    assert len(STANDARD) == 64

def test_standard_codons_match():
    assert set(STANDARD.keys()) == set(ALL_CODONS)

def test_standard_ggu_is_gly():
    assert STANDARD['GGU'] == 'Gly'

def test_standard_uaa_is_stop():
    assert STANDARD['UAA'] == 'Stop'

def test_standard_aug_is_met():
    assert STANDARD['AUG'] == 'Met'

def test_all_table_ids():
    ids = all_table_ids()
    assert 1 in ids
    assert 2 in ids
    assert len(ids) == 23

def test_get_code_table1_is_standard():
    assert get_code(1) == STANDARD

def test_get_code_vert_mito_uga_trp():
    code = get_code(2)
    assert code['UGA'] == 'Trp'
    assert code['AGA'] == 'Stop'

def test_get_code_yeast_mito_cun_thr():
    code = get_code(3)
    for codon in ('CUU', 'CUC', 'CUA', 'CUG'):
        assert code[codon] == 'Thr'

def test_every_code_has_64_entries():
    for tid in all_table_ids():
        assert len(get_code(tid)) == 64, f"Table {tid} missing codons"

def test_invalid_table_raises():
    import pytest
    with pytest.raises(KeyError):
        get_code(999)
```

**Step 2: Run to verify failure**

Run: `pytest tests/test_genetic_codes.py -v`

**Step 3: Implement genetic_codes.py**

Port `STANDARD` and the `CODES` dict from `all_codes.py` lines 73-172. Use the same `make_code(changes)` pattern.

```python
"""All 23 NCBI translation tables.

Reference: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
"""

STANDARD: dict[str, str] = {
    'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
    'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
    'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'Stop', 'UAG': 'Stop',
    'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'Stop', 'UGG': 'Trp',
    'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
    'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
    'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
    'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met',
    'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
    'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
    'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
    'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
    'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
    'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly',
}

# Variant codes defined as diffs from STANDARD.
# Copied verbatim from all_codes.py lines 103-172.
_CODES: dict[int, tuple[str, dict[str, str]]] = {
    1:  ("Standard", {}),
    2:  ("Vertebrate Mitochondrial", {'AGA':'Stop','AGG':'Stop','AUA':'Met','UGA':'Trp'}),
    3:  ("Yeast Mitochondrial", {'AUA':'Met','CUU':'Thr','CUC':'Thr','CUA':'Thr','CUG':'Thr','UGA':'Trp'}),
    4:  ("Mold/Protozoan/Coelenterate Mito", {'UGA':'Trp'}),
    5:  ("Invertebrate Mitochondrial", {'AGA':'Ser','AGG':'Ser','AUA':'Met','UGA':'Trp'}),
    6:  ("Ciliate Nuclear", {'UAA':'Gln','UAG':'Gln'}),
    9:  ("Echinoderm/Flatworm Mito", {'AAA':'Asn','AGA':'Ser','AGG':'Ser','UGA':'Trp'}),
    10: ("Euplotid Nuclear", {'UGA':'Cys'}),
    11: ("Bacterial/Archaeal/Plant Plastid", {}),
    12: ("Alternative Yeast Nuclear", {'CUG':'Ser'}),
    13: ("Ascidian Mitochondrial", {'AGA':'Gly','AGG':'Gly','AUA':'Met','UGA':'Trp'}),
    14: ("Alternative Flatworm Mito", {'AAA':'Asn','AGA':'Ser','AGG':'Ser','UAA':'Tyr','UGA':'Trp'}),
    16: ("Chlorophycean Mito", {'UAG':'Leu'}),
    21: ("Trematode Mitochondrial", {'AAA':'Asn','AGA':'Ser','AGG':'Ser','AUA':'Met','UGA':'Trp'}),
    22: ("Scenedesmus obliquus Mito", {'UAG':'Leu','UCA':'Stop'}),
    23: ("Thraustochytrium Mito", {'UUA':'Stop','UGA':'Trp'}),
    24: ("Rhabdopleuridae Mito", {'AGA':'Ser','AGG':'Lys','UGA':'Trp'}),
    25: ("Candidate Division SR1 / Gracilibacteria", {'UGA':'Gly'}),
    26: ("Pachysolen tannophilus Nuclear", {'CUG':'Ala'}),
    27: ("Karyorelictea Nuclear", {'UAG':'Gln','UAA':'Gln'}),
    29: ("Mesodinium Nuclear", {'UAA':'Tyr','UAG':'Tyr'}),
    30: ("Peritrich Nuclear", {'UAA':'Glu','UAG':'Glu'}),
    31: ("Blastocrithidia Nuclear", {'UGA':'Trp','UAG':'Glu','UAA':'Glu'}),
    33: ("Cephalodiscidae Mito", {'AGA':'Ser','AGG':'Lys','UAA':'Tyr','UGA':'Trp'}),
}


def get_code(table_id: int) -> dict[str, str]:
    """Return the codon->AA mapping for an NCBI translation table."""
    name, changes = _CODES[table_id]
    code = dict(STANDARD)
    code.update(changes)
    return code


def get_code_name(table_id: int) -> str:
    """Return the human-readable name of a translation table."""
    return _CODES[table_id][0]


def all_table_ids() -> list[int]:
    """Return sorted list of all available NCBI table IDs."""
    return sorted(_CODES.keys())
```

**Step 4: Run tests**

Run: `pytest tests/test_genetic_codes.py -v`

**Step 5: Commit**

```bash
git add src/codon_topo/core/genetic_codes.py tests/test_genetic_codes.py
git commit -m "feat: genetic codes module — all 23 NCBI translation tables"
```

---

### Task 3: Filtration module

**Files:**
- Create: `src/codon_topo/core/filtration.py`
- Create: `tests/test_filtration.py`

Extracts the two-fold and four-fold filtration checks from `verify_filtration.py`.

**Step 1: Write failing tests**

`tests/test_filtration.py`:
```python
import pytest
from codon_topo.core.filtration import (
    classify_degeneracy, check_twofold, check_fourfold, analyze_filtration,
)
from codon_topo.core.genetic_codes import STANDARD, get_code, all_table_ids

def test_classify_degeneracy_standard():
    deg = classify_degeneracy(STANDARD)
    assert sorted(deg[2]) == sorted(['Phe','Tyr','His','Gln','Asn','Lys','Asp','Glu','Cys'])
    assert sorted(deg[4]) == sorted(['Val','Pro','Thr','Ala','Gly'])
    assert sorted(deg[6]) == sorted(['Leu','Arg','Ser'])
    assert deg[1] == ['Met', 'Trp'] or set(deg[1]) == {'Met','Trp'}

def test_twofold_standard_all_pass():
    """Claim 1: all 9 two-fold AAs differ at exactly bit 5."""
    results = check_twofold(STANDARD)
    for aa, passed, differing in results:
        assert passed, f"{aa} failed: differs at {differing}"

def test_fourfold_standard_all_pass():
    """Claim 2: all 5 four-fold AAs have uniform prefix."""
    results = check_fourfold(STANDARD)
    for aa, passed in results:
        assert passed, f"{aa} failed four-fold check"

@pytest.mark.parametrize("table_id", [1,2,3,4,5,6,9,10,11,12,13,14,16,21,22,23,24,25,26,27,29,30,31,33])
def test_twofold_universal(table_id):
    """Two-fold filtration is 100% invariant across all 23 codes."""
    code = get_code(table_id)
    results = check_twofold(code)
    for aa, passed, differing in results:
        assert passed, f"Table {table_id}, {aa}: differs at {differing}"

def test_analyze_filtration_returns_summary():
    report = analyze_filtration(STANDARD)
    assert report['twofold_pass'] == 9
    assert report['twofold_fail'] == 0
    assert report['fourfold_pass'] == 5
    assert report['fourfold_fail'] == 0
```

**Step 2: Run to verify failure**

Run: `pytest tests/test_filtration.py -v`

**Step 3: Implement filtration.py**

```python
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
```

**Step 4: Run tests**

Run: `pytest tests/test_filtration.py -v`

**Step 5: Commit**

```bash
git add src/codon_topo/core/filtration.py tests/test_filtration.py
git commit -m "feat: filtration module — two-fold and four-fold degeneracy checks"
```

---

### Task 4: Homology module (persistent homology / connected components)

**Files:**
- Create: `src/codon_topo/core/homology.py`
- Create: `tests/test_homology.py`

Extracts `connected_components` and adds `persistent_homology` which computes beta_0 at each epsilon from 1..6.

**Step 1: Write failing tests**

`tests/test_homology.py`:
```python
import pytest
from codon_topo.core.homology import connected_components, persistent_homology, disconnection_catalogue
from codon_topo.core.encoding import codon_to_vector
from codon_topo.core.genetic_codes import STANDARD, get_code, all_table_ids


def test_single_codon_connected():
    v = [codon_to_vector('AUG')]
    assert connected_components(v, 1) == 1


def test_serine_standard_disconnected_eps1():
    """Claim 3: Serine uniquely disconnected at eps=1."""
    ser_codons = [c for c, aa in STANDARD.items() if aa == 'Ser']
    vectors = [codon_to_vector(c) for c in ser_codons]
    assert connected_components(vectors, 1) == 2


def test_serine_reconnects_at_eps4():
    ser_codons = [c for c, aa in STANDARD.items() if aa == 'Ser']
    vectors = [codon_to_vector(c) for c in ser_codons]
    assert connected_components(vectors, 3) == 2
    assert connected_components(vectors, 4) == 1


def test_glycine_connected():
    gly_codons = [c for c, aa in STANDARD.items() if aa == 'Gly']
    vectors = [codon_to_vector(c) for c in gly_codons]
    assert connected_components(vectors, 1) == 1


def test_persistent_homology_serine():
    ser_codons = [c for c, aa in STANDARD.items() if aa == 'Ser']
    vectors = [codon_to_vector(c) for c in ser_codons]
    ph = persistent_homology(vectors, max_eps=6)
    assert ph[1] == 2  # disconnected
    assert ph[4] == 1  # reconnected


@pytest.mark.parametrize("table_id", [1,2,3,4,5,6,9,10,11,12,13,14,16,21,22,23,24,25,26,27,29,30,31,33])
def test_serine_always_disconnected(table_id):
    """Universal Serine invariant: disconnected at eps<=4 in all codes."""
    code = get_code(table_id)
    ser_codons = [c for c, aa in code.items() if aa == 'Ser']
    if len(ser_codons) < 2:
        pytest.skip("Serine has <2 codons")
    vectors = [codon_to_vector(c) for c in ser_codons]
    assert connected_components(vectors, 1) > 1, f"Table {table_id}: Serine connected at eps=1"


def test_thr_yeast_mito_disconnected():
    """Novel finding: Thr disconnected in yeast mito, eps=2."""
    code = get_code(3)
    thr_codons = [c for c, aa in code.items() if aa == 'Thr']
    vectors = [codon_to_vector(c) for c in thr_codons]
    assert connected_components(vectors, 1) == 2
    assert connected_components(vectors, 2) == 1


def test_leu_chlorophycean_disconnected():
    """Novel finding: Leu disconnected in chlorophycean mito, eps=2."""
    code = get_code(16)
    leu_codons = [c for c, aa in code.items() if aa == 'Leu']
    vectors = [codon_to_vector(c) for c in leu_codons]
    assert connected_components(vectors, 1) == 2
    assert connected_components(vectors, 2) == 1


def test_ala_pachysolen_disconnected():
    """Novel finding: Ala disconnected in Pachysolen nuclear, eps=3."""
    code = get_code(26)
    ala_codons = [c for c, aa in code.items() if aa == 'Ala']
    vectors = [codon_to_vector(c) for c in ala_codons]
    assert connected_components(vectors, 1) == 2
    assert connected_components(vectors, 3) == 1


def test_ser_candida_three_components():
    """Novel finding: Serine has 3 components in Candida (table 12)."""
    code = get_code(12)
    ser_codons = [c for c, aa in code.items() if aa == 'Ser']
    vectors = [codon_to_vector(c) for c in ser_codons]
    assert connected_components(vectors, 1) == 3


def test_disconnection_catalogue_standard():
    cat = disconnection_catalogue(STANDARD)
    aa_names = [entry['aa'] for entry in cat]
    assert 'Ser' in aa_names
    ser = [e for e in cat if e['aa'] == 'Ser'][0]
    assert ser['reconnect_eps'] == 4
    assert ser['min_inter_distance'] == 4
    assert ser['n_components'] == 2
```

**Step 2: Run to verify failure**

Run: `pytest tests/test_homology.py -v`

**Step 3: Implement homology.py**

Port the union-find `connected_components` from preliminary scripts. Add `persistent_homology` and `disconnection_catalogue`.

```python
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
```

**Step 4: Run tests**

Run: `pytest tests/test_homology.py -v`

**Step 5: Commit**

```bash
git add src/codon_topo/core/homology.py tests/test_homology.py
git commit -m "feat: homology module — persistent homology and disconnection catalogue"
```

---

### Task 5: Embedding module

**Files:**
- Create: `src/codon_topo/core/embedding.py`
- Create: `tests/test_embedding.py`

Extracts the holomorphic embedding phi: GF(2)^6 -> C^3 from `verify_embedding.py`.

**Step 1: Write failing tests**

`tests/test_embedding.py`:
```python
import cmath
from codon_topo.core.embedding import (
    BASE_TO_COMPLEX, embed_codon,
)

I = complex(0, 1)

def test_base_to_complex():
    assert BASE_TO_COMPLEX['C'] == complex(1, 0)
    assert BASE_TO_COMPLEX['U'] == complex(0, 1)
    assert BASE_TO_COMPLEX['A'] == complex(-1, 0)
    assert BASE_TO_COMPLEX['G'] == complex(0, -1)

def test_embed_ggu():
    z = embed_codon('GGU')
    assert z == (complex(0,-1), complex(0,-1), complex(0,1))

def test_twofold_share_first_two_coords():
    """Property 1: two-fold synonymous codons agree in coords 1,2."""
    from codon_topo.core.genetic_codes import STANDARD
    from collections import defaultdict
    aa_codons = defaultdict(list)
    for c, aa in STANDARD.items():
        if aa != 'Stop':
            aa_codons[aa].append(c)
    for aa, codons in aa_codons.items():
        if len(codons) != 2:
            continue
        z1 = embed_codon(codons[0])
        z2 = embed_codon(codons[1])
        assert z1[0] == z2[0], f"{aa} coord 1 differs"
        assert z1[1] == z2[1], f"{aa} coord 2 differs"

def test_fourfold_exhaust_roots():
    """Property 2: four-fold AAs have third coord cycling all 4 roots."""
    from codon_topo.core.genetic_codes import STANDARD
    from collections import defaultdict
    aa_codons = defaultdict(list)
    for c, aa in STANDARD.items():
        if aa != 'Stop':
            aa_codons[aa].append(c)
    fourth_roots = {complex(1,0), complex(0,1), complex(-1,0), complex(0,-1)}
    for aa, codons in aa_codons.items():
        if len(codons) != 4:
            continue
        embeddings = [embed_codon(c) for c in codons]
        first_coords = set((z[0], z[1]) for z in embeddings)
        third_coords = set(z[2] for z in embeddings)
        assert len(first_coords) == 1, f"{aa}: first 2 coords not uniform"
        assert third_coords == fourth_roots, f"{aa}: third coord doesn't exhaust roots"
```

**Step 2: Run to verify failure**

Run: `pytest tests/test_embedding.py -v`

**Step 3: Implement embedding.py**

```python
"""Holomorphic embedding phi: GF(2)^6 -> C^3.

Maps each nucleotide base to a fourth root of unity (i^k),
then a codon (3 bases) maps to a point in C^3.
"""

BASE_TO_COMPLEX: dict[str, complex] = {
    'C': complex(1, 0),    # i^0 = 1
    'U': complex(0, 1),    # i^1 = i
    'A': complex(-1, 0),   # i^2 = -1
    'G': complex(0, -1),   # i^3 = -i
}


def embed_codon(codon: str) -> tuple[complex, complex, complex]:
    """Map a 3-letter codon to C^3 via fourth roots of unity."""
    return (
        BASE_TO_COMPLEX[codon[0]],
        BASE_TO_COMPLEX[codon[1]],
        BASE_TO_COMPLEX[codon[2]],
    )
```

**Step 4: Run tests**

Run: `pytest tests/test_embedding.py -v`

**Step 5: Commit**

```bash
git add src/codon_topo/core/embedding.py tests/test_embedding.py
git commit -m "feat: embedding module — holomorphic embedding phi: GF(2)^6 -> C^3"
```

---

### Task 6: Fano module

**Files:**
- Create: `src/codon_topo/core/fano.py`
- Create: `tests/test_fano.py`

Extracts Fano-line (XOR triple) computation from `verify_kras.py`.

**Step 1: Write failing tests**

`tests/test_fano.py`:
```python
from codon_topo.core.fano import is_fano_line, fano_partner, all_single_bit_fano_lines
from codon_topo.core.encoding import codon_to_vector

def test_ggu_guu_cac_is_fano_line():
    """KRAS: GGU XOR GUU XOR CAC = 0."""
    assert is_fano_line('GGU', 'GUU', 'CAC')

def test_fano_partner_ggu_guu():
    assert fano_partner('GGU', 'GUU') == 'CAC'

def test_non_fano_line():
    assert not is_fano_line('GGU', 'GUU', 'AAA')

def test_fano_line_symmetric():
    assert is_fano_line('GGU', 'GUU', 'CAC')
    assert is_fano_line('GUU', 'CAC', 'GGU')
    assert is_fano_line('CAC', 'GGU', 'GUU')

def test_all_single_bit_fano_from_ggu():
    """Every single-bit mutation from GGU produces a Fano line."""
    lines = all_single_bit_fano_lines('GGU')
    assert len(lines) == 6  # 6 bit positions
    # One of them should be GGU->GUU->CAC
    codons_and_partners = [(l['mutant_codon'], l['fano_partner']) for l in lines]
    assert ('GUU', 'CAC') in codons_and_partners
```

**Step 2: Run to verify failure**

Run: `pytest tests/test_fano.py -v`

**Step 3: Implement fano.py**

```python
"""Fano-line computation in GF(2)^6.

A Fano line is a triple of vectors (a, b, c) such that a XOR b XOR c = 0.
Equivalently, c = a XOR b.
"""
from codon_topo.core.encoding import (
    codon_to_vector, ALL_CODONS, DEFAULT_ENCODING,
)

# Reverse lookup: vector -> codon (default encoding)
_VECTOR_TO_CODON: dict[tuple[int, ...], str] = {
    codon_to_vector(c): c for c in ALL_CODONS
}


def _xor(a: tuple[int, ...], b: tuple[int, ...]) -> tuple[int, ...]:
    return tuple((x + y) % 2 for x, y in zip(a, b))


def is_fano_line(c1: str, c2: str, c3: str) -> bool:
    """Check if three codons form a Fano line (XOR = 0)."""
    v1, v2, v3 = codon_to_vector(c1), codon_to_vector(c2), codon_to_vector(c3)
    return all(x == 0 for x in _xor(_xor(v1, v2), v3))


def fano_partner(c1: str, c2: str) -> str:
    """Return the third codon completing the Fano line through c1 and c2."""
    v1, v2 = codon_to_vector(c1), codon_to_vector(c2)
    v3 = _xor(v1, v2)
    return _VECTOR_TO_CODON[v3]


def all_single_bit_fano_lines(codon: str) -> list[dict]:
    """For each single-bit mutation from codon, compute the Fano partner.

    Returns list of dicts: {bit_pos, mutant_codon, mutant_aa, fano_partner, partner_aa}.
    Note: amino acid lookups require a genetic code; this returns codons only.
    """
    v = codon_to_vector(codon)
    results = []
    for bit_pos in range(6):
        mutant = list(v)
        mutant[bit_pos] = 1 - mutant[bit_pos]
        mutant = tuple(mutant)
        partner = _xor(v, mutant)
        results.append({
            'bit_pos': bit_pos,
            'mutant_codon': _VECTOR_TO_CODON[mutant],
            'fano_partner': _VECTOR_TO_CODON[partner],
        })
    return results
```

**Step 4: Run tests**

Run: `pytest tests/test_fano.py -v`

**Step 5: Commit**

```bash
git add src/codon_topo/core/fano.py tests/test_fano.py
git commit -m "feat: fano module — XOR triples and Fano-line computation"
```

---

### Task 7: Regression test suite

**Files:**
- Create: `tests/test_regression.py`

Runs the full preliminary analysis programmatically and asserts exact match with known results from PRD Appendix 8 and `summary_table.py`.

**Step 1: Write the regression tests**

`tests/test_regression.py`:
```python
"""Regression tests: reproduce all preliminary results exactly.

These values come from the verified output of the preliminary scripts
and PRD Appendix 8. Any change that breaks these tests means the
refactored code disagrees with the original verification.
"""
import pytest
from codon_topo.core.encoding import codon_to_vector
from codon_topo.core.genetic_codes import STANDARD, get_code, all_table_ids
from codon_topo.core.filtration import check_twofold, check_fourfold, analyze_filtration
from codon_topo.core.homology import (
    connected_components, disconnection_catalogue, persistent_homology,
)
from codon_topo.core.embedding import embed_codon
from codon_topo.core.fano import is_fano_line


class TestFiltrationInvariants:
    """PRD Appendix 8, Claims 1-2."""

    @pytest.mark.parametrize("table_id", all_table_ids())
    def test_twofold_100pct(self, table_id):
        code = get_code(table_id)
        results = check_twofold(code)
        failures = [(aa, diff) for aa, ok, diff in results if not ok]
        assert failures == [], f"Table {table_id}: {failures}"

    def test_fourfold_standard(self):
        results = check_fourfold(STANDARD)
        assert all(ok for _, ok in results)

    @pytest.mark.parametrize("table_id", all_table_ids())
    def test_fourfold_breaks_only_on_stop_reassignment(self, table_id):
        """Four-fold breaks only when stops are reassigned to AAs."""
        code = get_code(table_id)
        results = check_fourfold(code)
        failures = [aa for aa, ok in results if not ok]
        if failures:
            # Verify these failures involve stop codon reassignments
            stop_reassigned = {
                c for c in ('UAA','UAG','UGA')
                if code[c] != 'Stop' and STANDARD[c] == 'Stop'
            }
            assert len(stop_reassigned) > 0, (
                f"Table {table_id}: four-fold failures {failures} "
                f"without stop reassignment"
            )


class TestSerineInvariant:
    """PRD Appendix 8, Claim 3."""

    @pytest.mark.parametrize("table_id", all_table_ids())
    def test_serine_disconnected(self, table_id):
        code = get_code(table_id)
        ser = [c for c, aa in code.items() if aa == 'Ser']
        if len(ser) < 2:
            pytest.skip("Serine has <2 codons")
        vectors = [codon_to_vector(c) for c in ser]
        assert connected_components(vectors, 1) > 1

    @pytest.mark.parametrize("table_id", all_table_ids())
    def test_serine_reconnects_at_4(self, table_id):
        code = get_code(table_id)
        ser = [c for c, aa in code.items() if aa == 'Ser']
        if len(ser) < 2:
            pytest.skip("Serine has <2 codons")
        vectors = [codon_to_vector(c) for c in ser]
        assert connected_components(vectors, 4) == 1

    def test_serine_min_inter_block_distance_4(self):
        cat = disconnection_catalogue(STANDARD)
        ser = [e for e in cat if e['aa'] == 'Ser'][0]
        assert ser['min_inter_distance'] == 4


class TestNovelDisconnections:
    """PRD Appendix 8.2."""

    def test_thr_yeast_mito(self):
        cat = disconnection_catalogue(get_code(3))
        thr = [e for e in cat if e['aa'] == 'Thr']
        assert len(thr) == 1
        assert thr[0]['reconnect_eps'] == 2

    def test_leu_chlorophycean(self):
        cat = disconnection_catalogue(get_code(16))
        leu = [e for e in cat if e['aa'] == 'Leu']
        assert len(leu) == 1
        assert leu[0]['reconnect_eps'] == 2

    def test_ala_pachysolen(self):
        cat = disconnection_catalogue(get_code(26))
        ala = [e for e in cat if e['aa'] == 'Ala']
        assert len(ala) == 1
        assert ala[0]['reconnect_eps'] == 3

    def test_ser_candida_three_components(self):
        cat = disconnection_catalogue(get_code(12))
        ser = [e for e in cat if e['aa'] == 'Ser']
        assert len(ser) == 1
        assert ser[0]['n_components'] == 3


class TestEmbedding:
    """PRD Appendix 8.1, embedding properties."""

    def test_twofold_coords_12_shared(self):
        from collections import defaultdict
        aa_codons = defaultdict(list)
        for c, aa in STANDARD.items():
            if aa != 'Stop':
                aa_codons[aa].append(c)
        for aa, codons in aa_codons.items():
            if len(codons) != 2:
                continue
            z1, z2 = embed_codon(codons[0]), embed_codon(codons[1])
            assert z1[0] == z2[0] and z1[1] == z2[1], f"{aa}"

    def test_fourfold_third_coord_exhausts_roots(self):
        from collections import defaultdict
        aa_codons = defaultdict(list)
        for c, aa in STANDARD.items():
            if aa != 'Stop':
                aa_codons[aa].append(c)
        roots = {complex(1,0), complex(0,1), complex(-1,0), complex(0,-1)}
        for aa, codons in aa_codons.items():
            if len(codons) != 4:
                continue
            thirds = set(embed_codon(c)[2] for c in codons)
            assert thirds == roots, f"{aa}"


class TestKRAS:
    """PRD Appendix 8.1, KRAS Fano line."""

    def test_ggu_guu_cac(self):
        assert is_fano_line('GGU', 'GUU', 'CAC')
```

**Step 2: Run to verify all pass**

Run: `pytest tests/test_regression.py -v`
Expected: All pass (since modules from Tasks 1-6 are already implemented).

**Step 3: Commit**

```bash
git add tests/test_regression.py
git commit -m "test: comprehensive regression suite against PRD Appendix 8"
```

---

### Task 8: Null Model A — Random Assignment

**Files:**
- Create: `src/codon_topo/analysis/null_models.py`
- Create: `tests/test_null_models.py`

**Step 1: Write failing tests**

`tests/test_null_models.py`:
```python
import pytest
from codon_topo.analysis.null_models import (
    null_model_a, null_model_b, null_model_c,
)

def test_null_model_a_returns_results():
    """Smoke test with small n_permutations."""
    result = null_model_a(n_permutations=100, seed=42)
    assert 'n_permutations' in result
    assert 'p_value_serine_unique' in result
    assert 'p_value_bit5_uniform' in result
    assert 0.0 <= result['p_value_serine_unique'] <= 1.0
    assert 0.0 <= result['p_value_bit5_uniform'] <= 1.0

def test_null_model_a_deterministic():
    r1 = null_model_a(n_permutations=50, seed=123)
    r2 = null_model_a(n_permutations=50, seed=123)
    assert r1 == r2

def test_null_model_b_returns_results():
    result = null_model_b(n_permutations=100, seed=42)
    assert 'p_value_serine_unique' in result

def test_null_model_c_returns_results():
    """Tests all 24 encodings. No permutations needed."""
    result = null_model_c()
    assert 'n_encodings' in result
    assert result['n_encodings'] == 24
    assert 'twofold_invariant_count' in result
    assert 'encoding_results' in result
```

**Step 2: Run to verify failure**

Run: `pytest tests/test_null_models.py -v`

**Step 3: Implement null_models.py**

```python
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
    """Extract degeneracy counts: how many AAs have each codon count."""
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
        'n_disconnected': len(disconnected),
        'exactly_one_disconnected': len(disconnected) == 1,
    }


def null_model_a(
    n_permutations: int = 100_000,
    seed: int | None = None,
) -> dict:
    """Null Model A: fix degeneracy structure, randomly assign codons.

    Keeps the same number of AAs with 1, 2, 3, 4, 6 codons but
    randomly assigns which codons go to which AA.
    """
    rng = random.Random(seed)

    # Get the degeneracy structure from the standard code
    aa_sizes = _degeneracy_structure(STANDARD)
    sense_codons = [c for c in ALL_CODONS if STANDARD[c] != 'Stop']
    stop_codons = [c for c in ALL_CODONS if STANDARD[c] == 'Stop']

    # Observed values
    obs = _score_code(STANDARD)

    count_serine_unique = 0
    count_bit5_uniform = 0

    for _ in range(n_permutations):
        shuffled = list(sense_codons)
        rng.shuffle(shuffled)

        # Assign codons to synthetic AAs with the same degeneracy sizes
        code = {c: 'Stop' for c in stop_codons}
        idx = 0
        for i, size in enumerate(aa_sizes):
            aa_name = f'AA{i}'
            for c in shuffled[idx:idx + size]:
                code[c] = aa_name
            idx += size

        score = _score_code(code)
        if score['exactly_one_disconnected']:
            count_serine_unique += 1
        if score['twofold_all_pass']:
            count_bit5_uniform += 1

    return {
        'n_permutations': n_permutations,
        'p_value_serine_unique': count_serine_unique / n_permutations,
        'p_value_bit5_uniform': count_bit5_uniform / n_permutations,
        'observed_exactly_one_disconnected': obs['exactly_one_disconnected'],
        'observed_twofold_all_pass': obs['twofold_all_pass'],
    }


def null_model_b(
    n_permutations: int = 100_000,
    seed: int | None = None,
) -> dict:
    """Null Model B: preserve block structure, shuffle block-to-AA assignment.

    A 'block' is a group of codons sharing the first two nucleotide positions
    (16 blocks of 4 codons each).
    """
    rng = random.Random(seed)

    # Build blocks: groups by first-2-base prefix
    blocks: dict[str, list[str]] = defaultdict(list)
    for codon in ALL_CODONS:
        blocks[codon[:2]].append(codon)
    block_list = list(blocks.values())  # 16 blocks

    # Standard: which blocks map to which AA pattern
    # Get the AA assignment per block in standard code
    std_block_aas = []
    for block_codons in block_list:
        aas = [STANDARD[c] for c in block_codons]
        std_block_aas.append(tuple(aas))

    obs = _score_code(STANDARD)
    count_serine_unique = 0

    for _ in range(n_permutations):
        # Shuffle which block gets which AA assignment pattern
        shuffled_patterns = list(std_block_aas)
        rng.shuffle(shuffled_patterns)

        code = {}
        for block_codons, pattern in zip(block_list, shuffled_patterns):
            for codon, aa in zip(block_codons, pattern):
                code[codon] = aa

        score = _score_code(code)
        if score['exactly_one_disconnected']:
            count_serine_unique += 1

    return {
        'n_permutations': n_permutations,
        'p_value_serine_unique': count_serine_unique / n_permutations,
    }


def null_model_c() -> dict:
    """Null Model C: test all 24 base-to-bit encodings.

    Determines whether filtration properties are encoding-dependent
    or encoding-invariant.
    """
    encodings = all_encodings()
    results = []
    twofold_invariant = 0

    for enc in encodings:
        score = _score_code(STANDARD, encoding=enc)
        results.append({
            'encoding': enc,
            'twofold_all_pass': score['twofold_all_pass'],
            'exactly_one_disconnected': score['exactly_one_disconnected'],
            'n_disconnected': score['n_disconnected'],
        })
        if score['twofold_all_pass']:
            twofold_invariant += 1

    return {
        'n_encodings': len(encodings),
        'twofold_invariant_count': twofold_invariant,
        'encoding_results': results,
    }
```

**Step 4: Run tests**

Run: `pytest tests/test_null_models.py -v`

**Step 5: Commit**

```bash
git add src/codon_topo/analysis/null_models.py tests/test_null_models.py
git commit -m "feat: null models A, B, C for statistical validation"
```

---

### Task 9: Package integration and coverage check

**Step 1: Re-export public API in `__init__.py`**

Update `src/codon_topo/__init__.py`:
```python
"""Codon Geometry Validation & Prediction Engine."""
from codon_topo.core.encoding import codon_to_vector, hamming_distance, ALL_CODONS
from codon_topo.core.genetic_codes import STANDARD, get_code, all_table_ids
from codon_topo.core.filtration import analyze_filtration
from codon_topo.core.homology import disconnection_catalogue, persistent_homology
from codon_topo.core.embedding import embed_codon
from codon_topo.core.fano import is_fano_line, fano_partner
```

**Step 2: Run full test suite with coverage**

Run: `pytest --cov=codon_topo --cov-report=term-missing -v`

Inspect any uncovered lines. Add targeted tests if coverage < 100% on core modules.

**Step 3: Run preliminary scripts as cross-check**

Run each preliminary script and visually confirm output hasn't changed. These are the oracle.

```bash
python3 codon_topo_preliminary/verify_filtration.py > /tmp/old_filtration.txt 2>&1
python3 codon_topo_preliminary/verify_kras.py > /tmp/old_kras.txt 2>&1
```

**Step 4: Commit**

```bash
git add src/codon_topo/__init__.py
git commit -m "feat: public API re-exports in package __init__"
```

---

### Task 10: Update CLAUDE.md with actual commands

**Step 1: Update CLAUDE.md**

Update the Build & Test Commands section to match the actual project structure (src layout, confirmed commands).

**Step 2: Commit**

```bash
git add CLAUDE.md
git commit -m "docs: update CLAUDE.md with verified build/test commands"
```
