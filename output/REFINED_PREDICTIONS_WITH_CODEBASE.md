# Refined Testable Predictions: Mapped to Codebase Modules

**Author**: Sergey Kornilov
**Date**: April 13, 2026
**Purpose**: Concrete testable predictions with specific refinements to existing code. Each prediction maps to the file(s) that need to be modified, the new code to add, and the literature priors.

---

## Executive Summary of Codebase State

```
CURRENT CODEBASE (280 tests passing, 96% coverage)
==================================================
src/codon_topo/core/
  encoding.py              3 fn   GF(2)^6 encoding + 24 permutations
  filtration.py            5 fn   Two-fold/four-fold checks
  homology.py              4 fn   Union-find components
  fano.py                  5 fn   XOR triples
  embedding.py             2 fn   Coordinate-wise map to C^3
  genetic_codes.py         4 fn   All 25 NCBI tables

src/codon_topo/analysis/
  null_models.py           6 fn   Models A, B, C
  reassignment_db.py       7 fn   Events + bit_position_bias()
  depth_calibration.py     4 fn   Epsilon vs age
  synbio_feasibility.py    3 fn   Structure-preservation scoring
  cosmic_query.py          7 fn   cBioPortal KRAS (negative)

MISSING (needed for refined predictions):
  * Grantham distance matrix data + loader
  * Hypercube coloring optimality module
  * tRNA database interface (GtRNAdb)
  * Position-weighted null models (D, E, F)
  * Wobble-compatibility scoring
  * Ambiguous intermediate simulation
```

---

## Literature Priors (Tavily-confirmed, April 13, 2026)

| Finding | Source | Implication |
|---------|--------|-------------|
| Freeland-Hurst 1998 ("code is one in a million") PMID 9732450 | Tavily confirmed | Direct prior for hypercube coloring theorem |
| Grantham (1974) matrix available as R package `grantham` | pattern.institute/grantham | Data source for Monte Carlo |
| CUG clade: "Evolutionary instability of CUG-Leu" Nature Commun 2018 | nature.com/articles/s41467-018-04374-7 | Critical for P2.2 (CUG test case) |
| Yeast mito: 2 tRNA copies (ancestral Thr + new from Leu) | PNAS 1200109109 | Direct evidence for P1.3 |
| GtRNAdb: active database at gtrnadb.ucsc.edu | Confirmed | Data source for tRNA counts |
| Ts/Tv in mitochondrial ≫ nuclear, 3rd position biased | Yang-Yoder 1999, PMC1207285 | Priors for P1.2 null model |
| Ambiguous intermediate theory (Schultz-Yarus 1994) | PMID 7666454 | Mechanism framework for P2.5 |
| Knight-Freeland-Landweber 2001 PMC3293468 | Confirmed | Review of code evolution |

---

## TIER 1: IMMEDIATE REFINEMENTS (this week, pre-meeting)

### P1.1 — Per-encoding min-distance emission (CRITICAL CORRECTION)

**File to modify**: `src/codon_topo/analysis/null_models.py`

**Current state**: `null_model_c()` returns `{twofold_all_pass, fourfold_all_pass, exactly_one_disconnected, n_disconnected, disconnected_aas}` per encoding. Does NOT emit per-AA minimum inter-block Hamming distance.

**Refinement needed**:

```python
# Add to null_models.py
def _score_code_detailed(code, encoding=None) -> dict:
    """Enhanced score that emits per-AA minimum inter-block Hamming distance."""
    vec_map = _get_vectors(encoding)
    groups = defaultdict(list)
    for codon, aa in code.items():
        if aa != "Stop":
            groups[aa].append((codon, vec_map[codon]))
    
    per_aa_info = {}
    for aa, codon_vecs in groups.items():
        vecs = [v for _, v in codon_vecs]
        components = connected_components(vecs, epsilon=1)
        if components > 1:
            # Compute inter-component min distance by brute force
            # (acceptable for <= 8 codons)
            blocks = _find_blocks(codon_vecs, epsilon=1)
            min_inter = min(
                hamming_distance(v1, v2)
                for b1, b2 in combinations(blocks, 2)
                for v1 in b1 for v2 in b2
            )
            per_aa_info[aa] = {
                "n_components": components,
                "min_inter_distance": min_inter,
                "reconnect_eps": _find_reconnect_eps(vecs)
            }
    return per_aa_info

def null_model_c_extended(code=None) -> dict:
    """Null Model C with per-encoding disconnection details."""
    # ... (extends existing null_model_c)
    for enc_idx, enc in enumerate(all_encodings()):
        details[enc_idx] = _score_code_detailed(code, enc)
    return {...}
```

**Tests to add** in `tests/test_null_models.py`:
```python
def test_serine_min_distance_varies_across_encodings():
    """Verify the kimi-k2.5 counterexample: Ser min distance 2 in some encodings."""
    result = null_model_c_extended()
    min_distances_for_ser = [
        enc_info.get("Ser", {}).get("min_inter_distance")
        for enc_info in result["per_encoding_details"]
    ]
    # Should find at least one encoding with min=2
    assert 2 in min_distances_for_ser
    # Should find encoding with min=4 (default)
    assert 4 in min_distances_for_ser
```

**Time**: 1-2 hours.

**Publication impact**: Prevents the overclaim. Makes the paper honest.

---

### P1.2 — Position-weighted null for bit-position bias

**File to modify**: `src/codon_topo/analysis/reassignment_db.py`

**Current state**: `bit_position_bias()` uses `scipy.stats.chisquare()` with uniform null (expected = total/6 per position).

**Refinement needed**:

```python
# Add to reassignment_db.py

# Position-specific Ts/Tv priors from Yang-Yoder 1999 (nuclear + mito aggregate)
# These are PRIORS - the actual values can be refined per lineage
TS_TV_PRIORS = {
    "nuclear": {
        # position 1: somewhat biased, modest Ts/Tv
        "pos1_ts_tv": 2.5,
        # position 2: most conserved (strongest selection)
        "pos2_ts_tv": 2.8,
        # position 3: most variable, wobble-driven
        "pos3_ts_tv": 4.5,
    },
    "mitochondrial": {
        "pos1_ts_tv": 8.0,
        "pos2_ts_tv": 6.0,
        "pos3_ts_tv": 15.0,
    }
}

def expected_bit_flips_under_neutral(
    events: list[ReassignmentEvent],
    compartment: str = "mitochondrial",
) -> list[float]:
    """Compute expected bit-flip count per bit position under Ts/Tv-weighted neutral null.
    
    Accounts for the empirical fact that transitions (Ts) are much more common
    than transversions (Tv), especially at 3rd position and in mitochondrial DNA.
    """
    priors = TS_TV_PRIORS[compartment]
    total_flips = sum(_count_flips_per_event(e) for e in events)
    # Rough weights: bit 0,1 from pos1; bit 2,3 from pos2; bit 4,5 from pos3
    weights = [
        priors["pos1_ts_tv"] / (priors["pos1_ts_tv"] + 1),  # bit 0 (Ts component)
        1 / (priors["pos1_ts_tv"] + 1),                      # bit 1 (Tv component)
        priors["pos2_ts_tv"] / (priors["pos2_ts_tv"] + 1),  # bit 2
        1 / (priors["pos2_ts_tv"] + 1),                      # bit 3
        priors["pos3_ts_tv"] / (priors["pos3_ts_tv"] + 1),  # bit 4
        1 / (priors["pos3_ts_tv"] + 1),                      # bit 5
    ]
    total_w = sum(weights)
    return [w / total_w * total_flips for w in weights]


def bit_position_bias_weighted(compartment: str = "mitochondrial") -> dict:
    """Position-weighted chi-square test for bit-position bias.
    
    Replaces the uniform null with a realistic Ts/Tv-weighted expectation.
    """
    db = build_reassignment_db()
    sense_events = [e for e in db if e.source_aa != "Stop" and e.target_aa != "Stop"]
    observed = _count_flips_per_bit(sense_events)  # list of 6 counts
    expected = expected_bit_flips_under_neutral(sense_events, compartment)
    
    from scipy.stats import chisquare
    stat, p = chisquare(observed, f_exp=expected)
    return {
        "bit_counts_observed": observed,
        "bit_counts_expected_weighted": expected,
        "chi2_statistic": float(stat),
        "chi2_p_value_weighted": float(p),
        "chi2_p_value_uniform": bit_position_bias()["chi2_p_value"],  # for comparison
        "compartment": compartment,
    }
```

**Tests to add**:
```python
def test_bit_position_bias_under_weighted_null():
    """Key test: does bias signal survive realistic null?"""
    uniform_result = bit_position_bias()
    weighted_mito = bit_position_bias_weighted("mitochondrial")
    weighted_nuc = bit_position_bias_weighted("nuclear")
    
    # Document the result regardless
    print(f"Uniform p = {uniform_result['chi2_p_value']}")
    print(f"Mito-weighted p = {weighted_mito['chi2_p_value_weighted']}")
    print(f"Nuclear-weighted p = {weighted_nuc['chi2_p_value_weighted']}")
    # Assertion depends on actual outcome - run it and see
```

**Time**: 2-3 hours.

**Publication impact**: DECIDES whether B4 is a real finding or a mutational-spectrum artifact.

---

### P1.3 — tRNA gene duplication correlation

**New file**: `src/codon_topo/analysis/trna_evidence.py`

**Current state**: No tRNA data in the codebase.

**Refinement needed**:

```python
"""tRNA evidence for variant-code disconnections.

Queries GtRNAdb and mitotRNAdb to test whether organisms with codon-graph 
disconnections show evidence of tRNA gene duplication or recruitment.

Confirmed prior (Tavily 2026-04-13):
- Yeast mito acquired tRNA^Thr derived from tRNA^His (PMC3113583)
- This is exactly what B3 predicts.
"""
from dataclasses import dataclass
import requests
from typing import Optional

@dataclass(frozen=True)
class TRNARepertoire:
    organism: str
    compartment: str  # "nuclear" or "mitochondrial"
    by_amino_acid: dict[str, int]  # AA -> gene count
    by_anticodon: dict[str, int]  # anticodon -> gene count
    source: str

# Pre-fetched/curated data for known disconnection cases
# These can be fetched from GtRNAdb gene lists (available as downloadable TSV)
CURATED_REPERTOIRES = {
    # Disconnection organisms
    "Saccharomyces_cerevisiae_mito": TRNARepertoire(
        organism="Saccharomyces cerevisiae",
        compartment="mitochondrial",
        by_amino_acid={"Thr": 2, "Leu": 0, ...},  # 2 tRNA^Thr (novel + ancestral)
        by_anticodon={"UGU": 1, "UAG": 1, ...},
        source="PNAS 1200109109 (Su et al. 2011)",
    ),
    "Saccharomyces_cerevisiae_nuclear": TRNARepertoire(
        organism="Saccharomyces cerevisiae",
        compartment="nuclear",
        by_amino_acid={"Thr": 11, "Leu": 10, ...},  # Standard distribution
        by_anticodon={...},
        source="GtRNAdb",
    ),
    # Controls (no disconnection): close relatives with standard code
    "Lachancea_thermotolerans_nuclear": TRNARepertoire(...),
    # ... all four disconnection organisms + controls
}


def test_trna_duplication_correlates_with_disconnection() -> dict:
    """Fisher's exact test: are disconnection organisms enriched for tRNA duplications?"""
    # For each disconnection organism:
    #   - Count tRNA gene copies for the reassigned AA
    #   - Compare to matched control (same clade, standard code)
    # 2x2 contingency: (disconnection, control) x (excess tRNAs, no excess)
    
    disconnection_organisms = [
        ("S.cerevisiae_mito", "Thr"),       # yeast mito
        ("S.obliquus_mito", "Leu"),         # chlorophycean
        ("P.tannophilus_nuclear", "Ala"),   # Pachysolen
        ("C.albicans_nuclear", "Ser"),      # Candida
    ]
    
    results = []
    for org, reassigned_aa in disconnection_organisms:
        excess = _compare_to_nearest_standard_relative(org, reassigned_aa)
        results.append({
            "organism": org,
            "reassigned_aa": reassigned_aa,
            "disconnection_aa_trna_count": excess["observed"],
            "control_aa_trna_count": excess["control"],
            "ratio": excess["observed"] / max(excess["control"], 1),
        })
    
    return {
        "organisms_tested": results,
        "hypothesis": "Disconnection organisms have extra tRNA genes for reassigned AA",
        "prior_evidence": "Confirmed for yeast mito Thr (PMC3113583)",
    }
```

**Tests to add**: `tests/test_trna_evidence.py` with fixture data.

**Time**: 4-6 hours (mostly data retrieval from GtRNAdb/mitotRNAdb).

**Publication impact**: If positive — this becomes the PRIMARY biological finding of the paper.

---

### P1.4 — Hypercube Coloring Optimality (the new direction)

**New file**: `src/codon_topo/analysis/coloring_optimality.py`

**Current state**: Doesn't exist.

**New data needed**: `src/codon_topo/data/grantham.json` (20×20 Grantham distance matrix).

**Refinement**:

```python
"""Hypercube coloring optimality: the new publishable direction.

Based on Freeland-Hurst 1998 (PMID 9732450: 'genetic code is one in a million')
and gemini3/glm-5 proposed theorem. Tests whether the standard code is a locally
optimal coloring of Q_6 by 22 AA classes (+ stop), minimizing physicochemical 
mismatch across 1-Hamming edges.
"""

import json
import random
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.stats import percentileofscore

from codon_topo.core.encoding import ALL_CODONS, codon_to_vector, hamming_distance
from codon_topo.core.genetic_codes import STANDARD

# Load Grantham (1974) distance matrix
# Data from: pattern.institute/grantham
def _load_grantham() -> dict[tuple[str, str], float]:
    data_path = Path(__file__).parent.parent / "data" / "grantham.json"
    with open(data_path) as f:
        raw = json.load(f)
    # Symmetric matrix
    out = {}
    for a, b_dict in raw.items():
        for b, d in b_dict.items():
            out[(a, b)] = float(d)
            out[(b, a)] = float(d)
    for a in raw:
        out[(a, a)] = 0.0
    return out


def hypercube_edge_mismatch_score(
    code: dict[str, str],
    distance: str = "grantham",
) -> float:
    """Compute total physicochemical mismatch across all Hamming-1 edges.
    
    F(code) = sum over {v,w} with d(v,w)=1: delta(color(v), color(w))
    where delta is physicochemical distance (Grantham) or 1 if different, 0 if same.
    """
    if distance == "grantham":
        dist = _load_grantham()
    else:
        dist = {("X", "Y"): 1.0 for X in set(code.values()) for Y in set(code.values())}
    
    total = 0.0
    codon_vecs = [(c, codon_to_vector(c)) for c in ALL_CODONS]
    for (c1, v1), (c2, v2) in combinations(codon_vecs, 2):
        if hamming_distance(v1, v2) == 1:
            aa1, aa2 = code[c1], code[c2]
            if aa1 != aa2:  # Only non-synonymous edges contribute
                # Handle "Stop" specially
                if aa1 == "Stop" or aa2 == "Stop":
                    total += 200.0  # Large penalty for stop-sense adjacency
                else:
                    total += dist.get((aa1, aa2), 200.0)
    return total


def monte_carlo_null(
    n_samples: int = 100_000,
    preserve_degeneracy: bool = True,
    seed: int | None = None,
) -> dict:
    """Generate null distribution of coloring scores with matched class sizes.
    
    preserve_degeneracy=True: same block-size distribution as standard code
        (3x6-fold, 5x4-fold, 1x3-fold, 9x2-fold, 2x1-fold)
    """
    rng = random.Random(seed)
    scores = []
    
    # Extract class sizes from STANDARD
    from collections import defaultdict
    aa_counts = defaultdict(int)
    for aa in STANDARD.values():
        aa_counts[aa] += 1
    class_sizes = sorted(aa_counts.values(), reverse=True)
    
    aa_labels = list(aa_counts.keys())
    sense_codons = [c for c in ALL_CODONS if STANDARD[c] != "Stop"]
    stop_codons = [c for c in ALL_CODONS if STANDARD[c] == "Stop"]
    
    for _ in range(n_samples):
        shuffled = list(sense_codons)
        rng.shuffle(shuffled)
        
        perm_code = {c: "Stop" for c in stop_codons}
        idx = 0
        for aa, size in zip(aa_labels, class_sizes):
            for c in shuffled[idx:idx + size]:
                perm_code[c] = aa
            idx += size
        
        scores.append(hypercube_edge_mismatch_score(perm_code))
    
    obs_score = hypercube_edge_mismatch_score(STANDARD)
    quantile = percentileofscore(scores, obs_score, kind="weak")
    
    return {
        "observed_score": obs_score,
        "null_mean": np.mean(scores),
        "null_std": np.std(scores),
        "quantile_of_observed": quantile,  # Lower is better for mismatch score
        "p_value": quantile / 100.0,
        "n_samples": n_samples,
    }


def affine_subspace_constrained_null(n_samples: int = 10_000) -> dict:
    """Stronger null: sample only colorings where synonymous blocks 
    are affine subspaces of dim <= 2.
    
    Tests whether within the algebraic constraint, the specific AA 
    assignment is optimal (addresses gemini3's proposed theorem).
    """
    # This is more involved - enumerate all valid affine-subspace 
    # decompositions of 64 sense codons with the right block sizes
    # ... [stub implementation]
    pass


def cross_table_optimality() -> dict:
    """Compare optimality quantile across all 25 NCBI tables.
    
    Do variant codes maintain similar optimality, or does reassignment 
    trade off error-correction for decoding flexibility?
    """
    from codon_topo.core.genetic_codes import get_code, all_table_ids
    results = {}
    for tid in all_table_ids():
        code = get_code(tid)
        results[tid] = {
            "score": hypercube_edge_mismatch_score(code),
            "name": get_code(tid).__class__.__name__,
        }
    return results
```

**Tests** `tests/test_coloring_optimality.py`:
```python
def test_hypercube_score_is_deterministic():
    s1 = hypercube_edge_mismatch_score(STANDARD)
    s2 = hypercube_edge_mismatch_score(STANDARD)
    assert s1 == s2

def test_standard_code_is_in_top_10_percent():
    """Replicate Freeland-Hurst 1998 finding at GF(2)^6 level."""
    result = monte_carlo_null(n_samples=10_000, seed=42)
    assert result["quantile_of_observed"] < 10.0  # Top 10% optimality
```

**Time**: 1-2 weeks for full implementation + paper-ready writeup.

**Publication impact**: HIGHEST — this is the central claim of the reframed paper.

---

## TIER 2: FOR PAPER 1 (next 3 months)

### P2.1 — Null Models D, E, F (proper nulls)

**File to modify**: `src/codon_topo/analysis/null_models.py`

**Refinements**:

```python
def null_model_d_wobble_preserving(n_permutations: int = 100_000, ...) -> dict:
    """Null D: Shuffle only within wobble-equivalence classes.
    
    Preserves the biological fact that 3rd-position degeneracy is the primary 
    source of synonymy. Tests whether filtration properties hold beyond this.
    """
    # Group codons by first-two-base prefix (16 blocks)
    # Within each block, shuffle AA assignment
    # ... this is actually Model B with exclude_stops=False, but make it 
    # explicit with different naming
    
def null_model_e_ts_tv_weighted(n_permutations: int = 100_000, ...) -> dict:
    """Null E: Sample codes weighted by position-specific Ts/Tv bias.
    
    Uses TS_TV_PRIORS from reassignment_db.py to make permutations 
    biologically realistic.
    """
    
def null_model_f_block_preserving_strict(n_permutations: int = 100_000, ...) -> dict:
    """Null F: Preserve BOTH block structure AND within-block AA identity 
    at third position, shuffle only first-two-position blocks.
    
    Hardest null - tests whether anything beyond wobble rules is special.
    """
```

**Time**: 1-2 weeks.

---

### P2.2 — CUG clade specific analysis

**File to modify**: `src/codon_topo/analysis/reassignment_db.py`

**Refinement**:

```python
def cug_clade_analysis() -> dict:
    """Test the CUG clade reassignment patterns.
    
    Prior (Nature Commun 2018): CUG reassignment is evolutionarily 
    unstable - different lineages reassign to Ser, Ala, or revert to Leu.
    
    Hypothesis: The algebraic "reach" of each destination AA predicts 
    which reassignment is evolutionarily stable.
    """
    # CUG's fresh location in each of the 24 encodings
    # For each possible destination (Ser, Ala, Leu), compute:
    #   - Resulting filtration status
    #   - Nearest wobble-compatible codon
    #   - Predicted stability
    cug_variants = {
        "CUG_to_Ser": {"table_id": 12, "organism": "Candida albicans"},
        "CUG_to_Ala": {"table_id": 26, "organism": "Pachysolen tannophilus"},
        "CUG_to_Leu": {"table_id": 1, "organism": "Standard (ancestral)"},
    }
    # ... compute metrics per variant
    return {...}
```

**Tests**: fold into `test_reassignment_db.py`.

**Time**: 3-5 days.

---

### P2.3 — Forward evolution simulation

**New file**: `src/codon_topo/analysis/code_evolution.py`

**Refinement**:

```python
"""Forward simulation of genetic code evolution under competing models.

Compares which model best reproduces the observed 25 NCBI tables:
  - Random reassignment
  - Wobble-weighted only
  - Algebra-preserving only
  - Wobble + algebra combined
  - tRNA-network-constrained (requires trna_evidence.py)
"""

def simulate_evolution(
    starting_code: dict[str, str] = STANDARD,
    n_steps: int = 50,
    model: str = "wobble",  # "random", "wobble", "algebra", "combined", "trna"
    n_replicates: int = 10_000,
    seed: int | None = None,
) -> dict:
    """Run forward simulations and return distribution of final codes."""
    # Use pre-computed transition probabilities per model
    # ... simulation logic
    
def compare_models_against_ncbi() -> dict:
    """Likelihood ratio test: which model best fits observed 25 tables?"""
    observed = [get_code(tid) for tid in all_table_ids()]
    for model in ["random", "wobble", "algebra", "combined", "trna"]:
        sim_dist = simulate_evolution(model=model)
        likelihood = _compute_likelihood(observed, sim_dist)
        # ...
    return {...}
```

**Time**: 3-4 weeks.

---

### P2.4 — Ambiguous intermediate simulation

**New file**: `src/codon_topo/analysis/ambiguous_intermediate.py`

**Refinement**:

```python
"""Ambiguous intermediate states in codon reassignment.

Based on Schultz-Yarus 1994 (PMID 7666454): codon reassignment proceeds 
through a state where a codon is read ambiguously by multiple tRNAs.

Simulates the ambiguity-to-clean-reassignment transition to predict 
which reassignments can occur.
"""

def simulate_ambiguous_intermediate(
    codon: str,
    current_aa: str,
    candidate_aa: str,
    trna_pool: dict,
) -> dict:
    """Model codon reassignment via ambiguous intermediate.
    
    Input: codon, current AA, candidate new AA, available tRNA pool.
    Output: probability of successful reassignment, expected ambiguity duration,
    fitness cost during ambiguous phase.
    """
    # ... Markov chain or similar model
    pass
```

**Time**: 2-3 weeks.

---

## TIER 3: FUTURE PAPERS / GRANT HORIZON

### P3.1 — Experimental orthogonal translation fidelity test
Requires Chin/Isaacs lab collaboration. Not a codebase refinement — but the codebase should PRODUCE the predictions:

**File to add**: `src/codon_topo/analysis/synbio_predictions.py`

```python
def generate_synbio_test_predictions() -> list[dict]:
    """Generate 10 matched-pair reassignments for experimental test.
    
    Uses P1.4 (coloring optimality) + P1.3 (tRNA evidence) to rank 
    candidate reassignments. Outputs JSON specification for wet lab.
    """
```

---

## REFINEMENT PRIORITY MATRIX (mapped to codebase)

```
  EASY (existing files, minor additions)     HARDER (new modules needed)
  ════════════════════════════════════       ══════════════════════════
  
T P1.1 null_models.py:                       P1.3 NEW trna_evidence.py
I   +null_model_c_extended()                 P1.4 NEW coloring_optimality.py
E P1.2 reassignment_db.py:                        + data/grantham.json
R   +bit_position_bias_weighted()            
                                             
1                                             
─────────────────────────────────────────────────────────────────────────
T P2.2 reassignment_db.py:                   P2.1 null_models.py:
I   +cug_clade_analysis()                     +null_model_d,e,f
E                                             P2.3 NEW code_evolution.py
R                                             P2.4 NEW ambiguous_intermediate.py
                                             
2                                             
─────────────────────────────────────────────────────────────────────────
T (no codebase work)                         P3.1 NEW synbio_predictions.py
I                                             
E                                             
R                                             
                                             
3                                             
```

---

## IMPLEMENTATION SEQUENCE (concrete next steps)

### Week 1 (Before meeting, April 14-16)
```bash
# Day 1: Fix the critical overclaim
1. Edit src/codon_topo/analysis/null_models.py
   → Add null_model_c_extended() with per-encoding min distance
   → Add reconnection epsilon tracking
2. Edit tests/test_null_models.py
   → Add test verifying counterexample
3. Run: python3.11 -m pytest tests/test_null_models.py -v
4. Update disconnection_catalogue.csv with per-encoding info

# Day 2: Position-weighted null for bit-position bias
1. Edit src/codon_topo/analysis/reassignment_db.py  
   → Add TS_TV_PRIORS
   → Add bit_position_bias_weighted()
2. Edit tests/test_reassignment_stats.py
3. Document result (p-value survival or not)

# Day 3-4: tRNA evidence for B3
1. Create src/codon_topo/analysis/trna_evidence.py
2. Curate data for 4 disconnection organisms from GtRNAdb
3. Implement Fisher's exact test
4. Create tests/test_trna_evidence.py

# Day 5: Start hypercube coloring
1. Create src/codon_topo/data/grantham.json from pattern.institute/grantham
2. Create src/codon_topo/analysis/coloring_optimality.py stub
3. Add hypercube_edge_mismatch_score()
4. Add test: standard code < mean of 1000 random codes
```

### Weeks 2-3
```
Complete hypercube coloring optimality module with full Monte Carlo.
Produce figure: quantile of standard code across NCBI tables.
```

### Weeks 4-8
```
Null models D, E, F.
CUG clade analysis.
Forward evolution simulation.
Draft paper.
```

---

## WHAT GETS CUT FROM THE PAPER

Based on this refinement plan, the following should be REMOVED from the paper scope:

1. ❌ "Holomorphic embedding" analysis — M5 terminology broken
2. ❌ "Fano line" language — M4 terminology broken  
3. ❌ PSL(2,7) / Clayworth Algebra — M7 already rejected
4. ❌ KRAS/Fano clinical predictions — B7 dead
5. ❌ Depth calibration with n=6 — B6 non-result
6. ❌ 83% structural robustness as-is — B5 encoding artifact
7. ⚠️ Bit-5 and prefix filtrations as DISCOVERIES — reframe as "expected consequences of chosen encoding"
8. ⚠️ Serine min distance 4 as INVARIANT — reframe as "ε=1 disconnection holds; min distance is encoding-dependent"

WHAT STAYS:
1. ✅ Systematic NCBI 25-table survey (methods)
2. ✅ Encoding-sensitivity analysis (Null Model C extended)
3. ✅ Bit-position bias UNDER PROPER NULL (P1.2)
4. ✅ tRNA duplication correlation (P1.3 — THE primary finding)
5. ✅ Hypercube coloring optimality (P1.4 — THE central theorem)
6. ✅ Clean KRAS/Fano negative (P as cautionary tale)
7. ✅ Open-source reproducible pipeline

---

## TAVILY RESEARCH REPORTS STATUS

Deep async research via `tvly.research()` submitted April 13, 2026 — all 5 queries returning `status: pending`. These reports complete asynchronously (minutes to hours). Regular search API (`tvly.search()`) used successfully for priors.

All key priors obtained via search:
- Freeland-Hurst 1998: PMID 9732450 ✓
- Grantham 1974: R package `grantham` at pattern.institute ✓
- Yeast mito Thr: PNAS 1200109109 (2 tRNA copies confirmed) ✓
- GtRNAdb at gtrnadb.ucsc.edu ✓
- Ts/Tv priors: Yang-Yoder 1999 PMC1207285 ✓
- Ambiguous intermediate: Schultz-Yarus 1994 PMID 7666454 ✓
- CUG instability: Nature Commun 2018 s41467-018-04374-7 ✓

Research reports can be polled later; not blocking the analysis.

---

*This document maps every refined prediction to specific files and functions in the codebase. Implementation estimates are based on existing patterns (280 tests, consistent structure). The "what gets cut" section is informed by consensus from 10 LLM models across 12 evaluation passes.*
