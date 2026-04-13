# Codon-Topo: Final Results Summary

**Date**: April 13, 2026 | **Version**: 0.2.0 | **Tests**: 367 passing | **Seed**: 135325

## Dataset

19 organisms total, 6 verified by tRNAscan-SE 2.0.12 on NCBI genome assemblies.
11 pairings (5 independent). 3-tier architecture: large-genome eukaryotes / small-genome eukaryotes / bacteria.

### tRNAscan-SE verified organisms

| Organism | Code | Assembly | Reassigned AA | Count | Total |
|----------|------|----------|--------------|-------|-------|
| T. thermophila | Table 6 | GCF_000189635.1 | Gln | **54** (15+39 supp) | 718 |
| P. tetraurelia | Table 6 | GCF_000165425.1 | Gln | **18** (7+11 supp) | 216 |
| O. trifallax | Table 6 | GCA_000295675.1 | Gln | **8** (2+6 supp) | 94 |
| S. coeruleus | Standard | GCA_001970955.1 | — | Gln=11 | 272 |
| I. multifiliis | Standard | GCF_000220395.1 | — | Gln=3 | 150 |
| Blastocrithidia | Table 31 | GCA_000436035.1 | Trp | **0** | 34 |

## Finding 1: Hypercube Coloring Optimality

- Observed score: 13,477 | Null: 14,954 +/- 628
- z-score: **-2.35** | Quantile: **0.6%** | p = **0.006** (conservative)
- 95% bootstrap CI for quantile: [0.2%, 1.2%]
- Robust across all rho in [0,1] (all p < 0.01)
- Score decomposition: Position 2 = 49.3%, Position 1 = 38.2%, Wobble = 12.5%

## Finding 2: Per-Table Optimality Preservation

- **24/25** tables significant at p < 0.05 (raw AND BH-FDR corrected)
- Mean quantile: 1.4% | Only Table 3 (yeast mito) fails at 7.4%
- FDR correction applied per crush/GLM-5 review recommendation

## Finding 3: Topology Avoidance

- Observed: 6/27 (22.2%) create new disconnections
- Possible: 931/1,280 (72.7%) create new disconnections
- Hypergeometric p = **4.8 x 10^-8** (correct finite-landscape null per codex review)
- Caveat: events share phylogenetic structure; p should be interpreted as depletion effect size

## Finding 4: tRNA Enrichment

- 11 pairings across 4 eukaryotic supergroups (Alveolata, Opisthokonta, Archaeplastida, Excavata)
- 6 tRNAscan-SE verified organisms from NCBI genome assemblies

| Test | n | Statistic | P-value |
|------|---|-----------|---------|
| Sign test | 11 | 8/11 elevated | 0.113 |
| Fisher+Stouffer (all) | 11 | Z=3.47 | **2.6 x 10^-4** |
| Fisher+Stouffer (independent) | 5 | Z=1.84 | **0.033** |
| AA-label exact+Stouffer | 11 | Z=4.64 | **1.7 x 10^-6** |

### Boundary conditions (informative negatives)

| Organism | Domain | Mechanism | tRNA duplication? |
|----------|--------|-----------|-------------------|
| Ciliates (large nuclear genomes) | Eukaryota | Gene duplication | **YES** — 39 suppressor tRNAs in Tetrahymena |
| Blastocrithidia (3 Mb genome) | Eukaryota | tRNA import | **NO** — 0 Trp tRNAs |
| Mycoplasma (0.8 Mb genome) | Bacteria | Anticodon modification | **NO** — 1 Trp tRNA reads both UGA+UGG |

**Interpretation**: tRNA gene duplication accompanies codon reassignment specifically in organisms with large nuclear genomes that can tolerate gene amplification. Streamlined genomes use alternative mechanisms (tRNA import, anticodon modification).

### Caveats (from 4-model code review)

1. **Stouffer on discrete p-values**: AA-label exact test has p_min=0.05 (1/20 AAs). Alternative interpretation: reassigned AA ranks #1 in 8/11 pairings, binomial p = 0.0002 under p0=0.05. [crush review]
2. **Non-independence**: Some pairings share control organisms. Independent-only Stouffer (n=5, p=0.033) is the conservative primary test. [codex review]
3. **Estimated counts**: Euplotes and Blepharisma tRNA counts are literature estimates, not tRNAscan-SE verified. [codex review]

## Clean Negatives

| Analysis | Result | Interpretation |
|----------|--------|---------------|
| Local mismatch cost | p=0.70 | Reassignment NOT driven by escaping bad neighborhoods |
| KRAS-Fano | p=1.0 | XOR structure has no clinical prediction power |
| Depth calibration | rho=0.0 | No monotonic epsilon-age relationship |
| Bit-position bias (deduplicated) | p=0.075 | Marginal after correcting for non-independence |
