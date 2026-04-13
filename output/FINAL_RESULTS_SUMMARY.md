# Codon-Topo: Final Results Summary

**Date**: April 13, 2026 | **Version**: 0.3.0 | **Tests**: 367 passing | **Seed**: 135325

## Dataset

21 organisms total, 10 verified by tRNAscan-SE 2.0.12 on NCBI genome assemblies.
17 pairings (6 independent). 3-tier architecture: large-genome eukaryotes / small-genome eukaryotes / bacteria.

### tRNAscan-SE verified organisms

| Organism | Code | Assembly | Reassigned AA | Count | Total |
|----------|------|----------|--------------|-------|-------|
| T. thermophila | Table 6 | GCF_000189635.1 | Gln | **54** (15+39 supp) | 718 |
| P. tetraurelia | Table 6 | GCF_000165425.1 | Gln | **18** (7+11 supp) | 216 |
| O. trifallax | Table 6 | GCA_000295675.1 | Gln | **8** (2+6 supp) | 94 |
| P. persalinus | Table 6 | GCA_001447515.1 | Gln | **20** (5+15 supp) | 262 |
| H. grandinella | Table 6 | GCA_006369765.1 | Gln | **9** (6+3 supp) | 130 |
| B. stoltei | Table 15 | GCA_965603825.1 | Trp | **6** (CCA:6) | 169 |
| F. salina | Standard | GCA_022984795.1 | — | Gln=3 | 89 |
| S. coeruleus | Standard | GCA_001970955.1 | — | Gln=11 | 272 |
| I. multifiliis | Standard | GCF_000220395.1 | — | Gln=3 | 150 |
| B. nonstop P57 | Table 31 | GCA_028554745.1 | Trp | **2** (CCA, stem-shortened) | 68 |

## Finding 1: Cross-Metric Coloring Optimality (NEW)

The standard code is significantly error-minimizing under **all four** physicochemical metrics:

| Metric | Quantile | p-value | Effect size (z) |
|--------|----------|---------|----------------|
| Grantham (1974) | 0.2% | **0.003** | 2.48 |
| Miyata (1979) | 0.0% | **0.001** | 3.30 |
| Polar Requirement (Woese 1973) | 0.3% | **0.004** | 2.86 |
| Kyte-Doolittle (1982) | 0.0% | **0.001** | 2.95 |

*Note: Freeland & Hurst (1998) used polar requirement only. Cross-metric robustness establishes this as a general structural property.*

- Robust across all rho in [0,1] (all p < 0.01, strengthens at rho=1)
- Stop penalty sensitivity: immaterial (tested 0, 150, 215, 300)
- Score decomposition: Position 2 = 49.3%, Position 1 = 38.2%, Wobble = 12.5%

## Finding 2: Per-Table Optimality Preservation

- **24/25** tables significant at p < 0.05 (raw AND BH-FDR corrected)
- Mean quantile: 1.4% | Only Table 3 (yeast mito) fails at 7.4%
- BH-FDR correction applied

## Finding 3: Topology Avoidance (SUPPORTED)

- Observed: 6/27 (22.2%) create new disconnections
- Possible: 931/1,280 (72.7%) create new disconnections
- Hypergeometric p = **4.8 x 10^-8** (finite-landscape effect size)
- Table-preserving permutation p = **1.0 x 10^-4** (n=10,000, seed=135325)
- **Phylogenetic sensitivity (NEW)**: 3.3x depletion survives all clade exclusions (all p < 10^-5)
- Lineage-collapsed (Sengupta et al. 2007): 27 unique events across 15+ independent lineages

## Finding 4: tRNA Enrichment

- 17 pairings across 4 eukaryotic supergroups (Alveolata, Opisthokonta, Archaeplastida, Excavata)
- 10 tRNAscan-SE verified organisms from NCBI genome assemblies

| Test | n | Statistic | P-value |
|------|---|-----------|---------|
| Fisher+Stouffer (all pairings) | 17 | Z=4.17 | **1.5 x 10^-5** |
| Fisher+Stouffer (6 independent) | 6 | Z=2.14 | **0.017** |
| **MIS worst-case (NEW)** | 6 | Z=1.73 | **0.042** |
| MIS median | 6 | Z=1.74 | **0.041** |
| MIS best-case | 6 | Z=2.14 | **0.016** |

*MIS = Maximal Independent Set enumeration (Bron-Kerbosch). 4/4 MIS significant at p<0.05, eliminating greedy selection bias.*

### Boundary conditions (informative negatives)

| Organism | Domain | Mechanism | tRNA duplication? |
|----------|--------|-----------|-------------------|
| Ciliates (large nuclear genomes) | Eukaryota | Gene duplication | **YES** — 39 suppressor tRNAs in Tetrahymena |
| B. nonstop P57 (24.7 Mb, GCA_028554745.1) | Eukaryota | Anticodon stem shortening 5bp→4bp + eRF1 Ser74Gly (Kachale 2023) | **NO** — 2 tRNA-Trp(CCA) read UGA+UGG |
| Mycoplasma (0.8 Mb genome) | Bacteria | Anticodon modification | **NO** — 1 Trp tRNA reads both UGA+UGG |

**Interpretation**: tRNA gene duplication accompanies codon reassignment in organisms with large nuclear genomes (ciliates, yeasts). Alternative mechanisms in other lineages: anticodon stem shortening of tRNA-Trp(CCA) from 5bp to 4bp (B. nonstop, Kachale 2023 Nature 613:751-758), anticodon modification (Mycoplasma). Multiple evolutionary routes to the same outcome.

### Caveats

1. **Non-independence addressed**: MIS enumeration replaces greedy selection. Worst-case p=0.042 across all 4 maximal independent sets.
2. **Estimated counts**: Euplotes tRNA counts are literature estimates, not tRNAscan-SE verified. B. japonicum replaced by B. stoltei (verified).
3. **Stouffer on discrete p-values**: AA-label exact test has p_min=0.05 (1/20 AAs) — inflated precision. MIS Fisher-Stouffer is the primary conservative test.

## Clean Negatives

| Analysis | Result | Interpretation |
|----------|--------|---------------|
| Local mismatch cost | p=0.70 | Reassignment NOT driven by escaping bad neighborhoods |
| KRAS-Fano | p=1.0 | XOR structure has no clinical prediction power |
| Depth calibration | rho=0.0 | No monotonic epsilon-age relationship |
| Bit-position bias (deduplicated) | p=0.075 | Marginal after correcting for non-independence |
| CUB correlation | p=0.99 | Codon usage not correlated with local mismatch cost |
