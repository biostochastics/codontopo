<div align="center">

  <img src="codon_topo.png" alt="CODON-TOPO logo" width="220" />

  # CODON-TOPO

  **Reproducibility pipeline for the *GF(2)⁶* genetic-code analysis**

  [![Version](https://img.shields.io/badge/version-0.6.1-blue)]()
  [![Status](https://img.shields.io/badge/paper-in%20press%2C%20BioSystems-blue)]()
  [![Tests](https://img.shields.io/badge/tests-440%20passing-success)]()
  [![Coverage](https://img.shields.io/badge/coverage-%E2%89%A596%25-brightgreen)]()
  [![Python](https://img.shields.io/badge/python-3.11%2B-yellow)]()
  [![License: CC BY-NC 4.0](https://img.shields.io/badge/license-CC%20BY--NC%204.0-lightgrey)](LICENSE)

</div>

---

## What this is

This repository accompanies Clayworth & Kornilov (2026), *Robust error-minimization in the genetic code across physicochemical metrics and variant codes: a graph-theoretic analysis in GF(2)⁶* (**BioSystems, in press**; article ref. BIO_BIOSYS-D-26-00689). It provides the Python package, data, and typesetting sources needed to regenerate every number, figure, and table in the paper from a single random seed. The current release accompanying the accepted manuscript is tag `v0.6.1`.

**Preprint:** [bioRxiv 2026.04.25.720843v1](https://www.biorxiv.org/content/10.64898/2026.04.25.720843v1)

## What the paper claims, in short

Encode each of the 64 codons as a 6-bit binary vector in *GF(2)⁶* (two bits per nucleotide position). Two independent questions can then be asked with the same algebraic scaffolding:

1. **Does the assignment of codons to amino acids minimize the physicochemical cost of point mutations?** Under a quartet-pattern shuffle null, the standard code lies in the extreme low-cost tail of 10,000 permutations across four independent distance metrics (Grantham *p* = 0.0062, Miyata *p* < 0.001, Woese polar requirement *p* = 0.003, Kyte–Doolittle *p* = 0.001), and the effect strengthens monotonically as the mutation graph is broadened from the hypercube *Q₆* toward the full single-nucleotide graph *H(3,4)*. Across the 27 NCBI translation tables, 26 preserve near-optimality after BH–FDR correction.
2. **Do natural codon reassignments avoid moves that disconnect an amino acid's codon family?** Under encoding-independent *H(3,4)* adjacency, only 6 of 28 observed reassignments are topology-breaking versus 66% of the 1,280 candidate moves (RR ≈ 0.32, permutation *p* ≤ 10⁻⁴). The *H(3,4)* result is stable by construction; the *Q₆* decomposition is representation-specific and fails to show depletion under 8 of 24 base-to-bit encodings, so *H(3,4)* is the primary test and *Q₆* is reported as a sensitivity. Under the primary *H(3,4)*/Δβ₀ > 0 definition, clade-exclusion sensitivity gives all seven regimes *p* < 10⁻³ (largest *p* ≈ 2 × 10⁻⁴); under the *Q₆*/new-disconnection sensitivity, all seven give *p* < 10⁻⁵. Event-level conditional-logit models show topology and local physicochemical cost are weakly correlated and jointly informative (ΔAICc(M1→M3) ≈ 110 under *Q₆*, 93 under *H(3,4)*).

The core contribution is the second axis: code evolution appears to be constrained by both physicochemical smoothness and codon-family topological integrity, and the two constraints are partly independent.

## Results in a little more detail

**Coloring optimality.** All five distance measures place the standard code in the extreme low-cost tail. The four code-independent measures give *z*-scores from 2.35 (Grantham) to 3.14 (Miyata), improvements of 9.9–19.9% over the null mean, and pass Bonferroni correction at α/5 = 0.01. Under the alignment-derived ProtSub matrix (a code-dependent sensitivity check), the standard code sits at the 0.04% quantile (*p* = 4 × 10⁻⁴), the most extreme percentile of any measure tested. Partial Spearman correlations across the five measures show that four of ten pairs are near-zero after adjustment (|ρ_partial| < 0.1; a descriptive point-estimate statement), so the convergence is not driven by a single latent factor. Decomposing the total mismatch score by codon position (Figure 5B) shows position 2 contributing 49.3%, position 1 contributing 38.2%, and the wobble position 12.5%, mirroring the biochemical hierarchy of mutational impact.

**ρ-robustness sweep.** Interpolating between *Q₆* (192 Hamming-1 edges) and *H(3,4)* (all 288 single-nucleotide edges) with a weight ρ ∈ {0, 0.25, 0.5, 0.75, 1}, every ρ yields *p* < 0.01 and *z* rises monotonically from 2.48 to 3.46. The optimality signal is not an artifact of ignoring within-nucleotide transversion edges.

**Per-table preservation across the 27 NCBI tables.** 26 of 27 tables sit in the top 5% of their own quartet-pattern shuffle null after BH–FDR correction; mean quantile 1.4%. Yeast mitochondrial code (table 3) is the sole marginal exception (*p*_BH = 0.075, six codon reassignments — the most extensively reassigned code). A standard-code-proximity audit (supplement §S8) makes explicit that every variant sits at the extreme low-*d_H* tail of its own null, so the "informative-distance" vs "near-standard" disaggregation rests on the absolute *d_H* ≥ 3 threshold rather than on a conditional-null comparison.

**Topology avoidance under *H(3,4)*.** Of the 28 de-duplicated observed reassignment events across 25 variant-code tables, 6 (21.4%) are topology-breaking versus 66.1% of the 1,280 candidate single-codon relabelings. Hypergeometric one-sided *p* = 1.28 × 10⁻⁶, table-preserving permutation *p* ≤ 10⁻⁴, risk ratio 0.32 (95% CI 0.16–0.66), depletion fold 3.1×. All four cells of the 2×2 adjacency × topology-breaking-definition audit (*Q₆* vs *H(3,4)* × new-disconnection vs Δβ₀>0) yield hypergeometric *p* < 10⁻⁵ with risk ratios 0.28–0.33.

**Conditional-logit decomposition.** Six models fit on 66 event-steps across 25 variant-code tables. The combined physicochemistry + topology model M3 is strongly favored: ΔAICc(M1→M3) = 110.4 under *Q₆*, 93.4 under encoding-independent *H(3,4)*; ΔAICc(M2→M3) = 91.2 (*Q₆*), 97.3 (*H(3,4)*). Adding a heuristic tRNA-distance proxy (M4) does not improve fit (LR = 0.12, *p* = 0.73). Spearman ρ between local physicochemical cost Δ_phys and topology disruption Δ_topo across the candidate landscape is 0.15 — the two features are almost independent predictors.

**tRNA-gene enrichment (exploratory).** 18 organism assemblies were scanned in this work with tRNAscan-SE 2.0.12. Fifteen of these populate the 24-pairing enrichment analysis as variant-code or standard-code-control repertoires; the remaining three (*Blastocrithidia nonstop*, *Mycoplasmoides genitalium*, *M. pneumoniae*) are mechanistic boundary cases outside the 24-pairing set. The full 24 pairings span 24 unique organisms and 48 slots; nine slots use literature/GtRNAdb/annotation counts, with primary sources including Bonitz *et al.* 1980 (*S. cerevisiae* mitochondrial, 24 tRNAs) and Kerscher *et al.* 2001 (*Y. lipolytica* mitochondrial, 27 tRNAs). Variant codes carry more tRNA genes for the reassigned amino acid than expected under an all-pairings Fisher–Stouffer combined test (*p* = 1.7 × 10⁻⁷; not adjusted for shared controls). Enumerating all 332 maximal-independent-pairing subsets of size 6 (via Bron–Kerbosch on the complement of the conflict graph) gives median Stouffer *p* = 0.037, best 0.012, worst 0.104; 264 of 332 (79.5%) sets fall below the 0.05 threshold. A pre-specified topology-breaking subset (*n* = 4) is null (Stouffer *p* = 0.387), so the pattern is best read as heterogeneous accommodation-of-decoding across variant-code lineages driven by UAR→Gln ciliates, UGA→Cys *Euplotes*, and *Blepharisma stoltei* UGA→Trp, rather than as compensation for codon-family disconnection. The mechanism landscape is not monolithic — *Blastocrithidia nonstop* achieves UGA→Trp via anticodon-stem shortening rather than gene duplication.

**Cross-study reanalysis of genome-recoding datasets.** Retrospective analysis of nine published genome-recoding experiments (>217,000 codon-level observations) shows results consistent with codon-family topology operating at a different biological layer from acute cellular fitness: Syn61 tolerated 18,218 boundary-crossing Ser→Ser swaps genome-wide, yet the same move type is 3.1-fold depleted in natural code evolution. The reanalysis does not establish causality but rules out "topology avoidance = engineering-scale lethality."

**Walsh–Hadamard / 2-adic spectral probe.** As a methodological bridge to the 2-adic codon algebra of Khrennikov & Kozyrev (2007) and Dragovich & Mišić (2019), we compute the Walsh spectral depth (block-indicator geometry) and a label-spectrum invariant *S* (wobble-free vs wobble-active partition of the 63 non-DC Walsh frequencies). Under a block-size matched null (*n* = 2,000), depth *z* = −17.74 and is mathematically constant across all 24 base-to-bit bijections; a stricter wobble-box-preserving null reveals that the depth is an algebraic invariant of the (wobble box × AA slot multiset) structure. The label-spectrum invariant *S* = 0.7514 is exact and encoding-independent.

**Falsifications and rejections.** The KRAS–Fano clinical prediction — that XOR triples in *GF(2)⁶* would predict co-mutation partners at KRAS G12 hotspots — was falsified on *n* = 1,670 MSK-IMPACT samples (*p* = 1.0 for all six G12 variants after Bonferroni correction). Three algebraic conjectures were pre-rejected on their own merits: a Serine minimum-Hamming-distance-4 invariant (encoding-dependent), a PSL(2,7) symmetry (literature-refuted), and a holomorphic embedding (finite domain).

The full status of all 15 claims — supported, exploratory, rejected, falsified, tautological — is registered in `src/codon_topo/reports/claim_hierarchy.py` and viewable with `codon-topo claims`.

## Reproduce the paper

One canonical command reproduces every displayed value in the manuscript and supplement from a clean clone at tag `v0.6.1`:

```bash
git clone https://github.com/biostochastics/codontopo && cd codontopo
git checkout v0.6.1
pip install -e ".[dev]"
bash scripts/build_publisher_release.sh
```

The wrapper `scripts/build_publisher_release.sh` runs the analysis pipeline (`codon-topo all` at seed 135325, *n* = 10 000), the provenance-emit and table-generation scripts, both R figure scripts, and the two Typst PDF compilations, in that order. The same block appears verbatim in the Elsevier response letter, `CHANGELOG.md`, main-text §5, and supplement §S24.

The single JSON file `output/manuscript_stats.json` holds every statistic cited in the paper; the Typst source reads from it, so within a single pipeline run the manuscript and its supplement render from shared versioned artifacts.

**Random seed:** `135325` for all Monte Carlo analyses.

## Command-line entry points

```bash
codon-topo all               # run everything and write manuscript_stats.json
codon-topo claims            # view the 15-claim hierarchy with p-values
codon-topo coloring          # coloring-optimality Monte Carlo (n = 10,000)
codon-topo per-table         # all 27 NCBI translation tables + BH-FDR
codon-topo topology-avoidance-k43  # encoding-independent H(3,4) topology test
codon-topo condlogit         # conditional-logit models M1-M4
codon-topo trna              # tRNA-gene enrichment across 18 verified genomes
codon-topo kras              # KRAS-Fano falsification test
codon-topo codonsafe         # cross-study reanalysis of 8 recoding datasets
```

`codon-topo --help` lists the 21 subcommands. All support `--json` for machine-readable output.

## Tests

```bash
python3.11 -m pytest tests/ -q                    # 454 tests, ~90 s
python3.11 -m pytest tests/test_regression.py -v  # the 105 numeric-regression tests
```

Use `python3.11 -m pytest` explicitly if your system default Python differs from where dev dependencies are installed.

## Layout

```
src/codon_topo/
  core/            # GF(2)^6 encoding, 27 NCBI tables, filtration, homology, Fano
  analysis/        # coloring optimality, topology avoidance, conditional logit, tRNA, KRAS, ...
  reports/         # claim_hierarchy.py (single source of truth for what the paper claims)
  visualization/   # CSV export + R ggplot2 figure scripts
tests/             # 454 unit + regression tests
scripts/           # analysis one-offs (stats generation, tRNAscan wrapper, cross-study reanalysis)
output/            # manuscript.typ, supplement.typ, manuscript_stats.json, figures/, PDFs
data/              # tRNAscan-SE results, genome assembly accessions, CodonSafe raw data
```

See [`CLAUDE.md`](CLAUDE.md) for contributor notes.

## Optional dependencies

- **R 4.5+** with `ggplot2`, `ggpubr`, `viridis`, `patchwork` — publication figures.
- **tRNAscan-SE 2.0.12** with Infernal 1.1.4 — regenerate the 18-genome tRNA scans (raw results are shipped in `data/trnascan_results/`).
- **Typst** and a TeX Live install — regenerate the manuscript PDFs. The LaTeX pipeline additionally needs `pandoc` 3.1+.
- **Tsour & Slavov (2026) SAAP table** — download `Supplementary_Data_8.High_confidence_SAAP_precursor_quant.xlsx` from [decode.slavovlab.net](https://decode.slavovlab.net/mass-spec/data) and pass it to `scripts/slavov_saap_analysis.py` to reproduce the cross-study cross-classification in §4.6 / §S17.

## License

Released under [CC BY-NC 4.0](LICENSE). Share and adapt with attribution; commercial use requires a separate license — contact the authors.

To cite, see [`CITATION.cff`](CITATION.cff) or use GitHub's "Cite this repository" button.
