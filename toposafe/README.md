# TopoSafe

**Cross-study retrospective meta-analysis of genome recoding experiments**

*Companion analysis to Clayworth & Kornilov (2026). Journal of Theoretical Biology.*

---

## What is TopoSafe?

TopoSafe tests whether codons classified as **topology-breaking** in the GF(2)^6 genetic code graph predict higher rates of recoding failure or fitness reduction, across 9 published genome recoding datasets (>248,000 codon positions).

### Datasets

| Study | Species | Codons | Topology variation |
|-------|---------|--------|-------------------|
| Napolitano 2016 (PNAS) | *E. coli* | 12,888 AGR→CGN | None (Arg connected) |
| Fredens 2019 (Nature) | *E. coli* | 18,218 Ser | All boundary-crossing |
| Ostrov 2016 (Science) | *E. coli* | 62,214 | Mixed (Ser/Leu/Arg) |
| Robertson 2025 (Science) | *E. coli* | 60,240 | Ser=100%, Ala=0% |
| Grome/Robertson 2025 (Nature) | *E. coli* | 1,195 | Stop codon only |
| Frumkin 2018 (PNAS) | *E. coli* | 60 | None (Arg connected) |
| Lajoie 2013 (Science) | *E. coli* | 321 | Stop codon only |
| Ding 2024 (Science) | Mammalian | Variable | Ser boundary (cross-organism) |
| Nyerges 2024 (bioRxiv) | *E. coli* | 62,007 | Multi-omics validation |

### Key Results

| Claim | Outcome | p-value / n |
|-------|---------|-------------|
| Local mismatch asymmetry (Syn57 Ser vs Ala) | **Positive** | p ≈ 0 (Mann-Whitney, n=60,240) |
| Boundary crossing predicts RNA-seq perturbation | **Null** | p = 0.40 |
| Local mismatch ~ SRZ correlation (Napolitano) | **Positive** | rho = -0.33, p ≈ 0 |
| Cross-organism Ser boundary confirmation | **Positive** | Ding 2024, mammalian cells |
| Topology predicts fix events in Syn57/Syn61 | **Null** | Fix events = regulatory/structural |

Full results in [`data/codonsafe/ANALYSIS_RESULTS.md`](../data/codonsafe/ANALYSIS_RESULTS.md).

---

## Quick Start

### 1. Install

```bash
# From the codon-topo repo root
pip install -e ".[codonsafe]"
```

### 2. Download raw data

```bash
bash toposafe/scripts/download_data.sh
```

> **Note**: Several datasets require institutional journal access (Nature, Science).
> See `data/codonsafe/DATA_MANIFEST.md` for full provenance and manual download instructions.

### 3. Run the analysis

```bash
bash toposafe/scripts/run.sh
```

Or step by step:

```bash
# Python pipeline (classification + statistics)
codon-topo codonsafe --output-dir=output/codonsafe

# R figures
Rscript src/codon_topo/visualization/R/codonsafe_figures.R
```

### 4. Outputs

| Path | Contents |
|------|----------|
| `output/codonsafe/analysis_stats.json` | All statistics (machine-readable) |
| `output/codonsafe/syn57_s9_*.csv` | Syn57 deviation classification tables |
| `output/codonsafe/ostrov_*.csv` | Ostrov segment analysis tables |
| `output/codonsafe/fig_cs1_*.pdf/png` | Figure CS1: Syn57 deviation directionality |
| `output/codonsafe/fig_cs2_*.pdf/png` | Figure CS2: Ostrov segment topology |
| `output/codonsafe/fig_css1_*.pdf/png` | Figure CS-S1: Sensitivity analysis |
| `output/codonsafe/fig_recoding_*.pdf/png` | Combined recoding reanalysis figure |
| `output/supplement_codonsafe.pdf` | Compiled TopoSafe supplement |

---

## Architecture

TopoSafe is built on the `codon-topo` core package:

```
codon-topo (core theory)
    └── codon_topo.analysis.codonsafe (TopoSafe analysis)
            ├── models.py          — StudyId, AnnotatedSwap, TopologyClassification
            ├── classify.py        — classify_swap_event(), build_reference_context()
            ├── aggregate.py       — cross-study aggregation
            ├── stats.py           — Fisher tests, Mann-Whitney, Spearman
            ├── normalize.py       — codon normalization utilities
            ├── run_analyses.py    — full pipeline driver
            └── loaders/
                    ├── fredens2019.py
                    ├── genbank_utils.py
                    ├── napolitano2016.py
                    └── ostrov2016.py
```

Visualization scripts: `src/codon_topo/visualization/R/codonsafe_figures.R`

---

## License

MIT. See main [`codon-topo`](../README.md) repository.
