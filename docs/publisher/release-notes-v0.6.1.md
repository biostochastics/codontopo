# codon-topo v0.6.1 — 2026-08-15

**Third pre-production release for BIOSYS-D-26-00689** (Clayworth & Kornilov, *Robust error-minimization in the genetic code across physicochemical metrics and variant codes: a graph-theoretic analysis in GF(2)⁶*, *BioSystems*, in press).

This release supersedes v0.6.0. It closes out a third external release-gate review that flagged three computational defects and one data-audit gap. **The four SUPPORTED scientific claims and the EXPLORATORY tRNA-enrichment classification are unchanged. Grantham conservative p under the quartet-pattern shuffle moves from 0.006099 to 0.006199 (still below the Bonferroni threshold α = 0.05/8 = 0.00625); MIS median Stouffer p tightens from 0.046 to ~0.037; the worst-case-MIS gate remains failed.**

## What changed

### Null-model preservation table corrected (Supplement §S2.1)
The table described the wrong preserved objects. The corrected table:
- The **quartet-pattern shuffle** does NOT fix the codon-family partition (the atomic 4-codon patterns move across quartet slots, so which specific codons decode as a named amino acid changes from draw to draw). It DOES fix each named amino acid's total codon count, because each mobile pattern carries its internal AA composition with it as a labeled unit. Stop-containing UA and UG quartets are held fixed and the remaining 14 patterns are permuted uniformly. State space corrected from ~$16!$ to $\approx 14! \approx 8.7 \times 10^{10}$.
- The **classical Haig–Hurst amino-acid permutation** fixes the codon-family partition (which codons decode together) and permutes the 20 AA labels uniformly. Per-named-AA codon count is NOT preserved (Trp's single-codon family becomes a 6-codon family if its label swaps with Leu's). State space $\approx 20!$.

Main §2.3.1 rewritten to match.

### Monte Carlo tail convention fixed (< → ≤)
`coloring_optimality.monte_carlo_null` (and `rho_robustness_sweep`, `encoding_sensitivity_of_optimality`, `per_table_proximity_audit`) counted null draws with `f < F_obs` rather than the standard lower-tail convention `f <= F_obs`, undercounting draws that tie the observed score. The Grantham null contains one draw exactly at $F = 13\,477$; conservative $p_{\text{cons}} = 61/10\,001 = 0.006099$ becomes $62/10\,001 = 0.006199$. All eight quartet-pattern-shuffle + Haig–Hurst *p*-values continue to pass Bonferroni at $\alpha = 0.05/8 = 0.00625$ (the previous claim "$\alpha = 0.05/5 = 0.01$" for eight tests was incoherent and is corrected). Diagnostic fields `n_at_or_below_observed`, `n_strictly_below_observed`, and `n_ties_at_observed` added to the JSON output.

### tRNA legacy repertoires audited against primary literature
- *Saccharomyces cerevisiae* mitochondrion (Bonitz et al. 1980 PNAS 77:3167–3170, PMC349575; RefSeq NC_001224.1): total is 24 tRNA genes. The v0.6.0 vector summed to 23 because Arg was set to 1; the standard S288C annotation has two Arg tRNAs (tRNA-Arg1 for AGR + tRNA-Arg2 for CGN). Fixed to **Arg = 2**.
- *Yarrowia lipolytica* mitochondrion (Kerscher et al. 2001 Comp Funct Genomics 2:80–90, PMC2447202, DOI 10.1002/cfg.72; EMBL AJ307410): the v0.6.0 vector was a **placeholder** (1 per amino acid, sum 20, table_id 3). Primary literature reports **27 functional tRNAs** under the **mould mitochondrial code (NCBI table 4, not table 3)**. Rewritten with the Kerscher distribution: Ile = 2 (I1, I2), Leu = 3 (L1, L2, L3), Lys = 2 (K1, K2), Met = 2 (initiator, elongator), Ser = 2 (S1, S2 — the ySer3 pseudogene in COB-I3 is not counted in the 27-tRNA total), Tyr = 2 (Y1, Y2), Arg = 1 (AGR reader only; no CGN-Arg tRNA is a hallmark of this genome), all others = 1. Table id 3 → 4.

Downstream effects:
- Fisher cell for the yeast-mitochondrial Thr pairing: (2/24) vs (1/27), one-sided $p = 0.455$ (was 0.554).
- MIS distribution: median Stouffer $p \approx 0.037$ (was 0.046), worst $\approx 0.104$ (was 0.123), topology-breaking subset $\approx 0.387$ (was 0.43); fraction of MIS below 0.05 rises from 57.2% to $\approx 79.5\%$.
- Pre-specified worst-case-MIS decision gate remains failed; tRNA-enrichment claim remains EXPLORATORY.

### H(3,4) clade-exclusion table added
The v0.6.0 prose claimed the topology-avoidance depletion was significant in every clade-exclusion regime, but only the $Q_6$ / new-disconnection seven-row table was persisted; the primary-cell $H(3,4)$ / $\Delta \beta_0 > 0$ analogue was not tabulated. New `topology_avoidance_k43_phylogenetic_sensitivity()` function emits `output/phylogenetic_sensitivity_k43.json` as part of the standard pipeline, and Supplement §S9 now displays both tables. All seven exclusions significant at $p < 10^{-3}$; excluding yeast mitochondrial gives $p \approx 4.2 \times 10^{-9}$, matching the direct calculation quoted in main §3.4. The stale "not tabulated in v0.5.2" phrase is removed.

### Cross-references and labels
- §2.1: Supplement §S3 → §S2 (encoding sweep target).
- §2.3.1: Supplement §S3.1 → §S2.1 (dual-null comparison target).
- §3.4: Supplement §S8 → §S9 (topology clade-exclusion target); the H(3,4) 6/28 rate is no longer mislabeled "new-disconnection" (it uses $\Delta \beta_0 > 0$).
- §3.5: Supplement §S13 → §S19.5, Figure S8 panel C (predictive-diagnostic target).
- §2.3.5: "four nested models" reframed as "four candidate models; likelihood-ratio tests only on nested pairs" (M1 and M2 are non-nested alternatives).
- Main Table 1: "FH per-table" and "FH weighted edges" replaced with explicit quartet-pattern-shuffle names.
- Tables 2 and 3 captions: "block-preserving null" → "quartet-pattern shuffle null".
- Figure S1 subtitle and Figure S5 x-axis label / title: "Freeland-Hurst block-preserving" → "quartet-pattern shuffle"; "transversion weight $\rho$" → "diagonal-edge weight $\rho$" (matches the main-text §2.1 caveat that $\rho$ is not a transition/transversion weight).
- Figure 3B caption clarified: 6 calibration points spanning 4 amino acids (universal Serine baseline + Thr, Leu, and Ala variant-code lineages; Serine and Leucine each contribute two points from different lineages).
- Supplement Table S6 caption "Full" → "Unrestricted" (matching the row label).

### Reproducibility
- Version bumped 0.6.0 → 0.6.1 across `pyproject.toml`, `CITATION.cff`, `llms.txt`, and the manuscript/supplement reproducibility stamps. `manuscript_stats._version` bumped to `0.6.1`.
- Full pipeline rerun at seed 135325 (same seed as v0.6.0) after all code fixes; all JSON artefacts and R-rendered figures regenerated.

## Scientific status (unchanged)

| Claim | Status | p-value |
|:--|:--|--:|
| Cross-metric coloring optimality (four metrics) | Supported | ≤ 0.0062 (Grantham QP after tie fix; all four metrics pass Bonferroni α = 0.00625) |
| Per-table preservation across 27 NCBI tables | Supported | BH–FDR; 11/12 informative + 10/11 distinct colorings retain top-5% |
| ρ-robustness | Supported | 0.0062 (Q₆) → 0.0003 (H(3,4)) |
| Topology avoidance (Q₆ + H(3,4)) | Supported | ≤ 10⁻⁴ (permutation) |
| tRNA enrichment (MIS median) | Exploratory | ≈ 0.037 (median); worst 0.104; topology-breaking subset 0.387 |
| KRAS–Fano clinical prediction | Falsified | 1.0 |

## Reproduction

```
git checkout v0.6.1
pip install -e ".[dev]"
codon-topo all --output-dir=./output --seed=135325 --n=10000
python3.11 scripts/patch_pt_keys.py
python3.11 scripts/emit_trna_provenance_table.py
Rscript src/codon_topo/visualization/R/all_figures.R output output/figures
Rscript src/codon_topo/visualization/R/strengthened_figures.R output output/figures
typst compile output/manuscript.typ output/manuscript.pdf
typst compile output/supplement.typ output/supplement.pdf
```

## Contact

Sergey Kornilov <sergey@biostochastics.com>
Biostochastics, LLC — Seattle, WA, USA
