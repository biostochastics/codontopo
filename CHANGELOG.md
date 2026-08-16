# CHANGELOG — CODON-TOPO

Release notes for the `codon-topo` pipeline and the accompanying
manuscript **BIOSYS-D-26-00689** (Clayworth & Kornilov, *BioSystems*, in
press). Each release corresponds to a manuscript revision or a
proofs-corrections pass.

---

## [v0.6.1] — 2026-08-15 — Pre-publication QA/QC of Figures, Tables, and Statistics

Consolidated close-out of a pre-publication QA/QC audit of every figure, table, and inferential quantity across the manuscript and the supplement. **All scientific conclusions of the accepted manuscript are unchanged**: four SUPPORTED claims (cross-metric coloring optimality, per-table preservation, ρ-robustness, topology-avoidance depletion), one EXPLORATORY (tRNA enrichment), one FALSIFIED (KRAS–Fano), and the remaining REJECTED/TAUTOLOGICAL entries as reported in the accepted version. Two substantive computational corrections (the Bron–Kerbosch MIS pivot in the tRNA-enrichment robustness bound and the Monte Carlo tie-inclusion in the coloring-optimality null) shifted several displayed values — the largest shift is the tRNA reclassification from SUGGESTIVE to EXPLORATORY (median MIS Stouffer *p* moves 0.046 → 0.037, worst 0.123 → 0.104, fraction of MIS < 0.05 rises 190/332 → 264/332, topology-breaking subset *p* moves 0.43 → 0.387); the Grantham quartet-pattern shuffle *p* moves 0.006099 → 0.006199 (tie-inclusion). No claim crosses a decision threshold. Every before/after value is enumerated in `docs/publisher/2026-08-14-corrections-ledger.md`.

### Statistics — corrections to inferential quantities

- **Null-model description.** Supplement §S2.1 preservation-property table corrected. The quartet-pattern shuffle *does not* fix the codon-family partition (atomic 4-codon quartet patterns permute across quartet slots, so which specific codons decode as a named amino acid varies from draw to draw); it *does* fix each named amino acid's total codon count. The classical Haig–Hurst amino-acid-permutation null *does* fix the codon-family partition; it *does not* fix per-named-AA codon count. State-space count for the quartet-pattern shuffle corrected from $≈ 16!$ to $≈ 14! ≈ 8.7 × 10^{10}$ — only 14 non-stop quartets are mobile in the standard code (UA and UG contain stops and are held fixed). Main §2.3.1 rewritten to match; sensitivity comparison against the classical HH-AA null (four physicochemical metrics × two nulls; all eight *p*-values pass Bonferroni at $α = 0.05/8 = 0.00625$) reported in Supplement §S2.1.

- **Monte Carlo tail convention.** The coloring-optimality Monte Carlo now counts null draws satisfying $f ≤ F_\text{obs}$ (standard lower-tail permutation-test convention, ties included) rather than strict $f < F_\text{obs}$. The Grantham null contains one draw exactly equal to $F_\text{obs} = 13\,477$; conservative $p_\text{cons}$ moves from $61/10\,001 = 0.006099$ to $62/10\,001 = 0.006199$. The four other code-independent metrics and ProtSub have no ties and are numerically unchanged. All five metrics continue to pass Bonferroni at $α = 0.05/5 = 0.01$; all eight quartet-pattern-shuffle + Haig–Hurst *p*-values pass at $α = 0.05/8 = 0.00625$.

- **tRNA reference vectors.** Two entries in the curated tRNA repertoire were reconciled with primary literature: *Saccharomyces cerevisiae* mitochondrion — Arg count corrected from 1 to 2 (tRNA-Arg1 for AGR + tRNA-Arg2 for CGN), reaching the Bonitz et al. 1980 (PMC349575) total of 24 tRNA genes; *Yarrowia lipolytica* mitochondrion — vector rewritten with the Kerscher et al. 2001 (PMC2447202, DOI 10.1002/cfg.72) 27-tRNA distribution under NCBI table 4 (mould mitochondrial code, not yeast mitochondrial) — Ile = 2, Leu = 3, Lys = 2, Met = 2, Ser = 2, Tyr = 2, others = 1; no CGN-Arg tRNA (a hallmark of this genome). The Fisher cell for the yeast-mitochondrial Thr pairing becomes (2/24) vs (1/27), one-sided *p* = 0.455 (was 0.554). Downstream MIS distribution: median Stouffer *p* ≈ 0.037 (was 0.046), worst *p* ≈ 0.104 (was 0.123), topology-breaking subset Stouffer *p* ≈ 0.387 (was 0.43); fraction of 332 MIS below 0.05 rises from 57.2 % to ≈ 79.5 %. The pre-specified worst-case-MIS decision gate remains failed; the tRNA-enrichment claim remains EXPLORATORY.

- **Maximal-independent-set (MIS) enumeration.** The MIS enumeration used for the tRNA-enrichment robustness bound (Supplement §S10.2) was recomputed with a corrected Bron–Kerbosch traversal on the complement of the conflict graph, yielding 332 MIS of size 6 rather than 2. The full Stouffer-*p* distribution across all 332 MIS (median, best, worst, and topology-breaking-subset) is now reported. The all-pairings Fisher–Stouffer result (*p* ≈ 1.7 × 10⁻⁷) is unchanged; the tRNA-enrichment claim, previously SUGGESTIVE, is reclassified as EXPLORATORY on the basis that the median-independent-subset signal is present but the strict worst-case bound is not met and the topology-breaking-only subset is null. Claim hierarchy remains 15 entries: 4 SUPPORTED / 5 EXPLORATORY / 1 FALSIFIED / 3 REJECTED / 2 TAUTOLOGICAL.

- **H(3,4) clade-exclusion sensitivity table.** The manuscript previously reported clade-exclusion robustness "in every regime" but only tabulated the $Q_6$ / new-disconnection seven-row table. The primary-cell $H(3,4)$ / $Δβ_0 > 0$ seven-row table is now generated and displayed alongside it in Supplement §S9. All seven exclusions are significant at *p* < 10⁻³; excluding yeast mitochondrial (which contributes 4 of the 6 lineage-collapsed *H(3,4)*-breakers) gives *p* ≈ 4.2 × 10⁻⁹, matching the direct calculation quoted in main §3.4.

- **Conditional-logit clade-exclusion — data plumbing.** The main-text §3.5 clade-exclusion summary now reads its ΔAICc aggregate keys (min, median, max across the seven regimes) from a lifted `manuscript_stats.condlogit.clade_exclusion` block rather than silently defaulting to 0 (which previously rendered "min = median = max = 0"). Fallback branches that could mask future data-lift gaps replaced with build-time asserts.

- **Cross-family multiple-comparison correction.** Supplement §S18 corrected: with eight analysis families, cross-family Bonferroni threshold is $α = 0.05/8 = 6.25 × 10⁻³$. The MIS median tRNA-enrichment $p = 0.037$ does **not** survive this threshold (nor does the worst-case $p = 0.104$), consistent with the EXPLORATORY classification. The four SUPPORTED claims all survive.

### Figures, tables, and prose — presentation and label corrections

- **Main Table 1 (claim hierarchy).** tRNA row moved from SUGGESTIVE to EXPLORATORY with updated *p*. "FH per-table" and "FH weighted edges" replaced with explicit "Per-table quartet-pattern shuffle + BH–FDR" and "Quartet-pattern shuffle (weighted edges)". Shaded-row caption reworded so FALSIFIED / REJECTED / TAUTOLOGICAL are described as non-inferential rather than "failed validation".

- **Main Tables 2 and 3.** Captions: "block-preserving null" → "quartet-pattern shuffle null" throughout.

- **Main Table 4.** *H(3,4)* topology-avoidance table now uses the primary-cell values (H(3,4), $Δβ_0 > 0$; N = 1,280, K = 846); the $2 × 2$ definition × adjacency audit reported inline with risk ratios 0.28–0.33 across all four cells.

- **Main Figure 3 (evolutionary evidence).** (A) bit-position bias distribution unchanged. (B) depth calibration — reproducible-seed jitter so coincident CUG-clade markers are visible; caption clarified: 6 calibration points spanning 4 amino acids (universal Serine baseline + Thr / Leu / Ala variant-code lineages; Serine and Leucine each contribute two points from different lineages). (C) topology avoidance shown for both $H(3,4)$ (primary, $Δβ_0 > 0$) and $Q_6$ (sensitivity, new-disconnection). (D) tRNA-rank panel now labels each bar with the control organism so all 24 pairings render distinctly.

- **Main Figure 4 (conditional logit).** Panel A shows six candidate conditional-logit models (M1–M4 under $Q_6$ + two $H(3,4)$-verification variants); likelihood-ratio tests applied only to nested pairs. AICc gaps itemised in the caption (M1 / M2 baselines 89.1–112.1 units above M3; $H(3,4)$ variant 16.9 units; M4 2.1 units). The main-text §2.3.5 description reframed as "four candidate models; likelihood-ratio tests only on nested pairs" (M1 and M2 are non-nested alternatives, both nested within M3; M3 nested within M4).

- **Main Figure 5 catalogue panel.** Named-vector label mapping restored; KRAS bar labelled "Falsified" (was rendering "Pending").

- **Supplement Figure S1 subtitle.** "Freeland–Hurst block-preserving random colorings" → "quartet-pattern shuffle random colorings".

- **Supplement Figure S5 title and *x*-axis.** "Transversion weight $ρ$" → "Diagonal-edge weight $ρ$" (matches the main-text §2.1 caveat that $ρ$ is not a transition/transversion weight).

- **Supplement Figure S7 label.** "disconnection × organism" → "variant-code × control".

- **Supplement Table S6 caption.** "Full" row → "Unrestricted" (matches the row label).

- **Supplement Table S8 caption.** Discloses the tRNAscan-SE two-pass distinction (Infernal-filtered `Supp` column versus first-pass Isotype/Anticodon totals in the Reassigned-AA parenthetical).

- **Supplement Table S10.1 — 24-row tRNA analysis-input provenance table.** New: for each of the 24 pairings, records variant + control organism, compartment, assembly accession, reassignment string, $Q_6$ and $H(3,4)$ topology status, exact 2×2 Fisher counts, per-organism source class (tRNAscan-SE 2.0.12 / GtRNAdb / literature / annotation), and per-pair one-sided Fisher-exact *p*-value.

- **Supplement §S12 — structural-preservation index.** Section renamed *Structural-preservation index (visualization-only)* and rewritten to describe the discrete four-indicator implementation actually used by Figure 5A ($S ∈ \{0.55, 0.75, 0.80, 1.00\}$). The "synthetic-biology feasibility" phrasing was dropped because no independently validated feasibility endpoint exists.

- **Supplement Figure FigB — coloring-optimality per-table panel.** Table 32 (Balanophoraceae plastid) and Table 28 (Condylostoma nuclear) rows restored; CSV exporter no longer drops these two tables.

- **Supplement §S21 (Walsh–Hadamard / 2-adic spectral panel) and §S22 (Slavov / Tsour SAAP grouped-bar figure).** New; provide descriptive bridges to the 2-adic-genetic-code and site-specific-amino-acid-polymorphism literatures respectively. Walsh conventions clarified (spectral depth definition; encoding invariance; block-size-null vs wobble-box-preserving-null distinction).

- **Main §2.6 tRNA decision gate.** Restored to the originally specified worst-case-MIS criterion (not met at *p* = 0.104); median MIS *p* reported descriptively as a summary statistic, not as a re-specified criterion.

- **Main §4.3 title.** *Exploratory tRNA accommodation across variant-code systems*; compensation-for-disconnection claims removed; the genome-size / mechanistic-repertoire correlation explicitly framed as a hypothesis motivated by three case examples (*Blastocrithidia nonstop* anticodon-stem shortening; *Mycoplasmoides* anticodon modification; ciliate tRNA duplication), not as a tested association.

- **Main §3.3 per-table optimality.** Reframed as "significantly low-cost *relative to its own quartet-shuffle null*"; every variant has $d_H ≤ 6$ from the standard code while every null draw has $d_H ≥ 30$, so the per-table quantile cannot identify a boundary between proximity-driven and independently-optimised codes. Both registry-entry (11 of 12 informative-distance variants after BH–FDR) and distinct-coloring (10 of 11) counts reported.

- **Main §2.3.5 conditional logit — restricted-candidate sensitivity.** No single *d* value designated the biological cut; three values ($d ∈ \{1, 2, 3\}$) reported as a bracket with $d = 1$ as the most defensible biological floor.

- **Main §2.3.5 conditional logit — IIA discussion.** IIA acknowledged as an untested structural assumption of the estimator even for the explanatory reading; the restricted-candidate refits framed as sensitivity analyses against candidate-set composition rather than as evidence that IIA is harmless.

- **Main §3.2 $ρ$-sweep.** Five pre-specified inferential grid points ($ρ ∈ \{0, 0.25, 0.5, 0.75, 1\}$) versus eleven plotted points (five inferential + six descriptive interpolation) explicitly distinguished.

- **Main §3.5 parametric-predictive rate.** Reported with explicit numerator / denominator / adjacency / definition ($5/66$ event-steps under $Q_6$ / $Δβ_0 > 0$; parametric predictive $p = 0.60$) and clearly distinguished from the lineage-collapsed $6/28$ rate under $H(3,4)$ / $Δβ_0 > 0$ quoted in main §3.4.

- **Supplement §S3 candidate-universe Table S3.** U3 no-op count corrected to 61 (only sense codons carry an AA identity); source-stop axis removed from the introductory framing.

- **Main §2.3.4 encoding-sweep claim.** Narrowed from a blanket "all encoding-dependent results tested across 24 encodings" to the exact list — coloring optimality + $Q_6$ topology-avoidance landscape depletion. M3-$H(3,4)$ renamed a *topology-verification* variant (its local physicochemical feature is still $Q_6$-based).

- **Supplement §S15 metric-correlation partials.** Descriptive language calibrated: four of ten metric-pair partial Spearman correlations are "approximately uncorrelated after adjustment" ($|ρ_\text{partial}| < 0.1$; a descriptive statement about point estimates, not a test of statistical independence). Main-text §3.1 count reconciled from "five" to "four" to match Table S19.

- **PSL(2,7) rejection.** Narrowed to "no 64-dimensional irreducible linear representation" rather than a general symmetry-exclusion claim.

- **Holomorphic-embedding rejection.** Grounded in the character-identity failure $χ(x + x) = χ(x)^2$ with $i^2 = -1 ≠ 1$, rather than in an appeal to domain type.

### Bibliography and reproducibility

- **Bibliography additions.** Kosiol 2004 (empirical codon-mutation matrix); Yurova Axelsson–Khrennikov 2024 (2-adic representation of the genetic code, second paper); Bonitz et al. 1980 (primary source for the *S. cerevisiae* mitochondrial tRNA vector, cited at §2.3.6); Kerscher et al. 2001 (primary source for the *Y. lipolytica* mitochondrial tRNA vector, cited at §2.3.6); Su et al. 2011; Atchley et al. 2005 in-text citation at the Serine / Factor 3 exploratory paragraph (§3.7).

- **Classical code-optimization citations added at Introduction §1.1** alongside Woese 1965 and Freeland–Hurst 1998: Novozhilov, Wolf & Koonin 2007; Novozhilov & Koonin 2009; Sella & Ardell 2006; Vetsigian, Woese & Goldenfeld 2006; Di Giulio 2005 (this-journal review of competing origin-of-the-code frameworks).

- **Bibliography validation pass.** Every entry in `references.bib` cross-checked for existence and metadata correctness against primary sources (DOI, journal, volume, pages, authors, diacritics): Kerscher 2001 `Br{\"a}ndt → Brandt` (spurious diacritic removed); Tsour 2026 gained volume 656, pages 506–515 (paper is no longer "in press"); Buschmann 2026 author capitalization `El-Hendi → el-Hendi` to match the primary source. All cited keys now resolve to bibliography entries and no orphan (uncited) entries remain in the .bib.

- **Cross-reference sweep.** §2.1 §S3 → §S2 (encoding sweep target); §2.3.1 §S3.1 → §S2.1 (dual-null comparison target); §3.4 §S8 → §S9 (topology clade-exclusion target); §3.5 §S13 → §S19.5, Figure S8 panel C (predictive-diagnostic target); §5.10 *H(3,4)* clade-exclusion paragraph quotes the matching primary-cell calculation ($H(3,4)$, $Δβ_0 > 0$, $N = 1280$, $K = 846$, $n = 24$, $x = 2$ → $p ≈ 4.2 × 10⁻⁹$), with the $Q_6$ / new-disconnection sensitivity variant ($≈ 3.6 × 10⁻¹¹$) called out separately; supplement pointers to §S13 / §S15 / §S16 / §S17 / §S18 updated to the current supplement numbering.

- **Reproducibility statement.** Downgraded from "bit-for-bit" to "within numerical / rendering tolerance"; transitive dependencies (BLAS routing, libgfortran ABI, Typst layout across versions) are not fully pinned.

- **Rendering.** Every R visualization script now reads statistics from `manuscript_stats.json` / per-analysis JSONs / T-CSVs rather than embedding hard-coded literals; Typst source-escape defects fixed; every JSON artefact parses under strict readers (R jsonlite, Typst `json()`).

### Reproduction

```
git clone https://github.com/biostochastics/codontopo && cd codontopo
git checkout v0.6.1
pip install -e ".[dev]"
bash scripts/build_publisher_release.sh
```

`scripts/build_publisher_release.sh` runs the analysis pipeline (`codon-topo all --seed=135325 --n=10000`), the provenance-emit and table-generation scripts, both R figure scripts, and the two Typst PDF compilations. Seed = 135325 across every Monte Carlo path.

---

## [v0.4.x] — earlier — Additional bridging analyses

Additional analyses that landed before the pre-publication QA/QC pass: ProtSub sensitivity (Supplement §S20), Walsh–Hadamard 2-adic spectral panel (§S21), non-power-representation discussion in §4.5, Slavov / Tsour SAAP figure (§S22), metric correlation analysis (§S20). See git log for per-commit detail.
