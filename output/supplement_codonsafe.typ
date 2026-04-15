// ============================================================
// Supplementary Material: Synthetic Genome Recoding Reanalysis
// Clayworth & Kornilov 2026
// ============================================================

#set page(paper: "a4", margin: (x: 2.5cm, y: 2.5cm), numbering: "S1")
#set text(font: "Libertinus Serif", size: 11pt)
#set par(justify: true, leading: 0.65em)
#set heading(numbering: "S1.1")
#show table.cell.where(y: 0): set text(weight: "bold")

#align(center)[
  #text(14pt, weight: "bold")[
    Supplementary Material: Cross-Study Reanalysis \
    of Synthetic Genome Recoding Outcomes
  ]
  #v(0.5em)
  #text(11pt)[Paul Clayworth and Sergey Kornilov]
  #v(1em)
]

= Dataset Inventory and Provenance <sec:supp-datasets>

Eight genome recoding datasets were fully analyzed, with one additional dataset (Nyerges et al., 2024) included for scope but pending quantitative table extraction from the published supplement.

#figure(
  table(
    columns: (auto, auto, auto, auto, auto, auto),
    align: (left, left, right, left, left, left),
    inset: 5pt,
    stroke: (x, y) => if y == 0 { (bottom: 0.7pt) } else { none },
    table.header([*Study*], [*Codons*], [*$n$*], [*Outcome*], [*Boundary var.*], [*Data source*]),
    [Robertson Syn57 (2025)], [4 Ser + 2 Ala + Stop], [60,240], [RNA-seq], [Ser=100%, Ala=0%], [Data S2/S8 GenBank],
    [Fredens Syn61 (2019)], [TCG/TCA$arrow.r$AGC/AGT], [18,218], [Binary fixes], [All cross], [Sup Data 3],
    [Ostrov (2016)], [7 codon types], [62,214], [Segment viability], [Mixed], [Table S4 + GenBank],
    [Napolitano (2016)], [AGR$arrow.r$CGN], [12,888], [SRZ covariates], [None], [Dataset S7],
    [Ochre (2025)], [TGA$arrow.r$TAA], [1,195], [Growth curves], [N/A (Stop)], [MOESM10],
    [Lajoie C321 (2013)], [UAG$arrow.r$UAA], [321], [Fitness data], [N/A (Stop)], [Table S26],
    [Frumkin (2018)], [CGU/CGC$arrow.r$CGG], [60], [Translation eff.], [None], [Dataset S1-S2],
    [Ding (2024)], [TCG (mammalian)], [Variable], [ncAA incorporation], [Yes], [Data S1-S3],
  ),
  caption: [
    Summary of analyzed datasets with data sources. Boundary variation indicates whether the dataset contains within-study variation in the connected-component boundary-crossing feature at $epsilon = 1$.
  ],
) <tbl:supp-datasets>

= Sensitivity Analysis: Deviation Filtering Threshold <sec:supp-sensitivity>

The filtering of genes with $>$10 CDS-level deviations removes structural rearrangements (predominantly _bioF_ with 311 changes). The null result for deviation directionality is robust across all tested thresholds.

#figure(
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto, auto),
    align: (center, right, right, right, right, right, right, right),
    inset: 5pt,
    stroke: (x, y) => if y == 0 { (bottom: 0.7pt) } else { none },
    table.header([*Cutoff*], [*$n$*], [*Removed*], [*Better*], [*Worse*], [*Same*], [*Binom. $p$*], [*Wilcoxon $p$*]),
    [$lt.eq 3$], [221], [506], [51], [45], [125], [0.62], [0.48],
    [$lt.eq 5$], [250], [477], [58], [56], [136], [0.56], [0.59],
    [$lt.eq 10$], [326], [401], [76], [87], [163], [0.83], [0.67],
    [$lt.eq 15$], [381], [346], [85], [103], [193], [0.90], [0.83],
    [$lt.eq 20$], [381], [346], [85], [103], [193], [0.90], [0.83],
    [$lt.eq 50$], [402], [325], [97], [112], [193], [0.88], [0.85],
    [All], [727], [0], [190], [326], [211], [1.00], [1.00],
  ),
  caption: [
    Sensitivity of the deviation directionality analysis to the per-gene filtering threshold. "Cutoff" = maximum deviations per gene retained. All binomial $p$-values are for the one-sided test (proportion better $> 0.5$); the two-sided test at the primary cutoff ($lt.eq 10$) gives $p = 0.43$, Clopper--Pearson 95% CI $[0.39, 0.55]$. The null result is robust across all thresholds.
  ],
) <tbl:supp-sensitivity>

= Ostrov Segment-Level Analysis <sec:supp-ostrov>

Among 87 segments in the Ostrov 57-codon genome design, 44 had post-deletion fitness data and 13 had lethal design exceptions. Bonferroni correction was applied for three simultaneous tests.

#figure(
  table(
    columns: (auto, auto, auto, auto),
    align: (left, right, right, right),
    inset: 5pt,
    stroke: (x, y) => if y == 0 { (bottom: 0.7pt) } else { none },
    table.header([*Test*], [*Statistic*], [*Raw $p$*], [*Bonf. $p$*]),
    [$rho$(essential recoded, DT)], [0.34], [0.022], [0.066],
    [$rho$(total recoded, DT)], [-0.10], [0.53], [1.00],
    [MW: problem vs normal (ess. recoded)], [$U$], [0.086], [0.26],
  ),
  caption: [
    Ostrov segment-level tests with Bonferroni correction for three simultaneous comparisons. DT = relative doubling time. MW = Mann--Whitney $U$ test. The essential-gene recoding load correlation is suggestive but not significant after correction.
  ],
) <tbl:supp-ostrov>

= GenBank Extraction Methods <sec:supp-genbank>

Codon positions were extracted by CDS-level comparison between parent and designed/final genomes using Biopython (Cock et al., 2009). CDS features were matched by `locus_tag` qualifier. The `/codon_start` qualifier was honored for correct reading-frame assignment; genes with `codon_start > 1` had the initial partial codon excluded from analysis.

For the Syn57 analysis, the MDS42 parent genome (Fredens et al. 2019, Supplementary Data 1; GenBank AP012306.1) was compared against the vS33A7 design genome (Robertson et al. 2025, Data file S2). The parent genome has 7,509 annotated features while the design genome has 110,025, resulting in 60,240 CDS-matched codon changes. The discrepancy from the paper-reported 101,302 reflects differences in CDS annotation completeness between the two genomes; the 60,240 positions represent the conservative intersection where both genomes have matching CDS annotations by locus tag.

For the design-to-final deviation analysis, the vS33A7 design genome (Data file S2) was compared against the verified Syn57 genome (Data file S8), yielding 727 CDS-level codon differences. Seven genes with $>$10 differences each (predominantly _bioF_ with 311 changes, likely a structural rearrangement or gene refactoring) were filtered as structural artifacts, leaving 326 genuine point deviations.

= Encoding Dependence of Hamming Distance <sec:supp-encoding>

Hamming distance in GF(2)#super[6] is encoding-dependent: the same pair of codons may have different Hamming distances under different base-to-bit bijections. By contrast, nucleotide edit distance (the number of nucleotide positions that differ) is encoding-invariant. For the Syn57 design deviations, we report Hamming distance under the default encoding ($C |-> 00$, $U |-> 01$, $A |-> 10$, $G |-> 11$). The observation that 83% of deviations are Hamming-1 (single-bit) moves under this encoding corresponds to deviations that change exactly one bit within one nucleotide position — a subset of single-nucleotide substitutions. All single-nucleotide substitutions are either Hamming-1 or Hamming-2 (within-nucleotide diagonal) in GF(2)#super[6]; the predominance of Hamming-1 moves indicates a preference for transitions or specific transversions under this encoding.

= Data Availability <sec:supp-data>

All analysis code is available in the `codon_topo.analysis.codonsafe` subpackage (version 0.3.0). The reproducible analysis pipeline can be run via:

```
python3.11 -m codon_topo.analysis.codonsafe.run_analyses
```

Raw supplementary data files were downloaded from publisher websites and are documented in `data/codonsafe/DATA_MANIFEST.md`. GenBank genome sequences are available under the accessions listed in the main text.
