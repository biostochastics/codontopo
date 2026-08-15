bioRxiv submission bundle — v0.5.0 (2026-08-14)
====================================================

Title:   Robust error-minimization in the genetic code across physicochemical
         metrics and variant codes: a graph-theoretic analysis in GF(2)^6
Authors: Paul Clayworth (Logocentricity Inc.),
         Sergey Kornilov (Biostochastics, LLC)
Corresp: sergey@biostochastics.com
Status:  In press at BioSystems (Elsevier), submission BIOSYS-D-26-00689
Release: v0.5.0-proofs-corrections-2026-08-14

Files in this bundle
--------------------

  manuscript.pdf                34 pages   Main text with figures embedded
                                           (Typst 0.14.2, 300 DPI, colour).

  supplementary_information.pdf 35 pages   Full supplement (§S1–§S22).
                                           +6 pages vs the R2 version — new
                                           Walsh spectral figure (§S16) and
                                           new Slavov SAAP figure (§S17)
                                           plus updated FigG_condlogit_diag.

  figures_highres/                          15 PNG + 15 TIFF at 300 DPI.
                                           Upload if bioRxiv asks for
                                           individual high-res figures.

  references.bib                            BibTeX source; informational.

  README.txt                                This file.

Upload checklist (bioRxiv)
--------------------------

  1. Manuscript file      -> manuscript.pdf
  2. Supplementary Info   -> supplementary_information.pdf
  3. Category             -> Bioinformatics (or Evolutionary Biology; pick primary)
  4. Subject area(s)      -> Genetics, Molecular Biology, Systems Biology
  5. Licence              -> CC BY 4.0 recommended
  6. Corresponding author -> Sergey Kornilov <sergey@biostochastics.com>
  7. Author list          -> Paul Clayworth (Logocentricity Inc., Austin, TX);
                             Sergey Kornilov (Biostochastics, LLC, Seattle, WA)
                             (Equal contribution — noted in the manuscript
                              frontmatter.)
  8. Prior/concurrent submission -> Accepted at Elsevier BioSystems
                             (BIOSYS-D-26-00689) as of the same date;
                             this bundle IS the version-of-record content.
  9. Funding              -> "None"
  10. Competing interests -> "None declared"

Changes vs any prior bioRxiv posting
------------------------------------

If this manuscript was previously posted to bioRxiv (e.g., the R2 version
or an earlier draft), the corrections below are all pre-production fixes
made during a self-audit pass after acceptance. None affects a scientific
claim, statistical significance level, or conclusion.

  * Figure 3 Panel B: coincident depth-calibration markers now jittered
    (Ala CUG-clade previously hidden under Ser CUG-clade at same coords).
  * Figure 3 Panel C: Q_6 subscript renders correctly (was literal text).
  * Figure 3 Panel D and supplement FigE: 24 pairing bars (was 18 with
    silently-stacked ranks).
  * Figure 4 Panel C: KRAS bar legend now "Falsified" (was "Pending").
  * Figure 5 Panel C: x-axis labels rotated, no overlap; p-values grammatical.
  * Supplement FigA / Figure 2 Panel A: empirical histogram (was Gaussian).
  * Supplement FigB / FigF / FigG_bit_position_bias / FigG_condlogit_diagnostics /
    FigC / FigS: various presentation and data-source corrections;
    see CHANGELOG.md in the release for the full list.
  * Two new supplement figures (Walsh §S16, Slavov §S17) added in response
    to reviewer requests that had shipped as tables only.
  * Per-table coloring MC rerun at authoritative n = 10,000 (was n = 1,000
    while the captions claimed 10,000). 26/27 significance count unchanged.
  * Data-plumbing fix: disconnection catalogue now correctly includes NCBI
    tables 28 (Condylostoma nuclear) and 32 (Balanophoraceae plastid),
    visible in Figure 1 Panel B.

Every corrected value is traceable to a specific commit SHA on the
public codontopo repository (release tag v0.5.0-proofs-corrections-
2026-08-14). See docs/publisher/2026-08-14-corrections-ledger.md for
per-item provenance.

Reproducibility
---------------

Both PDFs compiled from the same Typst 0.14.2 source that BioSystems'
production editor will typeset from (via LaTeX derivative). All
analyses are byte-reproducible from
  codon-topo all --output-dir=./output --seed=135325

Repository: https://github.com/biostochastics/codontopo
Release tag: v0.5.0-proofs-corrections-2026-08-14
