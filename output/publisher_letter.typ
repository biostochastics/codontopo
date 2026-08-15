// Cover letter to the BioSystems production editor for the
// pre-production corrections pass on BIOSYS-D-26-00689.
// Compile: typst compile output/publisher_letter.typ output/publisher_letter.pdf

#set page(
  paper: "us-letter",
  margin: (x: 1in, y: 1in),
)
#set text(font: "New Computer Modern", size: 11pt, lang: "en")
#set par(justify: true, leading: 0.65em, first-line-indent: 0pt, spacing: 1.1em)

#align(right)[
  #text(weight: "bold")[Sergey Kornilov, PhD] \
  Biostochastics, LLC \
  Seattle, WA, USA \
  #link("mailto:sergey@biostochastics.com")[sergey\@biostochastics.com]
]

#v(1em)

August 15, 2026

#v(0.5em)

Production Editor \
*BioSystems* \
Elsevier

#v(0.8em)

*Re: BIOSYS-D-26-00689 --- Pre-production corrections, second pass (v0.6.0)*  \
*"Robust error-minimization in the genetic code across physicochemical metrics and variant codes: a graph-theoretic analysis in GF(2)#super[6]"* \
(Clayworth & Kornilov, in press)

#v(0.4em)

Dear Production Editor,

Ahead of typesetting we ran a second internal QA/QC read of the manuscript, supplement, and pipeline, consolidating a small number of scientific corrections and a broader presentation sweep into a single v0.6.0 release. A per-item ledger with commit provenance is in `docs/publisher/2026-08-14-corrections-ledger.md` and a compact `CHANGELOG.md` covers both v0.5.1 and v0.6.0. The paper's four SUPPORTED scientific claims (cross-metric coloring optimality, per-table preservation across the 27 NCBI tables, $rho$-robustness, and topology-avoidance depletion) and the falsified KRAS--Fano result are unchanged in status and $p$-value. The tRNA-enrichment claim remains EXPLORATORY (as reclassified in v0.5.1 following the earlier MIS-enumeration correction).

*Substantive changes at v0.6.0.* Three items were flagged in the second QA/QC pass and are worth surfacing candidly rather than treating as routine.

+ *Null-model naming.* The coloring-optimality randomisation ensemble the code implements is a *quartet-pattern shuffle* (16 first-two-base codon quartets, each quartet's internal amino-acid pattern held atomic, patterns permuted across quartet slots with stop-containing quartets held fixed). The manuscript had labelled this "Freeland--Hurst block-preserving", which is close in spirit to the classical Haig--Hurst 1991 / Freeland--Hurst 1998 null but not identical: the classical construction preserves the 20 codon families (which codons belong to which family) and permutes the 20 amino-acid labels across those families ($20! approx 2.4 times 10^(18)$ possible codes), whereas the quartet-pattern shuffle preserves both the family membership and the quartet-slot structure ($16! approx 2.1 times 10^(13)$ possible codes). The two ensembles preserve non-nested structural properties. We have renamed the implemented null throughout the main text and supplement, and we now also report the classical Haig--Hurst AA-permutation null as a sensitivity companion (Supplement §S3.1, $n = 10{,}000$ per metric). All eight $p$-values (four physicochemical metrics under two nulls) pass Bonferroni at $alpha = 0.05/5 = 0.01$; for three of the four code-independent metrics the classical null gives a *more* extreme $p$-value than the quartet-pattern shuffle. The optimality claim is robust to the null choice.

+ *tRNA analysis population made explicit.* The 24 pairings that enter the Fisher--Stouffer enrichment analysis span 24 distinct organisms with mixed provenance: 15 of the 18 tRNAscan-SE-verified genomes populate the pairings, and the remaining 9 organism-slots use literature, GtRNAdb, or annotation counts. The three tRNAscan-verified genomes *not* in the pairings (_Blastocrithidia nonstop_ P57 and two _Mycoplasmoides_ species) are used only as mechanistic boundary cases in §4.3, because their reassignment routes (anticodon stem shortening; single-tRNA anticodon modification) act without changing tRNA gene counts and would be silent under a count-based Fisher test. Supplement §S10.1 now prints the complete 24-row input table with per-pairing $2 times 2$ counts, source class per organism, and topology status under both $Q_6$ and $H(3,4)$; the Fisher denominator convention (by-amino-acid sum over the *Std20* column of Table S8, excluding SeC / undetermined / suppressor / pseudogene rows) is stated explicitly.

+ *Data-plumbing correction to §3.5.* The conditional-logit clade-exclusion aggregate values (min / median / max $Delta"AICc"$(M1 $arrow.r$ M3)) were computed by the pipeline and written to `condlogit_clade_sensitivity.json`, but the manuscript template's default-value fallback (`.at("...", default: 0)`) meant that if the aggregate keys were absent from `manuscript_stats.json` the prose rendered zeros ("$Delta"AICc" >= 0$, median 0, max 0"). We now lift the aggregate keys explicitly, delete the fallback branches from the supplement, and add build-time `#assert` guards so a missing lift errors visibly rather than silently. The correct values (min $= 48.8$, median $= 100.7$, max $= 116.7$) are now rendered in §3.5.

*Presentation, framing, and reproducibility.* Numerous smaller items are covered in the ledger: §2.5 replaced with the discrete four-indicator structural-preservation index that Fig 5A actually uses (matching the §S12 body corrected in v0.5.1); §2.6 tRNA decision-gate wording restored to the originally specified worst-case-MIS criterion (which failed at $p = 0.123$; the median statistic is reported as a descriptive summary, not as a re-specified criterion); §4.3 retitled "Exploratory tRNA accommodation" and its compensation-for-disconnection language removed; §S18 cross-family multiple-testing statement corrected (the MIS-median $p = 0.046$ does not survive Bonferroni at $alpha = 0.00625$); ρ-sweep clarified (5 pre-specified inferential grid points, 6 additional plotted points as descriptive interpolation); candidate-universe Table S3 corrected (U3 no-op count 64 $arrow$ 61); Kosiol 2004 and the second Yurova Axelsson-Khrennikov 2024 paper now cited; reproducibility statement corrected to "within numerical/rendering tolerance", with `requirements.lock` and a `Dockerfile` shipped for the reference build. Several Typst source-escape rendering defects (literal `10{"^"}{-5}`, `1{","}280`, spacing after `e.g.`, doubled "Table Table S1", the Chan & Lowe DOI) are also fixed.

*Reproducibility.* A single `codon-topo all --seed=135325` regenerates every JSON and CSV artifact from the release tag; the Typst manuscript and supplement then render from those artifacts. Ruff format + check are clean across `src/` and `scripts/`. Release tag for what ships to production: `v0.6.0`.

*What ships.* The refreshed `manuscript.pdf`, `supplement.pdf`, the LaTeX-derivative source for Editorial Manager (`manuscript.tex` + `manuscript.bbl` + `references.bib` with BibTeX-safe diacritics), the figure TIFFs at 300 DPI LZW compression, this letter, the accompanying ledger, and the `requirements.lock` + `Dockerfile` for the reference build. All materials are internally consistent within the single-source-of-truth Typst render.

We appreciate your patience with these late-stage corrections and are happy to answer any questions before typesetting begins.

#v(1em)

Sincerely,

#v(2em)

Sergey Kornilov, PhD \
Corresponding author \
Biostochastics, LLC \
Seattle, WA, USA \
#link("mailto:sergey@biostochastics.com")[sergey\@biostochastics.com]

#v(0.4em)

#emph[On behalf of all authors.]

#v(1em)

#line(length: 100%, stroke: 0.4pt)

#v(0.4em)

*Attachments:*
- `docs/publisher/2026-08-14-corrections-ledger.md` --- full per-item ledger with commit provenance
- `CHANGELOG.md` --- release notes for v0.5.1 and v0.6.0
- Refreshed `manuscript.pdf`, `supplement.pdf`, `submission_em.zip`
- `requirements.lock`, `Dockerfile` --- pinned reference build environment
