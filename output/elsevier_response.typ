// ============================================================
// Response to Elsevier Production (BioSystems)
// Re: BIO 105919 / BIOSYS-D-26-00689
// Clayworth & Kornilov
// ============================================================
#set document(
  title: "Response to Elsevier BioSystems Production — BIO 105919 (BIOSYS-D-26-00689) — Corrected replacement files and QA/QC ledger",
  author: ("Sergey Kornilov", "Paul Clayworth"),
  keywords: (
    "BioSystems",
    "Elsevier production",
    "BIO 105919",
    "BIOSYS-D-26-00689",
    "response letter",
    "corrected PDF",
    "post-acceptance replacement",
    "QA/QC ledger",
  ),
  date: datetime(year: 2026, month: 8, day: 15),
)

#set page(paper: "a4", margin: (x: 2.5cm, y: 2.3cm))
#set text(font: "Libertinus Serif", size: 11pt, lang: "en")
#set par(justify: true, leading: 0.65em, spacing: 0.9em, first-line-indent: 0em)

// Header: sender + date
#grid(
  columns: (1fr, auto),
  align: (left, right),
  [
    #text(weight: "bold")[Sergey Kornilov, PhD] \
    Biostochastics, LLC \
    Seattle, WA, USA \
    #link("mailto:sergey@biostochastics.com")[sergey\@biostochastics.com] \
  ],
  [
    15 August 2026
  ],
)

#v(1.0em)

// Addressees
#text[
  Ms. S. Sreelakshmi \
  Data Administrator, Elsevier Production \
  #link("mailto:SreelakshmiS@elsevier.com")[SreelakshmiS\@elsevier.com] \
  #v(0.3em)
  Dr. Robert Prinz \
  Handling Editor, _BioSystems_ \
  Elsevier B.V.
]

#v(1.0em)

// Subject line
#text(weight: "bold")[
  Re: BIO 105919 — Article reference BIO_BIOSYS-D-26-00689 \
  "Robust error-minimization in the genetic code across physicochemical
  metrics and variant codes: a graph-theoretic analysis in GF(2)#super[6]"
]

#v(0.8em)

Dear Ms. Sreelakshmi and Dr. Prinz,

Thank you for accepting our article for publication in _BioSystems_ and for
your patience while we resolve the electronic-file issue flagged by
production. We are grateful for the care your team is taking with the
typesetting.

*Nature of the replacement.* We are supplying a corrected submission
bundle. The bundle addresses two categories of change simultaneously:

+ The immediate production issue: the previously supplied source PDF
  lacked figures and the reference list. The replacement `manuscript.tex`
  and `supplement.tex` embed every figure and render the complete
  reference list; both compile self-contained under pdflatex + elsarticle.

+ A pre-publication QA/QC pass carried out between acceptance and
  typesetting. Because that pass revised several displayed values,
  provenance statements, and framing sentences, the enclosed files
  contain post-acceptance author corrections in addition to the
  electronic-file repair. We therefore respectfully request the Handling
  Editor's approval of the revised accepted files. The full
  accepted-to-final diff is provided as `corrections_ledger.pdf`
  (attached); every user-visible change is enumerated there with its
  old value, new value, and justification.

*What did not change.* The scientific conclusions and the claim
hierarchy are unchanged. The four SUPPORTED headline findings
(cross-metric coloring optimality, per-table preservation across the 27
NCBI translation tables, ρ-robustness, and *H*(3, 4) topology-avoidance
depletion), the FALSIFIED KRAS–Fano prediction, the EXPLORATORY tRNA
enrichment classification, and the three REJECTED / two TAUTOLOGICAL
entries all retain their accepted status. The hierarchy total remains 15
claims.

*What was corrected.* The most consequential corrections are: (a) a
Monte Carlo tie-inclusion fix in the coloring-optimality null (Grantham
#emph[p] moves 0.006099 to 0.006199; other three metrics unchanged);
(b) a reconciliation of the #emph[S. cerevisiae] and #emph[Y. lipolytica]
mitochondrial tRNA vectors against primary literature (Bonitz #emph[et
al.] 1980 and Kerscher #emph[et al.] 2001, both added to the
bibliography); (c) a full Bron#sym.dash.en Kerbosch MIS enumeration for
the tRNA-enrichment robustness bound (the enumeration now returns 332
MIS of size 6 rather than 2, and the per-MIS Stouffer distribution
replaces the two-point best/worst summary); (d) a rewrite of the
main-text #sym.section 2.3.6 tRNA analysis population paragraph to
match the 24-row provenance table row-for-row; (e) a narrowing of the
primary #emph[H(3, 4)] clade-exclusion statement in #sym.section 3.4 to
$p < 10^(-3)$ (largest $p approx 2 times 10^(-4)$), with the #emph[Q]
#sub[6]/new-disconnection sensitivity retaining $p < 10^(-5)$; and (f)
a rewrite of Supplement #sym.section S18 as a descriptive
conservative-threshold sensitivity rather than a formal cross-family
Bonferroni correction (which is not well-defined across mixed
hypothesis-test / model-selection statistics). None of these
corrections change any claim's status, and the ledger states each
before/after value exactly.

*Reproducibility.* Both the corrected PDFs and every displayed statistic
are reproducible end-to-end from the public repository at tag `v0.6.1`
by one canonical command:

#block(
  inset: (left: 1.5em, right: 1.5em, y: 0.4em),
  fill: luma(245),
  radius: 3pt,
  breakable: false,
  [
    `git clone https://github.com/biostochastics/codontopo && cd codontopo` \
    `git checkout v0.6.1` \
    `pip install -e ".[dev]"` \
    `bash scripts/build_publisher_release.sh`
  ],
)

The wrapper script `build_publisher_release.sh` runs the analysis
pipeline (`codon-topo all` at seed 135325, *n* = 10 000), the
provenance-emit and table-generation scripts, both R figure scripts, and
the two Typst PDF compilations, in that order. The same command block
is quoted verbatim in `CHANGELOG.md`, in `README.md`, in the manuscript
Data-and-code availability paragraph, and in Supplement §S24, so no
ambiguity remains about which recipe reproduces the paper. Python
dependency versions are pinned in `requirements.lock` and a reference
build container is provided by `Dockerfile` at the repository root. A
SHA-256 manifest of the enclosed release artefacts is included as
`MANIFEST.sha256`.

*Further review.* If your copy-editors or the editorial office would
like us to prepare any additional derivatives — for instance a
single-column proof-ready PDF, a separate figure-only PDF, individual
TIFFs, or a merged manuscript-plus-supplement file — please let us know
and we will supply them promptly. We are equally happy to answer any
question that arises during typesetting or copy-editing, either by
email or by a brief call.

Thank you again for the opportunity to publish this work in _BioSystems_
and for your patience while we prepared the corrected files.

#v(0.6em)

With kind regards, on behalf of both authors,

#v(1.2em)

Sergey Kornilov \
_Corresponding author_ \
Biostochastics, LLC — Seattle, WA, USA

#v(0.6em)

#text(size: 9pt, fill: luma(90))[
  #text(weight: "bold")[Enclosures.] \
  1. `manuscript.pdf` — corrected self-contained PDF, all figures and references embedded. \
  2. `supplement.pdf` — supplementary material, self-contained. \
  3. `submission_bundle.zip` — corrected LaTeX source archive: `manuscript.tex`, `supplement.tex`, `elsarticle.cls`, `references.bib`, figure PNGs and TIFFs at 300 DPI, and the pdflatex-compiled PDFs. \
  4. `submission_em_v0.6.1/` — unpacked source of the same LaTeX bundle for reviewers preferring an unarchived layout. \
  5. `corrections_ledger.pdf` — pre-publication QA/QC accepted-to-final diff ledger. \
  6. `elsevier_response.pdf` — this letter. \
  7. `CHANGELOG.md` — consolidated release notes for the v0.6.1 release. \
  8. `2026-08-14-corrections-ledger.md` — Markdown source of the ledger. \
  9. `MANIFEST.sha256` — SHA-256 hashes of every enclosed file for integrity verification.
]
