// ============================================================
// Pre-publication QA/QC Ledger — BIOSYS-D-26-00689 (BIO 105919)
// Auto-generated from docs/publisher/2026-08-14-corrections-ledger.md
// via pandoc, then wrapped as Typst so all Unicode / Greek / math
// glyphs render without missing-character warnings.
// ============================================================
#set document(
  title: "BIOSYS-D-26-00689 (BIO 105919) — Pre-publication QA/QC Ledger",
  author: ("Paul Clayworth", "Sergey Kornilov"),
  keywords: (
    "BioSystems",
    "BIOSYS-D-26-00689",
    "BIO 105919",
    "pre-publication QA/QC",
    "corrections ledger",
    "figures",
    "tables",
    "statistics",
  ),
  date: datetime(year: 2026, month: 8, day: 15),
)

// Landscape A4 so that the wide correction-tables (ID / Change / Old value /
// New value / Justification / Commit) fit on one row without wrapping.
#set page(
  paper: "a4",
  flipped: true,
  margin: (x: 1.8cm, y: 1.8cm),
  numbering: "1",
  header: context {
    if counter(page).get().first() > 1 [
      #set text(size: 9pt, fill: luma(120))
      Pre-publication QA/QC Ledger — BIOSYS-D-26-00689 (BIO 105919)
      #h(1fr) #counter(page).display()
    ]
  },
)
#set text(font: "Libertinus Serif", size: 10pt, lang: "en")
#set par(justify: true, leading: 0.62em, spacing: 0.6em)
#set heading(numbering: none)
#show heading.where(level: 1): it => {
  v(1.2em, weak: false); text(size: 16pt, weight: "bold", it.body); v(0.5em, weak: false)
}
#show heading.where(level: 2): it => {
  v(0.9em, weak: false); text(size: 13pt, weight: "bold", it.body); v(0.3em, weak: false)
}
#show heading.where(level: 3): it => {
  v(0.6em, weak: false); text(size: 11pt, weight: "bold", it.body); v(0.25em, weak: false)
}
#show link: it => text(fill: rgb("#1a5490"), it)

// Small tabular figures inside tables
#show table: set text(number-width: "tabular", size: 9pt)
#show table.cell.where(y: 0): set text(weight: "bold")
#show figure.where(kind: table): set block(breakable: true)

// Body: include pandoc-generated content
#include "ledger_body.typ"
