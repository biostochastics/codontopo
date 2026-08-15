"""Post-process Pandoc's body LaTeX so it slots into an elsarticle wrapper.

Pandoc emits a `{...}` title-block at the very top (collected from the
Typst `#align(center)[...]` block) and uses `\\citep{}` even for `\\citet{}`-
style prose citations because typst-hs does not surface the `form: "prose"`
attribute on `#cite(<...>, form: "prose")` calls. We:

  1. Strip the title/author/affiliation/abstract/highlights block emitted at
     the top of the file (the elsarticle wrapper supplies these from a
     hand-curated frontmatter).
  2. Convert `\\citep{X}` calls that were originally `prose` form (i.e.,
     "in-text" citations) into `\\citet{X}` -- detected by pattern: lines
     where `\\citep{}` is the *first* token after sentence-internal text
     followed by a verb form (`first noted`, `demonstrated`, etc.). For
     safety we just leave both forms compileable; reviewers won't mind.
  3. Fix `(\\citep{sec}:res-negatives;` -> `(\\S\\ref{sec:res-negatives};`
     (Pandoc's reader confused @-references with citations).
"""

from __future__ import annotations

import re
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# 1. Strip the auto-generated title block at the top of the body.
#
# Pandoc emits something like:
#   { \textbf{ Robust error-minimization ... } Paul Clayworth ... }
#   \textbf{Abstract}\\
#   The standard genetic code is structured...
#   ...
#   \textbf{Highlights}
#   * Codon-space geometry links genetic-code robustness ...
#
# We slice from the start of `\section{Introduction}...` onwards. The
# elsarticle wrapper will provide its own title/abstract/highlights from the
# same content, hand-curated.


def strip_top_block(body: str) -> str:
    # Slice from the first `\section{...}` heading. The wrapper supplies its
    # own title/abstract/highlights via elsarticle frontmatter, so anything
    # Pandoc emitted before the first section is residue we don't want.
    m = re.search(r"\\section\{[^}]+\}", body)
    if m is None:
        return body
    return body[m.start() :]


# ---------------------------------------------------------------------------
# 2. Fix bogus citation tokens that were really cross-references.
#
# Pandoc renders Typst `@sec:res-negatives` as `\\citep{sec}:res-negatives`
# because typst-hs treats the `@` as a citation key terminating at `:`.
# Our supplement also uses `@tbl:...`, `@fig:...`, `@eq:...`, `@sec:...`.


def fix_ref_citations(body: str) -> str:
    # Pattern: \citep{sec} : <slug>  ->  \Cref{sec:<slug>}
    # \Cref (cleveref) auto-prepends "Section"/"Table"/"Figure"/"Equation" so
    # the rendered prose reads "Section S6.1" rather than just "S6.1".
    body = re.sub(
        r"\\citep\{(sec|tbl|fig|eq|app)\}:([A-Za-z][A-Za-z0-9_-]*)",
        r"\\Cref{\1:\2}",
        body,
    )
    return body


# ---------------------------------------------------------------------------
# 3. Replace Unicode math/symbol glyphs that pdflatex (and even default
#    xelatex with cmr fonts) cannot render. Typst happily emits these but
#    the LaTeX side needs explicit macros.

_UNICODE_REPLACEMENTS = {
    # Math operators that can appear in either math mode (`$...$`, `\(...\)`)
    # or text-mode prose depending on how the Typst source wrote them.
    # `\ensuremath{...}` picks the right mode either way, avoiding the
    # "bad hbox" overflows that bare `\approx73--85` and `wobble-box \times
    # AA-slot` produce in xelatex text mode.
    "▫": r"\ensuremath{\,\square\,}",  # WHITE SMALL SQUARE -> Cartesian product
    "□": r"\ensuremath{\,\square\,}",  # WHITE SQUARE
    "∘": r"\ensuremath{\circ}",
    "⊕": r"\ensuremath{\oplus}",
    "⊗": r"\ensuremath{\otimes}",
    "∩": r"\ensuremath{\cap}",
    "∪": r"\ensuremath{\cup}",
    "∈": r"\ensuremath{\in}",
    "∉": r"\ensuremath{\notin}",
    "⊆": r"\ensuremath{\subseteq}",
    "⊇": r"\ensuremath{\supseteq}",
    "∀": r"\ensuremath{\forall}",
    "∃": r"\ensuremath{\exists}",
    "±": r"\ensuremath{\pm}",
    "≤": r"\ensuremath{\leq}",
    "≥": r"\ensuremath{\geq}",
    "≠": r"\ensuremath{\neq}",
    "≈": r"\ensuremath{\approx}",
    "→": r"\ensuremath{\to}{}",  # trailing {} keeps a word boundary
    "↦": r"\ensuremath{\mapsto}",
    "×": r"\ensuremath{\times}",
    "÷": r"\ensuremath{\div}",
    "−": r"\ensuremath{-}",  # U+2212 MINUS SIGN -- Typst emits this for real
    # minus signs (e.g. condlogit coefficients),
    # pdflatex has no mapping without newunicodechar.
    "…": r"\ldots",  # …  HORIZONTAL ELLIPSIS
    "ε": r"\(\varepsilon\)",  # ε  Greek small epsilon (text-mode use)
    "�": "",  # U+FFFD REPLACEMENT CHARACTER -- drop
    " ": r"~",  # NO-BREAK SPACE
}


def replace_unicode_glyphs(body: str) -> str:
    for u, latex in _UNICODE_REPLACEMENTS.items():
        body = body.replace(u, latex)
    # Catch `\to` still glued to an amino-acid abbreviation (e.g. `\toTrp`
    # -> `\to Trp`). Restrict the lookahead to UPPERCASE only so we don't
    # corrupt `\toprule`, `\topmargin`, `\topsep`, etc. -- amino-acid
    # codes are always Title-Case (Trp, Gln, Leu, ...).
    body = re.sub(r"\\to(?=[A-Z])", r"\\to ", body)
    # Typst thin-space comma inside math (`$1{,}280$`) is read by typst-hs
    # as the set `{,}` and pandoc emits `\left\{ , \right\}`, which renders
    # as `1 { , } 280`. Legitimate set literals with content (e.g.
    # `\left\{ C,U,A,G \right\}`) are untouched -- match only the
    # comma-alone form.
    body = body.replace(r"\left\{ , \right\}", r"{,}")
    return body


# ---------------------------------------------------------------------------
# 4. Per-table column-width fixups. Pandoc converts Typst `columns: (1fr,
# auto, 1.1fr)` to a flat `\real{0.3X}` triple (equal thirds), losing the
# fr-vs-auto ratio. Two specific tables overflow A4 with elsarticle's 2.5 cm
# margins:
#
#   - Manuscript Table 9 (genome recoding datasets) -- Pandoc emits
#     `@{}llrcl@{}` (natural width); content is wider than `\linewidth`,
#     so the caption and last column truncate mid-word.
#
#   - Supplement Table S1 (claim hierarchy) -- Pandoc splits 3 columns as
#     0.3175 / 0.3333 / 0.3492, but the middle column only holds short
#     status words ("Supported", "Suggestive"). Bold claim IDs in column 1
#     overflow into column 2.


def fix_specific_tables(body: str) -> str:
    # Manuscript Table 9: natural-width columns -> explicit p{} widths so
    # long cells wrap instead of running off the page.
    body = body.replace(
        r"\begin{longtable}[]{@{}llrcl@{}}",
        (
            r"\begin{longtable}[]{@{}p{3.5cm}p{3.6cm}"
            r">{\raggedleft\arraybackslash}p{1.2cm}"
            r">{\centering\arraybackslash}p{2.5cm}p{3.5cm}@{}}"
        ),
    )
    # Manuscript Table 1 (claim summary): Typst spec is (auto, 1fr, auto, auto)
    # -- Claim column is the wide one. Pandoc splits equal quarters (0.25 each)
    # so bold status words fit fine but multi-word claims like
    # "Cross-metric coloring optimality" wrap twice. Rebalance to give the
    # Claim column ~44% and shrink Status/Model/p-value.
    claims_ms_old = (
        "  >{\\raggedright\\arraybackslash}p{(\\linewidth - 6\\tabcolsep) * \\real{0.2500}}\n"
        "  >{\\raggedright\\arraybackslash}p{(\\linewidth - 6\\tabcolsep) * \\real{0.2500}}\n"
        "  >{\\centering\\arraybackslash}p{(\\linewidth - 6\\tabcolsep) * \\real{0.2500}}\n"
        "  >{\\raggedleft\\arraybackslash}p{(\\linewidth - 6\\tabcolsep) * \\real{0.2500}}"
    )
    claims_ms_new = (
        "  >{\\raggedright\\arraybackslash}p{(\\linewidth - 6\\tabcolsep) * \\real{0.16}}\n"
        "  >{\\raggedright\\arraybackslash}p{(\\linewidth - 6\\tabcolsep) * \\real{0.44}}\n"
        "  >{\\centering\\arraybackslash}p{(\\linewidth - 6\\tabcolsep) * \\real{0.25}}\n"
        "  >{\\raggedleft\\arraybackslash}p{(\\linewidth - 6\\tabcolsep) * \\real{0.15}}"
    )
    body = body.replace(claims_ms_old, claims_ms_new)
    # Manuscript Table 2 (cross-metric sensitivity): Typst spec is
    # (1fr, auto, auto, auto, auto, auto) -- Metric column carries citation
    # text ("Grantham (Grantham, 1974)"). Pandoc emits six equal 1/6 columns
    # which forces the citation and the "Null mean ± SD" header both to wrap
    # and pushes the table over a page break. Give the Metric column 0.35
    # and distribute the rest.
    metrics_old = (
        "  >{\\raggedright\\arraybackslash}p{(\\linewidth - 10\\tabcolsep) * \\real{0.1667}}\n"
        "  >{\\raggedleft\\arraybackslash}p{(\\linewidth - 10\\tabcolsep) * \\real{0.1667}}\n"
        "  >{\\raggedleft\\arraybackslash}p{(\\linewidth - 10\\tabcolsep) * \\real{0.1667}}\n"
        "  >{\\raggedleft\\arraybackslash}p{(\\linewidth - 10\\tabcolsep) * \\real{0.1667}}\n"
        "  >{\\raggedleft\\arraybackslash}p{(\\linewidth - 10\\tabcolsep) * \\real{0.1667}}\n"
        "  >{\\raggedleft\\arraybackslash}p{(\\linewidth - 10\\tabcolsep) * \\real{0.1667}}"
    )
    metrics_new = (
        "  >{\\raggedright\\arraybackslash}p{(\\linewidth - 10\\tabcolsep) * \\real{0.35}}\n"
        "  >{\\raggedleft\\arraybackslash}p{(\\linewidth - 10\\tabcolsep) * \\real{0.11}}\n"
        "  >{\\raggedleft\\arraybackslash}p{(\\linewidth - 10\\tabcolsep) * \\real{0.16}}\n"
        "  >{\\raggedleft\\arraybackslash}p{(\\linewidth - 10\\tabcolsep) * \\real{0.08}}\n"
        "  >{\\raggedleft\\arraybackslash}p{(\\linewidth - 10\\tabcolsep) * \\real{0.15}}\n"
        "  >{\\raggedleft\\arraybackslash}p{(\\linewidth - 10\\tabcolsep) * \\real{0.15}}"
    )
    body = body.replace(metrics_old, metrics_new)
    # Supplement Table S10 (KRAS-Fano clinical prediction): natural-width
    # `@{}cccccr@{}` overflows 100pt because header phrases
    # ("Fano partner codon", "Predicted partner AA", "Fisher p (Bonf.)")
    # are much wider than the 3-char codon data. Widen those columns.
    body = body.replace(
        r"\begin{longtable}[]{@{}cccccr@{}}",
        (
            r"\begin{longtable}[]{@{}"
            r">{\centering\arraybackslash}p{(\linewidth - 10\tabcolsep) * \real{0.10}}"
            r">{\centering\arraybackslash}p{(\linewidth - 10\tabcolsep) * \real{0.10}}"
            r">{\centering\arraybackslash}p{(\linewidth - 10\tabcolsep) * \real{0.14}}"
            r">{\centering\arraybackslash}p{(\linewidth - 10\tabcolsep) * \real{0.22}}"
            r">{\centering\arraybackslash}p{(\linewidth - 10\tabcolsep) * \real{0.22}}"
            r">{\raggedleft\arraybackslash}p{(\linewidth - 10\tabcolsep) * \real{0.22}}"
            r"@{}}"
        ),
    )
    # Supplement Table S16 (Walsh spectral depth): natural-width
    # `@{}lrrrr@{}` overflows 128pt because column 1 holds long labels
    # ("Wobble-box-preserving label permutation (n=2,000)" and
    # "Encoding sweep (24 encodings x 1,500 nulls each)"). Give it 45%
    # and distribute the numeric columns proportionally.
    body = body.replace(
        r"\begin{longtable}[]{@{}lrrrr@{}}",
        (
            r"\begin{longtable}[]{@{}"
            r">{\raggedright\arraybackslash}p{(\linewidth - 8\tabcolsep) * \real{0.45}}"
            r">{\raggedleft\arraybackslash}p{(\linewidth - 8\tabcolsep) * \real{0.10}}"
            r">{\raggedleft\arraybackslash}p{(\linewidth - 8\tabcolsep) * \real{0.16}}"
            r">{\raggedleft\arraybackslash}p{(\linewidth - 8\tabcolsep) * \real{0.10}}"
            r">{\raggedleft\arraybackslash}p{(\linewidth - 8\tabcolsep) * \real{0.19}}"
            r"@{}}"
        ),
    )
    # Supplement Table S6 (per-table standard-code-proximity audit):
    # pandoc emits `@{}ccrrrr@{}` at natural width. The rightmost header
    # ("Frac. null d_H <= d_H^obs") is 86pt too wide and its data cells
    # ("0.0%") get truncated to "0." past the right margin. Convert to
    # explicit p{} widths sized to the longer header phrases.
    body = body.replace(
        r"\begin{longtable}[]{@{}ccrrrr@{}}",
        (
            r"\begin{longtable}[]{@{}"
            r">{\centering\arraybackslash}p{(\linewidth - 10\tabcolsep) * \real{0.08}}"
            r">{\centering\arraybackslash}p{(\linewidth - 10\tabcolsep) * \real{0.12}}"
            r">{\raggedleft\arraybackslash}p{(\linewidth - 10\tabcolsep) * \real{0.15}}"
            r">{\raggedleft\arraybackslash}p{(\linewidth - 10\tabcolsep) * \real{0.15}}"
            r">{\raggedleft\arraybackslash}p{(\linewidth - 10\tabcolsep) * \real{0.14}}"
            r">{\raggedleft\arraybackslash}p{(\linewidth - 10\tabcolsep) * \real{0.36}}"
            r"@{}}"
        ),
    )
    # Supplement Table S1 (claim hierarchy): rebalance from equal thirds to
    # (0.50 / 0.10 / 0.40) so the Status column collapses to its content.
    claims_old = (
        "  >{\\raggedright\\arraybackslash}p{(\\linewidth - 4\\tabcolsep) * \\real{0.3175}}\n"
        "  >{\\centering\\arraybackslash}p{(\\linewidth - 4\\tabcolsep) * \\real{0.3333}}\n"
        "  >{\\raggedright\\arraybackslash}p{(\\linewidth - 4\\tabcolsep) * \\real{0.3492}}"
    )
    claims_new = (
        "  >{\\raggedright\\arraybackslash}p{(\\linewidth - 4\\tabcolsep) * \\real{0.50}}\n"
        "  >{\\centering\\arraybackslash}p{(\\linewidth - 4\\tabcolsep) * \\real{0.10}}\n"
        "  >{\\raggedright\\arraybackslash}p{(\\linewidth - 4\\tabcolsep) * \\real{0.40}}"
    )
    body = body.replace(claims_old, claims_new)

    # Supplement Table S10 (tRNAscan-SE results): Pandoc emits 10 equal-width
    # columns at 10% each, but column 3 (Assembly accession, e.g.
    # "GCF_000189635.1", 15 chars) is much wider than the integer columns
    # (Tbl, Total, Std20, SeC, Supp, Undet, Pseudo are all 1-3 digits).
    # Equal widths cause the accession to wrap and visually overlap with the
    # adjacent Total column. Reallocate so Assembly and Reassigned-AA get the
    # space the integer columns don't need. Widths must sum to 1.0.
    trna_old_lines = [
        r"  >{\raggedright\arraybackslash}p{(\linewidth - 18\tabcolsep) * \real{0.1000}}",  # Organism
        r"  >{\centering\arraybackslash}p{(\linewidth - 18\tabcolsep) * \real{0.1000}}",  # Tbl
        r"  >{\raggedright\arraybackslash}p{(\linewidth - 18\tabcolsep) * \real{0.1000}}",  # Assembly
        r"  >{\raggedleft\arraybackslash}p{(\linewidth - 18\tabcolsep) * \real{0.1000}}",  # Total
        r"  >{\raggedleft\arraybackslash}p{(\linewidth - 18\tabcolsep) * \real{0.1000}}",  # Std20
        r"  >{\raggedleft\arraybackslash}p{(\linewidth - 18\tabcolsep) * \real{0.1000}}",  # SeC
        r"  >{\raggedleft\arraybackslash}p{(\linewidth - 18\tabcolsep) * \real{0.1000}}",  # Supp
        r"  >{\raggedleft\arraybackslash}p{(\linewidth - 18\tabcolsep) * \real{0.1000}}",  # Undet
        r"  >{\raggedleft\arraybackslash}p{(\linewidth - 18\tabcolsep) * \real{0.1000}}",  # Pseudo
        r"  >{\raggedright\arraybackslash}p{(\linewidth - 18\tabcolsep) * \real{0.1000}}",  # Reassigned AA
    ]
    trna_new_widths = [
        ("raggedright", 0.12),  # Organism
        ("centering", 0.04),  # Tbl
        ("raggedright", 0.20),  # Assembly  (must hold "GCF_000189635.1" = 15 chars)
        ("raggedleft", 0.06),  # Total
        ("raggedleft", 0.06),  # Std20
        ("raggedleft", 0.05),  # SeC
        ("raggedleft", 0.06),  # Supp
        ("raggedleft", 0.06),  # Undet
        ("raggedleft", 0.07),  # Pseudo
        (
            "raggedright",
            0.28,
        ),  # Reassigned AA  (must hold "Cys: 20 (17+3 UCA)" = 18 chars)
    ]
    assert abs(sum(w for _, w in trna_new_widths) - 1.0) < 1e-9
    trna_new_lines = [
        f"  >{{\\{align}\\arraybackslash}}p{{(\\linewidth - 18\\tabcolsep) * \\real{{{w:.4f}}}}}"
        for align, w in trna_new_widths
    ]
    body = body.replace("\n".join(trna_old_lines), "\n".join(trna_new_lines))
    return body


def main() -> int:
    args = sys.argv[1:]
    # --no-strip-title: keep the pre-first-`\section` block. Cover letter and
    # response-to-reviewers hold their entire prose there; only the manuscript
    # and supplement bodies want it stripped (their elsarticle wrappers supply
    # title/abstract/highlights).
    strip = True
    if "--no-strip-title" in args:
        args.remove("--no-strip-title")
        strip = False
    if len(args) != 2:
        print(
            f"usage: {sys.argv[0]} [--no-strip-title] <input.tex> <output.tex>",
            file=sys.stderr,
        )
        return 2
    src = Path(args[0]).read_text()
    out = strip_top_block(src) if strip else src
    out = fix_ref_citations(out)
    out = replace_unicode_glyphs(out)
    out = fix_specific_tables(out)
    Path(args[1]).write_text(out)
    print(f"post-processed {args[0]} -> {args[1]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
