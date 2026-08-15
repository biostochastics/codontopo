#!/usr/bin/env bash
# build_tex.sh -- one-shot Typst-to-LaTeX-to-PDF pipeline for JTB submission.
#
# Reads:
#   output/manuscript.typ
#   output/supplement.typ
#   output/cover_letter.typ
#   output/manuscript_stats.json   (drives all inline statistics)
#   output/references.bib
#
# Writes (under output-tex/):
#   *-frozen.typ                   -- Python-evaluated, set/show-stripped
#   *_body_raw.tex                 -- Pandoc output
#   *_body.tex                     -- post-processed (Unicode, citations, refs)
#   manuscript.pdf  supplement.pdf cover_letter.pdf
#
# Pre-reqs: pandoc 3.1+, xelatex (TeX Live), python 3.11+ in .venv-tex/
# Run from project root: ./scripts/build_tex.sh

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"

if [[ ! -d .venv-tex ]]; then
    echo "creating .venv-tex/..."
    python3.11 -m venv .venv-tex
fi
# shellcheck disable=SC1091
source .venv-tex/bin/activate

mkdir -p output-tex
cp -f output/references.bib output-tex/references.bib
# freeze_typst.py evaluates most #str/#calc.round references but leaves the
# top-level `#let stats = json("manuscript_stats.json")` in place, so pandoc
# still needs the JSON on disk next to the frozen typ. Refresh it every build
# to avoid drift from the canonical copy in output/.
cp -f output/manuscript_stats.json output-tex/manuscript_stats.json
# Editorial-Manager pdflatex wrappers live in templates/ so they survive
# a fresh clone; the XeLaTeX wrappers (main.tex, supplement.tex,
# cover_letter.tex, response_to_reviewers.tex) are still maintained
# manually in output-tex/.
cp -f templates/main_em.tex        output-tex/main_em.tex
cp -f templates/supplement_em.tex  output-tex/supplement_em.tex

# 1. Freeze each Typst file -- evaluates JSON-driven inline stats and helpers
echo "=== freezing Typst sources ==="
python scripts/freeze_typst.py \
    output/manuscript.typ \
    output/manuscript_stats.json \
    output-tex/manuscript-frozen.typ
python scripts/freeze_typst.py \
    output/supplement.typ \
    output/manuscript_stats.json \
    output-tex/supplement-frozen.typ
python scripts/freeze_typst.py \
    output/cover_letter.typ \
    output/manuscript_stats.json \
    output-tex/cover_letter-frozen.typ
python scripts/freeze_typst.py \
    output/response_to_reviewers.typ \
    output/manuscript_stats.json \
    output-tex/response_to_reviewers-frozen.typ

# 2. Pandoc Typst -> LaTeX (natbib for elsarticle compatibility)
echo "=== pandoc Typst -> LaTeX ==="
cd output-tex
pandoc manuscript-frozen.typ --from=typst --to=latex --natbib \
    -o manuscript_body_raw.tex
pandoc supplement-frozen.typ --from=typst --to=latex --natbib \
    -o supplement_body_raw.tex
pandoc cover_letter-frozen.typ --from=typst --to=latex \
    -o cover_letter_body_raw.tex
pandoc response_to_reviewers-frozen.typ --from=typst --to=latex \
    -o response_to_reviewers_body_raw.tex
cd "$PROJECT_ROOT"

# 3. Post-process LaTeX bodies: strip auto title block, fix \citep{sec}->\ref,
#    replace Unicode glyphs with macro forms, force `\to` word boundary,
#    fix Typst thin-space `{,}` mis-read as `\left\{,\right\}`.
echo "=== post-processing ==="
python scripts/postprocess_tex.py \
    output-tex/manuscript_body_raw.tex \
    output-tex/manuscript_body.tex
python scripts/postprocess_tex.py \
    output-tex/supplement_body_raw.tex \
    output-tex/supplement_body.tex
python scripts/postprocess_tex.py --no-strip-title \
    output-tex/cover_letter_body_raw.tex \
    output-tex/cover_letter_body.tex
python scripts/postprocess_tex.py --no-strip-title \
    output-tex/response_to_reviewers_body_raw.tex \
    output-tex/response_to_reviewers_body.tex

# 4. xelatex with bibtex.  -f forces continuation past warnings; the wrappers
#    use fontspec + newunicodechar so xelatex (not pdflatex) is required.
echo "=== compile LaTeX ==="
cd output-tex
for doc in main supplement cover_letter response_to_reviewers; do
    rm -f "${doc}".{aux,bbl,blg,log,out,fls,fdb_latexmk,toc,xdv,pdf}
done
# latexmk -f returns nonzero when forced past warnings (missing chars, etc).
# We accept that and check the PDF artifact afterwards instead of the exit code.
latexmk -xelatex -interaction=nonstopmode -f main.tex || true
latexmk -xelatex -interaction=nonstopmode -f supplement.tex || true
latexmk -xelatex -interaction=nonstopmode -f cover_letter.tex || true
latexmk -xelatex -interaction=nonstopmode -f response_to_reviewers.tex || true
cd "$PROJECT_ROOT"

# 3b. Derive Editorial-Manager-adapted body files.
#
#     EM constraints handled here:
#
#     1. Flat working directory -- strip `figures/` prefix.
#     2. TIFF-only uploads for Figure items, but pdflatex cannot
#        \includegraphics{...} a .tif file (Elsevier support docs +
#        thelatexlab.com/blog/submit-latex-to-elsevier-editorial-manager
#        both confirm .tif is stored but not consumable by \includegraphics).
#        We therefore rewrite every \pandocbounded{\includegraphics[opts]{FigN_...png}}
#        into \emfig{FigN_...}  -- the wrappers define \emfig to render
#        a placeholder box that names the TIFF, so the typesetter
#        knows which upload belongs at each position.
echo "=== deriving EM-adapted body files ==="
python3.11 << 'PY'
import re
from pathlib import Path

# Match \pandocbounded{ \includegraphics[opts]{figures/BASENAME.png} } (and
# variants without \pandocbounded or without the figures/ prefix).
PAT = re.compile(
    r"\\pandocbounded\{\s*"
    r"\\includegraphics(?:\[[^\]]*\])?\{(?:figures/)?([A-Za-z0-9_.-]+?)\.png\}"
    r"\s*\}"
)
PAT_BARE = re.compile(
    r"\\includegraphics(?:\[[^\]]*\])?\{(?:figures/)?([A-Za-z0-9_.-]+?)\.png\}"
)
for src_name, dst_name in [
    ("output-tex/manuscript_body.tex", "output-tex/manuscript_body_em.tex"),
    ("output-tex/supplement_body.tex", "output-tex/supplement_body_em.tex"),
]:
    text = Path(src_name).read_text()
    n1 = 0
    def sub(m):
        global n1
        n1 += 1
        # LaTeX text-mode `_` is special (subscript). Escape it here so
        # the filename renders correctly inside \texttt{...} in \emfig.
        escaped = m.group(1).replace("_", r"\_")
        return r"\emfig{" + escaped + "}"
    text = PAT.sub(sub, text)
    text = PAT_BARE.sub(sub, text)
    Path(dst_name).write_text(text)
    print(f"  {src_name}: {n1} \\includegraphics -> \\emfig substitutions")
PY

# 4b. Produce single-file self-contained submission .tex for editorial
#     systems (Elsevier Editorial Manager, arXiv-flat) that pick a "main"
#     .tex file and don't run a full latexmk with \input{...} resolution.
#     The submission_*.tex files below can be uploaded verbatim.
echo "=== flattening for editorial systems ==="
python scripts/flatten_tex.py \
    output-tex/main.tex \
    manuscript_body \
    output-tex/submission_manuscript.tex
python scripts/flatten_tex.py \
    output-tex/supplement.tex \
    supplement_body \
    output-tex/submission_supplement.tex

# 4c. Assemble two upload-ready staging directories + zip files, so the
#     author can hand either package to Editorial Manager without hunting
#     for the right files:
#
#       submission_flat/       -- Option A: one self-contained .tex per
#                                  document, works when EM picks a "main"
#                                  file and does not resolve \input.
#       submission_modular/    -- Option B: wrapper + body pair, works
#                                  when EM runs a full latexmk build.
#
#     Each directory bundles the .tex sources, references.bib, figures/,
#     and the freshly compiled PDFs for reviewer / auditor reference.
echo "=== assembling submission packages ==="
FLAT_DIR="output-tex/submission_flat"
MOD_DIR="output-tex/submission_modular"
rm -rf "$FLAT_DIR" "$MOD_DIR"
mkdir -p "$FLAT_DIR/figures" "$MOD_DIR/figures"

# Flat package: rename to bare document names for EM cleanliness.
cp output-tex/submission_manuscript.tex   "$FLAT_DIR/manuscript.tex"
cp output-tex/submission_supplement.tex   "$FLAT_DIR/supplement.tex"
python scripts/flatten_tex.py \
    output-tex/cover_letter.tex          cover_letter_body \
    "$FLAT_DIR/cover_letter.tex"
python scripts/flatten_tex.py \
    output-tex/response_to_reviewers.tex response_to_reviewers_body \
    "$FLAT_DIR/response_to_reviewers.tex"
cp output-tex/references.bib             "$FLAT_DIR/"
cp output-tex/main.pdf                   "$FLAT_DIR/manuscript.pdf"
cp output-tex/supplement.pdf             "$FLAT_DIR/supplement.pdf"
cp output-tex/cover_letter.pdf           "$FLAT_DIR/cover_letter.pdf"
cp output-tex/response_to_reviewers.pdf  "$FLAT_DIR/response_to_reviewers.pdf"
cp output/figures/*.png                  "$FLAT_DIR/figures/"

# Modular package: keep the wrapper + body pair as-is.
for f in main.tex manuscript_body.tex \
         supplement.tex supplement_body.tex \
         cover_letter.tex cover_letter_body.tex \
         response_to_reviewers.tex response_to_reviewers_body.tex \
         references.bib \
         main.pdf supplement.pdf cover_letter.pdf response_to_reviewers.pdf; do
    cp "output-tex/$f" "$MOD_DIR/"
done
# For EM clarity: name the main-PDF as `manuscript.pdf` alongside main.tex.
mv "$MOD_DIR/main.pdf" "$MOD_DIR/manuscript.pdf"
cp output/figures/*.png "$MOD_DIR/figures/"

# Zip both for one-click upload.
( cd output-tex && rm -f submission_flat.zip submission_modular.zip \
    && zip -q -r submission_flat.zip     submission_flat \
    && zip -q -r submission_modular.zip  submission_modular )

# 4d. Editorial-Manager-ready packages.
#
#   submission_em_pdflatex/   Plan A -- flat dir, pdflatex-compatible
#                             wrapper (no fontspec), figures adjacent
#                             (no subdir), local elsarticle.cls +
#                             pre-built .bbl bundled so EM's TeX Live
#                             2022 environment reliably builds it.
#
#   submission_em_pdf_first/  Plan B (Elsevier-recommended when PDF
#                             upload is allowed): the fresh XeLaTeX
#                             PDFs uploaded as "Manuscript", source
#                             archive zipped for "LaTeX source files".
#                             Bypasses EM's server-side compile.
echo "=== assembling EM (Editorial Manager) packages ==="
EM_A_DIR="output-tex/submission_em_pdflatex"
EM_B_DIR="output-tex/submission_em_pdf_first"
rm -rf "$EM_A_DIR" "$EM_B_DIR"
mkdir -p "$EM_A_DIR" "$EM_B_DIR"

# --- Plan A: pdflatex-compatible flat source ---
# Flatten wrapper + body into ONE file per document. Uploading both a
# wrapper and its \input body separately to EM risks it picking the
# body as "primary" and hitting the classic nullfont cascade because
# the body has no \documentclass. One file per doc eliminates that.
python scripts/flatten_tex.py \
    output-tex/main_em.tex       manuscript_body_em \
    output-tex/manuscript_em_flat.tex
python scripts/flatten_tex.py \
    output-tex/supplement_em.tex supplement_body_em \
    output-tex/supplement_em_flat.tex
cp output-tex/manuscript_em_flat.tex     "$EM_A_DIR/manuscript.tex"
cp output-tex/supplement_em_flat.tex     "$EM_A_DIR/supplement.tex"
cp output-tex/references.bib             "$EM_A_DIR/"
cp output-tex/main.bbl                   "$EM_A_DIR/manuscript.bbl"
cp output-tex/supplement.bbl             "$EM_A_DIR/supplement.bbl"
# Bundle our elsarticle.cls to sidestep EM's aries.cls-shim issue.
if [[ -n "$(kpsewhich elsarticle.cls 2>/dev/null || true)" ]]; then
    cp "$(kpsewhich elsarticle.cls)"     "$EM_A_DIR/"
fi
cp output/figures/*.png                  "$EM_A_DIR/"
cp templates/UPLOAD_INSTRUCTIONS_EM.txt  "$EM_A_DIR/UPLOAD_INSTRUCTIONS.txt"
# Ship the freshly-built XeLaTeX PDFs for reviewer reference.
cp output-tex/main.pdf                   "$EM_A_DIR/manuscript_reference.pdf"
cp output-tex/supplement.pdf             "$EM_A_DIR/supplement_reference.pdf"

# --- Plan B: PDF-first, source zipped ---
cp output-tex/main.pdf                   "$EM_B_DIR/manuscript.pdf"
cp output-tex/supplement.pdf             "$EM_B_DIR/supplement.pdf"
cp output-tex/cover_letter.pdf           "$EM_B_DIR/cover_letter.pdf"
cp output-tex/response_to_reviewers.pdf  "$EM_B_DIR/response_to_reviewers.pdf"
# Source bundle for the "LaTeX source files" upload slot.
( cd "$EM_A_DIR" && zip -q -r "$PROJECT_ROOT/$EM_B_DIR/manuscript_source.zip" . )

# Zip Plan A for one-click upload as "Manuscript" set.
( cd output-tex && rm -f submission_em_pdflatex.zip submission_em_pdf_first.zip \
    && zip -q -r submission_em_pdflatex.zip   submission_em_pdflatex \
    && zip -q -r submission_em_pdf_first.zip  submission_em_pdf_first )

# 4e. Ship an LZW-compressed TIFF set of every figure for the
#     camera-ready / acceptance stage (Elsevier requires TIFF or EPS
#     for the production build). Generated on-demand by ImageMagick
#     from output/figures/*.png so the raw R-produced PNGs remain the
#     canonical source. Requires `magick` (ImageMagick 7).
if command -v magick >/dev/null 2>&1; then
    echo "=== converting figures to TIFF (Elsevier production format) ==="
    mkdir -p output/figures_tiff
    for png in output/figures/*.png; do
        name=$(basename "$png" .png)
        tif="output/figures_tiff/${name}.tif"
        if [[ ! -f "$tif" ]] || [[ "$png" -nt "$tif" ]]; then
            magick "$png" -density 300 -units PixelsPerInch -compress LZW "$tif"
        fi
    done
    ( cd output && rm -f figures_tiff.zip && zip -q -r figures_tiff.zip figures_tiff )
    cp output/figures_tiff.zip output-tex/figures_tiff.zip
else
    echo "  note: skipping TIFF conversion (ImageMagick 'magick' not on PATH)"
fi

# 5. Report. Exit nonzero if any expected PDF is missing or empty -- `latexmk -f`
#    above can swallow real failures, so the artifact check is the gate.
echo
echo "=== build complete ==="
missing=0
for pdf in output-tex/main.pdf output-tex/supplement.pdf output-tex/cover_letter.pdf output-tex/response_to_reviewers.pdf; do
    if [[ -s "$pdf" ]]; then
        size=$(stat -f '%z' "$pdf")
        echo "  $pdf  ($size bytes)"
    else
        echo "  $pdf  MISSING or empty"
        missing=$((missing + 1))
    fi
done
if (( missing > 0 )); then
    echo
    echo "FAIL: $missing PDF(s) missing or empty -- check latexmk logs in output-tex/" >&2
    exit 1
fi
