#!/usr/bin/env bash
# build_em_bundle.sh -- Typst -> LaTeX -> PDF pipeline for the Elsevier
# Editorial Manager (EM) submission bundle.
#
# End-to-end automated: freeze inline Typst stats, pandoc-convert typst to
# latex, post-process, flatten wrapper + body, pdflatex + bibtex to PDF,
# regenerate TIFFs for Elsevier production, assemble the shipped
# submission bundle at output/publisher_deliverable/submission_em_<tag>/
# and zip it as submission_bundle.zip.
#
# Wired into scripts/build_publisher_release.sh so submission_bundle.zip
# stays in sync with the Typst-compiled reference PDFs after every
# analysis rerun.
#
# Requirements:
#   pandoc 3.1+, TeX Live (pdflatex, bibtex, kpsewhich), python 3.11.
#   Optional: ImageMagick 'magick' for TIFF regeneration.
#
# Usage:
#   bash scripts/build_em_bundle.sh                # tag = v0.6.1 (default)
#   BUNDLE_TAG=v0.7.0 bash scripts/build_em_bundle.sh

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"

BUNDLE_TAG="${BUNDLE_TAG:-v0.6.1}"
BUNDLE_DIR="output/publisher_deliverable/submission_em_${BUNDLE_TAG}"

# --- prerequisites ---
need() {
    if ! command -v "$1" >/dev/null 2>&1; then
        echo "error: $1 not on PATH -- install $2" >&2
        exit 1
    fi
}
need pandoc "pandoc 3.1+"
need pdflatex "TeX Live"
need bibtex "TeX Live"
need python3.11 "python 3.11"

# --- 0. staging directory + inputs ---
mkdir -p output-tex
cp -f output/references.bib        output-tex/references.bib
cp -f output/manuscript_stats.json output-tex/manuscript_stats.json
cp -f templates/main_em.tex        output-tex/main_em.tex
cp -f templates/supplement_em.tex  output-tex/supplement_em.tex

# --- 1. freeze inline Typst stats ---
echo "[1/6] freeze inline Typst stats"
python3.11 scripts/freeze_typst.py \
    output/manuscript.typ output/manuscript_stats.json \
    output-tex/manuscript-frozen.typ > /dev/null || {
    echo "  manuscript freeze failed" >&2; exit 1;
}
python3.11 scripts/freeze_typst.py \
    output/supplement.typ output/manuscript_stats.json \
    output-tex/supplement-frozen.typ > /dev/null || {
    echo "  supplement freeze failed" >&2; exit 1;
}

# --- 2. pandoc typst -> latex ---
echo "[2/6] pandoc typst -> latex body"
(
    cd output-tex
    pandoc manuscript-frozen.typ --from=typst --to=latex --natbib \
        -o manuscript_body_raw.tex
    pandoc supplement-frozen.typ --from=typst --to=latex --natbib \
        -o supplement_body_raw.tex
)

# --- 3. post-process (Unicode -> macros, \citep{sec} -> \ref, ...) ---
echo "[3/6] post-process latex body"
python3.11 scripts/postprocess_tex.py \
    output-tex/manuscript_body_raw.tex output-tex/manuscript_body.tex
python3.11 scripts/postprocess_tex.py \
    output-tex/supplement_body_raw.tex output-tex/supplement_body.tex

# --- 4. adapt for the flat EM directory (strip figures/ prefix, keep
#        \includegraphics so the reference PDF renders actual figures;
#        Elsevier's typesetter substitutes TIFFs in production) ---
echo "[4/6] strip figures/ prefix + copy figure PNGs"
python3.11 << 'PY'
import re
from pathlib import Path

STRIP = re.compile(r'(\\includegraphics(?:\[[^\]]*\])?\{)figures/([^}]+\})')

for src in ('output-tex/manuscript_body.tex', 'output-tex/supplement_body.tex'):
    p = Path(src)
    text = p.read_text()
    text2, n = STRIP.subn(r'\1\2', text)
    p.write_text(text2)
    print(f"    {src}: {n} figures/-prefix strips")
PY

# The wrappers input `manuscript_body_em.tex` / `supplement_body_em.tex`.
# Reuse the plain body under those names so the wrapper doesn't need to
# change.
cp output-tex/manuscript_body.tex output-tex/manuscript_body_em.tex
cp output-tex/supplement_body.tex output-tex/supplement_body_em.tex

# --- 5. flatten wrapper + body -> single self-contained .tex ---
echo "[5/6] flatten wrapper + body"
python3.11 scripts/flatten_tex.py \
    output-tex/main_em.tex       manuscript_body_em \
    output-tex/manuscript.tex
python3.11 scripts/flatten_tex.py \
    output-tex/supplement_em.tex supplement_body_em \
    output-tex/supplement.tex

# Copy PNGs so \includegraphics can find them in the compile CWD.
cp -f output/figures/*.png output-tex/

# --- 6. pdflatex + bibtex + pdflatex + pdflatex (three-pass for refs) ---
echo "[6/6] compile manuscript + supplement (pdflatex + bibtex)"
_compile_doc() {
    local doc="$1"
    rm -f "${doc}".{aux,bbl,blg,log,out,pdf,toc}
    pdflatex -interaction=nonstopmode -halt-on-error "${doc}.tex" \
        > "${doc}_pass1.log" 2>&1 || {
        echo "  ${doc}.tex pdflatex pass 1 failed:" >&2
        tail -30 "${doc}_pass1.log" >&2
        exit 1
    }
    bibtex "${doc}" > "${doc}_bibtex.log" 2>&1 || true
    pdflatex -interaction=nonstopmode "${doc}.tex" > "${doc}_pass2.log" 2>&1
    pdflatex -interaction=nonstopmode "${doc}.tex" > "${doc}_pass3.log" 2>&1
    if [[ ! -s "${doc}.pdf" ]]; then
        echo "  ${doc}.pdf missing after compile" >&2
        tail -30 "${doc}_pass3.log" >&2
        exit 1
    fi
    echo "    ${doc}.pdf ($(stat -f '%z' "${doc}.pdf") bytes)"
}
(
    cd output-tex
    _compile_doc manuscript
    _compile_doc supplement
)

# --- 7. TIFF regeneration for Elsevier production (optional) ---
if command -v magick > /dev/null 2>&1; then
    echo "[+] regenerate TIFFs (LZW, 300 DPI) for Elsevier production"
    mkdir -p output/figures_tiff
    for png in output/figures/*.png; do
        name=$(basename "$png" .png)
        tif="output/figures_tiff/${name}.tif"
        if [[ ! -f "$tif" ]] || [[ "$png" -nt "$tif" ]]; then
            magick "$png" -density 300 -units PixelsPerInch -compress LZW "$tif"
        fi
    done
else
    echo "  note: skipping TIFF regeneration (ImageMagick 'magick' not on PATH)"
fi

# --- 8. assemble submission bundle ---
echo "[+] assemble submission bundle at $BUNDLE_DIR"
rm -rf "$BUNDLE_DIR"
mkdir -p "$BUNDLE_DIR"
cp output-tex/manuscript.tex "$BUNDLE_DIR/"
cp output-tex/supplement.tex "$BUNDLE_DIR/"
cp output-tex/manuscript.pdf "$BUNDLE_DIR/"
cp output-tex/supplement.pdf "$BUNDLE_DIR/"
[[ -f output-tex/manuscript.bbl ]] && cp output-tex/manuscript.bbl "$BUNDLE_DIR/"
[[ -f output-tex/supplement.bbl ]] && cp output-tex/supplement.bbl "$BUNDLE_DIR/"
cp output-tex/references.bib "$BUNDLE_DIR/"
# bundle elsarticle.cls to sidestep EM class-shim issues
if [[ -n "$(kpsewhich elsarticle.cls 2>/dev/null || true)" ]]; then
    cp "$(kpsewhich elsarticle.cls)" "$BUNDLE_DIR/"
fi
cp output/figures/*.png "$BUNDLE_DIR/"
if [[ -d output/figures_tiff ]]; then
    cp output/figures_tiff/*.tif "$BUNDLE_DIR/"
fi

# Zip.
(
    cd output/publisher_deliverable
    rm -f submission_bundle.zip
    zip -q -r submission_bundle.zip "submission_em_${BUNDLE_TAG}"
)

echo
echo "=== EM bundle ready ==="
echo "  manuscript: output-tex/manuscript.pdf"
echo "  supplement: output-tex/supplement.pdf"
echo "  bundle dir: $BUNDLE_DIR/"
echo "  bundle zip: output/publisher_deliverable/submission_bundle.zip"
