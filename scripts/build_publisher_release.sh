#!/usr/bin/env bash
# Canonical publisher-release build.
#
# Reproduces every displayed number in the manuscript and supplement,
# regenerates every figure and table CSV, refreshes manuscript_stats.json,
# and compiles the Typst PDFs. Intended to be run from a clean clone at
# tag v0.6.1 (or head-of-main during development).
#
# Requires: Python 3.11 with the dev extra installed (pip install -e ".[dev]"),
#           R 4.4 with ggplot2 + ggpubr + jsonlite, and Typst 0.14.2.
#
# Usage: bash scripts/build_publisher_release.sh
#
# Every venue (Elsevier response letter, CHANGELOG, README, manuscript
# Data-and-code availability, supplement §S24) quotes this single
# recipe.

set -euo pipefail

# Repo root — the script assumes it is invoked from there.
cd "$(dirname "$0")/.."

echo "[1/8] codon-topo all: running analyses (seed 135325, n = 10000)…"
codon-topo all --output-dir=./output --seed=135325 --n=10000

echo "[2/8] patch_pt_keys.py: per-table + trna + condlogit key patches…"
python3.11 scripts/patch_pt_keys.py

echo "[3/8] emit_trna_provenance_table.py: 24-row provenance table…"
python3.11 scripts/emit_trna_provenance_table.py

echo "[4/8] generate_tables.py: T-CSV table exports…"
python3.11 scripts/generate_tables.py

echo "[5/8] all_figures.R: regenerate main + supplement figures (PNG + PDF + TIFF)…"
Rscript src/codon_topo/visualization/R/all_figures.R output output/figures

echo "[6/8] strengthened_figures.R: regenerate strengthened panels…"
Rscript src/codon_topo/visualization/R/strengthened_figures.R output/tables output/figures

echo "[7/8] typst compile output/manuscript.typ…"
typst compile output/manuscript.typ output/manuscript.pdf

echo "[8/8] typst compile output/supplement.typ…"
typst compile output/supplement.typ output/supplement.pdf

echo "done. Deliverables at output/manuscript.pdf and output/supplement.pdf."
