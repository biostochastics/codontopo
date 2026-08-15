#!/usr/bin/env bash
# =============================================================
# TopoSafe — Full analysis pipeline
#
# Runs the complete cross-study recoding reanalysis and
# generates all figures and output tables.
#
# Prerequisites:
#   1. Install codon-topo: pip install -e ".[codonsafe]"
#   2. Download raw data: bash toposafe/scripts/download_data.sh
#   3. R + ggplot2 + ggpubr for figures
#
# Run from repo root (codon-topo/):
#   bash toposafe/scripts/run.sh
# =============================================================

set -euo pipefail

echo "=== TopoSafe cross-study reanalysis ==="

# Step 1: Run Python analysis pipeline
echo "[1/2] Running CodonSafe Python analysis..."
codon-topo codonsafe --output-dir=output/codonsafe

# Step 2: Generate figures with R
echo "[2/2] Generating figures (R)..."
Rscript src/codon_topo/visualization/R/codonsafe_figures.R

echo ""
echo "=== Done ==="
echo "Outputs:"
echo "  Figures:  output/codonsafe/"
echo "  Stats:    output/codonsafe/analysis_stats.json"
echo "  Tables:   output/codonsafe/*.csv"
