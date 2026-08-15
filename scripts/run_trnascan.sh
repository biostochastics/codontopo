#!/usr/bin/env bash
# Run tRNAscan-SE 2.0.12 on all verified genome assemblies.
# Requires: tRNAscan-SE 2.0.12 with Infernal 1.1.4 in PATH.
#
# Usage: bash scripts/run_trnascan.sh
#
# Downloads assemblies from NCBI and scans each with tRNAscan-SE
# in eukaryotic mode (-E). Results go to data/trnascan_results/.

set -euo pipefail

OUTDIR="data/trnascan_results"
ASMDIR="data/assemblies"
mkdir -p "$OUTDIR" "$ASMDIR"

# Assembly accessions and short names
declare -A ASSEMBLIES=(
    ["t_thermophila"]="GCF_000189635.1"
    ["p_tetraurelia"]="GCF_000165425.1"
    ["o_trifallax"]="GCA_000295675.1"
    ["s_coeruleus"]="GCA_001970955.1"
    ["i_multifiliis"]="GCF_000220395.1"
    ["b_nonstop_P57"]="GCA_028554745.1"
    ["f_salina"]="GCA_022984795.1"
    ["p_persalinus"]="GCA_001447515.1"
    ["h_grandinella"]="GCA_006369765.1"
    ["b_stoltei"]="GCA_965603825.1"
)

for name in "${!ASSEMBLIES[@]}"; do
    acc="${ASSEMBLIES[$name]}"
    zipfile="$ASMDIR/${name}.zip"
    fastadir="$ASMDIR/${name}_data/ncbi_dataset/data/${acc}"

    echo "=== Processing $name ($acc) ==="

    # Download if not already present
    if [ ! -d "$fastadir" ]; then
        echo "  Downloading $acc..."
        curl -sL \
            "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/${acc}/download?include_annotation_type=GENOME_FASTA" \
            -o "$zipfile"
        unzip -o "$zipfile" -d "$ASMDIR/${name}_data" > /dev/null
    fi

    # Find the FASTA file
    fasta=$(find "$fastadir" -name "*.fna" | head -1)
    if [ -z "$fasta" ]; then
        echo "  ERROR: No FASTA found for $name"
        continue
    fi

    # Run tRNAscan-SE in eukaryotic mode
    echo "  Running tRNAscan-SE..."
    tRNAscan-SE -E \
        -o "$OUTDIR/${name}.out" \
        --stats "$OUTDIR/${name}.stats" \
        "$fasta" 2>&1 | grep -E "^(Status|tRNAs)" || true

    echo "  Done: $OUTDIR/${name}.out"
    echo ""
done

echo "All scans complete. Results in $OUTDIR/"
