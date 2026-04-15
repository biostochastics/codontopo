#!/usr/bin/env bash
# =============================================================
# TopoSafe — Raw data download script
#
# Downloads all supplementary data files needed for the
# cross-study recoding reanalysis from their original sources.
#
# Run from the repo root (codon-topo/):
#   bash toposafe/scripts/download_data.sh
#
# Outputs to: data/codonsafe/<study>/
# Total size: ~103 MB
# =============================================================

set -euo pipefail

DATADIR="data/codonsafe"

echo "=== TopoSafe data download ==="
echo "Target: $DATADIR"
echo ""

# ── Helper ────────────────────────────────────────────────────
download() {
    local url="$1"
    local dest="$2"
    mkdir -p "$(dirname "$dest")"
    if [ -f "$dest" ]; then
        echo "  [skip] $(basename "$dest") (exists)"
    else
        echo "  [get]  $(basename "$dest")"
        curl -sL --retry 3 -o "$dest" "$url"
    fi
}

# ── 1. Napolitano et al. 2016 (PNAS) ──────────────────────────
echo "1/8  Napolitano 2016 (PNAS: AGR recoding + CRAM)"
BASE="https://www.pnas.org/doi/suppl/10.1073/pnas.1605856113/suppl_file"
for n in 01 02 03 05 06 07 08; do
    download "${BASE}/pnas.1605856113.sd${n}.xlsx" \
        "${DATADIR}/napolitano2016/pnas.1605856113.sd${n}.xlsx"
done
# sd04 is .xls
download "${BASE}/pnas.1605856113.sd04.xls" \
    "${DATADIR}/napolitano2016/pnas.1605856113.sd04.xls"

# ── 2. Fredens et al. 2019 (Nature) ───────────────────────────
echo "2/8  Fredens 2019 (Nature: Syn61)"
echo "  NOTE: Fredens supplementary files require institutional access."
echo "  DOI: 10.1038/s41586-019-1192-5"
echo "  Download manually to: ${DATADIR}/fredens2019/"

# ── 3. Ostrov et al. 2016 (Science) ───────────────────────────
echo "3/8  Ostrov 2016 (Science: 57-codon design)"
BASE="https://www.science.org/doi/suppl/10.1126/science.aaf3639/suppl_file"
for f in table_s3.xlsx table_s4.xlsx table_s5.xlsx; do
    download "${BASE}/${f}" "${DATADIR}/ostrov2016/${f}"
done
echo "  NOTE: Ostrov GenBank (aaf3639_Ostrov_rE.coli-57.Design.gb, 14 MB)"
echo "  requires institutional access. DOI: 10.1126/science.aaf3639"

# ── 4. Robertson et al. 2025 (Science, Syn57) ─────────────────
echo "4/8  Robertson 2025 Syn57 (Science)"
echo "  NOTE: Robertson Syn57 data files require institutional access."
echo "  DOI: 10.1126/science.ady4368"
echo "  Download Data S1-S10 manually to: ${DATADIR}/robertson2025_syn57/"

# ── 5. Grome/Robertson et al. 2025 (Nature, Ochre) ────────────
echo "5/8  Grome/Robertson 2025 Ochre (Nature)"
echo "  NOTE: Ochre supplementary files require institutional access."
echo "  DOI: 10.1038/s41586-024-08501-x"
echo "  Download MOESM4-16 manually to: ${DATADIR}/robertson2025/"

# ── 6. Frumkin et al. 2018 (PNAS) ─────────────────────────────
echo "6/8  Frumkin 2018 (PNAS: CGG recoding)"
BASE="https://www.pnas.org/doi/suppl/10.1073/pnas.1719375115/suppl_file"
for n in 01 02; do
    download "${BASE}/pnas.1719375115.sd${n}.xlsx" \
        "${DATADIR}/frumkin2018/pnas.1719375115.sd${n}.xlsx"
done

# ── 7. Lajoie et al. 2013 (Science) ───────────────────────────
echo "7/8  Lajoie 2013 (Science: C321.ΔA)"
echo "  NOTE: Lajoie tables available via PMC (PMCID: PMC4924538)"
echo "  Manual download to: ${DATADIR}/lajoie2013/"

# ── 8. Ding et al. 2024 (Science) ─────────────────────────────
echo "8/8  Ding 2024 (Science: Mammalian recoding)"
echo "  NOTE: Ding supplementary files require institutional access."
echo "  DOI: 10.1126/science.adm8143"
echo "  Download Data S1-S3 manually to: ${DATADIR}/ding2024/"

echo ""
echo "=== Download complete ==="
echo "Files requiring manual download: see notes above."
echo "Full provenance: ${DATADIR}/DATA_MANIFEST.md"
