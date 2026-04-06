# CODON-TOPO Preliminary Verification Scripts

## Origin
These scripts were written during an initial independent verification session
of claims made in Clayworth TN-2026-11 ("Induced Geometry, Holomorphic
Embedding, and the Filtration Resolution of the Fano-Hexagram-Codon Conjecture").

## Scripts (in order of execution)

1. **verify_filtration.py** — Core verification of Claims 1-3:
   - Two-fold degeneracy (bit-5 filtration)
   - Four-fold degeneracy (prefix uniformity)
   - Persistent homology (Serine anomaly)

2. **verify_embedding.py** — Holomorphic embedding φ: GF(2)⁶ → ℂ³:
   - Base → fourth-root-of-unity mapping
   - Geometric properties of synonymous codon clusters
   - Serine block separation in ℂ³

3. **verify_kras.py** — KRAS Fano-line analysis:
   - GGU ⊕ GUU ⊕ CAC = 0 verification
   - All single-bit mutations from GGU with Fano partners
   - Hamming distances for known KRAS G12 mutations

4. **verify_mitochondrial.py** — Cross-code validation (6 codes):
   - Standard, vertebrate mito, invertebrate mito, yeast mito,
     echinoderm mito, ciliate nuclear
   - Serine persistence test
   - Threonine disconnection discovery (yeast mito)

5. **all_codes.py** — Comprehensive analysis (23 NCBI tables):
   - All known variant genetic codes
   - Complete disconnection catalogue
   - Filtration breakage analysis

6. **summary_table.py** — Interpretive summary with evolutionary context

## Key Results
- All combinatorial claims VERIFIED
- Serine disconnected (ε=4) in ALL 23 codes tested
- Novel disconnections found: Thr(ε=2), Leu(ε=2), Ala(ε=3)
- Two-fold filtration: 100% invariant, zero exceptions
- Four-fold filtration: breaks only on stop codon reassignments

## Dependencies
- Python 3.8+
- No external packages required (stdlib only)

## Status
Preliminary/exploratory. To be refactored into the codon_topo package
as specified in CODON_TOPO_PRD_v1.docx.
