"""
Final summary: the disconnection-depth catalogue with evolutionary dating.
"""

print("=" * 80)
print("DISCONNECTION CATALOGUE WITH EVOLUTIONARY CONTEXT")
print("=" * 80)

print("""
UNIVERSAL (present in ALL codes tested):
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  Serine    ε=4    d_min=4    23/23 codes    
  ──────────────────────────────────────────
  UCN block (4 codons) ←──── gap of 4 ────→ AGY block (2+ codons)
  
  Evolutionary date: Pre-LUCA (>3.5 Gya)
  Mechanism: Independent tRNA recruitment
  Evidence: Universal across all domains of life
  Framework prediction: Maximum reconnection distance = deepest event
  Status: CONSISTENT ✓


LINEAGE-SPECIFIC DISCONNECTIONS:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  1. Threonine (Yeast Mito, Table 3)
     ε=2    d_min=2    1/23 codes
     ──────────────────────────────────
     ACN block (4) ←── gap of 2 ──→ CUN block (4, ex-Leu)
     
     Evolutionary date: Within Saccharomycetales (~200-300 Mya)
     Mechanism: CUN reassignment from Leu → Thr via tRNA loss/gain
     Evidence: Restricted to budding yeast lineage
     Framework prediction: Small ε = shallow event
     Status: CONSISTENT ✓

  2. Leucine (Chlorophycean Mito, Table 16)
     ε=2    d_min=2    1/23 codes
     ──────────────────────────────────
     CUN+UUR block (6) ←── gap of 2 ──→ UAG singleton (1, ex-Stop)
     
     Evolutionary date: Within Chlorophyceae (green algae, ~500-700 Mya?)
     Mechanism: UAG stop codon reassigned to Leu
     Evidence: Restricted to certain green algal mitochondria
     Framework prediction: Small ε = shallow event
     Status: CONSISTENT ✓

  3. Leucine (Scenedesmus obliquus Mito, Table 22)
     ε=2    d_min=2    1/23 codes
     ──────────────────────────────────
     Same as Table 16 — UAG→Leu
     Species-specific variant of the Chlorophycean pattern
     Status: CONSISTENT ✓

  4. Alanine (Pachysolen tannophilus Nuclear, Table 26)
     ε=3    d_min=3    1/23 codes
     ──────────────────────────────────
     GCN block (4) ←── gap of 3 ──→ CUG singleton (1, ex-Leu)
     
     Evolutionary date: Within Saccharomycetales (~100-200 Mya?)
     Mechanism: CUG codon ambiguity/reassignment (related to the
                broader "CUG clade" phenomenon in yeasts)
     Evidence: Species-specific nuclear code change
     Framework prediction: Intermediate ε = intermediate depth?
     Status: INTERESTING — ε=3 between Ser(4) and Thr/Leu(2)
     Note: The CUG reassignment in yeasts is actively studied;
           Pachysolen assigns CUG→Ala, while Candida assigns CUG→Ser (Table 12)

  5. Serine (Alternative Yeast Nuclear, Table 12)
     ε=3    d_min=2    3 components!    1/23 codes
     ──────────────────────────────────
     AGY block (2) | CUG singleton (1) | UCN block (4)
     
     This is UNIQUE: Serine gains a THIRD disconnected block (CUG)
     Three components instead of two!
     Evolutionary date: Within Candida clade (~100-200 Mya)
     CUG→Ser is one of the most-studied codon reassignment events
     Status: HIGHLY INFORMATIVE — adds a third topology layer


PATTERN SUMMARY:
━━━━━━━━━━━━━━━

  ε_reconnect    Events                           Approximate age
  ────────────────────────────────────────────────────────────────
      4           Serine split (universal)          >3.5 Gya (pre-LUCA)
      3           Ala/CUG (Pachysolen)              ~100-200 Mya
                  Ser 3-component (Candida/Table12) ~100-200 Mya  
      2           Thr/CUN (yeast mito)              ~200-300 Mya
                  Leu/UAG (green algae mito)        ~500-700 Mya
""")

print("""
CRITICAL OBSERVATIONS:
━━━━━━━━━━━━━━━━━━━━━

1. TWO-FOLD FILTRATION: 100% preserved across ALL 23 codes.
   Not a single exception. This is the wobble floor — 
   it cannot be broken by codon reassignment because 
   reassignment operates at the block level, not within blocks.

2. FOUR-FOLD FILTRATION: Breaks ONLY when stop codons are 
   reassigned to an amino acid (Gln in ciliates, Glu in 
   peritrichs/Blastocrithidia). Stop codons sit in different 
   prefix blocks, so reassigning them violates prefix uniformity.
   This is EXACTLY the expected failure mode.

3. SERINE IS ALWAYS DISCONNECTED AT ε=4. No codon reassignment
   in any known code bridges the UCN-AGY gap. Even adding 
   AGA/AGG to Serine (invertebrate mito) doesn't help — 
   they join the AGY block, not bridge to UCN.

4. ALL non-Serine disconnections have ε ≤ 3. The framework 
   produces a clear hierarchy: ancient events → large ε, 
   recent events → small ε.

5. The ε=2 vs ε=3 distinction is less clean. Leu/UAG (ε=2)
   might be older than Ala/CUG (ε=3), despite having a smaller
   reconnection distance. This could mean:
   (a) The correlation is ordinal, not linear
   (b) The distance reflects structural placement, not just time
   (c) More data needed to calibrate
""")

