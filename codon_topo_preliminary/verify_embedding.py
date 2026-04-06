"""
Verification of the holomorphic embedding φ: GF(2)^6 → ℂ^3
and the claim about synonymous codon geometric clustering.
"""
import cmath
import math
from collections import defaultdict

# Standard genetic code (excluding stops)
GENETIC_CODE = {
    'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
    'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
    'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'Stop', 'UAG': 'Stop',
    'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'Stop', 'UGG': 'Trp',
    'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
    'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
    'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
    'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met',
    'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
    'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
    'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
    'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
    'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
    'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly',
}

# Paper's encoding: base → 2 bits → fourth root of unity
# C=(0,0) → i^0 = 1+0i
# U=(0,1) → i^1 = 0+i  
# A=(1,0) → i^2 = -1+0i
# G=(1,1) → i^3 = 0-i
BASE_TO_COMPLEX = {
    'C': complex(1, 0),    # i^0 = 1
    'U': complex(0, 1),    # i^1 = i
    'A': complex(-1, 0),   # i^2 = -1
    'G': complex(0, -1),   # i^3 = -i
}

# Verify against paper's table (page 4)
print("=" * 70)
print("VERIFICATION: Holomorphic Embedding φ: GF(2)^6 → ℂ³")
print("=" * 70)

print("\nBase encoding verification (vs paper's table on p.4):")
print(f"  C → {BASE_TO_COMPLEX['C']}  (paper: 1+0i)  ✓")
print(f"  U → {BASE_TO_COMPLEX['U']}  (paper: 0+i)   ✓")
print(f"  A → {BASE_TO_COMPLEX['A']} (paper: -1+0i)  ✓")
print(f"  G → {BASE_TO_COMPLEX['G']} (paper: 0-i)    ✓")

def embed_codon(codon):
    """Map codon to ℂ³ via fourth roots of unity"""
    return (BASE_TO_COMPLEX[codon[0]], 
            BASE_TO_COMPLEX[codon[1]], 
            BASE_TO_COMPLEX[codon[2]])

# Verify paper's sample values
print("\nSample codon embeddings (verify against paper p.4):")
test_codons = ['GGU', 'GUU', 'CAC', 'AUG']
test_labels = ['Gly, KRAS WT', 'Val, KRAS G12V', 'His, Fano partner', 'Met, start codon']
for codon, label in zip(test_codons, test_labels):
    z = embed_codon(codon)
    z_str = f"({z[0].real:+.0f}{z[0].imag:+.0f}i, {z[1].real:+.0f}{z[1].imag:+.0f}i, {z[2].real:+.0f}{z[2].imag:+.0f}i)"
    print(f"  {codon} ({label}): φ = {z_str}")

# ============================================================
# PROPERTY 1: 2-fold - codon images agree in coords 1,2; differ in coord 3
# ============================================================
print("\n" + "=" * 70)
print("PROPERTY 1 (2-fold): Synonymous codons agree in ℂ-coords 1,2;")
print("           differ only in coord 3 by a fixed complex number (1-i)")
print("=" * 70)

aa_groups = defaultdict(list)
for codon, aa in GENETIC_CODE.items():
    if aa != 'Stop':
        aa_groups[aa].append(codon)

twofold_aas = [aa for aa, codons in aa_groups.items() if len(codons) == 2]
prop1_pass = True

for aa in sorted(twofold_aas):
    codons = sorted(aa_groups[aa])
    z1 = embed_codon(codons[0])
    z2 = embed_codon(codons[1])
    
    coord1_same = (z1[0] == z2[0])
    coord2_same = (z1[1] == z2[1])
    coord3_diff = z1[2] - z2[2]
    
    ok = coord1_same and coord2_same
    if not ok:
        prop1_pass = False
    
    status = "✓" if ok else "✗"
    print(f"  {aa:4s}: {codons[0]}→{z1[2]:.1f}  {codons[1]}→{z2[2]:.1f}  "
          f"coords1,2 same={coord1_same and coord2_same}  "
          f"coord3 diff={coord3_diff:.2f} {status}")

print(f"\nPROPERTY 1 RESULT: {'VERIFIED ✓' if prop1_pass else 'FALSIFIED ✗'}")

# ============================================================
# PROPERTY 2: 4-fold - identical first two ℂ coords, third cycles through all 4 roots
# ============================================================
print("\n" + "=" * 70)
print("PROPERTY 2 (4-fold): Identical first two ℂ-coords;")
print("           third coord cycles through all 4 fourth roots of unity")
print("=" * 70)

fourfold_aas = [aa for aa, codons in aa_groups.items() if len(codons) == 4]
fourth_roots = {complex(1,0), complex(0,1), complex(-1,0), complex(0,-1)}
prop2_pass = True

for aa in sorted(fourfold_aas):
    codons = sorted(aa_groups[aa])
    embeddings = [embed_codon(c) for c in codons]
    
    first_coords = set((z[0], z[1]) for z in embeddings)
    third_coords = set(z[2] for z in embeddings)
    
    first_uniform = len(first_coords) == 1
    third_exhausts = third_coords == fourth_roots
    
    ok = first_uniform and third_exhausts
    if not ok:
        prop2_pass = False
    
    status = "✓" if ok else "✗"
    print(f"  {aa:4s}: first 2 coords unique={first_uniform}  "
          f"third exhausts roots={third_exhausts} {status}")

print(f"\nPROPERTY 2 RESULT: {'VERIFIED ✓' if prop2_pass else 'FALSIFIED ✗'}")

# ============================================================
# PROPERTY 3: Serine clusters geometrically separated in ℂ³
# ============================================================
print("\n" + "=" * 70)
print("PROPERTY 3 (Serine): UCN and AGY map to geometrically separated")
print("           regions of ℂ³, no shared coordinate in any dimension")
print("=" * 70)

ser_codons = sorted(aa_groups['Ser'])
ser_embeddings = {c: embed_codon(c) for c in ser_codons}

ucn = {c: z for c, z in ser_embeddings.items() if c.startswith('UC')}
agy = {c: z for c, z in ser_embeddings.items() if c.startswith('AG')}

print("\nUCN block:")
for c, z in sorted(ucn.items()):
    print(f"  {c} → ({z[0]}, {z[1]}, {z[2]})")

print("\nAGY block:")
for c, z in sorted(agy.items()):
    print(f"  {c} → ({z[0]}, {z[1]}, {z[2]})")

# Check coordinate overlap per dimension
for dim in range(3):
    ucn_vals = set(z[dim] for z in ucn.values())
    agy_vals = set(z[dim] for z in agy.values())
    overlap = ucn_vals & agy_vals
    print(f"\n  Dimension {dim+1}: UCN values={ucn_vals}, AGY values={agy_vals}")
    print(f"  Overlap: {overlap if overlap else 'NONE'}")

# Compute Euclidean distances in ℂ³ between blocks
print("\nInter-block distances in ℂ³:")
min_dist = float('inf')
for c1, z1 in ucn.items():
    for c2, z2 in agy.items():
        dist = math.sqrt(sum(abs(z1[k] - z2[k])**2 for k in range(3)))
        print(f"  d({c1}, {c2}) = {dist:.4f}")
        min_dist = min(min_dist, dist)

print(f"\n  Minimum inter-block distance: {min_dist:.4f}")

# Compare to Arg and Leu (also 6-fold, but connected)
print("\n--- COMPARISON: Arg and Leu block structure in ℂ³ ---")
for aa in ['Arg', 'Leu']:
    codons = sorted(aa_groups[aa])
    # Split into 4-codon and 2-codon blocks
    # Find the split
    prefixes = defaultdict(list)
    for c in codons:
        prefixes[c[:2]].append(c)
    
    print(f"\n{aa} blocks:")
    blocks = list(prefixes.values())
    for block_codons in blocks:
        embs = [embed_codon(c) for c in block_codons]
        print(f"  {block_codons}:")
        for c, z in zip(block_codons, embs):
            print(f"    {c} → ({z[0]}, {z[1]}, {z[2]})")
    
    # Check if blocks share any coordinate
    if len(blocks) == 2:
        b1_embs = [embed_codon(c) for c in blocks[0]]
        b2_embs = [embed_codon(c) for c in blocks[1]]
        for dim in range(3):
            v1 = set(z[dim] for z in b1_embs)
            v2 = set(z[dim] for z in b2_embs)
            overlap = v1 & v2
            shared = "SHARED" if overlap else "DISJOINT"
            print(f"  Dim {dim+1}: {shared} (overlap={overlap if overlap else '{}'})")

