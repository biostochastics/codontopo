"""
Quick check: KRAS G12V mutation in the geometric framework.
Also checking the paper's claim about His codon CAC being on a Fano line.
"""

# The KRAS codon 12 is normally GGU (Glycine)
# G12V mutation: GGU → GUU (Glycine → Valine)

# In binary:
# GGU = (1,1,1,1,0,1) → Gly
# GUU = (1,1,0,1,0,1) → Val
# These differ at bit positions 2,3 (0-indexed): bits 2 and 3 flip from (1,1) to (0,1)
# That's the second base: G→U

BASE_TO_BITS = {'C': (0, 0), 'U': (0, 1), 'A': (1, 0), 'G': (1, 1)}

def codon_to_bits(codon):
    bits = []
    for base in codon:
        bits.extend(BASE_TO_BITS[base])
    return tuple(bits)

def hamming_distance(a, b):
    return sum(x != y for x, y in zip(a, b))

print("=" * 70)
print("KRAS G12V GEOMETRY CHECK")
print("=" * 70)

ggu = codon_to_bits('GGU')
guu = codon_to_bits('GUU')
print(f"\nKRAS WT:  GGU (Gly) = {ggu}")
print(f"KRAS G12V: GUU (Val) = {guu}")
print(f"Hamming distance: {hamming_distance(ggu, guu)}")

differing = [i for i in range(6) if ggu[i] != guu[i]]
print(f"Differing bit positions: {differing}")
print(f"This is a Hamming-distance-1 mutation (single bit flip)")

# The paper claims this is a "traversal" in the binary hypercube
# and that CAC (His) closes a Fano line with GGU and GUU

cac = codon_to_bits('CAC')
print(f"\nCAC (His) = {cac}")
print(f"d(GGU, CAC) = {hamming_distance(ggu, cac)}")
print(f"d(GUU, CAC) = {hamming_distance(guu, cac)}")

# Check if GGU ⊕ GUU ⊕ CAC = 0 in GF(2)^6 (Fano line condition)
xor_result = tuple((ggu[i] + guu[i] + cac[i]) % 2 for i in range(6))
print(f"\nGGU ⊕ GUU ⊕ CAC = {xor_result}")
is_fano_line = all(x == 0 for x in xor_result)
print(f"Is Fano line (XOR = 0)? {is_fano_line}")

if is_fano_line:
    print("VERIFIED ✓ — GGU, GUU, CAC form a Fano line in GF(2)^6")
else:
    print("FALSIFIED ✗")

# What does this mean biologically?
# If GGU→GUU (G12V) is a single-bit mutation, then the "third point"
# on the Fano line is CAC (His). The paper's prediction seems to be
# that His mutations should be geometrically linked to G12V.

# Let's also check: what OTHER single-bit mutations from GGU are possible?
print("\n" + "=" * 70)
print("ALL SINGLE-BIT MUTATIONS FROM GGU (Gly)")
print("=" * 70)

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

BITS_TO_CODON = {}
for codon in GENETIC_CODE:
    BITS_TO_CODON[codon_to_bits(codon)] = codon

for bit_pos in range(6):
    mutant_bits = list(ggu)
    mutant_bits[bit_pos] = 1 - mutant_bits[bit_pos]
    mutant_bits = tuple(mutant_bits)
    mutant_codon = BITS_TO_CODON[mutant_bits]
    mutant_aa = GENETIC_CODE[mutant_codon]
    
    # Find the third point on the Fano line
    third_bits = tuple((ggu[i] + mutant_bits[i]) % 2 for i in range(6))
    # Actually for a Fano line through GGU and mutant, the third point
    # is GGU ⊕ mutant (in GF(2)^6, a+b+c=0 means c=a⊕b)
    third_codon = BITS_TO_CODON.get(third_bits, "N/A")
    third_aa = GENETIC_CODE.get(third_codon, "N/A")
    
    change = "synonymous" if mutant_aa == 'Gly' else f"→ {mutant_aa}"
    print(f"  bit {bit_pos}: GGU→{mutant_codon} ({change}), "
          f"Fano third point: {third_codon} ({third_aa})")

# ============================================================
# Now the critical question: what is the KRAS prediction exactly?
# ============================================================
print("\n" + "=" * 70)
print("THE KRAS G12V PREDICTION ANALYSIS")
print("=" * 70)

print("""
The paper claims: "His mutations should be over-represented in 
KRAS G12V tumours."

The geometric basis: GGU (Gly12, WT) and GUU (Val, G12V) form a 
Fano line with CAC (His). 

But this needs unpacking:
- G12V means codon 12 mutates from GGT→GTT (DNA) / GGU→GUU (RNA)
- The "His" on the Fano line is the codon CAC
- But KRAS codon 12 doesn't normally encode His
- So the prediction must be about His mutations ELSEWHERE in the 
  tumor genome, not at KRAS codon 12 specifically

This is where the paper's claim gets vague. "His mutations should 
be over-represented" could mean:
1. Co-occurring mutations at His codons elsewhere in the genome
2. Mutations at the His codon in KRAS itself (but codon 12 is Gly, not His)
3. Something about the Fano-plane constraint propagation

Without seeing the full liquid biopsy section, the prediction is 
under-specified.
""")

# Check: what are the actual known KRAS G12 mutations?
print("Known KRAS codon 12 mutations (for reference):")
print("  GGU → GUU: G12V (Val) — very common")
print("  GGU → GAU: G12D (Asp) — very common") 
print("  GGU → GCU: G12A (Ala) — moderate")
print("  GGU → AGU: G12S (Ser) — less common")
print("  GGU → UGU: G12C (Cys) — common in NSCLC")
print("  GGU → GGC: G12G (Gly) — synonymous, not oncogenic")

print("\nIn binary, checking which of these are single-bit mutations:")
kras_mutations = {'G12V': 'GUU', 'G12D': 'GAU', 'G12A': 'GCU', 
                  'G12S': 'AGU', 'G12C': 'UGU'}
for name, codon in kras_mutations.items():
    bits = codon_to_bits(codon)
    hd = hamming_distance(ggu, bits)
    print(f"  {name}: GGU→{codon} = {ggu}→{bits}, Hamming distance = {hd}")

