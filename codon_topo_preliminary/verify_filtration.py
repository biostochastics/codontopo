"""
Verification of Clayworth TN-2026-11 Filtration Claims
=======================================================
Testing the combinatorial assertions about the genetic code
when encoded as 6-bit vectors in GF(2)^6.
"""

# Standard genetic code table
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

# Binary encoding as described in the paper
# C=00, U=01, A=10, G=11
BASE_TO_BITS = {'C': (0, 0), 'U': (0, 1), 'A': (1, 0), 'G': (1, 1)}

def codon_to_bits(codon):
    """Convert 3-letter codon to 6-bit tuple in GF(2)^6"""
    bits = []
    for base in codon:
        bits.extend(BASE_TO_BITS[base])
    return tuple(bits)

def hamming_distance(a, b):
    """Hamming distance between two bit tuples"""
    return sum(x != y for x, y in zip(a, b))

# Build the full mapping
print("=" * 70)
print("VERIFICATION OF TN-2026-11 FILTRATION CLAIMS")
print("=" * 70)

# Group codons by amino acid
from collections import defaultdict
aa_to_codons = defaultdict(list)
aa_to_bits = defaultdict(list)

for codon, aa in GENETIC_CODE.items():
    if aa == 'Stop':
        continue
    aa_to_codons[aa].append(codon)
    aa_to_bits[aa].append(codon_to_bits(codon))

# Classify degeneracy
print("\n--- AMINO ACID DEGENERACY CLASSES ---")
degeneracy_classes = defaultdict(list)
for aa, codons in sorted(aa_to_codons.items()):
    deg = len(codons)
    degeneracy_classes[deg].append(aa)
    
for deg in sorted(degeneracy_classes.keys()):
    aas = degeneracy_classes[deg]
    print(f"  {deg}-fold: {', '.join(sorted(aas))} ({len(aas)} amino acids)")

print(f"\n  Total amino acids (excl. Stop): {len(aa_to_codons)}")

# ============================================================
# CLAIM 1: Two-fold degeneracy - differ at exactly bit index 5
# ============================================================
print("\n" + "=" * 70)
print("CLAIM 1: All 9 two-fold degenerate AAs differ at exactly one")
print("         binary position (bit index 5, the last base-position parity bit)")
print("=" * 70)

twofold = degeneracy_classes[2]
print(f"\nTwo-fold degenerate AAs: {sorted(twofold)}")
print(f"Count: {len(twofold)} (paper claims 9)")

claim1_pass = True
for aa in sorted(twofold):
    bits_list = aa_to_bits[aa]
    assert len(bits_list) == 2
    b1, b2 = bits_list
    differing_positions = [i for i in range(6) if b1[i] != b2[i]]
    codons = sorted(aa_to_codons[aa])
    status = "✓" if differing_positions == [5] else "✗"
    if differing_positions != [5]:
        claim1_pass = False
    print(f"  {aa:4s}: {codons[0]}={b1} {codons[1]}={b2} | differ at positions: {differing_positions} {status}")

print(f"\nCLAIM 1 RESULT: {'VERIFIED ✓' if claim1_pass else 'FALSIFIED ✗'}")

# ============================================================
# CLAIM 2: Four-fold degeneracy - identical first-4-bit prefixes
# ============================================================
print("\n" + "=" * 70)
print("CLAIM 2: All 5 four-fold degenerate AAs share identical first-4-bit")
print("         prefixes, with last 2 bits exhausting GF(2)^2 = {00,01,10,11}")
print("=" * 70)

fourfold = degeneracy_classes[4]
print(f"\nFour-fold degenerate AAs: {sorted(fourfold)}")
print(f"Count: {len(fourfold)} (paper claims 5)")

claim2_pass = True
for aa in sorted(fourfold):
    bits_list = aa_to_bits[aa]
    prefixes = set(b[:4] for b in bits_list)
    suffixes = set(b[4:] for b in bits_list)
    expected_suffixes = {(0,0), (0,1), (1,0), (1,1)}
    
    prefix_ok = len(prefixes) == 1
    suffix_ok = suffixes == expected_suffixes
    status = "✓" if (prefix_ok and suffix_ok) else "✗"
    if not (prefix_ok and suffix_ok):
        claim2_pass = False
    
    codons = sorted(aa_to_codons[aa])
    print(f"  {aa:4s}: codons={codons}")
    print(f"         bits={[b for b in bits_list]}")
    print(f"         prefixes={prefixes} (unique={len(prefixes)}) | suffixes={suffixes}")
    print(f"         prefix_uniform={prefix_ok}, suffixes_exhaust_GF2²={suffix_ok} {status}")

print(f"\nCLAIM 2 RESULT: {'VERIFIED ✓' if claim2_pass else 'FALSIFIED ✗'}")

# ============================================================
# What about 3-fold (Ile) and 6-fold degeneracy?
# ============================================================
print("\n" + "=" * 70)
print("ADDITIONAL CHECK: Non-standard degeneracy classes")
print("=" * 70)

for deg in [1, 3, 6]:
    if deg in degeneracy_classes:
        print(f"\n  {deg}-fold degenerate: {sorted(degeneracy_classes[deg])}")
        for aa in sorted(degeneracy_classes[deg]):
            bits_list = aa_to_bits[aa]
            codons = sorted(aa_to_codons[aa])
            print(f"    {aa}: {codons}")
            for i, (c, b) in enumerate(zip(codons, bits_list)):
                print(f"      {c} = {b}")
            
            # Check structure
            if len(bits_list) > 1:
                # Compute pairwise Hamming distances
                for i in range(len(bits_list)):
                    for j in range(i+1, len(bits_list)):
                        hd = hamming_distance(bits_list[i], bits_list[j])
                        print(f"      d({codons[i]}, {codons[j]}) = {hd}")


# ============================================================
# CLAIM 3: Persistent homology - Serine unique disconnection
# ============================================================
print("\n" + "=" * 70)
print("CLAIM 3: Persistent homology at ε=1 (Hamming distance)")
print("         Serine should be the ONLY amino acid with disconnected")
print("         codon graph at ε=1, reconnecting at ε=4")
print("=" * 70)

def connected_components(bits_list, epsilon):
    """Find connected components at Hamming distance threshold epsilon"""
    n = len(bits_list)
    # Union-Find
    parent = list(range(n))
    
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    
    def union(x, y):
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py
    
    for i in range(n):
        for j in range(i+1, n):
            if hamming_distance(bits_list[i], bits_list[j]) <= epsilon:
                union(i, j)
    
    components = len(set(find(i) for i in range(n)))
    return components

print(f"\n{'Amino Acid':<12} {'Deg':<5} {'β₀(ε=1)':<10} {'β₀(ε=2)':<10} {'β₀(ε=3)':<10} {'β₀(ε=4)':<10} {'Connected at ε=1?'}")
print("-" * 75)

disconnected_at_1 = []
for aa in sorted(aa_to_codons.keys()):
    bits_list = aa_to_bits[aa]
    deg = len(bits_list)
    if deg == 1:
        # Single codon = trivially connected
        b0_1 = b0_2 = b0_3 = b0_4 = 1
    else:
        b0_1 = connected_components(bits_list, 1)
        b0_2 = connected_components(bits_list, 2)
        b0_3 = connected_components(bits_list, 3)
        b0_4 = connected_components(bits_list, 4)
    
    connected = "YES" if b0_1 == 1 else f"NO ({b0_1} components)"
    print(f"  {aa:<10} {deg:<5} {b0_1:<10} {b0_2:<10} {b0_3:<10} {b0_4:<10} {connected}")
    
    if b0_1 > 1:
        disconnected_at_1.append((aa, b0_1))

print(f"\nAmino acids disconnected at ε=1: {disconnected_at_1}")
if len(disconnected_at_1) == 1 and disconnected_at_1[0][0] == 'Ser':
    print("CLAIM 3 RESULT: VERIFIED ✓ — Serine is uniquely disconnected at ε=1")
else:
    print("CLAIM 3 RESULT: FALSIFIED ✗")

# Now check Serine specifically
print("\n--- SERINE DETAIL ---")
ser_codons = sorted(aa_to_codons['Ser'])
ser_bits = [codon_to_bits(c) for c in ser_codons]
print("Serine codons and their binary representations:")
for c, b in zip(ser_codons, ser_bits):
    print(f"  {c} = {b}")

# Identify the two clusters
print("\nPairwise Hamming distances:")
for i in range(len(ser_bits)):
    for j in range(i+1, len(ser_bits)):
        hd = hamming_distance(ser_bits[i], ser_bits[j])
        print(f"  d({ser_codons[i]}, {ser_codons[j]}) = {hd}")

# Check paper's claim: UCN block (prefix 0100) and AGY block (prefix 1011)
print("\nCluster identification:")
ucn = [(c, b) for c, b in zip(ser_codons, ser_bits) if c.startswith('UC')]
agy = [(c, b) for c, b in zip(ser_codons, ser_bits) if c.startswith('AG')]
print(f"  UCN block: {[c for c,b in ucn]}")
print(f"    Prefixes (first 4 bits): {set(b[:4] for c,b in ucn)}")
print(f"  AGY block: {[c for c,b in agy]}")
print(f"    Prefixes (first 4 bits): {set(b[:4] for c,b in agy)}")

# Minimum inter-block distance
min_inter = min(hamming_distance(b1, b2) for _, b1 in ucn for _, b2 in agy)
print(f"\n  Minimum inter-block Hamming distance: {min_inter}")
print(f"  Paper claims: 4")
print(f"  {'VERIFIED ✓' if min_inter == 4 else 'FALSIFIED ✗'}")

# ============================================================
# CHECK: 6-fold AAs structure (Arg, Leu, Ser)
# ============================================================
print("\n" + "=" * 70)
print("ADDITIONAL: 6-fold degenerate AA block structure")
print("=" * 70)

for aa in ['Arg', 'Leu', 'Ser']:
    codons = sorted(aa_to_codons[aa])
    bits = [codon_to_bits(c) for c in codons]
    
    # Try to identify two blocks based on first 2 bits
    blocks = defaultdict(list)
    for c, b in zip(codons, bits):
        blocks[b[:2]].append((c, b))
    
    print(f"\n{aa}: {codons}")
    print(f"  Blocks by first-2-bit prefix:")
    for prefix, members in sorted(blocks.items()):
        member_codons = [c for c, b in members]
        print(f"    prefix {prefix}: {member_codons} ({len(members)} codons)")
    
    # Check connectivity at ε=1
    comp = connected_components(bits, 1)
    print(f"  Connected components at ε=1: {comp}")

