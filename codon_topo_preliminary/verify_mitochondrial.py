"""
Testing the Clayworth filtration framework against variant genetic codes.
If the framework captures real evolutionary topology, different codes
should produce different persistent homology signatures corresponding
to known reassignment events.
"""
from collections import defaultdict

BASE_TO_BITS = {'C': (0, 0), 'U': (0, 1), 'A': (1, 0), 'G': (1, 1)}

def codon_to_bits(codon):
    bits = []
    for base in codon:
        bits.extend(BASE_TO_BITS[base])
    return tuple(bits)

def hamming_distance(a, b):
    return sum(x != y for x, y in zip(a, b))

def connected_components(bits_list, epsilon):
    n = len(bits_list)
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
    return len(set(find(i) for i in range(n)))

def analyze_code(code_table, code_name):
    """Full filtration analysis of a genetic code."""
    print(f"\n{'=' * 70}")
    print(f"ANALYSIS: {code_name}")
    print(f"{'=' * 70}")
    
    # Group by amino acid (exclude stops)
    aa_groups = defaultdict(list)
    aa_bits = defaultdict(list)
    stop_codons = []
    
    for codon, aa in code_table.items():
        if aa == 'Stop':
            stop_codons.append(codon)
            continue
        aa_groups[aa].append(codon)
        aa_bits[aa].append(codon_to_bits(codon))
    
    # Degeneracy classes
    deg_classes = defaultdict(list)
    for aa, codons in aa_groups.items():
        deg_classes[len(codons)].append(aa)
    
    print(f"\nStop codons: {sorted(stop_codons)}")
    print(f"Amino acids: {len(aa_groups)}")
    print(f"\nDegeneracy distribution:")
    for deg in sorted(deg_classes.keys()):
        aas = sorted(deg_classes[deg])
        print(f"  {deg}-fold: {', '.join(aas)} ({len(aas)} AAs)")
    
    # === FILTRATION CHECK ===
    # Check 2-fold: do they all differ at exactly bit 5?
    twofold = deg_classes.get(2, [])
    if twofold:
        print(f"\n--- TWO-FOLD DEGENERACY CHECK (bit index 5) ---")
        tw_pass = 0
        tw_fail = 0
        for aa in sorted(twofold):
            bits_list = aa_bits[aa]
            b1, b2 = bits_list
            differing = [i for i in range(6) if b1[i] != b2[i]]
            ok = differing == [5]
            status = "✓" if ok else f"✗ (differs at {differing})"
            if ok:
                tw_pass += 1
            else:
                tw_fail += 1
            codons = sorted(aa_groups[aa])
            print(f"  {aa:4s}: {codons[0]}↔{codons[1]} differ at {differing} {status}")
        print(f"  Result: {tw_pass}/{tw_pass+tw_fail} pass standard filtration")
    
    # Check 4-fold: identical first-4-bit prefix, last 2 exhaust GF(2)^2?
    fourfold = deg_classes.get(4, [])
    if fourfold:
        print(f"\n--- FOUR-FOLD DEGENERACY CHECK (prefix uniformity) ---")
        ff_pass = 0
        ff_fail = 0
        for aa in sorted(fourfold):
            bits_list = aa_bits[aa]
            prefixes = set(b[:4] for b in bits_list)
            suffixes = set(b[4:] for b in bits_list)
            expected_suffixes = {(0,0), (0,1), (1,0), (1,1)}
            ok = len(prefixes) == 1 and suffixes == expected_suffixes
            status = "✓" if ok else "✗"
            if ok:
                ff_pass += 1
            else:
                ff_fail += 1
            print(f"  {aa:4s}: prefixes={prefixes} suffixes={suffixes} {status}")
        print(f"  Result: {ff_pass}/{ff_pass+ff_fail} pass standard filtration")
    
    # === PERSISTENT HOMOLOGY ===
    print(f"\n--- PERSISTENT HOMOLOGY (connected components) ---")
    print(f"  {'AA':<6} {'Deg':<5} {'β₀(ε=1)':<9} {'β₀(ε=2)':<9} {'β₀(ε=3)':<9} {'β₀(ε=4)':<9} {'Status'}")
    print(f"  {'-'*55}")
    
    disconnected = []
    for aa in sorted(aa_groups.keys()):
        bits_list = aa_bits[aa]
        deg = len(bits_list)
        if deg == 1:
            b0 = {1: 1, 2: 1, 3: 1, 4: 1}
        else:
            b0 = {}
            for eps in [1, 2, 3, 4]:
                b0[eps] = connected_components(bits_list, eps)
        
        status = "connected" if b0[1] == 1 else f"DISCONNECTED ({b0[1]} components)"
        print(f"  {aa:<6} {deg:<5} {b0[1]:<9} {b0[2]:<9} {b0[3]:<9} {b0[4]:<9} {status}")
        
        if b0[1] > 1:
            # Find the blocks
            codons = sorted(aa_groups[aa])
            bits = [codon_to_bits(c) for c in codons]
            
            # Cluster by connectivity at ε=1
            n = len(bits)
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
                    if hamming_distance(bits[i], bits[j]) <= 1:
                        union(i, j)
            
            clusters = defaultdict(list)
            for i in range(n):
                clusters[find(i)].append(codons[i])
            
            block_info = []
            for root, members in clusters.items():
                block_info.append(members)
            
            # Min inter-block distance
            min_inter = float('inf')
            for i in range(n):
                for j in range(i+1, n):
                    if find(i) != find(j):
                        d = hamming_distance(bits[i], bits[j])
                        min_inter = min(min_inter, d)
            
            reconnect_eps = None
            for eps in range(1, 7):
                if connected_components(bits, eps) == 1:
                    reconnect_eps = eps
                    break
            
            disconnected.append({
                'aa': aa,
                'blocks': block_info,
                'min_inter_distance': min_inter,
                'reconnect_at': reconnect_eps
            })
    
    if disconnected:
        print(f"\n--- DISCONNECTED AMINO ACIDS (DETAIL) ---")
        for d in disconnected:
            print(f"\n  {d['aa']}:")
            for i, block in enumerate(d['blocks']):
                print(f"    Block {i+1}: {block}")
            print(f"    Min inter-block Hamming distance: {d['min_inter_distance']}")
            print(f"    Reconnects at ε={d['reconnect_at']}")
    else:
        print(f"\n  No disconnected amino acids found!")
    
    return {
        'name': code_name,
        'n_amino_acids': len(aa_groups),
        'degeneracy': {k: sorted(v) for k, v in deg_classes.items()},
        'disconnected': disconnected
    }


# ================================================================
# GENETIC CODE DEFINITIONS
# ================================================================

# Standard code (NCBI translation table 1)
STANDARD = {
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

# Vertebrate Mitochondrial Code (NCBI translation table 2)
# Key differences from standard:
# - UGA: Stop → Trp
# - AUA: Ile → Met  
# - AGA: Arg → Stop
# - AGG: Arg → Stop
VERT_MITO = dict(STANDARD)
VERT_MITO['UGA'] = 'Trp'
VERT_MITO['AUA'] = 'Met'
VERT_MITO['AGA'] = 'Stop'
VERT_MITO['AGG'] = 'Stop'

# Invertebrate Mitochondrial Code (NCBI translation table 5)
# Key differences from standard:
# - UGA: Stop → Trp
# - AUA: Ile → Met
# - AGA: Arg → Ser
# - AGG: Arg → Ser
INVERT_MITO = dict(STANDARD)
INVERT_MITO['UGA'] = 'Trp'
INVERT_MITO['AUA'] = 'Met'
INVERT_MITO['AGA'] = 'Ser'
INVERT_MITO['AGG'] = 'Ser'

# Yeast Mitochondrial Code (NCBI translation table 3)
# - UGA: Stop → Trp
# - AUA: Ile → Met
# - CUU,CUC,CUA,CUG: Leu → Thr
YEAST_MITO = dict(STANDARD)
YEAST_MITO['UGA'] = 'Trp'
YEAST_MITO['AUA'] = 'Met'
YEAST_MITO['CUU'] = 'Thr'
YEAST_MITO['CUC'] = 'Thr'
YEAST_MITO['CUA'] = 'Thr'
YEAST_MITO['CUG'] = 'Thr'

# Echinoderm/Flatworm Mitochondrial Code (NCBI table 9)
# - UGA: Stop → Trp
# - AUA: Ile → Ile (stays)
# - AAA: Lys → Asn
# - AGA: Arg → Ser
# - AGG: Arg → Ser
ECHINO_MITO = dict(STANDARD)
ECHINO_MITO['UGA'] = 'Trp'
ECHINO_MITO['AAA'] = 'Asn'
ECHINO_MITO['AGA'] = 'Ser'
ECHINO_MITO['AGG'] = 'Ser'

# Ciliate Nuclear Code (NCBI table 6)
# - UAA: Stop → Gln
# - UAG: Stop → Gln
CILIATE = dict(STANDARD)
CILIATE['UAA'] = 'Gln'
CILIATE['UAG'] = 'Gln'

# ================================================================
# RUN ALL ANALYSES
# ================================================================

results = []
results.append(analyze_code(STANDARD, "Standard Genetic Code (Table 1)"))
results.append(analyze_code(VERT_MITO, "Vertebrate Mitochondrial (Table 2)"))
results.append(analyze_code(INVERT_MITO, "Invertebrate Mitochondrial (Table 5)"))
results.append(analyze_code(YEAST_MITO, "Yeast Mitochondrial (Table 3)"))
results.append(analyze_code(ECHINO_MITO, "Echinoderm Mitochondrial (Table 9)"))
results.append(analyze_code(CILIATE, "Ciliate Nuclear (Table 6)"))

# ================================================================
# COMPARATIVE SUMMARY
# ================================================================
print("\n\n" + "=" * 70)
print("COMPARATIVE SUMMARY: TOPOLOGICAL SIGNATURES ACROSS GENETIC CODES")
print("=" * 70)

for r in results:
    name = r['name']
    disconn = r['disconnected']
    if disconn:
        disc_str = ", ".join(f"{d['aa']}(ε={d['reconnect_at']})" for d in disconn)
    else:
        disc_str = "None"
    print(f"\n  {name}:")
    print(f"    Amino acids: {r['n_amino_acids']}")
    print(f"    Disconnected at ε=1: {disc_str}")

# Key evolutionary events and their topological consequences
print("\n\n" + "=" * 70)
print("EVOLUTIONARY INTERPRETATION")
print("=" * 70)
print("""
Key codon reassignments and their topological effects:

1. STANDARD → VERTEBRATE MITO:
   - AGA, AGG: Arg → Stop (removes 2 codons from Arg's AG-block)
   - AUA: Ile → Met (transfers codon between amino acids)
   - UGA: Stop → Trp (expands Trp)
   Expected: Arg loses its 2-codon AG-block entirely → becomes 4-fold (CGN only)
             Ser's AG-block remains but loses its Arg neighbors
             
2. STANDARD → INVERTEBRATE MITO:
   - AGA, AGG: Arg → Ser (transfers AG codons from Arg to Ser!)
   - This should RECONNECT Serine's two blocks? Let's check...
   
3. STANDARD → YEAST MITO:
   - CUN: Leu → Thr (the entire 4-codon CU-block of Leu becomes Thr)
   - This should split Leu into just UUA/UUG (2-fold)
   
4. STANDARD → ECHINODERM MITO:
   - AGA, AGG: Arg → Ser (same as invertebrate)
   - AAA: Lys → Asn (expands Asn)
""")

# Specific test: does invertebrate mito reconnect Serine?
print("=" * 70)
print("CRITICAL TEST: Does AGA/AGG → Ser reconnect Serine?")
print("=" * 70)

inv_ser_codons = [c for c, aa in INVERT_MITO.items() if aa == 'Ser']
inv_ser_bits = [codon_to_bits(c) for c in sorted(inv_ser_codons)]
print(f"\nInvertebrate Mito Serine codons: {sorted(inv_ser_codons)}")
for c in sorted(inv_ser_codons):
    print(f"  {c} = {codon_to_bits(c)}")

for eps in [1, 2, 3, 4]:
    comp = connected_components(inv_ser_bits, eps)
    print(f"  β₀(ε={eps}) = {comp}")

# Check: what's the minimum distance from AGA/AGG to the UCN and AGY blocks?
print("\nInter-block distances with new AGA/AGG codons:")
new_codons = ['AGA', 'AGG']
old_ucn = ['UCA', 'UCC', 'UCG', 'UCU']
old_agy = ['AGC', 'AGU']

for nc in new_codons:
    nc_bits = codon_to_bits(nc)
    for oc in old_agy + old_ucn:
        oc_bits = codon_to_bits(oc)
        d = hamming_distance(nc_bits, oc_bits)
        print(f"  d({nc}, {oc}) = {d}")

