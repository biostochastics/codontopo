"""
Comprehensive test: all NCBI translation tables.
For each, compute persistent homology and identify disconnections.
Then correlate reconnection distance with evolutionary depth.
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

def get_blocks(codons):
    """Cluster codons by connectivity at eps=1"""
    bits_list = [codon_to_bits(c) for c in codons]
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
            if hamming_distance(bits_list[i], bits_list[j]) <= 1:
                union(i, j)
    clusters = defaultdict(list)
    for i in range(n):
        clusters[find(i)].append(codons[i])
    return list(clusters.values())

# ================================================================
# ALL 64 CODONS
# ================================================================
ALL_CODONS = []
bases = ['U', 'C', 'A', 'G']
for b1 in bases:
    for b2 in bases:
        for b3 in bases:
            ALL_CODONS.append(b1 + b2 + b3)

# ================================================================
# STANDARD CODE (base for modifications)
# ================================================================
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

def make_code(changes):
    """Make a variant code from standard + dictionary of changes"""
    code = dict(STANDARD)
    for codon, aa in changes.items():
        code[codon] = aa
    return code

# ================================================================
# ALL NCBI TRANSLATION TABLES
# Sources: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
# ================================================================
CODES = {
    1: ("Standard", {}),
    2: ("Vertebrate Mitochondrial", {
        'AGA': 'Stop', 'AGG': 'Stop', 'AUA': 'Met', 'UGA': 'Trp'
    }),
    3: ("Yeast Mitochondrial", {
        'AUA': 'Met', 'CUU': 'Thr', 'CUC': 'Thr', 'CUA': 'Thr', 'CUG': 'Thr', 'UGA': 'Trp'
    }),
    4: ("Mold/Protozoan/Coelenterate Mito", {
        'UGA': 'Trp'
    }),
    5: ("Invertebrate Mitochondrial", {
        'AGA': 'Ser', 'AGG': 'Ser', 'AUA': 'Met', 'UGA': 'Trp'
    }),
    6: ("Ciliate Nuclear", {
        'UAA': 'Gln', 'UAG': 'Gln'
    }),
    9: ("Echinoderm/Flatworm Mito", {
        'AAA': 'Asn', 'AGA': 'Ser', 'AGG': 'Ser', 'UGA': 'Trp'
    }),
    10: ("Euplotid Nuclear", {
        'UGA': 'Cys'
    }),
    11: ("Bacterial/Archaeal/Plant Plastid", {}),  # Same as standard
    12: ("Alternative Yeast Nuclear", {
        'CUG': 'Ser'
    }),
    13: ("Ascidian Mitochondrial", {
        'AGA': 'Gly', 'AGG': 'Gly', 'AUA': 'Met', 'UGA': 'Trp'
    }),
    14: ("Alternative Flatworm Mito", {
        'AAA': 'Asn', 'AGA': 'Ser', 'AGG': 'Ser', 'UAA': 'Tyr', 'UGA': 'Trp'
    }),
    16: ("Chlorophycean Mito", {
        'UAG': 'Leu'
    }),
    21: ("Trematode Mitochondrial", {
        'AAA': 'Asn', 'AGA': 'Ser', 'AGG': 'Ser', 'AUA': 'Met', 'UGA': 'Trp'
    }),
    22: ("Scenedesmus obliquus Mito", {
        'UAG': 'Leu', 'UCA': 'Stop'
    }),
    23: ("Thraustochytrium Mito", {
        'UUA': 'Stop', 'UGA': 'Trp'  # Note: UUA→Stop is unusual
    }),
    24: ("Rhabdopleuridae Mito", {
        'AGA': 'Ser', 'AGG': 'Lys', 'UGA': 'Trp'
    }),
    25: ("Candidate Division SR1 / Gracilibacteria", {
        'UGA': 'Gly'
    }),
    26: ("Pachysolen tannophilus Nuclear", {
        'CUG': 'Ala'
    }),
    27: ("Karyorelictea Nuclear", {
        'UAG': 'Gln', 'UAA': 'Gln', # Same as table 6 but with UGA→Trp sometimes
    }),
    29: ("Mesodinium Nuclear", {
        'UAA': 'Tyr', 'UAG': 'Tyr'
    }),
    30: ("Peritrich Nuclear", {
        'UAA': 'Glu', 'UAG': 'Glu'
    }),
    31: ("Blastocrithidia Nuclear", {
        'UGA': 'Trp', 'UAG': 'Glu', 'UAA': 'Glu'  
    }),
    33: ("Cephalodiscidae Mito", {
        'AGA': 'Ser', 'AGG': 'Lys', 'UAA': 'Tyr', 'UGA': 'Trp'
    }),
}

# ================================================================
# ANALYZE ALL CODES
# ================================================================

all_disconnections = []

print("=" * 80)
print("COMPREHENSIVE ANALYSIS: ALL NCBI TRANSLATION TABLES")
print("=" * 80)

for table_id in sorted(CODES.keys()):
    name, changes = CODES[table_id]
    code = make_code(changes)
    
    # Group by AA
    aa_groups = defaultdict(list)
    for codon, aa in code.items():
        if aa != 'Stop':
            aa_groups[aa].append(codon)
    
    # Check for disconnections
    disconnections = []
    for aa in sorted(aa_groups.keys()):
        codons = sorted(aa_groups[aa])
        if len(codons) < 2:
            continue
        bits_list = [codon_to_bits(c) for c in codons]
        b0_1 = connected_components(bits_list, 1)
        if b0_1 > 1:
            # Find reconnection epsilon
            reconnect = None
            for eps in range(2, 7):
                if connected_components(bits_list, eps) == 1:
                    reconnect = eps
                    break
            
            blocks = get_blocks(codons)
            
            # Min inter-block distance
            min_inter = float('inf')
            for i in range(len(bits_list)):
                for j in range(i+1, len(bits_list)):
                    # Check if in different blocks
                    ci, cj = codons[i], codons[j]
                    in_same = False
                    for block in blocks:
                        if ci in block and cj in block:
                            in_same = True
                            break
                    if not in_same:
                        d = hamming_distance(bits_list[i], bits_list[j])
                        min_inter = min(min_inter, d)
            
            disconnections.append({
                'aa': aa,
                'deg': len(codons),
                'components': b0_1,
                'reconnect_eps': reconnect,
                'min_inter': min_inter,
                'blocks': blocks
            })
    
    # Check filtration preservation
    twofold_pass = 0
    twofold_fail = 0
    fourfold_pass = 0
    fourfold_fail = 0
    
    for aa, codons in aa_groups.items():
        bits_list = [codon_to_bits(c) for c in codons]
        if len(codons) == 2:
            differing = [i for i in range(6) if bits_list[0][i] != bits_list[1][i]]
            if differing == [5]:
                twofold_pass += 1
            else:
                twofold_fail += 1
        elif len(codons) == 4:
            prefixes = set(b[:4] for b in bits_list)
            suffixes = set(b[4:] for b in bits_list)
            if len(prefixes) == 1 and suffixes == {(0,0),(0,1),(1,0),(1,1)}:
                fourfold_pass += 1
            else:
                fourfold_fail += 1
    
    # Print summary
    disc_str = ", ".join(f"{d['aa']}(ε={d['reconnect_eps']}, d={d['min_inter']}, {d['components']}comp)" 
                         for d in disconnections) if disconnections else "None"
    
    changes_str = ", ".join(f"{c}→{a}" for c, a in sorted(changes.items())) if changes else "identical to standard"
    
    print(f"\n  Table {table_id}: {name}")
    print(f"    Changes: {changes_str}")
    print(f"    2-fold filtration: {twofold_pass}/{twofold_pass+twofold_fail} pass")
    print(f"    4-fold filtration: {fourfold_pass}/{fourfold_pass+fourfold_fail} pass")
    print(f"    Disconnected AAs: {disc_str}")
    
    for d in disconnections:
        all_disconnections.append({
            'table': table_id,
            'code_name': name,
            'aa': d['aa'],
            'deg': d['deg'],
            'components': d['components'],
            'reconnect_eps': d['reconnect_eps'],
            'min_inter': d['min_inter'],
            'blocks': d['blocks']
        })

# ================================================================
# DISCONNECTION CATALOGUE
# ================================================================
print("\n\n" + "=" * 80)
print("COMPLETE DISCONNECTION CATALOGUE")
print("=" * 80)
print(f"\n{'Table':<8} {'Code':<35} {'AA':<6} {'Deg':<5} {'Comp':<6} {'ε_reconnect':<13} {'d_min':<7} {'Blocks'}")
print("-" * 120)

for d in sorted(all_disconnections, key=lambda x: (x['aa'], x['reconnect_eps'], x['table'])):
    blocks_str = " | ".join(str(b) for b in d['blocks'])
    print(f"  {d['table']:<6} {d['code_name']:<35} {d['aa']:<6} {d['deg']:<5} {d['components']:<6} {d['reconnect_eps']:<13} {d['min_inter']:<7} {blocks_str}")

# ================================================================
# NOVEL DISCONNECTIONS (excluding Serine, which is universal)
# ================================================================
print("\n\n" + "=" * 80)
print("NOVEL DISCONNECTIONS (non-Serine)")
print("=" * 80)

novel = [d for d in all_disconnections if d['aa'] != 'Ser']
if novel:
    for d in novel:
        print(f"\n  Table {d['table']}: {d['code_name']}")
        print(f"    {d['aa']} ({d['deg']}-fold, {d['components']} components)")
        print(f"    Reconnects at ε={d['reconnect_eps']}, min inter-block distance={d['min_inter']}")
        print(f"    Blocks:")
        for i, block in enumerate(d['blocks']):
            bits = [codon_to_bits(c) for c in block]
            print(f"      Block {i+1}: {block}")
            for c, b in zip(block, bits):
                print(f"        {c} = {b}")
else:
    print("  None found beyond Serine.")

# ================================================================
# FILTRATION BREAKAGE ANALYSIS
# ================================================================
print("\n\n" + "=" * 80)
print("FILTRATION BREAKAGE: WHERE DOES THE STANDARD STRUCTURE FAIL?")
print("=" * 80)

for table_id in sorted(CODES.keys()):
    name, changes = CODES[table_id]
    if not changes:
        continue
    code = make_code(changes)
    
    aa_groups = defaultdict(list)
    for codon, aa in code.items():
        if aa != 'Stop':
            aa_groups[aa].append(codon)
    
    failures = []
    for aa, codons in aa_groups.items():
        bits_list = [codon_to_bits(c) for c in sorted(codons)]
        if len(codons) == 2:
            differing = [i for i in range(6) if bits_list[0][i] != bits_list[1][i]]
            if differing != [5]:
                failures.append(f"  2-fold {aa}: {sorted(codons)} differ at {differing} (not [5])")
        elif len(codons) == 4:
            prefixes = set(b[:4] for b in bits_list)
            suffixes = set(b[4:] for b in bits_list)
            if len(prefixes) != 1 or suffixes != {(0,0),(0,1),(1,0),(1,1)}:
                failures.append(f"  4-fold {aa}: prefixes={prefixes}, suffixes={suffixes}")
    
    if failures:
        print(f"\n  Table {table_id}: {name}")
        for f in failures:
            print(f"  {f}")

