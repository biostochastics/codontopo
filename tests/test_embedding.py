from codon_topo.core.embedding import BASE_TO_COMPLEX, embed_codon

I = complex(0, 1)

def test_base_to_complex():
    assert BASE_TO_COMPLEX['C'] == complex(1, 0)
    assert BASE_TO_COMPLEX['U'] == complex(0, 1)
    assert BASE_TO_COMPLEX['A'] == complex(-1, 0)
    assert BASE_TO_COMPLEX['G'] == complex(0, -1)

def test_embed_ggu():
    z = embed_codon('GGU')
    assert z == (complex(0,-1), complex(0,-1), complex(0,1))

def test_twofold_share_first_two_coords():
    """Property 1: two-fold synonymous codons agree in coords 1,2."""
    from codon_topo.core.genetic_codes import STANDARD
    from collections import defaultdict
    aa_codons = defaultdict(list)
    for c, aa in STANDARD.items():
        if aa != 'Stop':
            aa_codons[aa].append(c)
    for aa, codons in aa_codons.items():
        if len(codons) != 2:
            continue
        z1 = embed_codon(codons[0])
        z2 = embed_codon(codons[1])
        assert z1[0] == z2[0], f"{aa} coord 1 differs"
        assert z1[1] == z2[1], f"{aa} coord 2 differs"

def test_fourfold_exhaust_roots():
    """Property 2: four-fold AAs have third coord cycling all 4 roots."""
    from codon_topo.core.genetic_codes import STANDARD
    from collections import defaultdict
    aa_codons = defaultdict(list)
    for c, aa in STANDARD.items():
        if aa != 'Stop':
            aa_codons[aa].append(c)
    fourth_roots = {complex(1,0), complex(0,1), complex(-1,0), complex(0,-1)}
    for aa, codons in aa_codons.items():
        if len(codons) != 4:
            continue
        embeddings = [embed_codon(c) for c in codons]
        first_coords = set((z[0], z[1]) for z in embeddings)
        third_coords = set(z[2] for z in embeddings)
        assert len(first_coords) == 1, f"{aa}: first 2 coords not uniform"
        assert third_coords == fourth_roots, f"{aa}: third coord doesn't exhaust roots"
