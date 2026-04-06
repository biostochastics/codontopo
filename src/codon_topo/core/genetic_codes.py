"""NCBI translation tables.

Reference: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
"""

STANDARD: dict[str, str] = {
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

_CODES: dict[int, tuple[str, dict[str, str]]] = {
    1:  ("Standard", {}),
    2:  ("Vertebrate Mitochondrial", {'AGA':'Stop','AGG':'Stop','AUA':'Met','UGA':'Trp'}),
    3:  ("Yeast Mitochondrial", {'AUA':'Met','CUU':'Thr','CUC':'Thr','CUA':'Thr','CUG':'Thr','UGA':'Trp'}),
    4:  ("Mold/Protozoan/Coelenterate Mito", {'UGA':'Trp'}),
    5:  ("Invertebrate Mitochondrial", {'AGA':'Ser','AGG':'Ser','AUA':'Met','UGA':'Trp'}),
    6:  ("Ciliate Nuclear", {'UAA':'Gln','UAG':'Gln'}),
    9:  ("Echinoderm/Flatworm Mito", {'AAA':'Asn','AGA':'Ser','AGG':'Ser','UGA':'Trp'}),
    10: ("Euplotid Nuclear", {'UGA':'Cys'}),
    11: ("Bacterial/Archaeal/Plant Plastid", {}),
    12: ("Alternative Yeast Nuclear", {'CUG':'Ser'}),
    13: ("Ascidian Mitochondrial", {'AGA':'Gly','AGG':'Gly','AUA':'Met','UGA':'Trp'}),
    14: ("Alternative Flatworm Mito", {'AAA':'Asn','AGA':'Ser','AGG':'Ser','UAA':'Tyr','UGA':'Trp'}),
    16: ("Chlorophycean Mito", {'UAG':'Leu'}),
    21: ("Trematode Mitochondrial", {'AAA':'Asn','AGA':'Ser','AGG':'Ser','AUA':'Met','UGA':'Trp'}),
    22: ("Scenedesmus obliquus Mito", {'UAG':'Leu','UCA':'Stop'}),
    23: ("Thraustochytrium Mito", {'UUA':'Stop','UGA':'Trp'}),
    24: ("Rhabdopleuridae Mito", {'AGA':'Ser','AGG':'Lys','UGA':'Trp'}),
    25: ("Candidate Division SR1 / Gracilibacteria", {'UGA':'Gly'}),
    26: ("Pachysolen tannophilus Nuclear", {'CUG':'Ala'}),
    27: ("Karyorelictea Nuclear", {'UAG':'Gln','UAA':'Gln'}),
    29: ("Mesodinium Nuclear", {'UAA':'Tyr','UAG':'Tyr'}),
    30: ("Peritrich Nuclear", {'UAA':'Glu','UAG':'Glu'}),
    31: ("Blastocrithidia Nuclear", {'UGA':'Trp','UAG':'Glu','UAA':'Glu'}),
    33: ("Cephalodiscidae Mito", {'AGA':'Ser','AGG':'Lys','UAA':'Tyr','UGA':'Trp'}),
}


def get_code(table_id: int) -> dict[str, str]:
    """Return the codon->AA mapping for an NCBI translation table."""
    name, changes = _CODES[table_id]
    code = dict(STANDARD)
    code.update(changes)
    return code


def get_code_name(table_id: int) -> str:
    """Return the human-readable name of a translation table."""
    return _CODES[table_id][0]


def all_table_ids() -> list[int]:
    """Return sorted list of all available NCBI table IDs."""
    return sorted(_CODES.keys())
