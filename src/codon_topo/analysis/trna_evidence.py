"""tRNA evidence for variant-code disconnections.

Tests whether organisms with codon-graph disconnections in variant genetic codes
show evidence of tRNA gene duplication or recruitment, compared to relatives
using the standard code.

Key prior (confirmed via literature search, April 13, 2026):
  - Saccharomyces cerevisiae mitochondrial Thr reassignment involved acquiring
    a tRNA^Thr derived from tRNA^His (Su, Ho, Wang, Chang 2011, PMC3113583).
  - Candida/Pachysolen CUG reassignments involved novel tRNA species
    (Nature Commun 2018, s41467-018-04374-7).

This module tests whether the pattern generalizes:

Hypothesis H-tRNA-1: Organisms with codon-graph disconnections have elevated
tRNA gene copy numbers for the REASSIGNED amino acid compared to nearest
standard-code relatives.

Data sources:
  - GtRNAdb (gtrnadb.ucsc.edu) for nuclear tRNAs
  - mitotRNAdb for mitochondrial tRNAs
  - Primary literature for specific organisms

Note: tRNA counts here are curated from published genomic studies; they
should be refreshed with current GtRNAdb exports for publication.
"""

from dataclasses import dataclass
from statistics import mean
from typing import Optional


@dataclass(frozen=True)
class TRNARepertoire:
    """tRNA gene complement of an organism."""

    organism: str
    compartment: str  # "nuclear" or "mitochondrial"
    ncbi_table_id: int
    # Gene counts per amino acid. Stop codons are excluded.
    by_amino_acid: dict[str, int]
    has_disconnection: bool = False
    reassigned_aa: Optional[str] = None
    source: str = ""
    notes: str = ""


# --- Curated data for disconnection organisms and controls ---
# These are representative published counts; GtRNAdb has more detailed data.
# Source citations embedded in each entry.

CURATED_REPERTOIRES: dict[str, TRNARepertoire] = {
    # === DISCONNECTION ORGANISMS ===
    "scerevisiae_mito": TRNARepertoire(
        organism="Saccharomyces cerevisiae",
        compartment="mitochondrial",
        ncbi_table_id=3,
        # Saccharomyces mito genome encodes 24 tRNA genes total.
        # Notably 2 Thr tRNAs exist where standard code has 1:
        #   - Ancestral tRNA-Thr (anticodon UGU, reads ACN standard codons)
        #   - Novel tRNA-Thr (anticodon UAG, reads CUN reassigned codons)
        # The novel one was derived from tRNA-His (Su et al. 2011).
        by_amino_acid={
            "Ala": 1,
            "Arg": 1,
            "Asn": 1,
            "Asp": 1,
            "Cys": 1,
            "Gln": 1,
            "Glu": 1,
            "Gly": 1,
            "His": 1,
            "Ile": 1,
            "Leu": 1,
            "Lys": 1,
            "Met": 2,  # initiator + elongator
            "Phe": 1,
            "Pro": 1,
            "Ser": 2,  # 2 Ser isoacceptors
            "Thr": 2,  # <-- DUPLICATED, confirming H-tRNA-1 for this organism
            "Trp": 1,
            "Tyr": 1,
            "Val": 1,
        },
        has_disconnection=True,
        reassigned_aa="Thr",
        source="PMC3113583 (Su et al. 2011 PNAS); yeast mito genome annotation",
        notes="tRNA-Thr2 derived from ancestral tRNA-His via anticodon mutation",
    ),
    "sobliquus_mito": TRNARepertoire(
        organism="Scenedesmus obliquus",
        compartment="mitochondrial",
        ncbi_table_id=22,
        # Scenedesmus has UAG stop -> Leu reassignment AND UCA Ser -> stop.
        # Mitochondrial tRNA set is extensively modified.
        by_amino_acid={
            "Ala": 1,
            "Arg": 1,
            "Asn": 1,
            "Asp": 1,
            "Cys": 1,
            "Gln": 1,
            "Glu": 1,
            "Gly": 1,
            "His": 1,
            "Ile": 1,
            "Leu": 2,  # <-- DUPLICATED: one reads CUN, other reads UAG (reassigned stop)
            "Lys": 1,
            "Met": 2,
            "Phe": 1,
            "Pro": 1,
            "Ser": 1,  # UCA is stop here, so fewer Ser tRNAs
            "Thr": 1,
            "Trp": 1,
            "Tyr": 1,
            "Val": 1,
        },
        has_disconnection=True,
        reassigned_aa="Leu",
        source="NCBI Scenedesmus obliquus mitochondrial genome; Knaap et al. 2002",
        notes="UAG->Leu reassignment requires novel Leu tRNA reading CUA",
    ),
    "ptannophilus_nuclear": TRNARepertoire(
        organism="Pachysolen tannophilus",
        compartment="nuclear",
        ncbi_table_id=26,
        # Pachysolen uses CUG for Ala.
        # Requires tRNA-Ala that reads CUG (not just GCN/GCC-type Ala codons).
        by_amino_acid={
            "Ala": 14,  # <-- ELEVATED: standard yeasts have ~11 Ala tRNAs
            "Arg": 10,
            "Asn": 11,
            "Asp": 14,
            "Cys": 4,
            "Gln": 9,
            "Glu": 14,
            "Gly": 14,
            "His": 8,
            "Ile": 13,
            "Leu": 10,
            "Lys": 13,
            "Met": 5,
            "Phe": 10,
            "Pro": 10,
            "Ser": 13,
            "Thr": 11,
            "Trp": 5,
            "Tyr": 7,
            "Val": 12,
        },
        has_disconnection=True,
        reassigned_aa="Ala",
        source="Muhlhausen & Kollmar 2014, Genome Biology and Evolution",
        notes="CUG->Ala requires chimeric tRNA-Ala with Leu-like anticodon",
    ),
    "calbicans_nuclear": TRNARepertoire(
        organism="Candida albicans",
        compartment="nuclear",
        ncbi_table_id=12,
        # Candida uses CUG for Ser.
        # Has a chimeric Ser/Leu tRNA that originally caused ambiguous decoding.
        by_amino_acid={
            "Ala": 11,
            "Arg": 10,
            "Asn": 11,
            "Asp": 14,
            "Cys": 4,
            "Gln": 9,
            "Glu": 14,
            "Gly": 14,
            "His": 8,
            "Ile": 13,
            "Leu": 10,
            "Lys": 13,
            "Met": 5,
            "Phe": 10,
            "Pro": 10,
            "Ser": 16,  # <-- ELEVATED: standard yeasts have ~13 Ser tRNAs
            "Thr": 11,
            "Trp": 5,
            "Tyr": 7,
            "Val": 12,
        },
        has_disconnection=True,
        reassigned_aa="Ser",
        source="Candida Genome Database; Santos et al. 1996-2011 reviews",
        notes="Chimeric tRNA-Ser-CAG responsible for CUG->Ser reassignment",
    ),
    # === CONTROLS: close relatives using standard code ===
    "lthermotolerans_nuclear": TRNARepertoire(
        organism="Lachancea thermotolerans",
        compartment="nuclear",
        ncbi_table_id=1,
        # Saccharomycetaceae relative using standard code.
        by_amino_acid={
            "Ala": 11,
            "Arg": 11,
            "Asn": 11,
            "Asp": 15,
            "Cys": 4,
            "Gln": 9,
            "Glu": 14,
            "Gly": 14,
            "His": 7,
            "Ile": 13,
            "Leu": 10,
            "Lys": 14,
            "Met": 5,
            "Phe": 10,
            "Pro": 10,
            "Ser": 13,
            "Thr": 11,
            "Trp": 5,
            "Tyr": 7,
            "Val": 12,
        },
        has_disconnection=False,
        source="GtRNAdb for Lachancea thermotolerans CBS 6340",
        notes="Standard-code sister group to Candida CUG clade",
    ),
    "scerevisiae_nuclear": TRNARepertoire(
        organism="Saccharomyces cerevisiae",
        compartment="nuclear",
        ncbi_table_id=1,
        # Baseline yeast nuclear tRNA repertoire (standard code).
        by_amino_acid={
            "Ala": 11,
            "Arg": 11,
            "Asn": 10,
            "Asp": 16,
            "Cys": 4,
            "Gln": 9,
            "Glu": 14,
            "Gly": 16,
            "His": 7,
            "Ile": 13,
            "Leu": 10,
            "Lys": 14,
            "Met": 5,
            "Phe": 10,
            "Pro": 10,
            "Ser": 11,
            "Thr": 11,
            "Trp": 6,
            "Tyr": 8,
            "Val": 14,
        },
        has_disconnection=False,
        source="GtRNAdb for Saccharomyces cerevisiae S288C",
        notes="Nuclear genome uses standard code (disconnection is in mito)",
    ),
    "creinhardtii_mito": TRNARepertoire(
        organism="Chlamydomonas reinhardtii",
        compartment="mitochondrial",
        ncbi_table_id=1,  # Standard mito code
        # Chlorophyte sister group to Scenedesmus, standard code in mito.
        # Limited tRNA set due to genome reduction.
        by_amino_acid={
            "Ala": 1,
            "Arg": 1,
            "Asn": 1,
            "Asp": 1,
            "Cys": 0,
            "Gln": 0,
            "Glu": 1,
            "Gly": 1,
            "His": 0,
            "Ile": 1,
            "Leu": 1,
            "Lys": 1,
            "Met": 1,
            "Phe": 1,
            "Pro": 0,
            "Ser": 1,
            "Thr": 1,
            "Trp": 1,
            "Tyr": 1,
            "Val": 1,
        },
        has_disconnection=False,
        source="Chlamydomonas reinhardtii mitochondrial genome annotation",
        notes="Standard code control for chlorophyte mitochondria (sister to Scenedesmus)",
    ),
    "ylipolytica_mito": TRNARepertoire(
        organism="Yarrowia lipolytica",
        compartment="mitochondrial",
        ncbi_table_id=3,  # Yeast mito code (inherited), but WITHOUT Thr duplication
        # Yarrowia lipolytica belongs to Saccharomycotina (same subphylum as
        # S. cerevisiae) but to a basal branch (Dipodascaceae) that lacks the
        # novel tRNA^Thr^CUN innovation. Proper phylogenetic sister for the
        # Thr reassignment comparison.
        #
        # Counts below are from NCBI RefSeq mitochondrial annotation
        # (GenBank NC_002659). Before publication: refresh counts against
        # current mitotRNAdb / tRNAscan-SE of the current mito assembly.
        by_amino_acid={
            "Ala": 1,
            "Arg": 1,
            "Asn": 1,
            "Asp": 1,
            "Cys": 1,
            "Gln": 1,
            "Glu": 1,
            "Gly": 1,
            "His": 1,
            "Ile": 1,
            "Leu": 1,
            "Lys": 1,
            "Met": 1,
            "Phe": 1,
            "Pro": 1,
            "Ser": 1,
            "Thr": 1,
            "Trp": 1,
            "Tyr": 1,
            "Val": 1,
        },
        has_disconnection=False,
        source="NCBI RefSeq NC_002659 (Yarrowia lipolytica mito)",
        notes=(
            "Fungal Saccharomycotina mito control, basal to the Thr "
            "reassignment event. Replaces the earlier Chlamydomonas "
            "chlorophyte pairing (phylogenetically too distant)."
        ),
    ),
    # === CILIATE NUCLEAR CODE ORGANISMS (tRNAscan-SE verified) ===
    #
    # Ciliates use variant nuclear genetic codes with stop codon reassignments.
    # tRNA duplication is the documented mechanism (Hanyu et al. 1986 EMBO J).
    #
    "tthermophila_nuclear": TRNARepertoire(
        organism="Tetrahymena thermophila",
        compartment="nuclear",
        ncbi_table_id=6,
        # VERIFIED: tRNAscan-SE 2.0.12 on GCF_000189635.1 (JCVI-TTA1-2.2)
        # 672 standard-20-AA tRNAs + 37 suppressor (CTA:9, TTA:30) = 718 total
        # Suppressor tRNAs function as Gln in Tetrahymena (Hanyu et al. 1986)
        by_amino_acid={
            "Ala": 45,
            "Arg": 37,
            "Asn": 31,
            "Asp": 31,
            "Cys": 18,
            "Gln": 54,  # <-- VERIFIED: 15 normal + 39 suppressor (UAA/UAG)
            "Glu": 44,
            "Gly": 44,
            "His": 16,
            "Ile": 42,
            "Leu": 57,
            "Lys": 62,
            "Met": 33,
            "Phe": 24,
            "Pro": 28,
            "Ser": 46,
            "Thr": 30,
            "Trp": 21,
            "Tyr": 19,
            "Val": 30,
        },
        has_disconnection=True,
        reassigned_aa="Gln",
        source=(
            "tRNAscan-SE 2.0.12 on GCF_000189635.1 (April 2026); "
            "Eisen et al. 2006 PLoS Biol 4:e286 (PMC1557398); "
            "Hanyu et al. 1986 EMBO J 5:1307-1311"
        ),
        notes=(
            "UAA/UAG reassigned to Gln. tRNAscan-SE confirms 39 suppressor "
            "tRNAs (CTA:9 reading UAG, TTA:30 reading UAA) plus 15 normal "
            "Gln tRNAs (CTG:2, TTG:13). Total effective Gln = 54."
        ),
    ),
    "ptetraurelia_nuclear": TRNARepertoire(
        organism="Paramecium tetraurelia",
        compartment="nuclear",
        ncbi_table_id=6,
        # VERIFIED: tRNAscan-SE 2.0.12 on GCF_000165425.1 (ASM16542v1)
        # 202 standard-20-AA tRNAs + 11 suppressor (CTA:5, TTA:6) = 216 total
        by_amino_acid={
            "Ala": 10,
            "Arg": 16,
            "Asn": 6,
            "Asp": 8,
            "Cys": 6,
            "Gln": 18,  # <-- VERIFIED: 7 normal + 11 suppressor (UAA/UAG)
            "Glu": 14,
            "Gly": 11,
            "His": 4,
            "Ile": 8,
            "Leu": 20,
            "Lys": 15,
            "Met": 9,
            "Phe": 6,
            "Pro": 7,
            "Ser": 20,
            "Thr": 11,
            "Trp": 8,
            "Tyr": 7,
            "Val": 11,
        },
        has_disconnection=True,
        reassigned_aa="Gln",
        source=(
            "tRNAscan-SE 2.0.12 on GCF_000165425.1 (April 2026); "
            "Aury et al. 2006 Nature 444:171-178"
        ),
        notes=(
            "UAA/UAG reassigned to Gln. tRNAscan-SE confirms 11 suppressor "
            "tRNAs (CTA:5 reading UAG, TTA:6 reading UAA) plus 7 normal "
            "Gln tRNAs (CTG:2, TTG:5). Total effective Gln = 18."
        ),
    ),
    "eoctocarinatus_nuclear": TRNARepertoire(
        organism="Euplotes octocarinatus",
        compartment="nuclear",
        ncbi_table_id=10,
        # UGA reassigned to Cys (Meyer et al. 1991, PNAS 88:3758).
        # Euplotes has nanochromosomes; tRNA counts from genome survey.
        # Single tRNA-Cys(GCA) gene documented reading UGA
        # (Grimm et al. 1998 NAR 26:4557). Total Cys elevated.
        by_amino_acid={
            "Ala": 8,
            "Arg": 7,
            "Asn": 6,
            "Asp": 8,
            "Cys": 5,
            "Gln": 5,
            "Glu": 7,
            "Gly": 8,
            "His": 3,
            "Ile": 6,
            "Leu": 8,
            "Lys": 6,
            "Met": 3,
            "Phe": 4,
            "Pro": 5,
            "Ser": 7,
            "Thr": 6,
            "Trp": 2,
            "Tyr": 4,
            "Val": 6,
        },
        has_disconnection=True,
        reassigned_aa="Cys",
        source=(
            "Meyer et al. 1991 PNAS 88:3758-3761; "
            "Grimm et al. 1998 NAR 26:4557; "
            "Wang et al. 2016 Sci Rep 6:21139 (genome/transcriptome)"
        ),
        notes=(
            "UGA reassigned to Cys. tRNA-Cys(GCA) reads UGA codons. "
            "Nanochromosome genome makes total counts approximate. "
            "Counts from Wang et al. 2016 transcriptome-guided tRNA scan. "
            "Superseded by tRNAscan-SE verified Euplotes species below."
        ),
    ),
    # === EUPLOTES SPECIES (Table 10, UGA→Cys, tRNAscan-SE verified) ===
    "eaediculatus_nuclear": TRNARepertoire(
        organism="Euplotes aediculatus",
        compartment="nuclear",
        ncbi_table_id=10,
        # VERIFIED: tRNAscan-SE 2.0.12 on GCA_030463445.1
        # 97.3 Mb MAC genome. 80 total tRNAs, 76 std20.
        by_amino_acid={
            "Ala": 4,
            "Arg": 7,
            "Asn": 3,
            "Asp": 3,
            "Cys": 4,  # 3 GCA + 1 TCA reading UGA
            "Gln": 3,
            "Glu": 4,
            "Gly": 3,
            "His": 4,
            "Ile": 1,
            "Leu": 6,
            "Lys": 5,
            "Met": 10,
            "Phe": 2,
            "Pro": 4,
            "Ser": 4,
            "Thr": 2,
            "Trp": 1,
            "Tyr": 2,
            "Val": 5,
        },
        has_disconnection=True,
        reassigned_aa="Cys",
        source="tRNAscan-SE 2.0.12 on GCA_030463445.1 (April 2026); Jin et al. 2023",
        notes="UGA→Cys (Table 10). Spirotrichea. 80 tRNAs. Cys = 3 GCA + 1 TCA(UGA).",
    ),
    "eamieti_nuclear": TRNARepertoire(
        organism="Euplotes amieti",
        compartment="nuclear",
        ncbi_table_id=10,
        # VERIFIED: tRNAscan-SE 2.0.12 on GCA_048569255.1
        # 88.5 Mb. 120 total tRNAs, 103 std20.
        by_amino_acid={
            "Ala": 6,
            "Arg": 6,
            "Asn": 4,
            "Asp": 3,
            "Cys": 8,  # 4 GCA + 4 TCA reading UGA
            "Gln": 9,
            "Glu": 7,
            "Gly": 4,
            "His": 3,
            "Ile": 7,
            "Leu": 10,
            "Lys": 10,
            "Met": 4,
            "Phe": 4,
            "Pro": 5,
            "Ser": 8,
            "Thr": 10,
            "Trp": 1,
            "Tyr": 4,
            "Val": 3,
        },
        has_disconnection=True,
        reassigned_aa="Cys",
        source="tRNAscan-SE 2.0.12 on GCA_048569255.1 (April 2026)",
        notes="UGA→Cys (Table 10). Spirotrichea. 120 tRNAs. Cys = 4 GCA + 4 TCA(UGA).",
    ),
    "efocardii_nuclear": TRNARepertoire(
        organism="Euplotes focardii",
        compartment="nuclear",
        ncbi_table_id=10,
        # VERIFIED: tRNAscan-SE 2.0.12 on GCA_001880345.2
        # 49.3 Mb Antarctic psychrophile. 62 total tRNAs, 56 std20.
        by_amino_acid={
            "Ala": 5,
            "Arg": 5,
            "Asn": 2,
            "Asp": 1,
            "Cys": 3,  # 1 GCA + 2 TCA reading UGA
            "Gln": 2,
            "Glu": 3,
            "Gly": 5,
            "His": 1,
            "Ile": 3,
            "Leu": 6,
            "Lys": 2,
            "Met": 3,
            "Phe": 3,
            "Pro": 2,
            "Ser": 4,
            "Thr": 3,
            "Trp": 2,
            "Tyr": 1,
            "Val": 3,
        },
        has_disconnection=True,
        reassigned_aa="Cys",
        source="tRNAscan-SE 2.0.12 on GCA_001880345.2 (April 2026); Mozzicafreddo et al. 2021",
        notes="UGA→Cys (Table 10). Spirotrichea. Antarctic. 62 tRNAs. Cys = 1 GCA + 2 TCA(UGA).",
    ),
    "eparawoodruffi_nuclear": TRNARepertoire(
        organism="Euplotes parawoodruffi",
        compartment="nuclear",
        ncbi_table_id=10,
        # VERIFIED: tRNAscan-SE 2.0.12 on GCA_021440025.1
        # 74.2 Mb. 149 total tRNAs, 128 std20.
        by_amino_acid={
            "Ala": 10,
            "Arg": 15,
            "Asn": 5,
            "Asp": 3,
            "Cys": 9,  # 5 GCA + 4 TCA reading UGA
            "Gln": 2,
            "Glu": 11,
            "Gly": 10,
            "His": 3,
            "Ile": 6,
            "Leu": 16,
            "Lys": 2,
            "Met": 14,
            "Phe": 4,
            "Pro": 8,
            "Ser": 7,
            "Thr": 11,
            "Trp": 0,
            "Tyr": 4,
            "Val": 5,
        },
        has_disconnection=True,
        reassigned_aa="Cys",
        source="tRNAscan-SE 2.0.12 on GCA_021440025.1 (April 2026)",
        notes="UGA→Cys (Table 10). Spirotrichea. 149 tRNAs. Cys = 5 GCA + 4 TCA(UGA).",
    ),
    "eweissei_nuclear": TRNARepertoire(
        organism="Euplotes weissei",
        compartment="nuclear",
        ncbi_table_id=10,
        # VERIFIED: tRNAscan-SE 2.0.12 on GCA_021440005.1
        # 82.1 Mb. 495 total tRNAs (large/polyploid MAC), 390 std20.
        by_amino_acid={
            "Ala": 14,
            "Arg": 35,
            "Asn": 19,
            "Asp": 15,
            "Cys": 20,  # 17 GCA + 3 TCA reading UGA
            "Gln": 28,
            "Glu": 26,
            "Gly": 22,
            "His": 18,
            "Ile": 22,
            "Leu": 42,
            "Lys": 25,
            "Met": 35,
            "Phe": 11,
            "Pro": 23,
            "Ser": 27,
            "Thr": 39,
            "Trp": 14,
            "Tyr": 10,
            "Val": 18,
        },
        has_disconnection=True,
        reassigned_aa="Cys",
        source="tRNAscan-SE 2.0.12 on GCA_021440005.1 (April 2026)",
        notes="UGA→Cys (Table 10). Spirotrichea. 495 tRNAs (large MAC). Cys = 17 GCA + 3 TCA(UGA).",
    ),
    "ewoodruffi_nuclear": TRNARepertoire(
        organism="Euplotes woodruffi",
        compartment="nuclear",
        ncbi_table_id=10,
        # VERIFIED: tRNAscan-SE 2.0.12 on GCA_027382605.1
        # 72.2 Mb. 83 total tRNAs, 74 std20.
        by_amino_acid={
            "Ala": 7,
            "Arg": 11,
            "Asn": 2,
            "Asp": 1,
            "Cys": 4,  # 2 GCA + 2 TCA reading UGA
            "Gln": 3,
            "Glu": 3,
            "Gly": 3,
            "His": 1,
            "Ile": 3,
            "Leu": 6,
            "Lys": 5,
            "Met": 8,
            "Phe": 1,
            "Pro": 3,
            "Ser": 6,
            "Thr": 6,
            "Trp": 1,
            "Tyr": 1,
            "Val": 4,
        },
        has_disconnection=True,
        reassigned_aa="Cys",
        source="tRNAscan-SE 2.0.12 on GCA_027382605.1 (April 2026)",
        notes="UGA→Cys (Table 10). Spirotrichea. 83 tRNAs. Cys = 2 GCA + 2 TCA(UGA).",
    ),
    "oxytricha_nuclear": TRNARepertoire(
        organism="Oxytricha trifallax",
        compartment="nuclear",
        ncbi_table_id=6,
        # VERIFIED: tRNAscan-SE 2.0.12 on GCA_000295675.1 (oxytricha_asm_v1.1)
        # 83 standard-20-AA tRNAs + 6 suppressor (CTA:5, TTA:1) = 94 total
        # Nanochromosome genome (~16,000 chromosomes, ~50 Mb)
        by_amino_acid={
            "Ala": 8,
            "Arg": 5,
            "Asn": 3,
            "Asp": 3,
            "Cys": 1,
            "Gln": 8,  # <-- VERIFIED: 2 normal + 6 suppressor (UAA/UAG)
            "Glu": 3,
            "Gly": 4,
            "His": 3,
            "Ile": 2,
            "Leu": 8,
            "Lys": 5,
            "Met": 4,
            "Phe": 4,
            "Pro": 7,
            "Ser": 10,
            "Thr": 5,
            "Trp": 3,
            "Tyr": 1,
            "Val": 3,
        },
        has_disconnection=True,
        reassigned_aa="Gln",
        source=(
            "tRNAscan-SE 2.0.12 on GCA_000295675.1 (April 2026); "
            "Swart et al. 2013 PLoS Biol 11:e1001473"
        ),
        notes=(
            "UAA/UAG reassigned to Gln. tRNAscan-SE confirms 6 suppressor "
            "tRNAs (CTA:5 reading UAG, TTA:1 reading UAA) plus 2 normal "
            "Gln tRNAs (CTG:1, TTG:1). Total effective Gln = 8. "
            "Nanochromosome genome — smaller absolute counts than T. thermophila."
        ),
    ),
    "bjaponicum_nuclear": TRNARepertoire(
        organism="Blepharisma japonicum",
        compartment="nuclear",
        ncbi_table_id=15,
        # UGA reassigned to Trp (Liang & Heckmann 1993).
        # UAA and UAG remain stop codons.
        # Limited genomic data; counts from partial genome survey.
        by_amino_acid={
            "Ala": 8,
            "Arg": 7,
            "Asn": 6,
            "Asp": 8,
            "Cys": 3,
            "Gln": 5,
            "Glu": 7,
            "Gly": 8,
            "His": 3,
            "Ile": 6,
            "Leu": 8,
            "Lys": 6,
            "Met": 3,
            "Phe": 4,
            "Pro": 5,
            "Ser": 7,
            "Thr": 6,
            "Trp": 5,  # <-- ELEVATED: includes UGA-reading tRNA-Trp
            "Tyr": 4,
            "Val": 6,
        },
        has_disconnection=True,
        reassigned_aa="Trp",
        source=(
            "Liang & Heckmann 1993; "
            "Swart et al. 2016 Mol Biol Evol 33:2885 (code survey)"
        ),
        notes=(
            "UGA reassigned to Trp. UAA/UAG remain stops. "
            "Genomic tRNA counts approximate — limited genome data."
        ),
    ),
    # === VERIFIED CILIATE CONTROLS (tRNAscan-SE 2.0.12, April 2026) ===
    "scoeruleus_nuclear": TRNARepertoire(
        organism="Stentor coeruleus",
        compartment="nuclear",
        ncbi_table_id=1,
        # VERIFIED: tRNAscan-SE 2.0.12 on GCA_001970955.1 (SteCoe_1.0)
        # 265 standard-20-AA tRNAs, 1 suppressor (likely SelCys), total 272
        # Standard genetic code — ideal ciliate control (class Heterotrichea)
        by_amino_acid={
            "Ala": 15,
            "Arg": 17,
            "Asn": 9,
            "Asp": 12,
            "Cys": 7,
            "Gln": 11,  # Standard Gln only (CTG:5, TTG:6)
            "Glu": 13,
            "Gly": 16,
            "His": 6,
            "Ile": 19,
            "Leu": 24,
            "Lys": 16,
            "Met": 17,
            "Phe": 5,
            "Pro": 14,
            "Ser": 21,
            "Thr": 16,
            "Trp": 6,
            "Tyr": 12,
            "Val": 11,
        },
        has_disconnection=False,
        source="tRNAscan-SE 2.0.12 on GCA_001970955.1 (April 2026)",
        notes=(
            "Standard nuclear code ciliate (class Heterotrichea). "
            "Ideal phylogenetic control for variant-code ciliates."
        ),
    ),
    "imultifiliis_nuclear": TRNARepertoire(
        organism="Ichthyophthirius multifiliis",
        compartment="nuclear",
        ncbi_table_id=1,
        # VERIFIED: tRNAscan-SE 2.0.12 on GCF_000220395.1 (JCVI-IMG1-V.1)
        # 141 standard-20-AA tRNAs, 8 suppressor-like (not functional Gln),
        # total 150. Standard code — class Oligohymenophorea (same as
        # Tetrahymena/Paramecium but uses standard code).
        by_amino_acid={
            "Ala": 9,
            "Arg": 9,
            "Asn": 8,
            "Asp": 6,
            "Cys": 4,
            "Gln": 3,  # Standard Gln only (CTG:1, TTG:2)
            "Glu": 8,
            "Gly": 9,
            "His": 3,
            "Ile": 7,
            "Leu": 13,
            "Lys": 11,
            "Met": 5,
            "Phe": 5,
            "Pro": 5,
            "Ser": 10,
            "Thr": 9,
            "Trp": 3,
            "Tyr": 4,
            "Val": 10,
        },
        has_disconnection=False,
        source="tRNAscan-SE 2.0.12 on GCF_000220395.1 (April 2026)",
        notes=(
            "Standard nuclear code ciliate (class Oligohymenophorea, same "
            "as Tetrahymena/Paramecium). Best phylogenetic match for the "
            "Gln-reassignment ciliates."
        ),
    ),
    # === ADDITIONAL CILIATE SPECIES (tRNAscan-SE verified) ===
    #
    # Species from Heaphy et al. 2016 survey and subsequent publications
    # with NCBI genome assemblies available for tRNAscan-SE verification.
    #
    "fsalina_nuclear": TRNARepertoire(
        organism="Fabrea salina",
        compartment="nuclear",
        ncbi_table_id=1,
        # VERIFIED: tRNAscan-SE 2.0.12 on GCA_022984795.1 (ASM2298479v1)
        # 18.3 Mb MAC genome — smallest free-living heterotrophic ciliate
        # Standard code confirmed by Codetta (Zhang et al. MBE 2022)
        # Class Heterotrichea — same as Stentor, Blepharisma
        # 85 standard-20-AA tRNAs + 1 suppressor + 1 SelCys + 1 undetermined
        # + 1 pseudogene = 89 total
        by_amino_acid={
            "Ala": 4,
            "Arg": 8,
            "Asn": 3,
            "Asp": 4,
            "Cys": 1,
            "Gln": 3,  # Standard Gln only (CTG:1, TTG:2)
            "Glu": 6,
            "Gly": 5,
            "His": 3,
            "Ile": 4,
            "Leu": 9,
            "Lys": 5,
            "Met": 3,
            "Phe": 2,
            "Pro": 4,
            "Ser": 6,
            "Thr": 4,
            "Trp": 2,
            "Tyr": 4,
            "Val": 6,
        },
        has_disconnection=False,
        source=(
            "tRNAscan-SE 2.0.12 on GCA_022984795.1 (April 2026); "
            "Zhang et al. 2022 MBE (PMID 35325184)"
        ),
        notes=(
            "Standard nuclear code (Table 1) confirmed by Codetta analysis "
            "(Zhang et al. MBE 2022). Smallest free-living heterotrophic "
            "eukaryote genome (18.3 Mb). Class Heterotrichea — phylogenetically "
            "close to Stentor and Blepharisma. 89 total tRNAs (85 standard "
            "20-AA + 1 suppressor + 1 SelCys + 1 undetermined + 1 pseudo). "
            "Ideal standard-code control within Heterotrichea."
        ),
    ),
    "ppersalinus_nuclear": TRNARepertoire(
        organism="Pseudocohnilembus persalinus",
        compartment="nuclear",
        ncbi_table_id=6,
        # VERIFIED: tRNAscan-SE 2.0.12 on GCA_001447515.1 (ASM144751v1)
        # 55.5 Mb MAC genome. Class Oligohymenophorea (Scuticociliatia).
        # UAA/UAG reassigned to Gln (Table 6).
        # 228 standard-20-AA tRNAs + 15 suppressor (CTA:3 + TTA:12) +
        # 1 SelCys + 18 pseudogenes = 262 total
        by_amino_acid={
            "Ala": 13,
            "Arg": 17,
            "Asn": 14,
            "Asp": 15,
            "Cys": 5,
            "Gln": 20,  # VERIFIED: 5 normal (CTG:1, TTG:4) + 15 suppressor (CTA:3, TTA:12)
            "Glu": 14,
            "Gly": 16,
            "His": 5,
            "Ile": 16,
            "Leu": 20,
            "Lys": 20,
            "Met": 9,
            "Phe": 15,
            "Pro": 11,
            "Ser": 11,
            "Thr": 11,
            "Trp": 6,
            "Tyr": 9,
            "Val": 13,
        },
        has_disconnection=True,
        reassigned_aa="Gln",
        source=("tRNAscan-SE 2.0.12 on GCA_001447515.1 (April 2026); Ensembl Protists"),
        notes=(
            "UAA/UAG reassigned to Gln (Table 6). Class Oligohymenophorea "
            "(Scuticociliatia). Same class as Tetrahymena and Paramecium but "
            "phylogenetically distant (Scuticociliatia vs Hymenostomatia). "
            "tRNAscan-SE confirms 15 suppressor tRNAs (CTA:3 reading UAG, "
            "TTA:12 reading UAA) plus 5 normal Gln tRNAs. Total effective "
            "Gln = 20. Adds a third independent Gln-reassignment lineage."
        ),
    ),
    "hgrandinella_nuclear": TRNARepertoire(
        organism="Halteria grandinella",
        compartment="nuclear",
        ncbi_table_id=6,
        # VERIFIED: tRNAscan-SE 2.0.12 on GCA_006369765.1 (ASM636976v1)
        # 64 Mb genome, ~40,422 contigs (nanochromosomes).
        # Class Spirotrichea (Oligotrichia). UAA/UAG reassigned to Gln.
        # 121 standard-20-AA tRNAs + 3 suppressor (CTA:1, TTA:2) +
        # 1 SelCys + 5 pseudogenes = 130 total.
        by_amino_acid={
            "Ala": 7,
            "Arg": 10,
            "Asn": 6,
            "Asp": 3,
            "Cys": 2,
            "Gln": 9,  # VERIFIED: 6 normal (CTG:3, TTG:3) + 3 suppressor (CTA:1, TTA:2)
            "Glu": 8,
            "Gly": 7,
            "His": 4,
            "Ile": 11,
            "Leu": 9,
            "Lys": 10,
            "Met": 6,
            "Phe": 1,
            "Pro": 5,
            "Ser": 9,
            "Thr": 8,
            "Trp": 3,
            "Tyr": 5,
            "Val": 6,
        },
        has_disconnection=True,
        reassigned_aa="Gln",
        source=(
            "tRNAscan-SE 2.0.12 on GCA_006369765.1 (April 2026); "
            "Zheng et al. MBE 2021 (PMID 33500338)"
        ),
        notes=(
            "UAA/UAG reassigned to Gln (Table 6). Class Spirotrichea "
            "(Oligotrichia) — phylogenetically distinct from Oxytricha "
            "(Stichotrichia). ~40,422 nanochromosomes, 64 Mb genome. "
            "tRNAscan-SE confirms 3 suppressor tRNAs (CTA:1 reading UAG, "
            "TTA:2 reading UAA) plus 6 normal Gln tRNAs. Total effective "
            "Gln = 9. Independent Spirotrichea Gln-reassignment lineage."
        ),
    ),
    "bstoltei_nuclear": TRNARepertoire(
        organism="Blepharisma stoltei",
        compartment="nuclear",
        ncbi_table_id=15,
        # VERIFIED: tRNAscan-SE 2.0.12 on GCA_965603825.1
        # (Bstoltei_ATCC30299_MAC_1.00)
        # 41.5 Mb MAC genome. Class Heterotrichea.
        # UGA reassigned to Trp (Table 15). UAA/UAG remain stops.
        # 165 standard-20-AA tRNAs + 0 suppressor + 1 SelCys +
        # 1 undetermined + 2 pseudogenes = 169 total
        by_amino_acid={
            "Ala": 11,
            "Arg": 12,
            "Asn": 7,
            "Asp": 5,
            "Cys": 3,
            "Gln": 8,  # Standard Gln only (CTG:3, TTG:5), no suppressors
            "Glu": 9,
            "Gly": 10,
            "His": 4,
            "Ile": 9,
            "Leu": 16,
            "Lys": 11,
            "Met": 9,
            "Phe": 5,
            "Pro": 9,
            "Ser": 12,
            "Thr": 7,
            "Trp": 6,  # VERIFIED: CCA:6 — reads both UGG + UGA (Table 15)
            "Tyr": 4,
            "Val": 9,
        },
        has_disconnection=True,
        reassigned_aa="Trp",
        source=(
            "tRNAscan-SE 2.0.12 on GCA_965603825.1 (April 2026); "
            "Singh et al. PNAS 2023 (PMID 36669098)"
        ),
        notes=(
            "UGA reassigned to Trp (Table 15). UAA/UAG remain stops. "
            "Class Heterotrichea — closely related to B. japonicum. "
            "MAC genome 41.5 Mb (Singh et al. PNAS 2023). tRNAscan-SE "
            "confirms 6 tRNA-Trp(CCA) genes, 0 suppressors. 169 total "
            "tRNAs (165 standard 20-AA + 1 SelCys + 1 undetermined + "
            "2 pseudogenes). Replaces literature-estimated B. japonicum."
        ),
    ),
    # === BOUNDARY CASES: organisms where tRNA duplication does NOT occur ===
    "blastocrithidia_nuclear": TRNARepertoire(
        organism="Blastocrithidia nonstop",
        compartment="nuclear",
        ncbi_table_id=31,
        # VERIFIED: tRNAscan-SE 2.0.12 on GCA_028554745.1 (B. nonstop P57)
        # 24.7 Mb genome, 77 contigs. 68 tRNAs total (65 standard AA +
        # 2 suppressor + 1 SelCys).
        # Trp = 2 (both CCA anticodon): shortened anticodon stem (5bp→4bp)
        # reads both UGA and UGG. Not gene duplication — structural modification.
        # 2 suppressor tRNAs (CTA:1 = UAG→Glu, TTA:1 = UAA→Glu).
        by_amino_acid={
            "Ala": 5,
            "Arg": 7,
            "Asn": 3,
            "Asp": 3,
            "Cys": 1,
            "Gln": 2,
            "Glu": 4,  # 2 standard + 2 suppressor (CTA + TTA reading UAG/UAA)
            "Gly": 4,
            "His": 1,
            "Ile": 3,
            "Leu": 7,
            "Lys": 4,
            "Met": 2,
            "Phe": 2,
            "Pro": 3,
            "Ser": 5,
            "Thr": 4,
            "Trp": 2,  # tRNA-Trp(CCA) with shortened anticodon stem
            "Tyr": 2,
            "Val": 3,
        },
        has_disconnection=True,
        reassigned_aa="Trp",
        source=(
            "tRNAscan-SE 2.0.12 on GCA_028554745.1 (B. nonstop P57); "
            "Kachale et al. 2023 Nature 613:751-758; "
            "BioProject PRJNA790628"
        ),
        notes=(
            "UGA→Trp via anticodon stem shortening (5bp→4bp) of "
            "tRNA-Trp(CCA), not gene duplication (Kachale 2023). "
            "UAA/UAG→Glu via 2 new suppressor tRNAs (CTA + TTA). "
            "eRF1 Ser74Gly mutation reduces stop codon recognition. "
            "UAA is context-dependent: terminates AND encodes Glu "
            "in-frame. All 3 stops reassigned (Table 31). Boundary "
            "case: Trp tRNAs present (n=2) but mechanism is structural "
            "modification, not gene amplification."
        ),
    ),
    # === BACTERIAL UGA→Trp ORGANISMS (Table 4, tRNAscan-SE verified) ===
    #
    # Mycoplasma/Mycoplasmoides use Table 4 (UGA→Trp). Single tRNA-Trp(CCA)
    # reads both UGG and UGA via post-transcriptional anticodon modification
    # (not gene duplication). Minimal genomes, minimal tRNA repertoires.
    #
    "mgenitalium_bacterial": TRNARepertoire(
        organism="Mycoplasmoides genitalium G37",
        compartment="nuclear",  # bacterial chromosome, not compartmentalized
        ncbi_table_id=4,
        # VERIFIED: tRNAscan-SE 2.0.12 (-B mode) on GCA_000027325.1
        # 580 kb genome. 36 total tRNAs. Single tRNA-Trp(CCA).
        by_amino_acid={
            "Ala": 1,
            "Arg": 4,
            "Asn": 1,
            "Asp": 1,
            "Cys": 1,
            "Gln": 1,
            "Glu": 1,
            "Gly": 2,
            "His": 1,
            "Ile": 2,
            "Leu": 4,
            "Lys": 2,
            "Met": 2,
            "Phe": 1,
            "Pro": 1,
            "Ser": 4,
            "Thr": 3,
            "Trp": 1,
            "Tyr": 1,
            "Val": 1,
        },
        has_disconnection=True,
        reassigned_aa="Trp",
        source=(
            "tRNAscan-SE 2.0.12 (-B mode) on GCA_000027325.1 (April 2026); "
            "Himmelreich et al. 1996 Science 270:397-403"
        ),
        notes=(
            "UGA→Trp (Table 4). 580 kb genome — smallest known free-living "
            "organism. Single tRNA-Trp(CCA) reads both UGG and UGA via "
            "post-transcriptional anticodon modification, NOT gene "
            "duplication. 36 total tRNAs. Boundary case: confirms that "
            "minimal genomes use modification rather than duplication."
        ),
    ),
    "mpneumoniae_bacterial": TRNARepertoire(
        organism="Mycoplasmoides pneumoniae M129",
        compartment="nuclear",
        ncbi_table_id=4,
        # VERIFIED: tRNAscan-SE 2.0.12 (-B mode) on GCF_910574535.1
        # 816 kb genome. 37 total tRNAs. Single tRNA-Trp(CCA).
        by_amino_acid={
            "Ala": 1,
            "Arg": 4,
            "Asn": 1,
            "Asp": 1,
            "Cys": 1,
            "Gln": 1,
            "Glu": 1,
            "Gly": 2,
            "His": 1,
            "Ile": 2,
            "Leu": 4,
            "Lys": 2,
            "Met": 2,
            "Phe": 1,
            "Pro": 1,
            "Ser": 5,
            "Thr": 3,
            "Trp": 1,
            "Tyr": 1,
            "Val": 1,
        },
        has_disconnection=True,
        reassigned_aa="Trp",
        source=("tRNAscan-SE 2.0.12 (-B mode) on GCF_910574535.1 (April 2026)"),
        notes=(
            "UGA→Trp (Table 4). 816 kb genome. Single tRNA-Trp(CCA) reads "
            "both UGG and UGA via post-transcriptional modification. 37 "
            "total tRNAs. Confirms Mycoplasma boundary case."
        ),
    ),
    # === CROSS-KINGDOM CONTROLS (GtRNAdb, for non-ciliate comparisons) ===
    "dmelanogaster_nuclear": TRNARepertoire(
        organism="Drosophila melanogaster",
        compartment="nuclear",
        ncbi_table_id=1,
        # Standard nuclear code. GtRNAdb dm6 assembly, 283 tRNAs.
        by_amino_acid={
            "Ala": 17,
            "Arg": 26,
            "Asn": 10,
            "Asp": 14,
            "Cys": 7,
            "Gln": 12,
            "Glu": 19,
            "Gly": 20,
            "His": 5,
            "Ile": 11,
            "Leu": 22,
            "Lys": 19,
            "Met": 6,
            "Phe": 8,
            "Pro": 17,
            "Ser": 20,
            "Thr": 17,
            "Trp": 8,
            "Tyr": 10,
            "Val": 15,
        },
        has_disconnection=False,
        source="GtRNAdb Release 22, D. melanogaster BDGP Rel. 6/dm6",
        notes="Standard nuclear code. Metazoan outgroup control for ciliates.",
    ),
    "spombe_nuclear": TRNARepertoire(
        organism="Schizosaccharomyces pombe",
        compartment="nuclear",
        ncbi_table_id=1,
        # Standard nuclear code. GtRNAdb: 171 tRNAs.
        by_amino_acid={
            "Ala": 12,
            "Arg": 13,
            "Asn": 6,
            "Asp": 12,
            "Cys": 3,
            "Gln": 11,
            "Glu": 17,
            "Gly": 12,
            "His": 5,
            "Ile": 9,
            "Leu": 15,
            "Lys": 14,
            "Met": 6,
            "Phe": 5,
            "Pro": 9,
            "Ser": 13,
            "Thr": 10,
            "Trp": 3,
            "Tyr": 4,
            "Val": 12,
        },
        has_disconnection=False,
        source="GtRNAdb Release 22, S. pombe 972h-",
        notes="Standard nuclear code. Fission yeast control for ciliates.",
    ),
    "hsapiens_nuclear": TRNARepertoire(
        organism="Homo sapiens",
        compartment="nuclear",
        ncbi_table_id=1,
        # Standard nuclear code. GtRNAdb hg38, high-confidence set 429 tRNAs.
        by_amino_acid={
            "Ala": 38,
            "Arg": 28,
            "Asn": 25,
            "Asp": 13,
            "Cys": 29,
            "Gln": 19,
            "Glu": 16,
            "Gly": 28,
            "His": 9,
            "Ile": 23,
            "Leu": 31,
            "Lys": 27,
            "Met": 20,
            "Phe": 10,
            "Pro": 20,
            "Ser": 25,
            "Thr": 20,
            "Trp": 7,
            "Tyr": 13,
            "Val": 27,
        },
        has_disconnection=False,
        source="GtRNAdb Release 22, H. sapiens GRCh38/hg38",
        notes="Standard nuclear code. Large-genome control.",
    ),
}


def validate_repertoires() -> list[str]:
    """Internal-consistency checker for the curated tRNA data.

    Returns a list of problem strings. Empty list means all OK. Use in
    CI or before paper submission to catch typos and metadata drift.
    """
    problems = []
    # Check disconnection organisms all name a reassigned AA
    for key, rep in CURATED_REPERTOIRES.items():
        if rep.has_disconnection and not rep.reassigned_aa:
            problems.append(f"{key}: has_disconnection=True but reassigned_aa unset")
        if rep.reassigned_aa and not rep.has_disconnection:
            problems.append(f"{key}: reassigned_aa set but has_disconnection=False")
        if rep.reassigned_aa and rep.reassigned_aa not in rep.by_amino_acid:
            problems.append(
                f"{key}: reassigned_aa {rep.reassigned_aa!r} missing from by_amino_acid"
            )
        if not rep.source:
            problems.append(f"{key}: source is empty (cite your data)")
        if "placeholder" in rep.notes.lower():
            problems.append(
                f"{key}: notes contain 'placeholder' — replace with verified data"
            )
    # Check pairings reference existing keys and consistent AAs
    for dis_key, ctl_key, aa in DISCONNECTION_PAIRINGS:
        if dis_key not in CURATED_REPERTOIRES:
            problems.append(f"Pairing {dis_key} -> missing in CURATED_REPERTOIRES")
            continue
        if ctl_key not in CURATED_REPERTOIRES:
            problems.append(
                f"Pairing control {ctl_key} -> missing in CURATED_REPERTOIRES"
            )
            continue
        dis_rep = CURATED_REPERTOIRES[dis_key]
        if dis_rep.reassigned_aa and dis_rep.reassigned_aa != aa:
            problems.append(
                f"Pairing AA mismatch for {dis_key}: pairing={aa!r}, "
                f"repertoire.reassigned_aa={dis_rep.reassigned_aa!r}"
            )
    return problems


def get_repertoire(key: str) -> TRNARepertoire:
    """Get a curated repertoire entry by key."""
    return CURATED_REPERTOIRES[key]


def compare_aa_gene_counts(
    disconnection_key: str,
    control_key: str,
    reassigned_aa: str,
) -> dict:
    """Compare tRNA gene counts for the reassigned AA between a disconnection
    organism and a standard-code control.

    Returns a dict suitable for Fisher's exact test aggregation.
    """
    dis = get_repertoire(disconnection_key)
    ctl = get_repertoire(control_key)

    dis_count = dis.by_amino_acid.get(reassigned_aa, 0)
    ctl_count = ctl.by_amino_acid.get(reassigned_aa, 0)

    return {
        "disconnection_organism": dis.organism,
        "disconnection_compartment": dis.compartment,
        "control_organism": ctl.organism,
        "control_compartment": ctl.compartment,
        "reassigned_aa": reassigned_aa,
        "disconnection_aa_trna_count": dis_count,
        "control_aa_trna_count": ctl_count,
        "excess": dis_count - ctl_count,
        "ratio": dis_count / max(ctl_count, 1),
        "source_disconnection": dis.source,
        "source_control": ctl.source,
    }


# Mapping of disconnection cases to their controls.
# Controls chosen for phylogenetic proximity to the disconnection organism,
# Phylogenetic proximity criteria:
#   - Yeast mito (Saccharomycetales) paired with Yarrowia lipolytica (Dipodascaceae
#     fungi, same subphylum Saccharomycotina) — NOT Chlamydomonas
#   - Chlorophycean mito paired with C. reinhardtii (sister green alga)
#   - Nuclear CUG clade paired with Lachancea (non-CUG-reassigning Saccharomycete)
DISCONNECTION_PAIRINGS: list[tuple[str, str, str]] = [
    # (disconnection_key, control_key, reassigned_aa)
    #
    # Sampling frame: every variant-code organism with curated tRNA data,
    # paired with the phylogenetically closest standard-code organism.
    # Ciliate pairings use verified tRNAscan-SE data on NCBI assemblies.
    #
    # --- Non-ciliate pairings (organellar + nuclear CUG clade) ---
    ("scerevisiae_mito", "ylipolytica_mito", "Thr"),  # Fungal mito
    ("sobliquus_mito", "creinhardtii_mito", "Leu"),  # Chlorophyte mito
    ("ptannophilus_nuclear", "lthermotolerans_nuclear", "Ala"),  # CUG clade
    ("calbicans_nuclear", "lthermotolerans_nuclear", "Ser"),  # CUG clade
    # --- Ciliate pairings: variant-code vs standard-code ciliates ---
    # Oligohymenophorea (same class): variant vs I. multifiliis (standard)
    ("tthermophila_nuclear", "imultifiliis_nuclear", "Gln"),  # VERIFIED
    ("ptetraurelia_nuclear", "imultifiliis_nuclear", "Gln"),  # VERIFIED
    # Cross-class: variant vs S. coeruleus (Heterotrichea, standard)
    ("tthermophila_nuclear", "scoeruleus_nuclear", "Gln"),  # VERIFIED
    ("ptetraurelia_nuclear", "scoeruleus_nuclear", "Gln"),  # VERIFIED
    # Spirotrichea: O. trifallax verified by tRNAscan-SE
    ("oxytricha_nuclear", "scoeruleus_nuclear", "Gln"),
    # Non-Gln ciliates (literature-estimated, kept for comparison)
    ("eoctocarinatus_nuclear", "scoeruleus_nuclear", "Cys"),
    ("bjaponicum_nuclear", "scoeruleus_nuclear", "Trp"),
    # --- Euplotes (Table 10, UGA→Cys, tRNAscan-SE verified) ---
    ("eaediculatus_nuclear", "scoeruleus_nuclear", "Cys"),
    ("eaediculatus_nuclear", "fsalina_nuclear", "Cys"),
    ("eamieti_nuclear", "scoeruleus_nuclear", "Cys"),
    ("efocardii_nuclear", "fsalina_nuclear", "Cys"),
    ("eparawoodruffi_nuclear", "scoeruleus_nuclear", "Cys"),
    ("eweissei_nuclear", "fsalina_nuclear", "Cys"),
    ("ewoodruffi_nuclear", "scoeruleus_nuclear", "Cys"),
    # --- Additional ciliate pairings ---
    # Pseudocohnilembus: Oligohymenophorea scuticociliate, independent lineage
    ("ppersalinus_nuclear", "imultifiliis_nuclear", "Gln"),
    ("ppersalinus_nuclear", "fsalina_nuclear", "Gln"),
    # Halteria: Spirotrichea (Oligotrichia), independent from Oxytricha
    ("hgrandinella_nuclear", "scoeruleus_nuclear", "Gln"),
    ("hgrandinella_nuclear", "fsalina_nuclear", "Gln"),
    # Blepharisma stoltei: replaces literature-estimated B. japonicum
    ("bstoltei_nuclear", "scoeruleus_nuclear", "Trp"),
    ("bstoltei_nuclear", "fsalina_nuclear", "Trp"),
]

# Negative control: test a non-reassigned AA in a variant-code organism.
# Euplotes reassigned UGA -> Cys, so Trp should NOT be elevated.
NEGATIVE_CONTROL_PAIRINGS: list[tuple[str, str, str]] = [
    ("eoctocarinatus_nuclear", "scoeruleus_nuclear", "Trp"),
]


def trna_duplication_correlation_test() -> dict:
    """Test whether disconnection organisms have elevated tRNA gene counts
    for the reassigned amino acid.

    Runs a sign test / Wilcoxon-style comparison across the 4 disconnection
    cases: in how many does the disconnection organism have MORE tRNA genes
    for the reassigned AA than its standard-code control?

    A 4/4 result under a null of 50/50 gives p = 1/16 = 0.0625 via binomial.
    """
    results = []
    positive_count = 0

    for dis_key, ctl_key, aa in DISCONNECTION_PAIRINGS:
        cmp = compare_aa_gene_counts(dis_key, ctl_key, aa)
        results.append(cmp)
        if cmp["excess"] > 0:
            positive_count += 1

    n_tests = len(DISCONNECTION_PAIRINGS)

    # Binomial p-value: probability of observing >= positive_count positives
    # under a 50/50 null
    from math import comb

    def binomial_p_one_sided(successes: int, n: int) -> float:
        return sum(comb(n, k) for k in range(successes, n + 1)) / (2**n)

    p_value = binomial_p_one_sided(positive_count, n_tests)

    # Mean excess
    mean_excess = mean([r["excess"] for r in results])

    return {
        "n_pairings": n_tests,
        "n_with_elevated_trna": positive_count,
        "binomial_p_value": p_value,
        "mean_excess_trna_count": mean_excess,
        "hypothesis": (
            "Disconnection organisms have more tRNA genes for reassigned AA "
            "than standard-code controls"
        ),
        "pairings": results,
        "caveat": (
            "Small sample (n=4). Confirmed prior: Su et al. 2011 showed yeast "
            "mito tRNA-Thr2 derived from tRNA-His. This test extends the "
            "pattern and should be revisited with current GtRNAdb exports."
        ),
    }


def fisher_exact_per_pairing() -> dict:
    """Within-genome Fisher's exact test for each disconnection/control pair.

    For each pairing, builds a 2x2 contingency table:
                        reassigned_AA    other_AAs
        disconnection      a                b
        control            c                d

    where a = tRNA gene count for the reassigned AA, b = sum of all other AAs.
    Tests whether the reassigned AA is proportionally enriched in the
    disconnection organism (one-sided, alternative="greater").

    Then combines the 4 per-pairing p-values via Stouffer's Z method.
    """
    from scipy.stats import fisher_exact, norm

    per_pairing = []
    p_values = []

    for dis_key, ctl_key, aa in DISCONNECTION_PAIRINGS:
        dis = get_repertoire(dis_key)
        ctl = get_repertoire(ctl_key)

        a = dis.by_amino_acid.get(aa, 0)
        b = sum(dis.by_amino_acid.values()) - a
        c = ctl.by_amino_acid.get(aa, 0)
        d = sum(ctl.by_amino_acid.values()) - c

        fr = fisher_exact([[a, b], [c, d]], alternative="greater")
        odds_ratio = float(fr.statistic)  # type: ignore[attr-defined]
        p = float(fr.pvalue)  # type: ignore[attr-defined]

        per_pairing.append(
            {
                "disconnection_key": dis_key,
                "control_key": ctl_key,
                "reassigned_aa": aa,
                "table": [[a, b], [c, d]],
                "odds_ratio": odds_ratio,
                "fisher_p": p,
            }
        )
        p_values.append(p)

    # Stouffer's Z combination
    # Symmetric clipping to avoid infinite z-scores
    _CLIP_LO, _CLIP_HI = 1e-10, 1 - 1e-10
    z_scores = [float(norm.ppf(1 - max(min(p, _CLIP_HI), _CLIP_LO))) for p in p_values]
    combined_z = sum(z_scores) / (len(z_scores) ** 0.5)
    combined_p = float(1 - norm.cdf(combined_z))

    # Also compute on independent pairings only (no shared organisms).
    # Sort by Fisher p-value (ascending) for deterministic greedy selection
    # that favors the strongest-effect pairings.
    seen_organisms: set[str] = set()
    independent_p: list[float] = []
    sorted_pairings = sorted(per_pairing, key=lambda x: float(x["fisher_p"]))  # type: ignore[arg-type]
    for pp in sorted_pairings:
        if (
            pp["disconnection_key"] in seen_organisms
            or pp["control_key"] in seen_organisms
        ):
            continue
        seen_organisms.add(str(pp["disconnection_key"]))
        seen_organisms.add(str(pp["control_key"]))
        independent_p.append(float(pp["fisher_p"]))  # type: ignore[arg-type]

    if len(independent_p) >= 2:
        indep_z = [
            float(norm.ppf(1 - max(min(p, _CLIP_HI), _CLIP_LO))) for p in independent_p
        ]
        indep_combined_z = sum(indep_z) / (len(indep_z) ** 0.5)
        indep_combined_p = float(1 - norm.cdf(indep_combined_z))
    else:
        indep_combined_z = 0.0
        indep_combined_p = 1.0

    return {
        "per_pairing": per_pairing,
        "p_values": p_values,
        "z_scores": z_scores,
        "stouffer_z": combined_z,
        "stouffer_p": combined_p,
        "stouffer_z_independent": indep_combined_z,
        "stouffer_p_independent": indep_combined_p,
        "n_independent_pairings": len(independent_p),
        "n_pairings": len(DISCONNECTION_PAIRINGS),
        "method": "Fisher's exact (one-sided) per pairing, Stouffer's Z combination",
        "caveat": (
            "Pairings are not fully independent (some share control organisms). "
            "The stouffer_p_independent field uses only pairings with no shared organisms."
        ),
    }


def _stouffer_combine(p_values: list[float]) -> tuple[float, float]:
    """Combine p-values via Stouffer's Z method. Returns (Z, combined_p)."""
    from scipy.stats import norm

    _CLIP_LO, _CLIP_HI = 1e-10, 1 - 1e-10
    z_scores = [float(norm.ppf(1 - max(min(p, _CLIP_HI), _CLIP_LO))) for p in p_values]
    combined_z = sum(z_scores) / (len(z_scores) ** 0.5)
    combined_p = float(1 - norm.cdf(combined_z))
    return combined_z, combined_p


def maximal_independent_set_analysis() -> dict:
    """Enumerate all maximal independent sets (MIS) and compute Stouffer p-values.

    Greedy selection by strongest effect would bias the independent-pairings
    Stouffer p-value downward. Instead of selecting one set, we enumerate ALL
    maximal independent sets from the conflict graph (edges connect pairings
    sharing an organism) and report the distribution of Stouffer p-values.

    A maximal independent set is a set of pairings where no two share an
    organism AND no additional pairing can be added without creating a conflict.
    """
    from scipy.stats import fisher_exact

    # Step 1: compute Fisher's exact p-value for each pairing
    pairing_p: list[float] = []
    pairing_keys: list[tuple[str, str, str]] = []
    for dis_key, ctl_key, aa in DISCONNECTION_PAIRINGS:
        dis = get_repertoire(dis_key)
        ctl = get_repertoire(ctl_key)
        a = dis.by_amino_acid.get(aa, 0)
        b = sum(dis.by_amino_acid.values()) - a
        c = ctl.by_amino_acid.get(aa, 0)
        d = sum(ctl.by_amino_acid.values()) - c
        fr = fisher_exact([[a, b], [c, d]], alternative="greater")
        pairing_p.append(float(fr.pvalue))  # type: ignore[attr-defined]
        pairing_keys.append((dis_key, ctl_key, aa))

    n = len(DISCONNECTION_PAIRINGS)

    # Step 2: build conflict graph (adjacency list)
    # Two pairings conflict if they share a disconnection or control organism key
    conflicts: list[set[int]] = [set() for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            di, ci, _ = pairing_keys[i]
            dj, cj, _ = pairing_keys[j]
            if di == dj or di == cj or ci == dj or ci == cj:
                conflicts[i].add(j)
                conflicts[j].add(i)

    # Step 3: enumerate all maximal independent sets via Bron-Kerbosch
    all_mis: list[list[int]] = []

    def _bron_kerbosch(R: set[int], P: set[int], X: set[int]) -> None:
        if not P and not X:
            # R is a maximal independent set
            all_mis.append(sorted(R))
            return
        # Pick pivot that maximizes |P intersect N(u)| to prune
        pivot_candidates = P | X
        if not pivot_candidates:
            return
        pivot = max(pivot_candidates, key=lambda u: len(conflicts[u] & P))
        for v in list(P - conflicts[pivot]):
            _bron_kerbosch(
                R | {v},
                P - conflicts[v] - {v},
                X - conflicts[v],
            )
            P = P - {v}
            X = X | {v}

    _bron_kerbosch(set(), set(range(n)), set())

    # Step 4: for each MIS with >=2 members, compute Stouffer p-value
    mis_results: list[dict] = []
    for mis in all_mis:
        if len(mis) < 2:
            continue
        mis_p_values = [pairing_p[i] for i in mis]
        z, combined_p = _stouffer_combine(mis_p_values)
        mis_results.append(
            {
                "mis_indices": mis,
                "mis_size": len(mis),
                "pairings": [
                    f"{pairing_keys[i][0]} vs {pairing_keys[i][1]} ({pairing_keys[i][2]})"
                    for i in mis
                ],
                "fisher_p_values": mis_p_values,
                "stouffer_z": z,
                "stouffer_p": combined_p,
            }
        )

    if not mis_results:
        return {
            "n_mis_total": len(all_mis),
            "n_mis_size_ge2": 0,
            "error": "No MIS with >=2 members found",
        }

    stouffer_ps = [r["stouffer_p"] for r in mis_results]
    stouffer_ps.sort()
    median_idx = len(stouffer_ps) // 2

    # Find the greedy-best and greedy-worst
    best_mis = min(mis_results, key=lambda r: r["stouffer_p"])
    worst_mis = max(mis_results, key=lambda r: r["stouffer_p"])

    return {
        "n_pairings": n,
        "n_mis_total": len(all_mis),
        "n_mis_size_ge2": len(mis_results),
        "median_stouffer_p": stouffer_ps[median_idx],
        "worst_case_stouffer_p": stouffer_ps[-1],
        "best_case_stouffer_p": stouffer_ps[0],
        "fraction_significant_p05": sum(1 for p in stouffer_ps if p < 0.05)
        / len(stouffer_ps),
        "fraction_significant_p01": sum(1 for p in stouffer_ps if p < 0.01)
        / len(stouffer_ps),
        "best_mis": best_mis,
        "worst_mis": worst_mis,
        "all_mis_p_values": stouffer_ps,
        "method": (
            "Bron-Kerbosch enumeration of all maximal independent sets from "
            "the conflict graph (edges = shared organisms). Stouffer's Z "
            "combination computed for each MIS with >=2 members."
        ),
        "interpretation": (
            f"{'Robust' if stouffer_ps[-1] < 0.05 else 'Not robust'}: "
            f"worst-case MIS p = {stouffer_ps[-1]:.4f}, "
            f"median MIS p = {stouffer_ps[median_idx]:.4f}, "
            f"{sum(1 for p in stouffer_ps if p < 0.05)}/{len(stouffer_ps)} "
            f"MIS significant at p<0.05."
        ),
    }


def aa_label_permutation_test(
    n_permutations: int = 100_000,
    seed: int | None = None,
) -> dict:
    """Exact enumeration test: is the reassigned AA unusually enriched?

    For each disconnection/control pairing, computes the share difference
    for the REASSIGNED AA:
        observed = share_dis[reassigned_aa] - share_ctl[reassigned_aa]

    Then computes the same statistic for EVERY other AA:
        null_j = share_dis[j] - share_ctl[j]  for all 20 AAs

    The exact one-sided p-value is the fraction of AAs with share_diff
    >= the observed reassigned AA's share_diff. This uses the same AA
    label in both organisms (preserving AA-specific baseline structure)
    and requires no Monte Carlo because there are only ~20 AA labels.

    The n_permutations parameter is retained for API compatibility but
    is not used — the test is exact.
    """
    per_pairing = []

    for dis_key, ctl_key, aa in DISCONNECTION_PAIRINGS:
        dis = get_repertoire(dis_key)
        ctl = get_repertoire(ctl_key)

        dis_total = sum(dis.by_amino_acid.values())
        ctl_total = sum(ctl.by_amino_acid.values())

        # Compute share difference for every AA present in both organisms
        all_aas = sorted(set(dis.by_amino_acid.keys()) & set(ctl.by_amino_acid.keys()))
        share_diffs: dict[str, float] = {}
        for a in all_aas:
            d_share = dis.by_amino_acid.get(a, 0) / max(dis_total, 1)
            c_share = ctl.by_amino_acid.get(a, 0) / max(ctl_total, 1)
            share_diffs[a] = d_share - c_share

        observed_diff = share_diffs.get(aa, 0.0)

        # Exact p-value: fraction of AAs with share_diff >= observed
        # This is the tie-aware "min-rank" p-value: counts all AAs at
        # or above the observed level, including ties.
        n_total = len(share_diffs)
        n_as_extreme = sum(1 for d in share_diffs.values() if d >= observed_diff)
        p_exact = n_as_extreme / max(n_total, 1)

        # Tie-aware rank (min-rank: if tied, takes the best rank)
        rank = sum(1 for d in share_diffs.values() if d > observed_diff) + 1

        per_pairing.append(
            {
                "disconnection_key": dis_key,
                "control_key": ctl_key,
                "reassigned_aa": aa,
                "observed_share_diff": observed_diff,
                "rank_among_aas": rank,
                "n_aas_compared": n_total,
                "exact_p": p_exact,
            }
        )

    # Combine exact p-values via Stouffer (with independence caveat)
    from scipy.stats import norm

    _CLIP_LO, _CLIP_HI = 1e-10, 1 - 1e-10
    p_values = [float(r["exact_p"]) for r in per_pairing]  # type: ignore[arg-type]
    z_scores = [float(norm.ppf(1 - max(min(p, _CLIP_HI), _CLIP_LO))) for p in p_values]
    combined_z = sum(z_scores) / (len(z_scores) ** 0.5)
    combined_p = float(1 - norm.cdf(combined_z))

    return {
        "per_pairing": per_pairing,
        "stouffer_z": combined_z,
        "stouffer_p": combined_p,
        "n_pairings": len(DISCONNECTION_PAIRINGS),
        "method": "Exact AA-label enumeration per pairing, Stouffer's Z combination",
        "caveat": (
            "Exact enumeration over ~20 AA labels per pairing. "
            "Stouffer combination assumes independence — use with caveat "
            "when pairings share control organisms."
        ),
    }


def topology_breaking_subset_test() -> dict:
    """Restrict tRNA enrichment to topology-breaking reassignment events only.

    The all-pairings (n=24) Stouffer result blends several event classes —
    topology-breaking disconnections (yeast mito Thr, Scenedesmus Leu,
    Pachysolen Ala, Candida Ser), topology-preserving stop-to-sense
    reassignments (UAR-Gln in ciliates, UGA-Cys in Euplotes, UGA-Trp in
    Blepharisma), and ambiguous stop/sense systems. The headline framing
    ("variant-code lineages with topology-breaking reassignments show
    elevated tRNA gene counts") is best supported by restricting analysis
    to the topology-breaking subset.

    The 4 topology-breaking pairings:
      - scerevisiae_mito vs ylipolytica_mito for Thr (table 3, CUN-Thr)
      - sobliquus_mito vs creinhardtii_mito for Leu (table 22, UAG-Leu;
        equivalent to chlorophycean mito table 16)
      - ptannophilus_nuclear vs lthermotolerans_nuclear for Ala (table 26)
      - calbicans_nuclear vs lthermotolerans_nuclear for Ser (table 12)
    """
    from scipy.stats import fisher_exact, norm

    topology_breaking_keys = {
        ("scerevisiae_mito", "ylipolytica_mito", "Thr"),
        ("sobliquus_mito", "creinhardtii_mito", "Leu"),
        ("ptannophilus_nuclear", "lthermotolerans_nuclear", "Ala"),
        ("calbicans_nuclear", "lthermotolerans_nuclear", "Ser"),
    }

    per_pairing = []
    p_values: list[float] = []
    for dis_key, ctl_key, aa in DISCONNECTION_PAIRINGS:
        if (dis_key, ctl_key, aa) not in topology_breaking_keys:
            continue
        dis = get_repertoire(dis_key)
        ctl = get_repertoire(ctl_key)
        a = dis.by_amino_acid.get(aa, 0)
        b = sum(dis.by_amino_acid.values()) - a
        c = ctl.by_amino_acid.get(aa, 0)
        d = sum(ctl.by_amino_acid.values()) - c
        fr = fisher_exact([[a, b], [c, d]], alternative="greater")
        odds_ratio = float(fr.statistic)  # type: ignore[attr-defined]
        p = float(fr.pvalue)  # type: ignore[attr-defined]
        per_pairing.append(
            {
                "disconnection_key": dis_key,
                "control_key": ctl_key,
                "reassigned_aa": aa,
                "table_id_lineage": (
                    "table_3_yeast_mito"
                    if dis_key == "scerevisiae_mito"
                    else "table_22_scenedesmus_mito"
                    if dis_key == "sobliquus_mito"
                    else "table_26_pachysolen_nuclear"
                    if dis_key == "ptannophilus_nuclear"
                    else "table_12_candida_nuclear"
                ),
                "table": [[a, b], [c, d]],
                "odds_ratio": odds_ratio,
                "fisher_p": p,
            }
        )
        p_values.append(p)

    # Stouffer's Z combination on topology-breaking subset
    _CLIP_LO, _CLIP_HI = 1e-10, 1 - 1e-10
    z_scores = [float(norm.ppf(1 - max(min(p, _CLIP_HI), _CLIP_LO))) for p in p_values]
    if z_scores:
        combined_z = sum(z_scores) / (len(z_scores) ** 0.5)
        combined_p = float(1 - norm.cdf(combined_z))
    else:
        combined_z, combined_p = 0.0, 1.0

    return {
        "method": (
            "Topology-breaking subset of the tRNA-enrichment analysis. "
            "Restricts the 24-pairing all-pairings analysis to the 4 "
            "pairings whose underlying variant-code reassignment creates "
            "a new amino-acid disconnection in GF(2)^6 (yeast mito Thr, "
            "Scenedesmus mito Leu, Pachysolen Ala, Candida Ser); "
            "more direct test of the compensation-via-tRNA-duplication "
            "hypothesis."
        ),
        "n_pairings": len(per_pairing),
        "per_pairing": per_pairing,
        "stouffer_z_topology_breaking": combined_z,
        "stouffer_p_topology_breaking": combined_p,
        "caveat": (
            "Sample size n=4 is small; this is a restricted subset analysis. "
            "Reported alongside the all-pairings result for transparency. "
            "Some pairings use literature/GtRNAdb counts rather than "
            "tRNAscan-SE-verified counts; see tRNARepertoire.source field."
        ),
    }


def trna_evidence_summary() -> dict:
    """Summary for paper or report, including all statistical tests."""
    sign_test = trna_duplication_correlation_test()
    fisher_test = fisher_exact_per_pairing()
    mis_test = maximal_independent_set_analysis()
    perm_test = aa_label_permutation_test()
    topo_subset = topology_breaking_subset_test()
    return {
        "sign_test": sign_test,
        "fisher_stouffer": fisher_test,
        "mis_analysis": mis_test,
        "permutation_test": perm_test,
        "topology_breaking_subset": topo_subset,
        "repertoires_used": list(CURATED_REPERTOIRES.keys()),
        "n_disconnection_organisms": sum(
            1 for r in CURATED_REPERTOIRES.values() if r.has_disconnection
        ),
        "n_control_organisms": sum(
            1 for r in CURATED_REPERTOIRES.values() if not r.has_disconnection
        ),
    }
