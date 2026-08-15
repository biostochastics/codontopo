"""GenBank comparison utilities for extracting codon recoding positions.

Compares a designed/recoded genome against a parent genome to identify
every codon that was changed, along with its gene context, position,
and strand information.

Uses Biopython for robust GenBank parsing including complement strands,
join locations, and CDS qualifier extraction.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


@dataclass(frozen=True)
class CodonChange:
    """A single codon difference between parent and designed/final genome.

    Parameters
    ----------
    genome_pos_nt : int
        0-based position of the first nucleotide of the codon in the genome.
    source_codon_dna : str
        Parent codon (DNA alphabet, uppercase).
    target_codon_dna : str
        Designed/final codon (DNA alphabet, uppercase).
    gene : str or None
        Gene name from CDS qualifier.
    locus_tag : str or None
        Locus tag from CDS qualifier.
    strand : int
        1 for plus strand, -1 for minus strand.
    codon_index_in_cds : int
        0-based index of the codon within the CDS.
    is_essential : bool or None
        Whether the gene is essential (if known).
    source_aa : str
        Amino acid encoded by source codon.
    target_aa : str
        Amino acid encoded by target codon.
    is_synonymous : bool
        Whether source and target encode the same amino acid.
    """

    genome_pos_nt: int
    source_codon_dna: str
    target_codon_dna: str
    gene: str | None
    locus_tag: str | None
    strand: int
    codon_index_in_cds: int
    is_essential: bool | None
    source_aa: str
    target_aa: str
    is_synonymous: bool


def load_genbank(path: Path | str) -> SeqRecord:
    """Load a GenBank file, returning the first record."""
    return next(SeqIO.parse(str(path), "genbank"))


def _extract_codons_from_cds(
    record: SeqRecord,
    feature,
) -> list[tuple[int, str, int]]:
    """Extract (genome_pos, codon_dna, codon_index) for each codon in a CDS.

    Handles both simple and compound (join) locations, plus and minus strands.

    Returns list of (genome_pos_nt, codon_dna, codon_index_in_cds) tuples.
    genome_pos_nt is the 0-based position of the first nucleotide of the
    codon on the PLUS strand of the genome (even for minus-strand genes).
    """
    # Honor /codon_start qualifier (1-based; default is 1)
    codon_start = int(feature.qualifiers.get("codon_start", ["1"])[0]) - 1

    # Extract the CDS nucleotide sequence
    cds_seq = str(feature.extract(record.seq)).upper()
    if codon_start > 0:
        cds_seq = cds_seq[codon_start:]
    cds_len = len(cds_seq)

    # Get the genome positions that this CDS covers
    # For compound locations (joins), this handles all parts
    # Skip codon_start bases from the position list to match the trimmed sequence
    positions = []
    _all_positions: list[int] = []
    for part in feature.location.parts:
        start = int(part.start)
        end = int(part.end)
        strand = part.strand if part.strand is not None else 1
        if strand >= 0:
            _all_positions.extend(range(start, end))
        else:
            _all_positions.extend(range(end - 1, start - 1, -1))

    # Trim codon_start positions from the front
    positions = _all_positions[codon_start:]

    if len(positions) != cds_len:
        return []

    # Group into codons
    result = []
    for i in range(0, cds_len - 2, 3):
        codon = cds_seq[i : i + 3]
        if len(codon) != 3:
            break
        # genome_pos is the position of the first base of this codon
        # on the plus strand
        genome_pos = min(positions[i], positions[i + 1], positions[i + 2])
        codon_index = i // 3
        result.append((genome_pos, codon, codon_index))

    return result


def compare_genomes(
    parent: SeqRecord,
    design: SeqRecord,
    recode_map: dict[str, str] | None = None,
) -> list[CodonChange]:
    """Compare two GenBank records to find all codon differences in CDS regions.

    Parameters
    ----------
    parent : SeqRecord
        The parent/wild-type genome.
    design : SeqRecord
        The designed/recoded genome.
    recode_map : dict or None
        If provided, only report changes matching these source->target pairs
        (DNA alphabet, uppercase). If None, report all codon differences.

    Returns
    -------
    list of CodonChange
        Every codon that differs between parent and design within annotated CDS.
    """
    from Bio.Data.CodonTable import standard_dna_table

    codon_table = standard_dna_table.forward_table
    stop_codons = set(standard_dna_table.stop_codons)

    def translate_codon(codon: str) -> str:
        if codon in stop_codons:
            return "Stop"
        return codon_table.get(codon, "?")

    # Build CDS index from the parent genome
    # Match CDS features between parent and design by locus_tag or gene name
    parent_cds = {}
    for f in parent.features:
        if f.type != "CDS":
            continue
        quals = f.qualifiers
        key = quals.get("locus_tag", [None])[0] or quals.get("gene", [None])[0]
        if key:
            parent_cds[key] = f

    changes = []
    for f in design.features:
        if f.type != "CDS":
            continue
        quals = f.qualifiers
        key = quals.get("locus_tag", [None])[0] or quals.get("gene", [None])[0]
        if not key or key not in parent_cds:
            continue

        gene = quals.get("gene", [None])[0]
        locus_tag = quals.get("locus_tag", [None])[0]
        strand = f.location.strand if f.location.strand is not None else 1

        parent_codons = _extract_codons_from_cds(parent, parent_cds[key])
        design_codons = _extract_codons_from_cds(design, f)

        if len(parent_codons) != len(design_codons):
            continue

        for (p_pos, p_codon, p_idx), (d_pos, d_codon, d_idx) in zip(
            parent_codons, design_codons
        ):
            if p_codon == d_codon:
                continue

            if recode_map is not None:
                if p_codon not in recode_map:
                    continue
                if recode_map[p_codon] != d_codon:
                    continue

            p_aa = translate_codon(p_codon)
            d_aa = translate_codon(d_codon)

            changes.append(
                CodonChange(
                    genome_pos_nt=p_pos,
                    source_codon_dna=p_codon,
                    target_codon_dna=d_codon,
                    gene=gene,
                    locus_tag=locus_tag,
                    strand=strand,
                    codon_index_in_cds=p_idx,
                    is_essential=None,
                    source_aa=p_aa,
                    target_aa=d_aa,
                    is_synonymous=(p_aa == d_aa),
                )
            )

    return changes


def compare_genome_by_sequence(
    parent: SeqRecord,
    design: SeqRecord,
    recode_map: dict[str, str],
) -> list[CodonChange]:
    """Compare genomes by scanning the full sequence for recoded codons.

    Faster alternative to CDS-based comparison when the genomes have
    the same length and coordinates are directly comparable. Scans every
    position where the parent has a source codon from recode_map and
    the design has the corresponding target.

    Then annotates each hit with CDS context from the parent features.
    """
    from Bio.Data.CodonTable import standard_dna_table

    codon_table = standard_dna_table.forward_table
    stop_codons = set(standard_dna_table.stop_codons)

    def translate_codon(codon: str) -> str:
        if codon in stop_codons:
            return "Stop"
        return codon_table.get(codon, "?")

    parent_seq = str(parent.seq).upper()
    design_seq = str(design.seq).upper()

    if len(parent_seq) != len(design_seq):
        raise ValueError(
            f"Genome lengths differ: parent={len(parent_seq)}, "
            f"design={len(design_seq)}. Use compare_genomes() instead."
        )

    # Build a position-to-CDS index for annotation
    pos_to_cds: dict[int, dict] = {}
    for f in parent.features:
        if f.type != "CDS":
            continue
        quals = f.qualifiers
        gene = quals.get("gene", [None])[0]
        locus_tag = quals.get("locus_tag", [None])[0]
        strand = f.location.strand if f.location.strand is not None else 1
        codons = _extract_codons_from_cds(parent, f)
        for genome_pos, codon, codon_idx in codons:
            pos_to_cds[genome_pos] = {
                "gene": gene,
                "locus_tag": locus_tag,
                "strand": strand,
                "codon_index": codon_idx,
            }

    # Scan for recoded positions
    changes = []
    for pos in range(len(parent_seq) - 2):
        p_codon = parent_seq[pos : pos + 3]
        d_codon = design_seq[pos : pos + 3]

        if p_codon == d_codon:
            continue
        if p_codon not in recode_map:
            continue
        if recode_map[p_codon] != d_codon:
            continue

        cds_info = pos_to_cds.get(pos, {})
        p_aa = translate_codon(p_codon)
        d_aa = translate_codon(d_codon)

        changes.append(
            CodonChange(
                genome_pos_nt=pos,
                source_codon_dna=p_codon,
                target_codon_dna=d_codon,
                gene=cds_info.get("gene"),
                locus_tag=cds_info.get("locus_tag"),
                strand=cds_info.get("strand", 1),
                codon_index_in_cds=cds_info.get("codon_index", -1),
                is_essential=None,
                source_aa=p_aa,
                target_aa=d_aa,
                is_synonymous=(p_aa == d_aa),
            )
        )

    return changes


def changes_to_events(
    changes: list[CodonChange],
    study_id: str,
    table_id: int = 11,
    organism: str = "E_coli",
    strain: str | None = None,
    fragment_boundaries: dict[str, tuple[int, int]] | None = None,
) -> list:
    """Convert CodonChange list to CodonSwapEvent + RecodingOutcome pairs.

    Parameters
    ----------
    fragment_boundaries : dict or None
        Mapping of fragment_name -> (start_pos, end_pos) for assigning unit_ids.
    """
    from codon_topo.analysis.codonsafe.models import (
        CodonSwapEvent,
        StudyId,
    )
    from codon_topo.analysis.codonsafe.normalize import normalize_codon

    study = StudyId(study_id)

    results = []
    for ch in changes:
        try:
            source_rna = normalize_codon(ch.source_codon_dna)
            target_rna = normalize_codon(ch.target_codon_dna)
        except ValueError:
            continue

        unit_id = None
        if fragment_boundaries:
            for frag_name, (fstart, fend) in fragment_boundaries.items():
                if fstart <= ch.genome_pos_nt < fend:
                    unit_id = frag_name
                    break

        gene_str = ch.gene or ch.locus_tag or "unknown"
        event = CodonSwapEvent(
            study=study,
            event_id=f"{study_id}_{gene_str}_{ch.codon_index_in_cds}_{target_rna}",
            source_codon=source_rna,
            target_codon=target_rna,
            table_id=table_id,
            organism=organism,
            strain=strain,
            gene=ch.gene,
            locus_tag=ch.locus_tag,
            genome_pos_nt=ch.genome_pos_nt,
            codon_index_in_cds=ch.codon_index_in_cds,
            is_essential_gene=ch.is_essential,
            unit_id=unit_id,
        )
        results.append((event, None))

    return results
