#!/usr/bin/env python3.11
"""Emit the 24-row tRNA analysis-input provenance table.

For each of the 24 (disconnection, control, reassigned-AA) pairings in
DISCONNECTION_PAIRINGS, print columns:

  - variant organism / control organism
  - compartment
  - assembly accession (or "literature" / "annotation" / "GtRNAdb")
  - NCBI translation table id + reassignment (source AA -> target AA at
    the reassigned codon)
  - topology status under Q_6 (creates new AA-family disconnection?) and
    H(3,4) (same question, encoding-independent)
  - Fisher 2x2 focal count and total denominator for the variant and
    control
  - source label per organism (tRNAscan-SE 1st-pass / 2nd-pass /
    GtRNAdb / literature / annotation)
  - whether the organism is in the 18-assembly tRNAscan-SE-verified set
  - per-pair one-sided Fisher-exact p-value (target = greater)

Writes:
  output/tables/T-trna-24row-provenance.csv  (source of truth)
  output/manuscript_stats.trna_provenance     (24-row list; consumed by
                                               supplement §S10.1 table)
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

from scipy.stats import fisher_exact

from codon_topo.analysis.trna_evidence import (
    DISCONNECTION_PAIRINGS,
    get_repertoire,
)
from codon_topo.core.encoding import nucleotide_distance
from codon_topo.core.genetic_codes import STANDARD, get_code


OUT_ROOT = Path("output")
CSV_PATH = OUT_ROOT / "tables" / "T-trna-24row-provenance.csv"

# Table IDs for the 4 topology-breaking (Q6) pairings whose underlying
# reassignment creates a new AA-family disconnection at epsilon=1.
_Q6_CREATES_KEYS: dict[tuple[str, str, str], int] = {
    ("scerevisiae_mito", "ylipolytica_mito", "Thr"): 3,
    ("sobliquus_mito", "creinhardtii_mito", "Leu"): 22,
    ("ptannophilus_nuclear", "lthermotolerans_nuclear", "Ala"): 26,
    ("calbicans_nuclear", "lthermotolerans_nuclear", "Ser"): 12,
}


def _partition_components_h34(codons: list[str]) -> int:
    n = len(codons)
    if n == 0:
        return 0
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x: int, y: int) -> None:
        px, py = find(x), find(y)
        if px != py:
            parent[px] = py

    for i in range(n):
        for j in range(i + 1, n):
            if nucleotide_distance(codons[i], codons[j]) == 1:
                union(i, j)
    return len({find(i) for i in range(n)})


def _h34_impact(table_id: int, aa: str) -> str:
    variant = get_code(table_id)
    std = [c for c, a in STANDARD.items() if a == aa]
    var = [c for c, a in variant.items() if a == aa]
    n_std = _partition_components_h34(std)
    n_var = _partition_components_h34(var)
    if n_std == 1 and n_var > 1:
        return "creates"
    if n_std >= 2 and n_var > n_std:
        return "extends"
    if n_var < n_std:
        return "reduces"
    return "no-effect"


def _source_class(source_str: str) -> str:
    s = source_str.lower()
    if "trnascan-se" in s:
        return "tRNAscan-SE 2.0.12"
    if "gtrnadb" in s:
        return "GtRNAdb"
    if "annotation" in s or "refseq" in s or "candida genome database" in s:
        return "annotation"
    return "literature"


def _accession(source_str: str) -> str:
    for tok in source_str.replace("(", " ").replace(")", " ").split():
        if tok.startswith(("GCF_", "GCA_")):
            return tok.rstrip(";.")
    return ""


def _verified(source_str: str) -> str:
    return "yes" if _source_class(source_str) == "tRNAscan-SE 2.0.12" else "no"


def _reassignment_string(table_id: int, aa: str) -> str:
    if table_id == 1:
        return "(standard-code control)"
    variant = get_code(table_id)
    diffs = [
        f"{c}: {STANDARD[c]}→{variant[c]}"
        for c in sorted(variant)
        if STANDARD[c] != variant[c] and variant[c] == aa
    ]
    return "; ".join(diffs) if diffs else "(no direct reassignment for this AA)"


def _fisher_row(dis_key: str, ctl_key: str, aa: str) -> dict:
    dis = get_repertoire(dis_key)
    ctl = get_repertoire(ctl_key)
    a = dis.by_amino_acid.get(aa, 0)
    b = sum(dis.by_amino_acid.values()) - a
    c = ctl.by_amino_acid.get(aa, 0)
    d = sum(ctl.by_amino_acid.values()) - c
    fr = fisher_exact([[a, b], [c, d]], alternative="greater")
    return {
        "focal_variant": a,
        "denom_variant": a + b,
        "focal_control": c,
        "denom_control": c + d,
        "odds_ratio": float(fr.statistic),  # type: ignore[attr-defined]
        "p_greater": float(fr.pvalue),  # type: ignore[attr-defined]
    }


def build_rows() -> list[dict]:
    def _source_display(
        key: str, source_class: str, source_full: str, accession: str
    ) -> str:
        """Precompute the accession-or-citation string for Table S11b so the
        Typst side is a plain field lookup (no multi-branch #{ if...else }
        block that freeze_typst.py's block-splice resolver cannot evaluate)."""
        if accession:
            return accession
        if source_class == "literature" and key.startswith("scerevisiae"):
            return "Bonitz et al. 1980"
        if source_class == "literature" and key.startswith("ylipolytica"):
            return "Kerscher et al. 2001"
        return source_full

    rows: list[dict] = []
    for dis_key, ctl_key, aa in DISCONNECTION_PAIRINGS:
        dis = get_repertoire(dis_key)
        ctl = get_repertoire(ctl_key)
        table_id = dis.ncbi_table_id
        q6_creates = (dis_key, ctl_key, aa) in _Q6_CREATES_KEYS
        h34 = _h34_impact(table_id, aa) if q6_creates else "no-effect"
        fisher = _fisher_row(dis_key, ctl_key, aa)
        reass = _reassignment_string(table_id, aa)
        rows.append(
            {
                "variant_key": dis_key,
                "variant_short": dis_key.split("_")[0],
                "variant_organism": dis.organism,
                "variant_compartment": dis.compartment,
                "variant_table_id": dis.ncbi_table_id,
                "variant_source_class": _source_class(dis.source),
                "variant_accession": _accession(dis.source),
                "variant_source_full": dis.source,
                "variant_source_display": _source_display(
                    dis_key,
                    _source_class(dis.source),
                    dis.source,
                    _accession(dis.source),
                ),
                "variant_verified": _verified(dis.source),
                "control_key": ctl_key,
                "control_short": ctl_key.split("_")[0],
                "control_organism": ctl.organism,
                "control_compartment": ctl.compartment,
                "compartment_marker": (
                    dis.compartment[0].upper()
                    if dis.compartment == ctl.compartment
                    else f"{dis.compartment[0].upper()}/{ctl.compartment[0].upper()}"
                ),
                "control_table_id": ctl.ncbi_table_id,
                "control_source_class": _source_class(ctl.source),
                "control_accession": _accession(ctl.source),
                "control_source_full": ctl.source,
                "control_source_display": _source_display(
                    ctl_key,
                    _source_class(ctl.source),
                    ctl.source,
                    _accession(ctl.source),
                ),
                "control_verified": _verified(ctl.source),
                "reassigned_aa": aa,
                "reassignment": reass,
                "reassignment_short": reass if len(reass) <= 30 else reass[:29] + "…",
                "q6_creates_disconnection": q6_creates,
                "q6_break_marker": "Y" if q6_creates else "n",
                "h34_impact": h34,
                **fisher,
            }
        )
    return rows


def main() -> int:
    rows = build_rows()
    CSV_PATH.parent.mkdir(parents=True, exist_ok=True)
    with CSV_PATH.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    stats_path = OUT_ROOT / "manuscript_stats.json"
    if stats_path.exists():
        stats = json.loads(stats_path.read_text())
        stats.setdefault("trna", {})["provenance_24row"] = rows
        n_organisms_verified = len(
            {r["variant_key"] for r in rows if r["variant_verified"] == "yes"}
            | {r["control_key"] for r in rows if r["control_verified"] == "yes"}
        )
        stats["trna"]["n_provenance_pairings"] = len(rows)
        stats["trna"]["n_verified_organisms"] = n_organisms_verified
        stats_path.write_text(json.dumps(stats, indent=2, default=str, allow_nan=False))
    counts_by_class: dict[str, int] = {}
    for r in rows:
        for k in ("variant_source_class", "control_source_class"):
            counts_by_class[r[k]] = counts_by_class.get(r[k], 0) + 1
    print(
        f"emitted {len(rows)} pairings to {CSV_PATH}; "
        f"organism-source-class occurrences: {counts_by_class}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
