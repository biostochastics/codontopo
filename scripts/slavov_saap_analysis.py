# pyright: reportAttributeAccessIssue=false, reportArgumentType=false
# pandas Series / scipy.stats tuple-return types are false positives here.
"""Tsour/Slavov (2026) SAAP cross-analysis: minimum source-target codon
Hamming distance for each high-confidence amino-acid substitution event.

This is the analysis behind §4.6 of the manuscript and §S17 of the
supplement. The input is the 'Supplementary Data 8' file from Tsour et
al. (2026, Nature, doi:10.1038/s41586-026-10678-2), which contains 5,873
high-confidence SAAP-BP pairs at positional probability > 0.9. Download
it from the Slavov lab decode portal:

  https://decode.slavovlab.net/mass-spec/data
  → Supplementary_Data_8.High_confidence_SAAP_precursor_quant.xlsx

Usage:
  python scripts/slavov_saap_analysis.py \\
      --input /path/to/Supplementary_Data_8.High_confidence_SAAP_precursor_quant.xlsx \\
      --output output/slavov_saap_codon_distances.json
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from itertools import product
from pathlib import Path

from codon_topo.core.encoding import (
    DEFAULT_ENCODING,
    codon_to_vector,
    hamming_distance,
    nucleotide_distance,
)
from codon_topo.core.genetic_codes import STANDARD

_AA3_TO_1 = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "Stop": "*",
}

_STANDARD_AAS = "ACDEFGHIKLMNPQRSTVWY"


def _build_pair_distances() -> dict[tuple[str, str], tuple[int, int]]:
    codons_by_aa: dict[str, list[str]] = defaultdict(list)
    for codon, aa3 in STANDARD.items():
        codons_by_aa[_AA3_TO_1[aa3]].append(codon)

    pair_dists: dict[tuple[str, str], tuple[int, int]] = {}
    for src in _STANDARD_AAS:
        for tgt in _STANDARD_AAS:
            if src == tgt:
                continue
            pairs = list(product(codons_by_aa[src], codons_by_aa[tgt]))
            min_nt = min(nucleotide_distance(c1, c2) for c1, c2 in pairs)
            min_bit = min(
                hamming_distance(
                    codon_to_vector(c1, DEFAULT_ENCODING),
                    codon_to_vector(c2, DEFAULT_ENCODING),
                )
                for c1, c2 in pairs
            )
            pair_dists[(src, tgt)] = (min_nt, min_bit)
    return pair_dists


def _parse_aas(s: str) -> tuple[str, str] | None:
    parts = s.split(" to ")
    if len(parts) != 2:
        return None
    src, tgt = parts[0].strip(), parts[1].strip()
    if "/" in src or "/" in tgt or len(src) != 1 or len(tgt) != 1:
        return None
    if src not in _STANDARD_AAS or tgt not in _STANDARD_AAS:
        return None
    return (src, tgt)


def run_analysis(input_xlsx: Path, output_json: Path) -> dict:
    try:
        import pandas as pd
    except ImportError as e:  # pragma: no cover
        raise SystemExit(
            "pandas is required for the SAAP analysis: pip install pandas openpyxl"
        ) from e

    from scipy.stats import binomtest, fisher_exact, mannwhitneyu  # type: ignore
    import numpy as np

    df = pd.read_excel(input_xlsx, sheet_name="Sheet1")
    pair_dists = _build_pair_distances()

    df["parsed"] = df["AAS"].apply(_parse_aas)
    df_clean = df[df["parsed"].notnull()].copy()
    df_clean["d_nt"] = df_clean["parsed"].apply(lambda p: pair_dists[p][0])
    df_clean["d_bit"] = df_clean["parsed"].apply(lambda p: pair_dists[p][1])

    unique_obs = set(df_clean["parsed"].unique())
    unobs_pairs = [p for p in pair_dists if p not in unique_obs]
    obs_dnt = [pair_dists[p][0] for p in unique_obs]
    unobs_dnt = [pair_dists[p][0] for p in unobs_pairs]
    obs_dbit = [pair_dists[p][1] for p in unique_obs]
    unobs_dbit = [pair_dists[p][1] for p in unobs_pairs]

    u_nt, p_nt = mannwhitneyu(obs_dnt, unobs_dnt, alternative="less")
    u_bit, p_bit = mannwhitneyu(obs_dbit, unobs_dbit, alternative="less")

    obs_1nt = sum(1 for d in obs_dnt if d == 1)
    obs_other = sum(1 for d in obs_dnt if d != 1)
    unobs_1nt = sum(1 for d in unobs_dnt if d == 1)
    unobs_other = sum(1 for d in unobs_dnt if d != 1)
    odds, p_fisher = fisher_exact(
        [[obs_1nt, obs_other], [unobs_1nt, unobs_other]], alternative="greater"
    )

    events_1nt = int((df_clean["d_nt"] == 1).sum())
    events_total = len(df_clean)
    all_pair_1nt_frac = sum(1 for v in pair_dists.values() if v[0] == 1) / len(
        pair_dists
    )
    bt = binomtest(events_1nt, events_total, all_pair_1nt_frac, alternative="greater")

    # Full per-distance breakdown for the §S17 grouped-bar figure.
    # Distances observed in the standard code: 1, 2, 3 (max min-NT distance
    # between codons of two distinct amino acids is 3).
    distances = sorted({v[0] for v in pair_dists.values()})
    event_distance_distribution = {}
    baseline_distance_distribution = {}
    for d in distances:
        n_events_d = int((df_clean["d_nt"] == d).sum())
        n_pairs_d = int(sum(1 for v in pair_dists.values() if v[0] == d))
        event_distance_distribution[str(d)] = {
            "count": n_events_d,
            "fraction": round(float(n_events_d / events_total), 4),
        }
        baseline_distance_distribution[str(d)] = {
            "count": n_pairs_d,
            "fraction": round(float(n_pairs_d / len(pair_dists)), 4),
        }

    result = {
        "n_high_confidence_saap_events": int(len(df)),
        "n_parseable_events": int(events_total),
        "n_unique_observed_aas_pairs": int(df_clean["parsed"].nunique()),
        "event_single_nt_count": events_1nt,
        "event_single_nt_fraction": round(float(events_1nt / events_total), 4),
        "baseline_single_nt_fraction_all_380_pairs": round(float(all_pair_1nt_frac), 4),
        "event_enrichment_binomial_p_one_sided": float(bt.pvalue),
        "event_distance_distribution": event_distance_distribution,
        "baseline_distance_distribution": baseline_distance_distribution,
        "n_all_pairs_baseline": int(len(pair_dists)),
        "unique_pairs_fisher_odds_ratio": float(odds),
        "unique_pairs_fisher_p_one_sided": float(p_fisher),
        "mannwhitney_nt_observed_smaller": {
            "U": float(u_nt),
            "p_one_sided": float(p_nt),
            "median_obs": float(np.median(obs_dnt)),
            "median_unobs": float(np.median(unobs_dnt)),
        },
        "mannwhitney_bit_observed_smaller": {
            "U": float(u_bit),
            "p_one_sided": float(p_bit),
            "median_obs": float(np.median(obs_dbit)),
            "median_unobs": float(np.median(unobs_dbit)),
        },
        "data_source": (
            "Tsour et al. 2026 Nature, Supplementary Data 8 "
            "(high-confidence SAAP, positional probability > 0.9). "
            "doi:10.1038/s41586-026-10678-2"
        ),
        "methodology": (
            "For each amino-acid substitution observed in the Tsour & "
            "Slavov 2026 proteomic dataset, compute the minimum nucleotide-"
            "Hamming distance between any source codon (in the standard "
            "genetic code) and any target codon. Enrichment tests: "
            "(1) event-weighted binomial against all-pair baseline; "
            "(2) Fisher exact on single-NT pairs unique observed vs unique "
            "unobserved; (3) Mann-Whitney U on NT and bit distances."
        ),
    }
    output_json.write_text(json.dumps(result, indent=2))
    return result


def _main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Path to Tsour et al. 2026 Supplementary Data 8 XLSX",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("output/slavov_saap_codon_distances.json"),
        help="Path for the output JSON",
    )
    args = parser.parse_args()

    if not args.input.is_file():
        sys.stderr.write(
            f"Input file not found: {args.input}\n\n"
            f"Download from https://decode.slavovlab.net/mass-spec/data\n"
            f"  (file: Supplementary_Data_8.High_confidence_SAAP_precursor_quant.xlsx)\n"
        )
        return 2

    args.output.parent.mkdir(parents=True, exist_ok=True)
    res = run_analysis(args.input, args.output)
    print(
        f"Wrote {args.output}: "
        f"{res['event_single_nt_count']}/{res['n_parseable_events']} events "
        f"single-NT ({res['event_single_nt_fraction'] * 100:.1f}%) "
        f"vs baseline {res['baseline_single_nt_fraction_all_380_pairs'] * 100:.1f}%, "
        f"binomial p={res['event_enrichment_binomial_p_one_sided']:.2g}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(_main())
