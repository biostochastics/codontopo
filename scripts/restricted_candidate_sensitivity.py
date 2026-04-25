#!/usr/bin/env python3
"""Restricted-candidate-set sensitivity analysis for the conditional logit.

The full 1,280-move candidate set includes biologically catastrophic moves
(e.g., reassigning AUG-Met, multi-codon changes, reassignments to stop in
essential codons). Models with strongly negative coefficients on Δ_topo and
Δ_phys may be partly rediscovering that natural reassignments are not
biologically catastrophic, inflating apparent ΔAICc.

This script refits M1-M4 (Q_6 topology) on a restricted candidate set:
candidates whose target AA is already serviced by a codon at Hamming distance
≤d from the reassigned codon (delta_trna ≤ d), with d ∈ {1, 2, 3}. The
observed move is always retained regardless of its delta_trna, so likelihoods
remain comparable.

The qualitative explanatory claim ("topology adds value beyond physicochemistry")
should survive at every threshold; the unrestricted-set ΔAICc magnitudes are
upper bounds, with the d=2 filter giving a more biologically-calibrated effect
size.

Output: writes output/condlogit_restricted_candidate.json.
"""

from __future__ import annotations

import json
import sys
from copy import copy
from pathlib import Path

repo_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(repo_root / "src"))

from codon_topo.analysis.evolutionary_simulation import (  # noqa: E402
    StepChoiceSet,
    build_all_choice_sets,
    fit_all_models,
    likelihood_ratio_test,
)


def _filter_choice_set(
    cs: StepChoiceSet,
    max_trna: int,
) -> StepChoiceSet:
    """Return a copy of cs whose candidate list is filtered to delta_trna ≤ max_trna.

    The observed move is always retained, even if its delta_trna exceeds max_trna,
    so likelihood is well-defined.
    """
    kept = [m for m in cs.candidates if (m.delta_trna <= max_trna) or m.is_observed]
    new_cs = copy(cs)
    new_cs.candidates = kept
    return new_cs


def filter_all_choice_sets(
    all_cs: dict[int, list[list[StepChoiceSet]]],
    max_trna: int,
) -> dict[int, list[list[StepChoiceSet]]]:
    out: dict[int, list[list[StepChoiceSet]]] = {}
    for tid, orderings in all_cs.items():
        filtered_orderings = [
            [_filter_choice_set(cs, max_trna) for cs in order] for order in orderings
        ]
        out[tid] = filtered_orderings
    return out


def summarise_candidate_counts(
    all_cs: dict[int, list[list[StepChoiceSet]]],
) -> dict:
    counts = []
    obs_in_set = 0
    obs_total = 0
    for orderings in all_cs.values():
        for cs in orderings[0]:  # first ordering only for counting
            counts.append(len(cs.candidates))
            obs = cs.observed_move
            if obs is not None:
                obs_total += 1
                if any(m.is_observed for m in cs.candidates):
                    obs_in_set += 1
    return {
        "n_choice_sets": len(counts),
        "candidates_min": min(counts) if counts else 0,
        "candidates_max": max(counts) if counts else 0,
        "candidates_mean": float(sum(counts) / max(len(counts), 1)),
        "observed_in_filtered_set": obs_in_set,
        "observed_total": obs_total,
    }


def main() -> None:
    print("Building full candidate choice sets...")
    all_cs_full = build_all_choice_sets()
    full_summary = summarise_candidate_counts(all_cs_full)
    print(
        f"  Full: {full_summary['n_choice_sets']} choice sets, "
        f"~{full_summary['candidates_mean']:.0f} candidates each"
    )

    results: dict = {
        "_doc": (
            "Restricted-candidate-set sensitivity for the conditional logit. "
            "Candidates are filtered by delta_trna (Hamming distance to nearest "
            "existing target-AA codon). Observed moves are always retained."
        ),
        "full_set_summary": full_summary,
        "by_max_trna": {},
    }

    # Primary threshold: Hamming ≤2 from at least one existing target-AA
    # codon. Also report ≤1 (most stringent biological-plausibility filter)
    # and ≤3 (looser bound) to bracket the choice.
    for max_trna in (1, 2, 3):
        print(f"\n--- Refitting with delta_trna ≤ {max_trna} ---")
        all_cs_filt = filter_all_choice_sets(all_cs_full, max_trna)
        filt_summary = summarise_candidate_counts(all_cs_filt)
        print(f"  Filtered: ~{filt_summary['candidates_mean']:.0f} candidates each")
        print(
            f"  Observed events retained: "
            f"{filt_summary['observed_in_filtered_set']}/{filt_summary['observed_total']}"
        )

        fits = fit_all_models(all_cs_filt)

        m1 = fits.get("M1_phys")
        m2 = fits.get("M2_topo")
        m3 = fits.get("M3_phys_topo")
        m4 = fits.get("M4_full")
        m3k = fits.get("M3_phys_topo_k43")
        m2k = fits.get("M2_topo_k43")

        block: dict = {
            "max_trna": max_trna,
            "candidate_summary": filt_summary,
            "model_aicc": {name: float(f["aicc"]) for name, f in fits.items()},
            "model_log_likelihood": {
                name: float(f["log_likelihood"]) for name, f in fits.items()
            },
            "delta_aicc": {},
            "lr_tests": {},
        }

        if m1 and m3:
            block["delta_aicc"]["M1_to_M3"] = m1["aicc"] - m3["aicc"]
            block["lr_tests"]["M1_vs_M3"] = likelihood_ratio_test(m1, m3)
        if m2 and m3:
            block["delta_aicc"]["M2_to_M3"] = m2["aicc"] - m3["aicc"]
            block["lr_tests"]["M2_vs_M3"] = likelihood_ratio_test(m2, m3)
        if m3 and m4:
            block["delta_aicc"]["M3_to_M4"] = m3["aicc"] - m4["aicc"]
            block["lr_tests"]["M3_vs_M4"] = likelihood_ratio_test(m3, m4)
        if m1 and m3k:
            block["delta_aicc"]["M1_to_M3_k43"] = m1["aicc"] - m3k["aicc"]
            block["lr_tests"]["M1_vs_M3_k43"] = likelihood_ratio_test(m1, m3k)
        if m2k and m3k:
            block["delta_aicc"]["M2k43_to_M3k43"] = m2k["aicc"] - m3k["aicc"]
            block["lr_tests"]["M2_k43_vs_M3_k43"] = likelihood_ratio_test(m2k, m3k)

        # Effective candidate-set size summary for this filter.
        block["effective_n_candidates_mean"] = filt_summary["candidates_mean"]

        results["by_max_trna"][str(max_trna)] = block

        print(
            f"  ΔAICc(M1→M3) = {block['delta_aicc'].get('M1_to_M3', float('nan')):.1f}"
        )
        print(
            f"  ΔAICc(M2→M3) = {block['delta_aicc'].get('M2_to_M3', float('nan')):.1f}"
        )
        print(
            f"  ΔAICc(M3→M4) = {block['delta_aicc'].get('M3_to_M4', float('nan')):.1f}"
        )
        if "M1_to_M3_k43" in block["delta_aicc"]:
            print(f"  ΔAICc(M1→M3_H(3,4)) = {block['delta_aicc']['M1_to_M3_k43']:.1f}")

    out_path = repo_root / "output" / "condlogit_restricted_candidate.json"
    out_path.write_text(json.dumps(results, indent=2, default=str))
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
