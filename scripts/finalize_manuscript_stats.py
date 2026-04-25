#!/usr/bin/env python3.11
"""Splice new sensitivity blocks into manuscript_stats.json.

Used as a one-shot post-processor when `codon-topo all` was launched against
an older cli.py that did not yet include the new sensitivity keys
(topology_audit, topology_q6_encoding_sweep, condlogit.encoding_robustness,
topology_avoidance.k43_definitions_audit). Reads the per-analysis JSONs
and injects the new blocks into manuscript_stats.json without re-running
any heavy computation.

Idempotent: running multiple times produces the same result.
"""

from __future__ import annotations

import json
from pathlib import Path

OUTPUT = Path("/Users/biostochastics/clayworth/output")


def load_json(path: Path) -> dict:
    with path.open() as f:
        return json.load(f)


def write_json(path: Path, data: dict) -> None:
    with path.open("w") as f:
        json.dump(data, f, indent=2, default=str)


def main() -> None:
    stats = load_json(OUTPUT / "manuscript_stats.json")
    topo = load_json(OUTPUT / "topology_avoidance.json")
    evosim = load_json(OUTPUT / "evolutionary_simulation.json")
    clade_path = OUTPUT / "condlogit_clade_sensitivity.json"
    clade = load_json(clade_path) if clade_path.exists() else None
    restricted_path = OUTPUT / "condlogit_restricted_candidate.json"
    restricted = load_json(restricted_path) if restricted_path.exists() else None

    # 1. Inject topology_audit (2x2 sensitivity)
    if "definitions_audit" in topo:
        stats["topology_audit"] = topo["definitions_audit"]
        print("  [+] topology_audit added")

    # 2. Inject topology_q6_encoding_sweep
    if "Q6_encoding_sweep" in topo:
        stats["topology_q6_encoding_sweep"] = topo["Q6_encoding_sweep"]
        print("  [+] topology_q6_encoding_sweep added")

    # 3. Inject denominator_sensitivity
    if "denominator_sensitivity" in topo:
        stats["topology_denominator_sensitivity"] = topo["denominator_sensitivity"]
        print("  [+] topology_denominator_sensitivity added")

    # 4. Inject encoding_robustness from evolutionary_simulation.json into
    #    the condlogit block in manuscript_stats.json
    if "encoding_robustness" in evosim:
        if "condlogit" not in stats:
            stats["condlogit"] = {}
        stats["condlogit"]["encoding_robustness"] = evosim["encoding_robustness"]
        print("  [+] condlogit.encoding_robustness added")

    # 5. Inject the K_4^3 model fits if missing (for completeness)
    if "model_fits" in evosim:
        if "condlogit" not in stats:
            stats["condlogit"] = {}
        existing_fits = stats["condlogit"].get("model_fits", {})
        for name, f in evosim["model_fits"].items():
            if name not in existing_fits:
                existing_fits[name] = f
                print(f"  [+] condlogit.model_fits.{name} added")
        stats["condlogit"]["model_fits"] = existing_fits

    # 6. Inject lr_tests (extended) if available
    if "likelihood_ratio_tests" in evosim:
        if "condlogit" not in stats:
            stats["condlogit"] = {}
        existing_lr = stats["condlogit"].get("lr_tests", {})
        for name, lr in evosim["likelihood_ratio_tests"].items():
            if name not in existing_lr:
                existing_lr[name] = lr
                print(f"  [+] condlogit.lr_tests.{name} added")
        stats["condlogit"]["lr_tests"] = existing_lr

    # 7. Inject condlogit clade-exclusion sensitivity
    if clade is not None:
        if "condlogit" not in stats:
            stats["condlogit"] = {}
        stats["condlogit"]["clade_exclusion"] = clade
        print("  [+] condlogit.clade_exclusion added")

    # 7b. Inject condlogit restricted-candidate sensitivity
    if restricted is not None:
        if "condlogit" not in stats:
            stats["condlogit"] = {}
        stats["condlogit"]["restricted_candidate"] = restricted
        print("  [+] condlogit.restricted_candidate added")

    # 8. Inject phys_topo_rho from evosim diagnostics block (used by abstract/results)
    diag = evosim.get("diagnostics", {})
    ptc = diag.get("phys_topo_correlation", {})
    if "spearman_rho" in ptc:
        if "condlogit" not in stats:
            stats["condlogit"] = {}
        stats["condlogit"]["phys_topo_rho"] = ptc["spearman_rho"]
        stats["condlogit"]["phys_topo_correlation"] = ptc
        print(f"  [+] condlogit.phys_topo_rho = {ptc['spearman_rho']:.4f} added")

    # 9. Inject best-model posterior predictive
    pp = diag.get("posterior_predictive", {})
    if pp:
        if "condlogit" not in stats:
            stats["condlogit"] = {}
        stats["condlogit"]["posterior_predictive"] = pp
        print("  [+] condlogit.posterior_predictive added")

    # 10. Backfill improvement_pct in each metric block (cli.py
    # forwards an improvement_pct field that multi_metric_sensitivity
    # does not actually compute, leaving it as None and breaking the
    # manuscript prose). Compute from null_mean and observed.
    metrics = stats.get("metrics", {})
    for name, m in metrics.items():
        if (
            m.get("improvement_pct") is None
            and m.get("null_mean")
            and m.get("observed")
        ):
            m["improvement_pct"] = (
                100.0 * (m["null_mean"] - m["observed"]) / m["null_mean"]
            )
            print(
                f"  [+] metrics.{name}.improvement_pct = {m['improvement_pct']:.2f}% backfilled"
            )

    # 11. Backfill risk_ratio for topology_avoidance_q6 from topology_avoidance.json Q6 block
    tq6 = stats.get("topology_avoidance_q6", {})
    topo_q6 = topo.get("Q6", {})
    if (
        tq6.get("risk_ratio") is None
        and topo_q6.get("rate_observed")
        and topo_q6.get("rate_possible")
    ):
        rr_q6 = topo_q6["rate_observed"] / topo_q6["rate_possible"]
        tq6["risk_ratio"] = rr_q6
        print(f"  [+] topology_avoidance_q6.risk_ratio = {rr_q6:.3f} backfilled")
    if (
        tq6.get("depletion_fold") is None
        and topo_q6.get("rate_possible")
        and topo_q6.get("rate_observed")
    ):
        depletion = topo_q6["rate_possible"] / max(topo_q6["rate_observed"], 1e-10)
        tq6["depletion_fold"] = depletion
        print(
            f"  [+] topology_avoidance_q6.depletion_fold = {depletion:.2f}x backfilled"
        )

    # 12. Backfill per_table fields used by the manuscript prose:
    # marginal_table, marginal_p_bh, informative_*/near_standard_* counts.
    # The "marginal" exception is the table with the *highest* BH-corrected
    # p-value (closest to or above 0.05). Informative-distance tables have
    # >=3 codon reassignments from the standard code; near-standard have
    # <=2. The single source of truth for reassignment counts is
    # codon_topo.core.genetic_codes.get_changes.
    pt = stats.get("per_table", {})
    rows = pt.get("per_table", []) or []
    if rows and pt.get("marginal_table") is None:
        # Find row with highest p_bh
        marginal = max(
            rows,
            key=lambda r: r.get(
                "p_bh", r.get("p_value_bh", r.get("p_value_conservative_bh", 0.0))
            ),
        )
        marginal_p = marginal.get(
            "p_bh", marginal.get("p_value_bh", marginal.get("p_value_conservative_bh"))
        )
        pt["marginal_table"] = marginal["table_id"]
        pt["marginal_p_bh"] = round(float(marginal_p), 4)
        print(
            f"  [+] per_table.marginal_table = {pt['marginal_table']} "
            f"(p_bh = {pt['marginal_p_bh']}) backfilled"
        )
    if rows and pt.get("informative_total") is None:
        # Need reassignment counts to partition tables into informative
        # (>=3 reassignments) vs near-standard (<=2). Lazy import to keep
        # this script lightweight when the codon_topo package is not on
        # PYTHONPATH (e.g. the standalone post-processor).
        import sys

        sys.path.insert(0, str(OUTPUT.parent / "src"))
        try:
            from codon_topo.core.genetic_codes import get_changes
        except ImportError:
            print("  [!] codon_topo not importable; skipping informative split")
        else:
            inf_total = inf_sig = near_total = near_sig = 0
            INFORMATIVE_THRESHOLD = 3
            for r in rows:
                tid = r["table_id"]
                if tid == 1:
                    # Standard code is the reference; it's neither
                    # informative-distance nor a near-standard variant
                    # (it IS the standard). Excluded from the split.
                    continue
                n_changes = len(get_changes(tid))
                p_bh = r.get(
                    "p_bh", r.get("p_value_bh", r.get("p_value_conservative_bh", 1.0))
                )
                is_sig = float(p_bh) < 0.05
                if n_changes >= INFORMATIVE_THRESHOLD:
                    inf_total += 1
                    if is_sig:
                        inf_sig += 1
                else:
                    near_total += 1
                    if is_sig:
                        near_sig += 1
            pt["informative_total"] = inf_total
            pt["informative_significant"] = inf_sig
            pt["near_standard_total"] = near_total
            pt["near_standard_significant"] = near_sig
            print(
                f"  [+] per_table informative={inf_sig}/{inf_total}, "
                f"near-standard={near_sig}/{near_total} backfilled"
            )

    write_json(OUTPUT / "manuscript_stats.json", stats)
    print(f"\nUpdated {OUTPUT / 'manuscript_stats.json'}")


if __name__ == "__main__":
    main()
