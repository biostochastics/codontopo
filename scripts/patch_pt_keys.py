#!/usr/bin/env python3.11
"""Patch manuscript_stats.json's per_table section with disaggregation keys.

Adds:
- informative_significant, informative_total
- near_standard_significant, near_standard_total
- marginal_table, marginal_p_bh, marginal_quantile

Uses the same logic as cli.py:_per_table_bh_disaggregation but reads inputs
from the emitted JSON files rather than re-running the pipeline. Idempotent.
"""

import json
import math
from pathlib import Path

from scipy.stats import norm

OUT = Path("output")

stats = json.loads((OUT / "manuscript_stats.json").read_text())
pertable_full = json.loads((OUT / "coloring_optimality.json").read_text())["per_table"]

# reassignment counts per table
db_counts = stats.get("reassignment", {}).get("per_table_counts", [])
counts_by_tid = {r["table_id"]: r["n_events"] for r in db_counts}

informative_sig = 0
informative_total = 0
near_sig = 0
near_total = 0
marginal_row = None
for row in pertable_full["per_table"]:
    tid = row["table_id"]
    if tid == 1:
        continue
    n_events = counts_by_tid.get(tid, 0)
    is_sig = row["p_value_bh"] < 0.05
    if n_events >= 3:
        informative_total += 1
        if is_sig:
            informative_sig += 1
        if marginal_row is None or row["p_value_bh"] > marginal_row["p_value_bh"]:
            marginal_row = row
    else:
        near_total += 1
        if is_sig:
            near_sig += 1

if marginal_row is None:
    marginal_row = max(pertable_full["per_table"], key=lambda r: r["p_value_bh"])

new_pt = dict(stats["per_table"])
new_pt.update(
    {
        "informative_significant": informative_sig,
        "informative_total": informative_total,
        "near_standard_significant": near_sig,
        "near_standard_total": near_total,
        "marginal_table": marginal_row["table_id"],
        "marginal_p_bh": round(marginal_row["p_value_bh"], 3),
        "marginal_quantile": round(marginal_row["quantile"], 2),
        # Also expose the full per_table list for downstream tables
        "per_table": pertable_full["per_table"],
    }
)
stats["per_table"] = new_pt

# --- derive per-metric improvement_pct where it wasn't emitted ---
for name, m in stats.get("metrics", {}).items():
    if m.get("improvement_pct") is None and m.get("null_mean") not in (None, 0):
        obs = m.get("observed")
        if obs is not None:
            m["improvement_pct"] = 100.0 * (m["null_mean"] - obs) / m["null_mean"]

# --- derive Q6 topology risk_ratio and depletion_fold if missing ---
q6 = stats.get("topology_avoidance_q6", {})
if q6.get("risk_ratio") is None and q6.get("rate_possible"):
    q6["risk_ratio"] = (
        q6["rate_observed"] / q6["rate_possible"]
        if q6.get("rate_observed") is not None
        else None
    )
if q6.get("depletion_fold") is None and q6.get("rate_observed"):
    q6["depletion_fold"] = (
        q6["rate_possible"] / q6["rate_observed"] if q6.get("rate_possible") else None
    )
if q6.get("risk_ratio_ci_95") is None and q6.get("risk_ratio") is not None:
    # Log-normal CI using x/n, K/N with a Katz-style approximation
    import math

    x, n = q6["observed_breaks"], q6["observed_total"]
    K, N = q6["possible_breaks"], q6["possible_total"]
    if x > 0 and K > 0 and n - x >= 0 and N - K >= 0:
        se_log = (
            math.sqrt(1 / x - 1 / n + 1 / K - 1 / N)
            if (1 / x - 1 / n + 1 / K - 1 / N) > 0
            else 0
        )
        rr = q6["risk_ratio"]
        q6["risk_ratio_ci_95"] = [
            rr * math.exp(-1.96 * se_log),
            rr * math.exp(1.96 * se_log),
        ]

# --- also inject the phys/topo Spearman rho into condlogit ---
evosim = json.loads((OUT / "evolutionary_simulation.json").read_text())
phys_topo = (evosim.get("diagnostics") or {}).get("phys_topo_correlation") or {}
if phys_topo.get("spearman_rho") is not None:
    stats["condlogit"]["phys_topo_rho"] = phys_topo["spearman_rho"]

# --- add mis_median_z / mis_worst_z from inverse-normal of the p-values ---
trna_block = stats.get("trna", {})
for _key, _pkey in (
    ("mis_median_z", "mis_median_p"),
    ("mis_worst_z", "mis_worst_p"),
    ("mis_best_z", "mis_best_p"),
):
    _p = trna_block.get(_pkey)
    if _p is not None and _p not in (0, 1):
        trna_block[_key] = float(norm.ppf(1 - _p))
stats["trna"] = trna_block

# --- lift the H(3,4) phylogenetic sensitivity into stats.phylo_k43 ---
try:
    phylo_k43 = json.loads((OUT / "phylogenetic_sensitivity_k43.json").read_text())
    stats["phylo_k43"] = {
        "adjacency": phylo_k43.get("adjacency"),
        "definition": phylo_k43.get("definition"),
        "all_significant": phylo_k43.get("all_clade_exclusions_significant"),
        "lineage_collapsed": phylo_k43.get("lineage_collapsed"),
        "clade_exclusion": phylo_k43.get("clade_exclusion"),
    }
except FileNotFoundError:
    pass

# --- lift classical Haig-Hurst AA-permutation null (B2b sensitivity) ---
try:
    hh = json.loads((OUT / "coloring_optimality_hh_aa_null.json").read_text())
    stats.setdefault("coloring", {})["haig_hurst_aa_null"] = hh
except FileNotFoundError:
    pass

# --- lift clade-exclusion aggregates into condlogit for §3.5 prose ---
try:
    cce = json.loads((OUT / "condlogit_clade_sensitivity.json").read_text())
    stats["condlogit"]["clade_exclusion"] = {
        "delta_M1_M3_min": cce.get("delta_M1_M3_min"),
        "delta_M1_M3_median": cce.get("delta_M1_M3_median"),
        "delta_M1_M3_max": cce.get("delta_M1_M3_max"),
        "delta_M2_M3_min": cce.get("delta_M2_M3_min"),
        "delta_M2_M3_median": cce.get("delta_M2_M3_median"),
        "delta_M2_M3_max": cce.get("delta_M2_M3_max"),
        "all_M3_favored_over_M1": cce.get("all_M3_favored_over_M1"),
        "all_M3_favored_over_M2": cce.get("all_M3_favored_over_M2"),
        "n_regimes": len(cce.get("rows", [])),
        "rows": cce.get("rows", []),
    }
except FileNotFoundError:
    pass

# --- pull the per-table proximity audit into stats.coloring for §S8 ---
prox_audit = pertable_full  # already loaded above
prox_key = json.loads((OUT / "coloring_optimality.json").read_text()).get(
    "per_table_proximity_audit"
)
if prox_key and "per_table_proximity_audit" not in stats["coloring"]:
    stats["coloring"]["per_table_proximity_audit"] = prox_key


def _clean_nan(obj):
    if isinstance(obj, dict):
        return {k: _clean_nan(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_clean_nan(v) for v in obj]
    if isinstance(obj, float) and (math.isnan(obj) or math.isinf(obj)):
        return None
    return obj


stats = _clean_nan(stats)
(OUT / "manuscript_stats.json").write_text(
    json.dumps(stats, indent=2, default=str, allow_nan=False)
)
print(
    f"patched: informative_sig/total={informative_sig}/{informative_total}, "
    f"near_sig/total={near_sig}/{near_total}, "
    f"marginal_table={new_pt['marginal_table']} (p_bh={new_pt['marginal_p_bh']}), "
    f"phys_topo_rho={stats['condlogit'].get('phys_topo_rho')}"
)
