"""Generate all supporting tables for the paper."""

import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

OUT = Path(__file__).parent.parent / "output" / "tables"
OUT.mkdir(parents=True, exist_ok=True)


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"  {path.name}: {len(rows)} rows")


def write_json(path: Path, data) -> None:
    with open(path, "w") as f:
        json.dump(data, f, indent=2, default=str)
    print(f"  {path.name}")


# ── Table 1: Claim hierarchy ─────────────────────────────────────────
print("Table 1: Claim hierarchy")
from codon_topo.reports.claim_hierarchy import CLAIM_HIERARCHY

rows = []
for c in CLAIM_HIERARCHY:
    rows.append({
        "id": c.id,
        "status": c.status.value,
        "statement": c.statement,
        "p_value": c.evidence_p_value,
        "null_model": c.null_model,
        "sample_size": c.sample_size,
        "publishable": c.is_publishable,
    })
write_csv(OUT / "T1_claim_hierarchy.csv", rows)


# ── Table 2: Disconnection catalogue ─────────────────────────────────
print("Table 2: Disconnection catalogue")
from codon_topo.core.homology import disconnection_catalogue
from codon_topo.core.genetic_codes import all_table_ids, get_code, get_code_name

rows = []
for tid in all_table_ids():
    code = get_code(tid)
    cat = disconnection_catalogue(code)
    for entry in cat:
        rows.append({
            "table_id": tid,
            "table_name": get_code_name(tid),
            "aa": entry["aa"],
            "n_components": entry["n_components"],
            "min_inter_distance": entry["min_inter_distance"],
            "reconnect_eps": entry["reconnect_eps"],
            "blocks": str(entry["blocks"]),
        })
write_csv(OUT / "T2_disconnection_catalogue.csv", rows)


# ── Table 3: Coloring optimality Monte Carlo ─────────────────────────
print("Table 3: Coloring optimality (n=10000)")
from codon_topo.analysis.coloring_optimality import monte_carlo_null

mc = monte_carlo_null(n_samples=10_000, seed=135325)
write_json(OUT / "T3_coloring_optimality.json", mc)
# Also a summary row for the paper table
write_csv(OUT / "T3_coloring_optimality_summary.csv", [{
    "observed_score": mc["observed_score"],
    "null_mean": mc["null_mean"],
    "null_std": mc["null_std"],
    "quantile_pct": mc["quantile_of_observed"],
    "p_value_conservative": mc["p_value_conservative"],
    "n_samples": mc["n_samples"],
}])


# ── Table 4: Per-table optimality ────────────────────────────────────
print("Table 4: Per-table optimality (n=1000)")
from codon_topo.analysis.coloring_optimality import per_table_optimality

pt = per_table_optimality(n_samples=1000, seed=135325)
rows = []
for t in pt["per_table"]:
    rows.append({
        "table_id": t["table_id"],
        "table_name": get_code_name(t["table_id"]),
        "observed_score": t["observed_score"],
        "quantile_pct": t["quantile"],
        "p_value": t["p_value"],
    })
write_csv(OUT / "T4_per_table_optimality.csv", rows)


# ── Table 5: Rho robustness sweep ────────────────────────────────────
print("Table 5: Rho robustness sweep (n=1000)")
from codon_topo.analysis.coloring_optimality import rho_robustness_sweep

rho = rho_robustness_sweep(
    rho_values=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
    n_samples=1000,
    seed=135325,
)
rows = []
for r in rho["per_rho"]:
    rows.append({
        "rho": r["rho"],
        "observed_score": r["observed_score"],
        "null_mean": r["null_mean"],
        "null_std": r["null_std"],
        "quantile_pct": r["quantile"],
        "p_value": r["p_value"],
    })
write_csv(OUT / "T5_rho_robustness.csv", rows)


# ── Table 6: Score decomposition ─────────────────────────────────────
print("Table 6: Score decomposition")
from codon_topo.analysis.coloring_optimality import score_decomposition_by_position

dec = score_decomposition_by_position()
write_json(OUT / "T6_score_decomposition.json", dec)


# ── Table 7: tRNA enrichment ─────────────────────────────────────────
print("Table 7: tRNA enrichment tests")
from codon_topo.analysis.trna_evidence import (
    fisher_exact_per_pairing,
    aa_label_permutation_test,
    trna_duplication_correlation_test,
    CURATED_REPERTOIRES,
)

sign = trna_duplication_correlation_test()
fisher = fisher_exact_per_pairing()
exact = aa_label_permutation_test()

# Per-pairing table combining all tests
rows = []
for i, (fp, ep) in enumerate(zip(fisher["per_pairing"], exact["per_pairing"])):
    rows.append({
        "pairing": f"{fp['disconnection_key']} vs {fp['control_key']}",
        "reassigned_aa": fp["reassigned_aa"],
        "fisher_OR": round(fp["odds_ratio"], 3),
        "fisher_p": round(fp["fisher_p"], 4),
        "aa_rank": ep["rank_among_aas"],
        "n_aas": ep["n_aas_compared"],
        "exact_p": round(ep["exact_p"], 4),
        "share_diff": round(ep["observed_share_diff"], 4),
    })
write_csv(OUT / "T7_trna_per_pairing.csv", rows)

# Summary
write_csv(OUT / "T7_trna_summary.csv", [{
    "test": "Sign test",
    "n_pairings": sign["n_pairings"],
    "statistic": f"{sign['n_with_elevated_trna']}/{sign['n_pairings']}",
    "p_value": round(sign["binomial_p_value"], 6),
}, {
    "test": "Fisher+Stouffer (all)",
    "n_pairings": fisher["n_pairings"],
    "statistic": f"Z={fisher['stouffer_z']:.3f}",
    "p_value": round(fisher["stouffer_p"], 6),
}, {
    "test": "Fisher+Stouffer (independent)",
    "n_pairings": fisher["n_independent_pairings"],
    "statistic": f"Z={fisher['stouffer_z_independent']:.3f}",
    "p_value": round(fisher["stouffer_p_independent"], 6),
}, {
    "test": "AA-label exact+Stouffer",
    "n_pairings": exact["n_pairings"],
    "statistic": f"Z={exact['stouffer_z']:.3f}",
    "p_value": round(exact["stouffer_p"], 6),
}])


# ── Table 8: Bit-position bias ────────────────────────────────────────
print("Table 8: Bit-position bias")
from codon_topo.analysis.reassignment_db import (
    bit_position_bias_weighted,
    bit_bias_deduplicated,
    nucleotide_position_bias,
)

bpb = bit_position_bias_weighted("mitochondrial")
dedup = bit_bias_deduplicated()
nuc = nucleotide_position_bias()

write_csv(OUT / "T8_bit_bias.csv", [{
    "test": "6-bin uniform",
    "observed": str(bpb["bit_counts_observed"]),
    "chi2": bpb["chi2_statistic_uniform_reference"],
    "p_value": bpb["chi2_p_value_uniform_reference"],
    "n_events": bpb["n_events"],
}, {
    "test": "6-bin mito Ts/Tv weighted",
    "observed": str(bpb["bit_counts_observed"]),
    "chi2": bpb["chi2_statistic_weighted"],
    "p_value": bpb["chi2_p_value_weighted"],
    "n_events": bpb["n_events"],
}, {
    "test": "6-bin de-duplicated",
    "observed": str(dedup["unique_events_histogram"]),
    "chi2": dedup["chi2_unique"],
    "p_value": dedup["p_value_unique"],
    "n_events": dedup["n_unique"],
}, {
    "test": "3-bin nucleotide position",
    "observed": str(nuc["position_counts"]),
    "chi2": nuc["chi2_statistic"],
    "p_value": nuc["chi2_p_value"],
    "n_events": sum(nuc["position_counts"]),
}])


# ── Table 9: Topology avoidance ───────────────────────────────────────
print("Table 9: Topology avoidance")
from codon_topo.analysis.synbio_feasibility import topology_avoidance_test

ta = topology_avoidance_test()
write_csv(OUT / "T9_topology_avoidance.csv", [{
    "observed_creates_disc": ta["observed_creates_disc"],
    "observed_total": ta["observed_total"],
    "rate_observed_pct": round(100 * ta["rate_observed"], 1),
    "possible_creates_disc": ta["possible_creates_disc"],
    "possible_total": ta["possible_total"],
    "rate_possible_pct": round(100 * ta["rate_possible"], 1),
    "hypergeom_p": ta["hypergeom_p"],
    "fisher_p": ta["fisher_p"],
    "odds_ratio": ta["odds_ratio"],
}])


# ── Table 10: Reassignment database ──────────────────────────────────
print("Table 10: Reassignment database")
from codon_topo.analysis.reassignment_db import build_reassignment_db

db = build_reassignment_db()
rows = []
for e in db:
    rows.append({
        "table_id": e.table_id,
        "table_name": e.table_name,
        "codon": e.codon,
        "source_aa": e.source_aa,
        "target_aa": e.target_aa,
        "hamming_to_nearest_target": e.hamming_to_nearest_target,
    })
write_csv(OUT / "T10_reassignment_db.csv", rows)

print(f"\nDone. {len(list(OUT.glob('*')))} files written to {OUT}/")
