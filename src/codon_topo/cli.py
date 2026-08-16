"""Command-line interface for codon-topo."""

import json
from pathlib import Path

import click
import numpy as np

from codon_topo import DEFAULT_SEED


class _NumpyEncoder(json.JSONEncoder):
    """JSON encoder that converts numpy types to native Python types."""

    def default(self, o):
        if isinstance(o, np.bool_):
            return bool(o)
        if isinstance(o, np.integer):
            return int(o)
        if isinstance(o, np.floating):
            # NaN/Inf are not JSON-spec; emit null so downstream readers
            # (e.g. Typst's json()) parse cleanly.
            if not np.isfinite(o):
                return None
            return float(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        return super().default(o)


def _write_json_strict(path, data) -> None:
    """JSON dump with allow_nan=False so NaN/Inf raise rather than emit non-spec.

    Combined with _NumpyEncoder above, this guarantees that every JSON file
    written by the pipeline parses under strict (RFC 8259) JSON readers,
    including the Typst json() function used by manuscript.typ and
    supplement.typ.
    """
    with open(path, "w") as f:
        json.dump(data, f, indent=2, cls=_NumpyEncoder, default=str, allow_nan=False)


def _per_table_reassignment_counts(db: list) -> list[dict]:
    """Summary of raw reassignment events per NCBI translation table.

    Used by Supplement §S11. Returns a sorted list of dicts:
        [{"table_id": int, "table_name": str, "n_events": int}, ...]
    """
    from collections import defaultdict

    by_table: dict[int, dict] = defaultdict(
        lambda: {"table_id": 0, "table_name": "", "n_events": 0}
    )
    for e in db:
        row = by_table[e.table_id]
        row["table_id"] = e.table_id
        row["table_name"] = e.table_name
        row["n_events"] += 1
    return sorted(by_table.values(), key=lambda r: r["table_id"])


def _per_table_bh_disaggregation(
    pertable_result: dict, db: list, informative_min_events: int = 3
) -> dict:
    """Disaggregate per-table BH-FDR results into informative vs near-standard.

    A table is "informative-distance" when it has >= informative_min_events
    reassignments relative to the standard code; otherwise it is
    "near-standard". Also returns the marginal (largest p_value_bh) table for
    the sole exception language used in main-text §3.3.
    """
    counts_by_tid = {
        r["table_id"]: r["n_events"] for r in _per_table_reassignment_counts(db)
    }

    informative_sig = 0
    informative_total = 0
    near_sig = 0
    near_total = 0
    marginal_row = None
    for row in pertable_result["per_table"]:
        tid = row["table_id"]
        if tid == 1:
            continue  # standard code — trivially significant, excluded from variant framing
        n_events = counts_by_tid.get(tid, 0)
        is_sig = row["p_value_bh"] < 0.05
        if n_events >= informative_min_events:
            informative_total += 1
            if is_sig:
                informative_sig += 1
        else:
            near_total += 1
            if is_sig:
                near_sig += 1
        # Marginal = largest BH p-value among informative-distance tables that
        # are still below 0.10; falls back to the largest BH p-value overall
        # if no informative table sits in [0.05, 0.10].
        if n_events >= informative_min_events:
            if marginal_row is None or row["p_value_bh"] > marginal_row["p_value_bh"]:
                marginal_row = row

    if marginal_row is None:
        marginal_row = max(pertable_result["per_table"], key=lambda r: r["p_value_bh"])

    return {
        "informative_significant": informative_sig,
        "informative_total": informative_total,
        "near_standard_significant": near_sig,
        "near_standard_total": near_total,
        "marginal_table": marginal_row["table_id"],
        "marginal_p_bh": round(marginal_row["p_value_bh"], 3),
        "marginal_quantile": round(marginal_row["quantile"], 2),
    }


def _json_out(data: dict | list) -> None:
    click.echo(json.dumps(data, indent=2, cls=_NumpyEncoder, default=str))


def _try_rich_table(headers: list[str], rows: list[list[str]], title: str = "") -> None:
    try:
        from rich.console import Console
        from rich.table import Table

        table = Table(title=title)
        for h in headers:
            table.add_column(h)
        for row in rows:
            table.add_row(*[str(c) for c in row])
        Console().print(table)
    except ImportError:
        if title:
            click.echo(f"\n{title}")
            click.echo("=" * len(title))
        click.echo("\t".join(headers))
        for row in rows:
            click.echo("\t".join(str(c) for c in row))


@click.group()
@click.version_option(package_name="codon-topo")
def main() -> None:
    """Codon Geometry Validation & Prediction Engine.

    Analyze the algebraic structure of genetic codes encoded as 6-bit
    binary vectors in GF(2)^6.
    """


@main.command()
@click.option("--table", "table_id", default=1, help="NCBI translation table ID.")
@click.option("--all-tables", is_flag=True, help="Run across all 27 NCBI tables.")
@click.option("--json", "as_json", is_flag=True, help="Output as JSON.")
def filtration(table_id: int, all_tables: bool, as_json: bool) -> None:
    """Check two-fold and four-fold filtration properties."""
    from codon_topo.core.filtration import analyze_filtration
    from codon_topo.core.genetic_codes import all_table_ids, get_code

    tables = all_table_ids() if all_tables else [table_id]
    results = []
    for tid in tables:
        code = get_code(tid)
        r = analyze_filtration(code)
        r["table_id"] = tid
        results.append(r)

    if as_json:
        _json_out(results if all_tables else results[0])
        return

    headers = ["Table", "2-fold Pass", "2-fold Fail", "4-fold Pass", "4-fold Fail"]
    rows = [
        [
            str(r["table_id"]),
            str(r["twofold_pass"]),
            str(r["twofold_fail"]),
            str(r["fourfold_pass"]),
            str(r["fourfold_fail"]),
        ]
        for r in results
    ]
    _try_rich_table(headers, rows, title="Filtration Report")


@main.command()
@click.option("--table", "table_id", default=1, help="NCBI translation table ID.")
@click.option("--all-tables", is_flag=True, help="Run across all 27 NCBI tables.")
@click.option(
    "--extended", is_flag=True, help="Run null_model_c_extended across 24 encodings."
)
@click.option("--json", "as_json", is_flag=True, help="Output as JSON.")
def disconnections(
    table_id: int, all_tables: bool, extended: bool, as_json: bool
) -> None:
    """Compute the disconnection catalogue."""
    from codon_topo.core.homology import disconnection_catalogue
    from codon_topo.core.genetic_codes import all_table_ids, get_code

    if extended:
        from codon_topo.analysis.null_models import null_model_c_extended

        result = null_model_c_extended(code=get_code(table_id))
        if as_json:
            _json_out(result)
            return
        click.echo(f"Extended disconnection analysis (table {table_id}, 24 encodings)")
        for aa in result["universal_disconnected_aas"]:
            details = result["invariant_details"][aa]
            click.echo(
                f"  {aa}: distances={details['min_distance_values']}, "
                f"invariant={details['distance_is_invariant']}"
            )
        return

    tables = all_table_ids() if all_tables else [table_id]
    all_results = []
    for tid in tables:
        code = get_code(tid)
        cat = disconnection_catalogue(code)
        for entry in cat:
            entry["table_id"] = tid
        all_results.extend(cat)

    if as_json:
        _json_out(all_results)
        return

    if not all_results:
        click.echo("No disconnections found.")
        return

    headers = ["Table", "AA", "Components", "Min Inter-distance", "Reconnect ε"]
    rows = [
        [
            str(e["table_id"]),
            e["aa"],
            str(e["n_components"]),
            str(e["min_inter_distance"]),
            str(e.get("reconnect_eps", "—")),
        ]
        for e in all_results
    ]
    _try_rich_table(headers, rows, title="Disconnection Catalogue")


@main.command()
@click.option(
    "--null",
    "null_type",
    default="freeland_hurst",
    type=click.Choice(["freeland_hurst", "class_size", "haig_hurst_aa"]),
    help="Null model type (freeland_hurst = quartet-pattern shuffle; haig_hurst_aa = classical AA-permutation).",
)
@click.option("--n", "n_samples", default=10_000, help="Number of Monte Carlo samples.")
@click.option("--seed", default=DEFAULT_SEED, help="Random seed.")
@click.option("--no-stops", is_flag=True, help="Exclude stop-codon edges from scoring.")
@click.option("--json", "as_json", is_flag=True, help="Output as JSON.")
def coloring(
    null_type: str, n_samples: int, seed: int, no_stops: bool, as_json: bool
) -> None:
    """Run hypercube coloring Monte Carlo optimality test."""
    from codon_topo.analysis.coloring_optimality import monte_carlo_null

    result = monte_carlo_null(
        n_samples=n_samples,
        seed=seed,
        null_type=null_type,
        include_stops=not no_stops,
    )

    if as_json:
        _json_out(result)
        return

    click.echo(f"Hypercube Coloring Optimality ({null_type}, n={n_samples})")
    click.echo(f"  Observed score:   {result['observed_score']:.2f}")
    click.echo(
        f"  Null mean ± std:  {result['null_mean']:.2f} ± {result['null_std']:.2f}"
    )
    click.echo(f"  Quantile:         {result['quantile_of_observed']:.2f}%")
    click.echo(f"  P-value (cons):   {result['p_value_conservative']:.6f}")
    click.echo(f"  Interpretation:   {result['interpretation']}")


@main.command("bit-bias")
@click.option(
    "--compartment",
    default="mitochondrial",
    type=click.Choice(["uniform", "nuclear", "mitochondrial"]),
    help="Ts/Tv weight compartment.",
)
@click.option("--json", "as_json", is_flag=True, help="Output as JSON.")
def bit_bias(compartment: str, as_json: bool) -> None:
    """Run bit-position bias test under Ts/Tv-weighted null."""
    from codon_topo.analysis.reassignment_db import bit_position_bias_weighted

    result = bit_position_bias_weighted(compartment)

    if as_json:
        _json_out(result)
        return

    click.echo(f"Bit-Position Bias Test (compartment={compartment})")
    click.echo(f"  Observed counts:    {result['bit_counts_observed']}")
    click.echo(f"  Weighted p-value:   {result['chi2_p_value_weighted']:.4f}")
    click.echo(f"  Uniform p-value:    {result['chi2_p_value_uniform_reference']:.4f}")
    click.echo(f"  N events:           {result['n_events']}")


@main.command()
@click.option("--json", "as_json", is_flag=True, help="Output as JSON.")
def trna(as_json: bool) -> None:
    """Run tRNA duplication correlation test."""
    from codon_topo.analysis.trna_evidence import trna_duplication_correlation_test

    result = trna_duplication_correlation_test()

    if as_json:
        _json_out(result)
        return

    click.echo("tRNA Duplication Correlation Test")
    click.echo(f"  Pairings:           {result['n_pairings']}")
    click.echo(
        f"  Elevated tRNA:      {result['n_with_elevated_trna']}/{result['n_pairings']}"
    )
    click.echo(f"  Binomial p-value:   {result['binomial_p_value']:.4f}")
    click.echo(f"  Mean excess:        {result['mean_excess_trna_count']:.1f}")
    click.echo(f"  Caveat:             {result['caveat'][:80]}...")


@main.command()
@click.option(
    "--offline", is_flag=True, help="Skip cBioPortal API call (report structure only)."
)
@click.option("--json", "as_json", is_flag=True, help="Output as JSON.")
def kras(offline: bool, as_json: bool) -> None:
    """Run KRAS-Fano enrichment test (expected: null result)."""
    from codon_topo.analysis.cosmic_query import (
        fano_predictions_for_kras,
        ws4_gate_decision,
    )

    preds = fano_predictions_for_kras()

    if offline:
        result = {
            "mode": "offline",
            "fano_predictions": preds,
            "note": "Offline mode: no cBioPortal query. Use without --offline to test.",
        }
    else:
        from codon_topo.analysis.cosmic_query import CBioPortalClient

        client = CBioPortalClient()
        raw = client.get_kras_mutations()
        result = ws4_gate_decision(raw)
        result["fano_predictions"] = preds

    if as_json:
        _json_out(result)
        return

    click.echo("KRAS-Fano Enrichment Test")
    if offline:
        click.echo("  Mode: OFFLINE (predictions only)")
    else:
        click.echo(f"  Gate passed:  {result.get('pass', 'N/A')}")
        click.echo(f"  Significant:  {result.get('significant_variants', [])}")
    click.echo("  Fano predictions:")
    for var, info in preds.items():
        click.echo(
            f"    {var}: partner={info['fano_partner_codon']} ({info['fano_partner_aa']})"
        )


@main.command()
@click.option("--json", "as_json", is_flag=True, help="Output as JSON.")
def claims(as_json: bool) -> None:
    """Print the claim hierarchy table."""
    from codon_topo.reports.claim_hierarchy import (
        CLAIM_HIERARCHY,
        hierarchy_summary_table,
    )

    if as_json:
        out = [
            {
                "id": c.id,
                "status": c.status.value,
                "statement": c.statement,
                "p_value": c.evidence_p_value,
                "null_model": c.null_model,
                "sample_size": c.sample_size,
            }
            for c in CLAIM_HIERARCHY
        ]
        _json_out(out)
        return

    click.echo(hierarchy_summary_table())


@main.command("topology-avoidance")
@click.option("--json", "as_json", is_flag=True, help="Output as JSON.")
def topology_avoidance(as_json: bool) -> None:
    """Test whether natural reassignments avoid creating disconnections."""
    from codon_topo.analysis.synbio_feasibility import topology_avoidance_test

    result = topology_avoidance_test()

    if as_json:
        _json_out(result)
        return

    click.echo("Topology Avoidance Test")
    click.echo(
        f"  Observed: {result['observed_creates_disc']}/{result['observed_total']} "
        f"create new disconnections ({result['rate_observed']:.1%})"
    )
    click.echo(
        f"  Possible: {result['possible_creates_disc']}/{result['possible_total']} "
        f"create new disconnections ({result['rate_possible']:.1%})"
    )
    click.echo(f"  Fisher OR={result['odds_ratio']:.3f}, p={result['fisher_p']:.6f}")


@main.command("rho-sweep")
# Default n matches the coloring-optimality MC convention (CLAUDE.md: 10,000).
@click.option("--n", "n_samples", default=10_000, help="Monte Carlo samples per rho.")
@click.option("--seed", default=135325, help="Random seed.")
@click.option("--json", "as_json", is_flag=True, help="Output as JSON.")
def rho_sweep(n_samples: int, seed: int, as_json: bool) -> None:
    """Sweep the diagonal-edge weight rho to test optimality robustness."""
    from codon_topo.analysis.coloring_optimality import rho_robustness_sweep

    result = rho_robustness_sweep(n_samples=n_samples, seed=seed)

    if as_json:
        _json_out(result)
        return

    click.echo(f"Rho Robustness Sweep (n={n_samples} per rho)")
    for r in result["per_rho"]:
        click.echo(
            f"  rho={r['rho']:.2f}: quantile={r['quantile']:.1f}%, p={r['p_value']:.4f}"
        )
    click.echo(f"  All p<0.01: {result['all_significant_p01']}")


@main.command("decompose")
@click.option("--json", "as_json", is_flag=True, help="Output as JSON.")
def decompose(as_json: bool) -> None:
    """Decompose mismatch score by nucleotide position and AA pair."""
    from codon_topo.analysis.coloring_optimality import score_decomposition_by_position

    result = score_decomposition_by_position()

    if as_json:
        _json_out(result)
        return

    click.echo(f"Score Decomposition (total={result['total_score']:.0f})")
    for pos, val in result["by_nucleotide_position"].items():
        frac = result["position_fractions"][pos]
        click.echo(f"  {pos:12s}: {val:.0f} ({frac:.1%})")
    click.echo("Top AA pairs:")
    for entry in result["top_aa_pairs"][:5]:
        click.echo(
            f"  {entry['pair'][0]:4s}-{entry['pair'][1]:4s}: "
            f"{entry['score']:.0f} ({entry['fraction']:.1%})"
        )


@main.command("per-table")
# Default n matches the coloring-optimality MC convention (CLAUDE.md: 10,000).
@click.option("--n", "n_samples", default=10_000, help="Monte Carlo samples per table.")
@click.option("--seed", default=135325, help="Random seed.")
@click.option("--json", "as_json", is_flag=True, help="Output as JSON.")
def per_table(n_samples: int, seed: int, as_json: bool) -> None:
    """Test whether variant codes preserve coloring optimality."""
    from codon_topo.analysis.coloring_optimality import per_table_optimality

    result = per_table_optimality(n_samples=n_samples, seed=seed)

    if as_json:
        _json_out(result)
        return

    click.echo(f"Per-Table Optimality (n={n_samples} per table)")
    click.echo(
        f"  Significant (p<0.05 raw): {result['n_significant_p05_raw']}/{result['n_tables']}"
    )
    click.echo(
        f"  Significant (BH FDR<0.05): {result['n_significant_p05_bh']}/{result['n_tables']}"
    )
    click.echo(f"  Mean quantile: {result['mean_quantile']:.2f}%")
    headers = ["Table", "Quantile", "P-value"]
    rows = [
        [str(t["table_id"]), f"{t['quantile']:.1f}%", f"{t['p_value']:.4f}"]
        for t in result["per_table"]
    ]
    _try_rich_table(headers, rows, title="Per-Table Optimality")


@main.command("metric-sensitivity")
@click.option(
    "--n", "n_samples", default=1_000, help="Monte Carlo sample size per metric."
)
@click.option("--seed", default=DEFAULT_SEED, help="Random seed.")
@click.option("--json", "as_json", is_flag=True, help="JSON output.")
def metric_sensitivity(n_samples: int, seed: int, as_json: bool) -> None:
    """Test coloring optimality across multiple distance metrics."""
    from codon_topo.analysis.coloring_optimality import multi_metric_sensitivity

    results = multi_metric_sensitivity(n_samples=n_samples, seed=seed)
    if as_json:
        _json_out(results)
    else:
        _try_rich_table(
            ["Metric", "Quantile", "p-value", "Effect size (z)"],
            [
                [
                    r["metric"],
                    f"{r['quantile']:.1f}%",
                    f"{r['p_value_conservative']:.4f}",
                    f"{r['effect_size_z']:.2f}",
                ]
                for r in results["per_metric"]
            ],
            title="Multi-Metric Sensitivity Analysis",
        )
        click.echo(f"\nAll significant at p<0.01: {results['all_significant_p01']}")


@main.command("mis-analysis")
@click.option("--json", "as_json", is_flag=True, help="JSON output.")
def mis_analysis(as_json: bool) -> None:
    """Enumerate all maximal independent sets for tRNA enrichment."""
    from codon_topo.analysis.trna_evidence import maximal_independent_set_analysis

    results = maximal_independent_set_analysis()
    if as_json:
        _json_out(results)
    else:
        click.echo(f"MIS enumerated: {results['n_mis_size_ge2']}")
        click.echo(f"Median Stouffer p: {results['median_stouffer_p']:.4f}")
        click.echo(f"Worst-case Stouffer p: {results['worst_case_stouffer_p']:.4f}")
        click.echo(
            f"Fraction significant (p<0.05): {results['fraction_significant_p05']:.0%}"
        )
        click.echo(f"\n{results['interpretation']}")


@main.command("phylo-sensitivity")
@click.option("--json", "as_json", is_flag=True, help="JSON output.")
def phylo_sensitivity(as_json: bool) -> None:
    """Test topology avoidance robustness to phylogenetic clade exclusion."""
    from codon_topo.analysis.synbio_feasibility import (
        topology_avoidance_phylogenetic_sensitivity,
    )

    results = topology_avoidance_phylogenetic_sensitivity()
    if as_json:
        _json_out(results)
    else:
        lc = results["lineage_collapsed"]
        click.echo(
            f"Lineage-collapsed: {lc['n_events']} events, "
            f"{lc['depletion_fold']:.1f}x depletion, p={lc['hypergeom_p']:.2e}"
        )
        _try_rich_table(
            ["Excluded Clade", "n remaining", "p-value", "Significant"],
            [
                [
                    r["excluded_clade"],
                    str(r["n_events_remaining"]),
                    f"{r['hypergeom_p']:.2e}",
                    "YES" if r["significant_p05"] else "NO",
                ]
                for r in results["clade_exclusion"]
            ],
            title="Clade-Exclusion Sensitivity",
        )


@main.command("condlogit")
@click.option(
    "--max-orderings", default=720, help="Max orderings per table for order-averaging."
)
@click.option("--seed", default=DEFAULT_SEED, help="Random seed.")
@click.option("--json", "as_json", is_flag=True, help="Output as JSON.")
def condlogit(max_orderings: int, seed: int, as_json: bool) -> None:
    """Run the conditional logit model of reassignment choice (M1-M4)."""
    from codon_topo.analysis.evolutionary_simulation import (
        run_evolutionary_simulation_analysis,
    )

    result = run_evolutionary_simulation_analysis(
        max_orderings_per_table=max_orderings,
        seed=seed,
    )

    if as_json:
        _json_out(result)
        return

    click.echo("Conditional Logit Model Comparison")
    click.echo(f"  Tables: {result['n_tables']}, Events: {result['total_events']}")
    for name, aicc in result["aicc_ranking"]:
        fit = result["model_fits"][name]
        click.echo(
            f"  {name:16s}: AICc={aicc:.1f}  logL={fit['log_likelihood']:.1f}  "
            f"k={fit['n_params']}"
        )
    for test_name, lr in result["likelihood_ratio_tests"].items():
        click.echo(
            f"  LR {test_name}: stat={lr['lr_statistic']:.1f}, "
            f"df={lr['df']}, p={lr['p_value']:.2e}"
        )


@main.command("condlogit-restricted")
@click.option(
    "--max-orderings", default=720, help="Max orderings per table for order-averaging."
)
@click.option(
    "--max-trna",
    "max_trna_list",
    multiple=True,
    type=int,
    default=(1, 2, 3),
    help="Restrict candidates to delta_trna <= this value (repeatable).",
)
@click.option("--json", "as_json", is_flag=True, help="Output as JSON.")
def condlogit_restricted(
    max_orderings: int,
    max_trna_list: tuple[int, ...],
    as_json: bool,
) -> None:
    """Restricted-candidate sensitivity for the conditional logit.

    Refits M1-M4 on candidate sets filtered to moves whose target AA is
    already serviced by a codon at Hamming distance <= max-trna from the
    reassigned codon. The observed move is always retained.
    """
    from codon_topo.analysis.evolutionary_simulation import (
        run_restricted_candidate_sensitivity,
    )

    result = run_restricted_candidate_sensitivity(
        max_orderings_per_table=max_orderings,
        max_trna_thresholds=tuple(max_trna_list),
    )

    if as_json:
        _json_out(result)
        return

    click.echo("Restricted-Candidate Sensitivity (conditional logit)")
    click.echo(
        f"  Full set: ~{result['full_set_summary']['candidates_mean']:.0f} "
        f"candidates per choice set ({result['full_set_summary']['n_choice_sets']} sets)"
    )
    for k, block in result["by_max_trna"].items():
        cs = block["candidate_summary"]
        d = block["delta_aicc"]
        click.echo(
            f"  delta_trna<={k}: ~{cs['candidates_mean']:.0f} candidates, "
            f"obs retained {cs['observed_in_filtered_set']}/{cs['observed_total']}"
        )
        if "M1_to_M3" in d:
            click.echo(f"    ΔAICc(M1→M3 Q_6)    = {d['M1_to_M3']:.1f}")
        if "M2_to_M3" in d:
            click.echo(f"    ΔAICc(M2→M3 Q_6)    = {d['M2_to_M3']:.1f}")
        if "M3_to_M4" in d:
            click.echo(f"    ΔAICc(M3→M4)        = {d['M3_to_M4']:.1f}")
        if "M1_to_M3_k43" in d:
            click.echo(f"    ΔAICc(M1→M3 H(3,4)) = {d['M1_to_M3_k43']:.1f}")


@main.command("condlogit-heavy")
@click.option(
    "--output-dir",
    required=True,
    help="Directory to write condlogit_clade_sensitivity.json + "
    "condlogit_restricted_candidate.json.",
)
@click.option(
    "--max-orderings", default=720, help="Max orderings per table for order-averaging."
)
def condlogit_heavy(output_dir: str, max_orderings: int) -> None:
    """Run the two memory-heavy conditional-logit sensitivity analyses in a
    fresh Python process and write their JSON artefacts to ``--output-dir``.

    This subcommand is invoked by ``codon-topo all`` via ``subprocess.run``
    so that the ~750 MB feature-bundle allocations do not inherit the
    ~2.5 GB accumulated resident set of the main pipeline process. On
    macOS the inherited-VM fork of joblib Pool workers by the main
    pipeline triggers jetsam SIGKILL under memory pressure; running these
    two steps in a fresh interpreter drops the baseline from ~2.5 GB to
    ~50 MB and keeps the workers under the jetsam threshold.

    Writes:
      - ``condlogit_clade_sensitivity.json`` (7 clade-exclusion regimes)
      - ``condlogit_restricted_candidate.json`` (three delta_trna cuts)
    """
    from pathlib import Path as _Path

    from codon_topo.analysis.evolutionary_simulation import (
        run_clade_exclusion_sensitivity,
        run_restricted_candidate_sensitivity,
    )

    out = _Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    click.echo("[condlogit-heavy] Running clade-exclusion sensitivity (7 regimes)...")
    clade_result = run_clade_exclusion_sensitivity(
        max_orderings_per_table=max_orderings,
    )
    _write_json(out / "condlogit_clade_sensitivity.json", clade_result)
    click.echo(f"[condlogit-heavy] Wrote {out / 'condlogit_clade_sensitivity.json'}")

    click.echo(
        "[condlogit-heavy] Running restricted-candidate sensitivity "
        "(delta_trna<=1,2,3)..."
    )
    restricted_result = run_restricted_candidate_sensitivity(
        max_orderings_per_table=max_orderings,
    )
    _write_json(out / "condlogit_restricted_candidate.json", restricted_result)
    click.echo(f"[condlogit-heavy] Wrote {out / 'condlogit_restricted_candidate.json'}")


@main.command("topology-avoidance-k43")
@click.option("--seed", default=DEFAULT_SEED, help="Random seed.")
@click.option("--json", "as_json", is_flag=True, help="Output as JSON.")
def topology_avoidance_k43_cmd(seed: int, as_json: bool) -> None:
    """Topology avoidance test under K4^3 (encoding-independent) adjacency."""
    from codon_topo.analysis.synbio_feasibility import topology_avoidance_k43

    result = topology_avoidance_k43(seed=seed)

    if as_json:
        _json_out(result)
        return

    click.echo("Topology Avoidance (K4^3, encoding-independent)")
    click.echo(
        f"  Observed: {result['observed_breaks']}/{result['observed_total']} "
        f"topology-breaking ({result['rate_observed']:.1%})"
    )
    click.echo(
        f"  Possible: {result['possible_breaks']}/{result['possible_total']} "
        f"({result['rate_possible']:.1%})"
    )
    click.echo(f"  Depletion: {result['depletion_fold']:.1f}x")
    click.echo(
        f"  RR={result['risk_ratio']:.2f} "
        f"(95% CI [{result['risk_ratio_ci_95'][0]:.2f}, "
        f"{result['risk_ratio_ci_95'][1]:.2f}])"
    )
    click.echo(f"  Hypergeometric p={result['hypergeom_p']:.2e}")
    click.echo(f"  Permutation p={result['permutation_p']:.4f}")


@main.command("codonsafe")
@click.option("--json", "as_json", is_flag=True, help="Output as JSON.")
def codonsafe_cmd(as_json: bool) -> None:
    """Run the CodonSafe cross-study reanalysis of genome recoding datasets.

    Requires the [codonsafe] extra: pip install -e '.[codonsafe]'
    and raw data files in data/codonsafe/.
    """
    try:
        from codon_topo.analysis.codonsafe.run_analyses import main as cs_main
    except ImportError as e:
        click.echo(
            f"Error: {e}\nInstall codonsafe dependencies: pip install -e '.[codonsafe]'"
        )
        raise SystemExit(1)

    cs_main()

    if as_json:
        # Read and emit the combined stats JSON
        stats_path = Path("output/codonsafe/analysis_stats.json")
        if stats_path.exists():
            click.echo(stats_path.read_text())


@main.command("walsh-2adic")
@click.option("--n", "n_samples", default=2000, help="Null draws per test.")
@click.option("--seed", default=DEFAULT_SEED, help="Random seed.")
@click.option("--json", "as_json", is_flag=True, help="Output as JSON.")
def walsh_2adic_cmd(n_samples: int, seed: int, as_json: bool) -> None:
    """Walsh-Hadamard / 2-adic spectral probe of synonymous-block structure.

    Runs Test 2 (block-size null), Test 2b (wobble-box-preserving label
    permutation), and Test 3 (24-encoding invariance). Reports the spectral
    depth (encoding-invariant) and the wobble-free label-spectrum fraction
    S = E_high / E_total. See `src/codon_topo/analysis/walsh_2adic.py`.
    """
    from codon_topo.analysis.walsh_2adic import (
        compute_walsh_results,
        label_spectrum_null_test,
    )

    result = compute_walsh_results(
        n_null_main=n_samples,
        n_null_per_encoding=n_samples,
        n_null_label_perm=n_samples,
        seed=seed,
    )
    label = label_spectrum_null_test(n_null=n_samples, seed=seed)
    result["label_spectrum_S"] = label
    if as_json:
        _json_out(result)
        return
    t2 = result["test2_block_size_null"]
    t2b = result["test2b_label_permutation_null"]
    t3 = result["test3_encoding_invariance"]
    click.echo(
        f"Block-size null: depth={t2['observed_depth']}, z={t2['z_score']:.2f}, "
        f"frac_below={t2['fraction_null_at_or_below_observed']:.4f}"
    )
    click.echo(
        f"Wobble-box null: depth invariant at {t2b['observed_depth']} "
        f"(null std={t2b['null_std']:.4f})"
    )
    click.echo(
        f"24-encoding sweep: depths_unique={t3['observed_depths_unique']}, "
        f"all_z_negative={t3['all_z_negative']}"
    )
    click.echo(
        f"Label-spectrum S: observed={label['observed_S']:.6f}, "
        f"null_mean={label['null_mean']:.4f}, z={label['z_score']:.2f}"
    )


@main.command("protsub")
@click.option("--n", "n_samples", default=10_000, help="Monte Carlo sample size.")
@click.option("--seed", default=DEFAULT_SEED, help="Random seed.")
@click.option("--json", "as_json", is_flag=True, help="Output as JSON.")
def protsub_cmd(n_samples: int, seed: int, as_json: bool) -> None:
    """ProtSub matrix Monte Carlo (structure-aware code-dependent metric).

    ProtSub is alignment-derived (MSA of 2,320 Pfam domains) and therefore
    code-dependent in the sense of Di Giulio (2001). Reported as a
    structure-aware sensitivity check on top of the four primary
    physicochemical metrics.
    """
    from codon_topo.analysis.coloring_optimality import monte_carlo_null

    result = monte_carlo_null(
        n_samples=n_samples, seed=seed, null_type="freeland_hurst", metric="protsub"
    )
    if as_json:
        _json_out(result)
        return
    click.echo(
        f"ProtSub: observed={result['observed_score']:.1f}, "
        f"null_mean={result['null_mean']:.1f}±{result['null_std']:.1f}, "
        f"z={result['effect_size_z']:.2f}, "
        f"quantile={result['quantile_of_observed']:.2f}%, "
        f"p={result['p_value_conservative']:.4f}"
    )


@main.command("metric-corr")
@click.option(
    "--n-null", default=2000, help="Number of null codes for F-score correlations."
)
@click.option("--seed", default=DEFAULT_SEED, help="Random seed.")
@click.option("--json", "as_json", is_flag=True, help="Output as JSON.")
def metric_corr_cmd(n_null: int, seed: int, as_json: bool) -> None:
    """Pairwise and partial Spearman correlations of the five distance metrics.

    Computes correlations at two levels: (a) across the 190 unordered AA
    pairs of the 20 standard amino acids; (b) across n_null random codes
    drawn from the Freeland-Hurst null, with partial correlations
    controlling for the other metrics.
    """
    from codon_topo.analysis.coloring_optimality import (
        aa_pair_metric_correlations,
        code_level_F_score_correlations,
    )

    pairs = aa_pair_metric_correlations()
    codes = code_level_F_score_correlations(n_null=n_null, seed=seed)
    out = {"aa_pair_level": pairs, "code_level": codes}
    if as_json:
        _json_out(out)
        return
    click.echo(
        f"AA-pair Spearman ρ: off-diagonal "
        f"[{pairs['offdiagonal_min']:.2f}, {pairs['offdiagonal_max']:.2f}], "
        f"median {pairs['offdiagonal_median']:.2f}"
    )
    click.echo(
        f"Code-level F-score Spearman ρ: off-diagonal "
        f"[{codes['F_offdiagonal_min']:.2f}, {codes['F_offdiagonal_max']:.2f}], "
        f"median {codes['F_offdiagonal_median']:.2f}"
    )
    click.echo(
        f"Partial correlations: "
        f"{codes['n_partials_independent']} of {codes['n_partial_pairs_total']} "
        f"pairs independent (|partial ρ| < 0.1)"
    )


@main.command("all")
@click.option("--output-dir", default="./output", help="Directory for CSV/JSON output.")
@click.option("--seed", default=DEFAULT_SEED, help="Random seed for Monte Carlo.")
@click.option("--n", "n_samples", default=10_000, help="Monte Carlo sample size.")
def run_all(output_dir: str, seed: int, n_samples: int) -> None:
    """Run all analyses and write reports to output directory.

    Generates per-analysis JSON files plus a consolidated manuscript_stats.json
    that the Typst manuscript reads for all inline statistics.
    """
    from codon_topo.core.filtration import analyze_filtration
    from codon_topo.core.homology import disconnection_catalogue
    from codon_topo.core.genetic_codes import all_table_ids, get_code
    from codon_topo.analysis.coloring_optimality import (
        aa_pair_metric_correlations,
        code_level_F_score_correlations,
        monte_carlo_null,
        cross_table_optimality,
        multi_metric_sensitivity,
        per_table_optimality,
        rho_robustness_sweep,
        score_decomposition_by_position,
        stop_penalty_sensitivity,
        codon_usage_vs_local_mismatch,
        mechanistic_discriminant_test,
    )
    from codon_topo.analysis.walsh_2adic import (
        compute_walsh_results,
        label_spectrum_null_test,
    )
    from codon_topo.analysis.reassignment_db import (
        build_reassignment_db,
        hamming_path_lengths,
        bit_position_bias_weighted as _bit_bias_weighted,
        directionality_summary,
    )
    from codon_topo.analysis.trna_evidence import (
        trna_duplication_correlation_test,
        maximal_independent_set_analysis,
        fisher_exact_per_pairing,
    )
    from codon_topo.analysis.depth_calibration import (
        compute_correlation,
        depth_calibration_table,
    )
    from codon_topo.analysis.synbio_feasibility import (
        topology_avoidance_test as _topo_avoidance,
        topology_avoidance_k43,
        topology_avoidance_k43_phylogenetic_sensitivity,
        topology_avoidance_phylogenetic_sensitivity,
    )
    from codon_topo.analysis.cosmic_query import fano_predictions_for_kras
    from codon_topo.reports.claim_hierarchy import CLAIM_HIERARCHY
    from codon_topo.reports.catalogue import build_catalogue

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    n_steps = 24
    step = 0

    def _step(label: str) -> None:
        nonlocal step
        step += 1
        click.echo(f"  [{step}/{n_steps}] {label}")

    click.echo(f"Running all analyses (seed={seed}, n={n_samples})...")

    # 1. Filtration across all tables
    _step("Filtration...")
    filt_results = []
    for tid in all_table_ids():
        code = get_code(tid)
        r = analyze_filtration(code)
        r["table_id"] = tid
        filt_results.append(r)
    _write_json(out / "filtration.json", filt_results)

    # 2. Disconnection catalogue
    _step("Disconnections...")
    disc_results = []
    for tid in all_table_ids():
        code = get_code(tid)
        cat = disconnection_catalogue(code)
        for entry in cat:
            entry["table_id"] = tid
        disc_results.extend(cat)
    _write_json(out / "disconnections.json", disc_results)

    # 3. Coloring optimality + multi-metric + rho sweep + per-table + decomposition
    _step(f"Coloring optimality (n={n_samples})...")
    # include_null_samples: persist the full 10k null draws so figures
    # (Fig 2A, FigA) can plot the empirical histogram instead of a
    # Gaussian approximation of it. Only the primary MC block gets them;
    # per-table / rho / encoding-sweep would otherwise multiply the JSON
    # size by ~n_samples.
    coloring_result = monte_carlo_null(
        n_samples=n_samples, seed=seed, include_null_samples=True
    )
    cross_table = cross_table_optimality()

    _step("Multi-metric sensitivity (n=1000)...")
    metric_results = multi_metric_sensitivity(
        n_samples=min(n_samples, 10_000), seed=seed
    )

    _step("Classical Haig-Hurst AA-permutation null (B2b sensitivity)...")
    hh_aa_null = {
        "method": (
            "Classical Haig-Hurst 1991 / Freeland-Hurst 1998 AA-permutation "
            "null. Preserves stop positions and the unlabeled codon-family "
            "partition (which codons decode together); permutes the 20 AA "
            "labels uniformly across the 20 sense-codon families (20! ~ "
            "2.4e18 possible codes). Per-named-AA codon count is NOT "
            "preserved: e.g., Trp's 1-codon family can become a 6-codon "
            "family if its label is swapped with Leu's. Reported as a "
            "sensitivity companion to the primary quartet-pattern shuffle."
        ),
        "seed": seed,
        "n_samples": min(n_samples, 10_000),
        "per_metric": {},
    }
    for _metric in ("grantham", "miyata", "polar_requirement", "kyte_doolittle"):
        _r = monte_carlo_null(
            n_samples=min(n_samples, 10_000),
            seed=seed,
            null_type="haig_hurst_aa",
            metric=_metric,
        )
        hh_aa_null["per_metric"][_metric] = {
            "observed_score": _r["observed_score"],
            "null_mean": _r["null_mean"],
            "null_std": _r["null_std"],
            "null_min": _r["null_min"],
            "null_max": _r["null_max"],
            "z": _r["effect_size_z"],
            "quantile_of_observed": _r["quantile_of_observed"],
            "n_beaten_observed": _r["n_beaten_observed"],
            "p_value_raw": _r["p_value_raw"],
            "p_value_conservative": _r["p_value_conservative"],
        }
    _write_json(out / "coloring_optimality_hh_aa_null.json", hh_aa_null)

    _step("ProtSub Monte Carlo (code-dependent sensitivity check)...")
    protsub_result = monte_carlo_null(
        n_samples=min(n_samples, 10_000),
        seed=seed,
        null_type="freeland_hurst",
        metric="protsub",
    )
    _write_json(out / "protsub_optimality.json", protsub_result)

    _step("Metric correlation matrices (190 AA pairs + 2000 random codes)...")
    metric_corr_pairs = aa_pair_metric_correlations()
    _write_json(out / "metric_correlations.json", metric_corr_pairs)
    metric_corr_codes = code_level_F_score_correlations(n_null=2000, seed=seed)
    _write_json(out / "metric_F_correlations.json", metric_corr_codes)

    _step("Walsh-Hadamard / 2-adic spectral probe...")
    walsh_result = compute_walsh_results(
        n_null_main=min(n_samples, 2000),
        n_null_per_encoding=min(n_samples, 1500),
        n_null_label_perm=min(n_samples, 2000),
        seed=seed,
    )
    _write_json(out / "walsh_2adic.json", walsh_result)

    _step("Walsh label-spectrum (wobble-free fraction S)...")
    walsh_label = label_spectrum_null_test(n_null=min(n_samples, 10_000), seed=seed)
    _write_json(out / "walsh_label_spectrum.json", walsh_label)

    _step("Rho robustness sweep...")
    rho_result = rho_robustness_sweep(n_samples=min(n_samples, 10_000), seed=seed)

    _step("Per-table optimality...")
    # Match the coloring-optimality MC convention (CLAUDE.md: n=10,000).
    pertable_result = per_table_optimality(n_samples=min(n_samples, 10_000), seed=seed)

    _step("Per-table proximity audit (standard-code-proximity sensitivity)...")
    from codon_topo.analysis.coloring_optimality import (
        encoding_sensitivity_of_optimality,
        per_table_proximity_audit,
    )

    proximity_audit = per_table_proximity_audit(
        n_samples=min(n_samples, 1_000), seed=seed
    )

    _step("Coloring 24-encoding sensitivity sweep...")
    coloring_enc_sweep = encoding_sensitivity_of_optimality(
        n_samples=min(n_samples, 10_000),
        seed=seed,
        null_type="freeland_hurst",
    )

    decomp_result = score_decomposition_by_position()

    stop_results = stop_penalty_sensitivity(n_samples=1_000, seed=seed)
    _write_json(
        out / "coloring_optimality.json",
        {
            "monte_carlo": coloring_result,
            "cross_table": cross_table,
            "multi_metric": metric_results,
            "rho_sweep": rho_result,
            "per_table": pertable_result,
            "per_table_proximity_audit": proximity_audit,
            "encoding_sensitivity": coloring_enc_sweep,
            "decomposition": decomp_result,
            "stop_penalty": stop_results,
        },
    )

    # 7. Reassignment DB + bit bias
    _step("Reassignment analysis...")
    db = build_reassignment_db()
    paths = hamming_path_lengths()
    bit_bias_mito = _bit_bias_weighted("mitochondrial")
    dir_summary = directionality_summary()
    _write_json(
        out / "reassignment_analysis.json",
        {
            "n_events": len(db),
            "hamming_paths": paths,
            "bit_bias_mitochondrial": bit_bias_mito,
            "directionality": dir_summary,
        },
    )

    # 8. Topology avoidance (Q6 + K4^3 + definitions audit + 24-encoding sweep)
    _step("Topology avoidance (Q6 + K4^3)...")
    topo_q6 = _topo_avoidance()
    topo_k43 = topology_avoidance_k43(seed=seed)
    _step("Topology definitions audit (2x2 sensitivity)...")
    from codon_topo.analysis.synbio_feasibility import (
        topology_definitions_audit,
        topology_avoidance_q6_encoding_sweep,
        topology_denominator_sensitivity,
    )

    topo_audit = topology_definitions_audit()
    _step("Topology Q_6 24-encoding sweep...")
    topo_enc_sweep = topology_avoidance_q6_encoding_sweep()
    topo_denom = topology_denominator_sensitivity()
    _write_json(
        out / "topology_avoidance.json",
        {
            "Q6": topo_q6,
            "K43": topo_k43,
            "definitions_audit": topo_audit,
            "Q6_encoding_sweep": topo_enc_sweep,
            "denominator_sensitivity": topo_denom,
        },
    )

    # 9. tRNA evidence (with MIS analysis + topology-breaking subset)
    _step("tRNA evidence + MIS analysis + topology-breaking subset...")
    trna_result = trna_duplication_correlation_test()
    fisher_result = fisher_exact_per_pairing()
    mis_result = maximal_independent_set_analysis()
    from codon_topo.analysis.trna_evidence import topology_breaking_subset_test

    trna_topo_subset = topology_breaking_subset_test()
    _write_json(
        out / "trna_evidence.json",
        {
            "sign_test": trna_result,
            "fisher_stouffer": fisher_result,
            "mis_analysis": mis_result,
            "topology_breaking_subset": trna_topo_subset,
        },
    )

    # 10. Phylogenetic sensitivity (Q_6 / new-disconnection + H(3,4) / Delta beta_0)
    _step("Phylogenetic sensitivity (Q_6 + H(3,4))...")
    phylo_result = topology_avoidance_phylogenetic_sensitivity()
    _write_json(out / "phylogenetic_sensitivity.json", phylo_result)
    phylo_k43 = topology_avoidance_k43_phylogenetic_sensitivity(seed=seed)
    _write_json(out / "phylogenetic_sensitivity_k43.json", phylo_k43)

    # 11. Evolutionary simulation (conditional logit)
    _step("Conditional logit models (M1-M4)...")
    from codon_topo.analysis.evolutionary_simulation import (
        run_evolutionary_simulation_analysis,
    )

    evosim_result = run_evolutionary_simulation_analysis(seed=seed)
    _write_json(out / "evolutionary_simulation.json", evosim_result)

    # 11b+c. Conditional-logit clade-exclusion + restricted-candidate
    # sensitivities.
    #
    # These two analyses call fit_all_models() 7 and 3 times respectively
    # (~750 MB feature bundle each) and are the two heaviest steps in the
    # pipeline. Running them inline in the main pipeline process — which by
    # this point has accumulated ~2.5 GB resident from steps 1-20 — causes
    # macOS jetsam to SIGKILL the joblib Pool worker children under memory
    # pressure. We therefore invoke them as a fresh subprocess so they
    # start from a ~50 MB Python baseline; peak in-subprocess RAM is
    # ~800 MB, well under the jetsam threshold.
    _step(
        "Conditional logit sensitivity (clade-exclusion + restricted-candidate) "
        "in fresh subprocess..."
    )
    import subprocess as _sp

    _sp.run(
        ["codon-topo", "condlogit-heavy", "--output-dir", str(out)],
        check=True,
    )
    with open(out / "condlogit_clade_sensitivity.json") as _fh:
        condlogit_clade = json.load(_fh)
    with open(out / "condlogit_restricted_candidate.json") as _fh:
        condlogit_restricted_result = json.load(_fh)

    # 12. Depth calibration
    _step("Depth calibration...")
    corr = compute_correlation()
    cal_table = depth_calibration_table()
    _write_json(
        out / "depth_calibration.json",
        {
            "correlation": corr,
            "calibration_points": cal_table,
        },
    )

    # 13. KRAS-Fano predictions (offline)
    _step("KRAS-Fano predictions (offline)...")
    fano_preds = fano_predictions_for_kras()
    _write_json(out / "kras_fano.json", {"mode": "offline", "predictions": fano_preds})

    # 14. Claim hierarchy + catalogue + exploratory tests
    _step("Claims, catalogue & exploratory tests...")
    claims_data = [
        {
            "id": c.id,
            "status": c.status.value,
            "statement": c.statement,
            "p_value": c.evidence_p_value,
            "null_model": c.null_model,
            "sample_size": c.sample_size,
            "publishable": c.is_publishable,
        }
        for c in CLAIM_HIERARCHY
    ]
    catalogue = build_catalogue()
    catalogue_data = [
        {
            "id": p.id,
            "claim": p.claim,
            "workstream": p.workstream,
            "status": p.status,
            "evidence_strength": p.evidence_strength,
            "p_value": p.p_value,
        }
        for p in catalogue
    ]
    cub_result = codon_usage_vs_local_mismatch()
    mech_result = mechanistic_discriminant_test()
    _write_json(out / "claims.json", claims_data)
    _write_json(out / "catalogue.json", catalogue_data)
    _write_json(
        out / "exploratory_tests.json",
        {"cub_correlation": cub_result, "mechanistic_discriminant": mech_result},
    )

    # ================================================================
    # Generate consolidated manuscript_stats.json
    # Every number cited in manuscript.typ should come from here.
    # ================================================================
    click.echo("\n  Generating manuscript_stats.json...")

    # Build per-metric lookup from multi-metric results
    mm = {
        m["metric"].lower().replace(" ", "_").replace("-", "_"): m
        for m in metric_results["per_metric"]
    }

    # Best conditional logit model
    best_model = (
        evosim_result["aicc_ranking"][0][0]
        if evosim_result["aicc_ranking"]
        else "M3_phys_topo"
    )
    best_fit = evosim_result["model_fits"].get(best_model, {})
    lr_tests = evosim_result.get("likelihood_ratio_tests", {})

    manuscript_stats = {
        "_generated_by": "codon-topo all",
        "_seed": seed,
        "_n_samples": n_samples,
        "_version": "0.6.1",
        # Section 3.1: Cross-metric coloring optimality
        "coloring": {
            "observed_score": coloring_result["observed_score"],
            "per_table_proximity_audit": proximity_audit,
            "null_mean": coloring_result["null_mean"],
            "null_std": coloring_result["null_std"],
            "quantile": coloring_result["quantile_of_observed"],
            "p_value": coloring_result["p_value_conservative"],
            "n_samples": n_samples,
            # Classical Haig-Hurst AA-permutation null (Supplement §S2.1
            # sensitivity companion to the primary quartet-pattern shuffle).
            "haig_hurst_aa_null": hh_aa_null,
            # 24-encoding sensitivity summary (Supplement §S2)
            "encoding_sensitivity": {
                "n_encodings": coloring_enc_sweep["n_encodings"],
                "null_type": coloring_enc_sweep["null_type"],
                "quantile_min": coloring_enc_sweep["quantile_range"][0],
                "quantile_max": coloring_enc_sweep["quantile_range"][1],
                "mean_quantile": coloring_enc_sweep["mean_quantile"],
                "all_significant_p05": coloring_enc_sweep[
                    "all_encodings_significant_p05"
                ],
                "p_max": max(
                    r["p_value_conservative"]
                    for r in coloring_enc_sweep["per_encoding"]
                ),
                "p_min": min(
                    r["p_value_conservative"]
                    for r in coloring_enc_sweep["per_encoding"]
                ),
            },
        },
        "metrics": {
            name: {
                "observed": m.get("observed_score"),
                "null_mean": m.get("null_mean"),
                "null_std": m.get("null_std"),
                "z": m.get("effect_size_z"),
                "p": m.get("p_value_conservative"),
                "quantile": m.get("quantile"),
                "improvement_pct": (
                    m.get("improvement_pct")
                    if m.get("improvement_pct") is not None
                    else (
                        100.0 * (m["null_mean"] - m["observed_score"]) / m["null_mean"]
                        if m.get("null_mean") not in (None, 0)
                        and m.get("observed_score") is not None
                        else None
                    )
                ),
            }
            for name, m in mm.items()
        },
        # Section 3.2: Rho robustness
        "rho_sweep": {
            "per_rho": rho_result["per_rho"],
            "all_significant_p01": rho_result["all_significant_p01"],
        },
        # Section 3.3: Per-table optimality
        "per_table": {
            "n_significant_bh": pertable_result["n_significant_p05_bh"],
            "n_tables": pertable_result["n_tables"],
            "mean_quantile": pertable_result["mean_quantile"],
            "per_table": pertable_result["per_table"],
            **_per_table_bh_disaggregation(pertable_result, db),
        },
        # Section 3.1: Score decomposition
        "decomposition": {
            "position_fractions": decomp_result["position_fractions"],
            "total_score": decomp_result["total_score"],
        },
        # Section 3.4: Topology avoidance (Q6)
        "topology_avoidance_q6": (
            lambda tq: {
                "observed_breaks": tq["observed_creates_disc"],
                "observed_total": tq["observed_total"],
                "rate_observed": tq["rate_observed"],
                "possible_breaks": tq["possible_creates_disc"],
                "possible_total": tq["possible_total"],
                "rate_possible": tq["rate_possible"],
                "depletion_fold": _safe_depletion(
                    tq["rate_possible"], tq["rate_observed"]
                ),
                "risk_ratio": _safe_risk_ratio(
                    tq["rate_observed"], tq["rate_possible"]
                ),
                "risk_ratio_ci_95": _log_normal_ci_95(
                    tq["observed_creates_disc"],
                    tq["observed_total"],
                    tq["possible_creates_disc"],
                    tq["possible_total"],
                ),
                "hypergeom_p": tq["hypergeom_p"],
                "permutation_p": tq["permutation_p"],
            }
        )(topo_q6),
        # Section 3.4: Topology avoidance (K4^3, encoding-independent)
        "topology_avoidance_k43": {
            "observed_breaks": topo_k43["observed_breaks"],
            "observed_total": topo_k43["observed_total"],
            "rate_observed": topo_k43["rate_observed"],
            "possible_breaks": topo_k43["possible_breaks"],
            "possible_total": topo_k43["possible_total"],
            "rate_possible": topo_k43["rate_possible"],
            "depletion_fold": topo_k43["depletion_fold"],
            "risk_ratio": topo_k43["risk_ratio"],
            "risk_ratio_ci_95": topo_k43["risk_ratio_ci_95"],
            "hypergeom_p": topo_k43["hypergeom_p"],
            "permutation_p": topo_k43["permutation_p"],
        },
        # Section 3.4 supplementary: 2x2 audit + 24-encoding sweep + denominator sensitivity
        "topology_audit": topo_audit,
        "topology_q6_encoding_sweep": topo_enc_sweep,
        # Section 3.5: Conditional logit
        "condlogit": {
            "n_tables": evosim_result["n_tables"],
            "total_events": evosim_result["total_events"],
            "aicc_ranking": evosim_result["aicc_ranking"],
            "model_fits": evosim_result["model_fits"],
            "lr_tests": lr_tests,
            "encoding_robustness": evosim_result.get("encoding_robustness", {}),
            # Clade-exclusion sensitivity aggregates (used by main-text §3.5;
            # full per-regime rows live in condlogit_clade_sensitivity.json).
            "clade_exclusion": {
                "delta_M1_M3_min": condlogit_clade.get("delta_M1_M3_min"),
                "delta_M1_M3_median": condlogit_clade.get("delta_M1_M3_median"),
                "delta_M1_M3_max": condlogit_clade.get("delta_M1_M3_max"),
                "delta_M2_M3_min": condlogit_clade.get("delta_M2_M3_min"),
                "delta_M2_M3_median": condlogit_clade.get("delta_M2_M3_median"),
                "delta_M2_M3_max": condlogit_clade.get("delta_M2_M3_max"),
                "all_M3_favored_over_M1": condlogit_clade.get("all_M3_favored_over_M1"),
                "all_M3_favored_over_M2": condlogit_clade.get("all_M3_favored_over_M2"),
                "n_regimes": len(condlogit_clade.get("rows", [])),
                # Per-regime rows for the §S7 dynamic table.
                "rows": condlogit_clade.get("rows", []),
            },
            # Restricted-candidate-set sensitivity
            "restricted_candidate": condlogit_restricted_result,
            # Spearman rho between phys and topo features across the candidate
            # landscape (main-text abstract + conclusion + §4.2)
            "phys_topo_rho": (
                evosim_result.get("diagnostics", {})
                .get("phys_topo_correlation", {})
                .get("spearman_rho")
            ),
        },
        # Section 3.6: tRNA enrichment
        "trna": {
            "n_pairings": fisher_result.get("n_pairings"),
            "stouffer_z": fisher_result.get("stouffer_z"),
            "stouffer_p": fisher_result.get("stouffer_p"),
            "mis_worst_p": mis_result.get("worst_case_stouffer_p"),
            "mis_best_p": mis_result.get("best_case_stouffer_p"),
            "mis_median_p": mis_result.get("median_stouffer_p"),
            "n_mis": mis_result.get("n_mis_size_ge2"),
            "n_mis_total": mis_result.get("n_mis_total"),
            "mis_frac_significant_p05": mis_result.get("fraction_significant_p05"),
            "topology_breaking_p": trna_topo_subset.get("stouffer_p_topology_breaking"),
            "topology_breaking_z": trna_topo_subset.get("stouffer_z_topology_breaking"),
            "topology_breaking_n": trna_topo_subset.get("n_pairings"),
        },
        # Phylogenetic sensitivity (Supplement §S9)
        # Two cells reported: Q_6 / new-disconnection (legacy `phylo` key,
        # matches existing prose/tables) and H(3,4) / Delta beta_0 (`phylo_k43`,
        # the primary-cell H(3,4) companion added in R3-B3).
        "phylo": {
            "all_significant": phylo_result["all_clade_exclusions_significant"],
            "lineage_collapsed": phylo_result["lineage_collapsed"],
            "clade_exclusion": phylo_result["clade_exclusion"],
        },
        "phylo_k43": {
            "adjacency": phylo_k43["adjacency"],
            "definition": phylo_k43["definition"],
            "all_significant": phylo_k43["all_clade_exclusions_significant"],
            "lineage_collapsed": phylo_k43["lineage_collapsed"],
            "clade_exclusion": phylo_k43["clade_exclusion"],
        },
        # Reassignment DB (Supplement §S11)
        "reassignment": {
            "n_events": len(db),
            "n_dedup": topo_q6["observed_total"],
            # Per-table count of raw reassignment events (used by §S11 table)
            "per_table_counts": _per_table_reassignment_counts(db),
        },
        # ProtSub (structure-aware code-dependent sensitivity check; §3.1, §S15)
        "protsub": {
            "observed": protsub_result["observed_score"],
            "null_mean": protsub_result["null_mean"],
            "null_std": protsub_result["null_std"],
            "z": protsub_result["effect_size_z"],
            "p": protsub_result["p_value_conservative"],
            "quantile": protsub_result["quantile_of_observed"],
        },
        # Metric correlations (§S15)
        "metric_correlations": {
            "aa_pair_min": metric_corr_pairs["offdiagonal_min"],
            "aa_pair_median": metric_corr_pairs["offdiagonal_median"],
            "aa_pair_max": metric_corr_pairs["offdiagonal_max"],
            "code_level_min": metric_corr_codes["F_offdiagonal_min"],
            "code_level_median": metric_corr_codes["F_offdiagonal_median"],
            "code_level_max": metric_corr_codes["F_offdiagonal_max"],
            "n_partials_independent": metric_corr_codes["n_partials_independent"],
            "n_partial_pairs_total": metric_corr_codes["n_partial_pairs_total"],
        },
        # Walsh-Hadamard / 2-adic (§4.5, §S16)
        "walsh": {
            "spectral_depth": walsh_result["test2_block_size_null"]["observed_depth"],
            "block_size_z": walsh_result["test2_block_size_null"]["z_score"],
            "block_size_null_mean": walsh_result["test2_block_size_null"]["null_mean"],
            "block_size_null_std": walsh_result["test2_block_size_null"]["null_std"],
            "label_perm_invariant": walsh_result["test2b_label_permutation_null"][
                "is_algebraic_invariant"
            ],
            "encoding_invariant": walsh_result["test3_encoding_invariance"][
                "invariant"
            ],
            "encoding_z_min": walsh_result["test3_encoding_invariance"]["z_score_min"],
            "encoding_z_max": walsh_result["test3_encoding_invariance"]["z_score_max"],
            # Label-spectrum (§S16.4)
            "label_spectrum_S": walsh_label["observed_S"],
            "label_spectrum_null_mean": walsh_label["null_mean"],
            "label_spectrum_null_std": walsh_label["null_std"],
            "label_spectrum_z": walsh_label["z_score"],
            "label_spectrum_n_distinct_null": walsh_label["n_distinct_null_S_values"],
        },
        # Pipeline metadata
        "pipeline": {
            "n_claims": len(claims_data),
            "n_supported": sum(1 for c in claims_data if c["status"] == "supported"),
            "n_tests": "dynamic",  # filled by pytest count
        },
    }
    _write_json(out / "manuscript_stats.json", manuscript_stats)

    # Summary
    click.echo(f"\nDone. {len(list(out.glob('*.json')))} JSON files written to {out}/")
    for m in metric_results["per_metric"]:
        click.echo(
            f"  {m['metric']:20s}  quantile={m['quantile']:.1f}%  p={m['p_value_conservative']:.4f}"
        )
    click.echo(
        f"  ProtSub (sensitivity check): quantile="
        f"{protsub_result['quantile_of_observed']:.2f}%, "
        f"p={protsub_result['p_value_conservative']:.4f}"
    )
    click.echo(
        f"  Walsh spectral depth: {walsh_result['test2_block_size_null']['observed_depth']} "
        f"(z={walsh_result['test2_block_size_null']['z_score']:.2f}); "
        f"S = {walsh_label['observed_S']:.4f} (z={walsh_label['z_score']:.2f})"
    )
    click.echo(f"  MIS worst-case p: {mis_result.get('worst_case_stouffer_p', 'N/A')}")
    click.echo(
        f"  Phylo clade-exclusion: all significant = {phylo_result['all_clade_exclusions_significant']}"
    )
    click.echo(
        f"  Condlogit best model: {best_model} (AICc={best_fit.get('aicc', 'N/A')})"
    )
    click.echo("  manuscript_stats.json written for Typst dynamic rendering")


def _safe_depletion(rate_possible: float, rate_observed: float):
    return rate_possible / rate_observed if rate_observed > 0 else None


def _safe_risk_ratio(rate_observed: float, rate_possible: float):
    return rate_observed / rate_possible if rate_possible > 0 else None


def _log_normal_ci_95(x: int, n: int, K: int, N: int) -> list:
    """95% log-normal CI for the risk ratio (x/n)/(K/N)."""
    import math

    if x <= 0 or n <= 0 or K <= 0 or N <= 0 or n == x or N == K:
        return [None, None]
    log_rr = math.log((x / n) / (K / N))
    se_log_rr = math.sqrt((1 / x) - (1 / n) + (1 / K) - (1 / N))
    return [math.exp(log_rr - 1.96 * se_log_rr), math.exp(log_rr + 1.96 * se_log_rr)]


def _sanitize_nonfinite(obj):
    """Recursively convert NaN/Inf floats to None so the serialized JSON
    parses under strict readers (R jsonlite, Typst json(), etc.)."""
    import math

    if isinstance(obj, float):
        return None if not math.isfinite(obj) else obj
    if isinstance(obj, np.floating):
        return None if not np.isfinite(obj) else float(obj)
    if isinstance(obj, dict):
        return {k: _sanitize_nonfinite(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_sanitize_nonfinite(v) for v in obj]
    return obj


def _write_json(path: Path, data: dict | list) -> None:
    with open(path, "w") as f:
        json.dump(
            _sanitize_nonfinite(data),
            f,
            indent=2,
            cls=_NumpyEncoder,
            default=str,
            allow_nan=False,
        )
