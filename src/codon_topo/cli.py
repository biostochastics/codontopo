"""Command-line interface for codon-topo."""

import json
from pathlib import Path

import click
import numpy as np

from codon_topo import DEFAULT_SEED


class _NumpyEncoder(json.JSONEncoder):
    """JSON encoder that converts numpy types to native Python types."""

    def default(self, o):
        if isinstance(o, np.integer):
            return int(o)
        if isinstance(o, np.floating):
            return float(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        return super().default(o)


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
    type=click.Choice(["freeland_hurst", "class_size"]),
    help="Null model type.",
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
@click.option("--n", "n_samples", default=1_000, help="Monte Carlo samples per rho.")
@click.option("--seed", default=135325, help="Random seed.")
@click.option("--json", "as_json", is_flag=True, help="Output as JSON.")
def rho_sweep(n_samples: int, seed: int, as_json: bool) -> None:
    """Sweep transversion weight rho to test optimality robustness."""
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
@click.option("--n", "n_samples", default=1_000, help="Monte Carlo samples per table.")
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
        topology_avoidance_phylogenetic_sensitivity,
    )
    from codon_topo.analysis.cosmic_query import fano_predictions_for_kras
    from codon_topo.reports.claim_hierarchy import CLAIM_HIERARCHY
    from codon_topo.reports.catalogue import build_catalogue

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    n_steps = 14
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
    coloring_result = monte_carlo_null(n_samples=n_samples, seed=seed)
    cross_table = cross_table_optimality()

    _step("Multi-metric sensitivity (n=1000)...")
    metric_results = multi_metric_sensitivity(
        n_samples=min(n_samples, 10_000), seed=seed
    )

    _step("Rho robustness sweep...")
    rho_result = rho_robustness_sweep(n_samples=min(n_samples, 10_000), seed=seed)

    _step("Per-table optimality...")
    pertable_result = per_table_optimality(n_samples=min(n_samples, 1_000), seed=seed)

    _step("Per-table proximity audit (standard-code-proximity sensitivity)...")
    from codon_topo.analysis.coloring_optimality import per_table_proximity_audit

    proximity_audit = per_table_proximity_audit(
        n_samples=min(n_samples, 1_000), seed=seed
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

    # 10. Phylogenetic sensitivity
    _step("Phylogenetic sensitivity...")
    phylo_result = topology_avoidance_phylogenetic_sensitivity()
    _write_json(out / "phylogenetic_sensitivity.json", phylo_result)

    # 11. Evolutionary simulation (conditional logit)
    _step("Conditional logit models (M1-M4)...")
    from codon_topo.analysis.evolutionary_simulation import (
        run_clade_exclusion_sensitivity,
        run_evolutionary_simulation_analysis,
    )

    evosim_result = run_evolutionary_simulation_analysis(seed=seed)
    _write_json(out / "evolutionary_simulation.json", evosim_result)

    # 11b. Conditional-logit clade-exclusion sensitivity (Reviewer R1.C / R2.M1)
    _step("Conditional logit clade-exclusion sensitivity (7 regimes)...")
    condlogit_clade = run_clade_exclusion_sensitivity()
    _write_json(out / "condlogit_clade_sensitivity.json", condlogit_clade)

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
        "_version": "0.3.1",
        # Section 3.1: Cross-metric coloring optimality
        "coloring": {
            "observed_score": coloring_result["observed_score"],
            "null_mean": coloring_result["null_mean"],
            "null_std": coloring_result["null_std"],
            "quantile": coloring_result["quantile_of_observed"],
            "p_value": coloring_result["p_value_conservative"],
            "n_samples": n_samples,
        },
        "metrics": {
            name: {
                "observed": m.get("observed_score"),
                "null_mean": m.get("null_mean"),
                "null_std": m.get("null_std"),
                "z": m.get("effect_size_z"),
                "p": m.get("p_value_conservative"),
                "quantile": m.get("quantile"),
                "improvement_pct": m.get("improvement_pct"),
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
        },
        # Section 3.1: Score decomposition
        "decomposition": {
            "position_fractions": decomp_result["position_fractions"],
            "total_score": decomp_result["total_score"],
        },
        # Section 3.4: Topology avoidance (Q6)
        "topology_avoidance_q6": {
            "observed_breaks": topo_q6["observed_creates_disc"],
            "observed_total": topo_q6["observed_total"],
            "rate_observed": topo_q6["rate_observed"],
            "possible_breaks": topo_q6["possible_creates_disc"],
            "possible_total": topo_q6["possible_total"],
            "rate_possible": topo_q6["rate_possible"],
            "hypergeom_p": topo_q6["hypergeom_p"],
            "permutation_p": topo_q6["permutation_p"],
        },
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
        },
        # Section 3.6: tRNA enrichment
        "trna": {
            "n_pairings": fisher_result.get("n_pairings"),
            "stouffer_z": fisher_result.get("stouffer_z"),
            "stouffer_p": fisher_result.get("stouffer_p"),
            "mis_worst_p": mis_result.get("worst_case_stouffer_p"),
            "mis_best_p": mis_result.get("best_case_stouffer_p"),
            "n_mis": mis_result.get("n_mis_size_ge2"),
        },
        # Phylogenetic sensitivity
        "phylo": {
            "all_significant": phylo_result["all_clade_exclusions_significant"],
            "lineage_collapsed": phylo_result["lineage_collapsed"],
        },
        # Reassignment DB
        "reassignment": {
            "n_events": len(db),
            "n_dedup": topo_q6["observed_total"],
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
    click.echo(f"  MIS worst-case p: {mis_result.get('worst_case_stouffer_p', 'N/A')}")
    click.echo(
        f"  Phylo clade-exclusion: all significant = {phylo_result['all_clade_exclusions_significant']}"
    )
    click.echo(
        f"  Condlogit best model: {best_model} (AICc={best_fit.get('aicc', 'N/A')})"
    )
    click.echo("  manuscript_stats.json written for Typst dynamic rendering")


def _write_json(path: Path, data: dict | list) -> None:
    with open(path, "w") as f:
        json.dump(data, f, indent=2, cls=_NumpyEncoder, default=str)
