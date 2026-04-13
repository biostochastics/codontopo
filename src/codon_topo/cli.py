"""Command-line interface for codon-topo."""

import json
from pathlib import Path

import click

from codon_topo import DEFAULT_SEED


def _json_out(data: dict | list) -> None:
    click.echo(json.dumps(data, indent=2, default=str))


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
@click.option("--all-tables", is_flag=True, help="Run across all 25 NCBI tables.")
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
@click.option("--all-tables", is_flag=True, help="Run across all 25 NCBI tables.")
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


@main.command("all")
@click.option("--output-dir", default="./output", help="Directory for CSV/JSON output.")
@click.option("--seed", default=DEFAULT_SEED, help="Random seed for Monte Carlo.")
@click.option("--n", "n_samples", default=10_000, help="Monte Carlo sample size.")
def run_all(output_dir: str, seed: int, n_samples: int) -> None:
    """Run all analyses and write reports to output directory."""
    from codon_topo.core.filtration import analyze_filtration
    from codon_topo.core.homology import disconnection_catalogue
    from codon_topo.core.genetic_codes import all_table_ids, get_code
    from codon_topo.analysis.coloring_optimality import (
        monte_carlo_null,
        cross_table_optimality,
    )
    from codon_topo.analysis.reassignment_db import (
        build_reassignment_db,
        hamming_path_lengths,
        bit_position_bias_weighted as _bit_bias_weighted,
        directionality_summary,
    )
    from codon_topo.analysis.trna_evidence import trna_duplication_correlation_test
    from codon_topo.analysis.depth_calibration import (
        compute_correlation,
        depth_calibration_table,
    )
    from codon_topo.analysis.cosmic_query import fano_predictions_for_kras
    from codon_topo.reports.claim_hierarchy import CLAIM_HIERARCHY
    from codon_topo.reports.catalogue import build_catalogue

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    click.echo("Running all analyses...")

    # Filtration across all tables
    click.echo("  [1/8] Filtration...")
    filt_results = []
    for tid in all_table_ids():
        code = get_code(tid)
        r = analyze_filtration(code)
        r["table_id"] = tid
        filt_results.append(r)
    _write_json(out / "filtration.json", filt_results)

    # Disconnection catalogue
    click.echo("  [2/8] Disconnections...")
    disc_results = []
    for tid in all_table_ids():
        code = get_code(tid)
        cat = disconnection_catalogue(code)
        for entry in cat:
            entry["table_id"] = tid
        disc_results.extend(cat)
    _write_json(out / "disconnections.json", disc_results)

    # Coloring optimality
    click.echo(f"  [3/8] Coloring optimality (n={n_samples})...")
    coloring_result = monte_carlo_null(n_samples=n_samples, seed=seed)
    cross_table = cross_table_optimality()
    _write_json(
        out / "coloring_optimality.json",
        {
            "monte_carlo": coloring_result,
            "cross_table": cross_table,
        },
    )

    # Reassignment DB + bit bias
    click.echo("  [4/8] Reassignment analysis...")
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

    # tRNA evidence
    click.echo("  [5/8] tRNA evidence...")
    trna_result = trna_duplication_correlation_test()
    _write_json(out / "trna_evidence.json", trna_result)

    # Depth calibration
    click.echo("  [6/8] Depth calibration...")
    corr = compute_correlation()
    cal_table = depth_calibration_table()
    _write_json(
        out / "depth_calibration.json",
        {
            "correlation": corr,
            "calibration_points": cal_table,
        },
    )

    # KRAS-Fano predictions (offline — no API call in batch)
    click.echo("  [7/8] KRAS-Fano predictions (offline)...")
    fano_preds = fano_predictions_for_kras()
    _write_json(
        out / "kras_fano.json",
        {
            "mode": "offline",
            "predictions": fano_preds,
        },
    )

    # Claim hierarchy + catalogue
    click.echo("  [8/8] Claims & catalogue...")
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
    _write_json(out / "claims.json", claims_data)
    _write_json(out / "catalogue.json", catalogue_data)

    # Summary
    click.echo(f"\nDone. {len(list(out.glob('*.json')))} JSON files written to {out}/")
    click.echo(
        f"  Coloring: quantile={coloring_result['quantile_of_observed']:.2f}%, "
        f"p={coloring_result['p_value_conservative']:.6f}"
    )
    click.echo(
        f"  tRNA:     {trna_result['n_with_elevated_trna']}/4 elevated, "
        f"p={trna_result['binomial_p_value']:.4f}"
    )
    click.echo(
        f"  Depth:    rho={corr['spearman_rho']:.3f}, p={corr['spearman_p']:.3f}"
    )


def _write_json(path: Path, data: dict | list) -> None:
    with open(path, "w") as f:
        json.dump(data, f, indent=2, default=str)
