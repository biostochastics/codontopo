"""Aggregation and feature engineering for CodonSafe meta-analysis.

Converts lists of AnnotatedSwap records into pandas DataFrames suitable
for statistical modeling. Handles:
  - Flattening nested dataclasses into tabular form
  - Study-level aggregation for segment-level outcomes
  - Feature engineering for logistic/linear regression
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from codon_topo.analysis.codonsafe.models import AnnotatedSwap

if TYPE_CHECKING:
    import pandas as pd


def annotated_swaps_to_frame(records: list[AnnotatedSwap]) -> "pd.DataFrame":
    """Convert a list of AnnotatedSwap records to a pandas DataFrame.

    Each row corresponds to one (event, encoding) pair. Columns include:

    Identifiers:
      study, event_id, unit_id, encoding_id

    Outcome:
      outcome_type, success, fitness, fitness_unit

    Topology features:
      hamming, source_aa, target_aa, is_synonymous,
      changes_components_eps1, crosses_component_boundary_eps1,
      delta_F_grantham, delta_F_miyata, delta_F_polar_requirement,
      delta_F_kyte_doolittle,
      local_mismatch_source_grantham, local_mismatch_target_grantham

    Covariates:
      Expanded from event.covariates dict.
    """
    import pandas as pd

    rows = []
    for rec in records:
        row: dict = {
            # Identifiers
            "study": rec.event.study.value,
            "event_id": rec.event.event_id,
            "unit_id": rec.event.unit_id,
            "gene": rec.event.gene,
            "is_essential_gene": rec.event.is_essential_gene,
            "table_id": rec.event.table_id,
            "source_codon": rec.event.source_codon,
            "target_codon": rec.event.target_codon,
            # Encoding
            "encoding_id": rec.topo.encoding_id,
            # Topology
            "hamming": rec.topo.hamming,
            "source_aa": rec.topo.source_aa,
            "target_aa": rec.topo.target_aa,
            "is_synonymous": rec.topo.is_synonymous,
            "changes_components_eps1": rec.topo.changes_components_eps1,
            "crosses_component_boundary_eps1": rec.topo.crosses_component_boundary_eps1,
        }

        # Delta F per metric
        for metric, value in rec.topo.delta_edge_mismatch.items():
            row[f"delta_F_{metric}"] = value

        # Local mismatch per metric
        for metric, value in rec.topo.local_mismatch_source.items():
            row[f"local_mismatch_source_{metric}"] = value
        for metric, value in rec.topo.local_mismatch_target.items():
            row[f"local_mismatch_target_{metric}"] = value

        # Outcome
        if rec.outcome is not None:
            row["outcome_type"] = rec.outcome.outcome_type.value
            row["success"] = rec.outcome.success
            row["fitness"] = rec.outcome.fitness
            row["fitness_unit"] = rec.outcome.fitness_unit
        else:
            row["outcome_type"] = None
            row["success"] = None
            row["fitness"] = None
            row["fitness_unit"] = None

        # Covariates
        for k, v in rec.event.covariates.items():
            row[f"cov_{k}"] = v

        rows.append(row)

    return pd.DataFrame(rows)


def aggregate_to_units(
    df: "pd.DataFrame",
    *,
    unit_id_col: str = "unit_id",
    study_col: str = "study",
) -> "pd.DataFrame":
    """Aggregate per-event data to per-unit (segment) level.

    For studies like Fredens/Ostrov where the outcome is per-segment
    rather than per-codon, this aggregates topology features within
    each unit.

    Aggregated features:
      - n_events: number of swaps in the unit
      - frac_topology_breaking: fraction with changes_components_eps1=True
      - frac_ser_boundary: fraction with crosses_component_boundary_eps1=True
      - mean_delta_F_*: mean delta_F per metric
      - max_delta_F_*: max delta_F per metric
      - n_ser_crossings: count of Ser boundary crossings
    """

    agg_dict: dict[str, tuple[str, str]] = {
        "n_events": ("event_id", "count"),
    }

    # Topology breaking fraction
    if "changes_components_eps1" in df.columns:
        agg_dict["frac_topology_breaking"] = ("changes_components_eps1", "mean")

    # Ser boundary
    if "crosses_component_boundary_eps1" in df.columns:
        # Count non-None True values
        pass  # handled separately below

    # Delta F metrics
    delta_f_cols = [c for c in df.columns if c.startswith("delta_F_")]
    for col in delta_f_cols:
        metric = col.replace("delta_F_", "")
        agg_dict[f"mean_delta_F_{metric}"] = (col, "mean")
        agg_dict[f"max_delta_F_{metric}"] = (col, "max")

    grouped = df.groupby([study_col, unit_id_col])
    # groupby(...).agg(**kwargs) returns a DataFrame; the union type is overly broad
    result: pd.DataFrame = grouped.agg(**agg_dict)  # type: ignore[assignment]

    # Ser boundary crossings (handle None values)
    if "crosses_component_boundary_eps1" in df.columns:
        ser_counts = grouped["crosses_component_boundary_eps1"].apply(
            lambda s: s.dropna().sum()
        )
        ser_total = grouped["crosses_component_boundary_eps1"].apply(
            lambda s: s.dropna().count()
        )
        result["n_ser_crossings"] = ser_counts
        result["frac_ser_boundary"] = ser_counts / ser_total.replace(0, float("nan"))

    # Outcome (take first — should be uniform within a unit for segment-level)
    if "success" in df.columns:
        result["success"] = grouped["success"].first()
    if "fitness" in df.columns:
        result["fitness"] = grouped["fitness"].first()

    return result.reset_index()
