"""Statistical models for CodonSafe meta-analysis.

Primary analysis: logistic regression of recoding success/failure on
topology classification, adjusting for known confounders.

Secondary analysis: linear regression of continuous fitness on delta_F.

Sensitivity: encoding sweep (all 24 base-to-bit encodings),
leave-one-study-out cross-validation.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    import pandas as pd


@dataclass(frozen=True)
class ModelSpec:
    """Specification for a regression model.

    Parameters
    ----------
    formula : str
        Patsy/statsmodels formula string. Example:
        "success ~ changes_components_eps1 + C(study) + cov_mRNA_structure_deviation"
    family : str
        "binomial" for logistic regression, "gaussian" for linear.
    cluster_col : str or None
        Column for cluster-robust standard errors (e.g., "unit_id").
    study_fixed_effects : bool
        Whether to include C(study) in the formula.
    """

    formula: str
    family: str = "binomial"
    cluster_col: str | None = None
    study_fixed_effects: bool = True


@dataclass
class ModelResult:
    """Results from a fitted regression model.

    Parameters
    ----------
    spec : ModelSpec
        The model specification used.
    coefficients : dict
        Parameter name -> coefficient estimate.
    std_errors : dict
        Parameter name -> standard error.
    p_values : dict
        Parameter name -> p-value.
    odds_ratios : dict or None
        Parameter name -> exp(coefficient). Only for logistic models.
    ci_lower : dict
        Parameter name -> lower 95% CI bound.
    ci_upper : dict
        Parameter name -> upper 95% CI bound.
    n_obs : int
        Number of observations in the model.
    aic : float
        Akaike Information Criterion.
    bic : float
        Bayesian Information Criterion.
    pseudo_r2 : float or None
        McFadden pseudo-R^2 (logistic) or R^2 (linear).
    convergence : bool
        Whether the optimizer converged.
    summary_text : str
        Full text summary from statsmodels.
    """

    spec: ModelSpec
    coefficients: dict[str, float] = field(default_factory=dict)
    std_errors: dict[str, float] = field(default_factory=dict)
    p_values: dict[str, float] = field(default_factory=dict)
    odds_ratios: dict[str, float] | None = None
    ci_lower: dict[str, float] = field(default_factory=dict)
    ci_upper: dict[str, float] = field(default_factory=dict)
    n_obs: int = 0
    aic: float = 0.0
    bic: float = 0.0
    pseudo_r2: float | None = None
    convergence: bool = True
    summary_text: str = ""


def fit_model(df: "pd.DataFrame", spec: ModelSpec) -> ModelResult:
    """Fit a generalized linear model using statsmodels.

    Parameters
    ----------
    df : DataFrame
        Analysis table from annotated_swaps_to_frame().
    spec : ModelSpec
        Model specification.

    Returns
    -------
    ModelResult
        Fitted model results.
    """
    import numpy as np
    import statsmodels.api as sm
    import statsmodels.formula.api as smf

    if spec.family == "binomial":
        family = sm.families.Binomial()
    elif spec.family == "gaussian":
        family = sm.families.Gaussian()
    else:
        raise ValueError(f"Unknown family: {spec.family!r}")

    # Drop rows with missing outcome
    outcome_col = spec.formula.split("~")[0].strip()
    df_clean = df.dropna(subset=[outcome_col])

    model = smf.glm(spec.formula, data=df_clean, family=family)

    if spec.cluster_col and spec.cluster_col in df_clean.columns:
        fit = model.fit(
            cov_type="cluster",
            cov_kwds={"groups": df_clean[spec.cluster_col]},
        )
    else:
        fit = model.fit()

    params = fit.params
    ci = fit.conf_int()

    coefficients = {k: float(v) for k, v in params.items()}
    std_errors = {k: float(v) for k, v in fit.bse.items()}
    p_values = {k: float(v) for k, v in fit.pvalues.items()}
    ci_lower = {k: float(ci.loc[k, 0]) for k in params.index}
    ci_upper = {k: float(ci.loc[k, 1]) for k in params.index}

    odds_ratios = None
    if spec.family == "binomial":
        odds_ratios = {k: float(np.exp(v)) for k, v in params.items()}

    pseudo_r2 = None
    if spec.family == "binomial":
        null_model = smf.glm(f"{outcome_col} ~ 1", data=df_clean, family=family).fit()
        pseudo_r2 = float(1 - fit.llf / null_model.llf) if null_model.llf != 0 else None

    return ModelResult(
        spec=spec,
        coefficients=coefficients,
        std_errors=std_errors,
        p_values=p_values,
        odds_ratios=odds_ratios,
        ci_lower=ci_lower,
        ci_upper=ci_upper,
        n_obs=int(fit.nobs),
        aic=float(fit.aic),
        bic=float(fit.bic),
        pseudo_r2=pseudo_r2,
        convergence=fit.converged,
        summary_text=str(fit.summary()),
    )


def leave_one_study_out(
    df: "pd.DataFrame",
    spec: ModelSpec,
    study_col: str = "study",
) -> dict[str, Any]:
    """Leave-one-study-out cross-validation.

    For each study, fit the model on all other studies and predict
    on the held-out study. Reports per-study AUC (for binary outcomes)
    or RMSE (for continuous outcomes).

    Parameters
    ----------
    df : DataFrame
        Full analysis table.
    spec : ModelSpec
        Model specification.
    study_col : str
        Column identifying the study.

    Returns
    -------
    dict with keys:
      per_study: list of {study, n_test, auc_or_rmse, metric}
      mean_metric: float
      overall_metric_name: str
    """
    import numpy as np
    from sklearn.metrics import roc_auc_score

    studies = df[study_col].unique()
    outcome_col = spec.formula.split("~")[0].strip()
    is_binary = spec.family == "binomial"

    per_study = []
    for held_out in studies:
        train = df[df[study_col] != held_out].copy()
        test = df[df[study_col] == held_out].copy()

        if len(test) < 5 or len(train) < 10:
            continue

        # Remove study fixed effects for LOSO (can't have held-out level)
        rhs = spec.formula.split("~", 1)[1].strip()
        rhs_stripped = rhs.replace(f"C({study_col})", "").strip(" +")
        # Clean up dangling operators from removal
        import re

        rhs_stripped = re.sub(r"\+\s*\+", "+", rhs_stripped).strip(" +")
        if not rhs_stripped:
            continue  # study-only model, skip this fold

        loso_formula = f"{outcome_col} ~ {rhs_stripped}"
        loso_spec = ModelSpec(
            formula=loso_formula,
            family=spec.family,
            cluster_col=spec.cluster_col,
            study_fixed_effects=False,
        )

        try:
            fit_model(train, loso_spec)
        except Exception:
            continue

        # Predict on test set using the fitted model
        import statsmodels.api as sm
        import statsmodels.formula.api as smf

        if spec.family == "binomial":
            family = sm.families.Binomial()
        else:
            family = sm.families.Gaussian()

        model = smf.glm(loso_spec.formula, data=train, family=family).fit()
        try:
            test_clean = test.dropna(subset=[outcome_col])
            preds = model.predict(test_clean)
        except Exception:
            continue

        if len(test_clean) != len(preds):
            continue

        if is_binary:
            y_true = test_clean[outcome_col].astype(float)
            if y_true.nunique() < 2:
                continue
            metric_val = float(roc_auc_score(y_true, preds))
            metric_name = "auc"
        else:
            y_true = test_clean[outcome_col].astype(float)
            metric_val = float(np.sqrt(np.mean((y_true - preds) ** 2)))
            metric_name = "rmse"

        per_study.append(
            {
                "study": held_out,
                "n_test": len(test_clean),
                "metric_value": metric_val,
                "metric_name": metric_name,
            }
        )

    mean_val = (
        float(np.mean([r["metric_value"] for r in per_study]))
        if per_study
        else float("nan")
    )

    return {
        "per_study": per_study,
        "mean_metric": mean_val,
        "overall_metric_name": per_study[0]["metric_name"] if per_study else "unknown",
        "n_studies_evaluated": len(per_study),
    }


def encoding_sensitivity(
    df_builder: Any,
    *,
    spec: ModelSpec,
    table_id: int = 1,
) -> dict[str, Any]:
    """Run the primary model under all 24 encodings.

    Parameters
    ----------
    df_builder : callable
        Function that takes encoding_id (int) and returns a DataFrame
        ready for modeling. This is necessary because the topology
        features depend on the encoding.
    spec : ModelSpec
        Model specification.
    table_id : int
        NCBI translation table.

    Returns
    -------
    dict with keys:
      per_encoding: list of {encoding_id, coefficients, p_values, aic}
      topology_feature_robust: bool (True if significant at p<0.05 across all 24)
      topology_p_range: (min_p, max_p)
    """
    import numpy as np

    per_encoding = []
    topo_p_values = []

    for enc_id in range(24):
        df = df_builder(enc_id)
        try:
            result = fit_model(df, spec)
        except Exception:
            continue

        # Look for the topology coefficient
        topo_key = None
        for k in result.p_values:
            if "changes_components_eps1" in k or "topology_breaking" in k:
                topo_key = k
                break

        entry: dict[str, Any] = {
            "encoding_id": enc_id,
            "coefficients": result.coefficients,
            "p_values": result.p_values,
            "aic": result.aic,
            "n_obs": result.n_obs,
        }
        if topo_key:
            entry["topology_p"] = result.p_values[topo_key]
            entry["topology_coef"] = result.coefficients[topo_key]
            topo_p_values.append(result.p_values[topo_key])

        per_encoding.append(entry)

    robust = all(p < 0.05 for p in topo_p_values) if topo_p_values else False
    p_range = (
        (float(np.min(topo_p_values)), float(np.max(topo_p_values)))
        if topo_p_values
        else (float("nan"), float("nan"))
    )

    return {
        "per_encoding": per_encoding,
        "topology_feature_robust": robust,
        "topology_p_range": p_range,
        "n_encodings_fitted": len(per_encoding),
        "n_encodings_significant_p05": sum(1 for p in topo_p_values if p < 0.05),
    }
