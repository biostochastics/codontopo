"""Statistical utilities for confidence intervals and effect sizes.

Provides Beta posterior CIs for permutation p-values, risk ratio CIs,
quantile CIs, and other inferential summaries that GPT-5.2-pro flagged
as missing from the manuscript's statistical reporting.
"""

import math
from scipy import stats


def beta_posterior_ci(
    k: int,
    n: int,
    confidence: float = 0.95,
    prior_alpha: float = 0.5,
    prior_beta: float = 0.5,
) -> dict:
    """Bayesian credible interval for the true tail probability p.

    Given k exceedances in n permutations, models p ~ Beta(k + a, n - k + b)
    where (a, b) is the prior. Default: Jeffreys prior Beta(0.5, 0.5).

    Returns point estimate (posterior mean), credible interval, and the
    conservative frequentist estimate (k+1)/(n+1) for comparison.
    """
    post_a = k + prior_alpha
    post_b = (n - k) + prior_beta
    alpha = 1 - confidence

    posterior_mean = post_a / (post_a + post_b)
    ci_lo = float(stats.beta.ppf(alpha / 2, post_a, post_b))
    ci_hi = float(stats.beta.ppf(1 - alpha / 2, post_a, post_b))
    conservative_p = (k + 1) / (n + 1)

    return {
        "posterior_mean": posterior_mean,
        "ci_low": ci_lo,
        "ci_high": ci_hi,
        "confidence": confidence,
        "conservative_p": conservative_p,
        "k": k,
        "n": n,
    }


def risk_ratio_ci(a: int, n1: int, b: int, n2: int, confidence: float = 0.95) -> dict:
    """Risk ratio with log-normal CI.

    RR = (a/n1) / (b/n2) with 95% CI via log-normal approximation.
    a/n1 = observed rate, b/n2 = expected/possible rate.
    """
    p1 = a / max(n1, 1)
    p2 = b / max(n2, 1)

    if p2 == 0:
        return {"rr": 0.0, "ci_low": 0.0, "ci_high": 0.0, "p1": p1, "p2": p2}

    rr = p1 / p2
    z = stats.norm.ppf(1 - (1 - confidence) / 2)

    if a > 0 and n1 > a:
        se_log_rr = math.sqrt((1 / max(a, 1) - 1 / n1) + (1 / max(b, 1) - 1 / n2))
        ci_lo = math.exp(math.log(max(rr, 1e-20)) - z * se_log_rr)
        ci_hi = math.exp(math.log(max(rr, 1e-20)) + z * se_log_rr)
    else:
        se_log_rr = float("inf")
        ci_lo = 0.0
        ci_hi = float("inf")

    return {
        "rr": rr,
        "ci_low": ci_lo,
        "ci_high": ci_hi,
        "se_log_rr": se_log_rr,
        "p1": p1,
        "p2": p2,
        "confidence": confidence,
    }


def quantile_ci(k: int, n: int, confidence: float = 0.95) -> dict:
    """Confidence interval for a quantile estimate q = k/n.

    Uses the Clopper-Pearson (exact binomial) interval.
    """
    q_hat = k / max(n, 1)
    alpha = 1 - confidence

    if k == 0:
        ci_lo = 0.0
    else:
        ci_lo = float(stats.beta.ppf(alpha / 2, k, n - k + 1))

    if k == n:
        ci_hi = 1.0
    else:
        ci_hi = float(stats.beta.ppf(1 - alpha / 2, k + 1, n - k))

    return {
        "quantile": q_hat,
        "ci_low": ci_lo,
        "ci_high": ci_hi,
        "confidence": confidence,
        "k": k,
        "n": n,
    }


def percent_improvement(observed: float, null_mean: float) -> float:
    """Percent improvement of observed over null mean."""
    if null_mean == 0:
        return 0.0
    return 100.0 * (null_mean - observed) / null_mean
