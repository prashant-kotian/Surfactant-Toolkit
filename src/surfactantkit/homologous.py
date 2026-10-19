"""Klevens homologous-series CMC fit (Phase 2 worker C2)."""

from __future__ import annotations

import math
from dataclasses import dataclass

_R_GAS = 8.314462618


@dataclass
class HomologousSeriesResult:
    intercept_A: float          # log10(CMC) = A - B * n
    reduction_per_ch2: float    # B (positive: each CH2 lowers log10 CMC by B)
    slope_log10_per_ch2: float  # == -reduction_per_ch2 (the fitted slope)
    reduction_factor_per_ch2: float  # 10^-slope: multiplicative CMC step per CH2
    intercept_se: float
    slope_se: float
    r_squared: float
    rmse_log10: float
    n_points: int
    predicted: list[dict]       # [{"chain_length", "log10_cmc", "cmc", "ci_low", "ci_high", "extrapolation"}]
    unit: str
    converged: bool
    method: str
    params: dict
    caveat: str


def cmc_homologous_series(
    chain_lengths: list[int],
    cmc_values: list[float],
    unit: str = "mM",
    predict_at: list[int] | None = None,
    ci: float = 0.95,
) -> HomologousSeriesResult:
    """Klevens/log-linear homologous-series fit: log10(CMC) = A - B*n.

    Fits A and B (least squares) from a measured homologous series, reports
    slope SEs, R^2 and RMSE, and optionally extrapolates to unseen chain
    lengths with prediction intervals.  This is the standard empirical rule
    behind "each added CH2 cuts the CMC ~3x" -- a correlation, NOT a
    thermodynamic law: predictions far outside the fitted chain-length range
    are flagged as extrapolation per point.

    cmc_values: CMCs in `unit` ("mM", "uM", "M", "mol/L" -- label only; the
    log is taken of the number you pass, so keep one consistent unit).
    predict_at: chain lengths to predict (inside or outside the data range).
    """
    if not (0.5 < ci < 1.0):
        raise ValueError(f"ci must be in (0.5, 1.0), got {ci}")
    if unit not in ("mM", "uM", "M", "mol/L", "mol/m3"):
        raise ValueError(f"unit must be a CMC concentration label, got {unit!r}")
    n = len(chain_lengths)
    if n != len(cmc_values):
        raise ValueError("chain_lengths and cmc_values must be the same length")
    if n < 3:
        raise ValueError(
            "need at least 3 homologous points to fit a slope with any certainty"
        )
    if n != len(set(chain_lengths)):
        raise ValueError("chain_lengths must be distinct (duplicate n adds no leverage)")
    if any(c <= 0 for c in cmc_values):
        raise ValueError("cmc_values must be > 0 (log10 is undefined at 0)")
    if min(chain_lengths) == max(chain_lengths):
        raise ValueError("chain_lengths must span more than one value")

    ys = [math.log10(c) for c in cmc_values]
    xs = [float(v) for v in chain_lengths]
    xbar = sum(xs) / n
    ybar = sum(ys) / n
    sxx = sum((v - xbar) ** 2 for v in xs)
    if sxx <= 0:
        raise ValueError("chain_lengths carry no spread")
    slope = sum((xv - xbar) * (yv - ybar) for xv, yv in zip(xs, ys)) / sxx
    intercept = ybar - slope * xbar
    resid = [yv - (intercept + slope * xv) for xv, yv in zip(xs, ys)]
    dof = n - 2
    rss = sum(r * r for r in resid)
    s2 = rss / dof
    s = math.sqrt(s2)
    slope_se = math.sqrt(s2 / sxx)
    intercept_se = math.sqrt(s2 * (1.0 / n + xbar * xbar / sxx))
    ss_tot = sum((yv - ybar) ** 2 for yv in ys)
    r2 = 1.0 - rss / ss_tot if ss_tot > 0 else 1.0
    rmse = math.sqrt(rss / n)

    # t critical: scipy when available (calc extra), else normal approximation
    t_crit, t_source = _t_critical(ci, dof)

    predicted: list[dict] = []
    if predict_at:
        n_lo, n_hi = min(xs), max(xs)
        for np_ in predict_at:
            if np_ is None:
                continue
            xn = float(np_)
            y0 = intercept + slope * xn
            se_pred = s * math.sqrt(1.0 + 1.0 / n + (xn - xbar) ** 2 / sxx)
            predicted.append({
                "chain_length": int(np_),
                "log10_cmc": y0,
                "cmc": 10.0 ** y0,
                "ci_low": 10.0 ** (y0 - t_crit * se_pred),
                "ci_high": 10.0 ** (y0 + t_crit * se_pred),
                "extrapolation": bool(xn < n_lo or xn > n_hi),
            })

    caveat = (
        "Klevens is an empirical correlation over a narrow chain range -- "
        "prediction intervals include the fit's own residual scatter (RMSE "
        f"{rmse:.3f} log10) and extrapolated points are flagged; slopes far "
        "from the -0.28..-0.32 CH2 range known for ionic surfactant series "
        "usually mean a mixed or impure series"
    )
    return HomologousSeriesResult(
        intercept_A=float(intercept),
        reduction_per_ch2=float(-slope),
        slope_log10_per_ch2=float(slope),
        reduction_factor_per_ch2=float(10.0 ** (-slope)),
        intercept_se=float(intercept_se),
        slope_se=float(slope_se),
        r_squared=float(r2),
        rmse_log10=float(rmse),
        n_points=n,
        predicted=predicted,
        unit=unit,
        converged=True,
        method="klevens log-linear OLS",
        params={
            "ci": ci,
            "t_critical": t_crit,
            "t_source": t_source,
            "dof": dof,
            "chain_lengths": [int(v) for v in chain_lengths],
            "cmc_values": [float(v) for v in cmc_values],
        },
        caveat=caveat,
    )


def _t_critical(ci: float, dof: int) -> tuple[float, str]:
    """(t_{ci/2, dof}, source).  Falls back to the normal quantile with a
    disclosed source string when scipy is unavailable."""
    try:
        from scipy import stats
        return float(stats.t.ppf(0.5 + ci / 2.0, dof)), "student-t (scipy)"
    except ImportError:
        # normal quantile via a rational approximation (Acklam); disclosed
        p = 1.0 - (1.0 - ci) / 2.0
        return _normal_ppf(p), "normal-approximation (scipy not installed)"


def _normal_ppf(p: float) -> float:
    """Inverse standard-normal CDF (Acklam's rational approximation,
    |error| < 1.15e-9 -- adequate for interval display without scipy)."""
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00]
    plow, phigh = 0.02425, 1 - 0.02425
    if p < plow:
        q = math.sqrt(-2 * math.log(p))
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / \
            ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1)
    if p > phigh:
        q = math.sqrt(-2 * math.log(1 - p))
        return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / \
            ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1)
    q = p - 0.5
    r = q * q
    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / \
        (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1)
