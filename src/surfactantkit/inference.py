"""Uncertainty-quantification workers (Phase 2: I1, I3, I4, I5).

The rest of the toolkit returns point estimates; these workers attach the
uncertainty to them: Levenberg-Marquardt fits with parameter standard
errors (I3), Metropolis posterior sampling with real chain diagnostics
(I1), shared-parameter global fits across several datasets (I4), and
Monte Carlo error propagation through toolkit function chains (I5).

Everything here lazy-imports numpy/scipy and names the ``calc`` extra in
its ImportError, per the project's optional-dependency policy.  The core
library stays importable without them.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

_R_GAS = 8.314462618
_REQUIRED_CALC_HINT = (
    "inference workers require numpy+scipy -- install with "
    "'pip install surfactantkit[calc]'"
)


def _require_calc():
    try:
        import numpy as np
        import scipy  # noqa: F401
    except ImportError as exc:  # pragma: no cover - exercised only without extra
        raise ImportError(_REQUIRED_CALC_HINT) from exc
    return np


# ---------------------------------------------------------------------------
# Model registry (I3): named models so MCP/JSON callers can fit without
# shipping a callable.
# ---------------------------------------------------------------------------

def _model_linear(x, m, b):
    import numpy as np
    return m * np.asarray(x, dtype=float) + b


def _model_exponential_decay(x, amplitude, rate, offset):
    import numpy as np
    return amplitude * np.exp(-rate * np.asarray(x, dtype=float)) + offset


def _model_szyszkowski(x, K, *, gamma0_mN_m, gamma_max_mol_per_m2,
                       system_type, temperature_K=298.15, r_gas=_R_GAS):
    from .adsorption import szyszkowski_surface_tension
    return [szyszkowski_surface_tension(float(c), gamma0_mN_m, gamma_max_mol_per_m2,
                                        float(K), system_type, temperature_K, r_gas)
            for c in x]


def _model_frumkin(x, K, a, *, gamma0_mN_m, gamma_max_mol_per_m2,
                   system_type, temperature_K=298.15, r_gas=_R_GAS):
    from .adsorption import frumkin_surface_tension
    return [frumkin_surface_tension(float(c), gamma0_mN_m, gamma_max_mol_per_m2,
                                    float(K), float(a), system_type, temperature_K, r_gas)
            for c in x]


#: registry name -> (callable, default initial guess, bounds, fixed-kwarg names)
MODEL_REGISTRY = {
    "linear": (_model_linear, [1.0, 0.0], None, ()),
    "exponential_decay": (_model_exponential_decay, [1.0, 1.0, 0.0],
                          ([-math.inf, 0.0, -math.inf], [math.inf, math.inf, math.inf]), ()),
    "szyszkowski": (_model_szyszkowski, [1.0], ([0.0], [math.inf]),
                    ("gamma0_mN_m", "gamma_max_mol_per_m2", "system_type",
                     "temperature_K", "r_gas")),
    "frumkin": (_model_frumkin, [1.0, 0.0], ([0.0, -20.0], [math.inf, 20.0]),
                ("gamma0_mN_m", "gamma_max_mol_per_m2", "system_type",
                 "temperature_K", "r_gas")),
}


@dataclass
class NonlinearFitResult:
    """Least-squares fit with covariance-based parameter uncertainties."""

    param_names: list[str]
    param_values: list[float]
    param_sd: list[float]
    ci_low: list[float]
    ci_high: list[float]
    covariance: list[list[float]]
    r_squared: float
    adj_r_squared: float
    rmse: float
    max_abs_standardized_residual: float
    residual_lag1_autocorr: float
    n_points: int
    n_params: int
    dof: int
    converged: bool
    n_iterations: int
    method: str
    params: dict
    caveat: str


def nonlinear_fit_with_errors(
    x: list[float],
    y: list[float],
    model: str,
    p0: list[float] | None = None,
    fixed: dict | None = None,
    sigma: list[float] | None = None,
    bounds: tuple[list[float], list[float]] | None = None,
    param_names: list[str] | None = None,
    max_nfev: int = 4000,
) -> NonlinearFitResult:
    """Levenberg-Marquardt-class nonlinear least squares with parameter
    standard errors, 95% CIs, covariance, R^2 and residual diagnostics.

    This is the error-reporting replacement for the toolkit's golden-section
    K-fits: same model, but you also get sigma_j, the full covariance matrix,
    adjusted R^2, and standardized-residual / lag-1-autocorrelation checks
    (a drifting residual sequence means the model is wrong, not noisy).

    model: registry key -- 'linear', 'exponential_decay', 'szyszkowski',
    'frumkin' (see MODEL_REGISTRY).  Callers pass data, not code, so the
    same call works through MCP.
    fixed: fixed kwargs for models that need them (Szyszkowski/Frumkin need
    gamma0_mN_m, gamma_max_mol_per_m2, system_type; K and a are fit).
    sigma: optional per-point observation uncertainties (weights = 1/sigma).
    bounds: optional (lower, upper) parameter bounds; registry defaults apply
    when omitted (e.g. K >= 0).

    SEs come from (J^T W J)^-1 * s^2 with s^2 = RSS/(n-k) -- the standard
    nonlinear-regression covariance, which assumes the model is locally
    correct; it does not cover model misspecification.
    """
    np = _require_calc()
    from scipy import stats
    from scipy.optimize import least_squares

    if model not in MODEL_REGISTRY:
        raise ValueError(
            f"unknown model {model!r} -- choose one of {sorted(MODEL_REGISTRY)}"
        )
    func, default_p0, reg_bounds, fixed_names = MODEL_REGISTRY[model]
    if fixed is None:
        fixed = {}
    missing = [k for k in fixed_names if k in fixed_names and k not in fixed and k != "r_gas"
               and k != "temperature_K"]
    if missing and model in ("szyszkowski", "frumkin"):
        raise ValueError(
            f"model {model!r} requires fixed= kwargs {missing} "
            "(gamma0_mN_m, gamma_max_mol_per_m2, system_type are measured, not fitted)"
        )

    xs = np.asarray(x, dtype=float)
    ys = np.asarray(y, dtype=float)
    if xs.shape != ys.shape or xs.ndim != 1:
        raise ValueError("x and y must be same-length 1-D sequences")
    n = xs.size
    q = default_p0 if p0 is None else [float(v) for v in p0]
    k = len(q)
    if n <= k:
        raise ValueError(f"need more points than parameters: n={n}, k={k}")
    if sigma is not None:
        sig = np.asarray(sigma, dtype=float)
        if sig.shape != ys.shape:
            raise ValueError("sigma must match y in length")
        if np.any(sig <= 0):
            raise ValueError("sigma values must be positive")
        weights = 1.0 / sig
    else:
        weights = np.ones_like(ys)

    if bounds is None:
        bounds = reg_bounds
    if bounds is None:
        bounds = (-np.inf, np.inf)

    def _apply(p):
        try:
            out = np.asarray(func(xs, *p, **fixed), dtype=float)
        except TypeError:
            out = np.asarray(func(xs, *p), dtype=float)
        if out.shape != ys.shape:
            raise ValueError(
                f"model {model!r} returned shape {out.shape}, expected {ys.shape}"
            )
        return out

    def _resid(p):
        return (_apply(p) - ys) * weights

    # trf handles finite AND infinite bounds (lm does not accept bounds at all)
    sol = least_squares(_resid, q, bounds=bounds, max_nfev=max_nfev, method="trf")
    popt = sol.x
    resid = _resid(popt)
    rss = float(np.sum(resid ** 2))
    dof = n - k
    if dof <= 0:
        raise ValueError("no degrees of freedom left for error estimation")

    # covariance = pinv(J^T J) * s^2  (weights already inside jac)
    jac = sol.jac
    jtj = jac.T @ jac
    try:
        cov = np.linalg.inv(jtj) * (rss / dof)
    except np.linalg.LinAlgError:
        cov = np.linalg.pinv(jtj) * (rss / dof)
    sd = np.sqrt(np.maximum(np.diag(cov), 0.0))
    t_crit = float(stats.t.ppf(0.975, dof)) if k < n else 1.96

    _apply(popt)
    ss_tot = float(np.sum((ys - ys.mean()) ** 2))
    r2 = 1.0 - rss / ss_tot if ss_tot > 0 else 0.0
    adj = 1.0 - (1.0 - r2) * (n - 1) / dof if dof > 0 else 0.0
    rmse = math.sqrt(rss / n)
    s_res = math.sqrt(rss / dof)
    std_resid = resid / s_res if s_res > 0 else resid
    if n > 2 and float(np.sum(resid ** 2)) > 0:
        lag1 = float(np.corrcoef(resid[:-1], resid[1:])[0, 1])
        if not math.isfinite(lag1):
            lag1 = 0.0
    else:
        lag1 = 0.0

    if param_names is None:
        if model == "szyszkowski":
            param_names = ["K"]
        elif model == "frumkin":
            param_names = ["K", "a"]
        elif model == "linear":
            param_names = ["m", "b"]
        elif model == "exponential_decay":
            param_names = ["amplitude", "rate", "offset"]
        else:
            param_names = [f"p{i}" for i in range(k)]

    caveat = (
        "SEs from (J^T W J)^-1 * s^2 assume the model is locally correct -- "
        "they do not cover model misspecification; check the lag-1 residual "
        "autocorrelation before trusting them"
    )
    return NonlinearFitResult(
        param_names=list(param_names),
        param_values=[float(v) for v in popt],
        param_sd=[float(v) for v in sd],
        ci_low=[float(v - t_crit * s) for v, s in zip(popt, sd)],
        ci_high=[float(v + t_crit * s) for v, s in zip(popt, sd)],
        covariance=[[float(v) for v in row] for row in cov],
        r_squared=float(r2),
        adj_r_squared=float(adj),
        rmse=float(rmse),
        max_abs_standardized_residual=float(np.max(np.abs(std_resid))) if n else 0.0,
        residual_lag1_autocorr=float(lag1),
        n_points=n,
        n_params=k,
        dof=dof,
        converged=bool(sol.success),
        n_iterations=int(sol.nfev),
        method="levenberg-marquardt-trf (scipy least_squares)",
        params={
            "model": model,
            "fixed": {kk: vv for kk, vv in fixed.items()},
            "weighted_by_sigma": sigma is not None,
            "t_crit_95": t_crit,
            "n_function_evaluations": int(sol.nfev),
            "solver_message": str(sol.message),
        },
        caveat=caveat,
    )


# ---------------------------------------------------------------------------
# I5 -- Monte Carlo error propagation through a toolkit chain
# ---------------------------------------------------------------------------

@dataclass
class UncertaintyPropagationResult:
    output_mean: float
    output_sd: float
    output_median: float
    ci_low: float
    ci_high: float
    ci_level: float
    n_samples_requested: int
    n_samples_used: int
    n_failed: int
    seed: int
    converged: bool
    method: str
    params: dict
    caveat: str


def _sample_input(rng, value: float, sigma: float, dist: str):
    if dist == "normal":
        return float(rng.normal(value, sigma))
    if dist == "uniform":
        return float(rng.uniform(value - sigma, value + sigma))
    raise ValueError(f"unknown distribution {dist!r} -- use 'normal' or 'uniform'")


def propagate_uncertainty(
    function,
    inputs: dict[str, dict],
    n_samples: int = 10000,
    ci: float = 0.95,
    seed: int = 0,
) -> UncertaintyPropagationResult:
    """Monte Carlo error propagation through a function of named inputs.

    inputs: {name: {"value": center, "sigma": +/- uncertainty,
    "dist": "normal"|"uniform"}} (dist defaults to normal).  Each sample
    draws every input, calls ``function(**sampled)``, and collects the
    output distribution -> mean, sd, median, percentile interval.

    Linear chains (f = a*a + b*b style) get the exact propagated sd to MC
    precision; nonlinear chains get the true sampled distribution, which is
    the point -- symmetric +- errors on skewed outputs are a fabrication
    risk this tool removes.  Failed/non-finite draws are counted, not hidden
    (n_failed); >50% failures raises instead of returning nonsense.
    """
    np = _require_calc()
    if not isinstance(inputs, dict) or not inputs:
        raise ValueError("inputs must be a non-empty dict of {name: {value, sigma, ...}}")
    if n_samples < 100:
        raise ValueError("n_samples must be >= 100 for stable percentiles")
    if not (0.5 <= ci < 1.0):
        raise ValueError(f"ci must be in [0.5, 1.0), got {ci}")
    for name, spec in inputs.items():
        if "value" not in spec:
            raise ValueError(f"input {name!r} missing 'value'")
        sigma = float(spec.get("sigma", 0.0))
        if sigma < 0:
            raise ValueError(f"input {name!r}: sigma must be >= 0 (0 = fixed)")
        if sigma > 0:
            dist = spec.get("dist", "normal")
            if dist not in ("normal", "uniform"):
                raise ValueError(
                    f"input {name!r}: dist must be 'normal' or 'uniform', got {dist!r}"
                )

    rng = np.random.default_rng(seed)
    names = list(inputs)
    outputs: list[float] = []
    n_failed = 0
    for _ in range(n_samples):
        draw = {}
        for name in names:
            spec = inputs[name]
            sigma = float(spec.get("sigma", 0.0))
            if sigma == 0.0:
                draw[name] = float(spec["value"])
            else:
                draw[name] = _sample_input(rng, float(spec["value"]), sigma,
                                           spec.get("dist", "normal"))
        try:
            out = float(function(**draw))
        except Exception:  # noqa: BLE001 -- function is caller-supplied and
            # arbitrary; any exception type means this draw is unusable,
            # count as failed and keep sampling.
            n_failed += 1
            continue
        if not math.isfinite(out):
            n_failed += 1
            continue
        outputs.append(out)

    if not outputs:
        raise ValueError("every sample failed -- the chain does not accept these inputs")
    if n_failed > 0.5 * n_samples:
        raise ValueError(
            f"{n_failed}/{n_samples} samples failed -- chain rejects most inputs "
            "(check units/ranges of the sampled inputs)"
        )

    arr = np.asarray(outputs, dtype=float)
    alpha = 1.0 - ci
    return UncertaintyPropagationResult(
        output_mean=float(arr.mean()),
        output_sd=float(arr.std(ddof=1)),
        output_median=float(np.median(arr)),
        ci_low=float(np.quantile(arr, alpha / 2)),
        ci_high=float(np.quantile(arr, 1 - alpha / 2)),
        ci_level=ci,
        n_samples_requested=n_samples,
        n_samples_used=int(arr.size),
        n_failed=int(n_failed),
        seed=seed,
        converged=n_failed <= 0.1 * n_samples,
        method="monte-carlo-propagation",
        params={
            "inputs": {k: dict(v) for k, v in inputs.items()},
            "function": getattr(function, "__name__", "chain"),
        },
        caveat=(
            "input sigmas are user-claimed; the output interval reflects only "
            "those -- and the sampled output distribution's own skewness "
            "(compare sd vs the percentile interval width)"
        ),
    )


# Named chains for MCP/JSON callers (function objects can't cross the wire).
def _chain_gibbs_dg_from_cmc(cmc_mole_fraction: float, aggregation_number: float,
                             temperature_K: float = 298.15) -> float:
    from .thermodynamics import mass_action_free_energy_of_micellization
    return mass_action_free_energy_of_micellization(
        cmc_mole_fraction, aggregation_number, temperature_K)


def _chain_hydrodynamic_radius(diffusion_coefficient_cm2_per_s: float,
                               viscosity_mPas: float,
                               temperature_K: float = 298.15) -> float:
    from .dynamics import hydrodynamic_radius_stokes_einstein
    return hydrodynamic_radius_stokes_einstein(
        diffusion_coefficient_cm2_per_s, viscosity_mPas, temperature_K)


def _chain_debye_length(ionic_strength_M: float, temperature_K: float = 298.15) -> float:
    from .electrostatics import debye_length
    return debye_length(ionic_strength_M, temperature_K)


CHAIN_REGISTRY = {
    "gibbs_dg_from_cmc": (
        _chain_gibbs_dg_from_cmc,
        {
            "cmc_mole_fraction": {"value": 0.01, "sigma": 0.001},
            "aggregation_number": {"value": 60, "sigma": 5},
            "temperature_K": {"value": 298.15, "sigma": 0.5},
        },
    ),
    "hydrodynamic_radius": (
        _chain_hydrodynamic_radius,
        {
            "diffusion_coefficient_cm2_per_s": {"value": 5e-7, "sigma": 5e-8},
            "viscosity_mPas": {"value": 0.89, "sigma": 0.01},
            "temperature_K": {"value": 298.15, "sigma": 0.5},
        },
    ),
    "debye_length": (
        _chain_debye_length,
        {
            "ionic_strength_M": {"value": 0.01, "sigma": 0.001},
            "temperature_K": {"value": 298.15, "sigma": 0.5},
        },
    ),
}


def propagate_uncertainty_chain(
    chain: str,
    overrides: dict[str, dict] | None = None,
    n_samples: int = 10000,
    ci: float = 0.95,
    seed: int = 0,
) -> UncertaintyPropagationResult:
    """Named-chain variant of propagate_uncertainty for MCP/JSON callers.

    chain: key of CHAIN_REGISTRY.  Each entry ships representative inputs
    with realistic sigmas (documented, replaceable); ``overrides`` replaces
    any subset ({name: {"value": .., "sigma": ..}}).  These defaults are
    *demonstration* uncertainties -- a researcher should supply their own.
    """
    if chain not in CHAIN_REGISTRY:
        raise ValueError(f"unknown chain {chain!r} -- choose one of {sorted(CHAIN_REGISTRY)}")
    func, default_inputs = CHAIN_REGISTRY[chain]
    merged = {k: dict(v) for k, v in default_inputs.items()}
    if overrides:
        for k, v in overrides.items():
            if k not in merged:
                raise ValueError(
                    f"chain {chain!r} has no input {k!r}; inputs are {sorted(merged)}"
                )
            merged[k].update(v)
    return propagate_uncertainty(func, merged, n_samples=n_samples, ci=ci, seed=seed)


# ---------------------------------------------------------------------------
# I1 -- Metropolis posterior fit
# ---------------------------------------------------------------------------

@dataclass
class PosteriorResult:
    param_names: list[str]
    median: list[float]
    mean: list[float]
    sd: list[float]
    ci_low: list[float]
    ci_high: list[float]
    ci_level: float
    n_chains: int
    n_samples_per_chain: int
    burn_in: int
    acceptance_rate: list[float]
    r_hat: list[float]
    ess: list[float]
    converged: bool
    n_iterations: int
    method: str
    params: dict
    caveat: str


def mcmc_posterior_fit(
    x: list[float],
    y: list[float],
    model,
    p0: list[float],
    observation_sigma: float,
    param_names: list[str] | None = None,
    n_samples: int = 6000,
    burn_in: int = 1500,
    ci: float = 0.95,
    seed: int = 0,
    bounds: tuple[list[float], list[float]] | None = None,
) -> PosteriorResult:
    """Random-walk Metropolis posterior for any y = model(x, *params) model.

    4 independent chains (seeds seed..seed+3), diagonal-adaptive proposal
    tuned during burn-in to ~30% acceptance, then frozen.  Reports medians,
    credible intervals, per-chain acceptance, split-R-hat, and ESS -- the
    diagnostics that separate a real posterior from a made-up +- range.

    observation_sigma: known per-point observation noise (same units as y).
    This is the honest setup for "my instrument has this uncertainty"; the
    worker does NOT silently estimate sigma (use bootstrap_cmc / nonlinear
    fit residuals for that and state it).
    bounds: optional flat-prior box (uniform outside -> posterior = 0).
    """
    np = _require_calc()
    if not (0.5 <= ci < 1.0):
        raise ValueError(f"ci must be in [0.5, 1.0), got {ci}")
    if observation_sigma <= 0:
        raise ValueError("observation_sigma must be positive (known measurement noise)")
    if n_samples < 400:
        raise ValueError("n_samples must be >= 400 per chain")
    if burn_in < 100:
        raise ValueError("burn_in must be >= 100")
    if burn_in >= n_samples:
        raise ValueError("burn_in must be < n_samples")
    xs = np.asarray(x, dtype=float)
    ys = np.asarray(y, dtype=float)
    if xs.shape != ys.shape or xs.ndim != 1:
        raise ValueError("x and y must be same-length 1-D sequences")
    if xs.size < 4:
        raise ValueError("need at least 4 points to constrain a posterior")

    p0_arr = np.asarray(p0, dtype=float)
    k = p0_arr.size
    lo = -np.inf * np.ones(k) if bounds is None else np.asarray(bounds[0], float)
    hi = np.inf * np.ones(k) if bounds is None else np.asarray(bounds[1], float)
    if np.any(p0_arr <= lo) or np.any(p0_arr >= hi):
        raise ValueError("p0 must lie strictly inside bounds")

    if param_names is None:
        param_names = [f"p{i}" for i in range(k)]
    if len(param_names) != k:
        raise ValueError("param_names must match len(p0)")

    def log_post(p):
        if np.any(p < lo) or np.any(p > hi) or not np.all(np.isfinite(p)):
            return -np.inf
        try:
            yhat = np.asarray(model(xs, *p), dtype=float)
        except Exception:  # noqa: BLE001 -- model is caller-supplied and
            # arbitrary; any exception means this parameter draw is
            # infeasible, reject it (log-posterior -inf) rather than crash.
            return -np.inf
        if yhat.shape != ys.shape or not np.all(np.isfinite(yhat)):
            return -np.inf
        r = (ys - yhat) / observation_sigma
        return -0.5 * float(np.sum(r * r))

    # proposal scale: start at 1% of (bounded) parameter scale
    span = np.where(np.isfinite(hi - lo), hi - lo, np.maximum(np.abs(p0_arr), 1.0))
    step = 0.01 * span.copy()
    n_chains = 4
    all_draws = []
    accept_rates = []
    for chain in range(n_chains):
        rng = np.random.default_rng(seed + chain)
        p = p0_arr * (1.0 + 0.05 * chain) if chain else p0_arr.copy()
        p = np.clip(p, lo + 1e-12 * span, hi - 1e-12 * span)
        lp = log_post(p)
        if not math.isfinite(lp):
            # jitter start until inside the posterior support
            for _ in range(1000):
                p = p0_arr + rng.normal(0.0, 0.01 * span)
                lp = log_post(p)
                if math.isfinite(lp):
                    break
            if not math.isfinite(lp):
                raise ValueError("model returns non-finite likelihood at p0 -- check data/model")
        draws = np.empty((n_samples, k))
        accepted = 0
        adapt_window = max(50, burn_in // 10)
        for it in range(n_samples):
            prop = p + rng.normal(0.0, step)
            lp_prop = log_post(prop)
            if math.log(rng.random() + 1e-300) < lp_prop - lp:
                p, lp = prop, lp_prop
                accepted += 1
            if it < burn_in and it % adapt_window == 0 and it > 0:
                rate = accepted / (it + 1)
                step = step * (1.25 if rate > 0.35 else (0.8 if rate < 0.2 else 1.0))
            draws[it] = p
        all_draws.append(draws[burn_in:])
        accept_rates.append(accepted / n_samples)

    pooled = np.vstack(all_draws)  # (n_chains * kept, k)
    alpha = 1.0 - ci
    medians, means, sds, lows, highs, r_hats, ess_list = [], [], [], [], [], [], []
    for j in range(k):
        chains = [d[:, j] for d in all_draws]
        np.mean(chains)
        # split-R-hat: each chain -> first/second half -> 8 sequences
        split = []
        for c in chains:
            half = len(c) // 2
            split.append(c[:half])
            split.append(c[half:2 * half])
        seq_means = [float(np.mean(s)) for s in split]
        seq_vars = [float(np.var(s, ddof=1)) for s in split]
        M = len(split)
        W = float(np.mean(seq_vars))
        B = M * float(np.var(seq_means, ddof=1))
        n_seq = len(split[0])
        if W <= 0:
            r_hat = 1.0
        else:
            var_hat = (n_seq - 1) / n_seq * W + B / n_seq
            r_hat = math.sqrt(var_hat / W)
        # ESS via initial-positive-sequence autocorrelation on pooled draws
        z = pooled[:, j] - pooled[:, j].mean()
        v = float(np.dot(z, z)) / len(z)
        if v > 0:
            rho_sum = 0.0
            max_lag = min(len(z) - 1, 400)
            prev_rho = 1.0
            for lag in range(1, max_lag + 1):
                rho = float(np.dot(z[:-lag], z[lag:]) / (len(z) * v))
                if rho < 0:
                    break
                rho = min(rho, prev_rho)
                rho_sum += rho
                prev_rho = rho
            ess = len(z) / (1.0 + 2.0 * rho_sum)
        else:
            ess = float(len(z))
        ess = max(1.0, min(float(len(z)), ess))
        col = pooled[:, j]
        medians.append(float(np.median(col)))
        means.append(float(col.mean()))
        sds.append(float(col.std(ddof=1)))
        lows.append(float(np.quantile(col, alpha / 2)))
        highs.append(float(np.quantile(col, 1 - alpha / 2)))
        r_hats.append(float(r_hat))
        ess_list.append(float(ess))

    converged = all(r < 1.1 for r in r_hats) and all(e >= 100 for e in ess_list)
    return PosteriorResult(
        param_names=list(param_names),
        median=medians,
        mean=means,
        sd=sds,
        ci_low=lows,
        ci_high=highs,
        ci_level=ci,
        n_chains=n_chains,
        n_samples_per_chain=n_samples,
        burn_in=burn_in,
        acceptance_rate=[float(a) for a in accept_rates],
        r_hat=r_hats,
        ess=ess_list,
        converged=bool(converged),
        n_iterations=n_chains * n_samples,
        method="random-walk metropolis (4 chains, adaptive-diagonal proposal)",
        params={
            "observation_sigma": observation_sigma,
            "seed": seed,
            "proposal_steps_final": [float(s) for s in step],
            "draws": [[float(v) for v in row] for row in
                      pooled[:: max(1, len(pooled) // 200)]],  # ~200 corner points
            "model": getattr(model, "__name__", "model"),
        },
        caveat=(
            "credible intervals assume observation_sigma is correct and the "
            "model is true; split-R-hat < 1.1 and ESS >= 100 are reported per "
            "parameter -- do not quote an interval whose chain diagnostics failed"
        ),
    )


# ---------------------------------------------------------------------------
# I4 -- global multi-dataset fit (shared parameters)
# ---------------------------------------------------------------------------

@dataclass
class GlobalFitResult:
    shared_names: list[str]
    shared_values: list[float]
    shared_sd: list[float]
    shared_ci_low: list[float]
    shared_ci_high: list[float]
    local_values: list[list[float]]
    local_names: list[str]
    r_squared_per_dataset: list[float]
    r_squared_global: float
    chi2: float
    dof: int
    n_points_total: int
    converged: bool
    n_iterations: int
    method: str
    params: dict
    caveat: str


def global_multi_dataset_fit(
    datasets: list[dict],
    model,
    shared_p0: list[float],
    local_p0: list[float],
    shared_names: list[str] | None = None,
    local_names: list[str] | None = None,
    sigma: float | None = None,
    bounds: tuple[list[float], list[float]] | None = None,
    max_nfev: int = 6000,
) -> GlobalFitResult:
    """Fit ONE shared parameter vector across several datasets simultaneously.

    datasets: [{"x": [...], "y": [...]}, ...].  model(x, *shared, *local) is
    evaluated per dataset with the SAME shared vector and that dataset's own
    local vector -- the standard global-fit construction (one beta across a
    composition series with per-series offsets, one D across temperatures).

    Sharing shrinks the shared parameter's SE by the extra datasets' own
    information (nothing is assumed); local parameters absorb per-dataset
    differences.  sigma: common known observation noise for chi2 (optional).
    """
    np = _require_calc()
    from scipy.optimize import least_squares

    if not datasets:
        raise ValueError("datasets must be a non-empty list of {x, y} dicts")
    n_ds = len(datasets)
    sh0 = [float(v) for v in shared_p0]
    lo0 = [float(v) for v in local_p0]
    n_sh, n_lo = len(sh0), len(lo0)
    if n_sh == 0:
        raise ValueError("shared_p0 must have at least one shared parameter")
    xs_all, ys_all, ds_index = [], [], []
    for i, ds in enumerate(datasets):
        xi = np.asarray(ds["x"], dtype=float)
        yi = np.asarray(ds["y"], dtype=float)
        if xi.shape != yi.shape or xi.ndim != 1 or xi.size == 0:
            raise ValueError(f"dataset {i}: x and y must be same-length non-empty 1-D sequences")
        xs_all.append(xi)
        ys_all.append(yi)
        ds_index.extend([i] * xi.size)
    np.asarray(ds_index)
    y_concat = np.concatenate(ys_all)
    n_total = y_concat.size
    n_params = n_sh + n_lo * n_ds
    dof = n_total - n_params
    if dof <= 0:
        raise ValueError(
            f"no degrees of freedom: {n_total} points, {n_params} parameters "
            f"({n_sh} shared + {n_lo} local x {n_ds} datasets)"
        )

    if bounds is None:
        p_lo = -np.inf * np.ones(n_params)
        p_hi = np.inf * np.ones(n_params)
    else:
        p_lo = np.asarray(bounds[0], float)
        p_hi = np.asarray(bounds[1], float)
        if p_lo.size != n_params or p_hi.size != n_params:
            raise ValueError(f"bounds must have length {n_params} = n_shared + n_local*n_datasets")

    def _resid(theta):
        shared = theta[:n_sh]
        out = np.empty(n_total)
        offset = 0
        for i, (xi, yi) in enumerate(zip(xs_all, ys_all)):
            local = theta[n_sh + i * n_lo: n_sh + (i + 1) * n_lo]
            pred = np.asarray(model(xi, *shared, *local), dtype=float)
            out[offset:offset + xi.size] = pred - yi
            offset += xi.size
        return out

    theta0 = np.asarray(sh0 + lo0 * n_ds, dtype=float)
    sol = least_squares(_resid, theta0, bounds=(p_lo, p_hi), max_nfev=max_nfev)
    resid = _resid(sol.x)
    rss = float(np.sum(resid ** 2))
    jac = sol.jac
    jtj = jac.T @ jac
    try:
        cov = np.linalg.inv(jtj) * (rss / dof)
    except np.linalg.LinAlgError:
        cov = np.linalg.pinv(jtj) * (rss / dof)
    sd = np.sqrt(np.maximum(np.diag(cov), 0.0))
    from scipy import stats
    t_crit = float(stats.t.ppf(0.975, dof))

    y_mean = float(y_concat.mean())
    ss_tot = float(np.sum((y_concat - y_mean) ** 2))
    r2_global = 1.0 - rss / ss_tot if ss_tot > 0 else 0.0
    r2_per = []
    offset = 0
    for i, yi in enumerate(ys_all):
        r = resid[offset:offset + yi.size]
        offset += yi.size
        st = float(np.sum((yi - yi.mean()) ** 2))
        r2_per.append(1.0 - float(np.sum(r ** 2)) / st if st > 0 else 0.0)
    chi2 = float(np.sum((resid / sigma) ** 2)) if sigma else rss
    if sigma is not None and sigma <= 0:
        raise ValueError("sigma must be positive")

    if shared_names is None:
        shared_names = [f"shared_{j}" for j in range(n_sh)]
    if local_names is None:
        local_names = [f"local_{j}" for j in range(n_lo)]
    if len(shared_names) != n_sh or len(local_names) != n_lo:
        raise ValueError("name lists must match shared_p0/local_p0 lengths")

    shared_vals = [float(v) for v in sol.x[:n_sh]]
    local_vals = [[float(v) for v in sol.x[n_sh + i * n_lo: n_sh + (i + 1) * n_lo]]
                  for i in range(n_ds)]
    caveat = (
        "global-fit SEs assume one common observation noise and independent "
        "datasets; sharing parameters that the datasets do not actually share "
        "biases every result -- compare per-dataset residuals first"
    )
    return GlobalFitResult(
        shared_names=list(shared_names),
        shared_values=shared_vals,
        shared_sd=[float(v) for v in sd[:n_sh]],
        shared_ci_low=[float(v - t_crit * s) for v, s in zip(sol.x[:n_sh], sd[:n_sh])],
        shared_ci_high=[float(v + t_crit * s) for v, s in zip(sol.x[:n_sh], sd[:n_sh])],
        local_values=local_vals,
        local_names=list(local_names),
        r_squared_per_dataset=[float(v) for v in r2_per],
        r_squared_global=float(r2_global),
        chi2=float(chi2),
        dof=dof,
        n_points_total=n_total,
        converged=bool(sol.success),
        n_iterations=int(sol.nfev),
        method="global concatenated least_squares (trf)",
        params={
            "n_datasets": n_ds,
            "n_shared": n_sh,
            "n_local": n_lo,
            "sigma_assumed": sigma,
            "t_crit_95": t_crit,
            "solver_message": str(sol.message),
        },
        caveat=caveat,
    )
