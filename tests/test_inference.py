"""Numeric validation for surfactantkit.inference (I1 mcmc_posterior_fit,
I3 nonlinear_fit_with_errors, I4 global_multi_dataset_fit,
I5 propagate_uncertainty_chain) and homologous.cmc_homologous_series (C2).

Golds are independent of the code under test:
- linear least squares: parameter values and standard errors against the
  closed-form ordinary-least-squares solution (numpy lstsq + s^2 (X^T X)^-1);
- noise-free exact recovery for the nonlinear, global and homologous fits;
- uncertainty propagation against the first-order (delta-method) estimate
  for the Debye length, kappa^-1 proportional to I^-1/2;
- MCMC posterior median and sd against the OLS estimate and SE (flat prior,
  known sigma: the posterior is the OLS Gaussian), plus R-hat and
  seed determinism.
"""

import math

import numpy as np
import pytest

from surfactantkit.electrostatics import debye_length
from surfactantkit.homologous import cmc_homologous_series
from surfactantkit.inference import (
    _model_linear,
    global_multi_dataset_fit,
    mcmc_posterior_fit,
    nonlinear_fit_with_errors,
    propagate_uncertainty_chain,
)


def _noisy_line(m=2.0, b=1.0, n=12, noise=0.3, seed=5):
    rng = np.random.default_rng(seed)
    x = np.linspace(0.0, 5.0, n)
    y = m * x + b + rng.normal(0.0, noise, n)
    return x, y


def _ols(x, y):
    X = np.column_stack([x, np.ones_like(x)])
    coef, res, *_ = np.linalg.lstsq(X, y, rcond=None)
    dof = len(x) - 2
    s2 = float(res[0]) / dof
    se = np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X)))
    return coef, se


def test_nonlinear_fit_linear_matches_closed_form_ols():
    x, y = _noisy_line()
    r = nonlinear_fit_with_errors(x.tolist(), y.tolist(), "linear")
    coef, se = _ols(x, y)
    assert r.converged
    assert r.param_values == pytest.approx(coef.tolist(), rel=1e-6)
    assert r.param_sd == pytest.approx(se.tolist(), rel=1e-4)
    assert r.dof == len(x) - 2
    # the 95% interval brackets the estimate symmetrically
    for v, lo, hi in zip(r.param_values, r.ci_low, r.ci_high, strict=True):
        assert lo < v < hi
        assert (v - lo) == pytest.approx(hi - v, rel=1e-9)


def test_nonlinear_fit_exponential_decay_exact_recovery():
    x = np.linspace(0.0, 4.0, 15)
    y = 3.0 * np.exp(-1.7 * x) + 0.4
    r = nonlinear_fit_with_errors(x.tolist(), y.tolist(), "exponential_decay",
                                  p0=[1.0, 1.0, 0.0])
    assert r.converged
    assert r.param_values == pytest.approx([3.0, 1.7, 0.4], rel=1e-6)
    assert r.r_squared == pytest.approx(1.0, abs=1e-12)


def test_nonlinear_fit_unknown_model_refused():
    with pytest.raises(ValueError, match="unknown model"):
        nonlinear_fit_with_errors([1, 2, 3], [1, 2, 3], "cubic")


def test_global_fit_recovers_shared_slope_and_local_intercepts():
    x = np.linspace(0.0, 5.0, 8)
    datasets = [{"x": x.tolist(), "y": (2.0 * x + 1.0).tolist()},
                {"x": x.tolist(), "y": (2.0 * x + 5.0).tolist()}]
    r = global_multi_dataset_fit(datasets, _model_linear, shared_p0=[1.0],
                                 local_p0=[0.0])
    assert r.converged
    assert r.shared_values == pytest.approx([2.0], rel=1e-8)
    assert [v[0] for v in r.local_values] == pytest.approx([1.0, 5.0], rel=1e-8)
    assert r.r_squared_global == pytest.approx(1.0, abs=1e-12)


def test_propagation_debye_length_matches_delta_method():
    # kappa^-1 ~ I^-1/2, so sd_rel = 0.5 * sd_rel(I) = 0.05 for a 10 % input;
    # temperature fixed (sigma ~ 0) to isolate the ionic-strength term
    r = propagate_uncertainty_chain(
        "debye_length",
        overrides={"ionic_strength_M": {"value": 0.01, "sigma": 0.001},
                   "temperature_K": {"value": 298.15, "sigma": 1e-9}},
        n_samples=20000, seed=1)
    centre = debye_length(0.01, 298.15)
    assert r.output_median == pytest.approx(centre, rel=0.01)
    assert r.output_sd / centre == pytest.approx(0.05, rel=0.10)
    assert r.ci_low < centre < r.ci_high


def test_propagation_is_seed_deterministic_and_refuses_unknown_chain():
    a = propagate_uncertainty_chain("hydrodynamic_radius", n_samples=2000, seed=3)
    b = propagate_uncertainty_chain("hydrodynamic_radius", n_samples=2000, seed=3)
    assert a == b
    with pytest.raises(ValueError, match="unknown chain"):
        propagate_uncertainty_chain("cmc_to_area")


def test_mcmc_posterior_matches_ols_gaussian():
    x, y = _noisy_line(noise=0.3)
    coef, _ = _ols(x, y)
    # known-sigma posterior sd = sigma * sqrt(diag((X^T X)^-1))
    X = np.column_stack([x, np.ones_like(x)])
    sd_known = 0.3 * np.sqrt(np.diag(np.linalg.inv(X.T @ X)))
    r = mcmc_posterior_fit(x.tolist(), y.tolist(), _model_linear, p0=[1.0, 0.0],
                           observation_sigma=0.3, n_samples=6000, burn_in=1500,
                           seed=2)
    assert r.converged
    assert max(r.r_hat) < 1.1
    for k in range(2):
        assert r.median[k] == pytest.approx(coef[k], abs=0.25 * sd_known[k])
        assert r.sd[k] == pytest.approx(sd_known[k], rel=0.20)
    again = mcmc_posterior_fit(x.tolist(), y.tolist(), _model_linear,
                               p0=[1.0, 0.0], observation_sigma=0.3,
                               n_samples=6000, burn_in=1500, seed=2)
    assert again.median == r.median


def test_homologous_series_exact_recovery_and_extrapolation_flag():
    ns = [8, 10, 12, 14]
    cmcs = [10 ** (4.0 - 0.30 * n) for n in ns]
    r = cmc_homologous_series(ns, cmcs, unit="mM", predict_at=[12, 16])
    assert r.intercept_A == pytest.approx(4.0, abs=1e-9)
    assert r.reduction_per_ch2 == pytest.approx(0.30, abs=1e-9)
    assert r.reduction_factor_per_ch2 == pytest.approx(10 ** 0.30, rel=1e-9)
    by_n = {p["chain_length"]: p for p in r.predicted}
    assert by_n[16]["cmc"] == pytest.approx(10 ** (4.0 - 4.8), rel=1e-9)
    assert by_n[16]["extrapolation"] is True
    assert by_n[12]["extrapolation"] is False


def test_homologous_series_sodium_alkyl_sulfates_klevens_slope():
    # Mukerjee & Mysels (NSRDS-NBS 36, 1971) 25 C values, mM:
    # C8 130, C10 33, C12 8.2, C14 2.1 -> log10 drop ~0.30 per CH2
    r = cmc_homologous_series([8, 10, 12, 14], [130.0, 33.0, 8.2, 2.1])
    assert r.reduction_per_ch2 == pytest.approx(0.30, abs=0.01)
    assert r.r_squared > 0.999
    assert math.isfinite(r.slope_se)


def test_homologous_series_rejects_unknown_unit():
    with pytest.raises(ValueError, match="unit"):
        cmc_homologous_series([8, 10, 12], [130.0, 33.0, 8.2], unit="molal")
