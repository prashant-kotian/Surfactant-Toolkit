"""Tests for surfactantkit.curve_analysis.bootstrap_cmc (worker I2).

Gold: synthetic sigmoid with known midpoint, fixed seed -- point estimate
accurate, percentile interval covers truth, deterministic across runs,
input validation refuses under-determined data.
"""
import numpy as np
import pytest

from surfactantkit.curve_analysis import BootstrapCMCResult, bootstrap_cmc


def _synthetic_st(true_log10_cmc=2.1, n=16, noise=0.15, seed=3,
                   rising=False):
    rng = np.random.default_rng(seed)
    x = np.linspace(1.6, 2.6, n)
    b0, b2, b3 = 2.9, 0.3, 40.0
    core = b0 * (1.0 - np.tanh((x - true_log10_cmc) / b2))
    y = (b3 + core) if not rising else (b3 + 40.0 - core)
    y = y + rng.normal(0.0, noise, n)
    return x.tolist(), y.tolist()


def test_point_estimate_and_coverage():
    x, y = _synthetic_st()
    r = bootstrap_cmc(x, y, n_boot=400, seed=7)
    assert isinstance(r, BootstrapCMCResult)
    true_cmc = 10 ** 2.1
    assert r.cmc_point == pytest.approx(true_cmc, rel=0.05)
    assert r.cmc_ci_low < true_cmc < r.cmc_ci_high
    # interval must be tight enough to be informative but not degenerate
    assert 0.02 < (r.log10_cmc_ci_high - r.log10_cmc_ci_low) < 0.5
    assert r.converged and r.n_failed <= 40
    assert r.n_boot_used == 400
    assert r.r_squared > 0.95


def test_deterministic_for_same_seed():
    x, y = _synthetic_st()
    a = bootstrap_cmc(x, y, n_boot=200, seed=11)
    b = bootstrap_cmc(x, y, n_boot=200, seed=11)
    assert a == b
    c = bootstrap_cmc(x, y, n_boot=200, seed=12)
    assert c.cmc_ci_low != a.cmc_ci_low or c.cmc_ci_high != a.cmc_ci_high


def test_rising_curve_orientation():
    # conductivity-like sigmoid rising with x: midpoint must still recover
    rng = np.random.default_rng(5)
    x = np.linspace(-3.0, -1.0, 18)
    y = 200 + 600 * 0.5 * (1 + np.tanh((x - (-2.0)) / 0.15))
    y = y + rng.normal(0, 8, x.size)
    r = bootstrap_cmc(x.tolist(), y.tolist(), n_boot=300, seed=1)
    assert r.log10_cmc_point == pytest.approx(-2.0, abs=0.06)
    # Very low noise (r2 ~ 0.999) => extremely tight interval; the percentile
    # bootstrap does not correct the estimator's tiny finite-sample bias, so
    # allow coverage slack of 0.01 log units rather than assert strict
    # containment (strict containment IS asserted in the realistic-noise test).
    assert r.log10_cmc_ci_low <= -2.0 + 0.01
    assert r.log10_cmc_ci_high >= -2.0 - 0.01
    assert (r.log10_cmc_ci_high - r.log10_cmc_ci_low) < 0.05


def test_wider_ci_level_gives_wider_interval():
    x, y = _synthetic_st()
    r90 = bootstrap_cmc(x, y, n_boot=300, ci=0.90, seed=4)
    r99 = bootstrap_cmc(x, y, n_boot=300, ci=0.99, seed=4)
    w90 = r90.log10_cmc_ci_high - r90.log10_cmc_ci_low
    w99 = r99.log10_cmc_ci_high - r99.log10_cmc_ci_low
    assert w99 > w90


def test_validation_errors():
    x, y = _synthetic_st(n=5)
    with pytest.raises(ValueError, match=">= 8 points"):
        bootstrap_cmc(x, y, n_boot=100)
    x, y = _synthetic_st()
    with pytest.raises(ValueError, match="n_boot"):
        bootstrap_cmc(x, y, n_boot=10)
    with pytest.raises(ValueError, match="ci must be"):
        bootstrap_cmc(x, y, n_boot=100, ci=0.4)
    with pytest.raises(ValueError, match="distinct"):
        bootstrap_cmc([1.0] * 10, list(range(10)), n_boot=100)
    with pytest.raises(ValueError, match="finite"):
        bootstrap_cmc([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0],
                      [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, float("nan")],
                      n_boot=100)
    # flat curve -> no transition
    with pytest.raises(ValueError, match="no pre/post transition"):
        bootstrap_cmc(list(np.linspace(1.6, 2.6, 12)),
                      [40.0] * 12, n_boot=100)


def test_result_fields_and_caveat():
    x, y = _synthetic_st()
    r = bootstrap_cmc(x, y, n_boot=100, seed=1)
    d = r.__dict__
    for key in ("cmc_point", "cmc_ci_low", "cmc_ci_high", "ci_level",
                "n_boot_used", "n_failed", "converged", "r_squared",
                "method", "params", "caveat"):
        assert key in d
    assert "model-form error" in r.caveat
    assert r.ci_level == 0.95
