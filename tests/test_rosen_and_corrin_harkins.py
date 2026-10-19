"""Tests for the Rosen monolayer interaction parameter (mixed_micelle.py)
and the Corrin-Harkins salt-CMC relation."""

import math

import pytest

from surfactantkit.mixed_micelle import (
    corrin_harkins_predict_cmc,
    rosen_beta_sigma,
    rubingh_beta,
    solve_rosen_monolayer_x,
    solve_rubingh_x,
)


def test_rosen_monolayer_matches_rubingh_math_exactly():
    """Rosen's monolayer theory is mathematically identical to Rubingh's
    micellization theory -- the same inputs must give the same outputs,
    since solve_rosen_monolayer_x is a thin, documented alias, not a
    reimplementation. If this ever diverges, something broke."""
    alpha1, c_mix, c1, c2 = 0.5, 4.07, 11.50, 11.98  # reusing the cholate-SDS numbers
    x1_micelle = solve_rubingh_x(alpha1, c_mix, c1, c2)
    x1_monolayer = solve_rosen_monolayer_x(alpha1, c_mix, c1, c2)
    assert x1_monolayer == pytest.approx(x1_micelle)

    beta_micelle = rubingh_beta(x1_micelle, alpha1, c_mix, c1)
    beta_monolayer = rosen_beta_sigma(x1_monolayer, alpha1, c_mix, c1)
    assert beta_monolayer == pytest.approx(beta_micelle)


def test_corrin_harkins_round_trip():
    """Construct two (CMC, salt-concentration) points consistent with a
    chosen g and intercept, then verify the fit-and-predict function
    recovers a third point on the same line exactly -- mathematically
    guaranteed, no external literature dependency."""
    g_true = 0.6
    b_true = 1.0  # log10(CMC) = -g*log10(C) + b

    c1, c2, c_target = 0.01, 0.1, 0.05
    cmc1 = 10.0 ** (-g_true * math.log10(c1) + b_true)
    cmc2 = 10.0 ** (-g_true * math.log10(c2) + b_true)
    cmc_target_true = 10.0 ** (-g_true * math.log10(c_target) + b_true)

    cmc_predicted, g_fitted = corrin_harkins_predict_cmc(cmc1, c1, cmc2, c2, c_target)
    assert g_fitted == pytest.approx(g_true, abs=1e-9)
    assert cmc_predicted == pytest.approx(cmc_target_true, rel=1e-9)


def test_corrin_harkins_cmc_decreases_with_salt_for_positive_g():
    """Physical sanity check: for typical ionic surfactants, g > 0, and
    CMC should decrease as counterion (salt) concentration increases."""
    cmc_low_salt, _ = corrin_harkins_predict_cmc(
        cmc1=10.0, c_counterion1=0.01, cmc2=5.0, c_counterion2=0.1, c_counterion_target=0.02
    )
    cmc_high_salt, _ = corrin_harkins_predict_cmc(
        cmc1=10.0, c_counterion1=0.01, cmc2=5.0, c_counterion2=0.1, c_counterion_target=0.08
    )
    assert cmc_high_salt < cmc_low_salt


def test_corrin_harkins_rejects_equal_concentrations():
    with pytest.raises(ValueError):
        corrin_harkins_predict_cmc(10.0, 0.05, 5.0, 0.05, 0.02)


def test_corrin_harkins_rejects_nonpositive_input():
    with pytest.raises(ValueError):
        corrin_harkins_predict_cmc(0.0, 0.01, 5.0, 0.1, 0.02)
