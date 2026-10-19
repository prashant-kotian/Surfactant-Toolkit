"""Direct unit tests for kinetics.mass_action_equilibrium.

Regression driver: overflow fix (2026-10-05) -- K^(n-1)*x1^n evaluated as
the grouped x1*(K*x1)^(n-1). Physical premicellar-census inputs
(K = 1/CMC ~ 1e4 M^-1 for C12E8, n_max = 100) used to raise
OverflowError (K^99 = 1e396 > float64 max) even though every term of the
distribution is perfectly representable.
"""
import pytest

from surfactantkit.kinetics import mass_action_equilibrium


def test_large_K_n_max_no_overflow_and_mass_balance():
    c_total, K, n_max = 1e-3, 1e4, 100
    xn = mass_action_equilibrium(c_total, K, n_max)
    assert len(xn) == n_max
    assert all(v >= 0.0 for v in xn)
    mb = sum((i + 1) * v for i, v in enumerate(xn))
    assert mb == pytest.approx(c_total, rel=1e-9)


def test_grouped_form_matches_naive_form_where_defined():
    # same inputs below the naive form's overflow limit -> identical results
    c_total, K, n_max = 1e-4, 1e2, 40
    xn = mass_action_equilibrium(c_total, K, n_max)
    x1 = xn[0]
    for n in (2, 5, 40):
        assert xn[n - 1] == pytest.approx(
            (K ** (n - 1)) * (x1 ** n), rel=1e-12)


def test_fail_fast_contract_unchanged():
    with pytest.raises(ValueError, match="positive"):
        mass_action_equilibrium(0.0, 1.0, 10)
    with pytest.raises(ValueError, match="positive"):
        mass_action_equilibrium(1e-3, -1.0, 10)
    with pytest.raises(ValueError, match=">= 3"):
        mass_action_equilibrium(1e-3, 1.0, 2)
