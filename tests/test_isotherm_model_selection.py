"""Tests for surfactantkit.adsorption.select_isotherm_model -- the
Langmuir-vs-Frumkin model-selection step. Real gap this closes: both
szyszkowski_fit_K and frumkin_fit_K_and_a will each fit their OWN model to
any curve (Frumkin, having an extra free parameter, always achieves an
equal-or-lower raw RSS -- that alone proves nothing about which model is
actually correct). BIC penalizes the extra parameter, the same
model-selection principle and formula already used by
curve_analysis.cmc_from_surface_tension_curve's own 2-vs-3-segment choice.
Validated on noisy synthetic round trips (a real researcher's data always
has noise; this is the case that actually exercises the selection logic,
unlike an exact/noiseless curve where Frumkin trivially always wins)."""

import random

import pytest

from surfactantkit.adsorption import frumkin_surface_tension, select_isotherm_model


def _noisy(vals, sd, seed):
    rng = random.Random(seed)
    return [v + rng.gauss(0, sd) for v in vals]


def test_select_isotherm_model_picks_langmuir_for_true_a_zero():
    gamma0, gamma_max, system_type, T = 72.0, 3.5e-6, "nonionic", 298.15
    concentrations = [0.02, 0.05, 0.1, 0.2, 0.4, 0.7, 1.0, 1.4, 1.8, 2.2]
    exact = [frumkin_surface_tension(c, gamma0, gamma_max, 15.0, 0.0, system_type, T) for c in concentrations]
    tensions = _noisy(exact, sd=0.05, seed=0)

    result = select_isotherm_model(concentrations, tensions, gamma0, gamma_max, system_type, T)
    assert result.selected_model == "langmuir"
    assert result.frumkin_a == pytest.approx(0.0, abs=0.1)


def test_select_isotherm_model_picks_frumkin_for_real_interaction():
    gamma0, gamma_max, system_type, T = 72.0, 3.5e-6, "nonionic", 298.15
    concentrations = [0.02, 0.05, 0.1, 0.2, 0.4, 0.7, 1.0, 1.4, 1.8, 2.2]
    exact = [frumkin_surface_tension(c, gamma0, gamma_max, 15.0, 1.5, system_type, T) for c in concentrations]
    tensions = _noisy(exact, sd=0.05, seed=0)

    result = select_isotherm_model(concentrations, tensions, gamma0, gamma_max, system_type, T)
    assert result.selected_model == "frumkin"
    assert result.frumkin_a == pytest.approx(1.5, rel=0.15)
    assert result.delta_bic > 0  # positive delta_bic means Frumkin favored


def test_select_isotherm_model_returns_both_fits_r_squared():
    gamma0, gamma_max, system_type, T = 72.0, 3.5e-6, "nonionic", 298.15
    concentrations = [0.02, 0.05, 0.1, 0.2, 0.4, 0.7, 1.0, 1.4]
    tensions = [frumkin_surface_tension(c, gamma0, gamma_max, 10.0, 0.0, system_type, T) for c in concentrations]

    result = select_isotherm_model(concentrations, tensions, gamma0, gamma_max, system_type, T)
    assert result.langmuir_r_squared == pytest.approx(1.0, abs=1e-4)
    assert result.frumkin_r_squared == pytest.approx(1.0, abs=1e-4)
