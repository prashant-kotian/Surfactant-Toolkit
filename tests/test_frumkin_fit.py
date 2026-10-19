"""Tests for surfactantkit.adsorption.frumkin_fit_K_and_a -- closes a real
gap: frumkin_surface_tension only predicts forward from K and a already
known, and szyszkowski_fit_K only fits K alone under the implicit a=0
assumption. Neither could determine, from raw data alone, whether a
surfactant's adsorption is better described as plain Langmuir or genuinely
non-ideal (Frumkin). Validated via mathematically-guaranteed round trips
(construct a curve from chosen true K and a via frumkin_surface_tension
itself, fit both back), same disclosed pattern as test_szyszkowski_fit_k.py.
"""

import pytest

from surfactantkit.adsorption import frumkin_fit_K_and_a, frumkin_surface_tension


def test_frumkin_fit_recovers_true_k_and_a_round_trip_attractive():
    true_K, true_a = 25.0, 1.3
    gamma0, gamma_max, system_type, T = 72.0, 3.5e-6, "ionic_no_added_salt", 298.15
    concentrations = [0.05, 0.1, 0.2, 0.4, 0.8, 1.2, 1.6, 2.0]
    tensions = [frumkin_surface_tension(c, gamma0, gamma_max, true_K, true_a, system_type, T) for c in concentrations]

    result = frumkin_fit_K_and_a(concentrations, tensions, gamma0, gamma_max, system_type, T)
    assert result.K == pytest.approx(true_K, rel=1e-2)
    assert result.a == pytest.approx(true_a, rel=1e-2)
    assert result.r_squared == pytest.approx(1.0, abs=1e-6)
    assert not result.a_pinned_at_bound


def test_frumkin_fit_recovers_near_zero_a_langmuir_limit():
    true_K, true_a = 10.0, 0.02
    gamma0, gamma_max, system_type, T = 72.0, 3.5e-6, "ionic_no_added_salt", 298.15
    concentrations = [0.05, 0.1, 0.2, 0.4, 0.8, 1.2, 1.6, 2.0]
    tensions = [frumkin_surface_tension(c, gamma0, gamma_max, true_K, true_a, system_type, T) for c in concentrations]

    result = frumkin_fit_K_and_a(concentrations, tensions, gamma0, gamma_max, system_type, T)
    assert result.K == pytest.approx(true_K, rel=1e-2)
    assert result.a == pytest.approx(true_a, abs=0.05)


def test_frumkin_fit_flags_condensation_regime_warning():
    """|a| >= 2 is a real physical regime (see frumkin_theta's docstring),
    not a fitting artifact -- the fit must surface a warning, not hide it."""
    true_K, true_a = 15.0, 2.5
    gamma0, gamma_max, system_type, T = 72.0, 3.5e-6, "nonionic", 298.15
    concentrations = [0.05, 0.1, 0.2, 0.4, 0.8, 1.2, 1.6, 2.0]
    tensions = [frumkin_surface_tension(c, gamma0, gamma_max, true_K, true_a, system_type, T) for c in concentrations]

    result = frumkin_fit_K_and_a(concentrations, tensions, gamma0, gamma_max, system_type, T)
    assert any("condensation" in w for w in result.warnings)


def test_frumkin_fit_flags_pinned_at_bound_when_true_a_outside_bounds():
    true_K, true_a = 15.0, 5.0  # outside default a_bounds=(-4, 4)
    gamma0, gamma_max, system_type, T = 72.0, 3.5e-6, "nonionic", 298.15
    concentrations = [0.05, 0.1, 0.2, 0.4, 0.8, 1.2, 1.6, 2.0]
    tensions = [frumkin_surface_tension(c, gamma0, gamma_max, true_K, true_a, system_type, T) for c in concentrations]

    result = frumkin_fit_K_and_a(concentrations, tensions, gamma0, gamma_max, system_type, T)
    assert result.a_pinned_at_bound
    assert any("search bound" in w for w in result.warnings)


def test_frumkin_fit_rejects_bad_inputs():
    concentrations = [0.1, 0.2, 0.4, 0.8]
    tensions = [70.0, 65.0, 60.0, 55.0]
    with pytest.raises(ValueError):
        frumkin_fit_K_and_a([0.1, 0.2], tensions, 72.0, 4.0e-6, "nonionic")  # mismatched lengths
    with pytest.raises(ValueError):
        frumkin_fit_K_and_a([0.1, 0.2, 0.3], [70.0, 65.0, 60.0], 72.0, 4.0e-6, "nonionic")  # too few points (need >=4)
    with pytest.raises(ValueError):
        frumkin_fit_K_and_a([0.0, 0.2, 0.4, 0.8], tensions, 72.0, 4.0e-6, "nonionic")  # non-positive concentration
    with pytest.raises(ValueError):
        frumkin_fit_K_and_a(concentrations, tensions, 72.0, 0.0, "nonionic")  # non-positive gamma_max
    with pytest.raises(ValueError):
        frumkin_fit_K_and_a(concentrations, tensions, 72.0, 4.0e-6, "not_a_real_type")  # bad system_type
    with pytest.raises(ValueError):
        frumkin_fit_K_and_a(concentrations, tensions, 72.0, 4.0e-6, "nonionic", a_bounds=(2.0, -2.0))  # inverted bounds
