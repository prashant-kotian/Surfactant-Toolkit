"""Tests for surfactantkit.adsorption.szyszkowski_fit_K -- closes a real,
confirmed gap: szyszkowski_surface_tension's own docstring referenced this
function by name ("K... fit from real data, e.g. via szyszkowski_fit_K")
but it did not exist anywhere in the codebase (confirmed via grep before
building it). Real, standard two-step workflow: Gamma_max comes from the
pre-CMC log-slope independently (gibbs_gamma_max); this fits the single
remaining free parameter K via genuine nonlinear least squares against
the Szyszkowski equation itself. Validated via mathematically-guaranteed
round trips (construct a curve from a chosen true K via
szyszkowski_surface_tension itself, fit it back) -- no external raw
numeric literature table was sourced for this specific fit during this
pass, same disclosed pattern as other regression tools built this
session with no clean external numeric example.
"""

import pytest

from surfactantkit.adsorption import szyszkowski_fit_K, szyszkowski_surface_tension


def test_szyszkowski_fit_k_recovers_true_k_round_trip():
    true_K = 0.35  # 1/mM
    gamma0, gamma_max, system_type, T = 72.0, 4.0e-6, "nonionic", 298.15
    concentrations = [0.5, 1.0, 2.0, 4.0, 8.0, 14.0]
    tensions = [szyszkowski_surface_tension(c, gamma0, gamma_max, true_K, system_type, T) for c in concentrations]

    result = szyszkowski_fit_K(concentrations, tensions, gamma0, gamma_max, system_type, T)
    assert result.K == pytest.approx(true_K, rel=1e-4)
    assert result.r_squared == pytest.approx(1.0, abs=1e-6)
    assert result.n_points == 6


def test_szyszkowski_fit_k_recovers_true_k_for_ionic_no_added_salt():
    true_K = 1.8
    gamma0, gamma_max, system_type, T = 72.0, 3.2e-6, "ionic_no_added_salt", 298.15
    concentrations = [0.2, 0.5, 1.0, 2.0, 3.5, 5.0]
    tensions = [szyszkowski_surface_tension(c, gamma0, gamma_max, true_K, system_type, T) for c in concentrations]

    result = szyszkowski_fit_K(concentrations, tensions, gamma0, gamma_max, system_type, T)
    assert result.K == pytest.approx(true_K, rel=1e-4)
    assert result.r_squared == pytest.approx(1.0, abs=1e-6)


def test_szyszkowski_fit_k_fitted_tensions_match_predict_function():
    true_K = 0.6
    gamma0, gamma_max, system_type, T = 70.0, 4.5e-6, "nonionic", 298.15
    concentrations = [1.0, 3.0, 6.0, 10.0]
    tensions = [szyszkowski_surface_tension(c, gamma0, gamma_max, true_K, system_type, T) for c in concentrations]

    result = szyszkowski_fit_K(concentrations, tensions, gamma0, gamma_max, system_type, T)
    for c, fitted in zip(concentrations, result.fitted_surface_tensions_mN_per_m):
        expected = szyszkowski_surface_tension(c, gamma0, gamma_max, result.K, system_type, T)
        assert fitted == pytest.approx(expected)


def test_szyszkowski_fit_k_poor_fit_gives_lower_r_squared():
    """A curve that does NOT follow the Szyszkowski shape at all (e.g.
    surface tension increasing with concentration) should fit much worse
    than a genuine Szyszkowski curve -- r_squared is a real diagnostic,
    not always ~1."""
    gamma0, gamma_max, system_type, T = 72.0, 4.0e-6, "nonionic", 298.15
    concentrations = [0.5, 1.0, 2.0, 4.0, 8.0, 14.0]
    bad_tensions = [72.0 + 0.5 * c for c in concentrations]  # rising, unphysical for this model

    result = szyszkowski_fit_K(concentrations, bad_tensions, gamma0, gamma_max, system_type, T)
    assert result.r_squared < 0.5


def test_szyszkowski_fit_k_rejects_bad_inputs():
    concentrations = [1.0, 2.0, 3.0]
    tensions = [70.0, 65.0, 60.0]
    with pytest.raises(ValueError):
        szyszkowski_fit_K([1.0, 2.0], tensions, 72.0, 4.0e-6, "nonionic")  # mismatched lengths
    with pytest.raises(ValueError):
        szyszkowski_fit_K([1.0, 2.0], [70.0, 65.0], 72.0, 4.0e-6, "nonionic")  # too few points
    with pytest.raises(ValueError):
        szyszkowski_fit_K([0.0, 2.0, 3.0], tensions, 72.0, 4.0e-6, "nonionic")  # non-positive concentration
    with pytest.raises(ValueError):
        szyszkowski_fit_K(concentrations, tensions, 72.0, 0.0, "nonionic")  # non-positive gamma_max
    with pytest.raises(ValueError):
        szyszkowski_fit_K(concentrations, tensions, 72.0, 4.0e-6, "not_a_real_type")  # bad system_type
