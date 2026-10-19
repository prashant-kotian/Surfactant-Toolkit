"""Tests for surfactantkit.electrostatics.planar_pb_profile (worker P2).

Golds: (a) exact 1:1 Gouy-Chapman analytic potential profile;
(b) Gouy-Chapman closed-form surface charge; (c) Debye-Huckel linear limit
at small potential; (d) sinh scaling of sigma with psi0; (e) consistency
with the library's own debye_length().
"""
import math

import pytest

from surfactantkit.electrostatics import (
    WATER_REL_PERMITTIVITY_25C,
    debye_length,
    planar_pb_profile,
)


def test_matches_exact_analytic_profile():
    r = planar_pb_profile(50.0, 0.01, n_points=101, x_max_nm=3 * 3.04)
    assert r.converged
    assert r.psi_mV[0] == pytest.approx(50.0, abs=1e-6)
    # exact-solution agreement over the whole profile
    assert r.max_abs_err_vs_analytic_mV < 1e-6


def test_sigma_matches_gouy_chapman_closed_form():
    r = planar_pb_profile(50.0, 0.01)
    assert r.sigma_rel_error < 1e-6
    # independent numeric magnitude sanity: ~1e-2 C/m^2 for 50 mV @ 10 mM
    assert 1e-3 < abs(r.sigma_surface_C_m2) < 0.2


def test_debye_huckel_linear_limit():
    # at 0.1 mV the linear theory sigma = eps*psi0/lambda_D must hold
    r = planar_pb_profile(0.1, 0.01)
    ld_m = r.debye_length_nm * 1e-9
    sigma_lin = (WATER_REL_PERMITTIVITY_25C * 8.8541878128e-12) * 1e-4 / ld_m
    assert r.sigma_surface_C_m2 == pytest.approx(sigma_lin, rel=1e-3)


def test_sigma_sinh_scaling_with_potential():
    # sigma(100 mV)/sigma(50 mV) == sinh(u100/2)/sinh(u50/2)
    r50 = planar_pb_profile(50.0, 0.01)
    r100 = planar_pb_profile(100.0, 0.01)
    kT_eV = 8.617333262e-5 * 298.15
    u50 = 0.050 / kT_eV
    u100 = 0.100 / kT_eV
    ratio_num = abs(r100.sigma_surface_C_m2 / r50.sigma_surface_C_m2)
    ratio_theory = math.sinh(u100 / 2) / math.sinh(u50 / 2)
    assert ratio_num == pytest.approx(ratio_theory, rel=1e-5)


def test_debye_length_consistent_with_library():
    r = planar_pb_profile(50.0, 0.01, temperature_K=310.15)
    assert r.debye_length_nm == pytest.approx(
        debye_length(0.01, 310.15), rel=1e-12)
    assert r.kappa_nm_inv == pytest.approx(1.0 / r.debye_length_nm, rel=1e-12)


def test_profile_decays_and_charge_sign_follows_potential():
    r = planar_pb_profile(50.0, 0.01, n_points=51, x_max_nm=16.0)
    psi = r.psi_mV
    assert all(psi[i] >= psi[i + 1] - 1e-9 for i in range(len(psi) - 1))
    # lambda_D ~ 3.04 nm at 10 mM; at 16 nm (~5.3 lambda_D) the double layer
    # must be strongly screened
    assert psi[-1] < psi[0] / 50.0
    assert r.charge_density_C_m3[0] < 0  # countercharge near positive wall
    neg = planar_pb_profile(-50.0, 0.01, n_points=51)
    assert neg.charge_density_C_m3[0] > 0


def test_zero_potential_is_trivial():
    r = planar_pb_profile(0.0, 0.01, n_points=21)
    assert all(abs(p) < 1e-12 for p in r.psi_mV)
    assert r.sigma_surface_C_m2 == pytest.approx(0.0, abs=1e-15)


def test_multivalent_valence_runs():
    # z:z = 2:2 -- analytic z:z profile must still hold (same scaled equation)
    r = planar_pb_profile(25.0, 0.01, valence=2, n_points=61, x_max_nm=5)
    assert r.converged
    assert r.max_abs_err_vs_analytic_mV < 1e-5


def test_validation_errors():
    with pytest.raises(ValueError, match="positive"):
        planar_pb_profile(50.0, 0.0)
    with pytest.raises(ValueError, match="temperature"):
        planar_pb_profile(50.0, 0.01, temperature_K=0)
    with pytest.raises(ValueError, match="n_points"):
        planar_pb_profile(50.0, 0.01, n_points=5)
    with pytest.raises(ValueError, match="finite"):
        planar_pb_profile(float("nan"), 0.01)
    with pytest.raises(ValueError, match="valence"):
        planar_pb_profile(50.0, 0.01, valence=0)
    with pytest.raises(ValueError, match="x_max_nm"):
        planar_pb_profile(50.0, 0.01, x_max_nm=-1.0)


def test_result_metadata():
    r = planar_pb_profile(50.0, 0.01)
    assert r.method and r.params["ionic_strength_M"] == 0.01
    assert "mean-field" in r.caveat
    assert r.n_nodes > 100
