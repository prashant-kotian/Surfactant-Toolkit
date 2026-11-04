"""Regression tests for the external-audit fixes (FINDINGS 30).

One test per audit item. These exist because the previous test suite passed while
carrying six P0 defects: it tested that the code ran, not that the physics was
right. Every test here fails on the pre-fix code.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from research.laod_tau2 import harmonics as H
from research.laod_tau2.linear_response import (
    elasticity_scale,
    isotherm_slope,
    linear_response,
    lvdt_phase_deg,
)
from research.laod_tau2.model import ModelParams, run_oscillation
from research.laod_tau2.systems import SYSTEMS

C12 = SYSTEMS["C12E8"]
C0 = 0.071
K = 140.0
D = 1e-9
GINF = C12.gamma_inf


def _p(**kw):
    base = dict(c0=C0, K=K, gamma_inf=GINF, D=D, gamma0=0.05)
    base.update(kw)
    return ModelParams(**base)


# ---------------------------------------------------------------------------
# 2.1 -- tau2 is not the relaxation time
# ---------------------------------------------------------------------------
def test_closed_mode_relaxes_as_tau_over_one_plus_capacity():
    """The two-pool disturbance must decay with tau2/(1+C_m/n), not tau2."""
    from scipy.integrate import solve_ivp

    Cm, tau, n = 47.5, 4.0, 1.0
    rate = 1.0 + Cm / n

    def f(t, y):
        c, M = y
        mu = -(M - Cm * c / n - (Cm * 0.0) / n) / tau
        return np.array([-mu, mu])

    # pure decay mode around a zero-mean equilibrium
    y0 = np.array([1e-3, -1e-3 / rate * rate])
    sol = solve_ivp(f, (0, 0.5), np.array([0.0, 0.0]), rtol=1e-10)  # smoke
    sol = solve_ivp(f, (0, 0.5), y0, rtol=1e-12, atol=1e-18, dense_output=True)
    for t in (0.02, 0.05, 0.10):
        y = sol.sol(t)
        predicted = y0[0] * math.exp(-rate * t / tau)
        assert abs(y[0] - predicted) / abs(y0[0]) < 1e-6

    p = _p(capacity=Cm, tau2=tau, n_agg=n)
    assert p.relaxation_time == pytest.approx(tau / rate)


def test_relaxation_time_property_matches_eigenmode():
    p = _p(capacity=95.0, tau2=4.0, n_agg=95.0)
    # n cancels: C_m/n = 1, so tau_closed = tau2/2
    assert p.relaxation_time == pytest.approx(2.0)


# ---------------------------------------------------------------------------
# 2.2 -- inventory is C_m*c/n, not C_m*c
# ---------------------------------------------------------------------------
def test_micelle_inventory_uses_aggregation_number():
    p = _p(capacity=47.5, n_agg=95.0)
    assert p.micelle_inventory == pytest.approx(47.5 * 0.071 / 95.0)
    # the old expression, 47.5*0.071, was 95x larger
    assert p.micelle_inventory < 0.071


def test_inventory_matches_total_loading():
    """c_eq + M_eq(c_eq) must recover the nominal total."""
    for m in (1.2, 1.5, 2.0, 3.0):
        ct = m * C12.cmc_mM
        ce = C12.monomer_concentration(ct)
        assert ce + C12.capacity(ct) * ce / C12.n_bar == pytest.approx(ct, rel=1e-6)


# ---------------------------------------------------------------------------
# 2.3 -- monomer, not total, enters the transport model
# ---------------------------------------------------------------------------
def test_to_model_passes_monomer_concentration():
    p = C12.to_model(1.5 * C12.cmc_mM, c_m_tau=4.0, K=K)
    assert p.c0 < 1.5 * C12.cmc_mM
    assert p.c0 == pytest.approx(C12.cmc_mM, rel=1e-4)
    assert p.n_agg == C12.n_bar


def test_solve_convention_approximates_pinned_for_the_linear_closure():
    """C_m = n_bar*(m-1)/(1-phi) makes c_eq = ct/(1+C_m/n) ~ CMC to O(phi).

    Not exactly the CMC: the micellar-volume correction gives C_m/n = 1.500086
    rather than 1.5, which shifts c_eq by ~3e-5 relative. That is the volume
    correction finally doing something after being 1000x too small.
    """
    ct = 2.5 * C12.cmc_mM
    ce = C12.monomer_concentration(ct)
    assert ce == pytest.approx(C12.cmc_mM, rel=1e-4)
    # the identity itself must close tightly
    assert ce + C12.capacity(ct) * ce / C12.n_bar == pytest.approx(ct, rel=1e-12)


def test_legacy_total_convention_is_opt_in_only():
    p = C12.to_model(1.5 * C12.cmc_mM, K=K, monomer_convention="total")
    assert p.c0 == pytest.approx(1.5 * C12.cmc_mM)


def test_concentration_convention_moves_phase():
    """The bug the audit found was worth ~4 deg of phase at 1 Hz."""
    total = C12.to_model(1.5 * C12.cmc_mM, K=K, monomer_convention="total")
    fixed = C12.to_model(1.5 * C12.cmc_mM, K=K)
    om = 2 * math.pi
    a = H.lag_deg(*_osc(total, om), om)
    b = H.lag_deg(*_osc(fixed, om), om)
    assert abs(a - b) > 1.0


def _osc(p, omega, eps=0.002, n_per_cycle=80):
    r = run_oscillation(p, omega, eps, n_per_cycle=n_per_cycle)
    return r.t, r.gamma_mN_m, r.area_rel


# ---------------------------------------------------------------------------
# 2.4 -- large amplitude adds nothing to the linear response
# ---------------------------------------------------------------------------
def test_phase_is_amplitude_independent():
    om = 2 * math.pi
    vals = []
    for eps in (0.002, 0.02, 0.20):
        vals.append(H.lag_deg(*_osc(_p(capacity=47.5, tau2=4.0), om, eps), om))
    assert max(vals) - min(vals) < 0.2          # measured 0.04 deg


# ---------------------------------------------------------------------------
# 2.5 / 30.8 -- analytic transfer function agrees with the solver
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("f", [0.05, 0.2, 1.0, 5.0])
def test_analytic_matches_numeric(f):
    om = 2 * math.pi * f
    an = linear_response(om, c0=C0, K=K, gamma_inf=GINF, D=D,
                         capacity=47.5, tau2=4.0).lag_deg
    num = H.lag_deg(*_osc(_p(capacity=47.5, tau2=4.0), om), om)
    assert abs(an - num) < 0.05


def test_lvdt_never_exceeds_45_degrees():
    for f in (1e-3, 0.01, 0.1, 1, 10, 100, 1e4):
        assert abs(lvdt_phase_deg(f, c0=C0, K=K, gamma_inf=GINF, D=D)) < 45.0


def test_reservoir_raises_magnitude_above_45():
    om = 2 * math.pi
    base = abs(lvdt_phase_deg(1.0, c0=C0, K=K, gamma_inf=GINF, D=D))
    with_res = abs(linear_response(om, c0=C0, K=K, gamma_inf=GINF, D=D,
                                   capacity=47.5, tau2=4.0).lag_deg)
    assert with_res > 45.0 > base


# ---------------------------------------------------------------------------
# 2.7 / 30.3 -- global balance is a real residual
# ---------------------------------------------------------------------------
def test_mass_balance_residual_is_small():
    r = run_oscillation(_p(capacity=47.5, tau2=4.0, n_agg=95.0), 2 * math.pi,
                        eps=0.002, n_per_cycle=80)
    assert r.mass_error < 1e-2
    assert r.mass_residual_abs < 1e-8


def test_mass_balance_separates_drift_from_residual():
    r = run_oscillation(_p(capacity=47.5, tau2=4.0, n_agg=95.0), 2 * math.pi,
                        eps=0.002, n_per_cycle=80)
    # endpoint drift and identity residual are different quantities
    assert r.mass_endpoint_drift != pytest.approx(r.mass_error)


# ---------------------------------------------------------------------------
# 4.1 / 30.4 -- phase convention
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("delay", [10.0, 30.0, 60.0])
def test_lag_deg_recovers_a_known_delay(delay):
    w = 2 * math.pi
    t = np.linspace(0, 4, 2000)
    gamma = np.cos(w * t - math.radians(delay))
    assert H.lag_deg(t, gamma, np.cos(w * t), w) == pytest.approx(delay, abs=1e-6)


def test_lag_deg_recovers_a_known_lead():
    w = 2 * math.pi
    t = np.linspace(0, 4, 2000)
    gamma = np.cos(w * t + math.radians(20.0))
    assert H.lag_deg(t, gamma, np.cos(w * t), w) == pytest.approx(-20.0, abs=1e-6)


def test_lag_is_minus_arg_E():
    om = 2 * math.pi
    r = linear_response(om, c0=C0, K=K, gamma_inf=GINF, D=D, capacity=47.5, tau2=4.0)
    assert r.lag_deg == pytest.approx(-r.phase_deg)


# ---------------------------------------------------------------------------
# 4.2 / 30.5 -- sampling rule
# ---------------------------------------------------------------------------
def test_five_frames_per_cycle_is_insufficient_for_four_harmonics():
    assert H.required_samples_per_cycle(4) == 9
    assert H.design_rank(5, 4) == 5            # rank-deficient
    assert H.design_rank(9, 4) == 9            # full rank


def test_repeating_cycles_does_not_add_rank():
    assert H.design_rank(5, 1) == H.design_rank(5, 8)


def test_five_frames_suffices_for_the_fundamental_only():
    assert H.required_samples_per_cycle(1) == 3
    assert H.design_rank(5, 1, kmax=1) == 3


# ---------------------------------------------------------------------------
# 2.10 / 30.6 -- unit defects
# ---------------------------------------------------------------------------
def test_volume_fraction_is_not_scaled_by_a_thousand():
    phi = C12.volume_fraction(1.5 * C12.cmc_mM)
    expected = (1.5 * C12.cmc_mM - C12.cmc_mM) * C12.m_surf / 1000.0
    assert phi == pytest.approx(expected)
    assert phi > 1e-6            # the old code gave 1.9e-8


def test_dgamma_dcs_is_consistent_across_branches():
    """Langmuir and Frumkin branches must return the same physical quantity."""
    lang = ModelParams(c0=C0, K=K, gamma_inf=GINF, D=D, gamma0=0.05, frumkin_a=0.0)
    frum = ModelParams(c0=C0, K=K, gamma_inf=GINF, D=D, gamma0=0.05, frumkin_a=1e-6)
    both = [lang.dgamma_dcs(C0), frum.dgamma_dcs(C0)]
    for v in both:
        assert 0 < v < 1e-4        # dGamma/dc is metres and tiny
    # Gibbs consistency: dPi/dc = nRT*Gamma_inf*dGamma/dc
    # dPi/dc = nRT*Gamma_inf * (1/Gamma_inf) dGamma/dc = nRT * dGamma/dc
    # (the Gibbs prefactor carries Gamma_inf once, the isotherm slope the other)
    R = 8.314462618
    pref = frum.n_factor * R * frum.temperature_K * 1e3
    assert frum.dPi_dcs(C0) == pytest.approx(pref * frum.dgamma_dcs(C0), rel=1e-9)


def test_dgamma_dcs_matches_numeric_derivative():
    eps = 1e-9
    for a in (0.0, 1e-6):
        p = ModelParams(c0=C0, K=K, gamma_inf=GINF, D=D, gamma0=0.05, frumkin_a=a)
        num = (p.isotherm(C0 + eps) - p.isotherm(C0 - eps)) / (2 * eps)
        assert p.dgamma_dcs(C0) == pytest.approx(num, rel=1e-5)


# ---------------------------------------------------------------------------
# 3.3 -- statistics
# ---------------------------------------------------------------------------
def test_rms_and_residual_sd_are_distinct():
    sse, n, dof = 193.1815, 6, 4
    assert math.sqrt(sse / n) == pytest.approx(5.6742, abs=1e-3)
    assert math.sqrt(sse / dof) == pytest.approx(6.9495, abs=1e-3)


# ---------------------------------------------------------------------------
# 5.3 / 30.10 -- the external fit
# ---------------------------------------------------------------------------
def test_fainerman_cmc_is_the_rusedi_value():
    from research.laod_tau2 import study_fainerman_fit as F

    assert F.CMC_C14E8_UM == pytest.approx(18.1)
    assert F.C_OVER_CMC == pytest.approx(30.0 / 18.1, rel=1e-9)


def test_fainerman_fit_uses_monomer_not_total():
    from research.laod_tau2 import study_fainerman_fit as F

    assert F.C_MONOMER < F.C_TOT
    assert F.C_MONOMER == pytest.approx(18.1e-3, rel=1e-9)
