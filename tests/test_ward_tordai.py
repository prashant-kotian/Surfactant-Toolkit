"""Tests for surfactantkit.adsorption.ward_tordai_adsorption (worker P1).

Golds: (a) fixed surface concentration vs the exact solution
Gamma = 2 (c0 - cs) sqrt(D t / pi) for all t, with c0 - cs converted to
mol/m^3 (x1000 from mol/L) so Gamma is mol/m^2; (b) Langmuir transient
obeys the short-time clean-surface law 2 c0 sqrt(D t / pi) (same unit
conversion -- regression-locked, the missing factor 1000 starved the
Langmuir supply by 1000x); (c) long-time limit approaches the Langmuir
equilibrium monotonically from below.
"""
import math

import pytest

from surfactantkit.adsorption import WardTordaiResult, ward_tordai_adsorption

C0, D = 1e-3, 1e-9


def test_fixed_sink_matches_exact_for_all_times():
    r = ward_tordai_adsorption(C0, D, t_end_s=1.0, mode="fixed_cs",
                               cs_fixed_mol_per_L=0.0)
    # mol/L -> mol/m^3 (x1000): Gamma is mol/m^2 (SI).
    exact = 2 * (C0 * 1000.0) * math.sqrt(D * 1.0 / math.pi)
    assert r.gamma_final_mol_per_m2 == pytest.approx(exact, rel=1e-9)
    # the whole series must sit on the exact curve
    for t, g in zip(r.t_s, r.gamma_mol_per_m2):
        assert g == pytest.approx(2 * (C0 * 1000.0) * math.sqrt(D * t / math.pi),
                                  rel=1e-9)


def test_fixed_half_concentration_exact():
    r = ward_tordai_adsorption(C0, D, t_end_s=0.25, mode="fixed_cs",
                               cs_fixed_mol_per_L=C0 / 2)
    # mol/L -> mol/m^3 (x1000): Gamma is mol/m^2 (SI).
    exact = 2 * (C0 / 2 * 1000.0) * math.sqrt(D * 0.25 / math.pi)
    assert r.gamma_final_mol_per_m2 == pytest.approx(exact, rel=1e-9)


def test_langmuir_short_time_clean_surface_law():
    r = ward_tordai_adsorption(c0_mol_per_L=C0, D_m2_per_s=D, t_end_s=1.0,
                               n_times=40)
    assert 0.97 <= r.early_time_ratio <= 1.03


def test_langmuir_approaches_equilibrium_monotonically():
    r = ward_tordai_adsorption(c0_mol_per_L=C0, t_end_s=200_000.0,
                               n_times=80, gamma_inf_mol_per_m2=3e-6,
                               K_d_mol_per_L=1e-3)
    ge = r.gamma_equilibrium_mol_per_m2
    assert ge is not None
    frac = r.gamma_final_mol_per_m2 / ge
    assert 0.94 < frac < 1.0            # slow diffusive tail from below
    assert r.cs_final_mol_per_L < C0     # surface still depleting bulk
    assert r.converged


def test_langmuir_si_supply_not_starved() -> None:
    """P1 regression: without the x1000 mol/L -> mol/m^3 conversion the
    diffusion supply feeding the Langmuir boundary is 1000x too weak, so
    the subsurface equilibrates far below c0 (cs* ~ 0.94 c0 instead of
    ~ c0 for these parameters) and the transient is quantitatively wrong."""
    C0, D, KD, GINF = 1e-3, 1e-9, 1e-3, 3e-6
    r = ward_tordai_adsorption(
        c0_mol_per_L=C0, D_m2_per_s=D, t_end_s=200_000.0, n_times=80,
        gamma_inf_mol_per_m2=GINF, K_d_mol_per_L=KD,
    )
    # SI supply: subsurface essentially recovered (starved bug: 0.940 c0)
    assert r.cs_final_mol_per_L > 0.99 * C0
    # stored gamma satisfies the Langmuir isotherm at the final c_s
    iso = GINF * r.cs_final_mol_per_L / (KD + r.cs_final_mol_per_L)
    assert r.gamma_final_mol_per_m2 == pytest.approx(iso, rel=1e-6)
    # approaches equilibrium from below, and the clean-surface law holds
    # at the first node (quadrature + unit-conversion diagnostic)
    ge = r.gamma_equilibrium_mol_per_m2
    assert ge is not None
    frac = r.gamma_final_mol_per_m2 / ge
    assert 0.99 < frac < 1.0
    assert 0.97 <= r.early_time_ratio <= 1.03
    assert r.converged


def test_series_invariants():
    r = ward_tordai_adsorption(c0_mol_per_L=C0, t_end_s=500.0, n_times=50)
    g, cs = r.gamma_mol_per_m2, r.cs_mol_per_L
    assert all(g[i] <= g[i + 1] * (1 + 1e-12) for i in range(len(g) - 1))
    assert all(cs[i] <= cs[i + 1] * (1 + 1e-12) for i in range(len(cs) - 1))
    assert all(c <= C0 * (1 + 1e-12) for c in cs)
    assert r.gamma_equilibrium_mol_per_m2 is not None
    assert g[-1] < r.gamma_equilibrium_mol_per_m2
    assert len(r.t_s) == 50


def test_resolution_convergence():
    ra = ward_tordai_adsorption(c0_mol_per_L=C0, t_end_s=5.0, n_times=40,
                                n_sqrt=300)
    rb = ward_tordai_adsorption(c0_mol_per_L=C0, t_end_s=5.0, n_times=40,
                                n_sqrt=900)
    rel = abs(ra.gamma_final_mol_per_m2 - rb.gamma_final_mol_per_m2) \
        / rb.gamma_final_mol_per_m2
    assert rel < 1e-4


def test_validation_errors():
    with pytest.raises(ValueError, match="c0_mol_per_L"):
        ward_tordai_adsorption(0.0)
    with pytest.raises(ValueError, match="D_m2_per_s"):
        ward_tordai_adsorption(C0, 0.0)
    with pytest.raises(ValueError, match="t_end_s"):
        ward_tordai_adsorption(C0, D, 0.0)
    with pytest.raises(ValueError, match="n_times"):
        ward_tordai_adsorption(C0, D, 1.0, n_times=3)
    with pytest.raises(ValueError, match="n_sqrt"):
        ward_tordai_adsorption(C0, D, 1.0, n_sqrt=5)
    with pytest.raises(ValueError, match="mode"):
        ward_tordai_adsorption(C0, D, 1.0, mode="frumkin")
    with pytest.raises(ValueError, match="required"):
        ward_tordai_adsorption(C0, D, 1.0, mode="fixed_cs")
    with pytest.raises(ValueError, match=r"\[0, c0"):
        ward_tordai_adsorption(C0, D, 1.0, mode="fixed_cs",
                               cs_fixed_mol_per_L=2 * C0)
    with pytest.raises(ValueError, match="gamma_inf"):
        ward_tordai_adsorption(C0, D, 1.0, gamma_inf_mol_per_m2=0.0)


def test_result_metadata():
    r = ward_tordai_adsorption(c0_mol_per_L=C0, t_end_s=10.0)
    assert isinstance(r, WardTordaiResult)
    assert "Ward-Tordai" in r.method
    assert r.params["mode"] == "langmuir"
    assert "convective" in r.caveat


def test_langmuir_dilute_limit_matches_exact_henry_solution():
    """Kc << 1 makes the Langmuir boundary linear (Gamma = K_H c_s, K_H = Gamma_inf/K_d),
    for which diffusion-controlled uptake is exact:
    Gamma(t) = K_H c0 [1 - exp(t/tau) erfc(sqrt(t/tau))], tau = K_H^2 / D."""
    from scipy.special import erfcx
    c0, d, ginf, kd = 1e-8, 1e-9, 3e-6, 1e-3
    kh = ginf / kd / 1000.0            # m (c in mol/m^3)
    tau = kh * kh / d
    r = ward_tordai_adsorption(c0, d, t_end_s=5.0, n_times=16, mode="langmuir",
                               gamma_inf_mol_per_m2=ginf, K_d_mol_per_L=kd)
    for t, g in zip(r.t_s, r.gamma_mol_per_m2):
        exact = kh * c0 * 1000.0 * (1.0 - erfcx(math.sqrt(t / tau)))
        assert g == pytest.approx(exact, rel=2e-4)
