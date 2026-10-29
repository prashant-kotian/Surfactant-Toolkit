"""Tests for surfactantkit.dynamic_tension (worker O3).

Gold (fixed seed=11, sigma=0.05, c0=5e-3 M, t=0.1..50 s, 20 pts):
- forward gamma(t) strictly decreasing toward gamma_equilibrium;
- fixed-D joint fit recovers the IDENTIFIABLE lumped quantities
  (Gamma_eq ratio 0.9979, adsorption rate koff*(1+Kc) ratio 1.0021)
  with r2 0.999653, K-Gamma_max correlation -1.0000 and the
  identifiability_warning raised (the ridge is real and disclosed);
- 4-param fit r2 0.999666, lumped ratios 0.9985/0.9885, converged.
K and Gamma_max individually are NOT asserted (single-transient ridge --
the warning says quote the lumped product instead); CI containment of
truth was deliberately dropped for the same reason.
"""
import itertools

import pytest

from surfactantkit.dynamic_tension import (
    DynamicSurfaceTensionFitResult,
    DynamicSurfaceTensionResult,
    dynamic_surface_tension_forward,
    dynamic_surface_tension_ode,
)

TRUTH = {"D_m2_per_s": 6e-11, "K_per_M": 500.0,
             "gamma_max_mol_per_m2": 4e-6, "k_off_per_s": 0.5}
C0 = 5e-3
T = [float(x) for x in
     __import__("numpy").logspace(-1, 1.7, 20)]
NOISE = 0.05


def _synthetic(seed=11):
    import numpy as np
    fwd = dynamic_surface_tension_forward(T, C0, **TRUTH)
    rng = np.random.default_rng(seed)
    y = [float(min(71.9, max(0.1, g + rng.normal(0, NOISE))))
         for g in fwd.gamma_mN_per_m]
    return fwd, y


def _lumped_ratios(K, k_off, Gm):
    eq = Gm * K * C0 / (1.0 + K * C0)
    eq_t = TRUTH["gamma_max_mol_per_m2"] * TRUTH["K_per_M"] * C0 \
        / (1.0 + TRUTH["K_per_M"] * C0)
    rate = k_off * (1.0 + K * C0)
    rate_t = TRUTH["k_off_per_s"] * (1.0 + TRUTH["K_per_M"] * C0)
    return eq / eq_t, rate / rate_t


def test_forward_monotone_to_equilibrium():
    fwd, _ = _synthetic()
    assert isinstance(fwd, DynamicSurfaceTensionResult)
    g = fwd.gamma_mN_per_m
    assert all(b < a for a, b in itertools.pairwise(g))          # tension drops
    assert 0 < g[-1] < g[0] < 72.0
    assert g[-1] == pytest.approx(fwd.gamma_equilibrium_mN_per_m,
                                  abs=0.05)
    assert fwd.surface_excess_mol_per_m2[-1] > 0
    assert fwd.converged
    for k in ("converged", "method", "params", "caveat"):
        assert k in fwd.as_dict()


def test_fixed_D_fit_recovers_lumped_quantities():
    _, y = _synthetic()
    r = dynamic_surface_tension_ode(T, y, C0, sigma_mN_per_m=NOISE,
                                    fit_D=False)
    assert isinstance(r, DynamicSurfaceTensionFitResult)
    assert r.converged
    assert r.r_squared > 0.999
    a_ratio, k_ratio = _lumped_ratios(r.K_per_M, r.params["k_off_per_s_fit"],
                                      r.gamma_max_mol_per_m2)
    assert a_ratio == pytest.approx(1.0, abs=0.03)   # gold 0.9979
    assert k_ratio == pytest.approx(1.0, abs=0.10)   # gold 1.0021
    p = r.params
    assert p["log_K_vs_log_Gmax_correlation"] < -0.99   # gold -1.0000
    assert p["identifiability_warning"], "ridge must be disclosed"
    assert p["fit_D"] is False
    assert p["k_off_per_s_fit"] > 0
    for k in ("converged", "n_iterations", "method", "params", "caveat"):
        assert k in r.as_dict()


def test_four_param_fit():
    _, y = _synthetic()
    r = dynamic_surface_tension_ode(T, y, C0, sigma_mN_per_m=NOISE,
                                    fit_D=True)
    assert r.converged and r.r_squared > 0.999
    a_ratio, k_ratio = _lumped_ratios(r.K_per_M, r.params["k_off_per_s_fit"],
                                      r.gamma_max_mol_per_m2)
    assert a_ratio == pytest.approx(1.0, abs=0.05)   # gold 0.9985
    assert k_ratio == pytest.approx(1.0, abs=0.10)   # gold 0.9885
    # D is nearly unidentified from this transient (caveat says so);
    # only reject nonsense outside a 2-decade band of the truth
    assert 6e-13 < r.D_m2_per_s < 6e-9
    assert r.params["identifiability_warning"] is not None


def test_deterministic_same_input():
    _, y = _synthetic()
    a = dynamic_surface_tension_ode(T, y, C0, sigma_mN_per_m=NOISE,
                                    fit_D=False)
    b = dynamic_surface_tension_ode(T, y, C0, sigma_mN_per_m=NOISE,
                                    fit_D=False)
    assert a.K_per_M == b.K_per_M
    assert a.r_squared == b.r_squared
    assert a.as_dict() == b.as_dict()


def test_forward_errors():
    with pytest.raises(ValueError, match="positive"):
        dynamic_surface_tension_forward([1.0, -0.5], C0, **TRUTH)
    with pytest.raises(ValueError, match="sorted"):
        dynamic_surface_tension_forward([2.0, 1.0], C0, **TRUTH)
    with pytest.raises(ValueError, match="positive"):
        dynamic_surface_tension_forward(T, C0, 6e-11, -1.0, 4e-6)
    with pytest.raises(ValueError, match="Gibbs"):
        dynamic_surface_tension_forward(T, C0, 6e-11, 500.0, 4e-6,
                                        n_factor=1.5)
    with pytest.raises(ValueError, match="n_cells"):
        dynamic_surface_tension_forward(T, C0, 6e-11, 500.0, 4e-6,
                                        n_cells=10)
    with pytest.raises(ValueError, match="k_off"):
        dynamic_surface_tension_forward(T, C0, 6e-11, 500.0, 4e-6,
                                        k_off_per_s=0.0)


def test_fit_errors():
    _, y = _synthetic()
    with pytest.raises(ValueError, match="same length"):
        dynamic_surface_tension_ode(T, y[:-1], C0)
    with pytest.raises(ValueError, match="at least 6"):
        dynamic_surface_tension_ode(T[:5], y[:5], C0)
    with pytest.raises(ValueError, match=r"\(0, gamma0\)"):
        dynamic_surface_tension_ode(T, [0.0] * len(T), C0)
    with pytest.raises(ValueError, match="positive"):
        dynamic_surface_tension_ode(T, y, 0.0)
    with pytest.raises(ValueError, match="positive"):
        dynamic_surface_tension_ode(T, y, C0, sigma_mN_per_m=-0.1)
    with pytest.raises(ValueError, match="p0"):
        dynamic_surface_tension_ode(T, y, C0, p0=(1.0, 2.0))
    with pytest.raises(ValueError, match="positive"):
        dynamic_surface_tension_ode(T, y, C0, p0=(-1.0, 500.0, 4e-6, 0.5))


# ---------------------------------------------------------------------------
# equation of state, stability and cross-solver agreement
# ---------------------------------------------------------------------------

def test_long_time_limit_is_the_szyszkowski_equation():
    """gamma(t -> long) and gamma_equilibrium must equal the Langmuir isotherm's own
    equation of state (adsorption.szyszkowski_surface_tension), not gamma0 - nRT Gamma."""
    from surfactantkit.adsorption import szyszkowski_surface_tension
    for c0, K, gm in ((1e-4, 1e3, 3e-6), (1e-2, 1e3, 3e-6), (0.02, 1500.0, 4e-6)):
        r = dynamic_surface_tension_forward([200.0], c0, 1e-9, K, gm, k_off_per_s=200.0)
        ref = szyszkowski_surface_tension(c0, 72.0, gm, K, "nonionic")
        assert r.gamma_equilibrium_mN_per_m == pytest.approx(ref, rel=1e-12)
        assert r.gamma_mN_per_m[-1] == pytest.approx(ref, abs=0.05)


def test_uptake_is_monotone_with_fast_kinetics_at_low_concentration():
    """Stiff surface kinetics (k_off*dt >> 1) at low c0: the coupled half-cell/surface
    step must not overshoot or empty the surface between output times."""
    ts = [5.0 * (k / 16) ** 2 for k in range(1, 17)]
    for c0 in (1e-4, 1e-3, 1e-2):
        g = dynamic_surface_tension_forward(ts, c0, 1e-9, 1e3, 3e-6,
                                            k_off_per_s=200.0).surface_excess_mol_per_m2
        assert all(b >= a - 1e-18 for a, b in itertools.pairwise(g)), (c0, g)
        assert min(g) > 0.0


def test_agrees_with_ward_tordai_in_the_diffusion_limited_regime():
    """Fast kinetics -> diffusion control: the CN solver and the independent Ward-Tordai
    quadrature (adsorption.py) must give the same surface excess within 5 %."""
    from surfactantkit.adsorption import ward_tordai_adsorption
    for c0, tol in ((1e-4, 0.05), (1e-3, 0.03), (1e-2, 0.005)):
        wt = ward_tordai_adsorption(c0, 1e-9, t_end_s=5.0, n_times=16, mode="langmuir",
                                    gamma_inf_mol_per_m2=3e-6, K_d_mol_per_L=1e-3)
        pde = dynamic_surface_tension_forward(wt.t_s, c0, 1e-9, 1e3, 3e-6, k_off_per_s=1e4)
        for a, b in zip(pde.surface_excess_mol_per_m2, wt.gamma_mol_per_m2):
            assert a == pytest.approx(b, rel=tol)


# ---- numerical accuracy of the forward solver (stretched grid, coupled implicit step)

_ACC_CASES = [
    ([0.1, 0.5, 1.0, 3.0, 10.0, 30.0, 100.0], 0.02, 1e-9, 1500.0, 4e-6, 0.5),
    ([0.1, 1.0, 10.0, 100.0, 200.0], 0.001, 1e-9, 1e3, 3e-6, 200.0),
    ([0.01, 0.1, 1.0, 10.0, 100.0, 1000.0], 5e-4, 6e-11, 500.0, 4e-6, 0.5),
]


@pytest.mark.parametrize("ts, c0, D, K, gm, koff", _ACC_CASES)
def test_forward_converged_at_default_settings(ts, c0, D, K, gm, koff):
    """Default grid and step agree with a 4x finer grid and 4x smaller steps."""
    base = dynamic_surface_tension_forward(ts, c0, D, K, gm, k_off_per_s=koff)
    fine = dynamic_surface_tension_forward(ts, c0, D, K, gm, k_off_per_s=koff,
                                           n_cells=640, rel_step=0.0025)
    for a, b in zip(base.gamma_mN_per_m, fine.gamma_mN_per_m):
        assert a == pytest.approx(b, abs=2e-3)


@pytest.mark.parametrize("ts, c0, D, K, gm, koff", _ACC_CASES)
def test_forward_independent_of_requested_times(ts, c0, D, K, gm, koff):
    """gamma at the first output time does not depend on which later times
    are also requested (the grid used to be sized from the last one only)."""
    long = dynamic_surface_tension_forward(ts, c0, D, K, gm, k_off_per_s=koff)
    short = dynamic_surface_tension_forward(ts[:1], c0, D, K, gm, k_off_per_s=koff)
    assert long.gamma_mN_per_m[0] == pytest.approx(short.gamma_mN_per_m[0], abs=2e-3)


def test_forward_kinetic_limit_matches_analytic_langmuir_uptake():
    """Fast diffusion: the subsurface stays at c0, so the surface follows
    dG/dt = k_a c0 (Gmax - G) - k_off G exactly:
    G(t) = G_eq (1 - exp(-(k_a c0 + k_off) t)), G_eq = Gmax K c0 / (1 + K c0)."""
    import math
    c0, K, gm, koff = 1e-3, 800.0, 3e-6, 0.4
    ts = [0.05, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]
    r = dynamic_surface_tension_forward(ts, c0, 1e-3, K, gm, k_off_per_s=koff)
    rate = K * koff * c0 + koff
    g_eq = gm * K * c0 / (1.0 + K * c0)
    for t, g in zip(ts, r.surface_excess_mol_per_m2):
        assert g == pytest.approx(g_eq * (1.0 - math.exp(-rate * t)), rel=1e-4)
