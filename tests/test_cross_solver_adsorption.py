"""Cross-solver validation: Ward-Tordai Volterra integral vs the
finite-volume PDE with Langmuir kinetics (dynamic_tension).

Two fully independent implementations of the same physics must agree:

* ``ward_tordai_adsorption`` -- product integration of the Ward-Tordai
  integral (c_s piecewise linear, exact kernel weights) against the
  *equilibrium* Langmuir boundary condition (mol/L grid, explicit
  L->m^3 x1000 conversion; P1 unit-bug regression is locked in
  tests/test_ward_tordai.py);
* ``dynamic_surface_tension_forward`` -- implicit 1-D diffusion on a
  stretched mol/m^3 grid coupled to *kinetic* Langmuir uptake; run with
  fast kinetics (k_off = 200 1/s) it relaxes onto the equilibrium-BC
  (Ward-Tordai) limit, and with k_off = 1e6 1/s it must reproduce it.

Any concentration-unit mismatch (the P1 bug class: mol/L vs mol/m^3)
would show up here as an orders-of-magnitude divergence, not a few
percent, so this test doubles as an independent unit conversion check
of both solvers.
"""
import math

import pytest

from surfactantkit.adsorption import ward_tordai_adsorption
from surfactantkit.dynamic_tension import dynamic_surface_tension_forward

C0, D, GINF, KD = 1e-3, 1e-9, 3e-6, 1e-3  # K = 1/KD = 1000 1/M
T_END, N_PTS = 5.0, 16
T_S = [T_END * (k / N_PTS) ** 2 for k in range(1, N_PTS + 1)]


def _run_both():
    wt = ward_tordai_adsorption(
        c0_mol_per_L=C0, D_m2_per_s=D, t_end_s=T_END, n_times=N_PTS,
        mode="langmuir", gamma_inf_mol_per_m2=GINF, K_d_mol_per_L=KD,
    )
    pde = dynamic_surface_tension_forward(
        T_S, C0, D, K_per_M=1.0 / KD, gamma_max_mol_per_m2=GINF,
        n_factor=1.0, k_off_per_s=200.0, n_cells=60, max_substep_s=2e-4,
    )
    return wt, pde


def test_surface_excess_agrees_across_solvers() -> None:
    wt, pde = _run_both()
    assert wt.converged and pde.converged
    g_wt = wt.gamma_mol_per_m2
    g_pde = pde.surface_excess_mol_per_m2
    assert len(g_wt) == len(g_pde) == N_PTS
    for a, b in zip(g_wt, g_pde):
        # worst observed ~2.3% (early boundary-layer discretisation);
        # a mol/L<->mol/m^3 slip would be ~1000x, far outside this bound
        assert abs(b - a) / a <= 0.04
    # tight agreement at long times (quasi-equilibrium tail)
    assert abs(g_pde[-1] - g_wt[-1]) / g_wt[-1] <= 0.01


def test_tension_curves_agree_across_solvers() -> None:
    wt, pde = _run_both()
    # map the Ward-Tordai excess to gamma(t) with the Szyszkowski equation of state of the
    # same Langmuir isotherm, as dynamic_surface_tension_forward does
    gamma0, R, T = 72.0, 8.314462618, 298.15
    gy_wt = [gamma0 + R * T * GINF * math.log(1.0 - g / GINF) * 1e3 for g in wt.gamma_mol_per_m2]
    gy_pde = pde.gamma_mN_per_m
    worst = max(abs(a - b) for a, b in zip(gy_wt, gy_pde))
    assert worst <= 0.15  # mN/m; observed ~0.08 on a ~68 mN/m curve
    assert math.isfinite(pde.gamma_equilibrium_mN_per_m)


def test_subsurface_recovers_on_both_solvers() -> None:
    wt, pde = _run_both()
    # subsurface monotonically recovers toward bulk in both solvers
    cs_wt, cs_pde = wt.cs_mol_per_L, pde.subsurface_M
    assert all(cs_wt[i] <= cs_wt[i + 1] * (1 + 1e-12)
               for i in range(len(cs_wt) - 1))
    assert all(cs_pde[i] <= cs_pde[i + 1] * (1 + 1e-9)
               for i in range(len(cs_pde) - 1))
    # final subsurface within 3% of each other (observed ~0.5%)
    assert abs(cs_pde[-1] - cs_wt[-1]) / cs_wt[-1] <= 0.03
    # both end near (but below) bulk -- supply-limited tail
    assert cs_wt[-1] < C0 and cs_pde[-1] < C0
    assert cs_wt[-1] > 0.95 * C0 and cs_pde[-1] > 0.95 * C0


def test_diffusion_controlled_limit_agrees_tightly() -> None:
    """With quasi-instant kinetics (k_off = 1e6 1/s) the PDE is the
    local-equilibrium problem the Ward-Tordai integral solves: the two
    independent methods must agree to ~1e-4, from the dilute (Henry) to the
    near-saturated case (c0 = 0.02 M, K_d = 1e-3 M: Gamma(1 s) = 2.8568
    umol/m^2, consistent with the long-time asymptote
    Gamma_eq - (dGamma/dc_s) Gamma / sqrt(pi D t))."""
    for c0, t_end in ((1e-5, 5.0), (1e-3, 5.0), (0.02, 1.0)):
        wt = ward_tordai_adsorption(c0, D, t_end_s=t_end, n_times=12, mode="langmuir",
                                    gamma_inf_mol_per_m2=GINF, K_d_mol_per_L=KD)
        pde = dynamic_surface_tension_forward(list(wt.t_s), c0, D, 1.0 / KD, GINF,
                                              k_off_per_s=1e6)
        for a, b in zip(wt.gamma_mol_per_m2[2:], pde.surface_excess_mol_per_m2[2:]):
            assert a == pytest.approx(b, rel=2e-4), (c0, a, b)
    assert round(wt.gamma_final_mol_per_m2 * 1e6, 4) == 2.8568
