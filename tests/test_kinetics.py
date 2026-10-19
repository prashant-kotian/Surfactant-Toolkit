"""Real numeric validation for kinetics.py (O1 micelle_relaxation_odes,
S2 kmc_micelle_kinetics, and the shared mass_action_equilibrium target).

This file did not exist before -- every other Phase-2 worker has its own
dedicated test file (test_dpd_cmc_s4.py, test_obrien_white_e1.py,
test_poisson_boltzmann.py, etc.) but kinetics.py had none, despite
ROADMAP.md/COMPUTATIONAL_WORKERS_RESEARCH.md claiming specific numeric
golds for it ("kinetic CMC 8.2355e-05", "eq ratios 0.92/0.95",
"tau_fit/tau_O1 = 0.40 band"). The module's own docstring states the
real identity this suite checks: O1, S2, and mass_action_equilibrium
describe the SAME physical system and must agree with each other.
"""
import math

import pytest

from surfactantkit.kinetics import (
    aniansson_wall_two_regime,
    kmc_micelle_kinetics,
    mass_action_equilibrium,
    micelle_relaxation_odes,
)

pytest.importorskip("numpy")
pytest.importorskip("scipy")


# ---------------------------------------------------------------------------
# mass_action_equilibrium -- the closed-form target everything else checks against
# ---------------------------------------------------------------------------

def test_mass_action_equilibrium_exact_mass_balance():
    x = mass_action_equilibrium(0.01, 5000.0, 40)
    total = sum((i + 1) * v for i, v in enumerate(x))
    assert abs(total - 0.01) / 0.01 < 1e-9
    assert all(v >= 0 for v in x)
    assert x[0] > 0.0


def test_mass_action_equilibrium_monotone_decreasing_in_n():
    x = mass_action_equilibrium(0.01, 5000.0, 40)
    for i in range(1, len(x)):
        assert x[i] <= x[i - 1] + 1e-300


@pytest.mark.parametrize("bad_kwargs", [
    {"c_total_M": 0.0, "K_per_M": 5000.0, "n_max": 40},
    {"c_total_M": 0.01, "K_per_M": 0.0, "n_max": 40},
    {"c_total_M": 0.01, "K_per_M": 5000.0, "n_max": 2},
])
def test_mass_action_equilibrium_error_paths(bad_kwargs):
    with pytest.raises(ValueError):
        mass_action_equilibrium(**bad_kwargs)


# ---------------------------------------------------------------------------
# O1 -- micelle_relaxation_odes
# ---------------------------------------------------------------------------

def test_o1_equilibrium_matches_closed_form_exactly():
    """The core identity the module docstring promises: O1's own
    independent bisection-derived equilibrium must equal
    mass_action_equilibrium's target for the SAME (c_total, K, n_max)."""
    c_total, K_final, n_max = 0.01, 5000.0, 40
    target = mass_action_equilibrium(c_total, K_final, n_max)
    r = micelle_relaxation_odes(c_total, 3000.0, K_final, 1e8, n_max=n_max, n_times=200)
    assert r.equilibrium_concentrations_M == pytest.approx(target, rel=1e-9)
    assert r.kinetic_cmc_M == pytest.approx(target[0], rel=1e-9)


def test_o1_mass_conservation_and_tail_fit_quality():
    r = micelle_relaxation_odes(0.01, 3000.0, 5000.0, 1e8, n_max=40, n_times=200)
    assert r.converged is True
    assert r.max_relative_mass_error < 1e-10
    assert r.tail_exp_fit_r_squared > 0.999
    assert r.tau_dominant_s > 0.0
    assert r.params["jacobian_fd_rel_error"] < 1e-4


def test_o1_relaxation_times_are_sorted_longest_first():
    r = micelle_relaxation_odes(0.01, 3000.0, 5000.0, 1e8, n_max=40, n_times=200)
    times = r.relaxation_times_s
    assert len(times) >= 1
    assert all(times[i] >= times[i + 1] for i in range(len(times) - 1))
    assert times[0] == pytest.approx(r.tau_dominant_s, rel=1e-9)


@pytest.mark.parametrize("bad_kwargs", [
    {"c_total_M": 0.0, "K_initial_per_M": 3000.0, "K_final_per_M": 5000.0, "k_on_per_M_s": 1e8},
    {"c_total_M": 0.01, "K_initial_per_M": 0.0, "K_final_per_M": 5000.0, "k_on_per_M_s": 1e8},
    {"c_total_M": 0.01, "K_initial_per_M": 3000.0, "K_final_per_M": 5000.0, "k_on_per_M_s": 0.0},
    {"c_total_M": 0.01, "K_initial_per_M": 3000.0, "K_final_per_M": 5000.0, "k_on_per_M_s": 1e8, "n_max": 2},
])
def test_o1_error_paths(bad_kwargs):
    with pytest.raises(ValueError):
        micelle_relaxation_odes(**bad_kwargs)


# ---------------------------------------------------------------------------
# S2 -- kmc_micelle_kinetics
# ---------------------------------------------------------------------------

def test_s2_matches_mass_action_statistically():
    """Stochastic equilibrium must land within a real statistical band of
    the closed-form target for the realization's effective mass -- not
    exact equality (that's the whole point of it being stochastic)."""
    k = kmc_micelle_kinetics(0.01, 5000.0, 1e8, n_max=40, n_molecules=3000,
                              n_samples=240, seed=0)
    rel_diff = abs(k.equilibrium_c1_M - k.mass_action_c1_M) / k.mass_action_c1_M
    assert rel_diff < 0.02
    nagg_rel_diff = abs(k.equilibrium_mean_nagg - k.mass_action_mean_nagg) / k.mass_action_mean_nagg
    assert nagg_rel_diff < 0.1


def test_s2_exact_integer_mass_conservation():
    k = kmc_micelle_kinetics(0.01, 5000.0, 1e8, n_max=40, n_molecules=3000,
                              n_samples=240, seed=0)
    assert k.mass_conservation_max_error == 0
    assert k.converged is True


def test_s2_seed_determinism():
    k1 = kmc_micelle_kinetics(0.01, 5000.0, 1e8, n_max=40, n_molecules=3000,
                               n_samples=240, seed=0)
    k2 = kmc_micelle_kinetics(0.01, 5000.0, 1e8, n_max=40, n_molecules=3000,
                               n_samples=240, seed=0)
    assert k1.trajectory_c1_M == k2.trajectory_c1_M
    assert k1.trajectory_mean_nagg == k2.trajectory_mean_nagg
    assert k1.n_events == k2.n_events


def test_s2_different_seeds_differ():
    k1 = kmc_micelle_kinetics(0.01, 5000.0, 1e8, n_max=40, n_molecules=3000,
                               n_samples=240, seed=0)
    k2 = kmc_micelle_kinetics(0.01, 5000.0, 1e8, n_max=40, n_molecules=3000,
                               n_samples=240, seed=1)
    assert k1.trajectory_c1_M != k2.trajectory_c1_M


def test_s2_vs_o1_relaxation_time_same_order_of_magnitude():
    """tau_fit_s (KMC, micelle-count-biased) and tau_dominant_s (O1, true
    slowest linear mode) are explicitly DIFFERENT projections per both
    docstrings -- this is a loose cross-model sanity check, not a tight
    physical equivalence claim."""
    c_total, K, k_on, n_max = 0.01, 5000.0, 1e8, 40
    o1 = micelle_relaxation_odes(c_total, 3000.0, K, k_on, n_max=n_max, n_times=200)
    s2 = kmc_micelle_kinetics(c_total, K, k_on, n_max=n_max, n_molecules=3000,
                               n_samples=240, seed=0)
    assert not math.isnan(s2.tau_fit_s)
    ratio = s2.tau_fit_s / o1.tau_dominant_s
    assert 0.1 < ratio < 10.0


@pytest.mark.parametrize("bad_kwargs", [
    {"c_total_M": 0.0, "K_per_M": 5000.0, "k_on_per_M_s": 1e8},
    {"c_total_M": 0.01, "K_per_M": 0.0, "k_on_per_M_s": 1e8},
    {"c_total_M": 0.01, "K_per_M": 5000.0, "k_on_per_M_s": 0.0},
    {"c_total_M": 0.01, "K_per_M": 5000.0, "k_on_per_M_s": 1e8, "n_max": 2},
    {"c_total_M": 0.01, "K_per_M": 5000.0, "k_on_per_M_s": 1e8, "n_molecules": 100},
    {"c_total_M": 0.01, "K_per_M": 5000.0, "k_on_per_M_s": 1e8, "n_samples": 5},
])
def test_s2_error_paths(bad_kwargs):
    with pytest.raises(ValueError):
        kmc_micelle_kinetics(**bad_kwargs)


# ---------------------------------------------------------------------------
# O1b -- aniansson_wall_two_regime
#
# Gold case: real SDS (NaLS) multi-technique kinetics data from Lang,
# Tondre, Zana, Bauer, Hoffmann & Ulbricht, J. Phys. Chem. 79 (1975) 276,
# Table I at C=cmc -- tau1=29us (shock tube), tau2=2.3ms (P-jump); real
# SDS n_bar=62 at cmc=8.2mM, 25C. sigma defaults to Aniansson & Wall's own
# stated typical value 0.20*n_bar=12.4 (not independently measured here).
# Values below are hand-derived closed-form solutions, confirmed by an
# independent manual calculation and an exact round-trip before this
# function existed.
# ---------------------------------------------------------------------------

_SDS_TAU1_S = 29e-6
_SDS_TAU2_S = 2.3e-3
_SDS_CMC_M = 0.0081
_SDS_NBAR = 62.0


def test_aniansson_wall_matches_hand_derived_sds_gold():
    r = aniansson_wall_two_regime(
        tau1_s=_SDS_TAU1_S, tau2_s=_SDS_TAU2_S, cmc_M=_SDS_CMC_M, n_bar=_SDS_NBAR,
    )
    assert r.sigma == pytest.approx(12.4, rel=1e-9)
    assert r.a_dimensionless == pytest.approx(0.0, abs=1e-12)
    assert r.k_minus_per_s == pytest.approx(5302068.97, rel=1e-4)
    assert r.R_s_per_M == pytest.approx(1091.51, rel=1e-3)
    assert r.k_on_implied_per_M_per_s == pytest.approx(654576415.5, rel=1e-4)
    assert r.average_micellar_lifetime_s == pytest.approx(0.1426, rel=1e-3)
    assert r.converged is True


def test_aniansson_wall_implied_k_on_is_diffusion_limited():
    """Unforced physical sanity check: for a small 1:1 ionic surfactant,
    k_on_implied should land near known diffusion-limited association
    rates (~1e8-1e9 M^-1 s^-1) -- not a tautology, since k_minus and R are
    solved purely from tau1/tau2/n_bar/sigma/cmc with no k_on target."""
    r = aniansson_wall_two_regime(
        tau1_s=_SDS_TAU1_S, tau2_s=_SDS_TAU2_S, cmc_M=_SDS_CMC_M, n_bar=_SDS_NBAR,
    )
    assert 1e8 < r.k_on_implied_per_M_per_s < 1e9


def test_aniansson_wall_round_trip_reproduces_tau1_tau2_exactly():
    """Plugging the solved k_minus/R back into the forward tau1/tau2
    formulas must reproduce the original inputs exactly (closed-form
    algebraic inversion, not a fit)."""
    r = aniansson_wall_two_regime(
        tau1_s=_SDS_TAU1_S, tau2_s=_SDS_TAU2_S, cmc_M=_SDS_CMC_M, n_bar=_SDS_NBAR,
    )
    factor = r.params["shared_factor"]
    tau1_reproduced = (r.sigma ** 2 / r.k_minus_per_s) / factor
    tau2_reproduced = (r.R_s_per_M * _SDS_CMC_M) / (_SDS_NBAR ** 2 * factor)
    assert tau1_reproduced == pytest.approx(_SDS_TAU1_S, rel=1e-9)
    assert tau2_reproduced == pytest.approx(_SDS_TAU2_S, rel=1e-9)


def test_aniansson_wall_default_sigma_is_aniansson_wall_typical_value():
    r = aniansson_wall_two_regime(
        tau1_s=_SDS_TAU1_S, tau2_s=_SDS_TAU2_S, cmc_M=_SDS_CMC_M, n_bar=_SDS_NBAR,
    )
    assert r.sigma == pytest.approx(0.20 * _SDS_NBAR, rel=1e-12)


def test_aniansson_wall_explicit_sigma_overrides_default():
    r_default = aniansson_wall_two_regime(
        tau1_s=_SDS_TAU1_S, tau2_s=_SDS_TAU2_S, cmc_M=_SDS_CMC_M, n_bar=_SDS_NBAR,
    )
    r_explicit = aniansson_wall_two_regime(
        tau1_s=_SDS_TAU1_S, tau2_s=_SDS_TAU2_S, cmc_M=_SDS_CMC_M, n_bar=_SDS_NBAR,
        sigma=10.0,
    )
    assert r_explicit.sigma == pytest.approx(10.0, rel=1e-12)
    assert r_explicit.k_minus_per_s != pytest.approx(r_default.k_minus_per_s, rel=1e-6)


def test_aniansson_wall_c_total_defaults_to_cmc():
    r_default = aniansson_wall_two_regime(
        tau1_s=_SDS_TAU1_S, tau2_s=_SDS_TAU2_S, cmc_M=_SDS_CMC_M, n_bar=_SDS_NBAR,
    )
    r_explicit = aniansson_wall_two_regime(
        tau1_s=_SDS_TAU1_S, tau2_s=_SDS_TAU2_S, cmc_M=_SDS_CMC_M, n_bar=_SDS_NBAR,
        c_total_M=_SDS_CMC_M,
    )
    assert r_default.k_minus_per_s == pytest.approx(r_explicit.k_minus_per_s, rel=1e-12)
    assert r_default.a_dimensionless == pytest.approx(0.0, abs=1e-12)


def test_aniansson_wall_c_total_above_cmc_changes_shared_factor():
    r_at_cmc = aniansson_wall_two_regime(
        tau1_s=_SDS_TAU1_S, tau2_s=_SDS_TAU2_S, cmc_M=_SDS_CMC_M, n_bar=_SDS_NBAR,
    )
    r_above_cmc = aniansson_wall_two_regime(
        tau1_s=_SDS_TAU1_S, tau2_s=_SDS_TAU2_S, cmc_M=_SDS_CMC_M, n_bar=_SDS_NBAR,
        c_total_M=2 * _SDS_CMC_M,
    )
    assert r_above_cmc.a_dimensionless == pytest.approx(1.0, rel=1e-9)
    assert r_above_cmc.params["shared_factor"] > r_at_cmc.params["shared_factor"]
    assert r_above_cmc.k_minus_per_s != pytest.approx(r_at_cmc.k_minus_per_s, rel=1e-6)


@pytest.mark.parametrize("bad_kwargs", [
    {"tau1_s": 0.0, "tau2_s": _SDS_TAU2_S, "cmc_M": _SDS_CMC_M, "n_bar": _SDS_NBAR},
    {"tau1_s": _SDS_TAU1_S, "tau2_s": 0.0, "cmc_M": _SDS_CMC_M, "n_bar": _SDS_NBAR},
    {"tau1_s": _SDS_TAU1_S, "tau2_s": _SDS_TAU2_S, "cmc_M": 0.0, "n_bar": _SDS_NBAR},
    {"tau1_s": _SDS_TAU1_S, "tau2_s": _SDS_TAU2_S, "cmc_M": _SDS_CMC_M, "n_bar": 1.0},
    {"tau1_s": _SDS_TAU1_S, "tau2_s": _SDS_TAU2_S, "cmc_M": _SDS_CMC_M, "n_bar": _SDS_NBAR,
         "sigma": 0.0},
    {"tau1_s": _SDS_TAU1_S, "tau2_s": _SDS_TAU2_S, "cmc_M": _SDS_CMC_M, "n_bar": _SDS_NBAR,
         "c_total_M": 0.0},
    {"tau1_s": _SDS_TAU1_S, "tau2_s": _SDS_TAU2_S, "cmc_M": _SDS_CMC_M, "n_bar": _SDS_NBAR,
         "sigma": 100.0, "c_total_M": 0.1 * _SDS_CMC_M},
])
def test_aniansson_wall_error_paths(bad_kwargs):
    with pytest.raises(ValueError):
        aniansson_wall_two_regime(**bad_kwargs)
