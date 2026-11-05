"""Tests for the identifiability module.

These exist for one reason: the previous project's headline claims were
retracted because a sign convention and an algebra slip went unnoticed. Every
non-obvious step in this module is therefore pinned by a test against an
independent reference, not against itself.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from research.identifiability.models import (  # noqa: E402
    B_log_normal, B_single, B_two_mode, capacity_from_cap, lag_cap_deg, lag_deg,
    max_cap_single_relaxation, modulus_from_B,
)
from research.laod_tau2.linear_response import linear_response  # noqa: E402

D = 3.8e-10
H = 3.275302547717482e-06
E0 = 60.0


# --- 1. the model must reproduce the independently validated reference -----

def test_lvdt_matches_calibrated_reference():
    """B = 1 must reproduce laod_tau2.linear_response to machine precision.

    This is the single most important test in the file: it pins both the
    algebra and the sign convention against a function that was itself
    verified against a PDE solver in the earlier work.
    """
    for f in [1e-4, 1e-3, 1e-2, 0.1, 1.0, 10.0, 1e3]:
        omega = 2 * math.pi * f
        ref = abs(linear_response(omega, c0=0.071, K=140.0,
                                  gamma_inf=2.8e-6, D=D, capacity=0.0).lag_deg)
        got = lag_deg(omega, E0, D, H, 1.0 + 0j)
        assert abs(ref - got) < 1e-9, (f, ref, got)


def test_B_one_recovers_lvdt_everywhere():
    """B_single with Cm = 0 must be identical to LVDT, for any tau."""
    for tau in [1e-3, 1.0, 1e3]:
        for f in [1e-3, 1.0]:
            omega = 2 * math.pi * f
            assert lag_deg(omega, E0, D, H, B_single(omega, 0.0, tau)) == pytest.approx(
                lag_deg(omega, E0, D, H, 1.0 + 0j), abs=1e-12)


# --- 2. the 45 degree ceiling and the generalisation -----------------------

def test_lvdt_ceiling_is_45():
    """LVDT lag approaches 45 deg as frequency falls and never exceeds it."""
    freqs = [1e-1, 1e-2, 1e-3, 1e-4, 1e-6, 1e-8]
    vals = [lag_deg(2 * math.pi * f, E0, D, H, 1.0 + 0j) for f in freqs]
    assert all(v < 45.0 + 1e-9 for v in vals)
    # convergence to the ceiling is slow (like 45 - c sqrt(f)), so the limit
    # is checked at 1e-8 Hz rather than at a frequency someone might measure
    assert vals[-1] == pytest.approx(45.0, abs=1e-2)
    assert vals == sorted(vals), "lag must increase monotonically as f falls"


def test_cap_formula_is_exact():
    """|delta| -> 45 + |arg B|/2 as f -> 0, for several admittances.

    Verified against the numerical limit, so this is a test of the derivation
    rather than a restatement of it.
    """
    cases = {
        "single Cm=2 tau=8": B_single(2 * math.pi * 1e-7, 2.0, 8.0),
        "two-mode": B_two_mode(2 * math.pi * 1e-7, 1.5, 1.0, 50.0),
        "log-normal": B_log_normal(2 * math.pi * 1e-7, 2.0, 5.0, 1.5),
    }
    for name, B in cases.items():
        numeric = lag_deg(2 * math.pi * 1e-9, E0, D, H, complex(B))
        predicted = lag_cap_deg(0.0, complex(B))
        assert abs(numeric - predicted) < 5e-3, (name, numeric, predicted)


def test_cap_rises_towards_90_and_is_not_bound_by_67_5():
    """CORRECTED. The round-27 claim that a passive relaxation cannot exceed
    67.5 deg was an artefact of scanning tau at one frequency.

    The maximum of |arg B| over frequency is at w tau = sqrt(C+1):
        tan(|arg B|_max) = C / (2 sqrt(C+1)),
    which has no finite bound below 45 deg. This test checks the analytic
    optimum against a two-dimensional numerical search, so it fails if either
    the formula or the reasoning is wrong.
    """
    for C in [0.5, 2.0, 20.0, 200.0]:
        analytic = max_cap_single_relaxation(C)
        best = 0.0
        for tau in np.logspace(-6, 6, 120):
            for f in np.logspace(-8, 2, 60):
                B = B_single(2 * math.pi * f, C, tau)
                best = max(best, lag_cap_deg(0.0, B))
        assert abs(best - analytic) < 0.05, (C, best, analytic)
    assert max_cap_single_relaxation(1e6) > 89.0


def test_capacity_from_cap_inverts_the_ceiling():
    """The useful direction: a measured ceiling turns into a bound on C."""
    for C in [0.1, 0.5, 2.0, 10.0, 100.0]:
        d = max_cap_single_relaxation(C)
        assert capacity_from_cap(d) == pytest.approx(C, rel=1e-9)
    assert capacity_from_cap(45.0) == 0.0
    assert capacity_from_cap(50.0) == pytest.approx(0.420, abs=1e-3)


def test_fainerman_ceiling_bounds_capacity():
    """Fainerman 2008 reports a maximum phase angle of about 50 deg for C14EO8
    (text: 'the phase angle increases up to a value of 50 deg', and at
    30 umol/L the digitised maximum is 52 deg). That bounds C to ~0.4-0.6
    without fitting any curve -- the argument is a ceiling, not a regression.
    """
    assert 0.3 < capacity_from_cap(50.0) < 0.5
    assert 0.55 < capacity_from_cap(52.0) < 0.75


def test_published_data_lies_inside_the_bound():
    """The C14E8 data (max 52 deg) must satisfy whatever bound the theory gives.
    If a test here ever rejects real published data, the theory is wrong, not
    the data."""
    obs = [46.5, 50.5, 52.0, 48.5, 41.0, 32.5]
    assert max(obs) < 90.0 and max(obs) <= max_cap_single_relaxation(1e6)


# --- 3. structural facts used in the interpretation ------------------------

def test_E0_cancels_from_the_phase():
    """Phase is independent of E0. This is why phase-only data cannot
    identify the Gibbs elasticity, and the fit therefore excludes it."""
    for E in [1.0, 60.0, 400.0]:
        a = lag_deg(2 * math.pi * 0.05, 1.0, D, H, B_single(1.0, 1.0, 1.0))
        b = lag_deg(2 * math.pi * 0.05, E, D, H, B_single(1.0, 1.0, 1.0))
        assert a == pytest.approx(b, abs=1e-12)


def test_two_mode_and_log_normal_can_degenerate_to_single():
    """The extra freedom in the richer models can be switched off, so they can
    never fit worse than the simpler model. Hence a chi-square improvement
    alone is not evidence for the extra mechanism."""
    omega = 2 * math.pi * 0.05
    single = B_single(omega, 1.2, 7.0)
    assert B_two_mode(omega, 1.2, 7.0, 7.0) == pytest.approx(single, abs=1e-12)


def test_modulus_is_finite_and_positive_over_the_band():
    """Guard against the overflow/NaN that an unbounded tau produces."""
    for tau in np.logspace(-4, 3, 20):
        E = modulus_from_B(2 * math.pi * 0.1, E0, D, H, B_single(2 * math.pi * 0.1, 5.0, tau))
        assert np.isfinite(E) and abs(E) > 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])


# --- 4. the degeneracy that survives the full complex observable -----------

def test_only_h_over_sqrt_D_is_identifiable():
    """Scaling h and D together so that h/sqrt(D) is fixed leaves E unchanged.

    Consequence: no amount of phase *or amplitude* precision separates the
    isotherm slope from the diffusion coefficient. An independent value of one
    is required. This is a structural degeneracy, not a noise limit, and it is
    the reason the amplitude test reports h/sqrt(D) = 1.0000 while h and D
    individually wander by factors of several.
    """
    for f in [0.005, 0.05, 0.5]:
        w = 2 * math.pi * f
        B = complex(B_single(w, 1.0, 5.0))
        h0, D0 = H, D
        k = 4.0                     # scale h by 4, D by 16 -> h/sqrt(D) fixed
        E_ref = modulus_from_B(w, E0, D0, h0, B)
        E_scaled = modulus_from_B(w, E0, D0 * k ** 2, h0 * k, B)
        assert abs(E_ref - E_scaled) < 1e-12 * abs(E_ref)


def test_phase_is_blind_to_the_amplitude_of_E():
    """The phase carries no information about E0 at all, for any frequency.

    This is why the phase-only fit returns E0 = 350 mN/m from data generated
    with E0 = 72: it is not a bad fit, it is an unconstrained parameter.
    """
    vals = set()
    for E_scale in [0.5, 1.0, 7.0, 50.0]:
        vals.add(round(lag_deg(2 * math.pi * 0.03, E0 * E_scale, D, H,
                               B_single(2 * math.pi * 0.03, 2.0, 9.0)), 10))
    assert len(vals) == 1


# --- 5. design engine: invariants that must hold ---------------------------

def test_gibbs_elasticity_rises_with_coverage():
    """E0 = RT Gamma/(1-theta) must increase monotonically with concentration.

    An earlier version of the engine returned RT*h*c0 = RT*Gamma*theta*(1-theta),
    which falls with coverage and inverted every design conclusion. Langmuir has
    dPi/dlnc = RT*theta/(1-theta) -> 0 at low coverage and divergent at saturation,
    so E0 must rise, steeply, and without bound until the isotherm breaks down.
    """
    from research.identifiability.design_engine import elasticity_scale, isotherm_slope
    gi, K = 2.8e-6, 140.0
    cs = np.logspace(-6, -1, 20)          # mol/m^3, i.e. 1 uM to 100 mM
    E = [elasticity_scale(gi, K, c) for c in cs]
    assert all(E[i] < E[i + 1] for i in range(len(E) - 1)), "E0 must increase with c"


def test_discrimination_is_zero_for_identical_models():
    """A model cannot be distinguished from itself, at any noise level."""
    from research.identifiability.design_engine import discrimination
    fr = np.logspace(-3, 0, 6)
    m = dict(E0=50.0, C=0.5, tau=10.0, D=D)
    assert discrimination(fr, m, m, [H]) == pytest.approx(0.0, abs=1e-12)


def test_discrimination_increases_with_relaxation_strength():
    """Bigger surface relaxation must be easier to detect, monotonically."""
    from research.identifiability.design_engine import discrimination
    fr = np.logspace(-3, 0, 8)
    ds = [discrimination(fr, dict(E0=50.0, C=C, tau=10.0, D=D),
                         dict(E0=50.0, C=0.0, tau=10.0, D=D), [H])
          for C in (0.05, 0.2, 0.6, 2.0)]
    assert all(ds[i] < ds[i + 1] for i in range(len(ds) - 1))


def test_worst_case_nuisance_is_a_minimum_not_an_average():
    """The metric must be the minimum over the nuisance prior, so widening the
    prior can only reduce it. This is what makes the number honest: a design is
    not allowed to look good only at the nominal isotherm slope."""
    from research.identifiability.design_engine import discrimination
    fr = np.logspace(-3, 0, 8)
    t = dict(E0=50.0, C=0.6, tau=10.0, D=D)
    r = dict(E0=50.0, C=0.0, tau=10.0, D=D)
    narrow = discrimination(fr, t, r, [H])
    wide = discrimination(fr, t, r, [H * f for f in (0.3, 0.5, 1.0, 2.0, 3.3)])
    assert wide <= narrow + 1e-12


# --- 6. structure -> design chain ----------------------------------------- -

def test_cmc_from_structure_is_physical():
    """The CMC must fall steeply with chain length and land in a plausible range.

    A C12 polyoxyethylene ether CMC near 1e-4 mol/L is right; 1e-7 would mean a
    unit error of three orders, which happened twice in this project (mol/L vs
    mol/m3, and a stray factor in an f-string that printed a 500x-CMC
    concentration as 59 mM when it was 59 uM).
    """
    from research.identifiability.from_structure import cmc_from_structure
    c12 = cmc_from_structure(12)["cmc_mol_per_L"]
    c16 = cmc_from_structure(16)["cmc_mol_per_L"]
    assert 1e-5 < c12 < 1e-3, c12
    assert c16 < c12, "CMC must decrease with chain length"
    # Roughly a decade per two carbons for ethoxylates. The Nagarajan
    # prediction sits at the low end of the experimental scatter for the long
    # chains (C16E8 is measured around 1e-6 mol/L), so the window is set to
    # the predictor's own validated range rather than to a literature value.
    assert 1e-7 < c16 < 5e-6, c16
    assert c12 / c16 > 20.0, (c12, c16)


def test_structure_design_returns_monotone_sweep():
    """Discrimination must be non-decreasing along the concentration sweep, and
    must saturate rather than diverge."""
    from research.identifiability.from_structure import design_for_surfactant
    s = design_for_surfactant(14)
    assert "error" not in s
    Ds = [r["D"] for r in s["sweep"]]
    assert all(Ds[i] <= Ds[i + 1] + 1e-9 for i in range(len(Ds) - 1))
    assert max(Ds) < 100.0
