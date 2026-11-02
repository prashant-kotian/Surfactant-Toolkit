"""Anchor tests: the forward model must reproduce known limits, or nothing
downstream means anything.

Run:  PYTHONPATH=src:. python -m pytest research/laod_tau2/tests -q

Four anchors, in order of how much they constrain:
  1. LVDT (small-amplitude diffusion-controlled) -- the model must reproduce
     the Lucassen-van den Tempel closed form from the full PDE, with no
     fitting. This is the anchor that caught the first structural bug
     (writing c_s instead of Gamma as the state gave 1/E0 elasticity).
  2. Linear limit: higher harmonics must vanish as eps -> 0. If they do not,
     they are a numerical artefact.
  3. Fast-reservoir limit (omega*tau2 -> 0): with capacity C_m the response
     must equal a micelle-free system with D_eff = D*(1+C_m).
  4. Slow-reservoir limit (omega*tau2 -> infinity): micelles frozen, so the
     response must equal the micelle-free system at the same D.
Anchors 3 and 4 bracket every intermediate micellar response.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from research.laod_tau2 import harmonics as H
from research.laod_tau2 import lvdt
from research.laod_tau2.model import ModelParams, run_oscillation

BASE = dict(c0=1.0, K=5.0, gamma_inf=3.5e-6, D=4.0e-10, gamma0=72.0)


def _amp(p, omega, eps, **kw):
    r = run_oscillation(p, omega, eps, n_cycles=kw.pop("n_cycles", 10), **kw)
    g = np.asarray(r.gamma_mN_m)
    return g.max() - g.min()


@pytest.mark.parametrize("omega", [0.1, 1.0, 10.0])
def test_lvdt_reproduced_from_pde(omega):
    """Full reaction-diffusion model vs the LVDT closed form.

    NOTE: the nominal LVDT derivation neglects the kappa*(c10-c) term only in
    the absence of micelles, so this test is run with capacity = 0.
    """
    p = ModelParams(**BASE)
    st = lvdt.LangmuirState(**BASE)
    ep, epp = st.lvdt(omega)
    predicted = 2.0 * 0.02 * abs(complex(ep, epp))
    measured = _amp(p, omega, 0.02)
    assert measured == pytest.approx(predicted, rel=0.02), (measured, predicted)


def test_insoluble_limit_approaches_gibbs_elasticity():
    """D -> 0 (no exchange within a cycle): amplitude -> 2*eps*E0, monotonically.

    The residual gap at finite D is real (the interface is still marginally
    diffusion-relaxed), so the test asserts the monotone approach and a
    tolerance that the approach actually meets, rather than a number chosen
    to pass.
    """
    st = lvdt.LangmuirState(**BASE)
    target = 2.0 * 0.02 * st.E0
    amps = [_amp(ModelParams(**{**BASE, "D": d}), 10.0, 0.02) for d in (1e-11, 1e-12, 1e-13)]
    assert amps[0] < amps[1] < amps[2], amps
    assert amps[2] == pytest.approx(target, rel=0.15), (amps[2], target)


def test_harmonics_vanish_in_linear_limit():
    """I2/I1 must scale down with eps: the 2nd harmonic is a nonlinear effect."""
    p = ModelParams(**BASE)
    ratios = []
    for eps in (0.01, 0.02, 0.04):
        r = run_oscillation(p, 1.0, eps, n_cycles=10)
        hg = H.harmonics(np.asarray(r.t), r.gamma_mN_m, 1.0)
        ratios.append(hg[2] / hg[1])
    assert ratios[0] < 0.02
    assert ratios[0] < ratios[1] < ratios[2]


def test_fast_reservoir_equals_effective_diffusivity():
    """omega*tau2 -> 0: response == micelle-free system with D*(1+C_m)."""
    cap, tau2, omega = 2.0, 1e-6, 1.0
    with_micelles = ModelParams(**BASE, capacity=cap, tau2=tau2)
    effective = ModelParams(**{**BASE, "D": BASE["D"] * (1.0 + cap)})
    a = _amp(with_micelles, omega, 0.10)
    b = _amp(effective, omega, 0.10)
    assert a == pytest.approx(b, rel=0.05), (a, b)


def test_slow_reservoir_equals_frozen_micelles():
    """omega*tau2 -> infinity: response == micelle-free system at the same D."""
    with_micelles = ModelParams(**BASE, capacity=2.0, tau2=1e6)
    plain = ModelParams(**BASE)
    a = _amp(with_micelles, 1.0, 0.10)
    b = _amp(plain, 1.0, 0.10)
    assert a == pytest.approx(b, rel=0.02), (a, b)


def test_equilibrium_is_stationary():
    """eps = 0: nothing should evolve (guards against a spurious drift term)."""
    p = ModelParams(**BASE, capacity=1.0, tau2=0.1)
    r = run_oscillation(p, 1.0, 0.0, n_cycles=6)
    g = np.asarray(r.gamma_mN_m)
    assert np.ptp(g) < 1e-4
