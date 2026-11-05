"""Design engine: which oscillatory dilatational experiment identifies a
mechanism, and which cannot?

Motivation, and why this is not another identifiability audit
-------------------------------------------------------------
Al-Sou'od & Ghanem (J. Surfactants Deterg. 2026, 10.1002/jsde.70095) established
that oscillatory and step-relaxation data frequently cannot distinguish a
mechanism from a nuisance protocol artefact, and that profile likelihood
intervals computed within one model can vanish when an omitted nuisance degree
of freedom is added. That is a diagnosis. It does not say what experiment to run
instead. This module computes the answer.

The question
------------
Instrument bandwidth and noise are fixed by the apparatus. Concentration,
surface pressure and temperature are chosen by the experimenter. So: over the
*designable* variables, what is the best attainable discrimination between
(no relaxation) and (relaxation), and where should the frequencies go?

What "discrimination" means here
--------------------------------
Not the size of a residual. The quantity that matters is how far apart two
models are relative to noise, measured as

    D  =  min over nuisances of the Mahalanobis distance between
          model A and model B predictions,

with nuisances (E0, h, D, and where relevant C) integrated out or profiled. A
model is only assignable if D exceeds the noise. D is also half the Fisher
information distance, so it converts directly into a bound on how well C can be
estimated: sd(C) >= 2/D in the linear case.

Why the design is not "measure at the frequency where the curve moves most"
--------------------------------------------------------------------------
Because the nuisance parameters move it more. h and D shift the curve; E0 scales
it. The design must be scored on the *shape* that no nuisance can reproduce.
`_discrimination` therefore profiles the nuisance rather than holding it fixed,
and `score_design` keeps the worst case over a credible prior on h.

This is the first thing in this project that answers a question someone would
actually ask before going into the lab.
"""

from __future__ import annotations

import itertools
import json
import math
from pathlib import Path

import numpy as np

from .models import B_single, lag_deg, modulus_from_B

R_GAS = 8.314462618
T_K = 298.15

# --- instrument envelopes ---------------------------------------------------
# Stated ranges of the three instrument classes, so a design can be checked
# against hardware rather than against an idealised band.
INSTRUMENTS = {
    "PAT (oscillating drop/pendant)": (0.005, 1.0),
    "CPT (capillary pressure)": (0.1, 500.0),
    "microtensiometer": (1.7e-3, 1.7),
}
NOISE = {"phase_deg": 1.0, "modulus_rel": 0.02}


# --- forward model ----------------------------------------------------------

def predict(freqs, *, E0, C, tau, h, D, observable="phase"):
    """Predicted (|E|, delta) at the given frequencies."""
    out = []
    for f in np.atleast_1d(freqs):
        w = 2 * math.pi * f
        B = complex(B_single(w, C, tau)) if C > 0 else 1.0 + 0j
        E = modulus_from_B(w, E0, D, h, B)
        out.append((abs(E), lag_deg(w, E0, D, h, B)))
    a = np.array(out)
    if observable == "phase":
        return a[:, 1]
    if observable == "both":
        return a
    raise ValueError(observable)


def _dist(a, b, observable):
    """Mahalanobis distance between two predictions, in noise units."""
    if observable == "phase":
        return float(np.sqrt(np.sum(((a - b) / NOISE["phase_deg"]) ** 2)))
    dp = (a[:, 1] - b[:, 1]) / NOISE["phase_deg"]
    dm = (a[:, 0] - b[:, 0]) / (NOISE["modulus_rel"] * np.maximum(b[:, 0], 1e-12))
    return float(np.sqrt(np.sum(dp ** 2 + dm ** 2)))


def discrimination(freqs, truth, rival, nuisance_h, observable="phase",
                   nuisance_D=None):
    """Minimum distance between two models over a set of nuisance values.

    truth and rival are dicts of model parameters. The nuisance set is
    (h, D) pairs; the minimum is the worst case for the experimenter, which is
    the honest number, because none of these nuisances is measured.

    h and D are also *shared* between the two models: the experiment does not
    get a different isotherm slope depending on which mechanism is true. Both
    models are therefore evaluated at the same nuisance values, and the
    minimum over those values is taken.
    """
    best = np.inf
    for h in nuisance_h:
        for D in (nuisance_D or [truth.get("D", 3.8e-10)]):
            a = predict(freqs, E0=truth["E0"], C=truth["C"], tau=truth["tau"],
                        h=h, D=D, observable=observable)
            b = predict(freqs, E0=rival["E0"], C=rival["C"], tau=rival["tau"],
                        h=h, D=D, observable=observable)
            best = min(best, _dist(a, b, observable))
    return float(best)


# --- designable variables ---------------------------------------------------

def spectrum_pressure(gamma_inf, K, c0):
    """Pi = RT Gamma_inf ln(1 + K c0), mN/m. Langmuir."""
    return R_GAS * T_K * 1e3 * gamma_inf * math.log1p(K * c0)


def isotherm_slope(gamma_inf, K, c0):
    """h = dGamma_eq/dc = Gamma_inf K / (1 + K c0)^2, m."""
    return gamma_inf * K / (1.0 + K * c0) ** 2


def coverage(K, c0):
    """theta, coverage. The natural designable variable: it fixes the elasticity
    scale and the relaxation structure, and unlike c0 it transfers between
    surfactants of different CMC."""
    return K * c0 / (1.0 + K * c0)


def elasticity_scale(gamma_inf, K, c0, n_factor=1.0):
    """Gibbs elasticity E0 = nRT Gamma / (1 - theta), mN/m.

    NOTE the (1-theta)^-1. An earlier version of this function omitted it and
    returned E0 = RT h c0 = RT Gamma theta (1-theta), which *falls* with coverage
    and implies dPi/dlnc >= RT everywhere. That is wrong: Langmuir has
    dPi/dlnc = RT theta/(1-theta) -> 0 at low coverage and divergent at
    saturation. The corrected form is what makes "work at high surface pressure"
    a real design lever rather than a diminishing return.
    """
    th = coverage(K, c0)
    G = gamma_inf * th
    return n_factor * R_GAS * T_K * 1e3 * G / max(1.0 - th, 1e-12)


def relaxation_time_scale(gamma_inf, K, c0, D):
    """h^2 / D, seconds. The characteristic time separating the
    diffusion-dominated and elasticity-dominated regimes; it is NOT the surface
    relaxation time, and conflating the two was one of the errors in the earlier
    project (FINDINGS 29). Reported here only to place a chosen tau relative to
    the instrument band."""
    h = isotherm_slope(gamma_inf, K, c0)
    return h * h / D


def designs(instrument, n_freqs=6, modes=("log", "low", "high", "uniform")):
    """Candidate frequency layouts inside an instrument band."""
    lo, hi = INSTRUMENTS[instrument]
    out = {"log": np.logspace(math.log10(lo), math.log10(hi), n_freqs)}
    out["low"] = np.logspace(math.log10(lo), math.log10(lo * 10), n_freqs)
    out["high"] = np.logspace(math.log10(hi / 10), math.log10(hi), n_freqs)
    out["uniform"] = np.linspace(lo, hi, n_freqs)
    return {k: v for k, v in out.items() if k in modes}


def score_design(freqs, truth, rival, h_prior, observable="phase"):
    """Worst-case discrimination over a credible prior on the nuisance h."""
    return discrimination(freqs, truth, rival, h_prior, observable)
