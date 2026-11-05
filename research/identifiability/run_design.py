"""The design study: where should the frequencies go, and at what concentration?

Produces the figures and tables for the paper's central section. Everything is
computed; nothing is transcribed.

Three questions, in order of how much they are worth:

Q1  Given an instrument band and a fixed number of frequencies, which frequency
    layout maximises the worst-case discrimination between "relaxation" and "no
    relaxation", once the nuisance isotherm slope is allowed to vary over its
    credible range?
Q2  Given freedom to choose the concentration, where should the experiment be
    run? The designable variable is coverage theta, because it sets both the
    elasticity scale E0 = RT Gamma/(1-theta) and the characteristic time
    h^2/D, and it transfers between surfactants with different CMCs.
Q3  Does adding |E| to the phase rescue a design that phase alone cannot?

The output is the paper's answer to "what experiment should I run", which is
what the 2026 identifiability literature poses and does not answer.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from .design_engine import (
    INSTRUMENTS, NOISE, coverage, designs, discrimination, elasticity_scale,
    isotherm_slope, predict, relaxation_time_scale, spectrum_pressure,
)

# Units: mol/m^3 throughout for concentration. 7 umol/L = 7e-3 mol/m^3.
# K and gamma_inf are the values that reproduce this project's calibrated
# reference point (h = 3.2753e-6 m, |delta| = 40.08 deg at 0.1 Hz, 71 umol/L).
# Fainerman's b2 = 1e5 is in L/mol-like units and must NOT be used directly;
# taking it at face value puts theta ~ 1 at 5 umol/L, which is unphysical.
GAMMA_INF = 2.8e-6      # mol/m^2
K_ADS = 140.0           # m^3/mol
D_MONO = 4.0e-10        # m^2/s, Fainerman's fitted monomer D
CMC = 7.0e-3            # mol/m^3 (= 7 umol/L)

# The rival is "no surface relaxation" (C = 0, LVDT). The truth is a small
# relaxation of the size Fainerman's own ceiling implies.
C_TRUE = 0.6
TAU_TRUE = 12.0         # s, order of the digitised C14E8 relaxation
H_PRIOR_FACTOR = np.array([0.5, 0.7, 1.0, 1.5, 2.0])


def system_at(c0):
    """Material properties at a chosen bulk concentration."""
    h = isotherm_slope(GAMMA_INF, K_ADS, c0)
    return dict(E0=elasticity_scale(GAMMA_INF, K_ADS, c0), h=h,
                D=D_MONO, theta=coverage(K_ADS, c0),
                Pi=spectrum_pressure(GAMMA_INF, K_ADS, c0),
                tau_diff=relaxation_time_scale(GAMMA_INF, K_ADS, c0, D_MONO))


def truth_rival(c0, C=C_TRUE, tau=TAU_TRUE):
    s = system_at(c0)
    truth = dict(E0=s["E0"], C=C, tau=tau, D=s["D"], h_ref=s["h"])
    rival = dict(E0=s["E0"], C=0.0, tau=tau, D=s["D"], h_ref=s["h"])
    return truth, rival, s


def q1_frequency_layout(c0, observable="phase", n_freqs=6):
    truth, rival, s = truth_rival(c0)
    h_prior = [s["h"] * f for f in H_PRIOR_FACTOR]
    rows = []
    for inst in INSTRUMENTS:
        for name, fr in designs(inst, n_freqs, ("log", "low", "high", "uniform")).items():
            d = discrimination(fr, truth, rival, h_prior, observable)
            rows.append(dict(instrument=inst, layout=name, freqs=fr.tolist(), D=d))
    return rows


def q2_concentration(observable="phase", instrument="PAT (oscillating drop/pendant)",
                     layout="log", n_freqs=6):
    lo, hi = INSTRUMENTS[instrument]
    fr = designs(instrument, n_freqs, (layout,))[layout]
    rows = []
    for c_ratio in [1, 2, 5, 10, 20, 50, 100, 200, 500]:
        c0 = c_ratio * CMC
        truth, rival, s = truth_rival(c0)
        h_prior = [s["h"] * f for f in H_PRIOR_FACTOR]
        d = discrimination(fr, truth, rival, h_prior, observable)
        rows.append(dict(c_over_cmc=c_ratio, c0_uM=c0 * 1e6, theta=s["theta"],
                         Pi=s["Pi"], E0=s["E0"], tau_diff=s["tau_diff"], D=d))
    return rows


def q3_observable(c0, instrument="PAT (oscillating drop/pendant)", n_freqs=6):
    truth, rival, s = truth_rival(c0)
    h_prior = [s["h"] * f for f in H_PRIOR_FACTOR]
    fr = designs(instrument, n_freqs, ("log",))["log"]
    return {obs: discrimination(fr, truth, rival, h_prior, obs)
            for obs in ("phase", "both")}


def tau_sweep(c0, observables=("phase", "both")):
    """For which surface relaxation times is the experiment informative at all?
    A design cannot be blamed for a time constant outside its band."""
    truth, rival, s = truth_rival(c0)
    h_prior = [s["h"] * f for f in H_PRIOR_FACTOR]
    out = {}
    taus = np.logspace(-3, 4, 29)
    for obs in observables:
        vals = []
        for tau in taus:
            t, r, _ = truth_rival(c0, tau=tau)
            fr = designs("PAT (oscillating drop/pendant)", 6, ("log",))["log"]
            vals.append(discrimination(fr, t, r, h_prior, obs))
        out[obs] = dict(tau=taus.tolist(), D=vals)
    return out


def run(save=True):
    res = dict(
        constants=dict(gamma_inf=GAMMA_INF, K=K_ADS, D=D_MONO, cmc_uM=CMC * 1e6,
                       C_true=C_TRUE, tau_true=TAU_TRUE, noise=NOISE,
                       h_prior=list(H_PRIOR_FACTOR)),
        q1=q1_frequency_layout(20 * CMC),
        q2=q2_concentration(),
        q3=q3_observable(100 * CMC),
        tau_sweep=tau_sweep(100 * CMC),
    )
    if save:
        Path("research/identifiability/design_results.json").write_text(
            json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    r = run()
    print("Design study. Distance D is in noise units (1 deg phase, 2% modulus),")
    print("minimised over a factor-2 prior on the nuisance isotherm slope.")
    print("D > 3 is conventionally 'detected'; D > 5 is comfortably assignable.\n")

    print("Q1  frequency layout, c/c_mc = 20")
    print(f"{'instrument':32s}{'layout':10s}{'D':>7s}")
    for row in sorted(r["q1"], key=lambda x: (x["instrument"], -x["D"])):
        print(f"{row['instrument']:32s}{row['layout']:10s}{row['D']:7.2f}")

    print("\nQ2  concentration, PAT log layout")
    print(f"{'c/c_mc':>7s}{'theta':>7s}{'Pi':>7s}{'E0':>8s}{'h^2/D':>10s}{'D':>7s}")
    for row in r["q2"]:
        print(f"{row['c_over_cmc']:7.0f}{row['theta']:7.3f}{row['Pi']:7.2f}"
              f"{row['E0']:8.1f}{row['tau_diff']:10.2e}{row['D']:7.2f}")

    print("\nQ3  phase alone vs |E| + phase, c/c_mc = 100")
    for k, v in r["q3"].items():
        print(f"  {k:7s} D = {v:.2f}")

    print("\ntau sweep (c/c_mc = 100): informative range at D > 3")
    for obs, d in r["tau_sweep"].items():
        arr = np.array(d["D"]); tau = np.array(d["tau"])
        ok = tau[arr > 3]
        rng = f"{ok.min():.2e} - {ok.max():.2e} s" if len(ok) else "none"
        print(f"  {obs:7s} {rng}")
