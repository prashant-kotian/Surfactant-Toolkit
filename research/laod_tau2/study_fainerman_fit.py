"""Fit our coupled-reservoir model to Fainerman et al.'s MEASURED phase angles.

WHY
---
Every quantitative claim in this project has been validated against its own forward
model. Fainerman, Petkov & Miller (Langmuir 2008, 24, 6447, Figure 5) publish measured
phase angle vs frequency for C14EO8 above the cmc. That is the only external micellar
dilatational dataset available, and until now we have only compared single points to it.

WHAT IS FITTED
--------------
The published curve shape phi(f) for the highest concentration (30 umol/L = 4.3 x cmc)
over 0.005-0.2 Hz, against two parameters: capacity C_m and relaxation time tau2.
Everything else is fixed (D, Gamma_inf, E0 anchored, strain 0.20).

DIGITISATION AND ITS LIMITS
---------------------------
The points below were read off a 700 dpi rendering of Figure 5 by eye, not by a
calibrated digitiser, so they carry roughly +/- 1 degree. Six points, two parameters,
so four degrees of freedom. This is a shape test, not a precision measurement: it can
tell us whether the model has the right functional form and roughly the right
parameters, and it cannot validate the model to better than a few degrees. Stated in
the output as well, because a fit like this is easy to over-read.

Also noted: the data are for C14EO8 while most of our work is C12E8. Both are CnE8
nonionics with similar CMC magnitude. Only the *shape* transfers.

Outputs research/laod_tau2/fainerman_fit.json
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from . import harmonics as H
from .model import ModelParams, run_oscillation
from .systems import SYSTEMS

R_GAS = 8.314462618
S = SYSTEMS["C12E8"]
T_REF = 298.15
STRAIN = 0.20                 # Fainerman's oscillations were small-amplitude
C_CMC_RATIO = 30.0 / 7.0      # 30 umol/L against cmc = 7 umol/L
C_TOT = C_CMC_RATIO * S.cmc_mM

# digitised from Figure 5, filled squares = 30 umol/L
FREQS = (0.005, 0.010, 0.020, 0.050, 0.100, 0.200)
PHI_OBS = (46.5, 50.5, 52.0, 48.5, 41.0, 32.5)
PHI_ERR = 1.0                 # digitisation uncertainty, degrees


def K_for_E0(E0: float) -> float:
    pref = S.n_ions * R_GAS * T_REF * S.gamma_inf * 1e3
    r = E0 / pref
    th = r / (1.0 + r)
    return th / (1.0 - th) / C_TOT


K_ANCHOR = K_for_E0(45.3)


def phi_model(cap: float, tau2: float, f: float) -> float:
    p = ModelParams(c0=C_TOT, K=K_ANCHOR, gamma_inf=S.gamma_inf, D=S.D,
                    n_factor=float(S.n_ions), capacity=cap, tau2=tau2,
                    temperature_K=T_REF)
    om = 2 * math.pi * f
    r = run_oscillation(p, om, STRAIN, n_cycles=18, n_return_cycles=8,
                        n_per_cycle=160)
    t = np.asarray(r.t)
    return H.lag_deg(t, r.gamma_mN_m, r.area_rel, om)


def sse(cap: float, tau2: float) -> float:
    return float(sum((phi_model(cap, tau2, f) - obs) ** 2
                     for f, obs in zip(FREQS, PHI_OBS)))


def main() -> dict:
    caps = (1.0, 3.0, 10.0, 30.0, 60.0, 120.0, 250.0, 313.0)
    taus = (0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 20.0)

    print(f"C14EO8 at {C_CMC_RATIO:.1f} x cmc; closure would give C_m = "
          f"{S.capacity(C_TOT):.0f}\n")
    hdr = f"{'C_m':>7} |" + "".join(f"{t:>7g}" for t in taus)
    print(hdr)

    grid = {}
    best = None
    for cap in caps:
        row = []
        for t2 in taus:
            v = sse(cap, t2)
            grid[f"{cap}_{t2}"] = v
            row.append(v)
            if best is None or v < best[2]:
                best = (cap, t2, v)
        print(f"{cap:7.1f} |" + "".join(f"{v:7.0f}" for v in row), flush=True)

    cap_b, tau_b, sse_b = best
    dof = len(FREQS) - 2
    rmse = math.sqrt(sse_b / dof)

    print(f"\nbest grid point: C_m = {cap_b:g}, tau2 = {tau_b:g} s")
    print(f"SSE = {sse_b:.1f}, RMSE = {rmse:.2f} deg over {len(FREQS)} points")
    print(f"{'f(Hz)':>8} {'obs':>7} {'model':>8} {'resid':>8}")
    for f, obs in zip(FREQS, PHI_OBS):
        m = phi_model(cap_b, tau_b, f)
        print(f"{f:8.3f} {obs:7.1f} {m:8.1f} {m - obs:8.1f}")

    closure = S.capacity(C_TOT)
    print(f"\nclosure predicts C_m = {closure:.0f}; best fit is {cap_b:g}. "
          f"ratio {closure / cap_b:.0f}x")

    return {"fitted": {"capacity": cap_b, "tau2": tau_b, "sse": sse_b,
                       "rmse_deg": rmse, "dof": dof},
            "closure_capacity": closure,
            "data": {"c_over_cmc": C_CMC_RATIO, "freqs": list(FREQS),
                     "phi_obs": list(PHI_OBS), "phi_err": PHI_ERR,
                     "source": "Fainerman 2008 Fig 5, 30 umol/L, visual digitisation"},
            "grid": {"caps": list(caps), "taus": list(taus), "sse": grid}}


if __name__ == "__main__":
    res = main()
    path = Path(__file__).with_name("fainerman_fit.json")
    path.write_text(json.dumps(res, indent=2, default=float))
    print(f"\nwrote {path}")
