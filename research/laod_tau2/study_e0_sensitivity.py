"""Does the C12E8 operability window survive a literature-anchored elasticity?

PROBLEM CAUGHT IN SELF-REVIEW (2026-11-03). Every quantitative result so far
(C12E8 window, recovery uncertainty, Protocol B) was computed with
K = 140 m^3/mol, which sits at theta = 0.937 and gives E0 = 103.5 mN/m at
1.5 x CMC. That is the SAME saturated regime that section 12 flags as a defect
in the supplied CTAB/SDS constants. If I condemn those numbers, I have to
recheck mine.

Czakaj's own Table 4 supplies an independent anchor: E0 = 45.3 mN/m at 0.5 mM
and 73.9 at 0.3 mM for the LAE system (CMC 1.0 mM). So this system genuinely
does reach high elasticities -- it is not a ~25 mN/m surfactant. The question
is whether the window and the tau2 leverage are specific to E0 = 103 or
survive across the plausible range.

Sweep E0 = 25, 45 (Czakaj's own), 75, and 103.5 (K = 140) at 1.5 x CMC,
omega = 6.28, and report the window edges and the tau2 leverage at each.

Outputs research/laod_tau2/e0_sensitivity.json.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from . import harmonics as H
from .lae_data import NOISE_FLOOR_RATIO
from .model import ModelParams, run_oscillation
from .systems import SYSTEMS

R_GAS = 8.314462618
S = SYSTEMS["C12E8"]
NOISE = NOISE_FLOOR_RATIO[1]
STRAIN = 0.20
T_REF = 298.15
TAU2 = 4.0


def K_for_E0(E0, c):
    """Langmuir K giving E0 at concentration c (mol/m^3)."""
    pref = S.n_ions * R_GAS * T_REF * S.gamma_inf * 1e3
    r = E0 / pref
    th = r / (1.0 + r)
    return th / (1.0 - th) / c, th


def probe(p, omega):
    r = run_oscillation(p, omega, STRAIN, n_cycles=16)
    t = np.asarray(r.t)
    hg = H.harmonics(t, r.gamma_mN_m, omega)
    return hg[2], H.lag_deg(t, r.gamma_mN_m, r.area_rel, omega)


def main():
    out = {"tau2": TAU2, "omega": 6.28, "noise": NOISE, "rows": []}
    for E0 in (25.0, 45.3, 73.9, 103.5):
        # anchor at 1.5 x CMC; hold K fixed as concentration varies
        c_anchor = 1.5 * S.cmc_mM
        K, th_anchor = K_for_E0(E0, c_anchor)
        window = []
        for mult in (1.0, 1.2, 1.5, 2.0, 3.0, 5.0, 10.0):
            c = mult * S.cmc_mM
            p = ModelParams(c0=c, K=K, gamma_inf=S.gamma_inf, D=S.D,
                            n_factor=float(S.n_ions), capacity=S.capacity(c),
                            tau2=TAU2, temperature_K=T_REF)
            i2, lag = probe(p, 6.28)
            plo = ModelParams(**{**p.__dict__, "tau2": TAU2 / 2.0})
            phi = ModelParams(**{**p.__dict__, "tau2": TAU2 * 2.0})
            lag_dbl = probe(phi, 6.28)[1]
            lag_half = probe(plo, 6.28)[1]
            window.append({"mult": mult, "I2": i2, "lag": lag,
                           "dlag_per_4x": abs(lag_dbl - lag_half),
                           "detectable": i2 > NOISE,
                           "theta": p.theta0})
        best = max(window, key=lambda r: r["dlag_per_4x"] if r["detectable"] else -1)
        usable = [r["mult"] for r in window
                  if r["detectable"] and r["dlag_per_4x"] >= 3.0]
        out["rows"].append({
            "E0": E0, "K": K, "theta_1.5x": th_anchor,
            "window": window,
            "best_mult": best["mult"], "best_dlag": best["dlag_per_4x"],
            "usable_mult": usable,
        })
    return out


if __name__ == "__main__":
    res = main()
    path = Path(__file__).with_name("e0_sensitivity.json")
    path.write_text(json.dumps(res, indent=2, default=float))
    print(f"C12E8, tau2 = {TAU2} s, omega = 6.28, noise = {NOISE}\n")
    for d in res["rows"]:
        print(f"--- E0 = {d['E0']:.1f} mN/m  (K = {d['K']:.3f}, "
              f"theta@1.5x = {d['theta_1.5x']:.3f}) ---")
        print(f"{'mult':>6} {'theta':>6} {'I2':>8} {'lag':>7} {'dlag/4x':>8}  flag")
        for r in d["window"]:
            flag = ("det" if r["detectable"] else "BELOW NOISE")
            if r["detectable"] and r["dlag_per_4x"] < 3.0:
                flag += " / insensitive"
            print(f"{r['mult']:6.1f} {r['theta']:6.3f} {r['I2']:8.4f} "
                  f"{r['lag']:7.2f} {r['dlag_per_4x']:8.2f}  {flag}")
        print(f"best concentration = {d['best_mult']}x CMC "
              f"(dlag {d['best_dlag']:.2f} deg); usable = {d['usable_mult']}\n")
    print(f"wrote {path}")
