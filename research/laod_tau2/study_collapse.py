"""Two-axis collapse: is delta a function of (omega*tau2, C_m) alone?

FINDINGS 19.1 showed delta collapses on omega*tau2 to within 3 degrees across a
fourfold range in E0, for a single capacity. That was one column of a two-dimensional
problem. This module tests the full claim:

    delta = f(omega*tau2, C_m)      (approximately, for a given strain)

If it holds, the design rule stops being a three-system heuristic and becomes a
master chart: an experimenter looks up (their tau2, their capacity, their frequency)
and reads off the expected lag and hence the expected sensitivity.

Grid: omega*tau2 in {0.25 ... 64} x C_m in {2, 5, 10, 20, 47.5, 85, 200}.
Each cell is run at three E0 anchors (25, 45.3, 103.5 mN/m) so that the E0
independence is tested AT EVERY CAPACITY rather than at one.

Reported per cell: delta at each anchor, the spread across anchors, and the local
lever arm d(delta)/d ln tau2 (which is what determines precision).

Outputs research/laod_tau2/collapse_map.json
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
STRAIN = 0.20
C_TOT = 1.5 * S.cmc_mM
TAU2 = 4.0
OMEGA_TAU = (0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0)
CAPACITIES = (2.0, 5.0, 10.0, 20.0, 47.5, 85.0, 200.0)
E0_ANCHORS = (25.0, 45.3, 103.5)


def K_for_E0(E0: float) -> float:
    """Langmuir K placing E0 at the working concentration (a = 0)."""
    pref = S.n_ions * R_GAS * T_REF * S.gamma_inf * 1e3
    r = E0 / pref
    th = r / (1.0 + r)
    return th / (1.0 - th) / C_TOT


def lag(K: float, cap: float, tau2: float, omega: float) -> float:
    p = ModelParams(c0=C_TOT, K=K, gamma_inf=S.gamma_inf, D=S.D, n_factor=float(S.n_ions),
                    capacity=cap, tau2=tau2, temperature_K=T_REF)
    r = run_oscillation(p, omega, STRAIN, n_cycles=14, n_return_cycles=6)
    t = np.asarray(r.t)
    return H.lag_deg(t, r.gamma_mN_m, r.area_rel, omega)


def cell(cap: float, omega_tau: float, Ks: dict) -> dict:
    """delta at each E0 anchor, plus the lever arm at the middle anchor."""
    omega = omega_tau / TAU2
    vals = {f"{e:.0f}": lag(K, cap, TAU2, omega) for e, K in Ks.items()}
    arr = list(vals.values())
    # lever arm: +/- a factor of 1.25 in tau2 at the middle anchor
    K_mid = Ks[E0_ANCHORS[1]]
    lo = lag(K_mid, cap, TAU2 / 1.25, omega)
    hi = lag(K_mid, cap, TAU2 * 1.25, omega)
    slope = (hi - lo) / math.log(1.25 ** 2)
    return {"delta": vals, "spread": max(arr) - min(arr), "lever": abs(slope)}


def main() -> dict:
    Ks = {e: K_for_E0(e) for e in E0_ANCHORS}
    out = {"strain": STRAIN, "tau2": TAU2, "c_over_cmc": 1.5,
           "omega_tau": list(OMEGA_TAU), "capacities": list(CAPACITIES),
           "e0_anchors": list(E0_ANCHORS), "cells": {}}
    for cap in CAPACITIES:
        row = {}
        for wt in OMEGA_TAU:
            row[f"{wt:g}"] = cell(cap, wt, Ks)
        out["cells"][f"{cap:g}"] = row
        print(f"  C_m = {cap:6.1f} done", flush=True)
    return out


if __name__ == "__main__":
    res = main()
    path = Path(__file__).with_name("collapse_map.json")
    path.write_text(json.dumps(res, indent=2, default=float))

    print("\n=== delta at omega*tau2, by capacity (E0 = 45.3 anchor) ===")
    hdr = f"{'C_m':>7} |" + "".join(f"{w:>7g}" for w in OMEGA_TAU)
    print(hdr); print("-" * len(hdr))
    for cap in CAPACITIES:
        row = res["cells"][f"{cap:g}"]
        print(f"{cap:7.1f} |" + "".join(f"{row[f'{w:g}']['delta']['45']:7.1f}" for w in OMEGA_TAU))

    print("\n=== lever arm |d delta / d ln tau2|, degrees per e-fold ===")
    print(hdr); print("-" * len(hdr))
    for cap in CAPACITIES:
        row = res["cells"][f"{cap:g}"]
        print(f"{cap:7.1f} |" + "".join(f"{row[f'{w:g}']['lever']:7.2f}" for w in OMEGA_TAU))

    print("\n=== E0 spread across anchors, degrees ===")
    print(hdr); print("-" * len(hdr))
    worst = 0.0
    for cap in CAPACITIES:
        row = res["cells"][f"{cap:g}"]
        print(f"{cap:7.1f} |" + "".join(f"{row[f'{w:g}']['spread']:7.2f}" for w in OMEGA_TAU))
        worst = max(worst, max(row[f"{w:g}"]["spread"] for w in OMEGA_TAU))
    print(f"\nworst E0 spread anywhere on the grid: {worst:.2f} deg")
    print(f"wrote {path}")
