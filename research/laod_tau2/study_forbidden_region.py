"""The forbidden region: where the inversion sensitivity vanishes.

The two-axis collapse map (collapse_map.json) revealed something more useful than
the collapse itself. The lever arm |d delta / d ln tau2| is not monotone in Omega:
it passes through ZERO wherever delta has a maximum in Omega. A fine scan at
C_m = 20 confirms it: delta rises to 74.11 deg at Omega ~ 3.5, and the lever arm is
0.17 deg/e-fold at Omega = 4. Past the maximum the sign reverses.

That matters operationally. An experimenter who happens to sit on that locus gets a
precise, high-signal, and completely uninformative measurement: the lag barely moves
when tau2 changes, and on the far side it moves the wrong way, so a least-squares
inversion can settle on a wrong-branch solution while reporting a small error bar.

This module locates the zero-sensitivity locus across capacity, i.e. the forbidden
operating region of the design chart.

Method: for each C_m, scan Omega on a fine grid, compute the lever arm, and report
the Omega at which it is minimised together with the branch structure (whether the
lag is still rising or has turned over).

Outputs research/laod_tau2/forbidden_region.json
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from . import harmonics as H
from .model import ModelParams, run_oscillation
from .study_collapse import C_TOT, S, T_REF, K_for_E0

TAU2 = 4.0
STRAIN = 0.20
CAPACITIES = (2.0, 5.0, 10.0, 20.0, 47.5, 85.0, 200.0)
# Fine grid around the expected extremum, spanning the useful band.
OMEGA_TAU = (1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 8.0, 10.0, 12.0, 16.0)
K = K_for_E0(45.3)


def lag(cap: float, tau2: float, omega: float) -> float:
    p = ModelParams(c0=C_TOT, K=K, gamma_inf=S.gamma_inf, D=S.D, n_factor=float(S.n_ions),
                    capacity=cap, tau2=tau2, temperature_K=T_REF)
    r = run_oscillation(p, omega, STRAIN, n_cycles=14, n_return_cycles=6)
    t = np.asarray(r.t)
    return H.lag_deg(t, r.gamma_mN_m, r.area_rel, omega)


def scan(cap: float) -> list[dict]:
    rows = []
    for wt in OMEGA_TAU:
        om = wt / TAU2
        d = lag(cap, TAU2, om)
        lo = lag(cap, TAU2 / 1.25, om)
        hi = lag(cap, TAU2 * 1.25, om)
        # signed slope: positive means delta rises as tau2 rises
        signed = (hi - lo) / math.log(1.25 ** 2)
        rows.append({"omega_tau": wt, "omega": om, "delta": d,
                     "lever_signed": signed, "lever": abs(signed)})
    return rows


def find_null(rows: list[dict]) -> dict:
    """Omega at minimum |lever|, and whether delta has turned over there."""
    i = min(range(len(rows)), key=lambda j: rows[j]["lever"])
    r = rows[i]
    peak = max(rows, key=lambda x: x["delta"])
    return {
        "omega_tau_at_null": r["omega_tau"],
        "delta_at_null": r["delta"],
        "min_lever": r["lever"],
        "delta_peak_omega_tau": peak["omega_tau"],
        "delta_peak": peak["delta"],
        "sign_at_null": 1 if r["lever_signed"] > 0 else -1,
        "sign_before": 1 if rows[max(i - 1, 0)]["lever_signed"] > 0 else -1,
        "sign_after": 1 if rows[min(i + 1, len(rows) - 1)]["lever_signed"] > 0 else -1,
    }


def main() -> dict:
    out = {"strain": STRAIN, "tau2": TAU2, "e0": 45.3,
           "omega_tau_grid": list(OMEGA_TAU), "capacities": {}}
    for cap in CAPACITIES:
        rows = scan(cap)
        out["capacities"][f"{cap:g}"] = {"rows": rows, "null": find_null(rows)}
        print(f"  C_m = {cap:6.1f}  null at Omega ~ {out['capacities'][f'{cap:g}']['null']['omega_tau_at_null']:5.1f}"
              f"  (min lever {out['capacities'][f'{cap:g}']['null']['min_lever']:.2f} deg/e-fold)", flush=True)
    return out


if __name__ == "__main__":
    res = main()
    path = Path(__file__).with_name("forbidden_region.json")
    path.write_text(json.dumps(res, indent=2, default=float))

    print("\n=== signed lever arm, deg per e-fold (negative = lag falls as tau2 rises) ===")
    hdr = f"{'C_m':>7} |" + "".join(f"{w:>7g}" for w in OMEGA_TAU)
    print(hdr); print("-" * len(hdr))
    for cap in CAPACITIES:
        rows = res["capacities"][f"{cap:g}"]["rows"]
        print(f"{cap:7.1f} |" + "".join(f"{r['lever_signed']:>7.1f}" for r in rows))

    print("\n=== null locus ===")
    print(f"{'C_m':>7} {'Omega_null':>11} {'min lever':>10} {'sign flips?':>12}")
    for cap in CAPACITIES:
        n = res["capacities"][f"{cap:g}"]["null"]
        flip = "YES" if n["sign_before"] != n["sign_after"] else "no"
        print(f"{cap:7.1f} {n['omega_tau_at_null']:>11.1f} {n['min_lever']:>10.2f} {flip:>12}")
    print(f"\nwrote {path}")
