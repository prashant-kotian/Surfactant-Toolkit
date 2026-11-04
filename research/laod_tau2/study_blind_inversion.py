"""Blind inversion trial on the exact solver - the reliability question.

WHY THIS EXISTS
---------------
Every recovery study in this repository starts the fit near the truth and reports
the width of the cost minimum. That measures curvature, not reliability.

FINDINGS 20.2 then showed that delta is NON-MONOTONE in tau2: there is a locus
where d(delta)/d ln tau2 passes through zero and changes sign. On that locus the
cost surface has a flat direction, so a bounded local least-squares can return a
confident-looking tau2 that is wrong, and nothing in the fit diagnostics says so.

This module asks the question that matters for a user of the method:
    given a dataset, what does the protocol return, and is it right?

WHY THERE IS NO SIMILARITY TABLE HERE
-------------------------------------
FINDINGS 21 retracted the claim that delta is a function of omega*tau2 alone: at
matched Omega the lag changes by 11 degrees depending on how the product is split.
So this module tabulates delta(tau2) separately AT EACH MEASURED FREQUENCY and at
several capacities. Interpolation is then exact along the only axis the fitter
moves, and no similarity assumption enters anywhere.

The table is still an interpolation, so it is validated against the solver at
random off-grid points before any trial is run (see `validate_table`).

PROTOCOLS
---------
A : capacity known from the n_bar closure, fit tau2.  (the paper's headline mode)
B : capacity also free, fit (tau2, C_m).              (the exploratory mode)

Each is run twice on identical data:
  local      - bounded least squares from a generic start (tau2 = 0.3 s)
  multistart - the same fit from 7 starts, keep the lowest cost

Outputs research/laod_tau2/blind_inversion.json
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy.interpolate import RectBivariateSpline, CubicSpline
from scipy.optimize import least_squares

from . import harmonics as H
from .model import ModelParams, run_oscillation
from .systems import SYSTEMS

R_GAS = 8.314462618
S = SYSTEMS["C12E8"]
T_REF = 298.15
STRAIN = 0.20
C_TOT = 1.5 * S.cmc_mM
SIGMA_G = 0.05          # mN/m
OMEGAS = (3.14, 6.28, 10.0)
CAPS = (5.0, 10.0, 20.0, 47.5, 85.0)
TAU_GRID = np.exp(np.linspace(math.log(0.05), math.log(50.0), 40))
NPC, NCY, NRET = 160, 14, 8

# truths: tau2 chosen to straddle the null locus, which for C_m=47.5 lies near
# Omega ~ 6, i.e. tau2 ~ 1 s at 6.28 rad/s.  (FINDINGS 20.2)
TRUTHS = [
    (4.00, 47.5),   # the paper's operating point
    (1.00, 47.5),   # ON the null locus at 6.28 rad/s -> the dangerous case
    (0.50, 47.5),   # below the null (positive branch)
    (2.00, 47.5),   # just above the null
    (25.0, 47.5),   # far past the peak
    (1.00, 20.0),   # null at a different capacity
    (0.50, 10.0),
]
N_REPEATS = 10


def K_for_E0(E0: float) -> float:
    pref = S.n_ions * R_GAS * T_REF * S.gamma_inf * 1e3
    r = E0 / pref
    th = r / (1.0 + r)
    return th / (1.0 - th) / C_TOT


K_ANCHOR = K_for_E0(45.3)


def solver_lag(tau2: float, cap: float, omega: float, sigma: float = 0.0,
               seed: int = 0) -> float:
    p = ModelParams(c0=C_TOT, K=K_ANCHOR, gamma_inf=S.gamma_inf, D=S.D,
                    n_factor=float(S.n_ions), capacity=cap, tau2=tau2,
                    temperature_K=T_REF)
    r = run_oscillation(p, omega, STRAIN, n_cycles=NCY, n_return_cycles=NRET,
                        n_per_cycle=NPC)
    t = np.asarray(r.t)
    g = H.add_noise(t, r.gamma_mN_m, sigma, seed=seed) if sigma else np.asarray(r.gamma_mN_m)
    return H.lag_deg(t, g, r.area_rel, omega)


def build_table(force: bool = False) -> dict:
    """delta(tau2) at each (measured frequency, capacity) on TAU_GRID."""
    cache = Path(__file__).with_name("blind_table_v2.npz")
    if cache.exists() and not force:
        d = np.load(cache)
        return {om: {cap: CubicSpline(TAU_GRID, d[f"{om}_{cap}"]) for cap in CAPS}
                for om in OMEGAS}
    data, out = {}, {}
    for om in OMEGAS:
        out[om] = {}
        for cap in CAPS:
            z = np.array([solver_lag(t2, cap, om) for t2 in TAU_GRID])
            data[f"{om}_{cap}"] = z
            out[om][cap] = CubicSpline(TAU_GRID, z)
            print(f"  om={om:5.2f} cap={cap:6.1f}  z range "
                  f"{z.min():.1f}-{z.max():.1f}", flush=True)
    np.savez(cache, **data)
    return out


_T = None


def table() -> dict:
    global _T
    if _T is None:
        _T = build_table()
    return _T


def tab_lag(tau2: float, cap: float, om: float) -> float:
    """Interpolate in tau2 (exact axis); interpolate linearly across capacities."""
    t = table()[om]
    caps = np.array(CAPS, float)
    vals = np.array([t[c](min(max(tau2, TAU_GRID[0]), TAU_GRID[-1])) for c in CAPS])
    if cap <= caps[0]:
        return float(vals[0])
    if cap >= caps[-1]:
        return float(vals[-1])
    return float(np.interp(math.log(cap), np.log(caps), vals))


def validate_table(n: int = 6, seed: int = 999) -> dict:
    """Check the table against the solver at random off-grid points."""
    rng = np.random.default_rng(seed)
    errs = []
    for _ in range(n):
        t2 = float(np.exp(rng.uniform(math.log(0.3), math.log(20.0))))
        om = float(rng.choice(OMEGAS))
        cap = 47.5
        exact = solver_lag(t2, cap, om)
        errs.append(tab_lag(t2, cap, om) - exact)
    errs = np.abs(errs)
    return {"max_abs_err_deg": float(errs.max()), "mean_abs_err_deg": float(errs.mean())}


def fit(obs: dict, cap_fixed: float | None, start: float) -> dict:
    def res(x):
        tau2 = math.exp(x[0])
        cap = cap_fixed if cap_fixed is not None else math.exp(x[1])
        return [tab_lag(tau2, cap, om) - obs[om] for om in OMEGAS]

    x0 = [math.log(start)] if cap_fixed is not None else [math.log(start), math.log(30.0)]
    lo = [math.log(TAU_GRID[0])] if cap_fixed is not None else [math.log(TAU_GRID[0]), math.log(CAPS[0])]
    hi = [math.log(TAU_GRID[-1])] if cap_fixed is not None else [math.log(TAU_GRID[-1]), math.log(CAPS[-1])]
    sol = least_squares(res, x0=x0, bounds=(lo, hi), xtol=1e-12, ftol=1e-12,
                        diff_step=1e-4, max_nfev=400)
    dof = max(len(sol.fun) - len(x0), 1)
    s2 = max(2 * sol.cost / dof, 1e-12)
    try:
        cov = np.linalg.inv(sol.jac.T @ sol.jac) * s2
        se = float(np.sqrt(np.diag(cov))[0])
    except np.linalg.LinAlgError:
        se = float("nan")
    r = {"tau2": float(math.exp(sol.x[0])), "relSE": se, "cost": float(2 * sol.cost),
         "at_bound": bool(min(abs(sol.x[0] - lo[0]), abs(sol.x[0] - hi[0])) < 1e-3)}
    if cap_fixed is None:
        r["capacity"] = float(math.exp(sol.x[1]))
    return r


def trial(t2: float, cap: float, seed: int, protocol: str, multi: bool) -> dict:
    obs = {om: solver_lag(t2, cap, om, SIGMA_G, seed=seed * 37 + int(om * 11))
           for om in OMEGAS}
    cap_arg = cap if protocol == "A" else None
    if multi:
        best = None
        for st in (0.06, 0.15, 0.4, 1.0, 3.0, 10.0, 30.0):
            r = fit(obs, cap_arg, st)
            if best is None or r["cost"] < best["cost"]:
                best = r
        r = best
    else:
        r = fit(obs, cap_arg, 0.3)
    err = r["tau2"] / t2 - 1.0
    se = r["relSE"] if np.isfinite(r["relSE"]) else 9.0
    return {"truth": t2, "cap_truth": cap, "est": r["tau2"], "relerr": err,
            "relSE": se, "wrong": abs(err) > 0.5,
            "confident_wrong": abs(err) > 0.5 and se < 0.1,
            "at_bound": r["at_bound"], "reported_cap": r.get("capacity"),
            "om_null": r["tau2"] * 6.28}


def main() -> dict:
    build_table()
    val = validate_table()
    print(f"\ntable validation: max |err| = {val['max_abs_err_deg']:.3f} deg, "
          f"mean = {val['mean_abs_err_deg']:.3f} deg\n", flush=True)
    out = {"validation": val,
           "design": {"omegas": list(OMEGAS), "sigma_g": SIGMA_G, "capacities": list(CAPS),
                      "n_repeats": N_REPEATS, "cycles": NCY, "tau_grid": [float(TAU_GRID[0]), float(TAU_GRID[-1]), len(TAU_GRID)]},
           "protocols": {}}
    for protocol in ("A", "B"):
        for multi in (False, True):
            key = f"{protocol}_{'multistart' if multi else 'local'}"
            rows = [trial(t2, cap, s, protocol, multi)
                    for (t2, cap) in TRUTHS for s in range(N_REPEATS)]
            nw = sum(r["wrong"] for r in rows)
            nc = sum(r["confident_wrong"] for r in rows)
            ok = [abs(r["relerr"]) for r in rows if not r["wrong"]]
            out["protocols"][key] = {
                "n": len(rows), "n_wrong": nw, "n_confident_wrong": nc,
                "at_bound": sum(r["at_bound"] for r in rows),
                "median_abserr_on_success": float(np.median(ok)) if ok else None,
                "rows": rows}
            med = f"{np.median(ok)*100:.2f}%" if ok else "n/a"
            print(f"{key:18} n={len(rows):3d} wrong={nw:3d} "
                  f"confident-wrong={nc:3d} at-bound={out['protocols'][key]['at_bound']:3d} "
                  f"median|err|ok={med}", flush=True)
    return out


if __name__ == "__main__":
    res = main()
    path = Path(__file__).with_name("blind_inversion.json")
    path.write_text(json.dumps(res, indent=2, default=float))
    print(f"\nwrote {path}")
