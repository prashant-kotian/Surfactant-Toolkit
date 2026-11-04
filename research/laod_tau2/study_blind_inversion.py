"""Blind inversion trial: run the published protocol as a naive user would.

WHY THIS EXISTS
---------------
The recovery studies so far start the fit near the truth and report the width of
the cost minimum. That measures curvature, not reliability. FINDINGS 20.2 then
showed delta is NON-MONOTONE in tau2: there is a locus where the sensitivity is
exactly zero and the sign flips either side. On that locus a bounded local
least-squares will happily return a confident-looking tau2 that is wrong, because
the cost surface has a flat direction there.

So the question is not "how precise is the estimator near the truth" but
    "given a dataset, what does the protocol return, and is it right?"

DESIGN
------
1. A forward table delta(Omega, C_m) is built ONCE from the solver. This exploits
   the collapse established in FINDINGS 19.1/20.1: delta is (approximately) a
   function of Omega = omega*tau2 and capacity alone, with E0 entering weakly for
   Omega <~ 8. Using the table rather than the solver makes an otherwise
   intractable fit feasible, at the cost of the (bounded) E0 dependence; that is
   acceptable here because the question is which BRANCH the estimator lands on,
   not the third significant figure.
2. Truths are chosen across the parameter space, including points whose data span
   the null locus.
3. Synthetic data are generated from the SOLVER (not the table), so the table is
   never the thing being fitted to itself.
4. The published protocol is run blind: bounded local least squares from a
   generic start. Then a multi-start search on the same data.
5. A trial counts as a failure if the estimate is wrong by more than 50% --
   including the dangerous case of a wrong branch returned with a small
   reported uncertainty.

Outputs research/laod_tau2/blind_inversion.json
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy.interpolate import RectBivariateSpline
from scipy.optimize import least_squares

from . import harmonics as H
from .model import ModelParams, run_oscillation
from .systems import SYSTEMS

R_GAS = 8.314462618
S = SYSTEMS["C12E8"]
T_REF = 298.15
STRAIN = 0.20
C_TOT = 1.5 * S.cmc_mM
SIGMA_G = 0.05
N_RETURN_CYCLES = 8
OMEGAS = (3.14, 6.28, 10.0)      # a 3-point sweep; each is a separate measurement

# Forward table grid in (Omega, C_m)
_G_OM = np.exp(np.linspace(math.log(0.05), math.log(120.0), 64))
_G_CAP = np.exp(np.linspace(math.log(0.5), math.log(600.0), 32))

TRUTHS = [
    (4.0, 47.5),     # the paper's operating point
    (4.0, 20.0),
    (1.0, 47.5),
    (12.0, 47.5),
    (0.5, 10.0),
    (40.0, 85.0),
]
N_REPEATS = 8


def K_for_E0(E0: float) -> float:
    pref = S.n_ions * R_GAS * T_REF * S.gamma_inf * 1e3
    r = E0 / pref
    th = r / (1.0 + r)
    return th / (1.0 - th) / C_TOT


K_ANCHOR = K_for_E0(45.3)


def solver_lag(tau2: float, cap: float, omega: float, sigma: float = 0.0,
               seed: int = 0) -> float:
    """Lag from the actual ODE solve (used to GENERATE data, never to fit)."""
    p = ModelParams(c0=C_TOT, K=K_ANCHOR, gamma_inf=S.gamma_inf, D=S.D,
                    n_factor=float(S.n_ions), capacity=cap, tau2=tau2,
                    temperature_K=T_REF)
    r = run_oscillation(p, omega, STRAIN, n_cycles=12, n_return_cycles=N_RETURN_CYCLES,
                        n_per_cycle=160)
    t = np.asarray(r.t)
    g = H.add_noise(t, r.gamma_mN_m, sigma, seed=seed) if sigma else np.asarray(r.gamma_mN_m)
    return H.lag_deg(t, g, r.area_rel, omega)


def build_table(force: bool = False) -> dict:
    cache = Path(__file__).with_name("blind_lag_table.npz")
    if cache.exists() and not force:
        d = np.load(cache)
        z = d["z"]
    else:
        z = np.empty((_G_OM.size, _G_CAP.size))
        for i, om in enumerate(_G_OM):
            for j, cap in enumerate(_G_CAP):
                p = ModelParams(c0=C_TOT, K=K_ANCHOR, gamma_inf=S.gamma_inf, D=S.D,
                                n_factor=float(S.n_ions), capacity=cap, tau2=om,
                                temperature_K=T_REF)
                r = run_oscillation(p, 1.0, STRAIN, n_cycles=10,
                                    n_return_cycles=6, n_per_cycle=120)
                t = np.asarray(r.t)
                # omega = 1 rad/s, so tau2 == Omega
                z[i, j] = H.lag_deg(t, r.gamma_mN_m, r.area_rel, 1.0)
            print(f"  table row {i + 1}/{_G_OM.size}", flush=True)
        np.savez(cache, z=z, om=_G_OM, cap=_G_CAP)
    return {"spline": RectBivariateSpline(np.log(_G_OM), np.log(_G_CAP), z, kx=3, ky=3)}


_T = None


def tab_lag(tau2: float, cap: float, omega: float) -> float:
    global _T
    if _T is None:
        _T = build_table()
    om = min(max(omega * tau2, _G_OM[0]), _G_OM[-1])
    c = min(max(cap, _G_CAP[0]), _G_CAP[-1])
    return float(_T["spline"](math.log(om), math.log(c), grid=False))


def fit(obs: dict, cap_fixed: float | None, start: float) -> dict:
    def res(x):
        tau2 = math.exp(x[0])
        cap = cap_fixed if cap_fixed is not None else math.exp(x[1])
        return [tab_lag(tau2, cap, om) - obs[om] for om in OMEGAS]

    x0 = [math.log(start)] if cap_fixed is not None else [math.log(start), math.log(20.0)]
    lo = [math.log(1e-2)] if cap_fixed is not None else [math.log(1e-2), math.log(0.6)]
    hi = [math.log(1e3)] if cap_fixed is not None else [math.log(1e3), math.log(500.0)]
    sol = least_squares(res, x0=x0, bounds=(lo, hi), xtol=1e-10, ftol=1e-10,
                        diff_step=1e-3, max_nfev=200)
    dof = max(len(sol.fun) - len(x0), 1)
    s2 = max(2 * sol.cost / dof, 1e-12)
    try:
        cov = np.linalg.inv(sol.jac.T @ sol.jac) * s2
        se = float(np.sqrt(np.diag(cov))[0])
    except np.linalg.LinAlgError:
        se = float("nan")
    out = {"tau2": float(math.exp(sol.x[0])), "relSE": se, "cost": float(2 * sol.cost)}
    if cap_fixed is None:
        out["capacity"] = float(math.exp(sol.x[1]))
    return out


def trial(t2: float, cap: float, seed: int, protocol: str, multi: bool) -> dict:
    obs = {om: solver_lag(t2, cap, om, SIGMA_G, seed=seed * 37 + int(om * 11))
           for om in OMEGAS}
    cap_arg = cap if protocol == "A" else None
    if multi:
        best = None
        for st in (0.02, 0.1, 0.5, 2.0, 8.0, 30.0, 100.0):
            r = fit(obs, cap_arg, st)
            if best is None or r["cost"] < best["cost"]:
                best = r
        r = best
    else:
        r = fit(obs, cap_arg, 0.3)          # generic start, as the protocol specifies
    err = r["tau2"] / t2 - 1.0
    return {"truth": t2, "cap_truth": cap, "est": r["tau2"], "relerr": err,
            "relSE": r["relSE"], "wrong": abs(err) > 0.5,
            "confident_wrong": abs(err) > 0.5 and (r["relSE"] or 9) < 0.1,
            "reported_cap": r.get("capacity")}


def main() -> dict:
    build_table()
    out = {"design": {"omegas": list(OMEGAS), "sigma_g": SIGMA_G,
                      "n_repeats": N_REPEATS, "cycles": N_RETURN_CYCLES}, "protocols": {}}
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
                "median_abserr_on_success": float(np.median(ok)) if ok else None,
                "rows": rows}
            med = f"{np.median(ok)*100:.2f}%" if ok else "n/a"
            print(f"{key:20} n={len(rows):3d} wrong={nw:3d} confident-wrong={nc:3d} "
                  f"median|err|={med}", flush=True)
    return out


if __name__ == "__main__":
    res = main()
    path = Path(__file__).with_name("blind_inversion.json")
    path.write_text(json.dumps(res, indent=2, default=float))
    print(f"\nwrote {path}")
