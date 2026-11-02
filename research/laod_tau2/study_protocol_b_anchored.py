"""Protocol B at a physically anchored isotherm (closes the last gap in FINDINGS 13).

Protocol B was verified with K = 140 m^3/mol, which section 13.5 identified as the
saturated regime (theta = 0.937, E0 = 103.5 mN/m). The question this answers:
does a concentration series still break the (tau2, capacity) degeneracy when the
isotherm is anchored to a literature E0 of 25 mN/m?

Method is unchanged from the original Protocol B run: synthetic observations at
the model's own phase lag for a known tau2 and a known capacity from the closure,
then a joint two-parameter least squares on (tau2, kappa) where kappa is a scale
on n_bar. Correlation between the two parameters is the quantity of interest --
it is the degeneracy. The answer is compared across anchors for C12E8, and then
for CTAB, to see whether the earlier finding ("two concentrations suffice") is
an artefact of the saturated isotherm.

Outputs research/laod_tau2/protocol_b_anchored.json.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from . import harmonics as H
from .model import ModelParams, run_oscillation
from .systems import SYSTEMS

R_GAS = 8.314462618
EPS = 0.20
TAU2 = 4.0
OMEGAS = (3.0, 5.0, 6.28)
SERIES = [(1.5,), (1.3, 1.5), (1.3, 1.5, 2.0), (1.2, 1.5, 2.0, 3.0)]
LAG_SD = 0.18          # deg, from the sampling study (FINDINGS 11)


def K_for_E0(name: str, E0: float, a: float = 0.0) -> float:
    """Langmuir/Frumkin K so that E0 is reached at 1.5 x CMC."""
    s = SYSTEMS[name]
    pref = s.n_ions * R_GAS * 298.15 * s.gamma_inf * 1e3
    c = 1.5 * s.cmc_mM
    # solve E0 = pref*(th/(1-th) - 2 a th^2)
    lo, hi = 1e-4, 1.0 - 1e-6

    def f(th):
        return pref * (th / (1.0 - th) - 2.0 * a * th * th) - E0

    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    th = 0.5 * (lo + hi)
    return th / (1.0 - th) * math.exp(-2.0 * a * th) / c


def lag(name: str, K: float, a: float, ct: float, tau2: float, cap: float, om: float) -> float:
    s = SYSTEMS[name]
    p = ModelParams(c0=ct, K=K, gamma_inf=s.gamma_inf, D=s.D, n_factor=float(s.n_ions),
                    capacity=cap, tau2=tau2, frumkin_a=a, temperature_K=298.15)
    r = run_oscillation(p, om, EPS, n_cycles=16)
    t = np.asarray(r.t)
    return H.lag_deg(t, r.gamma_mN_m, r.area_rel, om)


def joint(name: str, K: float, a: float, mults, tau2_true: float):
    s = SYSTEMS[name]
    obs = {(m, om): lag(name, K, a, m * s.cmc_mM, tau2_true,
                        s.capacity(m * s.cmc_mM), om)
           for m in mults for om in OMEGAS}

    def res(x):
        tau2 = math.exp(x[0])
        kappa = math.exp(x[1])
        out = []
        for m in mults:
            cap = (m - 1.0) * s.n_bar * kappa
            for om in OMEGAS:
                out.append(lag(name, K, a, m * s.cmc_mM, tau2, cap, om) - obs[(m, om)])
        return out

    hi_tau = max(4.0 * tau2_true, 0.2)
    sol = least_squares(res, x0=[math.log(tau2_true), math.log(1.0)],
                        bounds=([math.log(tau2_true / 4), math.log(0.05)],
                                [math.log(hi_tau * 2), math.log(5.0)]),
                        xtol=1e-9, ftol=1e-9, diff_step=1e-3)
    dof = max(len(res(sol.x)) - 2, 1)
    s2 = max(2 * sol.cost / dof, LAG_SD ** 2)
    try:
        cov = np.linalg.inv(sol.jac.T @ sol.jac) * s2
        se = np.sqrt(np.diag(cov))
        corr = cov[0, 1] / math.sqrt(cov[0, 0] * cov[1, 1])
    except np.linalg.LinAlgError:
        return None
    return {"tau2": math.exp(sol.x[0]), "relSE": float(se[0]), "corr": float(corr),
            "kappa": math.exp(sol.x[1])}


def main():
    out = {"tau2": TAU2, "eps": EPS, "lag_sd": LAG_SD, "cases": []}

    # C12E8 across E0 anchors -- does the degeneracy-breaking survive?
    for E0 in (25.0, 45.3, 103.5):
        K = K_for_E0("C12E8", E0)
        for mults in SERIES:
            r = joint("C12E8", K, 0.0, mults, TAU2)
            out["cases"].append({"system": "C12E8", "E0": E0, "K": K,
                                 "series": list(mults), "result": r})

    # CTAB at its anchored E0, if the instrument could reach the frequency
    Kc = K_for_E0("CTAB", 30.0, 0.7)
    for mults in SERIES:
        r = joint("CTAB", Kc, 0.7, mults, 0.065)
        out["cases"].append({"system": "CTAB", "E0": 30.0, "K": Kc,
                             "series": list(mults), "result": r})
    return out


if __name__ == "__main__":
    res = main()
    path = Path(__file__).with_name("protocol_b_anchored.json")
    path.write_text(json.dumps(res, indent=2, default=float))
    print(f"Protocol B at anchored isotherms (tau2 = {res['tau2']} s for C12E8)\n")
    print(f"{'system':7} {'E0':>6} {'series':>20} {'tau2_fit':>9} {'relSE':>7} {'corr':>7}")
    for c in res["cases"]:
        r = c["result"]
        tag = ", ".join(f"{m}x" for m in c["series"])
        if r is None:
            print(f"{c['system']:7} {c['E0']:6.1f} {tag:>20} {'SINGULAR':>9}")
        else:
            print(f"{c['system']:7} {c['E0']:6.1f} {tag:>20} {r['tau2']:9.3f} "
                  f"{r['relSE']:7.3f} {r['corr']:7.3f}")
    print(f"\nwrote {path}")
