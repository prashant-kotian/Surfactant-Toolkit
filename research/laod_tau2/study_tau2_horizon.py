"""Where does the inverse-tau2 method actually stop working?

study_ionic_generality.py says CTAB (tau2 = 65 ms) and SDS (38 ms) fail on a
PAT-1M: the second harmonic never clears the noise floor inside the 0.005-1 Hz
band, and the frequencies where it would clear it are far outside the band and
leave fewer than 5 frames per cycle at 25 fps.

Two questions before that is believed:

  1. Is the failure an artefact of the assumed limiting elasticity? We fix
     E0 = 30 mN/m (CTAB), which fixes theta via the isotherm. Sweep the whole
     literature range 20-40 mN/m AND the Frumkin interaction parameter, and see
     whether any admissible combination clears the noise floor in-band.
  2. Is the dividing line tau2 itself? Sweep tau2 continuously on the C12E8
     isotherm to find the threshold at which the method becomes operable.

Outputs research/laod_tau2/tau2_horizon.json.
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
STRAIN = 0.20
T_REF = 298.15
NOISE = NOISE_FLOOR_RATIO[1]      # 0.05, the pessimistic end

# Operability thresholds, fixed BEFORE looking at the sweep.
I2_MIN = NOISE                     # nonlinearity must be measurable
DLAG_MIN = 3.0                     # deg of lag bought by a factor of FOUR in
                                   # tau2 (tau2/2 -> tau2*2). Below this the
                                   # inversion is noise-dominated. Read the
                                   # per-doubling figure as ~DLAG/2 landing in
                                   # the same ballpark to within a factor <2:
                                   # the response is NOT linear in log(tau2),
                                   # it saturates as tau2 grows (see FINDINGS 13).
FRAMES_MIN = 5.0                   # FINDINGS section 11
FPS = 25.0                         # PAT-1M oscillatory frame rate


def theta_for_E0(E0, gamma_inf, n_factor, a, t_ref=T_REF):
    pref = n_factor * R_GAS * t_ref * gamma_inf * 1e3

    def f(th):
        return pref * (th / (1.0 - th) - 2.0 * a * th * th) - E0

    lo, hi = 1e-4, 1.0 - 1e-6
    if f(lo) > 0 or f(hi) < 0:
        return None
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def params_from_E0(name, mult, tau2, E0, a=None):
    s = SYSTEMS[name]
    a = s.isotherm_a if a is None else a
    c = mult * s.cmc_mM
    th = theta_for_E0(E0, s.gamma_inf, float(s.n_ions), a)
    if th is None:
        return None
    K = th / (1.0 - th) * math.exp(-2.0 * a * th) / c
    return ModelParams(c0=c, K=K, gamma_inf=s.gamma_inf, D=s.D,
                       n_factor=float(s.n_ions), capacity=s.capacity(c),
                       tau2=tau2, frumkin_a=a, temperature_K=T_REF)


def probe(p, omega, eps=STRAIN):
    r = run_oscillation(p, omega, eps, n_cycles=16)
    t = np.asarray(r.t)
    hg = H.harmonics(t, r.gamma_mN_m, omega)
    return hg[2], H.lag_deg(t, r.gamma_mN_m, r.area_rel, omega)


def indicator(p, omega):
    """I2 and degrees-of-lag bought by a factor of 4 in tau2."""
    i2, lag = probe(p, omega)
    _, lo = probe(ModelParams(**{**p.__dict__, "tau2": p.tau2 / 2.0}), omega)
    _, hi = probe(ModelParams(**{**p.__dict__, "tau2": p.tau2 * 2.0}), omega)
    fpc = FPS / (omega / (2 * math.pi))
    ok = (i2 > I2_MIN) and (abs(hi - lo) >= DLAG_MIN) and (fpc >= FRAMES_MIN)
    return {"I2": i2, "dlag_per_4x": abs(hi - lo), "frames_per_cycle": fpc,
            "usable": ok}


def main():
    out = {"strain": STRAIN, "noise": NOISE, "i2_min": I2_MIN,
           "dlag_min": DLAG_MIN, "frames_min": FRAMES_MIN}

    # ---- 1. Is CTAB's failure an artefact of the assumed E0? ----------------
    # In-band only: omega = 3.14 (0.5 Hz) and 6.28 (1 Hz), the top of the PAT band.
    art = {}
    for name, tau2, e0s in (("CTAB", 0.065, (20.0, 25.0, 30.0, 35.0, 40.0)),
                            ("SDS", 0.038, (15.0, 22.0, 30.0))):
        rows = []
        for e0 in e0s:
            for a in (0.0, 0.35, 0.7, 1.0):
                for omega in (3.14, 6.28):
                    p = params_from_E0(name, 1.5, tau2, e0, a)
                    if p is None:
                        rows.append({"E0": e0, "a": a, "omega": omega,
                                     "feasible": False})
                        continue
                    d = indicator(p, omega)
                    rows.append({"E0": e0, "a": a, "omega": omega, "feasible": True,
                                 "theta": p.theta0, **d})
        feas = [r for r in rows if r.get("feasible")]
        best = max(feas, key=lambda r: r["I2"])
        art[name] = {"rows": rows, "n_feasible": len(feas),
                     "best_I2": best["I2"], "best_case": best,
                     "any_usable": any(r["usable"] for r in feas)}
    out["isotherm_artefact_check"] = art

    # ---- 2. The tau2 horizon on the C12E8 isotherm --------------------------
    # Scan tau2 from 10 ms to 8 s at the measured optimum (1.5 x CMC), at two
    # frequencies: the top of the PAT band and half of it.
    horizon = {}
    for omega in (3.14, 6.28):
        rows = []
        for tau2 in (0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 4.0, 8.0):
            p = params_from_E0("C12E8", 1.5, tau2, 30.0, 0.0)
            if p is None:
                continue
            rows.append({"tau2": tau2, **indicator(p, omega)})
        usable = [r["tau2"] for r in rows if r["usable"]]
        horizon[str(omega)] = {
            "rows": rows,
            "tau2_min_usable": min(usable) if usable else None,
        }
    out["tau2_horizon"] = horizon
    return out


if __name__ == "__main__":
    res = main()
    path = Path(__file__).with_name("tau2_horizon.json")
    path.write_text(json.dumps(res, indent=2, default=float))

    print("=== 1. Does CTAB/SDS fail because of the assumed E0? (in-band only) ===")
    for name, d in res["isotherm_artefact_check"].items():
        b = d["best_case"]
        print(f"{name}: {d['n_feasible']} feasible (E0, a, omega) combos; "
              f"best I2 = {d['best_I2']:.4f} at E0={b['E0']}, a={b['a']}, "
              f"omega={b['omega']}; any usable = {d['any_usable']}")
        print(f"   at that point theta = {b['theta']:.3f}, "
              f"dlag/4x = {b['dlag_per_4x']:.2f} deg")

    print("\n=== 2. tau2 horizon (C12E8 isotherm, 1.5x CMC, E0 = 30) ===")
    for om, d in res["tau2_horizon"].items():
        f = float(om) / (2 * math.pi)
        print(f"-- omega = {om} ({f:.2f} Hz) -- tau2_min_usable = {d['tau2_min_usable']}")
        print(f"{'tau2(s)':>9} {'I2':>9} {'dlag/4x':>8} {'fpc':>6}  usable")
        for r in d["rows"]:
            print(f"{r['tau2']:9.3f} {r['I2']:9.4f} {r['dlag_per_4x']:8.2f} "
                  f"{r['frames_per_cycle']:6.1f}  {'YES' if r['usable'] else ''}")
    print(f"\nwrote {path}")
