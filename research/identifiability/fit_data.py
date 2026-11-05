"""Fit the mechanism models to published C14E8 phase data; report what the
data can and cannot identify.

Data: 6 points of |delta| versus frequency, digitised from Fainerman et al.
(Langmuir 2008, Fig. 5, 30 umol/L C14E8), estimated 1 deg precision. Six
points is typical for oscillatory dilatational work, so the identifiability
answer obtained here is a statement about the field's information content,
not about an idealised experiment.

Three structural facts drive the design:

1. E0 does not appear in arg(E). In
       E = E0 / (1 + sqrt(B) exp(-i pi/4) / u)
   E0 is a real multiplicative factor, so the phase is completely blind to it.
   Phase-only data therefore cannot identify the Gibbs elasticity at all --
   yet phase-only papers report "E". Here E0 is excluded from the parameter
   vector and the consequence is reported rather than papered over.

2. h (the isotherm slope) and D enter only through u = h sqrt(w/D), and the
   observable depends on them through a single shape parameter. They cannot be
   separated from each other by phase data either. h is therefore treated as a
   fixed but *uncertain* nuisance and the whole fit is repeated across its
   plausible range. Reporting one number would be false precision.

3. Model comparison and parameter identifiability are different questions and
   are answered separately: AICc for the first, profile likelihood for the
   second. A model can win on AICc while its parameters remain one-sided or
   flat, and that combination is the situation of interest.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from .models import B_single, B_two_mode, B_log_normal, lag_deg

CHI2_95 = 3.841458820694124

# --- medium -----------------------------------------------------------------
D_M2_S = 3.8e-10
# Calibrated so the B = 1 model reproduces the project's validated LVDT
# baseline (|delta| = 40.08 deg at 0.1 Hz), i.e. h = gamma_inf K / (1+K c0)^2
# with gamma_inf = 2.8e-6 mol/m2, K = 140 m3/mol, c0 = 0.071 mol/m3.
H_NOMINAL = 3.275302547717482e-06
H_RANGE = (H_NOMINAL / 3.0, H_NOMINAL * 3.0)


def load_data(path="research/laod_tau2/fainerman_fit.json"):
    d = json.loads(Path(path).read_text())["data"]
    return (np.array(d["freqs"], float), np.array(d["phi_obs"], float),
            float(d["phi_err"]))


# --- models: parameter vectors contain NO E0 -------------------------------
#    M1  B = 1                                  k = 0
#    M2  B = 1 + Cm/(1+iw tau)                  k = 2
#    M3  B = 1 + Cm[w/(1+iwt1)+(1-w)/(1+iwt2)]  k = 3 (w fixed at 0.5)
#    M4  B = log-normal distribution of tau     k = 3

E0_PLACEHOLDER = 60.0    # mN/m; irrelevant to the phase, kept for |E| use


def predict(kind, theta, omega, D, h):
    out = []
    for w in np.atleast_1d(omega):
        if kind == "M1_lvdt":
            B = 1.0 + 0j
        elif kind == "M2_single":
            Cm, tau = theta
            B = B_single(w, Cm, tau)
        elif kind == "M3_two_mode":
            Cm, t1, t2 = theta
            B = B_two_mode(w, Cm, t1, t2, weight=0.5)
        elif kind == "M4_log_normal":
            Cm, tc, wd = theta
            B = B_log_normal(w, Cm, tc, wd)
        else:
            raise ValueError(kind)
        out.append(lag_deg(w, E0_PLACEHOLDER, D, h, complex(B)))
    return np.array(out)


# bounds in INTERNAL coordinates: Cm linear, time constants in log10
SPEC = {
    "M1_lvdt": dict(lb=np.array([]), ub=np.array([]), log=()),
    "M2_single": dict(lb=np.array([0.02, -4.0]), ub=np.array([60.0, 3.0]),
                      log=(1,)),
    "M3_two_mode": dict(lb=np.array([0.02, -4.0, -4.0]),
                        ub=np.array([60.0, 3.0, 3.0]), log=(1, 2)),
    "M4_log_normal": dict(lb=np.array([0.02, -4.0, 0.15]),
                          ub=np.array([60.0, 3.0, 3.0]), log=(1,)),
}
NK = {"M1_lvdt": 0, "M2_single": 2, "M3_two_mode": 3, "M4_log_normal": 3}


def to_natural(kind, p):
    lg = SPEC[kind]["log"]
    return np.array([10.0 ** p[i] if i in lg else p[i] for i in range(len(p))])


def _res(kind, p, omega, phi, sigma, D, h):
    return (predict(kind, to_natural(kind, p), omega, D, h) - phi) / sigma


def fit(kind, omega, phi, sigma, D=D_M2_S, h=H_NOMINAL, n_starts=30, seed=7):
    lo, hi = SPEC[kind]["lb"], SPEC[kind]["ub"]
    if len(lo) == 0:
        return float(np.sum(((predict(kind, np.array([]), omega, D, h) - phi) / sigma) ** 2)), np.array([])
    rng = np.random.default_rng(seed)
    starts = [np.clip(np.array([1.0, math.log10(1.0), math.log10(10.0)])[:len(lo)], lo, hi)]
    for _ in range(n_starts):
        starts.append(lo + rng.random(len(lo)) * (hi - lo))
    best = None
    for s in starts:
        try:
            r = least_squares(lambda pp: _res(kind, pp, omega, phi, sigma, D, h),
                              np.clip(s, lo, hi), bounds=(lo, hi), method="trf",
                              xtol=1e-13, ftol=1e-13, gtol=1e-13, max_nfev=5000)
        except Exception:
            continue
        c = float(np.sum(r.fun ** 2))
        if best is None or c < best[0]:
            best = (c, to_natural(kind, r.x))
    return best


def profile(kind, idx, grid, omega, phi, sigma, D=D_M2_S, h=H_NOMINAL,
            n_starts=10, seed=11):
    lo, hi = SPEC[kind]["lb"], SPEC[kind]["ub"]
    lg = SPEC[kind]["log"]
    keep = [i for i in range(len(lo)) if i != idx]
    rng = np.random.default_rng(seed)
    out = []
    for g in grid:
        gp = math.log10(g) if idx in lg else g
        if gp < lo[idx] or gp > hi[idx]:
            out.append(np.inf)
            continue
        if not keep:
            out.append(float(np.sum(_res(kind, np.array([gp]), omega, phi, sigma, D, h) ** 2)))
            continue

        def r_rest(pk, gp=gp):
            p = np.empty(len(lo))
            p[idx] = gp
            p[keep] = pk
            return _res(kind, p, omega, phi, sigma, D, h)

        lo_k, hi_k = lo[keep], hi[keep]
        best = None
        for _ in range(n_starts):
            s = lo_k + rng.random(len(lo_k)) * (hi_k - lo_k)
            try:
                r = least_squares(r_rest, s, bounds=(lo_k, hi_k), method="trf",
                                  xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=3000)
            except Exception:
                continue
            c = float(np.sum(r.fun ** 2))
            if best is None or c < best:
                best = c
        out.append(best if best is not None else np.inf)
    return np.array(out)


def classify(grid, chi2, thresh=CHI2_95):
    d = chi2 - np.nanmin(chi2)
    inside = grid[d <= thresh]
    if len(inside) == 0:
        return dict(shape="narrow", interval=None, lower_open=False, upper_open=False)
    lo_open, hi_open = bool(d[0] <= thresh), bool(d[-1] <= thresh)
    shape = ("flat" if lo_open and hi_open
             else "one_sided" if lo_open or hi_open else "bounded")
    return dict(shape=shape, interval=(float(inside.min()), float(inside.max())),
                lower_open=lo_open, upper_open=hi_open)


def aicc(chi2, n, k):
    if n - k - 1 <= 0:
        return np.inf
    return n * math.log(max(chi2, 1e-300) / n) + 2 * k + 2 * k * (k + 1) / (n - k - 1)


def compare(omega, phi, sigma, D=D_M2_S, h=H_NOMINAL):
    out = {}
    for kind in SPEC:
        b = fit(kind, omega, phi, sigma, D, h)
        if b is None:
            out[kind] = None
            continue
        chi2, th = b
        k = NK[kind]
        out[kind] = dict(theta=[float(x) for x in th], chi2=chi2, k=k,
                         aicc=aicc(chi2, len(phi), k),
                         rms_deg=math.sqrt(chi2 / len(phi)) * sigma)
    ok = [k for k in out if out[k] is not None]
    if not ok:
        return out, None, {}
    w = {}
    if len(ok) == len(SPEC):
        a = np.array([out[k]["aicc"] for k in SPEC])
        a = a - a.min()
        wv = np.exp(-0.5 * a)
        wv /= wv.sum()
        w = {k: float(x) for k, x in zip(SPEC, wv)}
    return out, min(ok, key=lambda x: out[x]["aicc"]), w


def run(save=True):
    f, phi, sigma = load_data()
    omega = 2 * np.pi * f
    fits, best, w = compare(omega, phi, sigma)

    # --- nuisance sensitivity: repeat the whole comparison over h ----------
    sens = {}
    for hf in (0.33, 0.5, 1.0, 2.0, 3.0):
        h = H_NOMINAL * hf
        _f, _b, _w = compare(omega, phi, sigma, h=h)
        sens[str(hf)] = dict(best=_b,
                             chi2={k: (_f[k]["chi2"] if _f[k] else None) for k in SPEC},
                             weights=_w)

    prof = {}
    for kind, idx, tag in (("M2_single", 1, "M2_tau"),):
        grid = np.logspace(-4, 3, 45)
        c = profile(kind, idx, grid, omega, phi, sigma)
        prof[tag] = dict(grid=grid.tolist(), chi2=c.tolist(), **classify(grid, c))
    # Cm is the parameter that bounds the lag; profile it too
    grid_cm = np.logspace(-2, 2, 41)
    c = profile("M2_single", 0, grid_cm, omega, phi, sigma)
    prof["M2_Cm"] = dict(grid=grid_cm.tolist(), chi2=c.tolist(), **classify(grid_cm, c))

    out = dict(data=dict(freqs=f.tolist(), phase_obs=phi.tolist(), sigma=sigma,
                         n=len(f), D=D_M2_S, h=H_NOMINAL,
                         source="Fainerman 2008 Fig 5, 30 uM C14E8"),
               fits=fits, best=best, akaike_weights=w,
               nuisance_sensitivity_h=sens, profiles=prof,
               note=("E0 is structurally unidentifiable from phase alone; it does not "
                     "enter arg(E). Parameters reported are Cm and time constants only."))
    if save:
        Path("research/identifiability/fit_results.json").write_text(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    r = run()
    print(f"n = {r['data']['n']}, sigma = {r['data']['sigma']} deg, "
          f"h = {r['data']['h']:.3e} m")
    print(f"\n{'model':16s}{'k':>3s}{'chi2':>9s}{'RMS deg':>9s}{'AICc':>9s}  params")
    for k, v in r["fits"].items():
        if v is None:
            print(f"{k:16s} FAILED"); continue
        print(f"{k:16s}{v['k']:3d}{v['chi2']:9.2f}{v['rms_deg']:9.2f}{v['aicc']:9.2f}  "
              f"{[round(x, 3) for x in v['theta']]}")
    print("\nbest:", r["best"])
    print("Akaike weights:", {k: round(v, 3) for k, v in r["akaike_weights"].items()})
    print("\nnuisance sensitivity to h (h/h_nominal -> best model, chi2 by model):")
    for k, v in r["nuisance_sensitivity_h"].items():
        cs = {m: (round(c, 1) if c is not None else None) for m, c in v["chi2"].items()}
        print(f"  {k:>5s}  best={v['best']:14s} {cs}")
    print("\nprofiles:")
    for k, v in r["profiles"].items():
        print(f"  {k:9s} shape={v['shape']:10s} interval={v['interval']}")
