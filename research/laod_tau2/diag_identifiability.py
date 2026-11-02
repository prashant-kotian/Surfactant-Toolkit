"""Diagnostics for the failed kill criterion: WHERE does tau2 become identifiable?

The first recoverability run failed (relative SE up to 1e3). A failed fit can
mean one of three very different things:

  (i)   the observable carries no information about tau2 in the accessible
        frequency window  -> the idea is dead for that class of surfactant;
  (ii)  tau2 and capacity are degenerate -> fixable by estimating capacity
        from independent equilibrium data (CMC, aggregation number), leaving
        a one-parameter inversion;
  (iii) the fitting protocol is bad -> fixable.

This script separates them:

  diag_leverage()      noise-free leverage of each observable with respect to
                       tau2 across the instrument window. If |dI2/dln tau2| is
                       below the noise floor everywhere, it is case (i).
  recover_1param()     inversion with capacity FIXED at a value estimated
                       independently -- the realistic protocol. If this
                       succeeds where the 2-parameter fit failed, it is (ii).
"""

from __future__ import annotations

import math

import numpy as np

from . import harmonics as H
from .model import ModelParams, run_oscillation
from .study_recoverability import BASE, Features, measure

# Frequency window of a real oscillating-drop tensiometer (rad/s).
WINDOW = (2.0 * math.pi * 0.01, 2.0 * math.pi * 1.0)
OMEGAS = [0.05, 0.1, 0.3, 0.5, 1.0, 2.0, 3.0, 6.28]
NOISE_RATIO = 0.01     # 1% of I1 is the FT-rheology detection floor
NOISE_LAG = 0.5        # degrees


def _feat(tau2, cap, omega, eps=0.20):
    p = ModelParams(**BASE, tau2=tau2, capacity=cap)
    r = run_oscillation(p, omega, eps, n_cycles=20)
    t = np.asarray(r.t)
    hg = H.harmonics(t, r.gamma_mN_m, omega)
    return Features(omega=omega, eps=eps, I2_over_I1=hg[2] / hg[1],
                    I3_over_I1=hg[3] / hg[1], lag_deg=H.lag_deg(t, r.gamma_mN_m, r.area_rel, omega))


def diag_leverage(capacity=0.2, eps=0.20):
    """Per-omega sensitivity of each observable to tau2, vs the noise floor."""
    print("\n=== DIAGNOSTIC 1: does tau2 leave a measurable trace in the instrument window? ===")
    print("leverage = |change in observable| for a 2x change in tau2; usable if > noise floor")
    print(f"{'omega':>6} {'tau2':>7} {'dI2/I1':>8} {'dI3/I1':>8} {'dlag[deg]':>9} {'usable?':>8}")
    taus = [1e-3, 1e-2, 1e-1, 1.0, 10.0]
    rows = []
    for tau2 in taus:
        best = 0.0
        for om in OMEGAS:
            a = _feat(tau2, capacity, om, eps)
            b = _feat(2 * tau2, capacity, om, eps)
            dI2 = abs(b.I2_over_I1 - a.I2_over_I1)
            dI3 = abs(b.I3_over_I1 - a.I3_over_I1)
            dlag = abs(b.lag_deg - a.lag_deg)
            usable = (dI2 > NOISE_RATIO * 0.5) or (dI3 > NOISE_RATIO * 0.5) or (dlag > NOISE_LAG)
            best = max(best, max(dI2 / 0.01, dI3 / 0.01, dlag / NOISE_LAG))
            rows.append((tau2, om, dI2, dI3, dlag, usable))
    for tau2, om, dI2, dI3, dlag, usable in rows:
        print(f"{om:6.2f} {tau2:7.0e} {dI2:8.4f} {dI3:8.4f} {dlag:9.2f} {'YES' if usable else 'no':>8}")
    print("\nsummary: for each tau2, the largest leverage/noise ratio over the window")
    for tau2 in taus:
        sel = [max(r[2] / 0.01, r[3] / 0.01, r[4] / NOISE_LAG) for r in rows if r[0] == tau2]
        print(f"  tau2 = {tau2:6.0e} s : max leverage/noise = {max(sel):7.2f} "
              f"({'IDENTIFIABLE' if max(sel) > 2 else 'NOT identifiable'})")


def recover_1param(truth_tau2, truth_cap, grid, sigma=0.05, seed=7):
    """Invert for tau2 alone, with capacity taken from independent data."""
    from scipy.optimize import minimize_scalar

    obs = [measure(ModelParams(**BASE, tau2=truth_tau2, capacity=truth_cap), om, eps, 0.0, 0)
           for (om, eps) in grid]
    scale = np.array([0.005, 0.005, 2.0])

    def cost(logtau):
        tau = math.exp(logtau)
        res = []
        for (om, eps), o in zip(grid, obs):
            pred = measure(ModelParams(**BASE, tau2=tau, capacity=truth_cap), om, eps, 0.0, 0)
            res += [(pred.I2_over_I1 - o.I2_over_I1) / scale[0],
                    (pred.I3_over_I1 - o.I3_over_I1) / scale[1],
                    (pred.lag_deg - o.lag_deg) / scale[2]]
        return float(np.sum(np.square(res)))

    # add noise to the "observed" set instead of the predictions
    obs = [measure(ModelParams(**BASE, tau2=truth_tau2, capacity=truth_cap), om, eps, sigma, seed + i)
           for i, (om, eps) in enumerate(grid)]

    def cost_noisy(logtau):
        tau = math.exp(logtau)
        res = []
        for (om, eps), o in zip(grid, obs):
            pred = measure(ModelParams(**BASE, tau2=tau, capacity=truth_cap), om, eps, 0.0, 0)
            res += [(pred.I2_over_I1 - o.I2_over_I1) / scale[0],
                    (pred.I3_over_I1 - o.I3_over_I1) / scale[1],
                    (pred.lag_deg - o.lag_deg) / scale[2]]
        return float(np.sum(np.square(res)))

    res = minimize_scalar(cost_noisy, bounds=(math.log(1e-4), math.log(10.0)),
                          method="bounded", options={"xatol": 1e-4})
    fit = math.exp(res.x)
    clean = cost(math.log(fit))
    return {"tau2_true": truth_tau2, "tau2_fit": fit,
            "ratio": fit / truth_tau2, "chi2": clean, "n_obs": 3 * len(grid)}


def main():
    diag_leverage()
    print("\n=== DIAGNOSTIC 2: one-parameter inversion (capacity known independently) ===")
    grid = [(w, 0.20) for w in (0.1, 0.3, 1.0, 3.0, 6.28)]
    print(f"{'tau2_true':>10} {'tau2_fit':>10} {'fit/true':>9} {'chi2':>8} {'n_obs':>6}")
    for tau2 in (1e-2, 1e-1, 1.0):
        r = recover_1param(tau2, 0.2, grid)
        print(f"{r['tau2_true']:10.0e} {r['tau2_fit']:10.3e} {r['ratio']:9.2f} {r['chi2']:8.1f} {r['n_obs']:6d}")


if __name__ == "__main__":
    main()
