"""KILL-CRITERION STUDY: can tau2 actually be recovered from LAOD harmonics?

This is the cheapest test that can end the project, so it is run first.

Question
--------
Given (a) realistic tensiometer noise, (b) the frequency window a real
oscillating-drop instrument can reach, and (c) the strain amplitudes those
instruments allow, is the micellar relaxation time tau2 recoverable from the
harmonic content of gamma(t) -- and how precisely?

Design
------
1. Generate a "measurement" with the forward model at a known tau2.
2. Add realistic noise (default sigma = 0.05 mN/m, the resolution of a good
   pendant-drop tensiometer).
3. Fit (tau2, capacity) back out by weighted least squares over the harmonic
   observables (I2/I1, I3/I1, phase lag) across a frequency sweep.
4. Report the recovered tau2 with an uncertainty, and whether the fit could
   distinguish the true capacity from a wrong one.

If the recovered tau2 has a relative uncertainty of order 1 or worse, or if
the two parameters are degenerate, the method fails the kill criterion and
there is no paper in this form.

Usage
-----
    python -m research.laod_tau2.study_recoverability
"""

from __future__ import annotations

import json
import math
import sys
from dataclasses import asdict, dataclass

import numpy as np

from . import harmonics as H
from .model import ModelParams, run_oscillation

# A representative nonionic system above its CMC (C12E8-like, literature
# orders of magnitude; not a fit to any one paper).
BASE = dict(c0=1.0, K=5.0, gamma_inf=3.5e-6, D=4.0e-10, gamma0=72.0)

# Instruments: (name, highest reliable omega rad/s, achievable eps)
INSTRUMENTS = {
    "oscillating_drop_1Hz": (2.0 * math.pi * 1.0, 0.10),
    "oscillating_drop_0.1Hz": (2.0 * math.pi * 0.1, 0.30),
    "capillary_pressure_10Hz": (2.0 * math.pi * 10.0, 0.05),
}


@dataclass
class Features:
    omega: float
    eps: float
    I2_over_I1: float
    I3_over_I1: float
    lag_deg: float


def measure(p: ModelParams, omega: float, eps: float, sigma: float, seed: int) -> Features:
    """One synthetic 'measurement' at a single (omega, eps)."""
    r = run_oscillation(p, omega, eps, n_cycles=20)
    t = np.asarray(r.t)
    g = H.add_noise(t, r.gamma_mN_m, sigma, seed=seed)
    a = np.asarray(r.area_rel)
    hg = H.harmonics(t, g, omega)
    return Features(
        omega=omega,
        eps=eps,
        I2_over_I1=hg[2] / hg[1],
        I3_over_I1=hg[3] / hg[1],
        lag_deg=H.lag_deg(t, g, a, omega),
    )


def predict(tau2: float, capacity: float, grid, sigma: float = 0.0, seed: int = 0):
    out = []
    for i, (omega, eps) in enumerate(grid):
        p = ModelParams(**BASE, tau2=tau2, capacity=capacity)
        out.append(measure(p, omega, eps, sigma, seed + i))
    return out


def _residual(obs, pred):
    """Pooled residual on the three observables, scaled by plausible noise."""
    scale = np.array([0.005, 0.005, 2.0])  # absolute floors on I2/I1, I3/I1, deg
    res = []
    for o, q in zip(obs, pred):
        res.append((q.I2_over_I1 - o.I2_over_I1) / scale[0])
        res.append((q.I3_over_I1 - o.I3_over_I1) / scale[1])
        res.append((q.lag_deg - o.lag_deg) / scale[2])
    return np.asarray(res)


def recover(obs, grid, p0=(1e-2, 0.2), bounds=((1e-4, 10.0), (1e-3, 5.0))):
    """Invert the harmonic observables for (tau2, capacity)."""
    from scipy.optimize import least_squares

    def fun(x):
        tau2, cap = x
        return _residual(obs, predict(tau2, cap, grid))

    sol = least_squares(
        fun, x0=p0, bounds=([bounds[0][0], bounds[1][0]], [bounds[0][1], bounds[1][1]]),
        xtol=1e-10, ftol=1e-10, diff_step=1e-3,
    )
    dof = max(len(sol.fun) - 2, 1)
    s2 = 2.0 * sol.cost / dof
    try:
        cov = np.linalg.inv(sol.jac.T @ sol.jac) * max(s2, 1.0)
        se = np.sqrt(np.diag(cov))
        corr = cov[0, 1] / math.sqrt(cov[0, 0] * cov[1, 1])
    except np.linalg.LinAlgError:
        se = np.array([np.nan, np.nan])
        corr = np.nan
    return {
        "tau2": float(sol.x[0]),
        "capacity": float(sol.x[1]),
        "tau2_rel_se": float(se[0] / sol.x[0]) if sol.x[0] else float("inf"),
        "capacity_rel_se": float(se[1] / sol.x[1]) if sol.x[1] else float("inf"),
        "corr_tau2_capacity": float(corr),
        "cost": float(sol.cost),
    }


def main() -> int:
    sigma = 0.05
    grid_1hz = [(w, 0.10) for w in (0.05, 0.1, 0.3, 1.0, 2.0, 2.0 * math.pi * 1.0)]
    grid_deep = [(w, 0.30) for w in (0.05, 0.1, 0.3, 1.0, 2.0, 2.0 * math.pi * 1.0)]

    print("=" * 78)
    print("KILL-CRITERION STUDY: tau2 recovery from LAOD harmonics")
    print(f"noise sigma = {sigma} mN/m (realistic pendant-drop resolution)")
    print("=" * 78)

    truth = [(1e-2, 0.2), (1e-1, 0.2), (1.0, 0.2)]
    results = {}
    for grid, label in ((grid_1hz, "eps=0.10 (1 Hz instrument, safe amplitude)"),
                        (grid_deep, "eps=0.30 (low-frequency instrument, deep strain)")):
        print(f"\n--- observation window: {label} ---")
        print(f"{'tau2_true':>10} {'tau2_fit':>10} {'rel_se':>8} {'cap_fit':>8} {'corr':>7}")
        for tau2, cap in truth:
            obs = predict(tau2, cap, grid, sigma=sigma, seed=7)
            fit = recover(obs, grid)
            results[f"{label}|{tau2}"] = fit
            print(
                f"{tau2:10.0e} {fit['tau2']:10.3e} {fit['tau2_rel_se']:8.2f} "
                f"{fit['capacity']:8.3f} {fit['corr_tau2_capacity']:7.2f}"
            )

    with open("research/laod_tau2/recoverability_results.json", "w") as f:
        json.dump(results, f, indent=2)
    print("\nwrote research/laod_tau2/recoverability_results.json")

    worst = max(r["tau2_rel_se"] for r in results.values() if r["tau2"] == r["tau2"])
    ok = worst < 0.5
    print(f"\nKILL CRITERION (worst relative SE < 0.5): {'PASS' if ok else 'FAIL'}  (worst = {worst:.2f})")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
