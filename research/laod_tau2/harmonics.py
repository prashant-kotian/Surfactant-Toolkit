"""Demodulation of the periodic LAOD response, plus noise and convergence tooling.

Everything downstream of the forward model depends on getting the harmonic
content right, so this module is deliberately paranoid:

  * harmonics() fits a truncated Fourier series by linear least squares over
    an integer number of cycles (no spectral leakage; no window needed).
  * convergence() checks the reported ratios are stable under model
    refinement (grid, step, cycle count) -- a harmonic ratio that moves with
    resolution is a numerical artefact, not physics.
  * add_noise() injects realistic tensiometer noise so that a sensitivity
    study can answer "could this be measured?" rather than "does the model
    produce a number?".

Conventions: I1, I2, I3 are AMPLITUDES (not peak-to-peak); ratios are
I2/I1 and I3/I1, which is the FT-rheology convention (a 1% ratio is the
usual detection floor).
"""

from __future__ import annotations

import math

import numpy as np

DEFAULT_HARMONICS = (1, 2, 3, 4)


def _design(t: np.ndarray, omega: float, kmax: int) -> np.ndarray:
    th = omega * (t - t[0])
    cols = [np.ones_like(th)]
    for k in range(1, kmax + 1):
        cols += [np.cos(k * th), np.sin(k * th)]
    return np.vstack(cols).T


def harmonics(t, y, omega: float, kmax: int = 4) -> dict:
    """Least-squares Fourier decomposition: {0: DC, 1: I1, 2: I2, ...}."""
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    n_cycles = omega * (t[-1] - t[0]) / (2.0 * math.pi)
    if abs(n_cycles - round(n_cycles)) > 0.05:
        raise ValueError(
            f"harmonic fit needs an integer number of cycles; got {n_cycles:.3f}"
        )
    coef, *_ = np.linalg.lstsq(_design(t, omega, kmax), y, rcond=None)
    out = {0: float(coef[0])}
    for k in range(1, kmax + 1):
        out[k] = float(math.hypot(coef[2 * k - 1], coef[2 * k]))
    return out


def first_harmonic_phase(t, y, omega: float) -> float:
    """Phase of the fundamental, radians, in [0, 2*pi).

    y(t) = A*cos(theta - phi) with theta = omega*(t - t0); phi is returned.
    """
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    coef, *_ = np.linalg.lstsq(_design(t, omega, 1)[:, :3], y, rcond=None)
    return float(math.atan2(coef[2], coef[1]))


def lag_deg(t, gamma, area, omega: float) -> float:
    """Lag of gamma behind the area oscillation, degrees (positive = lag).

    The White Paper's prototype never reports this quantity, but it is the
    observable that disentangles a micellar reservoir from EOS nonlinearity
    (see docs), so it is measured everywhere downstream.
    """
    dphi = first_harmonic_phase(t, area, omega) - first_harmonic_phase(t, gamma, omega)
    return (math.degrees(dphi) + 360.0) % 360.0


def add_noise(t, y, sigma: float, seed: int = 0):
    """Additive white measurement noise (mN/m)."""
    rng = np.random.default_rng(seed)
    return np.asarray(y, dtype=float) + rng.normal(0.0, sigma, size=len(t))


def convergence(runner, *, n_cycles_list=(12, 18, 24), n_cells_list=(64, 96, 144)) -> dict:
    """Sweep resolution and report the spread in I2/I1 and I3/I1.

    `runner` is a callable (n_cycles, n_cells) -> (t, gamma, area, omega).
    A ratio whose spread across the sweep exceeds a few percent of its value
    is not converged and must not be quoted.
    """
    rows = []
    for nc in n_cycles_list:
        for ncel in n_cells_list:
            t, gamma, area, omega = runner(nc, ncel)
            H = harmonics(t, gamma, omega)
            rows.append(
                {
                    "n_cycles": nc,
                    "n_cells": ncel,
                    "I2_over_I1": H[2] / H[1],
                    "I3_over_I1": H[3] / H[1],
                    "lag_deg": lag_deg(t, gamma, area, omega),
                }
            )
    ref = rows[-1]
    worst = max(abs(r["I2_over_I1"] - ref["I2_over_I1"]) / max(ref["I2_over_I1"], 1e-30) for r in rows)
    return {"rows": rows, "worst_relative_spread_I2": worst}
