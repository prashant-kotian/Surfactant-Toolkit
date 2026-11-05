"""EXPLORATORY PROBE -- not a result.

Question: given a realistic frequency window and phase-noise level, WHICH
relaxation mechanisms can dilatational rheology actually distinguish?

The literature assigns mechanisms routinely ("reorientation is negligible",
"the dominant mechanism is diffusion") usually by fitting one model and
inspecting residuals. The prior question -- is the assignment identifiable at
all from this observable, in this band, at this noise? -- is rarely asked.

Three mechanisms, all of which shift phase in the same direction:
  M1  pure diffusion (Lucassen-van den Tempel)
  M2  diffusion + finite adsorption kinetics (single barrier time tau_a)
  M3  diffusion + slow surface reorientation (extra relaxation tau_r)

Method: generate synthetic data from one mechanism, fit the others by
1-parameter scan, and report chi-square. If a wrong mechanism fits within noise,
the assignment is not identifiable from that observable.

Purpose: to check the proposed paper has a real question with a real answer
BEFORE committing. No numbers here belong in a manuscript.
"""

from __future__ import annotations

import math

import numpy as np

R = 8.314462618
T = 298.15
D = 3.8e-10          # m^2/s
GINF = 2.8e-6        # mol/m^2
C0 = 0.071           # mol/m^3
E0_TARGET = 40.0     # mN/m, anchor so coverage is realistic


def h_from_E0(E0=E0_TARGET):
    """Isotherm slope that gives this Gibbs elasticity at c0."""
    return E0 / (R * T * 1e3 * C0)


H = h_from_E0()

FREQS = np.array([0.005, 0.007, 0.01, 0.014, 0.02, 0.03, 0.05,
                  0.07, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 1.0])
OM = 2 * math.pi * FREQS
SIG_E = 0.02         # relative, on |E|
SIG_PHI = 1.0        # degrees


def E_complex(omega, *, tau_a=None, tau_r=None, E0=None):
    """Complex modulus.

    Base: LVDT,  E = E0 * i w / (i w + D q / h),  q = sqrt(i w / D).

    tau_a (adsorption barrier): the flux cannot keep up with the isotherm, so
    the effective h is reduced by 1/(1 + i w tau_a) -- a serial resistance
    between subsurface and surface.

    tau_r (surface reorientation): the surface excess relaxes internally with
    time tau_r, so the elastic response is divided by (1 + i w tau_r).
    """
    q = np.sqrt(1j * omega / D)
    if q.real < 0:
        q = -q
    h_eff = H
    if tau_a is not None:
        h_eff = H / (1.0 + 1j * omega * tau_a)
    e0 = E0_TARGET if E0 is None else E0
    E = e0 * 1j * omega / (1j * omega + D * q / h_eff)
    if tau_r is not None:
        E = E / (1.0 + 1j * omega * tau_r)
    return E


def obs(omega, **kw):
    E = E_complex(omega, **kw)
    return abs(E), -math.degrees(float(np.angle(E)))


def chi2(truth, fit_kw, free_E0=False, E0_fit=None):
    """Weighted chi-square of the fitted model against synthetic truth.

    IMPORTANT: the fitted model must be a GENUINELY DIFFERENT model, and any
    parameter that would not be known in practice (E0, hence h) must be free.
    An earlier version of this script scanned the wrong mechanism's parameter
    over a grid that contained the true value, so it reported chi2 = 0 and
    "indistinguishable" for every case. That is the same design error that
    produced the retracted conclusions in FINDINGS 27.
    """
    r = []
    for om in OM:
        e_t, p_t = obs(om, **truth)
        e_f, p_f = obs(om, **(fit_kw or {}), E0=E0_fit if free_E0 else None)
        r.append((e_f - e_t) / (SIG_E * e_t))
        r.append((p_f - p_t) / SIG_PHI)
    return float(np.sum(np.square(r)))


def best_with_free_E0(truth, key, grid, e0_grid=None):
    """Fit the WRONG mechanism, profiling over its own E0.

    Returns (chi2, tau, E0). A wrong mechanism that reaches low chi2 after
    re-scaling E0 is not identifiable from this observable.
    """
    if e0_grid is None:
        e0_grid = E0_TARGET * np.logspace(-1.0, 1.0, 40)
    out = []
    for e0 in e0_grid:
        for g in grid:
            out.append((chi2(truth, {key: g}, free_E0=True, E0_fit=e0), g, e0))
    out.sort()
    return out[0]


if __name__ == "__main__":
    n_points = 2 * len(FREQS)
    reject = n_points + 3 * math.sqrt(2 * n_points)
    print(f"E0 = {E0_TARGET} mN/m, window {FREQS[0]}-{FREQS[-1]} Hz, "
          f"{len(FREQS)} freqs, {n_points} observables")
    print(f"noise: {SIG_E:.0%} on |E|, {SIG_PHI} deg on phase")
    print(f"chi2 ~ {n_points} = good fit; chi2 > {reject:.0f} = wrong mechanism excluded\n")

    grid_a = np.logspace(-4, 0, 40)
    grid_r = np.logspace(-3, 1, 40)

    print("CAN PURE DIFFUSION (M1, E0 free) FAKE A MORE COMPLEX MECHANISM?")
    for name, truth, key, grid in [
        ("barrier tau_a=10 ms", {"tau_a": 1e-2}, "tau_a", grid_a),
        ("barrier tau_a=100 ms", {"tau_a": 1e-1}, "tau_a", grid_a),
        ("reorient tau_r=100 ms", {"tau_r": 1e-1}, "tau_r", grid_r),
        ("reorient tau_r=1 s", {"tau_r": 1.0}, "tau_r", grid_r),
    ]:
        c, v, e0 = best_with_free_E0(truth, key, grid)
        verdict = "CANNOT be distinguished" if c < reject else "distinguished"
        print(f"  data={name:22s} fit=M1: chi2={c:8.1f} at E0={e0:5.1f}  -> {verdict}")

    print("\nCAN A BARRIER FAKE REORIENTATION (both against M1-dominated data)?")
    c, v, e0 = best_with_free_E0({"tau_r": 1e-1}, "tau_a", grid_a)
    print(f"  data=reorient 100 ms, fit=barrier: chi2={c:8.1f}  -> "
          f"{'yes' if c < reject else 'no'}")
    c, v, e0 = best_with_free_E0({"tau_a": 1e-1}, "tau_r", grid_r)
    print(f"  data=barrier 100 ms, fit=reorient: chi2={c:8.1f}  -> "
          f"{'yes' if c < reject else 'no'}")

    print("\nPRECISION IF IT IS DISTINGUISHABLE (95% profile, E0 profiled):")
    for truth, key, grid, label in [
        ({"tau_a": 1e-1}, "tau_a", grid_a, "barrier 100 ms"),
        ({"tau_r": 1.0}, "tau_r", grid_r, "reorient 1 s"),
    ]:
        c0, _, _ = best_with_free_E0(truth, key, grid)
        inside = [g for g in grid if chi2(truth, {key: g}, free_E0=True,
                                          E0_fit=E0_TARGET) <= c0 + 3.84]
        if inside:
            print(f"  {label:14s}: interval [{min(inside):.2e}, {max(inside):.2e}]"
                  f" -> width {100*(max(inside)-min(inside)):.0f}% of the value")
