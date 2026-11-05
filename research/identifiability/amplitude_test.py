"""Does adding the modulus to the phase recover what phase alone loses?

The negative results in this module are:

  * E0 does not enter arg(E), so phase-only data cannot identify the Gibbs
    elasticity at all;
  * h and D enter the observable only through u = h sqrt(w/D), so phase-only
    data cannot separate the isotherm slope from the diffusion coefficient;
  * the response depends on the relaxation only through arg B, so any two
    mechanisms with the same arg B are indistinguishable in phase.

All three are statements about *phase alone*. Modern instruments report |E|
as well, and |E| carries E0 as a multiplicative scale. The obvious question,
which decides whether this work is a complaint or a methods recommendation,
is whether the joint observable (|E|, delta) is identifiable where the phase
alone is not.

This file answers that by construction: simulate a known truth, fit with
delta only and with (|E|, delta), and report which parameters come back.
"""

from __future__ import annotations

import cmath
import math

import numpy as np
from scipy.optimize import least_squares

from .models import B_single, lag_deg, modulus_from_B

# truth: a concentrated nonionic system with a modest surface relaxation
TRUE = dict(E0=72.0, Cm=0.8, tau=12.0, D=3.8e-10, h=3.275302547717482e-06)
FREQS_PAT = np.logspace(math.log10(0.005), math.log10(1.0), 12)
SIG_PHI = 1.0        # deg
SIG_E_REL = 0.02     # relative on |E|


def observe(freqs, kind, p, sigma_phi=SIG_PHI, sigma_e=SIG_E_REL, truth=TRUE):
    """Synthetic data from the truth, with the parameter vector p ignored.

    Kept explicit so that the 'truth' cannot accidentally be used as a start.
    """
    out = []
    for f in freqs:
        w = 2 * math.pi * f
        B = complex(B_single(w, truth["Cm"], truth["tau"]))
        E = modulus_from_B(w, truth["E0"], truth["D"], truth["h"], B)
        out.append([abs(E), -math.degrees(cmath.phase(E))])
    return np.array(out)


def model_obs(freqs, kind, p):
    """p = (E0, Cm, log10 tau, log10 h[, log10 D]).

    Four parameters when D is held at its true value, five when D is allowed
    to float. The five-parameter case is what shows the h-D degeneracy.
    """
    if len(p) == 4:
        E0, Cm, lt, lh = p
        lD = math.log10(TRUE["D"])
    else:
        E0, Cm, lt, lh, lD = p
    tau, h = 10.0 ** lt, 10.0 ** lh
    out = []
    for f in freqs:
        w = 2 * math.pi * f
        E = modulus_from_B(w, E0, 10.0 ** lD, h, complex(B_single(w, Cm, tau)))
        out.append([abs(E), -math.degrees(cmath.phase(E))])
    return np.array(out)


def residuals(freqs, p, data, kind):
    m = model_obs(freqs, kind, p)
    if kind == "phase":
        return (m[:, 1] - data[:, 1]) / SIG_PHI
    return np.concatenate([(m[:, 0] - data[:, 0]) / (SIG_E_REL * data[:, 0]),
                           (m[:, 1] - data[:, 1]) / SIG_PHI])


START = dict(E0=40.0, Cm=2.0, lt=math.log10(1.0), lh=math.log10(3.3e-6))


def run(freqs=FREQS_PAT, free_D=False, n_starts=60, seed=3):
    """Fit the truth with phase only, and with (|E|, phase).

    free_D=True adds log10 D as a fifth parameter, which is the case that
    exposes the h/sqrt(D) degeneracy.
    """
    data = observe(freqs, "both", None)
    if free_D:
        lo = np.array([1.0, 0.02, -4.0, -8.0, -12.0])
        hi = np.array([400.0, 60.0, 3.0, -3.0, -8.0])
    else:
        lo = np.array([1.0, 0.02, -4.0, -7.0])
        hi = np.array([400.0, 60.0, 3.0, -4.0])
    res = {}
    for kind in ("phase", "both"):
        best = None
        rng = np.random.default_rng(seed)
        s0 = np.array([START["E0"], START["Cm"], START["lt"], START["lh"]])
        if free_D:
            s0 = np.append(s0, math.log10(1e-10))
        starts = [np.clip(s0, lo, hi)]
        for _ in range(n_starts):
            starts.append(lo + rng.random(len(lo)) * (hi - lo))
        for s in starts:
            try:
                r = least_squares(lambda p: residuals(freqs, p, data, kind),
                                  np.clip(s, lo, hi), bounds=(lo, hi),
                                  method="trf", xtol=1e-13, ftol=1e-13,
                                  gtol=1e-13, max_nfev=6000)
            except Exception:
                continue
            c = float(np.sum(r.fun ** 2))
            if best is None or c < best[0]:
                best = (c, r.x)
        chi2, p = best
        d = dict(chi2=chi2, E0=float(p[0]), Cm=float(p[1]),
                 tau=float(10 ** p[2]), h=float(10 ** p[3]),
                 D=float(10 ** p[4]) if free_D else TRUE["D"])
        d["h_over_sqrt_D"] = d["h"] / math.sqrt(d["D"])
        res[kind] = d
    return res


def report(r, truth=TRUE):
    print("Synthetic truth:", {k: f"{v:.3g}" for k, v in truth.items()})
    print(f"Band {FREQS_PAT[0]:.3f}-{FREQS_PAT[-1]:.1f} Hz, "
          f"{len(FREQS_PAT)} frequencies, noise {SIG_PHI} deg and "
          f"{SIG_E_REL:.0%} on |E|\n")
    print(f"{'observable':12s}{'chi2':>9s}{'E0':>9s}{'Cm':>8s}{'tau':>9s}"
          f"{'h':>12s}{'D':>12s}")
    for k, v in r.items():
        print(f"{k:12s}{v['chi2']:9.1f}{v['E0']:9.2f}{v['Cm']:8.3f}{v['tau']:9.2f}"
              f"{v['h']:12.3e}{v['D']:12.3e}")
    print("\nrecovery ratio (fitted / true):")
    for k, v in r.items():
        print(f"  {k:6s} E0 {v['E0']/truth['E0']:7.2f}  Cm {v['Cm']/truth['Cm']:6.2f}  "
              f"tau {v['tau']/truth['tau']:7.2f}  h {v['h']/truth['h']:7.2f}  "
              f"D {v['D']/truth['D']:7.2f}   h/sqrtD "
              f"{v['h_over_sqrt_D']/(truth['h']/math.sqrt(truth['D'])):7.4f}")


if __name__ == "__main__":
    r = run()
    print("Synthetic truth:", {k: (f"{v:.3g}") for k, v in TRUE.items()})
    print(f"Band {FREQS_PAT[0]:.3f}-{FREQS_PAT[-1]:.1f} Hz, "
          f"{len(FREQS_PAT)} frequencies, noise {SIG_PHI} deg and {SIG_E_REL:.0%} on |E|\n")
    print(f"{'observable':12s}{'chi2':>10s}{'E0':>9s}{'Cm':>9s}{'tau':>10s}{'h':>12s}")
    for k, v in r.items():
        print(f"{k:12s}{v['chi2']:10.1f}{v['E0']:9.2f}{v['Cm']:9.3f}"
              f"{v['tau']:10.2f}{v['h']:12.3e}")
    print("\nrecovery (ratio fitted/true):")
    for k, v in r.items():
        print(f"  {k:6s} E0 {v['E0']/TRUE['E0']:6.2f}  Cm {v['Cm']/TRUE['Cm']:6.2f}  "
              f"tau {v['tau']/TRUE['tau']:6.2f}  h {v['h']/TRUE['h']:6.2f}")
