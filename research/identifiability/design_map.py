"""Does the information in a frequency window suffice to identify a mechanism?

This asks the design question, not the fitting question: for a given
instrument band, concentration and surface relaxation time, how much phase
excursion does a relaxation produce above the instantaneous-equilibrium
(LVDT) baseline, and is that excursion larger than the noise?

The quantity that matters is

    signal(tau, band) = max_f |delta(f; tau)| - min_f |delta(f; tau)|

because a *constant* offset in phase is absorbed by the unknown isotherm slope
h: h and D enter only through u = h sqrt(w/D), and a shift in h slides the
whole curve along the frequency axis. Only the shape identifies the mechanism.
This is the same reason a single-frequency measurement can never identify a
mechanism, no matter how precise it is.

The output is a map: for each (tau, band) pair, how many sigma the mechanism
raises. Where that number is below 1, the literature's mechanistic assignment
is not supported by its own data, and no amount of fitting changes that.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from .models import B_single, lag_deg

D = 3.8e-10
H = 3.275302547717482e-06
E0 = 60.0

# Instrument bands, from the manufacturers' stated ranges and typical practice
BANDS = {
    "microtensiometer": np.logspace(np.log10(1.7e-3), np.log10(1.7), 24),
    "PAT (oscillating drop)": np.logspace(np.log10(0.005), np.log10(1.0), 24),
    "CPT (capillary pressure)": np.logspace(np.log10(0.1), np.log10(500.0), 24),
    "extended (PAT + CPT)": np.logspace(np.log10(0.005), np.log10(500.0), 40),
}


def excursion(tau, freqs, Cm, sigma=1.0, D=D, h=H, n_scale=120):
    """Peak phase excursion that CANNOT be absorbed by rescaling h, in sigmas.

    h and D enter only through u = h sqrt(w/D), so changing h slides the whole
    curve along the frequency axis. A relaxation that merely shifts the curve
    is therefore not evidence for a mechanism: it is evidence for a different
    isotherm slope. The first version of this function ignored that and
    reported large excursions even for tau far outside the band, where B is
    effectively real and constant and the only effect is a shift. That is the
    same nuisance-conflation error that invalidated the earlier project's
    fitting, so it is corrected here explicitly: the best-fitting rescaling of
    h is subtracted before the excursion is measured.
    """
    w = 2 * np.pi * freqs
    y = np.array([lag_deg(x, E0, D, h, complex(B_single(x, Cm, tau))) for x in w])
    best = np.inf
    for sc in np.logspace(-1.5, 1.5, n_scale):
        base = np.array([lag_deg(x, E0, D, h * sc, 1.0 + 0j) for x in w])
        best = min(best, float(np.max(np.abs(y - base))))
    return best / sigma


def maps(sigma=1.0, Cm_grid=(0.05, 0.2, 1.0, 5.0, 20.0),
         tau_grid=np.logspace(-4, 3, 29), save=True):
    out = {}
    for band, freqs in BANDS.items():
        rows = {}
        for Cm in Cm_grid:
            rows[str(Cm)] = [excursion(t, freqs, Cm, sigma) for t in tau_grid]
        out[band] = dict(freqs=freqs.tolist(), Cm=rows)
    res = dict(sigma=sigma, tau_grid=tau_grid.tolist(), bands=out,
               definition=("min over h-rescalings of max_f |delta(tau)-delta(Cm=0)| / sigma; "
                           "the nuisance isotherm slope is profiled out so that a mere "
                           "frequency shift does not count as a mechanism signature"))
    if save:
        Path("research/identifiability/discriminability_map.json").write_text(
            json.dumps(res, indent=1))
    return res


def summarise(res, threshold=3.0):
    """Which (band, Cm) combinations reach `threshold` sigma anywhere in tau,
    and over what range of tau."""
    tau = np.array(res["tau_grid"])
    lines = []
    for band, d in res["bands"].items():
        for Cm, vals in d["Cm"].items():
            v = np.array(vals)
            ok = tau[v >= threshold]
            if len(ok):
                lines.append((band, float(Cm), float(ok.min()), float(ok.max()),
                              float(v.max())))
            else:
                lines.append((band, float(Cm), None, None, float(v.max())))
    return lines


if __name__ == "__main__":
    res = maps()
    print("Phase excursion above the LVDT baseline, in units of the phase noise.")
    print("3 sigma is the conventional detection threshold.\n")
    print(f"{'band':26s}{'Cm':>7s}{'detectable tau range (s)':>28s}{'peak':>8s}")
    for band, Cm, lo, hi, peak in summarise(res):
        rng = f"{lo:.1e} - {hi:.1e}" if lo is not None else "none"
        print(f"{band:26s}{Cm:7.2f}{rng:>28s}{peak:8.1f}")
    print("\nFor reference: the C14E8 fit gives Cm = 0.77, tau = 11.5 s.")
