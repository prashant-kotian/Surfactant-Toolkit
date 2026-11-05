"""Does the design metric predict what published experiments could resolve?

This is the only validation available without new experiments, and it is worth
doing properly because it is the difference between a plausible design score
and a predictive one.

The prediction to test
----------------------
The design metric D (worst-case model distance in noise units) rises steeply
with surface pressure. If that is real, then a published study should be able
to resolve a relaxation at high surface pressure and not at low, *using the
same instrument and the same surfactant*.

Fainerman et al. (Langmuir 2008, 24, 6447) provide exactly that comparison:
C14EO8, one instrument (buoyant bubble PAT), 0.005-0.2 Hz, concentrations from
5 umol/L to 500 umol/L. Their reported outcomes:

  * c0 = 5 uM to cmc (7 uM): theory (eq 2, plain LVDT) describes the data --
    no additional mechanism needed. Surface pressure 2.3-3.0 mN/m.
  * c0 > cmc at 12-30 uM: LVDT alone fits, micelle effects small. Pi = 4-7 mN/m.
  * c0 = 50-500 uM (7-70x cmc): eq 3 is needed and agrees "excellently"; the
    phase angle exceeds the 45 deg ceiling and is attributed to fast micelle
    dissociation. Pi = 10-33 mN/m.

So the literature itself separates cleanly by surface pressure, not by
concentration per se. The design metric predicts the same separation. That is a
consistency check on real data, not a fit.

A second, independent check: the earlier work in this repository fitted four
mechanism models to the digitised phase data at 30 uM (Pi = 7.4 mN/m) and found
LVDT and the relaxation model equal on AICc (Delta AICc = 0.19, weights 0.52 /
0.47). The design metric at that surface pressure predicts non-discrimination.
Two independent analyses of the same dataset, one a full model comparison and
one a design calculation, agree.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from .design_engine import (
    NOISE, coverage, designs, discrimination, elasticity_scale, isotherm_slope,
    spectrum_pressure,
)

# Same constants as run_design.py; see the note there on units.
GAMMA_INF = 2.8e-6
K_ADS = 140.0
D_MONO = 4.0e-10

# Fainerman's reported outcomes, in his own terms.
OUTCOMES = [
    (5.0, "LVDT sufficient (below cmc)"),
    (7.0, "LVDT sufficient (at cmc)"),
    (12.0, "LVDT sufficient (small micelle effect)"),
    (15.0, "LVDT sufficient"),
    (20.0, "LVDT sufficient"),
    (30.0, "transitional; digitised data cannot separate models"),
    (50.0, "eq 3 required, excellent agreement"),
    (70.0, "eq 3 required"),
    (100.0, "eq 3 required"),
    (200.0, "eq 3 required"),
    (500.0, "eq 3 required"),
]


def design_score(c0_uM, instrument="PAT (oscillating drop/pendant)",
                 C=0.6, tau=12.0):
    c0 = c0_uM * 1e-3                 # umol/L -> mol/m^3
    h = isotherm_slope(GAMMA_INF, K_ADS, c0)
    E0 = elasticity_scale(GAMMA_INF, K_ADS, c0)
    freqs = designs(instrument, 6, ("log",))["log"]
    prior = [h * f for f in (0.5, 0.7, 1.0, 1.5, 2.0)]
    truth = dict(E0=E0, C=C, tau=tau, D=D_MONO)
    rival = dict(E0=E0, C=0.0, tau=tau, D=D_MONO)
    d = discrimination(freqs, truth, rival, prior, "phase")
    return dict(c0_uM=c0_uM, theta=coverage(K_ADS, c0),
                Pi=spectrum_pressure(GAMMA_INF, K_ADS, c0), E0=E0, D=d)


def run(save=True):
    rows = []
    for c0, outcome in OUTCOMES:
        r = design_score(c0)
        r["reported"] = outcome
        r["resolved"] = outcome.startswith("eq 3")
        rows.append(r)
    if save:
        Path("research/identifiability/validation_results.json").write_text(
            json.dumps(rows, indent=1))
    return rows


if __name__ == "__main__":
    rows = run()
    print("Fainerman 2008 C14EO8, one instrument, concentrations 5-500 uM.")
    print("'resolved' means the authors needed a mechanism beyond plain LVDT.\n")
    print(f"{'c0 uM':>7s}{'theta':>7s}{'Pi mN/m':>9s}{'E0':>7s}{'D':>7s}  reported")
    for r in rows:
        flag = "  <-- resolves" if r["resolved"] else ""
        print(f"{r['c0_uM']:7.0f}{r['theta']:7.3f}{r['Pi']:9.2f}{r['E0']:7.1f}"
              f"{r['D']:7.2f}  {r['reported']}{flag}")
    res = [r["D"] for r in rows if r["resolved"]]
    un = [r["D"] for r in rows if not r["resolved"] and r["c0_uM"] < 30]
    print(f"\nmax D where the authors did NOT need a mechanism: {max(un):.2f}")
    print(f"min D where they DID:                              {min(res):.2f}")
    print(f"separation achieved: {'yes' if min(res) > max(un) else 'NO -- overlapping'}")
