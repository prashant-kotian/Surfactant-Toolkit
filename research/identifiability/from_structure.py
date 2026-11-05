"""Drive the design calculation from molecular structure via SurfactantKit.

Why this matters for the paper
------------------------------
Everything in `run_design.py` runs at a *chosen* concentration ratio c/c_mc,
with the CMC supplied as a constant. That is adequate for a worked example and
inadequate for a claim about surfactants in general, because c/c_mc is not the
designable variable an experimenter sets — they set a concentration in mol/L,
and the CMC differs by orders of magnitude between homologues. A design rule
stated in c/c_mc silently assumes the CMC is known and adjustable, which it is
not: for a given surfactant the CMC is fixed by its structure.

SurfactantKit's Nagarajan molecular-thermodynamic predictor supplies the CMC
from the headgroup and tail length with no fitted parameters (validated against
seven experimental CMCs to 0.15 log10 units, per the module's own tests). That
converts the design question from "at what c/c_mc" into:

    for this surfactant, is any achievable concentration able to resolve a
    surface relaxation with this instrument?

which is the question an experimentalist actually faces when handed a bottle.

Concentration ceiling
---------------------
The upper limit on concentration is not free either: it is the solubility or
cloud point. Above it the system is not a single isotropic phase and the
dilational measurement means something different. We cap the design sweep at
`MAX_C_RATIO` times the CMC, a deliberately generous 500x, and report when the
optimum sits at that ceiling — which is itself the useful answer ("run as
concentrated as the phase diagram permits").

This module is the load-bearing use of the software. The physics modules
elsewhere in SurfactantKit are untouched; this only reads their output.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from .design_engine import (
    designs, discrimination, elasticity_scale, isotherm_slope, spectrum_pressure,
)

# A deliberately generous ceiling on achievable concentration, in multiples of
# the CMC. Real systems are often limited far below this by cloud point or
# solubility; the design sweep reports when the optimum sits on this cap.
MAX_C_RATIO = 500.0

# gamma_inf and K for the Gibbs isotherm.
#   gamma_inf is taken as a saturated-monolayer value typical of polyoxyethylene
#   alkyl ethers, 2.8e-6 mol/m2 (~5.9 A2 per molecule is too small; ~ 50 A2 is
#   the right order), and K is the adsorption constant in m3/mol.
#   These two are the ONLY calibrated numbers in the chain, and they are flagged
#   because the CMC from Nagarajan is independent of both.
GAMMA_INF = 2.8e-6
K_ADS = 140.0
D_MONO = 4.0e-10


def cmc_from_structure(n_carbons: int, headgroup: str = "polyoxyethylene",
                       temperature_K: float = 298.15) -> dict:
    """CMC via SurfactantKit's Nagarajan predictor. mol/L out."""
    from surfactantkit.nagarajan import nagarajan_cmc_predictor
    r = nagarajan_cmc_predictor(n_carbons=n_carbons, headgroup=headgroup,
                                temperature_K=temperature_K)
    return dict(cmc_mol_per_L=float(r.cmc_mol_per_L),
                log10_cmc=float(r.log10_cmc_mol_per_L),
                method=r.method, converged=bool(r.converged))


def design_for_surfactant(n_carbons: int, headgroup="polyoxyethylene",
                          instrument="PAT (oscillating drop/pendant)",
                          C=0.6, tau=12.0, n_freqs=6, layout="log",
                          c_ratios=None) -> dict:
    """Best discrimination achievable for this molecule in this band."""
    cmc = cmc_from_structure(n_carbons, headgroup)
    # UNIT TRAP, hit twice in this project: mol/L -> mol/m3 is a factor of 1000,
    # so multiply by 1e3. A C12 CMC of 1.2e-4 mol/L is 0.12 mol/m3, i.e. 118 uM.
    c_mc = cmc["cmc_mol_per_L"] * 1e3
    assert 1e-6 < c_mc < 100.0, (
        f"CMC {c_mc} mol/m3 outside any physical range -- almost certainly a "
        f"mol/L vs mol/m3 conversion error")
    if not (c_mc > 0 and math.isfinite(c_mc)):
        return dict(n_carbons=n_carbons, error="CMC prediction failed", cmc=cmc)

    freqs = designs(instrument, n_freqs, (layout,))[layout]
    if c_ratios is None:
        c_ratios = np.logspace(0, math.log10(MAX_C_RATIO), 25)

    rows = []
    for r in c_ratios:
        c0 = r * c_mc
        h = isotherm_slope(GAMMA_INF, K_ADS, c0)
        E0 = elasticity_scale(GAMMA_INF, K_ADS, c0)
        prior = [h * f for f in (0.5, 0.7, 1.0, 1.5, 2.0)]
        D = discrimination(freqs, dict(E0=E0, C=C, tau=tau),
                           dict(E0=E0, C=0.0, tau=tau),
                           [p for p in prior if p > 0], "phase")
        rows.append(dict(c_ratio=float(r), c0_mol_m3=float(c0),
                         Pi=float(spectrum_pressure(GAMMA_INF, K_ADS, c0)),
                         D=float(D)))

    best = max(rows, key=lambda x: x["D"])
    at_ceiling = best["c_ratio"] >= 0.98 * MAX_C_RATIO
    return dict(n_carbons=n_carbons, headgroup=headgroup, instrument=instrument,
                cmc=cmc, cmc_mol_m3=c_mc, sweep=rows, best=best,
                at_ceiling=bool(at_ceiling),
                detectable=bool(best["D"] > 3.0))


def run(save=True):
    """Homologue series: does the design verdict change with chain length?"""
    out = []
    for n in (10, 12, 14, 16):
        out.append(design_for_surfactant(n))
    res = dict(constants=dict(gamma_inf=GAMMA_INF, K=K_ADS, D=D_MONO,
                              max_c_ratio=MAX_C_RATIO, C=0.6, tau=12.0),
               series=out)
    if save:
        Path("research/identifiability/structure_design.json").write_text(
            json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    r = run()
    print("Design computed from molecular structure (Nagarajan CMC), PAT band,")
    print(f"relaxation C=0.6 tau=12 s, concentration swept to {MAX_C_RATIO:.0f}x CMC.\n")
    print(f"{'nC':>4s}{'log CMC':>9s}{'CMC (uM)':>10s}{'best c/c_mc':>13s}"
          f"{'c0 (uM)':>10s}{'Pi':>7s}{'D':>7s}  verdict")
    for s in r["series"]:
        if "error" in s:
            print(f"{s['n_carbons']:4d}  {s['error']}")
            continue
        b = s["best"]
        v = "detectable" if s["detectable"] else "NOT detectable"
        if s["at_ceiling"]:
            v += " (at concentration ceiling)"
        print(f"{s['n_carbons']:4d}{s['cmc']['log10_cmc']:9.2f}"
              f"{s['cmc']['cmc_mol_per_L']*1e6:10.1f}{b['c_ratio']:13.1f}"
              f"{b['c0_mol_m3']*1e3:10.1f}{b['Pi']:7.1f}{b['D']:7.2f}  {v}")
