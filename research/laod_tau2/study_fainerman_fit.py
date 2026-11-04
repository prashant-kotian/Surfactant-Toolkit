"""Fit our coupled-reservoir model to Fainerman et al.'s MEASURED phase angles.

WHY
---
Every quantitative claim in this project has been validated against its own forward
model. Fainerman, Petkov & Miller (Langmuir 2008, 24, 6447, Figure 5) publish measured
phase angle vs frequency for C14EO8 above the cmc. That is the only external micellar
dilatational dataset available.

REWRITTEN 2026-11-05 (FINDINGS 30.10) after the external audit found three defects in
the previous version:

  1. **Wrong CMC.** The script hardcoded CMC = 7 umol/L for C12E8, which has no
     traced source, and normalised the C14EO8 data by it. Rusdi et al. (Bull. Chem.
     Soc. Jpn. 2005, 78, 604) Table 1 gives **C14E8 CMC = 1.81e-5 mol/dm3 =
     18.1 umol/L at 298.2 K**. Fainerman's 30 umol/L is therefore **1.66 x cmc**,
     not 4.29 x. The whole concentration axis was off by a factor of 2.6.

  2. **Homolog substitution.** The fit used C12E8 parameters (Gamma_inf, D) for
     C14EO8 data. That is an exploratory stress test, not a measurement, and it is
     now labelled as such. C14E8 values from Rusdi are used where available.

  3. **The archive did not match the notebook.** fainerman_fit.json recorded
     C_m = 1.0 with E0 fixed at 45.3, while FINDINGS 27 reported C_m = 0.5 with
     E0 = 35 from an uncommitted ad-hoc script. Both are now produced by this
     script, which scans E0 *and* capacity, and both are written to the archive.

Also changed: the model is now the analytic linear transfer function
(`linear_response.py`) rather than the PDE, because the audit showed they agree to
0.01 deg and the analytic form makes the E0 scan cheap. Fainerman's oscillation was
5 % area, i.e. eps = 0.05, not the 0.20 this project design-pointed at; the data are
also small-amplitude, so the linear response is the appropriate comparison.

STATISTICS -- corrected
-----------------------
The previous output labelled sqrt(SSE/dof) "RMSE". That is a residual standard
deviation after two fitted parameters, not the RMS residual over the observations.
Both are reported now. With sigma = 1 deg the chi-square is ~193 on 4 dof, which is a
severe lack of fit, not "barely above noise".

Outputs research/laod_tau2/fainerman_fit.json
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from .linear_response import linear_response

# ---------------------------------------------------------------------------
# Measured data: Fainerman, Petkov & Miller, Langmuir 2008, 24, 6447, Fig. 5
# filled squares = 30 umol/L. Digitised by eye from a 700 dpi render.
# C14E8 cmc = 18.1 umol/L at 298.2 K (Rusdi et al. 2005, Table 1).
# => 30/18.1 = 1.66 x cmc. Digitisation uncertainty ~1 deg.
# ---------------------------------------------------------------------------
FREQS = (0.005, 0.010, 0.020, 0.050, 0.100, 0.200)
PHI_OBS = (46.5, 50.5, 52.0, 48.5, 41.0, 32.5)
PHI_ERR = 1.0

CMC_C14E8_UM = 18.1             # Rusdi 2005 Table 1, 298.2 K
C_DATA_UM = 30.0                # Fainerman 2008 Fig. 5
C_OVER_CMC = C_DATA_UM / CMC_C14E8_UM

# C14E8 parameters. Gamma_inf and D remain C12E8 values: no traced C14E8 pair was
# found. This is a stated weakness, not a silent one.
GAMMA_INF = 2.80e-6
D_MONO = 3.8e-10
T_REF = 298.15
STRAIN = 0.05                   # Fainerman used 5 % area oscillations

C_TOT = C_DATA_UM * 1e-3        # umol/L -> mol/m^3, TOTAL loading
# The transport variable is MONOMER. Above the CMC the monomer is pinned near the
# cmc, so the model sees ~18.1 umol/L, not 30. Passing the total here was the
# concentration error fixed in FINDINGS 30.1.
C_MONOMER = min(C_DATA_UM, CMC_C14E8_UM) * 1e-3


def K_for_E0(E0: float, c0: float) -> float:
    """Langmuir K that yields elasticity E0 at monomer concentration c0."""
    R = 8.314462618
    pref = R * T_REF * GAMMA_INF * 1e3
    r = E0 / pref
    th = r / (1.0 + r)
    return th / (1.0 - th) / c0


def phi_model(cap: float, tau2: float, f: float, E0: float = 45.3,
              n_agg: float = 1.0) -> float:
    """delta in degrees from the analytic transfer function."""
    K = K_for_E0(E0, C_MONOMER)
    r = linear_response(2.0 * math.pi * f, c0=C_MONOMER, K=K, gamma_inf=GAMMA_INF,
                        D=D_MONO, capacity=cap, tau2=tau2, n_agg=n_agg,
                        temperature_K=T_REF)
    # The model unambiguously gives gamma LEADING the area oscillation (delta < 0);
    # verified against the PDE solver and by direct peak timing. Fainerman's sign
    # convention cannot be recovered from a digitised magnitude-only figure, so the
    # fit is on the MAGNITUDE of the phase difference. This is a stated ambiguity,
    # not a sign chosen to improve the fit: no sign choice alters |delta|.
    return abs(r.lag_deg)


def sse(cap: float, tau2: float, E0: float, n_agg: float = 1.0) -> float:
    return float(sum((phi_model(cap, tau2, f, E0, n_agg) - obs) ** 2
                     for f, obs in zip(FREQS, PHI_OBS)))


def closure_capacity(cmc_ratio: float, n_bar: float = 95.0) -> float:
    """C_m = n_bar*(m-1) at the data concentration."""
    return n_bar * (cmc_ratio - 1.0)


def main() -> dict:
    caps = (0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 8.0, 10.0, 20.0, 30.0, 62.5, 100.0)
    taus = (0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0, 16.0, 24.0, 40.0)
    e0s = (30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 150.0)

    closure = closure_capacity(C_OVER_CMC)
    print(f"C14E8 at {C_DATA_UM:g} umol/L = {C_OVER_CMC:.2f} x cmc "
          f"(cmc {CMC_C14E8_UM:g} umol/L, Rusdi 2005)")
    print(f"closure would give C_m = {closure:.1f} at that concentration\n")

    best = None
    grid = {}
    # joint grid over (E0, C_m, tau2)
    for E0 in e0s:
        for cap in caps:
            for t2 in taus:
                v = sse(cap, t2, E0)
                grid[f"{E0}_{cap}_{t2}"] = v
                if best is None or v < best[3]:
                    best = (E0, cap, t2, v)

    E0b, capb, taub, sseb = best

    # PROFILE over capacity: for each C_m, the best achievable SSE after
    # optimising the nuisances (E0, tau2). This is the quantity to interpret;
    # the joint argmin is not, because E0 and C_m are strongly degenerate.
    profile = {}
    for cap in caps:
        prof_best = None
        for E0 in e0s:
            for t2 in taus:
                v = grid[f"{E0}_{cap}_{t2}"]
                if prof_best is None or v < prof_best[0]:
                    prof_best = (v, E0, t2)
        profile[cap] = {"sse": prof_best[0], "E0": prof_best[1],
                        "tau2": prof_best[2]}
    n = len(FREQS)
    dof = n - 2
    rms_resid = math.sqrt(sseb / n)          # RMS residual over observations
    resid_sd = math.sqrt(sseb / dof)         # residual SD (2 params fitted)
    chi2_red = sseb / dof / PHI_ERR ** 2

    print(f"best: E0 = {E0b:g}, C_m = {capb:g}, tau2 = {taub:g} s")
    print(f"SSE = {sseb:.2f}  RMS residual = {rms_resid:.2f} deg  "
          f"residual SD = {resid_sd:.2f} deg  reduced chi2 = {chi2_red:.0f}\n")
    print(f"{'f(Hz)':>8} {'obs':>7} {'model':>8} {'resid':>8}")
    rows = []
    for f, obs in zip(FREQS, PHI_OBS):
        m = phi_model(capb, taub, f, E0b)
        _ = m
        rows.append({"f": f, "obs": obs, "model": m, "resid": m - obs})
        print(f"{f:8.3f} {obs:7.1f} {m:8.1f} {m - obs:+8.1f}")

    # boundary diagnostics: a bound optimum is not an estimate
    at_cap_bound = capb in (min(caps), max(caps))
    at_tau_bound = taub in (min(taus), max(taus))
    at_e0_bound = E0b in (min(e0s), max(e0s))
    print(f"\nboundary hits: capacity={at_cap_bound} tau2={at_tau_bound} "
          f"E0={at_e0_bound}")

    # profile over capacity, with a chi-square confidence interval on 1 parameter
    sse_min = min(v["sse"] for v in profile.values())
    print(f"\nprofile over C_m (nuisances E0, tau2 re-optimised at each):")
    print(f"{'C_m':>7} {'bestSSE':>9} {'dSSE':>8} {'E0':>6} {'tau2':>6}")
    inside68, inside95 = [], []
    for cap in caps:
        pr = profile[cap]
        d = pr["sse"] - sse_min
        flag = ""
        if d <= 1.0:
            inside68.append(cap); flag = "  <1sigma"
        if d <= 3.84:
            inside95.append(cap)
        print(f"{cap:7.1f} {pr['sse']:9.2f} {d:8.2f} {pr['E0']:6g} {pr['tau2']:6g}{flag}")
    print(f"\nC_m inside 1-sigma profile: {inside68}")
    print(f"C_m inside 95% profile:     {inside95}")
    print(f"closure requires C_m = {closure:.1f}")
    print(f"closure {closure:.1f} vs fit {capb:g} -> ratio {closure / capb:.0f}x")

    return {
        "fitted": {
            "E0_mN_m": E0b, "capacity": capb, "tau2_s": taub,
            "sse": sseb, "n_obs": n, "dof": dof,
            "rms_residual_deg": rms_resid,
            "residual_sd_deg": resid_sd,
            "reduced_chi2_sigma1deg": chi2_red,
            "at_boundary": {"capacity": at_cap_bound, "tau2": at_tau_bound,
                            "E0": at_e0_bound},
        },
        "closure": {"capacity_at_data_concentration": closure,
                    "ratio_closure_over_fit": closure / capb},
        "data": {
            "system": "C14EO8",
            "concentration_umol_per_L": C_DATA_UM,
            "cmc_umol_per_L": CMC_C14E8_UM,
            "c_over_cmc": C_OVER_CMC,
            "cmc_source": "Rusdi, Moroi, Hlaing & Matsuoka, Bull. Chem. Soc. Jpn. "
                          "2005, 78, 604, Table 1 (298.2 K)",
            "freqs": list(FREQS), "phi_obs": list(PHI_OBS), "phi_err": PHI_ERR,
            "source": "Fainerman 2008 Fig 5, 30 umol/L, visual digitisation",
            "digitisation": "by eye, 700 dpi render, ~1 deg",
            "gamma_inf_and_D": "C14E8 values not traced; C12E8 substituted "
                               "explicitly (stated weakness)",
            "strain": STRAIN,
        },
        "residuals": rows,
        "profile_over_capacity": {str(k): v for k, v in profile.items()},
        "profile_inside_1sigma": inside68,
        "profile_inside_95": inside95,
        "grid": {"E0s": list(e0s), "caps": list(caps), "taus": list(taus),
                 "sse": grid},
        "notes": [
            "Corrected CMC 2026-11-05: previous run used 7 umol/L, giving 4.29x cmc.",
            "RMS residual and residual SD are distinct; both reported.",
            "Boundary optimum does not identify any single cause.",
        ],
    }


if __name__ == "__main__":
    res = main()
    path = Path(__file__).with_name("fainerman_fit.json")
    path.write_text(json.dumps(res, indent=2, default=float))
    print(f"\nwrote {path}")
