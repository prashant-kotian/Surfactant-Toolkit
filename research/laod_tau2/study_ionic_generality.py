"""Does the inverse-tau2 method transfer to the ionic surfactants (CTAB, SDS)?

Section 12 established that the isotherm constants originally supplied for
CTAB/SDS (K_F = 18.5 and 1.25 m^3/mol) put the monolayer at theta = 0.94-0.96,
i.e. E0 = 198-270 mN/m, which is unphysical: measured limiting elasticities for
these systems are 20-40 mN/m.

This module fixes that the honest way round. Rather than guessing K and
discovering the coverage is absurd, we INVERT: choose K such that the resulting
Langmuir/Frumkin limiting elasticity matches the literature value for each
surfactant. Every step of that inversion is recorded here so it can be audited.

    E0 = n * RT * Gamma_inf * theta / (1 - theta)
    => theta/(1-theta) = E0 / (n RT Gamma_inf)          [Langmuir]
    => theta = r/(1+r), r = E0/(n RT Gamma_inf)
    => K = theta / (c(1-theta))                         [Langmuir, at conc. c]

For a Frumkin isotherm the same E0 is reached at a slightly different theta
because the lateral-interaction term contributes -2a*theta to the chemical
potential; we solve numerically for theta rather than assume the Langmuir value.

Two things are then asked, per surfactant, across frequency:

  1. DETECTABILITY: does the second-harmonic amplitude clear the 0.025-0.05
     noise floor of a 0.1 mN/m instrument reading?
  2. INVERTIBILITY: how many degrees of phase lag does a factor of 2 in tau2
     buy? (The lag, not the harmonic ratio, is the operative observable --
     see FINDINGS section 10.)

Outputs research/laod_tau2/ionic_generality.json.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from . import harmonics as H
from .lae_data import NOISE_FLOOR_RATIO
from .model import ModelParams, run_oscillation
from .systems import SYSTEMS

R_GAS = 8.314462618

# Literature limiting elasticities, mN/m. Range, and the point value used.
# CTAB: adsorbed layer at the CMC; published Gibbs elasticities for CnTAB
# homologues fall in 20-40 mN/m. SDS: similar, often reported lower (~15-25)
# because of the high CMC and counterion screening.
E0_LITERATURE = {
    "CTAB": {"lo": 20.0, "hi": 40.0, "use": 30.0},
    "SDS": {"lo": 15.0, "hi": 30.0, "use": 22.0},
}

STRAIN = 0.20          # LAOD drive, generous compared with Czakaj's 0.065
T_REF = 298.15


def theta_for_E0(E0_mN_m: float, gamma_inf: float, n_factor: float, a: float) -> float:
    """Coverage giving E0 = n RT Gamma_inf [theta/(1-theta) - 2 a theta^2] (Frumkin).

    The Frumkin generalisation: Pi = -nRT Gamma_inf [ln(1-th) + a th^2] so
    E0 = -dPi/dln Gamma = nRT Gamma_inf [th/(1-th) - 2 a th^2].
    Solve for theta by bisection on (1e-4, 1-1e-4).
    """
    pref = n_factor * R_GAS * T_REF * gamma_inf * 1e3      # mN/m

    def f(th: float) -> float:
        return pref * (th / (1.0 - th) - 2.0 * a * th * th) - E0_mN_m

    lo, hi = 1e-4, 1.0 - 1e-6
    if f(lo) > 0:
        raise ValueError("E0 below the available range")
    if f(hi) < 0:
        raise ValueError("E0 above the available range (divergent regime)")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def K_from_theta(theta: float, c_mol_m3: float, a: float) -> float:
    """Invert the Frumkin isotherm K c = th/(1-th) * exp(-2 a th)."""
    return theta / (1.0 - theta) * math.exp(-2.0 * a * theta) / c_mol_m3


def make_params(name: str, mult: float, tau2: float) -> tuple[ModelParams, dict]:
    """Corrected-isotherm ModelParams for `name` at `mult` x CMC."""
    s = SYSTEMS[name]
    c = mult * s.cmc_mM                                  # mol/m^3 == mM
    e0 = E0_LITERATURE[name]["use"]
    th = theta_for_E0(e0, s.gamma_inf, float(s.n_ions), s.isotherm_a)
    K = K_from_theta(th, c, s.isotherm_a)
    p = ModelParams(
        c0=c, K=K, gamma_inf=s.gamma_inf, D=s.D, n_factor=float(s.n_ions),
        capacity=s.capacity(c), tau2=tau2, frumkin_a=s.isotherm_a,
        temperature_K=T_REF,
    )
    info = {"theta": th, "K_m3_per_mol": K, "E0_target": e0,
            "capacity": p.capacity, "c_mM": c}
    return p, info


def probe(p: ModelParams, omega: float, eps: float = STRAIN) -> dict:
    """I1, I2, phase lag at one frequency."""
    r = run_oscillation(p, omega, eps, n_cycles=16)
    t = np.asarray(r.t)
    hg = H.harmonics(t, r.gamma_mN_m, omega)
    return {
        "I1": hg[1], "I2": hg[2],
        "lag_deg": H.lag_deg(t, r.gamma_mN_m, r.area_rel, omega),
    }


def leverage(name: str, mult: float, tau2: float, omega: float) -> dict:
    """Phase-lag sensitivity to tau2: degrees per factor-of-2 change."""
    base_p, info = make_params(name, mult, tau2)
    lo = probe(ModelParams(**{**base_p.__dict__, "tau2": tau2 / 2.0}), omega)
    hi = probe(ModelParams(**{**base_p.__dict__, "tau2": tau2 * 2.0}), omega)
    mid = probe(base_p, omega)
    return {
        "omega": omega, "frames_per_cycle_at_25fps": 25.0 / (omega / (2 * math.pi)),
        "I1": mid["I1"], "I2": mid["I2"], "lag_deg": mid["lag_deg"],
        "lag_tau2_half": lo["lag_deg"], "lag_tau2_double": hi["lag_deg"],
        "dlag_per_4x": abs(hi["lag_deg"] - lo["lag_deg"]),
        "I2_over_noise": mid["I2"] / NOISE_FLOOR_RATIO[1],
        "detectable": mid["I2"] > NOISE_FLOOR_RATIO[1],
        "theta": info["theta"], "K": info["K_m3_per_mol"], "capacity": info["capacity"],
    }


def main() -> dict:
    # tau2 point values from the benchmark matrix (Antigravity tau2_validation.md)
    TAU2 = {"CTAB": 0.065, "SDS": 0.038}
    MULT = 1.5
    out = {"strain": STRAIN, "mult_cmc": MULT,
           "noise_floor_ratio": NOISE_FLOOR_RATIO[1], "systems": {}}

    for name, tau2 in TAU2.items():
        rows = []
        # PAT-1M oscillatory band is 0.005-1 Hz => omega <= 6.28 rad/s.
        # Go beyond it to locate the optimum, so we can state where the
        # instrument boundary actually falls.
        for omega in (0.3, 0.6, 1.26, 2.0, 3.0, 4.5, 6.28, 10.0, 15.7, 25.0, 40.0, 62.8):
            rows.append(leverage(name, MULT, tau2, omega))
        best_i2 = max(rows, key=lambda r: r["I2"])
        best_lag = max(rows, key=lambda r: r["dlag_per_4x"])
        out["systems"][name] = {
            "tau2_s": tau2, "theta": rows[0]["theta"], "K": rows[0]["K"],
            "capacity": rows[0]["capacity"],
            "rows": rows,
            "peak_I2": {"omega": best_i2["omega"], "I2": best_i2["I2"],
                        "detectable": best_i2["detectable"]},
            "peak_leverage": {"omega": best_lag["omega"],
                              "dlag_per_4x": best_lag["dlag_per_4x"]},
            "usable_omega": [r["omega"] for r in rows
                             if r["detectable"] and r["dlag_per_4x"] >= 3.0],
        }
    return out


if __name__ == "__main__":
    res = main()
    path = Path(__file__).with_name("ionic_generality.json")
    path.write_text(json.dumps(res, indent=2))
    for name, d in res["systems"].items():
        print(f"\n=== {name}  tau2 = {d['tau2_s']*1e3:.0f} ms, "
              f"theta = {d['theta']:.3f}, K = {d['K']:.3f}, Cm = {d['capacity']:.0f} ===")
        print(f"{'omega':>7} {'f(Hz)':>7} {'fpc':>6} {'I2':>9} {'I2/noise':>9} "
              f"{'lag':>7} {'dlag/4x':>8}  verdict")
        for r in d["rows"]:
            f = r["omega"] / (2 * math.pi)
            ok = r["detectable"] and r["dlag_per_4x"] >= 3.0
            print(f"{r['omega']:7.2f} {f:7.3f} {r['frames_per_cycle_at_25fps']:6.1f} "
                  f"{r['I2']:9.4f} {r['I2_over_noise']:9.2f} {r['lag_deg']:7.2f} "
                  f"{r['dlag_per_4x']:8.2f}  {'USABLE' if ok else ''}")
        print(f"peak I2 = {d['peak_I2']['I2']:.4f} at omega = {d['peak_I2']['omega']}")
        print(f"usable omega = {d['usable_omega']}")
    print(f"\nwrote {path}")
