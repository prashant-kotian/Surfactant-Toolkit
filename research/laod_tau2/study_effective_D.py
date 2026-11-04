"""Is an explicit micellar reservoir NECESSARY, or does Fainerman's effective D* suffice?

THE QUESTION
------------
Fainerman, Petkov & Miller, Langmuir 2008, 24, 6447, eqs (3)-(4) replace the monomer
diffusion coefficient D by an effective value

    D* = D (1 + beta)(1 + R beta),   beta = (c0 - cmc)/cmc,   R = D_micelle/D ~ 0.25

and keep the Lucassen-van den Tempel functional form with c0 set to cmc. On that basis
they reproduce the viscoelasticity modulus |E| for C14EO8 at every concentration up to
70 x cmc, using ONE extra parameter (beta, fixed by the concentration).

Our model carries the micellar reservoir as an explicit species with a relaxation time
tau2. If D* is adequate, our extra state variable buys nothing and the paper's central
"coupled reservoir" contribution is unnecessary. If D* fails somewhere, we need to say
exactly where and by how much.

WHAT IS TESTED
--------------
Three forward models at matched bulk chemistry:

  A  coupled   - our explicit reservoir (capacity C_m, relaxation tau2)
  B  LvT       - plain diffusion-controlled, micelles absent (D* = D)
  C  D*        - Fainerman/Joos effective diffusion coefficient

B and C are analytic. For A we solve the existing model.

TWO COMPARISONS, because they answer different questions:

  1. ACROSS CONCENTRATION at fixed frequency (what Fainerman actually did). Sweep
     c0/cmc and see whether C tracks A and B does not. This tests whether the
     effective-D shortcut reproduces the reservoir's effect.

  2. THE LIMIT TEST. As Omega = omega*tau2 -> 0 the reservoir tracks quasi-statically
     and the system should become indistinguishable from a pure diffusion problem with
     a shifted D. Fainerman's eq 4 is precisely that limit. So we predict C -> A as
     Omega -> 0 and C diverging from A as Omega grows. If that holds, eq 4 is the
     Omega -> 0 asymptote of our model rather than a competitor to it - which is a
     defensible and useful statement either way.

HOW D* ENTERS THE ANALYTIC FORMS
--------------------------------
With xi = sqrt(omega_D / (2 omega)) and omega_D = D (dGamma/dc)^-2:

    |E| = E0 (1 + 2 xi + 2 xi^2)^(-1/2),    phi = arctan[ xi / (1 + xi) ]

For model C, D -> D* so xi -> xi / sqrt((1+beta)(1+R beta)). For model B, beta = 0.

Mapping between our C_m and Fainerman's beta: the reservoir can supply
C_m * c0 of micellar surfactant, and beta = (c0-cmc)/cmc is the excess over cmc in the
same units, so the identification is beta = C_m * (c0/cmc) / ((c0/cmc) - 1) ... but
that presumes a stoichiometry we have not established. We therefore treat beta as the
free parameter it is in their paper and report the D* prediction across a beta range,
flagging the mapping as an open modelling question rather than inventing one.

Outputs research/laod_tau2/effective_D_test.json
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from . import harmonics as H
from .model import ModelParams, run_oscillation
from .systems import SYSTEMS

R_GAS = 8.314462618
S = SYSTEMS["C12E8"]
T_REF = 298.15
STRAIN = 0.20
C_TOT = 1.5 * S.cmc_mM
R_MIC = 0.25                    # Fainerman's R = D_m/D
CAPACITY = 47.5
TAU2 = 4.0


def K_for_E0(E0: float) -> float:
    pref = S.n_ions * R_GAS * T_REF * S.gamma_inf * 1e3
    r = E0 / pref
    th = r / (1.0 + r)
    return th / (1.0 - th) / C_TOT


K_ANCHOR = K_for_E0(45.3)


# --------------------------------------------------------------------------
# analytic diffusion-controlled predictions
# --------------------------------------------------------------------------
def lvt(D: float, omega: float, E0: float, dGdc: float) -> tuple[float, float]:
    """Lucassen-van den Tempel |E| and phase, for diffusion coefficient D."""
    # Fainerman eq 1: omega_D = D (dGamma/dc)^-2, i.e. D DIVIDED by (dG/dc)^2.
    omega_D = D / dGdc ** 2
    xi = math.sqrt(omega_D / (2.0 * omega))
    mod = E0 / math.sqrt(1.0 + 2.0 * xi + 2.0 * xi ** 2)
    phi = math.degrees(math.atan(xi / (1.0 + xi)))
    return mod, phi


def d_star(D: float, beta: float) -> float:
    """Fainerman/Joos effective diffusion coefficient, eq (4)."""
    return D * (1.0 + beta) * (1.0 + R_MIC * beta)


def coupled(omega: float) -> tuple[float, float]:
    """Our explicit-reservoir model: |E| and lag from the solver."""
    p = ModelParams(c0=C_TOT, K=K_ANCHOR, gamma_inf=S.gamma_inf, D=S.D,
                    n_factor=float(S.n_ions), capacity=CAPACITY, tau2=TAU2,
                    temperature_K=T_REF)
    r = run_oscillation(p, omega, STRAIN, n_cycles=20, n_return_cycles=10,
                        n_per_cycle=200)
    t = np.asarray(r.t)
    g = np.asarray(r.gamma_mN_m)
    phi_lag = H.lag_deg(t, g, r.area_rel, omega)
    a = np.asarray(r.area_rel)
    amp = 0.5 * (g.max() - g.min())
    eps = 0.5 * (a.max() - a.min())
    mod = amp / eps if eps > 0 else float("nan")
    return mod, phi_lag


def equilibrium_derivative(c_mM: float) -> tuple[float, float]:
    """E0 and dGamma/dc at the working concentration for the anchored isotherm."""
    p = ModelParams(c0=C_TOT, K=K_ANCHOR, gamma_inf=S.gamma_inf, D=S.D,
                    n_factor=float(S.n_ions), capacity=CAPACITY, tau2=TAU2,
                    temperature_K=T_REF)
    th = p._theta_eq(c_mM * 1.0)          # c in the model's own units (mol/m^3 == mM)
    dcsdG = 1.0 / (K_ANCHOR * S.gamma_inf * (1.0 - th) ** 2)
    E0 = R_GAS * T_REF * S.gamma_inf * 1e3 * th / (1.0 - th) * S.n_ions
    return E0, 1.0 / dcsdG


def main() -> dict:
    E0, dGdc = equilibrium_derivative(C_TOT)
    out = {"setup": {"c0_over_cmc": 1.5, "capacity": CAPACITY, "tau2": TAU2,
                     "R_micelle": R_MIC, "E0_mN_m": E0, "dGamma_dc": dGdc,
                     "strain": STRAIN},
           "across_frequency": [], "across_beta": []}

    # ---- comparison 1: frequency sweep, coupled vs LvT vs D* at several beta
    print(f"E0={E0:.2f} mN/m  dGamma/dc={dGdc:.4e} m  omega_D={S.D/dGdc**2:.4f} Hz")
    print(f"{'f(Hz)':>7} {'cpl lag':>9} {'LvT':>9} {'D* b=0.5':>9} {'D* b=5':>9} "
          f"{'cpl |E|':>8}")
    for f in (0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 1.6, 2.5):
        om = 2 * math.pi * f
        mod_c, lag_c = coupled(om)
        _, phi_b = lvt(S.D, om, E0, dGdc)
        _, phi_b05 = lvt(d_star(S.D, 0.5), om, E0, dGdc)
        _, phi_b5 = lvt(d_star(S.D, 5.0), om, E0, dGdc)
        out["across_frequency"].append(
            {"f_Hz": f, "coupled_lag": lag_c, "coupled_mod": mod_c,
             "LvT_lag": phi_b, "Dstar_b0.5_lag": phi_b05, "Dstar_b5_lag": phi_b5})
        print(f"{f:7.3f} {lag_c:9.2f} {phi_b:9.2f} {phi_b05:9.2f} {phi_b5:9.2f}"
              f"   {mod_c:7.2f}")

    # ---- comparison 2: how far can D* reach? sweep beta at fixed frequencies
    print("\nCan D* ever exceed 45 deg? (it cannot, by construction - phi -> 45 as w -> 0)")
    for beta in (0.0, 0.5, 2.0, 10.0, 100.0):
        best = max(lvt(d_star(S.D, beta), 2 * math.pi * 0.005, E0, dGdc)[1],
                   lvt(d_star(S.D, beta), 2 * math.pi * 0.0001, E0, dGdc)[1])
        out["across_beta"].append({"beta": beta, "max_phi_over_f": best,
                                   "coupled_lag_at_1Hz": next(
                                       r["coupled_lag"] for r in out["across_frequency"]
                                       if r["f_Hz"] == 1.0)})
        print(f"  beta={beta:7.1f}  max phi at f->0: {best:6.2f} deg")

    return out


if __name__ == "__main__":
    res = main()
    path = Path(__file__).with_name("effective_D_test.json")
    path.write_text(json.dumps(res, indent=2, default=float))
    print(f"\nwrote {path}")
