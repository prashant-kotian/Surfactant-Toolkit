"""Analytic limits used as anchors for the numerical model.

These are closed forms, derived from the published theory, and they are the
only reason to trust anything the PDE model says. If `model.py` cannot
reproduce these in the limits where they hold, it is wrong.

1. Lucassen & van den Tempel (1972), diffusion-controlled adsorption of a
   single surfactant with a Langmuir isotherm, small oscillation amplitude:

       E(omega) = E0 / (1 + zeta + i*zeta),          zeta = (dc/dGamma)*sqrt(D/(2*omega))

       E' = E0*(1 + zeta) / (1 + 2*zeta + 2*zeta^2)
       E''= E0*zeta       / (1 + 2*zeta + 2*zeta^2)

   with E0 = -dgamma/dlnGamma = n*R*T*Gamma*/(1 - Gamma/Gamma_inf) evaluated at
   the equilibrium coverage, and (dc/dGamma)_eq = (1 + K*c)^2/(Gamma_inf*K) for
   the Langmuir isotherm.

2. Fast-reservoir limit (omega*tau2 -> 0). When the micelles re-equilibrate
   with the monomers faster than any diffusion or area-change timescale, the
   two species diffuse as one with an effective diffusivity

       D_eff = D*(1 + C_m)

   where C_m is the micellar storage capacity (d(micellar surfactant)/
   d(monomer) at fixed total). LVDT then applies unchanged with D -> D_eff.
   This is the standard buffered-diffusion result.

3. Slow-reservoir limit (omega*tau2 -> infinity). The micelles cannot respond
   within a cycle, so they are inert and the response is LVDT with D and no
   capacity: D_eff = D.

Limits 2 and 3 bracket every micellar response, so a model that does not sit
between them for intermediate omega*tau2 is not physical.

Derivation note (E0 from the Langmuir-Szyszkowski equation of state):
    gamma = gamma0 + n*R*T*Gamma_inf*ln(1 - Gamma/Gamma_inf)
    dgamma/dGamma = -n*R*T/(1 - Gamma/Gamma_inf)
    E0 = -Gamma*dgamma/dGamma = n*R*T*Gamma_inf*theta/(1 - theta)
Consistent with adsorption.gibbs_marangoni_elasticity_langmuir, which is the
library function this is cross-checked against in the tests.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

R_GAS = 8.314462618


@dataclass(frozen=True)
class LangmuirState:
    """Equilibrium state of a Langmuir monolayer at bulk concentration c0."""

    c0: float          # mol/m^3
    K: float           # m^3/mol   (Langmuir adsorption constant)
    gamma_inf: float   # mol/m^2   (maximum surface excess)
    D: float           # m^2/s
    gamma0: float      # mN/m      (solvent surface tension)
    n_factor: float = 1.0
    temperature_K: float = 298.15

    @property
    def theta(self) -> float:
        """Fractional surface coverage at equilibrium."""
        return self.K * self.c0 / (1.0 + self.K * self.c0)

    @property
    def gamma_eq(self) -> float:
        """Equilibrium surface excess, mol/m^2."""
        return self.gamma_inf * self.theta

    @property
    def sigma_eq(self) -> float:
        """Equilibrium surface tension, mN/m."""
        return self.gamma0 + self.n_factor * R_GAS * self.temperature_K * self.gamma_inf * 1e3 * math.log(
            1.0 - self.theta
        )

    @property
    def E0(self) -> float:
        """Gibbs elasticity -dgamma/dln(Gamma), mN/m."""
        return (
            self.n_factor
            * R_GAS
            * self.temperature_K
            * self.gamma_inf
            * 1e3
            * self.theta
            / (1.0 - self.theta)
        )

    @property
    def dc_dgamma(self) -> float:
        """(dc/dGamma) at equilibrium for the Langmuir isotherm, m^-1."""
        return (1.0 + self.K * self.c0) ** 2 / (self.gamma_inf * self.K)

    def zeta(self, omega: float, D_eff: float | None = None) -> float:
        D = self.D if D_eff is None else D_eff
        return self.dc_dgamma * math.sqrt(D / (2.0 * omega))

    def lvdt(self, omega: float, D_eff: float | None = None) -> tuple[float, float]:
        """(E', E'') in mN/m from the Lucassen-van den Tempel closed form."""
        z = self.zeta(omega, D_eff)
        den = 1.0 + 2.0 * z + 2.0 * z * z
        return self.E0 * (1.0 + z) / den, self.E0 * z / den

    def lvdt_lag_deg(self, omega: float, D_eff: float | None = None) -> float:
        """Phase lag of gamma behind the imposed area oscillation, degrees."""
        ep, epp = self.lvdt(omega, D_eff)
        return math.degrees(math.atan2(epp, ep))

    def lvdt_harmonics(self, omega: float, eps: float, D_eff: float | None = None):
        """Predicted fundamental amplitude and second/third harmonic ratios.

        In the linear (small-eps) limit the response is a pure sinusoid: the
        first harmonic has amplitude eps*|E| and all higher harmonics vanish.
        Returned as (I1, I2_over_I1, I3_over_I1) with I2, I3 exactly zero --
        this is the falsifiable statement that higher harmonics are a
        *nonlinear* effect and cannot be produced by linear viscoelasticity.
        """
        ep, epp = self.lvdt(omega, D_eff)
        return eps * math.hypot(ep, epp), 0.0, 0.0

    def lvdt_crossover_frequency(self, D_eff: float | None = None) -> float:
        """omega at which zeta = 1/2, where E'' peaks. Diagnostic only."""
        z_target = 0.5
        return (self.dc_dgamma**2) * (self.D if D_eff is None else D_eff) / (2.0 * z_target**2)


def fast_reservoir_D_eff(D: float, capacity: float) -> float:
    """D_eff = D*(1 + C_m): the omega*tau2 -> 0 limit."""
    if capacity < 0:
        raise ValueError("capacity must be >= 0")
    return D * (1.0 + capacity)


def slow_reservoir_D_eff(D: float, capacity: float) -> float:
    """D_eff = D: the omega*tau2 -> infinity limit (micelles frozen out)."""
    if capacity < 0:
        raise ValueError("capacity must be >= 0")
    return D
