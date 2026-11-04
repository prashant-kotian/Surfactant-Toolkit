"""Analytic linear transfer function for the immobile-reservoir model.

Added 2026-11-05 in response to external audit item 2.5. The auditor derived
this and showed it reproduces the PDE solver's small-strain phase to ~0.01 deg.
It was verified independently here before being adopted (FINDINGS 30.8).

Derivation
----------
Linearise about the static state with time dependence exp(+i omega t).

    bulk monomer      dc/dt = D d2c/dx2 - dM/dt
    reservoir         dM/dt = -(M - C_m c/n) / tau
    isotherm slope    h     = dGamma_eq/dc
    surface           i omega G = D dc/dx|0            (adsorbed amount small
                                                       compared with the flux)
    elasticity        E     = dPi/dlnGamma * (dG/dG)   -> report |E|, arg E

Eliminating M in favour of c gives a single diffusion equation for the *total*
perturbation with an effective retardation:

    B(omega) = 1 + C_m / (n * (1 + i omega tau))

so with q^2 = i omega B / D and the decaying root Re(q) > 0,

    E(omega) = E0 * i omega / (i omega + D q / h).

Limiting cases (all checked in tests):
  * C_m -> 0            -> B = 1, recovers Lucassen-van den Tempel.
  * omega tau -> 0      -> B = 1 + C_m/n (equilibrium reservoir; maximum lag).
  * omega tau -> inf    -> B = 1 (reservoir frozen out).
  * |arg E| -> 45 deg as D q / h dominates, for ANY B, because arg(i omega /
    (i omega + real)) -> 45 deg. Hence delta > 45 deg cannot arise from the
    diffusion term alone once the interface is in local equilibrium -- it needs
    E0 or h to be complex, or the geometry to be non-planar.

Convention
----------
`phase_deg` returns delta = phi_gamma - phi_area in degrees with positive =
lag, matching harmonics.lag_deg. No wrapping into [0, 360).
"""

from __future__ import annotations

import cmath
import math
from dataclasses import dataclass


@dataclass
class LinearResponse:
    omega: float
    modulus: float        # |E|, mN/m
    phase_deg: float      # arg(E), the standard "phase angle", positive = gamma LEADS
    B: complex            # reservoir retardation factor
    q: complex            # bulk wavenumber, 1/m

    @property
    def lag_deg(self) -> float:
        """delta = phi_gamma - phi_area, positive = gamma LAGS = -arg(E).

        The relationship arg(E) = -delta follows from E = |E|e^{i arg E} and
        gamma = |gamma|e^{i(omega t - delta)} with A = |A|e^{i omega t}:
        E = dgamma/dlnA carries the phase -delta directly.

        The whole project quoted "delta > 45 deg" while its code computed
        arg(E). Those have opposite signs. Every delta in the paper therefore
        has the wrong sign under its own stated convention (FINDINGS 30.9).
        """
        return -self.phase_deg


def isotherm_slope(gamma_inf: float, K: float, c0: float,
                   frumkin_a: float = 0.0) -> float:
    """h = dGamma_eq/dc at c0, m (mol/m^2 per mol/m^3). Langmuir or Frumkin."""
    if frumkin_a == 0.0:
        return gamma_inf * K / (1.0 + K * c0) ** 2
    th = K * c0 / (1.0 + K * c0)
    dth = th * (1.0 - th) / (c0 * (1.0 - 2.0 * frumkin_a * th * (1.0 - th)))
    return gamma_inf * dth


def elasticity_scale(gamma_inf: float, K: float, c0: float, n_factor: float = 1.0,
                     temperature_K: float = 298.15, frumkin_a: float = 0.0) -> float:
    """E0 = nRT Gamma_inf * dtheta/dlnc, mN/m. The Gibbs elasticity."""
    R = 8.314462618
    h = isotherm_slope(gamma_inf, K, c0, frumkin_a)
    return n_factor * R * temperature_K * 1e3 * h * c0


def linear_response(omega: float, *, c0: float, K: float, gamma_inf: float, D: float,
                    capacity: float = 0.0, tau2: float = 1.0, n_agg: float = 1.0,
                    n_factor: float = 1.0, temperature_K: float = 298.15,
                    frumkin_a: float = 0.0) -> LinearResponse:
    """Analytic small-strain response. See module docstring for the derivation."""
    if omega <= 0:
        raise ValueError("omega must be positive")
    if tau2 <= 0 and capacity != 0.0:
        raise ValueError("tau2 must be positive when a reservoir is present")
    if n_agg <= 0:
        raise ValueError("n_agg must be positive")
    if c0 <= 0 or K <= 0 or D <= 0:
        raise ValueError("c0, K and D must be positive")

    h = isotherm_slope(gamma_inf, K, c0, frumkin_a)
    if h <= 0:
        raise ValueError("isotherm slope must be positive (choose the stable branch)")

    B = 1.0 + capacity / (n_agg * (1.0 + 1j * omega * tau2)) if capacity else 1.0 + 0j
    q = cmath.sqrt(1j * omega * B / D)
    if q.real < 0:                       # decaying branch
        q = -q
    E0 = elasticity_scale(gamma_inf, K, c0, n_factor, temperature_K, frumkin_a)
    E = E0 * 1j * omega / (1j * omega + D * q / h)
    return LinearResponse(omega=omega, modulus=abs(E), phase_deg=math.degrees(cmath.phase(E)),
                          B=B, q=q)


def phase_deg(f: float, **kw) -> float:
    """Convenience: phase lag delta in degrees at frequency f Hz."""
    return linear_response(2.0 * math.pi * f, **kw).phase_deg


def lvdt_phase_deg(f: float, *, c0: float, K: float, gamma_inf: float, D: float,
                   **kw) -> float:
    """Lucassen-van den Tempel limit (C_m = 0). Never exceeds 45 deg in magnitude."""
    return phase_deg(f, c0=c0, K=K, gamma_inf=gamma_inf, D=D, capacity=0.0, **_lvdt_kw(kw))


def _lvdt_kw(kw):
    keep = {k: v for k, v in kw.items()
            if k in ("n_factor", "temperature_K", "frumkin_a")}
    return keep
