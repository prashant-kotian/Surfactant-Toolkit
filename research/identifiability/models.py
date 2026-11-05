"""Mechanism models for interfacial dilatational rheology.

Convention: E(omega) = E' + iE'' is the complex dilatational modulus, and the
reported phase is the LAG delta = phi_gamma - phi_area = -arg(E), positive when
the surface tension trails the area. This matches research/laod_tau2 and the
project's stated convention (FINDINGS 30.9).

Every model here is written through a single object: the SURFACE RELAXATION
ADMITTANCE B(omega), defined by

    B(omega) = 1 + sum_k  C_k / (1 + i omega tau_k)

for a discrete set of relaxations, or the corresponding integral for a
distribution. Physically B is the factor by which the layer's response to a
compression is retarded relative to instantaneous local equilibrium: a
perturbation of the bulk concentration drives an excess that takes time to
re-equilibrate, so the effective diffusion coefficient seen by the wave is
B times the free one. B = 1 is the Lucassen-van den Tempel limit.

Writing all mechanisms this way is not a modelling convenience, it is what
makes them comparable. The central result of this module (see `lag_cap_deg`)
is that the observable depends on B only through its argument, which is why
mechanisms that differ wildly in their underlying physics can be
indistinguishable in the phase.

Nothing here is fitted to anything; the fits live in fit_data.py.
"""

from __future__ import annotations

import cmath
import math

import numpy as np

R_GAS = 8.314462618


# ---------------------------------------------------------------- admittances

def B_single(omega, capacity, tau):
    """One Debye-like relaxation: a single exchange time."""
    return 1.0 + capacity / (1.0 + 1j * omega * tau)


def B_two_mode(omega, capacity, tau_1, tau_2, weight=0.5):
    """Two relaxations in parallel. Reduces to B_single when the taus meet."""
    return 1.0 + capacity * (weight / (1.0 + 1j * omega * tau_1)
                             + (1.0 - weight) / (1.0 + 1j * omega * tau_2))


def B_log_normal(omega, capacity, tau_c, width, n=256, span=6.0):
    """Log-normal distribution of relaxation times, median tau_c, log width.

    This is the continuum model of the 'distribution of adsorption
    coefficients' and of surface structural heterogeneity. It exists here to
    represent the claim that a *spread* of times, rather than any specific
    time, explains the data.
    """
    z = np.linspace(-span / 2, span / 2, n)
    tau = tau_c * np.exp(width * z)
    w = np.exp(-0.5 * (z / max(width, 1e-9)) ** 2)
    w = w / np.trapezoid(w, z)
    integrand = 1.0 / (1.0 + 1j * omega * tau) * w
    return 1.0 + capacity * np.trapezoid(integrand, z)


def B_maxwell_like(omega, capacity, tau, exponent):
    """Fractional (Cole-Cole) relaxation: (1 + (i w tau)^a)^-1.

    Included because fractional-power-law relaxation is what a broad,
    self-similar distribution of times looks like, and because it is the
    model that surface-rheology papers reach for when a single time fails.
    """
    return 1.0 + capacity / (1.0 + (1j * omega * tau) ** exponent)


# ---------------------------------------------------------------- modulus

def modulus_from_B(omega, E0, D, h, B):
    """E = E0 iw / (iw + D q / h) with q^2 = i w B / D.

    With u = h sqrt(w/D) > 0 this is algebraically identical to

        E = E0 / (1 + sqrt(B) exp(-i pi/4) / u).

    Verify by dividing numerator and denominator by iw:
        D q / (h i w) = sqrt(i w B D) / (h i w)
                      = sqrt(B) sqrt(i w D) / (h i w)
                      = sqrt(B) sqrt(D) e^{-i pi/4} / (h sqrt(w))
                      = sqrt(B) e^{-i pi/4} / u.
    The test suite checks the two forms against each other and against the
    PDE-calibrated `laod_tau2.linear_response`.
    """
    u = h * math.sqrt(omega / D)
    return E0 / (1.0 + cmath.sqrt(B) * cmath.exp(-1j * math.pi / 4) / u)


def lag_deg(omega, E0, D, h, B):
    """|delta| in degrees, where delta = -arg E, i.e. the magnitude of the
    phase angle between surface tension and area.

    Squaring the sign convention here would be a mistake that has already cost
    this project once (FINDINGS 30.9), so it is pinned by a test against the
    calibrated reference rather than by inspection.
    """
    return abs(math.degrees(cmath.phase(modulus_from_B(omega, E0, D, h, B))))


def lag_cap_deg(omega, B):
    """Low-frequency limit of |delta| for a given admittance, in degrees.

    As u -> 0, arg(1 + Z) -> arg Z = arg B / 2 - 45 deg, so arg E -> 45 - arg B/2
    and hence

        |delta|_max = 45 deg + |arg B| / 2.

    Consequences:

    * B = 1 gives 45 deg exactly -- the classical LVDT ceiling. That ceiling is
      therefore not a property of diffusion as such but of *instantaneous*
      surface re-equilibration.
    * A passive relaxation has Im B <= 0 (a positive-weighted sum of
      1/(1+iw tau)), so 0 <= |arg B| < 45 deg and every passive mechanism is
      confined to |delta| in [0, 67.5] deg.
    * A measured |delta| above 67.5 deg falsifies the entire class at once --
      barrier, reorientation, exchange, or any distribution of times -- and
      points outside it: surface viscoelasticity, a non-passive layer, or
      breakdown of the planar / local-equilibrium assumption.

    The factor of one half is what makes the bound tight: the exponent of B
    enters the observable halved, so the phase is a compressed image of the
    surface relaxation spectrum.
    """
    return 45.0 + abs(math.degrees(cmath.phase(B))) / 2.0
