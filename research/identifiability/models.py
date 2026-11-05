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
    and hence, at a given omega,

        |delta| -> 45 deg + |arg B(omega)| / 2.

    B = 1 gives 45 deg -- the classical LVDT ceiling -- so that ceiling is a
    property of *instantaneous* surface re-equilibration, not of diffusion.

    CORRECTION (round 29). An earlier version of this module asserted that a
    passive relaxation is confined to 45-67.5 deg, because Im B <= 0 implies
    |arg B| < 45 in the limit of a single dominant time constant. That bound is
    wrong. The maximum of |arg B| over frequency, for B = 1 + C/(1+ix), occurs at
    x = sqrt(C+1) and equals atan(C / (2 sqrt(C+1))), which increases without
    bound towards 45 deg as C grows. The 67.5 figure came from scanning tau at a
    single frequency, so the optimum was never reached -- a search artefact
    masquerading as a theorem.

    The correct statement is the inversion, implemented in `capacity_from_cap`:
    a measured ceiling bounds C, and a ceiling near 90 deg requires a large C.
    Fainerman et al. (Langmuir 2008, 24, 6447) make exactly this point in words,
    citing Joos, for the fast micelle-dissociation case. So the ceiling is
    prior art as an observation; what is not is the closed-form bound below.
    """
    return 45.0 + abs(math.degrees(cmath.phase(B))) / 2.0


def max_cap_single_relaxation(C):
    """Largest |delta| achievable by B = 1 + C/(1+iw tau), over all w and tau.

    |arg B| is maximised at w tau = sqrt(C+1), giving
        tan(|arg B|_max) = C / (2 sqrt(C+1)).
    """
    if C <= 0:
        return 45.0
    return 45.0 + math.degrees(math.atan(C / (2.0 * math.sqrt(C + 1.0)))) / 2.0


def capacity_from_cap(delta_max_deg):
    """Inverse of `max_cap_single_relaxation`: the C implied by a measured ceiling.

    Solves T = tan(2(delta - 45)) for C, with T = C / (2 sqrt(C+1)), giving
        C = 2T^2 + 2T sqrt(T^2 + 1).
    A measurement of the ceiling therefore lower-bounds the surface relaxation
    strength with no fitting at all.
    """
    if delta_max_deg <= 45.0:
        return 0.0
    T = math.tan(math.radians(2.0 * (delta_max_deg - 45.0)))
    return 2.0 * T * T + 2.0 * T * math.sqrt(T * T + 1.0)


