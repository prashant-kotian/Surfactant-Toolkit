"""Translational diffusion and hydrodynamic radius (dynamic light
scattering)."""

from __future__ import annotations

import math

K_B = 1.380649e-23  # J/K, exact SI 2019


def hydrodynamic_radius_stokes_einstein(diffusion_coefficient_cm2_per_s: float, viscosity_mPas: float, temperature_K: float = 298.15) -> float:
    """Hydrodynamic radius (nm) from the Stokes-Einstein equation:
    R_h = k_B*T / (6*pi*eta*D).

    diffusion_coefficient_cm2_per_s: translational diffusion coefficient
    in the commonly-reported DLS unit, cm^2/s (converted internally to
    m^2/s -- do not pass an m^2/s value here, this is a classic unit trap).
    viscosity_mPas: solvent viscosity in mPa.s (=cP; water is ~0.89 at 25C).
    """
    if diffusion_coefficient_cm2_per_s <= 0:
        raise ValueError("diffusion_coefficient_cm2_per_s must be positive")
    if viscosity_mPas <= 0:
        raise ValueError("viscosity_mPas must be positive")
    if temperature_K <= 0:
        raise ValueError("temperature_K must be positive")

    D_SI = diffusion_coefficient_cm2_per_s * 1e-4  # cm^2/s -> m^2/s
    eta_SI = viscosity_mPas * 1e-3  # mPa.s -> Pa.s
    r_h_m = (K_B * temperature_K) / (6.0 * math.pi * eta_SI * D_SI)
    return r_h_m * 1e9  # m -> nm


def perrin_friction_factor(axial_ratio: float) -> float:
    """Perrin translational friction factor F_P for an ellipsoid of
    revolution (spheroid) relative to a SPHERE OF EQUAL VOLUME -- Perrin,
    J. Phys. Radium 5 (1934) 497-511 (closed form as commonly tabulated,
    e.g. Wikipedia's "Perrin friction factors", cross-checked internally
    below since no external worked numeric table was sourced this pass).

    axial_ratio = a/b, where a is the semi-axis along the axis of
    revolution and b is the (equal) equatorial semi-axis. axial_ratio > 1
    is PROLATE (elongated, rod-like -- the relevant case for wormlike/
    rodlike micelles); axial_ratio < 1 is OBLATE (flattened, disk-like).
    At axial_ratio == 1 (a sphere), F_P == 1 exactly (no shape
    correction) -- provable directly from the formula (xi -> 0 gives
    S -> 2, F_P -> 2*1/2 = 1) and confirmed in tests, the one point this
    formula is unambiguously checkable without an external source.

        xi = sqrt(|axial_ratio^2 - 1|) / axial_ratio
        S  = 2*atanh(xi)/xi   (prolate, axial_ratio > 1)
        S  = 2*atan(xi)/xi    (oblate,  axial_ratio < 1)
        F_P = 2 * axial_ratio^(2/3) / S

    Real, sourced alternative-method use: DYNAMIC LIGHT SCATTERING and
    the plain Stokes-Einstein equation (hydrodynamic_radius_stokes_einstein)
    assume a SPHERE -- for a rodlike or wormlike micelle (axial_ratio far
    from 1), the naive Stokes-Einstein calculation returns an inflated
    APPARENT hydrodynamic radius, not the true equal-volume-sphere size.
    See hydrodynamic_radius_perrin_corrected, which uses this factor to
    recover the true equivalent-sphere radius given an assumed/known
    shape (e.g. from the critical packing parameter's predicted
    morphology).
    """
    if axial_ratio <= 0:
        raise ValueError("axial_ratio must be positive")
    if axial_ratio == 1.0:
        return 1.0
    xi = math.sqrt(abs(axial_ratio ** 2 - 1.0)) / axial_ratio
    if axial_ratio > 1.0:
        s = 2.0 * math.atanh(xi) / xi
    else:
        s = 2.0 * math.atan(xi) / xi
    return 2.0 * (axial_ratio ** (2.0 / 3.0)) / s


def hydrodynamic_radius_perrin_corrected(diffusion_coefficient_cm2_per_s: float, viscosity_mPas: float, axial_ratio: float, temperature_K: float = 298.15) -> float:
    """TRUE equivalent-sphere hydrodynamic radius (nm) for a non-
    spherical (rodlike/wormlike or disklike) micelle, correcting the
    naive Stokes-Einstein result for shape via the Perrin friction
    factor: R_eff = R_h,apparent / F_P(axial_ratio), where R_h,apparent
    is what hydrodynamic_radius_stokes_einstein returns (the naive,
    sphere-assuming DLS result). At axial_ratio=1 this returns exactly
    the same value as hydrodynamic_radius_stokes_einstein (F_P=1, no
    correction) -- confirmed in tests.

    axial_ratio must be known or assumed independently (e.g. from the
    critical packing parameter's predicted morphology, or from a
    separate SANS/SAXS shape measurement) -- this function does NOT fit
    or infer axial_ratio from the diffusion data alone (that inverse
    problem is underdetermined from D alone without an independent
    volume estimate; do not guess axial_ratio)."""
    if axial_ratio <= 0:
        raise ValueError("axial_ratio must be positive")
    r_h_apparent = hydrodynamic_radius_stokes_einstein(diffusion_coefficient_cm2_per_s, viscosity_mPas, temperature_K)
    return r_h_apparent / perrin_friction_factor(axial_ratio)
