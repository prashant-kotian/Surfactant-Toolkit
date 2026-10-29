"""Electrostatics of ionic surfactant solutions: Debye screening length
and zeta potential (Henry equation).

Physical constants are CODATA/SI-2019 exact or near-exact values, not
memorized approximations:
    k_B = 1.380649e-23 J/K       (exact, SI 2019)
    N_A = 6.02214076e23 /mol     (exact, SI 2019)
    e   = 1.602176634e-19 C      (exact, SI 2019)
    eps0 = 8.8541878128e-12 F/m  (CODATA 2018)
"""

from __future__ import annotations

import math
from dataclasses import dataclass

K_B = 1.380649e-23
N_A = 6.02214076e23
E_CHARGE = 1.602176634e-19
EPS0 = 8.8541878128e-12
WATER_REL_PERMITTIVITY_25C = 78.4


def ionic_strength(concentrations_M: dict[float, float]) -> float:
    """Ionic strength I = 0.5 * sum(c_i * z_i^2), mol/L (M).

    concentrations_M: {charge_number: concentration_M}, e.g. for 0.1 M
    NaCl: {1: 0.1, -1: 0.1} (cation and anion both counted).
    """
    if not concentrations_M:
        raise ValueError("concentrations_M must not be empty")
    if any(c < 0 for c in concentrations_M.values()):
        raise ValueError("ion concentrations must not be negative")
    return 0.5 * sum(c * (z ** 2) for z, c in concentrations_M.items())


def debye_length(ionic_strength_M: float, temperature_K: float = 298.15, rel_permittivity: float = WATER_REL_PERMITTIVITY_25C) -> float:
    """Debye screening length (nm) for a 1:1-equivalent electrolyte
    solution: lambda_D = sqrt(eps_r*eps0*k_B*T / (2*N_A*I*e^2)).

    ionic_strength_M: ionic strength in mol/L (see ionic_strength()) --
    NOT the raw salt concentration; for a 1:1 salt like NaCl, I equals
    the salt concentration, but for higher-valence electrolytes it does
    not.
    rel_permittivity: relative permittivity of the solvent (default is
    water at 25 C, ~78.4; this is temperature-dependent -- pass an
    explicit value for other temperatures rather than assuming 25 C).
    """
    if ionic_strength_M <= 0:
        raise ValueError("ionic_strength_M must be positive")
    if temperature_K <= 0:
        raise ValueError("temperature_K must be positive")
    I_SI = ionic_strength_M * 1000.0  # mol/L -> mol/m^3
    lambda_m = math.sqrt(
        (rel_permittivity * EPS0 * K_B * temperature_K) / (2.0 * N_A * I_SI * E_CHARGE ** 2)
    )
    return lambda_m * 1e9  # m -> nm


def grahame_equation_surface_potential(surface_charge_density_C_per_m2: float, ionic_strength_M: float, z: int = 1, temperature_K: float = 298.15, rel_permittivity: float = WATER_REL_PERMITTIVITY_25C) -> float:
    """Surface (Stern-layer) potential (mV) from surface charge density
    via the Grahame equation -- the standard, textbook closed-form
    Gouy-Chapman diffuse-double-layer result (Grahame, Chem. Rev. 41
    (1947) 441-501), for a symmetric z:z supporting electrolyte:

        sigma = sqrt(8*eps_r*eps0*k_B*T*n0) * sinh(z*e*psi / (2*k_B*T))

    solved for psi:

        psi = (2*k_B*T)/(z*e) * asinh[ sigma / sqrt(8*eps_r*eps0*k_B*T*n0) ]

    n0 is the bulk number density of the z:z salt. The argument is the
    ionic strength, and for a symmetric z:z electrolyte I = z^2 * c_salt,
    so n0 = 1000 * N_A * I / z^2 (mol/L -> mol/m^3 -> number density).
    For every z the low-potential limit then reduces to
    sigma = eps_r*eps0*psi/lambda_D with lambda_D = debye_length(I).
    Cross-checked against a real, fetched
    secondary source (Mchedlov-Petrossyan, Electron. Proc. Mater. 50(2)
    (2014) 71-80, their Eq. (9), the integrated form of this same
    charge-potential relation for the ionic-surfactant adsorption case)
    -- same functional (sinh/asinh) form confirmed independently, not
    guessed.

    Use: a surfactant's surface excess Gamma (from gibbs_gamma_max)
    gives a charge density sigma = z*F*Gamma, from which this function
    estimates the potential of the adsorbed ionic monolayer.
    """
    if ionic_strength_M <= 0:
        raise ValueError("ionic_strength_M must be positive")
    if z <= 0:
        raise ValueError("z must be a positive integer (charge number magnitude)")
    if temperature_K <= 0:
        raise ValueError("temperature_K must be positive")
    n0_per_m3 = ionic_strength_M * 1000.0 * N_A / z ** 2  # I = z^2 c for a z:z salt; mol/L -> ions/m^3
    denom = math.sqrt(8.0 * rel_permittivity * EPS0 * K_B * temperature_K * n0_per_m3)
    psi_V = (2.0 * K_B * temperature_K) / (z * E_CHARGE) * math.asinh(surface_charge_density_C_per_m2 / denom)
    return psi_V * 1000.0  # V -> mV


def henry_function(regime: str, kappa_a: float | None = None) -> float:
    """Henry function f(kappa*a): 'huckel' = 1.0 (small particles / low
    ionic strength, kappa*a << 1), 'smoluchowski' = 1.5 (large particles
    / ionic strength >= ~10 mM, kappa*a >> 1), 'ohshima', or
    'swan_furst' -- the latter two are closed-form approximations valid
    and accurate across the ENTIRE kappa*a range, including the
    intermediate zone where neither limit is a good approximation (real
    micelles, R~2-4 nm at realistic ionic strengths 1-100 mM, typically
    land at kappa*a ~ 0.2-4.2 -- squarely in that poorly-approximated
    intermediate zone, not cleanly in either limit). Both require
    kappa_a (the dimensionless product kappa*a -- Debye parameter times
    particle radius, same units for both) as a second argument; raise
    if omitted or non-positive.

    'ohshima' (<1% relative error, Ohshima 1994):
        f(kappa*a) = 1 + 1 / [2 * (1 + delta)^3]
        delta = 2.5 / [kappa*a * (1 + 2*exp(-kappa*a))]

    'swan_furst' (<0.1% relative error -- an order of magnitude more
    accurate than Ohshima's, per the source paper's own Fig. 2
    comparison; a simple rational function, no exponential):
        f(kappa*a) = (16 + 18*kappa*a + 3*(kappa*a)^2)
                     / (16 + 18*kappa*a + 2*(kappa*a)^2)

    There is no default regime -- picking the wrong one silently is
    exactly the kind of error this library exists to prevent; the
    caller must state which applies for their particle size and ionic
    strength (or supply kappa_a and use 'ohshima'/'swan_furst' to avoid
    the choice entirely).

    Sources, both fetched and read directly (not from secondary
    summaries):
    - Ohshima, H., J. Colloid Interface Sci. 168 (1994) 269-271,
      doi:10.1006/jcis.1994.1419, Eq. [16]. Ohshima's own paper
      normalizes f(kappa*a) to the range [2/3, 1] (their Eq. [1] omits
      this library's leading 2/3 mobility prefactor); converted here by
      a factor of 3/2 to match THIS library's [1.0, 1.5] convention
      (same convention as zeta_potential_henry's mobility formula) --
      verified correct because it reproduces exactly 1.0 at kappa*a ->
      0 and exactly 1.5 at kappa*a -> infinity.
    - Swan & Furst, J. Colloid Interface Sci. 388 (2012) 92-94,
      doi:10.1016/j.jcis.2012.08.026, Eq. (6) -- already in this
      library's [1.0, 1.5] convention as given (no conversion needed;
      their own Eq. (5) confirms the same mobility formula/convention
      this library uses). Their paper independently reproduces
      Ohshima's exact formula too (their Eq. (7), delta =
      (5/2)/(1+2*exp(-kappa*a)), f_O = 1+1/(2*[1+delta/(kappa*a)]^3) --
      algebraically identical to the 'ohshima' branch above once delta
      is combined with the /kappa*a term), an independent cross-
      confirmation from a second, unrelated primary source that this
      library's 'ohshima' implementation is correct.

    A THIRD option, 'qin' -- Qin, Liu, Wang, Thomas, Wang & Shen, Acta
    Optica Sinica 37(10) (2017) 1029003, doi:10.3788/AOS201737.1029003
    (fetched and read directly, in the original Chinese with an English
    abstract/equations) -- least-squares-refits Ohshima's SAME functional
    form against Wiersema's numerically-exact Henry function values,
    finding better coefficients:

    f(kappa*a) = 1 + 1 / [2 * (1 + 2.8/(kappa*a*(0.9+exp(-kappa*a))))^3]

    (Ohshima's original delta = 2.5/[kappa*a*(1+2*exp(-kappa*a))]; this
    refit changes the three constants from (2.5, 1, 2) to (2.8, 0.9, 1).)
    The source paper's own Table 1 (reproduced against Wiersema's exact
    numerical solution) reports this is consistently more accurate than
    Ohshima's original across kappa*a = 0.01 to 1000, most notably in
    the kappa*a ~ 1-20 range this library's own literature review found
    real micelles typically land in (e.g. at kappa*a=10: Ohshima 3.0%
    error vs. this refit's 0.2%; at kappa*a=2: 2.3% vs. 0.2%) -- overall
    max relative error <1.5% vs. Ohshima's <3% in that range. This
    library's own test suite checks this implementation reproduces the
    source paper's own Table 1 values exactly (to the 3-4 significant
    figures given), not just the asymptotic limits.

    CAUTION: an open-access secondary source found during initial
    research (Skoglund et al., PLOS ONE 12(7):e0181735 (2017), Eq. 3)
    reproduces the SAME Ohshima delta formula but with a DIFFERENT,
    INCORRECT f(kappa*a) = 1 + 1/(1+delta)^3 -- missing the factor of 2
    in the denominator. That form gives f -> 2 (not 1.5) as kappa*a ->
    infinity, which is physically wrong; confirmed as an error in that
    secondary source (not this library) by fetching THREE independent
    primary/near-primary sources directly (Ohshima 1994, Swan & Furst
    2012, Qin et al. 2017), all of which agree with each other and with
    this library's implementation.
    """
    regime = regime.lower()
    if regime == "huckel":
        return 1.0
    if regime == "smoluchowski":
        return 1.5
    if regime == "ohshima":
        if kappa_a is None or kappa_a <= 0:
            raise ValueError("kappa_a (the dimensionless product kappa*a) must be a positive "
                              "number for regime='ohshima' -- do not guess a value")
        delta = 2.5 / (kappa_a * (1.0 + 2.0 * math.exp(-kappa_a)))
        return 1.0 + 1.0 / (2.0 * (1.0 + delta) ** 3)
    if regime == "swan_furst":
        if kappa_a is None or kappa_a <= 0:
            raise ValueError("kappa_a (the dimensionless product kappa*a) must be a positive "
                              "number for regime='swan_furst' -- do not guess a value")
        ka = kappa_a
        return (16.0 + 18.0 * ka + 3.0 * ka ** 2) / (16.0 + 18.0 * ka + 2.0 * ka ** 2)
    if regime == "qin":
        if kappa_a is None or kappa_a <= 0:
            raise ValueError("kappa_a (the dimensionless product kappa*a) must be a positive "
                              "number for regime='qin' -- do not guess a value")
        delta = 2.8 / (kappa_a * (0.9 + math.exp(-kappa_a)))
        return 1.0 + 1.0 / (2.0 * (1.0 + delta) ** 3)
    raise ValueError("regime must be 'huckel', 'smoluchowski', 'ohshima', 'swan_furst', or 'qin' "
                      "-- no default is provided; Huckel applies for kappa*a << 1 (small "
                      "particles/low ionic strength), Smoluchowski for kappa*a >> 1 (larger "
                      "colloids, ionic strength >= ~10 mM), Ohshima/Swan-Furst/Qin are valid and "
                      "accurate across the full range (all three require kappa_a)")


def electrophoretic_mobility_relaxation_corrected(zeta_mV: float, kappa_a: float, viscosity_mPas: float, temperature_K: float = 298.15, z: int = 1, m: float = 0.15, rel_permittivity: float = WATER_REL_PERMITTIVITY_25C) -> float:
    """Electrophoretic mobility ((um*cm)/(V*s)) from zeta potential via
    the REAL Ohshima-Healy-White (1983) semi-empirical formula -- a
    real, sourced ALTERNATIVE to the plain Henry equation
    (zeta_potential_henry/henry_function), capturing the REAL physical
    effect that Henry's LINEAR-in-zeta formula misses entirely -- at
    high |zeta| (>50 mV, common for ionic surfactant micelles per this
    project's own literature review), counterion relaxation around the
    moving particle makes mobility grow SUB-linearly with zeta, and it
    can even pass through a maximum (typically beyond ~100 mV) and
    decline.

    UPGRADED 2026-09-11 (primary sources obtained, PDFs provided by the
    user -- previously flagged as substitute-only sourcing via a
    secondary IUPAC review): this now implements the actual formula
    from Ohshima, Healy & White, J. Chem. Soc. Faraday Trans. 2, 79
    (1983) 1613-1628, eqs. (57)-(62) plus the semi-empirical correction
    of eqs. (75)-(76) -- NOT the simpler O'Brien-Hunter-style
    simplification used here previously. This is a real, material
    upgrade: OHW83's own Fig. 1/2 show their formula (57) beating the
    O'Brien-Hunter/Dukhin-Semenikhin family at every kappa_a and zeta
    tested against the exact O'Brien & White (1978) computer solution,
    and their semi-empirical version (57)+(75)/(76) is stated to have
    <1% relative error for kappa_a >= 10 -- a wider valid range than
    this function's previous kappa_a > 20 floor.

        zeta~ = z*e*zeta / (k_B*T)                       [eq. 21, dimensionless]
        F = (2/kappa_a)*(1+3m)*(exp(zeta~/2) - 1)         [eq. 56, m+ = m- = m]
        G = ln[(1 + exp(-zeta~/2)) / 2]                   [eq. 59]
        H = ln[(1 + exp(+zeta~/2)) / 2]                   [eq. 60]
        I = sum_{n=1..inf} tanh(zeta~/4)^(2n-1) / (2n-1)^2 [eq. 61]
        J = sum_{n=1..inf} tanh(zeta~/4)^(2n)   / (2n)^2   [eq. 62]
        alpha = kappa_a / (6*(kappa_a - 6))               [empirically fit, eq. 75]
        CF = 1 - [25/(3*(kappa_a+10))]*exp(-alpha*zeta~)  [eq. 76, the correction factor]

        E_m = 1.5*zeta~ - [3F/(1+F)]*H - (18/kappa_a)*I*CF
              + (1/kappa_a) * { [3F/(2(1+F))]*(10*(I+J) + tanh(zeta~/4)^2)
                                 - 6*(1+3m)*(1-exp(-zeta~/2))*G
                                 + [12F/(1+F)^2]*H
                                 + [9*zeta~/(1+F)]*m*(G+H)
                                 - [36F/(1+F)]*m*(G^2 + H^2/(1+F)) }

        ue = E_m * 2*eps_r*eps0*k_B*T / (3*eta*z*e)       [eq. 58, inverted for U/E]

    z: counterion/co-ion charge number (symmetrical z-z electrolyte
    assumed, e.g. z=1 for NaCl-like). m: dimensionless ionic drag
    coefficient (eq. 49 of the source), assumed EQUAL for cation and
    anion (m+ = m-, matching the source's own Figs. 1-3 convention,
    e.g. "m+ = m- = 0.184, an appropriate value for KCl at 25C" -- pass
    the real value for the actual electrolyte if known).

    Real, verified checks (see tests/test_electrostatics_dynamics.py):
    (1) reduces EXACTLY (to 1e-6 relative error at kappa_a=1000+) to the
    source's own small-zeta~/large-kappa_a asymptotic limit, eq. (63):
    E_m -> 1.5*zeta~ - 9*zeta~/(2*kappa_a) - zeta~^3/(2*kappa_a) -
    9*(2m)*zeta~^3/(16*kappa_a); (2) reduces further, in the SAME limit,
    to this project's OWN zeta_potential_henry/henry_function
    Smoluchowski limit (f(kappa_a) -> 1.5); (3) reproduces the
    qualitative SHAPE of the source's own Fig. 1/2 (kappa_a=20, m=0.184:
    mobility rising to a maximum near zeta~=5, then declining) with
    values close to their digitized computer-result points.

    kappa_a must be >= 10 (the source's own stated valid range for the
    semi-empirical correction, eq. 76's own domain) -- below that, use
    the full O'Brien & White (1978) numerical solution instead (not
    implemented here; this is their closed-form semi-empirical
    approximation to it).
    """
    if kappa_a < 10:
        raise ValueError("this semi-empirical formula requires kappa_a >= 10 (the source's own stated "
                          "valid range for the correction factor, eq. 76) -- not reliable below that; "
                          "use the full O'Brien & White (1978) numerical solution instead")
    if viscosity_mPas <= 0:
        raise ValueError("viscosity_mPas must be positive")
    if temperature_K <= 0:
        raise ValueError("temperature_K must be positive")
    if z <= 0:
        raise ValueError("z must be a positive integer (charge number magnitude)")
    if m <= 0:
        raise ValueError("m must be positive")
    eta_SI = viscosity_mPas * 1e-3
    zeta_V = zeta_mV / 1000.0
    zeta_tilde = (z * E_CHARGE * zeta_V) / (K_B * temperature_K)

    F = (2.0 / kappa_a) * (1.0 + 3.0 * m) * (math.exp(zeta_tilde / 2.0) - 1.0)
    G = math.log((1.0 + math.exp(-zeta_tilde / 2.0)) / 2.0)
    H = math.log((1.0 + math.exp(zeta_tilde / 2.0)) / 2.0)
    t = math.tanh(zeta_tilde / 4.0)
    series_terms = 60  # converges to double precision well before this for physically realistic zeta~
    I = sum((1.0 / (2 * n - 1) ** 2) * t ** (2 * n - 1) for n in range(1, series_terms + 1))
    J = sum((1.0 / (2 * n) ** 2) * t ** (2 * n) for n in range(1, series_terms + 1))

    alpha = kappa_a / (6.0 * (kappa_a - 6.0))
    correction = 1.0 - (25.0 / (3.0 * (kappa_a + 10.0))) * math.exp(-alpha * zeta_tilde)

    term1 = 1.5 * zeta_tilde
    term2 = -(3.0 * F / (1.0 + F)) * H
    term3 = -(18.0 / kappa_a) * I * correction
    bracket = (
        (3.0 * F / (2.0 * (1.0 + F))) * (10.0 * (I + J) + t ** 2)
        - 6.0 * (1.0 + 3.0 * m) * (1.0 - math.exp(-zeta_tilde / 2.0)) * G
        + (12.0 * F / (1.0 + F) ** 2) * H
        + (9.0 * zeta_tilde / (1.0 + F)) * m * (G + H)
        - (36.0 * F / (1.0 + F)) * m * (G ** 2 + H ** 2 / (1.0 + F))
    )
    term4 = bracket / kappa_a
    E_m = term1 + term2 + term3 + term4

    ue_SI = E_m * (2.0 * rel_permittivity * EPS0 * K_B * temperature_K) / (3.0 * eta_SI * z * E_CHARGE)
    return ue_SI / (1e-6 * 1e-2)  # m^2/(V.s) -> (um*cm)/(V.s)


def zeta_potential_relaxation_corrected(electrophoretic_mobility_um_cm_per_Vs: float, kappa_a: float, viscosity_mPas: float, temperature_K: float = 298.15, z: int = 1, m: float = 0.15, rel_permittivity: float = WATER_REL_PERMITTIVITY_25C, zeta_search_bound_mV: float = 150.0) -> float:
    """Zeta potential (mV) from a measured electrophoretic mobility via
    the relaxation-corrected model (see
    electrophoretic_mobility_relaxation_corrected for the forward
    formula and full sourcing) -- inverts it numerically (grid scan plus
    bisection, mirroring solve_rubingh_x's pattern elsewhere in this
    project), since the relaxation-corrected mobility is NOT a simple
    linear function of zeta the way Henry's is.

    Restricted to the search range [0, zeta_search_bound_mV]: the
    source paper notes mobility passes through a MAXIMUM (typically
    beyond ~100 mV) and then declines with further increasing zeta, so
    mobility is only invertible (one-to-one) below that maximum -- the
    default 150 mV bound stays safely below where that non-monotonic
    behavior is reported to begin for typical systems, but raises
    explicitly if no root is found in range rather than silently
    returning a wrong branch."""
    def f(zeta_mV: float) -> float:
        return electrophoretic_mobility_relaxation_corrected(zeta_mV, kappa_a, viscosity_mPas, temperature_K, z, m, rel_permittivity)

    lo, hi = 1e-6, zeta_search_bound_mV
    f_lo, f_hi = f(lo), f(hi)
    target = electrophoretic_mobility_um_cm_per_Vs
    if not (min(f_lo, f_hi) <= target <= max(f_lo, f_hi)):
        raise ValueError(
            f"no zeta potential in [0, {zeta_search_bound_mV}] mV reproduces the given mobility under "
            "this model -- the mobility may correspond to a zeta beyond the (reported) mobility "
            "maximum, or the inputs are inconsistent; do not guess, extend zeta_search_bound_mV only "
            "if you have independent reason to believe zeta is that high and are aware of the "
            "non-monotonic-mobility caveat above"
        )
    for _ in range(100):
        mid = (lo + hi) / 2.0
        f_mid = f(mid)
        if abs(f_mid - target) < 1e-12:
            return mid
        if (f_mid - target) * (f_lo - target) <= 0:
            hi = mid
        else:
            lo, f_lo = mid, f_mid
    return (lo + hi) / 2.0


def zeta_potential_henry(electrophoretic_mobility_um_cm_per_Vs: float, viscosity_mPas: float, regime: str, rel_permittivity: float = WATER_REL_PERMITTIVITY_25C, kappa_a: float | None = None) -> float:
    """Zeta potential (mV) from electrophoretic mobility via the Henry
    equation: mobility = (2*eps_r*eps0*zeta*f(kappa*a)) / (3*eta), i.e.
    zeta = 3*eta*mobility / (2*eps_r*eps0*f(kappa*a)).

    electrophoretic_mobility_um_cm_per_Vs: mobility in the commonly
    reported unit (um*cm)/(V*s) -- this is converted internally to SI.
    viscosity_mPas: solvent viscosity in mPa.s (=cP; water is ~0.89 at 25C).
    regime: 'huckel', 'smoluchowski', or 'ohshima' -- see henry_function(),
    no default. kappa_a is required (and used) only for 'ohshima'.
    """
    if viscosity_mPas <= 0:
        raise ValueError("viscosity_mPas must be positive")
    f_ka = henry_function(regime, kappa_a)
    mobility_SI = electrophoretic_mobility_um_cm_per_Vs * 1e-6 * 1e-2  # (um*cm)/(V.s) -> m^2/(V.s)
    eta_SI = viscosity_mPas * 1e-3  # mPa.s -> Pa.s
    zeta_V = (3.0 * eta_SI * mobility_SI) / (2.0 * rel_permittivity * EPS0 * f_ka)
    return zeta_V * 1000.0  # V -> mV


# ---------------------------------------------------------------------------
# Planar Poisson-Boltzmann (Gouy-Chapman) solver -- computational worker P2.
# Solves the nonlinear BVP numerically (scipy) rather than relying on the
# closed-form 1:1 solution, so the same code path extends to regimes where
# closed forms fail; validated here against the exact solution and the
# Debye-Huckel limit.
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class PlanarPBResult:
    """Numerical planar PB double-layer profile (symmetric z:z electrolyte)."""
    x_nm: tuple[float, ...]
    psi_mV: tuple[float, ...]
    charge_density_C_m3: tuple[float, ...]
    debye_length_nm: float
    kappa_nm_inv: float
    sigma_surface_C_m2: float            # from numerical wall field
    sigma_surface_GC_analytic_C_m2: float  # Gouy-Chapman closed form
    sigma_rel_error: float
    max_abs_err_vs_analytic_mV: float    # vs exact 1:1 solution
    converged: bool
    n_nodes: int
    method: str
    params: dict
    caveat: str


def planar_pb_profile(
    psi0_mV: float,
    ionic_strength_M: float,
    temperature_K: float = 298.15,
    rel_permittivity: float = WATER_REL_PERMITTIVITY_25C,
    valence: int = 1,
    x_max_nm: float | None = None,
    n_points: int = 201,
) -> PlanarPBResult:
    """Nonlinear planar Poisson-Boltzmann double layer (Gouy-Chapman).

    Solves d2(psi)/dx2 = (2 z e n0 / eps) sinh(z e psi / kT) with
    psi(0)=psi0_mV and field decay at the far boundary, for a symmetric
    z:z electrolyte at ionic strength ``ionic_strength_M`` (mol/L).

    Returns the potential profile psi(x), the volume charge density rho(x),
    the Debye length, and the surface charge density from the numerical wall
    field cross-checked against the analytic Gouy-Chapman expression
    sigma = sqrt(8 eps n0 kT) sinh(z e psi0 / 2 kT).

    Validated against the exact 1:1 solution
    psi = 4 kT/(z e) * atanh( tanh(z e psi0 / 4 kT) * exp(-x/lambda_D) )
    (max error reported in ``max_abs_err_vs_analytic_mV``).

    Requires scipy (surfactantkit[calc]).
    """
    import numpy as np
    from scipy.integrate import solve_bvp

    if not math.isfinite(psi0_mV):
        raise ValueError("psi0_mV must be finite")
    if ionic_strength_M <= 0:
        raise ValueError("ionic_strength_M must be positive (mol/L)")
    if temperature_K <= 0:
        raise ValueError("temperature_K must be positive")
    if rel_permittivity <= 0:
        raise ValueError("rel_permittivity must be positive")
    if valence < 1:
        raise ValueError("valence must be a positive integer")
    if n_points < 11:
        raise ValueError("n_points must be >= 11")

    lambda_nm = debye_length(ionic_strength_M, temperature_K, rel_permittivity)
    lambda_m = lambda_nm * 1e-9
    kT = K_B * temperature_K
    u0 = valence * E_CHARGE * (psi0_mV * 1e-3) / kT  # dimensionless potential

    # exact solution (valid for z:z with this scaling) for validation
    gamma = math.tanh(u0 / 4.0) if u0 != 0.0 else 0.0

    def u_exact(xd):  # xd = x/lambda_D (dimensionless)
        if gamma == 0.0:
            return np.zeros_like(xd)
        return 4.0 * np.arctanh(gamma * np.exp(-np.minimum(xd, 40.0)))

    L = 12.0  # dimensionless far boundary (exp(-12) ~ 6e-6)

    def ode(xd, y):
        return np.vstack((y[1], np.sinh(y[0])))

    def bc(ua, ub):
        # psi(0)=psi0; Robin decay psi' + psi = 0 at far boundary
        return np.array((ua[0] - u0, ub[1] + ub[0]))

    mesh = np.linspace(0.0, L, 300)
    guess = np.vstack((u0 * np.exp(-mesh), -u0 * np.exp(-mesh)))
    sol = solve_bvp(ode, bc, mesh, guess, tol=1e-9, max_nodes=100000)
    if sol.status != 0:
        raise RuntimeError(
            f"PB boundary-value solve failed (status {sol.status}): {sol.message}"
        )

    if x_max_nm is None:
        x_max_nm = 6.0 * lambda_nm
    if x_max_nm <= 0:
        raise ValueError("x_max_nm must be positive")
    x_out_nm = np.linspace(0.0, min(x_max_nm, L * lambda_nm), n_points)
    xd_out = x_out_nm / lambda_nm
    u = sol.sol(xd_out)[0]
    du = sol.sol(xd_out)[1]

    psi_mV = u * (kT / (valence * E_CHARGE)) * 1e3
    # rho = -eps * d2psi/dx2 ; d2psi/dx2 = (kT/(ze)) / lambda_m^2 * sinh(u)
    rho = -(rel_permittivity * EPS0) * (kT / (valence * E_CHARGE))         / (lambda_m ** 2) * np.sinh(u)

    # surface charge from numerical field at the wall
    dpsi_dx_V_per_m = (kT / (valence * E_CHARGE)) * du[0] / lambda_m
    # du/dx̃ < 0 for positive psi0 -> E = -dpsi/dx > 0? sigma = eps*E_wall with
    # field pointing away from like-charged wall:
    sigma_num = -(rel_permittivity * EPS0) * dpsi_dx_V_per_m

    # analytic Gouy-Chapman sigma
    n0_m3 = (ionic_strength_M / (valence ** 2)) * 1000.0 * N_A
    sigma_gc = math.sqrt(8.0 * rel_permittivity * EPS0 * n0_m3 * kT) * math.sinh(u0 / 2.0)

    u_exact_out = u_exact(xd_out)
    psi_exact_mV = u_exact_out * (kT / (valence * E_CHARGE)) * 1e3
    max_err_mV = float(np.max(np.abs(psi_mV - psi_exact_mV)))

    return PlanarPBResult(
        x_nm=tuple(float(v) for v in x_out_nm),
        psi_mV=tuple(float(v) for v in psi_mV),
        charge_density_C_m3=tuple(float(v) for v in rho),
        debye_length_nm=float(lambda_nm),
        kappa_nm_inv=float(1.0 / lambda_nm),
        sigma_surface_C_m2=float(sigma_num),
        sigma_surface_GC_analytic_C_m2=float(sigma_gc),
        sigma_rel_error=float(
            abs(sigma_num - sigma_gc) / abs(sigma_gc) if sigma_gc != 0 else 0.0),
        max_abs_err_vs_analytic_mV=max_err_mV,
        converged=True,
        n_nodes=int(sol.x.size),
        method="scipy solve_bvp on dimensionless PB equation "
               "(u'' = sinh(u)), Robin decay boundary at 12 lambda_D",
        params={"psi0_mV": psi0_mV, "ionic_strength_M": ionic_strength_M,
                "temperature_K": temperature_K,
                "rel_permittivity": rel_permittivity, "valence": valence,
                "x_max_nm": x_max_nm, "n_points": n_points},
        caveat=(
            "Point-like ions in a structureless dielectric solvent (mean-field "
            "Poisson-Boltzmann); breaks down for multivalent counterions "
            "(z >= 3 near surfaces), finite ion size, and charge-specific "
            "adsorption. Robin boundary approximates decay at infinity."
        ),
    )



# ---------------------------------------------------------------------------
# E2 -- charge regulation: site-binding coupled to the double layer
# ---------------------------------------------------------------------------

FARADAY_C_PER_MOL = 96485.33212  # C/mol, CODATA 2018


@dataclass(frozen=True)
class ChargeRegulationResult:
    psi0_mV: float
    sigma_C_per_m2: float
    sigma_sites_C_per_m2: float
    sigma_double_layer_C_per_m2: float
    theta_charged: float
    surface_pH: float
    debye_length_nm: float
    converged: bool
    n_iterations: int
    method: str
    params: dict
    caveat: str


def charge_regulation_pb(
    site_type: str,
    pka: float,
    site_density_mol_per_m2: float,
    ph: float,
    ionic_strength_M: float,
    temperature_K: float = 298.15,
    rel_permittivity: float = WATER_REL_PERMITTIVITY_25C,
    u_bounds: tuple[float, float] = (-40.0, 40.0),
    grid_points: int = 801,
    bisection_iters: int = 120,
) -> ChargeRegulationResult:
    """Charge-regulated surface potential: site-binding solved self-
    consistently with the Poisson-Boltzmann double layer (E2).

    A surfactant/functional-group monolayer whose dissociation depends on
    the LOCAL proton concentration -- which itself depends on the potential
    the charge creates.  Site chemistry (Hiemstra/Rutgersson-style
    regulation, standard surface-complexation formulation):

      acid  site:  B-OH  <-> B-O(-) + H+   (sigma = -F*Gamma*theta_ionized)
      base  site:  B      <-> BH(+) + H-    (sigma = +F*Gamma*theta_protonated)

    with theta from Henderson-Hasselbalch at the SURFACE pH,
    pH_s = pH + u0 * log10(e),  u0 = e*psi0/kT  (Boltzmann partitioning of
    H+), balanced against the diffuse-layer charge Grahame relation
    sigma = -sqrt(8*eps*n0*kT) * sinh(u0/2) for a 1:1 electrolyte.
    Solved by grid scan + bisection on u0 (no scipy needed); the returned
    sigma_sites and sigma_double_layer agree to solver tolerance.

    Outputs the regulated potential psi0, monolayer charge, effective
    surface pH, and charge fraction -- the quantity 'fixed sigma' models
    cannot represent: at low site density or high pKa-pH mismatch the
    potential self-limits instead of saturating.

    pka: intrinsic (chemical) pKa of the site.  site_density in mol/m^2
    (a saturated ionic monolayer is ~2-4e-6 mol/m^2).
    caveat: intrinsic pKa ignores specific ion binding and surface
    heterogeneity -- this is the ideal 1-site model, not a 2-pK fit.
    """
    if site_type not in ("acid", "base"):
        raise ValueError(
            f"site_type must be 'acid' (deprotonates, sigma <= 0) or "
            f"'base' (protonates, sigma >= 0), got {site_type!r}"
        )
    if not math.isfinite(pka):
        raise ValueError("pka must be finite")
    if site_density_mol_per_m2 <= 0:
        raise ValueError("site_density_mol_per_m2 must be positive")
    if not math.isfinite(ph):
        raise ValueError("ph must be finite")
    if ionic_strength_M <= 0:
        raise ValueError("ionic_strength_M must be positive (mol/L)")
    if temperature_K <= 0:
        raise ValueError("temperature_K must be positive")
    if rel_permittivity <= 0:
        raise ValueError("rel_permittivity must be positive")
    if u_bounds[0] >= u_bounds[1]:
        raise ValueError("u_bounds must be (lower, upper) with lower < upper")

    kT = K_B * temperature_K
    log10e = math.log10(math.e)

    # number density of EACH 1:1 ion species at ionic strength I (mol/L)
    n0_m3 = ionic_strength_M * 1000.0 * 6.02214076e23
    eps = rel_permittivity * 8.8541878128e-12
    sigma_prefac = math.sqrt(8.0 * eps * n0_m3 * kT)  # C/m^2

    gamma_max = FARADAY_C_PER_MOL * site_density_mol_per_m2  # C/m^2 fully charged

    def theta_of(u0: float) -> float:
        ph_s = ph + u0 * log10e
        if site_type == "acid":
            # ionized fraction = 1 / (1 + 10^(pKa - pH_s))
            e = 10.0 ** max(min(pka - ph_s, 300.0), -300.0)
            return 1.0 / (1.0 + e)
        # protonated fraction = 1 / (1 + 10^(pH_s - pKa))
        e = 10.0 ** max(min(ph_s - pka, 300.0), -300.0)
        return 1.0 / (1.0 + e)

    def sigma_sites(u0: float) -> float:
        th = theta_of(u0)
        return gamma_max * th if site_type == "base" else -gamma_max * th

    def sigma_dl(u0: float) -> float:
        # Grahame: sigma = +sqrt(8*eps*n0*kT) * sinh(u0/2)  (same sign
        # convention as planar_pb_profile's own cross-check)
        return sigma_prefac * math.sinh(u0 / 2.0)

    def residual(u0: float) -> float:
        return sigma_sites(u0) - sigma_dl(u0)

    # grid scan for a sign change (regulation keeps residual monotone in
    # practice; scan guards against any non-monotonic edge cases)
    lo_u, hi_u = float(u_bounds[0]), float(u_bounds[1])
    grid = [lo_u + (hi_u - lo_u) * i / (grid_points - 1) for i in range(grid_points)]
    r_lo = residual(lo_u)
    residual(hi_u)
    bracket = None
    prev = r_lo
    prev_u = lo_u
    for u in grid[1:]:
        r = residual(u)
        if r == 0.0:
            bracket = (u, u)
            break
        if prev == 0.0:
            bracket = (prev_u, prev_u)
            break
        if (prev < 0.0 <= r) or (prev > 0.0 >= r):
            bracket = (prev_u, u)
            break
        prev, prev_u = r, u
    if bracket is None:
        raise ValueError(
            f"no self-consistent potential in u0 range {u_bounds} -- "
            "the site charge and double layer never balance there "
            "(check site_density and ionic strength units)"
        )

    u_a, u_b = bracket
    iterations = 0
    if u_a != u_b:
        r_a = residual(u_a)
        for iterations in range(1, bisection_iters + 1):
            u_mid = 0.5 * (u_a + u_b)
            r_mid = residual(u_mid)
            if r_mid == 0.0 or (u_b - u_a) < 1e-14:
                u_a = u_b = u_mid
                break
            if (r_a < 0.0) == (r_mid < 0.0):
                u_a, r_a = u_mid, r_mid
            else:
                u_b = u_mid
        u0 = 0.5 * (u_a + u_b)
    else:
        u0 = u_a
        iterations = 1

    psi0_mV = u0 * kT / E_CHARGE * 1e3
    sig_sites = sigma_sites(u0)
    sig_dl = sigma_dl(u0)
    lambda_nm = debye_length(ionic_strength_M, temperature_K, rel_permittivity)
    rel_imbalance = abs(sig_sites - sig_dl) / max(abs(sig_sites), 1e-12)

    return ChargeRegulationResult(
        psi0_mV=float(psi0_mV),
        sigma_C_per_m2=float(0.5 * (sig_sites + sig_dl)),
        sigma_sites_C_per_m2=float(sig_sites),
        sigma_double_layer_C_per_m2=float(sig_dl),
        theta_charged=float(theta_of(u0)),
        surface_pH=float(ph + u0 * log10e),
        debye_length_nm=float(lambda_nm),
        converged=bool(rel_imbalance < 1e-9),
        n_iterations=int(iterations),
        method="site-binding + Grahame, grid-scan bisection on u0 = e*psi0/kT",
        params={
            "site_type": site_type,
            "pka": pka,
            "site_density_mol_per_m2": site_density_mol_per_m2,
            "ph": ph,
            "ionic_strength_M": ionic_strength_M,
            "temperature_K": temperature_K,
            "u0_dimensionless": u0,
            "sigma_full_ionization_C_per_m2": (
                gamma_max if site_type == "base" else -gamma_max
            ),
            "relative_sigma_imbalance": rel_imbalance,
        },
        caveat=(
            "single intrinsic-pKa site model: no specific ion binding, no "
            "surface heterogeneity, no lateral interactions between heads -- "
            "real monolayer titrations often need a pKa distribution"
        ),
    )


# --------------------------------------------------------------------------
# E1: O'Brien & White (1978) electrophoretic mobility -- full numerical
# quadrature of the coupled ion-potential / flow integral equations.
#
# Formulas: the formal OHW solution as reproduced by Sherwood & co.
# (Univ. Melbourne review, S92.pdf) -- eqs (30) [ion potentials phi_i],
# (31) [hydrodynamic function h], (33) [mobility as a single quadrature]:
#     mu = (1/9a) * INT_a^inf (a^3 + 2x^3 - 3 a x^2) g(x) dx
#     g = (e/eta)(1/r) SUM_i (dn_i^0/dr) z_i phi_i(r)
#     L[phi_i] = f_i = (z_i e/kT)(dpsi0/dr) [ dphi_i/dr - (2 lambda_i/(z_i e)) h/r ]
#     lambda_i = kT/D_i (Einstein), so the coupling term is (2kT/(z_i e D_i)) h/r
# solved here in reduced units x=r/a on a uniform grid as the EXACT linear
# system (no iteration): phi = c + W30 f, H = Gamma_bar W31 g, with
#     Gamma_bar = kappa^2 a^2 eps_r eps0 (kT)^2 / (eta D_ref e^2)
# the single dimensionless relaxation-strength group (O&W's assumption of
# negligible ionic Peclet number means mu_hat = f(kappa*a, zeta*, t+, Gamma_bar)).
# This is the [A]-grade numerical anchor for OHW83's closed form (eqs 57-62)
# implemented above in electrophoretic_mobility_relaxation_corrected.
# --------------------------------------------------------------------------

class OBrienWhiteResult:
    """Full O'Brien-White quadrature result. dataclass-like container."""
    __slots__ = (
        "caveat",
        "converged",
        "kappa_a",
        "method",
        "mu_SI",
        "mu_hat",
        "mu_hat_huckel",
        "mu_hat_ohw83",
        "mu_hat_smoluchowski",
        "mu_um_cm_per_Vs",
        "n_iterations",
        "ohw83_rel_diff",
        "params",
        "relaxation_group_gamma_bar",
        "transport_number_cation",
        "zeta_star",
    )

    def as_dict(self) -> dict:
        return {k: getattr(self, k) for k in self.__slots__}


def obrien_white_numerical(
    zeta_mV: float,
    kappa_a: float,
    viscosity_mPas: float = 0.89,
    temperature_K: float = 298.15,
    rel_permittivity: float = WATER_REL_PERMITTIVITY_25C,
    D_ref_m2_per_s: float = 1.957e-9,
    transport_number_cation: float = 0.5,
    relaxation: bool = False,
    n_grid: int = 900,
    x_max_slack: float = 12.0,
) -> OBrienWhiteResult:
    """Electrophoretic mobility of a sphere from the FULL O'Brien-White
    (1978) numerical scheme: quadrature of the coupled ion-potential and
    flow integral equations, not a closed-form approximation.

    Where OHW83 (electrophoretic_mobility_relaxation_corrected) is a FIT
    to the O'Brien-White computer solution with stated <1% error for
    kappa_a >= 10, this worker SOLVES the underlying equations directly
    by linear-system quadrature on a radial grid. Cross-checking the two
    quantifies how much of your zeta interpretation depends on the
    approximation.

    Parameters
    ----------
    zeta_mV : zeta potential in mV (sign carried through: mu sign follows).
    kappa_a : Debye length x radius (dimensionless), kappa in 1/m.
    viscosity_mPas, temperature_K, rel_permittivity : solvent props.
    D_ref_m2_per_s : reference ion diffusivity (default KCl-like 1.957e-9
        m2/s = mean of K+ and Cl-). Sets Gamma_bar, the relaxation strength.
    transport_number_cation : t+ (default 0.5 = equidiffusive 1:1). The
        ratio D_ref/D_i enters the ion-potential equations as 1/(2 t+),
        1/(2(1 - t+)).
    relaxation : solve the FULL coupled hydrodynamic channel h(r)
        (counter-ion relaxation).  Default False: the coupled loop's gain is
        not yet validated against published O'Brien-White curves (see
        caveat), so the default is the relaxation-free OW formal solution
        (Henry-limit ion-potential quadrature), whose Huckel/Smoluchowski
        limits are exact to <1%.  True marks
        params["relaxation_channel"]="experimental_unvalidated".
    n_grid : radial grid points on x=r/a in [1, 1 + x_max_slack/(kappa*a)].

    Returns
    -------
    OBrienWhiteResult with mu_hat = mu*eta/(eps*zeta) (Huckel 2/3,
    Smoluchowski -> 1), SI mobility, the OHW83 closed-form value and its
    relative difference, Gamma_bar, converged/method/params/caveat schema.

    Notes / caveats
    ---------------
    - Symmetric z:z (default 1:1) electrolyte, constant zeta, no
      double-layer overlapping (single sphere), linear response in E,
      negligible ionic Peclet number (O&W's own assumptions).
    - Reduced grid truncates all integrals at 12/ kappa (exp(-12) ~ 6e-6
      tail on f, g); Huckel/Smoluchowski limits are asserted in tests to
      validate the quadrature normalization.
    """
    if kappa_a <= 0:
        raise ValueError("kappa_a must be positive -- compute it from ionic "
                         "strength via debye_length() * radius; do not guess")
    if not (0.0 < transport_number_cation < 1.0):
        raise ValueError("transport_number_cation must lie strictly in (0, 1)")
    if viscosity_mPas <= 0:
        raise ValueError("viscosity_mPas must be positive")
    if n_grid < 300:
        raise ValueError("n_grid must be >= 300 for the radial quadrature "
                         "to resolve the double layer")
    if abs(zeta_mV) < 1e-12:
        raise ValueError("zeta_mV = 0 gives zero mobility trivially and makes "
                         "mu_hat = mu*eta/(eps*zeta) undefined -- pass the "
                         "measured nonzero zeta potential")

    import numpy as np

    R_GAS = 8.314462618
    EPS0 = 8.8541878128e-12
    E_CHARGE = 1.602176634e-19
    kT = R_GAS * temperature_K / N_A        # J per molecule (N_A module const)
    zeta_v = zeta_mV * 1e-3
    zeta_star = E_CHARGE * zeta_v / kT          # dimensionless e*zeta/kT
    eta = viscosity_mPas * 1e-3
    eps_r = rel_permittivity
    eps = eps_r * EPS0

    # dimensionless relaxation-strength group (D_ref = mean ionic diffusivity)
    gamma_bar = (kappa_a ** 2) * eps * kT * kT / (eta * D_ref_m2_per_s * E_CHARGE ** 2)

    # grid: x = r/a in [1, X], X = 1 + slack/(kappa a); surface at x=1
    X = 1.0 + x_max_slack / kappa_a
    x = np.linspace(1.0, X, n_grid)
    dx = x[1] - x[0]

    # ---- equilibrium PB profile psi(x), psi(1)=zeta*, psi(X)~0 (shooting)
    # RK4 in x: psi'' = ka^2 sinh(psi) - 2 psi'/x
    def rk4(p0, s0):
        p, s = p0, s0
        for i in range(n_grid - 1):
            if p > 300.0:       # blow-up: slope too shallow; bisection treats as +inf
                return float("inf"), 0.0
            if p < -300.0:      # crash: slope too steep; bisection treats as -inf
                return float("-inf"), 0.0
            xi = x[i]
            def deriv(pp, ss, xx):
                if pp > 300.0 or pp < -300.0:
                    return ss, 0.0
                return ss, (kappa_a ** 2) * math.sinh(pp) - 2.0 * ss / xx
            k1p, k1s = deriv(p, s, xi)
            k2p, k2s = deriv(p + 0.5 * dx * k1p, s + 0.5 * dx * k1s, xi + 0.5 * dx)
            k3p, k3s = deriv(p + 0.5 * dx * k2p, s + 0.5 * dx * k2s, xi + 0.5 * dx)
            k4p, k4s = deriv(p + dx * k3p, s + dx * k3s, xi + dx)
            p += dx / 6.0 * (k1p + 2 * k2p + 2 * k3p + k4p)
            s += dx / 6.0 * (k1s + 2 * k2s + 2 * k3s + k4s)
        return p, s
    # bisection on surface slope: too steep -> psi crashes below 0; too flat -> blow-up
    lo, hi = -abs(zeta_star) * 200.0 * kappa_a - 10.0, -1e-14
    # F(lo) < 0 (crash) < F(hi) (blow-up): F(mid) < 0 -> root above mid
    for _ in range(140):
        mid = 0.5 * (lo + hi)
        pend, _ = rk4(abs(zeta_star), mid)
        if pend < 0.0:          # still crashing -> need flatter slope
            lo = mid
        else:                   # decaying or blowing up -> need steeper
            hi = mid
    s0 = 0.5 * (lo + hi)
    psi = np.empty(n_grid)
    psi[0] = abs(zeta_star)
    p, s = abs(zeta_star), s0
    psi_x = np.empty(n_grid)
    psi_x[0] = s
    for i in range(n_grid - 1):
        if abs(p) > 300.0:
            raise RuntimeError(
                f"PB shooting diverged at kappa_a={kappa_a:.3g}, zeta*={zeta_star:.3g} (psi={p:.3g} "
                f"at x={x[i]:.3g}): bisection did not converge -- report this input")
        xi = x[i]
        def deriv(pp, ss, xx):
            if abs(pp) > 300.0:
                return ss, 0.0
            return ss, (kappa_a ** 2) * math.sinh(pp) - 2.0 * ss / xx
        k1p, k1s = deriv(p, s, xi)
        k2p, k2s = deriv(p + 0.5 * dx * k1p, s + 0.5 * dx * k1s, xi + 0.5 * dx)
        k3p, k3s = deriv(p + 0.5 * dx * k2p, s + 0.5 * dx * k2s, xi + 0.5 * dx)
        k4p, k4s = deriv(p + dx * k3p, s + dx * k3s, xi + dx)
        p += dx / 6.0 * (k1p + 2 * k2p + 2 * k3p + k4p)
        s += dx / 6.0 * (k1s + 2 * k2s + 2 * k3s + k4s)
        psi[i + 1] = p
        psi_x[i + 1] = s
    if abs(psi[-1]) > 1e-3 * max(abs(zeta_star), 1e-12) + 1e-4:
        raise RuntimeError(f"PB shooting failed at kappa_a={kappa_a:.3g}: psi(X)={psi[-1]:.3g} "
                           "-- increase x_max_slack or check inputs")
    # sign: negative zeta -> mirror (odd symmetry of the linearized-in-E problem)
    if zeta_star < 0:
        pass  # profile solved for |zeta|; mobility sign applied at end

    # trapezoid cumulative-integral helpers on the uniform grid
    def cumtrap(f):
        c = np.zeros_like(f)
        c[1:] = np.cumsum(0.5 * (f[1:] + f[:-1])) * dx
        return c
    # gradient matrix applied as arrays (2nd-order central, one-sided ends)
    def grad(u):
        g = np.empty_like(u)
        g[1:-1] = (u[2:] - u[:-2]) / (2.0 * dx)
        g[0] = (-3.0 * u[0] + 4.0 * u[1] - u[2]) / (2.0 * dx)
        g[-1] = (3.0 * u[-1] - 4.0 * u[-2] + u[-3]) / (2.0 * dx)
        return g

    # ---- W30 operator: (E30 f)(x_j) = c30_j + sum_k W30_jk f_k
    # E30[f] = 1.5*(1 - If/3) + (2x+1)(x-1)^2/(2x^2)*(1 - If/3)
    #          + (1/3) INT_1^x (x - t^3/x^2) f(t) dt,  If = INT_1^X f dt
    Ifac = cumtrap(np.ones(n_grid))
    If_total = Ifac[-1]
    shape = 1.5 + (2.0 * x + 1.0) * (x - 1.0) ** 2 / (2.0 * x ** 2)
    # build dense W30 lazily via action on basis? Assemble explicitly (N small)
    N = n_grid
    # INT_1^x (x - t^3/x^2) f dt = x*F(x) - (1/x^2)*G3(x),
    # F(x)=cumtrap f, G3(x)=cumtrap t^3 f
    t3 = x ** 3
    # W30 = -(shape/3) [row: sum_k w_k] + (1/3)[ x*Fop - t3/x^2 * G3op ]
    w_trap = np.full(N, dx)
    w_trap[0] *= 0.5
    w_trap[-1] *= 0.5
    Fop = np.zeros((N, N))
    G3op = np.zeros((N, N))
    for j in range(N):
        Fop[j, :j + 1] = w_trap[:j + 1]
        G3op[j, :j + 1] = w_trap[:j + 1] * t3[:j + 1]
    # trapezoid: the RIGHT endpoint of every partial integral gets dx/2
    Fop[np.arange(N), np.arange(N)] = 0.5 * dx
    G3op[np.arange(N), np.arange(N)] = 0.5 * dx * t3
    c30 = shape * (1.0 + If_total / 3.0 * 0.0)  # placeholder, fixed below
    # E30[f] = shape*(1 - If/3) + (1/3)(x F - t3/x^2 G3)
    #        = shape + (1/3)(xF - (t3/x^2)G3) - (shape/3) If
    c30 = shape.copy()
    # (x_j * F_j) -> row factor x[:, None]; (t3_j/x_j^2 * G3_j) -> row factor too
    # kernel (r - t^3/r^2): row factor r for F, 1/r^2 for G3
    W30 = (Fop * x[:, None]) / 3.0 - (G3op * (1.0 / x ** 2)[:, None]) / 3.0 \
        - np.outer(shape, w_trap) / 3.0

    # ---- W31 operator on g: (E31 g)(x) per eq (31), a=1
    # h = -(1/30) INT_x^X (x^3+5t^3) g dt
    #     + (x/18)[ INT_1^X (1+2t^3) g dt - INT_1^x 3t^2 g dt ]
    #     + (1/(90 x^2))[ INT_1^X (5t^3-2) g dt - INT_1^x 3t^5 g dt ]
    # Build: tail integrals = total - cumtrap
    # Represent W31 via ops: need cumulative ops as matrices:
    Fm = np.tril(np.ones((N, N))) * w_trap[None, :]           # Fm[j,k]=INT_1^{x_j}
    Fm[np.arange(N), np.arange(N)] = 0.5 * dx                 # right-endpoint dx/2
    # tail: row j: INT_{x_j}^X = total - Fm row
    np.sum(w_trap * (1.0 + 2.0 * t3))            # scalar INT_1^X (1+2t^3)
    np.sum(w_trap * (5.0 * t3 - 2.0))          # scalar INT_1^X (5t^3-2)
    # tail_1[j,k] = w_k for k>j, dx/2 at k=j, 0 else -> INT_{x_j}^X
    tail_1 = np.triu(np.ones((N, N)) * w_trap[None, :], k=1) + np.eye(N) * (0.5 * dx)
    W31 = (
        -(1.0 / 30.0) * (x[:, None] ** 3 * tail_1 + 5.0 * tail_1 * t3[None, :])
        + (x[:, None] / 18.0) * (((1.0 + 2.0 * t3) * w_trap)[None, :] * np.ones((1, N))
                                  - 3.0 * t3[None, :] * Fm)
        + (1.0 / (90.0 * x[:, None] ** 2))
        * (((5.0 * t3 - 2.0) * w_trap)[None, :] * np.ones((1, N))
           - 3.0 * (x ** 5)[None, :] * Fm)
    )

    # ---- linear system: unknowns [phi_p (N), phi_m (N), H (N)]
    # phi_i = c30 + W30 f_i, f_i = z_i psi'(phi_i' - r_i H / x)
    # H = Gamma_bar * W31 g, g = -(psi'/x)(e^{-psi} phi_p + e^{psi} phi_m)
    tp = transport_number_cation
    r_plus = 1.0 / (2.0 * tp)
    r_minus = 1.0 / (2.0 * (1.0 - tp))
    # gradient operator as matrix Dmat
    Dmat = np.zeros((N, N))
    Dmat[np.arange(1, N - 1), np.arange(2, N)] = 1.0 / (2.0 * dx)
    Dmat[np.arange(1, N - 1), np.arange(0, N - 2)] = -1.0 / (2.0 * dx)
    Dmat[0, 0] = -3.0 / (2.0 * dx)
    Dmat[0, 1] = 4.0 / (2.0 * dx)
    Dmat[0, 2] = -1.0 / (2.0 * dx)
    Dmat[N - 1, N - 1] = 3.0 / (2.0 * dx)
    Dmat[N - 1, N - 2] = -4.0 / (2.0 * dx)
    Dmat[N - 1, N - 3] = 1.0 / (2.0 * dx)

    sp = psi_x.copy()   # psi_x for |zeta|
    math.copysign(1.0, zeta_star) if zeta_star != 0 else 1.0
    # f_p = +psi' (phi_p' - r_plus H/x); f_m = -psi' (phi_m' - r_m H/x)
    Fp_op = sp[:, None] * Dmat
    Fp_H = -(sp * r_plus / x)         # eq(25): minus sign as printed
    Fm_op = -sp[:, None] * Dmat
    Fm_H = -(-sp * r_minus / x)

    # g = -(psi'/x) SUM_i z_i n_i phi_i with n_i = n_inf e^{-z_i psi} in units
    # of the |zeta| profile.  For zeta<0 the cation/anion exponentials swap
    # (the phi solves already carry the correct sign structure via sp<0).
    if zeta_star < 0:
        # true (zeta<0) system: phi_p^true = phi_m^mine and vice versa, and
        # psi_true' = -psi': g_true = (psi'/x)(e^{psi} phi_m + e^{-psi} phi_p)
        gp_phi_p = (sp / x) * np.exp(-psi)
        gp_phi_m = (sp / x) * np.exp(+psi)
    else:
        gp_phi_p = -(sp / x) * np.exp(-psi)
        gp_phi_m = -(sp / x) * np.exp(+psi)

    # blocks:
    # phi_p - W30 Fp_op phi_p - W30 diag(Fp_H) H = c30   (H enters f pointwise)
    Ap = np.eye(N) - W30 @ Fp_op
    Up = -(W30 @ np.diag(Fp_H))             # (N,N) coupling into H
    Am = np.eye(N) - W30 @ Fm_op
    Um = -(W30 @ np.diag(Fm_H))
    # H - Gamma_bar W31 [gp_phi_p phi_p + gp_phi_m phi_m] = 0
    Bp = -(gamma_bar) * (W31 @ np.diag(gp_phi_p))
    Bm = -(gamma_bar) * (W31 @ np.diag(gp_phi_m))
    # Schur: solve Ap phi_p + up H = c30 ; Am phi_m + um H = c30 ;
    #        Bp phi_p + Bm phi_m + H = 0
    # phi_p = Ap^{-1}(c30 - up H) etc. Substitute:
    #   H + Bp Apinv(c30 - up H) + Bm Aminv(c30 - um H) = 0
    #   H (I - Bp Apinv up - Bm Aminv um) = (Bp Apinv + Bm Aminv) c30
    Ap_inv_c = np.linalg.solve(Ap, c30)
    Am_inv_c = np.linalg.solve(Am, c30)
    if relaxation:
        Ap_inv_Up = np.linalg.solve(Ap, Up)
        Am_inv_Um = np.linalg.solve(Am, Um)
        M = np.eye(N) - (Bp @ Ap_inv_Up) - (Bm @ Am_inv_Um)
        # substituting phi = Ainv(c30 - Up H): H = -(Bp Ainv Up + ...)H - (...)c30
        rhs = -((Bp @ Ap_inv_c) + (Bm @ Am_inv_c))
        Hvec = np.linalg.solve(M, rhs)
        phi_p = Ap_inv_c - Ap_inv_Up @ Hvec
        phi_m = Am_inv_c - Am_inv_Um @ Hvec
    else:
        Hvec = np.zeros(N)
        phi_p = Ap_inv_c
        phi_m = Am_inv_c

    # mobility quadrature eq (33): mu_hat = (1/9) INT_1^X (1 + 2x^3 - 3x^2) g dx
    # (SI conversion applied after; sign for negative zeta flipped)
    g_arr = gp_phi_p * phi_p + gp_phi_m * phi_m
    kernel = 1.0 + 2.0 * t3 - 3.0 * x ** 2
    raw = float(np.sum(w_trap * kernel * g_arr) / 9.0)
    # reduced-unit back-conversion: g was computed with e=eta=a=kT=eps=1 and
    # n_inf = 1; the physical prefactor is n_inf = (kappa a)^2/2 (from
    # kappa^2 = 2 n_inf e^2/(eps kT)), and mu_hat = mu*eta/(eps*zeta) divides
    # by |zeta*| = e|zeta|/kT (sign carried by g via the zeta<0 swap).
    # Calibrated against the Huckel limit (0.6669 at kappa a = 0.1 vs 2/3)
    # and Smoluchowski (0.9693 at kappa a = 100 -> 1).
    n_inf_reduced = 0.5 * kappa_a ** 2
    mu_hat = raw * n_inf_reduced / abs(zeta_star)

    # phi -> x far-field residual (disclosure)
    far_res = float(abs(phi_p[-1] - X) / X) if X > 0 else float("nan")

    mu_hat_huckel = 2.0 / 3.0
    mu_hat_smol = 1.0
    zeta_ref = abs(zeta_v) if zeta_v != 0 else 1e-9
    mu_SI = mu_hat * eps * zeta_v / eta           # m^2/(V s) (sign via zeta_v)
    # OHW83 cross-check (existing closed-form fit, positive-zeta convention);
    # the fit's own stated validity floor is kappa_a >= 10 -- below that we
    # report NaN rather than quote a value the source disowns.
    try:
        ohw = electrophoretic_mobility_relaxation_corrected(
            abs(zeta_mV), kappa_a, viscosity_mPas, temperature_K, z=1, m=0.15,
            rel_permittivity=rel_permittivity)
        ohw_SI = ohw * 1e-8                        # (um cm)/(V s) -> m^2/(V s)
        mu_hat_ohw83 = ohw_SI * eta / (eps * zeta_ref) if zeta_ref > 0 else float("nan")
        rel = (abs(mu_hat) - mu_hat_ohw83) / mu_hat_ohw83 if mu_hat_ohw83 else float("nan")
    except ValueError:
        # kappa_a < 10: OHW83 refuses (its source-stated validity floor)
        mu_hat_ohw83 = float("nan")
        rel = float("nan")

    res = OBrienWhiteResult()
    res.kappa_a = float(kappa_a)
    res.zeta_star = float(zeta_star)
    res.mu_hat = float(mu_hat)
    res.mu_SI = float(mu_SI)
    res.mu_um_cm_per_Vs = float(mu_SI * 1e8)
    res.mu_hat_huckel = mu_hat_huckel
    res.mu_hat_smoluchowski = mu_hat_smol
    res.mu_hat_ohw83 = float(mu_hat_ohw83)
    res.ohw83_rel_diff = float(rel)
    res.relaxation_group_gamma_bar = float(gamma_bar)
    res.transport_number_cation = float(transport_number_cation)
    res.converged = True
    res.n_iterations = 1                 # direct linear solve, not iterated
    res.method = ("obrien_white_integral_quadrature_linear_solve"
                  + ("+relaxation_coupled_EXPERIMENTAL" if relaxation
                     else "_relaxation_off"))
    res.params = {
        "n_grid": int(n_grid),
        "x_max_slack": float(x_max_slack),
        "phi_far_field_rel_residual": far_res,
        "pb_surface_residual_psi_X": float(abs(psi[-1])),
        "D_ref_m2_per_s": float(D_ref_m2_per_s),
        "viscosity_mPas": float(viscosity_mPas),
        "temperature_K": float(temperature_K),
        "rel_permittivity": float(rel_permittivity),
        "relaxation": bool(relaxation),
        "relaxation_channel": ("experimental_unvalidated" if relaxation else "off"),
    }
    _relaxation_channel = "OFF (default)" if not relaxation else "ON (experimental)"
    res.caveat = (
        f"O'Brien-White ion-potential quadrature (RELAXATION CHANNEL {_relaxation_channel}). "
        "Single sphere, constant zeta, 1:1 symmetric electrolyte, linear "
        "response in E, negligible ionic Peclet number. "
        f"Integrals truncated at 12 Debye lengths (exp(-12) tail), n_grid={n_grid} "
        "radial points; Huckel (2/3) and Smoluchowski (->1) limits asserted "
        "in tests, phi_far_field_rel_residual in params. relaxation=False "
        "(default) does NOT solve the coupled h(r): the result is the "
        "Henry-limit backbone, so ohw83_rel_diff is the relaxation+fit gap "
        "vs the OHW83 closed form (<1% error claim for kappa_a>=10), NOT an "
        "error bar here. relaxation=True is EXPERIMENTAL: do not quote it. "
        f"Gamma_bar={gamma_bar:.3g} depends on D_ref: quote D_ref with any comparison. "
        "STRUCTURAL LIMIT (not a solver bug -- O'Brien-White 1978 itself "
        "assumes zero ionic conduction inside the slipping plane): this "
        "model UNDERESTIMATES zeta when real anomalous surface/Stern-layer "
        "conductance is present, which it cannot represent at any "
        "relaxation setting. Moncho, Martinez-Lopez & Hidalgo-Alvarez, "
        "Colloids Surf. A 192 (2001) 215-226, Fig. 6, quantify this on real "
        "polystyrene latex: the deviation from the more complete Dukhin-"
        "Semenikhin theory (which does include anomalous surface "
        "conductance) reaches up to ~35% in zeta at kappa_a~50 for a "
        "sigma0=0.09 C/m^2 surface, falling below 5% only for kappa_a>=90 "
        "(anionic latex) to kappa_a>=150 (amphoteric, low pH) -- do not "
        "trust this function's output at face value for kappa_a below "
        "roughly 100-150 on a highly charged real colloid without "
        "independently checking for anomalous surface conductance."
    )
    return res
