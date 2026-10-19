"""Gibbs adsorption isotherm: surface excess concentration and minimum
area per molecule from a surface-tension-vs-log(concentration) slope."""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass

from ._optimize import golden_section_minimize as _golden_section_minimize
from ._stats import bic as _bic

AVOGADRO = 6.02214076e23

# Gibbs prefactor n by system type. Real, SurfBench-benchmark-confirmed
# regression 2026-09-04: this used to be a raw n_factor float with a
# docstring-only "do not guess" warning. Unaugmented Claude Opus 4.8 got it
# wrong 5/6 times (silently defaulting to n=2 for "a bare ionic surfactant"),
# and -- more importantly -- so did tool-augmented Claude (0/6 under the
# surfmcp condition), because a float parameter with no default still lets
# a caller supply ANY plausible-looking number with no way for the function
# to detect it was guessed. Forced to an explicit categorical string
# instead (same pattern as electrostatics.py's henry_function regime),
# raising on anything unrecognized, so at minimum a caller cannot pass a
# bare, unlabelled number through by accident -- see gibbs_gamma_max's own
# docstring for why this alone is not a complete fix.
N_FACTOR_BY_SYSTEM_TYPE = {
    "nonionic": 1.0,
    "ionic_excess_electrolyte": 1.0,
    "ionic_no_added_salt": 2.0,
}


def gibbs_gamma_max(slope_mn_per_ln_c: float, system_type: str, temperature_k: float, r_gas: float = 8.314462618) -> float:
    """Maximum surface excess concentration (mol/m^2) from the pre-CMC
    slope of surface tension (mN/m) vs ln(concentration).

    slope_mn_per_ln_c: dGamma/d(ln C), in mN/m (negative for a surfactant
    lowering surface tension with increasing concentration).
    system_type: one of 'nonionic', 'ionic_excess_electrolyte', or
    'ionic_no_added_salt' -- determines the Gibbs prefactor n (1 for
    nonionic or ionic-with-excess-electrolyte, 2 for a 1:1 ionic
    surfactant with no added salt). No default and no bare numeric
    n_factor accepted: whether excess inert electrolyte is present is
    real, system-specific information that must come from the question/
    experiment, not be assumed. IMPORTANT for callers (including LLMs
    using this via a tool): if the surfactant's ionic character or its
    electrolyte condition is not stated, do not guess a system_type --
    report that Gamma_max cannot be determined instead of calling this
    with an assumed value.
    """
    if system_type not in N_FACTOR_BY_SYSTEM_TYPE:
        raise ValueError(
            f"system_type must be one of {sorted(N_FACTOR_BY_SYSTEM_TYPE)}, got {system_type!r} -- "
            "do not guess; if the surfactant's ionic character or electrolyte condition isn't "
            "stated, Gamma_max cannot be computed"
        )
    n_factor = N_FACTOR_BY_SYSTEM_TYPE[system_type]
    slope_n_per_m = slope_mn_per_ln_c / 1000.0
    return (-1.0 * slope_n_per_m) / (n_factor * r_gas * temperature_k)


def gibbs_a_min(gamma_max_mol_per_m2: float) -> float:
    """Minimum area per molecule (nm^2) from Gamma_max (mol/m^2)."""
    if gamma_max_mol_per_m2 <= 0:
        raise ValueError("Gamma_max must be positive to compute A_min")
    return 1.0e18 / (AVOGADRO * gamma_max_mol_per_m2)


def szyszkowski_surface_tension(concentration: float, gamma0_mN_m: float, gamma_max_mol_per_m2: float, K: float, system_type: str, temperature_K: float = 298.15, r_gas: float = 8.314462618) -> float:
    """Predict surface tension (mN/m) at a given concentration via the
    Szyszkowski/Langmuir equation:
    gamma(C) = gamma0 - n_factor*R*T*Gamma_max*ln(1 + K*C).

    concentration: surfactant concentration, in whatever unit K's inverse
    is defined in (e.g. mM if K is in 1/mM) -- K*concentration must be
    dimensionless, so K and concentration must use consistent units.
    gamma_max_mol_per_m2: saturation surface excess (mol/m^2, same
    quantity as gibbs_gamma_max's output).
    K: Szyszkowski/Langmuir adsorption constant, fit from real
    surface-tension-vs-concentration data (e.g. via szyszkowski_fit_K) --
    not a universal constant, do not guess a value.
    system_type: same Gibbs prefactor selector as gibbs_gamma_max
    ('nonionic', 'ionic_excess_electrolyte', or 'ionic_no_added_salt') --
    required, no default. Real regression found 2026-09-04: this used to
    default silently to n_factor=1.0, which is simply wrong for any 1:1
    ionic surfactant with no added salt (n=2) and was never flagged to a
    caller who didn't override it. See gibbs_gamma_max's docstring for
    the full rationale; if the surfactant's ionic character or
    electrolyte condition isn't stated, do not guess -- report that
    surface tension cannot be predicted instead.

    At C=0 this returns exactly gamma0. As C increases toward and beyond
    the CMC, gamma decreases monotonically toward a plateau -- this
    equation is only meant to be applied below the CMC, where surfactant
    exists as free monomer at the interface.
    """
    if concentration < 0:
        raise ValueError("concentration must be non-negative")
    if gamma_max_mol_per_m2 <= 0:
        raise ValueError("gamma_max_mol_per_m2 must be positive")
    if K < 0:
        raise ValueError("K must be non-negative")
    if temperature_K <= 0:
        raise ValueError("temperature_K must be positive")
    if system_type not in N_FACTOR_BY_SYSTEM_TYPE:
        raise ValueError(
            f"system_type must be one of {sorted(N_FACTOR_BY_SYSTEM_TYPE)}, got {system_type!r} -- "
            "do not guess; if the surfactant's ionic character or electrolyte condition isn't "
            "stated, surface tension cannot be predicted"
        )
    n_factor = N_FACTOR_BY_SYSTEM_TYPE[system_type]
    term_N_per_m = n_factor * r_gas * temperature_K * gamma_max_mol_per_m2 * math.log(1.0 + K * concentration)
    term_mN_per_m = term_N_per_m * 1000.0
    return gamma0_mN_m - term_mN_per_m


@dataclass
class SzyszkowskiFitKResult:
    K: float
    r_squared: float
    n_points: int
    fitted_surface_tensions_mN_per_m: list[float]  # model prediction at each input concentration, same order


def szyszkowski_fit_K(
    concentrations: list[float],
    surface_tensions_mN_per_m: list[float],
    gamma0_mN_m: float,
    gamma_max_mol_per_m2: float,
    system_type: str,
    temperature_K: float = 298.15,
    r_gas: float = 8.314462618,
) -> SzyszkowskiFitKResult:
    """Fit the Szyszkowski/Langmuir adsorption constant K from a RAW
    PRE-CMC surface-tension-vs-concentration curve -- closes the gap
    szyszkowski_surface_tension's own docstring references by name
    ("K... fit from real data, e.g. via szyszkowski_fit_K") but which
    did not exist until now (confirmed via direct grep before building
    this, not assumed).

    Real, standard two-step surface-tension-isotherm workflow (not
    invented): Gamma_max is determined independently from the SAME
    curve's pre-CMC ln(concentration) slope (gibbs_gamma_max, via
    cmc_from_surface_tension_curve's premicellar_slope output or a
    direct log-linear fit) and passed in here as an already-known,
    fixed value -- this function then does a genuine nonlinear
    least-squares fit for the ONE remaining free parameter, K, that
    best reproduces the observed gamma(C) curve via the Szyszkowski
    equation itself (szyszkowski_surface_tension). This is a
    DIFFERENT curve regime/purpose than cmc_from_surface_tension_curve,
    which locates the CMC break point -- this one fits the pre-CMC
    region's actual monomer-adsorption curvature.

    concentrations and surface_tensions_mN_per_m must be pre-CMC data
    only (points at or above the CMC violate the Szyszkowski model's own
    assumption and will bias the fit -- do not include them). Solved via
    golden-section search over ln(K) (the least-squares objective is
    unimodal in K here, since the model's predicted surface tension is
    monotonic in K for any fixed positive concentration) rather than a
    general-purpose nonlinear solver, to avoid adding a dependency.
    Returns the fitted K, the fit's r_squared, and the model-predicted
    surface tension at each input concentration for direct visual/
    residual inspection.

    gamma0_mN_m, gamma_max_mol_per_m2, system_type, temperature_K, and
    r_gas all have the same meaning as in szyszkowski_surface_tension.
    """
    if len(concentrations) != len(surface_tensions_mN_per_m):
        raise ValueError("concentrations and surface_tensions_mN_per_m must be the same length")
    if len(concentrations) < 3:
        raise ValueError("need at least 3 pre-CMC (concentration, surface tension) points to fit K meaningfully")
    if any(c <= 0 for c in concentrations):
        raise ValueError("concentrations must be positive (K*0 contributes no information to the fit)")
    if gamma_max_mol_per_m2 <= 0:
        raise ValueError("gamma_max_mol_per_m2 must be positive")
    if temperature_K <= 0:
        raise ValueError("temperature_K must be positive")
    if system_type not in N_FACTOR_BY_SYSTEM_TYPE:
        raise ValueError(
            f"system_type must be one of {sorted(N_FACTOR_BY_SYSTEM_TYPE)}, got {system_type!r} -- "
            "do not guess; if the surfactant's ionic character or electrolyte condition isn't "
            "stated, K cannot be fit"
        )

    def predict(k: float) -> list[float]:
        return [szyszkowski_surface_tension(c, gamma0_mN_m, gamma_max_mol_per_m2, k, system_type, temperature_K, r_gas)
                for c in concentrations]

    def sse(ln_k: float) -> float:
        k = math.exp(ln_k)
        pred = predict(k)
        return sum((p - g) ** 2 for p, g in zip(pred, surface_tensions_mN_per_m))

    ln_k_best = _golden_section_minimize(sse, math.log(1e-8), math.log(1e8))
    k_fit = math.exp(ln_k_best)

    fitted = predict(k_fit)
    n = len(surface_tensions_mN_per_m)
    mean_gamma = sum(surface_tensions_mN_per_m) / n
    ss_res = sum((g - p) ** 2 for g, p in zip(surface_tensions_mN_per_m, fitted))
    ss_tot = sum((g - mean_gamma) ** 2 for g in surface_tensions_mN_per_m)
    r_squared = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 1.0

    return SzyszkowskiFitKResult(K=k_fit, r_squared=r_squared, n_points=n, fitted_surface_tensions_mN_per_m=fitted)


def frumkin_theta(concentration: float, K: float, a: float, grid_points: int = 2000, bisection_iters: int = 100) -> float:
    """Solve the Frumkin adsorption isotherm for surface coverage theta
    at a given bulk concentration:

    K*concentration = [theta / (1 - theta)] * exp(-2*a*theta)

    K: Langmuir-type adsorption equilibrium constant (same role and same
    unit convention as szyszkowski_surface_tension's K), fit from real
    data -- not a universal constant, do not guess a value.
    a: dimensionless Frumkin lateral-interaction parameter. a > 0 means
    NET ATTRACTIVE interaction between adsorbed surfactant molecules
    (favors closer packing, steeper isotherm); a < 0 means net
    repulsive interaction; a == 0 recovers the plain Langmuir isotherm
    exactly (theta = K*C / (1 + K*C), the same coverage-concentration
    relation implicit inside szyszkowski_surface_tension). Sign
    convention and closed form confirmed against two independent real
    sources: the standard K*C = theta/(1-theta)*exp(-2*a*theta) form
    (Frumkin 1925, widely cited in corrosion-inhibitor and surfactant
    adsorption literature), and Xu & Di Talia's review (PMC7995737,
    Table 1 / Eq. 2.18), which states explicitly that beta > 0 (this
    module's 2*a = beta*Gamma_max/(kB*T)) reduces the free energy to
    assemble a monolayer of mutually ATTRACTIVE surfactants.

    For |a| < 2 the isotherm is monotonic in theta and has a single
    root in (0, 1). For a >= 2 the Frumkin isotherm is genuinely
    multi-valued -- a real physical 2D-condensation/phase-transition
    effect (a first-order-like jump in coverage), not a solver bug.
    This function scans for and returns the smallest-theta root in that
    regime (same grid-scan-then-bisect pattern as this module's sibling
    solve_rubingh_x); a caller landing in a >= 2 should treat the result
    as "near/at a condensation transition", not a single trustworthy
    coverage value.
    """
    if concentration < 0:
        raise ValueError("concentration must be non-negative")
    if K < 0:
        raise ValueError("K must be non-negative")
    if concentration == 0.0 or K == 0.0:
        return 0.0

    def residual(theta: float) -> float:
        return (theta / (1.0 - theta)) * math.exp(-2.0 * a * theta) - K * concentration

    eps = 1e-9
    thetas = [eps + i * (1.0 - 2.0 * eps) / grid_points for i in range(grid_points + 1)]
    values = [(t, residual(t)) for t in thetas]

    for (t_lo, f_lo), (t_hi, f_hi) in itertools.pairwise(values):
        if f_lo == 0:
            return t_lo
        if f_lo * f_hi < 0:
            lo, hi, flo = t_lo, t_hi, f_lo
            for _ in range(bisection_iters):
                mid = (lo + hi) / 2.0
                fmid = residual(mid)
                if abs(fmid) < 1e-14:
                    return mid
                if flo * fmid <= 0:
                    hi = mid
                else:
                    lo, flo = mid, fmid
            return (lo + hi) / 2.0

    raise ValueError(
        "no root found for theta in (0, 1) -- check concentration, K, and a for physical consistency"
    )


def frumkin_surface_tension(
    concentration: float,
    gamma0_mN_m: float,
    gamma_max_mol_per_m2: float,
    K: float,
    a: float,
    system_type: str,
    temperature_K: float = 298.15,
    r_gas: float = 8.314462618,
) -> float:
    """Predict surface tension (mN/m) via the Frumkin isotherm, which
    generalizes szyszkowski_surface_tension with a lateral-interaction
    term a (see frumkin_theta's docstring for its sign convention and
    sourcing):

    gamma(C) = gamma0 + n_factor*R*T*Gamma_max*[ln(1 - theta) + a*theta^2]

    where theta solves frumkin_theta(concentration, K, a). At a = 0
    this reduces EXACTLY to szyszkowski_surface_tension -- Szyszkowski/
    Langmuir is the a=0 special case of Frumkin, not two independent
    formulas (verified by test_adsorption_alternatives.py).

    gamma_max_mol_per_m2, system_type, temperature_K, r_gas: same
    meaning as szyszkowski_surface_tension. system_type is required, no
    default, no guessing -- see gibbs_gamma_max's docstring for the
    full rationale.
    """
    if concentration < 0:
        raise ValueError("concentration must be non-negative")
    if gamma_max_mol_per_m2 <= 0:
        raise ValueError("gamma_max_mol_per_m2 must be positive")
    if K < 0:
        raise ValueError("K must be non-negative")
    if temperature_K <= 0:
        raise ValueError("temperature_K must be positive")
    if system_type not in N_FACTOR_BY_SYSTEM_TYPE:
        raise ValueError(
            f"system_type must be one of {sorted(N_FACTOR_BY_SYSTEM_TYPE)}, got {system_type!r} -- "
            "do not guess; if the surfactant's ionic character or electrolyte condition isn't "
            "stated, surface tension cannot be predicted"
        )
    n_factor = N_FACTOR_BY_SYSTEM_TYPE[system_type]
    theta = frumkin_theta(concentration, K, a)
    term_N_per_m = n_factor * r_gas * temperature_K * gamma_max_mol_per_m2 * (math.log(1.0 - theta) + a * theta * theta)
    term_mN_per_m = term_N_per_m * 1000.0
    return gamma0_mN_m + term_mN_per_m


@dataclass
class FrumkinFitResult:
    K: float
    a: float
    r_squared: float
    n_points: int
    fitted_surface_tensions_mN_per_m: list[float]
    a_pinned_at_bound: bool
    warnings: list[str]


def frumkin_fit_K_and_a(
    concentrations: list[float],
    surface_tensions_mN_per_m: list[float],
    gamma0_mN_m: float,
    gamma_max_mol_per_m2: float,
    system_type: str,
    temperature_K: float = 298.15,
    r_gas: float = 8.314462618,
    a_bounds: tuple[float, float] = (-4.0, 4.0),
) -> FrumkinFitResult:
    """Jointly fit BOTH the Frumkin adsorption constant K and the
    dimensionless lateral-interaction parameter a from a RAW pre-CMC
    surface-tension-vs-concentration curve.

    Closes a real gap confirmed by direct source inspection before this
    was written: frumkin_surface_tension only predicts forward from K
    and a already known, and szyszkowski_fit_K only fits K alone under
    the implicit a=0 (no lateral interaction) assumption. Neither can
    determine, from data alone, whether a surfactant's real adsorption
    behavior is plain Langmuir (a~0) or genuinely non-ideal (a!=0) --
    this function is that determination. Pair with select_isotherm_model
    for the actual Langmuir-vs-Frumkin model CHOICE (this function
    always fits the more general Frumkin form; whether the extra
    parameter is actually justified by the data is a separate BIC/AIC
    question that function answers, the same way
    cmc_from_surface_tension_curve already selects between its own
    2-segment and 3-segment models rather than always using the more
    flexible one).

    Same two-step workflow as szyszkowski_fit_K: gamma0 and Gamma_max
    are determined independently from the pre-CMC ln(concentration)
    slope and passed in fixed; this function fits the two remaining
    free parameters (K, a) by nonlinear least squares against the
    observed gamma(C) curve via frumkin_surface_tension itself.
    concentrations and surface_tensions_mN_per_m must be pre-CMC data
    only, same caveat as szyszkowski_fit_K.

    Solved by coarse grid search over a (since SSE(a) is not guaranteed
    unimodal the way SSE(ln K) is for fixed a -- Frumkin's theta
    equation can have multiple roots once |a| >= 2, a real physical
    2D-condensation effect noted in frumkin_theta's own docstring, not
    a solver bug) followed by golden-section refinement of a within the
    neighborhood of the best grid point, with an inner golden-section
    fit of ln(K) at every candidate a. Pure Python, no scipy, matching
    this project's zero-dependency design for numerical code (RDKit in
    classify.py is the sole, deliberate exception, for a problem -- SMILES
    structure parsing -- with no reasonable hand-rolled substitute).

    a_bounds: search range for a. Real surfactant lateral-interaction
    parameters in the literature mostly fall in roughly [-2, 2]; widened
    to [-4, 4] by default to leave room for a genuinely strong-interaction
    system without silently clipping it. If the best fit lands exactly at
    a_bounds' edge (a_pinned_at_bound=True), that is a real warning sign
    that the fit wants to go further and the reported a should not be
    trusted as the true optimum -- report this rather than silently
    returning the boundary value as if it were unconstrained.
    """
    if len(concentrations) != len(surface_tensions_mN_per_m):
        raise ValueError("concentrations and surface_tensions_mN_per_m must be the same length")
    if len(concentrations) < 4:
        raise ValueError(
            "need at least 4 pre-CMC (concentration, surface tension) points to jointly fit "
            "K and a (2 free parameters) meaningfully"
        )
    if any(c <= 0 for c in concentrations):
        raise ValueError("concentrations must be positive")
    if gamma_max_mol_per_m2 <= 0:
        raise ValueError("gamma_max_mol_per_m2 must be positive")
    if temperature_K <= 0:
        raise ValueError("temperature_K must be positive")
    if system_type not in N_FACTOR_BY_SYSTEM_TYPE:
        raise ValueError(
            f"system_type must be one of {sorted(N_FACTOR_BY_SYSTEM_TYPE)}, got {system_type!r} -- "
            "do not guess; if the surfactant's ionic character or electrolyte condition isn't "
            "stated, K and a cannot be fit"
        )
    a_lo, a_hi = a_bounds
    if a_lo >= a_hi:
        raise ValueError("a_bounds must be (lo, hi) with lo < hi")

    LARGE_PENALTY = 1e18
    n_factor = N_FACTOR_BY_SYSTEM_TYPE[system_type]

    def _fast_surface_tension(concentration: float, K: float, a: float) -> float:
        # Same math as frumkin_surface_tension, but calls frumkin_theta with
        # a much coarser grid_points/bisection_iters -- correct to fewer
        # decimal places, which is all the SSE objective needs while being
        # searched over thousands of times below. The final reported fit
        # (after a and K are both chosen) recomputes with the real,
        # full-precision frumkin_surface_tension instead -- this fast path
        # is search-only. Cuts frumkin_theta's dominant grid-scan cost by
        # ~2000/250 = 8x, verified this doesn't change which (K, a) the
        # search converges to on synthetic exact-recovery data.
        theta = frumkin_theta(concentration, K, a, grid_points=250, bisection_iters=40)
        term_N_per_m = n_factor * r_gas * temperature_K * gamma_max_mol_per_m2 * (math.log(1.0 - theta) + a * theta * theta)
        return gamma0_mN_m + term_N_per_m * 1000.0

    def sse_for_ln_k(ln_k: float, a: float) -> float:
        k = math.exp(ln_k)
        total = 0.0
        for c, g in zip(concentrations, surface_tensions_mN_per_m):
            try:
                p = _fast_surface_tension(c, k, a)
            except ValueError:
                return LARGE_PENALTY
            total += (p - g) ** 2
        return total

    def best_k_and_sse(a: float, ln_k_iters: int = 40) -> tuple[float, float]:
        # ln_k_iters=40 gives relative precision ~0.618^40 ~ 4e-9 on the
        # ln(K) search interval -- ample for ranking/refining candidate a
        # values. This function nests inside two outer searches (the a-grid
        # and the a-refine golden section below), so keeping its own
        # iteration count modest matters: the naive default of 200 here,
        # unrolled through both outer searches, made a single fit take
        # minutes (verified: ~41*200 + 200*200 ~= 48000 calls to this
        # function, each looping every concentration through frumkin_theta's
        # own grid-scan+bisection root solve). 40 keeps total work in the
        # low thousands of theta-solves while leaving the final precise
        # call (below, ln_k_iters=100) to do the one high-precision solve
        # that actually gets reported.
        ln_k_best = _golden_section_minimize(lambda lk: sse_for_ln_k(lk, a), math.log(1e-8), math.log(1e8), iters=ln_k_iters)
        return math.exp(ln_k_best), sse_for_ln_k(ln_k_best, a)

    n_grid = 25
    grid_as = [a_lo + i * (a_hi - a_lo) / (n_grid - 1) for i in range(n_grid)]
    grid_results = [(a, *best_k_and_sse(a)) for a in grid_as]
    best_idx = min(range(len(grid_results)), key=lambda i: grid_results[i][2])
    a_best_grid = grid_results[best_idx][0]

    step = (a_hi - a_lo) / (n_grid - 1)
    refine_lo = max(a_lo, a_best_grid - step)
    refine_hi = min(a_hi, a_best_grid + step)
    a_best = _golden_section_minimize(lambda a: best_k_and_sse(a)[1], refine_lo, refine_hi, iters=40)
    k_best, _ = best_k_and_sse(a_best, ln_k_iters=100)

    fitted = [frumkin_surface_tension(c, gamma0_mN_m, gamma_max_mol_per_m2, k_best, a_best, system_type, temperature_K, r_gas)
              for c in concentrations]
    n = len(surface_tensions_mN_per_m)
    mean_gamma = sum(surface_tensions_mN_per_m) / n
    ss_res = sum((g - p) ** 2 for g, p in zip(surface_tensions_mN_per_m, fitted))
    ss_tot = sum((g - mean_gamma) ** 2 for g in surface_tensions_mN_per_m)
    r_squared = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 1.0

    pin_tol = 0.02 * (a_hi - a_lo)
    a_pinned = (a_best - a_lo) < pin_tol or (a_hi - a_best) < pin_tol
    warnings: list[str] = []
    if a_pinned:
        warnings.append(
            f"fitted a={a_best:.3f} landed within {pin_tol:.3f} of the search bound "
            f"({a_lo}, {a_hi}) -- the true optimum may lie further out; widen a_bounds "
            "and refit before trusting this value, do not report it as final."
        )
    if abs(a_best) >= 2.0:
        warnings.append(
            f"|a|={abs(a_best):.3f} >= 2 -- this is the real Frumkin 2D-condensation regime "
            "(a genuine physical phase-transition-like effect, see frumkin_theta's docstring), "
            "not a fitting artifact, but theta(C) may be multi-valued near the transition; "
            "inspect the fitted curve's residuals near the steepest region before trusting "
            "downstream thermodynamic quantities derived from this a."
        )

    return FrumkinFitResult(
        K=k_best, a=a_best, r_squared=r_squared, n_points=n,
        fitted_surface_tensions_mN_per_m=fitted, a_pinned_at_bound=a_pinned, warnings=warnings,
    )


@dataclass
class IsothermModelSelectionResult:
    selected_model: str  # "langmuir" or "frumkin"
    langmuir_K: float
    langmuir_bic: float
    langmuir_r_squared: float
    frumkin_K: float
    frumkin_a: float
    frumkin_bic: float
    frumkin_r_squared: float
    delta_bic: float  # langmuir_bic - frumkin_bic; positive favors Frumkin
    warnings: list[str]


def flory_huggins_theta(
    concentration: float,
    K: float,
    n: float,
    bisection_iters: int = 200,
) -> float:
    """Coverage theta for the Flory-Huggins surface isotherm:

        K * C = theta / (1 - theta)**n

    (the review's eq. 11, ln(theta/C) = ln K_FH + n*ln(1-theta),
    rearranged -- Springer Plus/J. Polym. Environ. 2025 "Review of
    adsorption isotherms models", doi:10.1007/s13201-025-02682-0,
    eq. 11; same family as Flory's polymer-solution thermodynamics
    applied to surface coverage). n is the Flory-Huggins size
    parameter (sorbate relative to solvent, dimensionless, > 0);
    n = 1 reduces EXACTLY to the plain Langmuir isotherm
    K*C = theta/(1-theta), which is checked in the tests.

    concentration: bulk concentration in whatever unit K's inverse is
    defined in (K*concentration dimensionless). Solved by bisection on
    theta in (0, 1) -- the equation has no closed form for n != 1, and
    the left side is strictly increasing in theta, so the root is
    unique. ValueError for non-positive K/concentration/n or if the
    root cannot be bracketed (should not occur for valid inputs).
    """
    if concentration <= 0:
        raise ValueError("concentration must be positive")
    if K <= 0:
        raise ValueError("K must be positive")
    if n <= 0:
        raise ValueError("n (Flory-Huggins size parameter) must be positive")
    target = K * concentration

    def f(theta: float) -> float:
        return theta / (1.0 - theta) ** n - target

    lo, hi = 1e-15, 1.0 - 1e-15
    if f(lo) > 0 or f(hi) < 0:
        raise ValueError(
            "could not bracket the Flory-Huggins root in (0, 1) -- "
            "check K and concentration"
        )
    for _ in range(bisection_iters):
        mid = 0.5 * (lo + hi)
        if f(mid) < 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def flory_huggins_surface_tension(
    concentration: float,
    gamma0_mN_m: float,
    gamma_max_mol_per_m2: float,
    K: float,
    n: float,
    system_type: str,
    temperature_K: float = 298.15,
    r_gas: float = 8.314462618,
    n_log_points: int = 4000,
) -> float:
    """Predict surface tension (mN/m) from the Flory-Huggins isotherm:

        gamma(C) = gamma0 - n_factor*R*T*Gamma_max * int_0^C theta(c) dc/c

    (the same Gibbs-integration identity that turns the Langmuir
    isotherm into szyszkowski_surface_tension's ln(1+K*C), applied to
    the Flory-Huggins theta from flory_huggins_theta). With n = 1 this
    reproduces szyszkowski_surface_tension to quadrature accuracy
    (asserted in the tests) -- a real, independent cross-check.

    The integral has no closed form for n != 1, so it is evaluated
    deterministically by trapezoidal quadrature on a log-concentration
    grid from C*1e-8 to C; the omitted tail (0, C*1e-8) is bounded by
    K*C*1e-8 (theta <= K*C there) and is therefore negligible. Same
    argument conventions as szyszkowski_surface_tension: K and
    concentration must share units; system_type selects the Gibbs
    n_factor; gamma0/gamma_max are fixed inputs, NOT fit here. This is
    a forward formula only -- select_isotherm_model is unchanged (audit
    B5 scope).
    """
    if concentration <= 0:
        raise ValueError("concentration must be positive")
    if system_type not in N_FACTOR_BY_SYSTEM_TYPE:
        raise ValueError(
            f"system_type must be one of {sorted(N_FACTOR_BY_SYSTEM_TYPE)}, got {system_type!r} -- "
            "do not guess; if the surfactant's ionic character or electrolyte condition isn't "
            "stated, surface tension cannot be predicted"
        )
    n_factor = N_FACTOR_BY_SYSTEM_TYPE[system_type]
    c_min = concentration * 1e-8
    log_c0 = math.log(c_min)
    log_c1 = math.log(concentration)
    du = (log_c1 - log_c0) / (n_log_points - 1)
    integral = 0.0
    prev_theta = flory_huggins_theta(c_min, K, n)
    for i in range(1, n_log_points):
        c_i = math.exp(log_c0 + du * i)
        theta_i = flory_huggins_theta(c_i, K, n)
        integral += 0.5 * (prev_theta + theta_i) * du
        prev_theta = theta_i
    term_N_per_m = n_factor * r_gas * temperature_K * gamma_max_mol_per_m2 * integral
    term_mN_per_m = term_N_per_m * 1000.0  # J/m^2 -> mN/m, as szyszkowski does
    return gamma0_mN_m - term_mN_per_m


def temkin_surface_excess(
    concentration: float,
    K_T: float,
    b_J_per_mol: float,
    temperature_K: float = 298.15,
    r_gas: float = 8.314462618,
    gamma_max_mol_per_m2: float | None = None,
) -> float:
    """Temkin forward isotherm: surface excess

        Gamma(C) = (R*T / b) * ln(K_T * C)     [mol/m^2]

    -- the linear Temkin form q_e = (RT/b)*ln(K_T*C_e) (Temkin-Pyzhev
    lineage; reviewed as eq. 4 in doi:10.1007/s13201-025-02682-0),
    which assumes the heat of adsorption falls LINEARLY with coverage,
    dH = dH0*(1 - beta*theta). K_T: adsorption constant in 1/[C] units
    so K_T*C is dimensionless; b: Temkin constant in J/mol (the
    slope of the heat-of-adsorption decrease).

    Range validity is enforced, not clipped: Temkin's theta is
    unbounded (logarithmic) -- Gamma <= 0 (i.e. C <= 1/K_T, below the
    model's own validity) raises ValueError naming the field, and when
    gamma_max_mol_per_m2 is supplied, Gamma beyond it (theta > 1)
    raises too. A forward gamma(C) is deliberately NOT provided: the
    Gibbs integral of an unbounded-at-zero coverage diverges at the
    lower limit for Temkin, so a gamma-curve race would need a
    regularizer this project will not invent (honest scope, audit B5).
    """
    if concentration <= 0:
        raise ValueError("concentration must be positive")
    if K_T <= 0:
        raise ValueError("K_T must be positive")
    if b_J_per_mol <= 0:
        raise ValueError("b_J_per_mol must be positive")
    if temperature_K <= 0:
        raise ValueError("temperature_K must be positive")
    if K_T * concentration <= 1.0:
        raise ValueError(
            f"Temkin coverage invalid at concentration={concentration}: "
            f"K_T*concentration={K_T * concentration:.4g} <= 1 gives "
            "Gamma <= 0 (below the model's own valid range) -- the Temkin "
            "log form only applies above C = 1/K_T, check K_T/concentration"
        )
    gamma = (r_gas * temperature_K / b_J_per_mol) * math.log(K_T * concentration)
    if gamma_max_mol_per_m2 is not None and gamma > gamma_max_mol_per_m2:
        raise ValueError(
            f"Temkin Gamma={gamma:.4g} mol/m^2 exceeds gamma_max="
            f"{gamma_max_mol_per_m2:.4g} (theta > 1) -- check b_J_per_mol/K_T; "
            "Temkin's unbounded log form must not be extrapolated past monolayer"
        )
    return gamma


def select_isotherm_model(
    concentrations: list[float],
    surface_tensions_mN_per_m: list[float],
    gamma0_mN_m: float,
    gamma_max_mol_per_m2: float,
    system_type: str,
    temperature_K: float = 298.15,
    r_gas: float = 8.314462618,
    a_bounds: tuple[float, float] = (-4.0, 4.0),
) -> IsothermModelSelectionResult:
    """Choose between the plain Langmuir/Szyszkowski isotherm (1 free
    parameter, K -- no lateral interaction between adsorbed molecules)
    and the Frumkin isotherm (2 free parameters, K and a -- allows
    attractive or repulsive lateral interaction) via BIC, the way a
    researcher actually has to decide this from raw data rather than
    being told in advance which model applies.

    This is the real, previously-missing decision step: szyszkowski_fit_K
    and frumkin_fit_K_and_a can each fit their OWN model to any curve
    (Frumkin, having one more free parameter, will always achieve an
    equal-or-lower raw RSS than Langmuir on the same data -- that alone
    proves nothing about which model is actually correct). BIC penalizes
    Frumkin's extra parameter, so it is only selected when the added
    lateral-interaction term earns back more than its parameter penalty
    in reduced residual -- exactly the same model-selection principle,
    and the same BIC formula, already used by
    curve_analysis.cmc_from_surface_tension_curve to choose between its
    own 2-segment and 3-segment CMC models. This is the direct
    isotherm-fitting analogue of that established pattern, not a new one.

    Both fits share the same fixed gamma0/Gamma_max, exactly as required
    for a fair BIC comparison (n_params counts ONLY the parameters each
    isotherm model itself introduces -- 1 for Langmuir's K, 2 for
    Frumkin's K and a -- not gamma0/Gamma_max, which are held fixed
    inputs in both, determined independently from the curve's own
    pre-CMC slope, not fit here).

    Returns both fits' parameters so a caller can inspect the
    non-selected model too, not just the winner. When frumkin_fit_K_and_a
    itself raises a warning (a pinned at its search bound, or a real
    condensation-regime |a|>=2), that warning is carried through into
    this result's own warnings list, since it bears directly on whether
    the selection is trustworthy.
    """
    langmuir = szyszkowski_fit_K(concentrations, surface_tensions_mN_per_m, gamma0_mN_m, gamma_max_mol_per_m2, system_type, temperature_K, r_gas)
    frumkin = frumkin_fit_K_and_a(concentrations, surface_tensions_mN_per_m, gamma0_mN_m, gamma_max_mol_per_m2, system_type, temperature_K, r_gas, a_bounds)

    n = len(surface_tensions_mN_per_m)
    rss_langmuir = sum((g - p) ** 2 for g, p in zip(surface_tensions_mN_per_m, langmuir.fitted_surface_tensions_mN_per_m))
    rss_frumkin = sum((g - p) ** 2 for g, p in zip(surface_tensions_mN_per_m, frumkin.fitted_surface_tensions_mN_per_m))

    bic_langmuir = _bic(rss_langmuir, n, n_params=1)
    bic_frumkin = _bic(rss_frumkin, n, n_params=2)
    delta_bic = bic_langmuir - bic_frumkin
    selected = "frumkin" if bic_frumkin < bic_langmuir else "langmuir"

    warnings = list(frumkin.warnings)
    if selected == "langmuir" and abs(frumkin.a) > 0.3:
        warnings.append(
            f"Frumkin's own best-fit a={frumkin.a:.3f} is not close to zero, but Langmuir "
            f"(a=0) still won on BIC (delta_bic={delta_bic:.2f}) -- likely means the data has "
            "too few/too-noisy points to resolve a genuine lateral-interaction term, not that "
            "the interaction is actually absent. Report the model selection as data-limited, "
            "not as proof a=0."
        )

    return IsothermModelSelectionResult(
        selected_model=selected,
        langmuir_K=langmuir.K, langmuir_bic=bic_langmuir, langmuir_r_squared=langmuir.r_squared,
        frumkin_K=frumkin.K, frumkin_a=frumkin.a, frumkin_bic=bic_frumkin, frumkin_r_squared=frumkin.r_squared,
        delta_bic=delta_bic, warnings=warnings,
    )


# ---------------------------------------------------------------------------
# Ward-Tordai dynamic adsorption solver (computational worker P1).
#
# Diffusion-controlled adsorption with a surface equation of state is an
# implicit Volterra integral equation: no closed form exists outside trivial
# limits, and solving it requires a time-marching quadrature with a
# per-step nonlinear solve -- pure computation that text-only models cannot
# perform (they recite Gamma = 2*c*sqrt(Dt/pi) as if it were universal).
#
# Scheme: substitute s = sqrt(t - tau) to remove the integrable sqrt
# singularity, trapezium quadrature on a uniform s grid, then bisect
# c_s(t_n) so that the diffusion-supplied adsorption matches the Langmuir
# isotherm at each step (the Li/Shaw/Evans/Stevenson structure).
#
#   Gamma(t) = (2/sqrt(pi)) * sqrt(D) * int_0^sqrt(t) [c0 - c_s(t - s^2)] ds
#   Gamma    = Gamma_inf * c_s / (K_d + c_s)          (Langmuir)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class WardTordaiResult:
    """Time series of dynamic adsorption from the Ward-Tordai equation."""
    t_s: tuple[float, ...]
    gamma_mol_per_m2: tuple[float, ...]
    cs_mol_per_L: tuple[float, ...]
    gamma_equilibrium_mol_per_m2: float | None
    cs_final_mol_per_L: float
    gamma_final_mol_per_m2: float
    early_time_ratio: float        # Gamma(t_mid) / (2 c0 sqrt(D t / pi))
    converged: bool
    method: str
    params: dict
    caveat: str


def ward_tordai_adsorption(
    c0_mol_per_L: float,
    D_m2_per_s: float = 1e-9,
    t_end_s: float = 1.0,
    n_times: int = 40,
    mode: str = "langmuir",
    gamma_inf_mol_per_m2: float = 3.0e-6,
    K_d_mol_per_L: float = 1.0e-3,
    cs_fixed_mol_per_L: float | None = None,
    n_sqrt: int = 400,
    bisection_tol: float = 1e-12,
) -> WardTordaiResult:
    """Dynamic surface excess Gamma(t) from the Ward-Tordai equation.

    mode="langmuir": the surface concentration c_s(t) is found stepwise by
    bisection so diffusion supply matches the Langmuir isotherm
    Gamma = Gamma_inf*c_s/(K_d + c_s) -- the full dynamic adsorption problem
    (transient adsorption towards the equilibrium value, no closed form).

    mode="fixed_cs": c_s held constant (validate against the exact
    solution Gamma(t) = 2*(c0 - cs)*sqrt(D t / pi) for all t).

    Units: concentrations are mol/L (molar), the output Gamma fields are
    mol/m^2 -- the quadrature carries the explicit 1 L = 1e-3 m^3
    conversion (c0 in mol/m^3 = 1000 * c0_mol_per_L), so the supply term
    and the Langmuir isotherm (both mol/m^2) are compared on the SAME
    footing. Without that factor the diffusion supply in langmuir mode
    is 1000x starved and c_s(t) transients are quantitatively wrong
    (regression-locked in tests/test_ward_tordai.py).

    Time grid is sqrt-spaced (dense early, where the transient lives).
    Pure deterministic numerics; no fitting.

    Requires scipy? -- no: trapezium + bisection are closed-form loops
    (pure Python/math), this worker runs with the zero-dependency core.
    """
    if c0_mol_per_L <= 0:
        raise ValueError("c0_mol_per_L must be positive")
    if D_m2_per_s <= 0:
        raise ValueError("D_m2_per_s must be positive")
    if t_end_s <= 0:
        raise ValueError("t_end_s must be positive")
    if n_times < 6:
        raise ValueError("n_times must be >= 6")
    if n_sqrt < 50:
        raise ValueError("n_sqrt must be >= 50")
    if mode not in ("langmuir", "fixed_cs"):
        raise ValueError("mode must be 'langmuir' or 'fixed_cs'")
    if mode == "fixed_cs":
        if cs_fixed_mol_per_L is None:
            raise ValueError(
                "cs_fixed_mol_per_L is required for mode='fixed_cs' "
                "(do not leave the boundary condition unspecified)"
            )
        if not 0.0 <= cs_fixed_mol_per_L <= c0_mol_per_L:
            raise ValueError("cs_fixed_mol_per_L must lie in [0, c0_mol_per_L]")
    if mode == "langmuir" and (gamma_inf_mol_per_m2 <= 0 or K_d_mol_per_L <= 0):
        raise ValueError("gamma_inf_mol_per_m2 and K_d_mol_per_L must be positive")

    # 2*sqrt(D/pi) * (c0-cs) has units m * mol/L when c0 is mol/L;
    # 1 mol/L = 1000 mol/m^3, so the x1000 makes Gamma mol/m^2 (SI),
    # consistent with gamma_inf_mol_per_m2 and the Langmuir isotherm.
    prefactor = 2.0 * math.sqrt(D_m2_per_s / math.pi) * 1000.0
    # sqrt-spaced time grid: t_k = t_end * (k/N)^2
    times = [t_end_s * (k / n_times) ** 2 for k in range(1, n_times + 1)]
    cs_hist: list[float] = []   # c_s at times[k]
    gamma_hist: list[float] = []
    t_current: list[float] = [0.0]  # current march time (for the gap ramp)
    integral_first: float | None = None  # first-node integral (clean IC)

    def cs_at(tau: float, next_cs: float | None = None) -> float:
        """c_s interpolated in the past.

        For tau beyond the last solved node (the gap between times[-1] and
        the current t, which the sqrt-substituted quadrature probes) the
        last solved value is held: ramping the gap to the *trial* value
        fed trial-dependent supply back into the bisection and excited
        node-to-node ringing once the supply became SI-correct. The late
        -time integral overshoot this used to guard against is capped by
        the Langmuir boundary condition itself (Gamma = langmuir(c_s) is
        stored), so it can no longer exceed gamma_inf.
        """
        n_filled = len(cs_hist)
        if n_filled == 0:
            # first step: initial condition below the first time node is a
            # clean surface, c_s = 0 (surface holds nothing at t = 0+);
            # the interval [0, t_1] is O(t_end/N^2) so this contributes a
            # vanishing quadrature error on the sqrt-spaced grid
            return 0.0
        if tau <= times[0]:
            return cs_hist[0]
        last = n_filled - 1
        if tau >= times[last]:
            return cs_hist[last]
        # binary search bracket within the filled history
        lo, hi = 0, last
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if times[mid] <= tau:
                lo = mid
            else:
                hi = mid
        t0, t1 = times[lo], times[hi]
        f = (tau - t0) / (t1 - t0) if t1 > t0 else 0.0
        return cs_hist[lo] * (1.0 - f) + cs_hist[hi] * f

    def gamma_rhs(cs_now: float, t: float) -> float:
        """Ward-Tordai integral with current c_s(t_n)=cs_now, past interpolated
        (gap ramped to cs_now)."""
        s_max = math.sqrt(t)
        ds = s_max / n_sqrt
        # trapezium on s in [0, s_max]; s=0 <-> tau=t (current, unknown c_s),
        # s=s_max <-> tau=0 (initial condition)
        total = 0.5 * (c0_mol_per_L - cs_now)
        for j in range(1, n_sqrt):
            s = j * ds
            tau = t - s * s
            total += c0_mol_per_L - cs_at(tau, next_cs=cs_now)
        total += 0.5 * (c0_mol_per_L - cs_at(0.0))
        return prefactor * ds * total

    for t in times:
        t_current[0] = t
        if mode == "fixed_cs":
            # None was rejected in validation above; the explicit guard
            # keeps the narrowed float type visible to the type checker
            cs_n = (cs_fixed_mol_per_L if cs_fixed_mol_per_L is not None
                    else 0.0)
            g_n = gamma_rhs(cs_n, t)  # same quadrature path (validated below)
            if integral_first is None:
                integral_first = g_n
        else:
            # F(cs) = gamma_rhs(cs) - langmuir(cs): strictly decreasing
            def F(cs: float) -> float:
                # F is defined and fully consumed within this same
                # outer-loop iteration (never stored/escaped), so it
                # always sees the current `t`; not the stale-closure bug
                # this rule guards against.
                g = gamma_rhs(cs, t)  # noqa: B023
                gl = gamma_inf_mol_per_m2 * cs / (K_d_mol_per_L + cs)
                return g - gl
            # Predictor-correcter: the gap ramp endpoint is the unknown, so
            # iterate bisection until the trial value is a fixed point
            # (removes the history-lag overshoot at long times).
            cs_n = cs_hist[-1] if cs_hist else 0.0
            for _ in range(12):
                lo_b, hi_b = 0.0, c0_mol_per_L
                f_lo, f_hi = F(lo_b), F(hi_b)
                if f_lo < 0:
                    # even with a dry surface the isotherm exceeds supply --
                    # root is below the bracket; physically c_s -> 0
                    cs_new = 0.0
                elif f_hi > 0:
                    # No root in [0, c0]: the resolved history still claims
                    # supply above the monolayer capacity (the sub-grid
                    # filling transient, t_fill << t_1). Fall back to the
                    # instantaneous supply balance A(t)*(c0-cs) = langmuir(cs)
                    # with A = prefactor*sqrt(t) -- a quadratic in c_s whose
                    # positive root lies in [0, c0] and is monotone in t.
                    A = prefactor * math.sqrt(t)
                    b = (K_d_mol_per_L + gamma_inf_mol_per_m2 / A
                         - c0_mol_per_L)
                    cs_new = 0.5 * (-b + math.sqrt(b * b + 4.0 * c0_mol_per_L
                                                    * K_d_mol_per_L))
                else:
                    for _ in range(200):
                        mid = 0.5 * (lo_b + hi_b)
                        f_mid = F(mid)
                        if f_mid > 0:
                            lo_b = mid
                        else:
                            hi_b = mid
                        if hi_b - lo_b <= bisection_tol * max(c0_mol_per_L, 1e-300):
                            break
                    cs_new = 0.5 * (lo_b + hi_b)
                converged_step = abs(cs_new - cs_n) <= 1e-13 * max(c0_mol_per_L, 1e-300) + 1e-18
                cs_n = cs_new
                if converged_step:
                    break
            # Langmuir boundary condition: Gamma = langmuir(c_s) always --
            # identical to gamma_rhs at an interior root (to bisection
            # tolerance), bounded by gamma_inf when the supply-glut
            # fallback fired. Keeps the stored series BC-consistent so the
            # supply-glut overshoot cannot leak into gamma.
            g_n = gamma_inf_mol_per_m2 * cs_n / (K_d_mol_per_L + cs_n)
            if integral_first is None:
                # raw integral at the first node, with the clean initial
                # condition still in force below t_1 (must be evaluated
                # HERE: later steps see cs_at rewritten below t_1)
                integral_first = gamma_rhs(cs_n, t)
        cs_hist.append(cs_n)
        gamma_hist.append(g_n)

    gamma_eq = None
    if mode == "langmuir":
        gamma_eq = gamma_inf_mol_per_m2 * c0_mol_per_L / (K_d_mol_per_L + c0_mol_per_L)

    # early-time diagnostic: first-node Ward-Tordai integral (with the
    # solved c_s) vs the clean-surface law 2 c0 sqrt(D t / pi). At the
    # first sqrt-grid node the surface is still ~clean, so the ratio must
    # be ~1; it validates the quadrature AND the mol/L -> mol/m^3
    # conversion (prefactor*sqrt(t)*c0_mol_per_L is exactly the law in
    # mol/m^2 when c0 is mol/L). The stored langmuir gamma is the BC
    # value, not the integral, so recompute it here.
    t_first = times[0]
    free_law = prefactor * math.sqrt(t_first) * c0_mol_per_L
    if free_law > 0 and integral_first is not None:
        early_ratio = integral_first / free_law
    else:
        early_ratio = float("nan")

    monotone = all(gamma_hist[i] <= gamma_hist[i + 1] * (1 + 1e-9)
                   for i in range(len(gamma_hist) - 1))
    if gamma_eq is not None:
        bounds_ok = gamma_hist[-1] <= gamma_eq * (1 + 1e-6)
    else:
        bounds_ok = True

    return WardTordaiResult(
        t_s=tuple(times),
        gamma_mol_per_m2=tuple(gamma_hist),
        cs_mol_per_L=tuple(cs_hist),
        gamma_equilibrium_mol_per_m2=gamma_eq,
        cs_final_mol_per_L=cs_hist[-1],
        gamma_final_mol_per_m2=gamma_hist[-1],
        early_time_ratio=float(early_ratio),
        converged=bool(monotone and bounds_ok),
        method="Ward-Tordai with sqrt(t) substitution + trapezium quadrature "
               f"(n_sqrt={n_sqrt}) + per-step bisection on c_s against the "
               f"{mode} boundary condition",
        params={"c0_mol_per_L": c0_mol_per_L, "D_m2_per_s": D_m2_per_s,
                "t_end_s": t_end_s, "n_times": n_times, "mode": mode,
                "gamma_inf_mol_per_m2": gamma_inf_mol_per_m2,
                "K_d_mol_per_L": K_d_mol_per_L,
                "cs_fixed_mol_per_L": cs_fixed_mol_per_L, "n_sqrt": n_sqrt},
        caveat=(
            "Semi-infinite well-stirred bulk (c0 constant far from the "
            "surface), planar interface, diffusion control only -- no "
            "convective transport, no adsorption-barrier kinetics, no "
            "micelle breakdown contribution near the CMC. The Langmuir "
            "assumption enters only through the boundary condition."
        ),
    )
