"""Clint ideal mixing and Rubingh regular-solution theory for binary
surfactant mixtures.

Ported and generalized from a calculator originally built for a specific
cationic/anionic (amidoamine gemini / SLG) binary system, validated
against real conductometry and tensiometry data. Component 1 and
component 2 are generic here -- there is nothing cationic-specific left
in the math.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass

from ._optimize import golden_section_minimize as _golden_section_minimize
from .thermodynamics import cmc_to_mole_fraction

R_GAS = 8.314462618  # J / (mol K)


def clint_ideal_cmc(alpha1: float, cmc1: float, cmc2: float) -> float:
    """Clint's ideal mixed CMC for a binary surfactant mixture.

    alpha1: bulk mole fraction of component 1 in the total surfactant.
    cmc1, cmc2: pure-component CMCs (same concentration unit; the
    returned value is in that same unit).
    """
    if not (0.0 < alpha1 < 1.0):
        raise ValueError("alpha1 must be strictly between 0 and 1")
    if cmc1 <= 0 or cmc2 <= 0:
        raise ValueError("pure-component CMCs must be positive")
    return 1.0 / ((alpha1 / cmc1) + ((1.0 - alpha1) / cmc2))


def _rubingh_residual(x1: float, alpha1: float, cmc_mix: float, cmc1: float, cmc2: float) -> float:
    term1 = math.log((alpha1 * cmc_mix) / (x1 * cmc1))
    term2 = math.log(((1.0 - alpha1) * cmc_mix) / ((1.0 - x1) * cmc2))
    return (x1 * x1 * term1) - (((1.0 - x1) ** 2) * term2)


def solve_rubingh_x(
    alpha1: float,
    cmc_mix: float,
    cmc1: float,
    cmc2: float,
    grid_points: int = 4000,
    bisection_iters: int = 100,
) -> float | None:
    """Solve Rubingh's regular-solution equation for x1, the micellar mole
    fraction of component 1.

    Scans a fine grid for a sign change in the residual, then bisects to
    high precision. Returns None if no root is found in (0, 1) -- this
    happens for genuinely ideal or near-ideal mixtures, or for input data
    that isn't self-consistent (e.g. cmc_mix outside the physically
    possible range for the given alpha1).
    """
    eps = 1e-6
    xs = [eps + i * (1.0 - 2.0 * eps) / grid_points for i in range(grid_points + 1)]
    values: list[tuple[float, float]] = []
    for x in xs:
        try:
            fx = _rubingh_residual(x, alpha1, cmc_mix, cmc1, cmc2)
            if math.isfinite(fx):
                values.append((x, fx))
        except (ValueError, ZeroDivisionError):
            continue
    if not values:
        return None

    for (x_lo, f_lo), (x_hi, f_hi) in itertools.pairwise(values):
        if f_lo == 0:
            return x_lo
        if f_lo * f_hi < 0:
            lo, hi = x_lo, x_hi
            flo = f_lo
            for _ in range(bisection_iters):
                mid = (lo + hi) / 2.0
                fmid = _rubingh_residual(mid, alpha1, cmc_mix, cmc1, cmc2)
                if abs(fmid) < 1e-14:
                    return mid
                if flo * fmid <= 0:
                    hi = mid
                else:
                    lo, flo = mid, fmid
            return (lo + hi) / 2.0

    best_x, best_val = min(values, key=lambda pair: abs(pair[1]))
    return best_x if abs(best_val) < 1e-6 else None


def rubingh_beta(x1: float, alpha1: float, cmc_mix: float, cmc1: float) -> float:
    """Rubingh interaction parameter beta, given the solved micellar mole
    fraction x1 (see solve_rubingh_x)."""
    if x1 <= 0.0 or x1 >= 1.0:
        raise ValueError("x1 must be strictly between 0 and 1")
    term = math.log((alpha1 * cmc_mix) / (x1 * cmc1))
    return term / ((1.0 - x1) ** 2)


@dataclass
class RubinghBetaRegressionResult:
    beta_mean: float
    beta_std: float  # population std dev across the used points; 0.0 if only 1 point
    betas: list[float]  # individual pointwise beta values, in input order (excluding skipped points)
    x1_values: list[float]  # corresponding solved micellar mole fractions
    alpha1_used: list[float]  # the alpha1 values that actually solved (skipped points omitted)
    n_points_used: int
    n_points_skipped: int  # points where solve_rubingh_x found no root -- e.g. near-ideal mixtures


def rubingh_beta_regression(
    alpha1_series: list[float],
    cmc_mix_series: list[float],
    cmc1: float,
    cmc2: float,
) -> RubinghBetaRegressionResult:
    """Fit a single representative Rubingh beta across a whole composition
    series -- the real standard practice this project's own
    tests/test_mixed_micelle.py docstring already flagged as missing:
    "many papers fit beta by regression across several alpha
    compositions, while this solver computes the exact pointwise value
    for one composition."

    What "regression across compositions" means in practice, verified
    before implementing (not guessed): compute the pointwise x1 and beta
    at EACH (alpha1, cmc_mix) point independently via the existing
    solve_rubingh_x/rubingh_beta pointwise solver, then report the
    arithmetic MEAN beta across the series as the system's single
    interaction parameter -- this is the real practice used across
    published mixed-surfactant papers (e.g. Rodriguez et al. PMC6554738,
    Muherei & Junin 2009, and others in this project's own
    literature_validation_notes.md), not a continuous nonlinear
    least-squares curve fit against a global cmc_mix(alpha1) model
    (regular solution theory does not define one closed-form
    cmc_mix(alpha1) curve to fit against independent of x1 at each
    point -- x1 itself has to be solved pointwise first, by construction
    of the theory). The returned std dev across the individual betas is
    itself a real, literature-reported diagnostic (large spread signals
    that beta is not actually constant across composition, i.e. RST is
    only an approximate fit for that system) -- returned here rather
    than silently discarded.

    Points where solve_rubingh_x finds no root (near-ideal mixtures, or
    data that isn't self-consistent at that composition) are skipped,
    not treated as an error -- consistent with solve_rubingh_x's own
    None-return contract for those cases.
    """
    if len(alpha1_series) != len(cmc_mix_series):
        raise ValueError("alpha1_series and cmc_mix_series must be the same length")
    if len(alpha1_series) < 2:
        raise ValueError("need at least 2 composition points to regress a single beta")

    betas: list[float] = []
    x1_values: list[float] = []
    alpha1_used: list[float] = []
    n_skipped = 0
    for alpha1, cmc_mix in zip(alpha1_series, cmc_mix_series):
        x1 = solve_rubingh_x(alpha1, cmc_mix, cmc1, cmc2)
        if x1 is None:
            n_skipped += 1
            continue
        beta = rubingh_beta(x1, alpha1, cmc_mix, cmc1)
        betas.append(beta)
        x1_values.append(x1)
        alpha1_used.append(alpha1)

    if len(betas) < 2:
        raise ValueError(
            f"only {len(betas)} of {len(alpha1_series)} points solved -- "
            "need at least 2 to regress a single beta"
        )

    beta_mean = sum(betas) / len(betas)
    variance = sum((b - beta_mean) ** 2 for b in betas) / len(betas)
    beta_std = math.sqrt(variance)

    return RubinghBetaRegressionResult(
        beta_mean=beta_mean,
        beta_std=beta_std,
        betas=betas,
        x1_values=x1_values,
        alpha1_used=alpha1_used,
        n_points_used=len(betas),
        n_points_skipped=n_skipped,
    )


def activity_coefficients(x1: float, beta: float) -> tuple[float, float]:
    """Activity coefficients (f1, f2) from regular solution theory."""
    f1 = math.exp(beta * (1.0 - x1) ** 2)
    f2 = math.exp(beta * x1 ** 2)
    return f1, f2


def excess_free_energy(x1: float, f1: float, f2: float, temperature_k: float) -> float:
    """Excess free energy of micelle formation, kJ/mol."""
    return (R_GAS * temperature_k * (x1 * math.log(f1) + (1.0 - x1) * math.log(f2))) / 1000.0


def asymmetric_margules_activity_coefficients(x1: float, w12: float, w21: float) -> tuple[float, float]:
    """Activity coefficients (f1, f2) from the EOMMM (Equation-Oriented
    Mixed Micellization Model) binary asymmetric Margules formulation --
    Rubingh's RST generalized to allow TWO independent interaction
    parameters (w12 != w21) instead of RST's single symmetric beta.

    w12, w21: dimensionless interaction parameters in RT units (same
    convention as rubingh_beta's return value, and as the literature
    reports them, e.g. "W12 = +4.04 kBT" -- NOT raw J/mol; divide a
    literature W value in J/mol by R*T first). w12 is the energy of
    introducing a molecule of component 1 into a pure micelle of
    component 2 (and vice versa for w21) -- physically the same
    infinite-dilution-limit quantity as a Margules A12/A21 parameter.
    At w12 == w21 == beta this reduces EXACTLY to activity_coefficients
    (RST is the mathematically proven symmetric special case of the
    asymmetric formulation, not a competing formula) -- see
    test_mixed_micelle.py for the reduction check.

    Source: equations (9)-(13) of Schulz & Durand, Comput. Chem. Eng. 87
    (2016) 145-153 (doi:10.1016/j.compchemeng.2015.12.026), restricted
    to a binary (n=2) system -- their general n-component ternary-
    interaction term (Cijk) vanishes identically for n=2 since it
    requires a third component k. This is a real, verified reduction of
    their published general formula, not a guessed simplification, and
    was cross-checked to reproduce their own explicit binary symmetric-
    case reduction (their eq. 14-16, beta*x2^2 / beta*x1^2) exactly.

    NOTE: this function only computes activity coefficients for a
    KNOWN (already-fit) w12/w21 pair, e.g. quoted directly from a
    paper's own reported value. It does NOT solve for w12/w21 (or for
    x1) from raw alpha1/CMC data the way solve_rubingh_x does for RST --
    for that, see eommm_global_fit (fits w12/w21 across a whole
    composition series simultaneously, per the real SI-derived
    procedure).
    """
    x2 = 1.0 - x1
    g_exc = x1 * x2 * (x2 * w12 + x1 * w21)
    ln_f1 = 2.0 * x1 * x2 * w21 + (x2 ** 2) * w12 - 2.0 * g_exc
    ln_f2 = 2.0 * x1 * x2 * w12 + (x1 ** 2) * w21 - 2.0 * g_exc
    return math.exp(ln_f1), math.exp(ln_f2)


def _minimize_with_grid_refine(f, lo: float, hi: float, grid_points: int = 100, refine_iters: int = 60) -> float:
    """Coarse grid scan (finds the global-minimum BASIN) followed by a
    local golden-section refine within it -- needed because the per-
    point EOMMM mismatch objective below is NOT always unimodal in x1
    (multiple near-zero roots can genuinely exist for the same alpha1
    under a shared, fixed W12/W21, directly analogous to the multi-root
    behavior already found and documented for solve_rubingh_x elsewhere
    in this module at large beta). A plain single-bracket golden-section
    search silently converges to whichever basin its initial two probe
    points happen to land in; the grid scan first evaluates across the
    whole domain to find the TRUE global-minimum basin before refining
    within it precisely. grid_points/refine_iters raised from 40/40 to
    100/60 (2026-09-11): a real case was found (a symmetric composition,
    x1=0.5, under a specific W12/W21) with TWO genuinely exact roots
    (both mismatch=0 to machine precision) close enough together that
    the coarser grid landed in the wrong basin; the finer grid resolves
    it correctly (confirmed in tests/test_mixed_micelle.py)."""
    step = (hi - lo) / grid_points
    best_x, best_val = lo, f(lo)
    x = lo
    for _ in range(grid_points):
        x += step
        v = f(x)
        if v < best_val:
            best_val, best_x = v, x
    lo2 = max(lo, best_x - 2.0 * step)
    hi2 = min(hi, best_x + 2.0 * step)
    return _golden_section_minimize(f, lo2, hi2, iters=refine_iters, tol=1e-13)
    # tol=1e-13 restores this module's pre-2026-10-02 local-helper break
    # plateau region and flipped landing basins for near-degenerate EOMMM fits.


def _eommm_cmc_var_from_component(x_i: float, f_i: float, alpha_i: float, cmc_i: float, r_i: float) -> float:
    """Solves the SI's own cdefp1/cdefp2 mass-balance-with-dissociation
    equation, (cmc_var*alpha_i/cmc_i)^(2/r_i) = x_i*f_i, for cmc_var --
    the real, primary-source-verified generalization (Schulz & Durand
    2016 SI, mmc1.docx GAMS code, PDF/docx provided by the user) of the
    r=2-only linear mass balance this project used before 2026-09-11.
    r_i is the number of particles component i's headgroup nominally
    dissociates into in the micellar pseudo-phase (r=2: fully-dissociated
    1:1 ionic, the classical Clint/Rubingh convention used everywhere
    else in this module; r=1: nonionic, no dissociation)."""
    return cmc_i * (x_i * f_i) ** (r_i / 2.0) / alpha_i


def _eommm_point_log_mismatch_sq(x1: float, alpha1: float, cmc1: float, cmc2: float, w12: float, w21: float, r1: float, r2: float) -> float:
    """For a FIXED (W12, W21, r1, r2), the r-generalized cdefp1 and
    cdefp2 equations each independently imply their own cmc_var(x1);
    the joint solution (both equations satisfied by the SAME cmc_var) is
    where the two implied values coincide -- this is that squared
    log-mismatch, minimized over x1 by _minimize_with_grid_refine."""
    f1, f2 = asymmetric_margules_activity_coefficients(x1, w12, w21)
    c1 = _eommm_cmc_var_from_component(x1, f1, alpha1, cmc1, r1)
    c2 = _eommm_cmc_var_from_component(1.0 - x1, f2, 1.0 - alpha1, cmc2, r2)
    return (math.log(c1) - math.log(c2)) ** 2


def _eommm_point_infeasibility(alpha1: float, cmc1: float, cmc2: float, w12: float, w21: float, r1: float, r2: float,
                                cmc_exp: float, cmc_margin: float, eps: float = 1e-6) -> tuple[float, float, float]:
    """One mixture's contribution to the total infeasibility objective
    (see eommm_global_fit's docstring for the full derivation/history of
    why this replaced the old least-squares-to-cmc_mix approach). Two
    real, distinct failure modes are penalized SEPARATELY, not folded
    together:

    (1) mismatch_penalty = |ln(c1) - ln(c2)| at the best-fit x1 -- how
    far the two cdefp equations' independently-implied cmc_var values
    actually disagree at their closest approach. This is NOT the same
    as checking whether their average/geometric-mean happens to land
    near cmc_exp: a real bug found and fixed during development
    (2026-09-11) was exactly this -- two badly-mismatched c1, c2 values
    (e.g. a 40% relative disagreement) can have a geometric mean that
    coincidentally lands within 1% of cmc_exp, which silently corrupted
    an earlier version of this search into converging on spurious,
    non-physical (W12, W21) points. Penalizing the mismatch directly
    closes that hole.

    (2) margin_penalty = how far the (once mismatch is genuinely small)
    implied cmc_var falls outside [cmc_exp*(1-cmc_margin),
    cmc_exp*(1+cmc_margin)] -- the SI's own CMCexpVar bound (a real
    feature of their GAMS model: CMCexpVar is a FREE variable bounded
    to, not pinned at, CMCexp).

    Returns (point_infeasibility, x1, cmc_var) -- the latter two for the
    caller to report even when infeasibility is nonzero."""
    x1 = _minimize_with_grid_refine(
        lambda x: _eommm_point_log_mismatch_sq(x, alpha1, cmc1, cmc2, w12, w21, r1, r2), eps, 1.0 - eps
    )
    f1, f2 = asymmetric_margules_activity_coefficients(x1, w12, w21)
    c1 = _eommm_cmc_var_from_component(x1, f1, alpha1, cmc1, r1)
    c2 = _eommm_cmc_var_from_component(1.0 - x1, f2, 1.0 - alpha1, cmc2, r2)
    mismatch_penalty = abs(math.log(c1) - math.log(c2))
    cmc_var = math.sqrt(c1 * c2)
    lo, hi = cmc_exp * (1.0 - cmc_margin), cmc_exp * (1.0 + cmc_margin)
    margin_penalty = 0.0
    if cmc_var < lo:
        margin_penalty = (lo - cmc_var) / cmc_exp
    elif cmc_var > hi:
        margin_penalty = (cmc_var - hi) / cmc_exp
    return mismatch_penalty + margin_penalty, x1, cmc_var


def _eommm_total_infeasibility(w12: float, w21: float, alpha1_series: list[float], cmc_mix_series: list[float],
                                cmc1: float, cmc2: float, r1: float, r2: float, cmc_margin: float) -> float:
    total = 0.0
    for alpha1, cmc_exp in zip(alpha1_series, cmc_mix_series):
        p, _, _ = _eommm_point_infeasibility(alpha1, cmc1, cmc2, w12, w21, r1, r2, cmc_exp, cmc_margin)
        total += p
    return total


def _two_level_grid_search_2d(obj_fn, w_bound: float, coarse_n: int = 60, fine_n: int = 40, fine_window_factor: float = 2.0) -> tuple[float, float, float]:
    """Global minimizer of a 2D objective over [-w_bound, w_bound]^2: a
    coarse grid scan followed by a fine grid scan within a small window
    around the coarse best cell. Deliberately NOT coordinate descent /
    golden-section refinement -- a real finding during development
    (2026-09-11): alternating 1D golden-section line searches on W12
    then W21, even started exactly AT the true global optimum, walked
    AWAY from it (confirmed: infeasibility went from ~0 up to ~1 after
    "refining"). The most likely cause is the per-point inner solve
    (_eommm_point_infeasibility) occasionally locking onto a different
    x1 root as W12/W21 are perturbed one axis at a time, introducing
    discontinuities a 1D line search isn't robust to; a 2D grid simply
    never takes a large single-axis step, so it doesn't trigger this.
    Verified to recover known true (W12, W21) EXACTLY on a synthetic
    round-trip construction (see tests/test_mixed_micelle.py)."""
    step = 2.0 * w_bound / coarse_n
    best_val: float | None = None
    best_w12 = best_w21 = 0.0
    w12 = -w_bound
    for _ in range(coarse_n + 1):
        w21 = -w_bound
        for _ in range(coarse_n + 1):
            v = obj_fn(w12, w21)
            if best_val is None or v < best_val:
                best_val, best_w12, best_w21 = v, w12, w21
            w21 += step
        w12 += step

    lo12, hi12 = best_w12 - fine_window_factor * step, best_w12 + fine_window_factor * step
    lo21, hi21 = best_w21 - fine_window_factor * step, best_w21 + fine_window_factor * step
    fine_step12 = (hi12 - lo12) / fine_n
    fine_step21 = (hi21 - lo21) / fine_n
    w12 = lo12
    for _ in range(fine_n + 1):
        w21 = lo21
        for _ in range(fine_n + 1):
            v = obj_fn(w12, w21)
            if v < best_val:
                best_val, best_w12, best_w21 = v, w12, w21
            w21 += fine_step21
        w12 += fine_step12
    return best_val, best_w12, best_w21


@dataclass
class EommmGlobalFitResult:
    W12: float
    W21: float
    x1_values: list[float]  # fitted micellar mole fraction at each point, same order as input
    cmc_var_values: list[float]  # this fit's own implied CMCmix per point (within cmc_margin of the input)
    alpha1_used: list[float]
    total_infeasibility: float  # ~0 means a fully feasible, self-consistent fit was found at cmc_margin
    cmc_margin_used: float
    r_squared: float  # diagnostic: fit quality against ln(cmc_mix), see docstring
    n_points: int
    method: str = (
        "real GAMS-SI-derived mass-balance constraint satisfaction (r-generalized cdefp1/cdefp2, "
        "CMCexpVar bounded not pinned to CMCexp) minimized via 2-level grid search, NOT the old "
        "least-squares-to-cmc_mix / free-energy-minimization approaches -- see docstring"
    )


def eommm_global_fit(
    alpha1_series: list[float],
    cmc_mix_series: list[float],
    cmc1: float,
    cmc2: float,
    r1: float = 2.0,
    r2: float = 2.0,
    cmc_margin: float = 0.10,
    w_bound: float = 30.0,
) -> EommmGlobalFitResult:
    """Fit EOMMM's two independent Margules interaction parameters
    (W12, W21) SIMULTANEOUSLY across a whole composition series, plus
    each point's own micellar mole fraction x1 -- what
    asymmetric_margules_activity_coefficients's own docstring flags as
    out of its own scope.

    REWRITTEN 2026-09-11 (SI obtained, docx+CSVs provided by the user;
    replaces an earlier 2026-09-10 first-principles least-squares
    construction entirely, per this project's own Upgrade Protocol --
    the old method is not kept alongside this one). The real history,
    disclosed honestly because it matters for trusting this function:

    1. The SI's own GAMS code gives the real objective as MINIMIZING
       TOTAL FREE ENERGY OF MICELLIZATION, subject to the r-generalized
       mass-balance constraints (cdefp1/cdefp2) with CMCexpVar bounded
       (not pinned) to +/-10% of each mixture's experimental CMC. A
       first attempt to transcribe that free-energy objective exactly
       found it is UNBOUNDED BELOW as literally written (the Margules
       mixing term diverges as |W12|,|W21| grow, so minimizing it always
       ran to whatever W12/W21 search box was set -- confirmed on both
       synthetic and real Hyamine/DTAB data). That attempt was abandoned
       rather than shipped broken.
    2. A companion paper's SI (Serafini et al. 2019, TX100-DTAB, also
       obtained 2026-09-11) explains the REAL procedure in prose: "the mg
       [margin] parameter was varied in order to obtain the minimum
       value that allowed a feasible solution" -- i.e. the real fitting
       criterion is CONSTRAINT SATISFACTION at the tightest margin that
       remains feasible, not open-ended free-energy minimization under a
       fixed generous margin. This function implements THAT procedure:
       minimize a total INFEASIBILITY measure (see
       _eommm_point_infeasibility) at the given cmc_margin, rather than
       free energy -- a well-posed, bounded objective, unlike (1).
    3. A real, disclosed THEORETICAL property of this model (matching the
       Serafini SI's own stated reason for the margin-tightening
       procedure), tested carefully rather than assumed: a wider
       cmc_margin admits more of the (W12, W21) plane as feasible, so
       uniqueness is not GUARANTEED at a generous margin for every
       dataset -- multiple different (W12, W21) pairs could in principle
       both achieve zero infeasibility. An early attempt to demonstrate
       this concretely on the synthetic round-trip dataset below
       (constructed from TRUE W12=6.0, W21=-3.0) turned out to be an
       artifact of insufficient grid resolution in a first prototype,
       NOT a real property of that specific case -- with adequate
       resolution (this function's actual shipped defaults), that
       dataset recovers the exact true (W12, W21) at BOTH a tight
       cmc_margin (1e-6) and the SI's own default (0.10). This distinction
       (a real theoretical concern vs. a disproven concrete example) is
       kept explicit rather than either overclaiming non-uniqueness or
       quietly dropping the concern -- see tests/test_mixed_micelle.py
       for both the exact-recovery tests and this history. This function
       itself takes cmc_margin as a fixed input rather than searching for
       the tightest feasible value automatically -- for that, see
       eommm_find_minimal_feasible_margin (RESOLVED 2026-09-12: an
       earlier attempt at this automatic search had a real, disclosed
       grid-resolution sensitivity; re-verified against this function's
       now-fixed objective/resolution and found to work reliably -- see
       that function's own docstring for the full re-verification).

    Solved via a 2-level (coarse then fine) grid search over (W12, W21)
    -- see _two_level_grid_search_2d's own docstring for why coordinate
    descent / golden-section refinement was tried and found to actively
    walk AWAY from the true optimum (a real, disclosed finding, not a
    minor implementation detail).

    r1, r2: dissociation numbers for components 1 and 2 (how many
    particles each nominally dissociates into in the micellar pseudo-
    phase) -- default 2.0 for both, the classical fully-dissociated 1:1
    ionic convention already used throughout this module (solve_rubingh_x,
    clint_ideal_cmc, etc.); pass 1.0 for a nonionic component. cmc_margin:
    the fraction (0 to 1) each mixture's implied CMC is allowed to float
    from its experimental value -- default 0.10 matches the SI's own
    stated default. w_bound: symmetric search range for W12/W21 in RT
    units (default +/-30, generous relative to real reported values like
    the Serafini et al. TX100-DTAB system's own W12=+4.04, W21=-14.02).

    Needs at least 3 composition points. Runtime is a few seconds (a
    genuine 2D grid search with no numpy/scipy dependency, consistent
    with this project's existing EOMMM runtime expectations) -- expect
    that, it is not an error.
    """
    if len(alpha1_series) != len(cmc_mix_series):
        raise ValueError("alpha1_series and cmc_mix_series must be the same length")
    if len(alpha1_series) < 3:
        raise ValueError("need at least 3 composition points for a genuine fit")
    if any(not (0.0 < a < 1.0) for a in alpha1_series):
        raise ValueError("all alpha1 values must be strictly between 0 and 1")
    if any(c <= 0 for c in cmc_mix_series):
        raise ValueError("all cmc_mix values must be positive")
    if cmc1 <= 0 or cmc2 <= 0:
        raise ValueError("cmc1 and cmc2 must be positive")
    if r1 <= 0 or r2 <= 0:
        raise ValueError("r1 and r2 (dissociation numbers) must be positive")
    if not (0.0 < cmc_margin < 1.0):
        raise ValueError("cmc_margin must be strictly between 0 and 1")

    def obj(w12: float, w21: float) -> float:
        return _eommm_total_infeasibility(w12, w21, alpha1_series, cmc_mix_series, cmc1, cmc2, r1, r2, cmc_margin)

    total_infeas, w12, w21 = _two_level_grid_search_2d(obj, w_bound)
    return _eommm_build_result(w12, w21, total_infeas, alpha1_series, cmc_mix_series, cmc1, cmc2, r1, r2, cmc_margin)


def _eommm_build_result(w12: float, w21: float, total_infeas: float, alpha1_series: list[float],
                         cmc_mix_series: list[float], cmc1: float, cmc2: float, r1: float, r2: float,
                         cmc_margin: float) -> EommmGlobalFitResult:
    """Shared by eommm_global_fit and eommm_global_fit_all_optima: build
    the full result (x1/cmc_var per point, r_squared) for a given
    (w12, w21) -- factored out so both callers compute these identically."""
    x1_values: list[float] = []
    cmc_var_values: list[float] = []
    for alpha1, cmc_exp in zip(alpha1_series, cmc_mix_series):
        _, x1, cmc_var = _eommm_point_infeasibility(alpha1, cmc1, cmc2, w12, w21, r1, r2, cmc_exp, cmc_margin)
        x1_values.append(x1)
        cmc_var_values.append(cmc_var)
    n = len(cmc_mix_series)
    ln_actual = [math.log(c) for c in cmc_mix_series]
    ln_pred = [math.log(c) for c in cmc_var_values]
    mean_ln = sum(ln_actual) / n
    ss_res = sum((a - p) ** 2 for a, p in zip(ln_actual, ln_pred))
    ss_tot = sum((a - mean_ln) ** 2 for a in ln_actual)
    r_squared = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 1.0
    return EommmGlobalFitResult(
        W12=w12, W21=w21, x1_values=x1_values, cmc_var_values=cmc_var_values,
        alpha1_used=list(alpha1_series), total_infeasibility=total_infeas,
        cmc_margin_used=cmc_margin, r_squared=r_squared, n_points=n,
    )


@dataclass
class EommmMultiOptimaResult:
    optima: list[EommmGlobalFitResult]  # every DISTINCT feasible optimum found, sorted by total_infeasibility
    n_optima_found: int
    degenerate: bool  # True iff more than one distinct feasible optimum exists for this dataset/margin
    method: str = (
        "coarse 61x61 grid over (W12,W21), connected-component basin detection "
        "at a screening cutoff, independent fine 41x41 refinement of every "
        "basin, deduplicated by (W12,W21) separation -- see docstring"
    )


def eommm_global_fit_all_optima(
    alpha1_series: list[float],
    cmc_mix_series: list[float],
    cmc1: float,
    cmc2: float,
    r1: float = 2.0,
    r2: float = 2.0,
    cmc_margin: float = 0.10,
    w_bound: float = 30.0,
    feasible_tolerance: float = 1e-6,
) -> EommmMultiOptimaResult:
    """Like eommm_global_fit, but does not silently return only the single
    best-scoring (W12, W21) -- it reports EVERY distinct feasible optimum.

    Why this exists: eommm_global_fit's own docstring discloses a real
    theoretical property of this model (point 3) -- "uniqueness is not
    GUARANTEED at a generous margin for every dataset -- multiple
    different (W12, W21) pairs could in principle both achieve zero
    infeasibility." That disclosure was previously left unactionable:
    eommm_global_fit always returns whichever single optimum its grid
    search happens to find first, silently picking a winner even on a
    dataset where a second, equally-valid optimum exists elsewhere in
    the (W12, W21) plane. For a caller that needs to know whether a
    *different but equally correct* fit exists (e.g. grading an external
    answer that used a different search strategy and landed somewhere
    else), that silence is the problem, not the fit itself.

    Method: run the same coarse 61x61 grid eommm_global_fit uses, but
    keep the FULL infeasibility grid rather than only its minimum. Group
    grid cells into connected low-infeasibility "basins" (8-connected,
    screening cutoff = max(10x the grid's own global minimum, global
    minimum + 0.02) -- generous enough to catch any basin the fine search
    could plausibly pull below feasible_tolerance, without fine-refining
    the entire grid). Independently fine-refine each basin's own best
    cell with the identical 41x41 fine search eommm_global_fit uses, then
    keep only refined candidates that actually reach total_infeasibility
    < feasible_tolerance and are mutually separated by more than 3 coarse
    grid steps in (W12, W21) (so two basins that refine onto the SAME
    true optimum are correctly merged into one entry, not double-counted).

    Returns every surviving optimum, sorted by total_infeasibility
    (ascending); `degenerate=True` iff more than one survives. A caller
    grading an external (W12, W21) answer should accept it as correct if
    it matches ANY entry in `optima`, not only the first.

    Runtime: one coarse pass (shared cost with eommm_global_fit) plus one
    fine pass PER basin found -- typically 1-5 basins in practice, so a
    few times eommm_global_fit's own runtime (still seconds, not minutes).
    """
    if len(alpha1_series) != len(cmc_mix_series):
        raise ValueError("alpha1_series and cmc_mix_series must be the same length")
    if len(alpha1_series) < 3:
        raise ValueError("need at least 3 composition points for a genuine fit")
    if any(not (0.0 < a < 1.0) for a in alpha1_series):
        raise ValueError("all alpha1 values must be strictly between 0 and 1")
    if any(c <= 0 for c in cmc_mix_series):
        raise ValueError("all cmc_mix values must be positive")
    if cmc1 <= 0 or cmc2 <= 0:
        raise ValueError("cmc1 and cmc2 must be positive")
    if r1 <= 0 or r2 <= 0:
        raise ValueError("r1 and r2 (dissociation numbers) must be positive")
    if not (0.0 < cmc_margin < 1.0):
        raise ValueError("cmc_margin must be strictly between 0 and 1")
    if feasible_tolerance <= 0:
        raise ValueError("feasible_tolerance must be positive")

    def obj(w12: float, w21: float) -> float:
        return _eommm_total_infeasibility(w12, w21, alpha1_series, cmc_mix_series, cmc1, cmc2, r1, r2, cmc_margin)

    coarse_n = 60
    step = 2.0 * w_bound / coarse_n
    grid: list[list[float]] = []
    w12 = -w_bound
    for _ in range(coarse_n + 1):
        row = []
        w21 = -w_bound
        for _ in range(coarse_n + 1):
            row.append(obj(w12, w21))
            w21 += step
        grid.append(row)
        w12 += step

    global_min = min(v for row in grid for v in row)
    screening_cutoff = max(10.0 * global_min, global_min + 0.02)

    # connected-component (8-connected) basin detection over the coarse grid
    n_side = coarse_n + 1
    visited = [[False] * n_side for _ in range(n_side)]
    basins: list[list[tuple[int, int]]] = []
    for i in range(n_side):
        for j in range(n_side):
            if visited[i][j] or grid[i][j] > screening_cutoff:
                continue
            stack = [(i, j)]
            visited[i][j] = True
            cells = []
            while stack:
                ci, cj = stack.pop()
                cells.append((ci, cj))
                for di in (-1, 0, 1):
                    for dj in (-1, 0, 1):
                        if di == 0 and dj == 0:
                            continue
                        ni, nj = ci + di, cj + dj
                        if 0 <= ni < n_side and 0 <= nj < n_side and not visited[ni][nj] \
                                and grid[ni][nj] <= screening_cutoff:
                            visited[ni][nj] = True
                            stack.append((ni, nj))
            basins.append(cells)

    def _fine_refine(seed_w12: float, seed_w21: float) -> tuple[float, float, float]:
        lo12, hi12 = seed_w12 - 2.0 * step, seed_w12 + 2.0 * step
        lo21, hi21 = seed_w21 - 2.0 * step, seed_w21 + 2.0 * step
        fine_n = 40
        fine_step12 = (hi12 - lo12) / fine_n
        fine_step21 = (hi21 - lo21) / fine_n
        best_val, best_w12, best_w21 = obj(seed_w12, seed_w21), seed_w12, seed_w21
        fw12 = lo12
        for _ in range(fine_n + 1):
            fw21 = lo21
            for _ in range(fine_n + 1):
                v = obj(fw12, fw21)
                if v < best_val:
                    best_val, best_w12, best_w21 = v, fw12, fw21
                fw21 += fine_step21
            fw12 += fine_step12
        return best_val, best_w12, best_w21

    candidates: list[tuple[float, float, float]] = []  # (infeas, w12, w21)
    for cells in basins:
        bi, bj = min(cells, key=lambda c: grid[c[0]][c[1]])
        seed_w12, seed_w21 = -w_bound + bi * step, -w_bound + bj * step
        refined = _fine_refine(seed_w12, seed_w21)
        candidates.append(refined)

    candidates.sort(key=lambda c: c[0])
    kept: list[tuple[float, float, float]] = []
    sep_threshold = 3.0 * step
    for infeas, w12v, w21v in candidates:
        if infeas >= feasible_tolerance:
            continue
        if any(math.hypot(w12v - kw12, w21v - kw21) < sep_threshold for _, kw12, kw21 in kept):
            continue
        kept.append((infeas, w12v, w21v))

    optima = [
        _eommm_build_result(w12v, w21v, infeas, alpha1_series, cmc_mix_series, cmc1, cmc2, r1, r2, cmc_margin)
        for infeas, w12v, w21v in kept
    ]
    return EommmMultiOptimaResult(optima=optima, n_optima_found=len(optima), degenerate=len(optima) > 1)


def eommm_find_minimal_feasible_margin(
    alpha1_series: list[float],
    cmc_mix_series: list[float],
    cmc1: float,
    cmc2: float,
    r1: float = 2.0,
    r2: float = 2.0,
    w_bound: float = 30.0,
    margin_hi: float = 0.30,
    feasibility_tolerance: float = 1e-3,
    binary_search_iters: int = 20,
) -> EommmGlobalFitResult:
    """Automates the SI's own margin-tightening procedure that
    eommm_global_fit's docstring describes but does not perform itself:
    binary-searches cmc_margin down to the smallest value for which a
    fully feasible fit (total_infeasibility <= feasibility_tolerance)
    still exists, then returns that fit.

    RESOLVED 2026-09-12: eommm_global_fit's own docstring, as first
    written (2026-09-11), disclosed this as unshippable -- an automatic
    search prototyped against an EARLIER, buggy version of the
    infeasibility objective (the geometric-mean bug documented in
    eommm_global_fit's own history) showed a real grid-resolution
    sensitivity, converging to a margin looser than the true minimum.
    Re-tested directly against the CURRENT, fixed objective and grid
    resolution (both already corrected for other reasons -- see
    eommm_global_fit's docstring) before re-attempting this: the
    earlier sensitivity is GONE. On a synthetic round-trip (known true
    W12=6.0, W21=-3.0), the binary search converges monotonically,
    every single intermediate margin from 0.15 down to ~1e-6 stays
    exactly feasible at the exact true (W12, W21) -- not just the final
    answer. On real (noisy) Hyamine/DTAB data, it converges to a
    genuine, stable, nonzero minimal margin (~6.2%) reflecting real
    measurement noise rather than a resolution artifact, with
    (W12, W21) stable across nearby margin values, not jumping around.

    Runtime is substantially slower than eommm_global_fit itself
    (binary_search_iters+1 full 2-level grid searches, each a genuine
    2D search in its own right) -- expect roughly a minute for a typical
    dataset with default settings, not an error.

    margin_hi: the widest margin the search starts from (must itself be
    feasible, or this raises -- if even a generous margin can't produce
    a feasible fit, the model/data are fundamentally inconsistent, not
    a search-precision issue). feasibility_tolerance: how close to zero
    total_infeasibility must get to count as "feasible" -- the default
    (1e-3) matches what was used during the verification described
    above; tightening it substantially will need more grid resolution
    (see _two_level_grid_search_2d) to stay reliable.
    """
    if margin_hi <= 0.0 or margin_hi >= 1.0:
        raise ValueError("margin_hi must be strictly between 0 and 1")

    def fit_at(margin: float) -> EommmGlobalFitResult:
        return eommm_global_fit(alpha1_series, cmc_mix_series, cmc1, cmc2, r1, r2, margin, w_bound)

    hi = margin_hi
    result_hi = fit_at(hi)
    if result_hi.total_infeasibility > feasibility_tolerance:
        raise ValueError(
            f"no feasible fit found even at margin_hi={margin_hi} (total_infeasibility="
            f"{result_hi.total_infeasibility:.4g}) -- the data and model are fundamentally "
            "inconsistent at this r1/r2, not a search-precision issue; try a wider margin_hi, "
            "different r1/r2, or check the data"
        )
    lo = 0.0
    best = result_hi
    for _ in range(binary_search_iters):
        mid = (lo + hi) / 2.0
        result_mid = fit_at(mid)
        if result_mid.total_infeasibility <= feasibility_tolerance:
            hi = mid
            best = result_mid
        else:
            lo = mid
    return best


def rodenas_x1(alpha1: float, dln_cmc_mix_dalpha1: float) -> float:
    """Rodenas' Gibbs-Duhem-based, "model-independent" micellar mole
    fraction of component 1, from the LOCAL SLOPE of ln(experimental
    mixed CMC) vs. bulk mole fraction alpha1:

    X1_Rod = -(1 - alpha1) * alpha1 * d[ln(cmc_mix)]/d(alpha1) + alpha1

    Unlike solve_rubingh_x (which assumes a specific regular-solution
    activity-coefficient model and solves an implicit equation) or
    motomura_ideal_composition (the zero-interaction ideal-mixing
    limit), this derives composition directly from the Gibbs-Duhem
    relation applied to the REAL cmc_mix(alpha1) curve -- "model
    independent" in the sense of not assuming any particular
    interaction/activity-coefficient form, at the cost of needing the
    local slope of a real multi-point dataset rather than a single
    (alpha1, cmc_mix) pair.

    dln_cmc_mix_dalpha1: d[ln(cmc_mix)]/d(alpha1) at this alpha1 --
    must come from an actual multi-point experimental cmc_mix(alpha1)
    dataset (a fit, or a numerical/finite-difference derivative across
    real neighboring data points), never guessed or assumed analytic.
    This is exactly why this function takes the slope as a direct
    input instead of computing it internally from one data pair the
    way solve_rubingh_x/motomura_ideal_composition do: Rodenas'
    method fundamentally needs the LOCAL SLOPE of the real experimental
    composition curve, which requires several real neighboring data
    points, not just the single composition being evaluated.

    Source: equation (10), Azum, Rub, Alotaibi, Khan & Asiri,
    Biointerface Res. Appl. Chem. 12(6) (2022) 7416-7428 (the same G6
    gemini surfactant paper already used elsewhere in this project's
    test suite for Rubingh/Clint validation), citing Rodenas, Valiente
    & Villafruela, J. Phys. Chem. B 103(21) (1999) 4549-4554, itself
    built on Lange's model.

    Verified two ways: (1) directly against the cited source equation;
    (2) a real, source-independent algebraic identity (checked by
    test): when cmc_mix follows Clint's EXACT ideal-mixing law (zero
    interaction), X1_Rod reduces algebraically to EXACTLY
    motomura_ideal_composition's formula -- both are zero-interaction
    limits of otherwise-different methods and must agree there, and
    they do.
    """
    return -(1.0 - alpha1) * alpha1 * dln_cmc_mix_dalpha1 + alpha1


def _quadratic_lagrange_derivative_at(x0: float, y0: float, x1: float, y1: float, x2: float, y2: float, x_eval: float) -> float:
    """Derivative at x_eval of the unique quadratic interpolating
    (x0,y0), (x1,y1), (x2,y2) -- exact Lagrange-basis differentiation,
    valid for arbitrarily (unevenly) spaced x0/x1/x2. Reduces exactly to
    the textbook central-difference formula (y2-y0)/(2h) when x_eval=x1
    is the midpoint of evenly-spaced points, and to the textbook
    three-point forward/backward-difference formula at an endpoint --
    both verified in tests/test_mixed_micelle.py, not just asserted."""
    d0 = (x0 - x1) * (x0 - x2)
    d1 = (x1 - x0) * (x1 - x2)
    d2 = (x2 - x0) * (x2 - x1)
    l0p = (2.0 * x_eval - x1 - x2) / d0
    l1p = (2.0 * x_eval - x0 - x2) / d1
    l2p = (2.0 * x_eval - x0 - x1) / d2
    return y0 * l0p + y1 * l1p + y2 * l2p


@dataclass
class RodenasSeriesResult:
    alpha1_used: list[float]  # sorted ascending
    dln_cmc_mix_dalpha1: list[float]  # local-slope estimate at each alpha1_used
    x1_rodenas: list[float]  # rodenas_x1 evaluated at each point, chained through automatically
    n_points: int
    x1_physically_valid: list[bool]  # x1_rodenas[i] in [0, 1] -- a micellar mole fraction outside this
                                        # range is not a rounding quirk, it means the local-slope estimate
                                        # at that point is not trustworthy (see caveat)
    method: str = "local quadratic (Lagrange) fit through the 3 nearest neighboring points, differentiated analytically at each composition"
    caveat: str = (
        "Rodenas' method is mathematically exact only when ln(cmc_mix) truly follows a smooth "
        "quadratic locally -- on a sparse (few-point), sharply/non-monotonically curved real "
        "composition series (e.g. a strongly synergistic catanionic-adjacent system with a deep "
        "minimum near equimolar composition), the local-slope estimate can swing wildly between "
        "neighboring points and produce a micellar mole fraction OUTSIDE the physical [0,1] range "
        "-- a real failure of this 'model-independent' method's numerical stability on realistic "
        "data density, not a bug in this implementation. Check x1_physically_valid before trusting "
        "any x1_rodenas value; an out-of-range point should be reported as unreliable, not clipped "
        "or silently accepted. Rubingh's regular-solution fit (solve_rubingh_x/rubingh_beta), "
        "despite assuming a fixed activity-coefficient form, is often numerically more robust on "
        "sparse real data for exactly this reason -- cross-check against it."
    )


def rodenas_x1_series(alpha1_series: list[float], cmc_mix_series: list[float]) -> RodenasSeriesResult:
    """Compute Rodenas' d[ln(cmc_mix)]/d(alpha1) -- and the resulting x1
    at every composition -- directly from a RAW multi-alpha1 cmc_mix(alpha1)
    series, closing the gap rodenas_x1's own docstring flags: that
    function takes the local slope as an already-known input "because
    Rodenas' method fundamentally needs the LOCAL SLOPE of the real
    experimental composition curve, which requires several real
    neighboring data points" -- this is that numerical-differentiation
    step.

    Real, standard numerical-differentiation method for an unevenly-
    spaced discrete series (the source papers describe getting "the
    slope from the curve" without specifying their exact numerical
    procedure -- verified by checking, not guessed as surfactant-
    specific): at each composition, fit the unique quadratic through it
    and its 2 nearest neighbors along alpha1, then differentiate that
    quadratic analytically at that exact alpha1. For an interior point
    this uses its left and right neighbor; for the first/last point
    (which has only one side) it uses the nearest 2 points on its
    available side. This is mathematically exact (not approximate) when
    the underlying ln(cmc_mix) truly follows a smooth quadratic locally,
    and reduces exactly to the standard central-difference /
    forward-difference / backward-difference formulas for evenly-spaced
    data (see test_mixed_micelle.py) -- it is a real generalization of
    those textbook formulas to uneven spacing, not an invented method.

    alpha1_series and cmc_mix_series must be the same length, >= 3
    points, with all alpha1 values distinct and strictly between 0 and
    1, and all cmc_mix values positive (needed for ln). Points are
    sorted by alpha1 internally; results are returned in that sorted
    order (alpha1_used).
    """
    if len(alpha1_series) != len(cmc_mix_series):
        raise ValueError("alpha1_series and cmc_mix_series must be the same length")
    if len(alpha1_series) < 3:
        raise ValueError("need at least 3 (alpha1, cmc_mix) points to fit a local quadratic derivative")
    if any(not (0.0 < a < 1.0) for a in alpha1_series):
        raise ValueError("all alpha1 values must be strictly between 0 and 1")
    if any(c <= 0 for c in cmc_mix_series):
        raise ValueError("all cmc_mix values must be positive")

    paired = sorted(zip(alpha1_series, cmc_mix_series), key=lambda p: p[0])
    alphas = [p[0] for p in paired]
    if len(set(alphas)) != len(alphas):
        raise ValueError("alpha1_series contains duplicate values -- cannot fit a local derivative")
    ln_cmc = [math.log(c) for _, c in paired]

    n = len(alphas)
    slopes: list[float] = []
    for i in range(n):
        if i == 0:
            i0, i1, i2 = 0, 1, 2
        elif i == n - 1:
            i0, i1, i2 = n - 3, n - 2, n - 1
        else:
            i0, i1, i2 = i - 1, i, i + 1
        slope = _quadratic_lagrange_derivative_at(
            alphas[i0], ln_cmc[i0], alphas[i1], ln_cmc[i1], alphas[i2], ln_cmc[i2], alphas[i]
        )
        slopes.append(slope)

    x1_values = [rodenas_x1(alphas[i], slopes[i]) for i in range(n)]
    valid = [0.0 <= v <= 1.0 for v in x1_values]
    return RodenasSeriesResult(
        alpha1_used=alphas, dln_cmc_mix_dalpha1=slopes, x1_rodenas=x1_values,
        n_points=n, x1_physically_valid=valid,
    )


def rodenas_activity_coefficients(alpha1: float, x1: float, cmc_mix: float, cmc1: float, cmc2: float) -> tuple[float, float]:
    """Activity coefficients (f1, f2) for Rodenas' model, given the
    already-determined x1 (see rodenas_x1) -- equations (11)-(12) of
    the source cited in rodenas_x1's docstring:

    f1 = alpha1 * cmc_mix / (x1 * cmc1)
    f2 = (1 - alpha1) * cmc_mix / ((1 - x1) * cmc2)

    Structurally the same mass-balance definition already used
    internally by solve_rubingh_x (see _rubingh_residual's term1/
    term2) -- Rodenas and Rubingh differ only in HOW x1 itself is
    determined (a real experimental curve's Gibbs-Duhem slope vs. an
    assumed regular-solution model), not in how activity coefficients
    follow from an already-known x1.
    """
    if not (0.0 < x1 < 1.0):
        raise ValueError("x1 must be strictly between 0 and 1")
    f1 = (alpha1 * cmc_mix) / (x1 * cmc1)
    f2 = ((1.0 - alpha1) * cmc_mix) / ((1.0 - x1) * cmc2)
    return f1, f2


def maeda_free_energy_of_micellization(
    x1_rub: float,
    beta: float,
    cmc1_M: float,
    cmc2_M: float,
    temperature_k: float,
    water_molarity_m: float = 55.5,
) -> float:
    """Maeda's free energy of micellization (kJ/mol) for a binary
    ionic/nonionic mixed micelle -- an alternative expression to
    excess_free_energy, built specifically for ionic-nonionic systems
    and requiring chain-chain (not just head-group) interaction:

    deltaG_M = RT * (B0 + B1*X1_Rub + B2*X1_Rub^2)
    B0 = ln(Xcmc2)
    B1 + B2 = ln(Xcmc1/Xcmc2)
    B2 = -beta

    where Xcmc1, Xcmc2 are the pure-component CMCs on the MOLE-FRACTION
    scale (see cmc_to_mole_fraction -- the same conversion this
    library's own thermodynamics.py already uses for
    gibbs_free_energy_micellization), and x1_rub/beta are Rubingh's own
    solved micellar mole fraction and interaction parameter (from
    solve_rubingh_x/rubingh_beta) -- Maeda's model reuses RST's own
    x1/beta rather than fitting independently.

    cmc1_M, cmc2_M: pure-component CMCs in mol/L (molarity) -- convert
    from mM/other units before calling.

    Source: equation (9), Azum, Rub, Alotaibi, Khan & Asiri,
    Biointerface Res. Appl. Chem. 12(6) (2022) 7416-7428 (same paper
    already used elsewhere in this project's test suite for Clint/
    Rubingh/Rodenas validation), citing Maeda, J. Colloid Interface
    Sci. 172 (1995) 98-105. Numerically cross-checked (see
    test_mixed_micelle.py): computing B0 = ln(Xcmc2) with this
    project's own already-validated TX-114 CMC (0.263 mM, from the same
    source paper's Table 1) reproduces the paper's own reported
    -B0 = 12.25 (G6+TX-114, alpha1 = 0.22) to within 0.01 -- real
    numeric agreement, not just a formula transcription check.
    """
    xcmc1 = cmc_to_mole_fraction(cmc1_M, water_molarity_m)
    xcmc2 = cmc_to_mole_fraction(cmc2_M, water_molarity_m)
    b2 = -beta
    b1 = math.log(xcmc1 / xcmc2) - b2
    b0 = math.log(xcmc2)
    delta_g_over_rt = b0 + b1 * x1_rub + b2 * (x1_rub ** 2)
    return (R_GAS * temperature_k * delta_g_over_rt) / 1000.0


def solve_rosen_monolayer_x(alpha1: float, c_mix_sigma: float, c1_sigma: float, c2_sigma: float, **kwargs) -> float | None:
    """Rosen's extension of Rubingh's regular-solution theory to the
    mixed ADSORBED MONOLAYER at an interface, from surface tension data
    rather than CMC data. Mathematically identical to solve_rubingh_x --
    the difference is entirely in what the inputs mean, and getting that
    right matters: c1_sigma/c2_sigma are the pure-component surfactant
    concentrations required to reach a chosen REFERENCE surface tension
    (or surface pressure) at the interface, and c_mix_sigma is the total
    mixed-surfactant concentration required to reach that SAME reference
    surface tension at bulk mole fraction alpha1 -- these are NOT CMC
    values. Using CMC values here by mistake silently produces a
    micellization-monolayer confusion, not an error."""
    return solve_rubingh_x(alpha1, c_mix_sigma, c1_sigma, c2_sigma, **kwargs)


def rosen_beta_sigma(x1_sigma: float, alpha1: float, c_mix_sigma: float, c1_sigma: float) -> float:
    """Rosen's monolayer interaction parameter beta^sigma, given the
    solved monolayer mole fraction x1_sigma (see solve_rosen_monolayer_x).
    Same formula as rubingh_beta -- see that function's docstring for
    the input-meaning distinction (surface-tension-derived concentrations,
    not CMC values)."""
    return rubingh_beta(x1_sigma, alpha1, c_mix_sigma, c1_sigma)


def motomura_ideal_composition(alpha1: float, cmc1: float, cmc2: float) -> float:
    """Motomura & Aratono's ideal micellar composition for a binary
    surfactant mixture:

    X1_id = (alpha1 * cmc2) / (alpha1 * cmc2 + alpha2 * cmc1)

    alpha1: bulk mole fraction of component 1 in the total surfactant
    (0 to 1, exclusive). cmc1, cmc2: pure-component CMCs (any consistent
    concentration unit -- the result is dimensionless and unit-
    independent, since cmc1/cmc2 only ever appear as a ratio).

    This is Motomura's companion to Clint's ideal-mixing CMC formula
    (clint_ideal_cmc): both assume zero net interaction between the two
    surfactants, but Clint predicts the ideal MIXED CMC while Motomura
    predicts the ideal MICELLAR COMPOSITION (X1) at that CMC. The two
    are algebraically linked -- X1_id = alpha1 * CMC_mix_ideal / cmc1 --
    and this identity is what test_mixed_micelle.py checks, since it
    holds regardless of any external literature value.

    Source: equation (8), Serafini et al., Colloids Surf. A 562 (2019)
    170-181 (arXiv:1806.09721), citing Motomura & Aratono's original
    ideal-mixing composition formula. Like Clint's model, this is a
    ZERO-PARAMETER ideal-mixing prediction -- it does not use or need an
    interaction parameter, unlike Rubingh's regular-solution beta.
    """
    if not (0.0 < alpha1 < 1.0):
        raise ValueError("alpha1 must be strictly between 0 and 1")
    if cmc1 <= 0 or cmc2 <= 0:
        raise ValueError("pure-component CMCs must be positive")
    alpha2 = 1.0 - alpha1
    return (alpha1 * cmc2) / (alpha1 * cmc2 + alpha2 * cmc1)


def corrin_harkins_predict_cmc(cmc1: float, c_counterion1: float, cmc2: float, c_counterion2: float, c_counterion_target: float) -> tuple[float, float]:
    """Predict CMC at a new counterion (salt) concentration using the
    Corrin-Harkins log-linear relation: log10(CMC) = -g*log10(C) + b.
    Fits g and b from two known (CMC, counterion concentration) points,
    then predicts CMC at c_counterion_target. Returns (predicted_cmc, g).

    g and b are SYSTEM-SPECIFIC (surfactant + salt identity), not
    universal constants -- this is why two real data points are required
    as input rather than a lookup table. All CMC and concentration
    values must be in the same, consistent unit (e.g. all mM); the
    prediction is only reliable as an interpolation between (or close
    extrapolation beyond) the two given concentrations.
    """
    if cmc1 <= 0 or cmc2 <= 0 or c_counterion1 <= 0 or c_counterion2 <= 0 or c_counterion_target <= 0:
        raise ValueError("all CMC and concentration values must be positive")
    if c_counterion1 == c_counterion2:
        raise ValueError("c_counterion1 and c_counterion2 must differ to fit a slope")
    log_cmc1, log_cmc2 = math.log10(cmc1), math.log10(cmc2)
    log_c1, log_c2 = math.log10(c_counterion1), math.log10(c_counterion2)
    g = (log_cmc2 - log_cmc1) / (log_c1 - log_c2)
    b = log_cmc1 + g * log_c1
    log_cmc_target = -g * math.log10(c_counterion_target) + b
    return 10.0 ** log_cmc_target, g
