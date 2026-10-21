"""Dynamic surface tension (DST): the time-dependent surface tension of a
FRESHLY CREATED interface, gamma(t), before equilibrium is reached --
distinct from every other module in this library, which describes
equilibrium (infinite-time) surface tension. DST governs any process
operating on short timescales where equilibrium gamma predicts poorly:
wetting speed, foam generation during agitation, emulsification during
high-shear mixing (Ritacco, Langevin, Diamant & Andelman, arXiv:1008.1230v2
-- an already-trusted real source in this project, used for GZ-17's own
gold CMC value).

Built 2026-09-25, closing a real, previously-total gap (grep confirmed
zero prior DST presence anywhere in this library) after a deep literature
research pass. The real underlying theory (Ward & Tordai, 1946) has NO
general closed-form solution -- an integral equation requiring the
adsorption isotherm as a closing relation, normally solved numerically.
This module does NOT attempt that numerical solver (a real, disclosed
scope limit, not a silently-skipped gap). Instead it implements the real,
closed-form, genuinely fittable ASYMPTOTIC regimes the same primary source
derives and validates against real DTAB tensiometry data:

  - Diffusion-Limited-Adsorption (DLA) regime: fit_dla_regime.
  - Kinetically-Limited-Adsorption (KLA) regime (a real, distinct, slower
    relaxation seen in real ionic-surfactant data when an adsorption
    barrier, e.g. electrostatic, is present): select_dst_kinetic_regime
    fits this ON TOP OF the DLA regime and decides via BIC whether the
    extra 2 parameters are genuinely earned by the data -- the direct
    structural analogue of adsorption.select_isotherm_model's
    Langmuir-vs-Frumkin decision and curve_analysis.cmc_from_surface_
    tension_curve's own 2-vs-3-segment decision. Same principle, new
    curve type.
  - Rosen-Hua equation: fit_rosen_hua, the standard practical 3-parameter
    empirical curve-shape fit (Hua & Rosen, J. Colloid Interface Sci.
    (1988); Rosen's "Surfactants and Interfacial Phenomena") for when a
    researcher just wants to characterize a measured gamma(t) curve's
    relaxation time/shape without invoking the electrostatic-barrier
    theory above -- genuinely complementary, not redundant.

Real, disclosed limitation carried through every function here: converting
a fitted DLA time constant into an actual diffusion coefficient requires
the equilibrium Langmuir isotherm parameters (theta_0_eq, i.e. Gamma at
this concentration; beta_a, the Langmuir adsorption parameter) from a
SEPARATELY-fitted equilibrium curve on the same system (adsorption.py's
szyszkowski_fit_K/select_isotherm_model) -- without those, only the raw
time constant itself is reported, not a diffusion coefficient. This module
does not guess a diffusion coefficient when that equilibrium data is
missing.
"""

from __future__ import annotations
import math
from dataclasses import dataclass, field

from ._optimize import golden_section_minimize as _golden_section_minimize
from ._stats import bic as _bic


@dataclass
class DlaFitResult:
    tau_d_s: float
    delta_gamma_eq_mN_m: float
    r_squared: float
    n_points: int
    diffusion_coefficient_m2_per_s: float | None = None
    gap: str | None = None


def fit_dla_regime(
    surface_age_s: list[float],
    gamma_mN_m: list[float],
    gamma0_mN_m: float,
    theta_0_eq: float | None = None,
    beta_a_per_mM: float | None = None,
    a_molecular_size_m: float | None = None,
) -> DlaFitResult:
    """Fit the diffusion-limited-adsorption (DLA) regime to real
    surface-age-vs-tension data:

        Delta_gamma(t) = Delta_gamma_eq * (1 - sqrt(tau_d/t))

    where Delta_gamma(t) = gamma0 - gamma(t) (the tension DROP from the
    pure-solvent/fresh-interface value). Real, closed-form asymptotic
    solution of the Ward-Tordai diffusion equation (Ritacco et al. 2010,
    arXiv:1008.1230v2, eq. 2), valid for t much greater than the molecular
    relaxation time. tau_d (diffusion time) and Delta_gamma_eq (equilibrium
    tension drop implied by this DST curve alone) are the two free
    parameters, fit here by golden-section search over tau_d with a
    closed-form ordinary-least-squares Delta_gamma_eq at each candidate
    (the model is LINEAR in Delta_gamma_eq for any fixed tau_d) -- the
    same nested closed-form-inner/golden-section-outer pattern this
    project already uses (adsorption.py's szyszkowski_fit_K,
    curve_analysis.py's _fit_szyszkowski).

    surface_age_s: real surface age (seconds) at each measurement -- NOTE
    for the maximum-bubble-pressure method specifically: raw bubble-
    detachment intervals typically need an apparatus-specific correction
    to the true surface age (t = t_raw / beta^2, calibrated against a
    known-diffusion-coefficient reference compound) -- this function
    assumes surface_age_s is ALREADY the true, corrected surface age, not
    a raw instrument reading; do not pass raw bubble intervals directly
    without that correction if using MBP data.
    gamma_mN_m: measured surface tension (mN/m) at each surface age.
    gamma0_mN_m: pure-solvent (t=0, freshest-possible-interface) surface
    tension.

    theta_0_eq, beta_a_per_mM, a_molecular_size_m: OPTIONAL. When all
    three are supplied (from a SEPARATELY-fitted equilibrium isotherm on
    the same system -- theta_0_eq is the equilibrium surface coverage at
    this bulk concentration, beta_a_per_mM is the Langmuir adsorption
    parameter, a_molecular_size_m is the molecule's effective size), the
    real diffusion coefficient is derived via the near-saturated-monolayer
    limit (Ritacco et al. 2010, eq. 4): D = a^2 / (beta_a^2 * tau_d). When
    any of the three is missing, diffusion_coefficient_m2_per_s is None
    and `gap` explicitly states why -- this function does NOT guess a
    diffusion coefficient from tau_d alone.
    """
    if len(surface_age_s) != len(gamma_mN_m):
        raise ValueError("surface_age_s and gamma_mN_m must be the same length")
    if len(surface_age_s) < 3:
        raise ValueError("need at least 3 (surface age, surface tension) points to fit tau_d and Delta_gamma_eq")
    if any(t <= 0 for t in surface_age_s):
        raise ValueError("surface_age_s values must all be positive (t=0 is the singular pure-solvent limit, not a fittable data point)")

    delta_gamma = [gamma0_mN_m - g for g in gamma_mN_m]
    t = list(surface_age_s)
    n = len(t)

    def delta_gamma_eq_and_rss(tau_d: float) -> tuple[float, float]:
        x = [1.0 - math.sqrt(tau_d / ti) for ti in t]
        sxx = sum(xi * xi for xi in x)
        if sxx == 0:
            return 0.0, sum(dg * dg for dg in delta_gamma)
        sxy = sum(xi * dg for xi, dg in zip(x, delta_gamma))
        deg = sxy / sxx
        rss = sum((dg - deg * xi) ** 2 for dg, xi in zip(delta_gamma, x))
        return deg, rss

    # RSS(tau_d) is NOT guaranteed unimodal across a wide search range (verified:
    # a real secondary local minimum can appear for tau_d orders of magnitude away
    # from the true value) -- golden-section alone can converge to the wrong
    # region. Same fix as adsorption.py's frumkin_fit_K_and_a: a coarse grid
    # search over ln(tau_d) first to find the right neighborhood, then
    # golden-section refinement within it.
    t_min, t_max = min(t), max(t)
    ln_lo, ln_hi = math.log(t_min * 1e-4), math.log(t_max * 10.0)
    n_grid = 60
    grid_taus = [math.exp(ln_lo + i * (ln_hi - ln_lo) / (n_grid - 1)) for i in range(n_grid)]
    grid_rss = [(td, delta_gamma_eq_and_rss(td)[1]) for td in grid_taus]
    best_idx = min(range(len(grid_rss)), key=lambda i: grid_rss[i][1])
    tau_d_grid_best = grid_rss[best_idx][0]

    step = (ln_hi - ln_lo) / (n_grid - 1)
    refine_lo = math.exp(max(ln_lo, math.log(tau_d_grid_best) - step))
    refine_hi = math.exp(min(ln_hi, math.log(tau_d_grid_best) + step))
    tau_d_best = _golden_section_minimize(lambda td: delta_gamma_eq_and_rss(td)[1], refine_lo, refine_hi)
    delta_gamma_eq, rss = delta_gamma_eq_and_rss(tau_d_best)

    mean_dg = sum(delta_gamma) / n
    ss_tot = sum((dg - mean_dg) ** 2 for dg in delta_gamma)
    r_squared = 1.0 - (rss / ss_tot) if ss_tot > 0 else 1.0

    diffusion_coefficient = None
    gap = None
    if theta_0_eq is not None and beta_a_per_mM is not None and a_molecular_size_m is not None:
        if beta_a_per_mM == 0:
            gap = "beta_a_per_mM is zero -- diffusion coefficient is undefined (division by zero in the near-saturation formula)"
        else:
            diffusion_coefficient = (a_molecular_size_m ** 2) / ((beta_a_per_mM ** 2) * tau_d_best)
    else:
        gap = (
            "diffusion coefficient not computed: requires theta_0_eq, beta_a_per_mM, and "
            "a_molecular_size_m from a SEPARATELY-fitted equilibrium isotherm on this same "
            "system (see adsorption.py's szyszkowski_fit_K/select_isotherm_model) -- only the "
            "raw time constant tau_d is reported without them, not guessed further"
        )

    return DlaFitResult(
        tau_d_s=tau_d_best,
        delta_gamma_eq_mN_m=delta_gamma_eq,
        r_squared=r_squared,
        n_points=n,
        diffusion_coefficient_m2_per_s=diffusion_coefficient,
        gap=gap,
    )


@dataclass
class KlaFitResult:
    tau_k_s: float
    amplitude_mN_m: float
    r_squared: float
    n_points_used: int
    n_points_excluded: int


def fit_kla_regime(surface_age_s: list[float], gamma_mN_m: list[float], gamma0_mN_m: float, delta_gamma_eq_mN_m: float) -> KlaFitResult:
    """Fit the kinetically-limited-adsorption (KLA) longer-time relaxation
    seen in real ionic-surfactant DST data when a real adsorption barrier
    (e.g. electrostatic, for an ionic headgroup) slows adsorption beyond
    what pure diffusion (fit_dla_regime) predicts (Ritacco et al. 2010,
    arXiv:1008.1230v2, eq. 5):

        Delta_gamma_eq - Delta_gamma(t) ~ exp(-t/tau_k)

    (tension approaches equilibrium FROM BELOW as t increases -- Delta_gamma
    grows toward Delta_gamma_eq, so the residual gap Delta_gamma_eq -
    Delta_gamma(t) decays toward zero; this is the physically correct sign
    convention, the opposite of a naive Delta_gamma(t)-Delta_gamma_eq
    reading). This is linear in t after taking the log of the residual gap
    -- fit here by ordinary least squares on ln(Delta_gamma_eq -
    Delta_gamma(t)) vs t, giving tau_k = -1/slope directly.

    delta_gamma_eq_mN_m: the equilibrium tension drop -- normally
    fit_dla_regime's own delta_gamma_eq_mN_m output, so KLA is fit against
    the SAME equilibrium value DLA converges to (real, disclosed physical
    consistency choice, not an independent free parameter here).

    Points where Delta_gamma(t) >= Delta_gamma_eq (already at or past the
    supplied equilibrium value -- no positive residual gap to take the log
    of) are excluded from the fit rather than raising, with the exclusion
    count reported in n_points_excluded so a caller can judge data quality;
    raises only if too few points remain to fit a line (need at least 2).
    """
    if len(surface_age_s) != len(gamma_mN_m):
        raise ValueError("surface_age_s and gamma_mN_m must be the same length")
    if any(t <= 0 for t in surface_age_s):
        raise ValueError("surface_age_s values must all be positive")

    delta_gamma = [gamma0_mN_m - g for g in gamma_mN_m]
    pairs = [(t, dg) for t, dg in zip(surface_age_s, delta_gamma) if dg < delta_gamma_eq_mN_m]
    n_excluded = len(surface_age_s) - len(pairs)
    if len(pairs) < 2:
        raise ValueError(
            f"only {len(pairs)} point(s) have Delta_gamma(t) < delta_gamma_eq_mN_m -- need at "
            "least 2 to fit the KLA exponential relaxation; most of this data may already be at "
            "or past the supplied equilibrium value"
        )

    t = [p[0] for p in pairs]
    ln_residual = [math.log(delta_gamma_eq_mN_m - p[1]) for p in pairs]
    n = len(t)
    mean_t = sum(t) / n
    mean_y = sum(ln_residual) / n
    sxx = sum((ti - mean_t) ** 2 for ti in t)
    if sxx == 0:
        raise ValueError("all included surface_age_s values are identical -- cannot fit a slope")
    sxy = sum((ti - mean_t) * (yi - mean_y) for ti, yi in zip(t, ln_residual))
    slope = sxy / sxx
    intercept = mean_y - slope * mean_t
    if slope >= 0:
        raise ValueError(
            f"fitted slope={slope:.6g} is non-negative -- the residual gap is not decaying "
            "with surface age in this data, no genuine KLA relaxation is present (or the sign "
            "convention/equilibrium value supplied is inconsistent with this data)"
        )

    fitted = [intercept + slope * ti for ti in t]
    ss_res = sum((yi - fi) ** 2 for yi, fi in zip(ln_residual, fitted))
    ss_tot = sum((yi - mean_y) ** 2 for yi in ln_residual)
    r_squared = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 1.0

    tau_k = -1.0 / slope
    amplitude = math.exp(intercept)

    return KlaFitResult(tau_k_s=tau_k, amplitude_mN_m=amplitude, r_squared=r_squared, n_points_used=n, n_points_excluded=n_excluded)


@dataclass
class DstRegimeSelectionResult:
    selected_regime: str  # "dla_only" or "dla_plus_kla"
    dla_fit: DlaFitResult  # ALWAYS the honest standalone fit_dla_regime baseline, never the joint co-fit
    kla_fit: KlaFitResult | None
    dla_only_bic: float
    combined_bic: float | None
    delta_bic: float | None  # dla_only_bic - combined_bic; positive favors the combined model
    interpretation: str
    # Only populated when selected_regime == "dla_plus_kla": the JOINT 4-parameter
    # fit's own co-fit DLA component -- genuinely different from dla_fit's standalone
    # values, since here DLA only has to explain whatever the KLA term doesn't. This
    # is what a caller wants when reporting "the diffusion time" of the ACCEPTED
    # model, as opposed to dla_fit's own (correctly-labeled) rejected-simpler-model baseline.
    combined_tau_d_s: float | None = None
    combined_delta_gamma_eq_mN_m: float | None = None


def select_dst_kinetic_regime(surface_age_s: list[float], gamma_mN_m: list[float], gamma0_mN_m: float) -> DstRegimeSelectionResult:
    """Decide, from real gamma(t) data alone, whether pure diffusion-
    limited adsorption (DLA) already explains the data or whether a
    genuine second, slower kinetically-limited-adsorption (KLA) relaxation
    is also present -- the direct DST-domain analogue of
    adsorption.select_isotherm_model's Langmuir-vs-Frumkin BIC decision
    and curve_analysis.cmc_from_surface_tension_curve's own 2-vs-3-segment
    BIC decision. Same principle: a more flexible model (here, DLA+KLA,
    4 total free parameters: tau_d, tau_k, and two linear amplitudes) will
    always fit equal-or-better than the simpler one (DLA alone, 2 free
    parameters) -- that alone proves nothing. BIC penalizes the extra
    parameters, so the combined model is only selected when the residual
    reduction genuinely earns back more than its parameter penalty.

    Real physical grounding for why this is a genuine judgment call, not
    an arbitrary curve-fitting exercise: Ritacco et al. 2010's real DTAB
    data shows tau_k running 2-3 orders of magnitude LARGER than tau_d at
    higher concentration -- a real, physically detectable second
    timescale attributed to the electrostatic adsorption barrier an ionic
    surfactant's headgroup creates, not present for a nonionic surfactant
    like their Triton X-100 comparison (which the same source confirms IS
    well-described by pure DLA alone even at very short surface ages).

    METHODOLOGY NOTE (2026-09-25): the combined model is fit JOINTLY (2D
    grid over (tau_d, tau_k) with a closed-form 2-parameter linear
    regression for the two amplitudes at each candidate pair), not by
    fitting DLA first and then fitting a residual-gap exponential
    separately. A sequential residual-fitting approach was tried first and
    REJECTED after direct testing on synthetic data with a known-present
    KLA signal: DLA's own flexible 2-parameter saturating shape partially
    absorbs the injected KLA perturbation into its own fitted tau_d/
    Delta_gamma_eq (both curves are monotonically saturating with a single
    principal timescale, so they are not cleanly separable by fitting one
    then subtracting), leaving a residual too small/differently-shaped to
    recover tau_k reliably -- verified this failure mode directly, not
    assumed. The joint fit avoids this degeneracy since both timescales
    compete on equal footing in one simultaneous optimization.

    The returned `dla_fit` is ALWAYS the honest standalone fit_dla_regime
    result (the simpler-model baseline used for the BIC comparison itself),
    in both the dla_only and dla_plus_kla branches -- never the joint
    fit's own co-fit tau_d/Delta_gamma_eq, which can differ substantially
    from the standalone values (the joint fit's DLA component only has to
    explain whatever the KLA term doesn't, so its tau_d is not directly
    comparable to a standalone DLA-only fit's tau_d). BUG FIX (2026-09-27):
    an earlier version of this function substituted the joint fit's own
    (tau_d, Delta_gamma_eq) into the returned `dla_fit` when reporting
    dla_plus_kla, WITH the combined 4-parameter model's r_squared attached
    to it -- misleadingly labeled as if it were a standalone DLA fit's own
    r_squared. Found via a real cross-check on real thesis DST data (GZ-69):
    an independent from-scratch brute-force 2D grid search confirmed the
    regime SELECTION itself (dla_plus_kla, tau_k~0.86s) was always correct,
    but the previously-reported `dla_fit.tau_d_s`/`r_squared` values were
    not a meaningful standalone-DLA characterization. The joint fit's own
    tau_d/Delta_gamma_eq are now exposed separately and honestly via
    `combined_tau_d_s`/`combined_delta_gamma_eq_mN_m` (populated only for
    dla_plus_kla) instead of being conflated with the standalone `dla_fit`.
    """
    dla = fit_dla_regime(surface_age_s, gamma_mN_m, gamma0_mN_m)
    n = len(surface_age_s)
    dla_rss = sum(
        (dla.delta_gamma_eq_mN_m * (1.0 - math.sqrt(dla.tau_d_s / t)) - (gamma0_mN_m - g)) ** 2
        for t, g in zip(surface_age_s, gamma_mN_m)
    )
    bic_dla = _bic(dla_rss, n, n_params=2)

    t = list(surface_age_s)
    delta_gamma = [gamma0_mN_m - g for g in gamma_mN_m]

    def fit_amplitudes(tau_d: float, tau_k: float) -> tuple[float, float, float]:
        """Closed-form 2-parameter (dg_eq, amp) linear regression for FIXED
        (tau_d, tau_k): Delta_gamma(t) = dg_eq*X1(t) + amp*X2(t), where
        X1(t)=1-sqrt(tau_d/t), X2(t)=-exp(-t/tau_k). Returns (dg_eq, amp, rss)."""
        x1 = [1.0 - math.sqrt(tau_d / ti) for ti in t]
        x2 = [-math.exp(-ti / tau_k) for ti in t]
        s11 = sum(a * a for a in x1)
        s22 = sum(b * b for b in x2)
        s12 = sum(a * b for a, b in zip(x1, x2))
        sy1 = sum(a * y for a, y in zip(x1, delta_gamma))
        sy2 = sum(b * y for b, y in zip(x2, delta_gamma))
        det = s11 * s22 - s12 * s12
        if abs(det) < 1e-30:
            return 0.0, 0.0, sum(y * y for y in delta_gamma)
        dg_eq = (sy1 * s22 - sy2 * s12) / det
        amp = (s11 * sy2 - s12 * sy1) / det
        rss = sum((y - dg_eq * a - amp * b) ** 2 for y, a, b in zip(delta_gamma, x1, x2))
        return dg_eq, amp, rss

    t_min, t_max = min(t), max(t)
    ln_lo, ln_hi = math.log(t_min * 1e-3), math.log(t_max * 10.0)
    n_grid = 25
    grid_vals = [math.exp(ln_lo + i * (ln_hi - ln_lo) / (n_grid - 1)) for i in range(n_grid)]

    best = None
    for td in grid_vals:
        for tk in grid_vals:
            if tk <= td:
                continue  # KLA is by definition the SLOWER (longer-timescale) relaxation
            _, _, rss = fit_amplitudes(td, tk)
            if best is None or rss < best[0]:
                best = (rss, td, tk)
    if best is None:
        # No grid point had tau_k > tau_d in this range (degenerate/tiny dataset) --
        # cannot even attempt a combined fit; DLA-only is the only honest answer.
        return DstRegimeSelectionResult(
            selected_regime="dla_only", dla_fit=dla, kla_fit=None,
            dla_only_bic=bic_dla, combined_bic=None, delta_bic=None,
            interpretation="Not enough distinct surface-age range to attempt a combined DLA+KLA fit -- reporting DLA alone.",
        )
    _, td_grid, tk_grid = best

    # Local refine: a second, finer grid around the coarse best point.
    step = (ln_hi - ln_lo) / (n_grid - 1)
    refine_lo_d, refine_hi_d = math.exp(math.log(td_grid) - step), math.exp(math.log(td_grid) + step)
    refine_lo_k, refine_hi_k = math.exp(math.log(tk_grid) - step), math.exp(math.log(tk_grid) + step)
    n_refine = 20
    refine_td = [refine_lo_d + i * (refine_hi_d - refine_lo_d) / (n_refine - 1) for i in range(n_refine)]
    refine_tk = [refine_lo_k + i * (refine_hi_k - refine_lo_k) / (n_refine - 1) for i in range(n_refine)]
    best_refine = None
    for td in refine_td:
        for tk in refine_tk:
            if tk <= td:
                continue
            dg_eq, amp, rss = fit_amplitudes(td, tk)
            if best_refine is None or rss < best_refine[0]:
                best_refine = (rss, td, tk, dg_eq, amp)
    combined_rss, tau_d_combined, tau_k_combined, dg_eq_combined, amp_combined = best_refine
    bic_combined = _bic(combined_rss, n, n_params=4)

    kla_result = KlaFitResult(
        tau_k_s=tau_k_combined, amplitude_mN_m=amp_combined,
        r_squared=1.0 - combined_rss / dla_rss if dla_rss > 0 else 0.0,
        n_points_used=n, n_points_excluded=0,
    )

    delta_bic = bic_dla - bic_combined
    if delta_bic > 0:
        return DstRegimeSelectionResult(
            selected_regime="dla_plus_kla",
            dla_fit=dla,
            kla_fit=kla_result,
            dla_only_bic=bic_dla,
            combined_bic=bic_combined,
            delta_bic=delta_bic,
            combined_tau_d_s=tau_d_combined,
            combined_delta_gamma_eq_mN_m=dg_eq_combined,
            interpretation=(
                f"A genuine second, slower relaxation (tau_k={tau_k_combined:.4g} s) is present beyond "
                f"pure diffusion-limited adsorption (tau_d={tau_d_combined:.4g} s) -- consistent with a "
                "real adsorption barrier (e.g. electrostatic, for an ionic surfactant). "
                f"tau_k/tau_d = {tau_k_combined/tau_d_combined:.3g}."
            ),
        )
    return DstRegimeSelectionResult(
        selected_regime="dla_only",
        dla_fit=dla,
        kla_fit=kla_result,
        dla_only_bic=bic_dla,
        combined_bic=bic_combined,
        delta_bic=delta_bic,
        interpretation=(
            "A combined DLA+KLA fit was attempted, but BIC does not justify its 2 extra "
            "parameters over pure DLA for this data -- likely too few/noisy points to resolve a "
            "genuine second regime, not proof one is absent."
        ),
    )


@dataclass
class RosenHuaFitResult:
    gamma_mesoequilibrium_mN_m: float
    t_star_s: float
    n_exponent: float
    r_squared: float
    n_points: int


def fit_rosen_hua(surface_age_s: list[float], gamma_mN_m: list[float], gamma0_mN_m: float) -> RosenHuaFitResult:
    """Fit the Rosen-Hua empirical dynamic-surface-tension equation (Hua &
    Rosen, J. Colloid Interface Sci. (1988); Rosen, "Surfactants and
    Interfacial Phenomena"), the standard practical curve-shape fit for a
    measured gamma(t) decay, independent of the Ward-Tordai-derived DLA/KLA
    theory above:

        gamma(t) = gamma_m + (gamma0 - gamma_m) / (1 + (t/t_star)^n)

    gamma_m: the DST curve's own apparent "mesoequilibrium" tension (a
    curve-fit asymptote, NOT necessarily the true thermodynamic equilibrium
    value -- a real, standard distinction in this method). t_star:
    characteristic relaxation time. n: dimensionless shape exponent.

    METHODOLOGY NOTE (2026-09-26): fit directly in REAL (gamma) space, not
    via the tempting log-linearized form (ln[(gamma0-gamma_t)/(gamma_t-
    gamma_m)] = n*ln(t/t_star), linear after taking logs, which this
    function's first implementation used for a fast closed-form fit). That
    first version was tested against a real independent answer (a real
    Claude/ChatGPT-generated fit on real GZ-69/71 benchmark data, built
    2026-09-25) and found to disagree sharply: the log-space fit reported
    r_squared~0.72 where a genuine real-space nonlinear fit on the SAME
    data reaches r_squared~0.985 (independently re-verified by direct
    residual calculation, not just trusting the other answer). ROOT CAUSE:
    minimizing SSE in log-space does not minimize the same objective as
    SSE in real gamma space -- points where (gamma_t - gamma_m) is small
    get hugely amplified by the log transform, badly distorting which fit
    the search actually prefers relative to what "r_squared" is documented
    (and expected by any real reader) to mean: goodness-of-fit on the
    actual physical quantity, gamma(t) itself.

    Fixed by fitting directly in real space: gamma(t) = gamma_m*(1-w(t)) +
    gamma0*w(t), where w(t) = 1/(1+(t/t_star)^n) -- for ANY fixed (n,
    t_star), this is linear in the single remaining parameter gamma_m
    (closed-form, ordinary least squares through the origin in the
    transformed (1-w, gamma_t - gamma0*w) coordinates). Searched via a 2D
    grid-then-refine over (n, t_star) -- the same joint-search pattern
    already used in this module for select_dst_kinetic_regime's (tau_d,
    tau_k) search, needed here for the same reason: this 2-parameter
    outer search surface is not guaranteed unimodal either. r_squared is
    now computed directly from real gamma(t) residuals, matching what it
    is documented to mean.
    """
    if len(surface_age_s) != len(gamma_mN_m):
        raise ValueError("surface_age_s and gamma_mN_m must be the same length")
    if len(surface_age_s) < 4:
        raise ValueError("need at least 4 (surface age, surface tension) points to fit gamma_m, t_star, and n meaningfully")
    if any(t <= 0 for t in surface_age_s):
        raise ValueError("surface_age_s values must all be positive")

    t = list(surface_age_s)
    gamma_t = list(gamma_mN_m)
    n_pts = len(t)

    def gamma_m_and_rss(n_exp: float, t_star: float) -> tuple[float, float]:
        w = [1.0 / (1.0 + (ti / t_star) ** n_exp) for ti in t]
        x = [1.0 - wi for wi in w]
        y = [gt - gamma0_mN_m * wi for gt, wi in zip(gamma_t, w)]
        sxx = sum(xi * xi for xi in x)
        if sxx == 0:
            return gamma0_mN_m, sum(yi * yi for yi in y)
        sxy = sum(xi * yi for xi, yi in zip(x, y))
        gamma_m = sxy / sxx
        rss = sum((gt - (gamma_m * xi + gamma0_mN_m * wi)) ** 2 for gt, xi, wi in zip(gamma_t, x, w))
        return gamma_m, rss

    t_min, t_max = min(t), max(t)
    ln_ts_lo, ln_ts_hi = math.log(t_min * 1e-3), math.log(t_max * 10.0)
    n_lo, n_hi = 0.05, 3.0

    n_grid = 30
    grid_ts = [math.exp(ln_ts_lo + i * (ln_ts_hi - ln_ts_lo) / (n_grid - 1)) for i in range(n_grid)]
    grid_n = [n_lo + i * (n_hi - n_lo) / (n_grid - 1) for i in range(n_grid)]

    best = None
    for ts in grid_ts:
        for ne in grid_n:
            try:
                _, rss = gamma_m_and_rss(ne, ts)
            except (OverflowError, ZeroDivisionError):
                continue
            if best is None or rss < best[0]:
                best = (rss, ne, ts)
    if best is None:
        raise ValueError("could not find any valid (n, t_star) candidate for this data")
    _, n_grid_best, ts_grid_best = best

    step_ts = (ln_ts_hi - ln_ts_lo) / (n_grid - 1)
    step_n = (n_hi - n_lo) / (n_grid - 1)
    refine_ts = [math.exp(max(ln_ts_lo, math.log(ts_grid_best) - step_ts) + i * (2 * step_ts) / 19) for i in range(20)]
    refine_n = [max(n_lo, n_grid_best - step_n) + i * (2 * step_n) / 19 for i in range(20)]

    best_refine = None
    for ts in refine_ts:
        for ne in refine_n:
            try:
                gm, rss = gamma_m_and_rss(ne, ts)
            except (OverflowError, ZeroDivisionError):
                continue
            if best_refine is None or rss < best_refine[0]:
                best_refine = (rss, ne, ts, gm)
    rss, n_exp, t_star, gamma_m_best = best_refine

    mean_gamma = sum(gamma_t) / n_pts
    ss_tot = sum((gt - mean_gamma) ** 2 for gt in gamma_t)
    r_squared = 1.0 - (rss / ss_tot) if ss_tot > 0 else 1.0

    return RosenHuaFitResult(
        gamma_mesoequilibrium_mN_m=gamma_m_best,
        t_star_s=t_star,
        n_exponent=n_exp,
        r_squared=r_squared,
        n_points=n_pts,
    )
