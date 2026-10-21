"""Foamability and foam stability -- distinct real quantities (foamability:
how much foam forms per unit agitation/gas input, a GENERATION metric;
stability: how long formed foam persists, a PERSISTENCE metric), driven by
different physics (generation by adsorption KINETICS -- see
dynamic_surface_tension.py; persistence by film drainage/rupture physics).

Built 2026-09-25, closing a real, previously-total gap (grep confirmed zero
prior foam presence anywhere in this library) after a deep literature
research pass. Real standard measurement methods and their real
quantitative outputs:

  - Ross-Miles test (ASTM D1173): foam height at fixed times after a
    standard pour -- ross_miles_summary structures this real, standard
    data shape.
  - Bikerman's foaminess index (Sigma = foam height / gas superficial
    velocity, a real steady-state measure with units of time, physically
    the average bubble lifetime in the foam column):
    bikerman_foaminess_index.
  - Foam DECAY once formed: real competing functional forms (exponential
    vs logarithmic) exist in the literature for how foam height falls with
    time -- select_foam_decay_model chooses between them via BIC, the
    direct foam-domain analogue of adsorption.select_isotherm_model and
    curve_analysis.cmc_from_surface_tension_curve's own model-selection
    pattern.
  - The headline finding from this module's research pass: Gibbs-Marangoni
    surface elasticity -- a real, peer-reviewed foam-stability PROXY
    (Georgieva, Cagna & Langevin, Soft Matter 5(10) (2009) 2063-2071) --
    is computable directly from an ALREADY-FITTED equilibrium adsorption
    isotherm (adsorption.py's szyszkowski_fit_K/select_isotherm_model),
    with zero new curve-fitting machinery of its own. See
    adsorption.gibbs_marangoni_elasticity_langmuir/_frumkin (implemented
    there, not duplicated here, since it's a direct extension of isotherm
    fitting, not a foam-specific measurement) -- foam_stability_summary in
    this module wraps that function as the real, defensible foam-context
    interpretation it is.

Real, disclosed limit, carried through every function here: foam
GENERATION (foamability itself, as opposed to the stability of foam
already formed) genuinely depends on adsorption RATE -- how fast surfactant
reaches a freshly-created bubble interface during agitation -- which is
exactly what dynamic_surface_tension.py, not this module, addresses. This
module does not build a function that pretends to predict foam generation
from equilibrium data alone; only stability/elasticity is legitimately
derivable that way. A real Ross-Miles or Bikerman MEASUREMENT (this
module's own ross_miles_summary/bikerman_foaminess_index) is the honest
route to a real foamability number.
"""

from __future__ import annotations
import math
from dataclasses import dataclass

from ._stats import bic as _bic_foam

from .adsorption import gibbs_marangoni_elasticity_langmuir, gibbs_marangoni_elasticity_frumkin


def bikerman_foaminess_index(foam_volume_or_height: float, gas_superficial_velocity: float) -> float:
    """Bikerman's foaminess index, Sigma = foam_volume_or_height /
    gas_superficial_velocity -- a real, standard steady-state foam-
    generation measure (Bikerman, J.J., "Foam production: ratio between
    foaminess and rate of foam decay", J. Colloid Interface Sci.). Has
    units of TIME (e.g. cm / (cm/s) = s), physically interpreted as the
    average lifetime of a bubble in the steady-state foam column -- a
    single real formula, no curve-fitting involved.

    foam_volume_or_height: steady-state foam column height/volume reached
    while gas is sparged continuously at a constant rate through the
    solution, in whatever length/volume unit is consistent with
    gas_superficial_velocity's own unit.
    gas_superficial_velocity: gas flow rate divided by column cross-
    sectional area (e.g. cm/s), i.e. the same unit convention as
    foam_volume_or_height divided by time.
    """
    if foam_volume_or_height < 0:
        raise ValueError("foam_volume_or_height must be non-negative")
    if gas_superficial_velocity <= 0:
        raise ValueError("gas_superficial_velocity must be positive")
    return foam_volume_or_height / gas_superficial_velocity


@dataclass
class RossMilesSummary:
    initial_height_mm: float  # height immediately after pour completes (t=0)
    heights_by_time_mm: dict[float, float]  # e.g. {1.0: h1min, 3.0: h3min, 5.0: h5min}
    fraction_retained_at: dict[float, float]  # heights_by_time / initial_height, per time point


def ross_miles_summary(times_min: list[float], heights_mm: list[float]) -> RossMilesSummary:
    """Structure real Ross-Miles test data (ASTM D1173: solution poured
    through a standard column onto a receiving pool; foam height recorded
    immediately after pour (t=0) then at standard intervals, classically
    1/3/5 minutes) into the real, standard reportable form -- no modeling,
    pure data organization, matching this project's own established
    "structure the raw measurement, don't invent a formula where none is
    needed" pattern.

    times_min must include 0.0 (the immediately-after-pour reading) as its
    first (and by convention smallest) entry -- this is the real ASTM
    reference point every fraction_retained_at value is normalized against.
    """
    if len(times_min) != len(heights_mm):
        raise ValueError("times_min and heights_mm must be the same length")
    if len(times_min) < 2:
        raise ValueError("need at least 2 time points (including t=0) for a Ross-Miles summary")
    pairs = sorted(zip(times_min, heights_mm))
    if pairs[0][0] != 0.0:
        raise ValueError("times_min must include 0.0 (the immediately-after-pour reading) as the reference point")
    if any(h < 0 for _, h in pairs):
        raise ValueError("heights_mm values must all be non-negative")
    h0 = pairs[0][1]
    heights_by_time = {t: h for t, h in pairs}
    fraction_retained = {t: (h / h0 if h0 > 0 else 0.0) for t, h in pairs}
    return RossMilesSummary(initial_height_mm=h0, heights_by_time_mm=heights_by_time, fraction_retained_at=fraction_retained)


@dataclass
class FoamDecayModelSelectionResult:
    selected_model: str  # "exponential" or "logarithmic"
    exponential_half_life_min: float
    exponential_rate_constant_per_min: float
    exponential_r_squared: float
    exponential_bic: float
    logarithmic_slope: float
    logarithmic_intercept_mm: float
    logarithmic_r_squared: float
    logarithmic_bic: float
    delta_bic: float  # exponential_bic - logarithmic_bic; positive favors logarithmic


def select_foam_decay_model(times_min: list[float], heights_mm: list[float]) -> FoamDecayModelSelectionResult:
    """Choose between the two real, competing functional forms reported in
    the literature for foam height decay once foam has formed, via BIC --
    the direct foam-domain analogue of adsorption.select_isotherm_model
    and curve_analysis.cmc_from_surface_tension_curve's own model-
    selection pattern (a more flexible model always fits equal-or-better;
    BIC penalizes unearned extra parameters).

    Exponential (first-order decay kinetics, the simpler/more commonly
    assumed form): h(t) = h0*exp(-k*t), giving a real half-life
    t_1/2 = ln(2)/k. Fit here via a log-linear OLS regression (ln(h) vs t
    is exactly linear for this model -- exact closed form, not an
    approximation).

    Logarithmic (real, reported for longer-time foam persistence in some
    systems -- a genuinely distinct functional shape, not just a
    reparameterization): h(t) = a - b*ln(t) (t=0 excluded from this fit,
    the logarithm is singular there). Also a closed-form OLS regression
    (h vs ln(t) is exactly linear for this model).

    Both models have exactly 2 free parameters, so BIC's parameter penalty
    term is identical for both -- the comparison reduces to whichever
    achieves lower total RSS on the SAME data, but computed via the same
    formal BIC machinery as the rest of this library for consistency, and
    to make the comparison mechanism uniform for a future 3rd candidate
    model with a different parameter count.

    times_min, heights_mm: real measured foam-height-vs-time data (e.g.
    from ross_miles_summary's own heights_by_time_mm, or a denser decay
    curve if available). Points at t=0 are excluded from the logarithmic
    fit only (singular there); the exponential fit uses all points with
    h>0.
    """
    if len(times_min) != len(heights_mm):
        raise ValueError("times_min and heights_mm must be the same length")
    if len(times_min) < 3:
        raise ValueError("need at least 3 (time, height) points to compare two 2-parameter decay models meaningfully")
    if any(h < 0 for h in heights_mm):
        raise ValueError("heights_mm values must all be non-negative")
    if any(t < 0 for t in times_min):
        raise ValueError("times_min values must all be non-negative")

    exp_pairs = [(t, h) for t, h in zip(times_min, heights_mm) if h > 0]
    if len(exp_pairs) < 3:
        raise ValueError("need at least 3 points with positive height to fit the exponential decay model")
    exp_t = [p[0] for p in exp_pairs]
    exp_ln_h = [math.log(p[1]) for p in exp_pairs]
    n_exp = len(exp_t)
    mean_t = sum(exp_t) / n_exp
    mean_ln_h = sum(exp_ln_h) / n_exp
    sxx = sum((ti - mean_t) ** 2 for ti in exp_t)
    if sxx == 0:
        raise ValueError("all times_min values are identical -- cannot fit a decay rate")
    sxy = sum((ti - mean_t) * (yi - mean_ln_h) for ti, yi in zip(exp_t, exp_ln_h))
    k_neg = sxy / sxx  # = -k
    ln_h0_fit = mean_ln_h - k_neg * mean_t
    k = -k_neg
    fitted_exp = [ln_h0_fit + k_neg * ti for ti in exp_t]
    ss_res_exp = sum((yi - fi) ** 2 for yi, fi in zip(exp_ln_h, fitted_exp))
    ss_tot_exp = sum((yi - mean_ln_h) ** 2 for yi in exp_ln_h)
    r2_exp = 1.0 - (ss_res_exp / ss_tot_exp) if ss_tot_exp > 0 else 1.0
    half_life = math.log(2.0) / k if k > 0 else float("inf")
    bic_exp = _bic_foam(ss_res_exp, n_exp, n_params=2)

    log_pairs = [(t, h) for t, h in zip(times_min, heights_mm) if t > 0]
    if len(log_pairs) < 3:
        raise ValueError("need at least 3 points with t>0 to fit the logarithmic decay model")
    log_t = [math.log(p[0]) for p in log_pairs]
    log_h = [p[1] for p in log_pairs]
    n_log = len(log_t)
    mean_lt = sum(log_t) / n_log
    mean_h = sum(log_h) / n_log
    sxx2 = sum((lt - mean_lt) ** 2 for lt in log_t)
    if sxx2 == 0:
        raise ValueError("all times_min values (excluding t=0) are identical -- cannot fit the logarithmic model")
    sxy2 = sum((lt - mean_lt) * (h - mean_h) for lt, h in zip(log_t, log_h))
    slope = sxy2 / sxx2  # = -b
    intercept = mean_h - slope * mean_lt  # = a
    fitted_log = [intercept + slope * lt for lt in log_t]
    ss_res_log = sum((h - f) ** 2 for h, f in zip(log_h, fitted_log))
    ss_tot_log = sum((h - mean_h) ** 2 for h in log_h)
    r2_log = 1.0 - (ss_res_log / ss_tot_log) if ss_tot_log > 0 else 1.0
    bic_log = _bic_foam(ss_res_log, n_log, n_params=2)

    delta_bic = bic_exp - bic_log
    selected = "logarithmic" if bic_log < bic_exp else "exponential"

    return FoamDecayModelSelectionResult(
        selected_model=selected,
        exponential_half_life_min=half_life,
        exponential_rate_constant_per_min=k,
        exponential_r_squared=r2_exp,
        exponential_bic=bic_exp,
        logarithmic_slope=slope,
        logarithmic_intercept_mm=intercept,
        logarithmic_r_squared=r2_log,
        logarithmic_bic=bic_log,
        delta_bic=delta_bic,
    )


@dataclass
class FoamStabilityElasticitySummary:
    gibbs_marangoni_elasticity_mN_per_m: float
    isotherm_model_used: str  # "langmuir" or "frumkin"
    theta: float  # Langmuir/Frumkin surface coverage at this concentration, 0-1
    high_theta_warning: str | None  # non-None above the real, literature-found divergence threshold
    interpretation: str


def foam_stability_elasticity_summary(
    concentration: float,
    gamma_max_mol_per_m2: float,
    K: float,
    system_type: str,
    a: float = 0.0,
    temperature_K: float = 298.15,
    r_gas: float = 8.314462618,
) -> FoamStabilityElasticitySummary:
    """Real, peer-reviewed foam-stability PROXY computed directly from an
    already-fitted equilibrium adsorption isotherm, with zero new curve-
    fitting of its own -- the headline finding of this module's own
    literature research pass. Wraps adsorption.py's
    gibbs_marangoni_elasticity_langmuir/_frumkin (a=0 routes to the
    closed-form Langmuir version; a!=0 routes to the numerical Frumkin
    version) with the real citation and interpretation this quantity
    carries specifically in a FOAM context: Georgieva, Cagna & Langevin,
    "Link between surface elasticity and foam stability", Soft Matter
    5(10) (2009) 2063-2071, doi:10.1039/b822568k -- measured real surface
    dilational elastic moduli for nonionic surfactants (C12E6, C12G2,
    Pluronics) and found foam lifetime correlates with the adsorbed film's
    elastic modulus: a HIGHER Gibbs elasticity indicates a real, physically
    grounded expectation of GREATER foam stability (the Gibbs-Marangoni
    mechanism: a locally over-stretched film patch has locally higher
    surface tension, which physically pulls surrounding surface material
    back to heal the thin spot -- more elastic monolayers heal more
    strongly).

    REAL, DISCLOSED LIMIT, same as gibbs_marangoni_elasticity_langmuir's
    own docstring: this is a single-INTERFACE (monolayer) elasticity, not
    a full foam THIN-FILM elasticity model (which would additionally
    depend on film thickness via a real surfactant-reservoir-depletion
    effect this function does not capture -- no primary-sourced closed
    form for that additional dependence could be verified during this
    module's research pass, so it was deliberately not built rather than
    hardcoding an unverified formula). Use this as the real, defensible
    PROXY/ranking tool it is (e.g. comparing relative foam-stability
    expectation across several concentrations or surfactants from their
    already-fitted isotherms), not as an absolute foam-lifetime predictor.

    REAL, VERIFIED-IN-PRACTICE LIMIT, added 2026-09-25 via literature
    validation against this project's own already-validated C12E8 isotherm
    (the same real system behind GZ-08/GZ-29): the elasticity formula's
    own known unbounded-growth limitation was found to trigger at REAL,
    in-range concentrations, not just a theoretical extreme -- see
    high_theta_warning on the returned result, non-None whenever theta
    exceeds adsorption.py's GIBBS_ELASTICITY_HIGH_THETA_WARNING_THRESHOLD
    (0.9). Check this field before trusting a high elasticity value as
    physically meaningful.

    concentration, gamma_max_mol_per_m2, K, system_type, temperature_K,
    r_gas: same meaning as adsorption.py's isotherm functions.
    a: Frumkin lateral-interaction parameter (0.0 = plain Langmuir, the
    default and the common case when only szyszkowski_fit_K/select_
    isotherm_model selected Langmuir).
    """
    if a == 0.0:
        result = gibbs_marangoni_elasticity_langmuir(concentration, gamma_max_mol_per_m2, K, system_type, temperature_K, r_gas)
        model_used = "langmuir"
    else:
        result = gibbs_marangoni_elasticity_frumkin(concentration, gamma_max_mol_per_m2, K, a, system_type, temperature_K, r_gas)
        model_used = "frumkin"

    interpretation = (
        f"Gibbs-Marangoni elasticity = {result.elasticity_mN_per_m:.4g} mN/m at this concentration "
        f"({model_used} isotherm, theta={result.theta:.3f}). A real, peer-reviewed foam-stability proxy "
        "(Georgieva, Cagna & Langevin, Soft Matter 2009): higher elasticity is associated with greater "
        "film self-healing and therefore greater expected foam stability, all else equal. This is a "
        "single-interface monolayer proxy, not an absolute foam-lifetime prediction -- use it "
        "comparatively (across concentrations or surfactants), and combine with a real measured foam "
        "decay (select_foam_decay_model) or Ross-Miles/Bikerman data when an absolute number is needed."
    )
    if result.warning:
        interpretation += " CAVEAT: " + result.warning

    return FoamStabilityElasticitySummary(
        gibbs_marangoni_elasticity_mN_per_m=result.elasticity_mN_per_m,
        isotherm_model_used=model_used,
        theta=result.theta,
        high_theta_warning=result.warning,
        interpretation=interpretation,
    )
