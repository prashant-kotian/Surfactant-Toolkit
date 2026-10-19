"""Model-selection orchestrator for binary surfactant mixtures: given a raw
(alpha1, cmc_mix) composition series plus the two pure-component CMCs, decide
which of the competing mixture theories in mixed_micelle.py actually fits this
system, rather than requiring the caller to already know which one to use.

This is the mixture-side analogue of adsorption.py's select_isotherm_model
(Langmuir vs. Frumkin via BIC) -- same idea (real competing-model fit +
information-criterion comparison, not a coin flip), applied to the mixture-
CMC-model family, and the same real gap this project's own
TOOLKIT_CAPABILITY_AUDIT.md flagged: mixed_micelle.py already has 9 real,
tested competing theories, but nothing ran them together and concluded which
one fits until this module.

Not every one of the 9 theories in mixed_micelle.py is a competing PARAMETRIC
predictor of cmc_mix(alpha1) from the same raw input, so not every one enters
the model-selection race here -- each exclusion is a real, structural reason,
disclosed explicitly rather than silently omitted:

  - Rosen's monolayer extension needs surface-tension-derived concentrations
    (c1_sigma/c2_sigma/c_mix_sigma at a reference tension), not CMC data --
    a genuinely different experimental input, not obtainable from a bare
    (alpha1_series, cmc_mix_series, cmc1, cmc2) call.
  - Corrin-Harkins needs a counterion-concentration series (a salt-addition
    experiment), not a composition series -- different experiment entirely.
  - Maeda's free energy of micellization is not an independent cmc_mix
    predictor at all -- it reuses Rubingh's own already-solved x1/beta as an
    input and returns a thermodynamic quantity (deltaG), not a competing
    cmc_mix(alpha1) curve. It doesn't compete with Rubingh, it extends it.
  - Rodenas' method is real and included, but not folded into the SAME BIC
    race as Clint/Rubingh/EOMMM: it is a model-INDEPENDENT decomposition of
    the given curve (a numerical derivative), not a parametric model that
    predicts cmc_mix from a fitted parameter -- there is no over/underfitting
    number of parameters to penalize, and no "residual against the data" in
    the same sense (it always exactly reproduces whatever local slope the
    data has). Reported alongside the BIC race as an independent cross-check,
    not a competing entry in it.

What actually competes, on a common footing (predicted cmc_mix(alpha1) vs.
actual, in ln-space, BIC-penalized by real parameter count):
  - Clint ideal mixing (0 parameters -- the null/no-interaction hypothesis)
  - Rubingh regular-solution theory, regressed beta (1 parameter)
  - EOMMM asymmetric Margules global fit (2 parameters, W12/W21)
"""
from __future__ import annotations

import itertools
import math
from dataclasses import dataclass, field

from ._stats import bic as _bic
from .mixed_micelle import (
    EommmGlobalFitResult,
    RodenasSeriesResult,
    RubinghBetaRegressionResult,
    clint_ideal_cmc,
    eommm_global_fit,
    eommm_global_fit_all_optima,
    rodenas_x1_series,
    rubingh_beta_regression,
)


def rubingh_predict_cmc_mix(
    alpha1: float,
    beta: float,
    cmc1: float,
    cmc2: float,
    grid_points: int = 4000,
    bisection_iters: int = 100,
) -> float | None:
    """Forward-predict the mixed CMC at a given bulk mole fraction alpha1,
    for a KNOWN (already-fit) Rubingh interaction parameter beta -- the
    missing forward direction this module needed and mixed_micelle.py did
    not have: solve_rubingh_x/rubingh_beta only go from a MEASURED cmc_mix
    to the implied (x1, beta), never the other way.

    Derivation (real, from mixed_micelle.py's own already-verified formulas,
    not a new assumption): at the true root x1 of the original Rubingh
    residual, both components' own mass-balance relations give the SAME
    cmc_mix:
        cmc_mix = x1 * cmc1 * f1(x1, beta) / alpha1                  (A)
        cmc_mix = (1 - x1) * cmc2 * f2(x1, beta) / (1 - alpha1)       (B)
    (f1, f2 from activity_coefficients -- see rubingh_beta's own definition,
    beta = term1/(1-x1)^2, term1 = log((alpha1*cmc_mix)/(x1*cmc1)), which is
    exactly (A) rearranged; (B) follows symmetrically from the residual's own
    x1^2*term1 = (1-x1)^2*term2 identity at the root). Given beta (rather
    than cmc_mix) as the known quantity, this leaves ONE equation in the one
    unknown x1: where (A) and (B) agree. Solved by the same grid-scan +
    bisection root-finding pattern already used throughout this module
    (solve_rubingh_x), just with the roles of cmc_mix and beta swapped.

    Verified against two independent existing formulas, not just the
    derivation above (see tests/test_mixture_model_selection.py):
      - at beta=0 (no interaction), the solved x1 reduces EXACTLY to
        motomura_ideal_composition's own zero-interaction formula, and the
        predicted cmc_mix reduces EXACTLY to clint_ideal_cmc -- both are the
        same zero-parameter null hypothesis approached two different ways,
        and they must agree.
      - round-trip: predicting cmc_mix at a given (alpha1, beta) and then
        re-solving beta from that predicted cmc_mix via the EXISTING
        solve_rubingh_x/rubingh_beta recovers the original beta.

    Returns None if no root is found in (0, 1) -- mirrors solve_rubingh_x's
    own contract for the same underlying reason (e.g. a beta/alpha1
    combination with no physically consistent composition).
    """
    def h(x1: float) -> float:
        f1 = math.exp(beta * (1.0 - x1) ** 2)
        f2 = math.exp(beta * x1 ** 2)
        return (x1 * cmc1 * f1) / alpha1 - ((1.0 - x1) * cmc2 * f2) / (1.0 - alpha1)

    eps = 1e-6
    xs = [eps + i * (1.0 - 2.0 * eps) / grid_points for i in range(grid_points + 1)]
    values: list[tuple[float, float]] = []
    for x in xs:
        try:
            fx = h(x)
            if math.isfinite(fx):
                values.append((x, fx))
        except (ValueError, OverflowError, ZeroDivisionError):
            continue
    if not values:
        return None

    root_x1 = None
    for (x_lo, f_lo), (x_hi, f_hi) in itertools.pairwise(values):
        if f_lo == 0:
            root_x1 = x_lo
            break
        if f_lo * f_hi < 0:
            lo, hi, flo = x_lo, x_hi, f_lo
            for _ in range(bisection_iters):
                mid = (lo + hi) / 2.0
                fmid = h(mid)
                if abs(fmid) < 1e-14:
                    root_x1 = mid
                    break
                if flo * fmid <= 0:
                    hi = mid
                else:
                    lo, flo = mid, fmid
            if root_x1 is None:
                root_x1 = (lo + hi) / 2.0
            break

    if root_x1 is None:
        best_x, best_val = min(values, key=lambda pair: abs(pair[1]))
        if abs(best_val) >= 1e-6:
            return None
        root_x1 = best_x

    f1 = math.exp(beta * (1.0 - root_x1) ** 2)
    return (root_x1 * cmc1 * f1) / alpha1


@dataclass
class MixtureModelCandidate:
    name: str
    n_params: int
    params: dict  # e.g. {"beta": ...} or {"W12": ..., "W21": ...}
    predicted_cmc_mix: list[float]  # same order as the input series
    rss_ln_space: float  # sum of squared residuals in ln(cmc_mix) space
    bic: float


@dataclass
class MixtureModelSelectionResult:
    alpha1_used: list[float]
    cmc_mix_used: list[float]
    candidates: list[MixtureModelCandidate]  # Clint, Rubingh, EOMMM -- in that order
    selected_model: str  # the candidate name with the lowest BIC
    rodenas_diagnostic: RodenasSeriesResult  # independent cross-check, NOT part of the BIC race
    excluded_models: list[str] = field(default_factory=list)  # Rosen/Corrin-Harkins/Maeda + why


def select_mixture_model(
    alpha1_series: list[float],
    cmc_mix_series: list[float],
    cmc1: float,
    cmc2: float,
) -> MixtureModelSelectionResult:
    """Given a raw binary-surfactant composition series (no model named),
    fit Clint, Rubingh, and EOMMM to it, compare them by BIC, and report
    which one actually fits -- the mixture-side analogue of
    adsorption.select_isotherm_model. See this module's own docstring for
    exactly why Rosen/Corrin-Harkins/Maeda are not entered into this race,
    and why Rodenas is reported separately rather than competing in it.

    Needs at least 3 composition points (EOMMM's own real minimum,
    inherited here since it is the strictest of the three competing
    models -- Clint and Rubingh alone would tolerate 2, but a 2-point
    fit cannot be genuinely compared against a 2-parameter model on equal
    footing).
    """
    if len(alpha1_series) != len(cmc_mix_series):
        raise ValueError("alpha1_series and cmc_mix_series must be the same length")
    if len(alpha1_series) < 3:
        raise ValueError("need at least 3 composition points for a genuine model comparison")

    n = len(alpha1_series)
    ln_actual = [math.log(c) for c in cmc_mix_series]

    def _rss(predicted: list[float | None]) -> float | None:
        if any(p is None for p in predicted):
            return None
        return sum((la - math.log(p)) ** 2 for la, p in zip(ln_actual, predicted))

    # --- Clint: 0 parameters, the null/no-interaction hypothesis ---
    clint_pred = [clint_ideal_cmc(a, cmc1, cmc2) for a in alpha1_series]
    clint_rss = _rss(clint_pred)
    clint_candidate = MixtureModelCandidate(
        name="clint_ideal", n_params=0, params={},
        predicted_cmc_mix=clint_pred, rss_ln_space=clint_rss,
        bic=_bic(clint_rss, n, 0),
    )

    # --- Rubingh: 1 parameter (beta, regressed across the series) ---
    rubingh_reg: RubinghBetaRegressionResult = rubingh_beta_regression(alpha1_series, cmc_mix_series, cmc1, cmc2)
    rubingh_pred = [rubingh_predict_cmc_mix(a, rubingh_reg.beta_mean, cmc1, cmc2) for a in alpha1_series]
    rubingh_rss = _rss(rubingh_pred)
    rubingh_candidate = MixtureModelCandidate(
        name="rubingh_regular_solution", n_params=1,
        params={"beta": rubingh_reg.beta_mean, "beta_std": rubingh_reg.beta_std},
        predicted_cmc_mix=rubingh_pred,
        rss_ln_space=rubingh_rss if rubingh_rss is not None else float("inf"),
        bic=_bic(rubingh_rss, n, 1) if rubingh_rss is not None else float("inf"),
    )

    # --- EOMMM: 2 parameters (W12, W21), already returns its own implied cmc_mix per point.
    # EOMMM's grid fit is documented to have NO uniqueness guarantee (see
    # eommm_global_fit's own docstring, point 3): near-degenerate feasible optima
    # with very different fit quality can coexist, and a single-start search
    # silently returns whichever basin its scan order hits. That made this
    # selection vulnerable to a 2026-10 refactor (helper consolidation) changing
    # scan order and landing on a worse-but-still-feasible optimum. Fix, per the
    # benchmark's own doctrine: rank the EOMMM candidate by BIC over ALL optima
    # via eommm_global_fit_all_optima, so the selection is scan-order invariant.
    multi = eommm_global_fit_all_optima(alpha1_series, cmc_mix_series, cmc1, cmc2)
    eommm_fits = list(multi.optima)
    _single = eommm_global_fit(alpha1_series, cmc_mix_series, cmc1, cmc2)
    # union of entry points, deduplicated: basin detection can merge two
    # optima separated by a near-zero-infeasibility plateau into one basin
    # and refine only one of them -- the single-start fit's own landing
    # point must stay in the pool for scan-order invariance.
    if all(abs(_single.W12 - o.W12) > 1e-6 or abs(_single.W21 - o.W21) > 1e-6
           for o in eommm_fits):
        eommm_fits.append(_single)
    if not eommm_fits:
        eommm_fits = [_single]
    eommm_fit: EommmGlobalFitResult | None = None
    eommm_rss: float | None = None
    eommm_bic = float("inf")
    for _fit in eommm_fits:
        _cand_rss = _rss(_fit.cmc_var_values)
        if _cand_rss is None:
            continue
        _cand_bic = _bic(_cand_rss, n, 2)
        if _cand_bic < eommm_bic:
            eommm_fit, eommm_rss, eommm_bic = _fit, _cand_rss, _cand_bic
    eommm_candidate = MixtureModelCandidate(
        name="eommm_asymmetric_margules", n_params=2,
        params={"W12": eommm_fit.W12, "W21": eommm_fit.W21} if eommm_fit is not None else {},
        predicted_cmc_mix=eommm_fit.cmc_var_values if eommm_fit is not None else [],
        rss_ln_space=eommm_rss if eommm_rss is not None else float("inf"),
        bic=eommm_bic,
    )

    candidates = [clint_candidate, rubingh_candidate, eommm_candidate]
    selected = min(candidates, key=lambda c: c.bic)

    rodenas = rodenas_x1_series(alpha1_series, cmc_mix_series)

    return MixtureModelSelectionResult(
        alpha1_used=list(alpha1_series), cmc_mix_used=list(cmc_mix_series),
        candidates=candidates, selected_model=selected.name, rodenas_diagnostic=rodenas,
        excluded_models=[
            ("rosen_monolayer: needs surface-tension-derived concentrations (c1_sigma/c2_sigma/c_mix_sigma "
            "at a reference tension), not CMC data -- a different experiment, not derivable from this input."),
            "corrin_harkins: needs a counterion-concentration (salt-addition) series, not a composition series.",
            ("maeda: not an independent cmc_mix predictor -- it reuses Rubingh's own x1/beta to compute a "
            "free energy, it does not compete with Rubingh for which model fits cmc_mix(alpha1)."),
        ],
    )
