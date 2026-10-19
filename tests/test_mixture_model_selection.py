"""Tests for surfactantkit.mixture_model_selection.

Same validation strategy already established elsewhere in this project for
purely-numerical solvers with no single external literature source to check
against directly (see test_mixed_micelle.py's own rubingh_beta_regression/
eommm_global_fit tests): mathematically-guaranteed round trips (construct a
series for a KNOWN true parameter, recover it) plus real reduction checks
(a general model must reduce EXACTLY to its own zero-interaction special
case). DTAB/SDS pure CMCs match the same real system already used throughout
test_mixed_micelle.py and Category B of the Paper 3 benchmark.
"""
import math

import pytest

from surfactantkit.mixed_micelle import clint_ideal_cmc
from surfactantkit.mixture_model_selection import (
    rubingh_predict_cmc_mix,
    select_mixture_model,
)

DTAB_PURE_CMC = 14.80  # mM
SDS_PURE_CMC = 8.00  # mM


def _series_for_true_beta(x1_values, beta_true, cmc1, cmc2):
    """Same construction as test_mixed_micelle.py's private helper of the
    same name -- self-contained copy here, matching this project's own
    per-module/per-test-file independence convention."""
    alpha1_series, cmc_mix_series = [], []
    for x1 in x1_values:
        f1 = math.exp(beta_true * (1.0 - x1) ** 2)
        f2 = math.exp(beta_true * x1 ** 2)
        cmc_mix = x1 * f1 * cmc1 + (1.0 - x1) * f2 * cmc2
        alpha1 = x1 * f1 * cmc1 / cmc_mix
        alpha1_series.append(alpha1)
        cmc_mix_series.append(cmc_mix)
    return alpha1_series, cmc_mix_series


# --- rubingh_predict_cmc_mix ---------------------------------------------


def test_rubingh_predict_cmc_mix_recovers_known_constructed_point():
    beta_true = -1.8
    x1_values = [0.3, 0.45, 0.6]
    alpha1_series, cmc_mix_series = _series_for_true_beta(x1_values, beta_true, DTAB_PURE_CMC, SDS_PURE_CMC)
    for alpha1, cmc_mix_true in zip(alpha1_series, cmc_mix_series):
        predicted = rubingh_predict_cmc_mix(alpha1, beta_true, DTAB_PURE_CMC, SDS_PURE_CMC)
        assert predicted == pytest.approx(cmc_mix_true, rel=1e-4)


def test_rubingh_predict_cmc_mix_antagonistic_case():
    beta_true = 1.8
    x1_values = [0.35, 0.5, 0.65]
    alpha1_series, cmc_mix_series = _series_for_true_beta(x1_values, beta_true, DTAB_PURE_CMC, SDS_PURE_CMC)
    for alpha1, cmc_mix_true in zip(alpha1_series, cmc_mix_series):
        predicted = rubingh_predict_cmc_mix(alpha1, beta_true, DTAB_PURE_CMC, SDS_PURE_CMC)
        assert predicted == pytest.approx(cmc_mix_true, rel=1e-4)


def test_rubingh_predict_cmc_mix_reduces_to_clint_at_beta_zero():
    """The required real consistency check: Rubingh's own zero-interaction
    limit (beta=0) must exactly match Clint's independent zero-parameter
    ideal-mixing formula -- both are the same null hypothesis, derived two
    different ways."""
    for alpha1 in [0.2, 0.35, 0.5, 0.65, 0.8]:
        predicted = rubingh_predict_cmc_mix(alpha1, 0.0, DTAB_PURE_CMC, SDS_PURE_CMC)
        expected = clint_ideal_cmc(alpha1, DTAB_PURE_CMC, SDS_PURE_CMC)
        assert predicted == pytest.approx(expected, rel=1e-6)


# --- select_mixture_model -------------------------------------------------


def test_select_mixture_model_validates_length_mismatch():
    with pytest.raises(ValueError, match="same length"):
        select_mixture_model([0.2, 0.5], [1.0, 2.0, 3.0], DTAB_PURE_CMC, SDS_PURE_CMC)


def test_select_mixture_model_validates_min_points():
    with pytest.raises(ValueError, match="at least 3"):
        select_mixture_model([0.2, 0.5], [1.0, 2.0], DTAB_PURE_CMC, SDS_PURE_CMC)


def test_select_mixture_model_picks_clint_for_a_genuinely_ideal_series():
    """A series constructed with beta_true=0 has zero real interaction --
    Clint's 0-parameter null hypothesis should win on BIC over Rubingh's
    1 param and EOMMM's 2, since neither buys a real fit improvement on
    genuinely ideal data."""
    x1_values = [0.25, 0.4, 0.5, 0.6, 0.75]
    alpha1_series, cmc_mix_series = _series_for_true_beta(x1_values, 0.0, DTAB_PURE_CMC, SDS_PURE_CMC)

    result = select_mixture_model(alpha1_series, cmc_mix_series, DTAB_PURE_CMC, SDS_PURE_CMC)
    assert result.selected_model == "clint_ideal"
    clint = next(c for c in result.candidates if c.name == "clint_ideal")
    assert clint.rss_ln_space == pytest.approx(0.0, abs=1e-6)


def test_select_mixture_model_picks_rubingh_for_a_real_regular_solution_series():
    """A series with a real, substantial, constant beta across composition
    should be correctly identified as Rubingh's regime, with beta recovered
    close to the true value -- the actual point of this orchestrator."""
    beta_true = -1.8
    x1_values = [0.25, 0.4, 0.5, 0.6, 0.75]
    alpha1_series, cmc_mix_series = _series_for_true_beta(x1_values, beta_true, DTAB_PURE_CMC, SDS_PURE_CMC)

    result = select_mixture_model(alpha1_series, cmc_mix_series, DTAB_PURE_CMC, SDS_PURE_CMC)
    assert result.selected_model == "rubingh_regular_solution"
    rubingh = next(c for c in result.candidates if c.name == "rubingh_regular_solution")
    assert rubingh.params["beta"] == pytest.approx(beta_true, abs=1e-2)
    # Clint (wrong model for this data) should fit noticeably worse
    clint = next(c for c in result.candidates if c.name == "clint_ideal")
    assert clint.rss_ln_space > rubingh.rss_ln_space


def test_select_mixture_model_reports_rodenas_diagnostic_alongside():
    x1_values = [0.25, 0.4, 0.5, 0.6, 0.75]
    alpha1_series, cmc_mix_series = _series_for_true_beta(x1_values, -1.0, DTAB_PURE_CMC, SDS_PURE_CMC)

    result = select_mixture_model(alpha1_series, cmc_mix_series, DTAB_PURE_CMC, SDS_PURE_CMC)
    assert result.rodenas_diagnostic.n_points == len(x1_values)
    assert all(math.isfinite(x) for x in result.rodenas_diagnostic.x1_rodenas)


def test_select_mixture_model_discloses_excluded_models():
    x1_values = [0.25, 0.4, 0.5, 0.6, 0.75]
    alpha1_series, cmc_mix_series = _series_for_true_beta(x1_values, -1.0, DTAB_PURE_CMC, SDS_PURE_CMC)

    result = select_mixture_model(alpha1_series, cmc_mix_series, DTAB_PURE_CMC, SDS_PURE_CMC)
    joined = " ".join(result.excluded_models)
    assert "rosen" in joined.lower()
    assert "corrin_harkins" in joined.lower() or "corrin-harkins" in joined.lower()
    assert "maeda" in joined.lower()


def test_select_mixture_model_all_three_candidates_present_in_order():
    x1_values = [0.25, 0.4, 0.5, 0.6, 0.75]
    alpha1_series, cmc_mix_series = _series_for_true_beta(x1_values, -1.0, DTAB_PURE_CMC, SDS_PURE_CMC)

    result = select_mixture_model(alpha1_series, cmc_mix_series, DTAB_PURE_CMC, SDS_PURE_CMC)
    names = [c.name for c in result.candidates]
    assert names == ["clint_ideal", "rubingh_regular_solution", "eommm_asymmetric_margules"]


def test_eommm_candidate_ranked_by_bic_over_all_optima_regression_2026_10():
    """Regression guard (2026-10-18): a forward-generated EOMMM-truth series
    (W12=-3.5, W21=-0.9, 0.5% jitter, seed 7 -- the hard50 H-06 bank series)
    has a second near-degenerate feasible basin near (W12=-4.0, W21=-1.2)
    with much worse fit quality.  A single-start eommm_global_fit could land
    there after the 2026-10 helper-consolidation refactor changed scan order,
    which made select_mixture_model wrongly pick Rubingh over the generating
    truth.  select_mixture_model must rank EOMMM candidates by BIC over ALL
    optima so the answer is scan-order invariant."""
    alpha1 = [0.10, 0.25, 0.40, 0.55, 0.70, 0.85]
    cmc_mix = [3.3751, 2.5588, 2.2482, 2.0694, 1.9853, 1.9065]
    sel = select_mixture_model(alpha1, cmc_mix, 1.8, 8.3)
    assert sel.selected_model == "eommm_asymmetric_margules"
    eom = next(c for c in sel.candidates if c.name == "eommm_asymmetric_margules")
    assert eom.params["W12"] == pytest.approx(-3.4, abs=0.15)
    assert eom.params["W21"] == pytest.approx(-0.9, abs=0.25)
    assert eom.bic == pytest.approx(-49.898065, abs=0.05)
