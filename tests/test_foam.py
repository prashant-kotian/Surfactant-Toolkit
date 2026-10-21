"""Tests for surfactantkit.foam -- built 2026-09-25. No real room-temperature
Ross-Miles/Bikerman validation dataset was found during this module's
literature research pass (only elevated-temperature/pressure EOR-context
data, not a fair like-for-like check) -- flagged as real remaining work in
ROADMAP.md, not silently treated as done. Bikerman's index and the decay-
model-selection machinery are tested via exact formula/exact-recovery
checks instead (the same discipline used in test_dynamic_surface_tension.py
for the same reason). The Gibbs-Marangoni elasticity summary is tested via
cross-check against adsorption.py's own already-validated functions.
"""
import math
import pytest

from surfactantkit.foam import (
    bikerman_foaminess_index, ross_miles_summary, select_foam_decay_model,
    foam_stability_elasticity_summary,
)
from surfactantkit.adsorption import gibbs_marangoni_elasticity_langmuir, gibbs_marangoni_elasticity_frumkin


def test_bikerman_foaminess_index_exact_formula():
    assert bikerman_foaminess_index(50.0, 2.0) == pytest.approx(25.0)
    assert bikerman_foaminess_index(0.0, 5.0) == 0.0


def test_bikerman_rejects_invalid_inputs():
    with pytest.raises(ValueError):
        bikerman_foaminess_index(-1.0, 2.0)
    with pytest.raises(ValueError):
        bikerman_foaminess_index(50.0, 0.0)


def test_ross_miles_summary_structures_real_data_shape():
    r = ross_miles_summary([0.0, 1.0, 3.0, 5.0], [180.0, 160.0, 140.0, 130.0])
    assert r.initial_height_mm == 180.0
    assert r.heights_by_time_mm[5.0] == 130.0
    assert r.fraction_retained_at[0.0] == 1.0
    assert r.fraction_retained_at[5.0] == pytest.approx(130.0 / 180.0)


def test_ross_miles_summary_requires_t_zero_reference_point():
    with pytest.raises(ValueError):
        ross_miles_summary([1.0, 3.0, 5.0], [160.0, 140.0, 130.0])


def test_select_foam_decay_model_recovers_exponential():
    h0, k_true = 200.0, 0.15
    times = [0.0, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0]
    heights = [h0 * math.exp(-k_true * t) for t in times]
    r = select_foam_decay_model(times, heights)
    assert r.selected_model == "exponential"
    assert r.exponential_rate_constant_per_min == pytest.approx(k_true, rel=1e-6)
    assert r.exponential_half_life_min == pytest.approx(math.log(2.0) / k_true, rel=1e-6)
    assert r.exponential_r_squared > 0.999


def test_select_foam_decay_model_recovers_logarithmic():
    a_true, b_true = 190.0, 30.0
    times = [0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0]
    heights = [a_true - b_true * math.log(t) for t in times]
    r = select_foam_decay_model(times, heights)
    assert r.selected_model == "logarithmic"
    assert r.logarithmic_slope == pytest.approx(-b_true, rel=1e-6)
    assert r.logarithmic_intercept_mm == pytest.approx(a_true, rel=1e-6)
    assert r.logarithmic_r_squared > 0.999


def test_select_foam_decay_model_rejects_too_few_points():
    with pytest.raises(ValueError):
        select_foam_decay_model([0.0, 1.0], [100.0, 90.0])


def test_gibbs_marangoni_elasticity_langmuir_matches_frumkin_at_a_zero():
    gmax, K, c, T = 3e-6, 500.0, 5.0, 298.15
    e_lang = gibbs_marangoni_elasticity_langmuir(c, gmax, K, "nonionic", T)
    e_frum = gibbs_marangoni_elasticity_frumkin(c, gmax, K, 0.0, "nonionic", T)
    assert e_lang.elasticity_mN_per_m == pytest.approx(e_frum.elasticity_mN_per_m, rel=1e-4)
    assert e_lang.theta == pytest.approx(e_frum.theta, rel=1e-4)


def test_gibbs_marangoni_elasticity_langmuir_increases_with_concentration():
    gmax, K, T = 3e-6, 500.0, 298.15
    e_low = gibbs_marangoni_elasticity_langmuir(1.0, gmax, K, "nonionic", T)
    e_high = gibbs_marangoni_elasticity_langmuir(10.0, gmax, K, "nonionic", T)
    assert e_high.elasticity_mN_per_m > e_low.elasticity_mN_per_m


def test_gibbs_marangoni_elasticity_dilute_limit_approaches_surface_pressure():
    # Real physical sanity check: at low K*C (dilute monolayer), E0 -> Pi
    # (the well-known ideal-2D-gas limiting result) -- verify E0/Pi -> 1
    # as concentration -> 0.
    gmax, K, T, R = 3e-6, 500.0, 298.15, 8.314462618
    c = 1e-6  # very dilute, K*c << 1
    e = gibbs_marangoni_elasticity_langmuir(c, gmax, K, "nonionic", T, R)
    pi = 1.0 * R * T * gmax * math.log(1.0 + K * c) * 1000.0  # n_factor=1 for nonionic, mN/m
    assert e.elasticity_mN_per_m / pi == pytest.approx(1.0, rel=1e-3)
    assert e.warning is None  # dilute, nowhere near the high-theta divergence regime


def test_gibbs_marangoni_elasticity_high_theta_warning_triggers_at_real_literature_concentration():
    # Real, verified-in-practice finding (2026-09-25 literature validation):
    # this project's own already-validated C12E8 isotherm (Ho, Lips & Velikov
    # 2023, PCCP, PMC10132450 -- the same real system behind GZ-08/GZ-29)
    # predicts an implausible ~1985 mN/m elasticity at its own real
    # concentration range's upper end (0.08 mM) -- larger than pure water's
    # own 72 mN/m surface tension. The warning must fire here, not just in
    # a synthetic extreme case.
    gamma_max, K, T = 1.935817e-06, 5173.38, 298.0
    c_mM = 0.08
    r = gibbs_marangoni_elasticity_langmuir(c_mM, gamma_max, K, "nonionic", T)
    assert r.theta > 0.9
    assert r.warning is not None
    assert r.elasticity_mN_per_m > 72.0  # exceeds pure water's own surface tension -- implausible


def test_foam_stability_elasticity_summary_wraps_langmuir_and_frumkin_correctly():
    gmax, K, c, T = 3e-6, 500.0, 5.0, 298.15
    s_lang = foam_stability_elasticity_summary(c, gmax, K, "nonionic", a=0.0, temperature_K=T)
    assert s_lang.isotherm_model_used == "langmuir"
    s_frum = foam_stability_elasticity_summary(c, gmax, K, "nonionic", a=1.0, temperature_K=T)
    assert s_frum.isotherm_model_used == "frumkin"
    expected = gibbs_marangoni_elasticity_langmuir(c, gmax, K, "nonionic", T)
    assert s_lang.gibbs_marangoni_elasticity_mN_per_m == pytest.approx(expected.elasticity_mN_per_m, rel=1e-9)
    assert s_lang.theta == pytest.approx(expected.theta, rel=1e-9)
    assert s_lang.high_theta_warning == expected.warning
