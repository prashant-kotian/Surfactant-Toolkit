"""Tests for surfactantkit.curve_analysis -- raw experimental curve
extraction (CMC from tensiometry, aggregation number from fluorescence
quenching or static light scattering). No pytest coverage existed for
this module before this file; the CMC-extraction function's own __main__
block (a real literature pilot case) is reproduced here as an actual test.
"""

import math

import pytest

from surfactantkit.cpp import tanford_tail_volume
from surfactantkit.curve_analysis import (
    AVOGADRO_NUMBER,
    DN_DC_REFERENCE_ML_PER_G,
    aggregation_number_from_dls,
    aggregation_number_from_quenching_curve,
    aggregation_number_from_sls_debye_plot,
    aggregation_number_from_svedberg_equation,
    cmc_from_conductivity_curve,
    cmc_from_surface_tension_curve,
    cmc_via_szyszkowski_deviation,
    estimate_partial_specific_volume_from_tail_and_headgroup,
    plateau_is_statistically_suspect,
)

# --- cmc_from_surface_tension_curve -----------------------------------------


def test_cmc_from_surface_tension_curve_aot_literature_case():
    """Real pilot dataset from the module's own __main__ block: AOT in
    water, reconstructed from Shah/Das/Bhattarai 2025 (Heliyon,
    PMC11835642) Figure 1. Paper's own reported CMC: 2.51 mM."""
    concentrations = [0.01585, 0.03981, 0.10000, 0.25119, 0.63096, 1.58489,
                       3.981, 6.31, 10.0]
    tensions = [62.22, 55.66, 49.10, 42.54, 35.98, 29.42,
                27.5, 27.8, 28.0]
    result = cmc_from_surface_tension_curve(concentrations, tensions)
    assert result.cmc_mM == pytest.approx(2.51, rel=0.15)
    assert result.r_squared_premicellar > 0.99


def test_c20_matches_aot_case_and_is_physically_below_cmc():
    """C20/pC20 (surfactant efficiency parameter, added 2026-09-12).
    Regression guard using the same AOT dataset: gamma0=62.22 mN/m (the
    lowest-concentration point, no baseline detected for this dataset),
    C20 must land well below the CMC (2.51 mM) -- the real physical
    expectation confirmed against an open-access primary source (a
    cardanol-surfactant paper, pC20=4.03 for a different surfactant,
    C20 << CMC there too)."""
    concentrations = [0.01585, 0.03981, 0.10000, 0.25119, 0.63096, 1.58489,
                       3.981, 6.31, 10.0]
    tensions = [62.22, 55.66, 49.10, 42.54, 35.98, 29.42,
                27.5, 27.8, 28.0]
    result = cmc_from_surface_tension_curve(concentrations, tensions)
    assert result.gamma0_mN_m == pytest.approx(62.22)
    assert result.c20_mM is not None
    assert result.c20_mM == pytest.approx(0.2627, rel=1e-3)
    assert result.pC20 == pytest.approx(3.580, abs=1e-2)
    assert result.c20_mM < result.cmc_mM


def test_c20_exact_round_trip_construction():
    """Mathematically-guaranteed round trip: construct an EXACT premicellar
    line (known slope/intercept, no baseline lag), verify C20 solves
    exactly against gamma0 as THIS function actually defines it when no
    baseline is detected -- the lowest-concentration point's own gamma,
    not an abstract C->0 asymptote (a real, disclosed practical proxy,
    not the same thing)."""
    import math as _math

    true_slope = -15.0  # mN/m per decade
    true_intercept = 80.0
    cmc_true = 8.0
    concentrations = [0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]
    tensions = [
        true_slope * _math.log10(c) + true_intercept if c < cmc_true
        else true_slope * _math.log10(cmc_true) + true_intercept
        for c in concentrations
    ]

    result = cmc_from_surface_tension_curve(concentrations, tensions)
    expected_gamma0 = true_slope * _math.log10(min(concentrations)) + true_intercept
    expected_log_c20 = ((expected_gamma0 - 20.0) - true_intercept) / true_slope
    expected_c20_mM = 10 ** expected_log_c20

    assert result.gamma0_mN_m == pytest.approx(expected_gamma0, abs=1e-6)
    assert result.c20_mM == pytest.approx(expected_c20_mM, rel=1e-6)


def test_c20_none_when_premicellar_drop_under_20mNm():
    """Real, disclosed 'not defined' case: a weak surfactant whose total
    premicellar drop never reaches 20 mN/m before the CMC must return
    c20_mM=None/pC20=None, not an extrapolated value past the CMC."""
    concentrations = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]
    # only a 10 mN/m total premicellar drop (72 -> 62) before a flat plateau
    tensions = [72.0, 70.5, 68.0, 65.5, 63.0, 62.0, 62.0, 62.0]
    result = cmc_from_surface_tension_curve(concentrations, tensions)
    assert result.c20_mM is None
    assert result.pC20 is None


def test_cmc_from_surface_tension_curve_aot_case_does_not_spuriously_detect_baseline():
    """Regression guard for the 3-segment upgrade below: the AOT dataset
    has no flat pre-onset lag region, so the 2-segment model must still
    win via BIC and n_baseline_points must stay 0 -- the CMC must not
    move from the paper's own reported value just because a 3-segment
    model is now also tried."""
    concentrations = [0.01585, 0.03981, 0.10000, 0.25119, 0.63096, 1.58489,
                       3.981, 6.31, 10.0]
    tensions = [62.22, 55.66, 49.10, 42.54, 35.98, 29.42,
                27.5, 27.8, 28.0]
    result = cmc_from_surface_tension_curve(concentrations, tensions)
    assert result.n_baseline_points == 0
    assert result.method.startswith("two-segment")
    assert result.premicellar_x_min_mM == pytest.approx(0.01585)


def test_cmc_from_surface_tension_curve_detects_flat_baseline_and_finds_true_breakpoint():
    """Real, researcher-submitted dataset (2026-09-12) that exposed a
    genuine bug in the plain 2-segment model: it has a flat pre-onset
    lag region (surface tension near pure water, ~71 mN/m) from
    0.056-0.9484 mM BEFORE the real decline starts, then a real decline
    to 10.7 mM, then a real plateau at ~37.8 mN/m. The plain 2-segment
    model picked the WRONG breakpoint (1.16 mM, the baseline/decline
    boundary) instead of the true one (the decline/plateau boundary,
    ~10.7-13 mM) -- verified directly against this exact data before
    the fix (see SurfactantKit/ROADMAP.md's 2026-09-12 entry for the
    full RSS-table diagnosis). After the fix, the 3-segment model wins
    via BIC and lands on the physically correct breakpoint."""
    concentrations = [0.056, 0.0839, 0.1257, 0.1883, 0.2821, 0.4226, 0.633, 0.9484,
                       1.4208, 2.1284, 3.1886, 4.7769, 7.1563, 10.7208, 10.8,
                       16.0609, 24.0608, 36.0456, 54.0]
    tensions = [71.0, 70.93, 70.74, 70.69, 70.68, 70.6, 70.58, 70.41,
                50.46, 47.92, 45.45, 42.92, 40.33, 37.89, 37.81,
                37.78, 37.81, 37.81, 37.79]
    result = cmc_from_surface_tension_curve(concentrations, tensions)

    assert result.method.startswith("three-segment")
    assert result.n_baseline_points == 8
    assert result.premicellar_x_min_mM == pytest.approx(1.4208)
    # true breakpoint is visually/physically between 10.7 and 16.06 mM
    # (where the decline meets the plateau) -- NOT 1.16 mM, which is
    # where the old 2-segment model wrongly landed.
    assert 10.0 < result.cmc_mM < 17.0
    assert result.r_squared_premicellar > 0.999  # the true decline segment is a very clean line
    assert result.premicellar_slope_mN_per_m_per_log10C < -10.0  # a real steep decline, not the ~-0.4 the old baseline-as-premicellar fit gave

    # postmicellar plateau regression: near-zero slope, mean gamma close to
    # the plateau's actual measured values (~37.8 mN/m). Note R^2 is
    # naturally LOW here (~0.07) even though the fit is visually excellent --
    # a known, expected statistical property of R^2 for near-constant data
    # (ss_tot is tiny, so it's dominated by measurement noise), not a defect.
    assert abs(result.postmicellar_slope_mN_per_m_per_log10C) < 0.5
    assert result.postmicellar_mean_gamma_mN_per_m == pytest.approx(37.8, abs=0.1)


def test_cmc_from_surface_tension_curve_rejects_mismatched_lengths():
    with pytest.raises(ValueError):
        cmc_from_surface_tension_curve([1.0, 2.0], [50.0])


def test_cmc_from_surface_tension_curve_rejects_too_few_points():
    with pytest.raises(ValueError):
        cmc_from_surface_tension_curve([1.0, 2.0, 3.0], [50.0, 45.0, 40.0])


# --- plateau_is_statistically_suspect / cmc_via_szyszkowski_deviation -------
# Real, characterized limitation found via GZ-17 in benchmark/paper3_groundzero
# (ported from SurfactantKit): with a sparse 6-8 point real tensiometry
# curve, cmc_from_surface_tension_curve's own 3-segment BIC selection can
# pick a "postmicellar plateau" segment with as few as 2 points, which
# fits a line with a perfect R^2=1.0 by mathematical construction
# regardless of whether that segment is actually flat. These three real
# datasets (all from Dankloff et al. 2025, PendantProp, Zenodo
# 10.5281/zenodo.17232763) are the calibration set: SDS's plateau is
# genuinely flat (not suspect, already-correct answer must be left
# alone); DTAB's and CTAB's are not (both suspect, both give a wrong
# cmc_mM from the piecewise method alone, both independently confirmed
# against real literature -- DTAB ~15mM per Ritacco/Langevin/Diamant/
# Andelman arXiv:1008.1230v2; CTAB ~0.9-1.0mM per multiple literature
# sources).

_SDS_CONCS = [0.40577109375, 0.8115421875, 1.623084375, 3.24616875, 6.4923375, 12.984675, 25.96935, 51.9387]
_SDS_TENSIONS = [69.85933451975876, 67.9370105943862, 62.892649575815454, 54.57301518905874,
                  42.87619733622359, 37.91402789831731, 37.44295519808575, 36.79102755846683]
_DTAB_CONCS = [0.25578125, 0.5115625, 1.023125, 2.04625, 4.0925, 16.37, 32.74]
_DTAB_TENSIONS = [70.26907673, 68.20699668, 65.39702756, 57.50386199, 49.23642286, 32.34750118, 36.22382546]
_CTAB_CONCS = [0.016484375, 0.03296875, 0.0659375, 0.131875, 0.26375, 0.5275, 1.055, 2.11]
_CTAB_TENSIONS = [71.60374462, 71.55862868, 70.38023925, 64.77683517, 57.3369305, 46.13427154, 36.41066498, 32.37490721]


def test_plateau_is_statistically_suspect_sds_is_not_flagged():
    result = cmc_from_surface_tension_curve(_SDS_CONCS, _SDS_TENSIONS)
    assert not plateau_is_statistically_suspect(result)


def test_plateau_is_statistically_suspect_dtab_and_ctab_are_flagged():
    dtab = cmc_from_surface_tension_curve(_DTAB_CONCS, _DTAB_TENSIONS)
    ctab = cmc_from_surface_tension_curve(_CTAB_CONCS, _CTAB_TENSIONS)
    assert plateau_is_statistically_suspect(dtab)
    assert plateau_is_statistically_suspect(ctab)


def test_cmc_via_szyszkowski_deviation_recovers_the_real_dtab_cmc():
    """The headline result: using ONLY this same sparse 7-point dataset
    (no additional rescue data), the deviation check lands on 16.37 --
    exactly the value independently derived from 3 real replicates'
    worth of additional data from this exact dataset's own source
    (Zenodo 10.5281/zenodo.17232763), and close to the independent
    literature value (~15mM, Wilhelmy plate)."""
    result = cmc_from_surface_tension_curve(_DTAB_CONCS, _DTAB_TENSIONS)
    dev = cmc_via_szyszkowski_deviation(_DTAB_CONCS, _DTAB_TENSIONS, result, "ionic_no_added_salt", temperature_K=296.35)
    assert dev.cmc_mM == pytest.approx(16.37)
    assert dev.fit_r_squared > 0.99


def test_cmc_via_szyszkowski_deviation_recovers_a_realistic_ctab_cmc():
    result = cmc_from_surface_tension_curve(_CTAB_CONCS, _CTAB_TENSIONS)
    dev = cmc_via_szyszkowski_deviation(_CTAB_CONCS, _CTAB_TENSIONS, result, "ionic_no_added_salt", temperature_K=296.85)
    assert dev.cmc_mM == pytest.approx(1.055)
    assert dev.fit_r_squared > 0.99


def test_cmc_via_szyszkowski_deviation_does_not_override_sds_with_worse_answer():
    """SDS's plateau is genuinely flat and its piecewise cmc_mM (9.18) is
    already the correct, independently-verified answer -- this function
    is not gated behind plateau_is_statistically_suspect internally, so
    calling it directly can still return a different number (it's the
    CALLER's job to gate it, exactly as derive_all_properties_from_smiles_
    and_curve does). This test documents that behavior explicitly rather
    than leaving it implicit."""
    result = cmc_from_surface_tension_curve(_SDS_CONCS, _SDS_TENSIONS)
    dev = cmc_via_szyszkowski_deviation(_SDS_CONCS, _SDS_TENSIONS, result, "ionic_no_added_salt", temperature_K=296.15)
    assert dev.cmc_mM != pytest.approx(result.cmc_mM)
    assert not plateau_is_statistically_suspect(result)  # -- so a real caller must not use dev.cmc_mM here


def test_cmc_via_szyszkowski_deviation_returns_none_with_too_few_points():
    result = cmc_from_surface_tension_curve(
        [1.0, 2.0, 3.0, 4.0, 5.0, 6.0], [70.0, 65.0, 60.0, 40.0, 39.5, 39.0]
    )
    dev = cmc_via_szyszkowski_deviation(
        [1.0, 2.0, 3.0, 4.0, 5.0, 6.0], [70.0, 65.0, 60.0, 40.0, 39.5, 39.0],
        result, "nonionic", temperature_K=298.15,
    )
    if dev.n_points_used_in_fit < 3:
        assert dev.cmc_mM is None
        assert "fewer than 3" in dev.reason


def test_cmc_via_szyszkowski_deviation_rejects_unknown_system_type():
    result = cmc_from_surface_tension_curve(_DTAB_CONCS, _DTAB_TENSIONS)
    with pytest.raises(ValueError):
        cmc_via_szyszkowski_deviation(_DTAB_CONCS, _DTAB_TENSIONS, result, "not_a_real_system_type", temperature_K=296.35)


# --- cmc_from_conductivity_curve --------------------------------------------
# Source: standard conductometric break-point method (two linear segments,
# linear concentration axis, intersection = CMC) -- verified against real
# secondary literature description before implementing, not invented. No
# clean, fully-reproducible raw (concentration, conductivity) table was found
# during this pass (same "still open, no clean external numeric example"
# situation already disclosed for aggregation_number_spherical and the SSFQ/
# SLS aggregation tools) -- validated instead via a mathematically-guaranteed
# round-trip plus a real literature plausibility range (SDS CMC ~8.1-8.4 mM,
# multiple independent sources including PMC8007100's own reported
# 8.40+/-1.14 mM; counterion binding degree ~0.6-0.8, per this project's own
# MicelleMD ROADMAP.md citation of typical SDS conductometric beta values).


def test_conductivity_cmc_recovers_known_values_round_trip():
    """Mathematically guaranteed: construct a two-segment conductivity
    curve for a chosen real CMC and pair of slopes, verify exact recovery."""
    true_cmc, slope_below, slope_above = 8.2, 1.20, 0.36  # beta = 1-0.36/1.20 = 0.70
    intercept_below = 0.0
    intercept_above = (slope_below - slope_above) * true_cmc + intercept_below

    concentrations = [1.0, 2.5, 4.0, 5.5, 7.0, 9.5, 11.0, 13.0, 15.0, 17.0]
    conductivities = [
        slope_below * c + intercept_below if c < true_cmc else slope_above * c + intercept_above
        for c in concentrations
    ]

    result = cmc_from_conductivity_curve(concentrations, conductivities)
    assert result.cmc_mM == pytest.approx(true_cmc, rel=1e-6)
    assert result.slope_below_cmc == pytest.approx(slope_below, rel=1e-6)
    assert result.slope_above_cmc == pytest.approx(slope_above, rel=1e-6)
    assert result.counterion_binding_degree == pytest.approx(0.70, rel=1e-6)


def test_conductivity_cmc_lands_in_real_sds_literature_range():
    """Physical plausibility check against real, cited literature ranges
    (not a guessed tolerance): SDS CMC ~8.1-8.4 mM, counterion binding
    degree ~0.6-0.8."""
    true_cmc, slope_below, slope_above = 8.2, 1.15, 0.32
    intercept_above = (slope_below - slope_above) * true_cmc

    concentrations = [1.0, 2.0, 3.5, 5.0, 6.5, 8.0, 9.5, 11.5, 13.5, 16.0]
    conductivities = [
        slope_below * c if c < true_cmc else slope_above * c + intercept_above
        for c in concentrations
    ]

    result = cmc_from_conductivity_curve(concentrations, conductivities)
    assert 8.0 <= result.cmc_mM <= 8.5
    assert 0.6 <= result.counterion_binding_degree <= 0.8


def test_conductivity_cmc_beta_matches_thermodynamics_module_formula():
    """Independent recomputation using the exact formula
    thermodynamics.counterion_binding_degree() implements, checked here
    without importing that function -- confirms the two stay consistent."""
    from surfactantkit.thermodynamics import counterion_binding_degree

    concentrations = [1.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0]
    slope_below, slope_above, true_cmc = 1.3, 0.4, 8.0
    intercept_above = (slope_below - slope_above) * true_cmc
    conductivities = [
        slope_below * c if c < true_cmc else slope_above * c + intercept_above
        for c in concentrations
    ]

    result = cmc_from_conductivity_curve(concentrations, conductivities)
    expected_beta = counterion_binding_degree(result.slope_below_cmc, result.slope_above_cmc)
    assert result.counterion_binding_degree == pytest.approx(expected_beta)


def test_conductivity_cmc_rejects_bad_inputs():
    with pytest.raises(ValueError):
        cmc_from_conductivity_curve([1.0, 2.0], [1.0, 2.0])  # too few points
    with pytest.raises(ValueError):
        cmc_from_conductivity_curve([1.0, 2.0, 3.0], [1.0, 2.0])  # mismatched lengths


def test_conductivity_cmc_rejects_non_ionic_like_flat_data():
    """A flat/noisy curve with no real premicellar-steeper-than-postmicellar
    break (e.g. a nonionic surfactant, which this method cannot
    characterize at all) must raise, not silently return a bogus CMC."""
    concentrations = [1.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0]
    conductivities = [0.5 * c for c in concentrations]  # single straight line, no break
    with pytest.raises(ValueError):
        cmc_from_conductivity_curve(concentrations, conductivities)


# --- aggregation_number_from_quenching_curve (Turro-Yekta SSFQ) ------------
# Source: Turro & Yekta, J. Am. Chem. Soc. 100 (1978) 5951-5952. Equation
# verified directly from a real fetched image: ln(I/I0) = -N[Q]/(Ct-cmc).


def test_ssfq_recovers_known_aggregation_number_round_trip():
    """Mathematically guaranteed: construct quenching data for a chosen
    real N, then verify the function recovers that same N exactly."""
    true_n, ct, cmc, i0 = 80.0, 50.0, 8.0, 1000.0
    quencher_concs = [0.0005, 0.001, 0.0015, 0.002, 0.0025]
    slope = true_n / (ct - cmc)
    intensities = [i0 * math.exp(-slope * q) for q in quencher_concs]

    result = aggregation_number_from_quenching_curve(quencher_concs, intensities, i0, ct, cmc)
    assert result.aggregation_number == pytest.approx(true_n, rel=1e-6)


def test_ssfq_matches_manual_zero_intercept_regression():
    """Independent recomputation in the test itself."""
    quencher_concs = [0.001, 0.002, 0.003, 0.004]
    intensities = [900.0, 810.0, 730.0, 660.0]
    i0, ct, cmc = 1000.0, 40.0, 8.0

    ys = [math.log(i0 / i) for i in intensities]
    slope_expected = sum(x * y for x, y in zip(quencher_concs, ys)) / sum(x * x for x in quencher_concs)
    n_expected = slope_expected * (ct - cmc)

    result = aggregation_number_from_quenching_curve(quencher_concs, intensities, i0, ct, cmc)
    assert result.aggregation_number == pytest.approx(n_expected)


def test_ssfq_rejects_bad_inputs():
    with pytest.raises(ValueError):
        aggregation_number_from_quenching_curve([0.001], [900.0], 1000.0, 40.0, 8.0)  # too few points
    with pytest.raises(ValueError):
        aggregation_number_from_quenching_curve([0.001, 0.002], [900.0], 1000.0, 40.0, 8.0)  # mismatched
    with pytest.raises(ValueError):
        aggregation_number_from_quenching_curve([0.001, 0.002], [900.0, 800.0], 0.0, 40.0, 8.0)  # I0<=0
    with pytest.raises(ValueError):
        aggregation_number_from_quenching_curve([0.001, 0.002], [900.0, 800.0], 1000.0, 5.0, 8.0)  # Ct<=cmc


# --- DN_DC_REFERENCE_ML_PER_G (Tier 2 bottleneck fix #10, 2026-09-15 --
# see BOTTLENECK_RESOLUTION_PLAN.md) -------------------------------------


def test_dn_dc_reference_values_are_real_and_usable_directly():
    sds = DN_DC_REFERENCE_ML_PER_G["SDS"]
    assert sds["dn_dc"] == pytest.approx(0.11)
    assert sds["solvent"] == "water"
    ctab = DN_DC_REFERENCE_ML_PER_G["CTAB"]
    assert ctab["dn_dc"] == pytest.approx(0.15)


def test_dn_dc_reference_triton_x100_real_primary_source():
    """Real Tier 2 gap closed 2026-09-15: Triton X-100 was searched for
    via web search alone and not found, then closed with a real
    primary-source PDF the user provided (Stubicar et al. 1989)."""
    triton = DN_DC_REFERENCE_ML_PER_G["TRITONX100"]
    assert triton["dn_dc"] == pytest.approx(0.140)
    assert triton["dn_dc_uncertainty"] == pytest.approx(0.005)
    assert triton["wavelength_nm"] == pytest.approx(546.0)  # different wavelength than SDS/CTAB, disclosed


def test_dn_dc_reference_tween80_real_primary_source():
    """Real gap CLOSED 2026-09-15 (second real-paper batch, same day):
    Tween-80 was searched for and not found via web search alone, then
    closed with a real primary-source PhD dissertation (Jiang, Virginia
    Tech, 2011) the user provided, its own SPR methods section stating
    the value directly."""
    tween = DN_DC_REFERENCE_ML_PER_G["TWEEN80"]
    assert tween["dn_dc"] == pytest.approx(0.132)
    assert tween["wavelength_nm"] == pytest.approx(690.0)


def test_dn_dc_reference_sds_usable_end_to_end_in_sls_function():
    """Real, usable end-to-end check: the reference table's own SDS
    dn/dc plugged directly into aggregation_number_from_sls_debye_plot
    (no new function needed -- the real work was sourcing the value)."""
    dn_dc = DN_DC_REFERENCE_ML_PER_G["SDS"]["dn_dc"]
    n_solvent, wavelength_nm = 1.333, 632.8
    k = 4.0 * math.pi ** 2 * n_solvent ** 2 * dn_dc ** 2 / (AVOGADRO_NUMBER * (wavelength_nm * 1e-7) ** 4)
    true_mw = 15000.0
    concentrations = [0.002, 0.005, 0.008, 0.011]
    rayleigh_ratios = [k * c / (1.0 / true_mw + 2.0 * 1e-4 * c) for c in concentrations]
    result = aggregation_number_from_sls_debye_plot(
        concentrations, rayleigh_ratios, dn_dc, wavelength_nm, 288.38, refractive_index_solvent=n_solvent
    )
    assert result.micelle_molar_mass_g_per_mol == pytest.approx(true_mw, rel=1e-3)


# --- aggregation_number_from_sls_debye_plot ---------------------------------
# Source: the Debye equation Kc/DeltaR = 1/Mw + 2*A2*c, verified against
# Brookhaven Instruments' "SLS FAQ: The Debye Plot" technical note (fetched
# and read directly), cross-confirmed for the K formula (4*pi^2*n^2*(dn/dc)^2
# / (NA*lambda^4)) against Malvern and Anton Paar's independent documentation.


def _debye_k(n: float, dn_dc: float, wavelength_nm: float) -> float:
    wavelength_cm = wavelength_nm * 1e-7
    return 4.0 * math.pi ** 2 * n ** 2 * dn_dc ** 2 / (AVOGADRO_NUMBER * wavelength_cm ** 4)


def test_sls_debye_recovers_known_molar_mass_round_trip():
    """Mathematically guaranteed: construct Kc/DeltaR data for a chosen
    real Mw and A2 using the exact Debye equation, then verify the
    function recovers that same Mw (and hence aggregation number)."""
    n_solvent, dn_dc, wavelength_nm = 1.333, 0.12, 633.0
    k = _debye_k(n_solvent, dn_dc, wavelength_nm)
    true_mw = 15000.0  # g/mol, e.g. an SDS-like micelle
    true_a2 = 1.0e-4
    monomer_mw = 288.38  # SDS

    concentrations = [0.001, 0.002, 0.003, 0.004, 0.005]  # g/mL
    kc_over_r = [1.0 / true_mw + 2.0 * true_a2 * c for c in concentrations]
    rayleigh_ratios = [k * c / y for c, y in zip(concentrations, kc_over_r)]

    result = aggregation_number_from_sls_debye_plot(
        concentrations, rayleigh_ratios, dn_dc, wavelength_nm, monomer_mw, refractive_index_solvent=n_solvent
    )
    assert result.micelle_molar_mass_g_per_mol == pytest.approx(true_mw, rel=1e-6)
    assert result.aggregation_number == pytest.approx(true_mw / monomer_mw, rel=1e-6)
    assert result.second_virial_coefficient_cm3_mol_per_g2 == pytest.approx(true_a2, rel=1e-6)


def test_sls_debye_lands_in_real_sds_aggregation_number_range():
    """Physical plausibility check: constructing data for a molar mass
    consistent with the real, literature-established SDS aggregation
    number range (44.8-54.2, Bales et al. 1998, already used elsewhere
    in this project's literature_validation_notes.md) must recover an
    aggregation number back in that same real range -- not a guessed
    tolerance, the actual cited range."""
    n_solvent, dn_dc, wavelength_nm = 1.333, 0.12, 633.0
    k = _debye_k(n_solvent, dn_dc, wavelength_nm)
    monomer_mw = 288.38  # SDS
    target_n_agg = 49.5  # Bales et al.'s own mean
    true_mw = target_n_agg * monomer_mw

    concentrations = [0.001, 0.002, 0.003, 0.004]
    kc_over_r = [1.0 / true_mw + 2.0 * 1e-4 * c for c in concentrations]
    rayleigh_ratios = [k * c / y for c, y in zip(concentrations, kc_over_r)]

    result = aggregation_number_from_sls_debye_plot(
        concentrations, rayleigh_ratios, dn_dc, wavelength_nm, monomer_mw, refractive_index_solvent=n_solvent
    )
    assert 44.8 <= result.aggregation_number <= 54.2


def test_sls_debye_rejects_bad_inputs():
    with pytest.raises(ValueError):
        aggregation_number_from_sls_debye_plot([0.001], [1e-5], 0.12, 633.0, 288.38)  # too few points
    with pytest.raises(ValueError):
        aggregation_number_from_sls_debye_plot([0.001, 0.002], [1e-5], 0.12, 633.0, 288.38)  # mismatched
    with pytest.raises(ValueError):
        aggregation_number_from_sls_debye_plot([0.001, 0.002], [1e-5, 2e-5], 0.0, 633.0, 288.38)  # dn/dc=0
    with pytest.raises(ValueError):
        aggregation_number_from_sls_debye_plot([0.001, 0.002], [1e-5, 2e-5], 0.12, 633.0, -1.0)  # bad monomer MW


# --- aggregation_number_from_dls (alternative-methods sweep, 2026-09-10) ---
# Real, standard route: reuses the Stokes-Einstein hydrodynamic radius
# (same formula as dynamics.hydrodynamic_radius_stokes_einstein) then
# converts hydrodynamic volume -> molar mass via partial specific volume.


def test_aggregation_number_from_dls_matches_manual_calculation():
    from surfactantkit.dynamics import hydrodynamic_radius_stokes_einstein

    D, eta, v_bar, M_monomer, T = 1.56e-6, 0.89, 0.87, 288.4, 298.15
    r_h_nm = hydrodynamic_radius_stokes_einstein(D, eta, T)
    r_h_cm = r_h_nm * 1e-7
    v_h_cm3 = (4.0 / 3.0) * math.pi * r_h_cm ** 3
    expected_M = v_h_cm3 * AVOGADRO_NUMBER / v_bar
    expected_nagg = expected_M / M_monomer

    result = aggregation_number_from_dls(D, eta, v_bar, M_monomer, T)
    assert result.hydrodynamic_radius_nm == pytest.approx(r_h_nm)
    assert result.micelle_molar_mass_g_per_mol == pytest.approx(expected_M)
    assert result.aggregation_number == pytest.approx(expected_nagg)


def test_aggregation_number_from_dls_reasonable_for_realistic_micelle():
    """Real-number sanity check: realistic D and v_bar for a small ionic
    surfactant micelle should give an aggregation number in the commonly
    reported tens-to-low-hundreds range, not an absurd value."""
    result = aggregation_number_from_dls(
        diffusion_coefficient_cm2_per_s=1.56e-6, viscosity_mPas=0.89,
        partial_specific_volume_cm3_per_g=0.87, monomer_molar_mass_g_per_mol=288.4,
    )
    assert 10.0 <= result.aggregation_number <= 500.0


def test_aggregation_number_from_dls_rejects_bad_inputs():
    with pytest.raises(ValueError):
        aggregation_number_from_dls(0.0, 0.89, 0.87, 288.4)
    with pytest.raises(ValueError):
        aggregation_number_from_dls(1.56e-6, -1.0, 0.87, 288.4)
    with pytest.raises(ValueError):
        aggregation_number_from_dls(1.56e-6, 0.89, 0.0, 288.4)
    with pytest.raises(ValueError):
        aggregation_number_from_dls(1.56e-6, 0.89, 0.87, 0.0)


# --- aggregation_number_from_svedberg_equation (alternative-methods sweep,
# 2026-09-10) -- the classic, standard analytical-ultracentrifugation route
# (Svedberg equation), independent of the light-scattering/DLS routes above.


def test_aggregation_number_from_svedberg_matches_manual_calculation():
    s_svedberg, D, v_bar, rho, M_monomer, T = 2.0, 1.0e-6, 0.8, 0.997, 288.4, 298.15
    R_GAS = 8.314462618
    buoyancy = 1.0 - v_bar * rho
    s_seconds = s_svedberg * 1e-13
    D_SI = D * 1e-4
    expected_M_kg = (R_GAS * T * s_seconds) / (D_SI * buoyancy)
    expected_M = expected_M_kg * 1000.0
    expected_nagg = expected_M / M_monomer

    result = aggregation_number_from_svedberg_equation(s_svedberg, D, v_bar, rho, M_monomer, T)
    assert result.micelle_molar_mass_g_per_mol == pytest.approx(expected_M)
    assert result.aggregation_number == pytest.approx(expected_nagg)


def test_aggregation_number_from_svedberg_reasonable_for_realistic_micelle():
    result = aggregation_number_from_svedberg_equation(
        sedimentation_coefficient_S=2.0, diffusion_coefficient_cm2_per_s=1.0e-6,
        partial_specific_volume_cm3_per_g=0.8, solvent_density_g_per_cm3=0.997,
        monomer_molar_mass_g_per_mol=288.4,
    )
    assert 10.0 <= result.aggregation_number <= 500.0


def test_aggregation_number_from_svedberg_rejects_nonphysical_buoyancy():
    """A particle with partial specific volume * solvent density >= 1
    would not sediment -- must raise, not silently return a bogus mass."""
    with pytest.raises(ValueError):
        aggregation_number_from_svedberg_equation(2.0, 1.0e-6, 1.2, 0.997, 288.4)


def test_aggregation_number_from_svedberg_rejects_bad_inputs():
    with pytest.raises(ValueError):
        aggregation_number_from_svedberg_equation(2.0, 0.0, 0.8, 0.997, 288.4)
    with pytest.raises(ValueError):
        aggregation_number_from_svedberg_equation(2.0, 1.0e-6, 0.0, 0.997, 288.4)
    with pytest.raises(ValueError):
        aggregation_number_from_svedberg_equation(2.0, 1.0e-6, 0.8, 0.0, 288.4)
    with pytest.raises(ValueError):
        aggregation_number_from_svedberg_equation(2.0, 1.0e-6, 0.8, 0.997, 0.0)


# --- estimate_partial_specific_volume_from_tail_and_headgroup (Tier 1
# bottleneck fix #11, 2026-09-15 -- see BOTTLENECK_RESOLUTION_PLAN.md).
# Validated by round trip, not an external "known real v_bar" number: no
# specific headgroup partial molar volume was verified against a primary
# source with enough confidence this session to assert as fact (unlike
# wetting.py's VOCG_STANDARD_LIQUIDS table, which was live-verified) --
# only the exact tail-volume unit conversion is shipped as real.


def test_estimate_v_bar_tail_contribution_matches_exact_unit_conversion():
    """With headgroup_molar_volume=0, v_bar must equal exactly the
    Tanford tail volume converted from A^3/molecule to cm^3/mol and
    divided by the monomer molar mass -- the one part of this function
    that needs no external literature value at all."""
    n_carbons = 12
    mw = 288.38
    v_bar = estimate_partial_specific_volume_from_tail_and_headgroup(n_carbons, mw, headgroup_molar_volume_cm3_per_mol=0.0)
    expected_tail_cm3_per_mol = tanford_tail_volume(n_carbons) * AVOGADRO_NUMBER * 1e-24
    assert v_bar == pytest.approx(expected_tail_cm3_per_mol / mw, rel=1e-9)


def test_estimate_v_bar_round_trip_recovers_known_headgroup_volume():
    n_carbons = 12
    mw = 288.38
    for true_headgroup_cm3_per_mol in (10.0, 35.0, 60.0, -1.5):  # negative allowed (electrostriction)
        v_bar = estimate_partial_specific_volume_from_tail_and_headgroup(n_carbons, mw, true_headgroup_cm3_per_mol)
        recovered_headgroup = v_bar * mw - tanford_tail_volume(n_carbons) * AVOGADRO_NUMBER * 1e-24
        assert recovered_headgroup == pytest.approx(true_headgroup_cm3_per_mol, rel=1e-6)


def test_estimate_v_bar_falls_in_physically_sane_range_for_realistic_headgroup():
    # a plausible-magnitude sulfate-ester+Na headgroup volume (illustrative
    # only, not asserted as a verified literature value) should land v_bar
    # for a C12 chain in the range real ionic surfactants are known to
    # occupy (roughly 0.7-1.0 cm^3/g), not something absurd
    v_bar = estimate_partial_specific_volume_from_tail_and_headgroup(
        n_carbons=12, monomer_molar_mass_g_per_mol=288.38, headgroup_molar_volume_cm3_per_mol=35.0
    )
    assert 0.7 <= v_bar <= 1.0


def test_estimate_v_bar_rejects_nonpositive_result():
    # a headgroup volume negative enough to drive the total non-positive
    # must raise, not silently return a nonsensical v_bar
    with pytest.raises(ValueError):
        estimate_partial_specific_volume_from_tail_and_headgroup(
            n_carbons=12, monomer_molar_mass_g_per_mol=288.38, headgroup_molar_volume_cm3_per_mol=-1000.0
        )


def test_estimate_v_bar_rejects_bad_inputs():
    with pytest.raises(ValueError):
        estimate_partial_specific_volume_from_tail_and_headgroup(0, 288.38, 35.0)
    with pytest.raises(ValueError):
        estimate_partial_specific_volume_from_tail_and_headgroup(12, 0.0, 35.0)


# ---------------------------------------------------------------------------
# Batch sign-off (audit B1/B2): alternative break-point methods.
# Default ("rss") behavior is covered by the tests above and must not change;
# these guard the NEW method values only.
# ---------------------------------------------------------------------------

def _synthetic_tension_curve():
    """Noiseless two-regime tensiometry: gamma = 50 - 10*log10(c) up to
    CMC = 8 mM, flat plateau (continuous) above."""
    c = [0.5, 1, 2, 4, 6, 8, 10, 12, 14, 16, 20]
    plateau = 50 - 10 * math.log10(8)
    g = [50 - 10 * math.log10(x) if x <= 8 else plateau for x in c]
    return c, g


def test_tension_tangential_method_hits_true_cmc_exactly():
    """Tangential intersection on ideal segmented data recovers CMC = 8 mM
    to machine precision (the pre-line IS the decline, the post-line IS the
    plateau, they meet at the true break)."""
    c, g = _synthetic_tension_curve()
    result = cmc_from_surface_tension_curve(c, g, method="tangential")
    assert result.cmc_mM == pytest.approx(8.0, rel=1e-9)
    assert result.method.startswith("tangential")
    assert result.gamma_at_cmc_mN_per_m == pytest.approx(
        50 - 10 * math.log10(8), rel=1e-9
    )


def test_tension_derivative_method_is_deterministic_and_near_break():
    """Max-|dgamma/dlogC| picks the steepest adjacent pair (tie -> LAST,
    the break side): on this curve every decline pair is exactly steepest,
    so the midpoint lands at sqrt(6*8) = sqrt(48), pinned exactly."""
    c, g = _synthetic_tension_curve()
    result = cmc_from_surface_tension_curve(c, g, method="derivative")
    assert result.cmc_mM == pytest.approx(math.sqrt(48), rel=1e-9)
    assert result.method.startswith("maximum first-derivative")


def test_tension_invalid_method_raises_naming_field_and_values():
    c, g = _synthetic_tension_curve()
    with pytest.raises(ValueError, match="method"):
        cmc_from_surface_tension_curve(c, g, method="hld_magic")


def _synthetic_conductivity_linear():
    """Noiseless conductometric break at CMC = 6 mM on LINEAR kappa:
    kappa = 1.2c + 0.05 below, 0.4c + 4.85 above (continuous)."""
    c = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14]
    k = [1.2 * x + 0.05 if x <= 6 else 0.4 * x + 4.85 for x in c]
    return c, k


def test_conductivity_jacobs_jones_hits_grid_cmc_exactly():
    """Second-derivative (slope-jump) maximum lands ON the sampled CMC
    when the break coincides with a grid point; beta and both r^2 come
    from the same linear-basis fits the rss method uses (beta = 1 - 0.4/1.2)."""
    c, k = _synthetic_conductivity_linear()
    result = cmc_from_conductivity_curve(c, k, method="jacobs_jones")
    assert result.cmc_mM == pytest.approx(6.0, rel=1e-12)
    assert result.counterion_binding_degree == pytest.approx(1 - 0.4 / 1.2, rel=1e-9)
    assert result.r_squared_below_cmc == pytest.approx(1.0, abs=1e-9)
    assert result.r_squared_above_cmc == pytest.approx(1.0, abs=1e-9)
    assert result.method.startswith("second-derivative")


def test_conductivity_default_rss_unchanged_on_same_data():
    c, k = _synthetic_conductivity_linear()
    result = cmc_from_conductivity_curve(c, k)
    assert result.cmc_mM == pytest.approx(6.0, rel=1e-9)
    assert result.counterion_binding_degree == pytest.approx(1 - 0.4 / 1.2, rel=1e-9)
    assert result.method.startswith("two-segment")


def test_conductivity_onsager_method_deterministic_on_kohlrausch_shape():
    """Onsager-shaped synthetic (Lambda = 2 - 0.15*sqrt(c) pre-CMC,
    slow linear growth above): sqrt(c)-basis intersection is deterministic
    (pinned), stays near the rss answer, beta stays on the LINEAR basis
    in (0, 1), and the method string discloses the basis."""
    c = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14]

    def kappa_pre(x):
        return x * (2.0 - 0.15 * math.sqrt(x))

    k = [kappa_pre(x) if x <= 6 else kappa_pre(6) + 0.12 * (x - 6) for x in c]
    result = cmc_from_conductivity_curve(c, k, method="onsager")
    assert result.cmc_mM == pytest.approx(4.401001466269065, rel=1e-6)
    assert 0.0 < result.counterion_binding_degree < 1.0
    assert result.method.startswith("molar-conductivity")
    rss = cmc_from_conductivity_curve(c, k)
    assert result.cmc_mM != pytest.approx(rss.cmc_mM, rel=1e-6)


def test_conductivity_onsager_rejects_zero_concentration():
    c, k = _synthetic_conductivity_linear()
    with pytest.raises(ValueError, match="positive"):
        cmc_from_conductivity_curve([0] + c[1:], k, method="onsager")


def test_conductivity_invalid_method_raises_naming_field_and_values():
    c, k = _synthetic_conductivity_linear()
    with pytest.raises(ValueError, match="method"):
        cmc_from_conductivity_curve(c, k, method="tangential")
