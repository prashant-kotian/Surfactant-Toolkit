"""Tests for surfactantkit.thermodynamics."""

import math

import pytest

from surfactantkit.thermodynamics import (
    COUNTERION_BINDING_DEGREE_REFERENCE,
    cmc_to_mole_fraction,
    counterion_binding_degree,
    critical_micellization_degree,
    entropy_micellization,
    gibbs_free_energy_micellization,
    mass_action_free_energy_of_micellization,
    vant_hoff_enthalpy,
    vant_hoff_multi_point_fit,
)

R_GAS = 8.314462618


def test_cmc_to_mole_fraction_dilute_approximation():
    # for a dilute CMC, X_cmc should be very close to CMC/55.5
    x = cmc_to_mole_fraction(0.008)  # 8 mM
    assert x == pytest.approx(0.008 / 55.5, rel=1e-3)


def test_cmc_to_mole_fraction_rejects_nonpositive():
    with pytest.raises(ValueError):
        cmc_to_mole_fraction(0.0)


def test_counterion_binding_degree_basic():
    # slope drops by half above the CMC -> beta = 0.5
    assert counterion_binding_degree(slope_below_cmc=2.0, slope_above_cmc=1.0) == pytest.approx(0.5)


def test_counterion_binding_degree_rejects_wrong_slope_order():
    with pytest.raises(ValueError):
        counterion_binding_degree(slope_below_cmc=1.0, slope_above_cmc=2.0)


def test_counterion_binding_degree_rejects_nonpositive():
    with pytest.raises(ValueError):
        counterion_binding_degree(0.0, 1.0)


def test_gibbs_free_energy_is_negative_for_spontaneous_micellization():
    """Sanity check: for a realistic ionic-surfactant CMC and counterion
    factor, deltaG_mic should be negative (spontaneous) and land in the
    commonly-reported -20 to -50 kJ/mol range for small-molecule
    surfactants."""
    x_cmc = cmc_to_mole_fraction(0.008)  # ~8 mM, SDS-like
    dg = gibbs_free_energy_micellization(x_cmc, 298.15, counterion_factor=1.6)
    assert dg < 0
    assert -60.0 < dg < -10.0


def test_gibbs_free_energy_matches_second_independent_literature_system():
    """Real SECOND independent literature validation for
    gibbs_free_energy_micellization (2026-09-10), beyond the single Fu
    et al. 2019 source referenced elsewhere in this project's literature
    validation notes: a different paper, different chemistry entirely
    (an alkylguanidinium cationic surfactant, not a conventional
    ammonium/sulfate headgroup) -- Bouchal, Hamel, Hesemann, In, Prelot
    & Zajac, Int. J. Mol. Sci. 17(2) (2016) 223, doi:10.3390/ijms17020223,
    dodecylguanidinium chloride (DDGC) in pure water, their own Table 2
    (conductivity): CMC = 6.2 mmol/kg, counterion binding degree beta =
    0.74, at 298 K, with their own reported deltaG_mic = -28.3 +/- 0.9
    kJ/mol.

    REVISED 2026-10-20 -- this is a documented DISAGREEMENT, not a
    validation. The paper's -28.3 kJ/mol is reproduced only by the factor
    (2 - beta) = 1.26 with beta = 0.74 taken as the BINDING degree. The
    standard pseudo-phase factor for a binding degree is (1 + beta) = 1.74
    (Zana, Langmuir 12 (1996) 1208), which gives about -39 kJ/mol. Either
    the paper's 0.74 is a binding degree used in the inverse convention,
    or its symbol means something else; this test used to treat the match
    as proof of (2 - beta), which is how the inverse convention entered
    this library. Both numbers are pinned so any future change is
    deliberate."""
    cmc_M = 6.2e-3  # mmol/kg ~= mM, dilute approximation
    beta = 0.74
    x_cmc = cmc_to_mole_fraction(cmc_M)
    standard = gibbs_free_energy_micellization(x_cmc, 298.0, 1.0 + beta)
    paper_convention = gibbs_free_energy_micellization(x_cmc, 298.0, 2.0 - beta)
    assert standard == pytest.approx(-39.2, abs=0.2)
    assert paper_convention == pytest.approx(-28.3, rel=0.01)   # the paper's own reported value


def test_gibbs_free_energy_counterion_factor_alpha_beta_identity_physical_limits():
    """(1 + beta) = (2 - alpha), alpha = 1 - beta (ionization vs. binding
    degree), checked at the physical limits: fully BOUND counterions
    (beta=1, alpha=0) leave solution with every surfactant ion -> factor 2;
    fully DISSOCIATED (beta=0, alpha=1) stay in solution -> factor 1, the
    nonionic value. Corrected 2026-10-20 (the inverse (2 - beta) form had
    these limits reversed)."""
    for alpha, beta, expected in [(0.0, 1.0, 2.0), (1.0, 0.0, 1.0)]:
        assert (1.0 + beta) == pytest.approx(expected)
        assert (2.0 - alpha) == pytest.approx(expected)
    # and at a real, in-between, sourced value (Bales 2001, SDS alpha=0.272):
    alpha_sds = 0.272
    beta_sds = 1.0 - alpha_sds
    assert (1.0 + beta_sds) == pytest.approx(2.0 - alpha_sds)
    assert (2.0 - alpha_sds) == pytest.approx(1.728)


def test_gibbs_free_energy_sds_with_bales_2001_counterion_factor():
    """Real worked example for SDS using a properly SOURCED
    counterion_factor (Bales, Messina, Vidal & Peric, J. Phys. Chem. B
    105 (2001) 6798-6804, alpha=0.272+/-0.017, PDF provided by the user
    2026-09-11) instead of an arbitrary illustrative value -- and the
    well-known salt-free SDS cmc0=0.0083 M already cited elsewhere in
    this project (Bales et al. 1998, Table 4, the same cmc0 used in that
    paper's own eq. 5 fit). No independently-published deltaG_mic number
    for this exact system was found to check the final result against,
    so this is a real, traceable calculation with a properly-sourced
    input. Corrected 2026-10-20 to the standard factor 2 - alpha = 1.728
    (was 1 + alpha = 1.272, the inverse convention, giving -27.77); the
    result now sits in the usual literature range for SDS (about -37 to
    -38 kJ/mol)."""
    counterion_factor_sds = 2.0 - 0.272  # Bales 2001's own alpha (ionization) for SDS
    x_cmc = cmc_to_mole_fraction(0.0083)
    dg = gibbs_free_energy_micellization(x_cmc, 298.15, counterion_factor_sds)
    assert dg < 0.0  # spontaneous, as it must be
    assert dg == pytest.approx(-37.73, abs=0.05)


# Tests for COUNTERION_BINDING_DEGREE_REFERENCE, added 2026-09-15 as
# part of the Tier 2 bottleneck-resolution work (real [MINE] path for
# item 1 -- see BOTTLENECK_RESOLUTION_PLAN.md).


def test_counterion_reference_alpha_beta_consistency_for_all_entries():
    for name, entry in COUNTERION_BINDING_DEGREE_REFERENCE.items():
        assert entry["alpha"] + entry["beta"] == pytest.approx(1.0), name
        assert 0.0 < entry["alpha"] < 1.0, name
        assert entry["alpha_confidence"] in ("high", "moderate", "low"), name


def test_counterion_reference_sds_matches_bales_2001_exactly():
    entry = COUNTERION_BINDING_DEGREE_REFERENCE["SDS"]
    assert entry["alpha"] == pytest.approx(0.272)
    assert entry["alpha_uncertainty"] == pytest.approx(0.017)
    assert entry["alpha_confidence"] == "high"


def test_counterion_reference_sds_usable_directly_in_gibbs_free_energy():
    """Real, usable end-to-end check: the reference table's own SDS
    entry plugged directly into gibbs_free_energy_micellization
    reproduces the exact same result as the already-verified hand-
    computed worked example above."""
    entry = COUNTERION_BINDING_DEGREE_REFERENCE["SDS"]
    counterion_factor = 1.0 + entry["beta"]          # = 2 - alpha
    assert counterion_factor == pytest.approx(2.0 - entry["alpha"])
    x_cmc = cmc_to_mole_fraction(0.0083)
    dg = gibbs_free_energy_micellization(x_cmc, 298.15, counterion_factor)
    assert dg == pytest.approx(-37.73, abs=0.05)


def test_counterion_reference_dtab_and_ctab_are_present_with_moderate_confidence():
    for name in ("DTAB", "CTAB"):
        entry = COUNTERION_BINDING_DEGREE_REFERENCE[name]
        assert entry["alpha_confidence"] == "moderate"
        assert 0.2 < entry["alpha"] < 0.35  # real, physically sensible range for these compounds


def test_counterion_reference_aot_real_primary_source_high_confidence():
    """Real Tier 2 gap closed 2026-09-15: AOT was searched for and not
    found via web search alone, then closed with a real primary-source
    PDF the user provided (Thapa et al. 2015, Table 1, Corrin-Harkins
    plot, NaCl medium, below c*)."""
    entry = COUNTERION_BINDING_DEGREE_REFERENCE["AOT"]
    assert entry["alpha"] == pytest.approx(0.61)
    assert entry["beta"] == pytest.approx(0.39)
    assert entry["alpha_confidence"] == "high"
    assert entry["alpha"] + entry["beta"] == pytest.approx(1.0)


def test_gibbs_free_energy_nonionic_vs_ionic_factor():
    # same CMC, but the ionic (counterion_factor > 1) case should be
    # more negative than the nonionic (factor=1) case, since 1 + beta > 1
    x_cmc = cmc_to_mole_fraction(0.008)
    dg_nonionic = gibbs_free_energy_micellization(x_cmc, 298.15, counterion_factor=1.0)
    dg_ionic = gibbs_free_energy_micellization(x_cmc, 298.15, counterion_factor=1.6)
    assert dg_ionic < dg_nonionic


def test_gibbs_free_energy_rejects_bad_mole_fraction():
    with pytest.raises(ValueError):
        gibbs_free_energy_micellization(0.0, 298.15)
    with pytest.raises(ValueError):
        gibbs_free_energy_micellization(1.0, 298.15)


def test_vant_hoff_enthalpy_round_trip():
    """Construct two (CMC, T) points consistent with a chosen deltaH,
    then verify the function recovers it -- mathematically guaranteed.

    Real bug found 2026-09-04: this test used to construct its synthetic
    data with the SAME (1/T1 - 1/T2) convention the function itself used,
    so it stayed green through a real sign bug in vant_hoff_enthalpy() --
    a round-trip test only proves internal self-consistency, not that the
    convention matches the true physics. Caught via an unaugmented Claude
    pilot answer disagreeing with the tool, then confirmed independently
    from the Gibbs-Helmholtz equation and against this project's own
    existing Fu et al. 2019 literature-validation entry (whose two-point
    estimates had the opposite sign from the paper's own local-derivative
    value, previously misread as pure two-point-vs-polynomial magnitude
    noise). Fixed here to use the corrected (1/T2 - 1/T1) convention,
    matching the now-fixed function; see test below for an independent,
    non-round-trip sign check against real physical reasoning.
    """
    R = 8.314462618
    T1, T2 = 293.15, 313.15
    x1 = 1.5e-4
    delta_h_true_kJ = -15.0  # typical small negative micellization enthalpy

    delta_inv_T = (1.0 / T2) - (1.0 / T1)
    ln_x2 = math.log(x1) + (delta_h_true_kJ * 1000.0 * delta_inv_T) / R
    x2 = math.exp(ln_x2)

    delta_h_recovered = vant_hoff_enthalpy(x1, T1, x2, T2)
    assert delta_h_recovered == pytest.approx(delta_h_true_kJ, abs=1e-6)


def test_vant_hoff_enthalpy_sign_matches_physical_reasoning():
    """Independent, non-round-trip sign check: if CMC DECREASES as
    temperature RISES, micellization is more favourable at higher T,
    which by the van't Hoff/Le Chatelier argument must be an ENDOTHERMIC
    (positive deltaH) process -- this is real physical reasoning, not
    derived from the function under test, so it cannot share a bug with
    it the way the round-trip test above could (and did)."""
    T1, T2 = 293.15, 313.15
    cmc1, cmc2 = 4.0e-5, 3.8e-5  # CMC2 < CMC1: CMC falls as T rises
    delta_h = vant_hoff_enthalpy(cmc1, T1, cmc2, T2)
    assert delta_h > 0, "CMC falling with rising T must give a positive (endothermic) deltaH"


def test_vant_hoff_enthalpy_rejects_equal_temperatures():
    with pytest.raises(ValueError):
        vant_hoff_enthalpy(1e-4, 298.15, 2e-4, 298.15)


def test_vant_hoff_enthalpy_default_counterion_factor_is_backward_compatible():
    """counterion_factor added (ported from SurfactantKit, real gap found
    via GZ-26 in benchmark/groundzero_50: this function had no way to
    report a self-consistent deltaH for an ionic surfactant whose deltaG
    already used gibbs_free_energy_micellization's counterion_factor).
    Default 1.0 must reproduce the pre-existing (nonionic) behavior
    exactly -- every call site written before this parameter existed
    must be unaffected."""
    T1, T2 = 293.15, 313.15
    cmc1, cmc2 = 4.0e-5, 3.8e-5
    assert vant_hoff_enthalpy(cmc1, T1, cmc2, T2) == vant_hoff_enthalpy(cmc1, T1, cmc2, T2, counterion_factor=1.0)


def test_vant_hoff_enthalpy_counterion_factor_scales_linearly():
    """Real consequence of the Gibbs-Helmholtz derivation in this
    function's own docstring: deltaH = counterion_factor * R *
    d(ln X_cmc)/d(1/T) -- counterion_factor must multiply the raw
    (factor=1) result exactly, not approximately, since it's a constant
    prefactor on a linear formula, not a fitted correction."""
    T1, T2 = 283.15, 318.15
    cmc1, cmc2 = 9.93e-3 / (9.93e-3 + 55.5), 11.47e-3 / (11.47e-3 + 55.5)
    raw = vant_hoff_enthalpy(cmc1, T1, cmc2, T2)
    factor = 2.0 - 0.268  # real alpha=0.268 for [C12mpy][Br], Fu et al. 2019
    corrected = vant_hoff_enthalpy(cmc1, T1, cmc2, T2, counterion_factor=factor)
    assert corrected == pytest.approx(raw * factor, rel=1e-9)
    # Reproduces the exact GZ-26 wide-span value an augmented Claude Opus
    # 4.8 answer derived by hand and self-corrected to, catching that the
    # gold reference's own deltaH had been left at the nonionic (factor=1)
    # tool default while its deltaG used this same 1.732 factor.
    assert corrected == pytest.approx(-5.3428, abs=1e-3)


# --- vant_hoff_multi_point_fit ----------------------------------------------
# Source: Kantonen, Henriksen & Gilson, BBA Gen. Subj. 1862 (2018) 692-704,
# Eq. 4 (verified via WebFetch before implementing, not guessed) --
# deltaG(Tj) = deltaH(Tr) - Tj*deltaS(Tr) + deltaCp(Tr)*[(Tj-Tr) - Tj*ln(Tj/Tr)].
# Validated via mathematically-guaranteed round trips: construct X_cmc(Tj)
# for chosen true deltaH/deltaS/deltaCp via this exact equation, verify the
# fit recovers them.


def _synthetic_series(true_dH_kJ, true_dS_J, true_dCp_J, tr, temps, counterion_factor=1.0):
    true_dS_kJ = true_dS_J / 1000.0
    true_dCp_kJ = true_dCp_J / 1000.0
    x_values = []
    for t in temps:
        g = (t - tr) - t * math.log(t / tr)
        dg = true_dH_kJ - t * true_dS_kJ + true_dCp_kJ * g
        x_values.append(math.exp(dg * 1000.0 / (counterion_factor * R_GAS * t)))
    return x_values


def test_vant_hoff_multi_point_fit_recovers_true_parameters_round_trip():
    true_dH_kJ, true_dS_J, true_dCp_J = -20.0, -50.0, -300.0
    tr = 298.15
    temps = [278.15, 288.15, 298.15, 308.15, 318.15]
    x_values = _synthetic_series(true_dH_kJ, true_dS_J, true_dCp_J, tr, temps)

    result = vant_hoff_multi_point_fit(x_values, temps, reference_temperature_K=tr)
    assert result.delta_H_kJ_per_mol == pytest.approx(true_dH_kJ, abs=1e-6)
    assert result.delta_S_J_per_mol_K == pytest.approx(true_dS_J, abs=1e-4)
    assert result.delta_Cp_J_per_mol_K == pytest.approx(true_dCp_J, abs=1e-3)
    assert result.r_squared == pytest.approx(1.0, abs=1e-9)
    assert result.n_points == 5
    assert result.reference_temperature_K == pytest.approx(tr)


def test_vant_hoff_multi_point_fit_default_reference_is_mean_temperature():
    true_dH_kJ, true_dS_J, true_dCp_J = -18.0, -40.0, -250.0
    tr = 300.0
    temps = [283.15, 293.15, 303.15, 313.15]
    x_values = _synthetic_series(true_dH_kJ, true_dS_J, true_dCp_J, tr, temps)

    result = vant_hoff_multi_point_fit(x_values, temps)  # no explicit reference
    assert result.reference_temperature_K == pytest.approx(sum(temps) / len(temps))
    # deltaCp and the fitted deltaG curve are reference-independent by construction:
    assert result.delta_Cp_J_per_mol_K == pytest.approx(true_dCp_J, abs=1e-3)
    assert result.r_squared == pytest.approx(1.0, abs=1e-9)


def test_vant_hoff_multi_point_fit_with_ionic_counterion_factor():
    true_dH_kJ, true_dS_J, true_dCp_J = -22.0, -55.0, -320.0
    tr = 298.15
    temps = [288.15, 298.15, 308.15, 318.15]
    factor = 1.6
    x_values = _synthetic_series(true_dH_kJ, true_dS_J, true_dCp_J, tr, temps, counterion_factor=factor)

    result = vant_hoff_multi_point_fit(x_values, temps, counterion_factor=factor, reference_temperature_K=tr)
    assert result.delta_H_kJ_per_mol == pytest.approx(true_dH_kJ, abs=1e-6)
    assert result.delta_Cp_J_per_mol_K == pytest.approx(true_dCp_J, abs=1e-3)


def test_vant_hoff_multi_point_fit_rejects_bad_inputs():
    with pytest.raises(ValueError):
        vant_hoff_multi_point_fit([1e-4, 2e-4], [293.15, 303.15])  # too few points
    with pytest.raises(ValueError):
        vant_hoff_multi_point_fit([1e-4, 2e-4, 3e-4], [293.15, 303.15])  # mismatched lengths
    with pytest.raises(ValueError):
        vant_hoff_multi_point_fit([0.0, 2e-4, 3e-4], [293.15, 303.15, 313.15])  # bad mole fraction
    with pytest.raises(ValueError):
        vant_hoff_multi_point_fit([1e-4, 2e-4, 3e-4], [293.15, 0.0, 313.15])  # bad temperature
    with pytest.raises(ValueError):
        vant_hoff_multi_point_fit([1e-4, 2e-4, 3e-4], [298.15, 298.15, 298.15])  # degenerate (all same T)


# --- mass_action_free_energy_of_micellization (alternative-methods sweep,
# 2026-09-10) -- eq. 21 of Boneva-Aroca et al., RSC Advances 13 (2023) 9387
# (open access, quoted verbatim via WebFetch, not guessed), the same
# well-known classical single-equilibrium mass-action result this project's
# own literature review flagged via Rusanov 2014 (paywalled primary source).


def test_mass_action_matches_cited_formula_directly():
    """Independent recomputation in the test itself, not copy-pasted
    from the implementation."""
    x_cmc, n, T = 1.5e-4, 60, 298.15
    R = 8.314462618
    expected_kJ = (R * T * ((1.0 + 1.0 / n) * math.log(x_cmc) - (1.0 / n) * math.log(n))) / 1000.0
    assert mass_action_free_energy_of_micellization(x_cmc, n, T) == pytest.approx(expected_kJ)


def test_mass_action_converges_to_pseudo_phase_separation_at_large_n():
    """The real, verified property tying the two models together: as n
    grows, mass_action_free_energy_of_micellization must converge to
    gibbs_free_energy_micellization(counterion_factor=1) -- the paper's
    own stated large-n limit, not just an incidental numerical match."""
    x_cmc, T = 2.0e-4, 298.15
    pseudo_phase = gibbs_free_energy_micellization(x_cmc, T, counterion_factor=1.0)
    for n in (10, 100, 1000, 100000):
        mass_action = mass_action_free_energy_of_micellization(x_cmc, n, T)
        diff = abs(mass_action - pseudo_phase)
        if n == 10:
            prev_diff = diff
        else:
            assert diff < prev_diff
            prev_diff = diff
    assert mass_action_free_energy_of_micellization(x_cmc, 100000, T) == pytest.approx(pseudo_phase, abs=1e-3)


def test_mass_action_diverges_meaningfully_from_pseudo_phase_at_small_n():
    """For a small aggregation number, the correction should be real and
    non-negligible (not silently ~0), matching the source paper's own
    finding that the correction remains ~0.1 (in RT units) even at
    n~100 -- here at a much smaller, more divergent n."""
    x_cmc, n, T = 2.0e-4, 15, 298.15
    pseudo_phase = gibbs_free_energy_micellization(x_cmc, T, counterion_factor=1.0)
    mass_action = mass_action_free_energy_of_micellization(x_cmc, n, T)
    assert abs(mass_action - pseudo_phase) > 0.5  # kJ/mol, a real, non-trivial difference


def test_mass_action_rejects_bad_inputs():
    with pytest.raises(ValueError):
        mass_action_free_energy_of_micellization(0.0, 50, 298.15)
    with pytest.raises(ValueError):
        mass_action_free_energy_of_micellization(1.5e-4, 1, 298.15)  # n must be > 1
    with pytest.raises(ValueError):
        mass_action_free_energy_of_micellization(1.5e-4, 50, 0.0)


# --- critical_micellization_degree (new capability, 2026-09-11) -----------
# Source: Rusanov, "The Mass-Action-Law Theory of Micellization Revisited,"
# Langmuir 30(48) (2014) 14443-14451, eqs. (15)-(16), PDF provided by the user.


def test_cmd_d2alpha_dc2_matches_papers_own_worked_example():
    # paper states alpha_m = 0.061 at n = 100 for this definition
    assert critical_micellization_degree(100, "d2alpha_dc2") == pytest.approx(0.061, abs=0.001)


def test_cmd_d2alpha_dlnc2_matches_papers_own_worked_example():
    # paper states alpha_m = 0.091 at n = 100 for this definition
    assert critical_micellization_degree(100, "d2alpha_dlnc2") == pytest.approx(0.091, abs=0.001)


def test_cmd_dlnc2_definition_always_exceeds_dc2_definition_for_n_gt_1():
    # real, checkable consistency property across n, not just at n=100
    for n in (5, 10, 30, 50, 100, 500, 2000):
        dc2 = critical_micellization_degree(n, "d2alpha_dc2")
        dlnc2 = critical_micellization_degree(n, "d2alpha_dlnc2")
        assert dlnc2 > dc2


def test_cmd_vanishes_as_n_grows():
    small_n = critical_micellization_degree(10, "d2alpha_dlnc2")
    large_n = critical_micellization_degree(10000, "d2alpha_dlnc2")
    assert large_n < small_n
    assert large_n == pytest.approx(0.0, abs=0.02)


def test_cmd_rejects_bad_inputs():
    with pytest.raises(ValueError):
        critical_micellization_degree(1, "d2alpha_dlnc2")  # n must be > 1
    with pytest.raises(ValueError):
        critical_micellization_degree(1.5, "d2alpha_dc2")  # unphysical (negative) for n < 2
    with pytest.raises(ValueError):
        critical_micellization_degree(50, "not_a_real_definition")


def test_entropy_micellization_completes_triad_by_construction():
    # deltaS = (deltaH - deltaG) / T -- purely definitional, must hold exactly
    dh, dg, T = -15.0, -35.0, 298.15
    ds = entropy_micellization(dg, dh, T)
    assert ds == pytest.approx(((dh - dg) / T) * 1000.0)


def test_entropy_micellization_rejects_bad_temperature():
    with pytest.raises(ValueError):
        entropy_micellization(-35.0, -15.0, 0.0)


# ---------------------------------------------------------------------------
# Batch sign-off (audit B8): optional deltaCp term in entropy_micellization.
# ---------------------------------------------------------------------------

def test_entropy_micellization_delta_cp_term_matches_kirchhoff():
    """deltaS(T) = deltaS(Tr) + deltaCp*ln(T/Tr); exact values pinned so
    the wiring can't drift silently."""
    ds = entropy_micellization(
        -10.0, -5.0, 350.0,
        delta_cp_J_mol_K=-100.0,
        reference_temperature_K=298.15,
    )
    expected = (5.0 / 298.15) * 1000.0 - 100.0 * math.log(350.0 / 298.15)
    assert ds == pytest.approx(expected, rel=1e-12)
    assert ds == pytest.approx(0.7364382819682085, rel=1e-9)


def test_entropy_micellization_reference_equals_temperature_is_classic():
    base = entropy_micellization(-10.0, -5.0, 318.15)
    with_cp = entropy_micellization(
        -10.0, -5.0, 318.15,
        delta_cp_J_mol_K=-100.0,
        reference_temperature_K=318.15,
    )
    assert with_cp == pytest.approx(base, rel=1e-12)


def test_entropy_micellization_half_specified_delta_cp_fails_fast():
    with pytest.raises(ValueError, match="delta_cp_J_mol_K.*reference_temperature_K"):
        entropy_micellization(-10.0, -5.0, 318.15, delta_cp_J_mol_K=-100.0)
    with pytest.raises(ValueError, match="delta_cp_J_mol_K.*reference_temperature_K"):
        entropy_micellization(-10.0, -5.0, 318.15, reference_temperature_K=298.15)
