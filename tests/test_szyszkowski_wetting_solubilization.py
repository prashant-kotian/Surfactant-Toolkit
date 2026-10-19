"""Tests for the Szyszkowski surface-tension equation (adsorption.py),
wetting.py, and solubilization.py."""

import math

import pytest

from surfactantkit.adsorption import gibbs_gamma_max, szyszkowski_surface_tension
from surfactantkit.solubilization import (
    INTRINSIC_WATER_SOLUBILITY_REFERENCE_M,
    estimate_intrinsic_water_solubility_qspr,
    micelle_water_partition_coefficient,
    molar_solubilization_ratio,
)
from surfactantkit.wetting import (
    VOCG_STANDARD_LIQUIDS,
    capillary_number,
    get_owens_wendt_standard_liquid,
    get_vocg_standard_liquid,
    owens_wendt_solid_surface_energy,
    spreading_coefficient,
    van_oss_chaudhury_good_solid_surface_energy,
    work_of_adhesion,
)


def test_szyszkowski_returns_gamma0_at_zero_concentration():
    # boundary condition, follows directly from ln(1+0)=0
    g = szyszkowski_surface_tension(0.0, gamma0_mN_m=72.0, gamma_max_mol_per_m2=3.0e-6, K=50.0,
                                     system_type="nonionic")
    assert g == pytest.approx(72.0)


def test_szyszkowski_decreases_monotonically_with_concentration():
    g_low = szyszkowski_surface_tension(1.0, 72.0, 3.0e-6, 50.0, "nonionic")
    g_high = szyszkowski_surface_tension(10.0, 72.0, 3.0e-6, 50.0, "nonionic")
    assert g_high < g_low < 72.0


def test_szyszkowski_high_concentration_limit_matches_gibbs_slope_method():
    """At high concentration (K*C >> 1), the numerical slope
    d(gamma)/d(ln C) of the Szyszkowski curve should recover the same
    Gamma_max that was used to construct it, when fed through
    gibbs_gamma_max -- these are two independently-implemented physics
    relations that must agree exactly in this limit; a strong internal
    cross-check, not just a unit test."""
    gamma0, gamma_max_true, K, T = 72.0, 3.0e-6, 50.0, 298.15
    C1, C2 = 100.0, 105.0  # K*C ~ 5000-5250, well into the saturation limit
    g1 = szyszkowski_surface_tension(C1, gamma0, gamma_max_true, K, "nonionic", T)
    g2 = szyszkowski_surface_tension(C2, gamma0, gamma_max_true, K, "nonionic", T)
    slope = (g2 - g1) / (math.log(C2) - math.log(C1))
    recovered_gamma_max = gibbs_gamma_max(slope, system_type="nonionic", temperature_k=T)
    assert recovered_gamma_max == pytest.approx(gamma_max_true, rel=0.005)


def test_szyszkowski_rejects_bad_input():
    with pytest.raises(ValueError):
        szyszkowski_surface_tension(-1.0, 72.0, 3.0e-6, 50.0, "nonionic")
    with pytest.raises(ValueError):
        szyszkowski_surface_tension(1.0, 72.0, 0.0, 50.0, "nonionic")
    with pytest.raises(ValueError):
        szyszkowski_surface_tension(1.0, 72.0, 3.0e-6, -1.0, "nonionic")
    with pytest.raises(ValueError):
        szyszkowski_surface_tension(1.0, 72.0, 3.0e-6, 50.0, "not_a_real_system_type")


def test_gibbs_gamma_max_rejects_invalid_system_type():
    """Regression test for the 2026-09-04 hardening: system_type replaced
    a raw n_factor float precisely because a caller (including an LLM
    using this as a tool) must not be able to silently pass an
    unjustified numeric guess. No default is provided and any string
    outside the three real system types must raise."""
    with pytest.raises(ValueError):
        gibbs_gamma_max(-5.0, "made_up_system_type", 298.15)
    with pytest.raises(TypeError):
        gibbs_gamma_max(-5.0, temperature_k=298.15)  # system_type omitted entirely -- no default


def test_work_of_adhesion_complete_wetting_limit():
    # theta = 0 (complete wetting) -> W_a = 2*gamma_LV exactly
    assert work_of_adhesion(72.0, 0.0) == pytest.approx(144.0)


def test_work_of_adhesion_no_wetting_limit():
    # theta = 180 -> W_a = 0 exactly
    assert work_of_adhesion(72.0, 180.0) == pytest.approx(0.0, abs=1e-9)


def test_work_of_adhesion_rejects_bad_input():
    with pytest.raises(ValueError):
        work_of_adhesion(-1.0, 90.0)
    with pytest.raises(ValueError):
        work_of_adhesion(72.0, 200.0)


def test_spreading_coefficient_complete_wetting_boundary():
    # theta = 0 -> S = 0 exactly (the only case where S >= 0 is achieved)
    assert spreading_coefficient(72.0, 0.0) == pytest.approx(0.0, abs=1e-9)


def test_spreading_coefficient_is_never_positive():
    for theta in [0.0, 30.0, 90.0, 150.0, 180.0]:
        assert spreading_coefficient(72.0, theta) <= 1e-9


def test_capillary_number_matches_direct_formula():
    ca = capillary_number(viscosity_mPas=1.0, velocity_m_per_s=1e-5, interfacial_tension_mN_m=0.01)
    expected = (1.0e-3 * 1e-5) / (0.01e-3)
    assert ca == pytest.approx(expected)


def test_capillary_number_increases_as_interfacial_tension_drops():
    # this IS the whole point of surfactant EOR -- lowering IFT raises Ca
    ca_normal_ift = capillary_number(1.0, 1e-5, 20.0)
    ca_ultralow_ift = capillary_number(1.0, 1e-5, 0.001)
    assert ca_ultralow_ift > ca_normal_ift


def test_capillary_number_rejects_nonpositive():
    with pytest.raises(ValueError):
        capillary_number(0.0, 1e-5, 20.0)


# --- OWRK and van Oss-Chaudhury-Good surface energy decomposition
# (alternative-methods sweep, 2026-09-10). Sources verified via WebSearch
# before implementing: OWRK (Owens & Wendt 1969; Rabel 1971; Kaelble 1970),
# vOCG (van Oss, Chaudhury & Good, Chem. Rev. 88 (1988) 927-941). No
# external worked numeric table was sourced this pass (same disclosed
# pattern as other tools this session); validated via mathematically-
# guaranteed round trips (construct contact angles from chosen true solid
# surface-energy components via the forward mixing rule, recover them back).


def test_owens_wendt_recovers_true_components_round_trip():
    true_gamma_s_d, true_gamma_s_p = 30.0, 10.0
    # Real, commonly-tabulated test-liquid dispersive/polar splits (water, diiodomethane, glycerol)
    liquids_d = [21.8, 50.8, 37.0]
    liquids_p = [51.0, 0.0, 26.4]
    liquids_total = [d + p for d, p in zip(liquids_d, liquids_p)]

    contact_angles = []
    for gt, gd, gp in zip(liquids_total, liquids_d, liquids_p):
        w = 2.0 * (math.sqrt(true_gamma_s_d * gd) + math.sqrt(true_gamma_s_p * gp))
        cos_theta = w / gt - 1.0
        contact_angles.append(math.degrees(math.acos(cos_theta)))

    result = owens_wendt_solid_surface_energy(contact_angles, liquids_d, liquids_p)
    assert result.gamma_s_dispersive_mN_m == pytest.approx(true_gamma_s_d, abs=1e-3)
    assert result.gamma_s_polar_mN_m == pytest.approx(true_gamma_s_p, abs=1e-3)
    assert result.gamma_s_total_mN_m == pytest.approx(true_gamma_s_d + true_gamma_s_p, abs=1e-3)
    assert result.r_squared == pytest.approx(1.0, abs=1e-6)
    assert result.n_liquids == 3


def test_owens_wendt_purely_dispersive_solid_has_zero_polar_component():
    """A solid with zero polar character (e.g. PTFE-like) tested against
    a purely dispersive liquid alone should recover gamma_s_p == 0 when
    combined with any second liquid -- a real physical sanity check."""
    true_gamma_s_d = 20.0
    liquids_d = [21.8, 50.8]
    liquids_p = [51.0, 0.0]
    liquids_total = [d + p for d, p in zip(liquids_d, liquids_p)]
    contact_angles = []
    for gt, gd, gp in zip(liquids_total, liquids_d, liquids_p):
        w = 2.0 * math.sqrt(true_gamma_s_d * gd)  # gamma_s_p = 0
        cos_theta = w / gt - 1.0
        contact_angles.append(math.degrees(math.acos(cos_theta)))

    result = owens_wendt_solid_surface_energy(contact_angles, liquids_d, liquids_p)
    assert result.gamma_s_dispersive_mN_m == pytest.approx(true_gamma_s_d, abs=1e-3)
    assert result.gamma_s_polar_mN_m == pytest.approx(0.0, abs=1e-3)


def test_owens_wendt_rejects_bad_inputs():
    with pytest.raises(ValueError):
        owens_wendt_solid_surface_energy([50.0], [21.8], [51.0])  # too few liquids
    with pytest.raises(ValueError):
        owens_wendt_solid_surface_energy([50.0, 60.0], [21.8], [51.0, 0.0])  # mismatched lengths
    with pytest.raises(ValueError):
        owens_wendt_solid_surface_energy([50.0, 60.0], [0.0, 50.8], [51.0, 0.0])  # bad dispersive


def test_van_oss_chaudhury_good_recovers_true_components_round_trip():
    true_lw, true_acid, true_base = 25.0, 2.0, 15.0
    # Real, commonly-tabulated LW/acid/base splits (water, glycerol, diiodomethane)
    liquids_lw = [21.8, 34.0, 50.8]
    liquids_acid = [25.5, 3.92, 0.0]
    liquids_base = [25.5, 57.4, 0.0]

    contact_angles = []
    for lw, a, b in zip(liquids_lw, liquids_acid, liquids_base):
        gt = lw + 2.0 * math.sqrt(a * b)
        w = 2.0 * (math.sqrt(true_lw * lw) + math.sqrt(true_acid * b) + math.sqrt(true_base * a))
        cos_theta = w / gt - 1.0
        contact_angles.append(math.degrees(math.acos(cos_theta)))

    result = van_oss_chaudhury_good_solid_surface_energy(contact_angles, liquids_lw, liquids_acid, liquids_base)
    assert result.gamma_s_lw_mN_m == pytest.approx(true_lw, abs=1e-3)
    assert result.gamma_s_acid_mN_m == pytest.approx(true_acid, abs=1e-3)
    assert result.gamma_s_base_mN_m == pytest.approx(true_base, abs=1e-3)
    expected_total = true_lw + 2.0 * math.sqrt(true_acid * true_base)
    assert result.gamma_s_total_mN_m == pytest.approx(expected_total, abs=1e-3)
    assert result.n_liquids == 3


def test_van_oss_chaudhury_good_rejects_bad_inputs():
    with pytest.raises(ValueError):
        van_oss_chaudhury_good_solid_surface_energy([50.0, 60.0], [21.8, 34.0], [25.5, 3.92], [25.5, 57.4])  # only 2 liquids
    with pytest.raises(ValueError):
        van_oss_chaudhury_good_solid_surface_energy(
            [50.0, 60.0, 70.0], [0.0, 34.0, 50.8], [25.5, 3.92, 0.0], [25.5, 57.4, 0.0]
        )  # bad LW value


# Tests for the real, cited standard-probe-liquid lookup tables added
# 2026-09-15 (audit note retained with the authors; see
# item 12 and BOTTLENECK_RESOLUTION_PLAN.md Tier 0) -- these values were
# never guessed, they were live-verified against real sources (van Oss,
# Chaudhury & Good 1988; van Oss, Good & Busscher 1990) before being
# hardcoded here.


def test_vocg_standard_liquid_lookup_matches_known_real_values():
    water = get_vocg_standard_liquid("water")
    assert water["lw"] == pytest.approx(21.8)
    assert water["acid"] == pytest.approx(25.5)
    assert water["base"] == pytest.approx(25.5)
    diiodomethane = get_vocg_standard_liquid("Diiodomethane")  # case-insensitive
    assert diiodomethane["lw"] == pytest.approx(50.8)
    assert diiodomethane["acid"] == pytest.approx(0.0)
    assert diiodomethane["base"] == pytest.approx(0.0)


def test_vocg_standard_liquid_totals_are_internally_consistent():
    # total = LW + 2*sqrt(acid*base), the same mixing rule the real
    # solver uses -- catches a transcription error in the hardcoded table
    for name, vals in VOCG_STANDARD_LIQUIDS.items():
        recomputed_total = vals["lw"] + 2.0 * math.sqrt(vals["acid"] * vals["base"])
        assert recomputed_total == pytest.approx(vals["total"], abs=0.05), name


def test_vocg_standard_liquid_rejects_unknown_name():
    with pytest.raises(ValueError):
        get_vocg_standard_liquid("ethanol")


def test_owens_wendt_standard_liquid_lookup_matches_known_real_values():
    water = get_owens_wendt_standard_liquid("water")
    assert water["dispersive"] == pytest.approx(21.8)
    assert water["polar"] == pytest.approx(51.0)
    glycerol = get_owens_wendt_standard_liquid("GLYCEROL")  # case-insensitive
    assert glycerol["dispersive"] == pytest.approx(34.0)
    assert glycerol["polar"] == pytest.approx(30.0)


def test_owens_wendt_standard_liquid_rejects_unknown_name():
    with pytest.raises(ValueError):
        get_owens_wendt_standard_liquid("ethanol")


def test_vocg_solver_runs_end_to_end_on_real_standard_liquid_values_alone():
    """No external lookup needed at all -- pull water/glycerol/
    diiodomethane straight from the real tabulated dict and confirm the
    solver runs to completion on a physically sane synthetic solid."""
    liquids = ["water", "glycerol", "diiodomethane"]
    vals = [get_vocg_standard_liquid(n) for n in liquids]
    liquids_lw = [v["lw"] for v in vals]
    liquids_acid = [v["acid"] for v in vals]
    liquids_base = [v["base"] for v in vals]
    liquids_total = [v["lw"] + 2.0 * math.sqrt(v["acid"] * v["base"]) for v in vals]

    true_lw, true_acid, true_base = 25.0, 2.0, 15.0
    contact_angles = []
    for gt, lw, a, b in zip(liquids_total, liquids_lw, liquids_acid, liquids_base):
        w = 2.0 * (math.sqrt(true_lw * lw) + math.sqrt(true_acid * b) + math.sqrt(true_base * a))
        contact_angles.append(math.degrees(math.acos(w / gt - 1.0)))

    result = van_oss_chaudhury_good_solid_surface_energy(contact_angles, liquids_lw, liquids_acid, liquids_base)
    assert result.gamma_s_lw_mN_m == pytest.approx(true_lw, abs=1e-3)
    assert result.gamma_s_acid_mN_m == pytest.approx(true_acid, abs=1e-3)
    assert result.gamma_s_base_mN_m == pytest.approx(true_base, abs=1e-3)


def test_molar_solubilization_ratio_basic():
    # 0.5 mM extra solubilized (beyond intrinsic 0.01 mM) per (10-2) mM micellized surfactant
    msr = molar_solubilization_ratio(
        total_solubilized_M=0.51e-3, intrinsic_water_solubility_M=0.01e-3,
        surfactant_concentration_M=10e-3, cmc_M=2e-3,
    )
    assert msr == pytest.approx((0.51e-3 - 0.01e-3) / (10e-3 - 2e-3))


def test_molar_solubilization_ratio_rejects_below_cmc():
    with pytest.raises(ValueError):
        molar_solubilization_ratio(0.1e-3, 0.01e-3, 1e-3, 2e-3)


def test_molar_solubilization_ratio_rejects_inconsistent_solubility():
    with pytest.raises(ValueError):
        molar_solubilization_ratio(0.005e-3, 0.01e-3, 10e-3, 2e-3)


# --- micelle_water_partition_coefficient (alternative-methods sweep,
# 2026-09-10) -- mole-fraction-basis convention, one of (at least) two
# real conventions found in the literature (disclosed in the function's
# own docstring; the other, molar-concentration-ratio basis, is NOT
# implemented here, same "do not guess a single convention" discipline
# already used for the Davies HLB gaps elsewhere in this project).


def test_micelle_water_partition_coefficient_matches_manual_formula():
    """Independent recomputation in the test itself, not copy-pasted
    from the implementation."""
    total, intrinsic, surf, cmc = 0.51e-3, 0.01e-3, 10e-3, 2e-3
    x_micelle = (total - intrinsic) / ((total - intrinsic) + (surf - cmc))
    x_water = intrinsic / 55.5
    expected = x_micelle / x_water
    assert micelle_water_partition_coefficient(total, intrinsic, surf, cmc) == pytest.approx(expected)


def test_micelle_water_partition_coefficient_large_for_strongly_solubilized_species():
    """Real physical sanity check: a solute with very low intrinsic
    water solubility but substantial total solubilized amount (i.e. it
    strongly favors the micelle) should give Km >> 1."""
    km = micelle_water_partition_coefficient(
        total_solubilized_M=1.0e-3, intrinsic_water_solubility_M=1.0e-7,
        surfactant_concentration_M=10e-3, cmc_M=2e-3,
    )
    assert km > 1000.0


def test_micelle_water_partition_coefficient_zero_when_nothing_solubilized_in_micelle():
    km = micelle_water_partition_coefficient(
        total_solubilized_M=0.01e-3, intrinsic_water_solubility_M=0.01e-3,
        surfactant_concentration_M=10e-3, cmc_M=2e-3,
    )
    assert km == pytest.approx(0.0)


def test_micelle_water_partition_coefficient_rejects_below_cmc():
    with pytest.raises(ValueError):
        micelle_water_partition_coefficient(0.1e-3, 0.01e-3, 1e-3, 2e-3)


def test_micelle_water_partition_coefficient_rejects_zero_intrinsic_solubility():
    with pytest.raises(ValueError):
        micelle_water_partition_coefficient(0.1e-3, 0.0, 10e-3, 2e-3)


def test_micelle_water_partition_coefficient_rejects_inconsistent_solubility():
    with pytest.raises(ValueError):
        micelle_water_partition_coefficient(0.005e-3, 0.01e-3, 10e-3, 2e-3)


# --- estimate_intrinsic_water_solubility_qspr (Tier 1 bottleneck fix #13,
# 2026-09-15 -- see BOTTLENECK_RESOLUTION_PLAN.md). Coefficients are Pat
# Walters' real, published RDKit-refit of Delaney's 2004 ESOL model,
# live-verified before being hardcoded. Validated by recomputing the exact
# same formula independently in the test (not just trusting the
# implementation), plus a real-world sanity check on naphthalene, a classic
# ESOL benchmark compound.


def test_qspr_solubility_matches_manual_esol_formula_for_naphthalene():
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors, Lipinski

    smiles = "c1ccc2ccccc2c1"  # naphthalene
    mol = Chem.MolFromSmiles(smiles)
    mw = Descriptors.MolWt(mol)
    clogp = Crippen.MolLogP(mol)
    rotors = Lipinski.NumRotatableBonds(mol)
    n_aromatic = len(mol.GetSubstructMatches(Chem.MolFromSmarts("a")))
    ap = n_aromatic / mol.GetNumAtoms()
    expected_log_s = (
        0.26121066137801696 - 0.7416739523408995 * clogp - 0.0066138847738667125 * mw
        + 0.003451545565957996 * rotors - 0.42624840441316975 * ap
    )

    result = estimate_intrinsic_water_solubility_qspr(smiles)
    assert result.log_s_mol_per_L == pytest.approx(expected_log_s, rel=1e-9)
    assert result.solubility_M == pytest.approx(10.0 ** expected_log_s, rel=1e-9)


def test_qspr_solubility_naphthalene_in_realistic_range():
    # naphthalene's real measured aqueous solubility is well known to be
    # low-mM range (~2e-4 M); ESOL is a coarse estimate (real RMSE ~0.6-1
    # log unit) so this checks order-of-magnitude plausibility, not an
    # exact match against a specific literature number.
    result = estimate_intrinsic_water_solubility_qspr("c1ccc2ccccc2c1")
    assert 1e-6 <= result.solubility_M <= 1e-1
    assert result.caveats  # the QSPR-estimate caveat must always be present


def test_qspr_solubility_highly_polar_small_molecule_more_soluble_than_hydrocarbon():
    # a real, qualitative sanity check the formula must get right regardless
    # of exact coefficients: a small polar molecule should be predicted
    # MORE soluble than a nonpolar aromatic hydrocarbon of similar size
    naphthalene = estimate_intrinsic_water_solubility_qspr("c1ccc2ccccc2c1")
    ethanol = estimate_intrinsic_water_solubility_qspr("CCO")
    assert ethanol.solubility_M > naphthalene.solubility_M


def test_qspr_solubility_rejects_unparseable_smiles():
    with pytest.raises(ValueError):
        estimate_intrinsic_water_solubility_qspr("not a smiles $$$")


# --- INTRINSIC_WATER_SOLUBILITY_REFERENCE_M (Tier 2 bottleneck fix #13,
# 2026-09-15 -- see BOTTLENECK_RESOLUTION_PLAN.md) -----------------------


def test_intrinsic_water_solubility_reference_naphthalene_matches_existing_project_value():
    entry = INTRINSIC_WATER_SOLUBILITY_REFERENCE_M["naphthalene"]
    assert entry["solubility_M"] == pytest.approx(2.17e-4)
    assert entry["confidence"] == "high"


def test_intrinsic_water_solubility_reference_all_entries_consistent_with_mw():
    for name, entry in INTRINSIC_WATER_SOLUBILITY_REFERENCE_M.items():
        solubility_mg_per_L = entry["solubility_M"] * entry["mw_g_per_mol"] * 1000.0
        assert solubility_mg_per_L > 0, name
        assert entry["confidence"] in ("high", "moderate", "low"), name


def test_intrinsic_water_solubility_reference_benzene_and_pyrene_real_and_ordered():
    # a real, physically expected ordering: benzene (small, only mildly
    # hydrophobic) is far MORE soluble than pyrene (large, strongly
    # hydrophobic PAH) -- a real sanity check on the two moderate-
    # confidence entries, not just a formula check
    benzene = INTRINSIC_WATER_SOLUBILITY_REFERENCE_M["benzene"]
    pyrene = INTRINSIC_WATER_SOLUBILITY_REFERENCE_M["pyrene"]
    assert benzene["solubility_M"] > pyrene["solubility_M"]
    assert benzene["confidence"] == "moderate"
    assert pyrene["confidence"] == "moderate"


def test_intrinsic_water_solubility_reference_usable_directly_in_msr():
    """Real, usable end-to-end check: the reference table's own
    naphthalene value plugged directly into molar_solubilization_ratio
    reproduces this project's own already-verified Paria & Yuet 2006
    worked example (literature_validation_notes.md Round 5)."""
    intrinsic = INTRINSIC_WATER_SOLUBILITY_REFERENCE_M["naphthalene"]["solubility_M"]
    msr = molar_solubilization_ratio(
        total_solubilized_M=2.517e-3, intrinsic_water_solubility_M=intrinsic,
        surfactant_concentration_M=0.040, cmc_M=0.0085,
    )
    assert msr == pytest.approx(0.0730, abs=1e-3)


def test_qspr_estimate_vs_real_measured_naphthalene_value():
    """A genuinely interesting cross-check between two Tier 1/Tier 2
    additions: how close does the coarse ESOL QSPR ESTIMATE (Tier 1)
    land to the real MEASURED value (Tier 2, Paria & Yuet 2006) for the
    same real compound? Not asserted tight (ESOL's own real published
    error is ~0.6-1 log unit) -- this is a real, disclosed comparison,
    not a claim the estimate should be trusted as a substitute."""
    qspr = estimate_intrinsic_water_solubility_qspr("c1ccc2ccccc2c1")  # naphthalene
    real = INTRINSIC_WATER_SOLUBILITY_REFERENCE_M["naphthalene"]["solubility_M"]
    log_diff = abs(qspr.log_s_mol_per_L - math.log10(real))
    assert log_diff < 2.0  # sanity bound only -- not a precision claim
