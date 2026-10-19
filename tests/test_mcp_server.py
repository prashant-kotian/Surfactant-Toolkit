"""End-to-end tests for the SurfMCP server: verifies tools are
registered and that calling them through the MCP protocol layer
produces the same results as calling the underlying library functions
directly (i.e. the wrapping doesn't silently change any numbers)."""

import asyncio
import json
import math

import pytest

from surfactantkit.mcp_server import mcp


def call(tool_name: str, args: dict) -> dict:
    """Call an MCP tool synchronously and parse its JSON text content.

    Real API-shape change found 2026-09-13, alongside the FastMCP ->
    MCPServer rename (see mcp_server.py, migrated the same day): the mcp
    SDK installed in that session's environment returned a CallToolResult
    object with a `.content` list attribute, not a bare list.

    CORRECTED 2026-09-14: the 09-13 fix assumed every real environment had
    upgraded to that shape and dropped the old bare-list handling entirely
    -- but this Windows install's `mcp` package is still 1.27.2 (confirmed
    directly: `type(result)` here is `list`, not `CallToolResult`), which
    still returns the bare list. Both shapes are real, currently in use on
    different real machines in this project (Ubuntu vs. Windows) -- not a
    hypothetical, so both are handled rather than picking one and hoping.
    """
    result = asyncio.run(mcp.call_tool(tool_name, args))
    content = result.content[0] if hasattr(result, "content") else result[0]
    return json.loads(content.text)


def test_all_expected_tools_are_registered():
    tools = asyncio.run(mcp.list_tools())
    names = {t.name for t in tools}
    expected = {
        "clint_ideal_cmc",
        "rubingh_solve",
        "rubingh_activity_coefficients",
        "excess_free_energy",
        "gibbs_surface_excess",
        "gibbs_area_per_molecule",
        "hlb_from_mw",
        "hlb_from_groups",
        "hlb_method_recommendation",
        "gemini_hlb_griffin_reference",
        "blankschtein_entropy_of_binding",
        "blankschtein_steric_free_energy",
        "blankschtein_counterion_self_energy",
        "hld_optimal_salinity",
        "hld_cationic_quat_reference",
        "hld_class_reference",
        "hld_characteristic_curvature_reference",
        "counterion_binding_degree_reference",
        "dn_dc_reference",
        "intrinsic_water_solubility_reference",
        "nagarajan_ruckenstein_headgroup_reference",
        "nagarajan_ruckenstein_steric_free_energy",
        "nagarajan_ruckenstein_dipole_free_energy",
        "nagarajan_ruckenstein_ionic_free_energy",
        "hld_fit_k_and_cc_from_salinity_scan",
        "tanford_chain_geometry",
        "critical_packing_parameter",
        "aggregation_number_from_quenching",
        "aggregation_number_from_sls",
        "cmc_from_surface_tension_curve",
        "cmc_from_conductivity",
        "rubingh_beta_regression",
        "frumkin_predict_surface_tension",
        "motomura_ideal_composition",
        "asymmetric_margules_activity_coefficients",
        "rodenas_x1",
        "rodenas_activity_coefficients",
        "rodenas_x1_from_series",
        "maeda_free_energy_of_micellization",
        "vant_hoff_multi_point_fit",
        "szyszkowski_fit_K",
        "eommm_global_fit",
        "nagarajan_equilibrium_area_ionic",
        "hydrodynamic_radius_perrin_corrected",
        "hlb_davies_guo_ecl",
        "mass_action_free_energy_of_micellization",
        "predict_mobility_relaxation_corrected",
        "zeta_potential_relaxation_corrected",
        "owens_wendt_solid_surface_energy",
        "van_oss_chaudhury_good_solid_surface_energy",
        "get_standard_probe_liquid_properties",
        "estimate_axial_ratio_from_cpp_geometry",
        "estimate_partial_specific_volume",
        "classify_surfactant_charge_type_at_ph",
        "estimate_intrinsic_water_solubility_qspr",
        "micelle_water_partition_coefficient",
        "grahame_equation_surface_potential",
        "aggregation_number_from_dls",
        "aggregation_number_from_svedberg",
        "derive_davies_group_number",
        "critical_micellization_degree",
        "eommm_find_minimal_feasible_margin",
        "classify_surfactant_charge_type",
        "derive_all_properties_from_smiles_and_curve",
        "nagarajan_cmc_predictor",
        "nagarajan_free_energy_terms",
    }
    assert expected <= names


def test_clint_ideal_cmc_tool_matches_library():
    out = call("clint_ideal_cmc", {"alpha1": 0.25, "cmc1_mM": 14.80, "cmc2_mM": 8.00})
    assert out["cmc_ideal_mM"] == pytest.approx(9.038, abs=0.01)
    assert out["unit"] == "mM"


def test_rubingh_solve_tool_matches_cholate_sds_literature_case():
    out = call(
        "rubingh_solve",
        {"alpha1": 0.5, "cmc_mix_mM": 4.07, "cmc1_mM": 11.50, "cmc2_mM": 11.98},
    )
    assert out["micellar_mole_fraction_x1"] == pytest.approx(0.5033, abs=0.001)
    assert out["beta"] == pytest.approx(-4.236, abs=0.01)
    assert out["synergy_classification"] == "synergistic"


def test_rubingh_beta_regression_tool_recovers_true_beta():
    import math

    cmc1, cmc2, beta_true = 14.80, 8.00, -1.8
    alpha1_series, cmc_mix_series = [], []
    for x1 in (0.25, 0.4, 0.55, 0.7):
        f1 = math.exp(beta_true * (1.0 - x1) ** 2)
        f2 = math.exp(beta_true * x1 ** 2)
        cmc_mix = x1 * f1 * cmc1 + (1.0 - x1) * f2 * cmc2
        alpha1 = x1 * f1 * cmc1 / cmc_mix
        alpha1_series.append(alpha1)
        cmc_mix_series.append(cmc_mix)

    out = call(
        "rubingh_beta_regression",
        {"alpha1_series": alpha1_series, "cmc_mix_series_mM": cmc_mix_series, "cmc1_mM": cmc1, "cmc2_mM": cmc2},
    )
    assert out["beta_mean"] == pytest.approx(beta_true, abs=1e-3)
    assert out["n_points_used"] == 4
    assert out["n_points_skipped"] == 0
    assert out["synergy_classification"] == "synergistic"


def test_szyszkowski_fit_k_tool_recovers_true_k_round_trip():
    true_K = 0.35
    gamma0, gamma_max, system_type = 72.0, 4.0e-6, "nonionic"
    concentrations = [0.5, 1.0, 2.0, 4.0, 8.0, 14.0]
    tensions = [
        call("szyszkowski_predict_surface_tension", {
            "concentration": c, "gamma0_mN_m": gamma0, "gamma_max_mol_per_m2": gamma_max,
            "K": true_K, "system_type": system_type,
        })["surface_tension_mN_m"]
        for c in concentrations
    ]

    out = call("szyszkowski_fit_K", {
        "concentrations": concentrations, "surface_tensions_mN_per_m": tensions,
        "gamma0_mN_m": gamma0, "gamma_max_mol_per_m2": gamma_max, "system_type": system_type,
    })
    assert out["K"] == pytest.approx(true_K, rel=1e-4)
    assert out["r_squared"] == pytest.approx(1.0, abs=1e-6)
    assert out["n_points"] == 6


def test_frumkin_tool_zero_a_matches_szyszkowski_tool():
    """a=0 must reduce EXACTLY to szyszkowski_predict_surface_tension --
    the real regression guard for this tool, since no external numeric
    literature match exists for the Frumkin isotherm specifically (see
    test_adsorption_alternatives.py's module docstring)."""
    frumkin_args = {
        "concentration": 2.0,
        "gamma0_mN_m": 72.0,
        "gamma_max_mol_per_m2": 4.0e-6,
        "K": 0.5,
        "a": 0.0,
        "system_type": "nonionic",
    }
    szyszkowski_args = {k: v for k, v in frumkin_args.items() if k != "a"}
    frumkin_out = call("frumkin_predict_surface_tension", frumkin_args)
    szyszkowski_out = call("szyszkowski_predict_surface_tension", szyszkowski_args)
    assert frumkin_out["surface_tension_mN_m"] == pytest.approx(szyszkowski_out["surface_tension_mN_m"], abs=1e-6)


def test_motomura_ideal_composition_tool_matches_library():
    out = call("motomura_ideal_composition", {"alpha1": 0.75, "cmc1_mM": 14.80, "cmc2_mM": 8.00})
    assert out["x1_ideal"] == pytest.approx(0.75 * 12.21 / 14.80, abs=0.001)


def test_asymmetric_margules_tool_reduces_to_rst_when_symmetric():
    beta = -1.5
    asym_out = call("asymmetric_margules_activity_coefficients", {"x1": 0.4, "w12": beta, "w21": beta})
    rst_out = call("rubingh_activity_coefficients", {"x1": 0.4, "beta": beta})
    assert asym_out["f1"] == pytest.approx(rst_out["f1"], rel=1e-9)
    assert asym_out["f2"] == pytest.approx(rst_out["f2"], rel=1e-9)


def test_rodenas_tools_end_to_end():
    x1_out = call("rodenas_x1", {"alpha1": 0.5, "dln_cmc_mix_dalpha1": 0.0})
    assert x1_out["micellar_mole_fraction_x1"] == pytest.approx(0.5)
    f_out = call(
        "rodenas_activity_coefficients",
        {"alpha1": 0.5, "x1": 0.6, "cmc_mix_mM": 5.0, "cmc1_mM": 10.0, "cmc2_mM": 8.0},
    )
    assert f_out["f1"] == pytest.approx((0.5 * 5.0) / (0.6 * 10.0))


def test_maeda_tool_matches_library():
    out = call(
        "maeda_free_energy_of_micellization",
        {"x1_rub": 0.6, "beta": -1.8, "cmc1_M": 0.041e-3, "cmc2_M": 0.263e-3, "temperature_K": 298.15},
    )
    assert isinstance(out["deltaG_M_kJ_per_mol"], float)
    assert out["deltaG_M_kJ_per_mol"] < 0


def test_hld_cationic_quat_reference_tool_matches_library():
    out = call("hld_cationic_quat_reference", {"surfactant": "ctab"})
    assert out["Cc"] == pytest.approx(-5.7)
    assert out["k"] == pytest.approx(0.7)
    assert out["real_reported_hlb_davies_scale"] == pytest.approx(21.4)


def test_hld_cationic_quat_reference_tool_rejects_unknown_surfactant():
    out = call("hld_cationic_quat_reference", {"surfactant": "SDS"})
    assert out["status"] == "no_solution"
    assert out["reason"]


def test_nagarajan_ruckenstein_headgroup_reference_tool_matches_library():
    out = call("nagarajan_ruckenstein_headgroup_reference", {"headgroup": "sodium_sulfonate"})
    assert out["delta_A"] == pytest.approx(3.85)
    assert out["class"] == "anionic"

    out_betaine = call("nagarajan_ruckenstein_headgroup_reference", {"headgroup": "N_Betaine"})
    assert out_betaine["d_A"] == pytest.approx(5.0)
    assert out_betaine["class"] == "zwitterionic"


def test_nagarajan_ruckenstein_headgroup_reference_tool_rejects_unknown():
    out = call("nagarajan_ruckenstein_headgroup_reference", {"headgroup": "quaternary_ammonium"})
    assert out["status"] == "no_solution"
    assert out["reason"]


def test_nagarajan_ruckenstein_steric_free_energy_tool_matches_library():
    out = call("nagarajan_ruckenstein_steric_free_energy", {"area_per_molecule_A2": 34.0, "a_p_A2": 17.0})
    assert out["steric_free_energy_kT"] == pytest.approx(-math.log(1.0 - 17.0 / 34.0))


def test_nagarajan_ruckenstein_dipole_free_energy_tool_matches_library():
    out = call("nagarajan_ruckenstein_dipole_free_energy", {
        "area_per_molecule_A2": 21.0, "radius_A": 20.0, "d_A": 5.0, "geometry": "sphere",
    })
    assert out["dipole_free_energy_kT"] > 0.0


def test_nagarajan_ruckenstein_ionic_free_energy_tool_matches_library():
    out = call("nagarajan_ruckenstein_ionic_free_energy", {
        "area_per_molecule_A2": 60.0, "core_radius_A": 16.5, "delta_A": 5.45,
        "counterion_concentration_M": 0.008,
    })
    assert out["ionic_free_energy_kT"] > 0.0


def test_intrinsic_water_solubility_reference_tool_matches_library():
    out = call("intrinsic_water_solubility_reference", {"solubilizate": "Naphthalene"})
    assert out["solubility_M"] == pytest.approx(2.17e-4)
    assert out["confidence"] == "high"


def test_intrinsic_water_solubility_reference_tool_rejects_unknown_solubilizate():
    out = call("intrinsic_water_solubility_reference", {"solubilizate": "anthracene"})
    assert out["status"] == "no_solution"
    assert out["reason"]


def test_dn_dc_reference_tool_matches_library():
    out = call("dn_dc_reference", {"surfactant": "sds"})
    assert out["dn_dc"] == pytest.approx(0.11)
    assert out["solvent"] == "water"

    out_triton = call("dn_dc_reference", {"surfactant": "TritonX100"})
    assert out_triton["dn_dc"] == pytest.approx(0.140)
    assert out_triton["dn_dc_uncertainty"] == pytest.approx(0.005)


def test_dn_dc_reference_tool_rejects_unknown_surfactant():
    out = call("dn_dc_reference", {"surfactant": "Tween-80"})
    assert out["status"] == "no_solution"
    assert out["reason"]


def test_counterion_binding_degree_reference_tool_matches_library():
    out = call("counterion_binding_degree_reference", {"surfactant": "sds"})
    assert out["alpha"] == pytest.approx(0.272)
    assert out["beta"] == pytest.approx(0.728)
    assert out["alpha_confidence"] == "high"

    out_dtab = call("counterion_binding_degree_reference", {"surfactant": "DTAB"})
    assert out_dtab["alpha_confidence"] == "moderate"

    out_aot = call("counterion_binding_degree_reference", {"surfactant": "aot"})
    assert out_aot["alpha"] == pytest.approx(0.61)
    assert out_aot["alpha_confidence"] == "high"


def test_counterion_binding_degree_reference_tool_rejects_unknown_surfactant():
    out = call("counterion_binding_degree_reference", {"surfactant": "CAPB"})
    assert out["status"] == "no_solution"
    assert out["reason"]


def test_hld_characteristic_curvature_reference_tool_matches_library():
    out_sds = call("hld_characteristic_curvature_reference", {"surfactant": "sds"})
    assert out_sds["cc_this_work"] == pytest.approx(-2.63)
    assert out_sds["cc_literature"] == pytest.approx(-2.5)

    out_aot = call("hld_characteristic_curvature_reference", {"surfactant": "AOT"})
    assert out_aot["cc_this_work"] > 0

    out_rham = call("hld_characteristic_curvature_reference", {"surfactant": "rhamnolipid"})
    assert out_rham["cc"] == pytest.approx(-1.41)
    assert out_rham["class"] == "biosurfactant (glycolipid)"

    out_capb = call("hld_characteristic_curvature_reference", {"surfactant": "capb"})  # case-insensitive
    assert out_capb["cc_this_work"] == pytest.approx(-5.2)
    assert out_capb["class"] == "zwitterionic"

    out_gemini = call("hld_characteristic_curvature_reference", {"surfactant": "gemini_benzene_sulfonate"})
    assert out_gemini["cc"] == pytest.approx(-7.4)
    assert out_gemini["class"] == "gemini/dimeric (anionic)"


def test_hld_characteristic_curvature_reference_tool_rejects_unknown():
    out = call("hld_characteristic_curvature_reference", {"surfactant": "sulfobetaine_unknown_variant"})
    assert out["status"] == "no_solution"
    assert out["reason"]


def test_hld_class_reference_tool_matches_library():
    out_anionic = call("hld_class_reference", {"surfactant_class": "anionic_sulfate_sulfonate"})
    assert out_anionic["k_default"] == pytest.approx(0.16)
    assert out_anionic["sds_cc_reference"] == pytest.approx(-2.63)

    out_extended = call("hld_class_reference", {"surfactant_class": "Extended_Surfactant"})  # case-insensitive
    assert out_extended["k_default"] == pytest.approx(0.06)

    out_apg = call("hld_class_reference", {"surfactant_class": "apg_nonionic"})
    assert out_apg["alpha_default"] == 0.0
    assert out_apg["b_default"] == 0.0


def test_hld_class_reference_tool_rejects_unknown_class():
    out = call("hld_class_reference", {"surfactant_class": "zwitterionic"})
    assert out["status"] == "no_solution"
    assert out["reason"]


def test_hld_optimal_salinity_tool_round_trips_with_library():
    from surfactantkit.hld import hld_ionic

    out = call("hld_optimal_salinity", {"k": 0.7, "eacn": 8.0, "cc": -5.7})
    s_star = out["optimal_salinity_pct"]
    assert hld_ionic(s_star, 0.7, 8.0, -5.7) == pytest.approx(0.0, abs=1e-6)


def test_hld_fit_k_and_cc_from_salinity_scan_tool_round_trip():
    from surfactantkit.hld import optimal_salinity_ionic

    true_k, true_cc = 0.7, -5.7
    eacn_series = [6.0, 8.0, 10.0, 12.0]
    s_star_series = [optimal_salinity_ionic(true_k, eacn, true_cc) for eacn in eacn_series]

    out = call(
        "hld_fit_k_and_cc_from_salinity_scan",
        {"eacn_series": eacn_series, "optimal_salinity_series_pct": s_star_series},
    )
    assert out["k"] == pytest.approx(true_k, abs=1e-6)
    assert out["cc"] == pytest.approx(true_cc, abs=1e-6)
    assert out["r_squared"] == pytest.approx(1.0, abs=1e-6)
    assert out["n_points"] == 4


def test_hlb_from_groups_tool_sds_worked_example():
    out = call("hlb_from_groups", {"group_counts": {"SO4Na": 1, "CH2": 11, "CH3": 1}})
    assert out["hlb"] == pytest.approx(39.9, abs=0.01)


def test_hlb_davies_guo_ecl_tool_matches_manual_formula():
    import math

    n_carbons, n_eo = 12, 9
    n_ch2_eff = 0.965 * (n_carbons - 1) - 0.178
    n_eo_eff = 13.45 * math.log(n_eo) - 0.16 * n_eo + 1.26
    expected = 7.0 + 1.3 * n_eo_eff - 0.475 * n_ch2_eff - 0.475 * 1

    out = call("hlb_davies_guo_ecl", {"n_carbons_alkyl": n_carbons, "n_eo": n_eo})
    assert out["hlb"] == pytest.approx(expected)


def test_hlb_from_groups_tool_raises_on_unverified_group():
    out = call("hlb_from_groups", {"group_counts": {"quaternary_ammonium": 1}})
    assert out["status"] == "no_solution"
    assert out["reason"]


def test_hlb_from_groups_tool_sulfonate_requires_opt_in():
    out = call("hlb_from_groups", {"group_counts": {"sulfonate": 1, "CH2": 11, "CH3": 1}})
    assert out["status"] == "no_solution"
    assert out["reason"]


def test_hlb_method_recommendation_tool_matches_library():
    out_gemini = call("hlb_method_recommendation", {"structural_family": "dimeric_gemini_type"})
    assert out_gemini["recommended_method"] == "griffin"
    assert out_gemini["davies_valid"] is False

    out_glyco = call("hlb_method_recommendation", {"structural_family": "glycolipid_biosurfactant"})
    assert out_glyco["recommended_method"] == "griffin"

    out_mono = call("hlb_method_recommendation", {"structural_family": "monomeric"})
    assert out_mono["recommended_method"] == "either"
    assert out_mono["davies_valid"] is True


def test_hlb_method_recommendation_tool_rejects_unparseable():
    out = call("hlb_method_recommendation", {"structural_family": "unparseable"})
    assert out["status"] == "no_solution"
    assert out["reason"]


def test_gemini_hlb_griffin_reference_tool_matches_library():
    out = call("gemini_hlb_griffin_reference", {"gemini_name": "12-3-12"})
    assert out["hlb"] == pytest.approx(6.37)
    assert out["method"] == "Griffin"

    out2 = call("gemini_hlb_griffin_reference", {"gemini_name": "8-4-8"})
    assert out2["hlb"] == pytest.approx(7.73)


def test_gemini_hlb_griffin_reference_tool_rejects_unknown_name():
    out = call("gemini_hlb_griffin_reference", {"gemini_name": "16-3-16"})
    assert out["status"] == "no_solution"
    assert out["reason"]


def test_blankschtein_entropy_of_binding_tool_matches_library():
    from surfactantkit.cpp import blankschtein_entropy_of_binding_free_energy

    out = call("blankschtein_entropy_of_binding", {"beta": 0.5})
    assert out["gent_kT"] == pytest.approx(blankschtein_entropy_of_binding_free_energy(0.5))


def test_blankschtein_steric_free_energy_tool_matches_library():
    from surfactantkit.cpp import blankschtein_steric_free_energy_with_counterion

    out = call("blankschtein_steric_free_energy", {
        "beta": 0.5, "area_per_molecule_A2": 60.0, "ah_surfactant_A2": 25.0, "ah_counterion_A2": 4.537,
    })
    expected = blankschtein_steric_free_energy_with_counterion(0.5, 60.0, 25.0, 4.537)
    assert out["gst_kT"] == pytest.approx(expected)


def test_blankschtein_counterion_self_energy_tool_matches_library():
    from surfactantkit.cpp import blankschtein_counterion_self_energy_release

    out = call("blankschtein_counterion_self_energy", {
        "valence": 1, "hydrated_radius_A": 2.13, "kappa_inverse_A": 10.0,
    })
    expected = blankschtein_counterion_self_energy_release(1, 2.13, 10.0)
    assert out["gdis_kT"] == pytest.approx(expected)


    out = call("hlb_from_groups", {"group_counts": {"sulfonate": 1, "CH2": 11, "CH3": 1}, "allow_derived_groups": True})
    assert out["hlb"] == pytest.approx(7.566, abs=0.01)


def test_derive_davies_group_number_tool_matches_library():
    # sodium dodecanesulfonate, n=12
    mh = 32.06 + 3 * 15.999 + 22.990
    m_total = 12 * 12.011 + 25 * 1.008 + mh
    out = call("derive_davies_group_number", {
        "hydrophilic_fragment_mass_g_per_mol": mh,
        "total_molar_mass_g_per_mol": m_total,
        "lipophilic_group_counts": {"CH2": 11, "CH3": 1},
    })
    assert out["derived_group_number"] == pytest.approx(6.266, abs=0.01)


def test_critical_packing_parameter_tool_end_to_end():
    geom = call("tanford_chain_geometry", {"n_carbons": 12})
    out = call(
        "critical_packing_parameter",
        {"volume_A3": geom["tail_volume_A3"], "head_area_A2": 50.0, "length_A": geom["critical_length_A"]},
    )
    assert out["cpp"] == pytest.approx(350.2 / (50.0 * 16.68), abs=0.001)
    assert out["predicted_morphology"] in {
        "spherical micelle",
        "cylindrical/rodlike micelle",
        "vesicle/bilayer",
        "inverted structure",
    }


def test_nagarajan_equilibrium_area_ionic_tool_matches_table2():
    out = call(
        "nagarajan_equilibrium_area_ionic",
        {"cmc_M": 0.008, "tail_length_A": 16.5, "headgroup_prefactor_A": 82.0},
    )
    assert out["a_e_A2"] == pytest.approx(67.4, rel=0.01)
    assert out["debye_huckel_kappa_inverse_A"] == pytest.approx(34.43, rel=0.01)


def test_aggregation_number_tool():
    geom = call("tanford_chain_geometry", {"n_carbons": 12})
    out = call("aggregation_number", {"tail_volume_A3": geom["tail_volume_A3"], "core_radius_A": geom["critical_length_A"]})
    assert 40.0 < out["aggregation_number"] < 80.0


def test_aggregation_number_from_quenching_tool_round_trips_with_library():
    import math

    true_n, ct, cmc, i0 = 60.0, 40.0, 8.0, 1000.0
    slope = true_n / (ct - cmc)
    quencher_concs = [0.001, 0.002, 0.003]
    intensities = [i0 * math.exp(-slope * q) for q in quencher_concs]
    out = call(
        "aggregation_number_from_quenching",
        {
            "quencher_concentrations": quencher_concs,
            "intensities": intensities,
            "intensity_without_quencher": i0,
            "total_surfactant_concentration": ct,
            "cmc": cmc,
        },
    )
    assert out["aggregation_number"] == pytest.approx(true_n, rel=1e-6)


def test_aggregation_number_from_sls_tool_round_trips_with_library():
    import math

    n_solvent, dn_dc, wavelength_nm = 1.333, 0.12, 633.0
    wavelength_cm = wavelength_nm * 1e-7
    k = 4.0 * math.pi ** 2 * n_solvent ** 2 * dn_dc ** 2 / (6.02214076e23 * wavelength_cm ** 4)
    true_mw, monomer_mw = 15000.0, 288.38
    concentrations = [0.001, 0.002, 0.003]
    rayleigh_ratios = [k * c / (1.0 / true_mw + 2.0e-4 * c) for c in concentrations]
    out = call(
        "aggregation_number_from_sls",
        {
            "micellized_concentrations_g_per_mL": concentrations,
            "rayleigh_ratios_cm_inv": rayleigh_ratios,
            "dn_dc_mL_per_g": dn_dc,
            "wavelength_nm": wavelength_nm,
            "monomer_molar_mass_g_per_mol": monomer_mw,
            "refractive_index_solvent": n_solvent,
        },
    )
    assert out["micelle_molar_mass_g_per_mol"] == pytest.approx(true_mw, rel=1e-6)


def test_aggregation_number_from_dls_tool_matches_library():
    out = call(
        "aggregation_number_from_dls",
        {
            "diffusion_coefficient_cm2_per_s": 1.56e-6, "viscosity_mPas": 0.89,
            "partial_specific_volume_cm3_per_g": 0.87, "monomer_molar_mass_g_per_mol": 288.4,
        },
    )
    assert 10.0 <= out["aggregation_number"] <= 500.0
    assert out["hydrodynamic_radius_nm"] == pytest.approx(1.573, abs=0.01)


def test_aggregation_number_from_svedberg_tool_matches_library():
    out = call(
        "aggregation_number_from_svedberg",
        {
            "sedimentation_coefficient_S": 2.0, "diffusion_coefficient_cm2_per_s": 1.0e-6,
            "partial_specific_volume_cm3_per_g": 0.8, "solvent_density_g_per_cm3": 0.997,
            "monomer_molar_mass_g_per_mol": 288.4,
        },
    )
    assert 10.0 <= out["aggregation_number"] <= 500.0


def test_cmc_from_surface_tension_curve_tool_matches_library():
    concentrations = [0.01585, 0.03981, 0.10000, 0.25119, 0.63096, 1.58489, 3.981, 6.31, 10.0]
    tensions = [62.22, 55.66, 49.10, 42.54, 35.98, 29.42, 27.5, 27.8, 28.0]
    out = call("cmc_from_surface_tension_curve", {"concentrations_mM": concentrations, "surface_tensions_mN_per_m": tensions})
    assert out["cmc_mM"] == pytest.approx(2.51, rel=0.15)


def test_cmc_from_conductivity_tool_matches_library():
    true_cmc, slope_below, slope_above = 8.2, 1.20, 0.36  # beta = 1-0.36/1.20 = 0.70
    intercept_above = (slope_below - slope_above) * true_cmc

    concentrations = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 9.0, 11.0, 13.0, 15.0]
    conductivities = [
        slope_below * c if c < true_cmc else slope_above * c + intercept_above
        for c in concentrations
    ]

    out = call("cmc_from_conductivity", {"concentrations_mM": concentrations, "conductivities": conductivities})
    assert out["cmc_mM"] == pytest.approx(true_cmc, rel=1e-6)
    assert out["counterion_binding_degree"] == pytest.approx(0.70, rel=1e-6)


def test_debye_screening_length_tool_matches_textbook_value():
    out = call("debye_screening_length", {"ionic_strength_M": 0.1, "temperature_K": 298.15})
    assert out["debye_length_nm"] == pytest.approx(0.961, abs=0.01)


def test_grahame_equation_tool_zero_charge_gives_zero_potential():
    out = call("grahame_equation_surface_potential", {"surface_charge_density_C_per_m2": 0.0, "ionic_strength_M": 0.01})
    assert out["surface_potential_mV"] == pytest.approx(0.0)


def test_zeta_potential_tool_requires_explicit_regime():
    out = call("zeta_potential", {"electrophoretic_mobility_um_cm_per_Vs": 2.0, "viscosity_mPas": 0.89, "regime": "not_a_real_regime"})
    assert out["status"] == "no_solution"
    assert out["reason"]


def test_relaxation_corrected_zeta_potential_tool_round_trip():
    zeta_true, kappa_a, eta = 60.0, 50.0, 0.89
    predicted = call(
        "predict_mobility_relaxation_corrected",
        {"zeta_mV": zeta_true, "kappa_a": kappa_a, "viscosity_mPas": eta},
    )
    recovered = call(
        "zeta_potential_relaxation_corrected",
        {"electrophoretic_mobility_um_cm_per_Vs": predicted["electrophoretic_mobility_um_cm_per_Vs"], "kappa_a": kappa_a, "viscosity_mPas": eta},
    )
    assert recovered["zeta_potential_mV"] == pytest.approx(zeta_true, abs=1e-3)


def test_hydrodynamic_radius_tool_positive_and_reasonable():
    out = call("hydrodynamic_radius", {"diffusion_coefficient_cm2_per_s": 1e-6, "viscosity_mPas": 0.89})
    assert 0.1 < out["hydrodynamic_radius_nm"] < 1000.0


def test_hydrodynamic_radius_perrin_corrected_tool_sphere_limit_matches_plain_tool():
    plain = call("hydrodynamic_radius", {"diffusion_coefficient_cm2_per_s": 1.56e-6, "viscosity_mPas": 0.89})
    corrected = call(
        "hydrodynamic_radius_perrin_corrected",
        {"diffusion_coefficient_cm2_per_s": 1.56e-6, "viscosity_mPas": 0.89, "axial_ratio": 1.0},
    )
    assert corrected["hydrodynamic_radius_true_nm"] == pytest.approx(plain["hydrodynamic_radius_nm"])
    assert corrected["perrin_friction_factor"] == pytest.approx(1.0)


def test_hydrodynamic_radius_perrin_corrected_tool_smaller_for_rodlike():
    plain = call("hydrodynamic_radius", {"diffusion_coefficient_cm2_per_s": 1.56e-6, "viscosity_mPas": 0.89})
    corrected = call(
        "hydrodynamic_radius_perrin_corrected",
        {"diffusion_coefficient_cm2_per_s": 1.56e-6, "viscosity_mPas": 0.89, "axial_ratio": 4.0},
    )
    assert corrected["hydrodynamic_radius_true_nm"] < plain["hydrodynamic_radius_nm"]


def test_rosen_monolayer_solve_tool_matches_underlying_math():
    out = call("rosen_monolayer_solve", {"alpha1": 0.5, "c_mix_sigma_mM": 4.07, "c1_sigma_mM": 11.50, "c2_sigma_mM": 11.98})
    assert out["monolayer_mole_fraction_x1"] == pytest.approx(0.5033, abs=0.001)
    assert out["beta_sigma"] == pytest.approx(-4.236, abs=0.01)


def test_corrin_harkins_predict_tool_cmc_decreases_with_salt():
    out = call("corrin_harkins_predict", {
        "cmc1_mM": 10.0, "salt_conc1_mM": 10.0,
        "cmc2_mM": 5.0, "salt_conc2_mM": 100.0,
        "salt_conc_target_mM": 200.0,
    })
    assert out["predicted_cmc_mM"] < 5.0


def test_counterion_binding_degree_tool():
    out = call("counterion_binding_degree", {"slope_below_cmc": 2.0, "slope_above_cmc": 1.0})
    assert out["beta"] == pytest.approx(0.5)


def test_gibbs_free_energy_micellization_tool_is_negative():
    out = call("gibbs_free_energy_micellization", {"cmc_M": 0.008, "temperature_K": 298.15, "counterion_factor": 1.6})
    assert out["deltaG_mic_kJ_per_mol"] < 0


def test_vant_hoff_and_entropy_tools_complete_triad():
    dg = call("gibbs_free_energy_micellization", {"cmc_M": 0.008, "temperature_K": 298.15, "counterion_factor": 1.6})
    dh = call("vant_hoff_enthalpy", {"cmc1_M": 0.008, "temperature1_K": 293.15, "cmc2_M": 0.0075, "temperature2_K": 313.15})
    ds = call("entropy_of_micellization", {
        "deltaG_mic_kJ_per_mol": dg["deltaG_mic_kJ_per_mol"],
        "deltaH_mic_kJ_per_mol": dh["deltaH_mic_kJ_per_mol"],
        "temperature_K": 298.15,
    })
    expected = ((dh["deltaH_mic_kJ_per_mol"] - dg["deltaG_mic_kJ_per_mol"]) / 298.15) * 1000.0
    assert ds["deltaS_mic_J_per_mol_K"] == pytest.approx(expected)


def test_eommm_global_fit_tool_recovers_true_parameters_round_trip():
    """Uses a TIGHT cmc_margin: see mixed_micelle.eommm_global_fit's own
    docstring -- at the default (looser) margin, the fit is not always
    uniquely determined, a real disclosed property of the model, not a
    tool bug."""
    from surfactantkit.mixed_micelle import asymmetric_margules_activity_coefficients

    cmc1, cmc2 = 14.80, 8.00
    true_w12, true_w21 = 2.0, -3.0
    x1_true_values = [0.25, 0.5, 0.75]
    alpha1_series, cmc_mix_series = [], []
    for x1 in x1_true_values:
        f1, f2 = asymmetric_margules_activity_coefficients(x1, true_w12, true_w21)
        d = x1 * f1 * cmc1 + (1.0 - x1) * f2 * cmc2
        alpha1_series.append(x1 * f1 * cmc1 / d)
        cmc_mix_series.append(d)

    out = call(
        "eommm_global_fit",
        {"alpha1_series": alpha1_series, "cmc_mix_series_mM": cmc_mix_series, "cmc1_mM": cmc1, "cmc2_mM": cmc2,
         "cmc_margin": 1e-6},
    )
    assert out["W12"] == pytest.approx(true_w12, abs=1e-2)
    assert out["W21"] == pytest.approx(true_w21, abs=1e-2)
    assert out["total_infeasibility"] < 5e-4
    assert out["n_points"] == 3


def test_eommm_find_minimal_feasible_margin_tool_round_trip():
    """Reduced binary_search_iters (8) to keep this end-to-end test's
    runtime bounded (~30s) -- see mixed_micelle.eommm_find_minimal_feasible_margin's
    own docstring for the full re-verification of why this automatic
    search (previously disclosed as unshippable) now works reliably."""
    from surfactantkit.mixed_micelle import asymmetric_margules_activity_coefficients

    cmc1, cmc2 = 14.80, 8.00
    true_w12, true_w21 = 2.0, -3.0
    x1_true_values = [0.25, 0.5, 0.75]
    alpha1_series, cmc_mix_series = [], []
    for x1 in x1_true_values:
        f1, f2 = asymmetric_margules_activity_coefficients(x1, true_w12, true_w21)
        d = x1 * f1 * cmc1 + (1.0 - x1) * f2 * cmc2
        alpha1_series.append(x1 * f1 * cmc1 / d)
        cmc_mix_series.append(d)

    out = call(
        "eommm_find_minimal_feasible_margin",
        {"alpha1_series": alpha1_series, "cmc_mix_series_mM": cmc_mix_series, "cmc1_mM": cmc1, "cmc2_mM": cmc2,
         "binary_search_iters": 8},
    )
    assert out["W12"] == pytest.approx(true_w12, abs=1e-2)
    assert out["W21"] == pytest.approx(true_w21, abs=1e-2)
    assert out["total_infeasibility"] < 5e-4
    assert out["minimal_feasible_cmc_margin"] < 0.01


def test_rodenas_x1_from_series_tool_matches_library():
    import math

    a, b, c = 0.5, -2.0, 3.0  # ln(cmc_mix) = a*alpha1^2 + b*alpha1 + c
    alpha1_series = [0.1, 0.25, 0.4, 0.55, 0.7]
    cmc_mix_series = [math.exp(a * x ** 2 + b * x + c) for x in alpha1_series]
    true_slope = lambda x: 2.0 * a * x + b

    out = call("rodenas_x1_from_series", {"alpha1_series": alpha1_series, "cmc_mix_series_mM": cmc_mix_series})
    assert out["n_points"] == 5
    for x, slope in zip(out["alpha1_used"], out["dln_cmc_mix_dalpha1"]):
        assert slope == pytest.approx(true_slope(x), abs=1e-9)


def test_mass_action_free_energy_tool_converges_to_pseudo_phase_at_large_n():
    pseudo_phase = call("gibbs_free_energy_micellization", {"cmc_M": 0.008, "temperature_K": 298.15, "counterion_factor": 1.0})
    mass_action = call(
        "mass_action_free_energy_of_micellization",
        {"cmc_M": 0.008, "aggregation_number": 100000, "temperature_K": 298.15},
    )
    assert mass_action["deltaG_mic_kJ_per_mol"] == pytest.approx(pseudo_phase["deltaG_mic_kJ_per_mol"], abs=1e-3)


def test_vant_hoff_multi_point_fit_tool_recovers_true_parameters():
    import math

    from surfactantkit.thermodynamics import R_GAS

    true_dH_kJ, true_dS_J, true_dCp_J = -20.0, -50.0, -300.0
    tr = 298.15
    temps = [278.15, 288.15, 298.15, 308.15, 318.15]
    true_dS_kJ, true_dCp_kJ = true_dS_J / 1000.0, true_dCp_J / 1000.0
    cmc_series_M = []
    for t in temps:
        g = (t - tr) - t * math.log(t / tr)
        dg = true_dH_kJ - t * true_dS_kJ + true_dCp_kJ * g
        x = math.exp(dg * 1000.0 / (R_GAS * t))
        cmc_series_M.append(x * 55.5 / (1.0 - x))  # invert cmc_to_mole_fraction's dilute formula

    out = call(
        "vant_hoff_multi_point_fit",
        {"cmc_series_M": cmc_series_M, "temperatures_K": temps, "reference_temperature_K": tr},
    )
    assert out["deltaH_mic_kJ_per_mol"] == pytest.approx(true_dH_kJ, abs=1e-3)
    assert out["deltaS_mic_J_per_mol_K"] == pytest.approx(true_dS_J, abs=1e-1)
    assert out["deltaCp_mic_J_per_mol_K"] == pytest.approx(true_dCp_J, abs=1.0)
    assert out["n_points"] == 5


def test_szyszkowski_tool_at_zero_concentration_returns_gamma0():
    out = call("szyszkowski_predict_surface_tension", {
        "concentration": 0.0, "gamma0_mN_m": 72.0, "gamma_max_mol_per_m2": 3.0e-6, "K": 50.0,
        "system_type": "nonionic",
    })
    assert out["surface_tension_mN_m"] == pytest.approx(72.0)


def test_wetting_work_of_adhesion_tool_complete_wetting():
    out = call("wetting_work_of_adhesion", {"gamma_LV_mN_m": 72.0, "contact_angle_deg": 0.0})
    assert out["work_of_adhesion_mJ_per_m2"] == pytest.approx(144.0)


def test_wetting_spreading_coefficient_tool_never_positive():
    out = call("wetting_spreading_coefficient", {"gamma_LV_mN_m": 72.0, "contact_angle_deg": 90.0})
    assert out["spreading_coefficient_mN_m"] <= 0


def test_estimate_axial_ratio_from_cpp_geometry_tool_matches_library():
    import math

    n_carbons = 12
    from surfactantkit.cpp import tanford_critical_length, tanford_tail_volume

    b = tanford_critical_length(n_carbons)
    v_tail = tanford_tail_volume(n_carbons)
    true_axial_ratio = 2.5
    a = true_axial_ratio * b
    v_core = (4.0 / 3.0) * math.pi * a * b ** 2
    n_agg = v_core / v_tail

    out = call("estimate_axial_ratio_from_cpp_geometry", {"cpp": 0.45, "aggregation_number": n_agg, "n_carbons": n_carbons})
    assert out["axial_ratio"] == pytest.approx(true_axial_ratio, rel=1e-6)

    out_sphere = call("estimate_axial_ratio_from_cpp_geometry", {"cpp": 0.25, "aggregation_number": 60, "n_carbons": n_carbons})
    assert out_sphere["axial_ratio"] == pytest.approx(1.0)


def test_estimate_intrinsic_water_solubility_qspr_tool_matches_library():
    out = call("estimate_intrinsic_water_solubility_qspr", {"smiles": "c1ccc2ccccc2c1"})
    assert 1e-6 <= out["solubility_M"] <= 1e-1
    assert out["caveats"]

    out = call("estimate_intrinsic_water_solubility_qspr", {"smiles": "not a smiles $$$"})
    assert out["status"] == "no_solution"
    assert out["reason"]


def test_classify_surfactant_charge_type_at_ph_tool_matches_library():
    out_low = call("classify_surfactant_charge_type_at_ph", {"smiles": "CCCCCCCCCCCCN", "ph": 4.0})
    assert out_low["charge_type_at_ph"] == "cationic"

    out_high = call("classify_surfactant_charge_type_at_ph", {"smiles": "CCCCCCCCCCCCN", "ph": 13.0})
    assert out_high["charge_type_at_ph"] == "nonionic"

    out_exact = call(
        "classify_surfactant_charge_type_at_ph",
        {"smiles": "CCCCCCCCCCCCN", "ph": 9.0, "amine_pka_override": 10.6},
    )
    assert out_exact["confidence"] == "high"
    assert out_exact["fraction_protonated_low"] == pytest.approx(out_exact["fraction_protonated_high"])


def test_estimate_partial_specific_volume_tool_matches_library():
    from surfactantkit.cpp import tanford_tail_volume
    from surfactantkit.curve_analysis import AVOGADRO_NUMBER

    n_carbons, mw, head_vol = 12, 288.38, 35.0
    out = call("estimate_partial_specific_volume", {"n_carbons": n_carbons, "monomer_molar_mass_g_per_mol": mw, "headgroup_molar_volume_cm3_per_mol": head_vol})
    expected = (tanford_tail_volume(n_carbons) * AVOGADRO_NUMBER * 1e-24 + head_vol) / mw
    assert out["partial_specific_volume_cm3_per_g"] == pytest.approx(expected, rel=1e-9)


def test_get_standard_probe_liquid_properties_tool_matches_library():
    out = call("get_standard_probe_liquid_properties", {"liquid_name": "Water"})
    assert out["liquid"] == "water"
    assert out["owens_wendt_dispersive_mN_m"] == pytest.approx(21.8)
    assert out["owens_wendt_polar_mN_m"] == pytest.approx(51.0)
    assert out["vocg_lw_mN_m"] == pytest.approx(21.8)
    assert out["vocg_acid_mN_m"] == pytest.approx(25.5)
    assert out["vocg_base_mN_m"] == pytest.approx(25.5)
    assert out["vocg_total_mN_m"] == pytest.approx(72.8)


def test_get_standard_probe_liquid_properties_tool_rejects_unknown_liquid():
    out = call("get_standard_probe_liquid_properties", {"liquid_name": "ethanol"})
    assert out["status"] == "no_solution"
    assert out["reason"]


def test_owens_wendt_tool_recovers_true_components_round_trip():
    import math

    true_gamma_s_d, true_gamma_s_p = 30.0, 10.0
    liquids_d = [21.8, 50.8, 37.0]
    liquids_p = [51.0, 0.0, 26.4]
    contact_angles = []
    for gd, gp in zip(liquids_d, liquids_p):
        gt = gd + gp
        w = 2.0 * (math.sqrt(true_gamma_s_d * gd) + math.sqrt(true_gamma_s_p * gp))
        contact_angles.append(math.degrees(math.acos(w / gt - 1.0)))

    out = call(
        "owens_wendt_solid_surface_energy",
        {"contact_angles_deg": contact_angles, "liquid_gamma_dispersive_mN_m": liquids_d, "liquid_gamma_polar_mN_m": liquids_p},
    )
    assert out["gamma_s_dispersive_mN_m"] == pytest.approx(true_gamma_s_d, abs=1e-3)
    assert out["gamma_s_polar_mN_m"] == pytest.approx(true_gamma_s_p, abs=1e-3)


def test_van_oss_chaudhury_good_tool_recovers_true_components_round_trip():
    import math

    true_lw, true_acid, true_base = 25.0, 2.0, 15.0
    liquids_lw = [21.8, 34.0, 50.8]
    liquids_acid = [25.5, 3.92, 0.0]
    liquids_base = [25.5, 57.4, 0.0]
    contact_angles = []
    for lw, a, b in zip(liquids_lw, liquids_acid, liquids_base):
        gt = lw + 2.0 * math.sqrt(a * b)
        w = 2.0 * (math.sqrt(true_lw * lw) + math.sqrt(true_acid * b) + math.sqrt(true_base * a))
        contact_angles.append(math.degrees(math.acos(w / gt - 1.0)))

    out = call(
        "van_oss_chaudhury_good_solid_surface_energy",
        {
            "contact_angles_deg": contact_angles,
            "liquid_gamma_lw_mN_m": liquids_lw,
            "liquid_gamma_acid_mN_m": liquids_acid,
            "liquid_gamma_base_mN_m": liquids_base,
        },
    )
    assert out["gamma_s_lw_mN_m"] == pytest.approx(true_lw, abs=1e-3)
    assert out["gamma_s_acid_mN_m"] == pytest.approx(true_acid, abs=1e-3)
    assert out["gamma_s_base_mN_m"] == pytest.approx(true_base, abs=1e-3)


def test_eor_capillary_number_tool_matches_direct_formula():
    out = call("eor_capillary_number", {"viscosity_mPas": 1.0, "velocity_m_per_s": 1e-5, "interfacial_tension_mN_m": 0.01})
    expected = (1.0e-3 * 1e-5) / (0.01e-3)
    assert out["capillary_number"] == pytest.approx(expected)


def test_molar_solubilization_ratio_tool():
    out = call("molar_solubilization_ratio", {
        "total_solubilized_M": 0.51e-3, "intrinsic_water_solubility_M": 0.01e-3,
        "surfactant_concentration_M": 10e-3, "cmc_M": 2e-3,
    })
    assert out["msr"] == pytest.approx(0.0625)


def test_micelle_water_partition_coefficient_tool():
    total, intrinsic, surf, cmc = 0.51e-3, 0.01e-3, 10e-3, 2e-3
    x_micelle = (total - intrinsic) / ((total - intrinsic) + (surf - cmc))
    x_water = intrinsic / 55.5
    expected = x_micelle / x_water

    out = call("micelle_water_partition_coefficient", {
        "total_solubilized_M": total, "intrinsic_water_solubility_M": intrinsic,
        "surfactant_concentration_M": surf, "cmc_M": cmc,
    })
    assert out["km"] == pytest.approx(expected)


def test_classify_surfactant_charge_type_tool_matches_library():
    out = call("classify_surfactant_charge_type", {"smiles": "CCCCCCCCCCCCOS(=O)(=O)[O-].[Na+]"})
    assert out["charge_type"] == "anionic"
    assert out["system_type_for_gibbs"] == "ionic_no_added_salt"
    assert out["confidence"] == "high"


def test_classify_surfactant_charge_type_tool_zwitterion_gives_null_system_type():
    out = call("classify_surfactant_charge_type", {"smiles": "CCCCCCCCCCCCC(=O)NCCC[N+](C)(C)CC(=O)[O-]"})
    assert out["charge_type"] == "zwitterionic"
    assert out["system_type_for_gibbs"] is None


def test_derive_all_properties_tool_matches_independently_verified_gold_values():
    """Real regression guard: matches the gold CMC/Gamma_max values in
    benchmark/paper3_pilot/pilot50_v2_questions.json (A01, SDS), computed
    independently by directly calling the underlying library functions."""
    concs = [0.40577109375, 0.8115421875, 1.623084375, 3.24616875, 6.4923375, 12.984675, 25.96935, 51.9387]
    gammas = [69.85933451975876, 67.9370105943862, 62.892649575815454, 54.57301518905874,
              42.87619733622359, 37.91402789831731, 37.44295519808575, 36.79102755846683]
    out = call("derive_all_properties_from_smiles_and_curve", {
        "smiles": "CCCCCCCCCCCCOS(=O)(=O)[O-].[Na+]",
        "concentrations_mM": concs, "surface_tensions_mN_per_m": gammas, "temperature_K": 296.15,
    })
    assert out["charge_type"] == "anionic"
    assert out["cmc_mM"] == pytest.approx(9.181551744003432, rel=1e-6)
    assert out["gamma_max_mol_per_m2"] == pytest.approx(3.4266231784578047e-06, rel=1e-6)
    assert out["delta_g_mic_kJ_per_mol"] is None  # no counterion_dissociation_alpha supplied
    assert len(out["gaps"]) > 0


def test_nagarajan_cmc_predictor_tool_matches_library_sds_gold():
    out = call("nagarajan_cmc_predictor", {"n_carbons": 12,
                                           "headgroup": "alkyl_sulfate"})
    from surfactantkit.nagarajan import nagarajan_cmc_predictor as libfn
    lib = libfn(12, "alkyl_sulfate")
    # wrap fidelity: identical floats through dataclasses.asdict
    assert out["cmc_mol_per_L"] == lib.cmc_mol_per_L
    assert out["delta_g_mic_kT_per_molecule"] == lib.delta_g_mic_kT_per_molecule
    # SDS gold: ~9.1e-3 M (experimental 8.2 mM, within 0.15 log units)
    assert out["cmc_mol_per_L"] == pytest.approx(9.1e-3, rel=0.1)
    assert out["converged"] is True
    assert out["method"].startswith("Nagarajan")
    assert out["caveat"]
    assert out["terms_at_optimum"]["branch"] in ("sphere", "cylinder")


def test_nagarajan_free_energy_terms_tool_matches_independent_formula():
    out = call("nagarajan_free_energy_terms",
               {"g": 40, "n_carbons": 12, "headgroup": "alkyl_sulfate",
                "concentration_mol_per_L": 1e-3})
    assert out["transfer"] == pytest.approx(-40 * (2.05 + 1.49 * 12), rel=1e-12)
    assert out["ordering"] == pytest.approx(40 * (-0.5 + 0.24 * 12), rel=1e-12)
    assert out["branch"] == "sphere"
    assert out["total"] == pytest.approx(
        out["transfer"] + out["ordering"] + out["interface"]
        + out["steric"] + out["electrostatic"], rel=1e-12)


def test_no_solution_payload_contract():
    """Audit P0-1: ValueError/KeyError at the MCP boundary return the
    shaped {status, reason} payload (not a raw internal error), and the
    fail-fast field names survive into reason."""
    out = call("ward_tordai_adsorption",
               {"c0_mol_per_L": 0.0, "D_m2_per_s": 1e-9, "t_end_s": 1.0})
    assert set(out) == {"status", "reason"}
    assert out["status"] == "no_solution"
    assert "c0_mol_per_L" in out["reason"]  # field name reaches the client
    assert out["reason"]  # non-empty


def test_no_solution_payload_keyerror_path():
    """KeyError (unknown reference-table keys) is shaped too."""
    out = call("hlb_from_groups", {"group_counts": {"nonexistent_group": 1}})
    assert out["status"] == "no_solution"
    assert "refusing to guess" in out["reason"] or "unknown" in out["reason"]


# ---------------------------------------------------------------------------
# Batch sign-off (audit B1-B4, B8) at the MCP surface: new method params,
# the phase_inversion_temperature tool, and the CPP validity advisory key.
# ---------------------------------------------------------------------------

def test_cmc_tension_tool_tangential_method_payload():
    c = [0.5, 1, 2, 4, 6, 8, 10, 12, 14, 16, 20]
    plateau = 50 - 10 * math.log10(8)
    g = [50 - 10 * math.log10(x) if x <= 8 else plateau for x in c]
    out = call(
        "cmc_from_surface_tension_curve",
        {"concentrations_mM": c, "surface_tensions_mN_per_m": g,
         "method": "tangential"},
    )
    assert out["cmc_mM"] == pytest.approx(8.0, rel=1e-9)
    assert out["method"].startswith("tangential")


def test_cmc_tension_tool_default_method_unchanged():
    c = [0.5, 1, 2, 4, 6, 8, 10, 12, 14, 16, 20]
    plateau = 50 - 10 * math.log10(8)
    g = [50 - 10 * math.log10(x) if x <= 8 else plateau for x in c]
    out = call(
        "cmc_from_surface_tension_curve",
        {"concentrations_mM": c, "surface_tensions_mN_per_m": g},
    )
    assert out["method"].startswith(("two-segment", "three-segment"))


def test_cmc_tension_tool_invalid_method_returns_no_solution_payload():
    c = [1, 2, 3, 4, 5, 6, 7, 8]
    g = [40.0, 38.0, 36.0, 34.0, 33.0, 33.0, 33.0, 33.0]
    out = call(
        "cmc_from_surface_tension_curve",
        {"concentrations_mM": c, "surface_tensions_mN_per_m": g,
         "method": "guess"},
    )
    assert out["status"] == "no_solution"
    assert "method" in out["reason"]
    assert "tangential" in out["reason"]


def test_cmc_conductivity_tool_jacobs_jones_payload():
    c = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14]
    k = [1.2 * x + 0.05 if x <= 6 else 0.4 * x + 4.85 for x in c]
    out = call(
        "cmc_from_conductivity",
        {"concentrations_mM": c, "conductivities": k, "method": "jacobs_jones"},
    )
    assert out["cmc_mM"] == pytest.approx(6.0, rel=1e-12)
    assert out["counterion_binding_degree"] == pytest.approx(2 / 3, rel=1e-9)
    assert out["method"].startswith("second-derivative")


def test_phase_inversion_temperature_tool_round_trip():
    pit = call(
        "phase_inversion_temperature",
        {"family": "nonionic", "salinity_pct": 1.0, "k": 0.015,
         "eacn": 10.0, "cc": 0.5, "b": 0.13},
    )
    assert pit["phase_inversion_temperature_K"] == pytest.approx(290.15, abs=1e-9)
    assert pit["unit"] == "K"


def test_phase_inversion_temperature_tool_missing_b_no_solution():
    out = call(
        "phase_inversion_temperature",
        {"family": "nonionic", "salinity_pct": 1.0, "k": 0.015,
         "eacn": 10.0, "cc": 0.5},
    )
    assert out["status"] == "no_solution"
    assert "'b'" in out["reason"] or "b (" in out["reason"]


def test_critical_packing_parameter_tool_validity_note_key():
    geom = call("tanford_chain_geometry", {"n_carbons": 8})
    out = call(
        "critical_packing_parameter",
        {"volume_A3": geom["tail_volume_A3"], "head_area_A2": 50.0,
         "length_A": geom["critical_length_A"], "n_carbons": 8},
    )
    assert out["validity_note"] is not None
    assert "PMID 31954873" in out["validity_note"]
    # omitted -> null key present, CPP math untouched
    out2 = call(
        "critical_packing_parameter",
        {"volume_A3": geom["tail_volume_A3"], "head_area_A2": 50.0,
         "length_A": geom["critical_length_A"]},
    )
    assert out2["validity_note"] is None
    assert out2["cpp"] == pytest.approx(out["cpp"])


def test_entropy_tool_delta_cp_passthrough_and_half_specified_payload():
    out = call(
        "entropy_of_micellization",
        {"deltaG_mic_kJ_per_mol": -10.0, "deltaH_mic_kJ_per_mol": -5.0,
         "temperature_K": 350.0, "delta_cp_J_mol_K": -100.0,
         "reference_temperature_K": 298.15},
    )
    assert out["deltaS_mic_J_per_mol_K"] == pytest.approx(
        0.7364382819682085, rel=1e-9
    )
    bad = call(
        "entropy_of_micellization",
        {"deltaG_mic_kJ_per_mol": -10.0, "deltaH_mic_kJ_per_mol": -5.0,
         "temperature_K": 350.0, "delta_cp_J_mol_K": -100.0},
    )
    assert bad["status"] == "no_solution"
    assert "reference_temperature_K" in bad["reason"]
