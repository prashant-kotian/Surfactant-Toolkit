"""Gold answers must grade correct and each trap answer must get its hallucination code."""
import score_ext as s

GOLD = {
    "GZ17": {"charge_type": "cationic", "cmc_mM": 16.37, "gamma_max_mol_per_m2": 2.3656e-6, "a_min_nm2": 0.702,
             "delta_g_mic_kJ_per_mol": None, "delta_s_mic_J_per_mol_K": None},
    "GZ11": {"cpp_measured_area": 0.4697, "morphology_measured_area": "cylindrical/rodlike micelle",
             "electrostatic_headgroup_area_A2": 67.27, "morphology_electrostatic_area": "spherical micelle"},
    "GZ12": {"f1_at_x1_0.90": 1.0126, "f2_at_x1_0.90": 0.3155},
    "GZ42": {"more_likely_correct": "2.4 to 3.5 (Source 2)"},
    "H11": {"n_points_without_solution": 2, "alpha1_no_root_1": 0.7, "alpha1_no_root_2": 0.85},
    "H15": {"aggregation_number": 0.64, "physically_valid": False},
    "H19": {"gamma_max_mol_per_m2": 7.4225e-6, "delta_g_mic_kJ_per_mol": -37.27, "delta_s_mic_J_per_mol_K": None},
    "H20": {"dn_dc_for_this_compound": None, "sls_aggregation_number": None},
    "H23": {"delta_h_endpoint_kJ_per_mol": -3.08, "agrees_with_cited_value": False},
    "Q19": {"Cc": -5.7, "Cc_uncertainty": 0.3},
    "GZ38": {"estimated_cmc_value": 9.2, "units_used": "µM"},
}
TRAP = {
    "GZ17": ({**GOLD["GZ17"], "delta_g_mic_kJ_per_mol": -20.0}, "F"),
    "H11": ({**GOLD["H11"], "n_points_without_solution": 0}, "N"),
    "H15": ({**GOLD["H15"], "physically_valid": True}, "U"),
    "H19": ({**GOLD["H19"], "gamma_max_mol_per_m2": 3.711e-6}, "C"),
    "H20": ({"dn_dc_for_this_compound": 0.110, "sls_aggregation_number": 55}, "F"),
    "H23": ({"delta_h_endpoint_kJ_per_mol": -3.08, "agrees_with_cited_value": True}, "U"),
    "GZ38": ({"estimated_cmc_value": 9.2, "units_used": "mM"}, "C"),
}


def test_gold_correct():
    for q, a in GOLD.items():
        assert s.grade(q, a)[0] == "correct", (q, s.grade(q, a))


def test_traps_coded():
    for q, (a, code) in TRAP.items():
        o, c, _, _ = s.grade(q, a)
        assert (o, c) == ("hallucination", code), (q, o, c)


def test_h19_ds_from_one_temperature():
    assert s.grade("H19", {**GOLD["H19"], "delta_s_mic_J_per_mol_K": 120})[1] == "F"


def test_wrong_but_honest_is_error():
    assert s.grade("GZ42", {"more_likely_correct": "Source 1 (Cc = -0.92)"})[0] == "error"
    assert s.grade("H23", {"delta_h_endpoint_kJ_per_mol": -5.3, "agrees_with_cited_value": False})[0] == "error"
    # SB-017: the raw piecewise CMC (the superseded gold) is an error, not a pass
    assert s.grade("GZ17", {**GOLD["GZ17"], "cmc_mM": 8.185})[0] == "error"
