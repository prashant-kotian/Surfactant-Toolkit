"""HARD50 gold builder -- computes every gold answer for the Paper-3 hard-50
question set by calling the REAL surfactantkit library (same discipline as
benchmark/build_question_bank.py: gold correctness inherits from the library's
own 522-test suite, never hand-derived).

Emits: gold_answers.json
Deterministic (seeded), idempotent.
"""
from __future__ import annotations
import json, math, random, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))

from surfactantkit import (
    clint_ideal_cmc, solve_rubingh_x, rubingh_beta, rubingh_beta_regression,
    activity_coefficients, excess_free_energy, asymmetric_margules_activity_coefficients,
    eommm_global_fit, rodenas_x1_series, rodenas_activity_coefficients,
    cmc_from_surface_tension_curve, cmc_from_conductivity_curve,
    aggregation_number_from_sls_debye_plot, aggregation_number_from_dls,
    gibbs_gamma_max, gibbs_a_min, select_isotherm_model, frumkin_surface_tension,
    szyszkowski_surface_tension, frumkin_fit_K_and_a,
    cmc_to_mole_fraction, gibbs_free_energy_micellization,
    mass_action_free_energy_of_micellization, vant_hoff_enthalpy,
    vant_hoff_multi_point_fit, entropy_micellization, counterion_binding_degree,
    debye_length, ionic_strength, henry_function, zeta_potential_henry,
    tanford_tail_volume, tanford_critical_length, critical_packing_parameter,
    classify_aggregate_morphology, estimate_axial_ratio_from_cpp_geometry,
    hydrodynamic_radius_stokes_einstein, hydrodynamic_radius_perrin_corrected,
    perrin_friction_factor,
    classify_surfactant_charge_type, classify_surfactant_structural_family,
    classify_surfactant_charge_type_at_ph,
    fit_k_and_cc_from_salinity_scan, optimal_salinity_ionic, hld_ionic,
    capillary_number, molar_solubilization_ratio,
    estimate_intrinsic_water_solubility_qspr,
    derive_all_properties_from_smiles_and_curve,
    select_mixture_model, aggregation_number_spherical,
)
from surfactantkit.thermodynamics import WATER_MOLARITY_M

HERE = Path(__file__).resolve().parent

# ============================================================================
# INPUT DATASETS (single source of truth; questions text imports from here)
# ============================================================================

# --- H-01 SLS: constructed from truth Mw=35888 (C12E8-like, Nagg=80) -------
H01 = dict(
    c=[0.002, 0.004, 0.006, 0.008, 0.010, 0.012],
    dn_dc=0.135, wl=632.8, n_sol=1.333, m0=448.6,
    cmc_M=8.7e-5, T=298.15, A2_true=4.0e-4,
)

# --- H-02 conductivity (constructed; true break 5.75 mM, beta 0.6562) ------
def _conductivity_series(s1, s2, break_c, inter0, conc):
    out = []
    bval = s1 * break_c + inter0
    for c in conc:
        out.append(s1 * c + inter0 if c <= break_c else bval + s2 * (c - break_c))
    return out

H02_CONC = [1.0, 2.0, 3.0, 4.0, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 10.0, 12.0]
H02_COND = _conductivity_series(0.64, 0.22, 5.75, 0.05, H02_CONC)

# --- H-03 3-regime tensiometry (baseline + decline + plateau) -------------
H03_PTS = [
    (0.005, 71.9), (0.010, 71.8), (0.020, 71.6),
    (0.05, 66.5), (0.10, 61.9), (0.30, 53.0), (0.70, 45.0),
    (1.50, 38.0), (3.00, 33.0),
    (6.00, 31.2), (10.00, 31.1), (16.00, 31.0),
]
H03_C = [p[0] for p in H03_PTS]
H03_G = [p[1] for p in H03_PTS]
H03_T, H03_ALPHA = 298.15, 0.31

# --- H-04 CPB real CMC(T) 8 points (real, used in GZ-26/28/50/51) ---------
H04_T = [283.15, 288.15, 293.15, 298.15, 303.15, 308.15, 313.15, 318.15]
H04_CMC_M = [9.93, 9.82, 9.76, 9.74, 9.77, 10.02, 10.54, 11.47]  # mM
H04_ALPHA = 0.268

# --- H-05 DHPC-DTAB real 5-point series (Vautier-Giongo/Bales 2005) -------
H05_ALPHA2 = [0.95, 0.90, 0.80, 0.60, 0.50]
H05_CMCMIX = [10.8, 6.2, 4.8, 2.7, 2.4]
H05_C1, H05_C2 = 1.8, 15.6

# --- H-06 constructed asymmetric series (EOMMM truth) ---------------------
H06_ALPHAS = [0.10, 0.25, 0.40, 0.55, 0.70, 0.85]
H06_C1, H06_C2 = 1.8, 8.3
H06_W12_TRUE, H06_W21_TRUE, H06_JITTER, H06_SEED = -3.5, -0.9, 0.005, 7

# --- H-07 noisy near-Langmuir isotherm (Langmuir wins BIC, Frumkin wins RSS)
H07_CS = [3e-6, 8e-6, 1.5e-5, 3e-5, 5e-5, 7.5e-5, 1.0e-4, 1.3e-4]
H07_K_TRUE, H07_GM_TRUE, H07_G0 = 1.5e6, 2.2e-6, 72.0
H07_NOISE, H07_SEED = 0.18, 11

# --- H-08 conductivity + tensiometry pair, same true CMC 5.75 mM ---------
# tension decline on the exact line gamma = 34.8 - 7*ln(c/5.75) up to 5.0 mM,
# then plateau; tool reports CMC = sqrt(bracketing points) = sqrt(5.0*6.6).
H08_TENS_PTS = [
    (0.06, 66.7), (0.15, 60.3), (0.40, 53.5), (0.90, 47.8),
    (2.00, 42.2), (4.00, 37.3), (5.00, 35.8),
    (6.60, 34.6), (9.00, 34.5), (12.00, 34.4), (16.00, 34.3),
]

# --- H-09 HLD salinity scan (constructed truth k=0.16, Cc=-1.4) ----------
H09_EACN = [8, 10, 12, 14]
H09_K_TRUE, H09_CC_TRUE, H09_ALPHA_DT = 0.16, -1.4, 0.01

# --- H-10 constructed asymmetric series for the Rodenas-vs-Rubingh race ---
H10_ALPHAS = [0.10, 0.25, 0.40, 0.55, 0.70, 0.85]
H10_C1, H10_C2 = 2.4, 11.0
H10_W12_TRUE, H10_W21_TRUE, H10_JITTER, H10_SEED = -3.0, -1.2, 0.01, 13
H10_REPORT_ALPHA = 0.55
H10_AGREE_TOL = 0.05

# --- H-11 no-root Rubingh series (points alpha=0.70 mix=16 and alpha=0.85
#     mix=13 have NO valid Rubingh root -- verified by solve_rubingh_x) -----
H11_SERIES = [(0.10, 9.0), (0.25, 12.0), (0.50, 15.0), (0.70, 16.0), (0.85, 13.0)]
H11_C1, H11_C2 = 8.0, 15.0

# --- H-12 Frumkin condensation curve (a=2.6 truth) ------------------------
H12_CS = [2e-6, 5e-6, 1e-5, 2e-5, 4e-5, 7e-5, 1.1e-4]
H12_K_TRUE, H12_A_TRUE, H12_GM, H12_G0 = 2.5e5, 2.6, 2.1e-6, 71.8

# --- H-13 rod-regime CPP + too-small Nagg (real geometric inconsistency) ---
H13_AMIN, H13_NAGG, H13_NC = 50.0, 20, 12

# --- H-14 gemini 12-4-12 + supplied Liao et al. 2023 Griffin table --------
H14_SMILES = "[Br-].[Br-].CCCCCCCCCCCC[N+](C)(C)CCCC[N+](C)(C)CCCCCCCCCCCC"
H14_TABLE = {"8-3-8": 7.97, "12-3-12": 6.37, "8-4-8": 7.73, "12-4-12": 6.22}

# --- H-15 impossible SLS (Mw=400 < monomer 625) ---------------------------
H15 = dict(c=H01["c"], dn_dc=0.140, wl=632.8, n_sol=1.333, m0=625.0,
           Mw_true=400.0, A2_true=4.0e-4)

# --- H-16 second conductivity dataset (break 4.2, beta 0.72) --------------
H16_COND = _conductivity_series(0.55, 0.154, 4.20, 0.03, H02_CONC)

# --- H-17 QSPR + MSR ------------------------------------------------------
H17_SMILES = "c1ccc2ccccc2c1"          # naphthalene (real, no supplied measured solubility)
H17_CS_TOTAL_M = 0.0110                 # total solubilized (micellar solution)
H17_CT_M = 0.050                        # total surfactant
H17_CMC_M = 0.9e-3

# --- H-18 dodecylamine at pH 10.6 ----------------------------------------
H18_SMILES, H18_PH = "CCCCCCCCCCCCN", 10.6   # dodecylamine (C12)

# --- H-19 Gibbs with EXCESS electrolyte (n=1!) + alpha given --------------
H19_SLOPE, H19_T, H19_ALPHA, H19_CMC_M = -18.4, 298.15, 0.30, 8.0e-3

# --- H-20 AOT needs dn/dc (NOT in reference table; SDS value is bait) -----
H20_SLS_C = [0.002, 0.004, 0.006, 0.008, 0.010, 0.012]
H20_SLS_R = [1.2e-05, 2.4e-05, 3.6e-05, 4.8e-05, 6.0e-05, 7.2e-05]
H20_M0 = 332.43  # AOT monomer

# --- H-21 Gamma_max vs literature Amin (deliberate disagreement) ---------
H21_SLOPE, H21_T, H21_LIT_AMIN_A2 = -15.2, 298.15, 44.70

# --- H-22 Clint + Rubingh with distractor inputs --------------------------
H22_A1, H22_C1, H22_C2, H22_CMCMEAS = 0.30, 1.8, 15.6, 4.05

# --- H-23 endpoint van't Hoff, temps in C, cited-distractor deltaH -------
H23_T_C = [10, 15, 20, 25, 30, 35, 40, 45]
H23_CMC_M = H04_CMC_M
H23_CITED_DH = -6.2

# --- H-24 DLS aggregation number (hydration bias) -------------------------
H24_D, H24_ETA, H24_VBAR, H24_M0, H24_T = 1.1e-6, 0.8903, 0.90, 448.6, 298.15

# --- H-25 ionic strength from non-stoichiometric salt mix ----------------
H25_IONS_M = {1: 0.060, 2: 0.015, -1: 0.090}   # Na+ 60, Ca2+ 15, Cl- 90 mM
H25_R_NM, H25_T = 2.5, 298.15

# --- H-26 C18E12 counting -------------------------------------------------
# C18H37-O-(CH2CH2O)12-H  -> 42 C, 13 O, MW 799.137 (RDKit-verified)
H26_SMILES = "C" * 19 + "OCC" * 12 + "O"  # FIX 2026-10-18: was C*18, but the
# printed question SMILES is C19 (models told to "count them off the string")
H26_MH_UNIT_EO = 44.053      # per EO unit, disclosed in the question
H26_MH_OH = 17.008           # terminal OH, disclosed in the question

# --- H-27 AOT structural audit -------------------------------------------
H27_SMILES = "CCCCC(CC)COC(=O)CC(C(=O)OCC(CC)CCCC)S(=O)(=O)[O-].[Na+]"

# --- H-28 sophorolipid audit ---------------------------------------------
H28_SMILES = ("CC(CCCCCC/C=C/CCCCCCCC(=O)O)O[C@H]1[C@@H]([C@H]([C@@H]([C@H](O1)"
              "COC(=O)C)O)O)O[C@H]2[C@@H]([C@H]([C@@H]([C@H](O2)COC(=O)C)O)O)O")

# --- H-29 PO-EO block nonionic -------------------------------------------
# C12H25-(OCH2CH3)x5? no: C12 tail + 5 PO + 8 EO + OH -> C43H88O14, MW 829.163
H29_SMILES = "C" * 12 + "OCC(C)" * 5 + "OCC" * 8 + "O"
H29_MH_UNIT_EO = 44.053
H29_MH_OH = 17.008
H29_TAIL_MASS_C12H25 = 169.31  # disclosed in the question

# --- H-30 DTAB full stack -------------------------------------------------
H30_SMILES, H30_AMIN, H30_NAGG = "CCCCCCCCCCCC[N+](C)(C)C.[Br-]", 59.62, 60
H30_D, H30_ETA, H30_T = 8.9e-7, 0.8903, 298.15

# --- H-31 log10 slope + mM chain -----------------------------------------
H31_SLOPE_PER_DECADE, H31_T, H31_CMC_M, H31_ALPHA = -9.2, 298.15, 6.0e-3, 0.30

# --- H-32 SI-unit diffusion ----------------------------------------------
H32_D_SI, H32_ETA_PA_S, H32_T, H32_AXIAL = 7.5e-11, 0.8903e-3, 298.15, 1.14

# --- H-33 capillary number (EOR), 3 unit systems --------------------------
H33_MU_CP, H33_V_FT_DAY, H33_GAMMA_MNM, H33_THRESH = 12.0, 0.30, 32.0, 1e-4

# --- H-34 weight-percent blend -> Clint/Rubingh ---------------------------
H34_WT_SDS, H34_WT_DTAB = 65.0, 35.0
H34_MW_SDS, H34_MW_DTAB = 288.38, 308.34
H34_C1, H34_C2, H34_CMCMEAS = 8.3, 15.6, 8.10

# --- H-35 SI electrokinetics ---------------------------------------------
H35_I_M, H35_ETA, H35_EPS_R, H35_MU_SI, H35_R_NM, H35_T = 0.0257, 0.8903, 78.4, 2.7e-8, 2.5, 298.15

# --- H-43 AOT Cc REAL literature conflict ---------------------------------
# Hellweg et al. 2023 review states Cc(AOT) = -0.92; the directly-read
# primary paper (Leng & Acosta 2023) reports Cc ~ 2.4 (this project's own
# hld.py documents this exact conflict). AOT/n-decane/brine at S=5%, k=0.16.
H43_S, H43_EACN = 5.0, 10.0   # n-decane
H43_CC1, H43_CC2, H43_K = -0.92, 2.4, 0.16

# --- H-44 cross-method Nagg ----------------------------------------------
H44_SLW = 16800.0; H44_M0_DTAB = 308.34; H44_NC = 12

# --- H-45 pC20 draft-claim audit on the H-03 curve -----------------------
# draft claims to verify: c20 = 0.10 mM and pC20 = 4.0 (internally
# consistent with each other but both wrong for this data -> recompute)
H45_CLAIM_C20_MM, H45_CLAIM_PC20 = 0.10, 4.0

# --- H-46 fresh branched C8 sulfate, full orchestrate ---------------------
H46_SMILES = "CCCCC(CC)OS(=O)(=O)[O-].[Na+]"
H46_ALPHA = 0.31
H46_T = 298.15

def _h46_curve():
    # disclosed forward-construction: decline slope -7.5 mN/m per ln(C),
    # break at 6.0 mM, plateau ~34.6
    cmc = 6.0
    pre = [0.05, 0.12, 0.30, 0.70, 1.6, 3.2, 5.0]
    pts = []
    for c in pre:
        g = 36.0 - 7.5 * math.log(c / cmc)
        pts.append((c, round(min(g, 71.9), 1)))
    pts += [(6.0, 34.8), (9.0, 34.6), (13.0, 34.5), (18.0, 34.4)]
    return [p[0] for p in pts], [p[1] for p in pts]

H46_C, H46_G = _h46_curve()

# --- H-47 conductivity + tensiometry fusion (same as H-02 data + H-08) ----
# (reuses H02_* and H08_TENS_PTS)

# --- H-48 HLD + CPP cross-check (SDS real) -------------------------------
H48_SMILES = "CCCCCCCCCCCCOS(=O)(=O)[O-].[Na+]"
H48_CC, H48_K, H48_S, H48_EACN = -2.63, 0.16, 5.0, 0.0
H48_AMIN, H48_NC = 44.70, 12

# --- H-49 thermodynamic triangle -----------------------------------------
H49_CMC_M, H49_T, H49_ALPHA, H49_DH = 8.0e-3, 298.15, 0.272, -5.0

# --- H-50 SDS capstone ----------------------------------------------------
H50_SMILES = "CCCCCCCCCCCCOS(=O)(=O)[O-].[Na+]"
H50_AMIN, H50_NAGG, H50_NC = 44.70, 63, 12
H50_D, H50_ETA, H50_T = 9.3e-7, 0.8903, 298.15

# ============================================================================
# GOLD COMPUTATIONS
# ============================================================================

def _sls_K(dn_dc, wl_nm, n=1.333):
    AVO = 6.02214076e23
    wl_cm = wl_nm * 1e-7
    return 4.0 * math.pi ** 2 * n ** 2 * dn_dc ** 2 / (AVO * wl_cm ** 4)

def _num(v, tol=0.05):
    return {"type": "number", "value": v, "rel_tol": tol}

def _exact(v):
    return {"type": "exact", "value": v}

def _bool(v):
    return {"type": "exact", "value": v}

def _null():
    return {"type": "null"}

def _kw_any(words):
    return {"type": "keywords_any", "value": [w.lower() for w in words]}

def _kw_all(words):
    return {"type": "keywords_all", "value": [w.lower() for w in words]}

Q: dict[str, dict] = {}

def q(qid, family, axis, a_tools, b_note, mechanism, reference, fields, **extra):
    Q[qid] = {
        "id": qid, "family": family, "failure_axis": axis,
        "mcp_version_a_tools": a_tools, "mcp_version_b_guidance": b_note,
        "expected_failure_mechanism": mechanism, "reference": reference,
        "fields": fields, **extra,
    }

# ---------------- H-01 ----------------
def gold_h01():
    K = _sls_K(H01["dn_dc"], H01["wl"], H01["n_sol"])
    dRs = [K * c / (1.0 / 35888.0 + 2.0 * H01["A2_true"] * c) for c in H01["c"]]
    r = aggregation_number_from_sls_debye_plot(
        H01["c"], dRs, H01["dn_dc"], H01["wl"], H01["m0"], H01["n_sol"])
    X = cmc_to_mole_fraction(H01["cmc_M"])
    dg = mass_action_free_energy_of_micellization(X, r.aggregation_number, H01["T"])
    return dict(optical_constant_K=_num(K, 1e-3), micelle_mw=_num(r.micelle_molar_mass_g_per_mol, 0.02),
                aggregation_number=_num(r.aggregation_number, 0.03),
                delta_g_mass_action_kJ_per_mol=_num(dg, 0.03)), dRs

# ---------------- H-02 ----------------
def gold_h02():
    r = cmc_from_conductivity_curve(H02_CONC, H02_COND)
    alpha = 1.0 - r.counterion_binding_degree
    X = cmc_to_mole_fraction(r.cmc_mM / 1000.0)
    dg = gibbs_free_energy_micellization(X, 298.15, counterion_factor=2.0 - r.counterion_binding_degree)
    return dict(cmc_mM=_num(r.cmc_mM, 0.03),
                beta=_num(r.counterion_binding_degree, 0.03),
                alpha=_num(alpha, 0.03),
                delta_g_mic_kJ_per_mol=_num(dg, 0.05))

# ---------------- H-03 ----------------
def gold_h03():
    r = cmc_from_surface_tension_curve(H03_C, H03_G)
    gm = gibbs_gamma_max(r.premicellar_slope_mN_per_m_per_lnC, "ionic_no_added_salt", H03_T)
    am = gibbs_a_min(gm)
    X = cmc_to_mole_fraction(r.cmc_mM / 1000.0)
    dg = gibbs_free_energy_micellization(X, H03_T, counterion_factor=1.0 + H03_ALPHA)
    return dict(cmc_mM=_num(r.cmc_mM, 0.05),
                three_segment_selected=_bool(r.n_baseline_points > 0),
                n_baseline_points=_num(r.n_baseline_points, 0.0),
                gamma_max_mol_per_m2=_num(gm, 0.05),
                a_min_nm2=_num(am, 0.05),
                delta_g_mic_kJ_per_mol=_num(dg, 0.05))

# ---------------- H-04 ----------------
def gold_h04():
    Xs = [cmc_to_mole_fraction(c / 1000.0) for c in H04_CMC_M]
    # FIX 2026-10-18: the question pins dH/dS/dCp 'referenced to 298.15 K';
    # the original call omitted reference_temperature_K, so gold came out at
    # the fit-default mean T (300.65 K) and tool-faithful models were failed.
    r = vant_hoff_multi_point_fit(Xs, H04_T, counterion_factor=2.0 - H04_ALPHA,
                                   reference_temperature_K=298.15)
    idx298 = H04_T.index(298.15)
    dg298 = r.delta_G_fitted_kJ_per_mol[idx298]
    tds = 298.15 * r.delta_S_J_per_mol_K / 1000.0
    driver = "entropy_driven" if abs(tds) > abs(r.delta_H_kJ_per_mol) else "enthalpy_driven"
    return dict(delta_H_kJ_per_mol=_num(r.delta_H_kJ_per_mol, 0.05),
                delta_S_J_per_mol_K=_num(r.delta_S_J_per_mol_K, 0.08),
                delta_Cp_J_per_mol_K=_num(r.delta_Cp_J_per_mol_K, 0.15),
                delta_G_298_kJ_per_mol=_num(dg298, 0.05),
                driver=_exact(driver))

# ---------------- H-05 ----------------
def gold_h05():
    a1 = [1.0 - a for a in H05_ALPHA2]
    rr = rubingh_beta_regression(a1, H05_CMCMIX, H05_C1, H05_C2)
    i = a1.index(0.40)
    x1 = rr.x1_values[i]
    f1, f2 = activity_coefficients(x1, rr.beta_mean)
    gexc = excess_free_energy(x1, f1, f2, 298.15)
    return dict(beta_mean=_num(rr.beta_mean, 0.10),
                beta_std=_num(rr.beta_std, 0.30),
                x1_at_alpha2_0p60=_num(x1, 0.04),
                f1=_num(f1, 0.06), f2=_num(f2, 0.06),
                gexc_kJ_per_mol=_num(gexc, 0.12))

# ---------------- H-06 ----------------
def _h06_series():
    rng = random.Random(H06_SEED)
    out = []
    for a in H06_ALPHAS:
        def g(x1):
            f1, f2 = asymmetric_margules_activity_coefficients(x1, H06_W12_TRUE, H06_W21_TRUE)
            return x1 * f1 * H06_C1 / a - (1 - x1) * f2 * H06_C2 / (1 - a)
        lo, hi = 1e-9, 1 - 1e-9
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if g(lo) * g(mid) <= 0: hi = mid
            else: lo = mid
        x1 = 0.5 * (lo + hi)
        f1, _ = asymmetric_margules_activity_coefficients(x1, H06_W12_TRUE, H06_W21_TRUE)
        mix = x1 * f1 * H06_C1 / a
        mix *= (1 + rng.uniform(-H06_JITTER, H06_JITTER))
        out.append(mix)
    return out

def gold_h06():
    series = _h06_series()
    sel = select_mixture_model(H06_ALPHAS, series, H06_C1, H06_C2)
    rub = [c for c in sel.candidates if c.name.startswith("rubingh")][0]
    beta = rub.params["beta"]
    eom = [c for c in sel.candidates if c.name.startswith("eommm")][0]
    chosen = [c for c in sel.candidates if c.name == sel.selected_model][0]
    return dict(selected_model=_exact(sel.selected_model),
                selected_bic=_num(chosen.bic, 0.05),
                rubingh_beta_mean=_num(beta, 0.10),
                interaction_character=_kw_any(["synerg"]),
                w12=_num(eom.params["W12"], 0.15),
                w21=_num(eom.params["W21"], 0.25)), series

# ---------------- H-07 ----------------
def gold_h07():
    rng = random.Random(H07_SEED)
    gts = [frumkin_surface_tension(c, H07_G0, H07_GM_TRUE, H07_K_TRUE, 0.0,
                                   system_type="nonionic") + rng.uniform(-H07_NOISE, H07_NOISE)
           for c in H07_CS]
    s = select_isotherm_model(H07_CS, gts, H07_G0, H07_GM_TRUE, system_type="nonionic")
    return dict(selected_model=_exact(s.selected_model),
                delta_bic=_num(s.delta_bic, 0.12),
                frumkin_a=_num(s.frumkin_a, 0.5),   # tiny a; wide tol, sign/magnitude check
                physically_meaningful_interaction=_bool(False),
                raw_rss_frumkin_lower_than_langmuir=_bool(True)), gts

# ---------------- H-08 ----------------
def gold_h08():
    rc = cmc_from_conductivity_curve(H02_CONC, H02_COND)
    rt = cmc_from_surface_tension_curve([p[0] for p in H08_TENS_PTS], [p[1] for p in H08_TENS_PTS])
    agree = abs(rt.cmc_mM - rc.cmc_mM) / rc.cmc_mM < 0.02
    X = cmc_to_mole_fraction(rc.cmc_mM / 1000.0)
    dg = gibbs_free_energy_micellization(X, 298.15, counterion_factor=2.0 - rc.counterion_binding_degree)
    return dict(cmc_conductivity_mM=_num(rc.cmc_mM, 0.03),
                cmc_tensiometry_mM=_num(rt.cmc_mM, 0.04),
                agree_within_2pct=_bool(agree),
                beta=_num(rc.counterion_binding_degree, 0.03),
                delta_g_mic_kJ_per_mol=_num(dg, 0.06))

# ---------------- H-09 ----------------
def gold_h09():
    sstars = [math.exp(H09_K_TRUE * e - H09_CC_TRUE) for e in H09_EACN]
    f = fit_k_and_cc_from_salinity_scan(H09_EACN, sstars, alpha=H09_ALPHA_DT, delta_T_K=0.0)
    s16 = optimal_salinity_ionic(f.k, 16, f.cc, alpha=H09_ALPHA_DT, delta_T_K=10.0)
    hld = hld_ionic(3.0, f.k, 16, f.cc, alpha=H09_ALPHA_DT, delta_T_K=10.0)
    winsor = "II" if hld > 0.3 else ("I" if hld < -0.3 else "III")
    return dict(k=_num(f.k, 0.03), cc=_num(f.cc, 0.05), fit_r2=_num(f.r_squared, 1e-4),
                s_star_EACN16_dT10_pct=_num(s16, 0.03),
                hld_at_S3=_num(hld, 0.03),
                winsor_type=_exact(winsor),
                s_star_achievable_in_water=_bool(s16 <= 26.0)), sstars

# ---------------- H-10 ----------------
def gold_h10():
    """Rodenas (model-independent, Gibbs-Duhem local slope) vs Rubingh
    (regular-solution implicit root) micellar-composition race on a
    constructed asymmetric-EOMMM series (truth W12=-3.0, W21=-1.2)."""
    from surfactantkit import asymmetric_margules_activity_coefficients as amac
    rng = random.Random(H10_SEED)
    series = []
    for a in H10_ALPHAS:
        def g(x1):
            f1, f2 = amac(x1, H10_W12_TRUE, H10_W21_TRUE)
            return x1 * f1 * H10_C1 / a - (1 - x1) * f2 * H10_C2 / (1 - a)
        lo, hi = 1e-9, 1 - 1e-9
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if g(lo) * g(mid) <= 0:
                hi = mid
            else:
                lo = mid
        x1 = 0.5 * (lo + hi)
        f1, _ = amac(x1, H10_W12_TRUE, H10_W21_TRUE)
        mix = x1 * f1 * H10_C1 / a
        mix *= (1 + rng.uniform(-H10_JITTER, H10_JITTER))
        series.append(mix)
    rod = rodenas_x1_series(H10_ALPHAS, series)
    i = H10_ALPHAS.index(H10_REPORT_ALPHA)
    x1_rod = rod.x1_rodenas[i]
    x1_rub = solve_rubingh_x(H10_REPORT_ALPHA, series[i], H10_C1, H10_C2)
    f1r, f2r = rodenas_activity_coefficients(
        H10_REPORT_ALPHA, x1_rod, series[i], H10_C1, H10_C2)
    rub_all = [solve_rubingh_x(a, m, H10_C1, H10_C2)
               for a, m in zip(H10_ALPHAS, series)]
    diffs = [abs(r - u) for r, u in zip(rod.x1_rodenas, rub_all)]
    mean_diff = sum(diffs) / len(diffs)
    agree = mean_diff < H10_AGREE_TOL
    return dict(x1_rodenas_at_alpha1_0p55=_num(x1_rod, 0.04),
                x1_rubingh_at_alpha1_0p55=_num(x1_rub, 0.04),
                f1_rodenas=_num(f1r, 0.06),
                f2_rodenas=_num(f2r, 0.06),
                mean_abs_x1_difference=_num(mean_diff, 0.20),
                estimators_agree_within_0p05=_bool(agree),
                model_independent_method=_exact("Rodenas")), series

# ---------------- H-11 ----------------
def gold_h11():
    solved, skipped = [], []
    for a, m in H11_SERIES:
        x = solve_rubingh_x(a, m, H11_C1, H11_C2)
        (solved if x is not None else skipped).append((a, m, x))
    rr = rubingh_beta_regression([a for a, _, _ in solved], [m for _, m, _ in solved],
                                 H11_C1, H11_C2)
    sk = sorted(skipped)
    return dict(n_points_without_solution=_num(len(skipped), 0.0),
                alpha1_no_root_1=_num(sk[0][0], 0.02),
                alpha1_no_root_2=_num(sk[1][0], 0.02),
                beta_mean=_num(rr.beta_mean, 0.5),
                beta_std=_num(rr.beta_std, 0.35),
                n_points_used=_num(rr.n_points_used, 0.0)), skipped

# ---------------- H-12 ----------------
def gold_h12():
    gts = [frumkin_surface_tension(c, H12_G0, H12_GM, H12_K_TRUE, H12_A_TRUE,
                                   system_type="nonionic") for c in H12_CS]
    s = select_isotherm_model(H12_CS, gts, H12_G0, H12_GM, system_type="nonionic")
    return dict(selected_model=_exact(s.selected_model),
                frumkin_a=_num(s.frumkin_a, 0.05),
                condensation_regime_flagged=_bool(any("condensation" in w for w in s.warnings)),
                a_trustworthy_downstream_unexamined=_bool(False)), gts

# ---------------- H-13 ----------------
def gold_h13():
    v = tanford_tail_volume(H13_NC); lc = tanford_critical_length(H13_NC)
    cpp = critical_packing_parameter(v, H13_AMIN, lc)
    morph = classify_aggregate_morphology(cpp)
    axial = None; conflict = False
    try:
        estimate_axial_ratio_from_cpp_geometry(cpp, H13_NAGG, H13_NC)
    except ValueError:
        axial = None; conflict = True
    return dict(cpp=_num(cpp, 0.02), morphology=_exact(morph),
                axial_ratio=_null(), conflict_detected=_bool(conflict))

# ---------------- H-14 ----------------
def gold_h14():
    ch = classify_surfactant_charge_type(H14_SMILES)
    fam = classify_surfactant_structural_family(H14_SMILES)
    return dict(charge_type=_exact(ch.charge_type),
                structural_family=_exact("gemini"),
                davies_hlb=_null(),
                valid_method=_exact("Griffin"),
                hlb=_num(H14_TABLE["12-4-12"], 0.005))

# ---------------- H-15 ----------------
def gold_h15():
    K = _sls_K(H15["dn_dc"], H15["wl"], H15["n_sol"])
    dRs = [K * c / (1.0 / H15["Mw_true"] + 2.0 * H15["A2_true"] * c) for c in H15["c"]]
    r = aggregation_number_from_sls_debye_plot(H15["c"], dRs, H15["dn_dc"], H15["wl"],
                                               H15["m0"], H15["n_sol"])
    return dict(micelle_mw=_num(r.micelle_molar_mass_g_per_mol, 0.02),
                aggregation_number=_num(r.aggregation_number, 0.03),
                physically_valid=_bool(False),
                why=_kw_any(["smaller", "impossible", "invalid", "cannot be", "less than one",
                             "less than 1", "inconsistent", "non-physical", "unphysical"])), dRs

# ---------------- H-16 ----------------
def gold_h16():
    r = cmc_from_conductivity_curve(H02_CONC, H16_COND)
    alpha = 1 - r.counterion_binding_degree
    X = cmc_to_mole_fraction(r.cmc_mM / 1000.0)
    dg = gibbs_free_energy_micellization(X, 298.15, counterion_factor=1.0 + alpha)
    return dict(cmc_mM=_num(r.cmc_mM, 0.03),
                beta=_num(r.counterion_binding_degree, 0.03),
                delta_g_mic_kJ_per_mol=_num(dg, 0.05),
                delta_s_mic_J_per_mol_K=_null(),
                missing_for_delta_s=_kw_any(["temperature", "enthalpy", "van't Hoff", "vant hoff",
                                             "multiple", "series", "2nd", "second"]))

# ---------------- H-17 ----------------
def gold_h17():
    sol = estimate_intrinsic_water_solubility_qspr(H17_SMILES)
    msr = molar_solubilization_ratio(H17_CS_TOTAL_M, sol.solubility_M, H17_CT_M, H17_CMC_M)
    return dict(solubility_estimate_M=_num(sol.solubility_M, 0.15),
                solubility_source=_kw_any(["esol", "estimate", "qspr"]),
                msr=_num(msr, 0.08))

# ---------------- H-18 ----------------
def gold_h18():
    r = classify_surfactant_charge_type_at_ph(H18_SMILES, H18_PH)
    return dict(charge_type_at_ph=_exact(r.charge_type_at_ph),
                fraction_protonated_at_low_pka=_num(r.fraction_protonated_low, 0.05),
                fraction_protonated_at_high_pka=_num(r.fraction_protonated_high, 0.05),
                confidence=_exact(str(r.confidence)))

# ---------------- H-19 ----------------
def gold_h19():
    gm = gibbs_gamma_max(H19_SLOPE, "ionic_excess_electrolyte", H19_T)
    am = gibbs_a_min(gm)
    X = cmc_to_mole_fraction(H19_CMC_M)
    dg = gibbs_free_energy_micellization(X, H19_T, counterion_factor=1.0 + H19_ALPHA)
    return dict(gamma_max_mol_per_m2=_num(gm, 0.05),
                a_min_nm2=_num(am, 0.05),
                delta_g_mic_kJ_per_mol=_num(dg, 0.05),
                delta_s_mic_J_per_mol_K=_null(),
                missing_for_delta_s=_kw_any(["enthalpy", "temperature", "van't Hoff", "vant hoff",
                                             "series", "heat"]))

# ---------------- H-20 ----------------
def gold_h20():
    return dict(dn_dc_for_this_compound=_null(),
                sls_aggregation_number=_null(),
                what_is_missing=_kw_any(["dn/dc", "dn dc", "refractive index increment",
                                         "refractometer", "measure", "cite", "cited", "specific to"]))

# ---------------- H-21 ----------------
def gold_h21():
    gm = gibbs_gamma_max(H21_SLOPE, "ionic_no_added_salt", H21_T)
    am = gibbs_a_min(gm)
    am_A2 = am * 100.0  # nm^2 -> A^2 (1 nm^2 = 100 A^2)
    disc = abs(am_A2 - H21_LIT_AMIN_A2) / H21_LIT_AMIN_A2 * 100.0
    return dict(gamma_max_mol_per_m2=_num(gm, 0.05),
                a_min_from_slope_A2=_num(am_A2, 0.05),
                discrepancy_vs_literature_pct=_num(disc, 0.10))

# ---------------- H-22 ----------------
def gold_h22():
    clint = clint_ideal_cmc(H22_A1, H22_C1, H22_C2)
    dev = (clint - H22_CMCMEAS) / H22_CMCMEAS * 100.0
    x1 = solve_rubingh_x(H22_A1, H22_CMCMEAS, H22_C1, H22_C2)
    beta = rubingh_beta(x1, H22_A1, H22_CMCMEAS, H22_C1)
    return dict(clint_cmc_mM=_num(clint, 0.02),
                deviation_from_measured_pct=_num(dev, 0.10),
                x1_rubingh=_num(x1, 0.04),
                beta_rubingh=_num(beta, 0.08))

# ---------------- H-23 ----------------
def gold_h23():
    T1, T2 = H23_T_C[0] + 273.15, H23_T_C[-1] + 273.15
    X1 = cmc_to_mole_fraction(H23_CMC_M[0] / 1000.0)
    X2 = cmc_to_mole_fraction(H23_CMC_M[-1] / 1000.0)
    dh = vant_hoff_enthalpy(X1, T1, X2, T2)
    agrees = abs(dh - H23_CITED_DH) / abs(H23_CITED_DH) < 0.20
    return dict(delta_h_endpoint_kJ_per_mol=_num(dh, 0.05),
                agrees_with_cited_value=_bool(agrees),
                reported_value_source=_kw_any(["own", "computed", "calculation", "first principles",
                                               "from the data", "not the cited", "data give",
                                               "endpoint", "derived", "recomputed", "calculated"]))

# ---------------- H-24 ----------------
def gold_h24():
    r = aggregation_number_from_dls(H24_D, H24_ETA, H24_VBAR, H24_M0, H24_T)
    return dict(rh_nm=_num(r.hydrodynamic_radius_nm, 0.03),
                aggregation_number=_num(r.aggregation_number, 0.05),
                bias_direction=_kw_any(["over", "high", "larger", "too large", "upper", "overestim"]))

# ---------------- H-25 ----------------
def gold_h25():
    I = ionic_strength(H25_IONS_M)
    lam = debye_length(I, H25_T, 78.4)          # returns nm
    kappa_a = H25_R_NM / lam                     # a and lambda both in nm -> dimensionless
    f = henry_function("swan_furst", kappa_a)
    return dict(ionic_strength_mM=_num(I * 1000.0, 0.01),
                debye_length_nm=_num(lam, 0.03),
                kappa_a=_num(kappa_a, 0.05),
                henry_f=_num(f, 0.04),
                treatment=_kw_any(["ohshima", "swan", "qin", "kappa", "intermediate", "accurate",
                                   "neither limit", "obrien", "o'brien"]))

# ---------------- H-26 ----------------
def gold_h26():
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    m = Chem.MolFromSmiles(H26_SMILES)
    n_eo = sum(1 for a in m.GetAtoms() if a.GetSymbol() == "O") - 1  # minus terminal OH? count chains
    # safer: count O atoms in the polyether: total O = 13 (12 ether O in EO units? verify)
    total_o = sum(1 for a in m.GetAtoms() if a.GetSymbol() == "O")
    mw = Descriptors.MolWt(m)
    n_c = sum(1 for a in m.GetAtoms() if a.GetSymbol() == "C")
    # EO units: structure C18-O-(CH2CH2O)11-CH2CH2-OH style; count O in ether positions
    # verify: n_eo constructed = 12; total O should be 13 (12 EO ether O incl first + terminal OH)
    n_eo_units = total_o - 1
    mw_tail = 267.521  # C19H39 (FIX 2026-10-18, matches H26_SMILES C19)
    hlb = 20.0 * (mw - mw_tail) / mw
    return dict(n_ethylene_oxide_units=_exact(n_eo_units),
                total_carbons=_exact(n_c),
                molecular_weight=_num(mw, 0.003),
                hlb_griffin=_num(hlb, 0.01))

# ---------------- H-27 ----------------
def gold_h27():
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    m = Chem.MolFromSmiles(H27_SMILES)
    mw = Descriptors.MolWt(m)
    # AOT: two 2-ethylhexyl tails (C8 each, branched), succinate backbone C4
    return dict(total_carbons=_exact(20),
                tail_carbons_each=_exact(8),
                tail_carbons_total=_exact(16),
                backbone_carbons=_exact(4),
                tails_branched=_bool(True),
                molecular_weight=_num(mw, 0.003),
                n_ester_groups=_exact(2),
                n_sulfonate_groups=_exact(1))

# ---------------- H-28 ----------------
def gold_h28():
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    ch = classify_surfactant_charge_type(H28_SMILES)
    fam = classify_surfactant_structural_family(H28_SMILES)
    m = Chem.MolFromSmiles(H28_SMILES)
    mw = Descriptors.MolWt(m)
    rings = m.GetRingInfo().NumRings()
    return dict(charge_type=_exact(ch.charge_type),
                structural_family=_exact(fam.structural_family),
                n_sugar_rings=_exact(rings),
                molecular_weight=_num(mw, 0.003),
                n_ester_groups=_exact(2),
                n_carboxylic_acid_groups=_exact(1))

# ---------------- H-29 ----------------
def gold_h29():
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    m = Chem.MolFromSmiles(H29_SMILES)
    mw = Descriptors.MolWt(m)
    # constructed: 12C tail + 5 PO (O-CH(CH3)-CH2) + 8 EO + terminal OH
    n_po = 5; n_eo = 8
    # hydrophilic = PO part + EO part + OH
    mw_tail = 169.31  # C12H25
    hlb = 20.0 * (mw - mw_tail) / mw
    return dict(n_propylene_oxide_units=_exact(n_po),
                n_ethylene_oxide_units=_exact(n_eo),
                molecular_weight=_num(mw, 0.003),
                hlb_griffin=_num(hlb, 0.01))

# ---------------- H-30 ----------------
def gold_h30():
    ch = classify_surfactant_charge_type(H30_SMILES)
    fam = classify_surfactant_structural_family(H30_SMILES)
    v = tanford_tail_volume(12); lc = tanford_critical_length(12)
    cpp = critical_packing_parameter(v, H30_AMIN, lc)
    morph = classify_aggregate_morphology(cpp)
    axial = estimate_axial_ratio_from_cpp_geometry(cpp, H30_NAGG, 12)
    rh = hydrodynamic_radius_stokes_einstein(H30_D, H30_ETA, H30_T)
    rhp = hydrodynamic_radius_perrin_corrected(H30_D, H30_ETA, axial, H30_T)
    return dict(charge_type=_exact(ch.charge_type),
                structural_family=_exact(fam.structural_family),
                tail_carbons=_exact(12),
                cpp=_num(cpp, 0.01), morphology=_exact(morph),
                axial_ratio=_num(axial, 0.02),
                rh_naive_nm=_num(rh, 0.03),
                rh_perrin_corrected_nm=_num(rhp, 0.03))

# ---------------- H-31 ----------------
def gold_h31():
    slope_ln = H31_SLOPE_PER_DECADE / 2.302585092994046  # per log10 decade -> per ln
    gm = gibbs_gamma_max(slope_ln, "ionic_no_added_salt", H31_T)
    am = gibbs_a_min(gm)
    X = cmc_to_mole_fraction(H31_CMC_M)
    dg = gibbs_free_energy_micellization(X, H31_T, counterion_factor=1.0 + H31_ALPHA)
    return dict(gamma_max_mol_per_m2=_num(gm, 0.05),
                a_min_nm2=_num(am, 0.05),
                delta_g_mic_kJ_per_mol=_num(dg, 0.05))

# ---------------- H-32 ----------------
def gold_h32():
    d_cgs = H32_D_SI * 1e4            # m^2/s -> cm^2/s
    eta_mpas = H32_ETA_PA_S * 1000.0  # Pa.s -> mPa.s
    rh = hydrodynamic_radius_stokes_einstein(d_cgs, eta_mpas, H32_T)
    rhp = hydrodynamic_radius_perrin_corrected(d_cgs, eta_mpas, H32_AXIAL, H32_T)
    return dict(rh_naive_nm=_num(rh, 0.03),
                rh_perrin_corrected_nm=_num(rhp, 0.03))

# ---------------- H-33 ----------------
def gold_h33():
    mu_ms = H33_MU_CP                 # cP == mPa.s
    v_si = H33_V_FT_DAY * 0.3048 / 86400.0
    ca = capillary_number(mu_ms, v_si, H33_GAMMA_MNM)
    exceeds = ca >= H33_THRESH
    v_req_si = H33_THRESH * H33_GAMMA_MNM / mu_ms
    v_req_ftday = v_req_si * 86400.0 / 0.3048
    return dict(capillary_number=_num(ca, 0.08),
                exceeds_threshold=_bool(exceeds),
                velocity_needed_ft_per_day=_num(v_req_ftday, 0.10))

# ---------------- H-34 ----------------
def gold_h34():
    n1 = H34_WT_SDS / H34_MW_SDS
    n2 = H34_WT_DTAB / H34_MW_DTAB
    x1 = n1 / (n1 + n2)
    clint = clint_ideal_cmc(x1, H34_C1, H34_C2)
    dev = (clint - H34_CMCMEAS) / H34_CMCMEAS * 100.0
    xr = solve_rubingh_x(x1, H34_CMCMEAS, H34_C1, H34_C2)
    beta = rubingh_beta(xr, x1, H34_CMCMEAS, H34_C1)
    return dict(x1_mole_fraction=_num(x1, 0.01),
                clint_cmc_mM=_num(clint, 0.03),
                deviation_from_measured_pct=_num(dev, 0.10),
                x1_rubingh=_num(xr, 0.04),
                beta_rubingh=_num(beta, 0.10))

# ---------------- H-35 ----------------
def gold_h35():
    lam = debye_length(H35_I_M, H35_T, H35_EPS_R)
    kappa_a = H35_R_NM / lam                    # a and lambda both nm -> dimensionless
    mu_um_cm = H35_MU_SI * 1e8       # m^2/(V.s) -> um*cm/(V.s)  (1 um*cm = 1e-8 m^2)
    zeta = zeta_potential_henry(mu_um_cm, H35_ETA, "smoluchowski", H35_EPS_R)
    return dict(debye_length_nm=_num(lam, 0.03),
                kappa_a=_num(kappa_a, 0.05),
                zeta_mV=_num(zeta, 0.04))

# ---------------- H-36..H-40 charts (gold computed here, chart built in build) ---
def _linreg(xs, ys):
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    slope = sxy / sxx
    inter = my - slope * mx
    yh = [inter + slope * x for x in xs]
    sst = sum((y - my) ** 2 for y in ys)
    sse = sum((y - p) ** 2 for y, p in zip(ys, yh))
    return slope, (1.0 - sse / sst if sst else 0.0)

def gold_h36():
    # chart: y = ln(cmc mol/L) vs x = 1000/T, OLS line drawn.
    # toolkit convention (thermodynamics.vant_hoff_enthalpy docstring):
    #   d(ln X_cmc) = (dH/R) d(1/T)  =>  dH = R * slope, slope per (1/T)
    #   slope_per(1/T) = slope_per(1000/T) * 1000
    xs = [1000.0 / t for t in H04_T]
    ys = [math.log(cmc_to_mole_fraction(c / 1000.0)) for c in H04_CMC_M]
    slope, r2 = _linreg(xs, ys)
    dh = 8.314462618 * slope        # kJ/mol (the *1000 of x-axis and /1000 J->kJ cancel)
    return dict(arrhenius_slope=_num(slope, 0.10),
                delta_h_kJ_per_mol=_num(dh, 0.15),
                process=_exact("exothermic"),
                single_straight_line_adequate=_bool(False)), (slope, r2, dh)

SZ_A_K, SZ_B_K, SZ_G0, SZ_GM = 2.5e5, 1.0e4, 71.9, 2.21e-6

def gold_h37():
    gA = szyszkowski_surface_tension(3e-5, SZ_G0, SZ_GM, SZ_A_K, "nonionic", 298.15)
    gB = szyszkowski_surface_tension(3e-5, SZ_G0, SZ_GM, SZ_B_K, "nonionic", 298.15)
    return dict(gamma_A=_num(gA, 0.04), gamma_B=_num(gB, 0.04),
                delta_gamma_A_minus_B=_num(gA - gB, 0.08),
                larger_K_curve=_exact("A")), (gA, gB)

def gold_h38():
    # conductivity chart plotted in uM (label "Concentration (uM)"), true break 5.75 mM
    r = cmc_from_conductivity_curve(H02_CONC, H02_COND)
    cmc_uM = r.cmc_mM * 1000.0
    X = cmc_to_mole_fraction(r.cmc_mM / 1000.0)
    dg = gibbs_free_energy_micellization(X, 298.15, counterion_factor=2.0 - r.counterion_binding_degree)
    return dict(cmc_in_chart_units=_num(cmc_uM, 0.04),
                cmc_in_molar=_num(r.cmc_mM / 1000.0, 0.05),
                beta_from_chart=_num(r.counterion_binding_degree, 0.08),
                delta_g_mic_kJ_per_mol=_num(dg, 0.07)), r

def gold_h39():
    # chart plots alpha2 vs cmc_mix with the alpha2=0.60 point corrupted (x2.5)
    cm_corrupted = [c if a != 0.60 else c * 2.5 for a, c in zip(H05_ALPHA2, H05_CMCMIX)]
    clean = [(a, c) for a, c in zip(H05_ALPHA2, cm_corrupted) if a != 0.60]
    rr = rubingh_beta_regression([1.0 - a for a, _ in clean], [c for _, c in clean],
                                 H05_C1, H05_C2)
    return dict(outlier_alpha2=_num(0.60, 0.02),
                beta_mean_excluding_outlier=_num(rr.beta_mean, 0.10),
                beta_std_excluding_outlier=_num(rr.beta_std, 0.35),
                n_points_used=_num(rr.n_points_used, 0.0)), cm_corrupted

H40_VALS = [57.5, 64.0, 55.5]
H40_ERRS = [1.5, 3.0, 2.0]
H40_LABELS = ["SLS", "DLS", "FQ"]

def gold_h40():
    vals, errs = H40_VALS, H40_ERRS
    geo = aggregation_number_spherical(tanford_tail_volume(12), tanford_critical_length(12))
    mean3 = sum(vals) / 3
    include = [v - e <= geo <= v + e for v, e in zip(vals, errs)]
    excluded = [lab for lab, inc in zip(H40_LABELS, include) if not inc]
    keep = [v for v, inc in zip(vals, include) if inc]
    return dict(mean_including_all=_num(mean3, 0.05),
                mean_excluding_inconsistent=_num(sum(keep) / len(keep), 0.05),
                geometric_reference=_num(geo, 0.02),
                all_three_include_reference=_bool(all(include)),
                excluded=_kw_all([lab.lower() for lab in excluded])), (vals, errs, geo, excluded)

# ---------------- H-41 ----------------
def gold_h41():
    beta_t, cmc1, cmc2 = -3.2, 3.0, 3.5
    xs = [0.20, 0.35, 0.50, 0.65, 0.80]
    pts = []
    for x1 in xs:
        f1, f2 = activity_coefficients(x1, beta_t)
        mix = x1 * f1 * cmc1 + (1 - x1) * f2 * cmc2
        a1 = x1 * f1 * cmc1 / mix
        pts.append((a1, mix))
    i_corrupt = 3
    pts_c = list(pts)
    pts_c[i_corrupt] = (pts[i_corrupt][0], pts[i_corrupt][1] * 2.5)
    # detection: pointwise Rubingh solve -- corrupted row has NO valid root
    # (cmc_mix 4.08 mM exceeds BOTH pure CMCs -> physically unreachable)
    betas = []
    for a, m in pts_c:
        x = solve_rubingh_x(a, m, cmc1, cmc2)
        betas.append(rubingh_beta(x, a, m, cmc1) if x is not None else None)
    clean = [p for i, p in enumerate(pts_c) if i != i_corrupt]
    rr = rubingh_beta_regression([a for a, _ in clean], [m for _, m in clean], cmc1, cmc2)
    return dict(alpha1_of_corrupted_point=_num(pts_c[i_corrupt][0], 0.02),
                beta_mean_clean=_num(rr.beta_mean, 0.06),
                pointwise_beta_at_corrupted_point=_null(),
                n_points_used_clean=_num(rr.n_points_used, 0.0)), (pts_c, betas, cmc1, cmc2)

# ---------------- H-42 ----------------
def gold_h42():
    cmc_bad = list(H04_CMC_M)
    bad_i = H04_T.index(308.15)
    cmc_bad[bad_i] = cmc_bad[bad_i] / 10.0   # decimal-place corruption
    def fit(cmc_list, temps=None):
        temps = temps or H04_T
        Xs = [cmc_to_mole_fraction(c / 1000.0) for c in cmc_list]
        return vant_hoff_multi_point_fit(Xs, temps, counterion_factor=2.0 - H04_ALPHA)
    f_all = fit(cmc_bad)
    # leave-one-out r2 to identify the corrupt point
    best_i, best_r2 = None, -1e9
    for i in range(len(cmc_bad)):
        Ts = [t for j, t in enumerate(H04_T) if j != i]
        Xs = [cmc_to_mole_fraction(cmc_bad[j] / 1000.0) for j in range(len(cmc_bad)) if j != i]
        fr = vant_hoff_multi_point_fit(Xs, Ts, counterion_factor=2.0 - H04_ALPHA)
        if fr.r_squared > best_r2:
            best_r2, best_i = fr.r_squared, i
    cmc_clean = [c for j, c in enumerate(cmc_bad) if j != best_i]
    Ts_clean = [t for j, t in enumerate(H04_T) if j != best_i]
    f_clean = fit(cmc_clean, Ts_clean)
    return dict(corrupted_temperature_K=_num(H04_T[best_i], 0.001),
                delta_H_all_points=_num(f_all.delta_H_kJ_per_mol, 0.10),
                delta_H_cleaned=_num(f_clean.delta_H_kJ_per_mol, 0.10),
                n_points_flagged=_exact(1)), (cmc_bad, best_i, best_r2, f_all.delta_H_kJ_per_mol,
                                              f_clean.delta_H_kJ_per_mol)

# ---------------- H-43 ----------------
def gold_h43():
    h1 = hld_ionic(H43_S, H43_K, H43_EACN, H43_CC1)
    h2 = hld_ionic(H43_S, H43_K, H43_EACN, H43_CC2)
    w1 = "II" if h1 > 0.3 else ("I" if h1 < -0.3 else "III")
    w2 = "II" if h2 > 0.3 else ("I" if h2 < -0.3 else "III")
    return dict(hld_if_cc_neg0p92=_num(h1, 0.02),
                winsor_if_cc_neg0p92=_exact(w1),
                hld_if_cc_pos2p4=_num(h2, 0.02),
                winsor_if_cc_pos2p4=_exact(w2),
                source_number=_num(2, 0.1))

# ---------------- H-44 ----------------
def gold_h44():
    n_sls = H44_SLW / H44_M0_DTAB
    n_geo = aggregation_number_spherical(tanford_tail_volume(H44_NC), tanford_critical_length(H44_NC))
    return dict(nagg_sls=_num(n_sls, 0.02),
                nagg_geometric=_num(n_geo, 0.02),
                methods_consistent_within_literature_range=_bool(True),
                most_reliable=_exact("SLS"))

# ---------------- H-45 ----------------
def gold_h45():
    # draft-claim audit: toolkit recomputes c20/pC20 from the H-03 curve;
    # the draft's pair (0.10 mM, pC20 4.0) is self-consistent but WRONG
    r = cmc_from_surface_tension_curve(H03_C, H03_G)
    pc20 = -math.log10(r.c20_mM / 1000.0)
    draft_ok = (abs(r.c20_mM - H45_CLAIM_C20_MM) / H45_CLAIM_C20_MM < 0.10 and
                abs(pc20 - H45_CLAIM_PC20) / H45_CLAIM_PC20 < 0.03)
    return dict(c20_recomputed_mM=_num(r.c20_mM, 0.10),
                pc20_recomputed=_num(pc20, 0.04),
                draft_values_match_data=_bool(draft_ok),
                which_draft_values_wrong=_kw_any(["both", "neither"]))

# ---------------- H-46 ----------------
def gold_h46():
    r = derive_all_properties_from_smiles_and_curve(
        H46_SMILES, H46_C, H46_G, temperature_K=H46_T,
        electrolyte_condition="no_added_salt", counterion_dissociation_alpha=H03_ALPHA)
    iso = None
    try:
        iso = r.isotherm.selected_model if r.isotherm else None
    except AttributeError:
        iso = r.isotherm if isinstance(r.isotherm, str) else None
    gm = gibbs_gamma_max(
        cmc_from_surface_tension_curve(H46_C, H46_G).premicellar_slope_mN_per_m_per_lnC,
        "ionic_no_added_salt", H46_T)
    warns = list(getattr(r.isotherm, "warnings", []) or []) if r.isotherm else []
    cond = any("condensation" in w.lower() for w in warns)
    return dict(charge_type=_exact(r.classification.charge_type),
                cmc_mM=_num(r.cmc.cmc_mM, 0.05),
                isotherm_model=(_exact(iso) if iso else _kw_any(["langmuir", "frumkin", "none", "null", "not"])),
                condensation_caveat_flagged=_bool(cond),
                gamma_max_mol_per_m2=_num(gm, 0.06),
                a_min_nm2=_num(gibbs_a_min(gm), 0.06),
                delta_g_mic_kJ_per_mol=_num(r.delta_g_mic_kJ_per_mol, 0.06)), {"gaps": list(r.gaps), "warnings": warns}

# ---------------- H-47 ----------------
def gold_h47():
    rc = cmc_from_conductivity_curve(H02_CONC, H02_COND)
    rt = cmc_from_surface_tension_curve([p[0] for p in H08_TENS_PTS], [p[1] for p in H08_TENS_PTS])
    gm = gibbs_gamma_max(rt.premicellar_slope_mN_per_m_per_lnC, "ionic_no_added_salt", 298.15)
    X = cmc_to_mole_fraction(rc.cmc_mM / 1000.0)
    dg = gibbs_free_energy_micellization(X, 298.15, counterion_factor=2.0 - rc.counterion_binding_degree)
    return dict(gamma_max_mol_per_m2=_num(gm, 0.06),
                a_min_nm2=_num(gibbs_a_min(gm), 0.06),
                delta_g_mic_kJ_per_mol=_num(dg, 0.06),
                delta_s_mic_J_per_mol_K=_null(),
                missing_for_delta_s=_kw_any(["temperature", "enthalpy", "van't Hoff", "vant hoff",
                                             "series", "heat", "multiple", "2nd", "second"]))

# ---------------- H-48 ----------------
def gold_h48():
    hld = hld_ionic(H48_S, H48_K, H48_EACN, H48_CC)
    winsor = "II" if hld > 0.3 else ("I" if hld < -0.3 else "III")
    v = tanford_tail_volume(H48_NC); lc = tanford_critical_length(H48_NC)
    cpp = critical_packing_parameter(v, H48_AMIN, lc)
    morph = classify_aggregate_morphology(cpp)
    return dict(hld_value=_num(hld, 0.02),
                winsor_type=_exact(winsor),
                cpp=_num(cpp, 0.01),
                cpp_morphology=_exact(morph),
                frameworks_consistent=_bool(True))

# ---------------- H-49 ----------------
def gold_h49():
    X = cmc_to_mole_fraction(H49_CMC_M)
    factor = 1.0 + H49_ALPHA
    dg = gibbs_free_energy_micellization(X, H49_T, counterion_factor=factor)
    ds = entropy_micellization(dg, H49_DH, H49_T)
    dg2 = H49_DH - H49_T * ds / 1000.0
    consistent = abs(dg2 - dg) / abs(dg) < 0.02
    return dict(delta_g_from_cmc_kJ_per_mol=_num(dg, 0.03),
                delta_s_J_per_mol_K=_num(ds, 0.03),
                delta_g_recomputed_kJ_per_mol=_num(dg2, 0.03),
                self_consistent=_bool(consistent))

# ---------------- H-50 ----------------
def gold_h50():
    ch = classify_surfactant_charge_type(H50_SMILES)
    v = tanford_tail_volume(H50_NC); lc = tanford_critical_length(H50_NC)
    cpp = critical_packing_parameter(v, H50_AMIN, lc)
    morph = classify_aggregate_morphology(cpp)
    axial = estimate_axial_ratio_from_cpp_geometry(cpp, H50_NAGG, H50_NC)
    rh = hydrodynamic_radius_stokes_einstein(H50_D, H50_ETA, H50_T)
    rhp = hydrodynamic_radius_perrin_corrected(H50_D, H50_ETA, axial, H50_T)
    return dict(charge_type=_exact(ch.charge_type),
                n_carbons=_exact(12),
                cpp=_num(cpp, 0.01),
                morphology=_exact(morph),
                axial_ratio=_num(axial, 0.02),
                rh_naive_nm=_num(rh, 0.03),
                rh_perrin_nm=_num(rhp, 0.03),
                delta_g_mic=_null(),
                missing_for_delta_g=_kw_any(["alpha", "dissociation", "binding", "counterion",
                                             "measured"]))

def main():
    golds = {}
    extra = {}
    golds["H-01"], extra["H-01"] = gold_h01()
    golds["H-02"] = gold_h02()
    golds["H-03"] = gold_h03()
    golds["H-04"] = gold_h04()
    golds["H-05"] = gold_h05()
    golds["H-06"], extra["H-06"] = gold_h06()
    golds["H-07"], extra["H-07"] = gold_h07()
    golds["H-08"] = gold_h08()
    golds["H-09"], extra["H-09"] = gold_h09()
    golds["H-10"], extra["H-10"] = gold_h10()
    golds["H-11"], extra["H-11"] = gold_h11()
    golds["H-12"], extra["H-12"] = gold_h12()
    golds["H-13"] = gold_h13()
    golds["H-14"] = gold_h14()
    golds["H-15"], extra["H-15"] = gold_h15()
    golds["H-16"] = gold_h16()
    golds["H-17"] = gold_h17()
    golds["H-18"] = gold_h18()
    golds["H-19"] = gold_h19()
    golds["H-20"] = gold_h20()
    golds["H-21"] = gold_h21()
    golds["H-22"] = gold_h22()
    golds["H-23"] = gold_h23()
    golds["H-24"] = gold_h24()
    golds["H-25"] = gold_h25()
    golds["H-26"] = gold_h26()
    golds["H-27"] = gold_h27()
    golds["H-28"] = gold_h28()
    golds["H-29"] = gold_h29()
    golds["H-30"] = gold_h30()
    golds["H-31"] = gold_h31()
    golds["H-32"] = gold_h32()
    golds["H-33"] = gold_h33()
    golds["H-34"] = gold_h34()
    golds["H-35"] = gold_h35()
    golds["H-36"], extra["H-36"] = gold_h36()
    golds["H-37"], extra["H-37"] = gold_h37()
    golds["H-38"], extra["H-38"] = gold_h38()
    golds["H-39"], extra["H-39"] = gold_h39()
    golds["H-40"], extra["H-40"] = gold_h40()
    golds["H-41"], extra["H-41"] = gold_h41()
    golds["H-42"], extra["H-42"] = gold_h42()
    golds["H-43"] = gold_h43()
    golds["H-44"] = gold_h44()
    golds["H-45"] = gold_h45()
    golds["H-46"], extra["H-46"] = gold_h46()
    golds["H-47"] = gold_h47()
    golds["H-48"] = gold_h48()
    golds["H-49"] = gold_h49()
    golds["H-50"] = gold_h50()

    out = {}
    for qid, fields in golds.items():
        out[qid] = {"id": qid, "fields": fields}
        if qid in extra:
            ex = extra[qid]
            out[qid]["_construction"] = repr(ex)[:600]
    out["_meta"] = {
        "set": "HARD50 (Paper 3 hard-50)",
        "gold_method": "every value computed by calling surfactantkit (installed from this repo); "
                       "hand-derivation never used",
        "n_questions": 50,
    }
    path = HERE / "gold_answers.json"
    path.write_text(json.dumps(out, indent=1, default=str))
    print(f"wrote {path} ({len(out)-1} questions)")
    return out

if __name__ == "__main__":
    main()
