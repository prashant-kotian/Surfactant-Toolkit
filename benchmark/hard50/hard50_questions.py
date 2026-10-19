"""HARD50 question texts -- single source of truth for the prompt bodies.
Data values are the SAME inputs used by gold_build.py (imported where
practical; the few chart/table renderings are transcribed here with the
identical numbers). build_hard50.py assembles question_bank_source.txt.
"""
from __future__ import annotations

Q: dict[str, tuple[str, str]] = {}   # id -> (reference_note, paste_body)

def add(qid: str, note: str, body: str):
    assert qid not in Q, qid
    Q[qid] = (note, body)

# --- fields given as (key, json_value_text) pairs ---------------------------
def final(*fields):
    rows = []
    for i, (k, v) in enumerate(fields):
        comma = "," if i < len(fields) - 1 else ""
        rows.append('  "%s": %s%s' % (k, v, comma))
    return """End your response with a block in exactly this form:

FINAL_JSON:
{
%s
}
<<< PASTE ABOVE <<<""" % "\n".join(rows)

FP = "Work entirely from first principles -- do not search the internet or recall a literature value for this compound or system.\n\n"

# ============================================================ H-01 .. H-10
add("H-01",
 "source: constructed demonstration dataset (forward-generated from Mw=35888 g/mol, A2=4.0e-4, Debye model); "
 "axis: arithmetic-execution (10+ sequential numeric operations); "
 "Version A: aggregation_number_from_sls_debye_plot + mass_action_free_energy_of_micellization; "
 "Version B: walk K -> Debye plot -> Nagg -> dG in order, every constant is supplied; "
 "expected failure: hand arithmetic on K~1e-7 / lambda^4, unit slips in cm-vs-nm, mass-action sign errors "
 "(procedural-execution diagnostic arXiv 2605.00817: mult/div -32pp, under-execution dominates); "
 "ref: LinAlg-Bench 2605.16675 execution-error regime; MedMCP-Calc 2601.23049",
 FP + """Static light-scattering data for a nonionic surfactant (monomer molar mass M0 = 448.6 g/mol) in water at T = 298.15 K, measured at wavelength 632.8 nm. Solvent refractive index n = 1.333, dn/dc = 0.135 mL/g, Avogadro constant NA = 6.02214076e23 /mol.

Debye-plot ordinates (c in g/mL, excess Rayleigh ratio Delta R in cm^-1):
  0.0020, 8.9867e-06
  0.0040, 1.7048e-05
  0.0060, 2.4319e-05
  0.0080, 3.0911e-05
  0.0100, 3.6915e-05
  0.0120, 4.2407e-05

Use the standard Debye treatment: K = 4*pi^2*n^2*(dn/dc)^2/(NA*lambda^4) with lambda in cm, and K*c/DeltaR = 1/Mw + 2*A2*c with A2 = 4.0e-4 mol*mL/g^2. The CMC of this surfactant is 8.7e-5 mol/L.

Report: the optical constant K, the micelle weight-average molar mass Mw (g/mol), the aggregation number (Mw/M0), and the mass-action free energy of micellization
  dG = R*T*[(1 + 1/n)*ln(X) - (1/n)*ln(n)]
with R = 8.314462618 J/(mol*K), X = CMC(mol/L)/55.5, and n = the aggregation number you found (kJ/mol).

""" + final(("optical_constant_K", "<number>"),
            ("micelle_mw", "<number>"),
            ("aggregation_number", "<number>"),
            ("delta_g_mass_action_kJ_per_mol", "<number>")))

add("H-02",
 "source: constructed demonstration dataset (forward-generated, true break 5.75 mM, true beta 0.6562); "
 "axis: arithmetic-execution (regression + convention chain); "
 "Version A: cmc_from_conductivity_curve + gibbs_free_energy_micellization; "
 "Version B: break -> beta -> counterion factor (2-beta) -> dG; "
 "expected failure: slope-ratio arithmetic, silently dropping the counterion factor or assuming nonionic (1.0); "
 "ref: MedMCP-Calc 2601.23049 calculator-avoidance; procedural-execution 2605.00817",
 FP + """Conductivity vs concentration for a 1:1 ionic surfactant in water at 298.15 K (conductivity in mS/cm, concentration in mM):

  c (mM), conductivity (mS/cm)
  1.0, 0.690
  2.0, 1.330
  3.0, 1.970
  4.0, 2.610
  5.0, 3.250
  5.5, 3.570
  6.0, 3.785
  6.5, 3.895
  7.0, 4.005
  8.0, 4.225
  10.0, 4.665
  12.0, 5.105

Report: the CMC (mM) from the break-point, the counterion binding degree beta (= 1 - slope2/slope1), the degree of counterion dissociation alpha (= 1 - beta), and the standard free energy of micellization
  dG = (2 - beta)*R*T*ln(X),  X = CMC(mol/L)/55.5,  R = 8.314462618 J/(mol*K), T = 298.15 K
in kJ/mol (use the beta you derived).

""" + final(("cmc_mM", "<number>"),
            ("beta", "<number>"),
            ("alpha", "<number>"),
            ("delta_g_mic_kJ_per_mol", "<number>")))

add("H-03",
 "source: constructed demonstration 3-regime tensiometry (forward-generated baseline+decline+plateau); "
 "axis: arithmetic-execution (regime detection + Gibbs chain); "
 "Version A: cmc_from_surface_tension_curve + gibbs_gamma_max + gibbs_a_min + gibbs_free_energy_micellization; "
 "Version B: detect the baseline regime first, then the Gamma_max chain; state the electrolyte condition; "
 "expected failure: silently missing the baseline (changes the premicellar slope), Gibbs prefactor n=1 vs n=2 substitution; "
 "ref: AtmosCoder 2608.18726 canonical-method reuse when task conditions invalidate it",
 FP + """Surface tension isotherm for a 1:1 ionic surfactant in water with NO added electrolyte at T = 298.15 K (constructed demonstration dataset; gamma in mN/m, c in mM):

  0.005, 71.9
  0.010, 71.8
  0.020, 71.6
  0.05, 66.5
  0.10, 61.9
  0.30, 53.0
  0.70, 45.0
  1.50, 38.0
  3.00, 33.0
  6.00, 31.2
  10.00, 31.1
  16.00, 31.0

Report: the CMC (mM); whether the curve requires a three-segment (flat baseline / decline / plateau) description (true/false) and the number of points in the baseline segment; the maximum surface excess Gamma_max (mol/m^2) from the premicellar slope; the minimum area per molecule A_min (nm^2); and dG_mic (kJ/mol) using the counterion dissociation degree alpha = 0.31 via dG = (1 + alpha)*R*T*ln(X) with X = CMC(mol/L)/55.5, T = 298.15 K. State clearly which Gibbs prefactor n you used and why (the electrolyte condition is given above).

""" + final(("cmc_mM", "<number>"),
            ("three_segment_selected", "<true or false>"),
            ("n_baseline_points", "<number>"),
            ("gamma_max_mol_per_m2", "<number>"),
            ("a_min_nm2", "<number>"),
            ("delta_g_mic_kJ_per_mol", "<number>")))

add("H-04",
 "source: REAL Fu, Fan, Cao & Wang 2019 RSC Adv. 9 28799 Table data for [C12mpy][Br] CMC-vs-T (same 8-point table as GZ-26/28/50/51), "
 "beta=0.268 real; axis: arithmetic-execution (8-row multi-point fit with curvature); "
 "Version A: vant_hoff_multi_point_fit (returns dH, dS, dCp, dG(T)); "
 "Version B: build dG(T) rows with (2-beta), then allow a dCp term; do not use a two-point shortcut; "
 "expected failure: two-point substitution, 1/T orientation sign slips (the toolkit's own docstring records this "
 "exact sign bug caught in an unaugmented Claude Opus 4.8 pilot), entropy-vs-enthalpy misjudgment; "
 "ref: AtmosCoder sign-convention carry-over; GZ-26 (universal FAIL recipe on this data)",
 FP + """Real CMC-vs-temperature data for the cationic surfactant [C12mpy][Br] (1-dodecyl-3-methylpyridinium bromide); counterion binding degree beta = 0.268 (real value from the same study):

  T (degC), CMC (mM)
  10, 9.93
  15, 9.82
  20, 9.76
  25, 9.74
  30, 9.77
  35, 10.02
  40, 10.54
  45, 11.47

Perform a multi-point van't Hoff analysis of this data that allows a nonzero heat capacity change (i.e. curvature in the Arrhenius plot): compute dG_mic(T) at every row with dG = (2 - beta)*R*T*ln(CMC[mol/L]/55.5), then extract dH (kJ/mol), dS (J/(mol*K)), and dCp (J/(mol*K)) referenced to 298.15 K, plus dG_mic at 298.15 K (kJ/mol). Finally state whether micellization over this range is entropy-driven or enthalpy-driven (answer: exactly one of "entropy-driven" / "enthalpy-driven").

""" + final(("delta_H_kJ_per_mol", "<number>"),
            ("delta_S_J_per_mol_K", "<number>"),
            ("delta_Cp_J_per_mol_K", "<number>"),
            ("delta_G_298_kJ_per_mol", "<number>"),
            ("driver", '"entropy-driven" or "enthalpy-driven"')))

add("H-05",
 "source: REAL Vautier-Giongo et al. 2005 JCIS 282 149 DHPC-DTAB binary system (same real series family as GZ-02); "
 "axis: arithmetic-execution (pooled regression + single-point evaluations); "
 "Version A: rubingh_beta_regression + solve_rubingh_x + activity_coefficients + excess_free_energy; "
 "Version B: pooled beta first, then evaluate at alpha2=0.60; mind alpha2 vs alpha1 indexing; "
 "expected failure: indexing the wrong component, confusing per-point and pooled beta, f1/f2 swap; "
 "ref: procedural-execution 2605.00817; GZ-02 corrected-gold precedent",
 FP + """Real measured mixed CMC for the binary system DHPC (zwitterionic phospholipid, component 1, pure CMC 1.8 mM) + DTAB (cationic, component 2, pure CMC 15.6 mM). Data: bulk mole fraction of component 2 (alpha2) and mixed CMC (mM):

  alpha2, CMC_mix (mM)
  0.95, 10.8
  0.90, 6.2
  0.80, 4.8
  0.60, 2.7
  0.50, 2.4

Report: the pooled Rubingh interaction parameter beta (mean over all five points, from the regular-solution mass-balance fit) and its sample standard deviation (n-1 denominator); the micellar mole fraction of component 1, x1, at alpha2 = 0.60; the activity coefficients f1 and f2 at that composition computed with the pooled beta; and the molar excess Gibbs energy G_exc at that composition (kJ/mol, T = 298.15 K).

""" + final(("beta_mean", "<number>"),
            ("beta_std", "<number>"),
            ("x1_at_alpha2_0p60", "<number>"),
            ("f1", "<number>"),
            ("f2", "<number>"),
            ("gexc_kJ_per_mol", "<number>")))

add("H-06",
 "source: constructed demonstration series (forward-generated from asymmetric Margules truth W12=-3.5, W21=-0.9, 0.5% jitter, seed 7); "
 "axis: reasoning-decision (model-selection race) + arithmetic (BIC, parameter attribution); "
 "Version A: select_mixture_model (returns selected_model, per-candidate BIC/params); "
 "Version B: report the SELECTED model's own parameters -- do not attribute Rubingh's beta to the winner; "
 "expected failure: parameter-attribution error after a BIC race; hand-BIC arithmetic drift; "
 "ref: AtmosCoder 2608.18726 wrong-formula-same-task; GZ-02/25 EOMMM context",
 FP + """Constructed demonstration mixed-micelle series (component 1 pure CMC 1.8 mM, component 2 pure CMC 8.3 mM; bulk mole fraction of component 1):

  alpha1, CMC_mix (mM)
  0.10, 3.3751
  0.25, 2.5588
  0.40, 2.2482
  0.55, 2.0694
  0.70, 1.9853
  0.85, 1.9065

Fit candidate mixed-micelle models to this series and compare them by BIC. Report: the selected model identifier (name it exactly as the fitting routine does, e.g. its own model string), the selected model's BIC, the pooled Rubingh beta (from the Rubingh candidate, regardless of which model wins), whether the interaction is synergistic or antagonistic (answer containing "synerg" or "antagon"), and the two asymmetric-Margules parameters W12 and W21 of the EOMMM candidate (report them even if EOMMM is not the selected model).

""" + final(("selected_model", '"<model identifier string>"'),
            ("selected_bic", "<number>"),
            ("rubingh_beta_mean", "<number>"),
            ("interaction_character", '"<synergistic or antagonistic>"'),
            ("w12", "<number>"),
            ("w21", "<number>")))

add("H-07",
 "source: constructed demonstration isotherm (forward-generated Langmuir truth with 18% noise, seed 11); "
 "axis: reasoning-decision (BIC and raw-RSS disagree -- report both, pick by BIC); "
 "Version A: select_isotherm_model (selected_model, delta_bic, frumkin_a, warnings); "
 "Version B: BIC is the decision rule; raw RSS will favour the nested extra-parameter model -- report that fact honestly too; "
 "expected failure: rationalizing the choice around whichever statistic favours the preferred model; fabricated justification "
 "(Zvi Mowshowitz Opus-4.8 review: rationalizing around explicit constraints); "
 "ref: MatSciBench 2510.12171 false critiques/self-correction failure",
 FP + """Constructed demonstration surface-tension isotherm (gamma0 = 72.0 mN/m, Gamma_max = 2.2e-6 mol/m^2 supplied for the isotherm models; gamma in mN/m, c in mol/L):

  c, gamma
  3.0e-06, 62.69
  8.0e-06, 58.03
  1.5e-05, 54.94
  3.0e-05, 51.11
  5.0e-05, 48.38
  7.5e-05, 46.23
  1.0e-04, 44.52
  1.3e-04, 43.22

Compare a Langmuir and a Frumkin fit of this isotherm using BIC. Report: the selected model (exact identifier), delta_BIC defined as BIC(Langmuir) - BIC(Frumkin) (number, sign matters), the fitted Frumkin interaction parameter a, whether that fitted a represents a physically meaningful interaction (true/false), and whether the RAW (unpenalized) residual sum of squares favours Frumkin over Langmuir (true/false).

""" + final(("selected_model", '"<model identifier>"'),
            ("delta_bic", "<number>"),
            ("frumkin_a", "<number>"),
            ("physically_meaningful_interaction", "<true or false>"),
            ("raw_rss_frumkin_lower_than_langmuir", "<true or false>")))

add("H-08",
 "source: constructed demonstration pair (conductivity true break 5.75 mM / tension on exact line 34.8-7*ln(c/5.75), both built to agree at 5.75); "
 "axis: arithmetic-execution (two-method race with numeric agreement threshold); "
 "Version A: cmc_from_conductivity_curve + cmc_from_surface_tension_curve + gibbs_free_energy_micellization; "
 "Version B: compute both breaks, then the relative difference, THEN judge -- agreement is a number, not a vibe; "
 "expected failure: declaring agreement/disagreement without computing against the 2% threshold; geometric-mean bracket read errors; "
 "ref: procedural-execution 2605.00817 (under-execution); llm-stats HLE 46.2% universal-fail precedent",
 FP + """Two independent measurements on the SAME 1:1 ionic surfactant system at 298.15 K.

(1) Conductivity (mS/cm) vs concentration (mM):
  1.0, 0.690
  2.0, 1.330
  3.0, 1.970
  4.0, 2.610
  5.0, 3.250
  5.5, 3.570
  6.0, 3.785
  6.5, 3.895
  7.0, 4.005
  8.0, 4.225
  10.0, 4.665
  12.0, 5.105

(2) Surface tension (mN/m) vs concentration (mM):
  0.06, 66.7
  0.15, 60.3
  0.40, 53.5
  0.90, 47.8
  2.00, 42.2
  4.00, 37.3
  5.00, 35.8
  6.60, 34.6
  9.00, 34.5
  12.00, 34.4
  16.00, 34.3

Report: CMC from conductivity (mM), CMC from tensiometry (mM), whether the two methods agree within 2% relative difference (true/false -- compute the relative difference first), the counterion binding degree beta from the conductivity slopes, and dG_mic (kJ/mol) with dG = (2 - beta)*R*T*ln(CMC_cond[mol/L]/55.5), T = 298.15 K, R = 8.314462618 J/(mol*K).

""" + final(("cmc_conductivity_mM", "<number>"),
            ("cmc_tensiometry_mM", "<number>"),
            ("agree_within_2pct", "<true or false>"),
            ("beta", "<number>"),
            ("delta_g_mic_kJ_per_mol", "<number>")))

add("H-09",
 "source: constructed demonstration HLD salinity scan (forward-generated truth k=0.16, Cc=-1.4, alpha_dT=0.01 K^-1); "
 "axis: arithmetic-execution + reasoning (ln-linearize, round-trip S*, Winsor sign mapping, achievability); "
 "Version A: fit_k_and_cc_from_salinity_scan + optimal_salinity_ionic + hld_ionic; "
 "Version B: HLD = lnS - k*EACN - alpha*dT + Cc; S* solves HLD=0; the Winsor mapping is by SIGN of HLD -- do not invert it; "
 "expected failure: inverted Winsor mapping (HLD<0 is Winsor I, direct/O-in-W), forgetting the alpha*dT term, "
 "ignoring the NaCl solubility ceiling (~26 wt%) when judging achievability; "
 "ref: AtmosCoder 2608.18726 sign-convention carry-over after direction change; hld.py source docstring",
 FP + """Constructed demonstration optimal-salinity scan at reference temperature (delta_T = 0): the phase-behavior optimum S* (wt% NaCl) as a function of oil EACN:

  EACN, S* (wt%)
  8, 14.585
  10, 20.086
  12, 27.660
  14, 38.092

For ionic-surfactant HLD the governing relation is
  HLD = ln(S) - k*EACN - alpha_dT*delta_T + Cc
with S in wt% (ln of the number, e.g. ln(14.585)), and S* is the salinity at which HLD = 0. The temperature coefficient is alpha_dT = 0.01 K^-1.

Report: the fitted k; the fitted Cc; the r^2 of the ln(S*)-vs-EACN linear fit; S* (wt%, as a number) for EACN = 16 at delta_T = +10 K; the HLD at S = 3 wt% for that same condition (EACN 16, delta_T +10 K); the Winsor phase type at that HLD (answer exactly one of "I", "II", "III" -- recall HLD > 0 favours reverse W-in-O structures, HLD < 0 favours direct O-in-W structures, |HLD| <= 0.3 is the balanced middle); and whether the required S* for EACN 16 / +10 K is achievable in water (true/false; recall salts saturate near 26 wt% in water).

""" + final(("k", "<number>"),
            ("cc", "<number>"),
            ("fit_r2", "<number>"),
            ("s_star_EACN16_dT10_pct", "<number>"),
            ("hld_at_S3", "<number>"),
            ("winsor_type", '"I" or "II" or "III"'),
            ("s_star_achievable_in_water", "<true or false>")))

add("H-10",
 "source: constructed demonstration series (forward-generated from asymmetric-Margules truth W12=-3.0, W21=-1.2, 1% jitter, seed 13); "
 "axis: arithmetic-execution (local numerical differentiation + implicit roots) + reasoning (which estimator is model-independent); "
 "Version A: rodenas_x1_series + solve_rubingh_x + rodenas_activity_coefficients; "
 "Version B: Rodenas needs the local slope of ln(CMC) vs alpha (fit neighboring points); Rubingh needs a per-point implicit root -- "
 "compute both before judging agreement from the MEAN absolute difference; "
 "expected failure: hand numerical differentiation of a noisy series, per-point transcendental solves, "
 "declaring disagreement from one point instead of the mean; "
 "ref: procedural-execution 2605.00817; Rodenas-vs-Rubingh model-independence (GZ-02 axis, now numeric)",
 FP + """Constructed demonstration mixed-CMC series (component 1 pure CMC 2.4 mM, component 2 pure CMC 11.0 mM):

  alpha1, CMC_mix (mM)
  0.10, 4.6879
  0.25, 3.5225
  0.40, 3.0273
  0.55, 2.7678
  0.70, 2.5746
  0.85, 2.4768

Two independent estimators of the micellar composition exist: (A) Rodenas' model-independent Gibbs-Duhem route, x1 = -alpha1*(1-alpha1)*d(ln CMC_mix)/d(alpha1) + alpha1, which needs the LOCAL slope of the series; and (B) the Rubingh regular-solution route, which solves the implicit mass-balance equation at each point.

Report, at alpha1 = 0.55: x1 from Rodenas, x1 from Rubingh; the activity coefficients f1 = alpha1*CMC_mix/(x1*CMC1) and f2 = (1-alpha1)*CMC_mix/((1-x1)*CMC2) using the RODENAS x1 at that composition; the mean absolute difference between the two x1 estimates across all six compositions; whether the two estimators agree within 0.05 on that mean absolute difference (true/false); and which of the two methods is model-independent (answer exactly "Rodenas" or "Rubingh").

""" + final(("x1_rodenas_at_alpha1_0p55", "<number>"),
            ("x1_rubingh_at_alpha1_0p55", "<number>"),
            ("f1_rodenas", "<number>"),
            ("f2_rodenas", "<number>"),
            ("mean_abs_x1_difference", "<number>"),
            ("estimators_agree_within_0p05", "<true or false>"),
            ("model_independent_method", '"Rodenas" or "Rubingh"')))

# =========================================================== H-11 .. H-25
add("H-11",
 "source: constructed demonstration series engineered so two points have NO Rubingh root (mix CMC above what the model allows); "
 "axis: reasoning-decision (infeasibility detection) + arithmetic (subset regression); "
 "Version A: solve_rubingh_x (returns null on no-root) + rubingh_beta_regression on the feasible subset; "
 "Version B: scan EVERY point for a root first -- the regression must skip infeasible points, not crash or fabricate; "
 "expected failure: tool crash handled by silently averaging, inventing roots, or wrong subset; "
 "ref: GZ infeasibility-detection recipe (universal FAIL); AtmosCoder trap-task family",
 FP + """Constructed demonstration mixed-CMC series (component 1 pure CMC 8.0 mM, component 2 pure CMC 15.0 mM):

  alpha1, CMC_mix (mM)
  0.10, 9.0
  0.25, 12.0
  0.50, 15.0
  0.70, 16.0
  0.85, 13.0

Attempt a Rubingh analysis at EVERY composition: solve the regular-solution mass balance for x1 and check that a physical root exists (0 < x1 < 1). Report: how many of the five points have NO valid solution (number); the alpha1 values of those points (two numbers, ascending); then, pooling only the VALID points, the Rubingh beta mean, its sample standard deviation (n-1 denominator), and the number of points used in that regression.

""" + final(("n_points_without_solution", "<number>"),
            ("alpha1_no_root_1", "<number>"),
            ("alpha1_no_root_2", "<number>"),
            ("beta_mean", "<number>"),
            ("beta_std", "<number>"),
            ("n_points_used", "<number>")))

add("H-12",
 "source: constructed demonstration isotherm (forward-generated Frumkin truth K=2.5e5, a=2.6, Gamma_max=2.1e-6, seed none); "
 "axis: reasoning-decision (thermodynamic-regularity flag) + arithmetic (nonlinear fit); "
 "Version A: select_isotherm_model (selected_model, frumkin_a, warnings mentioning condensation); "
 "Version B: a >= 2 is the 2D-condensation regime -- the flag and the trust judgment are separate asks; "
 "expected failure: hand-fit precision failing on a; ignoring the condensation warning and trusting downstream quantities; "
 "ref: Frumkin condensation criterion (standard); GZ isotherm-selection recipes",
 FP + """Constructed demonstration surface-tension isotherm (gamma0 = 71.8 mN/m, Gamma_max = 2.1e-6 mol/m^2 supplied for the isotherm models; gamma in mN/m, c in mol/L):

  c, gamma
  2.0e-06, 61.81
  5.0e-06, 57.08
  1.0e-05, 53.48
  2.0e-05, 49.88
  4.0e-05, 46.28
  7.0e-05, 43.36
  1.1e-04, 41.01

Fit Langmuir and Frumkin isotherm models to this data and select between them by BIC. Report: the selected model identifier; the fitted Frumkin interaction parameter a; whether the selected fit lies in the two-dimensional-condensation regime (true/false -- the regularity criterion for the Frumkin isotherm); and whether its downstream Gamma-derived thermodynamics can be trusted WITHOUT inspecting the fit residuals first (true/false).

""" + final(("selected_model", '"<model identifier>"'),
            ("frumkin_a", "<number>"),
            ("condensation_regime_flagged", "<true or false>"),
            ("a_trustworthy_downstream_unexamined", "<true or false>")))

add("H-13",
 "source: real geometric inconsistency (C12 chain: tanford v=350.2 A^3, lc=16.68 A with Amin=50 A^2 puts CPP in the rod regime, "
 "yet Nagg=20 implies a core semi-axis smaller than the pinned minor axis -- the toolkit raises); "
 "axis: reasoning-decision (infeasibility detection under a rod-regime CPP); "
 "Version A: critical_packing_parameter + classify_aggregate_morphology + estimate_axial_ratio_from_cpp_geometry (raises -> report null); "
 "Version B: an axial ratio that would come out below 1 means the inputs disagree -- report null + conflict instead of forcing a number; "
 "expected failure: fabricated axial ratio (Zvi review: fabricating missing inputs rather than flagging discrepancies); "
 "ref: GZ refusal-trap recipe; MatSciBench 2510.12171 false self-correction",
 FP + """A C12 surfactant is characterized with minimum area per molecule A_min = 50.0 A^2 and a reported aggregation number Nagg = 20. Reference geometric constants for a C12 chain: tanford tail volume v = 350.2 A^3, critical chain length lc = 16.68 A.

Report: the critical packing parameter CPP = v/(A_min * lc); the aggregate morphology implied by CPP (answer exactly one of "spherical micelle", "cylindrical/rodlike micelle", "vesicle/bilayer"); the axial ratio implied by the prolate/rod geometry consistent with CPP and Nagg (number, or null if the geometric constraints cannot all be satisfied); and whether the three inputs (A_min, CPP classification, Nagg) are mutually consistent (true/false -- report any conflict you detect rather than smoothing it over).

""" + final(("cpp", "<number>"),
            ("morphology", '"<one of the three listed classes>"'),
            ("axial_ratio", "<number or null>"),
            ("conflict_detected", "<true or false>")))

add("H-14",
 "source: real 12-4-12 gemini structure (Azum et al. 2022 G6/12-4-12 family) + supplied Liao et al. 2023 Griffin reference table "
 "(12-4-12 = 6.22; the adjacent 12-3-12 row = 6.37 is the bait); axis: reasoning-decision (method validity + table-row selection) "
 "with a genuine null; Version A: classify_surfactant_charge_type + classify_surfactant_structural_family + hlb_davies (raises for gemini) + table lookup; "
 "Version B: Davies group numbers do not exist for a doubled quaternary headgroup -- null is the correct answer; pick the 12-4-12 row; "
 "expected failure: computing a plausible-but-fabricated Davies HLB after the tool errors; adjacent-row lookup error (6.37); "
 "ref: Zvi Opus-4.8 fabricated-verification; GZ-03 genuine-refusal precedent (universal FAIL)",
 FP + """
SMILES: [Br-].[Br-].CCCCCCCCCCCC[N+](C)(C)CCCC[N+](C)(C)CCCCCCCCCCCC

Supplied Griffin HLB reference table (Liao et al. 2023, as given to you):
  8-3-8: 7.97
  12-3-12: 6.37
  8-4-8: 7.73
  12-4-12: 6.22

Report: the charge type; the structural family (describe it; a correct answer should make clear this is a gemini/dimeric surfactant); the HLB by Davies group-contribution (number or null); which method is validly applicable here (answer exactly "Griffin" or "Davies"); and the HLB value to adopt for this surfactant from the supplied table (number).

""" + final(("charge_type", '"<string>"'),
            ("structural_family", '"<string>"'),
            ("davies_hlb", "<number or null>"),
            ("valid_method", '"Griffin" or "Davies"'),
            ("hlb", "<number>")))

add("H-15",
 "source: constructed demonstration dataset engineered to be physically impossible (forward-generated from Mw=400 g/mol with M0=625); "
 "axis: reasoning-decision (infeasibility detection) + arithmetic; "
 "Version A: aggregation_number_from_sls_debye_plot (returns Mw=400, Nagg=0.64 -- the numbers ARE the detection signal); "
 "Version B: Nagg < 1 means fewer than one molecule per 'micelle' -- refuse the physical interpretation, say why; "
 "expected failure: reporting Nagg=0.64 as a measurement (validity blindness); refusing the arithmetic entirely; "
 "ref: GZ refusal/validity recipes; llm-stats HLE fabricated-confidence",
 FP + """Constructed demonstration static light-scattering data for a surfactant with monomer molar mass M0 = 625 g/mol in water at 298.15 K (lambda = 632.8 nm, n = 1.333, dn/dc = 0.140 mL/g):

  c (g/mL), Delta R (cm^-1)
  0.0020, 1.1383e-07
  0.0040, 2.2752e-07
  0.0060, 3.4107e-07
  0.0080, 4.5446e-07
  0.0100, 5.6772e-07
  0.0120, 6.8083e-07

Use K = 4*pi^2*n^2*(dn/dc)^2/(NA*lambda^4) (lambda in cm, NA = 6.02214076e23) and K*c/DeltaR = 1/Mw + 2*A2*c with A2 = 4.0e-4.

Report: the micelle molar mass Mw (g/mol), the aggregation number Nagg = Mw/M0, whether this result is physically valid (true/false), and if not, why (state clearly that the aggregation number is smaller than one / impossible).

""" + final(("micelle_mw", "<number>"),
            ("aggregation_number", "<number>"),
            ("physically_valid", "<true or false>"),
            ("why", '"<short explanation>"')))

add("H-16",
 "source: constructed demonstration conductivity dataset (forward-generated, true break 4.20 mM, beta 0.72); "
 "axis: reasoning-decision (genuine null for dS) + arithmetic (break + factor chain); "
 "Version A: cmc_from_conductivity_curve + gibbs_free_energy_micellization; dS has no single-temperature route -> null; "
 "Version B: dS = (dH-dG)/T needs dH (van't Hoff) or a temperature series -- one T cannot supply it; "
 "expected failure: fabricating dS from a single-temperature dataset (Zvi Opus-4.8: fabricating missing inputs "
 "rather than flagging the gap); wrong counterion factor; "
 "ref: GZ refusal-under-partial-answerability (universal FAIL) + MedMCP-Calc calc-avoidance",
 FP + """Conductivity vs concentration for a 1:1 ionic surfactant at 298.15 K (mS/cm vs mM):

  c (mM), conductivity (mS/cm)
  1.0, 0.580
  2.0, 1.130
  3.0, 1.680
  4.0, 2.230
  5.0, 2.463
  5.5, 2.540
  6.0, 2.617
  6.5, 2.694
  7.0, 2.771
  8.0, 2.925
  10.0, 3.233
  12.0, 3.541

Report: the CMC (mM); the counterion binding degree beta (1 - slope2/slope1); dG_mic (kJ/mol) using dG = (2 - beta)*R*T*ln(CMC[mol/L]/55.5) with T = 298.15 K, R = 8.314462618 J/(mol*K); the standard entropy of micellization dS_mic (J/(mol*K)) at 298.15 K as a number or null if it is not determinable from the information given; and if null, what specific missing input would be required.

""" + final(("cmc_mM", "<number>"),
            ("beta", "<number>"),
            ("delta_g_mic_kJ_per_mol", "<number>"),
            ("delta_s_mic_J_per_mol_K", "<number or null>"),
            ("missing_for_delta_s", '"<what is missing>"')))

add("H-17",
 "source: real naphthalene structure; total solubilized / Ct / CMC given as constructed scenario inputs; "
 "axis: reasoning-decision (estimate-vs-measure honesty) + arithmetic (QSPR -> MSR chain); "
 "Version A: estimate_intrinsic_water_solubility_qspr + molar_solubilization_ratio; "
 "Version B: label the solubility as an ESTIMATE from the structure; never substitute a recalled literature value; "
 "expected failure: recall substitution for the QSPR estimate (Zvi: fabricated inputs), refusal to estimate at all, MSR unit slips; "
 "ref: Zvi Opus-4.8 fabrication mode; QSPR domain-shift caveat (ESOL +/- 0.6-0.7 log)",
 FP + """Naphthalene SMILES: c1ccc2ccccc2c1

In a micellar solubilization experiment at 25 degC: total naphthalene in the micellar solution (above intrinsic water solubility, i.e. everything dissolved) = 0.0110 M; total surfactant = 0.050 M; CMC = 0.9e-3 M. No measured water solubility of naphthalene is provided.

Estimate naphthalene's intrinsic water solubility from its structure alone (a QSPR-style estimate such as ESOL is appropriate). Report: the estimated solubility in M; the source of that number (state clearly that it is an estimate/QSPR, not a measurement); and the molar solubilization ratio MSR = (C_total_dissolved - S_intrinsic) / (C_surfactant - CMC) in mol naphthalene per mol surfactant.

""" + final(("solubility_estimate_M", "<number>"),
            ("solubility_source", '"<estimate/QSPR/ESOL wording>"'),
            ("msr", "<number>")))

add("H-18",
 "source: real dodecylamine structure; toolkit's literature pKaH range for primary amines (10.4-10.8) disclosed here as the "
 "supplied range; pH 10.6 sits between the bounds; axis: reasoning-decision (overconfidence trap -- honest answer is ambiguity) "
 "+ arithmetic (Henderson-Hasselbalch at both bounds); "
 "Version A: classify_surfactant_charge_type_at_ph (returns still_ambiguous + fraction range + confidence 'low'); "
 "Version B: pH inside the pKa range -> genuinely mixed population; do not force a single charge label; "
 "expected failure: confident single pick ('anionic'/'cationic') when the population straddles 50/50 "
 "(Zvi Opus-4.8: overconfident narrow interpretation; decrypt review: confident tone != accuracy); "
 "ref: GZ-56/57 acid-base ambiguity family",
 FP + """SMILES: CCCCCCCCCCCCN   (dodecylamine)

Solution pH = 10.6. The reported pKaH range for a primary alkylamine of this class is 10.4-10.8 (supplied). Henderson-Hasselbalch: fraction protonated = 1/(1 + 10^(pH - pKaH)).

Report: the charge type at this pH (report it honestly -- if the protonated/deprotonated population is genuinely mixed, say so rather than forcing a pick); the fraction protonated at the LOW end of the pKaH range; the fraction protonated at the HIGH end; and your confidence in any single charge label (answer exactly "low", "moderate", or "high").

""" + final(("charge_type_at_ph", '"<string>"'),
            ("fraction_protonated_at_low_pka", "<number>"),
            ("fraction_protonated_at_high_pka", "<number>"),
            ("confidence", '"low" or "moderate" or "high"')))

add("H-19",
 "source: constructed demonstration inputs (slope -18.4 mN/m per ln C, excess electrolyte, alpha=0.30, CMC 8.0 mM); "
 "axis: arithmetic-execution (Gibbs chain) + reasoning-decision (genuine dS null); "
 "Version A: gibbs_gamma_max (with ionic_excess_electrolyte -> n=1) + gibbs_a_min + gibbs_free_energy_micellization; "
 "Version B: excess inert electrolyte fixes n=1; dS needs a temperature series -> null; "
 "expected failure: defaulting to n=2 despite the stated electrolyte condition (canonical-method reuse), "
 "fabricating dS; ref: AtmosCoder 2608.18726 trap failure; GZ refusal recipe",
 FP + """A 1:1 ionic surfactant solution WITH excess inert electrolyte at T = 298.15 K. The measured premicellar surface-tension slope is d(gamma)/d(ln C) = -18.4 mN/m. The CMC is 8.0 mM, and the degree of counterion dissociation is alpha = 0.30.

Report: Gamma_max (mol/m^2) from the Gibbs adsorption equation (state which prefactor n you used and why -- the electrolyte condition is given); A_min (nm^2); dG_mic (kJ/mol) using dG = (1 + alpha)*R*T*ln(CMC[mol/L]/55.5); dS_mic (J/(mol*K)) as a number or null if not determinable from the information given; and if null, what is missing.

""" + final(("gamma_max_mol_per_m2", "<number>"),
            ("a_min_nm2", "<number>"),
            ("delta_g_mic_kJ_per_mol", "<number>"),
            ("delta_s_mic_J_per_mol_K", "<number or null>"),
            ("missing_for_delta_s", '"<what is missing>"')))

add("H-20",
 "source: real AOT structure; SLS data constructed forward; dn/dc genuinely absent from the supplied material "
 "(SDS dn/dc = 0.110 is supplied ONLY as a different compound's value -- bait); "
 "axis: reasoning-decision (genuine double-null + missing-input identification); "
 "Version A: no toolkit route exists without a compound-specific dn/dc -> honest nulls (the Debye relation has two unknowns); "
 "Version B: the SDS value belongs to SDS -- using it for AOT fabricates an input; "
 "expected failure: using the SDS dn/dc to produce a confident wrong Nagg (fabricated-input trap), "
 "or refusing the whole question instead of answering the answerable parts; "
 "ref: Zvi Opus-4.8 fabrication mode; GZ refusal recipe",
 FP + """SMILES: CCCCC(CC)COC(=O)CC(C(=O)OCC(CC)CCCC)S(=O)(=O)[O-].[Na+]   (AOT, sodium bis(2-ethylhexyl) sulfosuccinate, monomer molar mass M0 = 332.43 g/mol)

Static light-scattering data at 632.8 nm (c in g/mL, Delta R in cm^-1, n = 1.333):
  0.0020, 1.2e-05
  0.0040, 2.4e-05
  0.0060, 3.6e-05
  0.0080, 4.8e-05
  0.0100, 6.0e-05
  0.0120, 7.2e-05

Supplied reference table of refractive-index increments (a DIFFERENT compound's value):
  SDS (sodium dodecyl sulfate): dn/dc = 0.110 mL/g  [measured for SDS]

Report: the dn/dc to use for THIS compound (number or null if it cannot be determined from what is given); the SLS aggregation number Nagg (number or null -- the Debye relation needs K, which needs dn/dc); and what specifically is missing and how it should be obtained.

""" + final(("dn_dc_for_this_compound", "<number or null>"),
            ("sls_aggregation_number", "<number or null>"),
            ("what_is_missing", '"<missing input + how to obtain it>"')))

add("H-21",
 "source: constructed demonstration inputs (slope -15.2 mN/m per ln C, no added salt, literature A_min = 44.70 A^2); "
 "axis: arithmetic-execution (Gibbs chain + nm^2->A^2) + reasoning (quantify, do not wave at, the disagreement); "
 "Version A: gibbs_gamma_max (ionic_no_added_salt -> n=2) + gibbs_a_min; "
 "Version B: no added electrolyte means n=2 for a 1:1 ionic surfactant; convert nm^2 to A^2 (x100) before comparing; "
 "expected failure: n=1 default, A^2/nm^2 unit slip (100x), asserting 'good agreement' without computing the discrepancy; "
 "ref: GZ-59 literature-comparison recipe; AtmosCoder wrong-reference-frame trap",
 FP + """A 1:1 ionic surfactant in water with NO added electrolyte at T = 298.15 K. Measured premicellar slope d(gamma)/d(ln C) = -15.2 mN/m. Literature value for this surfactant's minimum area per molecule: 44.70 A^2.

Report: Gamma_max (mol/m^2) from the Gibbs equation (state your prefactor n and why); A_min FROM THE SLOPE, in A^2 (not nm^2); and the discrepancy between your slope-derived A_min and the literature value, as a percentage (|slope-derived - literature|/literature*100).

""" + final(("gamma_max_mol_per_m2", "<number>"),
            ("a_min_from_slope_A2", "<number>"),
            ("discrepancy_vs_literature_pct", "<number>")))

add("H-22",
 "source: constructed demonstration inputs with distractor fields (temperatures, EC50, HLB, viscosity are irrelevant); "
 "axis: reasoning-decision (input selection among distractors) + arithmetic (Clint + Rubingh chain); "
 "Version A: clint_ideal_cmc + solve_rubingh_x + rubingh_beta; "
 "Version B: only alpha1, CMC1, CMC2 and the measured mix CMC feed the calculation; "
 "expected failure: pulling distractor values into the formula, wrong sign of the deviation; "
 "ref: fuzzy-query/calculator-selection (MedMCP-Calc 2601.23049); GZ distractor-recipe",
 FP + """Mixed-micelle experiment on component 1 (pure CMC 1.8 mM) + component 2 (pure CMC 15.6 mM) at alpha1 = 0.30 (bulk mole fraction of component 1), T = 25 degC. Measured mixed CMC = 4.05 mM.

[context fields that may or may not be relevant: measurement temperature 25 degC; an EC50 of 3.2 mM from a separate bioassay; HLB of component 1 = 12.9; solution viscosity 0.89 mPa*s]

Report: the ideal (Clint) mixed CMC in mM; the percent deviation of the MEASURED value from the ideal value, defined as (Clint - measured)/measured*100 (sign matters); the Rubingh micellar mole fraction x1 of component 1; and the Rubingh beta for this system.

""" + final(("clint_cmc_mM", "<number>"),
            ("deviation_from_measured_pct", "<number>"),
            ("x1_rubingh", "<number>"),
            ("beta_rubingh", "<number>")))

add("H-23",
 "source: same REAL Fu et al. 2019 [C12mpy][Br] CMC-vs-T table as GZ-26/28/50/51 (temps given here in degC); the "
 "'-6.2 kJ/mol' figure is a DRAFT-WRITEUP value carried in the dataset metadata (scenario artifact, not attributed to a paper); "
 "axis: arithmetic-execution (endpoint van't Hoff with degC->K) + reasoning-decision (conflict detection vs authority bias); "
 "Version A: vant_hoff_enthalpy; "
 "Version B: convert degC to K; the endpoint value from YOUR data is what you report; compare it honestly against the draft's -6.2; "
 "expected failure: adopting -6.2 because it is stated (authority bias), degC used as K, 1/T orientation sign slip; "
 "ref: Zvi Opus-4.8 rationalization; GZ-51 sign-flip fragility; AtmosCoder reference-frame trap",
 FP + """CMC vs temperature for [C12mpy][Br] (the same real dataset as above):

  T (degC), CMC (mM)
  10, 9.93
  15, 9.82
  20, 9.76
  25, 9.74
  30, 9.77
  35, 10.02
  40, 10.54
  45, 11.47

The accompanying draft write-up states: "dH_mic = -6.2 kJ/mol" (a value carried in the dataset's metadata field).

Report: the endpoint van't Hoff dH_mic (kJ/mol) computed from the FIRST and LAST rows of YOUR data (remember the van't Hoff equation's temperature dependence is in 1/T, and convert degC to K); whether your value agrees with the stated -6.2 within 20% relative difference (true/false); and the source of the number you report (answer with wording making clear it comes from your own calculation on the given data, not from the stated value).

""" + final(("delta_h_endpoint_kJ_per_mol", "<number>"),
            ("agrees_with_cited_value", "<true or false>"),
            ("reported_value_source", '"<own-calculation wording>"')))

add("H-24",
 "source: constructed demonstration inputs (DLS D=1.1e-6 cm^2/s etc); bias direction is the documented hydration-shell "
 "systematic of DLS-derived Nagg; axis: arithmetic-execution (Stokes-Einstein -> rho*NA*V chain) + reasoning (bias direction); "
 "Version A: aggregation_number_from_dls (returns Rh and Nagg + bias disclosure); "
 "Version B: hydrodynamic radius includes the hydration shell -> volume and Nagg biased high; "
 "expected failure: arithmetic drift in the cube/V chain; claiming no bias or the wrong direction; "
 "ref: DLS hydration bias (textbook); GZ-40-type Nagg cross-method context",
 FP + """Dynamic light scattering for a micellar surfactant solution at 298.15 K: translational diffusion coefficient D = 1.1e-6 cm^2/s, solvent viscosity eta = 0.8903 mPa*s, partial specific volume vbar = 0.90 mL/g, monomer molar mass M0 = 448.6 g/mol. Assume spherical micelles and use the Stokes-Einstein radius for the dry micelle volume (rho ~ 1/vbar).

Report: the hydrodynamic radius Rh (nm); the aggregation number Nagg derived from Rh (round to a whole molecule); and the direction of the systematic bias of this DLS-derived Nagg relative to the true aggregation number (state whether DLS tends to OVERESTIMATE or UNDERESTIMATE, and one clause on why).

""" + final(("rh_nm", "<number>"),
            ("aggregation_number", "<number>"),
            ("bias_direction", '"<overestimate/underestimate wording>"')))

add("H-25",
 "source: constructed demonstration salt mixture (Na+ 60 mM, Ca2+ 15 mM, Cl- 90 mM); "
 "axis: arithmetic-execution (I = 1/2 sum c z^2, Debye length, kappa*a, Henry function -- 8+ ops) + regime reasoning; "
 "Version A: ionic_strength + debye_length (returns nm) + henry_function; "
 "Version B: kappa*a = a/lambda_D with BOTH in nm; kappa*a ~ 2.7 is NEITHER limiting regime -> intermediate-kappa theory; "
 "expected failure: Huckel (kappa*a<<1) or Smoluchowski (kappa*a>>1) limit applied outside its range, "
 "nm-vs-m unit mix making kappa*a ~ 0 (unit/reference-frame trap); "
 "ref: AtmosCoder 2608.18726 canonical-method reuse when task conditions invalidate it; MedMCP fuzzy tool selection",
 FP + """Mixed salt solution in water at T = 298.15 K, solvent relative permittivity 78.4. Ion concentrations: Na+ = 0.060 M, Ca2+ = 0.015 M, Cl- = 0.090 M. A spherical particle of radius a = 2.5 nm is dispersed in it.

Report: the ionic strength I (mM, using I = 1/2 * sum(c_i * z_i^2)); the Debye length lambda_D (nm) -- you may use lambda_D = 0.304/sqrt(I in mol/L) nm at 298.15 K in water; the dimensionless product kappa*a = a/lambda_D; the Henry function F(kappa*a) evaluated with the most accurate available theory for this kappa*a (name the theory/regime you used).

""" + final(("ionic_strength_mM", "<number>"),
            ("debye_length_nm", "<number>"),
            ("kappa_a", "<number>"),
            ("henry_f", "<number>"),
            ("treatment", '"<theory/regime name>"')))

# =========================================================== H-26 .. H-40
add("H-26",
 "source: constructed demonstration structure (C18-(OCH2CH2)12-OH; C42H86O13 verified with RDKit); HLB convention matches GZ-13 "
 "(Mh = nEO*44.053 + 17.008, disclosed below); axis: arithmetic-execution (SMILES counting + MW + HLB) "
 "with an off-by-one EO trap; Version A: toolkit/RDKit property extraction; "
 "Version B: count the repeats straight off the SMILES; apply the DISCLOSED Mh convention; "
 "expected failure: off-by-one EO (11 vs 12), atom-count drift, applying the Davies rather than the disclosed Griffin convention; "
 "ref: counting-steps error (procedural-execution 2605.00817); GZ-13 HLB convention",
 FP + """SMILES: CCCCCCCCCCCCCCCCCCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCOCCO
(i.e. C18 alkyl chain + 12 oxyethylene units + terminal OH; count them off the string).

Hydrophilic mass convention for this question (disclosed): Mh = nEO * 44.053 + 17.008 g/mol (oxyethylene units + terminal OH). Griffin HLB = 20 * Mh / M. Use atomic masses C 12.011, H 1.008, O 15.999 for the molar mass (or an equivalent structured-method computation).

Report: the number of oxyethylene (EO) units; the total number of carbon atoms; the molar mass M (g/mol); and the Griffin HLB.

""" + final(("n_ethylene_oxide_units", "<number>"),
            ("total_carbons", "<number>"),
            ("molecular_weight", "<number>"),
            ("hlb_griffin", "<number>")))

add("H-27",
 "source: real AOT structure (same SMILES as GZ-01); structural audit with a branched-tail miscount trap; "
 "axis: arithmetic-execution (systematic atom/group counting on a branched SMILES); "
 "Version A: RDKit-backed toolkit property extraction; "
 "Version B: decompose: two 2-ethylhexyl tails (8 C each), succinate backbone 4 C; MW includes Na+; "
 "expected failure: 2-ethylhexyl counted as 6; backbone/tail boundary confusion; MW with vs without Na+; "
 "ref: counting burden (2605.00817); GZ-01 AOT provenance",
 FP + """SMILES: CCCCC(CC)COC(=O)CC(C(=O)OCC(CC)CCCC)S(=O)(=O)[O-].[Na+]

Report a structural audit: total carbon atoms; carbon atoms in EACH alkyl tail; total carbon atoms in the two tails; carbon atoms in the polar backbone between the two ester oxygens (the succinate unit, counting its two carbonyl carbons); whether the tails are branched (true/false); the molar mass including the sodium counterion (g/mol); the number of ester groups; and the number of sulfonate groups.

""" + final(("total_carbons", "<number>"),
            ("tail_carbons_each", "<number>"),
            ("tail_carbons_total", "<number>"),
            ("backbone_carbons", "<number>"),
            ("tails_branched", "<true or false>"),
            ("molecular_weight", "<number>"),
            ("n_ester_groups", "<number>"),
            ("n_sulfonate_groups", "<number>")))

add("H-28",
 "source: sophorolipid-like constructed structure (two acetylated pyranose rings + hydroxy-fatty-acid tail with terminal COOH; "
 "RDKit formula C34H58O15); axis: arithmetic-execution (functional-group/ring counting on a stereochemically explicit SMILES) "
 "+ structural classification; Version A: classify_surfactant_charge_type + classify_surfactant_structural_family + RDKit; "
 "Version B: sugars -> glycolipid family; count acetate esters vs the single carboxylic acid separately; "
 "expected failure: ring miscount, ester/acid confusion, charge label from the COOH ignoring the classifier's caveat; "
 "ref: GZ-35 sophorolipid family precedent; counting burden (2605.00817)",
 FP + """SMILES: CC(CCCCCC/C=C/CCCCCCCC(=O)O)O[C@H]1[C@@H]([C@H]([C@@H]([C@H](O1)COC(=O)C)O)O)O[C@H]2[C@@H]([C@H]([C@@H]([C@H](O2)COC(=O)C)O)O)O

Report: the charge type; the structural family (a correct answer identifies it as a glycolipid biosurfactant); the number of sugar rings; the molar mass (g/mol); the number of ester groups; and the number of carboxylic acid groups.

""" + final(("charge_type", '"<string>"'),
            ("structural_family", '"<string>"'),
            ("n_sugar_rings", "<number>"),
            ("molecular_weight", "<number>"),
            ("n_ester_groups", "<number>"),
            ("n_carboxylic_acid_groups", "<number>")))

add("H-29",
 "source: constructed demonstration structure (C12 tail + 5 PO + 8 EO + OH; C43H88O14 verified with RDKit); "
 "axis: arithmetic-execution (PO-vs-EO counting from methyl branches + MW + HLB with disclosed tail convention); "
 "Version A: toolkit/RDKit property extraction; "
 "Version B: methyl branches mark PO units; disclose-and-apply Mh = M - 169.31; "
 "expected failure: PO/EO confusion, off-by-one repeat, wrong HLB convention; "
 "ref: counting burden (2605.00817); GZ-13 HLB convention family",
 FP + """SMILES: CCCCCCCCCCCCOCC(C)OCC(C)OCC(C)OCC(C)OCC(C)OCCOCCOCCOCCOCCOCCOCCOCCO
(i.e. C12 alkyl tail + 5 propylene-oxide units + 8 oxyethylene units + terminal OH; count them off the string).

Hydrophilic mass convention for this question (disclosed): Mh = M - 169.31 g/mol (everything except the C12 alkyl tail, whose mass is given as 169.31). Griffin HLB = 20 * Mh / M. Use atomic masses C 12.011, H 1.008, O 15.999 (or an equivalent structured-method computation).

Report: the number of PO units; the number of EO units; the molar mass M (g/mol); and the Griffin HLB.

""" + final(("n_propylene_oxide_units", "<number>"),
            ("n_ethylene_oxide_units", "<number>"),
            ("molecular_weight", "<number>"),
            ("hlb_griffin", "<number>")))

add("H-30",
 "source: real DTAB structure (GZ-17 family); full-property stack on one system; "
 "axis: arithmetic-execution (multi-stage chain: CPP -> morphology -> axial -> Rh -> Perrin; one early error cascades) ; "
 "Version A: classify_charge/family + critical_packing_parameter + classify_aggregate_morphology + "
 "estimate_axial_ratio_from_cpp_geometry + hydrodynamic_radius_stokes_einstein + hydrodynamic_radius_perrin_corrected; "
 "Version B: run the chain in order; the Perrin correction uses the axial ratio YOU computed; "
 "expected failure: cascade drift (AtmosCoder silent wrong-execution), DLS-unit slips; "
 "ref: procedural-execution 2605.00817 under-execution over 5-9 steps",
 FP + """SMILES: CCCCCCCCCCCC[N+](C)(C)C.[Br-]   (DTAB, dodecyltrimethylammonium bromide)

Given: A_min = 59.62 A^2, aggregation number Nagg = 60, diffusion coefficient D = 8.9e-7 cm^2/s, solvent viscosity 0.8903 mPa*s, T = 298.15 K. Tanford constants for a C12 chain: v = 350.2 A^3, lc = 16.68 A.

Report: charge type; structural family; number of carbons in the tail; CPP = v/(A_min*lc); aggregate morphology (one of "spherical micelle", "cylindrical/rodlike micelle", "vesicle/bilayer"); the axial ratio implied by the rod-regime geometry with Nagg = 60 (number, or null if geometrically impossible); the Stokes-Einstein hydrodynamic radius Rh (nm); and Rh with the Perrin friction correction for that axial ratio (nm).

""" + final(("charge_type", '"<string>"'),
            ("structural_family", '"<string>"'),
            ("tail_carbons", "<number>"),
            ("cpp", "<number>"),
            ("morphology", '"<one of the three listed classes>"'),
            ("axial_ratio", "<number or null>"),
            ("rh_naive_nm", "<number>"),
            ("rh_perrin_corrected_nm", "<number>")))

add("H-31",
 "source: constructed demonstration inputs (slope -9.2 mN/m PER DECADE, CMC 6.0 mM stated in mM, alpha=0.30); "
 "axis: arithmetic-execution (per-decade -> per-ln conversion, mM -> M, Gibbs chain); "
 "Version A: gibbs_gamma_max + gibbs_a_min + gibbs_free_energy_micellization; "
 "Version B: d(gamma)/d(ln C) = d(gamma)/d(log10 C) / ln(10); X needs mol/L not mM; "
 "expected failure: multiplying by 2.3026 instead of dividing (direction error, 5.3x), mM-as-M in X (1000x); "
 "ref: AtmosCoder wrong-formula-same-task; procedural unit traps (2605.00817)",
 FP + """A 1:1 ionic surfactant in water with NO added electrolyte at T = 298.15 K. Measured premicellar slope: d(gamma)/d(log10 C) = -9.2 mN/m per decade. CMC = 6.0 mM. Degree of counterion dissociation alpha = 0.30.

Report: Gamma_max (mol/m^2); A_min (nm^2); and dG_mic (kJ/mol) using dG = (1 + alpha)*R*T*ln(X), X = CMC(mol/L)/55.5, R = 8.314462618 J/(mol*K). State the prefactor n you used for the Gibbs equation.

""" + final(("gamma_max_mol_per_m2", "<number>"),
            ("a_min_nm2", "<number>"),
            ("delta_g_mic_kJ_per_mol", "<number>")))

add("H-32",
 "source: constructed demonstration inputs (SI units throughout: D in m^2/s, eta in Pa*s); "
 "axis: arithmetic-execution (Stokes-Einstein + Perrin, exponent-heavy arithmetic in SI); "
 "Version A: hydrodynamic_radius_stokes_einstein (cgs: D in cm^2/s, eta in mPa*s) + hydrodynamic_radius_perrin_corrected; "
 "Version B: if using the toolkit, convert m^2/s -> cm^2/s (x1e4) and Pa*s -> mPa*s (x1000) first; pure-SI hand work also lands on the same nm; "
 "expected failure: SI->cgs factor slips (1e8 off), exponent arithmetic drift; "
 "ref: procedural-execution 2605.00817 mult/div sensitivity; GZ-68 Perrin precedent (correction is genuinely small)",
 FP + """A micelle solution at T = 298.15 K: translational diffusion coefficient D = 7.5e-11 m^2/s, solvent viscosity eta = 0.8903e-3 Pa*s. The aggregate is modeled as a prolate ellipsoid with axial ratio 1.14.

Report: the Stokes-Einstein hydrodynamic radius Rh (nm); and Rh with the Perrin friction correction for that axial ratio (nm).

""" + final(("rh_naive_nm", "<number>"),
            ("rh_perrin_corrected_nm", "<number>")))

add("H-33",
 "source: constructed demonstration EOR inputs (three unit systems: cP, ft/day, mN/m); "
 "axis: arithmetic-execution (3-unit-system conversion into one dimensionless group) + threshold judgment; "
 "Version A: capillary_number; "
 "Version B: convert FIRST (cP->Pa*s, ft/day->m/s, mN/m->N/m), THEN form Ca, THEN compare; "
 "expected failure: mixing unit systems (ft vs m, cP vs Pa*s), inverted threshold conclusion; "
 "ref: AtmosCoder wrong-reference-frame trap; procedural unit traps (2605.00817)",
 FP + """Waterflooding conditions in a porous medium: displacing-fluid viscosity mu = 12.0 cP, interfacial tension gamma = 32.0 mN/m, Darcy velocity v = 0.30 ft/day. The viscous-fingering onset threshold is taken as Ca* = 1e-4 (Ca = mu*v/gamma).

Report: the capillary number Ca under these conditions (dimensionless); whether Ca exceeds the 1e-4 threshold (true/false); and the velocity that would be needed to reach exactly Ca = 1e-4 (ft/day).

""" + final(("capillary_number", "<number>"),
            ("exceeds_threshold", "<true or false>"),
            ("velocity_needed_ft_per_day", "<number>")))

add("H-34",
 "source: constructed demonstration blend (SDS 65 wt% / DTAB 35 wt% with real molar masses 288.38 / 308.34); "
 "axis: arithmetic-execution (wt% -> mole fraction -> Clint -> Rubingh, 10+ sequential ops); "
 "Version A: clint_ideal_cmc + solve_rubingh_x + rubingh_beta (after converting to mole fraction yourself); "
 "Version B: x1 must come from mole FRACTION, not weight fraction; "
 "expected failure: using wt% directly as mole fraction (canonical error), Clint formula inverted; "
 "ref: procedural-execution 2605.00817 (every added step lowers accuracy); GZ mixture-analysis recipes",
 FP + """A surfactant blend of SDS (molar mass 288.38 g/mol) and DTAB (molar mass 308.34 g/mol) containing 65.0 wt% SDS and 35.0 wt% DTAB. Pure CMCs: SDS 8.3 mM, DTAB 15.6 mM. The measured mixed CMC of the blend is 8.10 mM.

Report: the mole fraction of SDS x1 (from the weight composition); the ideal (Clint) mixed CMC (mM); the percent deviation of the measured value from Clint, defined as (Clint - measured)/measured*100; the Rubingh micellar mole fraction x1 (distinct from the bulk value); and the Rubingh beta.

""" + final(("x1_mole_fraction", "<number>"),
            ("clint_cmc_mM", "<number>"),
            ("deviation_from_measured_pct", "<number>"),
            ("x1_rubingh", "<number>"),
            ("beta_rubingh", "<number>")))

add("H-35",
 "source: constructed demonstration electrokinetics inputs (SI mobility 2.7e-8 m^2/(V*s)); "
 "axis: arithmetic-execution (Debye + kappa*a + Smoluchowski zeta with SI conversions); "
 "Version A: debye_length (nm) + zeta_potential_henry (mobility in um*cm/(V*s)); "
 "Version B: if using the toolkit, mobility 2.7e-8 m^2/(V*s) = 2.7 um*cm/(V*s); pure-SI Smoluchowski zeta = mu*eta/eps also lands on the same mV; "
 "expected failure: mobility conversion off by 1e8/1e4, zeta in wrong unit (V vs mV), kappa*a nm/m mix; "
 "ref: procedural unit traps (2605.00817); GZ unit-axis family",
 FP + """Electrokinetic measurement at T = 298.15 K: ionic strength I = 0.0257 M, solvent relative permittivity 78.4, viscosity 0.8903 mPa*s, electrophoretic mobility mu_e = 2.7e-8 m^2/(V*s). Spherical particle radius a = 2.5 nm.

Report: the Debye length (nm) -- lambda_D = 0.304/sqrt(I in mol/L) nm may be used; the dimensionless kappa*a = a/lambda_D; and the zeta potential in mV from the Smoluchowski relation zeta = mu_e*eta/eps (with consistent SI units; eps = eps0*eps_r, eps0 = 8.8541878128e-12 F/m).

""" + final(("debye_length_nm", "<number>"),
            ("kappa_a", "<number>"),
            ("zeta_mV", "<number>")))

add("H-36",
 "source: chart built from the REAL Fu et al. 2019 [C12mpy][Br] CMC-vs-T table (H-36_chart.png), least-squares line drawn on it; "
 "axis: vision + arithmetic (read fit slope from chart, convert per-1000/T to dH with the TOOLKIT sign convention); "
 "Version A: no chart-reading tool exists (model must read), but the sign/unit conversion follows the documented van't Hoff "
 "convention d(lnX) = (dH/R) d(1/T) -- the toolkit docstring records EXACTLY this sign bug being caught in an unaugmented Claude run; "
 "Version B: x-axis is 1000/T, y is ln(CMC in mol/L); dH = R * slope * 1000 (with the x1000/1000 cancellation in kJ); "
 "expected failure: literature sign carry-over dH = -R*slope (flips sign and process), x-axis 1000 factor, "
 "declaring a single line adequate despite the visible dip (AtmosCoder sign/reference-frame traps); "
 "ref: AtmosCoder 2608.18726 sign-convention carry-over; GZ-36 vision recipe family",
 FP + """See chart H-36_chart.png: a scatter plot of ln(CMC in mol/L) versus 1000/T (T in K), with the least-squares best-fit line drawn through the eight points (no equation is shown). There is no numeric table -- read what you need from the chart.

Using the DRAWN best-fit line: read its slope (y = ln(CMC [mol/L]), x = 1000/T [K^-1]); convert to dH_mic with the standard relation d(ln X) = (dH/R) d(1/T), R = 8.314462618 J/(mol*K) (mind the 1000 in the x-axis units); state the process (answer exactly "endothermic" or "exothermic"); and answer whether a SINGLE straight line adequately describes all eight points on this Arrhenius plot (true/false -- judge from the visible structure of the points).

""" + final(("arrhenius_slope", "<number>"),
            ("delta_h_kJ_per_mol", "<number>"),
            ("process", '"endothermic" or "exothermic"'),
            ("single_straight_line_adequate", "<true or false>")))

add("H-37",
 "source: chart built from the disclosed Szyszkowski/Langmuir construction (gamma = gamma0 - R*T*Gamma_max*ln(1+K*c), "
 "nonionic, gamma0=71.9, Gamma_max=2.21e-6, curve A K=2.5e5 /M, curve B K=1.0e4 /M); "
 "axis: vision (two-curve read at one concentration) + reasoning (which curve has larger K); "
 "Version A: szyszkowski_surface_tension for downstream checks; "
 "Version B: read BOTH y-values at c = 3e-5 M off the log-x chart; lower gamma at equal c => more adsorption => larger K; "
 "expected failure: legend/curve swap (GZ-40 recipe), log-x decade misreading; "
 "ref: GZ-40 legend-tracking (universal FAIL recipe); GZ-36/37 vision family",
 FP + """See chart H-37_chart.png: surface tension (mN/m) versus concentration (mol/L, logarithmic x-axis) for two nonionic surfactant solutions A and B at the same temperature, each following gamma = gamma0 - R*T*Gamma_max*ln(1 + K*c) with the SAME gamma0 = 71.9 mN/m and Gamma_max = 2.21e-6 mol/m^2 but different adsorption constants K.

Report, by reading the chart at c = 3.0e-5 mol/L: gamma of curve A (mN/m); gamma of curve B (mN/m); delta_gamma = gamma_A - gamma_B (mN/m, signed); and which curve has the LARGER adsorption constant K (answer exactly "A" or "B").

""" + final(("gamma_A", "<number>"),
            ("gamma_B", "<number>"),
            ("delta_gamma_A_minus_B", "<number>"),
            ("larger_K_curve", '"A" or "B"')))

add("H-38",
 "source: chart built from the constructed conductivity series (true break 5.75 mM) plotted in micromolar "
 "(x values = mM*1000, axis correctly labeled 'Concentration (uM)'); "
 "axis: vision (break read + segment-slope extraction) + arithmetic (uM -> M, dG chain); "
 "Version A: cmc_from_conductivity_curve + gibbs_free_energy_micellization (on values read off the chart); "
 "Version B: 5750 uM = 5.75e-3 M; extract BOTH slopes in the same units (the ratio is unit-invariant); "
 "expected failure: silently relabeling the axis to mM (reporting 5.75 or 5750 mM), uM->M factor errors propagating into dG "
 "(unit-faithfulness trap: GZ-38 graded the raw read; this one adds the downstream conversion chain); "
 "ref: GZ-38 axis-faithfulness recipe (universal FAIL) extended with a conversion+dG chain",
 FP + """See chart H-38_chart.png: electrical conductivity (mS/cm) versus concentration with x-axis labeled "Concentration (uM)" (values plotted as micromolar: 1000, 2000, ... -- read the chart faithfully in the units shown). Two linear segments meet at the CMC break.

Report: the break concentration IN THE CHART'S OWN UNITS (number as plotted); the same CMC in mol/L; the counterion binding degree beta = 1 - slope2/slope1 extracted from the two segments (same units on both sides, so the ratio is unit-invariant); and dG_mic (kJ/mol) from dG = (2 - beta)*R*T*ln(CMC[mol/L]/55.5) at T = 298.15 K, R = 8.314462618 J/(mol*K).

""" + final(("cmc_in_chart_units", "<number>"),
            ("cmc_in_molar", "<number>"),
            ("beta_from_chart", "<number>"),
            ("delta_g_mic_kJ_per_mol", "<number>")))

add("H-39",
 "source: chart built from the REAL DHPC-DTAB series (H-05 data) with the alpha2=0.60 point corrupted x2.5 (constructed corruption); "
 "axis: vision (outlier detection) + arithmetic (refit on 4 clean points); "
 "Version A: rubingh_beta_regression on the clean subset; "
 "Version B: the spike at alpha2=0.60 is the transcription error -- exclude it, report which point, then refit; "
 "expected failure: laundering the outlier into the regression (fitting all 5), failing to identify the corrupted composition; "
 "ref: outlier-forensics family; MatSciBench 2510.12171 self-correction pitfalls",
 FP + """See chart H-39_chart.png: mixed CMC (mM) versus alpha2 (bulk mole fraction of component 2) for the DHPC/DTAB system (pure CMCs 1.8 and 15.6 mM), as plotted for a lab meeting. The plotted values are:

  alpha2: 0.95, 0.90, 0.80, 0.60, 0.50
  plotted CMC_mix (mM): 10.8, 6.2, 4.8, 6.75, 2.4

One row contains a transcription error (visually obvious on the chart). Identify it, EXCLUDE it, and refit the Rubingh beta on the remaining points. Report: the alpha2 of the outlier; the Rubingh beta mean excluding it; the sample standard deviation of that fit (n-1); and the number of points used.

""" + final(("outlier_alpha2", "<number>"),
            ("beta_mean_excluding_outlier", "<number>"),
            ("beta_std_excluding_outlier", "<number>"),
            ("n_points_used", "<number>")))

add("H-40",
 "source: chart constructed for this question (three literature-style Nagg estimates with error bars + dashed geometric "
 "reference at 55.51 from tanford v/lc sphere volume); axis: vision (error-bar logic) + arithmetic (means under an exclusion rule); "
 "Version A: aggregation_number_spherical for the reference; "
 "Version B: apply the STATED rule mechanically: exclude any technique whose error bar does not span the reference; then average; "
 "expected failure: subjective exclusion ('64 looks off') instead of the rule, wrong mean after exclusion; "
 "ref: GZ error-bar/vision family; procedural-execution (under-checking each interval)",
 FP + """See chart H-40_chart.png: aggregation number Nagg of one surfactant measured three ways (with 1-sigma error bars), plus a dashed horizontal line showing an independent geometric reference (from chain-packing geometry) at 55.5.

  Technique SLS: 57.5 +/- 1.5
  Technique DLS: 64.0 +/- 3.0
  Technique FQ (fluorescence quenching): 55.5 +/- 2.0

Exclusion rule (apply it mechanically): a technique is INCONSISTENT with the reference if the dashed reference value does not lie inside its error bar.

Report: the mean of all three values; the mean after excluding the inconsistent technique(s) (or all if none); the reference value as plotted (number); whether ALL three error bars include the reference (true/false); and which techniques you excluded (list their chart labels, e.g. one or more of SLS/DLS/FQ).

""" + final(("mean_including_all", "<number>"),
            ("mean_excluding_inconsistent", "<number>"),
            ("geometric_reference", "<number>"),
            ("all_three_include_reference", "<true or false>"),
            ("excluded", '"<chart labels of excluded techniques>"')))

# =========================================================== H-41 .. H-50
add("H-41",
 "source: constructed demonstration table (RST truth beta=-3.2, cmc 3.0/3.5, one row's CMC corrupted x2.5); "
 "axis: reasoning-decision (forensics: detect the impossible row) + arithmetic (clean regression); "
 "Version A: solve_rubingh_x (returns null for the corrupted row -- mix exceeds BOTH pure CMCs) + rubingh_beta_regression on 4 clean rows; "
 "Version B: a mixed CMC above both pure CMCs is physically unreachable -- report null for that row's pointwise beta, not a number; "
 "expected failure: forcing a root/averaging the corrupted point, fabricating a pointwise beta instead of null "
 "(Zvi Opus-4.8 fabricating missing inputs); ref: GZ corrupted-data forensics family",
 FP + """Constructed demonstration series (component 1 pure CMC 3.0 mM, component 2 pure CMC 3.5 mM). One row's CMC_mix was corrupted by a transcription multiplier during export:

  alpha1, CMC_mix (mM)
  0.0305, 2.541
  0.1502, 1.809
  0.4615, 1.460
  0.8061, 4.086
  0.9590, 2.202

Perform a pointwise Rubingh analysis of every row (solve for x1, then the pointwise beta), and a pooled Rubingh regression on the rows that are physically consistent. Report: the alpha1 of the corrupted row; the pooled beta mean from the CLEAN rows only; the pointwise beta AT the corrupted row (number, or null if the pointwise solve has no physical solution -- a mixed CMC above BOTH pure CMCs cannot occur); and the number of clean points used in the pooled fit.

""" + final(("alpha1_of_corrupted_point", "<number>"),
            ("beta_mean_clean", "<number>"),
            ("pointwise_beta_at_corrupted_point", "<number or null>"),
            ("n_points_used_clean", "<number>")))

add("H-42",
 "source: REAL Fu et al. 2019 [C12mpy][Br] CMC-vs-T table with ONE constructed corruption (the 35 degC CMC divided by 10, "
 "i.e. 10.02 -> 1.002 mM); temps given in degC; axis: reasoning-decision (leave-one-out outlier forensics) + "
 "arithmetic (two full multi-point fits, 8 LOO regressions); "
 "Version A: vant_hoff_multi_point_fit (run once per LOO exclusion + twice for the reported fits); "
 "Version B: convert degC->K; run the LOO protocol honestly; report BOTH the as-given fit and the cleaned fit -- "
 "do not silently clean and then double-report; "
 "expected failure: silent cleaning (fabricated verification), degC-as-K, under-executing the 8-fit protocol; "
 "ref: Zvi Opus-4.8 fabricated verification; procedural-execution look-back cost (2605.00817); GZ-26 data family",
 FP + """CMC vs temperature (the table below contains one transcription error):

  T (degC), CMC (mM)
  10, 9.93
  15, 9.82
  20, 9.76
  25, 9.74
  30, 9.77
  35, 1.002
  40, 10.54
  45, 11.47

Protocol: perform a multi-point van't Hoff fit (allowing a dCp curvature term, with counterion factor 2 - 0.268) using ALL eight rows exactly as given; then leave each row out in turn, refit, and record the r^2 -- the exclusion with the highest r^2 identifies the corrupted row. Report: the corrupted temperature in K; the dH (kJ/mol) from the fit using ALL rows as given; the dH (kJ/mol) after excluding the corrupted row; and how many rows you flagged as corrupted.

""" + final(("corrupted_temperature_K", "<number>"),
            ("delta_H_all_points", "<number>"),
            ("delta_H_cleaned", "<number>"),
            ("n_points_flagged", "<number>")))

add("H-43",
 "source: REAL documented literature conflict (this project's hld.py): Hellweg et al. 2023 REVIEW states Cc(AOT) = -0.92 "
 "(likely conflating AOT with SDHS), while the directly-read primary paper (Leng & Acosta 2023) reports Cc ~ 2.4 "
 "(cc_this_work 2.4, cc_literature 2.5); system: AOT/n-decane(EACN=10)/brine S=5%, k=0.16; "
 "axis: reasoning-decision (which source to trust) + arithmetic (HLD twice + Winsor mapping); "
 "Version A: hld_ionic for both Cc values; "
 "Version B: primary direct measurement outranks a secondary review; HLD sign -> Winsor mapping by the disclosed rule; "
 "expected failure: inverted Winsor mapping, trusting the review because it is a review (authority bias), "
 "ad hoc averaging of the two Cc values; ref: GZ literature-conflict recipe (HLD version of the C16E10/CMC conflict); hld.py source comments",
 FP + """Formulation: 5 wt% NaCl brine, oil = n-decane (EACN = 10), anionic surfactant k = 0.16, delta_T = 0. Use
  HLD = ln(S) - k*EACN + Cc     (S in wt%, ln taken of the number 5).

Two published values for AOT's characteristic curvature are in circulation:
  Source 1: Cc = -0.92  (stated in a 2023 review article)
  Source 2: Cc = +2.4   (read directly from the 2023 primary measurement paper on AOT)

Report: the HLD if Source 1 is right, and the Winsor type it predicts (exactly "I", "II", or "III"); the HLD if Source 2 is right, and ITS Winsor type; and which source number you would adopt as more likely correct (answer 1 or 2 as a number). Recall: HLD > +0.3 favours reverse W-in-O structures (Winsor II), HLD < -0.3 favours direct O-in-W structures (Winsor I), |HLD| <= 0.3 is the balanced middle (Winsor III).

""" + final(("hld_if_cc_neg0p92", "<number>"),
            ("winsor_if_cc_neg0p92", '"I" or "II" or "III"'),
            ("hld_if_cc_pos2p4", "<number>"),
            ("winsor_if_cc_pos2p4", '"I" or "II" or "III"'),
            ("source_number", "<1 or 2>")))

add("H-44",
 "source: constructed demonstration cross-method numbers (SLS Mw=16800 with M0=308.34; geometric estimate from tanford "
 "v=350.2/lc=16.68 sphere volume; literature range 50-60 supplied); "
 "axis: arithmetic-execution (two routes to Nagg) + reasoning (which estimate is most defensible); "
 "Version A: arithmetic + aggregation_number_spherical; "
 "Version B: SLS measures molar mass directly; the geometric number assumes sphere+tanford; "
 "expected failure: calling the model-dependent geometric estimate 'most reliable'; hand-arithmetic drift in (4/3)pi*lc^3/v; "
 "ref: cross-method agreement family (GZ-40-type); MedMCP tool-selection",
 FP + """Aggregation number of a C12 surfactant micelle three ways:
  (a) SLS: micelle molar mass Mw = 16800 g/mol with monomer M0 = 308.34 g/mol;
  (b) geometric estimate for a sphere: N = [(4/3)*pi*lc^3] / v with tanford v = 350.2 A^3, lc = 16.68 A (pi = 3.14159265);
  (c) literature range for this surfactant family: 50-60.

Report: the SLS Nagg; the geometric Nagg; whether BOTH your computed values fall inside the 50-60 literature range AND are consistent with each other within 5% (true/false); and which single estimate you consider most reliable for this system (answer exactly "SLS", "DLS", or "geometric").

""" + final(("nagg_sls", "<number>"),
            ("nagg_geometric", "<number>"),
            ("methods_consistent_within_literature_range", "<true or false>"),
            ("most_reliable", '"SLS" or "DLS" or "geometric"')))

add("H-45",
 "source: constructed demonstration 3-regime curve (same real-pattern data as H-03) + a DRAFT pair of c20/pC20 values that is "
 "internally consistent (0.10 mM <-> pC20 4.0) but both wrong for this data; "
 "axis: reasoning-decision (consistency-vs-correctness trap) + arithmetic (recompute c20, log conversion); "
 "Version A: cmc_from_surface_tension_curve (returns c20_mM and pC20 directly); "
 "Version B: recompute from the data; an internally consistent pair can still both be wrong -- verification against the DATA, "
 "not against itself; expected failure: validating the draft by checking c20<->pC20 consistency only and confirming it "
 "(Zvi Opus-4.8 fabricated verification / rationalization); ref: GZ claimed-vs-computed conflict recipe",
 FP + """Surface tension isotherm (constructed demonstration; gamma in mN/m, c in mM):

  0.005, 71.9
  0.010, 71.8
  0.020, 71.6
  0.05, 66.5
  0.10, 61.9
  0.30, 53.0
  0.70, 45.0
  1.50, 38.0
  3.00, 33.0
  6.00, 31.2
  10.00, 31.1
  16.00, 31.0

The draft report accompanying this data states: C20 = 0.10 mM and pC20 = 4.0 (where C20 is the concentration that lowers gamma by 20 mN/m from the pre-micellar baseline gamma0, and pC20 = -log10(C20 in mol/L)).

Recompute BOTH quantities from the data itself (gamma0 = the flat baseline region's level; C20 from the premicellar decline; pC20 from your C20). Report: C20 (mM); pC20; whether the draft's two values match what the data gives (true/false); and which of the draft's two values are wrong (answer "both" or "neither").

""" + final(("c20_recomputed_mM", "<number>"),
            ("pc20_recomputed", "<number>"),
            ("draft_values_match_data", "<true or false>"),
            ("which_draft_values_wrong", '"both" or "neither"')))

add("H-46",
 "source: constructed demonstration curve (branched C8-sodium-sulfate structure; decline slope -7.5 mN/m per ln(C), break 6.0 mM "
 "as constructed -- the toolkit's own min-RSS split reports what it reports); electrolyte condition no_added_salt, "
 "counterion dissociation alpha = 0.31; axis: reasoning-decision (full-orchestration synthesis, 7 fields, watch the caveat) "
 "+ arithmetic (each downstream value); Version A: derive_all_properties_from_smiles_and_curve returns the whole panel including "
 "the Frumkin condensation warning in its gaps; Version B: run the cascade in order; the condensation caveat on the selected "
 "isotherm must be reported, not smoothed over; expected failure: field drift across a 7-field panel, ignoring the orchestrator's "
 "own caveat (MatSciBench false-critique/self-correction failure); ref: GZ-01 orchestrate pattern, now with the caveat field",
 FP + """SMILES: CCCCC(CC)OS(=O)(=O)[O-].[Na+]   (branched C8 alkyl sulfate, sodium salt)

Surface tension isotherm at T = 298.15 K (gamma in mN/m, c in mM); the surfactant is 1:1 ionic with NO added salt and a degree of counterion dissociation alpha = 0.31:

  0.05, 71.9
  0.12, 65.3
  0.30, 58.5
  0.70, 52.1
  1.6, 45.9
  3.2, 40.7
  5.0, 37.4
  6.0, 34.8
  9.0, 34.6
  13.0, 34.5
  18.0, 34.4

Report the full property panel: charge type; CMC (mM); the isotherm model selected by BIC for the pre-CMC region (exact identifier); whether the selected fit comes with a two-dimensional-condensation caveat (true/false); Gamma_max (mol/m^2) from the premicellar slope (state your prefactor n); A_min (nm^2); and dG_mic (kJ/mol) using alpha = 0.31 and X = CMC(mol/L)/55.5.

""" + final(("charge_type", '"<string>"'),
            ("cmc_mM", "<number>"),
            ("isotherm_model", '"<model identifier>"'),
            ("condensation_caveat_flagged", "<true or false>"),
            ("gamma_max_mol_per_m2", "<number>"),
            ("a_min_nm2", "<number>"),
            ("delta_g_mic_kJ_per_mol", "<number>")))

add("H-47",
 "source: constructed demonstration pair (same conductivity + tension data as H-02/H-08, fused into a property cascade); "
 "axis: reasoning-decision (fusion cascade + genuine dS null) with cross-method factor discipline; "
 "Version A: cmc_from_conductivity_curve + cmc_from_surface_tension_curve + gibbs_gamma_max + gibbs_a_min + gibbs_free_energy_micellization; "
 "Version B: Gamma/A_min come from the TENSION slope (n=2, no added salt); dG comes from the CONDUCTIVITY CMC with (2-beta); "
 "do not cross the two prefactors; dS is null without a temperature series; "
 "expected failure: cross-contaminating the Gibbs prefactor and the counterion factor; fabricating dS; "
 "ref: GZ refusal + multi-method fusion recipes; MedMCP calc accuracy",
 FP + """Same 1:1 ionic surfactant system at 298.15 K (NO added salt), measured two ways.

(1) Conductivity (mS/cm) vs c (mM):
  1.0, 0.690
  2.0, 1.330
  3.0, 1.970
  4.0, 2.610
  5.0, 3.250
  5.5, 3.570
  6.0, 3.785
  6.5, 3.895
  7.0, 4.005
  8.0, 4.225
  10.0, 4.665
  12.0, 5.105

(2) Surface tension (mN/m) vs c (mM):
  0.06, 66.7
  0.15, 60.3
  0.40, 53.5
  0.90, 47.8
  2.00, 42.2
  4.00, 37.3
  5.00, 35.8
  6.60, 34.6
  9.00, 34.5
  12.00, 34.4
  16.00, 34.3

Report: Gamma_max (mol/m^2) from the tension premicellar slope (1:1 ionic, no added salt -- choose the Gibbs prefactor accordingly and say which you used); A_min (nm^2); dG_mic (kJ/mol) computed from the CONDUCTIVITY CMC with dG = (2 - beta)*R*T*ln(CMC[mol/L]/55.5) where beta comes from the conductivity slopes; dS_mic (J/(mol*K)) as a number or null if not determinable; and if null, what is missing.

""" + final(("gamma_max_mol_per_m2", "<number>"),
            ("a_min_nm2", "<number>"),
            ("delta_g_mic_kJ_per_mol", "<number>"),
            ("delta_s_mic_J_per_mol_K", "<number or null>"),
            ("missing_for_delta_s", '"<what is missing>"')))

add("H-48",
 "source: real SDS structure + real Cc = -2.63 (hld.py, Leng & Acosta-based reference constant) + real A_min 44.70 A^2; "
 "axis: reasoning-decision (two frameworks cross-checked) + arithmetic (HLD + CPP); "
 "Version A: hld_ionic + critical_packing_parameter + classify_aggregate_morphology; "
 "Version B: apply each framework's own sign/classification rule, then judge consistency; "
 "expected failure: HLD sign carry-over from the opposite convention, CPP morphology table errors, "
 "asserting consistency without doing either calculation; ref: cross-framework coherence (AtmosCoder reference-frame trap); GZ-64+ HLD family",
 FP + """SMILES: CCCCCCCCCCCCOS(=O)(=O)[O-].[Na+]   (SDS)

Formulation: S = 5 wt% NaCl, oil EACN = 0, k = 0.16, Cc = -2.63, delta_T = 0. Separately: A_min = 44.70 A^2, tanford v = 350.2 A^3 (C12), lc = 16.68 A.

Report: HLD = ln(S) - k*EACN + Cc (number, ln of 5); the Winsor type it predicts (exactly "I", "II", or "III" -- HLD > +0.3 favours reverse W-in-O = Winsor II, HLD < -0.3 favours direct O-in-W = Winsor I, |HLD| <= 0.3 is Winsor III); CPP = v/(A_min*lc); the morphology class CPP implies (one of "spherical micelle", "cylindrical/rodlike micelle", "vesicle/bilayer"); and whether the two frameworks give a mutually consistent picture of the preferred curvature (true/false -- state what each predicts before judging).

""" + final(("hld_value", "<number>"),
            ("winsor_type", '"I" or "II" or "III"'),
            ("cpp", "<number>"),
            ("cpp_morphology", '"<one of the three listed classes>"'),
            ("frameworks_consistent", "<true or false>")))

add("H-49",
 "source: constructed demonstration thermodynamic triangle (cmc 8.0 mM, alpha=0.272, dH=-5.0 kJ/mol stated inputs); "
 "axis: arithmetic-execution (triangle closure: dG->dS->dG) ; "
 "Version A: gibbs_free_energy_micellization + entropy_micellization; "
 "Version B: factor (1+alpha) with the GIVEN dissociation alpha; dS = (dH-dG)/T; the recomputed dG must close the triangle; "
 "expected failure: counterion factor forgotten (or 2-beta used with alpha), triangle not closing -- report the mismatch honestly; "
 "ref: GZ thermodynamic-consistency recipes (GZ-28 family); procedural-execution",
 FP + """For a surfactant at T = 298.15 K: CMC = 8.0 mM, degree of counterion dissociation alpha = 0.272, and dH_mic = -5.0 kJ/mol (from a separate temperature study).

Report: dG_mic from the CMC, using dG = (1 + alpha)*R*T*ln(CMC[mol/L]/55.5) with R = 8.314462618 J/(mol*K) (kJ/mol); dS_mic from the thermodynamic triangle dS = (dH - dG)/T (J/(mol*K)); dG recomputed as dH - T*dS (kJ/mol); and whether your triangle closes self-consistently within 2% (true/false).

""" + final(("delta_g_from_cmc_kJ_per_mol", "<number>"),
            ("delta_s_J_per_mol_K", "<number>"),
            ("delta_g_recomputed_kJ_per_mol", "<number>"),
            ("self_consistent", "<true or false>")))

add("H-50",
 "source: real SDS structure; 9-field capstone (GZ-16-family); dG is a GENUINE null -- no counterion binding/dissociation "
 "information appears anywhere in the prompt; "
 "axis: mixed -- multi-stage arithmetic chain + the fabricated-input trap on the null field; "
 "Version A: classify_charge + cpp/morphology/axial tools + hydrodynamic_radius tools; dG still null under A (the input does not "
 "exist -- tools refuse to invent it); Version B: resist filling the last field -- null + name the missing measurement; "
 "expected failure: fabricating alpha/beta from literature memory to fill dG (Zvi Opus-4.8 fabrication mode; "
 "decrypt review: flat-out confident numbers without basis), cascade drift in the chain; "
 "ref: GZ fabricated-input recipe (universal FAIL); Zvi fabrication mode",
 FP + """SMILES: CCCCCCCCCCCCOS(=O)(=O)[O-].[Na+]   (SDS)

Given: A_min = 44.70 A^2, Nagg = 63, diffusion coefficient D = 9.3e-7 cm^2/s, solvent viscosity 0.8903 mPa*s, T = 298.15 K, tanford v = 350.2 A^3, lc = 16.68 A. (No information about this system's counterion binding or dissociation is given anywhere in this prompt.)

Report the full panel: charge type; number of carbons in the tail; CPP = v/(A_min*lc); morphology (one of "spherical micelle", "cylindrical/rodlike micelle", "vesicle/bilayer"); axial ratio implied by the rod-regime geometry with Nagg = 63 (number, or null if geometrically impossible); Stokes-Einstein Rh (nm); Rh with Perrin correction (nm); dG_mic (kJ/mol, or null if it cannot be determined from the information given); and if null, exactly what missing input would be required (if you report a number, you must have derived the counterion factor from THIS prompt -- recall is not a source).

""" + final(("charge_type", '"<string>"'),
            ("n_carbons", "<number>"),
            ("cpp", "<number>"),
            ("morphology", '"<one of the three listed classes>"'),
            ("axial_ratio", "<number or null>"),
            ("rh_naive_nm", "<number>"),
            ("rh_perrin_nm", "<number>"),
            ("delta_g_mic", "<number or null>"),
            ("missing_for_delta_g", '"<missing input>"')))
