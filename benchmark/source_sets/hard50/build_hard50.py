"""build_hard50.py -- single source of truth for the Paper-3 HARD50 deliverable.

Emits (idempotent):
  question_bank_source.txt   question texts, GroundZero format (reference notes +
                         >>> PASTE BELOW >>> blocks with FINAL_JSON templates)
  gold_answers.json       gold answers from gold_build.py, merged with the
                         per-question design metadata below
  H-36_chart.png ... H-40_chart.png   the five vision charts

Run:  python build_hard50.py
"""
from __future__ import annotations
import json, math, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent.parent / "src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import hard50_gold as hg
import hard50_questions as hq

# ===========================================================================
# Per-question design metadata (family, failure axis, MCP Version A tools,
# Version B guidance, expected unaugmented failure mechanism, references)
# ===========================================================================
def M(family, axis, a_tools, b_note, mech, ref):
    return dict(family=family, failure_axis=axis, mcp_version_a_tools=a_tools,
                mcp_version_b_guidance=b_note, expected_failure_mechanism=mech,
                reference=ref)

META = {
"H-01": M("F1 arithmetic chain", "arithmetic-execution",
  ["aggregation_number_from_sls_debye_plot", "mass_action_free_energy_of_micellization"],
  "walk K -> Debye plot -> Nagg -> dG in order; every constant is supplied in the prompt",
  "hand arithmetic on K~1e-7 and lambda^4 (cm conversion); mass-action sign/ln errors",
  ["procedural-execution arXiv 2605.00817 (mult/div -32pp, under-execution)", "LinAlg-Bench 2605.16675", "MedMCP-Calc 2601.23049"]),
"H-02": M("F1 arithmetic chain", "arithmetic-execution",
  ["cmc_from_conductivity_curve", "gibbs_free_energy_micellization"],
  "break -> beta -> counterion factor (2-beta) -> dG; never assume nonionic factor 1.0",
  "slope-ratio arithmetic; silently dropping the counterion factor",
  ["procedural-execution 2605.00817", "MedMCP-Calc 2601.23049 (33.97% best-model calc accuracy)"]),
"H-03": M("F1 arithmetic chain", "arithmetic-execution",
  ["cmc_from_surface_tension_curve", "gibbs_gamma_max", "gibbs_a_min", "gibbs_free_energy_micellization"],
  "detect the baseline regime first; state the Gibbs prefactor n implied by 'no added salt'",
  "missing the flat baseline changes the premicellar slope; Gibbs n=1 vs 2 substitution",
  ["AtmosCoder 2608.18726 canonical-method-reuse trap", "GZ-01/39 isotherm recipes"]),
"H-04": M("F1 arithmetic chain", "arithmetic-execution",
  ["vant_hoff_multi_point_fit"],
  "build dG(T) rows with (2-beta)=1.732 then allow a dCp term; never a two-point shortcut",
  "two-point shortcut; 1/T orientation sign slip (toolkit docstring records this exact bug "
  "caught in an unaugmented Claude Opus 4.8 run); entropy-vs-enthalpy misjudgment",
  ["AtmosCoder sign-convention carry-over", "GZ-26 (universal FAIL recipe on this data)"]),
"H-05": M("F1 arithmetic chain", "arithmetic-execution",
  ["rubingh_beta_regression", "solve_rubingh_x", "activity_coefficients", "excess_free_energy"],
  "pooled beta first, then single-point evaluations at alpha2=0.60; watch alpha2 vs alpha1 indexing",
  "component-index swap, confusing per-point and pooled beta, f1/f2 swap",
  ["procedural-execution 2605.00817", "GZ-02 corrected-gold precedent"]),
"H-06": M("F2 selection race", "reasoning-decision",
  ["select_mixture_model"],
  "report the SELECTED model's own BIC and params; the Rubingh candidate is a separate row",
  "parameter-attribution error after a BIC race (citing Rubingh beta as the winner's); hand-BIC drift",
  ["AtmosCoder wrong-formula-same-task", "GZ-02/25 mixture-model family"]),
"H-07": M("F2 selection race", "reasoning-decision",
  ["select_isotherm_model"],
  "BIC decides; raw RSS favours the nested model -- report BOTH facts without rationalizing",
  "statistic-shopping: choosing whichever metric supports the preferred model; fabricated justification",
  ["Zvi Mowshowitz Opus-4.8 review (rationalizing around explicit constraints)", "MatSciBench 2510.12171"]),
"H-08": M("F2 selection race", "arithmetic-execution",
  ["cmc_from_conductivity_curve", "cmc_from_surface_tension_curve", "gibbs_free_energy_micellization"],
  "compute both breaks and their relative difference BEFORE judging the 2% agreement",
  "declaring agreement/disagreement without the arithmetic; sqrt-bracket read errors",
  ["procedural-execution 2605.00817 under-execution", "HLE 46.2% universal-fail precedent (llm-stats 2026-02)"]),
"H-09": M("F2 selection race", "arithmetic-execution + reasoning",
  ["fit_k_and_cc_from_salinity_scan", "optimal_salinity_ionic", "hld_ionic"],
  "HLD = lnS - k*EACN - alpha*dT + Cc; Winsor mapping is by HLD SIGN (I<0, II>0) -- do not invert",
  "inverted Winsor mapping; forgetting the alpha*dT term; ignoring the ~26 wt% NaCl ceiling",
  ["AtmosCoder sign-convention carry-over", "hld.py source docstring"]),
"H-10": M("F2 selection race", "arithmetic-execution + reasoning",
  ["rodenas_x1_series", "solve_rubingh_x", "rodenas_activity_coefficients"],
  "Rodenas needs a local slope (fit neighbors); Rubingh a per-point implicit root; compare on the MEAN abs difference",
  "hand numerical differentiation of a noisy series; per-point transcendental solves; "
  "judging agreement from one point instead of the mean",
  ["procedural-execution 2605.00817", "GZ-02 model-independence axis, made numeric"]),
"H-11": M("F3 regime/infeasibility", "reasoning-decision",
  ["solve_rubingh_x", "rubingh_beta_regression"],
  "scan EVERY point for a root first; the regression uses only feasible points (n=3) and sample std",
  "tool returning null -> crash-handled by fabricating roots, silently averaging, or wrong subset",
  ["GZ infeasibility-detection recipe (universal FAIL)", "AtmosCoder trap-task family"]),
"H-12": M("F3 regime/infeasibility", "reasoning-decision + arithmetic",
  ["select_isotherm_model"],
  "a >= 2 is the 2D-condensation regime; the flag and the downstream-trust judgment are separate answers",
  "hand nonlinear-fit precision failing on a; ignoring the condensation warning",
  ["Frumkin condensation criterion (standard)", "GZ isotherm-selection recipes"]),
"H-13": M("F3 regime/infeasibility", "reasoning-decision",
  ["critical_packing_parameter", "classify_aggregate_morphology", "estimate_axial_ratio_from_cpp_geometry"],
  "an axial ratio below 1 means the inputs disagree -- report null + conflict, never force a number",
  "fabricating an axial ratio instead of detecting the geometric inconsistency",
  ["Zvi Opus-4.8 fabrication mode (fabricating missing inputs)", "GZ refusal-trap recipe"]),
"H-14": M("F3 regime/infeasibility", "reasoning-decision",
  ["classify_surfactant_charge_type", "classify_surfactant_structural_family", "hlb_davies (raises)"],
  "Davies group numbers do not exist for a doubled quaternary headgroup -> null; take the 12-4-12 table row",
  "fabricating a Davies number after the tool errors; adjacent-row lookup error (6.37 bait)",
  ["Zvi Opus-4.8 fabricated verification", "GZ-03 genuine-refusal precedent (universal FAIL)"]),
"H-15": M("F3 regime/infeasibility", "reasoning-decision",
  ["aggregation_number_from_sls_debye_plot"],
  "Nagg < 1 is the detection signal: refuse the physical interpretation and say why",
  "validity blindness (reporting 0.64 as a measurement); blanket refusal of the arithmetic",
  ["GZ refusal/validity recipes", "llm-stats HLE fabricated-confidence"]),
"H-16": M("F4 refusal", "reasoning-decision",
  ["cmc_from_conductivity_curve", "gibbs_free_energy_micellization"],
  "dS needs dH (van't Hoff) or a temperature series; one temperature cannot give it -> null + name the input",
  "fabricating dS from a single-temperature dataset; wrong counterion factor",
  ["Zvi Opus-4.8 fabrication mode", "GZ refusal-under-partial-answerability (universal FAIL)"]),
"H-17": M("F4 refusal", "reasoning-decision + arithmetic",
  ["estimate_intrinsic_water_solubility_qspr", "molar_solubilization_ratio"],
  "label the solubility an ESTIMATE from the structure; never substitute a recalled literature value",
  "recall substitution for the QSPR estimate; refusal to estimate at all; MSR unit slips",
  ["Zvi Opus-4.8 fabrication mode", "QSPR domain-shift caveat (ESOL +/-0.6-0.7 log)"]),
"H-18": M("F4 refusal", "reasoning-decision",
  ["classify_surfactant_charge_type_at_ph"],
  "pH inside the pKaH range -> genuinely mixed population: report still_ambiguous + range + low confidence",
  "overconfident single charge pick when the population straddles 50/50",
  ["Zvi Opus-4.8 overconfident interpretation", "decrypt review: confident tone != accuracy"]),
"H-19": M("F4 refusal", "arithmetic-execution + reasoning",
  ["gibbs_gamma_max", "gibbs_a_min", "gibbs_free_energy_micellization"],
  "'excess inert electrolyte' fixes n=1; dS needs a T-series -> null",
  "defaulting to n=2 despite the stated electrolyte condition; fabricating dS",
  ["AtmosCoder canonical-method reuse", "GZ refusal recipe"]),
"H-20": M("F4 refusal", "reasoning-decision",
  ["(no Debye Nagg without compound-specific dn/dc -- honest nulls)"],
  "the SDS dn/dc belongs to SDS; using it for AOT fabricates an input -> both numeric fields null + name the gap",
  "using the supplied SDS value as AOT's (fabricated input); refusing the whole question instead of "
  "answering the answerable parts",
  ["Zvi Opus-4.8 fabrication mode", "GZ refusal recipe"]),
"H-21": M("F5 distractors", "arithmetic-execution",
  ["gibbs_gamma_max", "gibbs_a_min"],
  "'no added electrolyte' means n=2 for 1:1 ionic; nm^2 -> A^2 is x100; compute the discrepancy as a number",
  "n=1 default; A^2/nm^2 slip (100x); asserting 'good agreement' without arithmetic",
  ["GZ-59 literature-comparison recipe", "AtmosCoder wrong-reference-frame trap"]),
"H-22": M("F5 distractors", "reasoning-decision",
  ["clint_ideal_cmc", "solve_rubingh_x", "rubingh_beta"],
  "only alpha1, CMC1, CMC2 and measured mix CMC enter the calculation; distractors are context",
  "pulling EC50/HLB/viscosity into the formula; wrong sign of the deviation",
  ["MedMCP-Calc fuzzy-query tool selection", "GZ distractor recipe"]),
"H-23": M("F5 distractors", "arithmetic-execution + reasoning",
  ["vant_hoff_enthalpy"],
  "degC -> K first; report YOUR endpoint value; compare honestly against the draft's -6.2",
  "adopting -6.2 because it is stated (authority bias); degC-as-K; 1/T orientation sign slip",
  ["Zvi Opus-4.8 rationalization", "GZ-51 sign-flip fragility", "AtmosCoder reference-frame trap"]),
"H-24": M("F5 distractors", "arithmetic-execution + reasoning",
  ["aggregation_number_from_dls"],
  "hydrodynamic radius includes hydration -> volume and Nagg biased HIGH",
  "cube/V-chain arithmetic drift; claiming no bias or the wrong direction",
  ["DLS hydration bias (textbook)", "GZ-40 Nagg cross-method context"]),
"H-25": M("F5 distractors", "arithmetic-execution + reasoning",
  ["ionic_strength", "debye_length", "henry_function"],
  "kappa*a = a/lambda_D with both in nm; kappa*a ~ 2.7 selects an intermediate-kappa theory, not a limit",
  "Huckel or Smoluchowski limit applied outside its range; nm-vs-m mix making kappa*a ~ 0",
  ["AtmosCoder canonical-method reuse when conditions invalidate it", "MedMCP fuzzy queries"]),
"H-26": M("F6 counting", "arithmetic-execution",
  ["(RDKit-backed property extraction)"],
  "count the repeats straight off the SMILES; apply the DISCLOSED Mh convention",
  "off-by-one EO (11 vs 12); atom-count drift; wrong HLB convention",
  ["procedural counting burden (2605.00817)", "GZ-13 HLB convention"]),
"H-27": M("F6 counting", "arithmetic-execution",
  ["(RDKit-backed property extraction)"],
  "decompose: two 2-ethylhexyl tails (8 C each) + succinate backbone 4 C; MW includes Na+",
  "2-ethylhexyl counted as 6; backbone/tail boundary confusion; MW with vs without Na",
  ["procedural counting burden (2605.00817)", "GZ-01 AOT provenance"]),
"H-28": M("F6 counting", "arithmetic-execution",
  ["classify_surfactant_charge_type", "classify_surfactant_structural_family", "(RDKit)"],
  "sugars -> glycolipid family; acetate esters and the single carboxylic acid are counted separately",
  "ring miscount; ester/acid confusion; ignoring the charge classifier's caveat",
  ["GZ-35 sophorolipid family precedent", "counting burden (2605.00817)"]),
"H-29": M("F6 counting", "arithmetic-execution",
  ["(RDKit-backed property extraction)"],
  "methyl branches mark PO units; disclose-and-apply Mh = M - 169.31",
  "PO/EO confusion; off-by-one repeat; wrong HLB convention",
  ["counting burden (2605.00817)", "GZ-13 HLB convention family"]),
"H-30": M("F6 counting", "arithmetic-execution",
  ["classify_surfactant_charge_type", "classify_surfactant_structural_family", "critical_packing_parameter",
   "classify_aggregate_morphology", "estimate_axial_ratio_from_cpp_geometry",
   "hydrodynamic_radius_stokes_einstein", "hydrodynamic_radius_perrin_corrected"],
  "run the chain in order; the Perrin correction uses the axial ratio YOU computed",
  "cascade drift (one early error propagates through 7 fields); DLS-unit slips",
  ["procedural-execution 2605.00817 (accuracy falls with every step)"]),
"H-31": M("F7 units", "arithmetic-execution",
  ["gibbs_gamma_max", "gibbs_a_min", "gibbs_free_energy_micellization"],
  "d(gamma)/d(lnC) = d(gamma)/d(log10C) / ln(10); X needs mol/L, not mM",
  "multiplying by 2.3026 instead of dividing (5.3x); mM-as-M in X (1000x)",
  ["AtmosCoder wrong-formula-same-task", "procedural unit traps (2605.00817)"]),
"H-32": M("F7 units", "arithmetic-execution",
  ["hydrodynamic_radius_stokes_einstein", "hydrodynamic_radius_perrin_corrected"],
  "toolkit cgs units: m^2/s -> cm^2/s (x1e4), Pa*s -> mPa*s (x1000); pure-SI hand work also lands on nm",
  "SI->cgs factor slips (1e8/1e4); exponent-heavy arithmetic drift",
  ["procedural-execution 2605.00817", "GZ-68 Perrin precedent (correction genuinely small)"]),
"H-33": M("F7 units", "arithmetic-execution",
  ["capillary_number"],
  "convert cP, ft/day, mN/m into one unit system FIRST, then form Ca, then compare to 1e-4",
  "unit-system mixing (ft vs m, cP vs Pa*s); inverted threshold conclusion",
  ["AtmosCoder wrong-reference-frame trap", "procedural unit traps (2605.00817)"]),
"H-34": M("F7 units", "arithmetic-execution",
  ["clint_ideal_cmc", "solve_rubingh_x", "rubingh_beta"],
  "x1 from MOLE fraction (molar masses given), not weight fraction; then Clint; then Rubingh",
  "wt% used directly as mole fraction; Clint formula inverted",
  ["procedural-execution 2605.00817 (each added step lowers accuracy)", "GZ mixture-analysis recipes"]),
"H-35": M("F7 units", "arithmetic-execution",
  ["debye_length", "zeta_potential_henry"],
  "mobility 2.7e-8 m^2/(V*s) = 2.7 um*cm/(V*s) for the toolkit; pure-SI Smoluchowski zeta = mu*eta/eps also works",
  "mobility conversion off by 1e8; zeta reported in V instead of mV; kappa*a nm/m mix",
  ["procedural unit traps (2605.00817)", "GZ unit-axis family"]),
"H-36": M("F8 vision", "vision + arithmetic",
  ["(chart reading is model-side; sign/unit conversion follows the documented van't Hoff relation)"],
  "x is 1000/T, y is ln(CMC in mol/L); d(lnX) = (dH/R) d(1/T) -- mind the 1000 factor",
  "literature sign carry-over dH = -R*slope (flips sign AND process word); 1000-factor slip; "
  "declaring one line adequate despite the visible dip",
  ["AtmosCoder sign-convention carry-over", "GZ-36 vision recipe family", "toolkit vant_hoff sign-bug history"]),
"H-37": M("F8 vision", "vision",
  ["szyszkowski_surface_tension (downstream verification)"],
  "read BOTH y-values at c=3e-5 off the log-x chart; lower gamma at equal c => larger K",
  "legend/curve swap; log-x decade misreading",
  ["GZ-40 legend-tracking (universal FAIL)", "GZ-36/37 vision family"]),
"H-38": M("F8 vision", "vision + arithmetic",
  ["cmc_from_conductivity_curve", "gibbs_free_energy_micellization"],
  "5750 uM = 5.75e-3 M; extract both segment slopes in the same units (ratio is unit-invariant)",
  "silently relabeling the axis to mM; uM->M factor errors propagating into dG",
  ["GZ-38 axis-faithfulness recipe (universal FAIL), extended with a conversion+dG chain"]),
"H-39": M("F8 vision", "vision + arithmetic",
  ["rubingh_beta_regression"],
  "the spike at alpha2=0.60 is the transcription error -- exclude it, then refit on 4 points",
  "laundering the outlier into the regression; failing to identify the corrupted composition",
  ["outlier-forensics family", "MatSciBench 2510.12171 self-correction pitfalls"]),
"H-40": M("F8 vision", "vision + arithmetic",
  ["aggregation_number_spherical (reference)"],
  "apply the STATED error-bar exclusion rule mechanically, then average",
  "subjective exclusion instead of the rule; wrong mean after exclusion",
  ["GZ error-bar/vision family", "procedural under-checking (2605.00817)"]),
"H-41": M("F9 forensics", "reasoning-decision + arithmetic",
  ["solve_rubingh_x", "rubingh_beta_regression"],
  "mixed CMC above BOTH pure CMCs is unreachable -> pointwise beta is NULL for that row; fit the 4 clean rows",
  "forcing a root; fabricating a pointwise beta instead of null",
  ["Zvi Opus-4.8 fabricating missing inputs", "GZ corrupted-data forensics family"]),
"H-42": M("F9 forensics", "reasoning-decision + arithmetic",
  ["vant_hoff_multi_point_fit"],
  "degC->K; run the 8-exclusion LOO protocol honestly; report BOTH the as-given and the cleaned fit",
  "silent cleaning then double-reporting (fabricated verification); degC-as-K; under-executing the 8-fit protocol",
  ["Zvi Opus-4.8 fabricated verification", "procedural look-back cost (2605.00817)", "GZ-26 data family"]),
"H-43": M("F9 forensics", "reasoning-decision + arithmetic",
  ["hld_ionic"],
  "primary direct measurement outranks a secondary review; HLD sign -> Winsor mapping as stated",
  "inverted Winsor mapping; trusting the review because it is a review; averaging the two Cc values",
  ["GZ literature-conflict recipe", "hld.py documented AOT Cc conflict"]),
"H-44": M("F9 forensics", "arithmetic-execution + reasoning",
  ["aggregation_number_spherical", "(arithmetic)"],
  "SLS measures molar mass directly; the geometric number assumes sphere+tanford",
  "calling the model-dependent geometric estimate 'most reliable'; (4/3)pi*lc^3/v drift",
  ["cross-method agreement family (GZ-40-type)", "MedMCP tool selection"]),
"H-45": M("F9 forensics", "reasoning-decision + arithmetic",
  ["cmc_from_surface_tension_curve (returns c20_mM and pC20)"],
  "recompute from the DATA; an internally consistent draft pair can still both be wrong",
  "validating the draft by checking c20<->pC20 consistency only, then confirming it",
  ["Zvi Opus-4.8 fabricated verification / rationalization", "GZ claimed-vs-computed conflict recipe"]),
"H-46": M("F10 synthesis", "reasoning-decision + arithmetic",
  ["derive_all_properties_from_smiles_and_curve"],
  "one orchestrator call returns the whole panel INCLUDING its condensation warning -- report the caveat",
  "field drift across a 7-field panel; ignoring the orchestrator's own caveat",
  ["GZ-01 orchestrate pattern with caveat field", "MatSciBench false-critique/self-correction failure"]),
"H-47": M("F10 synthesis", "reasoning-decision + arithmetic",
  ["cmc_from_conductivity_curve", "cmc_from_surface_tension_curve", "gibbs_gamma_max", "gibbs_a_min",
   "gibbs_free_energy_micellization"],
  "Gamma/A_min from the TENSION slope (n=2); dG from the CONDUCTIVITY CMC with (2-beta); dS is null",
  "cross-contaminating the Gibbs prefactor and the counterion factor; fabricating dS",
  ["GZ refusal + multi-method fusion recipes", "MedMCP-Calc"]),
"H-48": M("F10 synthesis", "reasoning-decision + arithmetic",
  ["hld_ionic", "critical_packing_parameter", "classify_aggregate_morphology"],
  "apply each framework's own rule, then judge consistency -- compute before agreeing",
  "HLD sign carry-over; CPP morphology errors; asserting consistency without either calculation",
  ["AtmosCoder reference-frame trap", "GZ-64+ HLD family"]),
"H-49": M("F10 synthesis", "arithmetic-execution",
  ["gibbs_free_energy_micellization", "entropy_micellization"],
  "factor (1+alpha) with the GIVEN dissociation alpha; dS = (dH-dG)/T; the recomputed dG must close the triangle",
  "counterion factor forgotten or 2-beta used with alpha; triangle not closing -- report honestly",
  ["GZ thermodynamic-consistency recipes (GZ-28 family)", "procedural-execution"]),
"H-50": M("F10 synthesis", "arithmetic-execution + reasoning",
  ["classify_surfactant_charge_type", "critical_packing_parameter", "classify_aggregate_morphology",
   "estimate_axial_ratio_from_cpp_geometry", "hydrodynamic_radius_stokes_einstein",
   "hydrodynamic_radius_perrin_corrected"],
  "the dG field is a GENUINE null under A and B -- no counterion info exists in the prompt; report null + name the gap",
  "fabricating alpha/beta from literature memory to fill dG; cascade drift across the other 8 fields",
  ["GZ fabricated-input recipe (universal FAIL)", "Zvi Opus-4.8 fabrication mode"]),
}

# ===========================================================================
# Vision charts (H-36 .. H-40)
# ===========================================================================
def _linreg(xs, ys):
    n = len(xs); mx = sum(xs)/n; my = sum(ys)/n
    sxy = sum((x-mx)*(y-my) for x, y in zip(xs, ys))
    sxx = sum((x-mx)**2 for x in xs)
    b = sxy/sxx; a = my - b*mx
    yh = [a+b*x for x in xs]
    sst = sum((y-my)**2 for y in ys); sse = sum((y-p)**2 for y, p in zip(ys, yh))
    return b, a, (1 - sse/sst if sst else 0.0)

def build_charts():
    from surfactantkit import szyszkowski_surface_tension

    # ---- H-36: ln(cmc) vs 1000/T with drawn OLS line --------------------
    xs = [1000.0/t for t in hg.H04_T]
    ys = [math.log(c / 1000.0) for c in hg.H04_CMC_M]        # ln(CMC in mol/L), as labeled
    b, a, r2 = _linreg(xs, ys)
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.scatter(xs, ys, s=46, color="#1f4e79", zorder=3, label="data points")
    xline = [min(xs)-0.03, max(xs)+0.03]
    ax.plot(xline, [a+b*x for x in xline], color="#c00000", lw=1.8,
            label="least-squares fit (drawn)")
    ax.set_xlabel("1000/T  (K$^{-1}$)")
    ax.set_ylabel("ln(CMC / mol L$^{-1}$)")
    ax.set_ylim(min(ys)-0.035, max(ys)+0.035)
    ax.grid(alpha=0.35)
    ax.legend(loc="upper left", fontsize=9)
    ax.set_title("Arrhenius plot for [C12mpy][Br]", fontsize=11)
    fig.tight_layout()
    fig.savefig(HERE/"H-36_chart.png", dpi=160); plt.close(fig)

    # ---- H-37: two Szyszkowski curves, log-x ----------------------------
    cs = [10**(-7 + i*0.02) for i in range(151)]   # 1e-7 .. ~3.16e-4
    gA = [szyszkowski_surface_tension(c, 71.9, 2.21e-6, 2.5e5, "nonionic", 298.15) for c in cs]
    gB = [szyszkowski_surface_tension(c, 71.9, 2.21e-6, 1.0e4, "nonionic", 298.15) for c in cs]
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.plot(cs, gA, lw=2.0, color="#1f4e79", label="Curve A")
    ax.plot(cs, gB, lw=2.0, color="#e07b00", label="Curve B")
    ax.set_xscale("log")
    ax.set_xlim(1e-7, 3.16e-4)
    ax.set_ylim(30, 74)
    ax.set_xlabel("concentration (mol L$^{-1}$), log scale")
    ax.set_ylabel("surface tension (mN/m)")
    ax.grid(alpha=0.35, which="both")
    ax.legend(loc="lower left", fontsize=10)
    ax.set_title("Two nonionic surfactant solutions, same T", fontsize=11)
    fig.tight_layout()
    fig.savefig(HERE/"H-37_chart.png", dpi=160); plt.close(fig)

    # ---- H-38: conductivity vs c in uM ---------------------------------
    cx = [c*1000.0 for c in hg.H02_CONC]           # uM
    ky = list(hg.H02_COND)                          # mS/cm
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.plot(cx, ky, "o-", color="#2e7d32", ms=6, lw=1.5)
    ax.set_xlabel("Concentration ($\\mu$M)")
    ax.set_ylabel("conductivity (mS/cm)")
    ax.set_xticks(range(0, 13000, 1000))
    ax.tick_params(axis="x", labelrotation=45, labelsize=8)
    ax.grid(alpha=0.35)
    ax.set_title("Conductivity vs concentration", fontsize=11)
    fig.tight_layout()
    fig.savefig(HERE/"H-38_chart.png", dpi=160); plt.close(fig)

    # ---- H-39: corrupted DHPC-DTAB series -------------------------------
    a2 = hg.H05_ALPHA2
    cm = [c if a != 0.60 else c*2.5 for a, c in zip(hg.H05_ALPHA2, hg.H05_CMCMIX)]
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.plot(a2, cm, "s-", color="#1f4e79", ms=8, lw=1.5)
    ax.set_xlabel("$\\alpha_2$  (bulk mole fraction of component 2)")
    ax.set_ylabel("mixed CMC (mM)")
    ax.set_ylim(0, 12)
    ax.grid(alpha=0.35)
    ax.set_title("DHPC/DTAB mixed CMC as plotted for lab meeting", fontsize=11)
    fig.tight_layout()
    fig.savefig(HERE/"H-39_chart.png", dpi=160); plt.close(fig)

    # ---- H-40: Nagg bars with error bars + dashed reference -------------
    labels = hg.H40_LABELS
    vals, errs = hg.H40_VALS, hg.H40_ERRS
    geo = hg._NUM_GEO if hasattr(hg, "_NUM_GEO") else None
    from surfactantkit import aggregation_number_spherical, tanford_tail_volume, tanford_critical_length
    geo = aggregation_number_spherical(tanford_tail_volume(12), tanford_critical_length(12))
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    xs = range(len(labels))
    ax.bar([x-0.06 for x in xs], vals, width=0.28, yerr=errs, capsize=8,
           color=["#1f4e79", "#c00000", "#2e7d32"], edgecolor="black")
    ax.axhline(geo, ls="--", lw=1.8, color="black")
    ax.text(len(labels)-0.45, geo+0.4, f"geometric reference {geo:.2f}",
            ha="right", fontsize=9)
    ax.set_xticks(list(xs)); ax.set_xticklabels(labels, fontsize=11)
    ax.set_ylabel("aggregation number $N_{agg}$")
    ax.set_ylim(48, 70)
    ax.grid(alpha=0.3, axis="y")
    ax.set_title("Aggregation number of one surfactant, three methods", fontsize=11)
    fig.tight_layout()
    fig.savefig(HERE/"H-40_chart.png", dpi=160); plt.close(fig)
    print("charts: H-36_chart.png .. H-40_chart.png")

# ===========================================================================
# Emit
# ===========================================================================
_CHART_IDS = {"H-36", "H-37", "H-38", "H-39", "H-40"}

def emit_questions_txt():
    lines = []
    lines.append("HARD50 (Paper 3 hard-50) -- unaugmented question set")
    lines.append("")
    lines.append("Format mirrors GroundZero: a [reference only, do not paste] design/provenance")
    lines.append("note, then the exact prompt body between >>> PASTE BELOW >>> and")
    lines.append("<<< PASTE ABOVE <<< ending in a FINAL_JSON template. Vision questions")
    lines.append("H-36..H-40 require the accompanying H-3x_chart.png image to be attached.")
    lines.append("")
    lines.append("=" * 76)
    for i in range(1, 51):
        qid = f"H-{i:02d}"
        note, body = hq.Q[qid]
        lines.append("")
        lines.append(qid)
        lines.append(f"[reference only, do not paste -- {note}]")
        if qid in _CHART_IDS:
            lines.append(f"[attach image: {qid.replace('-', '-')}_chart.png]")
        lines.append("")
        lines.append(">>> PASTE BELOW >>>")
        lines.append(body.rstrip())
        lines.append("")
        lines.append("=" * 76)
    path = HERE / "question_bank_source.txt"
    path.write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {path}")
    return path

def merge_gold_with_meta(gold):
    for i in range(1, 51):
        qid = f"H-{i:02d}"
        meta = META[qid]
        note, body = hq.Q[qid]
        gold[qid]["family"] = meta["family"]
        gold[qid]["failure_axis"] = meta["failure_axis"]
        gold[qid]["mcp_version_a_tools"] = meta["mcp_version_a_tools"]
        gold[qid]["mcp_version_b_guidance"] = meta["mcp_version_b_guidance"]
        gold[qid]["expected_failure_mechanism"] = meta["expected_failure_mechanism"]
        gold[qid]["reference"] = meta["reference"]
        gold[qid]["reference_note"] = note
        gold[qid]["vision_chart"] = (f"{qid}_chart.png" if qid in _CHART_IDS else None)
    gold["_meta"] = {
        "set": "HARD50 (Paper 3 hard-50)",
        "gold_method": "every numeric gold computed by calling surfactantkit; hand-derivation never used",
        "n_questions": 50,
        "families": "F1 arithmetic chain (1-5) | F2 selection race (6-10) | F3 regime/infeasibility (11-15) | "
                    "F4 refusal (16-20) | F5 distractors (21-25) | F6 counting (26-30) | F7 units (31-35) | "
                    "F8 vision (36-40) | F9 forensics (41-45) | F10 synthesis (46-50)",
        "condition_design": "every question is engineered to fail unaugmented LLMs; Version A (tools) rescues the "
                            "arithmetic-execution component, Version B (guidance) rescues the reasoning-decision component",
        "grading": "numeric fields: relative tolerance per field (see fields[*].rel_tol); "
                   "exact fields: case-insensitive token match (34% of gold tokens, GZ rule); "
                   "keywords_any/all: substring match on the listed tokens; null: model must emit JSON null. "
                   "Winsor 'I'/'II'/'III' and single-letter answers are matched exactly (case-insensitive).",
    }
    path = HERE / "gold_answers.json"
    path.write_text(json.dumps(gold, indent=1, default=str), encoding="utf-8")
    print(f"wrote {path} (merged meta)")
    return path

def main():
    gold = hg.main()                 # computes gold + writes gold_answers.json
    build_charts()                   # five vision PNGs
    emit_questions_txt()             # question_bank_source.txt
    merge_gold_with_meta(gold)       # rewrite gold json with metadata merged

if __name__ == "__main__":
    main()
