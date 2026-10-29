"""Machine-readable validation-tier registry (audit TOOLKIT_AUDIT_2026-09-30
P1-4).

GENERATED from the audit's "### 1.1 Full technique -> tier table" rows:
every row's tool names become keys, the Tier cell becomes "tier", and the
Evidence cell becomes "evidence" (verbatim, so benchmark consumers can see
the citation / gold / test backing each number).  If the audit table
changes, regenerate this dict -- tests/test_validation_tiers.py parses the
same markdown and fails on any drift.

Tier definitions (transcribed from the audit's "### Tier definitions"):
    T1: Primary-source numeric anchor: a real paper's own reported number is reproduced (or the chain is digit-level matched). *Safe for \"match-the-literature\" gold claims.*
    T2: Source-verified formula, round-trip/proof-only output: the equation was verified against a real source before implementation, but **no independent published numeric example of the output exists** — validation is round-trip, limit identities, hand derivations, or plausibility ranges. *The \"difficult to generate literature-backed outputs\" set the audit asked for. Safe for computation golds; unsafe for \"the literature says\" claims.*
    T3: Derived/estimate by construction: explicitly non-literature (self-consistency calibration or statistical model), disclosed at the call site.
    T4: Documented dead end: the honest output is \"cannot compute\" — see §1.3.

Plus, per the table's own usage and its footnote "dagger = strong but
qualified (see cell)":
    "T1\u2020" (T1-dagger) = strong-but-qualified T1 -- primary anchor holds
        for the main quantity; the row's evidence cell carries the
        qualification
    "T1\u2020/T2" = qualified T1 for the primary case, T2 for the other
        variant named in the row (vant_hoff_enthalpy)
    "dependent" = follows from its inputs; validated only as an
        internal-consistency delta (no independent anchor)
    "chained"   = inherits the tiers of its callee tools

Wildcard keys ("tanford_*", "*_reference", ...) mirror grouped cells of the
audit table exactly; get_tool_tier() resolves names via fnmatch so both
exact names and wildcards work.  T4 documented dead ends (audit 1.3) are
intentionally NOT registered -- those techniques ship no tool at all.

Private API surface: exported from the package for the benchmark generator
and the tutorial; intentionally not MCP tools.
"""

from __future__ import annotations

from fnmatch import fnmatchcase
from typing import Final

VALIDATION_TIERS: Final[dict[str, dict[str, str]]] = {
    "clint_ideal_cmc": {"tier": "T1", "evidence": "Exact on 6+ published systems (PMC6554738, PMC4087020, …)"},
    "rubingh_solve": {"tier": "T1", "evidence": "~exact across 6 independent systems"},
    "rubingh_beta": {"tier": "T1†", "evidence": "Exact for fixed 1:1 (cholate-SDS); 5–20 % gap vs *regression* betas documented as method difference, not bug"},
    "rubingh_beta_regression": {"tier": "T2", "evidence": "Round-trip construction only; no paper publishes the same pointwise-mean procedure with raw data"},
    "motomura_ideal_composition": {"tier": "T2", "evidence": "Algebraic-identity tests only"},
    "rodenas_x1": {"tier": "T1†", "evidence": "Exact zero-intercept identity test; real-data series check exists"},
    "rodenas_x1_from_series": {"tier": "T2", "evidence": "Real data tested, but endpoints err **19–60 %** (paper uses a different global-fit differentiation) — documented limit, sparse-series output hard to anchor"},
    "asymmetric_margules": {"tier": "T1", "evidence": "Schulz & Durand SI CSVs digit-level matched (mmc2/3/6/8/10)"},
    "excess_free_energy": {"tier": "T1", "evidence": "Schulz & Durand SI CSVs digit-level matched (mmc2/3/6/8/10)"},
    "rodenas_activity_coefficients": {"tier": "T1", "evidence": "Schulz & Durand SI CSVs digit-level matched (mmc2/3/6/8/10)"},
    "eommm_global_fit": {"tier": "T1†", "evidence": "Formula confirmed exactly vs SI-fitted curves; real Hyamine/DTAB convergence; margin-dependent uniqueness disclosed"},
    "eommm_find_minimal_feasible_margin": {"tier": "T2", "evidence": "Round-trip + real-data convergence; the tightest-margin procedure itself has no published number to match"},
    "select_mixture_model": {"tier": "T2", "evidence": "BIC on synthetic constructions; scope choices (why Rosen/Maeda/Rodenas sit out) documented in the tool"},
    "maeda_free_energy_of_micellization": {"tier": "T1", "evidence": "Rub et al. 2021 (~6 %, attributed to paper rounding); component-indexing ambiguity honestly disclosed"},
    "rosen_monolayer_solve": {"tier": "T2", "evidence": "Verified mathematically identical to Rubingh (by design — a meaning distinction, not an independent anchor)"},
    "rosen_beta_sigma": {"tier": "T2", "evidence": "Verified mathematically identical to Rubingh (by design — a meaning distinction, not an independent anchor)"},
    "corrin_harkins_predict": {"tier": "T2", "evidence": "Round-trip + salt monotonicity + input guards only"},
    "gibbs_gamma_max": {"tier": "T1", "evidence": "4 real ionic systems ≤ 0.5 %"},
    "gibbs_a_min": {"tier": "T1", "evidence": "5 real systems ≤ 0.4 %"},
    "szyszkowski_surface_tension": {"tier": "T1", "evidence": "C12E8 (Ho/Lips/Velikov 2023) curve; Γ-round-trip 0.02 %"},
    "szyszkowski_fit_K": {"tier": "T2", "evidence": "Round-trip only — the one real candidate K table (Hreczuch 2024) is provably unusable (OCR/exponent contradiction, reconfirmed twice)"},
    "frumkin_theta": {"tier": "T2", "evidence": "Limits match Langmuir/Szyszkowski exactly; no published fitted K,a pair matched"},
    "frumkin_surface_tension": {"tier": "T2", "evidence": "Limits match Langmuir/Szyszkowski exactly; no published fitted K,a pair matched"},
    "frumkin_fit_K_and_a": {"tier": "T2", "evidence": "Limits match Langmuir/Szyszkowski exactly; no published fitted K,a pair matched"},
    "select_isotherm_model": {"tier": "T2", "evidence": "BIC logic on constructed data"},
    "hlb_griffin": {"tier": "T1", "evidence": "Ullmann's Table 13; −5…−7 % systematic offset explained (polydispersity)"},
    "hlb_davies": {"tier": "T1", "evidence": "Davies 1957 numbers; hard `KeyError` on unverified groups"},
    "hlb_from_groups": {"tier": "T1", "evidence": "Davies 1957 numbers; hard `KeyError` on unverified groups"},
    "hlb_davies_guo_ecl": {"tier": "T1", "evidence": "Primary PDF + patent-quoted coefficients"},
    "gemini_hlb_griffin_reference": {"tier": "T1", "evidence": "Curated literature tables **with per-entry confidence + citation in output** (good pattern — keep)"},
    "*_reference": {"tier": "T1", "evidence": "Curated literature tables **with per-entry confidence + citation in output** (good pattern — keep)"},
    "derive_davies_group_number*": {"tier": "T3", "evidence": "Self-consistency calibration by design, opt-in flag"},
    "DAVIES_DERIVED_HYDROPHILIC_GROUPS": {"tier": "T3", "evidence": "Self-consistency calibration by design, opt-in flag"},
    "tanford_*": {"tier": "T1", "evidence": "Bales 1998, Nagarajan Table 2, real SDS/DTAB paper ≤1.2 %, paper's own morphology calls"},
    "critical_packing_parameter": {"tier": "T1", "evidence": "Bales 1998, Nagarajan Table 2, real SDS/DTAB paper ≤1.2 %, paper's own morphology calls"},
    "classify_aggregate_morphology": {"tier": "T1", "evidence": "Bales 1998, Nagarajan Table 2, real SDS/DTAB paper ≤1.2 %, paper's own morphology calls"},
    "aggregation_number_spherical": {"tier": "T2", "evidence": "Bales worked-example self-consistency + range anchor (44.8–54.2); **still no independent (V,R)→N pair in literature**"},
    "nagarajan_equilibrium_area_ionic": {"tier": "T1", "evidence": "Nagarajan Table 2 within ~0.2 %"},
    "nagarajan_ruckenstein_*": {"tier": "T2", "evidence": "Hand-derivation + constants-completeness tests; no worked example of the components alone"},
    "blankschtein_*": {"tier": "T2", "evidence": "Hand derivations + cross-consistency with Nagarajan at β=0"},
    "estimate_axial_ratio_from_cpp_geometry": {"tier": "T2", "evidence": "Round-trip + regime guards; prolate-only (oblate explicitly out of scope)"},
    "cmc_from_surface_tension_curve": {"tier": "T1", "evidence": "AOT: 2.512 vs paper 2.51 mM; real researcher dataset regression; C20 verified vs AOT + primary source"},
    "cmc_from_conductivity_curve": {"tier": "T2", "evidence": "**No clean raw literature table exists** (PMC8007100 gives only summary CMC) — round-trip + SDS plausibility range + β cross-check only"},
    "aggregation_number_from_quenching_curve": {"tier": "T2", "evidence": "Equation verified from a real source image; **round-trip only** — notes explicitly call a real worked example a \"good next step\""},
    "aggregation_number_from_sls_debye_plot": {"tier": "T2", "evidence": "Formula verified vs 3 instrument docs; **round-trip only** (Brookhaven's own example uses instrument-internal counts, not re-derivable)"},
    "aggregation_number_from_dls": {"tier": "T2", "evidence": "Textbook form; realistic-range checks; hydration bias disclosed"},
    "aggregation_number_from_svedberg": {"tier": "T2", "evidence": "Textbook form; realistic-range checks"},
    "estimate_partial_specific_volume": {"tier": "T2", "evidence": "Physical-range checks; no external worked pair"},
    "hydrodynamic_radius_stokes_einstein": {"tier": "T1", "evidence": "Sutherland 2009 Taylor-dispersion D → 2.63 nm in range (+ documented DOSY negative control)"},
    "perrin_friction_factor": {"tier": "T2", "evidence": "Sphere limit exactly 1 + monotonicity; **no fetchable worked table** (two PDF attempts failed)"},
    "hydrodynamic_radius_perrin_corrected": {"tier": "T2", "evidence": "Sphere limit exactly 1 + monotonicity; **no fetchable worked table** (two PDF attempts failed)"},
    "ionic_strength": {"tier": "T1", "evidence": "Exact vs 2026 review (0.785 nm) + textbook (0.96 nm)"},
    "debye_length": {"tier": "T1", "evidence": "Exact vs 2026 review (0.785 nm) + textbook (0.96 nm)"},
    "grahame_equation_surface_potential": {"tier": "T2", "evidence": "Same functional form in secondary source + derivable low-potential limit; no blind numeric anchor"},
    "henry_function": {"tier": "T2", "evidence": "Source-verified asymptotics; Ohshima ~1 % vs exact — provenance solid, worked-example match not done"},
    "zeta_potential_henry": {"tier": "T1", "evidence": "Khademi exact round-trip + Serafini raw mobility <0.55 %"},
    "*_relaxation_corrected": {"tier": "T1†", "evidence": "Primary formula; eq-63 asymptote 1e-4 + Fig. 1/2 shape reproduction — not a blind table match"},
    "gibbs_free_energy_micellization": {"tier": "T1", "evidence": "Fu 2019: 0.03 %; Bales SDS worked example"},
    "mass_action_free_energy_of_micellization": {"tier": "T2", "evidence": "Rusanov primary PDF formula; convergence identity to pseudo-phase limit only (no published numeric example of this output reproduced)"},
    "critical_micellization_degree": {"tier": "T1", "evidence": "Rusanov worked example exact (0.061/0.091)"},
    "counterion_binding_degree_reference": {"tier": "T1", "evidence": "Bales EPR primary; confidence tiers returned"},
    "vant_hoff_enthalpy": {"tier": "T1†/T2", "evidence": "Sign bug found & fixed 2026-09-04; **magnitude = interval average only** — papers' polynomial-derivative ΔH will not match; magnitude figure honestly not re-issued post-fix"},
    "vant_hoff_multi_point_fit": {"tier": "T2", "evidence": "Kantonen Eq. 4 verified; round-trip + reference-invariance tests"},
    "entropy_of_micellization": {"tier": "dependent", "evidence": "Follows from inputs; observed failure mode = stochastic sign-substitution slips (documented, ~4.5 %)"},
    "hld_ionic": {"tier": "T1", "evidence": "Real reference tables cross-checked to papers' own reported pairs (tests assert \"all entries have real literature crosscheck\")"},
    "hld_nonionic": {"tier": "T1", "evidence": "Real reference tables cross-checked to papers' own reported pairs (tests assert \"all entries have real literature crosscheck\")"},
    "cc_mixing_rule": {"tier": "T1", "evidence": "Real reference tables cross-checked to papers' own reported pairs (tests assert \"all entries have real literature crosscheck\")"},
    "hld_optimal_salinity": {"tier": "T1", "evidence": "Real reference tables cross-checked to papers' own reported pairs (tests assert \"all entries have real literature crosscheck\")"},
    "fit_k_and_cc_from_salinity_scan": {"tier": "T2", "evidence": "Round-trip + manual-OLS only — **no multi-EACN raw dataset ever sourced** (disclosed since 2026-09-10)"},
    "molar_solubilization_ratio": {"tier": "T1", "evidence": "Paria & Yuet 2006 exact (0.073)"},
    "micelle_water_partition_coefficient": {"tier": "T2", "evidence": "Manual-formula + sanity; convention ambiguity disclosed (molar-ratio basis still unimplemented)"},
    "estimate_intrinsic_water_solubility_qspr": {"tier": "T3", "evidence": "ESOL/Pat Walters RDKit refit; 0.6–1 log-unit error stated in tool docstring"},
    "work_of_adhesion": {"tier": "T1", "evidence": "6 real systems ≤0.34 %"},
    "spreading_coefficient": {"tier": "T1", "evidence": "6 real systems ≤0.34 %"},
    "owens_wendt_*": {"tier": "T2", "evidence": "Formulas source-verified pre-build; round-trip + zero-polar check only — **no published worked OWRK/vOCG table matched**"},
    "van_oss_chaudhury_good_*": {"tier": "T2", "evidence": "Formulas source-verified pre-build; round-trip + zero-polar check only — **no published worked OWRK/vOCG table matched**"},
    "capillary_number": {"tier": "T2", "evidence": "Textbook dimensionless groups (trivial risk)"},
    "eor_capillary_number": {"tier": "T2", "evidence": "Textbook dimensionless groups (trivial risk)"},
    "classify_surfactant_charge_type*": {"tier": "T2", "evidence": "Rule/SMARTS tests only — **no external labeled-SMILES corpus** validates the rules as a whole; refusal-first design documented"},
    "classify_surfactant_structural_family": {"tier": "T2", "evidence": "Rule/SMARTS tests only — **no external labeled-SMILES corpus** validates the rules as a whole; refusal-first design documented"},
    "derive_all_properties_from_smiles_and_curve": {"tier": "chained", "evidence": "Inherits tiers of its callees; `gaps` list is a good pattern"},
}

TIER_DEFINITIONS: Final[dict[str, str]] = {
    "T1": "Primary-source numeric anchor: a real paper's own reported number is reproduced (or the chain is digit-level matched). *Safe for \"match-the-literature\" gold claims.*",
    "T2": "Source-verified formula, round-trip/proof-only output: the equation was verified against a real source before implementation, but **no independent published numeric example of the output exists** — validation is round-trip, limit identities, hand derivations, or plausibility ranges. *The \"difficult to generate literature-backed outputs\" set the audit asked for. Safe for computation golds; unsafe for \"the literature says\" claims.*",
    "T3": "Derived/estimate by construction: explicitly non-literature (self-consistency calibration or statistical model), disclosed at the call site.",
    "T4": "Documented dead end: the honest output is \"cannot compute\" — see §1.3.",
    "T1\u2020": "strong-but-qualified T1 -- the audit's dagger: primary anchor holds for the main quantity; see the row's evidence cell for the qualification (e.g. method-level rather than table-level match)",
    "T1\u2020/T2": "split qualification -- qualified T1 for the primary case, T2 for the other variant named in the row (vant_hoff_enthalpy)",
    "dependent": "follows from its inputs; validated only as an internal-consistency delta (no independent anchor)",
    "chained": "inherits the tiers of its callee tools (orchestrator output)",
}

#: MCP tool names (and library names cited in benchmark gold provenance)
#: that wrap a technique registered above under a different name.
TOOL_ALIASES: Final[dict[str, str]] = {
    "gibbs_surface_excess": "gibbs_gamma_max",
    "gibbs_area_per_molecule": "gibbs_a_min",
    "szyszkowski_predict_surface_tension": "szyszkowski_surface_tension",
    "frumkin_predict_surface_tension": "frumkin_surface_tension",
    "hlb_from_mw": "hlb_griffin",
    "hld_fit_k_and_cc_from_salinity_scan": "fit_k_and_cc_from_salinity_scan",
    "aggregation_number": "aggregation_number_spherical",
    "aggregation_number_from_quenching": "aggregation_number_from_quenching_curve",
    "aggregation_number_from_sls": "aggregation_number_from_sls_debye_plot",
    "cmc_from_conductivity": "cmc_from_conductivity_curve",
    "debye_screening_length": "debye_length",
    "zeta_potential": "zeta_potential_henry",
    "hydrodynamic_radius": "hydrodynamic_radius_stokes_einstein",
    "wetting_work_of_adhesion": "work_of_adhesion",
    "wetting_spreading_coefficient": "spreading_coefficient",
    "rubingh_activity_coefficients": "rubingh_solve",
    "activity_coefficients": "rubingh_solve",
    "solve_rubingh_x": "rubingh_solve",
    "asymmetric_margules_activity_coefficients": "asymmetric_margules",
    "entropy_micellization": "entropy_of_micellization",
    "optimal_salinity_ionic": "hld_optimal_salinity",
    "rodenas_x1_series": "rodenas_x1_from_series",
}

_T1D = "T1†"

#: Techniques shipped after the 2026-09-30 audit table, with the evidence
#: their own test files establish (same tier definitions as above).
POST_AUDIT_TIERS: Final[dict[str, dict[str, str]]] = {
    "eommm_global_fit_all_optima": {"tier": "T2", "evidence": "Returns eommm_global_fit's optimum on a known unique case; basin search finds several optima at a loose tolerance (tests/test_mixed_micelle.py)"},
    "counterion_binding_degree": {"tier": "T2", "evidence": "Slope-ratio arithmetic beta = 1 - S_above/S_below with input guards; the method's physical interpretation is questioned for some systems, as the tool states"},
    "hlb_method_recommendation": {"tier": "T2", "evidence": "Rule-based choice of HLB method citing Liao et al. 2023 for gemini/glycolipid families; every family, caveat and refusal tested (tests/test_hlb_cpp.py); no numeric output"},
    "get_standard_probe_liquid_properties": {"tier": _T1D, "evidence": "Standard tabulated dispersive/polar components (e.g. water 21.8/51.0 mN/m) returned verbatim; untabulated liquids refused"},
    "simulate_emulsion_stability": {"tier": "T2", "evidence": "Stokes velocity and creaming time equal the closed form v = drho g d^2/(18 eta) to 1e-12; monodisperse t50 = H/(2v) and the separated fraction with coalescence equal hand-integrated closed forms (tests/test_emulsion.py)"},
    "population_balance_emulsion": {"tier": "T2", "evidence": "Constant-kernel number decay equals the exact Smoluchowski solution; mass conserved for sum/product kernels and breakage; LSW ripening law exact (tests/test_emulsion.py)"},
    "surfactant_transport_1d": {"tier": "T2", "evidence": "Tracer breakthrough matches the Ogata-Banks analytic solution; retarded case matches the constant-R substitution (tests/test_transport.py)"},
    "analyze_micelle_trajectory": {"tier": "T2", "evidence": "Cluster census exact on seeded synthetic trajectories with known aggregates; analysis only (tests/test_trajectory.py)"},
    "conformer_volume_3d": {"tier": "T3", "evidence": "RDKit ETKDG + MMFF conformer estimate; anisotropy exact for collinear and tetrahedral points, volume linear in chain length; no literature volume matched"},
    "cloud_point_estimate": {"tier": _T1D, "evidence": "HLD closed form recomputed exactly; C12E5 inside the source's stated accuracy window (Zarate-Munoz 2015, about 29 C, HLD within 5 C); NaCl lowers the cloud point (tests/test_cloudpoint.py)"},
    "phase_inversion_temperature": {"tier": "T2", "evidence": "Algebraic inverse of hld_nonionic/hld_ionic; round-trips to HLD = 0 for both families; refuses guessed inputs (tests/test_hld.py)"},
    "rdkit_descriptor_panel": {"tier": _T1D, "evidence": "RDKit descriptors; aspirin TPSA 63.6 A^2 (Ertl definition) exact; Crippen logP within the method's stated accuracy for four compounds (tests/test_descriptors.py)"},
    "bootstrap_cmc": {"tier": "T2", "evidence": "Synthetic sigmoid with known CMC: point estimate within 5 %, percentile interval covers the truth, seed-deterministic (tests/test_bootstrap_cmc.py)"},
    "planar_pb_profile": {"tier": "T2", "evidence": "Exact 1:1 Gouy-Chapman profile (max error below 1e-6 mV) and Grahame surface charge; Debye-Hueckel limit at low potential (tests/test_poisson_boltzmann.py)"},
    "mc_aggregation_distribution": {"tier": "T2", "evidence": "Monte Carlo marginals match the exact Boltzmann sum over integer partitions within sampling error; seed-reproducible (tests/test_mc_aggregation.py)"},
    "dls_size_distribution": {"tier": "T2", "evidence": "Scattering vector matches the hand value (633 nm, 90 deg, water); monodisperse and bimodal synthetic correlations recover the input sizes and mass split (tests/test_dls_inversion.py)"},
    "ward_tordai_adsorption": {"tier": "T2", "evidence": "Fixed-sink case equals Gamma = 2 (c0 - cs) sqrt(D t / pi) at all times; dilute Langmuir case equals the exact linear-isotherm solution K_H c0 [1 - exp(t/tau) erfc(sqrt(t/tau))] to 2e-4; transient obeys the short-time law and approaches equilibrium from below (tests/test_ward_tordai.py)"},
    "nagarajan_cmc_predictor": {"tier": "T1", "evidence": "Seven experimental CMCs, the source's own validation class (Nagarajan 1985), within 0.15 log10 units; Klevens slope of the sulfate series (tests/test_nagarajan_cmc.py)"},
    "nagarajan_free_energy_terms": {"tier": "T2", "evidence": "Term-table identities recomputed independently; sphere/cylinder continuity; Debye kappa equals the source's closed form (tests/test_nagarajan_cmc.py)"},
    "mass_action_equilibrium": {"tier": "T2", "evidence": "Closed-form mass-action distribution; mass balance to 1e-9; grouped form equals the naive form where both are finite (tests/test_kinetics_mass_action.py)"},
    "micelle_relaxation_odes": {"tier": "T2", "evidence": "Steady state agrees with mass_action_equilibrium; relaxation times consistent in order of magnitude with kmc_micelle_kinetics (tests/test_kinetics.py)"},
    "kmc_micelle_kinetics": {"tier": "T2", "evidence": "Stochastic simulation of the same stepwise scheme; equilibrium consistent with mass_action_equilibrium; seed-reproducible (tests/test_kinetics.py)"},
    "aniansson_wall_two_regime": {"tier": "T2", "evidence": "Closed-form inversion of the Aniansson-Wall equations (Patist et al. 2002, eqs 3-4); round trip reproduces tau1 and tau2 exactly; hand-derived SDS case; implied k_on diffusion-limited (tests/test_kinetics.py)"},
    "charge_regulation_pb": {"tier": "T2", "evidence": "Weak-coupling limit recovers bare Henderson-Hasselbalch; strong coupling self-limits; site and double-layer charges balance (tests/test_charge_regulation.py)"},
    "obrien_white_numerical": {"tier": "T2", "evidence": "Hueckel limit 2/3 at small ka and monotone approach to Smoluchowski; ordering against the OHW83 approximation above ka = 10 (tests/test_obrien_white_e1.py)"},
    "dynamic_surface_tension_forward": {"tier": "T2", "evidence": "Long-time limit equals the Szyszkowski equation; kinetic limit equals the analytic Langmuir uptake; converged to 2e-3 mN/m against a 4x finer grid and step; surface excess within 2e-4 of the independent Ward-Tordai integral under diffusion control (tests/test_dynamic_tension_o3.py, tests/test_cross_solver_adsorption.py); no raw literature transient matched"},
    "dynamic_surface_tension_ode": {"tier": "T2", "evidence": "Joint fit recovers the identifiable lumped quantities (Gamma_eq, k_off(1 + Kc)) from synthetic data and flags the K-Gamma_max ridge (tests/test_dynamic_tension_o3.py)"},
    "eacn_group_contribution": {"tier": _T1D, "evidence": "n-alkanes EACN = N_C (Wade 1977, definitional); octylbenzene 8 (Queste et al. 2007); refuses branched, ring, heteroatom and multi-fragment structures (tests/test_eacn_c3.py)"},
    "eacn_linear_mixing": {"tier": "T2", "evidence": "Mole-fraction linear mixing rule; midpoint identity and monotone HLD feed (tests/test_eacn_c3.py)"},
    "brownian_self_assembly": {"tier": "T3", "evidence": "Reduced-unit Langevin model; qualitative CMC-like trend only (micellizes at 0.15 M, monomeric at 0.002 M); seed-deterministic (tests/test_selfassembly_s3.py)"},
    "brownian_cmc_scan": {"tier": "T3", "evidence": "Concentration scan of brownian_self_assembly; qualitative trend, not an absolute CMC"},
    "cmc_from_free_monomer": {"tier": "T2", "evidence": "Rummel et al. 2026 eqs 3-4: beta = 0 identity exact; forward/inverse round trip and independent hand evaluation of the published equation (tests/test_dpd_cmc_s4.py)"},
    "free_monomer_from_cmc": {"tier": "T2", "evidence": "Inverse of cmc_from_free_monomer; round trip exact (tests/test_dpd_cmc_s4.py)"},
    "anderson_cmc_analysis": {"tier": "T2", "evidence": "Eq-8 profile fit recovers an exactly generated trace to machine precision; jackknife sigma equals the closed-form leave-one-out value (tests/test_dpd_cmc_s4.py)"},
    "dpd_micelle_cmc": {"tier": "T3", "evidence": "In-house DPD (C8E4 model): force invariants and model constants verified; finite-size boxes, so the CMC is a qualitative estimate, as the tool states (tests/test_dpd_micelle_s4.py)"},
    "submit_dpd_micelle_cmc_job": {"tier": "T3", "evidence": "Asynchronous form of dpd_micelle_cmc"},
    "dpd_micelle_cmc_job_status": {"tier": "T3", "evidence": "Status/result lookup for submit_dpd_micelle_cmc_job"},
    "nonlinear_fit_with_errors": {"tier": "T2", "evidence": "Linear case equals closed-form OLS values and standard errors; noise-free exponential decay recovered exactly (tests/test_inference.py)"},
    "propagate_uncertainty_chain": {"tier": "T2", "evidence": "Debye-length chain matches the delta-method sd (I^-1/2 scaling) within 10 %; seed-deterministic (tests/test_inference.py)"},
    "mcmc_posterior_fit": {"tier": "T2", "evidence": "Flat-prior known-sigma linear case reproduces the OLS Gaussian (median, sd within 20 %), R-hat below 1.1; seed-deterministic (tests/test_inference.py)"},
    "global_multi_dataset_fit": {"tier": "T2", "evidence": "Noise-free shared slope and per-dataset intercepts recovered exactly (tests/test_inference.py)"},
    "cmc_homologous_series": {"tier": _T1D, "evidence": "Klevens slope 0.30 per CH2 from sodium alkyl sulfate CMCs (Mukerjee and Mysels 1971); exact recovery and extrapolation flag on a synthetic series (tests/test_inference.py)"},
}


def tools_in_tier(*tiers: str) -> list[str]:
    """Return the registry keys whose tier is in `tiers` (e.g. "T1", "T2").
    Wildcard keys are included when their tier matches, so callers get the
    literal audit-table keys the benchmark generator filters by.
    Unknown tier names raise ValueError (fail-fast, field names included).
    """
    unknown = [t for t in tiers if t not in TIER_DEFINITIONS]
    if unknown:
        raise ValueError(
            f"unknown tier(s) {unknown!r} -- valid tiers: {sorted(TIER_DEFINITIONS)}"
        )
    return [name for name, rec in VALIDATION_TIERS.items() if rec["tier"] in tiers]


def get_tool_tier(name: str) -> str | None:
    """Return the validation tier for a tool name: exact audit-table match,
    then an MCP alias of an audit-table technique, then the post-audit
    table, then audit-table wildcards such as "tanford_*".  Returns None
    when the name is not registered at all.
    """
    if name in VALIDATION_TIERS:
        return VALIDATION_TIERS[name]["tier"]
    if name in TOOL_ALIASES:
        return VALIDATION_TIERS[TOOL_ALIASES[name]]["tier"]
    if name in POST_AUDIT_TIERS:
        return POST_AUDIT_TIERS[name]["tier"]
    for pattern, rec in VALIDATION_TIERS.items():
        if "*" in pattern and fnmatchcase(name, pattern):
            return rec["tier"]
    return None
