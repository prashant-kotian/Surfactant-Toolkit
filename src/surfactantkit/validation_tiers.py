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
    "mass_action_free_energy_of_micellization": {"tier": "T1", "evidence": "Rusanov primary PDF; convergence identity to pseudo-phase limit"},
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
    """Return the validation tier for a tool name (exact match first, then
    audit-table wildcards such as "tanford_*").  Returns None when the name
    is not in the registry at all.
    """
    if name in VALIDATION_TIERS:
        return VALIDATION_TIERS[name]["tier"]
    for pattern, rec in VALIDATION_TIERS.items():
        if "*" in pattern and fnmatchcase(name, pattern):
            return rec["tier"]
    return None
