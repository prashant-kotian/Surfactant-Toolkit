# Surfactant Toolkit — validation-tier audit

> **Context.** Point-in-time audit record (2026-09-30), reproduced from the
> project's development repository as provenance for the validation-tier
> registry. Findings cite the codebase state at that date (counts, file names,
> local paths); all P0/P1 process items are resolved in this artifact. The
> tier table in §1.1 is machine-checked against the code registry by
> `tests/test_validation_tiers.py`; edit it only together with the registry.
**Date:** 2026-09-30 · **Scope:** `src/surfactantkit/`, repository tests and benchmark data only — no code, config, or doc files were modified; this file is the only addition.

---

## 0. Baseline verified today

| Fact | Value | How verified |
|---|---|---|
| MCP tools | **81** (`@mcp.tool` count) | `grep -c "@mcp.tool" src/surfactantkit/mcp_server.py` |
| Library modules | 14 theory modules + `mcp_server.py` (2,044 lines) | directory listing |
| Public library functions | 93 (excluding private) | AST scan |
| Tests | **522 passed, 0 failed, 121 s** | `PYTHONPATH=src python3 -m pytest tests/ -q` (sandbox needed `pip install rdkit mcp` first — see §C.5) |
| Version | 0.4.0, alpha, deps: `rdkit` only (`mcp` optional extra) | `pyproject.toml` |
| Last commit | 2026-09-21; `benchmark/hard50/` and `benchmark/extrahard50/` are **untracked** | `git status` |

---

## 1. Part A — Literature-backing audit: where outputs are hard to anchor

The benchmark's gold answers come from this toolkit, so every technique's output needs a clear answer to: *could this number ever be checked against a published value?* Three failure shapes exist:

- **A-1: no published worked example exists** for the technique's output (only the *equation* is source-verified; validation is round-trip/self-consistency);
- **A-2: published data exists but is structurally unusable** (image tables, mismatched conventions, no raw table);
- **A-3: the output is deliberately NOT literature-derived** (by-design estimate).

### Tier definitions
- **T1 — Primary-source numeric anchor:** a real paper's own reported number is reproduced (or the chain is digit-level matched). *Safe for "match-the-literature" gold claims.*
- **T2 — Source-verified formula, round-trip/proof-only output:** the equation was verified against a real source before implementation, but **no independent published numeric example of the output exists** — validation is round-trip, limit identities, hand derivations, or plausibility ranges. *The "difficult to generate literature-backed outputs" set the audit asked for. Safe for computation golds; unsafe for "the literature says" claims.*
- **T3 — Derived/estimate by construction:** explicitly non-literature (self-consistency calibration or statistical model), disclosed at the call site.
- **T4 — Documented dead end:** the honest output is "cannot compute" — see §1.3.

### 1.1 Full technique → tier table

| Technique (tool/function) | Tier | Evidence / why |
|---|---|---|
| `clint_ideal_cmc` | T1 | Exact on 6+ published systems (PMC6554738, PMC4087020, …) |
| `rubingh_solve` (x1) | T1 | ~exact across 6 independent systems |
| `rubingh_beta` (pointwise) | T1† | Exact for fixed 1:1 (cholate-SDS); 5–20 % gap vs *regression* betas documented as method difference, not bug |
| `rubingh_beta_regression` | **T2** | Round-trip construction only; no paper publishes the same pointwise-mean procedure with raw data |
| `motomura_ideal_composition` | **T2** | Algebraic-identity tests only |
| `rodenas_x1` (pointwise) | T1† | Exact zero-intercept identity test; real-data series check exists |
| `rodenas_x1_from_series` | **T2** | Real data tested, but endpoints err **19–60 %** (paper uses a different global-fit differentiation) — documented limit, sparse-series output hard to anchor |
| `asymmetric_margules`, `excess_free_energy`, `rodenas_activity_coefficients` | T1 | Schulz & Durand SI CSVs digit-level matched (mmc2/3/6/8/10) |
| `eommm_global_fit` | T1† | Formula confirmed exactly vs SI-fitted curves; real Hyamine/DTAB convergence; margin-dependent uniqueness disclosed |
| `eommm_find_minimal_feasible_margin` | **T2** | Round-trip + real-data convergence; the tightest-margin procedure itself has no published number to match |
| `select_mixture_model` | **T2** | BIC on synthetic constructions; scope choices (why Rosen/Maeda/Rodenas sit out) documented in the tool |
| `maeda_free_energy_of_micellization` | T1 | Rub et al. 2021 (~6 %, attributed to paper rounding); component-indexing ambiguity honestly disclosed |
| `rosen_monolayer_solve` / `rosen_beta_sigma` | T2 | Verified mathematically identical to Rubingh (by design — a meaning distinction, not an independent anchor) |
| `corrin_harkins_predict` | **T2** | Round-trip + salt monotonicity + input guards only |
| `gibbs_gamma_max` | T1 | 4 real ionic systems ≤ 0.5 % |
| `gibbs_a_min` | T1 | 5 real systems ≤ 0.4 % |
| `szyszkowski_surface_tension` | T1 | C12E8 (Ho/Lips/Velikov 2023) curve; Γ-round-trip 0.02 % |
| `szyszkowski_fit_K` | **T2** | Round-trip only — the one real candidate K table (Hreczuch 2024) is provably unusable (OCR/exponent contradiction, reconfirmed twice) |
| `frumkin_theta`/`frumkin_surface_tension`/`frumkin_fit_K_and_a` | **T2** | Limits match Langmuir/Szyszkowski exactly; no published fitted K,a pair matched |
| `select_isotherm_model` | **T2** | BIC logic on constructed data |
| `hlb_griffin` | T1 | Ullmann's Table 13; −5…−7 % systematic offset explained (polydispersity) |
| `hlb_davies`, `hlb_from_groups` | T1 | Davies 1957 numbers; hard `KeyError` on unverified groups |
| `hlb_davies_guo_ecl` | T1 | Primary PDF + patent-quoted coefficients |
| `gemini_hlb_griffin_reference`, `*_reference` tables (β, dn/dc, HLD Cc/k, probe liquids, water solubility) | T1 | Curated literature tables **with per-entry confidence + citation in output** (good pattern — keep) |
| `derive_davies_group_number*`, `DAVIES_DERIVED_HYDROPHILIC_GROUPS` | **T3** | Self-consistency calibration by design, opt-in flag |
| `tanford_*`, `critical_packing_parameter`, `classify_aggregate_morphology` | T1 | Bales 1998, Nagarajan Table 2, real SDS/DTAB paper ≤1.2 %, paper's own morphology calls |
| `aggregation_number_spherical` | **T2** | Bales worked-example self-consistency + range anchor (44.8–54.2); **still no independent (V,R)→N pair in literature** |
| `nagarajan_equilibrium_area_ionic` | T1 | Nagarajan Table 2 within ~0.2 % |
| `nagarajan_ruckenstein_*` free-energy components | **T2** | Hand-derivation + constants-completeness tests; no worked example of the components alone |
| `blankschtein_*` components | **T2** | Hand derivations + cross-consistency with Nagarajan at β=0 |
| `estimate_axial_ratio_from_cpp_geometry` | **T2** | Round-trip + regime guards; prolate-only (oblate explicitly out of scope) |
| `cmc_from_surface_tension_curve` (+C20/pC20) | T1 | AOT: 2.512 vs paper 2.51 mM; real researcher dataset regression; C20 verified vs AOT + primary source |
| `cmc_from_conductivity_curve` | **T2** | **No clean raw literature table exists** (PMC8007100 gives only summary CMC) — round-trip + SDS plausibility range + β cross-check only |
| `aggregation_number_from_quenching_curve` (SSFQ) | **T2** | Equation verified from a real source image; **round-trip only** — notes explicitly call a real worked example a "good next step" |
| `aggregation_number_from_sls_debye_plot` | **T2** | Formula verified vs 3 instrument docs; **round-trip only** (Brookhaven's own example uses instrument-internal counts, not re-derivable) |
| `aggregation_number_from_dls` | **T2** | Textbook form; realistic-range checks; hydration bias disclosed |
| `aggregation_number_from_svedberg` | **T2** | Textbook form; realistic-range checks |
| `estimate_partial_specific_volume` | **T2** | Physical-range checks; no external worked pair |
| `hydrodynamic_radius_stokes_einstein` | T1 | Sutherland 2009 Taylor-dispersion D → 2.63 nm in range (+ documented DOSY negative control) |
| `perrin_friction_factor`, `hydrodynamic_radius_perrin_corrected` | **T2** | Sphere limit exactly 1 + monotonicity; **no fetchable worked table** (two PDF attempts failed) |
| `ionic_strength`, `debye_length` | T1 | Exact vs 2026 review (0.785 nm) + textbook (0.96 nm) |
| `grahame_equation_surface_potential` | **T2** | Same functional form in secondary source + derivable low-potential limit; no blind numeric anchor |
| `henry_function` (huckel/smoluchowski/ohshima/swan_furst) | **T2** | Source-verified asymptotics; Ohshima ~1 % vs exact — provenance solid, worked-example match not done |
| `zeta_potential_henry` | T1 | Khademi exact round-trip + Serafini raw mobility <0.55 % |
| `*_relaxation_corrected` (OHW83) | T1† | Primary formula; eq-63 asymptote 1e-4 + Fig. 1/2 shape reproduction — not a blind table match |
| `gibbs_free_energy_micellization` | T1 | Fu 2019: 0.03 %; Bales SDS worked example |
| `mass_action_free_energy_of_micellization` | T2 | Rusanov primary PDF formula; convergence identity to pseudo-phase limit only (no published numeric example of this output reproduced) |
| `critical_micellization_degree` | T1 | Rusanov worked example exact (0.061/0.091) |
| `counterion_binding_degree_reference` | T1 | Bales EPR primary; confidence tiers returned |
| `vant_hoff_enthalpy` | T1†/**T2** | Sign bug found & fixed 2026-09-04; **magnitude = interval average only** — papers' polynomial-derivative ΔH will not match; magnitude figure honestly not re-issued post-fix |
| `vant_hoff_multi_point_fit` | **T2** | Kantonen Eq. 4 verified; round-trip + reference-invariance tests |
| `entropy_of_micellization` | dependent | Follows from inputs; observed failure mode = stochastic sign-substitution slips (documented, ~4.5 %) |
| `hld_ionic/nonionic`, `cc_mixing_rule`, `hld_optimal_salinity` | T1 | Real reference tables cross-checked to papers' own reported pairs (tests assert "all entries have real literature crosscheck") |
| `fit_k_and_cc_from_salinity_scan` | **T2** | Round-trip + manual-OLS only — **no multi-EACN raw dataset ever sourced** (disclosed since 2026-09-10) |
| `molar_solubilization_ratio` | T1 | Paria & Yuet 2006 exact (0.073) |
| `micelle_water_partition_coefficient` | **T2** | Manual-formula + sanity; convention ambiguity disclosed (molar-ratio basis still unimplemented) |
| `estimate_intrinsic_water_solubility_qspr` | **T3** | ESOL/Pat Walters RDKit refit; 0.6–1 log-unit error stated in tool docstring |
| `work_of_adhesion`, `spreading_coefficient` | T1 | 6 real systems ≤0.34 % |
| `owens_wendt_*`, `van_oss_chaudhury_good_*` | **T2** | Formulas source-verified pre-build; round-trip + zero-polar check only — **no published worked OWRK/vOCG table matched** |
| `capillary_number` / `eor_capillary_number` | **T2** | Textbook dimensionless groups (trivial risk) |
| `classify_surfactant_charge_type*`, `classify_surfactant_structural_family` | **T2** | Rule/SMARTS tests only — **no external labeled-SMILES corpus** validates the rules as a whole; refusal-first design documented |
| `derive_all_properties_from_smiles_and_curve` | chained | Inherits tiers of its callees; `gaps` list is a good pattern |

† = strong but qualified (see cell). **T2 count: ~28 of ~60 distinct techniques — this is the "no literature backing" set.**

### 1.2 Published data that exists but is unusable (A-2)
Confirmed dead-end data pulls, each with evidence — do not re-attempt casually:
- **Hreczuch 2024 Szyszkowski K table** — header/cell exponent contradiction, unusable (reconfirmed against raw XML).
- **Muherei & Junin 2009 Table 2A** — α→CMCid mapping is an embedded image; PDF redirects JS-gated/403.
- **mmc4.csv Δθ_cmc** — symbol defined only inside a figure; zero mentions in extractable text.
- **Rubingh β as published** — pointwise vs regression determination (5–20 % structural gap; the one fixed 1:1 system matches exactly).
- **HLD multi-EACN scans, conductivity raw tables, SSFQ/SLS worked examples** — searched repeatedly, never found clean open raw tables.

### 1.3 T4 — documented dead ends (honest output is "cannot compute")
1. Gouy-Chapman-corrected `n_factor` for Gibbs (self-consistent closure term unresolvable; tool **removed** rather than shipped wrong).
2. FLM reorientation/multistate isotherm (same closure class; deferred).
3. Viscosity → Nagg route (insensitive to N for compact spheres — circular; correctly rejected).
4. Davies ionic group numbers (no trustworthy literature revision exists; only the disclosed derived-opt-in path).
5. Km on molar-concentration basis; Oremusova-style `ln(CMC)` ΔG convention (conventions, not bugs — documented as unimplemented).
6. Schulz & Durand free-energy objective as literally transcribed (unbounded below; replaced by margin-tightening).
7. Full Blankschtein MT theory, QSPR/ML HLB (out of scope by explicit decision).

> **Benchmark implication:** gold questions should draw *literature-claim* stems from T1 techniques only; T2 is fine for pure computation but must never say "as reported in the literature"; T4 is exactly the refusal-trap material the EXTRA50 G7/G9 families already exploit.

---

## 2. Part B — New techniques that can ride on EXISTING tools (no new workers)

Ranked by value/risk. **Invariant for all:** defaults and existing outputs must not change — every gold bank in `paper3_*` pins current values.

| # | Addition | Rides on | What it adds | Source basis | Risk |
|---|---|---|---|---|---|
| B1 | **Alternative break-point estimators** (`method="rss"|"derivative"|"tangential"`) | `cmc_from_surface_tension_curve` | The other standard literature CMC extraction methods on the same raw curve — turns one technique into three for the same Excel template | Rosen (graphical tangent method); derivative-max is standard practice | Low if default `rss` preserved; cross-method disagreement becomes new benchmark material |
| B2 | **Conductometric CMC alternatives** (`method="rss"|"jacobs_jones"|"onsager"`) | `cmc_from_conductivity` | Jacobs–Jones second-derivative inflection and Onsager limiting-slope β are the two classic alternatives to segmented regression | Standard conductometry practice | Low (additive) |
| B3 | **PIT output** (phase-inversion temperature from the same HLD=0 equation, varying T) | `hld_optimal_salinity` / `hld_ionic` | Shinoda PIT falls straight out of the equation already implemented — a genuinely missing standard quantity | Shinoda/HLD-NAC standard | Low; needs one verified sign convention for α |
| B4 | **Validity flag for short/medium chains** | `critical_packing_parameter` | The 2019/2020 ACIS review says the space-filled-spherical-core model *cannot* apply to short-chain surfactants — an advisory output field, not new math | PMID 31954873 (already cited in METHOD_ALTERNATIVES review) | Very low; pure disclosure |
| B5 | **Third candidate in the isotherm BIC race** (e.g. Flory–Huggins/Temkin) | `select_isotherm_model` | Selector currently Langmuir-vs-Frumkin only; a third literature isotherm would first need its forward formula in `adsorption.py` (module-level, tool unchanged) | Real surfactant-adsorption alternatives per review | Medium (new formula needs T1/T2 sourcing first) |
| B6 | **Motomura as a BIC candidate** | `select_mixture_model` | Real papers already run Motomura/Rodenas/Maeda side-by-side; Motomura is 0-parameter like Clint → same RSS space | Review's own headline finding | Medium — must verify Motomura yields a CMC-mix prediction usable in the BIC race before committing |
| B7 | **Oblate-side axial ratio** (vesicle/bilayer regime) | `estimate_axial_ratio_from_cpp_geometry` | Currently rejects everything outside the rodlike regime; the oblate branch has the same Perrin math already shipped | Needs a sourced oblate-core construction (Israelachvili/Nagarajan) | Medium — do not guess the geometry |
| B8 | **Optional ΔCp term** in entropy/ΔG(T) reporting | `entropy_of_micellization` | Uses ΔCp that `vant_hoff_multi_point_fit` already returns — new optional input, same tool | Kantonen 2018 Eq. 4 already in library | Low |
| B9 | **Grow `DAVIES_DERIVED_HYDROPHILIC_GROUPS`** (more families as calibrations are done) | `hlb_from_groups` (`allow_derived_groups=True`) | Pure table growth behind the existing opt-in flag; imidazoline's structural gap is a template for when to *not* add | Same Griffin cross-calibration method | Low; keep derived/literature tables separate |
| B10 | **Extend confidence-tiered reference tables** (β, Cc/k, dn/dc, probe liquids) | existing `*_reference` tools | Tables are the most literature-anchored outputs in the kit; adding entries = adding T1 coverage | Per-entry sourcing as already practiced | Low |

**Explicitly out of scope per the constraint (would need a NEW worker):** standalone tools for `frumkin_fit_K_and_a` and `select_isotherm_model` — both are exported, tested, and used by the orchestrator, but unlike `szyszkowski_fit_K` they have **no MCP tool**; an MCP-only user can *predict* Frumkin but cannot *fit* it or run the BIC choice (except via the orchestrator). Recorded as a coverage asymmetry, not acted on.

**Do not re-attempt (T4):** Gouy-Chapman n_factor on `gibbs_gamma_max`, FLM isotherm, viscosity→Nagg, ionic mass-action extension until a verified counterion-binding equilibrium source is in hand.

---

## 3. Part C — Overall improvements (prioritized)

### P0 — real, verified today
1. **MCP error shaping.** `mcp_server.py` contains **zero `except ValueError` / `except Exception` handlers** across 81 tools — library exceptions surface as raw MCP internal errors. The project's own Gibbs `n_factor` experiment proved the mechanism: structured, readable "cannot compute + reason" responses moved SurfMCP 0/6 → 6/6. A thin decorator converting `ValueError`/`KeyError` into a consistent `{"status": "no_solution", "reason": …}` payload would generalize that win to *every* raising tool (`hlb_davies`, `gibbs_gamma_max` system_type, regime guards, …).
2. **Tool-description bloat vs MCPGAUGE.** Descriptions: **mean 866 chars, median 864, max 2,362 (`eommm_global_fit`), 17 tools > 1,200 chars, 48/81 > 800.** MCPGAUGE (arXiv 2508.12566) found naive MCP access *lowered* accuracy 9.5 % avg, attributed to description bloat/selection. Recommendation: trim the MCP-facing docstrings to short summaries + argument essentials (full reasoning stays in library docstrings), and keep SurfMCP Version A/B descriptions crisp — this is the single highest-leverage change for the augmented condition.
3. **README status line is stale.** Says "67 MCP tools, 438 tests … 2026-09-15"; actual is **81 tools / 522 tests**. It also still reports MSR as "no confirmed real-paper match yet — still open", contradicted by `literature_validation_notes.md` Round 5 (Paria & Yuet 2006, exact) and an existing test asserting 0.073. The README's own "re-verify rather than trust" warning was not applied to itself.

### P1 — structural
4. **Machine-readable validation-tier registry.** No file maps technique → tier → evidence (only prose across two long docs). §1.1 should become a small `VALIDATION_TIERS` dict (tool → tier → citation/test name) so the benchmark generator can filter gold stems by tier mechanically.
5. **Fresh-environment onboarding.** Bare `pytest tests/` fails with 17 collection errors (`ModuleNotFoundError: surfactantkit`); needs `PYTHONPATH=src` or documented `pip install -e ".[dev]"` plus `rdkit`/`mcp`. Cost an audit step today; will cost any new collaborator the same.
6. **Doc portability.** `ROADMAP.md` (stale since 2026-09-13, self-declared) points the reader to `H:\CodeProjects\PhD-Research\PROJECT_STATUS.md` — a path on another machine. Current-state truth is split across a Windows-only file, ROADMAP, README, and `literature_validation_notes.md`.
7. **Tool-name ↔ library-name parity.** Names drift: `hydrodynamic_radius` → `hydrodynamic_radius_stokes_einstein`, `entropy_of_micellization` → `entropy_micellization`, `gibbs_surface_excess` → `gibbs_gamma_max`, `aggregation_number_from_quenching` → `…_curve`. This cost real debugging time during EXTRA50 construction. A one-page alias table (or a parity test) would prevent repeat confusion for both humans and model agents.
8. **Commit the benchmark banks.** `benchmark/hard50/` and `paper3_extra50/` are untracked; GroundZero is committed. Gold JSONs should be under version control before any run mutates adjacent files.

### P2 — hygiene / nice-to-have
9. **Duplicated private helpers:** `_bic` (curve_analysis + mixture_model_selection), `_solve_3x3` (thermodynamics + wetting), `_golden_section_minimize` (adsorption + mixed_micelle) — candidates for one internal `_optimize`/`_stats` module (no behavior change).
10. **Uniform runtime disclosures:** EOMMM (seconds→minutes) and Frumkin fit (~2.8 s) document cost ad hoc; make "expected runtime" a standard docstring field so MCP clients can set timeouts.
11. **`mcp_server.py` size** (2,044 lines, 81 wrappers) — internal split by domain would reduce merge risk; wrappers themselves are thin and consistent (good).
12. **Classify has no external labeled-corpus test** — rules are tested against their own expectations; a small hand-labeled SMILES set (even 50 structures with charge+family) would put it on T1 footing.
13. **Commission the three missing worked examples** the notes themselves nominate: SSFQ raw worked example, SLS raw Rayleigh-ratio table, independent (V,R)→Nagg pair — each would upgrade a T2 technique to T1 without any code change.
14. **`vant_hoff_enthalpy` magnitude re-check** post-sign-fix: the sign bug fix never got a fresh paper-comparison percentage (the old number is explicitly marked stale in the notes).

### Noted as strengths (keep doing)
- Refusal-first design (`KeyError` on unknown Davies groups, required `system_type`, confidence tiers in reference outputs) — empirically proven on the benchmark (15/15 trap pass post-hardening).
- Dead ends are *removed or disclosed*, never shipped as plausible-looking numbers.
- Test suite breadth (507 test functions → 522 tests incl. parametrization), literature notes discipline, per-entry citations in reference tools, `gaps` list in the orchestrator.

---

## 4. Method & limits of this audit
- Evidence: direct code/doc/test reads, one full test run, static scans (tool counts, description lengths, exception handling, name parity). All numbers in §0 and §C are measured today, not recalled.
- Tier assignments in §1.1 synthesize `literature_validation_notes.md`, `ROADMAP.md`, and test names; I did not re-execute each validation against its source paper.
- Companion repo `surfactantkit-reports` was out of scope (ask was "our toolkit itself").
- Nothing was modified; this report is `TOOLKIT_AUDIT_2026-09-30.md` (new file only).
