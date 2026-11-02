# How Arena contributes to Surfactant-Dynamics-Sim

Prepared by Arena AI (Agent Mode) — 2026-11-02
Sources read: *Synergistic Integration of Scilab-Xcos Dynamic Modeling with SurfactantKit* (14 pp.)
and *Comprehensive Operational Roadmap … for Multiscale Dynamic Colloid Simulation* (11 pp.)
Repo audited: `prashant-kotian/Surfactant-Toolkit` @ v1.3.0 (20,388 lines in `src/surfactantkit`, 982 tests, 123 MCP tools)
Repo inspected: `prashant-kotian/surfactant-dynamics-sim` (scaffold: README, LICENSE, 2 white papers, 1 `.sce`, no CI, no tests, no data contract)

---

## 1. The real gap is narrower than the papers claim — and that is good news

Both documents argue that equilibrium toolkits "fail to capture transient phenomena". Auditing the shipped
library, most of the *physics* needed by all four proposed papers already exists and is already tested:

| Needed by the proposal | Already in SurfactantKit (v1.3.0) | File / function |
|---|---|---|
| Ward–Tordai dynamic adsorption | Yes — Volterra quadrature + bisection, T2-tier, exact Henry fixed-sink + erfc limits | `adsorption.ward_tordai_adsorption` |
| γ(t) forward model with Langmuir kinetics → FDM boundary layer, Richardson-extrapolated | Yes | `dynamic_tension.dynamic_surface_tension_forward` |
| Joint fit of (D, K, Γ∞, k_off) with covariance SEs | Yes | `dynamic_tension.dynamic_surface_tension_ode` |
| DLA / KLA regime discrimination, Rosen–Hua fit | Yes | `dynamic_surface_tension.fit_dla_regime`, `select_dst_kinetic_regime`, `fit_rosen_hua` |
| Equation of state (Langmuir–Szyszkowski, Frumkin, Flory–Huggins, Temkin) + BIC model selection | Yes | `adsorption.szyszkowski_*`, `frumkin_*`, `select_isotherm_model` |
| Gibbs/Marangoni dilational elasticity E₀ (Langmuir and Frumkin) | Yes | `adsorption.gibbs_marangoni_elasticity_langmuir/_frumkin` |
| Micellar relaxation τ₁, τ₂ → k⁻, R (Aniansson–Wall inversion) | Yes, T2 | `kinetics.aniansson_wall_two_regime`, `micelle_relaxation_odes`, `kmc_micelle_kinetics` |
| Non-ionic EDL (planar Poisson–Boltzmann, κ, ψ_s) | Yes | `electrostatics.planar_pb_profile`, `debye_length` |
| Ca = μv/γ, CPP, HLB/HLD, foam stability proxy | Yes | `wetting.capillary_number`, `selfassembly`, `hld`, `foam` |

What is genuinely **missing** — i.e. where new physics has to be written, tested and tiered:

1. **Disjoining-pressure isotherm Π(h)** and **thin-film Reynolds drainage** → `foam.py` explicitly stops at
   stability *proxies* and says so; there is no Π(h), no h(t), no h_crit. Proposal §3.3 cannot be built today.
2. **Oscillatory / LAOD solver.** No small- or large-amplitude area-perturbation solve, no harmonic
   demodulation (I₂ω/Iω), no LVDT closed-form check. Everything in the library is a *step* or *ramp* transient.
3. **Marangoni neck / dynamic-Ca(t) microfluidic model.** `capillary_number` takes a static γ.
4. **Interfacial (dilatational) rheology module.** E′/E″/loss angle are not exported as functions of ω and ε.

So the contribution is not "write new surfactant theory" — it is **a coupling layer, a validated dynamic
solver, three genuinely new modules, and a data contract**. That is a much stronger paper §3.4 story
("we reuse a tested molecular engine") than "we rewrote the thermodynamics in Scilab".

---

## 2. My role, and the limits I will not paper over

Your roadmap §4.2 assigns me the browser-sandbox role: literature mining, curve digitisation, JSON curation,
Drive/GitHub sync. I will own that, **with these honest constraints stated up front** so the plan doesn't
quietly fail in month 3:

- **Paywalled corpora are out of reach.** I cannot authenticate to ScienceDirect, ACS or RSC from this sandbox.
  What I *can* mine legally and reproducibly: PMC / Europe PMC open-access full text, arXiv, Zenodo/figshare
  datasets (e.g. the CC-BY-4.0 *DST-SurfDB* dynamic-tension database), publisher OA articles (many Langmuir /
  Colloids Surf. A papers are CC-BY), and anything you drop into the workspace or Google Drive.
  → **Curated-data channel: you (or Antigravity) supply paywalled PDFs; I extract, calibrate and QC them.**
- **No WebPlotDigitizer GUI.** I digitise headlessly in Python (axis calibration from two fiducials per axis,
  colour/threshold masking, monotone-curve tracing, per-curve residual QC against a log-time spline), and I ship
  the extraction script + calibration points + SHA-256 with every curve, so a reviewer can re-derive the number.
  Output is a digitisation *uncertainty* per point, not just an (x, y) pair — the fits must be weighted by it.
- **No GPU, no Scilab, no Intel Arc in this sandbox.** I write and review `.sce` / `.zcos`, and I build the
  *Python reference oracle* that Antigravity's CVode runs are graded against; the actual native Scilab execution
  stays on the host tier. Parity is enforced numerically, not by eyeballing plots.
- Network is allowlisted to GitHub / PyPI / npm only, so CI, packaging and PR automation are all easy here.

---

## 3. Workstreams (each with a definition of done)

**W0 — Contract & scaffold (first, 1 turn).** In `surfactant-dynamics-sim`: `data/schema/dst.schema.json`,
`laod.schema.json`, `film.schema.json`, `microfluidic.schema.json`; a `surfactantkit → solver` parameter-export
format (the single place where a thermodynamic parameter crosses into Scilab); schema validators + tests;
`.github/workflows/ci.yml` (pytest + schema lint + Sudbury-style "no unfilled provenance" check); pinned
`surfactantkit==1.3.x` so a paper's Figure 3 is reproducible from a tag.
*DoD: `pytest` green on a clean clone; one example record per schema; every record has DOI + licence + hash.*

**W1 — Analytical acceptance suite (Phase 1 of your plan, made falsifiable).**
Your Phase 1 says "verify CVode against exact analytic Ward–Tordai solutions". Ward–Tordai has **no** general
closed form — the library's own docstring says so. I will encode the three limits that *are* exact, and port the
existing T2 proofs as the co-simulation acceptance gate:
(i) fixed sink: Γ = 2(c₀ − c_s)√(Dt/π) at all t;
(ii) dilute Langmuir: Γ = K_H c₀[1 − e^{t/τ} erfc√(t/τ)] to 2×10⁻⁴;
(iii) LVDT small-amplitude: E′ = E₀/(1+2ζ+2ζ²)^½, E″ = E₀(1+ζ)/(1+2ζ+2ζ²)^½ with ζ = √(ω_D/2ω) — this is the
missing test for proposal §3.1's ε < 0.05 baseline, and it is what makes the *nonlinear* result publishable.
*DoD: a Scilab/Xcos run and the Python oracle agree to a stated tolerance over a parameter grid; the tolerance
is written into the repo as a test, not a claim.*

**W2 — `dynamic_sim/lvdt.py` + harmonic analyser.** Small-amplitude solver, then the large-amplitude (ε = 0.05…0.50)
version with Frumkin β and micellar source/sink; FFT demodulation of I₂ω, I₃ω; asymptotic checks (ε → 0 recovers
W2's LVDT; ω → ∞ and ω → 0 recover the elastic and diffusive limits).
*DoD: harmonic ratios converge under grid/strain refinement; ε < 0.02 case matches an independent solver to <1%.*

**W3 — `film.py` in SurfactantKit (closes the §3.3 gap).** Π(h) = Π_vdW + Π_el(h) (reuse `electrostatics.planar_pb`
for the constant-potential/constant-charge double layer) + Π_steric; Reynolds drainage dh/dt under imposed
oscillatory pressure; Gibbs elasticity feedback; h_crit and the subharmonic rupture map.
*DoD: Π(h) matches DLVO asymptotes in both limits; drainage reproduces the Scheludko/Stefan closed forms;
new functions get a validation tier (`test_validation_tiers.py` fails otherwise) — I expect T2/T3 and I will say
so at the call site rather than overclaiming.*

**W4 — Data mining pipeline (my §4.2 brief).** Harvest OA DST/LAOD/film records → digitise → schema-conform →
weighted refit of (D, K, Γ∞, k_a, k_d) → produce a coverage report ("50+ datasets" with the *actual* number, its
licence mix, and which Figure each point came from). Target systems from your Table 1: SDS, AOT, CTAB, DTAB,
Triton X-100, C₁₂E₈, β-casein, 12-s-12.
*DoD: a `data/curated/manifest.json` where every record is re-derivable from a public URL or your upload + hash.*

**W5 — Manuscript-grade evidence, not prose.** Figure-source TSVs straight from the pipeline (the Toolkit already
does this: `results/figs_data/*.tsv`), a LaTeX table of validated vs. new results, and a reproducibility appendix
in the style of `results/REPRODUCIBILITY.md`. I keep the "zero numbers typed by hand" discipline your existing
release already enforces.
*DoD: `make figures` regenerates every submitted figure from the committed data + pinned library version.*

**W6 — Review duty on Antigravity's `.sce`/`.zcos`.** Equation-by-equation audit of the Scilab state-space files,
sign/unit/EOS-consistency checks, and a numeric parity report (see §5 for three defects found in the current
prototype already).

---

## 4. Sprint 1 (what I would do in the next two turns if you say go)

1. W0 scaffold + CI in `surfactant-dynamics-sim`, as a PR from a branch off `main`.
2. W1 acceptance tests in Python (they are portable; Antigravity then points CVode at the same numbers).
3. Fix the three documentary/parameter defects in §5 below, in the white paper's equation block and in
   `scilab/solvers/dynamic_surface_tension_oscillatory.sce`.
4. Seed the data contract with **one** fully worked, openly licensed record (axis calibration, digitised points,
   uncertainty, fitted parameters, residual plot) so the schema is proven by an example, not aspiration.

---

## 5. Three defects I found now (fixing these is itself a contribution)

**5.1 The prototype's parameter set cannot reproduce an SDS transient.** Code Block 1 sets Γ∞ = 3.8×10⁻⁶ mol/m²,
so the *maximum* tension lowering the model can ever produce is nRTΓ∞ = **9.42 mN/m**. At C₀ = 5 mol/m³ real SDS
sits near γ ≈ 40–45 mN/m, i.e. a 28–33 mN/m drop. The parameter set therefore starts (x₀ = 0.8Γ∞ → γ = 57.6 mN/m)
already deep in a compressed monolayer and can only ever oscillate inside a 9.4 mN/m window. Also k_a/k_d =
1.5×10⁻¹/1.2×10⁻³ = **125 m³/mol ≈ 1.25×10⁵ L/mol**, one to two decades above published SDS Langmuir K.
→ Exported parameters must be produced *by* SurfactantKit (`szyszkowski_fit_K` / `dynamic_surface_tension_ode`),
never typed in; the prototype becomes a unit test of the export path.

**5.2 The two equation-of-state forms in §4 of the paper disagree.** Prose: γ = γ₀ − nRTΓ∞ ln(1 + K_L C).
Code: γ = γ₀ + nRTΓ∞ ln(1 − Γ/Γ∞). They are the same physics only when Γ/Γ∞ = K_L C/(1 + K_L C) is carried
along; printed side by side they look like a sign error, and a referee will stop reading there.
→ One canonical EOS block, one symbol table, prose and code generated from the same source.

**5.3 The Ward–Tordai asymptotics are stated on the wrong variable.** Both documents write "γ(t) ∝ t^½ (short time)
and γ(t) ∝ t^−½ (long time)". It is the **surface pressure** Π = γ₀ − γ that scales that way
(Π = 2nRTc₀√(Dt/π) at short time). As stated, Phase 1's acceptance criterion is untestable.

**5.4 (secondary) The two-compartment boundary layer is not mass-conserving.** With a *fixed* δ = √(πD/ω),
j_diff = (2D/δ)(c₀ − c_s) and dc_s/dt = (j_diff − j_ads)/δ gives a bulk reservoir that does not lose what the
surface gains — it is a one-way filter, not a boundary layer. For a fresh interface at millisecond times
(new droplets, tips) that is a physical error, and it is exactly what W1/W2 replace with the Ward–Tordai closure
already in SurfactantKit. Keep the compartment model only as the cheap export for LAOD steady-oscillation runs,
and document that restriction at the call site.

---

## 6. Decisions I need from you

1. Where do my deliverables land: PRs against `surfactant-dynamics-sim/main` (I have `gh` access), or should the
   coupling code live inside `SurfactantKit` (`dynamic_sim/` subpackage) so it ships in the same tested dist?
   (My recommendation: **new physics in SurfactantKit** — it inherits 982 tests, tiers and the MCP server;
   **Scilab artefacts, data and the paper in `surfactant-dynamics-sim`**.)
2. Data channel: OA-only for now, or will you/Drive supply the paywalled PDFs for SDS/CTAB/Triton X-100?
3. Which paper do we aim at first — §3.1 LAOD/JCIS or §3.4 software paper (CACE)? The software paper is the
   cheaper, faster, more certain submission; the LAOD paper is the one with the real novelty claim.
