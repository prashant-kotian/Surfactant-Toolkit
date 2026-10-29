# Problem Statements for Frontier-LLM Runs (Claude Opus 4.8 / ChatGPT 5.6 + MCP, free web)

Protocol (as agreed): paste ONLY the "PROBLEM" block below into each model, MCP server attached,
web access on. No method guidance, no tool names, no candidate techniques, no expected values —
models figure out everything. Output they must produce: full reasoning + final answer + sources.
We analyse afterwards against `APPLICATION_REFERENCE_RESULTS_2026-10-05.md` and the §4 epistemology
(`APPLICATION_AXIS_RESEARCH.md`).

Shared preamble used in every problem (protocol, not direction):

> You are a surfactant/colloidal scientist solving a research problem. Use any resources you
> have — the web, literature, and any computational tools/MCP servers connected — as you see
> fit; you have full authority over strategy. Rules of scholarship: cite the source of every
> external parameter you use; show your full reasoning and every calculation; report numerical
> answers with uncertainty where derivable; if a required value genuinely cannot be determined,
> say so and explain — never fabricate numbers, citations, or tool results. Output format:
> (1) reasoning, (2) final answers (explicit numbers), (3) sources.

---

## P1 — Optimal-salinity predictions for a co-surfactant SDS formulation (candidate CS-1)

**PROBLEM:**
A colleague works with the formulation **6 % (w/v) sodium dodecyl sulfate + 16 % (w/v) n-butanol
in NaCl brine against hydrocarbon oil phases, 1 : 1 brine:oil, at ambient temperature**. Some
n-alkane oils have already been studied experimentally by an independent group (their data are
in the published literature — find them). She now needs, before any laboratory work:

1. The optimal (Winsor-III, balanced) salinity in g NaCl/100 mL for **n-hexane, n-heptane,
   n-tetradecane, n-hexadecane, and ethyl palmitate** as the oil phase in exactly this
   formulation.
2. An uncertainty range for each prediction, from the uncertainties of your inputs.
3. An audit of which of these predictions is physically attainable in an aqueous NaCl system
   at 25 °C, and which is not.
4. Everything needed to falsify your predictions with a single phase-scan experiment per oil.

State clearly which inputs you found in the literature (with citations) and which you derived.
Do not assume any value you cannot source or derive.

---

## P2 — The aggregate-size distribution of a nonionic surfactant BELOW its CMC (candidate CS-2)

**PROBLEM:**
Textbooks say micelles exist only above the CMC; recent experiments suggest smaller aggregates
exist below it that interfacial-tension methods cannot see. For **octaethylene glycol monododecyl
ether (C12E8) in water at 25 °C**, determine numerically:

1. The mass fractions of surfactant present as monomers, dimers, trimers, and aggregates n ≥ 4
   at total concentrations **0.05, 0.10, 0.12, 0.25, 0.5, 1.0 and 2.0 × CMC**.
2. The free-monomer concentration across the same range.
3. A quantitative argument that your model and parameter choices are physically defensible
   (every parameter sourced or derived, with provenance), including at least one independent
   consistency check of the numbers (e.g., conservation, alternative calculation, or comparison
   with experimental evidence on sub-CMC aggregates).
4. The same census for **sodium dodecyl sulfate** under the same conditions, with an explicit
   statement of which parts of your model do and do not apply to an ionic surfactant.

If different defensible models give different answers, report the spread and say which you
trust and why. Never present a model-conditional number as a measurement.

---

## P3 — Dynamic surface tension of a surfactant solution from first principles (candidate CS-3)

**PROBLEM:**
For **aqueous sodium dodecyl sulfate at 10⁻⁴, 10⁻³ and 10⁻² M, 25 °C, against a clean air
interface** (initial surface tension 72.8 mN/m):

1. Predict the full dynamic surface tension curve γ(t) from 0.1 s to 1000 s at each
   concentration, and the dynamic surface excess Γ(t).
2. Source every physicochemical parameter you need (diffusion coefficient, adsorption
   parameters, etc.) from the literature — cite each — and propagate their uncertainty into
   the curves (at least a sensitivity band).
3. Validate your numerics with a method you trust: two independent calculations, convergence
   tests, or an analytic limiting case — show the comparison quantitatively.
4. Compare your curves against any published dynamic surface tension measurements for this
   system you can find; report agreement or disagreement honestly, including any regime where
   your prediction and the data diverge.

Deliver γ(t) at logarithmically spaced times (0.1, 1, 10, 100, 1000 s) for each concentration
as a table, with uncertainty, plus your validation evidence.

---

## P4 — Stability prediction for a published emulsion from its droplet-size data (candidate CS-4)

**PROBLEM:**
A 5 v/v % paraffin-oil-in-water emulsion stabilized by Tween surfactant was prepared with a
high-shear mixer (d₃₂ = 235 µm at 330 rpm); the authors report that the droplet-size distribution
showed no meaningful change over one month of storage (Pol. J. Chem. Tech. 26(1) 2024, 51;
find the paper). Using only the information in that paper plus literature-sourced physical
properties:

1. Predict the creaming/separation timeline of this emulsion in a 10 cm-tall quiescent column
   at 25 °C: report t₅₀ (time to 50 % separated) and the full separated-fraction curve.
2. Quantify how the assumed shape/width of the droplet-size distribution (at fixed d₃₂) changes
   that prediction — state how much error the common monodisperse approximation incurs here.
3. Use the paper's month-long observation of an unchanged droplet-size distribution to place a
   quantitative constraint (upper bound) on the coalescence rate constant of this emulsion, and
   show how the predicted timeline changes across that bound.
4. State which of your inputs are measured (from the paper), sourced (other literature), or
   assumed, and what experiment would falsify your prediction.

---

## P5 — Cloud points and phase-inversion temperatures of alkyl ethoxylates (candidate CS-5)

**PROBLEM:**
A formulator needs, before buying anything:

1. Predicted cloud points (°C) of pure **C12E5, C12E6, C12E8, C10E6 and C8E6** in water at
   25 °C, and of **C12E6 with 1 g NaCl/100 mL** — each with an uncertainty band traceable to
   your input parameters.
2. The phase-inversion temperature (°C) of **C12E5 with n-hexane and with n-octane** at
   0.01 g NaCl/100 mL.
3. A validation of your predictive scheme against at least one published cloud point or
   phase-inversion temperature you did NOT use as an input — report the residual in °C and
   whether it is inside the accuracy the source of your method claims.
4. A short statement of the domain where your predictions should not be trusted.

Show the sources of every constant; distinguish method constants from system inputs.

---

## P6 — Reproducing (and then correctly reporting) micellization free energies (candidate CS-6)

**PROBLEM:**
Published reports of the standard Gibbs free energy of micellization ΔG°mic for the same
surfactant differ by tens of kJ/mol purely because of thermodynamic convention choices.

1. For **dodecylguanidinium chloride in water at 298 K (CMC = 6.2 mmol/kg, counterion binding
   degree β = 0.74, reported ΔG°mic = −28.3 ± 0.9 kJ/mol)** — determine which thermodynamic
   convention reproduces the published value, and show the reproduction quantitatively
   (Bouchal et al., Int. J. Mol. Sci. 2016, 17, 223 — you may verify the reported values).
2. For **sodium dodecyl sulfate at 298 K (cmc₀ = 0.0083 M, degree of counterion dissociation
   α = 0.272 ± 0.017 from EPR)**, compute ΔG°mic under every standard convention you can
   substantiate from the literature, as a table.
3. Quantify the total spread across conventions for SDS, and identify which convention(s) the
   literature you surveyed predominantly uses — so practitioners can quote a number that is
   comparable across papers.
4. Give a one-sentence recommendation of which convention should be reported, with justification.

Every formula must be cited to a source you did not derive it from yourself.

---

## P7 — A blind first-principles prediction of an unmeasured (or held-out) CMC (candidate CS-7)

**PROBLEM:**
Predict the critical micelle concentration of **sodium hexadecyl sulfate (SHS, C16 alkyl
sulfate) in pure water at 25 °C** from first principles / semi-empirical theory:

1. FIRST write down your predicted CMC (mol/L) with an uncertainty, using no measured value
   for STS itself as an input.
2. THEN establish the credibility of your method: validate it against a set of compounds you
   choose from the literature (report predicted vs measured errors for each).
3. ONLY AFTER items 1–2, search the literature for any measured STS CMC under comparable
   conditions; report it as a held-out comparison with your honest error, clearly marked as
   post-prediction.
4. State what molecular-level assumptions your prediction rests on, and what would make it
   fail.

Chronology of items 1–3 matters and must be visible in your output.

---

# ANALYST APPENDIX — do NOT paste to models

## Reference outputs & checks (per problem)

| P | Reference | Key checks when grading |
|---|---|---|
| P1 | CS-A doc: S\* hexane 3.32 / heptane 3.89 / C14 11.87 / C16 16.32 / ethyl palmitate 3.89 [3.32,4.56]; butanol ΔCc ≈ +2.4; pure-SDS C16 unreachable (>36 g/100 mL) | Their fit must use Bai 2026 windows (findable: Processes 14:2693); k ≈ 0.15–0.16 expected; achievability vs NaCl solubility (~36 g/100 mL H₂O); watch: models using pure-SDS Cc instead of the blend's effective Cc → wrong by ~11× (the system-level trap) |
| P2 | CS-B tables (K=1/CMC derivation, C12E8 rows; SDS rows + ionic caveat) | mass balance ≤1e-15; provenance of K (Rusanov chain); does model handle ENSEMBLE/model-choice discussion or over-claim? Reference K-path is one defensible option — a different sourced K with justification is acceptable if provenance holds |
| P3 | CS-C tables (γ(t) WT vs PDE; worst 4.6 % / 0.08 mN/m) | re-execution with THEIR sourced parameters; cross-validation demanded; Γ uncertainty bands; data comparison honesty (published SDS pendant-drop data exist — findable) |
| P4 | CS-4 JSON: mono t₅₀ 15.1 s; σ=0.5 → 22.7 s (+50 %); bin span 4–103 s; k-bound from no-growth observation | ρ_oil=870 must be flagged ASSUMED (paper may give it — better if they find real); t₅₀-identical-across-σ anomaly: do they notice/explain? (excellent discriminator) |
| P5 | CS-5: C12E5 31.8 (window 25–34 vs Zarate ≈29), C12E6 48.47, C12E8 81.8, C10E6 57.8, C8E6 67.13, salty 46.30; PIT hexane 27.6 / octane 33.0 °C | held-out = Zarate-Munoz 2015 (findable); PIT slope direction check vs Shinoda; round-trip HLD consistency |
| P6 | CS-6: DDGC −28.41 (0.38 %); SDS table −15.11 … −43.67 (spread 28.6) | convention identified by reproduction; citations for every formula; no invented "standard" |
| P7 | CS-7: STS prediction 5.72×10⁻⁴ M (zero-fit); seven-gold worst 0.089 log₁₀ | chronology enforcement (predict-before-look); if they find the measured STS value (~0.5–0.6 mM typically reported — verify their citation) mark as post-hoc confirmation only if item 1 predates it in their transcript |

## Novelty status (answers: "are these results already on the web?")

- **P1** — measured windows for 3 oils: YES, findable (input). All five requested S\* values:
  NO — never computed in any source (new). The ΔCc-butanol quantification: NO.
- **P2** — qualitative premicellar evidence: YES (2025 NMR paper). These numeric mass-fraction
  tables at 7 concentrations for C12E8/SDS: NO — not published (new).
- **P3** — method-independent measured γ(t) for SDS: partially YES (anchors). Our
  cross-solver-validated curve set with input-uncertainty bands: NO (new computation).
- **P4** — the paper + DSD + month-stability: YES (input). t₅₀/polydispersity-error/
  coalescence-bound for that system: NO (new).
- **P5** — cloud points of common ethoxylates: partially YES (retrodiction anchors — by
  design). The complete prediction set incl. PIT pair with uncertainty bands: NO as such.
- **P6** — DDGC paper's own number: YES (anchor, by design). The cross-convention table for
  SDS with Bales inputs: NO as a published table (new).
- **P7** — measured STS CMC: likely YES (held-out, by design — chronology enforced). The
  zero-fit prediction with credibility protocol: NO (new).

## Run logistics

- Attach `surfactantkit-mcp` (121 tools) to each model; both models get IDENTICAL text.
- Save raw transcripts to `claude+gpt-application-output/` (never pushed — uploads rule).
- Analysis gate = re-execute all numbers + provenance audit + pillars §4 + cross-model diff,
  following `verify_augmented.py` patterns.
