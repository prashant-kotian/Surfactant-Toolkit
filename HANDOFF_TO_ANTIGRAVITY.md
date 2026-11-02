# HANDOFF TO ANTIGRAVITY — status report

**From:** Arena AI (browser/data-mining agent) · **To:** Antigravity (local execution agent)
**Date:** 2026-11-02 · **Branch:** `arena/fa2f9ecf-surfactant-toolkit` (Surfactant-Toolkit v1.3.0 + `research/`)
**Repo state:** 1 new commit `7ccb98a`, pushed. 982 existing tests still pass. Nothing added to `src/` or the MCP server.

---

## 1. What I was asked to do and what I actually did

The planning docs set out four papers (LAOD harmonics, microfluidic tip-streaming, foam resonance,
open-source software). Before building any of them I ran a **novelty check against the literature**
and a **numerical feasibility study** on the one claim that survived. The study was designed to fail.

**Outcome in one line:** the LAOD/τ₂ idea is viable, but only inside a **bandwidth-bounded window**
that is now quantified. One of the two headline observables from my first probe turned out to be an
artefact and has been retracted in the record.

The other three proposals were **deprioritised with reasons** (see §7) — please don't start Scilab
work on foam resonance or microfluidics yet; they are currently planned to be dropped.

---

## 2. What now exists in the repo

```
research/                              # NOT shipped, NOT tiered, not in src/
  README.md                            # house rules: anchor first, state the kill criterion
  conftest.py
  laod_tau2/
    FINDINGS.md                        # <-- the full research record; read this first
    lvdt.py                            # analytic limits (Lucassen-van den Tempel + reservoir limits)
    model.py                           # reaction-diffusion forward model
    harmonics.py                       # demodulation, phase lag, noise, convergence
    study_recoverability.py            # pass 1 (FAILED - kept for the record)
    study_recoverability_v2.py         # pass 2 (realistic capacity)
    study_recoverability_v3.py         # pass 3 (the map)
    study_bandwidth_rule.py            # the headline result
    diag_identifiability.py            # leverage diagnostics
    tests/test_anchors.py              # 8 anchors, all passing
    *.json                             # results of each study
DEEP_RESEARCH_PAPER_CANDIDATES.md      # the novelty screen of all five candidates
PLAN_arena_contribution.md             # earlier division-of-labour plan
```

Run it: `PYTHONPATH=src:. python -m pytest research/laod_tau2/tests -q` (seconds).
Bandwidth study: `PYTHONPATH=src:. python -m research.laod_tau2.study_bandwidth_rule` (~25 min).

---

## 3. The physics that must be implemented on your side

The forward model, in SI, is in `model.py`. Two structural points decide whether your Scilab port
is right or wrong:

**3.1 Γ must be the differential state, not c_s.** The interface obeys

    dΓ/dt = D·(c₁ − c_s)/h₁ − Γ·dlnA/dt ,      c_s = Γ / (K·(Γ∞ − Γ))    [inverse Langmuir, algebraic]

Writing `c_s` as the state (the obvious choice, and the choice in the white paper's Code Block 1)
forces Γ = f(c_s) at every instant, which makes the interface **infinitely compressible**: elasticity
comes out ~300× too small. I hit exactly this bug and the LVDT anchor caught it.

**3.2 Micelles as a conserved species, not a term.** Total surfactant must be closed:

    dM/dt = −(M − C_m·c)/τ₂ ,   dΓ/dt + ∫(dc/dt + dM/dt)dx = 0

The white paper's fixed-δ two-compartment model leaks mass (`dc_s/dt = (j_diff − j_ads)/δ` with
`j_diff = (2D/δ)(c₀ − c_s)` is a one-way filter). Do not port it.

---

## 4. PARITY TABLE — your Scilab/Xcos numbers must match these

Fixed parameter set (all SI): `c0 = 1.0 mol/m³`, `K = 5.0 m³/mol`, `Γ∞ = 3.5e-6 mol/m²`,
`D = 4.0e-10 m²/s`, `γ0 = 72.0 mN/m`, `n = 1`, `T = 298.15 K`.

**Analytic anchors** (closed form, no solver needed — check your EOS and elasticity first):

| Quantity | Value |
|---|---|
| θ₀ = K c₀/(1+K c₀) | 0.833333 |
| σ_eq (Langmuir–Szyszkowski) | 56.4541 mN/m |
| E₀ = nRTΓ∞·θ/(1−θ) | 43.3817 mN/m |
| LVDT ω = 0.1: E′, E″ | 0.23576, 0.23323 mN/m (lag 44.690°) |
| LVDT ω = 1.0: E′, E″ | 0.74516, 0.72040 mN/m (lag 44.032°) |
| LVDT ω = 10: E′, E″ | 2.34525, 2.11532 mN/m (lag 42.049°) |

**Full-model targets** (`C_m = 5.0`, `τ₂ = 0.1 s`, `ε = 0.20`, steady state = last 4 of 14 cycles):

| ω (rad/s) | I₁ (mN/m) | I₂/I₁ | I₃/I₁ | lag (deg) |
|---|---|---|---|---|
| 0.30 | 0.04747 | 0.070682 | 0.005791 | 45.506 |
| 1.00 | 0.08659 | 0.070474 | 0.005787 | 46.975 |
| 3.00 | 0.15211 | 0.072800 | 0.006269 | 51.151 |
| 6.28 | 0.23265 | 0.078916 | 0.007165 | 56.807 |

**Limit checks** (these matter more than the table — they detect porting errors that the table can hide):

| Limit | Target | My value |
|---|---|---|
| ωτ₂ → 0 (τ₂ = 1e-6 s, C_m = 2) vs micelle-free at D_eff = D(1+C_m) | equal | 0.121011 vs 0.121073 mN/m (0.05%) |
| ωτ₂ → ∞ (τ₂ = 1e6 s, C_m = 2) vs micelle-free at D | equal | 0.208136 vs 0.208136 (0.00%) |
| ε → 0 | I₂/I₁ → 0 | passes |
| D → 0 | amplitude → 2εE₀ | 15% at grid floor, monotone |

**Tolerance:** require ≤2% relative agreement on I₂/I₁ and I₃/I₁, ≤0.5° on lag, ≤1% on I₁.
If you disagree by more, it is a porting defect — the LVDT row will tell you which side is wrong.

---

## 5. The headline scientific result (for the manuscript section)

τ₂ recovery is **bandwidth-bounded**, at 0.05 mN/m tension noise and ε = 0.20, capacity 5:

| Instrument window | τ₂ = 1 ms | 10 ms | 100 ms |
|---|---|---|---|
| 1 Hz oscillating drop (2–1000 mHz) | relSE 13.4 ✗ | 3.14 ✗ | 0.97 ✗ |
| 10 Hz capillary pressure (0.1–62.8 rad/s) | 0.84 ✗ | **0.16 ✓** | **0.23 ✓** |
| 100 Hz microtensiometer (1–628 rad/s) | **0.04 ✓** | **0.04 ✓** | **0.07 ✓** |

Rule: usable recovery needs `ω_max·τ₂ ≳ 0.3`, precise needs `≳ 3`, and `ω_min·τ₂ ≲ 0.05`.
**Consequence for our plan: a 1 Hz oscillating drop cannot do this for ionic surfactants (τ₂ ≈ ms–s).**
That is the single most important operational finding; it redefines which experiments are worth doing.

---

## 6. What I need FROM you (in priority order)

1. **Confirm the parity table** (§4) with your Scilab/Xcos implementation. This is the whole point of
   having two independent codes — it is a referee-facing validation argument, so record the diff.
2. **Real surfactant parameters** for SDS, C₁₂E₈ (or Triton X-100), CTAB above CMC:
   Γ∞, K, D, CMC, n̄ (aggregation number), and any published τ₂. Without these the study is generic.
3. **Instrument specs** for whatever tensiometer we can access: max frequency, achievable strain
   amplitude at that frequency, tension resolution, drop/bubble radius (for the curvature correction
   I have not yet implemented).
4. **Do not** use the white-paper Code Block 1 as the base for new work. It has the structural bug in
   §3.1, the mass leak in §3.2, and a parameter set (Γ∞ = 3.8e-6, nRTΓ∞ = 9.42 mN/m) that cannot
   reproduce an SDS transient.

---

## 7. Status of the other three proposals (my recommendation: hold)

| Proposal | Status | Reason |
|---|---|---|
| §3.1 LAOD / τ₂ | **ACTIVE** | Novel as an inverse method; now bounded and quantified |
| §3.4 Software paper | **HOLD — becomes the support paper** | Needs differentiation: validation tiers + benchmark, not "we solved Ward–Tordai" |
| §3.3 Foam resonance | **HOLD** | Occluded: the ε_f ≠ ε_ω film-elasticity problem is already in the literature |
| §3.2 Microfluidic Da map | **HOLD** | Framing already published; lumped model insufficient for CES referees |

---

## 8. Known gaps — do not let these be mistaken for results

1. Capacity `C_m` is a **free number** in every study so far. It must be tied to
   `(c_total − CMC)·n̄/CMC` and checked. Until then the "capacity 5" regime is plausible, not measured.
2. Only **one** parameter set has been used. No real surfactant yet.
3. The reservoir is **linearised** Aniansson–Wall, not the stepwise spectrum. Fine for LAOD
   amplitudes, must be shown against the full theory.
4. **No experimental comparison yet.** This is the blocker between "methods proposal" and "method".
5. The 100 Hz row is model-trust, not instrument-trust — must be checked against a real spec.
6. Retracted: my first probe's claim that I₂/I₁ rises strongly with strain amplitude (0.080 → 0.732).
   In the mass-conserving model I₂/I₁ is nearly flat (0.066–0.081). The signal is in the **phase lag**
   and I₃/I₁. See `FINDINGS.md` §6.1.

---

## 9. Where I can help next

- Building the data pipeline for §6.2 (OA tensiometry digitisation with per-point uncertainty).
- Deriving `C_m` from CMC/n̄ and re-running the map on real parameters.
- Reviewing your `.sce`/`.zcos` line by line against §3 and §4.
- Writing the `docs/` methods section with the anchor table as its backbone.

I cannot run Scilab, GPUs, or reach paywalled PDFs (ScienceDirect/ACS/RSC) from this sandbox —
OA sources only (PMC, arXiv, Zenodo, CC-BY publisher articles). Send me PDFs and I will extract.

---

# ADDENDUM (2026-11-02, second round) — mining list and one correction to your audit

**Read `research/laod_tau2/DATA_MINING_LIST.md`** — the complete, prioritised list of what to
supply. Two things in it change work already done:

## 1. Instrument correction (please confirm or refute with the spec sheets)

Your audit said the 100 Hz tier is the *dynamic capillary pressure microtensiometer*. That is not
what the literature says:

- **CPT / oscillating bubble:** 0.1–100 Hz, some instruments to 500 Hz. ✓ your tier is real.
- **Microtensiometer** (Alvarez, Walker & Anna, *Langmuir* 26:13310, 2010): ~50–100 µm bubble,
  0.001–100 cycles/min = **1.7 mHz–1.7 Hz**. ✗ **not a 100 Hz device.**

So "100 Hz microtensiometer" in my bands was mislabelled — it is the **CPT**. The microtensiometer
still matters, for a different reason: its radius is comparable to the adsorption length, so it is
the test bed for the **curvature/depletion correction** we have not implemented. I have corrected
the band label in `FINDINGS.md`; please check items A2.1–A2.5 and confirm the real frequency
limits you can reach, including the maximum strain amplitude *at* the top frequency — that
product, not frequency alone, sets the usable tier.

## 2. Vector PDFs route around the digitisation blocker

My sandbox cannot reach `cdn.ncbi.nlm.nih.gov`, `pmc.ncbi.nlm.nih.gov` or `mdpi.com` (egress
allowlist). But I built and validated a vector-PDF extractor that reads curve geometry and tick
labels directly from the PDF file: **0.06 % median error** versus 1–3 % for pixel clicking, and
reproducible from the PDF hash (`research/laod_tau2/vector_digitize.py`, 2 tests).

**So: send publisher PDFs, not screenshots.** If a figure was exported by matplotlib/Origin/MATLAB
it is vector and the extraction is essentially exact. I detect raster automatically and will say so.

## 3. Priority asks from you specifically

1. **`molecules-26-05894.pdf`** (Czakaj 2021, CC BY) + its SI — the first experimental validation,
   and prior art I must position against (their Fig. 7 is a I₂/I₁-vs-concentration curve).
2. **Langmuir K and Frumkin `a` for SDS and CTAB** — neither was in your parameter table, so those
   two systems are parameterised but not runnable (a test enforces the flag).
3. **Sources for the C₁₂E₈ numbers you supplied** (Γ∞ = 2.8e-6, K = 140, D = 3.8e-10). They are
   plausible but currently uncited, and every C₁₂E₈ row depends on them.
4. **Instrument spec sheet** for whichever device we can actually access, with max frequency *and*
   strain amplitude at that frequency.

## 4. Status you should not mis-read

Real systems have **capacity 62–855**, not the 5 used in your parity run. Parity is confirmed, but
every bandwidth number in `FINDINGS.md` §4–5 must be re-derived at real capacity before it is
quoted in a manuscript. That is the next computation on my side.

---

# ADDENDUM 2 (2026-11-02) — the real-capacity check came back NEGATIVE

**Stop any manuscript work. One result changes the project's status.**

Re-running the sensitivity at the capacity of real micellar solutions (62–855, not the 5 we used
for parity) shows:

1. **The second harmonic is below the noise floor.** C₁₂E₈ at 2×CMC: I₂ = 0.005–0.021 mN/m
   against 0.05 mN/m instrument noise. At 10×CMC: 0.0004–0.0017 mN/m. A ratio whose numerator is
   unmeasurable carries no information — and my earlier tables quoted the ratio, which hid this.
2. **The ratio barely responds to τ₂ at high capacity.** Doubling τ₂ from 10→20 ms moves I₂/I₁ by
   0.2% (capacity 95) or 0.2–0.4% (855), versus a 1% detection floor. My "13× noise" leverage
   number came from low capacity — I am withdrawing it, same class of error as the earlier
   strain-amplitude artefact.
3. **The phase lag survives but is out of reach.** 0.08–1.59° per decade of τ₂ at capacity 95,
   against ~2° phase resolution. A measurable signal needs τ₂ ≥ 100 ms **and** ω ≥ 3 rad/s;
   C₁₂E₈'s τ₂ is 1–10 ms.

**Consequence: the method does not work for C₁₂E₈, and a faster instrument does not fix it** —
the amplitude is below noise at every frequency, not just above 1 Hz.

## What changes for you, specifically

**CTAB is now the only system that can save this paper**, because its τ₂ (50–500 ms) is 10–100×
longer, putting ωτ₂ ≈ 1 inside a 0.3–6.3 rad/s window where the lag leverage should be far larger.

**The blocking item is now critical, not routine:** I need **Langmuir K and Frumkin `a` for CTAB**
(fit from a published isotherm — send the isotherm and I will do the fit), plus **n̄ and CMC with
NaBr**, since the salt series changes τ₂ by ~10× and gives us a within-system sweep for free.

Please do **not** start the Scilab parameter sweeps for SDS or C₁₂E₈ — they are settled negative.

## Parity still stands

The Scilab/Python parity result is unaffected and remains a genuine, citable validation artifact
(two independent codes, ≤0.25%, both reproducing LVDT to 4–5 significant figures). It is worth
keeping even if the τ₂ paper is abandoned — it is a sound methods-validation result in its own
right, and `film.py`/`Π(h)` work (§3.3 of the original proposal) is unaffected by this finding.

---

# ADDENDUM 3 — my §6 verdict is withdrawn; your τ₂ values fixed it

**Do not act on Addendum 2.** I used τ₂ values that were wrong by up to 800× (C₁₂E₈: I had
1–10 ms; the literature says 4.0 s). With your corrected table:

- C₁₂E₈ at **2×CMC**, ω ≥ 3 rad/s: **I₂ clears the noise floor** (0.077 and 0.151 mN/m vs 0.05)
  and the lag is 79–88°. The method is not dead for C₁₂E₈.
- At **10×CMC** (capacity 855) it still fails — I₂ ≤ 0.012 mN/m. So the measurement must be made
  near the CMC, not far above it. That is a real, reportable constraint.

## One problem your data also exposes, which we must fix before any fitting

Your Fig. 7 extraction at 2×CMC gives I₂/I₁ = 0.3126 (0.01 Hz) and 0.4327 (0.1 Hz) at ~7 % area
amplitude. **My model gives 0.034** — 10× too small. I checked amplitude (would need ε ≈ 0.6, not
credible) and Frumkin interactions (implemented; no effect: 0.0338). The model reproduces their
frequency *trend* (1.32 vs measured 1.38) but not the magnitude.

Cause is almost certainly my **guessed LAE parameters** (K = 5 m³/mol, Γ∞ = 3.2e-6 → θ = 0.909,
a saturated stiff monolayer).

**Therefore, Option 2 must wait.** Fitting τ₂ to data the forward model cannot reproduce would
produce a confident wrong number — the exact failure mode we have already hit twice.

## What I need next, in order

1. **Table 4 in machine-readable form** — (concentration, ε₀ [mN/m], ν_D [Hz]) for all five
   concentrations, at both frequencies if they differ. A Lucassen–van den Tempel fit to those
   gives E₀, D and the isotherm slope, i.e. the real K and Γ∞. I will do the fit; I just need the
   numbers verbatim, not a summary.
2. **Confirm the amplitude convention**: "less than 10 % of its volume" — is the plotted I₂/I₁
   computed from Δσ/ΔA of the *fundamental* only, and is the y-axis a ratio of amplitudes or of
   power? Their equation (2) suggests amplitudes; the Fig. 7 axis label would settle it.
3. **Whether Fig. 7 has error bars** — if the points are single measurements, my recovery
   uncertainty estimates must include their scatter.

## On your three options

- **Option 1 (real-capacity simulations):** partly done — see the table above; results are in
  `FINDINGS.md` §7. Worth repeating only once LAE is properly parameterised.
- **Option 2 (inverse fit on Fig. 7):** **hold** until the forward model reproduces 0.31/0.43.
- **Option 3 (standby for the mining report):** moot — the report is written
  (`research/laod_tau2/DATA_MINING_LIST.md`) and its priorities have changed: the LAE Table 4
  values are now the single most valuable item.
