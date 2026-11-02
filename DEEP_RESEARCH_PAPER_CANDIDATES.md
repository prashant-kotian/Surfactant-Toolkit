# Deep research memo: what genuinely new work is available in SurfactantKit + Scilab/Xcos

Date: 2026-11-02 · Author: Arena AI · Status: **research only — no project work started**
Inputs: the two white papers; audit of `Surfactant-Toolkit` v1.3.0; literature checks; and a
numerical feasibility probe run from scratch in `/tmp/feas/laod_probe.py` (not committed).

---

## 0. Method, and what "genuine" means here

I judged each candidate on four tests, in this order — a candidate that fails 3 or 4 is dropped
regardless of how attractive it sounds:

- **T1 Existence** — is the effect real, or an artefact of the model?
- **T2 Novelty** — has someone already published the claim? (searched, not assumed)
- **T3 Separability** — can the signature be told apart from a confound that a referee will
  immediately propose? *This is the test that kills most "simulation shows X" papers.*
- **T4 Falsifiability** — is there a specific measurement that would prove it wrong?

I then ran the recommended candidate through a working numerical model to see whether the
signature survives contact with a realistic parameter set. That is the difference between this memo
and the earlier plan document.

---

## 1. The five candidate questions I screened

| # | Candidate question | T2 Novelty | Verdict |
|---|---|---|---|
| A | **Can τ₂ (micellar slow relaxation) be measured from the *harmonic content* of a large-amplitude oscillatory dilatational response — inverting a fast optical measurement for a quantity that currently needs pressure-jump / ultrasonic kit?** | Open in the frequency domain — see §2 | **RECOMMENDED** |
| B | A priori prediction of *dynamic* (not equilibrium) γ(t) from SMILES across a homologous series | Partly done for equilibrium (ethoxylate theory validated on 40 systems); dynamic side open | **Strong runner-up, second paper** |
| C | Frequency-dependent film elasticity governs oscillatory foam-film drainage | Occluded — the ε_f ≠ ε_ω problem is already stated in the literature | Dropped |
| D | Dynamic-Ca(t) regime map for microfluidic tip-streaming | Framed already; lumped model insufficient for CES | Dropped |
| E | Open-source validated toolkit + benchmark (software paper) | Real but crowded | **Keep — but as the *support* paper, not the discovery** |

---

## 2. Why candidate A is genuinely open

What is already published (so you do not waste a year rediscovering it):

- Aniansson–Wall gives τ₁ (µs, monomer exchange) and τ₂ (ms–min, whole-micelle
  formation/breakdown); the τ₂–foamability–DST correlation is Shah's classic result. τ₂ is
  normally obtained by **pressure-jump, temperature-jump or ultrasonic relaxation**.
- Morse, McCormick & Morse, *Phys. Rev. E* **105**, 034602/034603 (2022) built the rigorous
  nonlinear theory of micellar kinetics **and** applied it to adsorption at a quiescent interface.
  So the *forward theory* of "micellar reservoir feeding a relaxing interface" exists.
- Experimentally, nonlinearity in oscillatory surface tension **near and above the CMC has been
  observed** and attributed to micelles (e.g. ethyl lauroyl arginate, *Molecules* 26:5894, 2021);
  local maxima in E′/E″ near the CMC were noted as early as 2006.
- LAOD as a *method* was reviewed in full in 2025 (Risse et al., *Adv. Colloid Interface Sci.*
  345:103625) — general stress decomposition, Lissajous, Chebyshev.

What has **not** been done, as far as this search reaches:

> Nobody has treated the harmonic spectrum of a LAOD response as an **inverse problem for τ₂** —
> i.e. no paper extracts a micellar relaxation time from the *nonlinear* (harmonic) part of an
> oscillatory dilatational measurement, and validates it against an independently measured τ₂ on
> the same system. The existing literature stops at "the response is nonlinear near/above CMC and
> micelles are probably involved".

That is a **method claim**, not a simulation claim. It is the difference between "our model
reproduces the literature" (reject) and "we can now measure X without instrument Y, and here is
the error bar" (JCIS/Langmuir accept).

**The one-sentence novelty statement for the paper:**
*A frequency-domain inverse method that recovers the slow micellar relaxation time from the
higher-harmonic content of a large-amplitude oscillatory dilatational response, with a stated
identifiability criterion that separates micellar reservoir effects from equation-of-state
nonlinearity, and an error budget validated against independently measured τ₂.*

---

## 3. Numerical evidence from the probe

Model: interface mass balance with a dilatational term, kinetic Langmuir adsorption, a
boundary-layer subsurface compartment and a linearised micellar source (Aniansson–Wall τ₂ with
capacity C_m), γ from the Langmuir–Szyszkowski EOS; response demodulated by least-squares Fourier
fit over the last 4 of 40 cycles. System: Γ∞ = 5.0×10⁻⁶ mol/m², K = 2000 M⁻¹, D = 5×10⁻¹⁰ m²/s,
c₀ = 20 mM (above CMC; θ_eq = 0.976).

**Result 1 — the reservoir signature exists and is large.**

At ω = 1 rad/s, ε = 0.30: I₂/I₁ rises from **0.080** (no micelles) to **0.732** at ωτ₂ = 1, and the
fundamental phase lag goes from **−19°** to **−36…−86°**. This is far above the ~1% noise floor that
FT-rheology practice treats as detection threshold — a commercial oscillating-drop tensiometer can
see it.

**Result 2 — the signature is a fingerprint, not a monotone trend.** Sweeping ω at fixed ε = 0.30
with τ₂ = 0.1 s:

| ω (rad/s) | ωτ₂ | I₂/I₁ | I₃/I₁ | lag γ behind A |
|---|---|---|---|---|
| 0.05 | 0.005 | 0.153 | 0.024 | 89° |
| 0.10 | 0.01 | 0.153 | 0.023 | 88° |
| 0.30 | 0.03 | 0.144 | 0.019 | 84° |
| 1.0 | 0.10 | **0.054** | 0.009 | 68° |
| 3.0 | 0.30 | 0.306 | 0.184 | 41° |
| 10.0 | 1.0 | **0.732** | 0.524 | 36° |

There is a **minimum in I₂/I₁ near ωτ₂ ≈ 0.1 and a steep rise through ωτ₂ ≈ 1**. Without micelles,
I₂/I₁ is essentially frequency-flat (0.063 → 0.105 across ε = 0.01→0.40). A non-monotonic
I₂/I₁(ω) with a dip is therefore a qualitative, checkable fingerprint of a relaxing micellar
reservoir — and its location moves with τ₂, which is what makes it invertible.

**Result 3 — the confound test, and the reason this can be a method paper.** The obvious referee
objection: *"the interface is just nonlinear — you can get any I₂/I₁ you like by changing strain
amplitude."* Checked directly: the target (τ₂ = 0.1 s, C_m = 0.2, ε = 0.30, I₂/I₁ = 0.054, lag 68°)
cannot be matched by any strain amplitude without micelles — no-micelle runs give lag 18–22° at
every ε, and reaching I₂/I₁ = 0.054 requires ε ≈ 0.10, which yields lag 22°, not 68°. **The pair
(I₂/I₁, phase lag) is separable in (ε, τ₂) space.** That identifiability statement is the
intellectual core of the paper: it converts a nonlinear fit into a defensible measurement.

**Honest limitations of the probe (must be fixed before any claim is published):**
1. The micellar term is a *linearised* sink, not the full Aniansson–Wall stepwise spectrum.
2. δ = √(πD/ω_ref) is fixed, not mass-conserving — the 89° low-frequency lag is probably an
   artefact of this, and the low-ω branch (ωτ₂ → 0) needs the full Ward–Tordai closure.
3. No Frumkin lateral interaction, no polydispersity, no curvature/depletion (drop geometry).
4. c₀ is a single point; the LAOD literature lives at 0.1–10 × CMC and the C_m capacity must be
   tied to real aggregation numbers.
None of these threatens the qualitative finding; all of them change the numbers.

---

## 4. The paper, specified

**Working title.** *Extracting the slow micellar relaxation time from nonlinear oscillatory
dilatational rheology: a frequency-domain inverse method with an identifiability criterion.*

**Claim.** τ₂ can be recovered from a single LAOD frequency–amplitude sweep on a
tensiometer, without pressure-jump or ultrasonic equipment, with a stated error budget; and the
harmonic fingerprint (a dip in I₂/I₁ near ωτ₂ ≈ 0.1, a rise through ωτ₂ ≈ 1, phase lag tracking)
distinguishes a micellar reservoir from EOS nonlinearity.

**Why it is publishable.** It is a *new measurement capability* at zero extra instrument cost, aimed
at a quantity the foam/emulsion/formulation community already treats as the controlling variable for
dynamic performance. It is falsifiable, and it produces a design rule (which structures/packing give
slow-relaxing micelles).

**The decisive figure.** Measured I₂/I₁(ω) and phase lag (ω) for a micellar system above CMC versus
the same surfactant below CMC, with the inverted τ₂ overlaid against an independently measured τ₂.
Everything else in the paper supports that one figure.

**Validation chain (in order, each step falsifiable):**
1. ε → 0 reproduces the Lucassen–van den Tempel E′(ω), E″(ω) closed form (analytic anchor; the
   library has no such test today).
2. Below CMC the model must show *no* non-monotonic I₂/I₁ — i.e. the effect must vanish when the
   reservoir vanishes.
3. Forward comparison against Morse's micellar-adsorption theory in the small-perturbation limit.
4. τ₂ inversion against at least one system with a *literature* τ₂ from pressure-jump/ultrasonic
   (SDS has τ₂ ≈ 10⁻⁴–10¹ s concentration-dependent; Triton X-100 seconds; C_nEO_m minutes).
5. Sensitivity: what τ₂ range and precision are recoverable given realistic instrument noise
   (±0.05 mN/m) and the usable frequency window of an oscillating-drop tensiometer (≲1 Hz).

**Target journals, in order.** *J. Colloid Interface Sci.* (method + physical claim) →
*Langmuir* → *Colloids Surf. A*. If step 5 fails (τ₂ not recoverable in the achievable frequency
window), the fallback is a *Phys. Rev. E / Soft Matter* theory paper on the identifiability
criterion alone.

**Kill criteria — decide before starting, and honour them:**
- If step 1 or 2 fails → the model is wrong; stop.
- If step 5 shows the usable frequency window cannot reach ωτ₂ ≈ 1 for any real surfactant, the
  method is unmeasurable → convert to a theory paper or drop.
- If a literature search at submission time finds an existing harmonic-inversion paper → reframe as
  a validation/extension, or drop.

**Time to a submittable draft, single person with me doing the code and analysis: 4–6 months**,
dominated by measurement of steps 4–5, not by modelling.

---

## 5. Runner-up (the natural second paper)

**Structure → dynamic γ(t) prediction.** Equilibrium structure-based adsorption theory for
ethoxylates is already validated on 40 systems, but nobody closes the loop to *dynamic* tension and
foamability from SMILES. SurfactantKit already has every link (`classify` → `cpp`/`nagarajan` →
`szyszkowski_fit_K` → `dynamic_tension.*` → `kinetics.*`). The honest scope is a homologous series
(C_nE_m, n = 8–14, m = 4–10) where the trend, not the absolute value, is the claim. This is a
*demonstration* paper: lower risk, lower novelty, good companion to the software paper.

---

## 6. Where Scilab/Xcos is legitimately useful — said plainly

The recommended paper **does not need Xcos**, and pretending otherwise will hurt it. The honest
roles, in descending value:

1. **Numerical independence as a validation argument.** Two independently written codes (Python and
   Scilab/Xcos) agreeing to a stated tolerance on the same benchmark grid is a real, referee-visible
   argument. The library's own tests do this already (Volterra vs FDM duality); extending that
   discipline across languages is defensible.
2. **Event-driven problems only** — pinch-off and film rupture genuinely benefit from zero-crossing
   detection; the LAOD problem does not.
3. **Parameter sweeps as a `.sce` batch** with the model exported as an FMU for archival
   reproducibility.
Make it a *supporting* section of the software paper, not the core of the physics paper.

---

## 7. What I would do next, if you approve

1. Replace the probe's four limitations (§3) with the real physics, reusing SurfactantKit's
   `ward_tordai_adsorption` closure and `aniansson_wall_two_regime` — 1–2 sessions.
2. Run step 5 (the recoverability/sensitivity study) **first**, because it is the kill criterion and
   it is cheap. If τ₂ is not recoverable at realistic noise and frequency, we learn that in days
   rather than months.
3. Only then build the module, the tests, and the figure.

I have not written anything into either repository. Say the word and I will start at step 2.
