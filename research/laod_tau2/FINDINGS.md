# LAOD / τ₂ research record — findings so far

**Status:** kill-criterion study complete. **SUPERSEDED IN PART — see §6.** At capacity 5 the
method looked viable in a bandwidth-bounded window; at the capacity of real micellar solutions
the second-harmonic amplitude falls below the noise floor and the method fails for C12E8.
CTAB is the remaining test case.
**Date:** 2026-11-02 · Branch `arena/fa2f9ecf-surfactant-toolkit`
**Code:** `research/laod_tau2/` (research code, deliberately outside `src/`; not
shipped, not tiered). Anchors: 8 passing tests. Existing suite: 982 passed.

---

## 1. The question, restated operationally

Not "can harmonics be computed" — that is trivial. The question that decides the paper:

> Given a real instrument (bounded frequency window, finite strain amplitude,
> ~0.05 mN/m tension noise), can the micellar slow relaxation time τ₂ be
> **recovered with a useful uncertainty** from the harmonic content of a
> large-amplitude oscillatory dilatational response?

"Useful" was fixed in advance at **relative standard error < 30%**, and the study
was designed to fail. It partly did — which is the most valuable thing that
happened.

---

## 2. What the forward model is (and is not)

`model.py`: 1D monomer diffusion to a planar interface, finite-volume,
geometric grid; Langmuir isotherm at local equilibrium; surface mass balance
including the dilatational term `−Γ dlnA/dt`; and the micelles carried as an
**explicit conserved species** with a linearised Aniansson–Wall exchange

    M_eq(c) = C_m · c ,     dM/dt = −(M − M_eq(c))/τ₂

so total surfactant (bulk + micellar + adsorbed) is conserved *by construction*.
This is the deliberate replacement for the white paper's fixed-δ two-compartment
model, which leaks mass (documented in `model.py`'s module docstring).

Two reservoir parameters, and this distinction turned out to matter more than
anything else in the study:

- **capacity** `C_m = d(micellar)/d(monomer)` — how *much* surfactant is available
- **τ₂** — how *fast* it is released

## 3. Anchors (what makes the numbers trustworthy)

All in `tests/test_anchors.py`, all passing:

| Anchor | Result |
|---|---|
| LVDT closed form reproduced from the PDE, ω = 0.1/1/10 | 4–5 significant figures, no fitting |
| Insoluble limit D → 0 → 2εE₀ | monotone approach, 15% at the grid floor |
| Harmonics → 0 as ε → 0 | second harmonic scales with ε (nonlinear effect, not artefact) |
| Fast reservoir, ωτ₂ → 0 | equals micelle-free system with D_eff = D(1+C_m), 5% |
| Slow reservoir, ωτ₂ → ∞ | equals micelle-free system at same D, 2% |
| ε = 0 stationarity | no spurious drift |

The LVDT anchor was not decorative: it caught a **structural bug** in which I
wrote the subsurface concentration as the differential state instead of Γ. That
version satisfies Γ = f(c_s) at every instant, making the interface infinitely
compressible and reporting ~1/E₀ elasticity — a 300× error that looked like
"physics" until the anchor refuted it. Any group attempting this without the
analytic anchor will publish that bug.

---

## 4. Headline result: the bandwidth rule

`study_bandwidth_rule.py` / `bandwidth_rule.json`. Capacity 5, σ = 0.05 mN/m,
ε = 0.20, one-parameter inversion (capacity known from equilibrium data).

| Instrument window | τ₂ = 1 ms | τ₂ = 10 ms | τ₂ = 100 ms |
|---|---|---|---|
| 1 Hz oscillating drop (2–1000 mHz) | relSE 13.4 ✗ | 3.14 ✗ | 0.97 ✗ |
| 10 Hz capillary-pressure (0.1–62.8 rad/s) | 0.84 ✗ | **0.16 ✓** | **0.23 ✓** |
| 100 Hz oscillating-bubble CPT (1–628 rad/s) | **0.04 ✓** | **0.04 ✓** | **0.07 ✓** |

**Recovery requires the instrument window to bracket ωτ₂ ≈ 1.** Empirically:
usable (relSE < 0.3) needs `ω_max·τ₂ ≳ 0.3`, precise (relSE < 0.1) needs
`ω_max·τ₂ ≳ 3`, and the low-frequency branch must reach `ω_min·τ₂ ≲ 0.05`.

This is the paper's practical deliverable: **an instrument specification per
micellar system**, not a blanket "we can measure it". It also explains, in
advance, why the original white-paper plan (a 1 Hz oscillating drop, ionic
surfactants with τ₂ ≈ ms–s) would have failed — and that negative result is
worth stating in the paper, because it is the failure mode most groups will hit.

## 5. Where recovery is identifiable at all — the map

`study_recoverability_v3.py` / `recoverability_map.json`. Window 0.3–6.28 rad/s
(a deliberately *narrow* case, ωτ₂ never bracketing 1 from below).

| capacity | τ₂ = 1 ms | 10 ms | 100 ms | 300 ms | 1 s |
|---|---|---|---|---|---|
| 5.0 (one-parameter) | ✗ 1.42 | ✗ 0.58 | **✓ 0.19** | ✗ 0.37 | ✗ 14.9 |
| 1.4 (one-parameter) | ✗ 1.37 | ✗ 0.54 | **✓ 0.20** | ✗ 0.40 | ✗ 28.6 |

Two-parameter fits are consistently worse (relSE 0.65–12), i.e. **τ₂ and capacity
are partly degenerate and capacity must come from independent equilibrium data**
(CMC + aggregation number + mass balance). That is a stated protocol requirement,
not a defect — and it is exactly the kind of honest constraint referees reward.

Note τ₂ = 1 s fails here only because the window's *lower* end (0.3 rad/s) is
above 0.05/τ₂; the bandwidth study shows a 1 Hz window covers it fine when the
sweep starts at 0.01 rad/s. Both edges matter.

---

## 6. Two negative results worth recording

**6.1 Strain amplitude is a decoy, and a dangerous one.** A quick probe with a
non-conservative boundary-layer reservoir suggested I₂/I₁ rises steeply with
strain amplitude (0.080 → 0.732 across ε = 0.1→0.4). In the mass-conserving
model, I₂/I₁ at fixed ε is almost *insensitive* to both coverage and reservoir
strength: it sits at 0.066–0.081 for θ = 0.56–0.98, capacity 0.2–5, ω = 0.3–6.3.
**If we had published the probe's "large, distinctive signature", it would have
been an artefact.** The effect lives in the *phase lag* (0.14° → 6.7° per doubling
of τ₂) and in I₃/I₁, not in I₂/I₁.

**6.2 Coverage is not the lever; capacity is.** Sweeping θ = 0.56 → 0.98 at fixed
capacity moved nothing. Sweeping capacity 0.2 → 5 raised τ₂ leverage from 0.9× to
13× the noise floor. Pass 1 of the recoverability study failed *entirely* because
it used capacity 0.2 — a reservoir 25× too weak, which corresponds to "all the
surfactant is monomer", i.e. not a micellar solution at all. Lesson recorded: a
sensitivity study is only as meaningful as the parameter regime it is run in, and
the regime must be justified from equilibrium data before the run.

---

## 7. What this does and does not yet establish

**Established (model-internal, anchored):** the forward physics reaches the
analytic limits; within a bracketed bandwidth the inversion recovers τ₂ to 0.04–0.23
relative SE at realistic tension noise; the recovery is bandwidth-limited by a
quantified rule; τ₂/capacity degeneracy requires independent capacity input.

**Not yet established — required before any submission:**

1. **Real surfactant parameters.** Everything here uses one representative
   parameter set (Γ∞ = 3.5×10⁻⁶ mol/m², K = 5 m³/mol, D = 4×10⁻¹⁰ m²/s,
   capacity 5). Each of these must be sourced or measured per system, and the
   capacity must be derived from an actual aggregation number and CMC.
2. **The C_m ↔ micelle physics link.** Capacity is currently a free number.
   It must be tied to (c_total − CMC)·n̄ / CMC and checked for consistency.
3. **Nonlinear Aniansson–Wall.** The linearised reservoir is fine for LAOD
   amplitudes but must be shown to be, against the stepwise spectrum.
4. **An experimental comparison.** At least one published oscillatory
   dilatational dataset above the CMC, re-analysed, with the inverted τ₂
   compared to an independently measured one. Without this the paper is a
   methods proposal, not a method.
5. **Instrument reality check.** The 100 Hz microtensiometer result must be
   confirmed against a real instrument's specification (radius, acquisition
   rate, strain capability at that frequency), not assumed.

---

## 8. Honest assessment of the paper

The result is **not** "harmonics reveal everything" — it is narrower and better:
*a bounded, quantified regime in which a nonlinear surface-rheology measurement
can replace dedicated relaxation equipment, plus the condition under which it
cannot.* The honest title is about the bound, not the capability.

Realistic outcome: one solid methods paper (JCIS/Langmuir) **if** items 1–5 above
land. The failure mode is now identified precisely enough to avoid: an
under-reservoired parameter set and a bandwidth too narrow to reach ωτ₂ ≈ 1
produce a confident, wrong, unpublishable answer.

**Reproduce:** `PYTHONPATH=src:. python -m pytest research/laod_tau2/tests -q`,
then the three `study_*.py` modules (the bandwidth study takes ~25 min).

---

# 6. DECISIVE NEGATIVE RESULT (2026-11-02): at real capacity the signal is below the noise

This supersedes §5 above and should be read before anything else in this file.

The feasibility study used capacity 5. Real systems have **capacity 62–855**. Re-running at real
capacity changes the conclusion completely, and the reason is elementary once seen:

> A large capacity means the micelles buffer the interface. Γ barely changes, so γ barely
> changes, so the harmonic **amplitudes** collapse even where the **ratios** look informative.
> A ratio whose numerator is below the noise floor carries no information.

**Absolute amplitudes** (C₁₂E₈, ε = 0.20, noise floor 0.05 mN/m, `real_amplitude.json`):

| capacity | I₁ (mN/m) | I₂ (mN/m) | I₂/I₁ | I₂ above noise? |
|---|---|---|---|---|
| 95 (2×CMC) | 0.072–0.326 | 0.005–0.021 | 0.064–0.070 | **no** |
| 855 (10×CMC) | 0.005–0.024 | 0.0004–0.0017 | 0.071 | **no** |

Across a capacity sweep 0.2 → 855 at ω = 3, ε = 0.20: I₂ falls 0.0216 → 0.0009 mN/m — a factor of
24 — while I₂/I₁ moves only 0.0659 → 0.0713. **The ratio was hiding a collapse in absolute signal.**

**And the ratio carries almost no τ₂ information at high capacity.** Doubling τ₂ from 10 → 20 ms
changes I₂/I₁ by **0.2%** at capacity ≥ 5 (0.2–0.4% at 855), against a 1% detection floor. The
earlier "leverage 13× noise" figure came from *low* capacity, i.e. from a regime real micellar
solutions are not in. That number was an artefact of the parameter regime, not a property of the
method — the same class of error as §6.1.

**The phase lag survives, but only just, and only where we cannot look.** At capacity 95 the lag
changes 0.08–1.59° per decade of τ₂, versus an instrument's ~2° phase resolution. A measurable
3–5° change needs τ₂ ≳ 100 ms **and** ω ≥ 3 rad/s. C₁₂E₈'s τ₂ is 1–10 ms. So even the surviving
observable is outside the required window for this system.

## Verdict

**The method does not work for C₁₂E₈ with a 1 Hz pendant drop, and this is not a bandwidth
problem that a faster instrument fixes** — the second-harmonic amplitude is below the noise floor
at every frequency in the window. Two escape routes remain, and both must be tested before any
manuscript work resumes:

1. **CTAB.** Its τ₂ (50–500 ms) is 10–100× longer, so ωτ₂ ≈ 1 is reachable at ~2–6 rad/s, and
   the lag leverage should be 10–100× larger than C₁₂E₈'s. CTAB, not C₁₂E₈ or SDS, is the test
   case that decides whether this paper exists at all. It needs a Langmuir K (missing — see
   `DATA_MINING_LIST.md` B1).
2. **Higher strain amplitude.** I₂ scales steeply with ε, so ε = 0.30–0.40 (if the instrument and
   the monolayer tolerate it) may lift I₂ above the noise. Untested; the model's linearised
   reservoir is least reliable there, so this must be checked against the stepwise theory before
   being believed.

If CTAB also fails, the honest conclusion is that LAOD harmonic inversion is **not**
competitive with pressure-jump or ultrasonic relaxation for τ₂, and the paper should not be
written. That is a cheap, publishable-in-a-paragraph negative result, not a year lost.

**Test order for next session (cheapest decisive test first):**
1. CTAB at 2×CMC, τ₂ = 50 and 500 ms, ω = 0.3–6.28, ε = 0.20 → does I₂ clear 0.05 mN/m?
2. If yes: ε = 0.30, 0.40 → margin grows or saturates?
3. If no: stop, write the negative result.

---

# 7. §6 IS WITHDRAWN — I used the wrong τ₂ (2026-11-02, second correction)

Antigravity supplied the literature τ₂ values. Mine were wrong by up to three orders of magnitude:

| System | my τ₂ | literature τ₂ | error |
|---|---|---|---|
| C₁₂E₈ | 1–10 ms | **4.0 s** | ~800× too short |
| SDS | 2–40 ms | 50–250 mM; **38 ms peak at 200 mM** | range wrong |
| CTAB | 50–500 ms | **65 ms** (480 ms with NaBr) | range wrong |
| Triton X-100 | — | 3.5 s | not in my table |
| Brij 35 | — | 80 s | not in my table |

My §6 verdict ("C₁₂E₈ fails, method dead") rested on τ₂ = 1–10 ms, which put ωτ₂ ≈ 0.06 —
far outside the observable window. With τ₂ = 4.0 s, ωτ₂ = 25 and the picture changes:

**C₁₂E₈ at 2×CMC, ε = 0.20 (corrected τ₂ = 4.0 s):**

| ω (rad/s) | I₁ | I₂ | I₂/I₁ | lag (°) | I₂ above 0.05 noise |
|---|---|---|---|---|---|
| 0.30 | 0.090 | 0.0081 | 0.091 | 69.6 | no |
| 1.00 | 0.266 | 0.0261 | 0.098 | 81.4 | no |
| 3.00 | 0.785 | 0.0768 | 0.098 | 82.7 | **YES** |
| 6.28 | 1.618 | 0.1515 | 0.094 | 78.7 | **YES** |

So the §6 blanket claim is withdrawn. What survives of it, correctly stated:

- **The capacity–signal trade-off is real.** At 10×CMC (capacity 855) I₂ falls to
  0.0006–0.012 mN/m and does not clear the noise floor at any frequency. The measurement must be
  made close to the CMC, not far above it — the opposite of what one assumes.
- **The usable window is narrow but not empty:** ω ≥ 3 rad/s (≈ 0.5 Hz) at 2×CMC.
- The phase lag at ω ≥ 3 is large (79–88°) and moves with τ₂, which is now the more promising
  observable. Not yet quantified against noise at real capacity — that is the next run.

**Lesson (twice now, in the same direction):** every negative verdict I have produced came from a
guessed parameter, not from the model. Parameters must be sourced *before* a sensitivity study,
never after.

---

# 8. A serious problem: the model cannot reproduce Czakaj's measured harmonic ratios

Antigravity vector-digitised Czakaj's Fig. 7 (LAE, CC BY) at 2×CMC: **I₂/I₁ = 0.3126 at 0.01 Hz
and 0.4327 at 0.1 Hz**, with drop volume oscillations < 10 % (area amplitude ε ≈ 0.07).

My model at ε = 0.07, matched frequency: **I₂/I₁ = 0.034** — an order of magnitude *too small*.
Checked and excluded as the cause:

- **Amplitude:** their values would need ε ≈ 0.6 in my model (a 60 % area swing on a 10 µL drop is
  not credible), so this is not an amplitude misunderstanding.
- **Frumkin lateral interactions** (their Table 4 reports a ≈ 0.8): implementing the full Frumkin
  isotherm moved I₂/I₁ from 0.0341 to **0.0338** — no effect. Excluded.
- **Frequency dependence:** the model reproduces their *rise* with frequency (ratio 1.38 measured,
  ~1.32 modelled), so the qualitative mechanism — a slow reservoir catching up as ω increases —
  is right even though the magnitude is wrong.

**Most likely cause: guessed parameters.** I fabricated K = 5 m³/mol and Γ∞ = 3.2×10⁻⁶ m²/mol⁻¹
for LAE (giving θ = 0.909, a nearly saturated, stiff monolayer with little room to be nonlinear).
The correct values are recoverable from **their own Table 4** (ε₀ = 30.1…12.6 mN/m, ν_D =
0.0032…0.06 Hz): a Lucassen–van den Tempel fit to (ε₀, ν_D) at several concentrations gives E₀,
D and the isotherm slope directly. **That is the next computation, and it needs only data already
extracted.**

Until the model reproduces 0.31 / 0.43 for LAE, **no τ₂ inversion on real data is meaningful** —
which is precisely why Option 2 (inverse fitting on the Czakaj points) must not be run yet.

---

# 9. [WITHDRAWN — see §10] The LAE peak at the CMC: a conclusion I over-generalised

> **This section's headline claim is withdrawn.** LAE is a multicomponent, hydrolysing system
> (Czakaj et al.'s own title), so failing to fit it says nothing definitive about pure systems.
> The *observations* in 9.1/9.2 stand as facts about LAE; the *conclusion* does not. §10 has the
> corrected analysis.

Antigravity supplied Table 4 and the amplitude convention. Two checks, both easy, and together
they decide the question. This is the most important entry in this file.

## 9.1 The model cannot be parameterised to fit LAE, and the reason is diagnostic

Table 4 at c = 1.0 mM gives **E₀ = 12.6 mN/m, ν_D = 0.06 Hz**. Those are *two equations*, and the
state has more than two unknowns, so they do not pin (Γ∞, K, D). Sweeping the one free parameter
(coverage θ₀, with E₀ and ν_D held at their measured values) shows the model's I₂/I₁ at ε = 0.07,
0.1 Hz lies between **0.032 and 0.106 for every coverage from 0.1 to 0.95** — against a measured
**0.4327**. The model cannot reach the observation by any choice of parameters.

So my earlier "guessed parameters" diagnosis was wrong. It was not the substitution of K and Γ∞
that caused the factor-of-ten gap; the *structure* of the model caps I₂/I₁ near 0.1 at this
amplitude.

## 9.2 The measured peak sits at the CMC, where the micellar reservoir is empty

Czakaj's Fig. 7 at 0.1 Hz: 0.039 (0.2 mM), 0.019 (0.5), 0.222 (0.8), **0.433 (1.0 = CMC)**,
0.359 (1.5). The maximum is **at the CMC**, and it *declines* above it.

- At c = CMC the total concentration equals the CMC, so the micellar capacity is **zero**.
  A zero-capacity state cannot be the most nonlinear one if micellar buffering is the cause.
- Above the CMC the capacity grows with c, so a buffering mechanism predicts a **monotone rise**.
  The data show a peak and a decline.
- Quantitatively: the reservoir moves I₂/I₁ by ~8% across a 4000-fold capacity range in my model;
  the data move **23-fold** across a 3-fold concentration range. The mechanism is off by two
  orders of magnitude in *sensitivity*, not merely in scale.

Czakaj et al. attribute the effect to interfacial reorganisation and heterodimer (LAE–dodecanoate /
LAE–LAS) formation — a **compositional** nonlinearity, with the micelles as a secondary actor. Their
attribution is consistent with their own data; mine is not.

## 9.3 Consequence for the project

**The τ₂-from-harmonics premise is falsified at the mechanism level.** The observable is dominated
by composition-dependent interfacial reorganisation, not by the slow micellar relaxation time, so
inverting it for τ₂ would be attributing a compositional effect to a kinetic one — exactly the
error the method was supposed to avoid.

Options, in order of honesty:

1. **Drop the τ₂ paper.** Report the negative result: *"the second-harmonic ratio of oscillatory
   dilatational rheology near the CMC is dominated by interfacial compositional reorganisation,
   not by micellar relaxation; it is therefore not a τ₂ probe."* This is a real, useful,
   publishable-in-a-letter finding, and it is now supported by two independent arguments (§9.2).
2. **Re-scope to what the model CAN do:** the phase lag and I₂/I₁ as a *linear-viscoelastic
   consistency test* (Eq. 1 reproduced across their five concentrations) — a validation paper for
   SurfactantKit's LVDT/Frumkin machinery, not a discovery.
3. **Attack composition explicitly** — a multicomponent Frumkin model with an explicit
   heterodimer equilibrium. That is a different, larger project, and it is where the real physics
   of this system lives.

What is **not** an option: proceeding to Option 2 of Antigravity's list (inverse fitting τ₂ on
Fig. 7). The forward model cannot reproduce the data, and the data's shape contradicts the
mechanism. Any τ₂ extracted would be an artefact.

**Verified in scratch, worth repeating:** both attacks on LAE fail — the ratio is undefined for a
vanishing fundamental near the CMC (a second harmonic appears without a prior fundamental, which
itself indicates a discontinuity, i.e. reorganisation rather than smooth desorption).

---

# 10. §9 IS WITHDRAWN, AGAIN — and this time the method survives, with a quantified window

Antigravity's rebuttal was right on the substantive point: **LAE is an inherently multicomponent
system** (their own title), hydrolysing to LAS/dodecanoate, forming heterodimers. Testing a
single-component micellar model against it is a category error, and my §9 conclusion ("the τ₂
premise is falsified at the mechanism level") over-generalised from that one system to all systems.
Withdrawn.

But two of the claims built on it were assertions I had not tested. I tested them.

## 10.1 Detectable is not invertible — and I had only shown the first

I₂ clearing a noise floor does not mean it responds to τ₂. Both must hold **at the same frequency**.
At C₁₂E₈, 2×CMC, capacity 95, ε = 0.20, noise 0.05 mN/m and 2° phase:

| ω | ωτ₂ | I₂ | detect I₂? | Δlag per 2× τ₂ | sensitive? |
|---|---|---|---|---|---|
| 0.10 | 0.2–0.4 | 0.003 | no | 5.2–8.3° | yes (but unmeasurable) |
| 0.50 | 1–2 | 0.010 | no | 5.6–8.9° | yes (but unmeasurable) |
| 1.00 | 2–4 | 0.019 | no | 2.1–5.5° | yes (but unmeasurable) |
| 3.00 | 6–12 | 0.055 | yes | 0.2–2.8° | marginal |
| **6.28** | **12–25** | **0.113–0.152** | **yes** | **3.3–7.0°** | **yes** |

There is a genuine sweet spot at ω ≈ 6.28 rad/s, and it exists because the lag leverage peaks at
ωτ₂ ≈ 1 while the amplitude keeps growing with ω. **I₂/I₁ is nearly useless** (moves 0.004 on
doubling τ₂) — the phase lag is the operative observable, as suspected, but only at the top of the
instrument window. My §7 statement that the lag was "not yet quantified against noise" is now
resolved in the method's favour.

## 10.2 The concentration window is narrower than claimed, and lopsided

Sweeping concentration at ω = 6.28, τ₂ = 4 s:

| c/CMC | capacity | I₂ | lag | Δlag (4→8 s) | usable? |
|---|---|---|---|---|---|
| 1.0 | 0 | 0.129 | 32.4° | **0.00** | detect yes, **no sensitivity** |
| 1.2 | 19 | 0.052 | 50.4° | 7.5° | yes |
| 1.5 | 48 | 0.183 | 67.9° | 10.2° | yes (best) |
| 2.0 | 95 | 0.152 | 78.7° | 7.0° | yes |
| 3.0 | 190 | 0.080 | 84.3° | 3.4° | yes |
| 5.0 | 380 | 0.035 | 86.8° | 1.4° | **no** (below noise) |
| 10.0 | 855 | 0.012 | 88.0° | 0.3° | no |

**The window is ~1.2–3× CMC, centred near 1.5×, and it shuts abruptly between 3× and 5×.** Note the
shape is *not* monotone: I₂ is largest at 1.5×, not at the lowest capacity, because at exactly 1×
the monolayer is compressed enough that the reservoir is empty and the lag is nearly τ₂-independent.
Antigravity's "1.5–3×" is a good first guess; the measurement says 1.2–3× with an optimum at 1.5×,
and gives the reason it fails on both sides.

## 10.3 What survives, and what the paper actually is

**Survives:** a validated forward model; a two-code parity result; a quantified operational window in
both frequency and concentration; and a computed sensitivity of the phase lag to τ₂ for a real
nonionic system.

**Does not survive:** any claim that LAOD harmonics are a general τ₂ probe. The operable window is
narrow, and the observable is the phase lag rather than the harmonic ratio.

**Still untested, and the biggest remaining objection:** with an optimized design (ω = 6.28,
c = 1.5×CMC, ε = 0.20, the lag as observable), how precise is a recovered τ₂? That is a single
least-squares run and it decides whether the paper's headline is "τ₂ to ±30%" or "τ₂ only to a
factor of two". It must be run before the manuscript's abstract is written — everything else is
now in place for it.

**Also unresolved:** whether the 1.2–3× window is generic or specific to C₁₂E₈'s τ₂/capacity
combination. It should be repeatable in one session for SDS and CTAB once their isotherms exist.

## 10.4 On LAE, one correction to the rebuttal

The "Czakaj criterion" as phrased ("I₂/I₁ > 0.2 near the CMC implies compositional reorganisation")
is supported by the evidence only weakly, because **the model cannot produce I₂/I₁ > 0.11 for any
coverage, single- or multi-component, at ε = 0.07**. So the model's ceiling does most of the work in
that argument, and a robust criterion needs an independent handle — e.g. the *frequency* signature
(a compositional process has no reason to sit at ωτ₂ ≈ 1) or the *strain-amplitude* signature
(dimer expulsion should saturate differently from desorption). Recorded as a caveat so the paper
does not lean on a number the model itself produced.

---

# 11. RECOVERY UNCERTAINTY — the answer, and the constraint that actually matters

Design tested: C₁₂E₈, c = 1.5×CMC (0.106 mM, capacity 47.5), τ₂ = 4.0 s, ε = 0.20,
σ = 0.05 mN/m, sweep over the usable window. `study_recovery_uncertainty.py`,
`recovery_uncertainty.json`.

## 11.1 Headline numbers

| ω (rad/s) | lag (°) | Δlag per e-fold of τ₂ | relSE (Δ-method) |
|---|---|---|---|
| 3.00 | 77.2 | 6.0° | 0.011 |
| 4.00 | 74.5 | 8.9° | 0.006 |
| 5.00 | 71.6 | 11.2° | 0.004 |
| 6.28 | 67.9 | 13.3° | 0.003 |

Joint least squares: **1-parameter τ₂ = 4.000 s, relSE = 0.002**; 2-parameter
τ₂ = 4.000 s, relSE = 0.030, **capacity relSE = 0.00, correlation = 1.00**.

Two conclusions that hold regardless of the caveat below:
- **At one concentration the (τ₂, capacity) pair is perfectly degenerate (corr = 1.00).** A
  two-parameter inversion is meaningless from a single concentration; capacity must come from the
  equilibrium closure, or a concentration series must break it. This is the same degeneracy the
  earlier studies found, now quantified exactly.
- The phase lag is a strong observable: **13° per e-fold of τ₂** at ω = 6.28.

## 11.2 The constraint that actually governs the experiment: frames per cycle

The paragraph above assumes the waveform is sampled densely. A tensiometer records *frames*, and
the sub-0.1° uncertainties are achievable only if there are enough frames per cycle. Degrading the
sampling honestly:

| frames/cycle | lag estimate | scatter (sd) |
|---|---|---|
| **1** | **164°, 266°, 194°, 254°** (true ≈ 68°) | **~90°** |
| 2 | 82° | 43° |
| 5 | 67.9° | 0.33° |
| 20 | 67.9° | 0.18° |

**At one frame per cycle the phase estimate is not merely imprecise, it is wrong** — it depends on
where the samples happen to fall and returned 164–266° in four phase offsets. My first pass
reported relSE = 0.045 at ω = 6.28 by luck of the sampling phase; that number is withdrawn.

So the design requirement is **≥5 frames per cycle**, and the achievable precision follows:

| camera | usable ω for 5 frames/cycle | relSE at that ω |
|---|---|---|
| **5 fps** | 6.28 rad/s (1 Hz) | **≈ 0.025** |
| 1 fps | 1.26 rad/s (0.2 Hz) | ≈ 0.5–0.8 (unusable as a precision method) |

**This is the binding constraint on the whole project, and it is an instrument question, not a
theory question.** Everything else — concentration window, bandwidth, observable choice — is now
settled. The one thing nobody has established is the frame rate a PAT-1M uses *during oscillatory
measurements*, which is not the same as the 1 s profile interval quoted for its dynamic-tension
mode. If it is ≥5 fps at 1 Hz, the method delivers τ₂ to ~2.5%. If it is 1 fps, the method cannot
be done on that instrument at any frequency, and needs a faster camera or a capillary-pressure
tensiometer.

## 11.3 Status of the paper's claims

| Claim | Status |
|---|---|
| Forward model validated against LVDT, two independent codes | **solid** |
| Phase lag is the operative observable; I₂/I₁ is not | **solid** |
| Operational window ~1.2–3× CMC (optimum 1.5×) | **computed for C₁₂E₈; generality untested** |
| τ₂ recoverable to ~2.5% | **conditional on ≥5 frames/cycle** |
| (τ₂, capacity) separable from one concentration | **NO — |corr| ≈ 0.6; needs a series, and even then only to |corr| 0.1–0.6 (§13.8)** |
| Method transfers to SDS/CTAB | **untested; isotherms now supplied** |
| Czakaj criterion (compositional vs kinetic) | **needs an independent signature, not the model's 0.11 ceiling** |

The manuscript can be drafted around the first four rows. The fifth is a required addition:
**a concentration series is not optional for the two-parameter problem.**

---

# 12. Frame rate resolved; Protocol B verified; CTAB/SDS parameters need correcting

## 12.1 The sampling constraint is lifted

Antigravity confirms the PAT-1M runs at **25–50 fps** in oscillatory mode (the "1 profile/second"
in Czakaj's methods is their static DST mode). At 1 Hz that is 25–50 frames per cycle against my
≥5 requirement, so the scatter is ~0.18–0.33° and **τ₂ recoverable to ~2.5% stands**. The CPT
transducer at 10 kHz gives 100 samples/cycle even at 100 Hz, so the high-frequency tier is fine too.

This closes the last open hardware question in §11.2.

## 12.2 Protocol B verified — and better than expected

Joint fit of (τ₂, capacity-scale) across a concentration series, synthetic truth τ₂ = 4.0 s:

| series | τ₂ recovered | relSE | corr(τ₂, C_m) |
|---|---|---|---|
| 1.5× only | 4.000 | 0.012 | **−0.61** |
| 1.3×, 1.5× | 4.000 | 0.007 | −0.37 |
| 1.3×, 1.5×, 2.0× | 4.000 | 0.006 | −0.37 |
| 1.2×, 1.5×, 2.0×, 3.0× | 4.000 | 0.006 | −0.24 |

So a **two-concentration series already breaks the degeneracy** (corr −0.61 → −0.37); a third
buys little. Antigravity's Protocol B is correct, and the practical recommendation is
**two runs (1.3× and 1.5× CMC)**, not three. Note this differs from the single-concentration
result in §11.1, where capacity was *fixed by closure* and the correlation was 1.00 by construction
— the series helps only when C_m is genuinely free.

## 12.3 The supplied CTAB and SDS isotherm constants are wrong, and the tests now catch it

| system | K supplied | implied θ at CMC | implied E₀ | E₀ needed for θ≈0.7 |
|---|---|---|---|---|
| CTAB | 18.5 m³/mol | 0.945 | **270 mN/m** | K ≈ 1.7 m³/mol |
| SDS | 1.25 m³/mol | 0.911 | **198 mN/m** | K ≈ 0.16 m³/mol |

Real single-chain surfactants show E₀ ≈ 20–30 mN/m at the CMC, which corresponds to θ ≈ 0.6–0.8.
At θ = 0.95 the monolayer is saturated and has almost no room to be compressed, which is why I₂
collapsed to 0.003–0.011 mN/m in my CTAB/SDS runs even where the *lag* looked fine. **The
"CTAB is directly solvable on a PAT-1M" conclusion therefore does not yet follow** — with corrected
K the coverage moves to the compressible regime and the amplitudes will rise, but that must be
recomputed rather than assumed.

**A correction to my own reasoning, recorded:** I first wrote that "E₀ can never exceed nRTΓ∞".
That is false — for Langmuir, E₀ = nRTΓ∞·θ/(1−θ) diverges as θ→1 (the library's own
`gibbs_marangoni_elasticity_langmuir` docstring warns of exactly this). nRTΓ∞ is a *prefactor*, not
a ceiling. The correct statement is the one now in the test file: a *sensible* E₀ pins θ to
0.6–0.8, so any isotherm implying θ > 0.9 has too large a K.

## 12.4 What is ready to draft, and what is not

| Manuscript element | Status |
|---|---|
| Forward model + two-code parity + LVDT anchors | **ready** |
| Phase lag as the operative observable (13°/e-fold) | **ready** |
| Sampling requirement (≥5 frames/cycle; PAT-1M gives 25–50) | **ready** |
| τ₂ recovered to ~2.5% (Protocol A) | **ready** |
| Protocol B: two concentrations break the degeneracy | **ready** — new result |
| C₁₂E₈ concentration window 1.2–3× CMC | **ready (one system)** |
| CTAB / SDS applicability | **resolved in §13** — CTAB works, but only ≥16 Hz (CPT band), not on a PAT |
| Czakaj criterion | **needs an independent signature** (§10.4) |

The three "ready" rows that depend on numbers (window, 2.5%, Protocol B) are all computed for
**one surfactant**. The manuscript's weakest point is that every quantitative claim rests on
C₁₂E₈; the corrected CTAB run is the cheapest fix. **→ Done in §13, together with a
self-correction: those C₁₂E₈ rows were themselves run at a saturated isotherm (§13.5).**

---

# 13. Generality: the method transfers to CTAB — but not to a PAT-1M, and the reason is a frequency problem

Two things happened in this section. First, the CTAB/SDS blocker from §12 was cleared
honestly. Second, clearing it exposed a defect in **my own** headline studies, which I have to
record because it changes how every number above must be read.

## 13.1 Clearing the blocker by inverting, not guessing

§12 found that the supplied Frumkin constants put CTAB/SDS at θ = 0.94–0.96. The fix is not to
guess a better K. It is to note that **E₀ is the measured quantity** and K is derived:

    E₀ = n·RT·Γ∞·[θ/(1−θ) − 2aθ²]      (Frumkin generalisation of the Langmuir–Gibbs form)

so, given a literature E₀, θ is fixed and K follows from the isotherm:

    K = θ/(1−θ) · exp(−2aθ)/c

`study_ionic_generality.py` does this inversion and records every intermediate. Literature
anchors: CTAB 30 mN/m (range 20–40), SDS 22 mN/m (range 15–30). Result at 1.5× CMC:
**CTAB θ = 0.724, K = 0.690 m³/mol; SDS θ = 0.643, K = 0.052 m³/mol** — both physical.

## 13.2 The transfer fails on the PAT-1M, and it is not the isotherm

With corrected isotherms, CTAB's second harmonic **never reaches the 0.05 noise floor anywhere
inside the 0.005–1 Hz PAT band**. Best in-band I₂ = 0.0057 at ω = 6.28, i.e. 0.11× noise.
SDS is worse (0.0007, 0.01× noise).

I checked whether that was an artefact of my E₀ choice (`study_tau2_horizon.py`, part 1):
swept E₀ over **20–40 mN/m** (CTAB) and **15–30** (SDS) × Frumkin a = 0/0.35/0.7/1.0 × ω = 0.5/1 Hz
— 40 and 24 feasible combinations respectively. **Not one is usable.** Best achievable CTAB I₂ is
0.0068 at the extreme corner (E₀ = 40, a = 1.0, ω = 6.28). The failure is robust to the isotherm.

The cause is a frequency–timescale mismatch: **τ₂(CTAB) = 65 ms** ⇒ the response turns over near
ωτ₂ ≈ 1, i.e. **≈ 2.4 Hz**, and the PAT stops at 1 Hz. The instrument is simply too slow.

## 13.3 Where CTAB *does* work: it is an instrument problem, not a method problem

Extending to the tensiometer band (corrected isotherm, 1.5× CMC, ε = 0.20):

| f (Hz) | I₁ | I₂ | I₂/noise | lag (°) |
|---|---|---|---|---|
| 1.0 | 0.074 | 0.006 | 0.11 | 55.2 |
| 15.9 | 0.715 | **0.065** | **1.30** | 74.8 |
| 40.0 | 1.623 | 0.116 | 2.32 | 63.2 |
| 100.0 | 2.859 | 0.252 | 5.04 | 42.6 |

CTAB becomes detectable at **≈ 16 Hz** and comfortable at 25–100 Hz — the capillary-pressure
tensiometer band (0.1–500 Hz). SDS does **not** clear noise even at 100 Hz (I₂ = 0.033, 0.66×).

**So "CTAB solvable" is true, and "CTAB solvable on a PAT-1M" is false.** Both statements are now
supported rather than asserted. Two caveats stated plainly: (i) I₁ ≥ 0.7 means the *linear*
response is already ~1 mN/m and weakly nonlinear theory is strained, so the 100 Hz row is
indicative, not quantitative; (ii) the SDS negative at 100 Hz is **band-edge-limited** — the scan
stops at 100 Hz, not the method.

## 13.4 A hypothesised collapse that FAILED

I expected the response to be a function of ωτ₂ alone, which would have unified all three
surfactants into one design curve. Tested at matched ωτ₂ (0.5 … 64), it does **not** collapse:
at ωτ₂ = 4, I₁ = 0.482 (C₁₂E₈, K = 140) vs 0.205 (CTAB) — a factor 2.4 apart.

A second timescale is in play. Computing the diffusion time from the inverse-isotherm slope,
ω_D = D·(dc_s/dΓ)², gives τ_D = 6.3 ms (C₁₂E₈), 0.8 ms (CTAB), 0.02 ms (SDS) — against
τ₂ = 4000/65/38 ms. So **τ_D ≪ τ₂ by 3–4 orders of magnitude everywhere**, which is *why*
ωτ₂ looks like the controlling group locally, but the weak residual dependence on τ_D (through
the surface excess) is enough to break the collapse when E₀ and θ differ as much as they do here.
I am not claiming a resolved scaling law. The honest statement is: **no single-group collapse
found; τ₂ alone does not determine the response, so every parameter set needs its own run.**

## 13.5 SELF-CORRECTION: my own C₁₂E₈ studies used a saturated isotherm

While checking §13.2 I found the same defect I had just condemned in Antigravity's numbers, in
**my own**. Every headline result above (the 1.2–3× CMC window, the 2.5% recovery precision, the
Protocol B correlation reduction) was computed with **K = 140 m³/mol**, which gives
**θ = 0.937 and E₀ = 103.5 mN/m at 1.5× CMC** — the very regime §12 calls unphysical.

I could have quietly swapped in E₀ ≈ 25 and re-reported. That would have been outcome-shopping.
Instead `study_e0_sensitivity.py` re-ran the whole window at four anchor elasticities:

| Anchor E₀ (mN/m) | θ @1.5× | best conc. | Δlag per 4× τ₂ | usable window |
|---|---|---|---|---|
| 25.0 | 0.783 | **1.5×** | 19.4° | 1.2 / 1.5 / 2.0 / 3.0× |
| 45.3 (Czakaj's own) | 0.867 | **1.5×** | 18.8° | 1.2 / 1.5 / 2.0 / 3.0× |
| 73.9 (Czakaj's own) | 0.914 | **1.5×** | 18.0° | 1.5 / 2.0 / 3.0× |
| 103.5 (K = 140) | 0.937 | **1.5×** | 17.4° | 1.2 / 1.5 / 2.0 / 3.0× |

**The result is robust.** The optimum is 1.5× CMC at every anchor; the window edges hold; the
leverage varies by only 11% across a 4× range in E₀. So the conclusions are not artefacts of the
saturated isotherm — which is the outcome that lets me keep them. But the *reported* numbers must
now carry the anchor, and §10's Δlag figures gain roughly 10%.

There is independent support for the high-E₀ end: Czakaj's own Table 4 reports E₀ = **73.9 mN/m at
0.3 mM** and 45.3 at 0.5 mM for LAE. This is genuinely a high-elasticity system; a 20–30 mN/m
assumption for it would have been the error. **The general lesson stands, though: an isotherm
constant is not a free parameter. Anchor it to a measured E₀ or report the regime you are in.**

## 13.6 Correction to my own metric, caught before it reached the manuscript

The field I have been calling "Δlag per doubling of τ₂" is computed as
|lag(2τ₂) − lag(τ₂/2)| — that spans a factor **four**, i.e. two doublings. The labels in
`study_ionic_generality.py`, `study_tau2_horizon.py` and `study_e0_sensitivity.py` are renamed
`dlag_per_4x`. The true per-doubling figure, measured directly on C₁₂E₈ at 1.5× CMC:

| τ₂ span | Δlag |
|---|---|
| 1 → 2 s | 2.9° |
| 2 → 4 s | 7.3° |
| 4 → 8 s | 10.2° |

So the response is **not linear in log τ₂** — it steepens in the middle of the range. Quoting
"13°/e-fold" as a single number would have been wrong in both directions. The correct statement
is the table above, and that the *steepest* leverage (10.2°/doubling) sits near the operating
point τ₂ = 4–8 s.

## 13.7 What this changes for the manuscript

1. The generality claim becomes **conditional and therefore citable**: inversion works for
   C₁₂E₈-class nonionics (τ₂ ~ 1–8 s) on a PAT; it works for CTAB (65 ms) only on a CPT-class
   instrument at ≥ 16 Hz; SDS is unresolved.
2. The manuscript gains a **design rule**: operability requires the drive frequency to be near
   1/τ₂ *and* the instrument must reach it. That rule is what failed for CTAB, not the theory.
3. Every quantitative claim must now state its **E₀ anchor**. The manuscript can no longer say
   "K = 140" without saying "θ = 0.94, E₀ = 103 mN/m, and here is the 25 mN/m row."
4. The C₁₂E₈ concentration series (Protocol B) is still measured at one anchor (K = 140). Re-running
   Protocol B at E₀ = 25 is the remaining arithmetic, and it is the last quantitative gap.

## 13.8 Protocol B at the anchored isotherm: partly overturned

§13.7 item 4 was the last quantitative gap. It is now closed, and the answer is not the one I
reported in §12.

Re-running the joint (τ₂, capacity) fit at three C₁₂E₈ anchors (`study_protocol_b_anchored.py`):

| Anchor E₀ | 1 conc. | 2 conc. | 3 conc. | 4 conc. |
|---|---|---|---|---|
| 25.0 | +0.60 | +0.46 | +0.45 | −0.11 |
| 45.3 | +0.59 | −0.28 | −0.20 | +0.12 |
| 103.5 | −0.60 | +0.40 | +0.25 | +0.47 |
| CTAB 30 | +0.59 | +0.41 | −0.09 | **−0.63** |

(τ₂ is recovered to 0.5–1.6% relSE in every cell; the table is the parameter **correlation**, where
±1 is total degeneracy.)

**What survives:** a single concentration is always degenerate (|corr| ≈ 0.6); adding
concentrations reduces the *magnitude* in most cells; τ₂ itself is always recovered.

**What does not survive — two corrections to §12:**
1. **"Two concentrations suffice" is wrong**, or at least not robust. At E₀ = 25 the second
   concentration only moves |corr| from 0.60 to 0.46, and at E₀ = 103.5 it makes it *worse*
   (0.60 → 0.40 with the sign flipping, which is noise in a correlation this unconstrained).
   Reduction is not monotone in the number of concentrations, and **not monotone in likelihood**.
2. **The clean monotone −0.61 → −0.37 → −0.24 trend in §12 does not reproduce** at any other
   anchor. It was a single-anchor result that I over-read as a protocol property. This is the same
   error as §9 — generalising from one configuration — caught this time before the manuscript.

**Why the correlation is erratic:** with the lag as the only observable and both τ₂ and capacity
acting on it in the same direction, the two-parameter problem is weakly identified by construction.
Multi-concentration data helps, but the improvement is sensitive to where the concentrations sit
relative to the window edges, and 3 points at 3 frequencies is a thin design. The reliable claim is
narrower and more useful: **the phase lag determines a combination of τ₂ and capacity; breaking the
degeneracy needs a concentration series, and even then the correlation is 0.1–0.6, not ~0.**

**Consequence for the manuscript's headline uncertainty.** "τ₂ to ~2.5%" is a *one-parameter* result
(§11, Protocol A: capacity known from CMC and n̄). If capacity must also be fitted, the honest
figure is the relSE column here — **0.5–1.6%** on synthetic noise-free-in-model data with
σ_lag = 0.18°, which will degrade against real instrument noise. Both numbers must be quoted with
their conditioning: **Protocol A gives precision; Protocol B gives only partial separation.**

## 13.9 An open tension in the CTAB anchor that I could NOT resolve numerically

While writing the handoff I checked whether the CTAB anchor in §13.1 is actually consistent with
measured surface excesses, and found a tension I cannot close with this model.

Holding E₀ = 30 mN/m and sweeping the Frumkin interaction parameter a (which trades off against
coverage at fixed E₀):

| a | θ | Γ (µmol/m²) | I₂ @1 Hz | I₂ @16 Hz | I₂ @100 Hz |
|---|---|---|---|---|---|
| 0.7 | 0.724 | 2.32 | 0.0057 | 0.065 | 0.252 |
| 2.0 | 0.821 | 2.63 | 0.0073 | 0.070 | 0.690 |
| 4.0 | 0.892 | 2.85 | **0.0000** | **0.0000** | **0.0000** |
| 7.8 | 0.940 | 3.01 | **0.0000** | **0.0000** | **0.0000** |

Literature CTAB Γ(CMC) is ≈ 2.8–3.2 µmol/m², which needs **a ≈ 7.8** at this Γ∞. But at a ≥ 4 the
model returns an **identically zero** response — γ range = 0.0000 mN/m over the whole cycle, not a
small number but a stalled solution, with **`mass_error = nan`**. The implicit Frumkin solve goes
numerically stiff once the exp(2aθ) term dominates.

**So this is not evidence about CTAB; it is evidence that the model is invalid there.** I am
recording it as an unresolved numerical limit, not a result. Two consequences:

1. The CTAB conclusion in §13.3 is stated for **a = 0.7–2.0**, where the solver is trustworthy.
   Within that range I₂ at 16 Hz moves only 0.065 → 0.070 and at 100 Hz rises 0.25 → 0.69, so the
   *qualitative* conclusion (no PAT, fine on a CPT) is stable. It is **not** verified at the a ≈ 7.8
   that the measured Γ actually implies.
2. `mass_error = nan` is a live defect in `run_oscillation`'s bookkeeping, independent of this
   parameter question. It should be fixed before any mass-balance claim is made in the manuscript.

This is now a concrete data request rather than a modelling question: **a measured CTAB
surface-pressure isotherm (Π vs c) with its Γmax and a fitted independently** would settle whether
a ≈ 0.7 or a ≈ 8 is right, and if the latter, the model needs a stiffer-isotherm formulation
(logarithmic or van der Waals) before CTAB can be claimed at all.
