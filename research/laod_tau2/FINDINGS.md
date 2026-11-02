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

# 9. The LAE peak at the CMC rules out micellar buffering — and with it the paper's premise

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
