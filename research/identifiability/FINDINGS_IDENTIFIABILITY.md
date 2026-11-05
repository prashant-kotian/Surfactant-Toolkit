# Identifiability of relaxation mechanisms in interfacial dilatational rheology

**Status:** exploratory results, round 27. Not a manuscript. Nothing here has been
peer-reviewed, and the word "novel" is deliberately not used until the collision tests in
`NEXT_PAPER_GAP_ANALYSIS.md` §2.4 have been run properly.

**Data:** 6 points of |δ| versus frequency, C14E8 at 30 µmol/L, digitised from Fainerman
et al., *Langmuir* 2008, 24, 6447, Fig. 5. Precision ~1°.

---

## 1. What was done

Four mechanisms were written through a single object, the surface relaxation
admittance B(ω):

| model | B(ω) | parameters |
|---|---|---|
| M1 LVDT | 1 | none beyond E₀ |
| M2 single relaxation | 1 + C_m/(1+iωτ) | C_m, τ |
| M3 two relaxations (50/50) | 1 + C_m/2·[1/(1+iωτ₁) + 1/(1+iωτ₂)] | C_m, τ₁, τ₂ |
| M4 log-normal distribution of τ | ∫ w(z)/(1+iωτ(z)) dz | C_m, τ_c, width |

and the observable is the phase lag |δ| = |arg E| with

    E = E₀ / (1 + √B e^{-iπ/4} / u),    u = h √(ω/D).

This form is algebraically identical to the project's validated transfer function and
reproduces `research/laod_tau2/linear_response.py` to 1e-14 (test 1). Everything below is a
consequence of that form, which is why it was worth deriving rather than fitting numerically.

---

## 2. Result 1 — the 45° ceiling is about *instantaneous equilibrium*, not about diffusion

LVDT is usually quoted as "the dilatational phase of a diffusion-controlled layer cannot
exceed 45°". That is true but understates what is going on. Because E is a function of √B
multiplied by a complex constant, as u → 0

    arg E → arg √B − 45°    ⇒    |δ|_max = 45° + |arg B| / 2.

Verified numerically to 5e-3° across three admittances (test 3). Three consequences:

1. The 45° figure is the value for **B = 1**, i.e. a layer that re-equilibrates
   instantaneously with the subsurface. Diffusion sets the *shape* of |δ|(ω); the
   *ceiling* is set by surface relaxation.
2. A passive layer has Im B ≤ 0 (B is a positive-weighted sum of 1/(1+iωτ)), so
   0 ≤ |arg B| < 45° and **every passive relaxation mechanism is confined to
   45–67.5°**. Checked over C_m up to 60 and τ over 7 decades (test 4).
3. Therefore a measured |δ| **above 67.5° falsifies the entire mechanism class at once** —
   barrier, reorientation, exchange, any distribution of times — and points outside it
   (surface viscoelasticity, a non-passive layer, or failure of planar/local-equilibrium
   assumptions).

Point 3 is the useful one. The existing literature argues one mechanism against another,
which needs the data to be able to tell them apart. This bound is mechanism-agnostic and
needs no fitting: one measurement either violates it or does not.

The factor of ½ is what makes the bound tight. The exponent of B enters the observable
halved, so the phase is a **compressed** image of the surface relaxation spectrum — which is
the theoretical reason resolution is poor, and the reason to look for a better observable
rather than a better fit.

---

## 3. Result 2 — phase-only data cannot identify three of the four parameters

Directly from the form of E:

- **E₀ does not appear in arg E.** Test 7. Phase-only data therefore say nothing about the
  Gibbs elasticity, however precise. Yet phase-only papers report "E".
- **h and D enter only through u = h√(ω/D)**, so they are degenerate together.
- **The relaxation enters only through arg B**, so any two mechanisms with the same arg B
  are indistinguishable in phase.

Confirmed by simulation (`amplitude_test.py`), truth E₀ = 72 mN/m, C_m = 0.8, τ = 12 s,
12 frequencies 0.005–1 Hz, 1° noise:

| observable | E₀ recovered | C_m | τ | D |
|---|---|---|---|---|
| phase only | **349.9** (4.86× true) | 0.800 | 12.00 | — |
| \|E\| + phase | 72.00 | 0.800 | 12.00 | — |

The phase-only fit is not *wrong* — χ² = 0 in both cases. E₀ is simply unconstrained: 350
mN/m is thermodynamically impossible for a nonionic monolayer and the data do not care.

**But the fix is not simply "measure |E|".** With D free as a fifth parameter, both
observables recover only the combination h/√D, to 1.0000, while h and D individually
diverge by 17×:

| observable | h/h_true | D/D_true | (h/√D) ratio |
|---|---|---|---|
| phase | 1.08 | 1.17 | 1.0000 |
| \|E\| + phase | 0.46 | 0.21 | 1.0000 |

This is a **structural** degeneracy, provable analytically (test 10: scaling h by k and D by
k² leaves E bit-identical). No measurement precision removes it. It requires an independent
value of h, D, or Γ∞ — exactly the quantity that is rarely reported, and that surfaces and
interfaces papers routinely treat as known.

---

## 4. Result 3 — the models are not distinguishable on the actual data

AICc comparison on the C14E8 data at the calibrated isotherm slope:

| model | k | χ² | RMS (deg) | AICc | Akaike weight |
|---|---|---|---|---|---|
| **M1 LVDT** | 0 | 229.6 | 6.19 | **21.87** | **0.52** |
| M2 single | 2 | 62.5 | 3.23 | 22.06 | 0.47 |
| M3 two-relaxation | 3 | 62.5 | 3.23 | 32.06 | 0.003 |
| M4 log-normal | 3 | 62.5 | 3.23 | 32.06 | 0.003 |

Two things to read off:

- **M2, M3 and M4 give the same χ² to three figures.** The extra freedom is not used. M3
  converges to τ₁ = τ₂ = 11.5 s and M4 to a width pinned at its lower bound 0.15. A paper
  reporting "two relaxation times, τ₁ = 3 s and τ₂ = 40 s" from data of this quality would be
  reporting a local optimum in an unidentifiable direction.
- **LVDT and the relaxation model are 50/50.** ΔAICc = 0.19. The data cannot choose.

That second line is the one that matters for the field, and it is not flattering. This is
Fainerman's own figure — a paper arguing for a specific mechanism — and six points at 1°
precision cannot separate that mechanism from plain diffusion at conventional thresholds.

**The nuisance dependence is worse.** Repeating the whole comparison across a factor-of-3
uncertainty in h (the isotherm slope is itself estimated, and never reported with an
uncertainty):

| h/h_nominal | best model | χ²(M1) | χ²(M2) |
|---|---|---|---|
| 0.33 | M1 | 230.3 | 145.2 |
| 0.5 | M1 | 221.8 | 118.6 |
| 1.0 | M1 | 229.6 | 62.4 |
| 2.0 | **M2** | 351.0 | 12.2 |
| 3.0 | **M2** | 557.9 | 1.5 |

**The conclusion flips with a nuisance parameter that is not measured.** At the low end a
relaxation model is not just worse than LVDT, it is *itself* excluded (χ² = 145 for 6
points). At the high end it fits to χ² = 1.5. No amount of care in the fitting stage can
repair this: the mechanism assignment is controlled by the isotherm slope, and the isotherm
slope is controlled by the equilibrium isotherm, which is a different experiment.

This is a stronger and more specific statement than "report uncertainties". It says which
unreported quantity decides the answer.

---

## 5. Result 4 — only the *shape* of |δ|(ω) carries mechanism information

Fixed by the E₀ cancellation and the h–D degeneracy: a constant offset in phase is
absorbed by h, and a frequency shift is absorbed by D. What survives is the shape.

The first version of the design map in `design_map.py` ignored this: it measured the raw
excursion above the LVDT baseline and reported large signals for time constants far outside
the measurement band, where B is effectively real, the phase is merely offset, and the effect
is indistinguishable from a different isotherm slope. Profiling out the h-rescaling before
measuring the excursion changes the answer materially — the detectable τ range at C_m = 1 in
the PAT band shrinks from 1e-4–180 s to 0.56–180 s, and the low-C_m cases drop to
non-detection.

After the correction, at 3σ:

| band | C_m = 1 | C_m = 5 | C_m = 20 |
|---|---|---|---|
| microtensiometer (1.7 mHz–1.7 Hz) | τ ≈ 0.32–560 s | 0.06–1000 s | 0.03–1000 s |
| PAT oscillating drop (5 mHz–1 Hz) | 0.56–180 s | 0.10–1000 s | 0.06–1000 s |
| CPT (0.1–500 Hz) | 0.018–10 s | 0.002–32 s | 0.001–180 s |
| PAT + CPT | 0.018–180 s | 0.002–1000 s | 0.001–1000 s |

**C_m ≲ 0.2 is not detectable in any band at 1° noise.** The C14E8 fit gives C_m = 0.77,
which sits just above the threshold — and that fit is, per §4, also consistent with no
relaxation at all. That is not a coincidence, and it is not good news for the published
assignment.

---

## 6. What this says, and what it does not

**Supported:**
- The 45° → 45° + |arg B|/2 generalisation, exact, and the 67.5° class bound that follows.
- The identifiability statements in §3, which are analytic and verified.
- The AICc non-discrimination in §4, for this dataset.
- The h-sensitivity flip in §4, computed, not asserted.
- The phase-excursion map in §5, after profiling the nuisance.

**Not supported / not attempted:**
- Any claim that the field is *universally* wrong. Six points from one figure is one
  system. The result is that this dataset cannot discriminate; other datasets are better
  and some papers do report amplitude.
- Any claim of priority. The 67.5° bound is elementary enough that it is probably in a
  textbook or a review I have not read. It must be checked against Rusanov, Fainerman,
  Kretzschmar, Lucassen-Reynders and the Danov–Kralchevsky reviews before any novelty claim.
  **This is the single largest risk to the work, and it is unresolved.**
- The D substitution: C12E8's diffusion coefficient is used for C14E8. Stated in the code,
  and it affects the frequency scale, not the |arg B| statements.
- `h` is calibrated, not measured. Every number in §4's sensitivity table moves with it.

**The obvious next question**, which the current code cannot answer: is the *joint* design
(frequency band + amplitude + equilibrium isotherm with reported uncertainty) identifiable
where the *standard* design is not? That is answerable with what is already here, by
treating the isotherm parameters as fitted rather than calibrated — and it is the
constructive version of this paper rather than the pessimistic one.

---

## 7. Reproduce

```
python -m research.identifiability.fit_data        # §4  model comparison + profiles
python -m research.identifiability.design_map      # §5  discriminability map
python -m research.identifiability.amplitude_test  # §3  phase vs joint observable
pytest research/identifiability/tests -q           # 11 tests
```
