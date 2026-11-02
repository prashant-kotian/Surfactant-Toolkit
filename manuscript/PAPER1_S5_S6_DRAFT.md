# Paper 1 — Draft of §5 and §6

Companion to `PAPER1_S1_S2_DRAFT.md` and `PAPER1_S3_S4_DRAFT.md`. Numbers from
`recovery_uncertainty.json`, `ionic_generality.json`, `tau2_horizon.json` and the sampling study
recorded in FINDINGS §11.

---

## §5. Sampling Density and Experimental Protocol

### 5.1 The sampling requirement is a bias problem, not a noise problem

It is natural to treat temporal sampling as a signal-to-noise question and apply a Nyquist-style
rule: two samples per cycle at minimum, more is better. For the present purpose that framing is
dangerously wrong, because the quantity being extracted is a **phase**, and an under-sampled phase
estimate is not merely imprecise — it is biased, and the bias does not average out over repeated
measurements.

We quantified this directly. A waveform with a known phase lag of 67.9° was sampled and demodulated
at varying densities, with no added measurement noise, the only error being the sampling itself:

| frames / cycle | recovered δ | scatter |
|---|---|---|
| 1 | 164.3°, 266.2°, 194.1°, 254.1° (four sampling phases) | ≈ 90° |
| 2 | — | 43° |
| 5 | 67.9° | 0.33° |
| 20 | 67.9° | 0.18° |

At one frame per cycle the returned value is not a noisy measurement of 68°; it is a number that
bears no relation to it, and which varies by ±90° depending on where the single sample happens to
fall. We report this because it affected our own analysis: an earlier design curve of ours, built on
a 1 frame/s assumption, recommended an operating point at 1 frame per cycle with an apparent
relative error of 4.5%. That figure was an artefact of favourable sample placement and is withdrawn.

The mechanism is aliasing of the higher harmonics together with the lost phase reference; the
fundamental is not the only component present (§4.5), so a two-point demodulation is not merely
under-determined, it is inconsistent with the waveform it is meant to describe.

### 5.2 The criterion

**At least five frames per cycle.** At five frames the phase is recovered to 0.33° and the estimate
is unbiased; at twenty it improves only to 0.18°, so the returns saturate quickly once the
catastrophic regime is escaped. We recommend 5 as the specification and 10–25 as the practical
target, since the criterion must hold at the *highest* frequency in a sweep.

### 5.3 The PAT-1M satisfies this, and the confusion that says otherwise

The Sinterface PAT-1M in oscillatory mode acquires at 25–50 frames per second
[CIT: instrument, to be sourced]. At 1.0 Hz — the top of its oscillatory band and the operating point
of §3 — that is **25–50 frames per cycle**, comfortably above the criterion. At 0.2 Hz it is 125–250.

A frequently cited figure of approximately one profile per second refers to the instrument's *static*
drop-shape-tensiometry mode, in which a full profile is fitted to obtain a single equilibrium
tension. The two modes are reported in the same methods paragraphs and the numbers are easily
conflated; we flag it because reading the 1 s figure as a sampling rate would imply that the PAT-1M
violates the criterion by a factor of 25, and would wrongly rule out the only widely available
instrument that satisfies it.

### 5.4 Strain amplitude: why ε = 0.20

The second harmonic scales as the square of the strain, which we verified rather than assumed. Over
ε = 0.03–0.30 at 1.5 × CMC and ω = 6.28 rad/s, the local logarithmic slope of I₂ against ε is
2.00, 1.99, 1.95 and 1.84 at ε = 0.05, 0.10, 0.20 and 0.30 respectively — quadratic at small
amplitude with mild saturation at the top of the range. Going from Czakaj's ε ≈ 0.065 to ε = 0.20
raises I₂ from 0.0198 to 0.1827 mN/m, a factor of **9.21** against the 9.5 predicted by pure ε².

That factor is what lifts the second harmonic clear of the 0.05 mN/m noise floor. At ε = 0.065 the
C₁₂E₈ operating point gives I₂ = 0.0198, which is *below* the noise floor; the measurement would
fail the detectability criterion outright even though the underlying physics is unchanged.

Crucially, δ itself is unaffected: over the whole ε = 0.03–0.30 range the loss angle moves by
0.17° (67.97° to 67.80°). **The strain amplitude buys detectability and costs nothing in
invertibility.** An earlier draft of ours claimed the opposite — that operating at large amplitude
traded linear-theory validity against τ₂ sensitivity — and that was wrong (§4.5, corrected). The
loss angle here is an operational observable calibrated by the forward model, not a linear
material property, and it is stable across a sixfold amplitude range.

### 5.5 The protocol

For a new nonionic surfactant of unknown τ₂:

1. **Measure the CMC.** Do not take it from a table. The window of §3 is expressed in units of CMC
   and its lower edge is sharp; a 10% CMC error displaces the whole window by 10%.
2. **Determine n̄** (aggregation number) independently, to fix the capacity closure
   C_m = (c_tot − CMC)·n̄/CMC. This is the only external input the inversion requires.
3. **Set c = 1.5 × CMC.** This is the interior optimum of §3.2, invariant across a fourfold range in
   E₀.
4. **Set ε = 0.20**, subject to §5.6.
5. **Sweep frequency** across the instrument band, verifying that frames/cycle ≥ 5 at every point.
   Retain ≥ 12 cycles and analyse the last four.
6. **Estimate δ at each frequency** by least-squares Fourier decomposition over an integer number of
   cycles, then invert τ₂ by Protocol A (§7) using the Δ-method or a joint fit.

### 5.6 Precision budget, and what it costs

At the operating point (1.5 × CMC, ω = 6.28 rad/s, σ_γ = 0.05 mN/m, 60 Monte-Carlo realisations):

| step | value |
|---|---|
| Monte-Carlo scatter of recovered δ | **0.033°** |
| Monte-Carlo scatter at ω = 3.0 (for comparison) | 0.065° |
| Δδ/Δln τ₂ at 6.28 rad/s | 13.3° per e-fold |
| relSE by Δ-method | **0.25%** |
| relSE by joint 1-parameter fit | **0.23%** |
| relSE, 2-parameter fit (capacity also free) | 3.0%, with corr(τ₂, C_m) = 0.997 |

**We quote τ₂ to 0.2–0.3% under Protocol A**, where capacity is fixed by step 2 and only τ₂ is
estimated. The improvement from 0.065° at ω = 3.0 to 0.033° at ω = 6.28 is a factor of 2 for a
factor of 2 in frequency, which is why the top of the instrument band is the operating point.

Under Protocol B, where capacity is also free, the estimate degrades to 3% with a correlation of
0.997 between the two parameters — i.e. essentially perfect degeneracy at a single concentration.
No amount of averaging fixes that; only additional concentrations or an independent capacity do, and
§7 reports how partially the former succeeds.

**These are model-based figures with realistic amplitude noise.** They exclude error in CMC, n̄, and
E₀, all of which enter the forward model and none of which is a measurement uncertainty in the
inversion. Adding a 10% error in n̄ would dominate the budget; we have not propagated it, and it
should be treated as an open item rather than assumed small.

### 5.7 Honest hardware caveats

- **ε = 0.20 is three times the amplitude of the published LAOD protocol it replaces.** We have not
  verified drop stability at this amplitude on a 10 µL pendant drop with a 2 mm capillary; we have
  only the report that the instrument accommodates ε = 0.15–0.30 on standard capillaries
  [CIT: to be sourced]. **This is a required bench validation and should be the first experiment
  performed.** Nothing in the method survives if the drop detaches.
- **25–50 fps is quoted for oscillatory mode.** Should a particular instrument/software combination
  acquire more slowly, §5.2 is the binding constraint and the method fails at the top of the band
  first.
- **The 0.05 mN/m noise floor is our estimate**, propagated from a 0.1 mN/m instrument precision on
  a 2–4 mN/m fundamental. A better calibration of σ_γ moves both the window's upper edge and the
  achievable precision.

---

## §6. The Timescale-to-Bandwidth Horizon

### 6.1 The design rule

Sections 3–5 establish that the method requires a reservoir that responds on the timescale of the
drive. The relevant dimensionless group is Ω = ωτ₂, and operability requires the drive frequency to
sit in the band **Ω ≳ 1**, with the sensitivity maximum near Ω ≈ 3 (observed, §4.2) and useful
response out to Ω ≈ 30.

Because τ₂ differs by three orders of magnitude across common surfactants, this rule translates
directly into an instrument requirement. For each system, the operable frequency band is

    f ≈ Ω / (2π τ₂),      Ω ∈ [1, 30],                                    (13)

and the question of whether a given tensiometer can perform the inversion becomes a comparison of
(13) with its bandwidth:

| surfactant | τ₂ | operable band, eq. (13) | PAT-1M (0.005–1 Hz) | CPT (0.1–500 Hz) |
|---|---|---|---|---|
| C₁₂E₈ | 4.0 s | **0.04–1.19 Hz** | **covers it** (and 1.19 Hz exceeds the 1.0 Hz top of band by a factor of 1.2 — the optimum at 1.0 Hz is comfortably inside) | covers it |
| CTAB | 65 ms | **2.4–73 Hz** | **cannot reach** | **covers it** |
| SDS | 38 ms | **4.2–126 Hz** | **cannot reach** | **covers most** |

This is the paper's predictive content outside its own system, and it is a statement about
*instruments*, not about surfactants: the method is not limited to slow nonionics, it is limited to
surfactants whose relaxation time the available bandwidth can straddle.

### 6.2 Why the PAT-1M fails for the fast ionics — the mechanism

It is not that CTAB lacks a micellar reservoir. It has one, with a capacity comparable to C₁₂E₈ at
the same reduced concentration. The failure is that at 1 Hz the drive is **far below** the
relaxation rate: Ω = 2π × 1 × 0.065 = 0.41. The reservoir then tracks the area change almost
quasi-statically, the interface stays near compositional equilibrium, and — exactly as in the
over-buffered limit of §3.1 — the response goes flat.

The numbers show the suppression clearly. At 1.5 × CMC with E₀ anchored to 30 mN/m, CTAB's
effective modulus I₁/ε at 1 Hz is **0.37 mN/m**, against an insoluble-limit ceiling of 30. The
interface is 98% buffered and has almost nothing left to say. Its second harmonic is 0.0057 mN/m,
one-ninth of the noise floor. Sweep upward and the picture inverts:

| f (Hz) | Ω = ωτ₂ | I₁/ε (mN/m) | I₂ (mN/m) | δ (°) |
|---|---|---|---|---|
| 1.0 | 0.41 | 0.37 | 0.006 | 55.2 |
| 3.98 | 1.6 | 0.98 | 0.019 | 71.3 |
| 15.9 | 6.5 | 2.60 | **0.065** | 74.8 |
| 100 | 41 | 14.3 | 0.252 | 42.6 |

At Ω ≈ 6.5 the second harmonic reaches the noise floor, and the effective modulus has recovered by
a factor of 7. **CTAB would be measurable on a capillary-pressure tensiometer and is not measurable
on a PAT-1M.** That is the whole of the horizon result, and it is a one-line instrument specification
rather than a limitation of the theory.

### 6.3 A historical consequence, offered as a prediction

There is a well-known asymmetry in the LAOD literature: micellar exchange is reported for slow
nonionic systems far more often than for fast ionic ones. The usual explanation is experimental
difficulty with ionic systems — foam, evaporation, surface charge.

Equation (13) offers a second and more specific explanation. A drop-shape tensiometer is intrinsically
a low-frequency instrument (0.005–1 Hz), and a surfactant with τ₂ below ~0.1 s requires drive
frequencies above its reach. **Low-frequency drop-shape analysis is structurally incapable of
observing micellar exchange in fast surfactants**, independently of any chemical difficulty. We offer
this as a testable prediction rather than a claim about any particular study: the prediction is that
a fast ionic surfactant measured on a CPT at 10–100 Hz will show δ > 45°, whereas the same system on
a PAT-1M will not.

### 6.4 The 45° ceiling as an instrument-independent screen

The diagnostic of §4.1 gives §6 a second, cheaper application that requires no model at all. Because
no diffusion-only monolayer can exceed δ = 45° at any frequency, a measured δ above 45° proves
kinetic relaxation is present. The converse is weaker but useful: **a measured δ below 45° does not
prove the absence of a reservoir, only that the drive frequency was too low to see it.**

This asymmetry is worth stating carefully, because it is the trap §6.2 describes. At 1 Hz, CTAB
sits at δ = 55.2° — above the ceiling, so a reservoir would be detected — yet its I₂ is nine times
below the noise floor, so τ₂ could *not* be inverted. The ceiling tells you a reservoir is there;
it does not tell you the measurement carries enough information to quantify it. Detectability and
invertibility again diverge, and this time across instruments rather than across concentrations.

We checked that the ceiling is not an artefact of our own model by computing the LVDT loss angle
from Czakaj's own fitted (E₀, ν_D) pairs for the LAE system at 0.2–1.0 mM: every value falls below
33° at 0.01–1 Hz, as it must, since the LVDT form is diffusion-only by construction.

### 6.5 Limits of the horizon argument

- **Ω ≈ 3 is observed, not derived.** We did not solve for the optimum of the coupled
  surface/bulk problem, and we do not have a scaling law: at matched Ω, CTAB and C₁₂E₈ differ in
  fundamental amplitude by a factor of 2.4 (FINDINGS §13.4). Equation (13) is a design heuristic
  located by direct computation on three systems, not a collapse. Anyone applying it to a fourth
  surfactant should verify by running the forward model first.
- **The CTAB and SDS parameters are anchored, not measured.** E₀ = 30 and 22 mN/m were chosen from
  literature ranges, and the resulting coverages (θ = 0.72, 0.64) are lower than the saturated
  monolayers that measured surface excesses imply. We established (FINDINGS §15) that the anchored
  isotherm is internally consistent, but it is not independently validated, and CTAB coverage is an
  open question requiring a measured Π–c isotherm.
- **SDS does not clear the noise floor anywhere in the scanned range**, reaching only 0.0048 mN/m at
  10 Hz. That is a band-edge result — our scan stops at 100 Hz, where Ω ≈ 24, not a demonstration
  that SDS fails.
- **The PAT-1M and CPT bandwidths are taken from the instrument specifications** supplied by our
  collaborators; we have not operated either instrument.

---

## Items requiring your decision before §7

1. **§6 is a heuristic, not a law, and I have said so.** Equation (13) rests on three systems and
   has no scaling collapse behind it. If you want a stronger claim in the abstract, it will not be
   supportable.
2. **I need the instrument citation.** "Sinterface PAT-1M, 25–50 fps oscillatory, 0.005–1.0 Hz" and
   "CPT 0.1–500 Hz" are currently sourced to internal collaborator documents, not to a publication.
   §5.3 and §6.1 both depend on them. Please supply a citable reference or confirm that these
   should be described as manufacturer specifications.
3. **The ε = 0.20 drop-stability claim is hearsay in our draft.** You report that the PAT-1M
   accommodates ε = 0.15–0.30 without detachment on 1–2 mm capillaries. I have attributed it to you
   rather than to a source. This is the one claim in §5 that a reviewer will demand a figure for,
   and it is also the one that, if false, invalidates the protocol.
4. **A discrepancy I cannot resolve from here.** The horizon table puts CTAB's operable band at
   2.4–73 Hz, so a CPT at 10 Hz (Ω = 4.1) should clear the noise floor. Our direct scan does not
   reach the floor until ~16 Hz (Ω = 6.5). The design rule is therefore conservative by roughly a
   factor of two in frequency. I have written the rule as Ω ≳ 1 but quoted the measured threshold
   separately rather than force them to agree.
5. **§7 needs the Protocol A/B conditioning restated**, and I should warn in advance that Protocol B
   is a partial negative (FINDINGS §13.8): the τ₂–capacity correlation is 0.1–0.6 and is *not*
   monotone in the number of concentrations. It cannot be presented as a clean recipe.
