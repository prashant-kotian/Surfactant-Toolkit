# Paper 1 — Draft of §3 and §4

Companion to `PAPER1_S1_S2_DRAFT.md`. §1–§2 stand as revised (Langmuir title, loss angle δ pinned).
Every number below is drawn from `research/laod_tau2/e0_sensitivity.json`,
`ionic_generality.json`, `recovery_uncertainty.json` and the §11 sampling study; none is estimated.

---

## §3. The Concentration Inversion Window

### 3.1 Why a window must exist, and why it is not where one expects

Above the CMC the micellar reservoir grows with total concentration. The intuitive expectation is
monotone improvement: more reservoir, more buffering, stronger micellar signature. Both edges of the
observed window contradict that expectation, and the reason is that capacity enters the problem
twice with opposite effects.

- **It buffers the subsurface.** A larger C_m holds c_s closer to equilibrium during the oscillation
  cycle, which *is* the micellar effect and should grow with concentration.
- **It suppresses the nonlinearity that generates I₂.** The second harmonic arises from the
  curvature of the isotherm traversed during the cycle. By holding the interface near equilibrium,
  a large reservoir makes the trajectory quasi-linear, and I₂ collapses faster than the phase effect
  grows.

The first effect sets the lower edge, the second sets the upper edge. Consequently the method lives
in a *bounded* interval, and — the point that matters for replication — **the informative region is
not "near the CMC" in the loose sense** in which that phrase is used to motivate large-amplitude
experiments on micellar systems.

### 3.2 The window, measured

We fix τ₂ = 4.0 s (the C₁₂E₈ benchmark value), drive at ε = 0.20, and sweep concentration at
ω = 6.28 rad/s (1.0 Hz, the top of the PAT-1M oscillatory band). Two criteria are applied
simultaneously, as demanded by §1.3: I₂ above the 0.05 noise floor (detectability), and a change in
loss angle of at least 3° per factor-of-four in τ₂ (invertibility).

| c/CMC | θ | I₂ | δ (°) | Δδ per 4× τ₂ | detectable | verdict |
|---|---|---|---|---|---|---|
| 1.0 | 0.909 | 0.129 | 32.4 | **0.0** | yes | reservoir empty — **τ₂-invisible** |
| 1.2 | 0.923 | 0.052 | 50.4 | 17.4 | yes | usable |
| **1.5** | 0.937 | 0.183 | **67.9** | **17.4** | yes | **optimum** |
| 2.0 | 0.952 | 0.152 | 78.7 | 10.4 | yes | usable |
| 3.0 | 0.968 | 0.080 | 84.3 | 4.4 | yes | marginal upper edge |
| 5.0 | 0.980 | 0.035 | 86.8 | 1.2 | **no** | below noise floor |
| 10.0 | 0.990 | 0.012 | 88.0 | 0.4 | **no** | saturated |

**The window is 1.2–3.0 × CMC with an interior optimum at 1.5 × CMC.** It is not a plateau: Δδ varies
by a factor of 4 across the usable interval, so operating-point choice changes the achievable
precision by the same factor. The lower edge is sharp — at exactly 1.0 × CMC, Δδ = 0.00° while I₂
remains large (0.129), a textbook instance of the detectability/invertibility distinction of §4: the
nonlinearity is strongest where the reservoir is absent.

The upper edge is a genuine cliff, not a gentle roll-off: between 3.0 and 5.0 × CMC the second
harmonic falls from 0.080 to 0.035, crossing the noise floor, and Δδ falls from 4.4° to 1.2°.

### 3.3 The window is a property of the transport problem, not of the isotherm

The results of §3.2 are computed at a specific isotherm parameterisation (K = 140 m³/mol, which at
1.5 × CMC gives θ = 0.937 and E₀ = 103.5 mN/m). Since E₀ is itself uncertain by a factor of several
and §2.4 warns that an assumed E₀ pins the operating coverage, this is a potential confound. We
therefore re-ran the entire sweep at four anchor elasticities spanning 25–103.5 mN/m, each
corresponding to a different K and different coverage:

| Anchor E₀ (mN/m) | θ at 1.5× | K (m³/mol) | optimum | Δδ there (°) | usable range |
|---|---|---|---|---|---|
| 25.0 | 0.783 | 33.8 | **1.5×** | 19.4 | 1.2–3.0× |
| 45.3 | 0.867 | 61.3 | **1.5×** | 18.8 | 1.2–3.0× |
| 73.9 | 0.914 | 100.0 | **1.5×** | 18.0 | 1.5–3.0× |
| 103.5 | 0.937 | 140.0 | **1.5×** | 17.4 | 1.2–3.0× |

**The optimum is 1.5 × CMC at every anchor, and the leverage varies by only 11% across a fourfold
range in E₀.** We regard this as the single most important robustness result in the paper: the
window is set by the competition between reservoir buffering and isotherm curvature (§3.1), both of
which are properties of the transport problem, and neither of which depends on the particular
(K, Γ∞) pair used to place the operating point on the isotherm.

We state the corollary plainly, because it constrains how the result may be cited: **the window
edges are robust, the absolute numbers are not.** A paper quoting "Δδ = 17.4° at 1.5 × CMC" must
state its E₀ anchor. We therefore quote Δδ with its anchor throughout, and report the range
17.4–19.4° rather than a single figure.

### 3.4 What the window is not

- It is **not** a statement about where micellar exchange is physically strongest. Exchange is
  strongest as C_m → ∞. The window is where exchange is most *measurable*, which is a different and
  much narrower question.
- It is **not** transferable between surfactants without recomputation. §6 shows that the equivalent
  window for CTAB sits at frequencies the PAT-1M cannot reach; the *shape* of the argument
  transfers, the *numbers* do not.
- The upper edge depends on the noise floor, and therefore on instrument and averaging. With
  σ_γ = 0.05 mN/m we find the floor crossed between 3 and 5 × CMC. A tenfold improvement in noise
  would move that edge, and we have not computed where.

### 3.5 A cautionary note on the CMC itself

The window is expressed in units of CMC, so its location inherits any error in the CMC. For C₁₂E₈ we
use 0.071 mM. Because the lower edge is sharp (§3.2), a 10% error in the CMC displaces the entire
window by 10% and could place the nominal "1.5 × CMC" operating point outside the usable interval
at the low end. **The CMC must be measured, not taken from a table, before the protocol is applied
to a new system.** This is an unglamorous point that we expect to be the most common cause of
failure to reproduce the method.

---

## §4. Observables: Detectability versus Invertibility

### 4.1 The loss angle and the 45° ceiling

Section 2.6 establishes that δ is the lag of γ behind A, and that the tabulated LVDT modulus gives
δ = arctan[ξ/(1+ξ)] with ξ = √(ν_D/2ν), so that **any diffusion-only monolayer has δ < 45° at every
frequency**, approaching 45° as ω → 0 and 0° as ω → ∞.

This bound is the organising fact for the whole paper. It is instrument-independent, it requires no
fit, and it does not depend on E₀. A solution containing micelles routinely exceeds it. At
1.5 × CMC and ω = 6.28 rad/s we obtain δ = 67.9° where the same system with C_m = 0 gives 37.9°.

**A measured δ above 45° is therefore categorical evidence of a micellar (or other) reservoir**, and
no choice of D, E₀ or ν_D can produce it. The magnitude of the excess is what carries τ₂.

### 4.2 The reservoir raises the loss angle — and this is the opposite of the intuitive expectation

We initially expected the reservoir to make the interface *more elastic*, damping δ downwards from
its diffusion-limited value, on the reasoning that rapid monomer replenishment holds the subsurface
concentration constant and lets the surface excess follow the area more closely. Direct computation
shows the opposite, and the sign matters enough to state explicitly:

| f (Hz) | δ, C_m = 0 | δ, C_m = 47.5 | shift |
|---|---|---|---|
| 0.010 | 44.28 | 51.81 | **+7.5** |
| 0.050 | 43.40 | 69.62 | **+26.2** |
| 0.200 | 41.83 | 79.98 | **+38.2** |
| 0.500 | 39.98 | 76.86 | **+36.9** |
| 1.000 | 37.91 | 67.91 | **+30.0** |

The shift is **positive and large** — up to +38° — and is not monotone in frequency: at fixed
τ₂ = 4 s it peaks near 0.2 Hz and decreases on both sides. Nor is it monotone in τ₂: at
ω = 6.28 rad/s the shift is +14.4° at τ₂ = 0.05 s, rises to +38.8° at τ₂ = 0.5 s, and falls again to
+8.7° at τ₂ = 20 s.

The maximum therefore sits at intermediate ωτ₂ of order a few — approximately 3 on the available
scan, not at ωτ₂ = 1 as a single-relaxation-time argument would suggest. We note this without
explaining it: it is a property of the coupled surface/bulk problem rather than of the exchange
term alone, and we did not attempt to derive the optimum. What matters operationally is that the
sensitivity has a maximum in the τ₂ direction as well as in concentration, so a system whose τ₂ lies
far from the band centre is harder to invert for a second, independent reason.

The physical reading is that the reservoir is not primarily an elastic element. It is a **lagging
source term**: the micelles deliver surfactant *late* in the cycle, when the subsurface has been
depleted, which puts the γ response further behind the area drive and adds dissipation rather than
restoring elastic energy. Modelling it as an added elasticity gives the wrong sign for the effect.

This also resolves the ambiguity in the phrase "phase depression" that appears in the LAOD
literature: the *reservoir* raises δ; what *lowers* δ is the approach to the insoluble limit at high
frequency. The two effects act in the same frequency range and are easily confused.

### 4.3 Detectability: I₂

I₂ measures the departure from linearity and sets the noise-limited reach of the method. It is
driven by isotherm curvature and, as §3.1 shows, is *suppressed* by the reservoir. Its value as a
criterion is thus asymmetric: a small I₂ rules the measurement out, but a large I₂ does not rule it
in.

### 4.4 Invertibility: what actually carries τ₂

Three candidate observables, evaluated at the 1.5 × CMC optimum, as τ₂ is varied over a factor of
16:

| τ₂ (s) | I₁ | I₂ | **I₂/I₁** | **δ (°)** |
|---|---|---|---|---|
| 1.0 | 1.485 | 0.131 | 0.0882 | 78.12 |
| 2.0 | 2.057 | 0.168 | 0.0816 | 75.19 |
| 4.0 | 2.734 | 0.183 | 0.0669 | 67.92 |
| 8.0 | 3.299 | 0.157 | 0.0475 | 57.75 |
| 16.0 | 3.563 | 0.112 | 0.0315 | 48.78 |

**I₂/I₁ is not blind, but it is not resolvable.** It changes by 0.007–0.019 per doubling of τ₂. The
noise floor on the ratio, propagated from a 0.1 mN/m instrument reading on a 2–4 mN/m fundamental,
is 0.025–0.05. The τ₂-induced change is therefore *smaller than the measurement uncertainty*, and a
single-point inversion from the ratio cannot work. Any published inversion based on I₂/I₁ alone
should be regarded as under-determined.

**δ carries the signal.** It moves 7.3–10.2° per doubling of τ₂ in the region τ₂ = 2–8 s, against a
phase uncertainty of 0.18° achievable at adequate sampling (§5). The ratio of sensitivity to
uncertainty is 40–57, versus 0.2–0.8 for the harmonic ratio — a separation of roughly two orders of
magnitude. That is the quantitative content of the claim that the phase lag is the operative
observable.

**δ is not log-linear in τ₂.** The response is 2.9° for 1→2 s, 7.3° for 2→4 s and 10.2° for 4→8 s,
and flattens again beyond 8 s. Quoting a single "degrees per e-fold" would be wrong in both
directions; the local slope near the operating point is the correct figure.

### 4.5 The high-harmonic caveat

At ε = 0.20 the fundamental response is strongly nonlinear (I₁ = 1.5–3.6 mN/m), and the waveform
carries substantial power above the second harmonic. Two consequences:

1. δ, extracted from the fundamental only, is a *projection* of a nonlinear waveform, not a loss
   angle in the linear-response sense. Its τ₂ sensitivity is real and monotone (Table 4.4), but it
   should not be interpreted as the phase of a linear transfer function. We report it as an
   operational observable whose calibration is the model, not as a material property.
2. The strong higher-harmonic content means the observation is genuinely in the large-amplitude
   regime, which is the point — but it also means the analysis cannot be validated against linear
   LVDT theory beyond the C_m = 0 limit of §2.6, where we obtained 0.48° agreement.

An experimenter wishing to use δ as a *material* quantity would reduce ε until the harmonics are
negligible. We tested what that costs and did not find a penalty in δ: over ε = 0.03–0.30 the loss
angle is constant to 0.17° (67.97° at ε = 0.03, 67.92° at ε = 0.20). The reason is that in this
system the nonlinearity is small in the fundamental (I₂/I₁ ≈ 0.07) and the C_m ≠ 0 response is
already far from linear theory at *every* ε we can measure; the curved isotherm trajectory that
produces I₂ is not what produces the phase lag.

The trade-off that does bind is therefore **not** δ versus linearity, and we correct an earlier
draft that claimed it was. The binding constraint is signal-to-noise on I₂ for the *detectability*
criterion: I₂ scales as ε² while δ does not depend on ε, so ε is set by how far I₂ must be lifted
above the noise floor, and δ is then estimated from a waveform that happens to be nonlinear. This is
a more favourable situation than we first supposed, and §5 quantifies it.

### 4.6 Summary of the two criteria

| | Detectable | Invertible |
|---|---|---|
| Observable | I₂ | δ |
| Noise-limited by | σ_γ on the fundamental | sampling density, then σ_γ |
| Behaviour vs. reservoir | **decreases** with C_m | **increases** with C_m |
| Behaviour at 1.0 × CMC | maximal (I₂ = 0.129) | **zero** (Δδ = 0.00°) |
| Behaviour at 10 × CMC | below noise | saturated, τ₂-insensitive |
| Sensitivity vs. uncertainty | 0.2–0.8 | 40–57 |

The two criteria peak at *different* concentrations: detectability at the CMC, invertibility at
1.5 × CMC. The operating point must be chosen for invertibility, which means deliberately operating
where I₂ is not maximal. This is the practical content of §1.3, and it is why a "maximise the
nonlinear signal" design rule would select the one concentration at which the method cannot work.

---

## Items needing your decision before §5

1. **The "phase depression" narrative in your resolution was backwards, and I have reversed it.**
   You proposed that micellar exchange damps δ from 32.4° to 18.8°, and that this depression is the
   primary signature. Computation gives the opposite: the reservoir *raises* δ (37.9° → 67.9° at
   1.5 × CMC), and it exceeds the 45° diffusion ceiling. I have kept your *structure* — that the
   loss angle is the signature — and replaced the sign. Please confirm, since §6 will reuse this.
2. **The 18.8° vs 32.4° discrepancy was mine, not a convention conflict.** It came from comparing
   against a ν_D borrowed from the LAE system. With ν_D computed self-consistently the agreement is
   0.48°. No offset is applied anywhere. This removes the open item flagged in the §1–§2 draft.
3. **The LVDT ceiling is 45°, not 35°.** The arctan(1/√2) ≈ 35.3° figure I put in the earlier draft
   was wrong and has been deleted; δ = arctan[ξ/(1+ξ)] → 45°.
4. **Citations.** Thank you — those five are exactly what §1 needed. One caution before I insert
   them: I have not read Stadler 2016 or Risse 2025, so I will cite them only for claims I can
   verify from your description (instrument capability and review-level context respectively), and
   I will not attribute specific numbers to them. The Czakaj reference I have data for, but note the
   dataset in this repository is `Molecules` 26:5894 (2021), not `Colloids Surf. A` 627:127181.
   **Please confirm which Czakaj 2021 paper the harmonic data in §1 is attributed to** — I have the
   `Molecules` data extracted and the two may not be the same paper.
5. **ε = 0.20 throughout.** Czakaj's LAE measurements used ε ≈ 0.065. Our C₁₂E₈ numbers are at 0.20.
   §5 should state this mismatch explicitly, because it means the method as specified requires a
   strain amplitude three times larger than the published LAOD protocol it is meant to replace.
