# Direct Inversion of Micellar Relaxation Kinetics from Large-Amplitude Oscillatory Dilatational Surface Rheology: Theory, Bandwidth Limits, and Robust Parameter Identification for Model Nonionic Surfactants

**Paper 1 of 2 — the C₁₂E₈ foundation paper**

*Unified manuscript draft, assembled from `manuscript/PAPER1_*.md` at commit `dd48f4c` plus the
citations resolved in round 5. Body text is reproduced verbatim from the section drafts; no number
has been retyped. Editorial notes are confined to the Appendix and are not part of the paper.*

---

## Scope statement for co-authors

This paper is restricted to the nonionic surfactant C₁₂E₈ on a profile-analysis tensiometer
(PAT-1M, 0.005–1.0 Hz) with a Langmuir surface equation of state (a = 0). Every quantitative claim
in §3–§7 is computed for that system. The ionic surfactants (CTAB, SDS) appear only in §6, as a
predictive instrument-selection criterion derived from their relaxation timescales — not as
experimental validations. Multicomponent and Frumkin-interaction effects are deferred to Paper 2.

**Principal limitation, stated up front:** this is a methods paper validated against its own forward
model. It contains no measured C₁₂E₈ dataset, and the τ₂ = 4 s benchmark is a literature value
rather than one inverted from data. See §7.7.

---

## Abstract

The kinetics of micellar exchange at a surfactant-laden interface are carried by the slow
Aniansson–Wall relaxation time τ₂, a parameter that governs the response of foams and emulsions to
any non-steady area change. It is nevertheless almost invisible to conventional dilatational
rheology: in the linear regime the micellar reservoir competes with diffusion, and the two
contributions cannot be separated from a single frequency response. Here we show that the phase lag
of the surface-tension response in **large-amplitude** oscillatory dilatational rheology carries τ₂
quantitatively, and we establish the conditions under which that information can be recovered.

Using a mass-conserving reaction–diffusion model in which the micellar reservoir is an explicit
species with linearised Aniansson–Wall exchange, we identify a categorical signature that requires no
model: no diffusion-only monolayer can exceed a loss angle of 45° at any frequency, whereas a
micellar solution routinely reaches 68–80°. Because the reservoir acts as a *lagging source term*
that delivers surfactant late in the expansion cycle, it adds dissipation rather than elasticity — the
opposite of the intuitive expectation — and this sign is diagnostic.

We then separate **detectability** from **invertibility** and show they peak at different
concentrations. The second harmonic I₂ is maximal at the CMC, where the micellar reservoir is empty
and its sensitivity to τ₂ is exactly zero. The operative observable is instead the loss angle,
whose sensitivity exceeds its measurement uncertainty by two orders of magnitude where the harmonic
ratio's does not: across a doubling of τ₂ the ratio moves 0.007–0.019 against a noise floor of
0.025–0.05. For the model nonionic surfactant C₁₂E₈ on a profile-analysis tensiometer we measure a
concentration inversion window of **1.2–3.0 × CMC with an interior optimum at 1.5 × CMC**, and show
this window is invariant across a fourfold range in the assumed limiting elasticity — it is a
property of the transport problem, not of the isotherm parameterisation. A sampling density of
**≥ 5 frames per cycle** is required, below which the phase estimate is biased rather than noisy
(164–266° against a true 68° at one frame per cycle). At the operating point τ₂ is recovered with a
**statistical precision of 0.2–0.3%** for a known aggregation number.

Because the method requires the drive frequency to straddle the relaxation rate (ωτ₂ ≳ 1), it
converts directly into an instrument-selection criterion: the 0.005–1 Hz band of a drop-shape
tensiometer covers C₁₂E₈ (τ₂ = 4 s) but cannot reach CTAB (65 ms) or SDS (38 ms), which require a
capillary-pressure tensiometer above ≈ 16 Hz. This explains, without invoking chemical difficulty,
why low-frequency drop-shape analysis has not reported micellar exchange in fast ionic surfactants.

**The precision figure is statistical.** Accuracy is limited by systematic error in the aggregation
number, which propagates one-for-one: a 10% error in n̄ produces a 10% bias in τ₂, forty-four times
the statistical uncertainty. Both numbers must be quoted together, and we do so throughout.

---

## §1. Introduction

### 1.1 The measurement, and the parameter buried inside it

When the area of a surfactant-laden interface is driven sinusoidally, the surface tension responds
at the same frequency but out of phase. For a surfactant that is *only* transported by diffusion,
that phase lag is a nuisance: it is set by the frequency and by a single diffusion time, and the
standard Lucassen–van den Tempel (LVDT) modulus accounts for it completely [CIT: Lucassen & van den
Tempel]. The practical consequence is that decades of dilatational rheology have been built on a
one-parameter-per-frequency description, fitted to a diffusion coefficient, and validated by its
ability to collapse data across frequency.

Above the critical micelle concentration (CMC) this description is incomplete in a way that matters.
Micelles constitute a reservoir of surfactant that can dissociate and feed the interface within a
single oscillation cycle. The kinetics of that exchange are characterised by the slow relaxation
time τ₂ of the Aniansson–Wall stepwise association scheme [CIT: Aniansson & Wall]. τ₂ is not a
nuisance parameter: it controls how quickly a micellar solution can heal a depletion layer, it
governs foam and emulsion stability under dynamic area change, and it is the parameter that
distinguishes a "fast" from a "slow" surfactant in any application where the interfacial area is not
constant.

The difficulty is that τ₂ is essentially invisible to linear methods. In the linear regime the
micellar reservoir enters as a second relaxation process with a well-defined signature, but it does
so in competition with diffusion, and the two contributions are not separately identifiable from a
single frequency response. The classical resolution is to measure over a wide frequency range and
look for the resolvable separation of timescales — which requires an instrument spanning the
crossover, and which in practice has meant that τ₂ is reported for a handful of well-studied
systems and assumed for the rest.

### 1.2 Why large amplitude, and why nobody has used it for this

The second harmonic of the surface-tension response is generated only when the interface is driven
far enough from equilibrium that the response is genuinely nonlinear. Because it is a small signal
sitting on a large one, it has usually been treated as a quality-of-measurement indicator — a
warning that the strain amplitude was too high for linear analysis — rather than as a source of
information.

That framing has recently begun to change. Czakaj et al. [CIT: Czakaj 2021] measured the second
harmonic ratio of the surfactant LAE across concentrations spanning the CMC and found a pronounced,
non-monotonic structure: the ratio rises sharply as the CMC is approached, peaks essentially *at*
the CMC, and declines above it. They interpret this as evidence of interfacial compositional
reorganisation rather than of micellar kinetics. Whatever the correct interpretation, the result
establishes the central point: in the large-amplitude regime the second harmonic carries
concentration-dependent information that the fundamental does not.

What has not been established is whether that information can be **inverted** — whether a
large-amplitude oscillatory dilatational (LAOD) measurement, of the kind a standard tensiometer can
actually perform in a reasonable time, contains enough information to recover τ₂ quantitatively,
and under what conditions. That is the question this paper answers.

### 1.3 The distinction this paper turns on: detectable is not invertible

The natural first question is whether the nonlinearity is *measurable*. It is easy to answer, and
answering it alone leads to a systematic error of reasoning that we found ourselves making and
which, we suspect, is common: a signal that clears the noise floor is not thereby a signal that
determines a parameter.

We therefore separate two criteria and insist on both:

- **Detectability.** The second-harmonic amplitude I₂ exceeds the instrumental noise floor. This is
  necessary but weak — it says the nonlinearity exists, not that it encodes τ₂.
- **Invertibility.** The observable *changes* with τ₂ by an amount large compared with its
  uncertainty. A detectable nonlinearity whose sensitivity to τ₂ falls below its own uncertainty is
  useless for inversion, and the harmonic ratio is very nearly that: across a doubling of τ₂ at
  fixed concentration and frequency, I₂/I₁ moves by 0.007–0.019 against a noise floor of 0.025–0.05
  on the ratio itself. The signal is real but not resolvable — a weaker and more defensible claim
  than calling it blind.

The observable that does carry τ₂ is not the amplitude ratio but the **phase lag**. This is
counter-intuitive — phase is usually regarded as the noisier measurement — and it is the pivot on
which the whole method rests.

### 1.4 Why an operational window exists at all

Above the CMC the micellar reservoir grows with total concentration, so one expects the reservoir
effect to become more pronounced as concentration is raised. The opposite happens. The reservoir
capacity C_m enters the problem twice: it strengthens the buffering of the interfacial depletion
layer, and it *smears* the phase response toward the diffusion-limited limit by holding the
subsurface concentration near its equilibrium value. Beyond a few times the CMC the second harmonic
falls below the noise floor and the phase becomes τ₂-insensitive again. The method therefore lives
in a bounded window above the CMC, and — importantly for anyone who tries to reproduce it — the
window has two edges, not one.

At exactly the CMC there is no micellar reservoir at all, and the phase lag becomes strictly
independent of τ₂. The informative region is not "near the CMC" in the loose sense used to motivate
the LAE experiments; it is a specific interval.

### 1.5 A hardware constraint that is easy to miss

Recovering a phase to a useful precision from a sampled waveform requires adequate sampling *within
each cycle*. This sounds trivial and is not: an undersampled demodulation does not merely produce a
noisy phase, it produces a **wrong** one, and the error does not announce itself as scatter around
the true value. In our own analysis an estimate taken at one frame per cycle returned values of
164°, 266°, 194° and 254° against a true 68°. Any precision figure computed from such data is
meaningless, and a design curve built from it will recommend operating points that cannot be
realised.

We therefore derive and validate an explicit sampling requirement (≥ 5 frames per cycle) and show
that the PAT-1M in oscillatory mode satisfies it with margin, while noting that the same
instrument's *static* mode does not — a distinction that has caused confusion in the literature,
because the two modes are reported in the same methods sections.

### 1.6 Contributions

1. A mass-conserving forward model of LAOD response that carries the micellar reservoir as an
   explicit immobile species with a linearised Aniansson–Wall exchange, so that total surfactant is
   conserved by construction rather than asserted. The model removes a structural defect of the
   two-compartment boundary-layer formulation in which the bulk reservoir never loses what the
   surface gains.
2. The **detectable-versus-invertible** distinction (§1.3), together with the demonstration that
   the harmonic ratio's sensitivity to τ₂ sits below its own noise while the phase lag carries the
   information, and that in the LAOD regime the phase is dominated by higher harmonics rather than
   by the fundamental.
3. A measured **concentration inversion window** for C₁₂E₈, 1.2–3.0 × CMC with an interior optimum
   at 1.5 × CMC, and a demonstration that the optimum is invariant under a fourfold change in the
   assumed limiting elasticity E₀ — i.e. that the window is a property of the transport problem and
   not of the isotherm parameterisation.
4. An explicit **sampling-density criterion** (≥ 5 frames per cycle) with a demonstration that
   violating it biases rather than scatters the phase estimate.
5. **Two inversion protocols with stated conditioning**: a one-parameter fit (Protocol A), in which
   capacity is fixed by the independent CMC/aggregation-number closure and only τ₂ is estimated,
   and a multi-concentration regression (Protocol B), in which τ₂ and a capacity scale are fitted
   jointly. Protocol A gives τ₂ to ≈ 0.2 % on synthetic data and is the headline result.
   **Protocol B is reported as a partial negative:** it recovers τ₂ to 0.5–1.6 % but the τ₂–capacity
   correlation remains between 0.1 and 0.6, does *not* fall monotonically with the number of
   concentrations, and in one of three isotherm anchors gets worse when a second concentration is
   added. We state this because a protocol whose conditioning is misreported is worse than no
   protocol.

### 1.7 Scope and exclusions

We treat a single nonionic surfactant with a Langmuir equation of state. We do not treat lateral
interactions, ionic double-layer effects, multicomponent mixtures, adsorption barriers, or
curvature corrections to the pendant-drop geometry. The ionic surfactants are used only to
demonstrate the predictive content of the timescale-to-bandwidth criterion of §6. We make no claim
about the interpretation of the LAE second-harmonic anomaly; on the contrary, we show that the
model cannot reproduce it, and we treat that failure as evidence about LAE's compositional
complexity rather than as a falsification of the micellar mechanism.

---

## §2. Coupled Reaction–Diffusion Model with Langmuir Interfacial Kinetics

### 2.1 Geometry, state variables, and what is conserved

We consider a planar interface at x = 0 bounding a semi-infinite surfactant solution. The bulk
carries monomer at concentration c(x,t) and micellised surfactant at concentration M(x,t), both in
mol m⁻³. The interface carries surface excess Γ(t) in mol m⁻². The interface is driven by a
prescribed area oscillation

    A(t)/A₀ = 1 + ε sin(ωt),                                    (1)

with strain amplitude ε and angular frequency ω. The model returns γ(t), from which observables are
extracted in §2.6.

The essential structural choice is that **micelles are an explicit species rather than an
absorbing boundary condition**. This is what makes the mass balance a property of the discretisation
rather than of a closure assumption.

### 2.2 Bulk transport and micellar exchange

    ∂c/∂t  = D ∂²c/∂x² − ∂M/∂t,                                 (2)
    ∂M/∂t  = −[ (M − C_m c) − (M_eq⁰ − C_m c₀) ] / τ₂.          (3)

Equation (2) states that monomer is transported by Fickian diffusion and exchanged with the
micellar pool. Equation (3) is the linearised Aniansson–Wall relaxation of the micellar population
towards its local equilibrium value M_eq(c) = C_m c, with relaxation time τ₂. Two points about the
form of (3) deserve emphasis.

First, the exchange term is written in the **conservative** form. Because ∂c/∂t contains −∂M/∂t
with the same expression, the sum c + M is unaffected by exchange: whatever the micelles release the
monomer gains, and vice versa. This is not cosmetic. The two-compartment formulation in common use
prescribes a boundary-layer thickness δ = √(πD/ω) and writes the subsurface flux as
(2D/δ)(c₀ − c_s), which imposes a flux across a fixed slab and lets the bulk reservoir supply
surfactant without ever being depleted. For a freshly created interface, or at the strain
amplitudes considered here, the resulting mass error is first-order, not a small correction.

Second, the exchange is linearised about equilibrium with the *offset* made explicit in (3). Without
the (M_eq⁰ − C_m c₀) term the equation relaxes towards equilibrium rather than about it, which
produces a spurious steady drift superimposed on the oscillation.

The reservoir is characterised by exactly two parameters, and both are required:

- **capacity** C_m = dM/dc (dimensionless) — *how much* surfactant the micelles can supply within a
  cycle. This is obtained from equilibrium data, not fitted:
  C_m = (c_tot − CMC)·n̄/CMC with a micellar-volume correction, where n̄ is the mean aggregation
  number. For C₁₂E₈, n̄ = 95.
- **τ₂** — *how fast* they supply it.

A model containing capacity without τ₂, or τ₂ without capacity, cannot separate "strong reservoir"
from "fast reservoir"; this is the central deficiency of the prototype formulation this work
replaces.

Boundary conditions are no-flux for micelles at the interface (micellar diffusion is neglected; the
micelles are much larger and slower than monomers and matter only at the highest frequencies) and
Dirichlet c = c₀ at the far boundary, which is placed at L = 8√(D/ω_min) so that the outer boundary
never sees the oscillation.

### 2.3 Surface balance

    dΓ/dt = D ∂c/∂x|₀ − Γ d(lnA)/dt.                            (4)

The first term is adsorption from the subsurface; the second is the dilution of surface excess by
area expansion. The interface is taken to be in local equilibrium with the subsurface
concentration c_s ≡ c(0⁺,t), so that

    Γ = Γ(c_s)                                                 (5)

holds instantaneously. Equation (5) is algebraic, and how it is used numerically is the single most
important implementation decision in the model.

**Writing c_s as the differential state — the "obvious" choice — is wrong.** If c_s is integrated
alongside the bulk, the interface satisfies Γ = Γ(c_s) at every instant by construction and is
therefore infinitely compressible: the model then reports an elasticity of order 1/E₀ and no
nonlinearity. The correct formulation carries Γ as the differential state and obtains c_s
*algebraically* by inverting (5). All results in this paper use the inverted form.

### 2.4 Equation of state

For a Langmuir layer, with Γ∞ the saturation excess and K the adsorption equilibrium constant,

    θ = Γ/Γ∞,      K c_s = θ/(1 − θ),                           (6)
    γ = γ₀ + nRTΓ∞ ln(1 − θ),                                   (7)

with n = 1 for a nonionic surfactant. Linearising (7) about the static point gives the limiting
elasticity

    E₀ = nRTΓ∞ · θ/(1 − θ).                                     (8)

Two features of (8) must be stated explicitly because they are frequent sources of error. E₀
**diverges** as θ → 1; the prefactor nRTΓ∞ is not a ceiling on the elasticity but a scale. And,
conversely, an assumed E₀ pins the operating coverage: a typical measured E₀ of 20–40 mN/m for a
single-chain surfactant corresponds to θ between roughly 0.6 and 0.85, not to the θ ≈ 0.95 that
large fitted adsorption constants tend to imply. Parameter sets whose K places the monolayer in the
saturated regime propagate an unphysical elasticity into every downstream number, and we recommend
(8) be used as a routine consistency check on any reported (K, Γ∞) pair.

Because this paper is restricted to C₁₂E₈, interactions are set to zero and (6)–(7) are used
throughout. The generalised (Frumkin) boundary used in Paper 2 replaces ln(1−θ) by
ln(1−θ) + aθ² in (7), which changes (6) to the implicit isotherm

    K c_s = [θ/(1 − θ)] · exp(−2aθ).                            (9)

We record (9) here only to fix the sign convention: the interaction term enters the isotherm with a
**negative** exponent and enters the equilibrium condition with the opposite sign, and a residual
with the sign reversed will converge to a spurious root. For a > 2 the isotherm is multivalued and
the branch choice becomes a physical assumption; this bounds the model and is deferred to Paper 2.

### 2.5 The natural dimensionless group

The problem has two timescales: the micellar relaxation τ₂ and the experimental period 2π/ω. Their
ratio is the controlling group,

    Ω ≡ ωτ₂,                                                    (10)

and the transport problem is additionally characterised by the diffusion time implicit in the
inverse-isotherm slope, ω_D = D(∂c_s/∂Γ)². For C₁₂E₈ at 1.5 × CMC, τ_D = 1/ω_D ≈ 6 ms against
τ₂ = 4 s: the diffusion time is shorter by a factor of 634 (about 2.8 decades), so on the timescale
of micellar exchange the interface is quasi-equilibrium with respect to diffusion, and τ₂ is the
slow process that sets the phase lag. (We note that for the fast ionic surfactants of §6 the ratio
grows to three to four decades.)

We note explicitly, because we tested and refuted it, that Ω alone does **not** collapse the
response across surfactants: at matched Ω the fundamental amplitude differs by a factor of 2.4
between C₁₂E₈ and CTAB. No single-group scaling law is claimed.

### 2.6 Observables

Fourier decomposition of the steady-state response over an integer number of cycles,

    γ(t) = γ̄ + Σ_k [ a_k cos(kωt) + b_k sin(kωt) ],              (11)

defines the harmonics I_k = √(a_k² + b_k²). We report:

- **I₁**, the fundamental amplitude (the linear response, ≈ ε|E|);
- **I₂**, the second harmonic, the measure of nonlinearity, and the detectability statistic;
- **the harmonic ratio** I₂/I₁, the quantity most often reported in the literature and, for the
  purposes of inversion, effectively uninformative (§1.3);
- **the phase lag φ**, defined as the lag of γ behind A:

      φ = phase(A) − phase(γ),   positive when γ peaks *after* A.   (12)

**Sign convention (resolved).** We define the **loss angle** δ as the lag of γ behind A, positive
when γ peaks after A. In the linear limit a diffusion-controlled soluble monolayer must give δ → 0
as ω → ∞ (insoluble limit: the surface excess cannot follow the area, so γ responds elastically and
in phase) and δ → 45° as ω → 0. Under the tabulated LVDT modulus
ε(ω) = ε₀(1+ξ+iξ)/(1+2ξ+2ξ²) with ξ = √(ν_D/2ν), the loss angle is δ = arctan[ξ/(1+ξ)], which is
bounded by 45° as ξ → ∞. The model's harmonic definition reproduces this expression identically.

We verified the correspondence rather than assuming it. With C_m = 0 and ν_D computed
self-consistently from the isotherm, ν_D = D(∂c_s/∂Γ)²/2π = 0.500 Hz for C₁₂E₈ at 0.5 × CMC, model
and LVDT agree to **0.48° at worst over 0.002–20 Hz** (e.g. 18.78° vs 18.43° at 1 Hz). An earlier
apparent 32°-vs-18.8° discrepancy was an artefact of comparing against a ν_D imported from an
unrelated system; it is not a convention conflict, and no offset correction is applied anywhere in
this paper.

**This 45° ceiling is the central diagnostic of the paper.** As §4 shows, a solution containing
micelles routinely produces δ **above** 45° — up to 80° for C₁₂E₈ — which is impossible for any
diffusion-only monolayer at any frequency. The excess above 45° is therefore a *categorical*
signature of a micellar reservoir, and its magnitude is what carries τ₂.

### 2.7 Numerical method and what has been verified

Equations (2)–(4) are integrated as a differential-algebraic system in the state vector
[Γ, c₁…c_{N−1}, M₁…M_{N−1}], with c_s recovered algebraically from (5). Space is discretised by a
finite-volume scheme on a geometrically stretched grid, N = 96 cells, with cell widths expanding as
(exp(σξ)−1)/(exp(σ)−1) for σ = 4, so that resolution is concentrated at the interface where the
gradient is steepest. Time integration uses an implicit BDF scheme at rtol = 10⁻⁷, atol = 10⁻¹¹
with a step ceiling of T/160. The last four of N ≥ 12 cycles are retained as steady state, and
harmonics are extracted by least squares over an integer number of cycles.

Three verifications are performed and reported in the SI:

1. **Mass conservation.** The relative drift of ∫c dV + ∫M dV + Γ over a run is computed, not
   assumed. It is 2.2 × 10⁻³ at the operating point; the residual is attributable to
   exchange across the semi-infinite Dirichlet boundary, which the domain can legitimately exchange
   with, and is monitored for that reason.
2. **Grid and cycle convergence** against N ∈ {64, 96, 144} and N_cycles ∈ {12, 18, 24}.
3. **The no-micelle limit.** Setting C_m = 0 must recover diffusion-limited behaviour; the phase lag
   must then be τ₂-independent, which provides a null test for the inversion machinery.

### 2.8 Deliberate omissions

The model neglects: micellar diffusion; any adsorption or desorption barrier (the interface is at
local equilibrium); lateral interactions (a = 0 throughout this paper); pendant-drop curvature and
capillary depletion; and any heterogeneity in the adsorbed layer. Each omission is either immaterial
in the regime studied or is bounded by an explicit test, and each is a candidate mechanism for the
discrepancy between model and the LAE data discussed in §5. The model is not offered as a general
description of surfactant dynamics; it is offered as the minimal structure that contains τ₂ and
capacity as separately identifiable quantities.

---

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
  fundamental amplitude by a factor of 2.4 (Section 6.5). Equation (13) is a design heuristic
  located by direct computation on three systems, not a collapse. Anyone applying it to a fourth
  surfactant should verify by running the forward model first.
- **The CTAB and SDS parameters are anchored, not measured.** E₀ = 30 and 22 mN/m were chosen from
  literature ranges, and the resulting coverages (θ = 0.72, 0.64) are lower than the saturated
  monolayers that measured surface excesses imply. We established analytically that the anchored
  isotherm is internally consistent, but it is not independently validated, and CTAB coverage is an
  open question requiring a measured Π–c isotherm.
- **SDS does not clear the noise floor anywhere in the scanned range**, reaching only 0.0048 mN/m at
  10 Hz. That is a band-edge result — our scan stops at 100 Hz, where Ω ≈ 24, not a demonstration
  that SDS fails.
- **The PAT-1M and CPT bandwidths are taken from the instrument specifications** supplied by our
  collaborators; we have not operated either instrument.

---

---

## §7. Inversion Protocols and Concluding Remarks

### 7.1 The structural reason a single measurement is degenerate

Before presenting the protocols it is worth stating exactly why two of them are needed. At the
1.5 × CMC operating point the loss angle responds to the two reservoir parameters with
almost identical lever arms and opposite signs:

    ∂δ/∂ln C_m  = +14.1° per e-fold                                       (14)
    ∂δ/∂ln τ₂   = −13.3° per e-fold                                       (15)

A larger reservoir *raises* δ; a faster reservoir *lowers* it. A single-frequency measurement
therefore sees only a weighted difference, and the correlation between the two fitted parameters is
0.997 — not as a numerical accident, but as a direct consequence of (14)–(15). This is why the
second-harmonic ratio is useless for the two-parameter problem regardless of how precisely it is
measured, and why the problem is genuinely under-determined at one concentration.

The two protocols correspond to the two honest ways of breaking that symmetry: **supply C_m
independently**, or **supply additional concentrations**.

### 7.2 Protocol A — one-parameter inversion (recommended)

Fix capacity from an independent equilibrium closure, C_m = (c_tot − CMC)·n̄/CMC, and estimate τ₂
alone from the frequency dependence of δ.

*Inputs:* a measured CMC and aggregation number n̄. That is all.
*Design:* c = 1.5 × CMC, ε = 0.20, 5–25 frames per cycle, ≥ 12 cycles, analyse the last four,
sweep the full instrument band.
*Estimator:* joint least squares over frequencies, or the Δ-method on the widest-lever-arm point.

*Performance at C₁₂E₈, 1.5 × CMC, σ_γ = 0.05 mN/m:* Monte-Carlo scatter of δ is 0.033° at
ω = 6.28 rad/s; the lever arm is 13.3° per e-fold; the resulting relative standard error on τ₂ is
**0.23%** by joint fit and 0.25% by the Δ-method. τ₂ = 4.000 s is recovered exactly on synthetic
data.

*Conditioning:* one parameter, well identified, no correlation to report. This is the protocol whose
precision may be quoted.

### 7.3 Protocol B — multi-concentration regression (exploratory; documented partial negative)

Where n̄ is unavailable, fit τ₂ and a capacity scale jointly across a concentration series in
1.2–2.5 × CMC.

*Performance:* τ₂ is recovered to 0.5–1.6% in every configuration we tested, across three isotherm
anchors and for CTAB as well as C₁₂E₈. **The correlation between τ₂ and C_m does not fall below
0.1, does not decrease monotonically with the number of concentrations, and in one of three anchors
gets worse when a second concentration is added** (0.60 → 0.40, with a sign change).

*Why we report this rather than presenting the protocol as a recipe:* an earlier version of our own
analysis found a clean monotone improvement from one to four concentrations at a single isotherm
anchor (−0.61 → −0.37 → −0.24) and we nearly published it as a property of the method. Re-running at
other anchors showed it was an artefact of that anchor. A protocol whose conditioning is
misreported is worse than no protocol, and the correct statement is that **multi-concentration data
partially breaks the degeneracy, by an amount that depends on where the concentrations sit relative
to the window edges and cannot be predicted from the number of points alone.**

*Practical guidance:* four to five points spanning 1.2–2.5 × CMC is a reasonable design, but a
two-point fit should be regarded as under-determined and should not be reported as a joint
determination of τ₂ and C_m.

### 7.4 Systematic error dominates the precision budget

We emphasise what the headline figure does *not* include. Because δ responds to C_m at
+14.1°/e-fold (14) and Protocol A takes C_m from n̄, any error in n̄ propagates directly:

| n̄ error | C_m error | induced δ bias | **τ₂ bias** |
|---|---|---|---|
| 5% | 5% | 0.69° | **5.2%** |
| 10% | 10% | 1.34° | **10.1%** |
| 20% | 20% | 2.56° | **19.3%** |

**A 10% error in n̄ produces a 10% bias in τ₂ — forty-four times the 0.23% statistical
uncertainty.** No amount of signal averaging reduces it. The precision of the method is therefore
bounded by the accuracy of an equilibrium measurement made by an entirely different technique, and
the appropriate headline for this work is not "τ₂ to 0.2%" but "**τ₂ to 0.2% statistical, 1:1
systematic in n̄**".

We regard this as the single most important caveat in the paper, and it is the reason Protocol A is
described as requiring a *measured* n̄ rather than a tabulated one. The same argument applies to the
CMC, which sets both the concentration axis and the capacity closure (§3.5).

### 7.5 What the method buys, and what it costs

**Buys:** a quantitative, model-calibrated route to τ₂ from an instrument already present in most
interfacial laboratories, without a frequency sweep beyond 1 Hz and without an independent
micellar-kinetics measurement. The 45° ceiling is a free, instrument-independent presence test
requiring no fitting at all.

**Costs:** a strain amplitude three times larger than the standard LAOD protocol (ε = 0.20 against
0.065), which requires bench validation of drop stability that we have not performed; a measured n̄
and CMC; adequate sampling, which must be verified rather than assumed; and operation at a
concentration chosen for invertibility rather than for signal strength.

### 7.6 Conclusions

1. **The loss angle, not the harmonic ratio, is the operative observable.** The ratio's sensitivity
   to τ₂ lies below its own measurement uncertainty; the loss angle's exceeds it by two orders of
   magnitude. Inversions based on harmonic ratios alone should be regarded as under-determined.
2. **A 45° ceiling makes the micellar reservoir categorically detectable.** Any δ > 45° proves
   kinetic relaxation is present, independent of D, E₀ or any fit. The reservoir *raises* δ rather
   than depressing it, acting as a lagging source term rather than an elastic element.
3. **Detectability and invertibility peak at different concentrations,** and the operating point
   must be chosen for the latter. The inversion window for C₁₂E₈ is 1.2–3.0 × CMC, optimum
   1.5 × CMC, and is invariant across a fourfold range in E₀. At exactly 1.0 × CMC, I₂ is maximal
   and the τ₂ sensitivity is zero.
4. **Sampling density is a bias constraint, not a noise constraint.** Below five frames per cycle
   the phase is wrong, not merely imprecise. The PAT-1M at 25–50 fps satisfies this with margin;
   its static mode, sometimes quoted at one profile per second, does not.
5. **The method converts into an instrument-selection criterion.** Requiring ωτ₂ ≳ 1 means the
   0.005–1 Hz band covers C₁₂E₈ but not CTAB or SDS, which need a CPT above ~16 Hz. This is offered
   as the explanation for the historical absence of micellar-exchange reports from low-frequency
   drop-shape work on fast ionics.
6. **Precision is statistical; accuracy is systematic.** τ₂ is recovered to 0.2–0.3% for a known n̄,
   but n̄ errors propagate one-for-one and dominate by a factor of 44.

### 7.7 What this work does not establish

- It does not validate the model against a measured C₁₂E₈ LAOD dataset. Every quantitative figure
  is a forward-model result with injected realistic noise; the τ₂ = 4 s benchmark is a literature
  value, not one we inverted from data. **The single most valuable next experiment is a measured
  C₁₂E₈ frequency sweep at 1.5 × CMC and ε = 0.20 to test the predicted δ(ω) curve directly.**
- It does not explain the LAE second-harmonic anomaly. Our model cannot reproduce it (§3.4),
  and we treat that as evidence about LAE's compositional complexity rather than as a falsification
  of the micellar mechanism. No claim is made about compositional reorganisation.
- It does not establish the CTAB or SDS parameters. Those coverages are anchored to literature E₀
  values, not measured, and the CTAB monolayer capacity is in tension with measured surface excesses
  (Section 6.5).
- It does not derive the ωτ₂ ≈ 3 sensitivity maximum. That is an observation on three systems, and
  no scaling collapse exists — at matched ωτ₂, CTAB and C₁₂E₈ differ in amplitude by 2.4×.

---

## References

Items 1–5 are as supplied by our collaborators. **Two carry corrections and one is unverified**;
these are flagged rather than silently corrected.

1. Czakaj, A.; Jarek, E.; Krzan, M.; Warszyński, P. Ethyl Lauroyl Arginate, an Inherently
   Multicomponent Surfactant System. *Molecules* **2021**, *26*, 5894. doi:10.3390/molecules26195894.
   — *Source of the LAE harmonic dataset used in §1; the vector-digitised Fig. 7 and Table 4 in our
   repository belong to this paper.*
2. Lucassen, J.; van den Tempel, M. Dynamic measurements of dilational properties of a
   liquid interface. *Chem. Eng. Sci.* **1972**, *27*, 1283–1291.
3. Aniansson, E. A. G.; Wall, S. N. Kinetics of step-wise micelle association. *J. Phys. Chem.*
   **1974**, *78*, 1024–1030.
4. Stadler, K. et al. High-frequency capillary pressure tensiometer with piezo-actuator operating up
   to 500 Hz. *Meas. Sci. Technol.* **2016**, *27*, 045301.
5. Risse, S. et al. *Adv. Colloid Interface Sci.* **2025**, *325*, 103050.
6. Loglio, G. et al. *Colloids Surf. A* **2001**, *186*, 63–71. — cited for large-amplitude
   oscillatory drop stability.
7. Kotsmar, C. et al. *Colloids Surf. A* **2012**, *413*, 20–25. — cited for profile-analysis
   tensiometry at ΔA/A₀ ≈ 20–30%.
8. Miller, R. et al. Transient and oscillatory surface dilatational rheology: principles and
   instrumentation. *Adv. Colloid Interface Sci.* **2004**, *111*, 17–34.

### Reference flags

- **Bond number (refs 6, 7 as applied in §5.7).** The supplied statement that Bo < 0.15 for
  1.0–1.5 mm outer-diameter tips is not true at the stated conditions. With
  Bo = ΔρgR²/γ and R taken as the tip *radius* (the conservative and correct choice),
  Bo = 0.184 at OD 1.5 mm and γ = 30 mN/m; Bo = 0.15 corresponds to OD ≈ 1.35 mm. The assertion
  holds for OD ≤ 1.2 mm at γ = 30 mN/m, and for all listed diameters only if γ ≳ 37 mN/m.
  More importantly, **the Bond number does not govern oscillatory detachability**, which is
  set by contact-line pinning and the dynamic contact angle. We therefore cite refs 6 and 7 for
  the *empirical* large-amplitude stability result and have removed the Bond-number argument from
  the text.
- **Instrument acquisition rate (ref 8).** The 25–50 fps figure is a manufacturer sensor
  specification (PAL 25 fps; USB 3.0 50–100 fps), not a measurement reported in Miller et al. 2004.
  In the manuscript this is attributed to the instrument specification, with ref 8 cited for the
  instrumentation principles only.
- **Refs 4, 5, 7** are cited from the supplied descriptions; we have not verified their contents
  directly and make no claim beyond the property each is cited for.

---


---

# Appendix (editorial — not part of the manuscript)

This appendix replaces the per-section "items needing your decision" blocks that accompanied the
drafts. Nothing here belongs in the typeset paper; it is the resolved record so that the typeset
pass can be checked against it.

## A.1 Number provenance

Every quantitative claim in §3–§7 traces to a committed artifact:

| claim | source |
|---|---|
| window 1.2–3.0× CMC, optimum 1.5× | `e0_sensitivity.json` |
| E₀ invariance over 25–103.5 mN/m | `e0_sensitivity.json` |
| δ vs frequency, Cₘ on/off (§4.2) | computed in FINDINGS §15.3 |
| 45° ceiling, 0.48° LVDT agreement | FINDINGS §15.1–15.2 |
| I₂/I₁ 0.007–0.019 per doubling | computed from K = 140 sweep |
| ≥ 5 frames/cycle | FINDINGS §11 |
| ε² law, 9.21× amplification | FINDINGS §16.1 |
| 0.23% relSE, MC scatter 0.033° | `recovery_uncertainty.json` |
| corr = 0.997, 2-parameter relSE 3.0% | `recovery_uncertainty.json` |
| Protocol B non-monotonicity | `protocol_b_anchored.json` |
| ∂δ/∂ln Cₘ and ∂δ/∂ln τ₂ | FINDINGS §17.1 |
| n̄ propagation table | FINDINGS §17.2 |
| CTAB/SDS bands, 16 Hz crossing | `ionic_generality.json` + manual extension |
| Bond number refutation | FINDINGS §17.3 |

## A.2 Corrections made to supplied material

1. **Loss-angle sign.** The supplied narrative had the reservoir *depressing* δ (32.4° → 18.8°) and
   called that depression the signature. Computation gives the opposite: the reservoir *raises* δ to
   68–80°, above the 45° diffusion ceiling. The structure was kept; the sign was reversed. Endorsed.
2. **LVDT "convention conflict."** Not a conflict. With ν_D computed self-consistently the agreement
   is 0.48°. The earlier 18.8°-vs-32° gap was our own parameter borrowing.
3. **45° vs 35°.** The ceiling is 45°, not arctan(1/√2) ≈ 35.3°.
4. **"I₂/I₁ is blind"** → "unresolvable": it moves 0.007–0.019 against a 0.025–0.05 noise floor.
5. **"Two concentrations suffice"** — retracted; the correlation is not monotone in the number of
   concentrations.
6. **"Large amplitude trades linearity against τ₂ sensitivity"** — false; δ is constant to 0.17°
   over ε = 0.03–0.30.
7. **Bond number** — removed; Bo = 0.184 at the stated worst corner (OD 1.5 mm, γ = 30 mN/m), and
   Bo does not govern oscillatory detachment.
8. **Frumkin "numerical stiffness at a ≥ 4"** — false; it was a wrong-sign residual in the
   isotherm solver. Fixed.

## A.3 Multiplicity note (statistical)

The manuscript reports several derived figures from a family of model sweeps. Two deserve explicit
labelling because they were selected *after* inspection rather than pre-registered:

- **ωτ₂ ≈ 3 sensitivity maximum (§4.2, §6.1).** Located by scanning, not derived. No scaling
  collapse exists; at matched ωτ₂, CTAB and C₁₂E₈ differ in amplitude by 2.4×. Presented as an
  observation.
- **The 13.3° per e-fold lever arm at ω = 6.28 (§5.6, §7.2).** This is the largest of four
  frequencies scanned (5.99, 8.89, 11.16, 13.26 at ω = 3/4/5/6.28). Operating at the top of the
  instrument band because the lever arm is largest there is legitimate, but the reported precision
  is the best of the four, not a typical value. Report as "0.23% at ω = 6.28 rad/s" rather than
  "0.23%" without the qualifier.

## A.4 Provenance audit (final)

Verified via the GitHub API on 2025-11-03 against the user's stated commits:

| commit | repo | status |
|---|---|---|
| `f73fa80` | surfactant-dynamics-sim | **resolves** |
| `e9da2d0` | surfactant-dynamics-sim | **resolves** |
| `caa1f57` | Surfactant-Toolkit (public) | **does not exist** |
| `85e1218` | Surfactant-Toolkit (public) | **does not exist** |
| `27691b1` | Surfactant-Toolkit `main` | **is** the frozen public tip, as stated |

Both `85e1218` and `caa1f57` were described as internal-mirror commits. Neither is a public artifact,
so neither should be cited in the manuscript or in any preprint announcement. The manuscript text
lives on `arena/fa2f9ecf-surfactant-toolkit` at `a92c123`, `a88bd29`, `34aedd3`, `dd48f4c`.

**I cannot verify the current `.docx`/`.pdf`.** I have no access to the internal SurfactantKit
mirror. The synthesis documents I *can* see (`f73fa80`) were consistent with the corrected §4/§7 and
carried the two errors since fixed (Bond number; the "2.4 to 73 Hz" head frequency). If `e9da2d0`
applied those fixes, good; if `caa1f57` is where the fixes actually live, then the public-facing
record still needs them.

## A.5 What is not verified

- **Refs 4, 5, 7** (Stadler 2016, Risse 2025, Kotsmar 2012) are cited from supplied descriptions.
  We have not read them. They are cited only for the specific property named in the reference list.
- **The 25–50 fps acquisition rate** is a manufacturer sensor specification, not a figure from
  Miller et al. 2004. The text attributes it correctly; the reference is cited for instrumentation
  principles only.
- **No formatting spec exists in this repository.** `MANUSCRIPT_FORMATTING_SPEC.md` is not present
  on any branch of `Surfactant-Toolkit`, so this draft is plain Markdown with no journal-specific
  typesetting applied. The typesetting pass must supply the spec.
