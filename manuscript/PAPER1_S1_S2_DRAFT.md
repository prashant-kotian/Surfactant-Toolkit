# Direct Inversion of Micellar Relaxation Kinetics from Large-Amplitude Oscillatory Dilatational Surface Rheology: Theory, Bandwidth Limits, and Robust Parameter Identification for Model Nonionic Surfactants

**Paper 1 of 2 — the C₁₂E₈ foundation paper**

Draft of §1 and §2. Working document; not for circulation.

**Scope note for co-authors.** This paper is restricted to the nonionic surfactant C₁₂E₈ on a
profile-analysis tensiometer (PAT-1M, 0.005–1.0 Hz) with a Langmuir surface equation of state.
Every quantitative claim in §3–§7 is computed for that system. The ionic surfactants (CTAB, SDS)
appear only in §6, as a *predictive* instrument-selection criterion derived from their relaxation
timescales — not as experimental validations. Multicomponent and Frumkin-interaction effects are
deferred to Paper 2. See "Open items" at the end of this draft.

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
  uncertainty. A detectable nonlinearity that is insensitive to τ₂ is useless for inversion, and
  the harmonic-ratio observable turns out to be almost exactly that: across a factor of two in τ₂
  at fixed concentration and frequency, I₂/I₁ moves by roughly 0.004, which is an order of
  magnitude below the noise.

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
   the harmonic ratio is effectively τ₂-blind and the phase lag carries the information.
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

## §2. Model Formulation

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

**Sign convention (open item).** In the linear limit, a diffusion-controlled soluble monolayer
must give φ → 0 as ω → ∞ (insoluble limit: the surface excess cannot follow the area, so γ responds
elastically and in phase) and φ → 90° as ω → 0 (fully relaxed, dissipative response). The model
reproduces both trends, returning 8.4° at ω = 60 rad/s and 41° at ω = 0.031 rad/s in the
no-micelle limit. Note that this convention places both limits in (0°, 90°) and never permits a lag
above 90°.

The LVDT complex modulus as usually tabulated, ε(ω) = ε₀(1+ξ+iξ)/(1+2ξ+2ξ²) with ξ = √(ν_D/2ν),
is expressed as a complex amplitude whose argument is not the lag under our convention: we obtain
−arg ε ≈ +32° where the model gives +18.8° for nominally the same physical situation, and the sign
of the discrepancy is convention, not physics. Sign-level agreement must therefore be pinned to the
instrument's own reference channel before any overlay of model and literature is attempted. We flag
this as an open item rather than resolving it by assertion.

A useful diagnostic follows from the LVDT form itself: |arg ε| is **bounded** by arctan(1/√2) ≈
35.3° for every (E₀, ν_D, ν) — we verified a maximum of 32.4° over a wide parameter grid. A measured
phase lag far outside that range is therefore a categorical signature that the response is not
diffusion-controlled, which is the role the lag plays in §4.

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

## Open items for co-authors before §3 is drafted

1. **LVDT sign convention (§2.6).** Must be pinned to a reference. Until it is, comparisons to
   published LVDT moduli are magnitude-only. This is the highest-priority item: it affects §4 and
   any figure that overlays model and literature.
2. **Capacity scaling in the derivative.** The closure C_m = (c_tot − CMC)n̄/CMC is linear in
   concentration for fixed CMC. If the reservoir contribution scales this way, the inversion
   sensitivity to n̄ is strong enough that n̄ must be sourced rather than assumed. Confirm whether
   the derivative should include the term from 1 + n·c_mic/c₁ (see `systems.py` note 1).
3. **References.** All [CIT] markers need completion. I have not invented any citations; where a
   claim rests on a specific paper I have named the author, and where it rests on this
   repository's own computation I have said so.
4. **§2's title in the agreed outline** is "Coupled Nonlinear Reaction–Diffusion Dynamic Model
   (mass-conserving Aniansson–Wall + Frumkin boundary)". For Paper 1 the boundary is **Langmuir**
   (a = 0); Frumkin is Paper 2. Suggest retitling to "(mass-conserving Aniansson–Wall + Langmuir
   boundary)" so the title does not overclaim the scope.
5. **E₀ anchor statement.** Given §2.4, §3 must state the C₁₂E₈ anchor explicitly (K = 140 m³/mol
   ⇒ θ = 0.937, E₀ = 103.5 mN/m at 1.5 × CMC) and report the invitation-window result alongside its
   re-run at E₀ = 25 mN/m. The conclusion is invariant; the numbers are not identical.
