# Paper 1 — Abstract, §7, Conclusions, and References

Companion to the §1–§2, §3–§4 and §5–§6 drafts. This completes the manuscript body.

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
- It does not explain the LAE second-harmonic anomaly. Our model cannot reproduce it (FINDINGS §8),
  and we treat that as evidence about LAE's compositional complexity rather than as a falsification
  of the micellar mechanism. No claim is made about compositional reorganisation.
- It does not establish the CTAB or SDS parameters. Those coverages are anchored to literature E₀
  values, not measured, and the CTAB monolayer capacity is in tension with measured surface excesses
  (FINDINGS §13.9, §15).
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

## Note on the external synthesis document

The compiled deliverable in the `surfactant-dynamics-sim` repository (commit `f73fa80`) is
consistent with the corrected §4 and §7 positions: it carries the ε² verification (slopes
1.84–2.00, factor 9.21 vs 9.5 predicted) and the honest Protocol B statement (correlation 0.1–0.6,
non-monotone). Two items in it need correcting before it is treated as the record:

1. **The Bond-number sentence** carries the Bo < 0.15 error above and should be replaced with the
   empirical stability citation alone.
2. **The head-frequency claim for CTAB/SDS** ("2.4 to 73 Hz") is the design-rule estimate; the
   measured threshold at which I₂ clears the noise floor is ≈ 16 Hz. Both figures should appear, or
   the table should say "design rule" and quote Ω rather than Hz.

Two audit items, stated plainly because they affect provenance rather than science:

- **Commit `85e1218` in `Surfactant-Toolkit` does not exist.** The GitHub API returns
  "No commit found for SHA: 85e1218", and no such object exists on any local ref. The manuscript
  sections in this branch are at commits `a92c123`, `a88bd29` and `34aedd3`. Commit `f73fa80` in
  `surfactant-dynamics-sim` does resolve.
- **Neither `.docx` nor `.pdf` is present in this repository**, and the manuscript has not been
  compiled to Word here. The deliverables described above are in the other repository.
