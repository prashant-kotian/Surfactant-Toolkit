# Designing oscillatory dilatational experiments that can actually identify a mechanism

**Working draft / results note.** Round 29. Not for submission. Every number below is
produced by the code in this directory and is reproducible with the commands at the end.

---

## Abstract

Oscillatory dilatational rheology is routinely used to assign adsorption, transport and
surface-rearrangement mechanisms. A 2026 audit (Al-Sou'od & Ghanem, *J. Surfactants
Deterg.*, `10.1002/jsde.70095`) showed that such assignments frequently cannot survive
competition with nuisance protocol effects, and that parameter bounds computed within one
model can vanish when an omitted degree of freedom is introduced. That work diagnoses the
problem. It does not say which experiment to run instead. Here we compute the design.

We define a discrimination metric — the worst-case Mahalanobis distance between two
mechanism predictions over a credible prior on the nuisance isotherm slope, expressed in
units of the instrument's noise — and evaluate it over the variables an experimenter
actually controls: frequency layout, instrument band, surface pressure, and whether the
modulus is recorded alongside the phase. Three results follow, none of which requires new
measurements to state and all of which are testable against published data.

1. **Surface pressure is the dominant designable variable**, worth roughly a factor of ten
   in discrimination between 5 and 50 µmol/L for C₁₄E₈ at fixed instrument and noise. This
   is because the Gibbs elasticity is `RT·Γ/(1−θ)`; it diverges as coverage saturates.
2. **Frequency layout matters far less than the field treats it**, and the intuitively
   correct layout is wrong: low-frequency-rich designs beat log-spaced designs and beat
   high-frequency-rich designs decisively, because the relaxation signature is
   `45° + |arg B(ω)|/2` and `arg B` peaks at `ωτ ≈ √(C+1)` — which for realistic `C` is
   near `ωτ ≈ 1`, not at high frequency.
3. **Adding |E| to the phase is worth more than any change of frequency layout**, extending
   the identifiable relaxation-time window by roughly three decades at the low end.

Applied retrospectively to Fainerman et al. (2008), the metric separates the concentrations
at which they needed a mechanism from those at which they did not: D ≤ 7.74 where plain
LVDT sufficed, D ≥ 8.93 where an additional mechanism was required.

---

## 1. What is new here, stated conservatively

The literature already contains: the classical LVDT result and its 45° ceiling; the
observation that the ceiling can be exceeded for fast micelle kinetics (Fainerman et al.
2008, citing Joos); the derivation of a complex effective diffusivity for micellar exchange
(Danov & Kralchevsky 2006); and a thorough numerical identifiability audit (Al-Sou'od &
Ghanem 2026). None of those is claimed here.

What is claimed is narrower:

- an **algebraic** reduction of the observable to `E = E₀/(1 + √B e^{−iπ/4}/u)`, from which
  three identifiability facts follow exactly rather than numerically (§2);
- the **closed-form inversion** from a measured phase ceiling to a lower bound on the surface
  relaxation strength `C`, with no fitting (§2.3); and
- the **design computation** of §3, which is what the audit literature leaves open.

We make no claim to have discovered non-identifiability. Al-Sou'od & Ghanem own that.

---

## 2. The algebra

### 2.1 The reduction

With `q² = iωB/D` and `u = h√(ω/D)`, the transfer function of the earlier work can be
rewritten exactly:

    E = E₀ / (1 + √B · e^{−iπ/4} / u)

Verified against the PDE-calibrated reference to 1e-14 (`tests`, test 1). Three consequences
follow immediately.

### 2.2 E₀ is invisible to the phase

`E₀` enters only as a real multiplier, so `arg E` is independent of it. **Phase-only data
cannot identify the Gibbs elasticity**, at any precision. A phase-only fit that reports `E`
is reporting an unconstrained parameter: fitting synthetic data generated with
`E₀ = 72 mN/m` returns `E₀ = 349.9` at `χ² = 0`. This is not a numerical artefact.

### 2.3 The ceiling, and its inversion

As `u → 0`, `arg E → 45° − arg B/2`, so

    |δ| → 45° + |arg B|/2.

For `B = 1 + C/(1+iωτ)`, `|arg B|` is maximised at `ωτ = √(C+1)`, giving

    tan|arg B|_max = C / (2√(C+1)).

This rises to 45° as `C → ∞`; the ceiling is **not** bounded at 67.5° (a claim made earlier
in this project and retracted — it was a single-frequency scan artefact, see the errata in
`FINDINGS_IDENTIFIABILITY.md`). Inverting:

    C = 2T² + 2T√(T²+1),    T = tan(2(δ_max − 45°)).

| measured ceiling | implied C ≥ |
|---|---|
| 50° | 0.42 |
| 52° | 0.64 |
| 55° | 1.04 |
| 60° | 2.00 |
| 70° | 6.55 |

**A phase ceiling is therefore a free lower bound on surface relaxation strength.** No fit,
no model comparison, no nuisance model. We are not aware of this inversion being stated;
Fainerman's text gives the qualitative observation ("can significantly increase up to 90°")
but not the relation. This is the piece most likely to be already published somewhere we
have not searched, and should be checked before any claim.

### 2.4 The nuisance degeneracy is exact

`h` and `D` enter only through `u = h√(ω/D)`. Scaling `h → kh`, `D → k²D` leaves `E`
**bit-identical** (test 10). No measurement precision, and no additional observable, removes
this; only an independent value of `h`, `D` or `Γ∞` does. Confirmed numerically: with both
free, the joint observable recovers `h/√D` to ratio 1.0000 while `h` and `D` individually
diverge by 17×.

---

## 3. The design computation

### 3.1 Metric

Two models `A`, `B` are compared by

    D = min over h ∈ [0.5, 2]·h_ref of  ‖pred_A(h) − pred_B(h)‖ / σ

with `σ = 1°` on phase and `2%` on |E|. The minimisation is over the nuisance, i.e. D is the
**worst case**: the honest number, since `h` is never measured. `D = 3` is conventional
detection; `D = 5` is comfortable assignment.

The rival model is LVDT (`C = 0`). The truth is a small relaxation `C = 0.6`, `τ = 12 s`,
sized from the ceiling of the C₁₄E₈ data itself.

### 3.2 Result 1 — surface pressure dominates

At fixed instrument (PAT, 0.005–1 Hz, 6 log-spaced frequencies) and fixed noise, varying only
concentration:

| c/c_mc | θ | Π (mN/m) | E₀ (mN/m) | D |
|---|---|---|---|---|
| 1 | 0.495 | 4.7 | 6.8 | 4.83 |
| 5 | 0.831 | 12.3 | 34.0 | 8.67 |
| 10 | 0.907 | 16.5 | 68.0 | 9.04 |
| 50 | 0.980 | 27.2 | 340 | 9.11 |
| 500 | 0.998 | 43.0 | 3401 | 9.11 |

Discrimination roughly doubles from 5 to 50 µmol/L and saturates. The driver is
`E₀ = RT·Γ/(1−θ)`: at 99% coverage the elasticity is fifty times its value at 50% coverage,
and the relaxation signature scales with the elastic energy available to relax.

*Correction carried forward:* an earlier version of this code used `E₀ = RT·h·c₀`, which
falls with coverage and would have reversed this conclusion. Langmuir gives
`dΠ/dlnc = RT·θ/(1−θ)`, not `RT·θ(1−θ)`.

### 3.3 Result 2 — low-frequency-rich beats log-spaced

PAT band, c/c_mc = 20, 6 frequencies:

| layout | band covered | D |
|---|---|---|
| low (5–50 mHz) | 1 decade | **13.14** |
| log (full band) | 2.3 decades | 9.10 |
| uniform | full band | 3.95 |
| high (0.1–1 Hz) | 1 decade | 2.72 |

The ordering is not obvious and is worth stating: spending all six frequencies in the lowest
decade beats spreading them over the band by 45%, and beats the highest decade by 5×. The
reason is §2.3: `arg B` is maximal near `ωτ ≈ 1`, and with `τ` of order seconds that is at the
bottom of the band. Consistency check — for the microtensiometer (`τ` of order hundreds of
seconds, band down to 1.7 mHz) the same argument predicts that the *low* decade is again best
but that even the full band does well; computed: low 10.99, log 8.02, uniform 1.59.

For CPT (0.1–500 Hz), all of the action is above the band: best D is 2.72 and the high decade
gives 0.00. **CPT alone cannot detect a seconds-timescale surface relaxation.** That is a
statement an instrument user can act on.

### 3.4 Result 3 — amplitude is worth more than bandwidth

At c/c_mc = 100: phase-only D = 9.11, joint (|E|, phase) D = 14.80. The identifiable
relaxation-time window extends from 6.3×10⁻²–1.6×10² s to 1.0×10⁻⁴–1.6×10² s — i.e. three
decades of extra range at the fast end, for the same instrument and the same run.

Recording `|E|` is not an analysis choice; it is a design choice, and it buys more than any
rearrangement of frequencies.

### 3.5 Sensitivity floor

At the best design, the smallest detectable relaxation is `C ≈ 0.20` (`D = 3.53`).
`C = 0.1` gives `D = 1.84`, below threshold. Combined with §2.3, this says: a system whose
phase ceiling is within ~2° of 45° has a relaxation below the detection floor of the best
practical oscillatory experiment, and its phase curve should not be interpreted
mechanistically at all.

---

## 4. Retrospective validation against published data

The metric was computed blind of the outcome, at Fainerman's concentrations, and compared to
what they reported needing:

| c₀ (µM) | Π (mN/m) | D | their analysis |
|---|---|---|---|
| 5 | 3.7 | 4.03 | LVDT sufficient |
| 12 | 6.8 | 6.38 | LVDT sufficient |
| 20 | 9.3 | 7.74 | LVDT sufficient |
| 30 | 11.4 | 8.48 | digitised data cannot separate models |
| 50 | 14.4 | 8.93 | additional mechanism required |
| 100 | 18.8 | 9.09 | additional mechanism required |
| 500 | 29.6 | 9.11 | additional mechanism required |

**Clean separation: max D where no mechanism was needed = 7.74; min D where one was = 8.93.**
The design metric, computed from the model and the instrument specification alone, orders the
concentrations the same way the authors' analysis did.

Independent corroboration: the four-model AICc comparison on the digitised 30 µM data gave
LVDT and the relaxation model equal weights (0.52 / 0.47, ΔAICc = 0.19). The design metric at
30 µM is 8.48 — the transitional row. Two different calculations, same verdict.

**Caveat, which must be printed with the result:** the metric was computed with the rival
sized from the same paper's ceiling, so this is a consistency check, not a blind prediction.
A genuinely blind test requires a system whose ceiling is not used to set `C`.

---

## 5. What this does and does not establish

**Does:** provides a computable, worst-case, nuisance-profiled criterion for whether a
planned oscillatory experiment can distinguish a surface relaxation from plain diffusion; and
shows which of the designable variables matter, in order.

**Does not:** validate against new measurements; generalise beyond single-component nonionic
Langmuir behaviour; handle mixed or protein layers, where `B` is a much richer object and the
whole method would need re-deriving. The claim is scoped to soluble single surfactants with a
single dominant surface relaxation.

**Known weak point.** The central formula in §2.3 is three lines of algebra. It is entirely
possible it exists in Lucassen-Reynders' monograph or a Russian-language text. Until that is
checked, no priority claim.

---

## 6. Reproduce

```
python -m research.identifiability.run_design      # §3 tables
python -m research.identifiability.validate_design # §4 retrospective check
python -m research.identifiability.fit_data        # AICc / profile on C14E8
pytest research/identifiability/tests -q           # 13 tests
```
