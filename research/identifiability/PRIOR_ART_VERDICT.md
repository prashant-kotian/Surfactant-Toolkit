# Prior-art verdict on the identifiability direction

**Date:** 5 November 2026 · **Round 28** · Based on full texts now in hand, not
search snippets.

---

## 1. The paper that was sitting on the target

**Al-Sou'od & Ghanem, *J. Surfactants Deterg.* 2026, 1–16**, `10.1002/jsde.70095`,
"Protocol-Induced False Identifiability in Dynamic Surface Tension and Dilational
Measurements of Surfactant Interfaces."

This was not in my survey. It occupies most of what I claimed:

| I claimed | their paper |
|---|---|
| identifiability is rarely asked in this field | the entire subject |
| model comparison is the field's weak point | AICc, leave-one-segment/frequency-out, parametric bootstrap, Monte Carlo digitisation propagation |
| nuisance parameters drive conclusions | "protocol-induced" false identifiability; drift and segment-amplitude models beat the mechanism model by 5.68 and 16.84 AICc |
| profile likelihood is the right tool | used explicitly, including a boundary-case bootstrap |
| a reanalysis of published data would be valuable | two independent digitizations of Hansen 2008, 148 consensus points, plus a Poloxamer frequency-domain screen |

They also state the caveat I intended to make: "a parameter may be bounded within one
model while disappearing when an omitted nuisance degree of freedom is introduced."
And they name my second candidate's framing as prior art: Rosales-Anzola & Urbina-Villalba
(2025) on equation-of-state selection, Rosales-Anzola (2025) on thermodynamic-model
dependence.

**So "identifiability of interfacial relaxation mechanisms" as a topic is taken.** Proposing
it would repeat the failure mode that killed the last paper, one level up. I am not doing it.

---

## 2. What they did *not* do, and what that leaves

Three things, none of which they attempt:

### (a) Analytic structure

Their approach is narrative and numerical: a nested hierarchy M0/M1/M2, fitted, compared,
bootstrapped. There is no algebra. Concretely:

- **P1** — The 45° ceiling is never generalised. Risse et al.'s review states only the
  classical result ("in the absence of in-plane interactions, surfactants follow the
  Lucassen–Van den Tempel model"). Danov & Kralchevsky (*Adv. Colloid Interface Sci.* 2006,
  119, 17 — read in full) derive exactly the physics I use: under quasi-equilibrium
  micellisation the surface exchange acts as a complex `D_eff(ω)`, and they compute
  adsorption *kinetics*. They never form the complex dilational modulus, never mention the
  phase, and never state the ceiling. Zero occurrences of "dilational modulus", "complex
  modulus", "phase angle".
- **P2** — Danov & Kralchevsky handle polydispersity (15 uses of the word in part 1) as a
  *correction to the relaxation rate*, not as a continuum of relaxation times impinging on
  an observable. The distribution-of-times treatment with its 45–67.5° bound is absent.
- **P3** — The analytic identifiability statements are absent. Al-Sou'od use *profile
  likelihood numerically*; nobody states that `E₀` cancels from `arg E` identically, or that
  `h` and `D` enter only through `h/√D` so the degeneracy is exact and unremovable at any
  precision. Those are model-independent facts about the observable, not statements about
  one dataset — which is why they survive where a six-point fit does not.

### (b) The design question

Grep counts across Al-Sou'od: `D-optimal` 0, `optimal design` 0, `experimental design` 0,
`design of experiment` 0, `information matrix` 0. They *diagnose* non-identifiability. They
never ask which experiment would fix it. Their §5.6 closes with a list of *falsification*
experiments (better samples, synchronized A(t) and γ(t), varied step amplitudes) but no
computed design.

MBDoE is mature in the neighbouring literatures (`10.1007/s10450-025-00653-0` for isotherms,
Hunter–Reiner + G-criterion) and absent here. That gap is now *demonstrated* rather than
asserted, because the paper that did the diagnosis is in hand and its omissions are
checkable.

### (c) The ceiling as a screen

Consequence of P1. If every passive relaxation is confined to 45–67.5°, then a measured
|δ| above 67.5° falsifies the whole mechanism class in one number, with no fit, no model
comparison, no nuisance model. That is a screening instrument Al-Sou'od's workflow does not
contain — and it costs nothing to apply to historical data.

---

## 3. Honest assessment of what remains

**Strength of the remaining position:** the analytic results are exact, verified, and
checkable by a referee in an afternoon. They are not "our fit says" — they are "this
expression has this form". Al-Sou'od's contribution is valuable and largely empirical; mine
would be the theory that explains *why* their empirical finding occurs and bounds how far it
can go.

**Weakness, stated plainly:**

1. **Priority on P1 is plausible but unproven.** The derivation is three lines from
   Lucassen–van den Tempel. It is entirely possible it appears in Lucassen-Reynders'
   monograph, in Kretzschmar/Fainerman/Miller, or in a Russian-language text I cannot search.
   If it does, Results P1 and the ceiling become a remark, and only P3 and the design study
   carry the paper.
2. **A theory paper without data is a hard sell.** Al-Sou'od have 148 digitised points;
   I have 6 from one figure. The design study (b) is the answer — it is computation, not
   data — but it must be *useful*, not merely correct.
3. **Diminishing distinctiveness.** A year ago this direction would have been wide open. It
   is now one 2026 paper old, and the field is clearly moving toward identifiability
   auditing. The window is narrow and closing.

---

## 4. Recommendation

**Narrow to (b), with (a) as the supporting theory, and drop any framing of "identifiability
of mechanisms" as the contribution.**

Concretely the paper becomes: *given that oscillatory dilatational data frequently cannot
distinguish mechanisms (Al-Sou'od 2026, and independently demonstrated here on C14E8), what
minimal experiment can?* — answered by computing the design, with the analytic results
supplying the reason the naive designs fail.

That is defensible, constructive, does not require new bench data, does not compete with the
2026 paper, and uses the software. It is also the only version where being second is an
advantage rather than a problem, because it builds on their finding rather than racing it.

**Kill criterion:** if P1 turns out to be published, the paper rests on P3 + design alone.
That is thinner but still viable. If P3 is also published, stop.

---

## 5. What would settle it

- Lucassen-Reynders, *Colloidal Surfactants* (monograph) or her 2001 review on dilational
  rheology.
- Kretzschmar, Miller, Fainerman review series.
- Any Russian-language surface-rheology text covering the complex modulus of layers with a
  relaxation spectrum.
- Rusanov's surface thermodynamics.

If none of these contain `45° + |arg B|/2`, P1 stands.
