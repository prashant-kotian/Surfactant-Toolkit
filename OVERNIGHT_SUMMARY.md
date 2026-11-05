# Overnight summary — 5→6 November 2026

Read in this order. Everything is in `research/identifiability/`.

---

## The one thing that changed

**I retracted my own headline result** (`bc2c2be`). The claim that a passive surface
relaxation cannot produce a phase above 67.5° was wrong — it came from scanning the
relaxation time at a single frequency, so the optimum was never reached. The correct
ceiling rises towards 90°, and Fainerman 2008 says so in words.

It was caught only because you gave me the Fainerman PDF and I read his discussion of
phase angles above 45°. My own test "confirmed" the wrong bound because the test scanned
one dimension of a two-dimensional optimum. The test now checks the analytic optimum
against a 2D numerical search, which is what would have caught it the first time.

Left visible in the record rather than quietly deleted.

---

## What was built

| file | what it is |
|---|---|
| `PAPER_DESIGN_NOTE.md` | **start here** — the drafted results note, ordered by how strong each result is |
| `PRIOR_ART_VERDICT.md` | why the original framing was dropped (Al-Sou'od 2026 owns it) |
| `FINDINGS_IDENTIFIABILITY.md` | round-27 results, with errata marked in place |
| `design_engine.py` | the metric: worst-case model distance, nuisance-profiled |
| `run_design.py` | the three design studies |
| `validate_design.py` | retrospective check against Fainerman's reported outcomes |
| `from_structure.py` | design driven by SurfactantKit's Nagarajan CMC |
| `tests/` | 19 tests; the full research suite is 83 passed, 1 skipped |

---

## The four results, ranked by how much I trust them

**1. Surface pressure dominates everything else.** D roughly doubles from 5 to 50 µmol/L
and saturates. Driver: E₀ = RT·Γ/(1−θ), which diverges at saturation. *(High confidence —
but note I had the elasticity formula wrong first, and it inverted this conclusion until
caught.)*

**2. Low-frequency-rich beats log-spaced frequency layouts**, 13.14 vs 9.10, and beats
high-frequency-rich 2.72. Because `arg B` peaks near ωτ ≈ 1, which for second-scale
relaxations is the bottom of the band. CPT alone (0.1–500 Hz) scores 0.00 — it cannot see
a seconds-scale relaxation at all. *(High confidence; checkable by anyone.)*

**3. Recording |E| alongside phase is worth more than any change of frequency layout.**
D 9.11 → 14.80, three decades more identifiable relaxation-time range. *(High confidence.)*

**4. The retrospective check separates cleanly.** Max D where Fainerman needed no
mechanism: 7.74. Min D where he did: 8.93. *(Medium — see the caveat below.)*

**Caveat that must travel with result 4:** the rival model was sized from the same paper's
phase ceiling, so this is a consistency check, not a blind prediction. A real test needs a
system whose ceiling isn't used to set C.

---

## The structure-driven result, which surprised me

Using SurfactantKit's Nagarajan CMC, every surfactant's optimum sits at the concentration
ceiling, and `c/c_mc` turns out to measure the wrong thing. A given surface pressure is
reached at a *lower* multiple of the CMC for a longer chain, so the practical limit binds
earlier in those terms. The design collapses to one instruction:

> highest surface pressure the phase diagram permits, lowest frequency decade the
> instrument offers, record the modulus.

Weaker-sounding than a table of optima, more useful. And unreachable without the structure
step — with a fixed assumed CMC, all four surfactants look identical.

**Limitation, stated plainly:** Γ∞ and K are still calibrated, not molecular. Nagarajan
predicts the CMC, not the adsorption isotherm. So the absolute D values are an ordering,
not predictions. An isotherm module is the obvious next extension.

---

## Three unit errors, all mine

Caught and fixed overnight: elasticity scale using `RT·h·c₀` instead of `RT·Γ/(1−θ)` (this
one inverted a conclusion); CMC in mol/L vs mol/m³; and a stray factor of 1000 in a print
statement. Each now has a test or an assert. I'm flagging the count rather than hiding it —
it suggests the numerical layer needs more invariant tests before anything is submitted.

---

## What I need from you

**Decisions (blocking):**

1. **Framing.** The note is currently "here is a design criterion, validated retrospectively
   on one system." The alternative is to lead with the retraction + the analytic structure
   and treat the design as the constructive payoff. I lean to the second — being the paper
   that corrects a 45° folklore is more memorable — but it invites the "you got it wrong
   first" reading. Your call.

2. **Whether result 4 is enough validation** to write a manuscript, or whether you want me
   to find a second system with a published phase ceiling that wasn't used to build the
   model. I'd argue for the second system; it's the difference between a consistency check
   and a test.

**Papers that would close remaining risk:**

3. **Lucassen-Reynders** monograph or review, **Kretzschmar/Miller/Fainerman**, **Rusanov**.
   Still the four places where `45° + |arg B|/2` and the ceiling inversion might already be
   written down. Until checked, no priority claim anywhere.

4. **Any oscillatory dilatational dataset reporting |E| *and* phase as numbers** — not just
   curves. This is what would test result 3 on real data rather than simulation.

---

## Honest position

The design computation is solid and new as far as I can tell. The analytic structure is
correct but possibly not novel. The validation is one system and not blind. The whole thing
sits downstream of a 2026 paper that staked out the diagnosis, and it only works as "the
constructive sequel" — which is a real position, but a narrower one than where we started
this session.

I would not start writing until items 3 and 4 above are resolved.
