# Data acquisition status, and a correction to my own novelty claim

Date: 2026-11-02. Companion to `FINDINGS.md`. Two things happened that the next
person must know before they write a sentence of the manuscript.

---

## 1. Correction: Figure 7 of the LAE paper is prior art for one of my claims

**Czakaj, Jarek, Krzan, Warszyński, *Molecules* 26(19):5894 (2021)**, CC BY 4.0,
open access, PMC8512375.

They measured the oscillatory surface tension of ethyl lauroyl arginate (LAE) by
pendant drop (PAT-1M, 11 µL drop, 22 °C, oscillation <10% of volume, Fourier
transform of the τ variations, precision 0.1 mN/m) and **Figure 7 plots "the
ratio of second to first harmonics in the Fourier spectrum of surface tension
oscillations" against concentration at 0.01 and 0.1 Hz.** Their conclusion: the
nonlinear response appears near and above the CMC and "nonlinearity originates
from the presence of micelles".

This matters, and I want it stated plainly:

- **My earlier claim that "nobody has treated the harmonic spectrum of a LAOD
  response as a quantity of interest" is too strong.** The second-harmonic ratio
  above the CMC has been measured and published, in an open-access paper, with a
  correct qualitative attribution to micelles.
- **What remains genuinely open is narrower and must be stated as such:** they
  report I₂/I₁ as a *diagnostic of nonlinearity onset* versus concentration.
  They do not (a) invert it for a relaxation time, (b) sweep frequency finely
  enough to locate the ωτ₂ ≈ 1 signature, (c) separate the micellar reservoir
  from equation-of-state nonlinearity, or (d) report the phase lag alongside the
  harmonic ratio. Those four things are the paper. The novelty claim in
  `DEEP_RESEARCH_PAPER_CANDIDATES.md` §2 should be read as superseded by this
  paragraph.

Additional value in that paper for us: their **Table 4** gives (concentration,
ε₀, ν_D) — including ν_D = 0.0032–0.06 Hz — and Figure 6 fits E′, E″ to
Lucassen–van den Telt. That is an independent test case for our LVDT anchor and
for the sign of the nonlinearity.

---

## 2. Digitisation is blocked from this sandbox — a real infrastructure finding

Attempted, and failed, at the egress layer:

| Route | Result |
|---|---|
| `curl`/`urllib` → `cdn.ncbi.nlm.nih.gov/.../molecules-26-05894-g007.jpg` | `SSL_ERROR_SYSCALL` (host not on the allowlist) |
| `urllib` → `pmc.ncbi.nlm.nih.gov/.../pdf/molecules-26-05894.pdf` | TLS closed by egress |
| `urllib` → `mdpi.com/.../pdf` | TLS closed by egress |

**Consequence:** the plan in the operational roadmap that assigns Arena the job of
"downloading papers and extracting (x, y) coordinates from PDF figures
(WebPlotDigitizer routines)" **cannot be executed from this sandbox as configured**.
The allowlist is GitHub/PyPI/npm only.

Two workable options, and one is strictly better:

1. *(works today)* Antigravity or you download the PDF/figure and drop it in the
   workspace; I do the axis calibration, tracing, uncertainty and provenance
   bookkeeping, which is the labour-intensive part.
2. *(needs a platform change)* Add `cdn.ncbi.nlm.nih.gov`, `pmc.ncbi.nlm.nih.gov`
   and `mdpi.com` (or a general HTTPS allowlist) to the sandbox's egress policy.
   I would not spend effort on this for one paper.

Until one of those happens, **no digitised dataset exists**, and no claim in the
manuscript should imply otherwise.

---

## 3. What I did complete in this round

- **Capacity closure verified, not assumed.** `C_m = (c_tot−CMC)/CMC·n̄` is
  *exactly* the mass-action result in the monodisperse, pinned-monomer limit:
  `C_m = n̄·φ·c_tot/c₁` with `c₁ = CMC` reduces to the closure algebraically
  (`systems.py::verify_capacity_closure`, now a test). Two corrections added:
  a micellar volume-fraction term (≤10% at 10×CMC) and an explicit statement that
  the closure **overestimates near the CMC edge** (the exact derivative → 0 there),
  which is why every run uses c_tot ≥ 2×CMC.
  *Correction to my own earlier note:* I initially reported a 2000× discrepancy
  against the library's stepwise mass-action solver. That was **my** error — I fed
  the *adsorption* Langmuir K (140 m³/mol) into a *micellization* ladder, which are
  different constants. The closure is sound; the cross-check was invalid.
- **Real-system parameter sets** (`systems.py`): SDS, C₁₂E₈, CTAB with CMC, n̄, Γ∞,
  D, τ₁/τ₂ ranges from Antigravity. **Only C₁₂E₈ has a supplied Langmuir K.** SDS
  and CTAB are Frumkin systems with no fitted K — so those rows are parameterised
  but not yet runnable, and a test enforces that they are flagged as such.
- **The map for real systems** (Figure 1): with real τ₂ ranges, the conclusion
  inverts the plan — see below.
- 14 tests passing in `research/laod_tau2/tests/`.

---

## 4. The real-systems result changes the experimental recommendation

Applying the bandwidth rule (usable ω_max·τ₂ ≳ 0.3, precise ≳ 3) to the real
parameter table:

| System | τ₂ | 1 Hz pendant drop (PAT-1M) | 100 Hz capillary pressure |
|---|---|---|---|
| SDS, no salt | 2–40 ms | ✗ (ωτ₂ = 0.01–0.25) | borderline → usable at the upper end |
| SDS, 0.1 M NaCl | 1–5 s | ✗ (needs ω_min ≲ 0.05 s⁻¹) | ✓ |
| C₁₂E₈ | 1–10 ms | ✗ | ✓ (precise ≥5 ms) |
| CTAB, pure | 50–500 ms | ✓ at 500 ms, edge at 50 ms | ✓ |
| CTAB + NaBr | 0.5–5 s | ✓ | ✓ |

**CTAB is the only system a 1 Hz pendant drop can attack, and only at its slow
end.** For SDS and C₁₂E₈ the measurement requires a capillary-pressure
microtensiometer. That is the opposite of the original plan, which assumed the
PAT-1M data in the literature would do.

Practical consequence: the first experimental target should be **CTAB**, where
existing PAT-1M-class measurements can already bracket ωτ₂ ≈ 1, with SDS as the
demonstration that the method *fails* on the 1 Hz instrument — a negative control
that costs nothing and pre-empts the obvious referee objection.

---

## 5. Still open, in priority order

1. **Capacity is huge for real systems** (62–855), not the 5 used in the
   feasibility study. Every number in `FINDINGS.md` §4–5 must be re-derived at the
   real capacity before it is quoted. The strong-reservoir limit may also change
   the sensitivity qualitatively — that is the next computation, not a footnote.
2. **Langmuir/Frumkin K for SDS and CTAB** must be fitted from a published
   isotherm (I can do the fitting if the isotherm data is supplied).
3. **One experimental comparison** — LAE (CC BY, above) is the best first target
   because the data are open and the authors already report I₂/I₁.
4. The **100 Hz bandwidth figures** come from Antigravity's instrument audit and
   should be cited to the instrument specification, not to this repo.
