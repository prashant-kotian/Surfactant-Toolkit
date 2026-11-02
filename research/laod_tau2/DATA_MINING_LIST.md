# COMPLETE DATA-MINING LIST — LAOD / τ₂ paper

**For:** Prashant and Antigravity · **From:** Arena AI · **Date:** 2026-11-02
**Purpose:** everything needed to turn `research/laod_tau2/` from a feasibility study into a
submittable methods paper. Items marked **[A]** are for Antigravity (paywalled/installed access);
**[U]** is you (manual download); **[X]** is blocked for me by the sandbox egress allowlist.

---

## 0. READ THIS FIRST — two findings that change what to fetch

### 0.1 The instrument audit needs correcting, and it matters a lot

Antigravity's audit said the 100 Hz tier belongs to the *dynamic capillary pressure
microtensiometer (CPT)* and that profile-analysis tensiometers top out at ~1 Hz. The literature
says something different and more useful:

- Profile-analysis (pendant drop/bubble): **0.005–1 Hz**, occasionally to ~2 Hz. ✓ audit correct.
- **CPT / oscillating bubble: 0.1–100 Hz, some instruments to 500 Hz.** ✓ audit correct.
- **Microtensiometer (Alvarez, Walker & Anna, *Langmuir* 26:13310, 2010): ~50–100 µm bubble,
  oscillated 0.001–100 cycles/min = 1.7 mHz–1.7 Hz.** ✗ **NOT a 100 Hz instrument.**

So the 100 Hz column in my Table (and in `FINDINGS.md` §4) is **CPT, not the microtensiometer**,
and the microtensiometer's real value is different and still important: it is the
**curvature/depletion** test bed (its radius is comparable to the adsorption/diffusion length),
which is a physics correction I have not implemented. This is corrected in the repo; the
mining list below asks for the two instrument papers so the numbers are cited properly.

### 0.2 If you send **vector** PDFs, I do not need to digitise by eye

I built and validated a vector-PDF extractor (`research/laod_tau2/vector_digitize.py`,
2 tests): it reads the curve geometry and the tick labels **out of the PDF file**, giving
**0.06 % median error** against ground truth versus 1–3 % for pixel clicking. It is
also reproducible from the PDF hash.

**Therefore: please prefer the publisher PDF over a screenshot, and prefer a
figure exported by matplotlib/Origin/MATLAB over a scanned one.** For the Elsevier papers
this usually means the "PDF" of the article, not the "figure image". Where only a raster
exists, say so — I will digitise pixels, but every point will carry a larger uncertainty.

---

## 1. PRIORITY A — the paper cannot be submitted without these

### A1. LAE oscillatory dataset (the first experimental validation) — **[U]** or **[A]**
**Czakaj, Jarek, Krzan, Warszyński, *Molecules* 26(19):5894 (2021), CC BY 4.0.**
PMC8512375 · doi:10.3390/molecules26195894
- **Already retrieved:** no — blocked from my sandbox (`cdn.ncbi.nlm.nih.gov`,
  `pmc.ncbi.nlm.nih.gov`, `mdpi.com` all refuse TLS at the egress layer).
- **What I need:** the **publisher PDF** (`molecules-26-05894.pdf`), and if it exists, the
  **Supplementary Materials PDF**.
- **Which figures:** Fig. 5 (γ(t) and area oscillations at 0.2/0.5/0.8/1.0/1.5 mM, 0.01 and
  0.1 Hz), Fig. 6 (E′, E″ vs frequency with LVDT fits), **Fig. 7 (I₂/I₁ vs concentration —
  this is prior art for one of my claims, see `DATA_SEARCH_AND_CORRECTION.md` §1)**, Table 4
  (concentration, ε₀, ν_D).
- **Check on receipt:** open the PDF and zoom into Fig. 5 — if the curves stay sharp at 400 %
  they are vector; if they pixelate they are raster.
- **Why it matters:** these authors already published the second-harmonic ratio above the CMC.
  I must re-analyse their data to show what my method adds (inversion + phase lag), and to
  avoid claiming novelty I do not have.

### A2. Instrument specifications — **[A]**
Needed to replace "assume 100 Hz" with a citation.
1. **Alvarez NJ, Walker LM, Anna SL.** *A microtensiometer to probe the effect of radius of
   curvature on surfactant transport to a spherical interface.* **Langmuir 26(16):13310–13319
   (2010).** doi:10.1021/la101870m → radius, frequency range, resolution, strain capability.
2. **Alvarez NJ, Walker LM, Anna SL.** *Diffusion-limited adsorption to a spherical geometry:
   the impact of curvature and competitive time scales.* **Phys. Rev. E 82:011604 (2010).**
   → the curvature correction I have not implemented.
3. **Alvarez NJ, Walker LM, Anna SL.** *A criterion to assess the impact of confined volumes on
   surfactant transport to liquid–fluid interfaces.* **Soft Matter 8:8917 (2012).**
4. **Automated oscillating-bubble instrument, 1–500 Hz**, *Meas. Sci. Technol.* **27:065301
   (2016)** → the real upper bound for oscillating-bubble dilatational rheology.
5. **Capillary-pressure technique review** (Kovalchuk/Fainerman/Miller school, ~2012,
   *Adv. Colloid Interface Sci.* or equivalent) → the ≤100 Hz validity limit for the CPT model.

### A3. Independent τ₂ measurements for validation — **[A]**
The paper's核心 claim is that LAOD gives τ₂; it must be checked against τ₂ measured another way.
1. **Patist A, Kanicky JR, Shukla PK, Shah DO.** *Importance of micellar kinetics in relation to
   technological processes.* **J. Colloid Interface Sci. 245(1):1–15 (2002).**
   → the standard τ₂ compilation: SDS ms–s, Triton X-100 ~s, CₙEOₘ minutes; concentration and
   salt dependence; the SDS 200 mM maximum. **This is the backbone of our validation table.**
2. **Patist A, Axelberd T, Shah DO.** *Effect of long chain alcohols on micellar relaxation time
   and foaming properties of SDS solutions.* **J. Colloid Interface Sci. (1998).**
3. **A pressure-jump or ultrasonic-absorption study giving τ₂ (± error) for one of
   SDS / C₁₂E₈ / CTAB / Triton X-100** — ideally with the aggregation number from the same paper.
4. **Mysona JA, McCormick AV, Morse DC.** *Nonlinear dynamics in micellar surfactant solutions
   I & II.* **Phys. Rev. E 105:034602, 034603 (2022).** → the rigorous theory our linearised
   reservoir must be checked against. (I may reach the accepted manuscripts via OSTI; if not,
   please supply.)

---

## 2. PRIORITY B — needed to parameterise real systems

### B1. Equilibrium isotherms (SDS and CTAB have **no supplied Langmuir K** — this blocks them)
For each: surface tension vs concentration (a digitised table or a table in the paper), plus Γ∞,
K, and Frumkin `a` if given.
1. **SDS**: a Szyszkowski/Langmuir–Frumkin fit with reported Γ∞ and K at 22–25 °C, ideally with
   and without NaCl.
2. **CTAB**: same, plus a NaBr series (the audit says τ₂ changes by ~10× with NaBr — I want the
   isotherm at matched conditions).
3. **C₁₂E₈**: confirm Γ∞ = 2.8×10⁻⁶ mol/m² and **K = 140 m³/mol** — Antigravity supplied these
   without a source, and every C₁₂E₈ row depends on them. A published isotherm is needed.
4. **LAE**: the isotherm is in the A1 paper; just make sure Fig. 2/3 is readable.

### B2. Aggregation numbers (capacity derives from n̄, and capacity is the dominant lever)
1. **SDS**: n̄ vs concentration and vs NaCl (light scattering / SANS / fluorescence quenching).
2. **C₁₂E₈** and **CTAB**: n̄ at 2× and 10× CMC.
3. If n̄ is unavailable at our conditions, a **polydispersity index** (σ/n̄) is enough to bound
   the error — the Aniansson–Wall literature typically uses σ = 0.2 n̄.

### B3. Diffusion coefficients
PFG-NMR or Taylor-dispersion D for the monomer of each system. The audit's values
(5.5, 3.8, 4.0 × 10⁻¹⁰ m²/s) need a source; Stokes–Einstein estimates are not sufficient for a
methods paper.

---

## 3. PRIORITY C — context, positioning, and the referee's first three objections

1. **Risse K, Yang J, de Groot A, Giménez-Ribes G, Shen P, Drusch S, Hinderink EBA, Sagis LMC.**
   *Advances in large amplitude oscillatory dilatational surface rheology — a review.*
   **Adv. Colloid Interface Sci. 345:103625 (2025).** → forces us to state exactly what is new.
2. **Klaus S, Willenbacher N et al.** *Surface stress decomposition in large amplitude
   oscillatory interfacial dilatation of complex interfaces.* **J. Colloid Interface Sci. (2023).**
   → the 0th/2nd/4th harmonic decomposition; our phase-lag observable must be positioned against it.
3. **Lucassen J, van den Tempel M.** *Chem. Eng. Sci. 27:1283 (1972).* → the anchor. (Closed form
   already implemented; want the original for the citation and for the exact E′/E″ notation.)
4. **Morse et al. (as A3.4).**
5. **A review of micellar relaxation vs dynamic surface tension** (Patist/Shah school, or
   Eastoe & Dalton *Adv. Colloid Interface Sci.* 85:103, 2000) → for the DST bridge.

---

## 4. PRIORITY D — only if we decide to implement the curvature correction

Alvarez 2010b + 2012 (above), and the **"surfactant dynamics: hidden variables controlling fluid
flows"** JFM review (2020) which summarises the curvature/depletion literature. This is a
*second paper's* worth of work; fetch only when A1–A3 are in hand.

---

## 5. How to hand files to me (this matters — I cannot fetch them)

- Put files in the repo under **`research/laod_tau2/literature/`** (gitignored if large), or in
  the shared Drive folder and tell me the paths.
- **Naming:** `firstauthor_year_journal_short.pdf`, e.g. `czakaj_2021_molecules_LAE.pdf`.
- For each file, one line: **DOI + licence + whether the figures are vector or raster.**
  If you are unsure about vector/raster, just send it — I check automatically and report.
- **Do not crop or screenshot** unless the PDF is unavailable; a screenshot throws away 30× of
  the achievable precision.
- Supplements matter: harmonic ratios and raw tables often live only in the SI.

---

## 6. What I will produce from each item

| Item | Output |
|---|---|
| A1 (LAE) | `data/curated/LAE_Czakaj2021.json`: extracted γ(t) traces with per-point uncertainty, FFT-derived I₂/I₁ and phase lag, and a comparison against their Fig. 7 |
| A2 | `docs/instruments.md`: a citable table of window, radius, resolution, max strain, with the CPT/microtensiometer correction |
| A3 | `docs/tau2_validation.md`: recovered τ₂ vs independently measured τ₂, with error bars — **the decisive figure** |
| B1–B3 | `systems.py` entries with a source string per number; SDS and CTAB become runnable |
| C | the introduction's positioning paragraph, written so no reviewer can say "already known" |

---

## 7. Stop-list — do not spend effort on these

- WebPlotDigitizer licences or GUI workflows: superseded by the vector extractor.
- Scanning/OCR of tables: if a paper gives a table, copy the text, do not photograph it.
- The foam and microfluidic literatures: those proposals are on hold (see
  `DEEP_RESEARCH_PAPER_CANDIDATES.md` §7).
- Anything behind a paywall you cannot get: tell me it is unavailable rather than sending a
  figure lifted from a secondary source — unverifiable provenance is worse than a gap.
