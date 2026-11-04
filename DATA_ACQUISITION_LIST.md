# Experimental data acquisition list — what exists, what it can validate, how to get it

Compiled 4 Nov 2026 in response to: "search for supporting information where a complete experimental
dataset is given, so those data can validate our work."

**Sandbox limitation:** only github.com, api.github.com, codeload.github.com, registry.npmjs.org,
pypi.org and files.pythonhosted.org are reachable from this workspace. Zenodo, OSTI, ACS, Elsevier
and ResearchGate all return nothing. So this is a search-result survey, not a download. Items below
are tagged **[I DOWNLOAD]** where the user must fetch the file.

---

## 0. The honest headline first

There is **no open dataset that directly validates Paper 1's central claim** (τ₂ inverted from
large-amplitude oscillatory dilatational data of a *micellar* nonionic surfactant). That claim is
new, which is why it has not been measured. Anyone who says otherwise is wrong.

But there is a **much closer match than I expected**, and it cuts two ways:

> Fainerman, Petkov & Miller, *Langmuir* **2008**, 24, 6447–6452, doi:10.1021/la704058y.
> C₁₄EO₈, bubble profile analysis, concentrations from CMC (7 µM) up to **70 × CMC**.
> *"With increasing C₁₄EO₈ concentration the viscoelasticity modulus decreases and the phase angle
> increases. At the highest surfactant concentrations, the phase angle was more than 45°."*

That is **exactly our signature** — a micellar nonionic system whose loss angle exceeds 45° — measured
and published in 2008, eighteen years before this work. Two consequences, both must go in the paper:

1. **Validation:** the >45° behaviour is real and externally measured, not a model artefact. This is
   the strongest available corroboration and it costs nothing to cite.
2. **Novelty:** our framing must not imply that >45° in micellar nonionics is new. What is new is
   (a) using it as a *screen* with a stated model-class scope, (b) the τ₂ inversion, (c) the null
   locus. Dr. Kotian's own instruction to narrow novelty to "τ₂ inversion, not first-harmonic use"
   is independently confirmed by this paper.

Note also that Fainerman's interpretion is **Joos fast-micellar-kinetics**, whereas ours is a lagging
source term. Those are different physical pictures of the same data. That disagreement is itself
worth a paragraph, and is more interesting than either paper alone.

---

## 1. Tier 1 — data that can actually test the model

### 1.1 Fainerman, Petkov & Miller 2008 — C₁₄EO₈ micellar, up to 70× CMC

- **Citation:** Fainerman, V. B.; Petkov, J. T.; Miller, R. *Langmuir* **2008**, *24*(13), 6447–6452.
- **DOI:** 10.1021/la704058y · **PMID:** 18517233
- **What we need:** modulus **and phase angle vs frequency**, for each concentration (CMC → 70×).
- **Why it is the single best target:** it is the only paper found that reports *phase angle* (not just
  modulus) across a wide concentration range above the CMC, in a nonionic micellar system. Our model
  predicts phase angle directly, and our whole design argument is about phase, not modulus.
- **What it can validate:** the sign and magnitude of the δ increase with concentration, and whether
  δ crosses 45° where our model says it should. It cannot validate τ₂ inversion (no LAOD, and the
  frequency range is 0.005–0.2 Hz).
- **Caveat:** C₁₄EO₈ not C₁₂E₈; τ₂ differs, so only the *shapes* transfer, not the absolute τ₂.
- **[I DOWNLOAD]** — ACS paywall.

### 1.2 Fainerman et al. 2021 — C₁₄EO₈ modulus vs frequency, several concentrations (OPEN ACCESS)

- **Citation:** Fainerman, V. B.; Kovalchuk, V. I.; Aksenenko, E. V.; Ravera, F.; Liggieri, L.;
  Loglio, G.; Makievski, A. V.; Schneck, E.; Miller, R. *Colloids Interfaces* **2021**, *5*, 39.
- **DOI:** 10.3390/colloids5030039 · **PDF:** pdfs.semanticscholar.org/1bc0/d71c29cfc5a337caee5c0e321be2e6367bfa.pdf
- **What we need:** **Figure 12** — "experimental visco-elasticity modulus vs frequency for C₁₄EO₈
  solutions of different concentrations (0.2, 0.5, …)". Also **Figures 7–11** (same for C₁₀EO₈ and
  others) and **Table 1** (reorientation-model parameters incl. D).
- **Why:** open access, so digitally extractable without a paywall, and includes *frequency sweeps* at
  multiple concentrations — the exact structure of our Protocol B dataset. Table 1 also gives a
  literature D for C₁₄EO₈ (8.7×10⁻⁹ m²/s in the model fit; expected ~3.5×10⁻¹⁰), which we need for
  the third similarity group.
- **No download needed** — open access, I can digitise the figures directly if you confirm.

### 1.3 Mys, Fainerman, Aksenenko & Miller 2009 — Triton micellar, CMC to ~10× CMC

- **Citation:** Mys, A. V.; Fainerman, V. B.; Aksenenko, E. V.; Miller, R. "Adsorption layer
  characteristics of Triton surfactants: 4. Dynamic surface tension and dilational visco-elasticity
  of micellar solutions." *Colloids Surf. A* **2009**.
- **What we need:** viscoelasticity modulus **and phase angle vs frequency**, CMC → ~10× CMC, for
  Triton X-100 / X-114 / X-165 / X-305 / X-405.
- **Why:** covers the concentration window we operate in (1.2–3× CMC) and includes phase angle. Triton
  X-100 τ₂ ≈ 3.5 s is already in our benchmark matrix, so it is the closest dynamical analogue to
  C₁₂E₈'s 4 s.
- **Benefit:** five surfactants in one paper = a real multi-system test of the E₀-independence claim.
- **[I DOWNLOAD]** — Elsevier paywall.

### 1.4 Dynamic Surface Tension Database (DST-SurfDB) — OPEN, but probably not useful

- **Zenodo:** 10.5281/zenodo.21118257 · single JSON, 1.5 MB · 979 series, 233 surfactants, 1959–2024.
- **Fields:** SMILES, temperature, concentration, surface age, tension, method, DOI, Hua–Rosen fit
  parameters, and **`t_derivatives`** (dγ/d ln t).
- **Why it is tempting:** `t_derivatives` is the diffusive short-time slope, and we could check our D
  against independent adsorption kinetics.
- **Why it probably fails for us:** it is **dynamic surface tension of quiescent adsorption only** —
  no oscillatory, no phase angle, no modulus. It cannot test the LAOD claim at all. It could at best
  give a sanity check on D and K for C₁₂E₈ below/at the CMC.
- **Verdict:** low priority. **[I DOWNLOAD]** only if we want the D/K cross-check; 1.5 MB JSON.

### 1.5 Reichert et al. — C₁₂E₈ *large-amplitude* oscillatory, but submicellar

- **Citation:** Reichert, M. D.; Alvarez, N. J.; Brooks, C. F.; Grillet, A. M.; Mondy, L. A.;
  Anna, S. L.; Walker, L. M. "The Importance of Experimental Design on Measurement of Dynamic
  Interfacial Tension and Interfacial Rheology in Diffusion-Limited Surfactant Systems."
- **OSTI:** osti.gov/servlets/purl/1140848 (blocked from sandbox; try from a browser)
- **What it contains:** C₁₂E₈ **6.0 µM (submicellar, ~0.06× CMC)**, pendant drop R=1 mm and
  microtensiometer R=30 µm, E* vs ω, plus **12.5% amplitude** oscillations at 0.1–1 Hz.
- **Why it matters despite being submicellar:** it is the only C₁₂E₈ *large-amplitude* dilatational
  dataset found, and its §"startup effects" discussion — that large amplitude + poor startup gives
  erroneous elasticity — is the same failure mode our §3.4 nonlinearity budget worries about. It is a
  methods citation more than a validation source.
- **Cannot validate:** anything micellar. Explicitly a diffusion-limited (sub-CMC) study.
- **[I DOWNLOAD]** — check if OSTI full text is open.

---

## 2. Tier 2 — useful for parameter grounding, not validation

| Item | Identifier | Use | Access |
|---|---|---|---|
| SurfPro2 surfactant property DB | Zenodo 10.5281/zenodo.21456552; github.com/BigChemistry-RobotLab/SurfPro2 | CMC, cross-check our K/E₀ anchors | open |
| Surfactant adsorption to fluid interfaces (SDS/DTAB/Brij35) | Zenodo 10.5281/zenodo.18799416 | surface excess, Γ∞, min area/molecule | open |
| Interface Science Data Repository (ISDR) templates | Zenodo 10.5281/zenodo.19680999 | **deposit our own data in a citable format** | open |
| Czakaj et al. 2021 (LAE) | already in repo, `lae_data.py` | Paper 2, not this one | have it |

The ISDR templates are worth noting for a different reason: depositing our forward-model dataset and
analysis code there would answer the reviewer's items 10 and 11 (private/unauditable code, no
manifest) in one move, and give the paper a citable data DOI.

---

## 3. What I searched for and did NOT find

Recorded so nobody repeats the search:

- **No** open dataset of LAOD (large-amplitude) dilatational rheology on any micellar surfactant.
  The LAOD literature that exists (Lissajous / general stress decomposition) is on **proteins**
  (WPI, pea albumin, β-casein, Escin) and **nanoparticles**, not micellar surfactants.
- **No** open τ₂ values for C₁₂E₈ measured by any method — our 4 s is a literature benchmark we
  cannot trace to a primary measurement in this workspace.
- **No** Zenodo/figshare deposit of raw γ(t) traces from any oscillating-bubble dilatational study.
  The field does not yet share raw waveforms; it shares fitted E and φ values.
- Accessed-but-empty: zenodo.org/communities/isdr search for "dilatational rheology" (returned no
  records listing; the ISDR is mostly adsorption/isotherm data).

---

## 4. Recommended order, with expected payoff

1. **§1.2 (open access, no cost)** — digitise Fig. 12 frequency sweeps for C₁₄EO₈. Gives a genuine
   multi-concentration frequency dataset to run Protocol B against. Can be done in this workspace if
   you confirm the figure.
2. **§1.1 (download)** — the phase-angle-vs-concentration data. This is the only route to an
   *external* test of δ > 45°, and it simultaneously fixes the novelty framing.
3. **§1.3 (download)** — five Tritons, CMC→10× CMC with phase angle. Multi-system test of
   E₀-independence, the one surviving claim from the collapse work.
4. **§1.4/§2** — parameter cross-checks and data deposition. Low effort, small payoff.

**The decisive point:** none of these validate τ₂ inversion, because nobody has published LAOD on a
micellar surfactant. What they *can* do is validate the δ > 45° signature externally and test the
design window on other systems. That converts the paper from "model talking to itself" to "model,
plus independent confirmation of its signature, plus a validated design rule" — which is a real
improvement, but it is still not the measured C₁₂E₈ sweep that §7.7 names as the key missing
experiment.
