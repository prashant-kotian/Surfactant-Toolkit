# Literature survey and gap analysis — Paper 1

**Date:** 4 Nov 2026
**Purpose:** establish what is already known, what is already claimed, and therefore what this paper
can honestly claim. This document supersedes all earlier novelty assessments in this repository,
including §18.9 (review verdicts) and the novelty paragraph in `PAPER1_SUBMISSION.md`.

---

## 0. How this was done, and its limits

**Method.** Fourteen targeted queries to the Crossref REST API (`api.crossref.org/works`), plus
OpenAlex, on 4 Nov 2026. Queries covered five clusters: (i) micellar relaxation kinetics and
Aniansson–Wall theory; (ii) dilational viscoelasticity of micellar solutions; (iii) LAOD /
large-amplitude nonlinear interfacial rheology; (iv) instrument and tensiometry papers; (v) inverse
problems and identifiability in surface rheology. Full title, journal, year and DOI were retrieved
for each hit.

**Limits, stated plainly.** This is a *foundation*, not a complete systematic review. It is
query-driven, so it inherits the ranking behaviour of Crossref's relevance search; a paper not
returned for any of my fourteen phrasings is invisible to it. Citation-graph traversal (forward and
backward) has **not** been done and is the single biggest gap in this survey. Everything below
should be read as "these are the papers that must be read first", not "this is the field".

**Nothing here has been read yet.** All judgements are from titles, abstracts where Crossref
carries them, and metadata. Claims about what a paper argues are provisional until the PDF is in
hand. Two items (§3.1, §3.2) carry abstracts and are judged on those.

---

## 1. The five literatures this paper sits between

The work is interdisciplinary in a way that has not been faced properly. It borrows from five
separate communities, and **each one contains a result we currently claim or rely on**:

| # | Literature | What we take from it | Where we stand |
|---|---|---|---|
| L1 | Aniansson–Wall micelle kinetics | τ₂ as a two-time-scale relaxation | we use it, we do not advance it |
| L2 | Lucassen–van den Tempel + successors | the diffusion-controlled modulus and its phase | the 45° bound is implicit here (1972/1983) |
| L3 | Fainerman/Miller/Joos/Aksenenko school | micellar dilational data, Joos fast-kinetics model | **they have measured our signature** |
| L4 | Nonlinear/LAOD interfacial rheology | large-amplitude behaviour, harmonic decomposition | our ε = 0.20 operating point lives here |
| L5 | Inverse problems / identifiability | parameter recovery, model discrimination | **a 2026 paper pre-empts part of our claim** |

The novelty question is not "is this new to surfactant science" but "is this new to *all five*". So
far we have only ever checked L1–L3 partially, and never L4 or L5 at all. That is the root cause of
the loop the work has been stuck in.

---

## 2. MASTER DOI LIST — read these

Priority: **P0** = changes a claim, read before anything else. **P1** = must be read and cited.
**P2** = background/context. **DATA** = contains data we need.

### 2.1 P0 — these two may pre-empt us (found today, not previously known)

| DOI | Reference | Why P0 |
|---|---|---|
| `10.1002/jsde.70095` | Al-Sou'od, K. A.; Ghanem, R. "Protocol-Induced False Identifiability in Dynamic Surface Tension and Dilational Measurements of Surfactant Interfaces." *J. Surfactants Deterg.*, published online 23 Aug 2026. | **Directly adjacent to our §20 null locus and §22 blind-inversion result.** See §3.1. |
| `10.1021/acs.jpcb.5c07514` | Al-Sou'od, K. A. "Entropy-Regulated Dynamics of Surface Tension: Coupling Langmuir Kinetics with High-Precision Tensiometry." *J. Phys. Chem. B* **2026**, *130*(5), 1730–1742. | Langmuir kinetics + high-precision tensiometry + GUM/ISO uncertainty. Overlaps our §2 and our precision framing. See §3.2. |

### 2.2 P1 — the competing theory we have never cited

| DOI | Reference | Note |
|---|---|---|
| `10.1016/j.cis.2005.09.002` | Danov, K. D.; Kralchevsky, P. A. et al. "Mass transport in micellar surfactant solutions: 1. Relaxation of micelle concentration, aggregation number and polydispersity." *Adv. Colloid Interface Sci.* **2006**. | The rigorous theory. We must position against this. |
| `10.1016/j.cis.2005.09.003` | Danov, K. D.; Kralchevsky, P. A. et al. "Mass transport in micellar surfactant solutions: 2. Theoretical modeling of adsorption at a quiescent interface." *Adv. Colloid Interface Sci.* **2006**. | The adsorption side of the same theory. |
| `10.1006/jcis.1996.0537` | "Adsorption from Micellar Surfactant Solutions: Nonlinear Theory and Experiment." *J. Colloid Interface Sci.* **1996**. | Nonlinear theory — directly comparable to our large-amplitude treatment. |
| `10.1016/0166-6622(91)80280-2` | "The kinetics of the surface tension of micellar surfactant solutions." *Colloids Surf.* **1991**. | Earlier kinetics formulation. |
| `10.1063/1.2036968` | "Kinetics of nonionic surfactant adsorption at a fluid-fluid interface from a micellar solution." *J. Chem. Phys.* **2005**. | Nonionic specifically — our class. |
| `10.1016/j.colsurfa.2005.10.090` | "Micellar surfactant solutions: Dynamics of adsorption at fluid interfaces subjected to stationary expansion." *Colloids Surf. A* **2006**. | Expansion, not oscillation, but same physics. |

### 2.3 P1 — Aniansson–Wall, correctly cited

| DOI | Reference |
|---|---|
| `10.1021/j100603a016` | Aniansson, E. A. G.; Wall, S. N. *J. Phys. Chem.* **1974**, *78*, 1024–1030. **The correct DOI for our ref 3.** |
| `10.1021/j100575a019` | Aniansson, E. A. G.; Wall, S. N. *J. Phys. Chem.* **1975**, *79*, 857. |
| `10.1002/bbpc.19780820961` | Aniansson, E. A. G. "Theory of Micelle Formation Kinetics." *Ber. Bunsenges. Phys. Chem.* **1978**, *82*, 981–988. **The two-time-scale derivation** — this is the primary source for τ₁/τ₂, not the 1974 paper. |
| `10.1021/j100412a070` | "Kinetics of mixed and ionic micelles. A discussion of the slow process and a general formulation of the slow relaxation time." *J. Phys. Chem.* **1986**. **The closed form for τ₂.** |
| `10.1021/j150555a005` | "On the Statistical Mechanical Theory of Micelle Formation in Detergent Solutions." *J. Phys. Chem.* |

### 2.4 P1 — the 45° bound: where it actually comes from

| DOI | Reference | Note |
|---|---|---|
| `10.1016/0009-2509(72)80104-0` | Lucassen, J.; van den Tempel, M. *Chem. Eng. Sci.* **1972**, *27*, 1283–1291. | The diffusion-controlled modulus. |
| `10.1016/0001-8686(83)87004-3` | van den Tempel, M.; Lucassen-Reynders, E. H. "Relaxation processes at fluid interfaces." *Adv. Colloid Interface Sci.* **1983**, *18*, 281–301. | **The most likely primary statement of the 45° limit.** Read first for §5.6. |
| `10.1016/s0927-7757(01)00483-6` | Lucassen-Reynders, E. H.; Cagna, A.; Lucassen, J. "Gibbs elasticity, surface dilational modulus and diffusional relaxation in nonionic surfactant monolayers." *Colloids Surf. A* **2001**, *186*, 63–72. | **This is our ref 6, but our citation has the wrong title and year.** Nonionic. |
| `10.1016/0021-9797(72)90373-6` | Lucassen, J.; van den Tempel, M. "Longitudinal waves on visco-elastic surfaces." *J. Colloid Interface Sci.* **1972**, *41*, 491–498. |
| `10.1016/0021-9797(91)90260-f` | "A remark on the definition of the surface dilational modulus." *J. Colloid Interface Sci.* **1991**. | Definitional care — directly relevant to our δ convention. |

**Analytical note.** In LvT, tan φ = ξ/(1+ξ) with ξ = √(ω_D/2ω), so φ → 45° as ω → 0 and → 0°
as ω → ∞. **The 45° ceiling is therefore a *limit of a 1972 model*, derivable in two lines, not a
new bound.** Our paper must present it as "following from LvT" and claim only its use as an
inversion screen. This is consistent with what we already scoped, but now we can cite the actual
derivation rather than assert it.

### 2.5 P1 — micellar dilational data (L3)

| DOI | Reference | Tag |
|---|---|---|
| `10.1021/la704058y` | Fainerman, V. B.; Petkov, J. T.; Miller, R. "Surface Dilational Viscoelasticity of C₁₄EO₈ Micellar Solution…" *Langmuir* **2008**, *24*, 6447–6452. | **DATA — phase angle >45° up to 70× CMC.** |
| `10.3390/colloids5030039` | Fainerman et al. "A Multistate Adsorption Model for C₁₄EO₄ and C₁₄EO₈." *Colloids Interfaces* **2021**, *5*, 39. | **DATA — open access**, Fig. 12 frequency sweeps. |
| `10.1007/s11743-017-2016-y` | Fainerman et al. "Surface Tension and Adsorption Studies by Drop Profile Analysis Tensiometry." *J. Surfactants Deterg.* **2017**. | **Open access review**; good for isotherms + D values. |
| Mys, A. V.; Fainerman, V. B.; Aksenenko, E. V.; Miller, R. "Adsorption layer characteristics of Triton surfactants: 4. …micellar solutions." *Colloids Surf. A* **2009**. | **DATA — 5 Tritons, CMC→10×; DOI not captured, find via author search.** |
| `10.1016/s0927-7757(98)00373-2` | "Dilational viscoelasticity of surfactant monolayers." *Colloids Surf. A* **1998**. | background |
| `10.1016/j.cocis.2024.101849` | "Dynamic dilational viscoelasticity of surfactant layers at liquid–liquid interfaces." *Curr. Opin. Colloid Interface Sci.* **2024**. | current review |
| `10.1016/j.cocis.2021.101538` | "Salt effects on the dilational viscoelasticity of surfactant adsorption layers." **2022**. | ionic context |

### 2.6 P1 — LAOD / nonlinear (L4) — we have never engaged this literature

| DOI | Reference | Note |
|---|---|---|
| `10.1016/j.cis.2025.103625` | **Risse, K.; Yang, J.; de Groot, A.; Giménez-Ribes, G.; Shen, P.; Drusch, S.; Hinderink, E. B. A.; Sagis, L. M. C.** "Advances in large amplitude oscillatory dilatational surface rheology — A review." *Adv. Colloid Interface Sci.* **2025**, *345*, 103625. CC-BY. | **This is our ref 5. Our citation is WRONG: we have vol. 325, p. 103050. Correct is vol. 345, 103625.** Reviewer's "ref 5 incomplete" now resolved. |
| `10.1016/j.cocis.2014.09.003` | Sagis, L. M. C.; Fischer, P. "Nonlinear rheology of complex fluid–fluid interfaces." *Curr. Opin. Colloid Interface Sci.* **2014**. | the LAOD framework |
| `10.1122/1.2970095` | Ewoldt, R. H.; Hosoi, A. E.; McKinley, G. H. "New measures for characterizing nonlinear viscoelasticity in LAOS." *J. Rheol.* **2008**. | Chebyshev/harmonic measures |
| `10.1016/j.jcis.2023.??` | "Surface stress decomposition in large amplitude oscillatory interfacial dilatation of complex interfaces." *J. Colloid Interface Sci.* **2023**. | GSD method |
| `10.3390/colloids2040057` | "Effect of Amplitude on the Surface Dilational Visco-Elasticity of Protein Solutions." *Colloids Interfaces* **2018**. | **amplitude sweeps — directly relevant to our ε = 0.20 choice** |

### 2.7 P1 — instruments and the frequency horizon

| DOI | Reference | Note |
|---|---|---|
| `10.1088/0957-0233/27/6/065301` | "Automated system for measuring the surface dilational modulus of liquid–air interfaces." *Meas. Sci. Technol.* **2016**, *27*, 065301. | **Our ref 4 reads "27, 045301" — verify; the automated-system paper is 065301.** |
| `10.1006/jcis.2002.8665` | "Measurement of the Surface Dilational Viscoelasticity of Adsorbed Layers with a Capillary Pressure Tensiometer." *J. Colloid Interface Sci.* **2002**. | **CPT frequency range — our horizon argument depends on this** |
| `10.1021/la901391y` | "Optical Observation of High-Frequency Drop Oscillations … Capillary Pressure Tensiometry." *Langmuir* **2009**. | high-frequency limit |
| `10.1007/bf02870389` | FAST/STS-107: dilational modulus 0.01–0.32 Hz. *Colloid Polym. Sci.* **2006**. | independent frequency-range evidence |
| `10.1016/j.cis.2020.102295` | "Single bubble and drop techniques for characterizing foams and emulsions." **2020**. | instrument review |

### 2.8 P2 — surfactants, CMC, parameters

| DOI | Reference |
|---|---|
| `10.1246/bcsj.78.604` | "Micelle Formation and Surface Adsorption of Octaethylene Glycol Monoalkyl Ether (CₙE₈)." *Bull. Chem. Soc. Jpn.* **2005**. **C₁₂E₈ CMC source — we need this.** |
| `10.1134/s1061933x13040145` | "Diffusion-controlled adsorption kinetics … heptaethylene glycol monododecyl ether" (C₁₂E₇). **2013**. Nearest-neighbour to C₁₂E₈ kinetics. |
| `10.1021/jp970746s` / `10.1021/jp960148y` | Fainerman reorientation model, parts 1 and 2. **The isotherm model for ethoxylates.** |
| `10.1021/la980956b` | "Effect of the Reorientation of Oxyethylated Alcohol Molecules … on Equilibrium and Dynamic Surface Pressure." **1999**. |
| `10.1016/j.colsurfa.2015.06.015` | "Dynamic surface tension of C₁₀EO₈ at the aqueous solution/hexane vapour interface." **2015**. |
| `10.3390/molecules26195894` | Czakaj et al. — our existing LAE source (Paper 2). |

### 2.9 P2 — identifiability methodology (L5)

| DOI | Reference |
|---|---|
| `10.1021/jp049682t` | "Surface Dilational Modulus or Gibbs' Elasticity of Protein Adsorption Layers." **2004**. |
| `10.1016/j.cocis.2024.101877` | referenced by both Al-Sou'od papers; likely a measurement-practice review. Find and read. |
| `10.1021/acs.langmuir.1c00668` | cited by both Al-Sou'od papers; relevant to their framework. |
| `10.1021/acs.langmuir.2c03511` | ditto. |

---

## 3. Gap analysis

### 3.1 THREAT: `10.1002/jsde.70095` overlaps our newest result

Published **23 August 2026**. Abstract (verbatim from Crossref):

> "Dynamic surface tension and interfacial dilational measurements are used to assign adsorption,
> diffusion, and structural-relaxation mechanisms in surfactant systems. However, small unmeasured
> variations in area forcing or baseline stability **can generate identifiable relaxation modes**. We
> quantify this failure mode, termed **protocol-induced false identifiability** … With multi-protocol
> forcing, synthetic estimation recovered the four M2 dynamic parameters with **median absolute
> relative errors of 0.34%–2.02%** … Experimental falsification used two independent digitizations of
> a published polyoxyethylene (23) lauryl ether step-relaxation trace … Because the commercial
> ethoxylate and block-copolymer datasets lack batch-resolved compositional characterization
> sufficient for microscopic attribution, their reanalysis is restricted to model-discrimination
> screening rather than compound-specific mechanism assignment."

Read that alongside our own results:

| their finding | our finding |
|---|---|
| relaxation modes can be **artifacts** of protocol variation | δ sensitivity passes through zero on a physical null locus (§20.2) |
| model discrimination via AICc between nested M0/M1/M2 | model comparison by residual structure |
| median abs rel error 0.34–2.02% on synthetic recovery | 0.23–0.66% on synthetic recovery (§22) |
| "locked cross-protocol prediction" as validation standard | multi-start + relSE gate (§22.4) |
| reanalysis of *published* digitised data | we have only forward-model data |
| ethoxylate datasets lack batch-resolved characterisation | same limitation applies to us |

**What this means.** The *headline of our §22* — "a method can return an identifiable-looking answer
that is not real" — is **not new as of August 2026**. They got there first, with a more developed
methodological framework (nested models, AICc, leave-one-segment-out, locked prediction) and applied
it to real published data.

**What may still be ours.** Their failure mode is *protocol-induced* (area forcing / baseline drift).
Ours is *physical* — a genuine zero-sensitivity locus that exists in the forward map for any
diffusion-plus-reservoir system, independent of instrument error. That is a different mechanism. Our
§20.2 sign flip and its C_m-dependent location are not in their abstract. But **this distinction is
the entire novelty of our §20–22, and it now rests on a paper I have only read the abstract of.**

**Action:** this is the single most urgent read. If they have the zero-sensitivity locus, our §20–22
reduce to a worked example of their principle, and the paper's contribution collapses to τ₂ inversion
plus the LAOD angle.

### 3.2 THREAT: `10.1021/acs.jpcb.5c07514` overlaps our modelling and precision framing

Same author group, *J. Phys. Chem. B* **2026**, *130*, 1730–1742. Title alone: "Coupling Langmuir
Kinetics with High-Precision Tensiometry". Our §2 is "Coupled Reaction–Diffusion Model with Langmuir
Interfacial Kinetics". Their reference list contains GUM (JCGM 100:2008), ISO 5725-2, ISO/IEC 17025,
ISO 5725 — i.e. a **formal metrological uncertainty framework**. Our uncertainty treatment is a
Jacobian SE plus a coverage test we ran last week.

**What this means.** Two risks: (a) they may already couple Langmuir kinetics to precision claims, and
(b) their metrological rigour sets a bar our current uncertainty section does not meet. If a referee
knows this paper, "0.106% precision" without ISO 5725 repeatability/reproducibility treatment looks
naive.

### 3.3 GENUINE GAP: the 45° bound is stated nowhere I can find, but is derivable

I searched specifically for a statement of the 45° ceiling and found **no paper that states it as a
theorem**. It is a limit of the LvT expression (φ → 45° as ω/ω_D → 0). So:

- We **can** claim the screen, because nobody appears to have used it as an inversion test.
- We **cannot** claim the bound as new. It is a two-line consequence of a 1972 result and should be
  presented as such, with the derivation shown.
- The honest phrasing is: "implicit in [LvT 1972; van den Tempel & Lucassen-Reynders 1983], used here
  as a screening criterion, to our knowledge for the first time."

### 3.4 GENUINE GAP: τ₂ inversion from LAOD — still unclaimed

Nothing found, in any of the fourteen queries, that inverts a micellar relaxation time from
large-amplitude oscillatory dilatational data. The LAOD literature (L4) is entirely proteins and
particles; the micellar literature (L3) is entirely small-amplitude, linear, and analysed with
Joos/Danov–Kralchevsky equilibrium-kinetics models. **The intersection is empty.** This survives as
the paper's core claim.

But it is now a *narrow* claim: one parameter, one model class, synthetic data only.

### 3.5 GENUINE GAP: the combination is unusual

Nobody combines (a) explicit micellar reservoir + (b) large-amplitude forcing + (c) formal
identifiability analysis + (d) an instrument-horizon argument. Each element exists; the combination
does not appear to. That is a legitimate contribution, but it is a *methods* contribution and must be
pitched as one.

### 3.6 What the survey says we should NOT claim

| Claim in current manuscript | Status after survey |
|---|---|
| 45° ceiling novel | **No.** Derivable from LvT 1972; cite, derive, do not claim. |
| "no scaling collapse exists" | superseded by our own §21 (three-group) — reword |
| precision 0.106% as a headline result | **weakened.** Al-Sou'od uses GUM/ISO; our SE is a Jacobian estimate |
| false identifiability as a novel warning | **pre-empted** (§3.1) unless we can distinguish physical null locus |
| micellar nonionics exceed 45° | **pre-empted**, Fainerman 2008 |
| ωτ₂ ≈ 3 optimum | observation only, no derivation — keep as observation |
| τ₂ inversion from LAOD | **survives** |

---

## 4. Dependencies we have been treating as facts

| Fact used in paper | Actual provenance | Status |
|---|---|---|
| τ₂(C₁₂E₈) = 4 s | **unknown** — a number in our benchmark matrix with no traced source | **must be traced or dropped** |
| C₁₂E₈ CMC = 0.09–0.109 mM | secondary (0.109 mM from an Abbott webpage in one search) | needs a primary source: `10.1246/bcsj.78.604` |
| PAT-1M range 0.005–1.0 Hz | instrument specs, not literature | cite the instrument, not a paper |
| CPT range up to 500 Hz | our ref 4, whose citation appears wrong | **verify against `10.1088/0957-0233/27/6/065301`** |
| E₀ anchors 25/45.3/103.5 mN/m | Czakaj Table 4 — **different surfactant (LAE)** | **not C₁₂E₈ values.** Defensible as a range, not as C₁₂E₈. |
| Γ∞, K for C₁₂E₈ | our own fit | needs a literature check against Fainerman reorientation parameters |

**The τ₂ = 4 s problem is serious.** It is the number the entire paper is calibrated around, and after
five rounds of work we still cannot say where it came from. If it is unverifiable, §3 is built on a
guess and the reviewer's item 2 ("τ₂ = 4 s is a literature benchmark") is more damaging than we
thought.

---

## 5. Read order

**Wave 1 — decide whether the paper survives (read first, in this order):**
1. `10.1002/jsde.70095` — Al-Sou'od & Ghanem 2026. Does their false identifiability include our null locus?
2. `10.1021/acs.jpcb.5c07514` — Al-Sou'od 2026. How close is their Langmuir+precision framing to ours?
3. `10.1016/j.cis.2025.103625` — Risse et al. LAOD review. Where does our ε = 0.20 sit? Is anything like our inversion in there?

**Wave 2 — fix the theory base:**
4. `10.1016/0001-8686(83)87004-3` — van den Tempel & Lucassen-Reynders 1983. The 45° derivation.
5. `10.1016/j.cis.2005.09.002` and `10.1016/j.cis.2005.09.003` — Danov & Kralchevsky. The correct
   micellar mass-transport theory we should be building on.
6. `10.1002/bbpc.19780820961` — Aniansson 1978. The proper τ₁/τ₂ derivation.
7. `10.1021/j100412a070` — the closed form for τ₂.

**Wave 3 — data and parameters:**
8. `10.1021/la704058y` — Fainerman 2008. Phase angle vs concentration.
9. `10.3390/colloids5030039` — Fainerman 2021, open access. Frequency sweeps.
10. `10.1246/bcsj.78.604` — C₁₂E₈ CMC and isotherm.

**Wave 4 — background:** everything in §2.6 and §2.9.

---

## 6. Honest bottom line

The survey found **two 2026 papers by the same author** that occupy territory we believed was ours,
and one of them (the *J. Surfactant Deterg.* paper) is close enough to our §20–22 that the claim may
not survive contact. It also found that our 45° "bound" is a 1972 corollary, that our refs 4, 5 and 6
have wrong bibliographic details, and that τ₂ = 4 s — the number everything is calibrated to — has no
traced source.

What survives: **τ₂ inversion from large-amplitude oscillatory dilatational data is genuinely
unclaimed**, and the combination of reservoir + LAOD + identifiability + instrument horizon is
unusual. That is a narrower and more defensible paper than the one currently drafted.

None of this could have been established by more simulation. It required reading, which we had not
done. The user's diagnosis was correct.
