# Where the next paper should go — a gap analysis, with its own uncertainty stated

**Date:** 5 November 2026
**Status:** decision document. Not a manuscript, not a plan of record.
**Why this exists:** the LAOD/τ₂ paper was lost because novelty was inferred from keyword
absences rather than from reading the field. This document states what was actually searched,
what was not, and what the evidence supports.

---

## 0. What this survey is, and is not

**Done.** Keyword and phrase searches across six problem areas (listed below); reading of
abstracts and, where accessible, full text; four targeted novelty-collision checks on the
shortlisted candidate; one quantitative feasibility probe run in code (`probe_mechanisms.py`);
inspection of what SurfactantKit actually implements.

**Not done.** No systematic citation-graph traversal. No backward/forward snowballing from the
key primary papers (Lucassen 1972, Joos, Danov/Kralchevsky, Fainerman, Salager, Risse 2025).
No Web of Science / Scopus citation counts. No reading of paywalled methods sections. No
search of theses, patents, or conference abstracts, where negative and methodological results
often live.

**Therefore.** Every "no one has done X" below should be read as "I did not find X, and I looked
in these specific places". That is exactly the claim that failed last time. The mitigation is
that each candidate is paired with a **cheap falsification test** — a search or a 1-day check
that would kill it — and those are listed so they can be run before any more energy is spent.

---

## 1. The landscape: what is crowded, and should be avoided

| Area | State of the field | Verdict |
|---|---|---|
| CMC prediction | Very crowded. GNN/GCN (RMSE 0.30 → 0.21), graph transformers (R² 0.943–0.947, Jan 2026), COSMOmic, QSPR, group contributions, ML on 1,638-molecule sets. 2026 papers already claim state of the art. | **Avoid.** Marginal gains only, and the benchmark protocol problem is somebody else's paper. |
| HLD-NAC phase behaviour prediction | Active but consolidating; the 2025 book chapter and Abbott's guide cover formulation practice; ML-on-HLD-NAC already published (Appl. Sci. 2025, 15, 100). | Do not lead with it. Possible supporting role. |
| LAOD nonlinear analysis (GSD, Lissajous, Chebyshev) | Risse et al. 2025 review defines the agenda: expand GSD, revisit old interfaces. Heavily resourced, well-known groups. | **Avoid as a competition.** |
| Dilatational rheology model *derivation* | Mature to the point of exhaustion: LVDT, Joos, Jiang, Fainerman reorientation/compression, Danov–Kralchevsky micellar, proteins multi-conformation. New derivations are incremental. | Avoid. |
| Relaxation-spectrum inversion | **Mature in adjacent fields**, essentially absent in this one — see §2. | **This is the opening.** |

---

## 2. Recommended gap: identifiability of relaxation mechanisms in interfacial dilatational rheology

### 2.1 The observation

Interfacial dilatational rheology routinely assigns a mechanism. The asphaltene literature says
deviations from LVDT "largely originate from a distribution in adsorption coefficients";
reorientation models "excellently describe" TX45; EOR work concludes reorientation is
"negligible". These are fit-and-inspect claims.

The prior question — *given this frequency window, this noise level, and these observables,
which mechanisms are distinguishable at all, and how tightly is a relaxation time bounded?* — is
not asked. The polymer-rheology and electrochemical-impedance communities have a mature
apparatus for exactly this:

- relaxation-spectrum identification as a **Fredholm integral equation of the first kind**,
  regularised by Tikhonov, with L-curve / GCV for the parameter;
- **DRT** (distribution of relaxation times) in impedance, where deconvolution of overlapping
  processes is standard practice;
- open implementations (`pyReSpect`) with information-criterion selection of the number of modes;
- profile-likelihood structural and practical identifiability as standard in systems biology.

I found **no transfer of this apparatus to interfacial dilatational rheology.** The Risse 2025
review's stated forward agenda is GSD and revisiting old interfaces — not identifiability. That
is a gap between two adjacent literatures, which is the kind of gap that survives scrutiny
because it needs no new physics, only correct application.

### 2.2 Feasibility, tested rather than assumed

`research/discriminability/probe_mechanisms.py`, PAT window 0.005–1 Hz, 15 frequencies,
2% amplitude noise, 1° phase noise, E₀ profiled as a nuisance:

| synthetic truth | fitted model | χ² (30 observables; ≈30 = good, >53 = excluded) | verdict |
|---|---|---|---|
| barrier τ_a = 10 ms | pure diffusion | 104.5 | excluded |
| barrier τ_a = 100 ms | pure diffusion | 103.8 | excluded |
| reorientation τ_r = 100 ms | pure diffusion | 129.3 | excluded |
| reorientation τ_r = 1 s | pure diffusion | 84.9 | excluded |
| reorientation 100 ms | barrier | 1163.9 | excluded |
| barrier 100 ms | reorientation | 105.3 | excluded |

**There is real information in the observable.** Mechanism *class* is identifiable in this band
at realistic noise. But:

| quantity | 95% profile interval | width |
|---|---|---|
| barrier τ_a = 100 ms | [1.0×10⁻⁴, 6.2×10⁻¹] s | ~2 decades |

**The class is identifiable; the time constant is barely bounded.** That combination — you can
say which mechanism, but not how fast — is the honest, useful, and publishable result. It is
also the kind of result that could not have been written before someone actually computed it.

### 2.3 Why this is a *premium* paper rather than a technical note

The contribution is not "we applied Tikhonov regularisation". It is:

1. **A discriminability map**: for the standard instrument bands (PAT 0.005–1 Hz, CPT 0.1–500 Hz,
   microtensiometer 1.7 mHz–1.7 Hz), which (mechanism, timescale) pairs can be separated and
   which cannot. This immediately tells the field which published mechanistic assignments are
   supported by their own data and which are not.
2. **A reanalysis opportunity**: the assignments in the literature are attached to published
   curves. Applying the map to them is a concrete, testable, and potentially uncomfortable
   contribution — exactly what a high-impact methods paper does.
3. **An open implementation.** SurfactantKit already contains the pieces (`curve_analysis`,
   `inference` with UQ workers, `validation_tiers`, `mixture_model_selection`, `_optimize`,
   `_stats`). It is already a structure→property library with tiered validation; adding a
   rigorous identifiability layer is a natural extension, not a bolt-on.
4. **The correct statistical framing.** Profile likelihood with disconnected intervals rather
   than symmetric errors; nuisance parameters profiled, not fixed; coverage reported; boundary
   optima reported as bounds. The audit of Paper 1 taught exactly this lesson and it is a real
   differentiator, because most of this literature reports a point estimate and a symmetric
   standard error.

### 2.4 Falsification tests, to run before any further work

Each of these is under a day and would kill the candidate:

- **T1.** Search for "distribution of relaxation times" + "dilatational" / "interfacial rheology"
  in Scopus/WoS (not just web search). If a DRT-for-interfaces paper exists, the gap narrows to
  an incremental one.
- **T2.** Check whether any group has published a *negative* identifiability result for
  dilatational rheology — "these models cannot be distinguished" — including in thesis form.
  If so, the map already exists.
- **T3.** Verify that published mechanistic assignments really are fit-and-inspect, by opening
  the three most-cited ones and looking at their methods sections. If one already reports a
  profile likelihood, the novelty claim is weaker than stated here.
- **T4.** Confirm no paywalled 2024–2026 paper has done the transfer. This is the one I cannot
  do from here and is the largest single risk.

**If T1–T4 all come back clean, the candidate stands.** If T1 or T2 fails, this becomes a
technical note, not a premium paper, and should be dropped.

---

## 3. Second candidate: uncertainty propagation through HLD

**The gap.** HLD/HLD_N is used as a quantitative formulation tool, and the literature itself
documents that its inputs are shakier than its use implies. Salager's own 2020 review states
that multiple HLD definitions disagree, that the surfactant parameter is "not very accurate
unless very pure species are used", and that supplier data can be "erroneous for formulating
systems that require extreme precision". Yet HLD predictions are reported as point values.

**What I did not find:** a systematic propagation of parameter uncertainty (Cc, EACN, k, α_T)
through to formulation outcomes, with confidence intervals on the predicted optimum. Furth et al.
2024 did ML for EACN — that reduces one input's error, it does not propagate it.

**Fit to the software:** `hld.py`, `eacn.py`, `nagarajan.py`, `inference.py` already exist. This
is largely an assembly job.

**Honest weakness:** it is a "how reliable is your tool" paper, which formulators value and
high-impact journals may not. Prestige ceiling is lower than §2. Also, "nobody propagated
uncertainty" is a weaker claim than "nobody asked whether the question is answerable" — absence
of propagation is often just absence of interest.

---

## 4. Do not do these

- **Another CMC predictor.** Six independent 2025–2026 papers, state of the art already at
  R² 0.947. A seventh adds nothing.
- **Another LAOD nonlinearity paper.** The GSD programme is well resourced and has stated its
  own agenda.
- **Anything requiring new bench data.** Neither candidate above needs it. That is deliberate:
  the project has no laboratory, and every previous plan that assumed measurements died.
- **The audit's four-stage experimental programme** for the τ₂ paper. Two years, to confirm a
  negative already established.

---

## 5. Recommendation

**Go with §2, subject to T1–T4.** It is the only candidate where:
the question is prior to the field's current practice rather than competing with it;
the apparatus is mature elsewhere and absent here;
the claim is a *limitation* result, which cannot be wrong in the way the τ₂ paper was wrong —
you cannot over-claim that something is not identifiable if you have computed the profile;
and the software already exists in a form that makes the implementation honest.

The one thing that would sink it is T1/T2 finding that the map already exists. Run those first.

---

## 6. Sources consulted

Colloid/interfacial grand challenges (Front. Mater. 2025, 10.3389/fmats.2025.1563178);
AI in colloid science review (Curr. Opin. Colloid Interface Sci. 2025, 80, 101965);
Risse et al., LAOD review, Adv. Colloid Interface Sci. 2025, 345, 103625;
Guzmán et al., dilational rheology foundations and tools;
Fainerman/Rusanov/Miller reorientation and compression models;
Danov, Kralchevsky, Denkov, Ananthapadmanabhan, Lips, Adv. Colloid Interface Sci. 2006, 119, 1
and 17;
Salager et al., HLD_N for cosmetics/pharma, Cosmetics 2020, 7, 57;
Salager et al., surfactant classification confusion, 2022;
Abbott, HLD-NAC formulator's guide, 2025;
ML-assisted phase behaviour, Appl. Sci. 2025, 15, 100;
Furth et al. 2024 (EACN by ML);
Asphaltene mixture-effect dilational rheology, Langmuir 2017, 33, 1927;
Tikhonov/regularised relaxation-spectrum identification (Polymers 2023, 15, 958; Materials 2024,
17, 4870; Polymers 2024, 16, 2300; pyReSpect);
DRT in impedance spectroscopy (IOP 2025, 10.1149/1945-7111/adda7b);
profile-likelihood identifiability (Raue et al.; unidentifiable-combinations work);
CMC prediction state of the art (GCN, ChemRxiv; PharmHGT, 2026, S0927025626000674);
particle adsorption Thiele-modulus analysis, arXiv 2507.21026;
adsorption kinetics practical guidelines, ACS 2026.
