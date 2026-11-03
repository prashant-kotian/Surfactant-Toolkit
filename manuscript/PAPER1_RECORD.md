# Paper 1 — Editorial Record (NOT part of the manuscript)

This file holds the material that must never be typeset into the paper: provenance of every
number, the record of corrections, the statistical multiplicity note, and the commit audit. It is
kept in the repository so the paper remains auditable, and is deliberately a separate file so a
build step cannot pick it up by accident.

---

# Appendix (editorial — not part of the manuscript)

This appendix replaces the per-section "items needing your decision" blocks that accompanied the
drafts. Nothing here belongs in the typeset paper; it is the resolved record so that the typeset
pass can be checked against it.

## A.1 Number provenance

Every quantitative claim in §3–§7 traces to a committed artifact:

| claim | source |
|---|---|
| window 1.2–3.0× CMC, optimum 1.5× | `e0_sensitivity.json` |
| E₀ invariance over 25–103.5 mN/m | `e0_sensitivity.json` |
| δ vs frequency, Cₘ on/off (§4.2) | computed in FINDINGS §15.3 |
| 45° ceiling, 0.48° LVDT agreement | FINDINGS §15.1–15.2 |
| I₂/I₁ 0.007–0.019 per doubling | computed from K = 140 sweep |
| ≥ 5 frames/cycle | FINDINGS §11 |
| ε² law, 9.21× amplification | FINDINGS §16.1 |
| 0.23% relSE, MC scatter 0.033° | `recovery_uncertainty.json` |
| corr = 0.997, 2-parameter relSE 3.0% | `recovery_uncertainty.json` |
| Protocol B non-monotonicity | `protocol_b_anchored.json` |
| ∂δ/∂ln Cₘ and ∂δ/∂ln τ₂ | FINDINGS §17.1 |
| n̄ propagation table | FINDINGS §17.2 |
| CTAB/SDS bands, 16 Hz crossing | `ionic_generality.json` + manual extension |
| Bond number refutation | FINDINGS §17.3 |

## A.2 Corrections made to supplied material

1. **Loss-angle sign.** The supplied narrative had the reservoir *depressing* δ (32.4° → 18.8°) and
   called that depression the signature. Computation gives the opposite: the reservoir *raises* δ to
   68–80°, above the 45° diffusion ceiling. The structure was kept; the sign was reversed. Endorsed.
2. **LVDT "convention conflict."** Not a conflict. With ν_D computed self-consistently the agreement
   is 0.48°. The earlier 18.8°-vs-32° gap was our own parameter borrowing.
3. **45° vs 35°.** The ceiling is 45°, not arctan(1/√2) ≈ 35.3°.
4. **"I₂/I₁ is blind"** → "unresolvable": it moves 0.007–0.019 against a 0.025–0.05 noise floor.
5. **"Two concentrations suffice"** — retracted; the correlation is not monotone in the number of
   concentrations.
6. **"Large amplitude trades linearity against τ₂ sensitivity"** — false; δ is constant to 0.17°
   over ε = 0.03–0.30.
7. **Bond number** — removed; Bo = 0.184 at the stated worst corner (OD 1.5 mm, γ = 30 mN/m), and
   Bo does not govern oscillatory detachment.
8. **Frumkin "numerical stiffness at a ≥ 4"** — false; it was a wrong-sign residual in the
   isotherm solver. Fixed.

## A.3 Multiplicity note (statistical)

The manuscript reports several derived figures from a family of model sweeps. Two deserve explicit
labelling because they were selected *after* inspection rather than pre-registered:

- **ωτ₂ ≈ 3 sensitivity maximum (§4.2, §6.1).** Located by scanning, not derived. No scaling
  collapse exists; at matched ωτ₂, CTAB and C₁₂E₈ differ in amplitude by 2.4×. Presented as an
  observation.
- **The 13.3° per e-fold lever arm at ω = 6.28 (§5.6, §7.2).** This is the largest of four
  frequencies scanned (5.99, 8.89, 11.16, 13.26 at ω = 3/4/5/6.28). Operating at the top of the
  instrument band because the lever arm is largest there is legitimate, but the reported precision
  is the best of the four, not a typical value. Report as "0.23% at ω = 6.28 rad/s" rather than
  "0.23%" without the qualifier.

## A.4 Provenance audit (final)

Verified via the GitHub API on 2025-11-03 against the user's stated commits:

| commit | repo | status |
|---|---|---|
| `f73fa80` | surfactant-dynamics-sim | **resolves** |
| `e9da2d0` | surfactant-dynamics-sim | **resolves** |
| `caa1f57` | Surfactant-Toolkit (public) | **does not exist** |
| `85e1218` | Surfactant-Toolkit (public) | **does not exist** |
| `27691b1` | Surfactant-Toolkit `main` | **is** the frozen public tip, as stated |

Both `85e1218` and `caa1f57` were described as internal-mirror commits. Neither is a public artifact,
so neither should be cited in the manuscript or in any preprint announcement. The manuscript text
lives on `arena/fa2f9ecf-surfactant-toolkit` at `a92c123`, `a88bd29`, `34aedd3`, `dd48f4c`.

**I cannot verify the current `.docx`/`.pdf`.** I have no access to the internal SurfactantKit
mirror. The synthesis documents I *can* see (`f73fa80`) were consistent with the corrected §4/§7 and
carried the two errors since fixed (Bond number; the "2.4 to 73 Hz" head frequency). If `e9da2d0`
applied those fixes, good; if `caa1f57` is where the fixes actually live, then the public-facing
record still needs them.

## A.5 What is not verified

- **Refs 4, 5, 7** (Stadler 2016, Risse 2025, Kotsmar 2012) are cited from supplied descriptions.
  We have not read them. They are cited only for the specific property named in the reference list.
- **The 25–50 fps acquisition rate** is a manufacturer sensor specification, not a figure from
  Miller et al. 2004. The text attributes it correctly; the reference is cited for instrumentation
  principles only.
- **No formatting spec exists in this repository.** `MANUSCRIPT_FORMATTING_SPEC.md` is not present
  on any branch of `Surfactant-Toolkit`, so this draft is plain Markdown with no journal-specific
  typesetting applied. The typesetting pass must supply the spec.

---

## A.6 Typeset deliverable audit (2025-11-03)

I downloaded and inspected the typeset artifact at `surfactant-dynamics-sim` commit `0cb880f`
(`docs/Direct_Inversion_Micellar_Relaxation_Kinetics_LAOD.docx` and `.pdf`).

**Verified correct:**
- 60 pages, confirmed by opening the PDF (11,318 words).
- Body text is faithful to the manuscript; no number drifted in typesetting.
- Times New Roman throughout; running page numbers present; title block and corresponding-author
  line present; 7 sections plus references.
- TimesNewRomanPS-BoldItal is used (1,904 chars), so the italic restriction did not strip emphasis
  entirely.

**Defects, in order of severity:**

1. **The editorial appendix was typeset into the paper** (pages 56–60). `build_full_manuscript_docx.py`
   line 332 explicitly handles `# Appendix` and emits it under a sanitised heading, "Appendix:
   Computational Provenance, Sensitivity Audits, and Multiplicity Notes". As printed, the submission
   contains: an itemised list of the corrections made to supplied material (including "Endorsed"),
   the commit-SHA audit naming commits that do not resolve, and internal repository paths. None of
   this is publishable, and item A.1 also discloses internal tooling. **This is why the manuscript is
   now split into `PAPER1_SUBMISSION.md` and this file.**

2. **Markdown not rendered in the appendix**: 21 lines retained backticks and `**`, in both the
   .docx and the PDF.

3. **Italics are not applied to any body text.** Measured across all 60 pages: 0 italic characters
   from non-bold-italic Times New Roman. The spec says italics are "restricted exclusively to
   subsection headings", so 0.00% body italics is correct — but the consequence is that every
   `*term*`, species name and emphasis marker in the source was silently dropped rather than
   rendered. Subsection headings are present in the text stream, so the restriction was applied; it
   is worth confirming that headings are actually italic in the .docx run properties rather than
   merely unstyled.

4. **"60 pages" is partly a mechanical artefact.** With 2.0 double spacing on 10.4k words, 60 pages
   is expected. Page 60 ends mid-appendix on the sentence about the formatting spec. Once the
   appendix is removed the paper will be substantially shorter, and the page count should be
   re-derived rather than carried into any cover note.

**Not verifiable by me:** commits `caa1f57` and `45d7035` return "No commit found for SHA" from the
public API. They were described as living in an internal mirror; if so they cannot be cited publicly.
`f73fa80`, `e9da2d0` and `0cb880f` all resolve.


---

## A.7 The "Reference flags" defect (2025-11-03)

The 37-page typeset artifact at `f858592` corrected the appendix leak and the paragraph
reconstitution (795 fragments -> 260 paragraphs; verified: 37 pages, no editorial appendix).

**But the paper still ended on an editorial section, and this one was my fault.** `PAPER1_SUBMISSION.md`
carried a subsection headed "Reference flags" as the final section of the manuscript. As printed on
page 37 it told the reader that:

- a statement "supplied" to us was "not true", and that we had rebutted it;
- the manuscript "removed the Bond-number argument from the text";
- refs "4, 5, 7 are cited from the supplied descriptions; we have not verified their contents";
- a reference was cited "for the instrumentation principles only".

That is an internal audit printed as the last page of a journal submission. It was removed from the
submission copy, together with the annotations attached to references 1, 6 and 7 and the sentence
"Items 1-5 are as supplied by our collaborators". The substance is preserved here:

- **Bond number.** Bo = ΔρgR²/γ with R = tip radius gives 0.184 at OD 1.5 mm, γ = 30 mN/m; Bo = 0.15
  corresponds to OD ≈ 1.35 mm. The claim holds for OD ≤ 1.2 mm at γ = 30 mN/m, or all listed
  diameters only if γ ≳ 37 mN/m. Independently, Bo is static and does not govern oscillatory
  detachability. The paper therefore cites refs 6 and 7 empirically and contains no Bond argument.
- **Acquisition rate.** 25-50 fps is a manufacturer sensor specification, not a figure from Miller et
  al. 2004. The paper cites ref 8 for instrumentation principles only.
- **Refs 4, 5, 7** are cited from descriptions supplied to us and were not read directly. They are
  cited only for the property named against them, and this limitation should be checked before
  submission by someone with access to the papers.

**Root cause, so it does not recur:** I created the split by truncating `PAPER1_FULL_DRAFT.md` at the
`# Appendix` marker. Everything *after* that marker in the editorial sense was removed, but the
"Reference flags" subsection sat *before* it, at the end of the references. The split boundary was
in the wrong place. Any future split must be verified by reading the last page of the rendered PDF,
not by trusting the marker.

**Verification now required before submission:** re-render and confirm the paper terminates at
reference 8, with no trailing section.
