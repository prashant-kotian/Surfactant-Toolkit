# Data notes

## Benchmark composition

SurfBench-200 is one 200-question benchmark (ids `SB-001..SB-200`). It was
authored as four 50-question working groups — groundzero (definitional /
standard practice), hard50 (multi-step computation), extrahard (adversarial),
and orig_hard50 (toolkit-execution items carried in from the project's first
benchmark release) — and merged by `benchmark/build_question_set.py`. The
group names survive only as `source_sets/` directories and `authoring_alias`
columns for provenance; all reporting uses the merged ids and the question
classification (`skill type` and `subconcept` per question) it carries.
Question texts are frozen; versions are tracked through the gold files, not
silent edits. The orig_hard50 grading JSONs are committed frozen inputs from
the first-release harness; their raw transcripts are not re-committed here.

## Scoring

Field-level comparison against `gold_answers.json` (numeric tolerance by
magnitude class; categorical labels by token-set equality after stripping
justifications; NO-JSON detection per transcript). The rule set lives in the
per-source-set graders and is exercised by their self-tests.

## Exclusions (all symmetric, declared in `benchmark/build_results.py` and
## `benchmark/question_set.json`)

- groundzero n=43: six questions requiring chart/image input (GZ-15, GZ-36..40;
  chart files retained in the source directory for reference) and one retired
  ambiguous question (GZ-26).
- hard50 n=47: three figure-dependent questions (H-36..38).
- extrahard n=49 in every arm: X-19's question text carried an authoring typo
  until 2026-10-16. It is excluded from all arms, including the augmented
  cells run after the fix (which had included it at n=50 until 2026-10-19),
  so every unaugmented/augmented pair is scored on the same questions.
- orig_hard50 scores n=50; NOT_ANSWERED records (Claude unaugmented 1, ChatGPT
  unaugmented 4) reflect no answer content captured and are excluded from n
  as ungradable rather than scored as failures.
- 2026-09-21 Gemini unaugmented groundzero: 11 ungradable (no-content) records
  removed from n; an earlier alternative grading pass (22/50) is not
  reproducible from the transcripts in this repository and is therefore not
  reported.

## Classification

Labels are rule-derived from gold field names, declared failure axes, and
question text (see `build_question_set.py`), with 28 manual corrections in
`question_set_overrides.py` after full review of the generated assignment.
Small categories (n ≤ 3: rheology, solubilization, data_fitting,
analysis_synthesis, foaming_films, micelle_thermodynamics, surface_wetting)
are reported but should be read as descriptive only.

## Corrections

H-04 (reference-temperature mismatch), H-26 (surfactant identity), GZ-48
(factor-of-2 row mix-up): golds corrected on 2026-10-18 with per-field
`remediation_note`. The GZ-48 correction revealed a masked model error
(Antigravity unaugmented PASS->FAIL): reported rather than smoothed. See
`results/REPRODUCIBILITY.md`.

Re-grade completeness (2026-10-19): the DeepSeek augmented groundzero and
hard50 cells had NOT been re-graded after the gold corrections (their graders
sat outside the documented chain; the hard50 one also had a stale import and
answer path). Re-graded: GZ-48 FAIL->PASS, H-04 FAIL->PASS, H-26 PASS->FAIL
(DeepSeek's 42-carbon answer matched the superseded gold). `gold_state` is now
verified from the gold values stored in each cell's records rather than
inferred from file dates; every per-field cell verifies as post-fix, and the
count-only hard50 unaugmented cells are labelled "not field-verifiable". All
graders that produce a reported cell now run in REPRODUCE.md section 1, and
Claude's hard50 unaugmented transcript is committed in-repo
(`benchmark/source_sets/hard50/claude_unaugmented/`).

Claude's groundzero-unaugmented row (2026-10-19): the 10-question partial used
through 2026-10-18 was a different, genuinely truncated Sept-2026 capture; a
complete GZ-01..50 transcript was subsequently located and verified (all 50
answer blocks structurally intact, independently-derived rather than copied
from any gold version). Grading it surfaced that `grade_unaugmented.py` had no
mechanism to resolve long-prose judgment fields -- `results/results.json`'s
verdict conversion silently scored every unresolved field as a FAIL, which
also understated ChatGPT's and Gemini's groundzero-unaugmented rows (Gemini's
score did not move; no judgment field happened to swing its total). A
`JUDGMENT_RESOLUTIONS` table (31 entries, one rationale per field, same
discipline as the Antigravity graders) now resolves every pending field.
