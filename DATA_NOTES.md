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
- extrahard unaugmented cells and the Antigravity augmented cell of 2026-10-09
  score n=49: X-19's question text carried an authoring typo until 2026-10-16;
  cells executed after the fix include it at n=50 (Claude and ChatGPT
  augmented extrahard are full-set).
- orig_hard50 scores n=50; NOT_ANSWERED records (Claude unaugmented 1, ChatGPT
  unaugmented 4) reflect no answer content captured and are excluded from n
  as ungradable rather than scored as failures.
- Claude's groundzero unaugmented row is a 10-question partial (the only
  recovered Sept-2026 transcript); it is reported as `partial coverage` in the
  provenance ledger and excluded from main comparisons.
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
`remediation_note`; all cells re-graded programmatically. The GZ-48 correction
revealed a masked model error (Antigravity unaugmented PASS->FAIL): reported
rather than smoothed. See `results/REPRODUCIBILITY.md`.
