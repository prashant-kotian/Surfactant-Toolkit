# Data notes

## Banks
Three 50-question banks (GroundZero, Hard50, ExtraHard), each with two
conditions (unaugmented / augmented). Question texts are frozen; banks are
versioned through the gold files, not silent edits.

## Scoring
Field-level comparison against `gold_answers.json` (numeric tolerance by
magnitude class; categorical labels by token-set equality after stripping
justifications; NO-JSON detection per transcript). The rule set lives in the
per-bank graders and is exercised by their self-tests.

## Exclusions (all symmetric, declared in `benchmark/build_results.py`)
- GroundZero n=43: six questions requiring chart/image input (GZ-15, GZ-36..40;
  chart files retained in the bank directory for reference) and one retired
  ambiguous question (GZ-26).
- Hard50 n=47: three figure-dependent questions (H-36..38).
- ExtraHard unaugmented cells and the Antigravity augmented cell of 2026-10-09
  score n=49: X-19's question text carried an authoring typo until 2026-10-16;
  cells executed after the fix include it at n=50 (Claude and ChatGPT
  augmented ExtraHard are full-bank).
- Claude's GroundZero unaugmented row is a 10-question partial (the only
  recovered Sept-2026 transcript); it is reported as `partial coverage` in the
  provenance ledger and excluded from main comparisons.
- 2026-09-21 Gemini unaugmented GroundZero: 11 ungradable (no-content) records
  removed from n; an earlier alternative grading pass (22/50) is not
  reproducible from the transcripts in this repository and is therefore not
  reported.

## Corrections
H-04 (reference-temperature mismatch), H-26 (surfactant identity), GZ-48
(factor-of-2 row mix-up): golds corrected on 2026-10-18 with per-field
`remediation_note`; all cells re-graded programmatically. The GZ-48 correction
revealed a masked model error (Antigravity unaugmented PASS->FAIL): reported
rather than smoothed. See `results/REPRODUCIBILITY.md`.
