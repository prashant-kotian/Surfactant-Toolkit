# Reproducibility verification (2026-10-19)

The results chain was re-run end to end in this repository after the
single-benchmark restructuring, on Linux with Python 3.13
(`pip install rdkit "mcp[cli]" pytest`).

## Procedure

All commands from `REPRODUCE.md` §1, in order (graders included), followed by
`build_question_set.py` and `build_bench200.py`, then the full test suite.

## Checks performed

- The five grader scripts ran from the tree and rewrote the committed
  grading JSONs; parsed contents identical.
- `build_results.py`: 30 result cells, 1 declared data gap; the provenance
  ledger and per-question CSV regenerated from the rewritten JSONs.
- `results/per_question_results.csv`: 1346 verdict rows reproduced; all
  verdicts and denominators unchanged (cell identifiers use the source-set
  names `groundzero`/`hard50`/`extrahard` after the 2026-10-19 rename).
- `build_bench200.py`: `results/bench200.json` byte-stable across re-runs;
  each grand cell equals the sum of its per-source-set cells (e.g. Claude
  augmented = 36 + 45 + 50 + 50 = 181 of 190 gradable); CSV row count equals
  the sum of grand denominators (1791).
- `build_question_set.py`: 200 items, deterministic (rule table + curated
  overrides; no network, no model access).
- Test suite: 823 passed.

## Accepted non-determinism (documented, not chased)

- `source_sha` in the ledger hashes the grading-JSON bytes as last written;
  a fresh re-grade may reformat (not revalue) them.
- `gold_state` flags derive from source-file mtimes; a re-run after the
  2026-10-18 gold corrections reads `post-fix` by design.
- Exclusion policy and ungradeable-record rules live in `build_results.py`
  and `question_set.json` as data; changing them changes denominators by
  design, not by drift.

## Verdict

**PASS.** Reported numbers = this tree's output. The fourth source set's
orig_hard50 grading JSONs are frozen inputs (first-release harness); their
generation is outside this repository by design.
