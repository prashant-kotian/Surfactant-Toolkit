# Reproducibility verification

## Current record (2026-10-28, Windows, Python 3.14)

- SB-166, SB-167 and SB-197 on one target per question (DATA_NOTES.md): corrected Ward-Tordai and
  dynamic-tension solvers; new wordings and answers in `repair_2026-10-20/answers_2026-10-28/`;
  `orig_q16_q17_q47_grade.py` -> `orig_q16_q17_q47_results.json`. Unaugmented grand totals: Claude
  132, ChatGPT 129, DeepSeek 123, Gemini 105, Qwen 106 of 200; augmented unchanged.
- Chain re-run (build_results, build_tables, build_question_set, build_bench200,
  build_figure_data): all exit 0, all figure-data cross-checks pass.
- Application P3 recomputed (`benchmark/applications/compute_p3.py`).
- Test suite: 982 passed.

## Previous record (2026-10-27, Windows, Python 3.14)

- New: `benchmark/literature_checks/` (five published results recomputed from their printed values; `run_all.py` regenerates `case1..5_results.json`; screening record `screened_papers.csv`) with five tests.
- Toolkit text only: the SDS counterion-binding source is credited to B. L. Bales alone (`thermodynamics.py`); no value changes; chain outputs (TABLE200.md, per-question CSV) byte-identical.
- Test suite: 951 passed, in full runs under both mcp 1.27.2 and mcp 2.3.0.

## Previous record (2026-10-26 night, Windows, Python 3.14)

- SB-167 reference values per wording (augmented 9.62 from the release toolkit, unaugmented 33.29
  from the stated calculation; `repair_2026-10-20/q17_regrade.py`). Augmented arms re-run against the
  current tool (`repair_2026-10-20/answers_q17/<arm>/`): Claude, ChatGPT, DeepSeek and Gemini CLI all
  answered 9.62. All verdicts unchanged; no reported cell changes.
- Gemini web GZ-50 re-asked in a new chat (`answers_fill/gemini_unaug/`, PASS -> PASS); GZ-13 and
  GZ-33 returned no answer in new chats and keep their recorded answers (DATA_NOTES.md, Net-200 fill).
- Application reference values (`benchmark/applications/compute_*.py`) and the reproducibility sets
  (`benchmark/reproducibility/`) regenerate byte-identical outputs.
- Chain re-run (build_results, build_tables, build_question_set, build_bench200,
  build_figure_data): all exit 0, all figure-data cross-checks pass; grand totals as in TABLE200.md.
- Test suite: 946 passed, in full runs under both mcp 1.27.2 and mcp 2.3.0.

## Previous record (2026-10-26, Windows, Python 3.14)

- SB-197 reference value recomputed from the release toolkit (57.982) and every arm regraded
  (`repair_2026-10-20/q47_regrade.py` -> `q47_regrade_results.json`): unaugmented arms from their
  stored answers, augmented arms from re-runs against the current tool
  (`repair_2026-10-20/answers_q47/<arm>/`). All nine verdicts unchanged (unaugmented FAIL,
  augmented PASS); no reported cell changes.
- Chain re-run (build_results, build_tables, build_question_set, build_bench200,
  build_figure_data): all exit 0, all figure-data cross-checks pass.
- Test suite: 940 passed (mcp 1.27.2).

## Previous record (2026-10-25, Windows, Python 3.14)

- Toolkit changes of 2026-10-25 (dynamic tension equation of state and surface step, emulsion
  creaming model, Aniansson-Wall inversion, Davies -SO4Na; see DATA_NOTES.md). No reported cell
  changes: the SB-197 reference-value update was staged until every augmented re-run answer existed.
- Downstream chain re-run (build_bench200.py, build_figure_data.py): bench200.json, TABLE200.md and
  the per-question CSV are byte-identical to the previous record.
- Test suite: 940 passed, in full runs under both mcp 1.27.2 and mcp 2.3.0.

## Previous record (2026-10-24, Windows, Python 3.14)

This is the record for the current tree; the dated sections below are the
history of earlier states and their test counts.

- Downstream chain re-run (`build_results.py`, `build_tables.py`,
  `build_question_set.py`, `build_bench200.py`, `build_figure_data.py`), all
  exit 0, all figure-data cross-checks pass. The raw graders were not re-run
  in this pass: their JSONs are pre-correction artifacts (`REPRODUCE.md` §1).
- 27 result cells, 1 declared data gap; `per_question_results_200.csv` 1800
  rows; grand pass rates as in `TABLE200.md` (including the 2026-10-23
  GZ-17 gold, see `DATA_NOTES.md`).
- `source_sha` now hashes canonical JSON content, so a formatting-only
  rewrite of a grading JSON no longer changes the ledger.
- `CLASSIFICATION200.md` column headers now follow the model order of the
  data (one column per model); F7 coverage is derived from the gradable count
  and F9's date column is named `date_in_source_path`.
- Test suite: 935 passed, under both mcp 1.27.2 and mcp 2.3.0.

## Earlier record (2026-10-19)

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

## Re-run after the Claude groundzero correction (2026-10-19, Windows, Python 3.14)

Full §1 chain re-run (graders included) through `build_figure_data.py`:

- Only the groundzero unaugmented Claude (5/10 -> 27/43) and ChatGPT
  (22/43 -> 27/43) cells changed; every other cell identical (see DATA_NOTES.md).
- `results/per_question_results.csv`: 1379 rows (+33 = Claude's 10 -> 43
  gradable); grand CSV rows = sum of grand denominators (1824).
- Re-graded session JSONs differing only in formatting (zero parsed
  differences) were left at their committed bytes.
- Ledger `src` paths now always POSIX (`build_results.py` previously emitted
  `\` separators on Windows).
- `build_figure_data.py`: F8 failure anatomy no longer cross-attributes
  per-field details between models sharing one multi-model grading JSON; F8
  failures reconcile to F1 (456 = 456). All figure-data cross-checks pass.
- Test suite: 828 passed (823 + 5 F8 regression tests).

## Re-run after the exclusion and re-grade audit (2026-10-19, Windows, Python 3.14)

- X-19 now excluded from every arm (the three post-fix augmented extrahard
  cells drop from n=50 to 49); grand augmented: Claude 180/189, ChatGPT
  174/189, DeepSeek 173/186.
- DeepSeek augmented groundzero/hard50 re-graded against the corrected golds
  (GZ-48 FAIL->PASS, H-04 FAIL->PASS, H-26 PASS->FAIL); `gold_state` verified
  from stored gold values. All reported-cell graders now in REPRODUCE.md §1;
  Antigravity and Claude-hard50 cells re-graded from the repository with zero
  verdict differences.
- Grand CSV rows = sum of grand denominators (1821).
- Test suite: 832 passed (828 + 4 exclusion/gold-state regression tests).

## Re-run after the counterion-factor correction and chart/typo repair (2026-10-21, Windows, Python 3.13)

- Toolkit: ionic factor corrected to (1 + beta) = (2 - alpha) in thermodynamics, orchestrate and the MCP
  tool descriptions; Guo-Rong-Ying EO weight corrected (0.33). Ported from the development copy:
  `foam.py`, `dynamic_surface_tension.py`, Gibbs-Marangoni elasticity in `adsorption.py`.
- Re-runs and regrades are applied by `benchmark/repair_2026-10-20/corrections.py` and logged per cell
  (`results.json` `corrections`); grading of the re-runs: `benchmark/repair_2026-10-20/grade_repair.py`
  (repair set) and `grade_repair.py --dg` (counterion set). Answers are in `answers/` and `answers_dg/`
  with per-answer meta (model, tools/MCP state, chart attachment, question variant).
- Rebuild order: `grade_repair.py`, `grade_repair.py --dg`, `build_results.py`, `build_bench200.py`,
  `build_tables.py`, `build_figure_data.py` (all cross-checks pass; F8 uses the corrected field detail).
- Claude augmented counterion re-runs completed 2026-10-21 (connector enabled on the run account); no re-runs pending.

## Accepted non-determinism (documented, not chased)

- `source_sha` in the ledger hashes canonical grading-JSON content (since
  2026-10-24; before that, the bytes as last written, which a re-grade could
  reformat without revaluing).
- `gold_state` is verified from the gold values stored in each cell's grading
  records (not file dates); count-only sources without per-field gold read
  `not field-verifiable`.
- Exclusion policy and ungradeable-record rules live in `build_results.py`
  and `question_set.json` as data; changing them changes denominators by
  design, not by drift.

## Verdict

**PASS.** Reported numbers = this tree's output. The fourth source set's
orig_hard50 grading JSONs are frozen inputs (first-release harness); their
generation is outside this repository by design.
