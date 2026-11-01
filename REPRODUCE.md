# Reproduction guide

Environment: Python 3.11+ on Linux/macOS/Windows;
`pip install -e ".[mcp,calc,dev]"` (tested with mcp 1.27.2 and 2.3.0).
All commands run from the repository root.

## 1. Results chain (no network, no model access required)

```bash
export PYTHONPATH=$PWD/src
python3 benchmark/source_sets/groundzero/grade_unaugmented.py
python3 benchmark/source_sets/groundzero/grade_chat_runs.py
(cd benchmark/source_sets/groundzero && python3 grade_antigravity.py && python3 grade_augmented_deepseek.py)
python3 benchmark/source_sets/hard50/grade_chat_runs.py
(cd benchmark/source_sets/hard50 && python3 grade_antigravity.py && python3 grade_augmented_deepseek.py && python3 grade_claude50.py)
python3 benchmark/source_sets/extrahard/grade_unaugmented.py
python3 benchmark/source_sets/extrahard/grade_augmented_chat_runs.py
python3 benchmark/build_results.py
python3 benchmark/build_tables.py
python3 benchmark/build_question_set.py
python3 benchmark/build_bench200.py
python3 benchmark/build_figure_data.py   # results/figs_data/*.tsv for manuscript figures
```

Windows (PowerShell):

```powershell
$env:PYTHONPATH = "$PWD\src"
python benchmark\source_sets\groundzero\grade_unaugmented.py
python benchmark\source_sets\groundzero\grade_chat_runs.py
Push-Location benchmark\source_sets\groundzero; python grade_antigravity.py; python grade_augmented_deepseek.py; Pop-Location
python benchmark\source_sets\hard50\grade_chat_runs.py
Push-Location benchmark\source_sets\hard50; python grade_antigravity.py; python grade_augmented_deepseek.py; python grade_claude50.py; Pop-Location
python benchmark\source_sets\extrahard\grade_unaugmented.py
python benchmark\source_sets\extrahard\grade_augmented_chat_runs.py
python benchmark\build_results.py
python benchmark\build_tables.py
python benchmark\build_question_set.py
python benchmark\build_bench200.py
python benchmark\build_figure_data.py
```

The graders read only committed transcripts (`benchmark/source_sets/*/`,
`benchmark/transcripts/`) and source-set golds, and rewrite the committed
grading JSONs; the orig_hard50 grading JSONs are frozen inputs from the
first-release harness (their transcripts and the first-release grading scripts are in
`benchmark/transcripts/orig_hard50_first_release/`; `repair_2026-10-20/orig_hard50_grade.py`
reproduces every frozen verdict).
`build_results.py` assembles `results/results.json`, `build_tables.py` renders
the per-source-set diagnostic tables, `build_question_set.py` rebuilds the
benchmark classification, and `build_bench200.py` produces the headline
tables. Expected at this release: 27 result cells in `results.json` (three
source sets x nine model arms), 1 declared data gap, and
`results/per_question_results.csv` byte-identical to the committed file (1350
records). The fourth source set, `benchmark/source_sets/orig_hard50/`, has no
in-tree grader; its committed grading JSONs are consumed directly by
`build_bench200.py`.

The provenance ledger's `source_sha` hashes canonical JSON content (sorted
keys, compact separators), so rewriting the same data with different
formatting does not change it. Verdict data downstream of the corrections
layer must not change, and `bench200.json`, `TABLE200.md`,
`CLASSIFICATION200.md` and `results/figs_data/*.tsv` reproduce
byte-identically. Note on the raw grading JSONs: they are deliberately
preserved as PRE-correction audit artifacts. The machine-readable corrections
applied downstream (`benchmark/repair_2026-10-20/`, logged per cell in
`results.json`) change some verdicts when the raw graders are re-run against
corrected golds, so re-running the graders rewrites those files; the diff is
exactly the logged corrections (see `DATA_NOTES.md`) and must not be
committed. The grand-table cell totals must equal the sums of their
per-source-set cells; `build_bench200.py` derives them from the same JSONs.

### Results outside the benchmark

```bash
python3 benchmark/applications/compute_p2_p6.py      # application reference values (also compute_p3.py, compute_p4.py)
python3 benchmark/reproducibility/extract_runs.py     # R-01..R-18 answers -> runs.json
python3 benchmark/reproducibility/score_runs.py       # -> scored.json
python3 benchmark/reproducibility/sb_extension/score_ext.py   # eleven SurfBench-200 questions x5 -> scored_ext.json
python3 benchmark/literature_checks/run_all.py        # five published results recomputed -> case*_results.json
```

Each regenerates its committed output; see the README in each folder.

## 2. Test suite

```bash
python3 -m pytest -q        # 982 passed (tests the working tree under src/)
python3 -X utf8 benchmark/llm_checks/test_score.py   # scorer self-test: SELF-TEST: PASS
```

## 3. Executing new model runs (optional)

Augmented condition: start the toolkit MCP server (`src/surfactantkit/mcp_server`)
behind any MCP-capable client; the `benchmark/mcp_harness/` scripts show the
exact server configuration used for the reported runs. Pose the question text
unmodified: `benchmark/source_sets/*/question_bank*.txt` (for `orig_hard50`,
`question_bank_unaugmented.txt` in the unaugmented condition and the
single-question `Q##.md` files, which name the tool to call, in the augmented
condition), saving one
transcript file per question in the run's transcript directory (any layout —
new graders read the per-question `.txt` files with a `FINAL_JSON:` block as
the last element, which all runs follow).

Note: `results/results.json` records the toolkit commit used at grade time in
`toolkit_state`; for this self-contained repository that is this repository's
own HEAD (override with `$TOOLKIT_DIR`).
