# Reproduction guide

Environment: Python 3.11+ on Linux/macOS/Windows; `pip install rdkit "mcp[cli]"`.
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

The graders read only committed transcripts (`benchmark/source_sets/*/`,
`benchmark/transcripts/`) and source-set golds, and rewrite the committed
grading JSONs; the orig_hard50 grading JSONs are frozen inputs from the
first-release harness (their grader is not part of this repository).
`build_results.py` assembles `results/results.json`, `build_tables.py` renders
the per-source-set diagnostic tables, `build_question_set.py` rebuilds the
benchmark classification, and `build_bench200.py` produces the headline
tables. Expected at this release: 30 result cells from the three in-tree
graders (plus 9 orig_hard50 cells read as frozen inputs), 1 declared data gap,
and `results/per_question_results.csv` byte-identical to the committed file.
Re-runs may differ from the committed bytes only in JSON formatting (the
provenance ledger's `source_sha` column tracks this); verdict data must not
change. The grand-table cell totals must equal the sums of their
per-source-set cells; `build_bench200.py` derives them from the same JSONs.

## 2. Test suite

```bash
PYTHONPATH=src python3 -m pytest -q        # 832 passed
```

## 3. Executing new model runs (optional)

Augmented condition: start the toolkit MCP server (`src/surfactantkit/mcp_server`)
behind any MCP-capable client; the `benchmark/mcp_harness/` scripts show the
exact server configuration used for the reported runs. Unaugmented condition:
pose the text of `benchmark/source_sets/*/question_bank*.txt` (or, for
`orig_hard50`, the single-question `Q##.md` files) unmodified, saving one
transcript file per question in the run's transcript directory (any layout —
new graders read the per-question `.txt` files with a `FINAL_JSON:` block as
the last element, which all runs follow).

Note: `results/results.json` records the toolkit commit used at grade time in
`toolkit_state`; for this self-contained repository that is this repository's
own HEAD (override with `$TOOLKIT_DIR`).
