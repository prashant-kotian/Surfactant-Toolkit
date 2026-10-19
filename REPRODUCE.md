# Reproduction guide

Environment: Python 3.11+ on Linux/macOS/Windows; `pip install rdkit "mcp[cli]"`.
All commands run from the repository root.

## 1. Results chain (no network, no model access required)

```bash
export PYTHONPATH=$PWD/src
python3 benchmark/groundzero_50/grade_unaugmented.py
python3 benchmark/groundzero_50/grade_chat_runs.py
python3 benchmark/hard50/grade_chat_runs.py
python3 benchmark/extrahard50/grade_unaugmented.py
python3 benchmark/extrahard50/grade_augmented_chat_runs.py
python3 benchmark/build_results.py
python3 benchmark/build_tables.py
```

The graders read only committed transcripts (`benchmark/**/`,
`benchmark/transcripts/`) and bank golds, and rewrite the committed grading
JSONs; `build_results.py` then assembles `results/results.json` and
`build_tables.py` renders the manuscript tables. Expected at this release:
30 result cells, 1 declared data gap, and `results/per_question_results.csv`
byte-identical to the committed file. Re-runs may differ from the committed
bytes only in JSON formatting (the provenance ledger's `source_sha` column
tracks this); verdict data must not change.

## 2. Test suite

```bash
PYTHONPATH=src python3 -m pytest -q        # 823 passed
```

## 3. Executing new model runs (optional)

Augmented condition: start the toolkit MCP server (`src/surfactantkit/mcp_server`)
behind any MCP-capable client; the `benchmark/mcp_harness/` scripts show the
exact server configuration used for the reported runs. Unaugmented condition:
pose the text of `benchmark/*/question_bank*.txt` unmodified, saving one
transcript file per question in the run's transcript directory (any layout —
new graders read the per-question `.txt` files with a `FINAL_JSON:` block as
the last element, which all runs follow).

Note: `results/results.json` records the toolkit commit used at grade time in
`toolkit_state`; for this self-contained repository that is this repository's
own HEAD (override with `$TOOLKIT_DIR`).
