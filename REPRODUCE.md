# Reproduce the manuscript numbers

```bash
pip install rdkit "mcp[cli]"            # only external deps (see pyproject extras)
python3 benchmark/paper3_groundzero/grade_unaugmented_gz50.py
python3 benchmark/paper3_groundzero/grade_session_new_gz50.py
python3 benchmark/paper3_hard50/grade_session_new_h50.py
PYTHONPATH=src python3 benchmark/paper3_extra50/grade_unaugmented_5models.py
PYTHONPATH=src python3 benchmark/paper3_extra50/grade_augmented_claude_deepseek.py
python3 benchmark/results_master.py      # -> benchmark/RESULTS_MASTER.{json,md}
python3 benchmark/make_tables.py         # -> paper3/RESULTS_*
```
Expected at the freeze tag: `cells=30 gaps=1`; `RESULTS_PER_Q.csv`
byte-identical to the committed copy. The per-question graders are
deterministic given the committed transcripts and golds; re-running the
chain may rewrite JSON serialization (different `source_sha` values in the
provenance ledger) without changing any verdict — verdict-level equality is
the reproducibility criterion, and `toolkit_state` reports this repo's own
HEAD (self-contained layout; `$TOOLKIT_DIR` overrides).

Augmented-condition execution (optional, needs API access):
`benchmark/harness/` + `src/surfactantkit/mcp_server` (MCP config per
harness README). Test suite: `PYTHONPATH=src python3 -m pytest -q` → 823 passed.
