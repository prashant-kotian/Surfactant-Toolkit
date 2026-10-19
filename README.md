# Surfactant Toolkit

A computation toolkit for surfactant and interfacial chemistry, with the
benchmark used to evaluate large language models with and without tool access
(Paper 3). This repository is the complete submission artifact: the library,
its test suite, three question banks with recomputed reference answers (golds),
raw model-answer transcripts, the graders, and the results pipeline. Every
number reported in the manuscript is regenerated from this tree; none is
typed by hand.

## Layout

```
src/surfactantkit/      toolkit library (CMC, mixed micellization, thermodynamics,
                        model selection, adsorption, HLD, QSPR connectors) + MCP server
tests/                  test suite (823 tests)
benchmark/
  groundzero_50/        Question bank 1: definitional/standard-practice items (GZ-01..50)
  hard50/               Question bank 2: multi-step computation items (H-01..50)
  extrahard50/          Question bank 3: adversarial items (X-01..50)
  mcp_harness/          Model Context Protocol harness used for the augmented condition
  transcripts/          raw model answers, one file per run/question
  build_results.py      collects all graded cells into results/results.json (+ summary)
  build_tables.py       renders manuscript tables and SI artifacts from that JSON
results/                generated: results.json, tables.md/.tex, per_question_results.csv,
                        REPRODUCIBILITY.md
```

Each bank directory contains its `question_bank*.txt` inputs, `gold_answers.json`
(reference values, with `remediation_note` fields where golds were corrected),
the grader scripts, and per-run grading JSONs.

## Install and run

```bash
pip install rdkit "mcp[cli]"        # only external deps (see pyproject extras)
PYTHONPATH=src python3 -m pytest -q  # 823 passed
```

Reproduce the reported numbers (deterministic given the committed transcripts):

```bash
python3 benchmark/groundzero_50/grade_unaugmented.py
python3 benchmark/groundzero_50/grade_chat_runs.py
python3 benchmark/hard50/grade_chat_runs.py
PYTHONPATH=src python3 benchmark/extrahard50/grade_unaugmented.py
PYTHONPATH=src python3 benchmark/extrahard50/grade_augmented_chat_runs.py
python3 benchmark/build_results.py     # -> results/results.json, results_summary.md
python3 benchmark/build_tables.py      # -> results/tables.md, tables.tex, per_question_results.csv
```

`per_question_results.csv` must match the committed copy byte for byte; see
`results/REPRODUCIBILITY.md` for the recorded verification.

## Data and provenance policy

- Reference values are recomputed by the toolkit (`src/surfactantkit/`); gold
  corrections made during audit are recorded as machine-readable
  `remediation_note` fields on the affected entries, with cause.
- Questions excluded from scoring (undeliverable multimodal inputs, retired
  ambiguous wording, and items affected by a pre-fix question packet) are
  declared in `build_results.py` as data, applied identically to every model
  in the affected bank, with the denominator reported (n after exclusions).
- Answers lost to capture defects (no recoverable record) are treated as
  ungradable and excluded from denominators, not scored against the model;
  cells covering fewer than the full bank are flagged `partial coverage`.
- Every grading input carries a content hash in the provenance ledger
  (`results/results_summary.md`).

License: MIT. Citation: `CITATION.cff`.
