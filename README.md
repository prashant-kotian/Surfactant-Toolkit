# Surfactant Toolkit

A computation toolkit for surfactant and interfacial chemistry, with
**SurfBench-200**, the benchmark used to evaluate large language models with
and without tool access (Paper 3). This repository is the complete submission
artifact: the library, its test suite, the full benchmark with recomputed
reference answers (golds), raw model-answer transcripts, the graders, and the
results pipeline. Every number reported in the manuscript is regenerated from
this tree; none is typed by hand.

## The benchmark: SurfBench-200

SurfBench-200 is a single 200-question benchmark for toolkit-augmented vs
unaugmented reasoning in surfactant science. Each question `SB-001..SB-200`
carries two classification labels used by the result tables:

- **skill type** — `conceptual` (definitions, identification, conventions),
  `computational` (numeric execution through formulas and fits), and
  `analytical` (multi-step reasoning: model selection, cross-source audits,
  reconciliation of conflicting data);
- **sub-concept** — 16 topic keys (cmc measurement, mixed micellization,
  adsorption isotherms, emulsion/HLD, colloid electrostaticics, micelle
  structure, packing/morphology, kinetics/transport, thermodynamics,
  solubilization, and related areas).

Labels are assigned deterministically from each question's gold metadata and
text by `benchmark/build_question_set.py` (rule table plus a curated override
file); the resulting table is committed as `benchmark/question_set.json`.

For authoring convenience the questions were written as four 50-question
working groups; those names are internal provenance only
(`benchmark/source_sets/`) and all reporting uses the merged benchmark.

| Deliverable | Path |
|---|---|
| Headline pass rates, augmented vs unaugmented | `results/TABLE200.md` |
| Classification matrices + generated conclusions | `results/CLASSIFICATION200.md` |
| Per-question table (model × condition × question) | `results/per_question_results_200.csv` |
| Machine-readable grand/class data | `results/bench200.json` |
| Per-source-set diagnostic view and grading ledger | `results/tables.md`, `results/results_summary.md` |

To rename the benchmark, change the `BENCH_NAME` constant in
`benchmark/build_question_set.py` and re-run the pipeline.

## Layout

```
src/surfactantkit/      toolkit library (CMC, mixed micellization, thermodynamics,
                        model selection, adsorption, HLD, QSPR connectors) + MCP server
tests/                  test suite (832 tests)
benchmark/
  question_set.json     the benchmark: 200 ids, skill types, sub-concepts, exclusions
  build_question_set.py assembles question_set.json from the source-set golds
  build_bench200.py     grand + classification tables from the grading JSONs
  source_sets/          authoring provenance (four working groups, merged above)
    groundzero/         SB-001..050  (GZ-01..50): definitional / standard practice
    hard50/             SB-051..100  (H-01..50): multi-step computation
    extrahard/          SB-101..150  (X-01..50): adversarial
    orig_hard50/        SB-151..200  (Q01..Q50): toolkit-execution items from the
                        project's first benchmark release
  mcp_harness/          Model Context Protocol harness used for the augmented condition
  transcripts/          raw model answers, one file per run/question
  build_results.py      collects the graded cells into results/results.json (+ ledger)
  build_tables.py       renders the per-set diagnostic tables from that JSON
docs/                   validation-tier audit table (guarded by tests/test_validation_tiers.py)
results/                generated: TABLE200.md, CLASSIFICATION200.md, bench200.json,
                        figs_data/*.tsv (paste-ready figure sources),
                        per_question_results_200.csv, tables.md/.tex,
                        per_question_results.csv, REPRODUCIBILITY.md
```

Each source-set directory contains its question bank, `gold_answers.json`
(reference values, with `remediation_note` fields where golds were corrected),
the grader scripts, and per-run grading JSONs.

## Install and run

```bash
pip install rdkit "mcp[cli]"         # only external deps (see pyproject extras)
PYTHONPATH=src python3 -m pytest -q  # 832 passed
```

Reproduce the reported numbers (deterministic given the committed transcripts)
— see `REPRODUCE.md` for the exact chain:

```bash
python3 benchmark/build_results.py     # grading JSONs -> results/results.json + ledger
python3 benchmark/build_tables.py      # -> results/tables.md, tables.tex, per-question CSV
python3 benchmark/build_question_set.py  # golds -> benchmark/question_set.json
python3 benchmark/build_bench200.py    # -> TABLE200.md, CLASSIFICATION200.md, bench200.json
```

The four `build_*` steps consume the committed grading JSONs; REPRODUCE.md
lists the full chain including the graders that regenerate those JSONs.
`results/per_question_results.csv` must match its committed
copy byte for byte; see `results/REPRODUCIBILITY.md` for the recorded
verification.

## Data and provenance policy

- Reference values are recomputed by the toolkit (`src/surfactantkit/`); gold
  corrections made during audit are recorded as machine-readable
  `remediation_note` fields on the affected entries, with cause.
- Questions excluded from scoring (undeliverable multimodal inputs, retired
  ambiguous wording, and items affected by a pre-fix question packet) are
  declared in `build_results.py` and `question_set.json` as data, applied
  identically to every model and both conditions, with the denominator
  reported (n after exclusions).
- Answers lost to capture defects (no recoverable record) are treated as
  ungradable and excluded from denominators, not scored against the model;
  cells covering fewer than the full benchmark are flagged `partial coverage`.
- Every grading input carries a content hash in the provenance ledger
  (`results/results_summary.md`).

License: MIT. Citation: `CITATION.cff`.
