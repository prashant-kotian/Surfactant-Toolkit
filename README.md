# Surfactant Toolkit Paper-3 artifact: toolkit + benchmark + results chain

Complete, self-contained artifact for *Paper 3* (agentic LLM evaluation on
surfactant/interfacial-chemistry computation): the **SurfToolkit** library, its
test suite, three 50-question benchmark banks with toolkit-recomputed golds,
the per-question run transcripts (raw model outputs), the committed graders,
and a single machine-graded results pipeline that emits every number in the
manuscript.

```
src/surfactantkit/          the toolkit (library + MCP server; see docs)
tests/                      823-test suite   (pip install rdkit "mcp[cli]"; PYTHONPATH=src pytest -q)
benchmark/paper3_groundzero/    GroundZero-50 bank: gold, graders, transcripts (incl. image-question charts)
benchmark/paper3_hard50/        Hard50 bank: same layout
benchmark/paper3_extra50/       ExtraHard-50 bank: same layout
benchmark/harness/              MCP harness (augmented condition)
session_pull/                   chat-UI transcripts, one .txt per question (+ .RAW captures for recaptured items)
new extrahard set output-.../   ExtraHard unaugmented transcripts + grading JSON
benchmark/results_master.py     assembles ALL reported cells from committed JSONs, with sha provenance
benchmark/make_tables.py        generates manuscript tables + prose-section from RESULTS_MASTER.json
paper3/                         manuscript draft (skeleton + rendered section), SI tables, gauntlet record
```

**Policy that shapes this repo:** every manuscript number is generated — no
hand-typed values, report-only figures are refused and listed as gaps, and
exclusions are declared data applied symmetrically across a bank/condition.
Directory names with historical typos (`unagumented 18th sept`) are preserved
verbatim because committed graders read those exact paths; provenance beats
typo fixes here. The internal engineering history (audit briefs, agent
run-logs, scratch) stays in the authors' private work repo; see
`MANIFEST.md` for what was excluded, and `paper3/GAUNTLET.md` for the
fresh-clone reproducibility record this artifact passed.

Citation: `CITATION.cff`. License: MIT (`LICENSE`).
