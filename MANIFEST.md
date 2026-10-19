# MANIFEST — contents policy of this artifact

INCLUDED (submission-relevant, reviewer-runnable):
- Toolkit source + full test suite + pyproject (freeze tag `paper3-freeze-2026-10-18` of `surfactantkit`).
- All three question banks: question texts, gold JSONs (with per-field
  `remediation_note` entries where golds were corrected 2026-10-18), bank
  builders, image assets for excluded multimodal questions.
- Every raw transcript the results chain reads (chat UI .txt per question,
  re-captured .RAW files, agy CLI batch outputs, parsed answer JSONs).
- All grader scripts (each runs from a fresh clone) + the two committed
  generators (`results_master.py`, `make_tables.py`) + their outputs
  (RESULTS_MASTER.json/.md, paper3/RESULTS_TABLES.*, RESULTS_PER_Q.csv).
- Manuscript draft v1 (§4–§7) + freeze gauntlet record (paper3/GAUNTLET.md).

EXCLUDED (internal working record; available from authors on request):
- Agent/brief correspondence (BACKEND_BRIEF, ARENA_AGENT_BRIEF, MOVES,
  MORNING_REPORT), FIXLOG adjudication narrative, RUNS_REPORT (contains
  superseded counts).
- Dev-time scratch: root `scratch_build_*`, `run_agy_overnight_*` (machine
  paths for the benchmark box), `benchmark/development/` (MCP dev logs,
  merge-audit scripts), per-bank research notes (`*_RESEARCH.md`,
  `GENUINE_BOTTLENECK_AUDIT.md`, `BOTTLENECK_*`, DESIGN/GRADE_REPORT
  markdowns), `paste*` staging dirs, superseded duplicate partial transcript
  dirs, and the quarantined `grade_runs_x50.py` wrapper (known header-parse
  bug; never used for reported numbers — retained privately as a cautionary
  artifact).
- `benchmark/paper3_scale50/` (development category builders for a bank not
  part of the three reported here).

NOTHING in the reported-results path was removed: re-running every generator
in this repo reproduces `RESULTS_PER_Q.csv` byte-identically (see
REPRODUCE.md).
