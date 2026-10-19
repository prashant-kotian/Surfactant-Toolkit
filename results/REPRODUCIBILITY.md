# Reproducibility verification (fresh clone, 2026-10-18)

Verified by cloning the source repositories at the freeze tags to an
empty Linux directory (Python 3.13) and running the full chain.

## Commands (exact)
```
git clone <surfactantkit> sk && git clone <extra-surfactantkit-work> work
pip install rdkit "mcp[cli]"                 # pyproject deps for the test suite
export PYTHONPATH=$PWD/sk/src
cd work
python3 benchmark/paper3_groundzero/grade_unaugmented_gz50.py
python3 benchmark/paper3_groundzero/grade_session_new_gz50.py
python3 benchmark/paper3_hard50/grade_session_new_h50.py
python3 benchmark/paper3_extra50/grade_unaugmented_5models.py
python3 benchmark/paper3_extra50/grade_augmented_claude_deepseek.py
python3 benchmark/results_master.py
python3 benchmark/make_tables.py
cd ../sk && PYTHONPATH=src python3 -m pytest -q
```

## Results
- Every grader ran from the clone; regenerated grading JSONs are **JSON-level
  equal** to the committed artifacts (6/6 files compared by parsed content).
- `results_master.py` -> `cells=30 gaps=1`, counts identical
  (X50 unaug 32/33/28/21/27; aug 50/48/50/49; H50new 45/44 aug; agy GZ
  42/43 unaug post-unmasking).
- `RESULTS_PER_Q.csv` and `RESULTS_SECTION.md`: **byte-identical**.
  `RESULTS_TABLES.md` differs only in the provenance-ledger `source_sha`
  column, because the freshly rewritten grading JSONs differ from the
  committed (backend-machine-written) bytes while parsing equal —
  serialization, not data.
- Toolkit test suite on the fresh clone: **823 passed** (same run, same
  local dir: 823 passed).
- `toolkit_state` resolves from sibling `surfactantkit`/`SurfactantKit`/
  `surfactantkit-computation` or `$TOOLKIT_DIR` (hardened this commit).

## Known, accepted non-determinism (documented, not chased)
- `gold_state` fields derive from source-file mtimes of the machine that
  last ran a grader; a fresh re-grade on a later date legitimately flips
  flags to post-fix.
- Development scripts (`scratch_build_*.py`, `run_agy_overnight_windows.py`,
  `benchmark/development/merge_audit/`) intentionally retains the benchmark
  machine's `H:\CodeProjects` paths; **not part of the reproducibility
  chain** (excluded-at-freeze decision, per build plan).
- `ungradeable`-record rule and exclusions live in `results_master.py` as
  data; changing them changes denominators by design, not by drift.

## Verdict
**PASS.** Manuscript numbers = this clone's output. Freeze tag:
`paper3-freeze-2026-10-18` (work repo) + matching toolkit tag.
