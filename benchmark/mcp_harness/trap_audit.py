"""One-off targeted audit: run every no_solution and unit_trap question in the
full Tier-1 bank against claude-opus, unaugmented condition only. Purpose:
find out whether these trap types are a real, full-population sub-50% Claude
weakness (not just an artifact of the tiny n=2-3 samples seen in prior mixed
pilots). Reuses run_pilot.py's run_one() directly so grading/extraction/
provider logic is identical to the main harness -- no reimplementation."""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_pilot  # noqa: E402  (triggers _load_env() as a side effect)

QB_PATH = Path(__file__).resolve().parent.parent / "question_bank_tier1.json"
qb = json.loads(QB_PATH.read_text())

subset = [q for q in qb if q.get("trap_type") in ("no_solution", "unit_trap")]
print(f"running {len(subset)} questions x claude-opus x unaugmented")

rows = []
for q in subset:
    r = run_pilot.run_one("claude-opus", q, "unaugmented")
    flag = "OK " if r["correct"] else "XX "
    if r["error"]:
        flag = "ERR"
    print(f"  [{flag}] {q['id']:7s} {q['category']} {q['trap_type']:12s} "
          f"{r['seconds']}s {r['grade_reason']}")
    rows.append(r)

out = run_pilot.RESULTS_DIR / f"trap_audit_{__import__('datetime').datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
run_pilot.RESULTS_DIR.mkdir(exist_ok=True)
out.write_text(json.dumps({"models": ["claude-opus"], "conditions": ["unaugmented"],
                            "subset_ids": [q["id"] for q in subset], "rows": rows}, indent=2))

from collections import defaultdict
by_trap = defaultdict(lambda: [0, 0])
for r in rows:
    by_trap[r["trap_type"]][1] += 1
    if r["correct"]:
        by_trap[r["trap_type"]][0] += 1
print("\n=== accuracy by trap_type (claude-opus, unaugmented) ===")
for t, (c, n) in sorted(by_trap.items()):
    print(f"  {t:14s} {c}/{n} = {100*c/n:.0f}%")
print(f"\nwrote {out}")
