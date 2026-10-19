"""Select a small, deterministic, stratified pilot subset from the Tier 1 bank.

The pilot's job is to validate the whole pipeline (API call -> tool use ->
answer extraction -> grading) end to end, cheaply -- not to measure final
scores. So the subset is small but deliberately spans all 9 categories and all
trap types, since those are the parts most likely to expose a plumbing bug
(no-solution traps stress the grader; multi-tool chains stress tool-calling;
unit traps stress the unit-conversion path).
"""

from __future__ import annotations
import json
import random
from pathlib import Path

BANK = Path(__file__).resolve().parent.parent / "question_bank_tier1.json"


def select(n: int = 20, seed: int = 7) -> list[dict]:
    bank = json.loads(BANK.read_text())
    rng = random.Random(seed)
    chosen: dict[str, dict] = {}

    def take(pool, k):
        pool = [q for q in pool if q["id"] not in chosen]
        for q in rng.sample(pool, min(k, len(pool))):
            chosen[q["id"]] = q

    # guarantee coverage of the trap types most likely to break the pipeline
    take([q for q in bank if q.get("trap_type") == "no_solution"], 3)
    take([q for q in bank if q.get("trap_type") == "unit_trap"], 3)
    take([q for q in bank if q.get("trap_type") == "multi_tool_chain"], 4)
    # guarantee at least one question from each category A..I (only if that
    # category isn't already represented by a trap picked above)
    for cat in "ABCDEFGHI":
        if not any(q["category"] == cat for q in chosen.values()):
            take([q for q in bank if q["category"] == cat], 1)
    # fill the rest at random up to n; never trim guaranteed coverage back out
    take(bank, max(0, n - len(chosen)))

    return sorted(chosen.values(), key=lambda q: q["id"])


if __name__ == "__main__":
    from collections import Counter
    sub = select()
    print(f"pilot subset: {len(sub)} questions")
    print("by category:", dict(sorted(Counter(q["category"] for q in sub).items())))
    print("by trap:", dict(Counter(q.get("trap_type", "none") for q in sub)))
    print("by difficulty:", dict(Counter(q["difficulty"] for q in sub)))
    print("tools touched:", len({t for q in sub for t in q["tools_required"]}))
    for q in sub:
        print(f"  {q['id']} [{q.get('trap_type','none')}] {q['subcategory']}")
