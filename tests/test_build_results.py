"""Regression tests for benchmark/build_results.py: exclusion symmetry and verified gold state."""
import importlib.util
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "benchmark"))
_SPEC = importlib.util.spec_from_file_location("build_results", ROOT / "benchmark" / "build_results.py")
br = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(br)

CELLS = json.loads((ROOT / "results" / "results.json").read_text(encoding="utf-8"))["cells"]


def test_every_arm_excludes_the_same_questions():
    # an unaug/aug pair scored on different question sets is not a paired comparison
    for key, cell in CELLS.items():
        bank = key.split("|")[0]
        want = {q for _, qids, _ in br.EXCLUSIONS.get(bank, []) for q in qids}
        assert set(cell["excluded"]) == want, key
        assert not want & set(cell["per_q"]), key


def test_no_reported_cell_is_graded_against_superseded_gold():
    assert not [k for k, c in CELLS.items() if c["gold_state"] == "PRE-FIX gold"]


def test_verify_gold_detects_superseded_values():
    cur = br.current_gold("hard50", "H-26")
    stale = {"H-26": {"fields": {"total_carbons": {"ok": True, "gold": cur["total_carbons"] - 1}}}}
    fresh = {"H-26": {"fields": {"total_carbons": {"ok": True, "gold": cur["total_carbons"]}}}}
    assert br.verify_gold("hard50", stale) == "PRE-FIX gold"
    assert br.verify_gold("hard50", fresh) == "post-fix"
    assert br.verify_gold("hard50", {"passed": [], "failed": []}) == "not field-verifiable"


def test_cell_notes_belong_to_their_own_model():
    for key, cell in CELLS.items():
        note = cell.get("note") or ""
        if note.startswith("Claude:"):
            assert key.split("|")[1] == "Claude", key
