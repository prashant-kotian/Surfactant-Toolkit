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


def test_every_departure_from_the_standard_exclusions_is_logged():
    # Since the 2026-10-20 corrections, exclusions can differ between arms (re-run answers available
    # for some arms only, re-runs still pending). Every difference must carry a logged reason, and
    # an excluded question is never scored.
    for key, cell in CELLS.items():
        bank = key.split("|")[0]
        standard = {q for _, qids, _ in br.EXCLUSIONS.get(bank, []) for q in qids}
        log = cell.get("corrections", [])
        included = {c["qid"] for c in log if str(c.get("change", "")).startswith(("included", "replaced",
                                                                                    "regraded"))}
        excluded_logged = {c["qid"] for c in log if c.get("change") in ("excluded", "stays excluded")}
        excluded = set(cell["excluded"])
        assert not excluded & set(cell["per_q"]), key
        assert standard - excluded <= included, (key, standard - excluded - included)
        assert excluded - standard <= excluded_logged, (key, excluded - standard - excluded_logged)


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
