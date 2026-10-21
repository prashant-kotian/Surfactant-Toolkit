"""The 2026-10-20 corrections layer used by build_results.py: rules are arm-independent and logged."""
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("corrections", ROOT / "benchmark/repair_2026-10-20/corrections.py")
C = importlib.util.module_from_spec(spec)
spec.loader.exec_module(C)


def _rec(verdict, **fields):
    return {"verdict": verdict, "fields": {k: dict(v) for k, v in fields.items()}}


def test_choice_rule_fixes_both_directions():
    obj = {"GZ-20": _rec("FAIL", isotherm_model={"ok": False, "model": "langmuir"}, cmc_mM={"ok": True}),
           "GZ-29": _rec("PASS", selected_model={"ok": True, "model": "langmuir (BIC-selected over Frumkin)"})}
    out, log = C.choice_regrade(obj)
    assert out["GZ-20"]["verdict"] == "PASS" and out["GZ-29"]["verdict"] == "FAIL"
    assert obj["GZ-20"]["verdict"] == "FAIL"                        # input untouched
    assert {x["qid"] for x in log} == {"GZ-20", "GZ-29"}


def test_choice_rule_keeps_other_failures_and_reads_got_key():
    obj = {"GZ-01": _rec("FAIL", isotherm_model={"ok": False, "got": "Langmuir"}, cmc_mM={"ok": False})}
    out, _ = C.choice_regrade(obj)
    assert out["GZ-01"]["fields"]["isotherm_model"]["ok"] is True
    assert out["GZ-01"]["verdict"] == "FAIL"                        # another field still fails


def test_choice_rule_leaves_irreproducible_records_alone():
    obj = {"GZ-20": _rec("PASS", isotherm_model={"ok": False, "model": "frumkin"})}   # verdict contradicts flags
    out, log = C.choice_regrade(obj)
    assert out == obj and log[0]["change"] == "none"


def test_apply_excludes_stale_counterion_verdicts_without_rerun(monkeypatch):
    monkeypatch.setattr(C, "DG_RES", {})
    monkeypatch.setattr(C, "REPAIR_RES", {})
    v, ex, log = C.apply("groundzero", "Qwen", "unaug", {"GZ-16": True, "GZ-20": False, "GZ-01": True}, {}, [])
    assert "GZ-16" not in v and "GZ-16" in ex and "GZ-20" in ex and v["GZ-01"] is True


def test_apply_uses_reruns_and_original_chart_runs(monkeypatch):
    monkeypatch.setattr(C, "DG_RES", {"deepseek_aug": {"GZ-16": {"verdict": "FAIL"}}})
    monkeypatch.setattr(C, "REPAIR_RES", {"deepseek_aug": {"GZ-36": {"verdict": "PASS"}}})
    base_ex = ["GZ-15", "GZ-36", "GZ-26"]
    v, ex, _ = C.apply("groundzero", "DeepSeek", "aug", {"GZ-15": True, "GZ-16": True, "GZ-36": False}, {}, base_ex)
    assert v["GZ-15"] is True and "GZ-15" not in ex                 # chart delivered in the original run
    assert v["GZ-36"] is True and v["GZ-16"] is False                # re-runs replace the originals
    assert "GZ-26" in ex                                             # retired question untouched


def test_gemini_cli_extrahard_items_stay_excluded(monkeypatch):
    monkeypatch.setattr(C, "DG_RES", {})
    monkeypatch.setattr(C, "REPAIR_RES", {})
    v, ex, _ = C.apply("extrahard", "Gemini-CLI", "aug", {"X-19": True, "X-41": True, "X-47": False}, {}, ["X-19"])
    assert set(ex) == {"X-19", "X-41", "X-47"} and not ({"X-41", "X-47"} & set(v))


def test_h46_unaug_regraded_from_stored_answer_against_corrected_gold(monkeypatch):
    monkeypatch.setattr(C, "DG_RES", {})
    G = C._load_grader()
    gold = {k: fd["value"] for k, fd in G.H50.GOLD["H-46"]["fields"].items()}
    good = {"fields": {k: {"got": val} for k, val in gold.items()}}
    stale = {"fields": {k: {"got": (-29.953153909017875 if k.startswith("delta_g") else val)} for k, val in gold.items()}}
    v, ex, _ = C.apply("hard50", "Gemini-CLI", "unaug", {"H-46": True}, {"H-46": good}, [])
    assert v["H-46"] is True
    v, ex, _ = C.apply("hard50", "Gemini-CLI", "unaug", {"H-46": True}, {"H-46": stale}, [])
    assert v["H-46"] is False                                        # the inverse-factor value now fails
