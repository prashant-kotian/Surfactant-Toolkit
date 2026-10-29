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
    monkeypatch.setattr(C, "FILL_RES", {})
    v, ex, log = C.apply("groundzero", "Qwen", "unaug", {"GZ-16": True, "GZ-20": False, "GZ-01": True}, {}, [])
    assert "GZ-16" not in v and "GZ-16" in ex and "GZ-20" in ex and v["GZ-01"] is True


def test_apply_uses_reruns_and_original_chart_runs(monkeypatch):
    monkeypatch.setattr(C, "DG_RES", {"deepseek_aug": {"GZ-16": {"verdict": "FAIL"}}})
    monkeypatch.setattr(C, "REPAIR_RES", {"deepseek_aug": {"GZ-36": {"verdict": "PASS"}}})
    monkeypatch.setattr(C, "FILL_RES", {})                          # isolate from the live net-200 fill results
    base_ex = ["GZ-15", "GZ-36", "GZ-26"]
    v, ex, _ = C.apply("groundzero", "DeepSeek", "aug", {"GZ-15": True, "GZ-16": True, "GZ-36": False}, {}, base_ex)
    assert v["GZ-15"] is True and "GZ-15" not in ex                 # chart delivered in the original run
    assert v["GZ-36"] is True and v["GZ-16"] is False                # re-runs replace the originals
    assert "GZ-26" in ex                                             # retired question untouched


def test_gemini_cli_extrahard_items_stay_excluded(monkeypatch):
    monkeypatch.setattr(C, "DG_RES", {})
    monkeypatch.setattr(C, "REPAIR_RES", {})
    monkeypatch.setattr(C, "FILL_RES", {})
    v, ex, _ = C.apply("extrahard", "Gemini", "aug", {"X-19": True, "X-41": True, "X-47": False}, {}, ["X-19"])
    assert set(ex) == {"X-19", "X-41", "X-47"} and not ({"X-41", "X-47"} & set(v))


def test_h46_unaug_regraded_from_stored_answer_against_corrected_gold(monkeypatch):
    monkeypatch.setattr(C, "DG_RES", {})
    monkeypatch.setattr(C, "FILL_RES", {})
    G = C._load_grader()
    gold = {k: fd["value"] for k, fd in G.H50.GOLD["H-46"]["fields"].items()}
    good = {"fields": {k: {"got": val} for k, val in gold.items()}}
    stale = {"fields": {k: {"got": (-29.953153909017875 if k.startswith("delta_g") else val)} for k, val in gold.items()}}
    v, ex, _ = C.apply("hard50", "Claude", "unaug", {"H-46": True}, {"H-46": good}, [])
    assert v["H-46"] is True
    v, ex, _ = C.apply("hard50", "Claude", "unaug", {"H-46": True}, {"H-46": stale}, [])
    assert v["H-46"] is False                                        # the inverse-factor value now fails


def test_fill_overlay_lifts_gz26_and_fills_missing_records(monkeypatch):
    monkeypatch.setattr(C, "DG_RES", {})
    monkeypatch.setattr(C, "REPAIR_RES", {})
    monkeypatch.setattr(C, "FILL_RES", {"gemini_unaug": {"GZ-26": {"verdict": "FAIL"}, "GZ-01": {"verdict": "PASS"},
                                                         "GZ-05": {"verdict": "NO_JSON"}}})
    v, ex, log = C.apply("groundzero", "Gemini", "unaug", {"GZ-02": True}, {}, ["GZ-26"])
    assert v["GZ-26"] is False and "GZ-26" not in ex           # retired question now scored from its re-run
    assert v["GZ-01"] is True                                   # never-captured record filled
    assert "GZ-05" not in v                                     # ungradable re-run: still absent, not a fail


def test_fill_overlay_replaces_gemini_cli_extrahard_exclusions(monkeypatch):
    monkeypatch.setattr(C, "DG_RES", {})
    monkeypatch.setattr(C, "REPAIR_RES", {})
    monkeypatch.setattr(C, "FILL_RES", {"agy_aug": {"X-41": {"verdict": "PASS"}}})
    v, ex, log = C.apply("extrahard", "Gemini", "aug", {"X-41": False, "X-47": True}, {}, ["X-19"])
    assert v["X-41"] is True and "X-41" not in ex and {"X-19", "X-47"} <= set(ex)
    assert not [e for e in log if e["qid"] == "X-41" and e["change"] in ("excluded", "stays excluded")]


def test_orig_overlay_only_returns_gradable_q_items(monkeypatch):
    monkeypatch.setattr(C, "FILL_RES", {"chatgpt_unaug": {"Q24": {"verdict": "PASS"}, "Q25": {"verdict": "NO_JSON"},
                                                          "GZ-26": {"verdict": "FAIL"}}})
    assert C.orig_overlay("ChatGPT", "unaug") == {"Q24": "PASS"}

def test_orig_regrade_is_inert_until_written_then_overrides(monkeypatch):
    # nothing written yet (pending augmented re-runs) -> no override for any arm
    if not (ROOT / "benchmark/repair_2026-10-20/q47_regrade_results.json").exists():
        assert C.orig_regrade("Claude", "aug") == {}
    monkeypatch.setattr(C, "Q3_RES", {})   # the 2026-10-28 one-target-per-question results merge last
    monkeypatch.setattr(C, "Q17_RES", {})
    monkeypatch.setattr(C, "Q47_RES", {"agy_aug": {"Q47": {"verdict": "PASS", "detail": ""}},
                                       "claude_unaug": {"Q47": {"verdict": "FAIL",
                                                                "detail": "part1 got 38.056 vs gold 57.982"}}})
    assert C.orig_regrade("Gemini", "aug") == {"Q47": ("PASS", "")}
    assert C.orig_regrade("Claude", "unaug")["Q47"][0] == "FAIL"
    assert C.orig_regrade("ChatGPT", "aug") == {}


def test_orig_regrade_merges_q17_augmented_only(monkeypatch):
    monkeypatch.setattr(C, "Q3_RES", {})
    monkeypatch.setattr(C, "Q47_RES", {"claude_aug": {"Q47": {"verdict": "PASS", "detail": ""}}})
    monkeypatch.setattr(C, "Q17_RES", {"claude_aug": {"Q17": {"verdict": "PASS", "detail": ""}}})
    assert C.orig_regrade("Claude", "aug") == {"Q47": ("PASS", ""), "Q17": ("PASS", "")}
    assert "Q17" not in C.orig_regrade("Claude", "unaug")


def test_q17_reference_values_per_wording():
    # unaugmented wording states its own weights (ether O 1.3) -> frozen 33.29; augmented wording
    # is the release tool's output (EO unit 0.33)
    spec = importlib.util.spec_from_file_location("q17_regrade", ROOT / "benchmark/repair_2026-10-20/q17_regrade.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert mod.unaugmented_reference() == 33.29
    assert mod.reference_value() == 9.62


def test_q16_q17_q47_one_reference_per_question():
    # one target per question in both wordings (audit 2026-10-28), computed from the release toolkit
    spec = importlib.util.spec_from_file_location(
        "orig_q16_q17_q47_grade", ROOT / "benchmark/repair_2026-10-20/orig_q16_q17_q47_grade.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert mod.references() == {"Q16": (2.8568, 4), "Q17": (9.62, 2), "Q47": (57.93, 2)}


def test_orig_regrade_q16_q17_q47_results_take_precedence(monkeypatch):
    monkeypatch.setattr(C, "Q47_RES", {"claude_unaug": {"Q47": {"verdict": "FAIL", "detail": "old"}}})
    monkeypatch.setattr(C, "Q17_RES", {})
    monkeypatch.setattr(C, "Q3_RES", {"claude_unaug": {"Q47": {"verdict": "PASS", "detail": ""},
                                                        "Q16": {"verdict": "PASS", "detail": ""}}})
    assert C.orig_regrade("Claude", "unaug") == {"Q47": ("PASS", ""), "Q16": ("PASS", "")}
