"""The 2026-10-20 repair grader: gold answers pass, wrong answers fail, rules are arm-independent."""
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("grade_repair", ROOT / "benchmark/repair_2026-10-20/grade_repair.py")
G = importlib.util.module_from_spec(spec)
spec.loader.exec_module(G)


def _final(d: dict) -> str:
    return "FINAL_JSON:\n" + json.dumps(d)


def _gz_gold(qid: str) -> dict:
    return {k: v for k, v in G.GZ_GOLD[qid]["gold"].items() if k not in G.GZ.META_KEYS}


@pytest.mark.parametrize("qid", ["GZ-15", "GZ-37", "GZ-38", "GZ-40"])
def test_groundzero_gold_back_passes(qid):
    ans = _gz_gold(qid)
    if qid == "GZ-38":
        ans["units_used"] = "µM"            # what a correct answer writes
    assert G.grade_gz("any_arm", qid, _final(ans))["verdict"] == "PASS"


@pytest.mark.parametrize("qid", ["H-36", "H-37", "H-38"])
def test_hard50_gold_back_passes(qid):
    ans = {k: fd["value"] for k, fd in G.H50.GOLD[qid]["fields"].items()}
    assert G.grade_h(qid, _final(ans))["verdict"] == "PASS"


def test_x19_gold_back_passes_and_uncorrected_value_fails():
    gold = {k: fd["value"] for k, fd in G.X50.GOLD["X-19"]["fields"].items()} if hasattr(G.X50, "GOLD") else None
    if gold is None:
        gold = {"mu_SI": 4.1e-08, "zeta_plain_mV": 52.5665011077, "zeta_corr_mV": 69.4550175843,
                "shift_pct": 32.1279067861, "larger": "corrected"}
    assert G.grade_x("X-19", _final(gold))["verdict"] == "PASS"
    wrong = dict(gold, zeta_corr_mV=60.18, shift_pct=14.48)
    assert G.grade_x("X-19", _final(wrong))["verdict"] == "FAIL"


@pytest.mark.parametrize("unit,ok", [("µM", True), ("μM", True), ("uM", True), ("micromolar", True),
                                     ("mM", False), ("M", False), ("", False)])
def test_gz38_unit_rule(unit, ok):
    assert G.FIELD_RULES[("GZ-38", "units_used")][0](unit) is ok


@pytest.mark.parametrize("ans,ok", [("B", True), ("Surfactant B", True), ("surfactant b.", True),
                                    ("A", False), ("Surfactant A", False), (None, False)])
def test_gz40_series_rule(ans, ok):
    assert G.FIELD_RULES[("GZ-40", "which_shows_dip_then_rise")][0](ans) is ok


@pytest.mark.parametrize("v,ok", [(2.51, True), (2.3, True), (2.75, True), (2.2, False), (2.0, False), (1.5, False)])
def test_gz15_cmc_uses_ten_percent(v, ok):
    assert G.FIELD_RULES[("GZ-15", "estimated_cmc_mM")][0](v) is ok


def test_gz36_rule_rejects_positive():
    rule = G.FIELD_RULES[("GZ-36", "slope_sign")][0]
    assert rule("negative") and rule("Negative (ln CMC decreases)") and not rule("positive")


def test_unjudged_free_text_is_never_silently_scored():
    gold = _gz_gold("GZ-39")
    r = G.grade_gz("an_arm_without_judgments", "GZ-39", _final(gold))
    assert r["verdict"] == "NEEDS_JUDGMENT"


def _h50_gold_answer(qid):
    """A correct answer built from the gold: null fields answered null, keyword fields answered
    with a sentence naming the first expected keyword."""
    out = {}
    for k, fd in G.H50.GOLD[qid]["fields"].items():
        if fd["type"] == "null":
            out[k] = None
        elif fd["type"] == "keywords_any":
            out[k] = f"needs a {fd['value'][0]} series"
        else:
            out[k] = fd["value"]
    return out


@pytest.mark.parametrize("qid", ["H-02", "H-03", "H-08", "H-16", "H-19", "H-31", "H-46", "H-47", "H-49"])
def test_counterion_corrected_hard50_gold_passes_and_inverse_answer_fails(qid):
    gold = _h50_gold_answer(qid)
    assert G.grade_h(qid, _final(gold))["verdict"] == "PASS"
    dg_key = next(k for k in gold if k.startswith("delta_g"))
    inverse = dict(gold)
    # an answer that still uses the inverse factor lands ~20-25% too weak and must fail
    inverse[dg_key] = gold[dg_key] * 0.78
    assert G.grade_h(qid, _final(inverse))["verdict"] == "FAIL"


def test_counterion_corrected_extrahard_golds():
    for qid in ("X-41", "X-47"):
        gold = {k: fd["value"] for k, fd in G.X50.GOLD[qid]["fields"].items()}
        assert G.grade_x(qid, _final(gold))["verdict"] == "PASS"
    old_x41 = {"dG_correct": -29.2597899505, "rel_diff_pct": 20.0, "claim_agrees": False,
               "what_wrong": "formula", "direction": "less_negative"}           # pre-correction gold
    assert G.grade_x("X-41", _final(old_x41))["verdict"] == "FAIL"


def test_dg_question_set_maps_ids():
    assert G.QIDS_DG["GZ16"] == "GZ-16" and G.QIDS_DG["H46"] == "H-46" and G.QIDS_DG["X41"] == "X-41"


def test_curly_quotes_are_normalised():
    ans = _gz_gold("GZ-37")
    text = _final(ans).replace('"', "“")
    assert G.grade_gz("any_arm", "GZ-37", text)["verdict"] == "PASS"


@pytest.mark.parametrize("ans,ok", [("langmuir", True), ("Langmuir (BIC-selected over Frumkin)", True),
                                    ("Langmuir/Szyszkowski (BIC-selected; data-limited)", True), ("Szyszkowski", True),
                                    ("frumkin (BIC-selected over Langmuir)", False),
                                    ("Gibbs adsorption from linear surface tension vs ln(concentration)", False),
                                    (None, False)])
def test_gz20_isotherm_choice_rule(ans, ok):
    assert G.FIELD_RULES[("GZ-20", "isotherm_model")][0](ans) is ok