"""Grader for EXTRA50 (X-01..X-50) runs -- implements HARD50 DESIGN.md
section 7 rules verbatim, against gold_answers.json's own field-type vocabulary
(num/int/str/bool/null, "tol" not "rel_tol" -- see xb_common.py's num()/
integer()/exact()/boolean() builders).

Rules (HARD50 DESIGN.md Sec.7, applied here unchanged):
  - num:  relative tolerance (field's own "tol"); sign errors always fail
          (relative denominator -- see grade_field).
  - int:  structural count, exact match (no tolerance -- these are
          counted results: atoms, rings, n_flagged, etc.).
  - str:  case-insensitive; short/single-token golds (<=4 chars, or no
          internal whitespace at <=12 chars) require exact equality --
          the prompt prescribes the answer vocabulary for these. Longer
          golds use the GZ string_matches rule: >=34% of the gold's
          non-stopword tokens must appear as a substring in the answer.
  - bool: True/False equality; a string "true"/"false" (any case) is
          accepted and normalized first.
  - null: the model must emit JSON null -- a number, a string, or an
          excuse all fail. (No question in the current X-bank uses this
          type by design -- DESIGN.md Sec.0 item 2 -- but the rule is
          implemented for completeness/future questions.)
  - No partial credit at question level: ALL fields must pass.
  - Only fields in the gold dict are graded (== the FINAL_JSON contract
    fields, by xb_common.validate()'s own key-parity gate).

Usage:
  python grade_x50.py              # run self-tests only (no transcripts yet)
  python grade_x50.py --selftest   # same, explicit
"""
from __future__ import annotations

import json
import re
import unicodedata
from pathlib import Path

HERE = Path(__file__).resolve().parent
GOLD: dict = json.loads((HERE / "gold_answers.json").read_text(encoding="utf-8"))
EXPECTED_IDS = [f"X-{i:02d}" for i in range(1, 51)]

STOP = {"the", "and", "for", "with", "from", "that", "this", "not", "are", "was", "were",
        "has", "have", "had", "but", "its", "our", "your", "all", "any", "can", "may",
        "yes", "no", "to", "of", "in", "on", "at", "is", "be", "as", "it", "by", "or",
        "a", "an", "if", "then", "than", "true", "false", "null", "none", "not"}


# ---------------------------------------------------------------- normalize
def _norm(s) -> str:
    s = unicodedata.normalize("NFKC", str(s)).strip().lower()
    return re.sub(r"[\"'`]", "", s).strip()


def _tokens(s: str) -> list[str]:
    return [t for t in re.split(r"[^a-z0-9]+", s.lower()) if t and t not in STOP]


# ---------------------------------------------------------------- grading
def grade_field(fd: dict, ans) -> bool:
    """fd: one entry of a question's gold["fields"] dict. ans: the model's
    value for that key from its parsed FINAL_JSON (use "__MISSING__" sentinel
    if the key is absent)."""
    t = fd["type"]

    if t == "null":
        return ans is None

    if t == "num":
        if isinstance(ans, bool) or ans is None or ans == "__MISSING__":
            return False
        if isinstance(ans, str):
            try:
                ans = float(ans.strip().rstrip("%"))
            except Exception:
                return False
        if not isinstance(ans, (int, float)):
            return False
        g, tol = float(fd["value"]), float(fd.get("tol", 0.01))
        return abs(float(ans) - g) <= tol * abs(g)

    if t == "int":
        if isinstance(ans, bool) or ans is None or ans == "__MISSING__":
            return False
        try:
            return float(ans) == float(fd["value"])
        except Exception:
            return False

    if t == "bool":
        g = bool(fd["value"])
        if isinstance(ans, bool):
            return ans == g
        if isinstance(ans, str):
            na = _norm(ans)
            if na not in ("true", "false"):
                return False
            return (na == "true") == g
        return False

    # t == "str": GZ string_matches rule, verbatim from HARD50's grader.
    if ans is None or ans == "__MISSING__":
        return False
    gold = _norm(fd["value"])
    a = _norm(ans)
    if not gold:
        return False
    if len(gold) <= 4:
        # single short token (e.g. "I", "II", "SLS"): full equality required
        return a == gold
    if " " not in gold and len(gold) <= 12:
        # one short word gold (e.g. "rodenas", "exothermic"): must appear as
        # a whole token of the answer, not just a substring
        return gold == a or gold in re.split(r"[^a-z0-9]+", a)
    gt = _tokens(gold)
    if not gt:
        return False
    at = set(re.split(r"[^a-z0-9]+", a))
    # short tokens (<=3 chars, e.g. roman numerals "i"/"ii"/"iii") require a
    # whole-token match -- substring containment on a 1-3 char token produces
    # spurious hits against unrelated longer words (e.g. "i" inside
    # "definitely"). Longer tokens keep substring leniency for stemming.
    def _hit(tok: str) -> bool:
        if len(tok) <= 3:
            return tok in at
        return any(tok in w for w in at)
    hit = sum(1 for tok in gt if _hit(tok))
    return hit / len(gt) >= 0.34


def grade_question(qid: str, answer: dict) -> dict:
    """answer: the model's parsed FINAL_JSON dict for this question.
    Returns {"pass": bool, "fields": {key: bool}}."""
    gold = GOLD[qid]["fields"]
    fields = {}
    all_ok = True
    for key, fd in gold.items():
        ok = grade_field(fd, answer.get(key, "__MISSING__"))
        fields[key] = ok
        all_ok = all_ok and ok
    return {"pass": all_ok, "fields": fields}


def grade_all(answers: dict[str, dict]) -> dict:
    """answers: {qid: parsed FINAL_JSON dict}. Missing qids are scored as
    a full-field miss (not silently skipped)."""
    out = {}
    for qid in EXPECTED_IDS:
        ans = answers.get(qid, {})
        out[qid] = grade_question(qid, ans)
    return out


# ---------------------------------------------------------------- self-tests
def _gold_as_answer(qid: str) -> dict:
    """The gold's own field values, echoed back as a hypothetical perfect
    FINAL_JSON answer -- str(bool(...)) is fine since grade_field normalizes."""
    return {k: fd["value"] for k, fd in GOLD[qid]["fields"].items()}


def test_gold_back_in_passes_all_50():
    answers = {qid: _gold_as_answer(qid) for qid in EXPECTED_IDS}
    results = grade_all(answers)
    failed = [qid for qid, r in results.items() if not r["pass"]]
    assert not failed, f"gold-as-answer should be 50/50, but these failed: {failed}"
    assert len(results) == 50


def test_perturbed_numeric_field_fails():
    """Bump every numeric/int field in every question by 10x its tolerance
    (or a flat 50% for int/no-tol fields) -- each such question must now fail,
    and specifically on the field we perturbed."""
    checked = 0
    for qid in EXPECTED_IDS:
        gold_fields = GOLD[qid]["fields"]
        num_keys = [k for k, fd in gold_fields.items() if fd["type"] in ("num", "int")]
        if not num_keys:
            continue
        key = num_keys[0]
        fd = gold_fields[key]
        answer = _gold_as_answer(qid)
        true_val = float(fd["value"])
        tol = float(fd.get("tol", 0.0)) if fd["type"] == "num" else 0.0
        bump = max(tol, 0.05) * 10.0
        perturbed = true_val * (1.0 + bump) if true_val != 0 else bump
        answer[key] = perturbed
        result = grade_question(qid, answer)
        assert not result["fields"][key], f"{qid}.{key}: perturbed value unexpectedly passed"
        assert not result["pass"], f"{qid}: perturbed-field question unexpectedly passed overall"
        checked += 1
    assert checked >= 40, f"expected most of the 50 questions to have a numeric field, only checked {checked}"


def test_wrong_string_categorical_fails():
    """Every 'str'-type field fed an adversarial wrong token must fail --
    unless the wrong token happens to collide with the real gold string
    (checked explicitly and excluded)."""
    checked = 0
    for qid in EXPECTED_IDS:
        gold_fields = GOLD[qid]["fields"]
        str_keys = [k for k, fd in gold_fields.items() if fd["type"] == "str"]
        for key in str_keys:
            fd = gold_fields[key]
            wrong = "zzz_definitely_wrong_category_zzz"
            if _norm(wrong) == _norm(fd["value"]):
                continue  # pathological collision, skip rather than false-fail
            answer = _gold_as_answer(qid)
            answer[key] = wrong
            result = grade_question(qid, answer)
            assert not result["fields"][key], f"{qid}.{key}: wrong string unexpectedly passed"
            checked += 1
    assert checked >= 40, f"expected most questions to carry a str field, only checked {checked}"


def test_missing_key_fails_every_type():
    """Dropping a key entirely (as if the model never emitted it) must fail
    that field for every type, including bool/int (no silent default)."""
    for qid in EXPECTED_IDS:
        for key, fd in GOLD[qid]["fields"].items():
            assert not grade_field(fd, "__MISSING__"), f"{qid}.{key} ({fd['type']}): missing key unexpectedly passed"


def test_no_partial_credit():
    """A question with exactly one wrong field (all others perfect) must
    fail overall, never score partial credit."""
    checked = 0
    for qid in EXPECTED_IDS:
        gold_fields = GOLD[qid]["fields"]
        if len(gold_fields) < 2:
            continue
        answer = _gold_as_answer(qid)
        first_key = next(iter(gold_fields))
        fd = gold_fields[first_key]
        if fd["type"] in ("num", "int"):
            answer[first_key] = float(fd["value"]) * 5.0 + 1000.0
        elif fd["type"] == "bool":
            answer[first_key] = not bool(fd["value"])
        else:
            answer[first_key] = "zzz_wrong_zzz"
        result = grade_question(qid, answer)
        assert not result["pass"], f"{qid}: single wrong field should fail the whole question"
        n_ok = sum(1 for v in result["fields"].values() if v)
        assert n_ok == len(gold_fields) - 1, f"{qid}: expected exactly one field to fail"
        checked += 1
    assert checked >= 40


def run_self_tests() -> None:
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for t in tests:
        t()
        print(f"[ok] {t.__name__}")
    print(f"\nself-tests passed: {len(tests)}/{len(tests)}")


def main() -> None:
    assert sorted(GOLD) == EXPECTED_IDS, "gold_answers.json id set mismatch"
    run_self_tests()
    print("\nNo EXTRA50 model transcripts exist yet (bank not yet run) -- "
          "nothing further to grade. Re-run this file once answers.json-style "
          "parsed FINAL_JSON per model is available, via grade_all({qid: {...}}).")


if __name__ == "__main__":
    main()
