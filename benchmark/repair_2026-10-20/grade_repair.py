"""Grade the 2026-10-20 re-runs of the excluded questions.

The 9 chart questions (GZ-15, GZ-36..GZ-40, H-36..H-38) were excluded because
the original harnesses never delivered the chart image; X-19 because its
question body carried a kappa*a typo (1.45 vs 14.5) at run time. They were
re-run on 2026-10-20 with the chart attached / the corrected X-19 text, under
the same model settings and MCP conditions as the original arms, with every
condition verified on the live page before sending (see answers/*/*.meta.json).

Each answer is graded by its own benchmark's existing grader, unchanged:
  GZ-*  groundzero/grade_unaugmented.grade_field (+ the long-prose judgment rule)
  H-*   hard50/grade_unaugmented_legacy.grade_field
  X-19  extrahard/grade_x50.grade_question
Curly quotes are normalised first, exactly as groundzero/hard50 grade_chat_runs do.

Long free-text gold fields (> 80 chars) cannot be keyword-matched; as in the
groundzero grader they get one written resolution per (arm, question, field) in
JUDGMENTS below. Anything unresolved is reported as NEEDS_JUDGMENT, never
silently passed or failed.

Usage: python grade_repair.py            -> writes repair_grade_results.json
"""
from __future__ import annotations

import importlib.util
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SETS = HERE.parent / "source_sets"
ANSWERS = HERE / "answers"

ARMS = {   # folder -> (model label used in build_results, condition)
    "claude_unaug": ("Claude", "unaug"), "claude_aug": ("Claude", "aug"),
    "chatgpt_unaug": ("ChatGPT", "unaug"), "chatgpt_aug": ("ChatGPT", "aug"),
    "deepseek_unaug": ("DeepSeek", "unaug"), "deepseek_aug": ("DeepSeek", "aug"),
    "qwen_unaug": ("Qwen", "unaug"), "gemini_unaug": ("Gemini", "unaug"),
    "agy_unaug": ("Gemini-CLI", "unaug"), "agy_aug": ("Gemini-CLI", "aug"),   # counterion re-runs only
}
QIDS = {"GZ15": "GZ-15", "GZ36": "GZ-36", "GZ37": "GZ-37", "GZ38": "GZ-38", "GZ39": "GZ-39",
        "GZ40": "GZ-40", "H36": "H-36", "H37": "H-37", "H38": "H-38", "X19": "X-19"}
# second run set: questions re-asked after the 2026-10-20 counterion-factor correction
QIDS_DG = {q: f"{q[:2]}-{q[2:]}" if q.startswith("GZ") else f"{q[0]}-{q[1:]}"
           for q in ("GZ16", "GZ20", "H02", "H03", "H08", "H16", "H19", "H31", "H38", "H46", "H47", "H49",
                     "X41", "X47")}


def _load(name: str, path: Path):
    sys.path.insert(0, str(path.parent))
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    sys.path.pop(0)
    return mod


GZ = _load("gz_grader", SETS / "groundzero" / "grade_unaugmented.py")
H50 = _load("h50_grader", SETS / "hard50" / "grade_unaugmented_legacy.py")
X50 = _load("x50_grader", SETS / "extrahard" / "grade_x50.py")

with open(SETS / "groundzero" / "gold_answers.json", encoding="utf-8") as f:
    GZ_GOLD = json.load(f)
GZ_TEMPLATES = GZ.extract_templates()

def _norm(v) -> str:
    return " ".join(str(v or "").lower().replace("μ", "u").replace("µ", "u").split())


def _letter(v) -> str:
    s = _norm(v).replace("surfactant", "").replace("curve", "").strip(" .:")
    return s.upper()


# Field rules: applied identically to EVERY arm, before the shared comparator. Each exists
# because the shared keyword matcher cannot score this gold string correctly (it strips
# non-ASCII such as the micro sign, and drops one-letter tokens), and follows the gold's
# own grading_note. Defined from the gold text, not from any model's answer.
FIELD_RULES = {
    ("GZ-38", "units_used"): (
        lambda m: _norm(m).split(" ")[0] in ("um", "umol/l", "micromolar") or _norm(m).startswith("um "),
        "gold is 'uM (as labeled on the chart)': pass iff the reported unit is micromolar as labelled"),
    ("GZ-40", "which_shows_dip_then_rise"): (
        lambda m: _letter(m) == "B", "gold 'Surfactant B': pass iff the answer names series B"),
    ("GZ-40", "which_shows_plain_plateau"): (
        lambda m: _letter(m) == "A", "gold 'Surfactant A': pass iff the answer names series A"),
    ("GZ-36", "slope_sign"): (
        lambda m: _norm(m).startswith("negative"),
        "gold sign is 'negative' (the rest is explanation); same resolution as groundzero JUDGMENT_RESOLUTIONS"),
    ("GZ-20", "isotherm_model"): (
        lambda m: (re.sub(r"[^a-z]+", " ", _norm(m)).split() or [""])[0] in ("langmuir", "szyszkowski"),
        "gold 'langmuir (BIC-selected over Frumkin, ...)': the choice is the first word, the rest is justification; "
        "the shared token matcher failed a plain 'langmuir' and passed 'frumkin (BIC-selected over Langmuir)'; "
        "'Szyszkowski' (the Langmuir equation of state) counts as the Langmuir choice (author decision 2026-10-21)"),
    ("GZ-15", "estimated_cmc_mM"): (
        lambda m: GZ.is_close(m, GZ_GOLD["GZ-15"]["gold"]["estimated_cmc_mM"], rel=0.10),
        "author decision 2026-10-20: +/-10% (shared comparator formula, tolerance 0.10); the gold note's "
        "'wider tolerance' gave no number"),
}

# (arm, qid, field) -> (passed, rationale): written judgments for free-text fields,
# each made by reading the answer against the gold and its grading_note.
_GZ39_REASON = ("identifies micelle formation at/above the CMC (monomer activity pinned, surface tension "
                "plateaus) as what the pre-CMC isotherm cannot capture -- the gold's explanation")
_GZ39_REGION = ("grading_note credits placing the deviation in the HIGH-concentration region near/above the CMC; "
                "the answer does")
JUDGMENTS: dict[tuple[str, str, str], tuple[bool, str]] = {
    ("claude_unaug", "GZ-39", "deviation_reason"): (True, _GZ39_REASON),
    ("claude_aug", "GZ-39", "deviation_reason"): (True, _GZ39_REASON),
    ("chatgpt_unaug", "GZ-39", "deviation_reason"): (True, _GZ39_REASON),
    ("deepseek_unaug", "GZ-39", "deviation_reason"): (True, _GZ39_REASON),
    ("qwen_unaug", "GZ-39", "deviation_reason"): (True, _GZ39_REASON),
    ("gemini_unaug", "GZ-39", "deviation_reason"): (True, _GZ39_REASON),
    ("gemini_unaug", "GZ-39", "deviation_region"): (
        True, _GZ39_REGION + " ('concentrations above ~2-3 mM', the gold's own range)"),
    ("chatgpt_aug", "GZ-39", "deviation_reason"): (True, _GZ39_REASON),
    ("chatgpt_aug", "GZ-39", "deviation_region"): (
        True, _GZ39_REGION + " ('approximately 2-3 mM and above', the gold's own range)"),
    ("chatgpt_unaug", "GZ-39", "deviation_region"): (
        True, _GZ39_REGION + " ('above 1-2 mM, clear divergence by 3-4 mM and higher')"),
    ("deepseek_unaug", "GZ-39", "deviation_region"): (
        True, _GZ39_REGION + " ('2-3 mM and above', the gold's own range)"),
    ("claude_unaug", "GZ-39", "deviation_region"): (
        True, _GZ39_REGION + " ('above approximately 1.5-2 mM, the high-concentration tail')"),
    ("claude_aug", "GZ-39", "deviation_region"): (
        True, _GZ39_REGION + " ('above roughly 1-2 mM, at and beyond the CMC (~2.5 mM)')"),
    ("qwen_unaug", "GZ-39", "deviation_region"): (
        True, _GZ39_REGION + " ('high-concentration region above ~2 mM, at and above the CMC')"),
}
# Fields where the shared keyword matcher fails on wording a reader would accept are routed to
# JUDGMENTS too (not only > 80-char gold); the comparator itself is left unchanged.
JUDGED_FIELDS = {("GZ-39", "deviation_region")}

# GZ-15's grading_note asks for a "wider tolerance than the other numeric fields" on
# estimated_cmc_mM but gives no number; the author set +/-10% on 2026-10-20 (FIELD_RULES above).
# Disclose as an author-set tolerance. Answers vs gold are printed below for transparency.
GZ15_NOTE = "GZ-15 estimated_cmc_mM graded at +/-10% (author decision); gold note's 'wider tolerance' was unquantified"


def normalize_quotes(text: str) -> str:
    return text.replace("“", '"').replace("”", '"').replace("‘", "'").replace("’", "'")


def grade_gz(arm: str, qid: str, text: str) -> dict:
    ans = GZ.extract_last_json_block(normalize_quotes(text))
    if not isinstance(ans, dict):
        return {"verdict": "NO_JSON", "fields": {}}
    gold = GZ_GOLD[qid]["gold"]
    keys = (GZ_TEMPLATES.get(qid, set()) & set(gold)) - GZ.META_KEYS
    fields, all_ok, pending = {}, True, False
    for k in sorted(keys):
        rule = FIELD_RULES.get((qid, k))
        if rule:
            ok, kind = bool(rule[0](ans.get(k))), "field-rule"
        elif (qid, k) in JUDGED_FIELDS:
            ok, kind = None, "needs-judgment"
        else:
            ok, kind = GZ.grade_field(ans.get(k), gold[k])
        entry = {"ok": ok, "kind": kind, "gold": gold[k], "model": ans.get(k)}
        if rule:
            entry["rationale"] = rule[1]
        if kind == "needs-judgment":
            res = JUDGMENTS.get((arm, qid, k))
            if res is None:
                pending = True
            else:
                entry.update(ok=res[0], kind="judgment", rationale=res[1])
                ok = res[0]
        if ok is False:
            all_ok = False
        fields[k] = entry
    verdict = "FAIL" if not all_ok else ("NEEDS_JUDGMENT" if pending else "PASS")
    return {"verdict": verdict, "fields": fields}


def grade_h(qid: str, text: str) -> dict:
    ans = H50.extract_final_json(normalize_quotes(text))
    if not isinstance(ans, dict):
        return {"verdict": "NO_JSON", "fields": {}}
    fields, all_ok = {}, True
    for k, fd in H50.GOLD[qid]["fields"].items():
        got = ans.get(k, "__MISSING__")
        ok = bool(H50.grade_field(fd, got))
        fields[k] = {"ok": ok, "gold": fd.get("value"), "model": got}
        all_ok = all_ok and ok
    return {"verdict": "PASS" if all_ok else "FAIL", "fields": fields}


def grade_x(qid: str, text: str) -> dict:
    ans = GZ.extract_last_json_block(normalize_quotes(text))     # same balanced-brace FINAL_JSON parse
    if not isinstance(ans, dict):
        return {"verdict": "NO_JSON", "fields": {}}
    r = X50.grade_question(qid, ans)
    return {"verdict": "PASS" if r["pass"] else "FAIL",
            "fields": {k: {"ok": bool(v), "model": ans.get(k)} for k, v in r["fields"].items()}}


def grade_all(answers_dir: Path = ANSWERS, qids: dict = QIDS) -> dict:
    out = {}
    for arm in ARMS:
        d = answers_dir / arm
        if not d.is_dir():
            continue
        out[arm] = {}
        for short, qid in qids.items():
            p = d / f"{short}.txt"
            if not p.exists():
                continue                    # not run (or not yet): absent, not a fail
            text = p.read_text(encoding="utf-8", errors="replace")
            if qid.startswith("GZ"):
                r = grade_gz(arm, qid, text)
            elif qid.startswith("H"):
                r = grade_h(qid, text)
            else:
                r = grade_x(qid, text)
            meta = d / f"{short}.meta.json"
            r["meta"] = json.loads(meta.read_text(encoding="utf-8")) if meta.exists() else None
            out[arm][qid] = r
    return out


def main() -> int:
    sys.stdout.reconfigure(encoding="utf-8")
    dg = "--dg" in sys.argv                    # grade the counterion-correction re-runs instead
    qids = QIDS_DG if dg else QIDS
    res = grade_all(HERE / "answers_dg" if dg else ANSWERS, qids)
    print(f"{'arm':<16}" + "".join(f"{q:>7}" for q in qids.values()))
    pending = []
    for arm, qs in res.items():
        row = f"{arm:<16}"
        for qid in qids.values():
            v = qs.get(qid, {}).get("verdict", "-")
            row += f"{ {'PASS': 'P', 'FAIL': 'F', 'NO_JSON': 'NJ', 'NEEDS_JUDGMENT': 'J?', '-': '.'}[v]:>7}"
            for k, fv in qs.get(qid, {}).get("fields", {}).items():
                if fv.get("kind") == "needs-judgment":
                    pending.append((arm, qid, k, fv["gold"], fv["model"]))
        print(row)
    if pending:
        print(f"\n{len(pending)} field(s) need a written judgment (add to JUDGMENTS):")
        for arm, qid, k, g, m in pending:
            print(f"- ({arm!r}, {qid!r}, {k!r})\n    gold : {g}\n    model: {m}")
    print(f"\nNote: {GZ15_NOTE}. Answers given:")
    for arm, qs in res.items():
        f = qs.get("GZ-15", {}).get("fields", {}).get("estimated_cmc_mM")
        if f:
            dev = abs(float(f["model"]) - float(f["gold"])) / float(f["gold"]) * 100 if f["model"] is not None else None
            print(f"  {arm:<16} {f['model']} mM ({dev:.0f}% from 2.51)" if dev is not None else f"  {arm:<16} none")
    name = "dg_regrade_results.json" if dg else "repair_grade_results.json"
    (HERE / name).write_text(json.dumps(res, indent=1, default=str), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
