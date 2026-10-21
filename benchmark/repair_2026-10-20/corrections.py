"""Corrections layer applied by build_results.py (2026-10-20). Data and rules only; every change is logged.

1. Choice-field regrade (groundzero). Golds of the form '<choice> (<justification>)' were scored by the
   shared token matcher, which failed plain correct answers ('langmuir') and passed opposite ones
   ('frumkin (BIC-selected over Langmuir)'). Rule, identical for every arm: the field passes iff the
   answer's first word is the gold's choice. A record is only touched when its stored verdict is
   reproducible from its own field flags (verdict PASS iff no field is False); otherwise it is left
   unchanged and logged.

2. Chart/typo repair (10 questions excluded until now: GZ-15, GZ-36..40, H-36..38, X-19). Re-run on
   2026-10-20 with the chart attached / corrected X-19 text and graded by grade_repair.py. Arms whose
   original runs already delivered the charts (DeepSeek augmented except GZ-36; Gemini CLI
   augmented) keep their original verdicts. A question stays excluded for an arm only when it has
   no gradable re-run answer (reason logged).

3. Counterion-factor correction (14 questions). The question text of GZ-16, GZ-20, H-02, H-03, H-08,
   H-16, H-19, H-31, H-38, H-47, H-49, X-41, X-47 stated the inverse factor at run time; every arm
   was re-run on the corrected text (grade_repair.py --dg) and those verdicts replace the originals.
   H-46 (gold-only error): augmented arms re-run; unaugmented arms are regraded from their stored
   answers against the corrected gold where the record holds the answer, else re-run.
   No gradable re-run answer -> excluded for that arm (reason logged), never the stale verdict.
"""
from __future__ import annotations

import copy
import importlib.util
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

CHOICE_FIELDS = {("GZ-01", "isotherm_model"): "langmuir", ("GZ-20", "isotherm_model"): "langmuir",
                 ("GZ-29", "selected_model"): "frumkin", ("GZ-23", "valid_hlb_method"): "griffin",
                 ("GZ-35", "valid_hlb_method"): "griffin"}

# The Szyszkowski equation is the Langmuir isotherm written as a surface-tension equation of state, so
# naming it is the Langmuir choice (author decision 2026-10-21; one line to revert).
CHOICE_SYNONYMS = {"langmuir": {"langmuir", "szyszkowski"}}

ARM = {("Claude", "unaug"): "claude_unaug", ("Claude", "aug"): "claude_aug",
       ("ChatGPT", "unaug"): "chatgpt_unaug", ("ChatGPT", "aug"): "chatgpt_aug",
       ("DeepSeek", "unaug"): "deepseek_unaug", ("DeepSeek", "aug"): "deepseek_aug",
       ("Qwen", "unaug"): "qwen_unaug", ("Gemini", "unaug"): "gemini_unaug",
       ("Gemini", "aug"): "agy_aug"}   # Gemini augmented arm = Gemini CLI (Antigravity)

REPAIR = {"groundzero": ["GZ-15", "GZ-36", "GZ-37", "GZ-38", "GZ-39", "GZ-40"],
          "hard50": ["H-36", "H-37", "H-38"], "extrahard": ["X-19"]}
CHARTS_ORIGINALLY_DELIVERED = {
    ("DeepSeek", "aug"): {"GZ-15", "GZ-37", "GZ-38", "GZ-39", "GZ-40", "H-36", "H-37", "H-38"},
    ("Gemini", "aug"): {"GZ-15", "GZ-36", "GZ-37", "GZ-38", "GZ-39", "GZ-40", "H-36", "H-37", "H-38"},
}
NOT_RERUN = {(("Gemini", "aug"), "X-19"): "Gemini-CLI extrahard tier retired; X-19 not re-run",
             (("Gemini", "aug"), "X-41"): "Gemini-CLI extrahard tier retired; not re-run on corrected text",
             (("Gemini", "aug"), "X-47"): "Gemini-CLI extrahard tier retired; not re-run on corrected text"}

DG = {"groundzero": ["GZ-16", "GZ-20"],
      "hard50": ["H-02", "H-03", "H-08", "H-16", "H-19", "H-31", "H-38", "H-46", "H-47", "H-49"],
      "extrahard": ["X-41", "X-47"]}
GOLD_ONLY = {"H-46"}


def _load_grader():
    spec = importlib.util.spec_from_file_location("grade_repair", HERE / "grade_repair.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _results(name: str) -> dict:
    p = HERE / name
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else {}


REPAIR_RES = _results("repair_grade_results.json")
DG_RES = _results("dg_regrade_results.json")
FILL_RES = _results("fill_regrade_results.json")


def _head(v) -> str:
    w = re.sub(r"[^a-z ]", " ", str(v or "").lower()).split()
    return w[0] if w else ""


def choice_regrade(obj: dict) -> tuple[dict, list[dict]]:
    """Groundzero V-format records ({qid: {verdict, fields}}) -> corrected copy and change log."""
    out, log = copy.deepcopy(obj), []
    for (qid, field), choice in CHOICE_FIELDS.items():
        rec = out.get(qid)
        if not isinstance(rec, dict) or not isinstance(rec.get("fields"), dict):
            continue
        f = rec["fields"].get(field)
        if not isinstance(f, dict):
            continue
        verdict = str(rec.get("verdict", "")).upper()
        if verdict not in ("PASS", "FAIL"):
            continue
        flags = [x.get("ok") for x in rec["fields"].values() if isinstance(x, dict)]
        if (verdict == "PASS") != all(x is not False for x in flags):
            log.append(dict(qid=qid, field=field, change="none",
                            reason="stored verdict not reproducible from field flags; left unchanged"))
            continue
        new_ok = _head(f.get("model", f.get("got"))) in CHOICE_SYNONYMS.get(choice, {choice})
        if new_ok == f.get("ok"):
            continue
        old_ok = f.get("ok")
        f.update(ok=new_ok, kind="choice-rule",
                 rationale=f"choice field: passes iff the answer's first word is '{choice}'")
        flags = [x.get("ok") for x in rec["fields"].values() if isinstance(x, dict)]
        rec["verdict"] = "PASS" if all(x is not False for x in flags) else "FAIL"
        log.append(dict(qid=qid, field=field, field_ok=f"{old_ok}->{new_ok}",
                        answer=str(f.get("model", f.get("got")))[:80], verdict=f"{verdict}->{rec['verdict']}"))
    return out, log


def _rerun(res: dict, arm: str | None, qid: str):
    r = (res.get(arm) or {}).get(qid) if arm else None
    v = r.get("verdict") if isinstance(r, dict) else None
    return v if v in ("PASS", "FAIL") else None, v


def _stored_answer_regrade(qid: str, rec) -> bool | None:
    """Regrade one hard50 record from the answer values it stores, against the current gold."""
    if not isinstance(rec, dict) or not isinstance(rec.get("fields"), dict):
        return None
    G = _load_grader()
    gold = G.H50.GOLD[qid]["fields"]
    got = {k: (f.get("got", f.get("model")) if isinstance(f, dict) else None) for k, f in rec["fields"].items()}
    if set(gold) - set(got):
        return None
    return all(bool(G.H50.grade_field(fd, got[k])) for k, fd in gold.items())


def apply(bank: str, model: str, cond: str, verdicts: dict, obj, excluded: list[str]):
    """-> (verdicts, excluded, log). `verdicts` covers every question the source has; `excluded` is
    the standard exclusion list for this cell."""
    key, arm = (model, cond), ARM.get((model, cond))
    v, ex, log = dict(verdicts), list(excluded), []

    def lift(q):
        if q in ex:
            ex.remove(q)

    def exclude(q, why):
        v.pop(q, None)
        if q not in ex:
            ex.append(q)
        log.append(dict(qid=q, change="excluded", reason=why))

    for q in REPAIR.get(bank, []):
        if q in DG.get(bank, []):
            continue                                   # H-38: the counterion re-run supersedes
        if (key, q) in NOT_RERUN:
            log.append(dict(qid=q, change="stays excluded", reason=NOT_RERUN[(key, q)]))
            continue
        if q in CHARTS_ORIGINALLY_DELIVERED.get(key, set()):
            if q in verdicts:
                lift(q)
                log.append(dict(qid=q, change="included (original run delivered the chart)",
                                verdict="PASS" if verdicts[q] else "FAIL"))
            continue
        got, raw = _rerun(REPAIR_RES, arm, q)
        if got:
            lift(q)
            v[q] = got == "PASS"
            log.append(dict(qid=q, change="included (2026-10-20 re-run)", verdict=got))
        else:
            log.append(dict(qid=q, change="stays excluded",
                            reason=f"no gradable re-run answer ({raw or 'not run'})"))

    for q in DG.get(bank, []):
        if (key, q) in NOT_RERUN:
            exclude(q, NOT_RERUN[(key, q)])
            continue
        old = verdicts.get(q)
        if q in GOLD_ONLY and cond == "unaug":
            ok = _stored_answer_regrade(q, obj.get(q) if isinstance(obj, dict) else None)
            if ok is not None:
                lift(q)
                v[q] = ok
                log.append(dict(qid=q, change="regraded stored answer vs corrected gold",
                                verdict=f"{old}->{ok}"))
                continue
        got, raw = _rerun(DG_RES, arm, q)
        if got:
            lift(q)
            v[q] = got == "PASS"
            log.append(dict(qid=q, change="replaced by corrected-text re-run", verdict=f"{old}->{got == 'PASS'}"))
        else:
            exclude(q, f"awaiting gradable re-run on corrected text ({raw or 'not run yet'})")

    # 4. Net-200 fill (2026-10-21): GZ-26 re-run on the current text, answers never captured in the
    #    original runs, Gemini-CLI augmented X-19/X-41/X-47. A gradable fill answer replaces the
    #    exclusion or the missing record for that arm; anything else stays as above.
    for q, r in sorted((FILL_RES.get(arm) or {}).items()) if arm else []:
        if not q.startswith(BANK_PREFIX.get(bank, "?")) or r.get("verdict") not in ("PASS", "FAIL"):
            continue
        old_v = v.get(q)
        lift(q)
        v[q] = r["verdict"] == "PASS"
        log[:] = [e for e in log if not (e.get("qid") == q and e.get("change") in ("excluded", "stays excluded"))]
        log.append(dict(qid=q, change="included (net-200 fill re-run)",
                        verdict=f"{old_v}->{v[q]}" if old_v is not None else r["verdict"]))
    return v, ex, log


BANK_PREFIX = {"groundzero": "GZ-", "hard50": "H-", "extrahard": "X-"}


def orig_overlay(model: str, cond: str) -> dict:
    """orig_hard50 fill verdicts for an arm: {"Q07": "PASS"|"FAIL", ...} (re-run of NOT_ANSWERED items)."""
    arm = ARM.get((model, cond))
    return {q: r["verdict"] for q, r in (FILL_RES.get(arm) or {}).items()
            if q.startswith("Q") and r.get("verdict") in ("PASS", "FAIL")}


def detail(bank: str, model: str, cond: str, qid: str, obj) -> list | None:
    """Per-field (ok, gold, got) for a question whose verdict this layer set (re-run or regraded),
    so failure-anatomy figures describe the verdict actually scored. None = not set by this layer."""
    arm = ARM.get((model, cond))
    if arm and isinstance((FILL_RES.get(arm) or {}).get(qid), dict) and FILL_RES[arm][qid].get("verdict") in ("PASS", "FAIL"):
        return [(bool(f.get("ok")), f.get("gold"), f.get("model")) for f in FILL_RES[arm][qid].get("fields", {}).values()]
    for res, qset in ((DG_RES, DG.get(bank, [])), (REPAIR_RES, REPAIR.get(bank, []))):
        if qid in qset and arm:
            r = (res.get(arm) or {}).get(qid)
            if isinstance(r, dict) and r.get("verdict") in ("PASS", "FAIL"):
                return [(bool(f.get("ok")), f.get("gold"), f.get("model")) for f in r.get("fields", {}).values()]
    if qid in GOLD_ONLY and cond == "unaug" and isinstance(obj, dict) and isinstance(obj.get(qid), dict):
        G = _load_grader()
        rec = obj[qid].get("fields", {})
        out = []
        for k, fd in G.H50.GOLD[qid]["fields"].items():
            got = rec.get(k, {}).get("got", rec.get(k, {}).get("model")) if isinstance(rec.get(k), dict) else None
            out.append((bool(G.H50.grade_field(fd, got)), fd.get("value"), got))
        return out
    return None
