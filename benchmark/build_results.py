#!/usr/bin/env python3
"""build_results.py - Paper 3, Phase 0B: single machine-graded results master.

BUILD-TIME ONLY (never shipped to the model). Assembles every benchmark cell
from committed, re-runnable grader outputs; hand-typed report numbers are NOT
trusted - a cell with no machine-readable source in the repo becomes a `gap`.

Contract
--------
* CELLS registry: bank x model x condition -> the committed grading JSON and
  how to read it. Every such JSON is produced by a grader in this repo that
  runs against the corrected golds (GZ-48 / H-04 / H-26 fixes 2026-10-18).
* EXCLUSIONS are DATA, applied to every cell of every arm: images the harness
  cannot deliver (GZ-15,36-40; H-36,37,38), retired-ambiguous questions
  (GZ-26), and question-authoring windows (X-19). X-19 is excluded from all
  arms, including the chat-model augmented cells run against the fixed packet,
  so every unaugmented/augmented pair is scored on the same questions.
* `gold_state` is verified, not inferred: for the questions whose gold was
  corrected 2026-10-18 (GZ-48, H-04, H-26), every gold value stored in the
  cell's grading records must equal the current gold ("post-fix"); a mismatch
  is "PRE-FIX gold"; sources without per-field gold are "not field-verifiable".

Usage:  python3 benchmark/build_results.py
Writes benchmark/../results/results.json and ../results/results_summary.md.
"""
from __future__ import annotations

import json
import subprocess
from datetime import datetime, timezone
from hashlib import sha1
from pathlib import Path

HERE = Path(__file__).resolve().parent            # <repo>/benchmark
REPO = HERE.parent
B = HERE

GOLDS = {"groundzero": B / "source_sets/groundzero/gold_answers.json",
         "hard50": B / "source_sets/hard50/gold_answers.json",
         "extrahard": B / "source_sets/extrahard/gold_answers.json"}
GOLD_FIX_DATE = "2026-10-18"

# ---------------------------------------------------------------- registry --
# fmt: V = {qid:{verdict,fields}} | P = {passed,failed[,n_pass]} |
#      B = {qid:bool} | A = {qid:{pass,fields:{name:bool}}}
# Gemini: unaugmented = Gemini web (multi-model unaug cells); augmented = Gemini CLI (Antigravity, agy-cli
# harness). The Gemini CLI unaugmented runs were set aside 2026-10-21 (agy keeps coding/web-search
# tools, so they were not a true unaugmented baseline) and are not part of the benchmark results.
CELLS: list[tuple[str, dict]] = [
    # ---- GroundZero 50 -----------------------------------------------------
    ("groundzero", dict(model="Claude", cond="aug", src=B / "source_sets/groundzero/claude_augmented_grade_results.json",
                  fmt="V", harness="mcp")),
    ("groundzero", dict(model="ChatGPT", cond="aug", src=B / "source_sets/groundzero/chatgpt_augmented_grade_results.json",
                  fmt="V", harness="mcp")),
    ("groundzero", dict(model="Qwen", cond="unaug", src=B / "source_sets/groundzero/qwen_unaugmented_grade_results.json",
                  fmt="V", harness="mcp")),
    ("groundzero", dict(model="DeepSeek", cond="unaug", src=B / "source_sets/groundzero/deepseek_unaugmented_grade_results.json",
                  fmt="V", harness="mcp")),
    ("groundzero", dict(model="*", cond="unaug", src=B / "source_sets/groundzero/unaugmented_2026-09-21/grade_results.json",
                  fmt="V", harness="mcp", multi=True,
                  only=["chatgpt", "gemini", "claude"], note_for=["claude"],
                  gap_note="Claude: full GZ-01..50 transcript (located 2026-10-19; supersedes "
                           "a genuinely truncated GZ-01..10-only Sept-18 capture used through "
                           "2026-10-18). Full coverage, n=43 after standard exclusions, same as "
                           "ChatGPT/Gemini; long-prose fields resolved via the grader's "
                           "JUDGMENT_RESOLUTIONS table (one rationale per field)")),
    ("groundzero", dict(model="Gemini", cond="aug", src=B / "source_sets/groundzero/antigravity_augmented_grade_results.json",
                  fmt="V", harness="agy-cli")),
    ("groundzero", dict(model="DeepSeek", cond="aug", src=B / "source_sets/groundzero/deepseek_augmented_grade_results.json",
                  fmt="V", harness="mcp", si_only=True)),
    # ---- Hard50 (new) ------------------------------------------------------
    ("hard50", dict(model="Claude", cond="aug", src=B / "source_sets/hard50/claude_augmented_grade_results.json",
                    fmt="V", harness="mcp")),
    ("hard50", dict(model="ChatGPT", cond="aug", src=B / "source_sets/hard50/chatgpt_augmented_grade_results.json",
                    fmt="V", harness="mcp")),
    ("hard50", dict(model="Claude", cond="unaug", src=B / "source_sets/hard50/claude_unaugmented_grade_results.json", pick="results",
                    fmt="OK", harness="mcp")),
    ("hard50", dict(model="*", cond="unaug", src=B / "source_sets/hard50/unaugmented_grade_stats.json", pick="baseline",
                    fmt="P", harness="cli", multi=True, count_only=True)),
    ("hard50", dict(model="Gemini", cond="aug", src=B / "source_sets/hard50/antigravity_augmented_grade_results.json",
                    fmt="V", harness="agy-cli")),
    ("hard50", dict(model="DeepSeek", cond="aug", src=B / "source_sets/hard50/deepseek_augmented_grade_results.json",
                    fmt="OK", harness="mcp", si_only=True)),
    # ---- Extra50 ------------------------------------------------------------
    ("extrahard", dict(model="*", cond="unaug", src=REPO / "benchmark/transcripts/extrahard50_unaugmented/grade_results.json",
                 fmt="P", harness="mcp", multi=True)),
    ("extrahard", dict(model="Claude", cond="aug", src=B / "source_sets/extrahard/augmented_runs_2026-10-12/grade_results_augmented.json", pick="Claude",
                 fmt="P", harness="mcp")),
    ("extrahard", dict(model="DeepSeek", cond="aug", src=B / "source_sets/extrahard/augmented_runs_2026-10-12/grade_results_augmented.json", pick="DeepSeek",
                 fmt="P", harness="mcp")),
    ("extrahard", dict(model="ChatGPT", cond="aug", src=B / "source_sets/extrahard/augmented_runs_2026-10-12/grade_results_chatgpt.json",
                 fmt="P", harness="mcp")),
    ("extrahard", dict(model="Gemini", cond="aug", src=B / "source_sets/extrahard/antigravity_augmented_grade_results.json",
                 fmt="A", harness="agy-cli")),
]

EXCLUSIONS: dict[str, list[tuple[list[str], list[str], str]]] = {
    "groundzero": [
        (["*"], ["GZ-15", "GZ-36", "GZ-37", "GZ-38", "GZ-39", "GZ-40"],
         "question requires image/chart input the harnesses cannot deliver"),
        (["*"], ["GZ-26"], "ambiguous question wording; retired 2026-10-18"),
    ],
    "hard50": [
        (["*"], ["H-36", "H-37", "H-38"],
         "image/figure-based questions excluded from all cells"),
    ],
    "extrahard": [
        (["*"], ["X-19"],
         "question body carried a kappa*a typo (1.45 vs 14.5) at run time for "
         "the pre-fix packets (fixed in-bank 2026-10-16); excluded from every "
         "arm so unaugmented/augmented pairs share one question set"),
    ],
}

def _load_corrections():
    import importlib.util
    spec = importlib.util.spec_from_file_location("corrections", B / "repair_2026-10-20" / "corrections.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# 2026-10-20 corrections (choice-field regrade, chart/typo repair re-runs, counterion-factor re-runs):
# data and rules live in repair_2026-10-20/corrections.py; every change is logged in the cell.
CORR = _load_corrections()

NON_Q = {"n_questions", "n_graded", "n_pass", "n_no_json", "pass_pct",
         "n_pass_5models", "summary", "suspects", "suspects_5models",
         "arith_fail_pct", "reason_fail_pct", "arith_fields", "reason_fields",
         "qcause", "missing_parse"}


def load(p: Path):
    return json.loads(p.read_text(encoding="utf-8"))


def sha(p: Path) -> str:
    """Content hash. JSON is hashed in canonical form (sorted keys, compact
    separators) so a grader that rewrites the same data with different
    whitespace or key order leaves the ledger unchanged."""
    if p.suffix == ".json":
        canon = json.dumps(load(p), sort_keys=True, separators=(",", ":"), ensure_ascii=False)
        return sha1(canon.encode("utf-8")).hexdigest()[:12]
    return sha1(p.read_bytes()).hexdigest()[:12]


def gold_ids(bank: str) -> list[str]:
    g = load(GOLDS[bank])
    items = g["questions"] if isinstance(g, dict) and "questions" in g else g
    return [q for q in items if q not in NON_Q and not str(q).startswith("_")]


CORRECTED = {"groundzero": ["GZ-48", "GZ-17"], "hard50": ["H-04", "H-26"]}   # 2026-10-18 gold fixes; GZ-17 2026-10-23


def current_gold(bank: str, q: str) -> dict:
    g = load(GOLDS[bank])
    e = (g["questions"] if isinstance(g, dict) and "questions" in g else g)[q]
    if "gold" in e:                                     # groundzero layout
        return {k: v for k, v in e["gold"].items() if not isinstance(v, str)}
    return {k: f.get("value") for k, f in e.get("fields", {}).items()}   # hard50 layout


def same(a, b) -> bool:
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return abs(a - b) <= 1e-6 * max(1.0, abs(b))
    return a == b


def verify_gold(bank: str, obj) -> str:
    """post-fix / PRE-FIX gold / not field-verifiable, from the gold values the
    grader stored in this cell's records for the corrected questions."""
    checked = False
    for q in CORRECTED.get(bank, []):
        rec = obj.get(q) if isinstance(obj, dict) else None
        if not isinstance(rec, dict) or rec.get("verdict") == "NO_JSON":   # no answer: no stored gold to check
            continue
        flds = {k: f for k, f in rec.get("fields", {}).items()
                if isinstance(f, dict) and "gold" in f}
        if not flds:
            return "not field-verifiable"
        cur = current_gold(bank, q)
        for k, f in flds.items():
            if k in cur and not same(f["gold"], cur[k]):
                return "PRE-FIX gold"
        checked = True
    return "post-fix" if checked or bank not in CORRECTED else "not field-verifiable"


def normalize_cell(obj, fmt: str, ids: list[str]):
    """-> (verdicts {qid: bool}, count_only: bool). Unknown qids dropped."""
    if fmt == "V":
        # NO_JSON-with-no-records = nothing was gradable (transcript section
        # absent OR capture lost the answer). Per Phase-0 doctrine this is an
        # UNGRADEABLE record (like the ChatGPT "unparseable" literals), not a
        # model fail -> excluded from n; the cell then carries a partial flag.
        v = {q: (str(rec.get("verdict", "")).upper() == "PASS")
             for q, rec in obj.items()
             if isinstance(rec, dict) and "verdict" in rec
             and not (str(rec.get("verdict", "")).upper() == "NO_JSON"
                     and not rec.get("fields"))}
        return v, False
    if fmt == "OK":   # legacy: {qid:{fields:{name:{ok,...}}}} ; verdict = all fields ok
        v = {}
        for q, rec in obj.items():
            if isinstance(rec, dict) and "fields" in rec:
                flds = [f for f in rec["fields"].values() if isinstance(f, dict)]
                v[q] = bool(flds) and all(f.get("ok") for f in flds) and rec.get("found_block", True)
        return v, False
    if fmt == "P":
        passed = {str(q).strip() for q in obj.get("passed", [])}
        if "failed" in obj:
            v = {q: (q in passed) for q in ids}
            return v, True
        v = {q: (q in passed) for q in obj}
        return v, False
    if fmt == "B":
        return {q: bool(x) for q, x in obj.items() if q not in NON_Q}, False
    if fmt == "A":
        return {q: bool(rec.get("pass")) for q, rec in obj.items()
                if isinstance(rec, dict) and "pass" in rec}, False
    raise ValueError(fmt)


def main() -> int:
    now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%MZ")
    bank_ids = {b: gold_ids(b) for b in GOLDS}
    import os
    cands = [os.environ.get("TOOLKIT_DIR"), (REPO if (REPO / "src/surfactantkit").exists() else None), REPO.parent / "surfactantkit",
             REPO.parent / "SurfactantKit", REPO.parent / "surfactantkit-computation"]
    tk_sha = "unknown"
    for c0 in cands:
        if not c0:
            continue
        try:
            tk = subprocess.run(["git", "-C", str(c0), "rev-parse", "--short=12", "HEAD"],
                                capture_output=True, text=True, check=True).stdout.strip()
            dirty = bool(subprocess.run(["git", "-C", str(c0), "status", "--porcelain"],
                                        capture_output=True, text=True).stdout.strip())
            tk_sha = tk + ("-dirty" if dirty else "")
            break
        except Exception:
            continue

    cells, gaps = {}, []
    for bank, spec in CELLS:
        base = dict(bank=bank, cond=spec["cond"], harness=spec["harness"],
                    src=spec["src"].relative_to(REPO).as_posix())
        if not spec["src"].exists():
            gaps.append(dict(**base, reason="grading JSON not present in repo (report-only number)"))
            continue
        data = load(spec["src"])
        if spec.get("pick"):
            data = data[spec["pick"]]
        groups = {}
        if spec.get("multi"):
            for nm, obj in data.items():
                if nm in NON_Q or not isinstance(obj, dict):
                    continue
                if spec.get("only") and nm not in spec["only"]:
                    continue
                if spec.get("model") not in ("*", nm) and spec.get("model") != nm:
                    norm = {"chatgpt": "ChatGPT", "claude": "Claude", "deepseek": "DeepSeek",
                            "gemini": "Gemini", "qwen": "Qwen"}.get(nm)
                    if norm != spec.get("model"):
                        continue
                groups[nm] = obj
        else:
            groups[spec["model"]] = data

        for nm, obj in groups.items():
            model = {"chatgpt": "ChatGPT", "claude": "Claude", "deepseek": "DeepSeek",
                     "gemini": "Gemini", "qwen": "Qwen"}.get(nm, nm)
            ex = [q for pats, qids, _ in EXCLUSIONS.get(bank, [])
                  if {spec["cond"], f"{spec['cond']}:{model}", "*"} & set(pats)
                  for q in qids]
            ex = list(dict.fromkeys(ex))
            ex_reason = next((r for pats, qids, r in EXCLUSIONS.get(bank, [])
                              if {spec["cond"], f"{spec['cond']}:{model}", "*"} & set(pats)), None)
            corrections = []
            if bank == "groundzero" and spec["fmt"] == "V":
                obj, corrections = CORR.choice_regrade(obj)
            verdicts, count_only = normalize_cell(obj, spec["fmt"], bank_ids[bank])
            verdicts, ex, log = CORR.apply(bank, model, spec["cond"], verdicts, obj, ex)
            corrections += log
            ids = [q for q in bank_ids[bank] if q not in ex]
            graded = [q for q in ids if q in verdicts]
            if not graded:
                gaps.append(dict(**base, model=model, reason="no gradable records in source"))
                continue
            passed = [q for q in graded if verdicts[q]]
            cells[f"{bank}|{model}|{spec['cond']}"] = dict(
                n=len(graded), passed=len(passed),
                pct=round(100 * len(passed) / len(graded), 1),
                excluded=ex, exclude_reason=ex_reason if ex else None,
                per_q={q: bool(verdicts[q]) for q in graded},
                harness=spec["harness"], source=base["src"],
                source_sha=sha(spec["src"]),
                partial=(len(graded) < 50 - len(ex)
                         and spec["fmt"] == "V"),
                gold_state=verify_gold(bank, obj),
                count_only=count_only, si_only=bool(spec.get("si_only")),
                corrections=corrections,
                note=(spec.get("gap_note")
                      if nm in spec.get("note_for", [nm]) else None),
            )

    # documented gaps (no machine source)
    # GZ50 Claude unaug: RESOLVED 2026-10-19 -- a complete GZ-01..50 transcript
    # (benchmark/source_sets/groundzero/claude_unaugmented/) supersedes the
    # truncated GZ-01..10 Sept-18 capture used through 2026-10-18; a full cell
    # via the CELLS entry above, not a gap.

    gaps.append(dict(bank="groundzero", model="Gemini", cond="unaug",
                     reason="repo transcripts for the Sept-21 Gemini direct run: 11 no-content "
                            "records (ungradeable, excluded); long-prose fields resolved via "
                            "the grader's JUDGMENT_RESOLUTIONS table (none swung the total, "
                            "n=33). The 200Q report's 22/50 came from a fuller backend grading "
                            "pass; transcript confirmed absent on benchmark machine "
                            "(RUNS_REPORT 2026-10-18 search)."))

    out = dict(
        generated=now, contract="machine-graded sources only; report numbers land in gaps",
        toolkit_state=tk_sha,
        gold_sha={b: sha(GOLDS[b]) for b in GOLDS},
        gold_fix_date=GOLD_FIX_DATE,
        cells=cells,
        gaps=gaps,
    )
    (HERE / "../results/results.json").write_text(json.dumps(out, indent=1), encoding="utf-8")

    md = ["# Grading ledger, per source set (generated - do not hand-edit; rerun benchmark/build_results.py). Headline table: results/TABLE200.md", ""]
    for cond in ("unaug", "aug"):
        md += [f"## {cond}", "",
               "| bank | model | passed/n | % | gold | count-only | source |",
               "|---|---|---|---|---|---|---|"]
        for k, c in sorted(cells.items()):
            if not k.endswith("|" + cond):
                continue
            bank, model, _ = k.split("|")
            md.append(f"| {bank} | {model} | {c['passed']}/{c['n']} | {c['pct']} | "
                      f"{c['gold_state']} | {'yes' if c['count_only'] else ''} | {c['source']} |")
        md.append("")
    md += ["## gaps", ""]
    for g in gaps:
        md.append(f"- {g.get('bank')} {g.get('model', '')} {g.get('cond', '')}: {g['reason']}")
    (HERE / "../results/results_summary.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(f"cells={len(cells)} gaps={len(gaps)} toolkit={tk_sha}")
    for k, c in sorted(cells.items()):
        flag = " [SI-only]" if c["si_only"] else ""
        flag += " [PRE-FIX]" if c["gold_state"] == "PRE-FIX gold" else ""
        flag += " [count-only]" if c["count_only"] else ""
        print(f"  {k:42s} {c['passed']:3d}/{c['n']:3d} {c['pct']:5.1f}%{flag}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
