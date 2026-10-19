"""SurfBench pilot orchestrator.

Runs the selected question subset across chosen models and conditions, extracts
each model's final answer, grades it against the code-computed gold answer, and
writes a timestamped result file to benchmark/results/ (gitignored).

Two model tiers:
  * smoke tier  -- cheap, tool-capable models, thinking OFF: for validating the
                   pipeline end to end without flagship spend.
  * flagship tier -- the real benchmark models with extended thinking ON. Exact
                   model IDs are resolved/confirmed at real-run time; the ones
                   here are placeholders to be pinned before the scored run.

Usage:
  python run_pilot.py --models gemini-flash,gpt-4o-mini,qwen-plus --conditions unaugmented,surfmcp --limit 2
"""

from __future__ import annotations
import argparse
import datetime
import json
import os
import time
from pathlib import Path


def _load_env():
    """Load benchmark/harness/.env (gitignored) into os.environ if present, without
    overriding variables already set in the real environment. Keeps API keys out of
    every command line and out of the (public) repo."""
    envf = Path(__file__).resolve().parent / ".env"
    if not envf.exists():
        return
    for line in envf.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        k, v = line.split("=", 1)
        os.environ.setdefault(k.strip(), v.strip().strip('"').strip("'"))


_load_env()

import providers
import prompts
from extract import extract_final
from grading import grade
from select_pilot import select

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"

MODELS = {
    # --- smoke tier (cheap, thinking off) ---
    "gemini-flash": {"provider": "gemini", "model_id": "gemini-2.5-flash", "kw": {}},
    "gpt-4o-mini":  {"provider": "openai", "model_id": "gpt-4o-mini", "kw": {}},
    "qwen-plus":    {"provider": "qwen",   "model_id": "qwen-plus-2025-09-11", "kw": {}},
    # --- flagship tier (thinking ON) -- IDs confirmed live via each provider's
    # own API on 2026-09-02, not guessed. Verification notes:
    #   claude-opus-4-8   -- fixed id, no separate dated alias exists/needed
    #   gpt-5.6-sol       -- fixed id, exact match to the originally chosen model;
    #                        NOT date-suffixed like gpt-5/-5.1/.../-5.5 siblings,
    #                        appears to be OpenAI's naming convention for this tier
    #   qwen-plus-2025-09-11 -- real dated snapshot exists and was confirmed to
    #                        respond; "qwen3-plus" (the originally recalled name)
    #                        does NOT exist as a model -- was an imprecise mention
    #                        earlier in the project, not a real DashScope model id
    #   gemini-2.5-pro    -- NO dated snapshot alias is exposed on this account
    #                        (checked via the models.list endpoint); Gemini is
    #                        rolling-only here. Document this as a real, disclosed
    #                        reproducibility limitation in the methods section
    #                        rather than pretend a pinned date exists.
    # gemini-2.5-pro deprecated ("no longer available to new users", confirmed
    # via the live API's own error message 2026-09-04) -- gemini-3.1-pro-preview
    # is the real, current replacement, but this account's Gemini key is on the
    # free tier, which has a hard 0 quota for any gemini-3.1-pro variant
    # (confirmed via the live API's QuotaFailure detail, not guessed) -- needs a
    # paid-tier upgrade before this model works, same category of real,
    # account-side blocker as GPT-5.6-sol's credit_balance_exhausted.
    "gemini-pro":   {"provider": "gemini", "model_id": "gemini-3.1-pro-preview", "kw": {"thinking_budget": 8000}},
    "gpt-5.6-sol":  {"provider": "openai", "model_id": "gpt-5.6-sol", "kw": {"reasoning_effort": "high"}},
    # note: same underlying model_id as the smoke-tier "qwen-plus" entry above,
    # just with a real thinking budget for the scored run -- Qwen's actual
    # separate qwen-max model was never selected or tested, deliberately not
    # used here to avoid silently swapping in an untested model under this name
    "qwen-plus-flagship": {"provider": "qwen", "model_id": "qwen-plus-2025-09-11", "kw": {"thinking_budget": 4000}},
    "claude-opus":  {"provider": "anthropic", "model_id": "claude-opus-4-8", "kw": {"thinking_budget": 8000}},
}


def run_one(model_key, question, condition):
    spec = MODELS[model_key]
    system = prompts.build_system(condition)
    user = prompts.build_user(question)
    t0 = time.time()
    try:
        out = providers.run(spec["provider"], spec["model_id"], system, user, condition, **spec["kw"])
    except Exception as e:
        out = {"final_text": "", "thinking": "", "tool_calls": [], "raw": {},
               "error": f"{type(e).__name__}: {e}"}
    dt = time.time() - t0

    ans = extract_final(out.get("final_text", ""))
    g = grade(question["gold_answer"], question.get("tolerance", {}),
              question.get("grading_method", "numeric_tolerance"), ans, out.get("final_text", ""),
              question.get("trap_type", "none"))
    return {
        "id": question["id"], "category": question["category"],
        "subcategory": question["subcategory"], "trap_type": question.get("trap_type", "none"),
        "difficulty": question["difficulty"], "tools_required": question["tools_required"],
        "model": model_key, "condition": condition,
        "correct": g["correct"], "grade_reason": g["reason"],
        "extracted_answer": ans,
        "n_tool_calls": len(out.get("tool_calls", [])),
        "tools_used": [tc["name"] for tc in out.get("tool_calls", [])],
        "thinking_chars": len(out.get("thinking", "") or ""),
        "final_text": (out.get("final_text", "") or "")[:1200],
        "seconds": round(dt, 1), "error": out.get("error"),
        "raw": out.get("raw", {}),
    }


def summarize(rows):
    from collections import defaultdict
    acc = defaultdict(lambda: [0, 0])          # (model,cond) -> [correct,total]
    trap = defaultdict(lambda: [0, 0])         # (cond,trap) -> [correct,total]
    errs = 0
    for r in rows:
        if r["error"]:
            errs += 1
        acc[(r["model"], r["condition"])][0] += int(r["correct"])
        acc[(r["model"], r["condition"])][1] += 1
        trap[(r["condition"], r["trap_type"])][0] += int(r["correct"])
        trap[(r["condition"], r["trap_type"])][1] += 1

    print("\n=== accuracy by model x condition ===")
    for (m, c), (ok, tot) in sorted(acc.items()):
        print(f"  {m:14s} {c:12s}  {ok}/{tot}  ({100*ok/tot:.0f}%)")
    print("\n=== accuracy by condition x trap type ===")
    for (c, tt), (ok, tot) in sorted(trap.items()):
        print(f"  {c:12s} {tt:18s}  {ok}/{tot}")
    print(f"\nAPI errors: {errs}/{len(rows)}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", default="gemini-flash,gpt-4o-mini,qwen-plus")
    ap.add_argument("--conditions", default="unaugmented,surfmcp")
    ap.add_argument("--limit", type=int, default=2, help="number of questions from the subset")
    ap.add_argument("--n", type=int, default=20, help="pilot subset size to draw from")
    args = ap.parse_args()

    models = [m.strip() for m in args.models.split(",") if m.strip()]
    conditions = [c.strip() for c in args.conditions.split(",") if c.strip()]
    for m in models:
        if m not in MODELS:
            raise SystemExit(f"unknown model '{m}'. known: {list(MODELS)}")

    subset = select(n=args.n)[: args.limit]
    print(f"running {len(subset)} questions x {len(models)} models x {len(conditions)} conditions "
          f"= {len(subset)*len(models)*len(conditions)} API calls")

    rows = []
    for q in subset:
        for m in models:
            for c in conditions:
                r = run_one(m, q, c)
                flag = "OK " if r["correct"] else "XX "
                if r["error"]:
                    flag = "ERR"
                print(f"  [{flag}] {q['id']:7s} {m:13s} {c:12s} "
                      f"tools={r['n_tool_calls']} think={r['thinking_chars']}c "
                      f"{r['seconds']}s {'| '+r['error'] if r['error'] else ''}")
                rows.append(r)

    RESULTS_DIR.mkdir(exist_ok=True)
    stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    outfile = RESULTS_DIR / f"pilot_{stamp}.json"
    outfile.write_text(json.dumps({"models": models, "conditions": conditions,
                                   "subset_ids": [q["id"] for q in subset], "rows": rows}, indent=2))
    summarize(rows)
    print(f"\nwrote {outfile}")


if __name__ == "__main__":
    main()
