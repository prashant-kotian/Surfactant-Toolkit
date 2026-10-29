#!/usr/bin/env python3
"""SB-197 (orig_hard50 Q47) regrade after the dynamic-tension equation-of-state change.

Superseded for Q47 by orig_q16_q17_q47_grade.py (2026-10-28: one wording-consistent target per
question); kept, with its results file, as the record of the earlier regrade.

Q47 asks for gamma(t = 0.1 s) from dynamic_surface_tension_forward. Its reference value was
computed when the tool converted the surface excess with gamma0 - nRT*Gamma; the tool now uses
the Szyszkowski equation of the same Langmuir isotherm, and the reference value is recomputed
here from the release toolkit (no hand-typed number).

Unaugmented arms: graded from the answer stored in the frozen record. Augmented arms reported the
earlier tool output, so each is re-run against the current tool (answers_q47/<arm>/Q47.txt) and
graded on that answer. Grading rule as orig_hard50_grade.py (3 decimal places).

The script writes q47_regrade_results.json only when every augmented arm has a re-run answer;
until then it lists the pending arms and leaves the results untouched.

    python q47_regrade.py
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(HERE))
from orig_hard50_grade import GOLDS, grade_question  # noqa: E402
from surfactantkit.dynamic_tension import dynamic_surface_tension_forward  # noqa: E402

ORIG = ROOT / "benchmark" / "source_sets" / "orig_hard50"
ARM_OF = {("unaug", "Claude"): "claude_unaug", ("unaug", "ChatGPT"): "chatgpt_unaug",
          ("unaug", "DeepSeek"): "deepseek_unaug", ("unaug", "Qwen"): "qwen_unaug",
          ("unaug", "Gemini"): "gemini_unaug", ("aug", "Claude"): "claude_aug",
          ("aug", "ChatGPT"): "chatgpt_aug", ("aug", "DeepSeek"): "deepseek_aug",
          ("aug", "Antigravity"): "agy_aug"}


def reference_value() -> float:
    """Q47 as posed: t_s = [0.1, 0.5, 1, 3, 10, 30, 100], c0 = 0.02 M, D = 1e-9, K = 1500,
    Gamma_inf = 4e-6, defaults otherwise; first element of gamma, 3 decimals."""
    r = dynamic_surface_tension_forward([0.1, 0.5, 1.0, 3.0, 10.0, 30.0, 100.0], 0.02, 1.0e-9,
                                        1500.0, 4.0e-6)
    return round(r.gamma_mN_per_m[0], 3)


def main() -> int:
    gold = reference_value()
    GOLDS["Q47"] = [(gold, 3)]
    out, pending = {}, []
    for cond, fname in (("unaug", "unaugmented_grade_results.json"),
                        ("aug", "augmented_grade_results.json")):
        d = json.loads((ORIG / fname).read_text(encoding="utf-8"))
        for model, md in d.items():
            arm = ARM_OF[(cond, model)]
            rec = md["questions"].get("Q47", {})
            rerun = HERE / "answers_q47" / arm / "Q47.txt"
            if rerun.exists():
                m = re.search(r"ANSWER Q47:\s*(-?[\d.]+)", rerun.read_text(encoding="utf-8"))
                ans, source = ([float(m.group(1))] if m else None), f"answers_q47/{arm}/Q47.txt"
            elif cond == "unaug":
                m = re.search(r"got (-?[\d.]+)", rec.get("detail", ""))
                ans, source = ([float(m.group(1))] if m else None), "stored answer (frozen record)"
            else:
                pending.append(arm)
                continue
            verdict, detail = grade_question("Q47", ans)
            out[arm] = {"Q47": {"verdict": verdict, "detail": detail, "source": source,
                                "previous_verdict": rec.get("verdict"),
                                "fields": {"gamma_t0.1": {"ok": verdict == "PASS", "gold": gold,
                                                          "model": ans[0] if ans else None}}}}
    print(f"reference value (release toolkit): {gold}")
    for arm, r in out.items():
        q = r["Q47"]
        print(f"  {arm:15} {q['previous_verdict']:>5} -> {q['verdict']:<5} {q['source']}")
    if pending:
        print(f"NOT WRITTEN: augmented re-run pending for {pending}")
        return 2
    (HERE / "q47_regrade_results.json").write_text(
        json.dumps({"reference_value": gold, "arms": out}, indent=1) + "\n", encoding="utf-8")
    print("q47_regrade_results.json written")
    return 0


if __name__ == "__main__":
    sys.exit(main())
