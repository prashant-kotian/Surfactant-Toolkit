#!/usr/bin/env python3
"""SB-167 (orig_hard50 Q17) regrade after the Guo effective-chain-length EO weight change.

Superseded for Q17 by orig_q16_q17_q47_grade.py (2026-10-28: one wording-consistent target per
question); kept, with its results file, as the record of the earlier regrade.

The two wordings of Q17 define their answers differently:

- Augmented wording: "MCP tool: hlb_davies_guo_ecl. Inputs: n_carbons_alkyl = 12, n_eo = 6",
  report the tool's hlb. Its reference value is what the release tool returns; the tool now
  weights the effective EO count with the Davies -CH2CH2O- group number 0.33 (it used the ether-O
  number 1.3 before), so the reference value is recomputed here from the release toolkit and each
  augmented arm is re-run against the current tool (answers_q17/<arm>/Q17.txt).
- Unaugmented wording: states every input of the calculation itself (base 7, ether-oxygen weight
  1.3, CH2/CH3 weight 0.475, both Guo transforms). Its reference value is the result of that
  stated calculation, recomputed below and checked against the frozen reference (33.29). It does
  not depend on the toolkit, so the unaugmented verdicts are unchanged and not overridden.

The script writes q17_regrade_results.json only when every augmented arm has a re-run answer;
until then it lists the pending arms and leaves the results untouched.

    python q17_regrade.py
"""
from __future__ import annotations

import json
import math
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(HERE))
from orig_hard50_grade import GOLDS, grade_question  # noqa: E402
from surfactantkit.hlb import hlb_davies_guo_ecl  # noqa: E402

ORIG = ROOT / "benchmark" / "source_sets" / "orig_hard50"
AUG_ARMS = {"Claude": "claude_aug", "ChatGPT": "chatgpt_aug", "DeepSeek": "deepseek_aug",
            "Antigravity": "agy_aug"}


def reference_value() -> float:
    """Augmented Q17 as posed: hlb_davies_guo_ecl(12, 6), field hlb, 2 decimals."""
    return round(hlb_davies_guo_ecl(12, 6), 2)


def unaugmented_reference() -> float:
    """Unaugmented Q17 as worded: 7 + 1.3 * EO_eff - 0.475 * (CH2_eff + 1), with
    CH2_eff = 0.965 * 11 - 0.178 and EO_eff = 13.45 ln 6 - 0.16 * 6 + 1.26."""
    ch2 = 0.965 * 11 - 0.178
    eo = 13.45 * math.log(6) - 0.16 * 6 + 1.26
    return round(7 + 1.3 * eo - 0.475 * (ch2 + 1), 2)


def main() -> int:
    unaug = unaugmented_reference()
    if abs(unaug - GOLDS["Q17"][0][0]) > 1e-9:
        raise SystemExit(f"unaugmented wording gives {unaug}, frozen reference {GOLDS['Q17'][0][0]}")
    gold = reference_value()
    GOLDS["Q17"] = [(gold, 2)]
    d = json.loads((ORIG / "augmented_grade_results.json").read_text(encoding="utf-8"))
    out, pending = {}, []
    for model, md in d.items():
        arm = AUG_ARMS[model]
        rec = md["questions"].get("Q17", {})
        rerun = HERE / "answers_q17" / arm / "Q17.txt"
        if not rerun.exists():
            pending.append(arm)
            continue
        m = re.search(r"ANSWER Q17:\s*(-?[\d.]+)", rerun.read_text(encoding="utf-8"))
        ans = [float(m.group(1))] if m else None
        verdict, detail = grade_question("Q17", ans)
        out[arm] = {"Q17": {"verdict": verdict, "detail": detail, "source": f"answers_q17/{arm}/Q17.txt",
                            "previous_verdict": rec.get("verdict"),
                            "fields": {"hlb": {"ok": verdict == "PASS", "gold": gold,
                                               "model": ans[0] if ans else None}}}}
    print(f"reference value, augmented wording (release toolkit): {gold}")
    print(f"reference value, unaugmented wording (stated calculation): {unaug} (unchanged)")
    for arm, r in out.items():
        q = r["Q17"]
        print(f"  {arm:15} {q['previous_verdict']:>5} -> {q['verdict']:<5} {q['source']}")
    if pending:
        print(f"NOT WRITTEN: augmented re-run pending for {pending}")
        return 2
    (HERE / "q17_regrade_results.json").write_text(
        json.dumps({"reference_value": gold, "unaugmented_reference_value": unaug, "arms": out},
                   indent=1) + "\n", encoding="utf-8")
    print("q17_regrade_results.json written")
    return 0


if __name__ == "__main__":
    sys.exit(main())
