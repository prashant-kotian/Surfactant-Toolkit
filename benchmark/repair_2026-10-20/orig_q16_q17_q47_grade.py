#!/usr/bin/env python3
"""SB-166, SB-167 and SB-197 (orig_hard50 Q16, Q17, Q47) graded on one target per question
(audit 2026-10-28).

- Q16 (Ward-Tordai, Gamma at 1 s): wording unchanged. The reference value is recomputed from the
  release ward_tordai_adsorption, which now integrates the Ward-Tordai equation by product
  integration (2.8568 umol/m^2; the finite-difference solver and the long-time asymptote give the
  same). Unaugmented arms: stored answers regraded. Augmented arms: answers from the current tool
  (answers_2026-10-28/<arm>/Q16.txt).
- Q17 (Davies HLB with the Guo correction): the unaugmented wording now states the Davies group
  number of an EO unit (0.33), so both wordings have the target the tool gives (9.62).
  Unaugmented arms: answers to the new wording (answers_2026-10-28/<arm>/Q17.txt). Augmented
  arms: wording unchanged, answers in answers_q17/<arm>/Q17.txt.
- Q47 (dynamic surface tension at 0.1 s): both wordings now state every setting (gamma0, n, T,
  k_off, the equilibrium constant K with k_a = K k_off) and ask for 2 decimals; the reference value
  is the release dynamic_surface_tension_forward with those inputs. All arms answered the new
  wording (answers_2026-10-28/<arm>/Q47.txt).

Grading rule as orig_hard50_grade.py (one unit of the last requested decimal). The results file
is written only when every arm has a gradable answer; until then the pending ones are listed.

    python orig_q16_q17_q47_grade.py
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
from surfactantkit.adsorption import ward_tordai_adsorption  # noqa: E402
from surfactantkit.dynamic_tension import dynamic_surface_tension_forward  # noqa: E402
from surfactantkit.hlb import hlb_davies_guo_ecl  # noqa: E402

ORIG = ROOT / "benchmark" / "source_sets" / "orig_hard50"
NEW = HERE / "answers_2026-10-28"
ARM_OF = {("unaug", "Claude"): "claude_unaug", ("unaug", "ChatGPT"): "chatgpt_unaug",
          ("unaug", "DeepSeek"): "deepseek_unaug", ("unaug", "Qwen"): "qwen_unaug",
          ("unaug", "Gemini"): "gemini_unaug", ("aug", "Claude"): "claude_aug",
          ("aug", "ChatGPT"): "chatgpt_aug", ("aug", "DeepSeek"): "deepseek_aug",
          ("aug", "Antigravity"): "agy_aug"}


def references() -> dict:
    """Reference values from the release toolkit, rounded as each question asks."""
    q16 = ward_tordai_adsorption(0.02, 1e-9, t_end_s=1.0, mode="langmuir",
                                 gamma_inf_mol_per_m2=3.0e-6, K_d_mol_per_L=1.0e-3)
    q47 = dynamic_surface_tension_forward([0.1, 0.5, 1.0, 3.0, 10.0, 30.0, 100.0], 0.02, 1.0e-9,
                                          1500.0, 4.0e-6, gamma0_mN_per_m=72.0, n_factor=1.0,
                                          temperature_K=298.15, k_off_per_s=0.5)
    return {"Q16": (round(q16.gamma_final_mol_per_m2 * 1e6, 4), 4),
            "Q17": (round(hlb_davies_guo_ecl(12, 6), 2), 2),
            "Q47": (round(q47.gamma_mN_per_m[0], 2), 2)}


def parse(path: Path, qid: str):
    txt = path.read_text(encoding="utf-8")
    m = (re.search(rf"ANSWER {qid}:\s*(-?\d+(?:\.\d+)?)", txt)
         or re.search(r"answer\W{0,3}:\s*(-?\d+(?:\.\d+)?)", txt, re.I))
    return [float(m.group(1))] if m else None


def main() -> int:
    ref = references()
    for q, g in ref.items():
        GOLDS[q] = [g]
    out, pending = {}, []
    for cond, fname in (("unaug", "unaugmented_grade_results.json"),
                        ("aug", "augmented_grade_results.json")):
        d = json.loads((ORIG / fname).read_text(encoding="utf-8"))
        for model, md in d.items():
            arm = ARM_OF[(cond, model)]
            for qid in ("Q16", "Q17", "Q47"):
                rec = md["questions"].get(qid, {})
                if qid == "Q16" and cond == "unaug":
                    m = re.search(r"got (-?[\d.]+)", rec.get("detail", ""))
                    ans, source = ([float(m.group(1))] if m else None), "stored answer (frozen record)"
                elif qid == "Q17" and cond == "aug":
                    p = HERE / "answers_q17" / arm / "Q17.txt"
                    ans, source = (parse(p, qid) if p.exists() else None), f"answers_q17/{arm}/Q17.txt"
                else:
                    p = NEW / arm / f"{qid}.txt"
                    ans, source = (parse(p, qid) if p.exists() else None), f"answers_2026-10-28/{arm}/{qid}.txt"
                if ans is None:
                    pending.append(f"{arm}/{qid}")
                    continue
                verdict, detail = grade_question(qid, ans)
                out.setdefault(arm, {})[qid] = {
                    "verdict": verdict, "detail": detail, "source": source,
                    "previous_verdict": rec.get("verdict"),
                    "fields": {"answer": {"ok": verdict == "PASS", "gold": ref[qid][0], "model": ans[0]}}}
    print("reference values (release toolkit):", {q: g[0] for q, g in ref.items()})
    for arm in sorted(out):
        for qid, r in sorted(out[arm].items()):
            print(f"  {arm:15} {qid} {r['previous_verdict']:>5} -> {r['verdict']:<5} "
                  f"{r['fields']['answer']['model']}  ({r['source']})")
    if pending:
        print(f"NOT WRITTEN: answers pending for {pending}")
        return 2
    (HERE / "orig_q16_q17_q47_results.json").write_text(
        json.dumps({"reference_values": {q: g[0] for q, g in ref.items()}, "arms": out}, indent=1) + "\n",
        encoding="utf-8")
    print("orig_q16_q17_q47_results.json written")
    return 0


if __name__ == "__main__":
    sys.exit(main())
