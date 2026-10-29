# -*- coding: utf-8 -*-
"""
Grades the 5 LLMs' unaugmented HARD-60 answers against the real toolkit-
executed gold answers. See module docstring details in the two sibling
files (GOLDS here, ANSWERS in answers_data.py).
"""
from __future__ import annotations
import sys
from answers_data import MODELS

LENIENT_QS = {"Q19", "Q35", "Q41", "Q42", "Q44", "Q46", "Q52", "Q60"}
ADVERSARIAL_QS = {"Q53", "Q56"}

# gold entry: list of numeric (value, dp) tuples and/or plain strings (for
# category/rule answers compared loosely, case-insensitive substring)
GOLDS = {
    "Q01": [(2.6829, 4)],
    "Q02": [(0.8799, 4)],
    "Q03": [(-3.69997, 5)],
    "Q04": [(0.9275, 4)],
    "Q05": [(1.3019, 4), (3.6339, 4)],
    "Q06": [(1.7805, 4)],
    "Q07": [(0.7554, 4)],
    "Q08": [(-26.976, 3)],
    "Q09": [(0.6881, 4)],
    "Q10": [(0.4634, 4), (0.2307, 4)],
    "Q11": [(9.9403, 4), (0.1595, 4)],
    "Q12": [(52.8761, 4)],
    "Q13": [(1.1497, 4)],
    "Q14": [(0.9074, 4)],
    "Q15": [(0.0495, 4)],
    "Q16": [(2.8566, 4)],
    "Q17": [(33.29, 2)],
    "Q18": [(0.1568, 4), (1.5473, 4)],
    "Q19": [(-5.7, 1)],
    "Q20": [(298.61, 2)],
    "Q21": [(1.2456, 4)],
    "Q22": [(0.5006, 4), "vesicle/bilayer"],
    "Q23": [(1.383, 3)],
    "Q24": [(-2.0394, 4)],
    "Q25": [(20.219, 3)],
    "Q26": [(0.0244, 4)],
    "Q27": [(2.3307, 4)],
    "Q28": [(40.386, 3)],
    "Q29": [(2.6086, 4)],
    "Q30": [(-0.02066, 5)],
    "Q31": [(-23.148, 3)],
    "Q32": [(-23.332, 3)],
    "Q33": [(0.0836, 4)],
    "Q34": [(-23.55, 2)],
    "Q35": [(-5.284, 3)],
    "Q36": [(1.510, 3), (37.558, 3)],
    "Q37": [(0.601, 3), (0.4314, 4)],
    "Q38": ["frumkin", (6.75, 2)],
    "Q39": [(0.6410, 4), (0.861, 3)],
    "Q40": [(132.695, 3)],
    "Q41": [(3.8137, 4)],
    "Q42": [(3.512, 3), (3.239, 3), (3.821, 3)],
    "Q43": [(8.630, 3)],
    "Q44": [(0.0388, 4)],
    "Q45": [(73.732, 3)],
    "Q46": [(35.72, 2)],
    "Q47": [(64.441, 3)],
    "Q48": [(39.8200, 4), (3.5713, 4)],
    "Q49": [(41.2292, 4), (43.3510, 4)],
    "Q50": [(87.2662, 4)],
    "Q51": [(19.300, 3)],
    "Q52": [(0.2701, 4)],
    "Q53": ["False", (0.987, 3)],
    "Q54": [(1.2034, 4)],
    "Q55": [(23627.76, 2)],
    "Q56": ["False", (0.2731, 4)],
    "Q57": [(0.2728, 4), (0.0076, 4)],
    "Q58": [(-0.3019, 4), (0.5157, 4)],
    "Q59": [(8.0, 1)],
    "Q60": [(284.30, 2), (0.6535, 4)],
}

QORDER = [f"Q{i:02d}" for i in range(1, 61)]


def num_close(val, gold, dp, slack=1.0):
    tol = slack * (10 ** (-dp))
    return abs(val - gold) <= tol + 1e-9


def grade_question(qid, ans):
    """Returns (status, detail) where status in
    {PASS, FAIL, PARTIAL_ADVERSARIAL_CORRECT, NOT_ANSWERED, LENIENT_INFO}"""
    gold = GOLDS[qid]
    if ans is None:
        return "NOT_ANSWERED", ""

    if qid in ADVERSARIAL_QS:
        bool_gold = gold[0]
        num_gold, dp = gold[1]
        bool_ans = str(ans[0])
        bool_ok = bool_ans.strip().lower() == bool_gold.lower()
        num_ans = ans[1]
        num_ok = num_close(num_ans, num_gold, dp)
        if bool_ok and num_ok:
            return "PASS", ""
        if bool_ok and not num_ok:
            return "PARTIAL_ADVERSARIAL_CORRECT", f"bool correct (gold={bool_gold}); number off: got {num_ans}, gold {num_gold}"
        return "FAIL", f"bool WRONG: got {bool_ans}, gold {bool_gold}"

    if qid in LENIENT_QS:
        # informational only: report numeric deltas, not pass/fail
        parts = []
        ai = 0
        for g in gold:
            if isinstance(g, tuple):
                gv, dp = g
                if ai < len(ans) and isinstance(ans[ai], (int, float)):
                    parts.append(f"got {ans[ai]} vs gold {gv} (delta {ans[ai]-gv:+.4g})")
                else:
                    parts.append(f"got {ans[ai] if ai < len(ans) else '?'} vs gold {gv}")
                ai += 1
        return "LENIENT_INFO", "; ".join(parts)

    # strict grading, multi-part all-or-nothing
    ai = 0
    all_ok = True
    details = []
    for g in gold:
        if isinstance(g, tuple):
            gv, dp = g
            if ai >= len(ans) or not isinstance(ans[ai], (int, float)):
                all_ok = False
                details.append(f"part{ai+1} missing/non-numeric")
            else:
                ok = num_close(ans[ai], gv, dp)
                if not ok:
                    all_ok = False
                    details.append(f"part{ai+1} got {ans[ai]} vs gold {gv}")
            ai += 1
        else:
            # string/category part
            if ai >= len(ans):
                all_ok = False
                details.append(f"part{ai+1} missing (expected '{g}')")
            else:
                got = str(ans[ai]).strip().lower()
                ok = g.lower() in got or got in g.lower()
                if not ok:
                    all_ok = False
                    details.append(f"part{ai+1} got '{ans[ai]}' vs gold '{g}'")
            ai += 1
    return ("PASS" if all_ok else "FAIL"), "; ".join(details)


def main():
    results = {name: {} for name in MODELS}
    for name, answers in MODELS.items():
        for qid in QORDER:
            status, detail = grade_question(qid, answers.get(qid))
            results[name][qid] = (status, detail)

    # Summary counts
    print(f"{'Model':10s} {'PASS':>5s} {'FAIL':>5s} {'NOT_ANS':>8s} {'ADV_PARTIAL':>12s} {'STRICT_N':>9s} {'STRICT_%':>9s}")
    for name in MODELS:
        c = results[name]
        strict_qs = [q for q in QORDER if q not in LENIENT_QS]
        n_pass = sum(1 for q in strict_qs if c[q][0] == "PASS")
        n_fail = sum(1 for q in strict_qs if c[q][0] == "FAIL")
        n_na = sum(1 for q in strict_qs if c[q][0] == "NOT_ANSWERED")
        n_adv = sum(1 for q in strict_qs if c[q][0] == "PARTIAL_ADVERSARIAL_CORRECT")
        n_strict = len(strict_qs)
        pct = 100.0 * n_pass / n_strict
        print(f"{name:10s} {n_pass:5d} {n_fail:5d} {n_na:8d} {n_adv:12d} {n_strict:9d} {pct:8.1f}%")

    print()
    print("=== Per-question detail (strict-graded questions only) ===")
    for qid in QORDER:
        if qid in LENIENT_QS:
            continue
        line = f"{qid}: "
        cells = []
        for name in MODELS:
            status, detail = results[name][qid]
            tag = {"PASS": "P", "FAIL": "F", "NOT_ANSWERED": "NA", "PARTIAL_ADVERSARIAL_CORRECT": "ADV-OK"}[status]
            cells.append(f"{name}={tag}")
        print(line + " ".join(cells))

    print()
    print("=== Adversarial questions (Q53, Q56) detail ===")
    for qid in ADVERSARIAL_QS:
        print(f"-- {qid} --")
        for name in MODELS:
            status, detail = results[name][qid]
            print(f"  {name}: {status} ({detail})")

    print()
    print("=== Lenient questions (Q19,35,41,42,44,46,52,60) informational ===")
    for qid in sorted(LENIENT_QS):
        print(f"-- {qid} --")
        for name in MODELS:
            status, detail = results[name][qid]
            print(f"  {name}: {detail if detail else status}")

    return results


if __name__ == "__main__":
    main()
