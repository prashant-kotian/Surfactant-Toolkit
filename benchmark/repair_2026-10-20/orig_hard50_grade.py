"""Grader for orig_hard50 (Q01-Q50) re-run answers, net-200 fill (2026-10-21).

GOLDS and num_close are copied verbatim from the first-release harness that produced the frozen
benchmark/source_sets/orig_hard50/*_grade_results.json (surfactantkit repo, benchmark/blind_run/
grade_blind_hard60.py at origin/main ce2a951). Scoring follows that harness's orig-HARD50 subset rule
(grade_orig_hard50_subset.py): flat strict, every numeric part within 10^-dp of the gold, string
parts by containment; no lenient or adversarial exemptions.
"""
from __future__ import annotations
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


def grade_question(qid: str, ans) -> tuple[str, str]:
    gold = GOLDS[qid]
    if ans is None:
        return "NOT_ANSWERED", ""
    ok_all, details = True, []
    for i, g in enumerate(gold):
        if isinstance(g, tuple):
            gv, dp = g
            if i >= len(ans) or not isinstance(ans[i], (int, float)) or isinstance(ans[i], bool):
                ok_all = False
                details.append(f"part{i + 1} missing/non-numeric")
            elif not num_close(ans[i], gv, dp):
                ok_all = False
                details.append(f"part{i + 1} got {ans[i]} vs gold {gv}")
        else:
            got = str(ans[i]).strip().lower() if i < len(ans) else ""
            if not got or not (g.lower() in got or got in g.lower()):
                ok_all = False
                details.append(f"part{i + 1} got '{got}' vs gold '{g}'")
    return ("PASS" if ok_all else "FAIL"), "; ".join(details)