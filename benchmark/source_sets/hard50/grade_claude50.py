"""Grade the Claude HARD50 full-solutions transcript (5th model) with the SAME
DESIGN.md §7 rules as grade_hard50.py. Writes claude_unaugmented_grade_results.json and
prints a full analysis (baseline + audit-adjusted views).

Parser: per-question dividers `####…####` / `  H-NN` / `####…####` (49 present;
H-24 absent from the source file the model worked from).
"""
from __future__ import annotations
import json
import re
import sys
from pathlib import Path

import grade_unaugmented_legacy as gh

HERE = Path(__file__).resolve().parent
# Recovered transcript (2026-10-18), all 50 questions present (incl. H-24, which
# the original ####-divided source was missing); committed in-repo 2026-10-19 so
# this cell regrades from the repository alone.
CLAUDE_FILE = HERE / "claude_unaugmented" / "HARD50_Solved_Outputs_H01-H50.txt"

# divider-based split: 74 '#' line, id, 74 '#' line (legacy format, unused
# for the recovered file below but kept for reference)
DIV = r"^#+[ \t]*\n[ \t]*(H-\d+)[ \t]*\n#+[ \t]*$"


def parse_claude(text: str) -> dict[str, str]:
    """Recovered file uses trailing 'END OF QUESTION H-NN' markers, so the
    content for each ID is the text BEFORE its marker, not after."""
    parts = re.split(r"^END OF QUESTION (H-\d+)[ \t]*$", text, flags=re.M)
    out: dict[str, str] = {}
    for i in range(1, len(parts), 2):
        out[parts[i]] = parts[i - 1]
    return out


# ---------------- audit-adjusted view (same three fixes as analyze_grading.py) ----
AUDIT_ADJ = {
    "H-04": {"delta_H_kJ_per_mol": -2.846400782241967,
             "delta_S_J_per_mol_K": 115.12279612768563},
    "H-17": {"solubility_estimate_M_rel_tol": 0.5},
    "H-18": {"charge_type_at_ph_keywords": ["ambiguous", "mixed", "both",
                                            "coexist", "intermediate"]},
}


def regrade_field(qid, f, fd, got, adjusted):
    import re as _re
    if adjusted and qid == "H-04" and f in ("delta_H_kJ_per_mol", "delta_S_J_per_mol_K"):
        fd = dict(fd, value=AUDIT_ADJ["H-04"][f])
    if adjusted and qid == "H-17" and f == "solubility_estimate_M":
        fd = dict(fd, rel_tol=AUDIT_ADJ["H-17"]["solubility_estimate_M_rel_tol"])
    if adjusted and qid == "H-18" and f == "charge_type_at_ph":
        if got is None:
            return False
        hay = str(got).lower()
        return any(k in hay for k in AUDIT_ADJ["H-18"]["charge_type_at_ph_keywords"])
    t = fd["type"]
    if t == "number":
        if isinstance(got, bool) or got is None:
            return False
        if isinstance(got, str):
            try:
                got = float(got.strip())
            except Exception:
                return False
        g, tol = float(fd["value"]), fd.get("rel_tol", 0.01)
        return abs(float(got) - g) <= tol * abs(g)
    return gh.grade_field(fd, got)


def main():
    text = CLAUDE_FILE.read_text(encoding="utf-8", errors="replace")
    blocks = parse_claude(text)
    print(f"[parse] claude: {len(blocks)}/50 question blocks "
          f"(missing: {sorted(set(gh.GOLD) - set(blocks))})")

    results = gh.grade_model("claude", blocks)

    # audit-adjusted view
    adj = {}
    for qid, e in results.items():
        e2 = dict(e, fields={})
        all_ok = e.get("pass", False) if not e.get("fields") else True
        for f, d in e.get("fields", {}).items():
            ok = regrade_field(qid, f, gh.GOLD[qid]["fields"][f], d["got"], adjusted=True)
            e2["fields"][f] = dict(d, ok=ok)
            all_ok = all_ok and ok
        e2["pass"] = bool(e.get("fields")) and all_ok
        adj[qid] = e2

    summary = gh.summarize({"claude": results})
    summary_adj = gh.summarize({"claude": adj})

    # cross-check with the existing 4-model results (gold suspects across 5 models)
    # optional: the first-release 4-model results are not shipped in this repo
    old_p = HERE / "grade_hard50_results.json"
    merged = dict(json.load(open(old_p))["results"]) if old_p.exists() else {}
    merged["claude"] = results
    suspects = gh.cross_check(merged)

    json.dump({"summary": {k: {kk: vv for kk, vv in v.items() if kk != "axis_q"}
                           for k, v in summary.items()},
               "summary_audit_adj": {k: {kk: vv for kk, vv in v.items() if kk != "axis_q"}
                                     for k, v in summary_adj.items()},
               "results": results, "results_audit_adj": adj,
               "suspects_5models": suspects},
              open(HERE / "claude_unaugmented_grade_results.json", "w"), indent=1, default=str)

    s, sa = summary["claude"], summary_adj["claude"]
    print("\n================ CLAUDE SUMMARY ================")
    print(f"  baseline   pass: {len(s['passed'])}/{s['n']} = {s['pass_pct']:.1f}%  "
          f"(fail {s['n'] - len(s['passed'])})")
    print(f"  audit-adj  pass: {len(sa['passed'])}/{sa['n']} = {sa['pass_pct']:.1f}%  "
          f"(fail {sa['n'] - len(sa['passed'])})")
    print(f"  arith fields failed: {s['arith_fields'][0]}/{s['arith_fields'][1]} "
          f"= {s['arith_fail_pct']:.1f}%")
    print(f"  reason fields failed: {s['reason_fields'][0]}/{s['reason_fields'][1]} "
          f"= {s['reason_fail_pct']:.1f}%")
    if s["errors"]:
        print(f"  extraction errors: {s['errors']}")

    # per-question detail
    print("\n per-question (baseline | audit-adj):")
    causes = {"arith_only": 0, "reason_only": 0, "both": 0, "error": 0, "pass": 0}
    causes_adj = dict(causes)
    axis_fail = {}
    fam_fail = {}
    fam_tot = {}
    for qid in sorted(results):
        e, ea = results[qid], adj[qid]

        def cause(en):
            if en.get("error"):
                return "error"
            fa = any(not d["ok"] and d["bucket"] == "arithmetic" for d in en["fields"].values())
            fr = any(not d["ok"] and d["bucket"] == "reasoning" for d in en["fields"].values())
            if fa and fr:
                return "both"
            if fa:
                return "arith_only"
            if fr:
                return "reason_only"
            return "pass"

        c, ca = cause(e), cause(ea)
        causes[c] += 1
        causes_adj[ca] += 1
        meta = gh.META.get(qid, {})
        fam_tot[meta.get("family")] = fam_tot.get(meta.get("family"), 0) + 1
        if not ea["pass"]:
            fam_fail[meta.get("family")] = fam_fail.get(meta.get("family"), 0) + 1
            ax = meta.get("failure_axis")
            axis_fail[ax] = axis_fail.get(ax, 0) + 1
        fails = [f for f, d in e["fields"].items() if not d["ok"]]
        flag = "" if e["pass"] else ("*" if ea["pass"] else "")
        print(f"  {qid} {'PASS' if e['pass'] else 'FAIL'}{flag} [{c}"
              f"{'' if c == ca else ' -> ' + ca}]"
              + (f"  missed: {fails}" if fails else ""))

    print("\n causes baseline :", causes)
    print(" causes audit-adj:", causes_adj)
    print(" family fails (audit-adj):", {k: f"{v}/{fam_tot[k]}" for k, v in sorted(fam_fail.items(), key=lambda x: -x[1])})
    print(" axis fails (audit-adj):", dict(sorted(axis_fail.items(), key=lambda x: -x[1])))
    print(f"\n gold suspects across 5 models: {len(suspects)}")
    for x in suspects:
        print("  ", x["qid"], x["field"], "gold=", x["gold"], "agree=", x["agreement"],
              "passed_by=", x["passed_by"])
    print("\n-> claude_unaugmented_grade_results.json")


if __name__ == "__main__":
    main()
