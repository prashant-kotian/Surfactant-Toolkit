#!/usr/bin/env python3
"""build_bench200.py — grand tables for SurfBench-200 (single benchmark).

Consumes:
  benchmark/question_set.json      (ids, types, subconcepts; run build_question_set.py first)
  results/results.json             (30 per-set cells from the four authoring graders;
                                    source-set names are internal provenance only)
  benchmark/source_sets/orig_hard50/*_grade_results.json (fourth authoring set)

Emits:
  results/bench200.json                machine-readable grand + class data
  results/TABLE200.md                  headline: % pass, augmented vs unaugmented
  results/CLASSIFICATION200.md         skill-type and sub-concept matrices + conclusions
  results/per_question_results_200.csv one row per question x model x condition

All numbers derive from committed grading JSONs; nothing is typed by hand.
Exclusion policy is the one declared in build_results.py / question_set.json and
is applied identically in both conditions; orig NOT_ANSWERED records are
ungradeable (no answer content) and are excluded from n, flagged in output.
"""
from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent          # benchmark/
ROOT = HERE.parent
RESULTS = ROOT / "results"
SET_ALIAS = {"groundzero": "groundzero", "hard50": "hard50", "extrahard": "extrahard"}
MODEL_ORDER = ["Claude", "ChatGPT", "DeepSeek", "Qwen", "Gemini", "Gemini-CLI"]
NAME_FIX = {"Antigravity": "Gemini-CLI", "Gemini CLI": "Gemini-CLI", "GPT-4o": "ChatGPT"}


def norm(m: str) -> str:
    return NAME_FIX.get(m, m)


def main() -> int:
    qs = json.loads((HERE / "question_set.json").read_text(encoding="utf-8"))
    alias2sb = {it["alias"]: it for it in qs["items"]}
    sb_sorted = [it["sb"] for it in qs["items"]]

    # records[(model, cond)][sb] = (passed: bool, note: str|None)
    rec: dict = defaultdict(dict)
    notes: dict = defaultdict(list)

    rj = json.loads((RESULTS / "results.json").read_text(encoding="utf-8"))
    for key, cell in rj["cells"].items():
        setname, model, cond = key.split("|")
        model = norm(model)
        flag = ("partial transcript" if cell.get("partial") else
                "count-only grading" if cell.get("count_only") else
                "SI-only grading" if cell.get("si_only") else None)
        for alias, passed in cell["per_q"].items():
            it = alias2sb.get(alias)
            if it is None:
                raise SystemExit(f"question_set missing alias {alias}")
            rec[(model, cond)][it["sb"]] = (bool(passed), flag)
        if flag:
            notes[(model, cond)].append(f"{SET_ALIAS.get(setname, setname)}: {flag}")

    for cond, fname in (("unaug", "unaugmented_grade_results.json"), ("aug", "augmented_grade_results.json")):
        d = json.loads((HERE / "source_sets/orig_hard50" / fname).read_text(encoding="utf-8"))
        for model, md in d.items():
            model = norm(model)
            na = 0
            for alias, q in md["questions"].items():
                v = q["verdict"]
                it = alias2sb.get(alias)
                if it is None:
                    raise SystemExit(f"question_set missing alias {alias}")
                if v == "NOT_ANSWERED":
                    na += 1
                    continue
                rec[(model, cond)][it["sb"]] = (v == "PASS", None)
            if na:
                notes[(model, cond)].append(f"orig_hard50: {na} unanswered (ungradeable, excluded from n)")

    def stats(subset_sbs):
        out = {}
        for model in MODEL_ORDER:
            for cond in ("unaug", "aug"):
                cell = rec.get((model, cond), {})
                rows = [cell[s] for s in subset_sbs if s in cell]
                if not rows:
                    out[(model, cond)] = None
                    continue
                n = len(rows)
                p = sum(1 for x, _ in rows if x)
                out[(model, cond)] = (p, n, round(100.0 * p / n, 1))
        return out

    grand = stats(sb_sorted)
    per_set_stats = {}
    for setname, sname in (("groundzero", "groundzero"), ("hard50", "hard50"),
                           ("extrahard", "extrahard"), ("orig_hard50", "orig_hard50")):
        sbs = [it["sb"] for it in qs["items"] if it["source_set"] == sname]
        per_set_stats[setname] = (sbs, stats(sbs))
    by_type = {t: stats([it["sb"] for it in qs["items"] if it["type"] == t])
               for t in ("conceptual", "computational", "analytical")}
    subs = sorted({it["subconcept"] for it in qs["items"]})
    by_sub = {s: stats([it["sb"] for it in qs["items"] if it["subconcept"] == s])
              for s in subs}
    sub_n = {s: sum(1 for it in qs["items"] if it["subconcept"] == s) for s in subs}

    def fmt(st):
        return "—" if st is None else f"{st[0]}/{st[1]} ({st[2]}%)"

    # ---------- TABLE200.md ----------
    L = ["# SurfBench-200 — grand results",
         "",
         "Single benchmark, 200 questions (the four authoring sets are internal provenance only).",
         "`n` = gradable questions after the declared exclusion policy (see DATA_NOTES.md);",
         "`Aug` = toolkit/MCP-augmented run, `Una` = unaugmented baseline. Lift = Aug − Una (pp).",
         "",
         "| Model | Unaugmented | Augmented | Lift | Gradable n /200 | Flags |",
         "|---|---|---|---|---|---|"]
    for m in MODEL_ORDER:
        u, a = grand.get((m, "unaug")), grand.get((m, "aug"))
        lift = f"+{a[2] - u[2]:.1f} pp" if u and a else "—"
        nf = u[1] if u else (a[1] if a else 0)
        seen = {}
        for c in ("unaug", "aug"):
            for nt in notes.get((m, c), []):
                seen[nt] = seen.get(nt, 0) + 1
        fl = "; ".join(f"{k}{'' if v == 1 else ' (both arms)'}" for k, v in seen.items()) or "—"
        L.append(f"| {m} | {fmt(u)} | {fmt(a)} | {lift} | {nf} | {fl} |")
    L += ["", "## Per-authoring-set breakdown (provenance; not separate benchmarks)", "",
          "| Model | " + " | ".join(f"{s} U/A" for s in per_set_stats) + " |",
          "|---|" + "---|" * (len(per_set_stats) + 1)]
    for m in MODEL_ORDER:
        cells = []
        for setname, (sbs, st) in per_set_stats.items():
            u, a = st.get((m, "unaug")), st.get((m, "aug"))
            cells.append(f"{u[0]}/{u[1]} · {a[0]}/{a[1]}" if u and a else
                         (f"{u[0]}/{u[1]} · —" if u else (f"— · {a[0]}/{a[1]}" if a else "—")))
        L.append(f"| {m} | " + " | ".join(cells) + " |")
    (RESULTS / "TABLE200.md").write_text("\n".join(L) + "\n", encoding="utf-8")

    # ---------- CLASSIFICATION200.md ----------
    def matrix(title, table, nn):
        C = [f"## {title}", "",
             "| Category | n | Claude U→A | ChatGPT U→A | DeepSeek U→A | Qwen U→A | Gemini U→A | GCLI U→A | Lift (core 3) |",
             "|---|---|---|---|---|---|---|---|---|---|"]
        concl = []
        for cat, st in table.items():
            row = [f"| {cat} | {nn.get(cat, '')} |"]
            lifts = []
            for m in MODEL_ORDER:
                u, a = st.get((m, "unaug")), st.get((m, "aug"))
                fu = "—" if u is None else f"{u[2]}%"
                fa = "—" if a is None else f"{a[2]}%"
                if u and a:
                    row.append(f"{fu}→{fa} |")
                    if m in ("Claude", "ChatGPT", "DeepSeek"):
                        lifts.append(a[2] - u[2])
                elif u:
                    row.append(f"{fu}→— |")
                else:
                    row.append("—→" + fa + " |" if a else "— |")
            tot = ""
            if lifts:
                avg = sum(lifts) / len(lifts)
                tot = f"{avg:+.1f} pp"
                pu = sum(st[(m, "unaug")][0] for m in ("Claude", "ChatGPT", "DeepSeek") if st.get((m, "unaug")))
                pn = sum(st[(m, "unaug")][1] for m in ("Claude", "ChatGPT", "DeepSeek") if st.get((m, "unaug")))
                pa = sum(st[(m, "aug")][0] for m in ("Claude", "ChatGPT", "DeepSeek") if st.get((m, "aug")))
                pan = sum(st[(m, "aug")][1] for m in ("Claude", "ChatGPT", "DeepSeek") if st.get((m, "aug")))
                if pn and pan:
                    concl.append((cat, round(100.0 * pu / pn, 1), round(100.0 * pa / pan, 1), pn, avg))
            else:
                concl.append((cat, None, None, 0, None))
            C.append("".join(row) + f" {tot} |")
        C.append("")
        return C, concl

    C, tconcl = matrix("Skill type", by_type, {t: sum(1 for i in qs["items"] if i["type"] == t) for t in by_type})
    S, sconcl = matrix("Sub-concept", by_sub, sub_n)
    lines = ["# SurfBench-200 — classification matrix",
             "",
             "Pass % per category, unaugmented→augmented (core-3 lift = mean pp over",
             "Claude/ChatGPT/DeepSeek, the arms run in both conditions). Rows pool every",
             "gradable model×condition cell; `n` is benchmark questions in the category.",
             ""] + C + S

    scored = [c for c in sconcl if c[1] is not None and c[3] >= 30]
    worst_u = sorted(scored, key=lambda c: c[1])[:5]
    worst_a = sorted(scored, key=lambda c: c[2])[:5]
    biggest_gap = sorted([c for c in scored if c[4] is not None], key=lambda c: c[4])[-3:]
    least_help = sorted([c for c in scored if c[4] is not None], key=lambda c: c[4])[:3]
    lines += ["## Conclusions (auto-generated from the matrices above)", "",
              f"- Weakest sub-concepts unaugmented: " + "; ".join(f"{c[0]} ({c[1]}%)" for c in worst_u) + ".",
              f"- Weakest even augmented (persistent blind spots): " + "; ".join(f"{c[0]} ({c[2]}%)" for c in worst_a) + ".",
              f"- Largest augmentation benefit: " + "; ".join(f"{c[0]} ({c[4]:+.1f} pp)" for c in biggest_gap[::-1]) + ".",
              f"- Where the toolkit helps least: " + "; ".join(f"{c[0]} ({c[4]:+.1f} pp)" for c in least_help) + ".",
              f"- Skill-type split — " + "; ".join(
                  f"{t}: {st[('Claude','unaug')][2] if st.get(('Claude','unaug')) else '—'}%→{st[('Claude','aug')][2] if st.get(('Claude','aug')) else '—'}% (Claude)"
                  for t, st in by_type.items()) + ".",
              ""]
    (RESULTS / "CLASSIFICATION200.md").write_text("\n".join(lines), encoding="utf-8")

    # ---------- CSV ----------
    with open(RESULTS / "per_question_results_200.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["sb_id", "authoring_alias", "source_set", "type", "subconcept",
                    "model", "condition", "passed", "flag"])
        for it in qs["items"]:
            for m in MODEL_ORDER:
                for c in ("unaug", "aug"):
                    cell = rec.get((m, c), {}).get(it["sb"])
                    if cell is None:
                        continue
                    p, flag = cell
                    w.writerow([it["sb"], it["alias"], it["source_set"], it["type"],
                                it["subconcept"], m, c, int(p), flag or ""])

    # ---------- bench200.json ----------
    def pack(st):
        return {f"{m}|{c}": (list(v) if v else None) for (m, c), v in st.items()}
    json.dump({"benchmark": qs["benchmark"], "n": qs["n"],
                "excluded": qs["excluded_from_scoring"],
                "grand": pack(grand), "per_set": {k: pack(v[1]) for k, v in per_set_stats.items()},
                "by_type": {t: pack(v) for t, v in by_type.items()},
                "by_subconcept": {k: pack(v) | {"_n": sub_n[k]} for k, v in by_sub.items()}},
              open(RESULTS / "bench200.json", "w"), indent=1)

    for (m, c), v in sorted(grand.items(), key=lambda kv: (MODEL_ORDER.index(kv[0][0]), kv[0][1])):
        if v:
            print(f"{m:11} {c:6} {v[0]:3}/{v[1]:3}  {v[2]:5.1f}%")
    print("rows_csv:", sum(1 for _ in open(RESULTS / "per_question_results_200.csv")) - 1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
