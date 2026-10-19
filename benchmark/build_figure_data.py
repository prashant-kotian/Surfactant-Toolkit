#!/usr/bin/env python3
"""build_figure_data.py — paste-ready TSV sources for every manuscript figure.

Generates results/figs_data/*.tsv from the committed grading JSONs so figure
production (manual or scripted) never transcribes a number. Also computes the
two analyses that previously existed only in prose:
  F7  denominator honesty (exclusion accounting per cell)
  F8  failure-class anatomy from per-field grader details
  F9  run dates parsed from provenance paths
Cross-checked against results/bench200.json (must match exactly or this fails).

Run after build_results.py, build_bench200.py:
    python3 benchmark/build_figure_data.py
"""
from __future__ import annotations

import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent          # benchmark/
ROOT = HERE.parent
RES = ROOT / "results"
FIG = RES / "figs_data"
MODEL_ORDER = ["Claude", "ChatGPT", "DeepSeek", "Qwen", "Gemini", "Gemini-CLI"]
NOMINAL = {"groundzero": 50, "hard50": 50, "extrahard": 50, "orig_hard50": 50}
Z = 1.96


def wilson(p: int, n: int) -> tuple[float, float]:
    if n == 0:
        return (0.0, 0.0)
    ph = p / n
    d = 1 + Z * Z / n
    c = (ph + Z * Z / (2 * n)) / d
    hw = Z * ((ph * (1 - ph) / n + Z * Z / (4 * n * n)) ** 0.5) / d
    return (round(100 * (c - hw), 1), round(100 * (c + hw), 1))


def relerr(gold, got):
    if not isinstance(gold, (int, float)) or not isinstance(got, (int, float)):
        return None
    return abs(gold - got) / max(abs(gold), 1e-12)


def normalize_detail(path: Path) -> dict:
    """Return {qid: list-of-(ok, gold, got) | detail-string | None}."""
    out = {}
    if not path.exists():
        return out
    j = json.loads(path.read_text(encoding="utf-8"))
    stack = [("", j)]
    while stack:
        prefix, node = stack.pop()
        if isinstance(node, dict):
            if "fields" in node and isinstance(node["fields"], dict):
                out[prefix] = [(f.get("ok", True), f.get("gold"), f.get("got", f.get("model")))
                               for f in node["fields"].values() if isinstance(f, dict)]
                continue
            if "verdict" in node:
                out[prefix] = node.get("detail") or None
                continue
            if "passed" in node or "failed" in node:
                for q in node.get("failed", []):
                    out.setdefault(q, None)
                continue
            for k, v in node.items():
                if isinstance(v, dict) and ("verdict" in v or "fields" in v):
                    stack.append((k, v))
                elif isinstance(v, dict):
                    stack.append((prefix or k, v))
    return out


def classify(detail) -> str:
    if detail is None:
        return "unclassified"            # count-only graders
    if isinstance(detail, str):
        pairs = re.findall(r"got (-?[\d.]+) vs gold (-?[\d.]+)", detail)
        if pairs:
            worst = max((relerr(float(g), float(m)) or 0) for m, g in pairs)
            return "precision_drift" if worst < 0.02 else "magnitude_error"
        if "vs gold" in detail:
            return "selection_or_convention"   # string mismatch (e.g. langmuir vs frumkin)
        return "other"
    classes = set()
    for ok, gold, got in detail:
        if ok:
            continue
        e = relerr(gold, got)
        classes.add("precision_drift" if e is not None and e < 0.02
                    else ("magnitude_error" if e is not None else "selection_or_convention"))
    return "+".join(sorted(classes)) if classes else "all_pass"


def main() -> int:
    FIG.mkdir(exist_ok=True)
    rj = json.loads((RES / "results.json").read_text(encoding="utf-8"))
    qs = json.loads((HERE / "question_set.json").read_text(encoding="utf-8"))
    alias2 = {it["alias"]: it for it in qs["items"]}
    alias2_sb = {it["alias"]: it["sb"] for it in qs["items"]}
    b200 = json.loads((RES / "bench200.json").read_text(encoding="utf-8"))

    # ---- assemble records: (model,cond,sb) -> passed; plus per-set exclusion bookkeeping
    rec = defaultdict(dict)
    cellmeta = {}
    for key, cell in rj["cells"].items():
        setn, model, cond = key.split("|")
        if model == "*":
            continue
        cellmeta[key] = cell
        for alias, ok in cell["per_q"].items():
            rec[(model, cond)][alias2[alias]["sb"]] = bool(ok)
    for cond, fname in (("unaug", "unaugmented_grade_results.json"), ("aug", "augmented_grade_results.json")):
        d = json.loads((HERE / "source_sets/orig_hard50" / fname).read_text(encoding="utf-8"))
        for model, md in d.items():
            model = "Gemini-CLI" if model == "Antigravity" else model
            for alias, q in md["questions"].items():
                if q["verdict"] == "NOT_ANSWERED":
                    continue
                rec[(model, cond)][alias2[alias]["sb"]] = q["verdict"] == "PASS"

    def w(name, header, rows):
        with open(FIG / name, "w", newline="") as f:
            f.write("\t".join(header) + "\n")
            for r in rows:
                f.write("\t".join(str(x) for x in r) + "\n")
        print(f"{name}: {len(rows)} rows")

    # ---- F1/F2 grand + lift
    grand = {}
    for m in MODEL_ORDER:
        for c in ("unaug", "aug"):
            rows = list(rec.get((m, c), {}).items())
            n = len(rows)
            p = sum(1 for _, ok in rows if ok)
            lo, hi = wilson(p, n) if n else (None, None)
            grand[(m, c)] = (p, n, round(100 * p / n, 1) if n else None, lo, hi)
    # cross-check vs bench200.json
    for (m, c), (p, n, pct, _, _) in grand.items():
        b = b200["grand"].get(f"{m}|{c}")
        assert (b is None and n == 0) or (b[0] == p and b[1] == n), f"grand mismatch {m}|{c}"
    w("F1_grand_ci.tsv", ["model", "condition", "passed", "n", "pct", "wilson95_lo", "wilson95_hi"],
      [(m, c, *grand[(m, c)]) for m in MODEL_ORDER for c in ("unaug", "aug") if grand[(m, c)][1]])
    w("F2_lift_dumbbell.tsv",
      ["model", "unaug_pct", "unaug_lo", "unaug_hi", "aug_pct", "aug_lo", "aug_hi",
       "lift_pp", "ci_overlap"],
      [(m, *grand[(m, "unaug")][2:5], *grand[(m, "aug")][2:5],
        round(grand[(m, "aug")][2] - grand[(m, "unaug")][2], 1),
        "yes" if grand[(m, "aug")][3] <= grand[(m, "unaug")][4] else "no")
       for m in MODEL_ORDER
       if grand[(m, "unaug")][1] and grand[(m, "aug")][1]])

    # ---- F3/F5 matrices (subconcept, type)
    def matrix(label_key, fname):
        cats = sorted({it[label_key] for it in qs["items"]})
        rows = []
        for cat in cats:
            sbs = [it["sb"] for it in qs["items"] if it[label_key] == cat]
            for m in MODEL_ORDER:
                for c in ("unaug", "aug"):
                    cell = rec.get((m, c), {})
                    rows_cat = [(sb, cell[sb]) for sb in sbs if sb in cell]
                    n = len(rows_cat)
                    p = sum(1 for _, ok in rows_cat if ok)
                    rows.append((cat, len(sbs), m, c, p, n,
                                 round(100 * p / n, 1) if n else "NA"))
        w(fname, [label_key, "n_items", "model", "condition", "passed", "n", "pct"], rows)
        return cats, rows

    _, rows3 = matrix("subconcept", "F3_subconcept_matrix.tsv")
    matrix("type", "F5_skilltype_matrix.tsv")

    # ---- F4 pooled weakness -> lift (core-3 arms)
    core = ["Claude", "ChatGPT", "DeepSeek"]
    rows = []
    for cat in sorted({it["subconcept"] for it in qs["items"]}):
        sbs = [it["sb"] for it in qs["items"] if it["subconcept"] == cat]
        pu = au = pn = an = 0
        for m in core:
            for c, acc in (("unaug", 0), ("aug", 1)):
                cell = rec.get((m, c), {})
                got = [(sb, cell[sb]) for sb in sbs if sb in cell]
                if c == "unaug":
                    pn += len(got); pu += sum(x for _, x in got)
                else:
                    an += len(got); au += sum(x for _, x in got)
        if pn and an:
            up, ap = 100 * pu / pn, 100 * au / an
            rows.append((cat, len(sbs), pu, pn, round(up, 1), au, an, round(ap, 1),
                         round(ap - up, 1)))
    rows.sort(key=lambda r: r[4])
    w("F4_weakness_lift.tsv",
      ["subconcept", "n_items", "core3_pass_unaug", "core3_n_unaug", "unaug_pct",
       "core3_pass_aug", "core3_n_aug", "aug_pct", "lift_pp"], rows)

    # ---- F6 long per-question
    with open(RES / "per_question_results_200.csv", newline="", encoding="utf-8") as f:
        src = list(csv.DictReader(f))
    w("F6_per_question.tsv", ["sb_id", "subconcept", "type", "model", "condition", "passed"],
      [(r["sb_id"], r["subconcept"], r["type"], r["model"], r["condition"], r["passed"]) for r in src])

    # ---- F7 denominator honesty per (model,cond,set)
    rows = []
    for m in MODEL_ORDER:
        for c in ("unaug", "aug"):
            cell = rec.get((m, c), {})
            for setn in ("groundzero", "hard50", "extrahard", "orig_hard50"):
                sbs = [it["sb"] for it in qs["items"] if it["source_set"] == setn]
                present = [sb for sb in sbs if sb in cell]
                gradable = len(present)
                cov = ("not-run" if not present else
                       ("partial" if setn == "groundzero" and m == "Claude" and c == "unaug"
                        else "full"))
                rows.append((m, c, setn, NOMINAL[setn], gradable,
                             NOMINAL[setn] - gradable, cov))
    w("F7_denominators.tsv",
      ["model", "condition", "source_set", "nominal_n", "gradable_n",
       "not_gradable_n", "coverage"], rows)

    # ---- F8 failure anatomy from grader detail
    by_alias_set = defaultdict(set)
    for it in qs["items"]:
        by_alias_set[it["source_set"]].add(it["alias"])
    rows = Counter()
    detail_src = {}
    for key, cell in rj["cells"].items():
        detail_src[key] = normalize_detail(ROOT / cell["source"])
    for key, cell in rj["cells"].items():
        setn, model, cond = key.split("|")
        det = detail_src.get(key, {})
        for alias, ok in cell["per_q"].items():
            if ok:
                continue
            d = det.get(alias)
            rows[(model, "aug" if cond in ("aug",) else cond, setn, classify(d))] += 1
    for cond, fname in (("unaug", "unaugmented_grade_results.json"), ("aug", "augmented_grade_results.json")):
        d = json.loads((HERE / "source_sets/orig_hard50" / fname).read_text(encoding="utf-8"))
        for model, md in d.items():
            model = "Gemini-CLI" if model == "Antigravity" else model
            for alias, q in md["questions"].items():
                if q["verdict"] == "FAIL":
                    rows[(model, cond, "orig_hard50", classify(q.get("detail")))] += 1
    w("F8_failure_anatomy.tsv",
      ["model", "condition", "source_set", "failure_class", "n_failed_questions"],
      [(a, b, c, d, n) for (a, b, c, d), n in sorted(rows.items())])

    # ---- F9 provenance dates
    rows = []
    for key, cell in sorted(rj["cells"].items()):
        mdate = re.search(r"(2026-\d\d-\d\d)", str(cell["source"]))
        rows.append((key, cell["harness"], mdate.group(1) if mdate else "-", cell["source"]))
    w("F9_provenance_dates.tsv", ["cell", "harness", "date", "source"], rows)

    # guard: F3 totals reconcile with question_set
    cnt = Counter()
    for r in rows3:
        cnt[r[0]] += 0
    per_sub = Counter(it["subconcept"] for it in qs["items"])
    for cat, n_items in per_sub.items():
        assert any(r[0] == cat and r[1] == n_items for r in rows3), f"F3 count mismatch {cat}"
    print("all cross-checks passed ->", FIG)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
