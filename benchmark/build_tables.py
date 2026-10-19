#!/usr/bin/env python3
"""build_tables.py — Phase 2: generate every manuscript results artifact from
../results/results.json. Numbers are NEVER hand-typed: rerun after each re-grade.

Outputs (paper3/):
  tables.md / .tex   main tables + exclusion + provenance appendix
  per_question_results.csv          per-cell per-question pass matrix (SI)

Flag policy (PENDING = manual overlay; rest auto-derived from master):
  * si_only / count_only / PRE-FIX  -> automatic footnote per cell
  * entries in PENDING              -> §1b recovery, agy tier, known gaps
"""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent            # benchmark/
ROOT = HERE.parent
OUT = ROOT / "results"
MASTER = json.loads((HERE / "../results/results.json").read_text(encoding="utf-8"))
CELLS: dict = MASTER["cells"]

MODELS = ["Claude", "ChatGPT", "DeepSeek", "Gemini", "Qwen", "Gemini-CLI"]
DISPLAY = {"Gemini-CLI": "Antigravity CLI", "Gemini": "Gemini (web)"}
BANKS = ["GZ50", "H50new", "X50"]
BANKNAME = {"GZ50": "GroundZero-50", "H50new": "Hard50 (new)", "X50": "ExtraHard-50"}

# Manual overlay — remove entries as the corresponding backend steps land.
PENDING = {   # §1b RESOLVED at merge 4ccbd87 (full recaptures, re-graded 2026-10-18)
    "X50|Gemini-CLI|aug": "pro-high tier (retired) — disclosed, final",
    "X50|Gemini-CLI|unaug": "pro-high tier (retired) — disclosed, final",
}


def flags(key: str) -> list[str]:
    c = CELLS[key]
    f = []
    if c.get("si_only"):
        f.append("SI-only")
    if c.get("count_only"):
        f.append("count-only")
    if c.get("partial"):
        f.append(f"partial coverage (n={c['n']})")
    if c.get("gold_state") != "post-fix":
        f.append("PRE-FIX gold")
    if key in PENDING:
        f.append(PENDING[key])
    return f


def cell(bank: str, model: str, cond: str) -> dict | None:
    return CELLS.get(f"{bank}|{model}|{cond}")


def build_table(cond: str) -> tuple[list[str], list[str]]:
    """-> (header, rows) with one row per (bank, model); '-' when no cell."""
    hdr = ["Bank", "Model", "passed", "n", "%", "flags"]
    rows = []
    for bank in BANKS:
        for m in MODELS:
            c = cell(bank, m, cond)
            if c is None:
                continue
            fl = flags(f"{bank}|{m}|{cond}")
            rows.append([BANKNAME[bank], DISPLAY.get(m, m), str(c["passed"]), str(c["n"]),
                         f"{c['pct']:.1f}", "; ".join(fl)])
    return hdr, rows


def build_delta() -> tuple[list[str], list[str]]:
    hdr = ["Bank", "Model", "unaug %", "aug %", "Δ pts", "comparability"]
    rows = []
    for bank in BANKS:
        for m in MODELS:
            u, a = cell(bank, m, "unaug"), cell(bank, m, "aug")
            if not (u and a):
                continue
            d = a["pct"] - u["pct"]
            comp = []
            if u["harness"] != a["harness"]:
                comp.append("harness differs")
            for tag in (f"{bank}|{m}|unaug", f"{bank}|{m}|aug"):
                for fl in flags(tag):
                    if "tier" in fl or "PRE-FIX" in fl or "§1b" in fl:
                        comp.append(fl.split(" —")[0].split(" (")[0])
            rows.append([BANKNAME[bank], DISPLAY.get(m, m), f"{u['pct']:.1f}",
                         f"{a['pct']:.1f}", f"{d:+.1f}", "; ".join(sorted(set(comp)) or ["ok"])])
    return hdr, rows


def build_exclusions() -> tuple[list[str], list[str]]:
    hdr = ["Bank", "n bank", "n cells", "excluded", "reason"]
    rows = []
    for bank in BANKS:
        exs, reasons = set(), set()
        for k, c in CELLS.items():
            if k.startswith(bank + "|"):
                exs |= set(c["excluded"])
                if c.get("exclude_reason"):
                    reasons.add(c["exclude_reason"])
        ns = {c["n"] for k, c in CELLS.items() if k.startswith(bank + "|")}
        rows.append([BANKNAME[bank], "50", "/".join(str(x) for x in sorted(ns)),
                     ", ".join(sorted(exs)) or "—", " | ".join(sorted(reasons)) or "—"])
    return hdr, rows


def build_provenance() -> tuple[list[str], list[str]]:
    hdr = ["cell", "passed/n", "%", "harness", "source", "sha", "gold", "notes"]
    rows = []
    for k in sorted(CELLS):
        c = CELLS[k]
        rows.append([k.replace("|", " / "), f"{c['passed']}/{c['n']}", f"{c['pct']:.1f}",
                     c["harness"], c["source"], c["source_sha"], c["gold_state"],
                     "; ".join(flags(k))])
    return hdr, rows


def md_table(hdr, rows) -> str:
    out = ["| " + " | ".join(hdr) + " |", "|" + "---|" * len(hdr)]
    out += ["| " + " | ".join(r) + " |" for r in rows]
    return "\n".join(out)


def tex_table(hdr, rows, caption, label) -> str:
    def esc(s):
        return (s.replace("&", r"\&").replace("_", r"\_").replace("%", r"\%")
                 .replace("§", r"S\textsection "))
    lines = [rf"\begin{{table}}[t]\centering\small\caption{{{esc(caption)}}}\label{{{label}}}",
             r"\begin{tabular}{" + "l" * len(hdr) + "}", r"\toprule",
             " & ".join(esc(h) for h in hdr) + r" \\", r"\midrule"]
    for r in rows:
        lines.append(" & ".join(esc(x) for x in r) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    return "\n".join(lines)


def render() -> dict[str, str]:
    t = {}
    t["tbl:unaug"] = md_table(*build_table("unaug"))
    t["tbl:aug"] = md_table(*build_table("aug"))
    t["tbl:delta"] = md_table(*build_delta())
    t["tbl:exclusions"] = md_table(*build_exclusions())
    t["tbl:provenance"] = md_table(*build_provenance())
    tex = [r"% generated by benchmark/build_tables.py -- do not hand-edit", ""]
    for cap, lab, key in [("Unaugmented (pure-reasoning) condition", "tab:unaug", "unaug"),
                          ("Augmented (tool-harness) condition", "tab:aug", "aug"),
                          ("Augmentation lift", "tab:delta", "delta"),
                          ("Declared exclusions", "tab:excl", "exclusions"),
                          ("Provenance ledger (appendix)", "tab:prov", "provenance")]:
        fn = {"unaug": build_table, "aug": build_table, "delta": build_delta,
              "exclusions": build_exclusions, "provenance": build_provenance}[key]
        args = ("unaug",) if key in ("unaug",) else ("aug",) if key == "aug" else ()
        tex.append(tex_table(*fn(*args), cap, lab))
        tex += ["", ""]
    t["tex:all"] = "\n".join(tex)
    for k in sorted(CELLS):
        c = CELLS[k]
        t[f"cell:{k}"] = f"{c['passed']}/{c['n']} ({c['pct']:.1f}\\%)"
        t[f"cellpct:{k}"] = f"{c['pct']:.1f}"
        t[f"cellpass:{k}"] = str(c["passed"])
        t[f"celln:{k}"] = str(c["n"])
    for nm, sec in {"unaug": "Unaugmented", "aug": "Augmented", "delta": "Augmentation lift",
                    "exclusions": "Declared exclusions", "provenance": "Provenance ledger"}.items():
        t[f"ref:tbl:{nm}"] = f"tables.md — {sec}"
    gaps = MASTER["gaps"]
    t["meta:gaps"] = f"{len(gaps)} declared gap(s)"
    t["meta:toolkit"] = MASTER["toolkit_state"]
    return t


def main() -> int:
    OUT.mkdir(exist_ok=True)
    t = render()
    (OUT / "tables.md").write_text(
        "# Results tables (generated — rerun benchmark/build_tables.py)\n\n"
        "## Unaugmented\n\n" + t["tbl:unaug"] +
        "\n\n## Augmented\n\n" + t["tbl:aug"] +
        "\n\n## Augmentation lift\n\n" + t["tbl:delta"] +
        "\n\n## Declared exclusions\n\n" + t["tbl:exclusions"] +
        "\n\n## Provenance ledger\n\n" + t["tbl:provenance"] + "\n", encoding="utf-8")
    (OUT / "tables.tex").write_text(t["tex:all"], encoding="utf-8")
    with open(OUT / "per_question_results.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["cell", "question", "pass"])
        for k in sorted(CELLS):
            for q, p in CELLS[k]["per_q"].items():
                w.writerow([k, q, int(bool(p))])
    print(f"results/: tables+csv rendered from {len(CELLS)} cells "
          f"({len(PENDING)} pending flags active)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
