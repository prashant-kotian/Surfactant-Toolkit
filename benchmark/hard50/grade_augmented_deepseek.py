"""Grade the DeepSeek (Better DeepSeek extension + MCP) augmented HARD50 run.
Reuses grade_hard50.py's extract_final_json/grade_field (the same balanced-
brace extractor and field comparator used for the unaugmented baseline).

Usage: python grade_augmented_deepseek_hard50.py
"""
import json
from pathlib import Path

from grade_hard50 import GOLD, extract_final_json, grade_field, field_bucket

HERE = Path(__file__).resolve().parent
ANSWERS_DIR = HERE / "augmented_runs_deepseek_2026-10-15" / "hard50"

results = {}
for qid, gold in sorted(GOLD.items()):
    num = int(qid.split("-")[1])
    fpath = ANSWERS_DIR / f"H{num:02d}.txt"
    entry = {"found_block": fpath.exists(), "fields": {}, "pass": False}
    if not fpath.exists():
        entry["error"] = "no answer file (known limitation: Better DeepSeek "\
            "markdown-corrupts this SMILES before the tool call)" if num == 30 else "no answer file"
        results[qid] = entry
        continue
    block = fpath.read_text(encoding="utf-8")
    ans = extract_final_json(block)
    if ans is None:
        entry["error"] = "no FINAL_JSON"
        results[qid] = entry
        continue
    if "__unparseable__" in ans:
        entry["error"] = "unparseable FINAL_JSON"
        results[qid] = entry
        continue
    entry["answer"] = ans
    all_ok = True
    for f, fd in gold["fields"].items():
        ok = grade_field(fd, ans.get(f, "__MISSING__"))
        entry["fields"][f] = {
            "ok": ok, "gold": fd.get("value"), "got": ans.get(f, "__MISSING__"),
            "bucket": field_bucket(fd), "type": fd["type"],
        }
        all_ok = all_ok and ok
    entry["pass"] = all_ok
    results[qid] = entry

n_pass = sum(1 for r in results.values() if r["pass"])
n_total = len(results)
print(f"DeepSeek augmented HARD50: {n_pass}/{n_total} = {100*n_pass/n_total:.1f}%")
print()
for qid, r in sorted(results.items()):
    status = "PASS" if r["pass"] else "FAIL"
    detail = r.get("error", "")
    if not r["pass"] and "fields" in r and r["fields"]:
        bad = [f for f, fv in r["fields"].items() if not fv["ok"]]
        detail = f"bad fields: {bad}"
    print(f"  {qid}: {status}  {detail}")

with open(HERE / "deepseek_augmented_grade_results.json", "w", encoding="utf-8") as f:
    json.dump(results, f, indent=1, default=str)
