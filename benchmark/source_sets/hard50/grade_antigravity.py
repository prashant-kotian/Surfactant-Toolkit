"""Grade Antigravity (gemini-3.8-flash-high) HARD50 runs, both conditions.
Concatenates the 5 batch files per condition, splits on '=== H## ===' markers,
extracts each question's FINAL_JSON, and grades against gold_answers.json.
"""
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from grade_unaugmented_legacy import GOLD, extract_final_json, grade_field

HERE = Path(__file__).resolve().parent


def load_condition(raw_dir: Path) -> dict:
    """Concatenate batch1..5.txt, split on markers, return {qid: block_text}."""
    full_text = ""
    for i in range(1, 6):
        full_text += (raw_dir / f"batch{i}.txt").read_text(encoding="utf-8", errors="replace") + "\n"
    chunks = re.split(r"=== (H\d+) ===", full_text)
    blocks = {}
    for i in range(1, len(chunks), 2):
        qid_raw = chunks[i]
        content = chunks[i + 1]
        qid = f"H-{qid_raw[1:]}"
        # last marker wins if a qid appears more than once (shouldn't happen)
        blocks[qid] = content
    return blocks


def grade_condition(label: str, raw_dir: Path):
    blocks = load_condition(raw_dir)
    expected_ids = {f"H-{n:02d}" for n in range(1, 51)}
    missing = sorted(expected_ids - set(blocks.keys()), key=lambda x: int(x.split("-")[1]))
    results = {}
    for qid in sorted(expected_ids, key=lambda x: int(x.split("-")[1])):
        if qid not in blocks:
            results[qid] = {"verdict": "NO_JSON", "fields": {}, "reason": "question ID marker not found in transcript"}
            continue
        ans = extract_final_json(blocks[qid])
        if ans is None:
            results[qid] = {"verdict": "NO_JSON", "fields": {}, "reason": "no FINAL_JSON parsed"}
            continue
        gold = GOLD[qid]["fields"]
        field_results = {}
        all_ok = True
        for k, fd in gold.items():
            ok = grade_field(fd, ans.get(k, "__MISSING__"))
            field_results[k] = {"ok": ok, "gold": fd.get("value"), "got": ans.get(k, "__MISSING__")}
            all_ok = all_ok and ok
        results[qid] = {"verdict": "PASS" if all_ok else "FAIL", "fields": field_results}

    n_pass = sum(1 for r in results.values() if r["verdict"] == "PASS")
    print(f"\n=== HARD50 {label} (gemini-3.8-flash-high): {n_pass}/50 = {100*n_pass/50:.1f}% ===")
    if missing:
        print(f"  MISSING question IDs (not found in transcript at all): {missing}")
    for qid, r in results.items():
        if r["verdict"] != "PASS":
            bad = [k for k, fv in r.get("fields", {}).items() if not fv["ok"]]
            reason = r.get("reason", "")
            print(f"  {qid}: {r['verdict']}  {bad if bad else reason}")

    out_path = HERE / f"antigravity_{label}_grade_results.json"   # the file build_results.py reads
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=1, default=str)
    return results, n_pass


if __name__ == "__main__":
    grade_condition("augmented", HERE / "antigravity_aug_raw")
    grade_condition("unaugmented", HERE / "antigravity_unaug_raw")
