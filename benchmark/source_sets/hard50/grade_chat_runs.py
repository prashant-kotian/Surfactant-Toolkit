"""Grade the new per-question H50(new) runs from this session: Claude (aug),
ChatGPT (aug). Reuses grade_hard50's extractor/comparator. Handles Claude's
curly-quote JSON rendering."""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from grade_unaugmented_legacy import GOLD, extract_final_json, grade_field

HERE = Path(__file__).resolve().parent
SESSION = HERE.parent.parent / "transcripts"


def normalize_quotes(text):
    return text.replace("“", '"').replace("”", '"').replace("‘", "'").replace("’", "'")


def grade_condition(label, dir_name):
    results = {}
    for n in range(1, 51):
        qid = f"H-{n:02d}"
        p = SESSION / dir_name / f"H{n:02d}.txt"
        if not p.exists():
            results[qid] = {"verdict": "NO_JSON", "fields": {}}
            continue
        text = normalize_quotes(p.read_text(encoding="utf-8", errors="replace"))
        ans = extract_final_json(text)
        if ans is None:
            results[qid] = {"verdict": "NO_JSON", "fields": {}}
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
    n_fail = sum(1 for r in results.values() if r["verdict"] == "FAIL")
    n_nojson = sum(1 for r in results.values() if r["verdict"] == "NO_JSON")
    print(f"\n=== H50(new) {label}: PASS={n_pass} FAIL={n_fail} NO_JSON={n_nojson} / 50 ===")
    for qid, r in results.items():
        if r["verdict"] != "PASS":
            bad = [k for k, fv in r.get("fields", {}).items() if not fv["ok"]]
            print(f"  {qid}: {r['verdict']}  {bad}")
    out_path = HERE / f"grade_results_session_{label.replace(' ', '_')}.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=1, default=str)
    return n_pass, results


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    grade_condition("claude_augmented", "claude_augmented_hard50")
    grade_condition("chatgpt_augmented", "chatgpt_augmented_hard50")
