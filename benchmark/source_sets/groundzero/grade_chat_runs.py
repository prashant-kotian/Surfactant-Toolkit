"""Grade the new per-question GZ50 runs from this session: Claude (aug),
ChatGPT (aug), Qwen (unaug), DeepSeek (unaug). Reuses grade_unaugmented's
extractor/comparator. Handles Claude's curly-quote JSON rendering.
"""
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from grade_unaugmented import extract_last_json_block, extract_templates, grade_field, META_KEYS

HERE = Path(__file__).resolve().parent
SESSION = HERE.parent.parent / "transcripts"
GOLD_PATH = HERE / "gold_answers.json"

with open(GOLD_PATH, encoding="utf-8") as f:
    gold_all = json.load(f)
gold_all = {k: v for k, v in gold_all.items() if k.startswith("GZ-") and 1 <= int(k.split("-")[1]) <= 50}
templates = extract_templates()


def normalize_quotes(text):
    return text.replace("“", '"').replace("”", '"').replace("‘", "'").replace("’", "'")


def load_question_file(dir_name, qid_prefix, num):
    qid_file = f"{qid_prefix}{num:02d}.txt"
    p = SESSION / dir_name / qid_file
    if not p.exists():
        return None
    text = p.read_text(encoding="utf-8", errors="replace")
    text = normalize_quotes(text)
    return extract_last_json_block(text)


def grade_condition(label, dir_name, qid_prefix="GZ"):
    results = {}
    for n in range(1, 51):
        qid = f"GZ-{n:02d}"
        gold_entry = gold_all[qid]
        gold = gold_entry["gold"]
        template_keys = templates.get(qid, set())
        graded_keys = (template_keys & set(gold.keys())) - META_KEYS

        model_answer = load_question_file(dir_name, qid_prefix, n)
        if model_answer is None:
            results[qid] = {"verdict": "NO_JSON", "fields": {}}
            continue

        field_results = {}
        all_pass = True
        for k in graded_keys:
            gv = gold[k]
            mv = model_answer.get(k) if isinstance(model_answer, dict) else None
            ok, kind = grade_field(mv, gv)
            field_results[k] = {"ok": ok, "kind": kind, "gold": gv, "model": mv}
            if ok is False:
                all_pass = False
        verdict = "PASS" if all_pass else "FAIL"
        results[qid] = {"verdict": verdict, "fields": field_results}

    n_pass = sum(1 for r in results.values() if r["verdict"] == "PASS")
    n_fail = sum(1 for r in results.values() if r["verdict"] == "FAIL")
    n_nojson = sum(1 for r in results.values() if r["verdict"] == "NO_JSON")
    print(f"\n=== GZ50 {label}: PASS={n_pass} FAIL={n_fail} NO_JSON={n_nojson} / 50 ===")
    for qid, r in results.items():
        if r["verdict"] != "PASS":
            bad = [k for k, fv in r.get("fields", {}).items() if fv.get("ok") is False]
            needs_j = [k for k, fv in r.get("fields", {}).items() if fv.get("ok") is None]
            print(f"  {qid}: {r['verdict']}  fail_fields={bad}  needs_judgment={needs_j}")

    out_path = HERE / f"grade_results_session_{label.replace(' ', '_')}.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=1, default=str)
    return n_pass, results


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    grade_condition("claude_augmented", "claude_augmented_groundzero50")
    grade_condition("chatgpt_augmented", "chatgpt_augmented_groundzero50")
    grade_condition("qwen_unaugmented", "qwen_unaugmented_groundzero50")
    grade_condition("deepseek_unaugmented", "deepseek_unaugmented_groundzero50")
