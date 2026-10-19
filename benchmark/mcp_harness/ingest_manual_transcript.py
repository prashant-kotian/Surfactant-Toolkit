"""Parse a manually-collected transcript (see
PhD-Research/SurfBench_Manual_Runs/questions_export.txt Step 3 format) back
into the same {models, conditions, subset_ids, rows} schema run_pilot.py's
own results files use, so grading.py can be run on manually-collected
unaugmented answers exactly as it would on API-collected ones.

Transcript format expected (one .txt file per model):
    ### <question_id>
    <raw response text, including the FINAL_JSON: line>

    ### <next question_id>
    ...

Usage:
    python ingest_manual_transcript.py --model gpt-5.6-sol --transcript gpt_manual.txt --out results/manual_gpt_unaug.json
"""
from __future__ import annotations
import argparse
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from extract import extract_final
from grading import grade

_ID_RE = re.compile(r"^###\s+(\S+)\s*$", re.MULTILINE)


def parse_transcript(text: str) -> dict[str, str]:
    """Split on '### <id>' headers, return {id: response_text}."""
    matches = list(_ID_RE.finditer(text))
    if not matches:
        raise ValueError("No '### <id>' headers found -- check the transcript follows the "
                         "Step 3 format exactly (### on its own line, then the response text).")
    out = {}
    for i, m in enumerate(matches):
        qid = m.group(1)
        start = m.end()
        end = matches[i + 1].start() if i + 1 < len(matches) else len(text)
        response = text[start:end].strip()
        if qid in out:
            raise ValueError(f"Duplicate question id '{qid}' found in transcript -- "
                             f"each id should appear exactly once.")
        out[qid] = response
    return out


def load_question_bank() -> dict[str, dict]:
    bench = Path(__file__).resolve().parents[1]
    t1 = json.loads((bench / "question_bank_tier1.json").read_text())
    t2 = json.loads((bench / "question_bank_tier2.json").read_text())
    return {q["id"]: q for q in (t1 + t2)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True, help="logical model name, e.g. gpt-5.6-sol")
    ap.add_argument("--transcript", required=True, help="path to the manually-collected .txt transcript")
    ap.add_argument("--out", required=True, help="output path for the results JSON, e.g. results/manual_gpt_unaug.json")
    ap.add_argument("--condition", default="unaugmented")
    args = ap.parse_args()

    text = Path(args.transcript).read_text(encoding="utf-8")
    answers = parse_transcript(text)
    bank = load_question_bank()

    unknown = [qid for qid in answers if qid not in bank]
    if unknown:
        print(f"WARNING: {len(unknown)} question id(s) in the transcript don't match the "
              f"real question bank (typo? wrong id?): {unknown[:10]}{'...' if len(unknown) > 10 else ''}")

    rows = []
    n_correct = 0
    for qid, response_text in answers.items():
        q = bank.get(qid)
        if q is None:
            continue
        ans = extract_final(response_text)
        g = grade(q["gold_answer"], q.get("tolerance", {}), q.get("grading_method", "numeric_tolerance"),
                  ans, response_text, q.get("trap_type", "none"))
        n_correct += int(g["correct"])
        rows.append({
            "id": qid,
            "category": q["category"],
            "subcategory": q["subcategory"],
            "trap_type": q["trap_type"],
            "difficulty": q["difficulty"],
            "tools_required": q["tools_required"],
            "model": args.model,
            "condition": args.condition,
            "correct": g["correct"],
            "grade_reason": g["reason"],
            "extracted_answer": ans,
            "n_tool_calls": 0,
            "tools_used": [],
            "thinking_chars": 0,
            "final_text": response_text,
            "seconds": None,
            "error": None,
            "raw": {"source": "manual_transcript"},
        })

    out_data = {
        "models": [args.model],
        "conditions": [args.condition],
        "subset_ids": sorted(answers.keys()),
        "rows": rows,
    }
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(out_data, indent=2))
    print(f"Parsed and graded {len(rows)}/{len(answers)} answers into {args.out}.")
    print(f"Coverage: {len(rows)}/{len(bank)} of the full question bank.")
    if rows:
        print(f"Accuracy on this batch: {n_correct}/{len(rows)} = {100*n_correct/len(rows):.1f}%")


if __name__ == "__main__":
    main()
