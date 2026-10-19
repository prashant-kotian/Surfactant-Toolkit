"""Grade the Claude and DeepSeek MCP-augmented EXTRAHARD runs (2026-10-12),
one FINAL_JSON per X##.txt file, and compare against their already-graded
unaugmented pass rates.
"""
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).parent
AUG_DIR = HERE / "augmented_runs_2026-10-12"
UNAUG_RESULTS = HERE.parents[2] / "benchmark/transcripts/extrahard50_unaugmented" / "grade_results.json"  # repo-relative 2026-10-18

sys.path.insert(0, str(HERE))
import grade_x50 as G

FOLDERS = {
    "Claude": AUG_DIR / "claude",
    "DeepSeek": AUG_DIR / "deepseek",
}


def extract_final_json(text: str):
    idx = text.rfind("FINAL_JSON")
    if idx == -1:
        return None
    brace = text.find("{", idx)
    if brace == -1:
        return None
    depth = 0
    end = None
    for j in range(brace, len(text)):
        if text[j] == "{":
            depth += 1
        elif text[j] == "}":
            depth -= 1
            if depth == 0:
                end = j + 1
                break
    if end is None:
        return None
    try:
        return json.loads(text[brace:end])
    except json.JSONDecodeError:
        return None


def parse_folder(folder: Path) -> dict[str, dict]:
    answers = {}
    for qid_num in range(1, 51):
        qid = f"X-{qid_num:02d}"
        path = folder / f"X{qid_num:02d}.txt"
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8")
        parsed = extract_final_json(text)
        if parsed is not None:
            answers[qid] = parsed
    return answers


def main():
    unaug = json.loads(UNAUG_RESULTS.read_text(encoding="utf-8"))
    all_results = {}
    for model, folder in FOLDERS.items():
        answers = parse_folder(folder)
        missing_parse = [qid for qid in G.EXPECTED_IDS if qid not in answers]
        results = G.grade_all(answers)
        passed = [qid for qid, r in results.items() if r["pass"]]
        failed = [qid for qid, r in results.items() if not r["pass"]]
        all_results[model] = {"results": results, "passed": passed, "failed": failed, "missing_parse": missing_parse}
        print(f"=== {model} AUGMENTED: {len(passed)}/50 ({len(passed)*2}.0%) ===")
        if missing_parse:
            print(f"  could not parse FINAL_JSON for: {missing_parse}")
        for qid in failed:
            bad_fields = [k for k, v in results[qid]["fields"].items() if not v]
            print(f"  FAIL {qid}: bad fields = {bad_fields}")
        print()

    print("=== Family breakdown (G1..G10, X-01..50 in blocks of 5) ===")
    header = "Family".ljust(10) + "".join(m.ljust(22) for m in FOLDERS)
    print(header)
    for g in range(10):
        lo, hi = g * 5 + 1, g * 5 + 5
        fam_ids = [f"X-{i:02d}" for i in range(lo, hi + 1)]
        row = f"G{g+1}".ljust(10)
        for model in FOLDERS:
            p = sum(1 for qid in fam_ids if qid in all_results[model]["passed"])
            row += f"{p}/5".ljust(22)
        print(row)

    print()
    print("=== Augmented vs Unaugmented ===")
    for model in FOLDERS:
        n_aug = len(all_results[model]["passed"])
        n_unaug = len(unaug[model]["passed"])
        delta = n_aug - n_unaug
        sign = "+" if delta >= 0 else ""
        print(f"{model}: unaugmented {n_unaug}/50 ({n_unaug*2}.0%)  ->  augmented {n_aug}/50 ({n_aug*2}.0%)   [{sign}{delta}]")

    with open(AUG_DIR / "grade_results_augmented.json", "w", encoding="utf-8") as f:
        json.dump(
            {m: {"passed": v["passed"], "failed": v["failed"], "missing_parse": v["missing_parse"]}
             for m, v in all_results.items()},
            f, indent=1,
        )


if __name__ == "__main__":
    main()
