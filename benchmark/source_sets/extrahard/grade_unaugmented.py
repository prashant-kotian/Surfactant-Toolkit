"""Grade the 5 unaugmented EXTRAHARD runs (manual + overnight-automated collection)."""
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).parent
OUTDIR = HERE.parents[2] / "benchmark/transcripts/extrahard50_unaugmented"  # repo-relative since 2026-10-18

sys.path.insert(0, str(HERE))
import grade_x50 as G

FILES = {
    "ChatGPT": "chatgpt_unaugmented_transcript.txt",
    "Claude": "claude_unaugmented_transcript.txt",
    "DeepSeek": "deepseek_unaugmented_transcript.txt",
    "Qwen": "qwen_unaugmented_transcript.txt",
    "Gemini": "gemini_unaugmented_transcript.txt",
}


def parse_file(path: Path) -> dict[str, dict]:
    text = path.read_text(encoding="utf-8")
    # strict: "X-NN" with a real hyphen, standalone (not a math variable like "X=0.0044").
    # the one known typo ("X=08" in the ChatGPT file) is matched narrowly by exact value.
    headers = [(m.start(), m.end(), int(m.group(1))) for m in re.finditer(r"^X-(\d{2})\b", text, re.M | re.I)]
    headers += [(m.start(), m.end(), 8) for m in re.finditer(r"^X=08\b", text, re.M)]
    headers.sort(key=lambda h: h[0])
    answers = {}
    for i, (s, e, num) in enumerate(headers):
        nexts = headers[i + 1][0] if i + 1 < len(headers) else len(text)
        sec = text[e:nexts]
        idx = sec.rfind("FINAL_JSON")
        if idx == -1:
            continue
        brace = sec.find("{", idx)
        if brace == -1:
            continue
        depth = 0
        end = None
        for j in range(brace, len(sec)):
            if sec[j] == "{":
                depth += 1
            elif sec[j] == "}":
                depth -= 1
                if depth == 0:
                    end = j + 1
                    break
        if end is None:
            continue
        raw_json = sec[brace:end]
        try:
            parsed = json.loads(raw_json)
        except json.JSONDecodeError:
            continue
        qid = f"X-{num:02d}"
        # keep the LAST valid parse per qid (handles duplicate headers)
        answers[qid] = parsed
    return answers


def main():
    all_results = {}
    for model, fname in FILES.items():
        path = OUTDIR / fname
        answers = parse_file(path)
        missing_parse = [qid for qid in G.EXPECTED_IDS if qid not in answers]
        results = G.grade_all(answers)
        passed = [qid for qid, r in results.items() if r["pass"]]
        failed = [qid for qid, r in results.items() if not r["pass"]]
        all_results[model] = {
            "results": results,
            "passed": passed,
            "failed": failed,
            "missing_parse": missing_parse,
        }
        print(f"=== {model}: {len(passed)}/50 ===")
        if missing_parse:
            print(f"  could not parse FINAL_JSON for: {missing_parse}")
        for qid in failed:
            bad_fields = [k for k, v in results[qid]["fields"].items() if not v]
            print(f"  FAIL {qid}: bad fields = {bad_fields}")
        print()

    # family breakdown
    print("=== Family breakdown (G1..G10, X-01..50 in blocks of 5) ===")
    header = "Family".ljust(10) + "".join(m.ljust(10) for m in FILES)
    print(header)
    for g in range(10):
        lo, hi = g * 5 + 1, g * 5 + 5
        fam_ids = [f"X-{i:02d}" for i in range(lo, hi + 1)]
        row = f"G{g+1}".ljust(10)
        for model in FILES:
            p = sum(1 for qid in fam_ids if qid in all_results[model]["passed"])
            row += f"{p}/5".ljust(10)
        print(row)

    print()
    print("=== Overall ===")
    for model in FILES:
        n = len(all_results[model]["passed"])
        print(f"{model}: {n}/50 = {n*2}.0%")

    with open(OUTDIR / "grade_results.json", "w", encoding="utf-8") as f:
        json.dump(
            {m: {"passed": v["passed"], "failed": v["failed"], "missing_parse": v["missing_parse"]}
             for m, v in all_results.items()},
            f, indent=1,
        )


if __name__ == "__main__":
    main()
