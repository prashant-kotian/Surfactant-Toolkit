"""Real API run of GZ-01..50 (unaugmented) against Gemini.

Resumable per-question-file design (this project's own established pattern,
survived 3 separate mid-run kills on the original 42-question set) -- re-
running this script skips any question that already has a saved output file.

Qwen/DashScope is NOT run here: confirmed real, current blocker (2026-09-19),
not a transient error -- "The free quota has been exhausted. To continue
accessing the model on a paid basis, please complete your payment
information..." (same real, already-documented blocker from 2026-09-14/15).
"""
import base64
import json
import os
import re
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "harness"))
from dotenv import load_dotenv
load_dotenv(os.path.join(os.path.dirname(__file__), "..", "harness", ".env"))

from providers import run
from prompts import SYSTEM_BASE

HERE = os.path.dirname(__file__)
TXT_PATH = os.path.join(HERE, "question_bank_unaugmented.txt")
OUT_DIR = os.path.join(HERE, "gemini_cli_2026-09-19")
os.makedirs(OUT_DIR, exist_ok=True)

VISION_IMAGES = {
    "GZ-15": "GZ-15_chart.png",
    "GZ-36": "GZ-36_chart.png",
    "GZ-37": "GZ-37_chart.png",
    "GZ-38": "GZ-38_chart.png",
    "GZ-39": "GZ-39_chart.png",
    "GZ-40": "GZ-40_chart.png",
}

MODEL = "gemini-3.1-pro-preview"


def parse_questions(text: str) -> dict[str, str]:
    """Extract {GZ-NN: paste_block_text} for GZ-01..GZ-50 only."""
    pattern = re.compile(
        r"^(GZ-\d+)\s*$\n(?:\[.*?\]\s*\n)?\s*>>> PASTE BELOW >>>\n(.*?)\n<<< PASTE ABOVE <<<",
        re.S | re.M,
    )
    out = {}
    for m in pattern.finditer(text):
        qid, body = m.group(1), m.group(2).strip()
        num = int(qid.split("-")[1])
        if 1 <= num <= 50:
            out[qid] = body
    return out


def main():
    with open(TXT_PATH, encoding="utf-8") as f:
        text = f.read()
    questions = parse_questions(text)
    print(f"Parsed {len(questions)} questions (expected 50).")

    ordered_ids = sorted(questions.keys(), key=lambda k: int(k.split("-")[1]))
    for qid in ordered_ids:
        out_path = os.path.join(OUT_DIR, f"{qid}.json")
        if os.path.exists(out_path):
            print(f"{qid}: already done, skipping")
            continue

        user_text = questions[qid]
        kw = {}
        if qid in VISION_IMAGES:
            img_path = os.path.join(HERE, VISION_IMAGES[qid])
            with open(img_path, "rb") as imf:
                kw["image_base64"] = base64.b64encode(imf.read()).decode()

        print(f"{qid}: calling Gemini...", flush=True)
        t0 = time.time()
        try:
            result = run("gemini", MODEL, SYSTEM_BASE, user_text, "unaugmented", **kw)
        except Exception as e:
            print(f"{qid}: FAILED with {e!r} -- stopping run, rerun this script to resume from here", flush=True)
            break
        elapsed = time.time() - t0

        record = {
            "question_id": qid,
            "model": MODEL,
            "elapsed_s": round(elapsed, 1),
            "error": result.get("error"),
            "final_text": result.get("final_text"),
            "raw": result.get("raw"),
        }
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(record, f, indent=2, ensure_ascii=False)

        status = "ERROR: " + str(result.get("error")) if result.get("error") else "ok"
        print(f"{qid}: {status} ({elapsed:.1f}s)", flush=True)
        time.sleep(8)  # real rate-limit hit at GZ-33 on 2026-09-19 (HTTP 429) -- pace requests

    print("\nDone. Compiling combined output...")
    combined = []
    for qid in ordered_ids:
        out_path = os.path.join(OUT_DIR, f"{qid}.json")
        if os.path.exists(out_path):
            with open(out_path, encoding="utf-8") as f:
                combined.append(json.load(f))
    combined_path = os.path.join(HERE, "gemini_cli_answers_combined.json")
    with open(combined_path, "w", encoding="utf-8") as f:
        json.dump(combined, f, indent=2, ensure_ascii=False)
    print(f"Wrote {combined_path} ({len(combined)} questions).")


if __name__ == "__main__":
    main()
