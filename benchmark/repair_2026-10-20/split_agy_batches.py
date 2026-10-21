"""Split the Gemini-CLI (Antigravity) counterion-correction batch outputs into one answer per question.

Each batch reply lists answers under a markdown heading naming the question (see HEAD). Every block is written
to answers_dg/<arm>/<QID>.txt in the same shape the browser arms use, so grade_repair.py --dg grades
all arms with the same code. A question missing from a batch is reported, never filled in.

Usage: python split_agy_batches.py <arm> <batch.txt> [<batch.txt> ...]   (arm: agy_unaug | agy_aug)
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
# headings vary by batch: "### Question H02", "## Question GZ16 Analysis", "### Solution for Question GZ16"
HEAD = re.compile(r"^#+[^\n]*?\bQuestion\s+((?:GZ|H|X)\d+)\b[^\n]*$", re.M)


def split(text: str) -> dict[str, str]:
    heads = list(HEAD.finditer(text))
    out = {}
    for i, m in enumerate(heads):
        end = heads[i + 1].start() if i + 1 < len(heads) else len(text)
        out[m.group(1)] = text[m.start():end].strip() + "\n"
    return out


def main() -> int:
    arm, files = sys.argv[1], sys.argv[2:]
    dest = HERE / "answers_dg" / arm
    dest.mkdir(parents=True, exist_ok=True)
    for f in files:
        blocks = split(Path(f).read_text(encoding="utf-8", errors="replace"))
        for qid, body in blocks.items():
            if "FINAL_JSON" not in body:
                print(f"{arm} {qid}: no FINAL_JSON in block")
            (dest / f"{qid}.txt").write_text(body, encoding="utf-8")
        print(f"{arm}: {Path(f).name} -> {sorted(blocks)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
