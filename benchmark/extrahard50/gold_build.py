"""EXTRA50 runner: builds all 5 chunks (X-01..X-50), validates, writes
gold_answers.json (grading source of truth) + question_bank.txt (prompt file).

    PYTHONPATH=../../../src:. python gold_build.py
"""
from __future__ import annotations
import json
import re
from pathlib import Path

import xb_c1, xb_c2, xb_c3, xb_c4, xb_c5
from xb_common import validate, check_near_duplicates

HERE = Path(__file__).resolve().parent
EXPECTED = [f"X-{i:02d}" for i in range(1, 51)]


def main():
    Q = {}
    for m in (xb_c1, xb_c2, xb_c3, xb_c4, xb_c5):
        got = m.build()
        dup = set(got) & set(Q)
        assert not dup, f"question ids collide across chunks: {dup}"
        Q.update(got)

    errs = validate(Q)
    assert not errs, "\n".join(errs)
    assert sorted(Q) == EXPECTED, f"id set mismatch: {sorted(set(EXPECTED) - set(Q))}"

    # ---- gold_answers.json (records include prompts; grading source of truth) ----
    recs = {k: Q[k] for k in EXPECTED}
    (HERE / "gold_answers.json").write_text(json.dumps(recs, indent=1, sort_keys=True), encoding="utf-8")

    # ---- questions file (prompt-only, like HARD50/GZ txt) ----
    lines = ["=" * 76, "EXTRA50 - Question Bank (unaugmented LLM evaluation)", "=" * 76, ""]
    for k in EXPECTED:
        q = recs[k]
        lines += [
            "-" * 76,
            f"## {k} | family {q['family']} | arith {q['axes']['arith']} | "
            f"reason {q['axes']['reason']}",
            f"META levers={','.join(q['levers'])}  (tools/mechanism in gold_answers.json only -- never put tool names in a paste artifact)",
            "-" * 76,
            q["body"].rstrip(),
            "",
        ]
    (HERE / "question_bank.txt").write_text("\n".join(lines), encoding="utf-8")

    # ---- audits ----
    bodies = {k: recs[k]["body"] for k in EXPECTED}
    intra = check_near_duplicates(bodies)
    print("intra-bank near-duplicates:", intra or "none")

    h50 = HERE.parent / "hard50" / "question_bank.txt"
    if h50.exists():
        hb = {}
        for blk in re.split(r"^={16}$", h50.read_text(), flags=re.M):
            m = re.search(r"END OF QUESTION (H-\d\d)", blk)
            if m and len(blk) > 400:
                hb[m.group(1)] = blk
        assert len(hb) == 50, f"HARD50 parse got {len(hb)} blocks"
        sides = check_near_duplicates({**bodies, **hb})
        cross = []
        for c in sides:
            a, b = c.split(":")[0].split("/")
            if {a[:2], b[:2]} == {"X-", "H-"}:
                cross.append(c)
        print("vs-HARD50 cross near-duplicates:", cross or "none")

    nfields = sum(len(q["fields"]) for q in recs.values())
    nint = sum(1 for q in recs.values() for f in q["fields"].values() if f["type"] == "int")
    nnum = sum(1 for q in recs.values() for f in q["fields"].values() if f["type"] == "num")
    print(f"OK: {len(recs)} questions, {nfields} gold fields "
          f"({nnum} numeric, {nint} int, {nfields - nnum - nint} str/bool)")
    print("wrote gold_answers.json, question_bank.txt")


if __name__ == "__main__":
    main()
