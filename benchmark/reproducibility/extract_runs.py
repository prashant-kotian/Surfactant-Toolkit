#!/usr/bin/env python3
"""extract_runs.py -- pull every FINAL_JSON answer for the reproducibility set (GZ-51..68, five
repeats, unaugmented only) out of the raw transcripts into one table, runs.json.

Each question has a unique set of answer keys, so an answer is assigned to its question by its
keys, not by position in the file; templates echoed from the prompt ("<number>") are skipped.
For each answer the text the model wrote before it is kept for the qualitative checks.

    python extract_runs.py
"""
from __future__ import annotations

import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
GZ = HERE / "transcripts"
NEW = GZ / "new output" / "New folder" / "reproducibility type"

# question id <- key that only that question's answer block contains
SIGNATURE = [
    ("GZ-68", "axial_ratio"), ("GZ-66", "hlb_davies"), ("GZ-67", "hld_value"),
    ("GZ-51", "delta_h_kJ_per_mol"), ("GZ-52", "selected_mixture_model"),
    ("GZ-53", "delta_s_scenario_a_J_per_mol_K"), ("GZ-55", "physically_meaningful_interaction"),
    ("GZ-56", "charge_type_at_ph5"), ("GZ-57", "delta_g_mic_kJ_per_mol"),
    ("GZ-58", "f1_at_x1_0.15"), ("GZ-59", "tail_molar_volume_contribution_cm3_per_mol"),
    ("GZ-60", "estimated_cmc_mM"), ("GZ-61", "kappa_a"), ("GZ-62", "counterion_binding_degree"),
    ("GZ-63", "aggregation_number"), ("GZ-64", "best_estimate_cmc_mM"), ("GZ-65", "a_min_nm2"),
    ("GZ-54", "cpp"),  # last: GZ-66 and GZ-68 also carry cpp
]

# (model, protocol, attempt, file): batched = one chat answers a block of questions per attempt;
# fresh = every prompt in its own new session (attempt read from the [GZ-xx-rN] header)
SOURCES = [
    ("ChatGPT", "fresh", None, GZ / "unaugmented_runs_20260921" / "chatgpt_GZ51-68_repro_x5_unaugmented.txt"),
    ("Gemini", "fresh", None, GZ / "unaugmented_runs_20260921" / "gemini_GZ51-68_repro_x5_unaugmented.txt"),
    ("Gemini (batched)", "batched", 1, NEW / "gemini" / "GZ-51_to_68_Outputs-gemini-1st-attenpt.txt"),
    ("Gemini (batched)", "batched", 2, NEW / "gemini" / "gemini-gz-51-68-2nd-attempt.txt"),
    ("Gemini (batched)", "batched", 3, NEW / "gemini" / "gemini-GZ-69-78-rd--3rs-atempt.txt"),
    ("Gemini (batched)", "batched", 4, NEW / "gemini" / "gemini-GZ-51-68-4th-attempt.txt"),
    ("Gemini (batched)", "batched", 5, NEW / "gemini" / "gemini gz-51-68-5th attempt.txt"),
    ("Qwen", "batched", 1, NEW / "qween" / "Qwen_51-68-1st attempt.txt"),
    ("Qwen", "batched", 2, NEW / "qween" / "Qwen_GZ-51-68-2ND-ATTEMPT.txt"),
    ("Qwen", "batched", 3, NEW / "qween" / "Qwen_GZ-51-68-3rd-attempt.txt"),
    ("Qwen", "batched", 4, NEW / "qween" / "qween-gz-51-68output-4th-attempt.txt"),
    ("Qwen", "batched", 5, NEW / "qween" / "Qwen_gz51-68-5th attempt.txt"),
]
_C1 = GZ / "unagumented 18th sept" / "question repeat type"
_C2 = GZ / "unaggumented output 19th sept"
for n, f51, f61, f66 in (
        (1, "GZ_51-60_outputs 1st attempt.txt", "GZ_61-65_solutions 1st.txt", "SDS_analysis_GZ66-68 1st.txt"),
        (2, "GZ51-60_surfactant_answers 2nd attempt.txt", "surfactant_solutions 2nd.txt", "SDS_analysis_all_questions 2nd.txt"),
        (3, "GZ_51-60_solutions 3rd attempt.txt", "GZ_surfactant_solutions 3rd.txt", "SDS_analysis_GZ66-67-68 3rd.txt"),
        (4, "GZ_51-60_answers 4th attempt.txt", "GZ_61-65_solutions 4th.txt", "SDS_GZ66-68_compiled 4th.txt"),
        (5, "GZ_51-60_outputs_1 5th attempt.txt", "GZ_61-65_solutions 5th.txt", "SDS_analysis_GZ66-67-68 5th.txt")):
    SOURCES += [("Claude", "batched", n, _C1 / f51), ("Claude", "batched", n, _C2 / f61),
                ("Claude", "batched", n, _C2 / "New folder" / f66)]


def qid_of(obj: dict) -> str | None:
    for q, k in SIGNATURE:
        if k in obj:
            return q
    return None


def objects(text: str):
    """Yield (start, end, dict) for every top-level {...} that parses as JSON."""
    i = 0
    while (i := text.find("{", i)) != -1:
        depth, j = 0, i
        while j < len(text):
            if text[j] == "{":
                depth += 1
            elif text[j] == "}":
                depth -= 1
                if depth == 0:
                    break
            j += 1
        raw = text[i:j + 1]
        try:  # strict=False: answers that wrap a long string over several lines
            obj = json.loads(raw.replace("\u201c", '"').replace("\u201d", '"'), strict=False)
        except ValueError:
            obj = None
        if isinstance(obj, dict) and not any(isinstance(v, str) and v.startswith("<") for v in obj.values()):
            yield i, j + 1, obj
            i = j + 1
        else:
            i += 1


KV = re.compile(r"^[ \t]*([a-z][a-z0-9_.]*)[ \t]*=[ \t]*(-?[\d.]+(?:e-?\d+)?)[ \t]*$", re.M | re.I)


def kv_objects(text: str):
    """Fallback for answers written as 'key = value' lines (one Claude attempt): consecutive
    lines that name answer keys form one answer."""
    keys = {k for _, k in SIGNATURE} | {"zeta_potential_mV", "cmc_mM", "a_min_nm2",
                                        "gamma_max_mol_per_m2", "aggregation_number"}
    cur, start, last = {}, 0, -1
    for m in KV.finditer(text):
        if m.group(1) not in keys:
            continue
        if cur and text.count("\n", last, m.start()) > 3:
            yield start, last, cur
            cur = {}
        if not cur:
            start = m.start()
        cur[m.group(1)] = float(m.group(2))
        last = m.end()
    if cur:
        yield start, last, cur


def main() -> int:
    rows, log = [], []
    for model, proto, att, path in SOURCES:
        text = path.read_text(encoding="utf-8", errors="replace")
        if proto == "fresh":
            parts = re.split(r"PROMPT \d+ \[(GZ-\d\d)-r(\d)\]", text)
            segs = [(parts[k], int(parts[k + 1]), parts[k + 2]) for k in range(1, len(parts), 3)]
        else:
            segs = [(None, att, text)]
        for want, run, seg in segs:
            body = seg.split("RESPONSE", 1)[-1] if proto == "fresh" else seg
            prev = 0
            found = []
            objs = list(objects(body))
            if proto == "batched" and not any(qid_of(o) for _, _, o in objs):
                objs = list(kv_objects(body))
            for s, e, obj in objs:
                q = qid_of(obj)
                if q is None or (want and q != want):
                    continue
                found.append((q, body[prev:s], obj))
                prev = e
            if proto == "fresh":
                found = found[-1:]  # the model's own final block in that session
                if not found:
                    log.append(f"{model} {want}-r{run}: no answer block")
            for q, before, obj in found:
                rows.append({"model": model, "protocol": proto, "qid": q, "run": run,
                             "answer": obj, "text": before[-6000:], "file": path.name})
    # one answer per (model, qid, run): the last one in the file wins (restated blocks)
    best = {}
    for r in rows:
        best[(r["model"], r["qid"], r["run"])] = r
    out = sorted(best.values(), key=lambda r: (r["model"], r["qid"], r["run"]))
    (HERE / "runs.json").write_text(json.dumps(out, indent=1, ensure_ascii=True), encoding="utf-8")
    models = sorted({r["model"] for r in out})
    print("answers per model:", {m: sum(1 for r in out if r["model"] == m) for m in models})
    for m in models:
        miss = [f"{q}-r{n}" for q, _ in sorted(SIGNATURE) for n in range(1, 6)
                if (m, q, n) not in best]
        print(f"{m}: missing {len(miss)}", miss[:30])
    for l in log:
        print("LOG", l)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
