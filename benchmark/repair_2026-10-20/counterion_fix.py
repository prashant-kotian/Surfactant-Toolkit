"""Correct the counterion-factor convention in SurfBench questions and golds (2026-10-20).

Standard pseudo-phase result (Zana, Langmuir 12 (1996) 1208):
    deltaG_mic = (1 + beta) R T ln X = (2 - alpha) R T ln X
beta = counterion BINDING degree (1 - S2/S1), alpha = IONIZATION degree (S2/S1).
Questions written with (2 - beta) for a binding beta, or (1 + alpha) for an ionization
alpha, used the inverse. This script, question by question:
  * rewrites the stated formula in every question-bank file and in the gold record;
  * rescales each deltaG gold by new_factor / old_factor, after checking that the stored
    gold really was built with the old factor (refuses otherwise);
  * recomputes the fields that depend on deltaG (H-49 triangle, X-41 comparison).
H-04 keeps its numbers: Fu et al. 2019's 0.268 is the slope ratio (ionization), so
(2 - 0.268) was right; only its "binding degree beta" label is wrong and is renamed.
Idempotent: re-running on corrected files changes nothing.

Usage:  python counterion_fix.py          (dry run, prints the plan)
        python counterion_fix.py --apply
"""
from __future__ import annotations

import json
import math
import re
import sys
from pathlib import Path

SETS = Path(__file__).resolve().parents[1] / "source_sets"
R = 8.314462618
DRY = "--apply" not in sys.argv


def rtlnx(cmc_mM: float, T: float, x_style: str) -> float:
    c = cmc_mM / 1000.0
    x = c / 55.5 if x_style == "div" else c / (c + 55.5)
    return R * T * math.log(x) / 1000.0


# question -> (bank, value given, kind, CMC mM, T, X convention, gold fields scaled with deltaG)
#   kind "beta": stated (2 - beta) -> (1 + beta);  kind "alpha": stated (1 + alpha) -> (2 - alpha);
#   kind "alpha_unstated": formula not in the text, only the gold is wrong.
PLAN = {
    "H-02": ("hard50", 0.65625, "beta", 5.75, 298.15, "div", ["delta_g_mic_kJ_per_mol"]),
    "H-08": ("hard50", 0.65625, "beta", 5.75, 298.15, "div", ["delta_g_mic_kJ_per_mol"]),
    "H-16": ("hard50", 0.72, "beta", 4.2, 298.15, "div", ["delta_g_mic_kJ_per_mol"]),
    "H-38": ("hard50", 0.65625, "beta", 5.75, 298.15, "div", ["delta_g_mic_kJ_per_mol"]),
    "H-47": ("hard50", 0.65625, "beta", 5.75, 298.15, "div", ["delta_g_mic_kJ_per_mol"]),
    "H-03": ("hard50", 0.31, "alpha", 4.242640687119285, 298.15, "div", ["delta_g_mic_kJ_per_mol"]),
    "H-19": ("hard50", 0.30, "alpha", 8.0, 298.15, "div", ["delta_g_mic_kJ_per_mol"]),
    "H-31": ("hard50", 0.30, "alpha", 6.0, 298.15, "div", ["delta_g_mic_kJ_per_mol"]),
    "H-46": ("hard50", 0.31, "alpha_unstated", 5.477225575051662, 298.15, "div", ["delta_g_mic_kJ_per_mol"]),
    "H-49": ("hard50", 0.272, "alpha", 8.0, 298.15, "div", ["delta_g_from_cmc_kJ_per_mol"]),
    "X-41": ("extrahard", 0.75, "beta", 4.40, 298.15, "sum", ["dG_correct"]),
    "X-47": ("extrahard", 0.6551540732, "beta", 4.574362416, 298.15, "sum", ["deltaG"]),
    "GZ-16": ("groundzero", 0.272, "alpha", 9.1815, 296.15, "sum", ["delta_g_mic_kJ_per_mol"]),  # Dankloff SDS at 23 C
    "GZ-20": ("groundzero", 0.39, "alpha", 2.5119, 298.15, "sum", ["delta_g_mic_kJ_per_mol"]),
}

BANKS = {"hard50": ["question_bank.txt", "question_bank_source.txt"],
         "groundzero": ["question_bank_augmented.txt", "question_bank_unaugmented.txt"],
         "extrahard": ["question_bank.txt"]}


def factors(v: float, kind: str) -> tuple[float, float]:
    if kind == "beta":
        return 2.0 - v, 1.0 + v
    return 1.0 + v, 2.0 - v            # alpha / alpha_unstated


def fix_text(text: str, kind: str) -> str:
    if kind == "beta":
        return re.sub(r"\(\s*2\s*-\s*beta\s*\)", "(1 + beta)", text)
    if kind == "alpha":
        return re.sub(r"\(\s*1\s*\+\s*alpha\s*\)", "(2 - alpha)", text)
    return text


def gold_value(gold: dict, bank: str, qid: str, field: str):
    e = gold[qid]
    return e["gold"][field] if bank == "groundzero" else e["fields"][field]["value"]


def set_gold(gold: dict, bank: str, qid: str, field: str, val):
    if bank == "groundzero":
        gold[qid]["gold"][field] = val
    else:
        gold[qid]["fields"][field]["value"] = val


def question_spans(text: str, bank: str) -> dict[str, tuple[int, int]]:
    """Character span of each question in a bank file."""
    spans = {}
    if bank == "hard50" and "END OF QUESTION H-" in text:      # paste-style bank
        start = 0
        for m in re.finditer(r"END OF QUESTION (H-\d\d)", text):
            spans[m.group(1)] = (start, m.end())
            start = m.end()
    else:                                                       # header-line banks
        pat = {"hard50": r"(?m)^(H-\d\d)\s*$", "groundzero": r"(?m)^(GZ-\d\d)\s*$"}.get(bank, r"(?m)^Question (X-\d\d)\b")
        ms = list(re.finditer(pat, text))
        for i, m in enumerate(ms):
            spans[m.group(1)] = (m.start(), ms[i + 1].start() if i + 1 < len(ms) else len(text))
    return spans


def main() -> int:
    golds = {b: json.loads((SETS / b / "gold_answers.json").read_text(encoding="utf-8")) for b in BANKS}
    texts = {(b, f): (SETS / b / f).read_text(encoding="utf-8") for b, fs in BANKS.items() for f in fs}
    report = []
    for qid, (bank, v, kind, cmc, T, xs, fields) in PLAN.items():
        old_f, new_f = factors(v, kind)
        gold = golds[bank]
        base = None
        for field in fields:
            old = gold_value(gold, bank, qid, field)
            if base is None and T is not None:
                base = rtlnx(cmc, T, xs)
            implied = old / base if base else None
            if implied is not None and abs(implied - new_f) < 1e-3:
                report.append(f"{qid} {field}: already corrected ({old:.4f})")
                continue
            if implied is not None and abs(implied - old_f) > 2e-3:
                raise SystemExit(f"{qid} {field}: gold {old} implies factor {implied:.4f}, expected old {old_f:.4f} "
                                 "-- not built the way this script assumes; stopping.")
            new = old * new_f / old_f
            set_gold(gold, bank, qid, field, new)
            report.append(f"{qid} {field}: {old:.4f} -> {new:.4f} kJ/mol  (factor {old_f:.4f} -> {new_f:.4f})")
        # dependent fields
        if qid == "H-49":
            f = gold[qid]["fields"]
            dg = f["delta_g_from_cmc_kJ_per_mol"]["value"]
            ds = (-5.0 - dg) * 1000.0 / 298.15                       # dH = -5.0 kJ/mol stated
            f["delta_s_J_per_mol_K"]["value"] = ds
            f["delta_g_recomputed_kJ_per_mol"]["value"] = -5.0 - 298.15 * ds / 1000.0
            report.append(f"H-49 dS -> {ds:.3f} J/(mol K); recomputed dG -> {f['delta_g_recomputed_kJ_per_mol']['value']:.4f}")
        if qid == "X-41":
            f = gold[qid]["fields"]
            dg = f["dG_correct"]["value"]
            draft = dg / (1.0 + v)                  # exact no-factor value, as xb_c5.build_x41 computes it
            f["rel_diff_pct"]["value"] = abs(draft - dg) / abs(dg) * 100.0
            assert f["claim_agrees"]["value"] is False and f["what_wrong"]["value"] == "formula"
            assert f["direction"]["value"] == "less_negative"
            report.append(f"X-41 rel_diff_pct -> {f['rel_diff_pct']['value']:.3f} (claim/what_wrong/direction unchanged)")
        # texts: every bank file, only inside this question's span; plus the gold record's own strings
        if kind != "alpha_unstated":
            for (b, fname), text in list(texts.items()):
                if b != bank:
                    continue
                span = question_spans(text, b).get(qid)
                if not span:
                    continue
                s, e = span
                texts[(b, fname)] = text[:s] + fix_text(text[s:e], kind) + text[e:]
            gold[qid] = json.loads(fix_text(json.dumps(gold[qid], ensure_ascii=False), kind))
    # H-04: label only (0.268 is Fu's slope ratio = ionization degree; numbers unchanged)
    for (b, fname), text in list(texts.items()):
        if b == "hard50":
            span = question_spans(text, b).get("H-04")
            if span:
                s, e = span
                seg = text[s:e].replace("counterion binding degree beta = 0.268",
                                        "degree of counterion ionization alpha = 0.268 (S2/S1)")
                seg = re.sub(r"\(\s*2\s*-\s*beta\s*\)", "(2 - alpha)", seg)
                texts[(b, fname)] = text[:s] + seg + text[e:]
    h04 = json.dumps(golds["hard50"]["H-04"], ensure_ascii=False)
    h04 = h04.replace("(2-beta)=1.732", "(2-alpha)=1.732").replace("(2-beta)", "(2-alpha)").replace(
        "beta=0.268 real", "alpha=0.268 real (Fu's slope ratio, written 'beta' in the paper)")
    golds["hard50"]["H-04"] = json.loads(h04)
    report.append("H-04: label 'binding degree beta' -> 'ionization alpha' (values unchanged)")

    for line in report:
        print(line)
    changed = [k for k, t in texts.items() if t != (SETS / k[0] / k[1]).read_text(encoding="utf-8")]
    print("\nfiles with text changes:", [f"{b}/{f}" for b, f in changed])
    if DRY:
        print("\nDRY RUN -- nothing written. Re-run with --apply.")
        return 0
    for (b, f), t in texts.items():
        (SETS / b / f).write_text(t, encoding="utf-8")
    sys.path.insert(0, str(Path(__file__).parent))
    from _json_style import write                # each file keeps its own committed layout
    for b, g in golds.items():
        write(SETS / b / "gold_answers.json", b, g)
    print("\nAPPLIED.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
