#!/usr/bin/env python3
"""Score Paper 3 test-kit answers. Reads ONLY the FINAL_JSON block of answers/<run_id>/answer.txt.
Numeric keys default to |ans - truth| <= tol (accept_abs: either sign passes); kind "interval"
passes on lo <= ans <= hi (bands generated from the release JSONs); kind "sign" passes on
sign only (sign "+" -> ans > 0). String keys: exact match, case-insensitive, whitespace-normalized.
No phrase matching, no guessing: absent/malformed block -> parsed=no and every key reported as
not-scored. Outputs per-run lines + summary.tsv."""
import json, re, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
def extract_final_json(text):
    idx = text.rfind("FINAL_JSON")
    if idx < 0: return None
    m = re.search(r"\{", text[idx:])
    if not m: return None
    start = idx + m.start(); depth = 0
    for i in range(start, len(text)):
        if text[i] == "{": depth += 1
        elif text[i] == "}":
            depth -= 1
            if depth == 0:
                block = text[start:i+1]
                # brief-28 #1: curly double quotes inside the FINAL_JSON block are plain "
                fixed = block.replace("\u201c", '"').replace("\u201d", '"')
                for cand in (block, fixed):
                    try: return json.loads(cand)
                    except Exception: pass
                return None
    return None
YESNO_TOKEN = re.compile(r"^\s*(yes|no|cannot_determine|cannot determine)(?=[\s,\-–—]|$)", re.I)
def score_run(ansobj, keys):
    out = {}
    for k, spec in keys.items():
        a = ansobj.get(k, "<<absent>>")
        if a == "<<absent>>": out[k] = ("absent", False); continue
        kind = spec.get("kind", "tol")
        numeric = isinstance(a, (int, float)) and not isinstance(a, bool)
        if kind == "interval":
            lo, hi = spec["lo"], spec["hi"]
            ok = numeric and lo <= a <= hi
            out[k] = (f"{a} vs interval [{lo}, {hi}]", ok)
        elif kind == "sign":
            want = spec.get("sign", "+")
            ok = numeric and (a > 0 if want == "+" else a < 0)
            out[k] = (f"{a} vs sign {want}", ok)
        elif numeric and isinstance(spec["value"], (int, float)):
            v = spec["value"]
            ok = abs(a - v) <= spec["tol"] or (spec.get("accept_abs") and abs(abs(a) - abs(v)) <= spec["tol"])
            out[k] = (f"{a} vs {v}+/-{spec['tol']}", ok)
        else:
            # brief-28 #2: JSON booleans are answers to yes/no keys; #3: a leading yes/no/
            # cannot_determine token (then dash/comma/space) is scored on that token only.
            if isinstance(a, bool): a = "yes" if a else "no"
            s = str(a).strip()
            mt = YESNO_TOKEN.match(s)
            if mt: s = mt.group(1)
            norm = lambda q: str(q).strip().lower().replace("_", " ")
            ok = norm(s) == norm(spec["value"])
            out[k] = (f"{s[:24]} vs {spec['value']}", ok)
    return out
def main(answers="answers", summary="summary.tsv"):
    adir = HERE / answers
    rows = ["case\tmodel\tstage\tkeys_correct/keys_total\tprinted_result_follows\tparsed"]
    lines = []
    for rp in sorted(adir.glob("*/answer.txt")) if adir.is_dir() else []:
        run_id = rp.parent.name
        parts = run_id.split("-")
        case, stage, model = parts[0], parts[1], parts[2] if len(parts) > 2 else "?"
        tj = HERE / "truth" / f"{case}.json"
        if not tj.exists(): lines.append(f"{run_id}: unknown case {case}"); continue
        truth = json.loads(tj.read_text())
        keys = truth["stages"][stage] if stage in truth["stages"] else None
        if keys is None: lines.append(f"{run_id}: unknown stage {stage}"); continue
        obj = extract_final_json(rp.read_text(encoding="utf-8", errors="replace"))
        if obj is None:
            lines.append(f"{run_id}: parsed=no (missing/malformed FINAL_JSON)")
            rows.append(f"{case}\t{model}\t{stage}\t0/{len(keys)}\tnot-parsed\tno"); continue
        res = score_run(obj, keys)
        good = sum(1 for _, ok in res.values() if ok)
        prf = "n/a"
        if "printed_result_follows" in keys:
            prf = "OK" if res["printed_result_follows"][1] else f"WRONG ({res['printed_result_follows'][0]})"
        for k, (det, ok) in res.items():
            lines.append(f"{run_id}: {k} {'PASS' if ok else 'FAIL'} [{det}]")
        rows.append(f"{case}\t{model}\t{stage}\t{good}/{len(keys)}\t{prf}\tyes")
    (HERE / summary).write_text("\n".join(rows) + "\n")
    print("\n".join(lines) if lines else "(no answer files yet - nothing to score; this is a valid state)")
    print(f"\nsummary -> {summary} ({len(rows)-1} rows)")
if __name__ == "__main__":
    main(*sys.argv[1:])
