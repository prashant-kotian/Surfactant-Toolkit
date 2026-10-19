"""Grader for HARD50 unaugmented runs — implements DESIGN.md section 7 exactly.

Usage:
  python grade_hard50.py                 # grade all transcripts in UPLOADS
  python grade_hard50.py --json out.json # write detailed results (default: grade_hard50_results.json)

Field classification for the reason/arithmetic split (reported separately):
  ARITHMETIC (arithmetic + final answer) : type 'number', plus integer 'exact' (structural counts)
  REASONING (reasoning/decision)         : type 'exact' with str/bool value, keywords_any/all, null

Question passes only if ALL template fields pass (no partial credit).
"""
from __future__ import annotations
import argparse, json, re, sys, unicodedata
from pathlib import Path

HERE = Path(__file__).resolve().parent
UPLOADS = Path("/home/user/uploads")

GOLD = {k: v for k, v in json.load(open(HERE / "gold_answers.json")).items() if k.startswith("H-")}


def _load_meta() -> dict:
    """Extract the META dict literal from build_hard50.py without importing it
    (importing would pull in the full surfactantkit package)."""
    import ast
    src = (HERE / "build_hard50.py").read_text()
    tree = ast.parse(src)
    fn = None
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name == "M":
            fn = node
    assert fn is not None, "M() helper not found"
    # define M in a tiny namespace
    ns: dict = {}
    exec(compile(ast.Module(body=[fn], type_ignores=[]), "<M>", "exec"), ns)
    for node in tree.body:
        if isinstance(node, ast.Assign):
            for t in node.targets:
                if isinstance(t, ast.Name) and t.id == "META":
                    return eval(compile(ast.Expression(node.value), "<META>", "eval"), ns)
    raise RuntimeError("META not found")


META = _load_meta()

STOP = {"the", "and", "for", "with", "from", "that", "this", "not", "are", "was", "were",
        "has", "have", "had", "but", "its", "our", "your", "all", "any", "can", "may",
        "yes", "no", "to", "of", "in", "on", "at", "is", "be", "as", "it", "by", "or",
        "a", "an", "if", "then", "than", "true", "false", "null", "none", "not"}


# ---------------------------------------------------------------- extraction
def _blocks(text: str, pattern: str) -> dict[str, str]:
    parts = re.split(pattern, text, flags=re.M)
    out: dict[str, str] = {}
    for i in range(1, len(parts), 2):
        out[parts[i]] = parts[i + 1]   # later duplicate overwrites -> last copy wins
    return out


def parse_chatgpt(text):
    return _blocks(text, r"^QUESTION (H-\d+)[ \t]*$")


def parse_deepseek(text):
    return _blocks(text, r"^=+\n(H-\d+)\n=+[ \t]*$")


def parse_gemini(text):
    return _blocks(text, r"^\*\*Question (H-\d+)\*\*[ \t]*$")


def parse_qwen(text):
    # headers: "### H-07" and "### H-36 (chart read)"; file has a preamble copy of
    # H-01 before the first ### header (later copy wins -> the ### blocks are graded)
    return _blocks(text, r"^### (H-\d+)(?: \(.*?\))?[ \t]*$")


PARSERS = {"chatgpt": parse_chatgpt, "deepseek": parse_deepseek, "gemini": parse_gemini,
           "qwen": parse_qwen}

FILES = {
    "chatgpt": UPLOADS / "Chatgpt-unagumented H-01-to-H-50.txt",
    "deepseek": UPLOADS / "DeepSeek-unagumented-H-01-to-H-50.txt",
    "gemini": UPLOADS / "Gemini unagumented H-01-to-H-50.txt",
    "qwen": UPLOADS / "Qwen unagumented H-01-to-H-50.txt",
}


def _repair_json(s: str) -> str:
    s = re.sub(r"//.*?$", "", s, flags=re.M)
    s = s.replace("True", "true").replace("False", "false").replace("None", "null")
    s = re.sub(r"'", '"', s)
    s = re.sub(r",\s*([}\]])", r"\1", s)
    # unquoted keys
    s = re.sub(r'([{,]\s*)([A-Za-z_][A-Za-z0-9_]*)\s*:', r'\1"\2":', s)
    return s


def extract_final_json(block: str):
    """Return the LAST FINAL_JSON object in the block, or None."""
    idx = block.rfind("FINAL_JSON")
    if idx < 0:
        return None
    start = block.find("{", idx)
    if start < 0:
        return None
    depth, in_str, esc = 0, False, False
    for i in range(start, len(block)):
        c = block[i]
        if in_str:
            if esc:
                esc = False
            elif c == "\\":
                esc = True
            elif c == '"':
                in_str = False
            continue
        if c == '"':
            in_str = True
        elif c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0:
                raw = block[start:i + 1]
                for cand in (raw, _repair_json(raw)):
                    try:
                        obj = json.loads(cand)
                        if isinstance(obj, dict):
                            return obj
                    except Exception:
                        continue
                return {"__unparseable__": raw}
    return None


# ---------------------------------------------------------------- grading
def _norm(s: str) -> str:
    s = unicodedata.normalize("NFKC", str(s)).strip().lower()
    return re.sub(r"[\"'`]", "", s).strip()


def _tokens(s: str) -> list[str]:
    return [t for t in re.split(r"[^a-z0-9]+", s.lower()) if t and t not in STOP]


def grade_field(fd: dict, ans) -> bool:
    t = fd["type"]
    if t == "number":
        if isinstance(ans, bool) or ans is None:
            return False
        if isinstance(ans, str):
            try:
                ans = float(ans.strip().rstrip("%"))
            except Exception:
                return False
        if not isinstance(ans, (int, float)):
            return False
        g, tol = float(fd["value"]), fd.get("rel_tol", 0.01)
        return abs(float(ans) - g) <= tol * abs(g)
    if t == "null":
        return ans is None
    # exact / keywords
    if t in ("keywords_any", "keywords_all"):
        if ans is None:
            return False
        hay = _norm(ans)
        toks = fd["value"] if isinstance(fd["value"], list) else [fd["value"]]
        toks = [_norm(x) for x in toks]
        if t == "keywords_any":
            return any(k in hay for k in toks)
        return all(k in hay for k in toks)
    # exact
    g = fd["value"]
    if isinstance(g, bool):
        if isinstance(ans, str):
            return _norm(ans) in ("true", "false") and (ans.strip().lower() == "true") == g
        if isinstance(ans, bool):
            return ans == g
        return False
    if isinstance(g, int):
        if isinstance(ans, bool):
            return False
        try:
            return float(ans) == float(g)
        except Exception:
            return _norm(ans) == _norm(g)
    # string gold
    gold = _norm(g)
    if ans is None:
        return False
    a = _norm(ans)
    # short golds (single letter / roman / single short token): exact case-insensitive
    if len(gold) <= 4 or " " not in gold and len(gold) <= 12 and not re.search(r"\s", str(g).strip()):
        # single-token golds: full equality required (prompt prescribes vocabulary)
        if len(gold) <= 4:
            return a == gold
        # e.g. 'rodenas', 'exothermic', 'gemini', 'low': token must match a token of the answer
        return gold == a or gold in re.split(r"[^a-z0-9]+", a)
    gt = _tokens(gold)
    if not gt:
        return False
    at = set(re.split(r"[^a-z0-9]+", a))
    hit = sum(1 for tok in gt if any(tok in w for w in at))
    return hit / len(gt) >= 0.34


def field_bucket(fd: dict) -> str:
    if fd["type"] == "number":
        return "arithmetic"
    if fd["type"] == "exact" and isinstance(fd["value"], int) and not isinstance(fd["value"], bool):
        return "arithmetic"   # structural counts are counted results
    return "reasoning"


def grade_model(name: str, blocks: dict[str, str]):
    res = {}
    for qid, gold in sorted(GOLD.items()):
        block = blocks.get(qid)
        entry = {"found_block": block is not None, "fields": {}, "pass": False,
                 "axis": META.get(qid, {}).get("failure_axis"), "family": META.get(qid, {}).get("family")}
        if block is None:
            entry["error"] = "no block"
            res[qid] = entry
            continue
        ans = extract_final_json(block)
        if ans is None:
            entry["error"] = "no FINAL_JSON"
            res[qid] = entry
            continue
        if "__unparseable__" in ans:
            entry["error"] = "unparseable FINAL_JSON"
            entry["raw"] = ans["__unparseable__"][:300]
            res[qid] = entry
            continue
        entry["answer"] = ans
        all_ok = True
        for f, fd in gold["fields"].items():
            ok = grade_field(fd, ans.get(f, "__MISSING__"))
            entry["fields"][f] = {
                "ok": ok,
                "gold": fd.get("value"),
                "got": ans.get(f, "__MISSING__"),
                "bucket": field_bucket(fd),
                "type": fd["type"],
            }
            all_ok = all_ok and ok
        entry["pass"] = all_ok
        res[qid] = entry
    return res


def summarize(results: dict[str, dict]):
    out = {}
    for model, rq in results.items():
        qpass = [q for q, e in rq.items() if e["pass"]]
        n_q = len(rq)
        # field level
        stats = {"arithmetic": [0, 0], "reasoning": [0, 0]}   # [fail, total]
        axis_q = {}   # per question: did it fail an arith field / reason field
        for q, e in rq.items():
            fa = fr = False
            for f, d in e["fields"].items():
                stats[d["bucket"]][1] += 1
                if not d["ok"]:
                    stats[d["bucket"]][0] += 1
                    if d["bucket"] == "arithmetic":
                        fa = True
                    else:
                        fr = True
            axis_q[q] = (fa, fr, e["pass"], e.get("error"))
        out[model] = {
            "n": n_q,
            "passed": qpass,
            "pass_pct": 100.0 * len(qpass) / n_q,
            "arith_fail_pct": 100.0 * stats["arithmetic"][0] / max(1, stats["arithmetic"][1]),
            "reason_fail_pct": 100.0 * stats["reasoning"][0] / max(1, stats["reasoning"][1]),
            "arith_fields": stats["arithmetic"],
            "reason_fields": stats["reasoning"],
            "errors": [q for q, e in rq.items() if e.get("error")],
            "axis_q": axis_q,
        }
    return out


def cross_check(results: dict[str, dict]):
    """Flag fields where >=2 models agree with each other but all disagree with gold,
    and fields where NO model passes (candidates for gold review)."""
    suspects = []
    for qid, gold in GOLD.items():
        vals = {}
        for model, rq in results.items():
            e = rq.get(qid, {})
            if "fields" not in e:
                continue
            for f, d in e["fields"].items():
                vals.setdefault(f, []).append((model, d["got"], d["ok"], d["gold"], d["type"]))
        for f, lst in vals.items():
            failed = [x for x in lst if not x[2]]
            passed = [x for x in lst if x[2]]
            if not failed:
                continue
            # unanimous failure with unanimous same answer
            gv = failed[0][3]
            got_vals = [x[1] for x in failed]
            same = None
            if not passed and len(failed) >= 2:
                allnum = all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in got_vals)
                if allnum and isinstance(gv, (int, float)) and not isinstance(gv, bool) and float(gv) != 0:
                    nums = [float(v) for v in got_vals]
                    mean = sum(nums) / len(nums)
                    g = float(gv)
                    # mutually consistent cluster that disagrees with gold beyond gold's tolerance
                    if max(nums) - min(nums) <= 0.02 * abs(mean) and abs(mean - g) > abs(g) * 1e-9:
                        same = "numeric-cluster"
                elif allnum and isinstance(gv, (int, float)) and float(gv) == 0:
                    nums = [float(v) for v in got_vals]
                    if max(nums) - min(nums) <= 1e-9:
                        same = "numeric-cluster-zero-gold"
                else:
                    norms = {_norm(v) if not isinstance(v, (dict, list)) else str(v) for v in got_vals}
                    if len(norms) == 1 and "__MISSING__" not in norms and _norm(gv) not in norms:
                        same = "string-unanimous"
            if same:
                suspects.append({"qid": qid, "field": f, "gold": gv, "type": failed[0][4],
                                 "model_values": {m: v for m, v, *_ in failed},
                                 "passed_by": [x[0] for x in passed],
                                 "agreement": same})
    return suspects


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", default=str(HERE / "grade_hard50_results.json"))
    args = ap.parse_args()

    results = {}
    for model, path in FILES.items():
        if not path.exists():
            print(f"[skip] {model}: file missing", file=sys.stderr)
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        blocks = PARSERS[model](text)
        print(f"[parse] {model}: {len(blocks)}/50 question blocks")
        results[model] = grade_model(model, blocks)

    summary = summarize(results)
    suspects = cross_check(results)

    json.dump({"summary": {k: {kk: vv for kk, vv in v.items() if kk != 'axis_q'} for k, v in summary.items()},
               "results": results, "suspects": suspects},
              open(args.json, "w"), indent=1, default=str)

    print("\n================ SUMMARY ================")
    for model, s in summary.items():
        print(f"\n{model.upper()}:")
        print(f"  overall pass      : {len(s['passed'])}/{s['n']} = {s['pass_pct']:.1f}%")
        print(f"  passed questions  : {', '.join(s['passed']) or '(none)'}")
        print(f"  arithmetic fields : {s['arith_fields'][0]}/{s['arith_fields'][1]} failed = {s['arith_fail_pct']:.1f}%")
        print(f"  reasoning fields  : {s['reason_fields'][0]}/{s['reason_fields'][1]} failed = {s['reason_fail_pct']:.1f}%")
        if s["errors"]:
            print(f"  extraction errors : {s['errors']}")
    print(f"\n================ GOLD SUSPECTS ({len(suspects)}) ================")
    for s in suspects:
        print(f"  {s['qid']} {s['field']}: gold={s['gold']!r} agreement={s['agreement']} "
              f"values={s['model_values']} passed_by={s['passed_by']}")
    print(f"\ndetailed results -> {args.json}")


if __name__ == "__main__":
    main()
