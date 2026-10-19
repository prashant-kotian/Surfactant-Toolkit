"""Grade ChatGPT/Claude/Gemini unaugmented GZ-01..50 answers against gold_answers.json.

Real fix vs the first attempt: the fields a model is actually expected to report
are defined by each question's OWN FINAL_JSON template (in question_bank_augmented.txt),
not by every key in gold['gold'] -- that dict mixes real answer-fields with
grading-instruction prose ('grading_note', 'reasoning') and reference/context
data the model was never asked to echo back (e.g. a literature value handed to
it in the question). Only template-field keys are graded.

Numeric fields: tolerance-based. Short string fields with a clean gold value:
lenient core-keyword match. Fields whose only gold information is a prose
grading_note (no clean string to match) are NOT auto-graded -- they're flagged
NEEDS_JUDGMENT and printed for a human/LLM read, since that's a real semantic
judgment call, not something keyword matching can honestly do.
"""
import json
import os
import re
import sys

BENCH = os.path.dirname(os.path.abspath(__file__))  # portable since 2026-10-18 (was a hardcoded Windows dev path)
GOLD_PATH = os.path.join(BENCH, "gold_answers.json")
TEMPLATE_PATH = os.path.join(BENCH, "question_bank_augmented.txt")

SOURCES = {
    "chatgpt": os.path.join(BENCH, "unaugmented_2026-09-21", "chatgpt_GZ01-50_unaugmented.txt"),
    "gemini": os.path.join(BENCH, "unaugmented_2026-09-21", "gemini_GZ01-50_unaugmented.txt"),
    "claude": os.path.join(BENCH, "claude_unaugmented", "answers.txt"),
}

META_KEYS = {"grading_note", "reasoning", "note"}


def extract_last_json_block(text: str):
    idx = text.rfind("FINAL_JSON")
    if idx == -1:
        return None
    tail = text[idx:]
    start = tail.find("{")
    if start == -1:
        return None
    depth = 0
    for i, ch in enumerate(tail[start:], start):
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
            if depth == 0:
                blob = tail[start:i + 1]
                try:
                    return json.loads(blob)
                except Exception:
                    cleaned = re.sub(r"```(?:json)?", "", blob)
                    try:
                        return json.loads(cleaned)
                    except Exception:
                        return None
    return None


def extract_templates() -> dict:
    """Question ID -> set of keys in its own FINAL_JSON template."""
    with open(TEMPLATE_PATH, "r", encoding="utf-8-sig") as f:
        text = f.read()
    id_re = re.compile(r"^(GZ-\d\d)$", re.MULTILINE)
    matches = list(id_re.finditer(text))
    templates = {}
    for i, m in enumerate(matches):
        label = m.group(1)
        start = m.start()
        end = matches[i + 1].start() if i + 1 < len(matches) else len(text)
        block = text[start:end]
        tm = re.search(r"FINAL_JSON:\s*\n\{(.*?)\n\}", block, re.S)
        if not tm:
            templates[label] = set()
            continue
        keys = re.findall(r'"([a-zA-Z0-9_.]+)"\s*:', tm.group(1))
        templates[label] = set(keys)
    return templates


def parse_writer_format(path: str) -> dict:
    with open(path, "r", encoding="utf-8-sig", errors="replace") as f:
        text = f.read()
    blocks = re.split(r"={10,}\nPROMPT \d+ \[(GZ-\d\d)\]\n={10,}\n", text)
    out = {}
    for i in range(1, len(blocks), 2):
        qid = blocks[i]
        content = blocks[i + 1] if i + 1 < len(blocks) else ""
        resp_idx = content.find("RESPONSE\n")
        response_text = content[resp_idx:] if resp_idx != -1 else content
        out[qid] = extract_last_json_block(response_text)
    return out


def parse_compiled_format(path: str) -> dict:
    with open(path, "r", encoding="utf-8-sig", errors="replace") as f:
        text = f.read()
    parts = re.split(r"[-=]{10,}\n(GZ-\d\d)\s", text)
    out = {}
    for i in range(1, len(parts), 2):
        qid = parts[i]
        content = parts[i + 1] if i + 1 < len(parts) else ""
        out[qid] = extract_last_json_block(content)
    return out


def is_close(a, b, rel=0.20, abs_tol=1e-6) -> bool:
    try:
        af, bf = float(a), float(b)
    except (TypeError, ValueError):
        return False
    if abs(af - bf) <= abs_tol:
        return True
    denom = max(abs(af), abs(bf), 1e-9)
    return abs(af - bf) / denom <= rel


_STOPWORDS = {"the", "a", "an", "of", "for", "and", "or", "with", "no", "not",
              "over", "on", "to", "by", "at", "this", "that", "is", "are"}


def string_matches(model_val, gold_val) -> bool:
    """Lenient: PASS if the model's text contains a meaningful fraction of
    gold's distinguishing (non-stopword) tokens -- catches paraphrase, doesn't
    require exact wording."""
    if model_val is None:
        return False
    m = re.sub(r"[^a-z0-9 ]", " ", str(model_val).lower())
    g = re.sub(r"[^a-z0-9 ]", " ", str(gold_val).lower())
    g_tokens = [t for t in g.split() if t and t not in _STOPWORDS and len(t) > 2]
    if not g_tokens:
        return str(model_val).strip().lower() == str(gold_val).strip().lower()
    hits = sum(1 for t in g_tokens if t in m)
    return hits / len(g_tokens) >= 0.34  # lowered from 0.5 -- real paraphrases rarely hit half


def grade_field(model_val, gold_val):
    if gold_val is None:
        return (model_val is None), "null-check"
    if isinstance(gold_val, bool):
        return (str(model_val).strip().lower() == str(gold_val).strip().lower()), "bool"
    if isinstance(gold_val, (int, float)):
        return is_close(model_val, gold_val), "numeric"
    if isinstance(gold_val, str):
        if len(gold_val) > 80:  # long prose gold value -- not a clean string to match
            return None, "needs-judgment"
        return string_matches(model_val, gold_val), "string"
    if isinstance(gold_val, dict):
        return True, "dict-skipped"
    return False, "unknown"


# Resolutions for fields whose gold is long free-text prose (no clean string
# to auto-match) -- same discipline as the Antigravity graders: one rationale
# per field, nothing left unscored. Ported from extra-surfactantkit-work
# commit 03b7c3a after the full Claude GZ-01..50 transcript was located.
JUDGMENT_RESOLUTIONS = {
    ("chatgpt", "GZ-04", "charge_confidence"): (True, "hedged but correctly ties confidence to the real pH-dependence gold is testing for"),
    ("gemini", "GZ-04", "charge_confidence"): (False, "field is null -- no answer given"),
    ("claude", "GZ-04", "charge_confidence"): (True, "'medium' is a synonym of gold's 'moderate'"),
    ("chatgpt", "GZ-14", "lateral_interaction_present"): (True, "gold states this is genuinely ambiguous by design; False is one of the two defensible reads"),
    ("gemini", "GZ-14", "lateral_interaction_present"): (False, "field is null -- no answer given"),
    ("claude", "GZ-14", "lateral_interaction_present"): (True, "gold states this is genuinely ambiguous by design; False is one of the two defensible reads"),
    ("chatgpt", "GZ-25", "interaction_character"): (False, "literal unfilled '<string>' placeholder, no real content"),
    ("gemini", "GZ-25", "interaction_character"): (False, "forces the blanket 'synergistic' label gold explicitly says is the wrong framing"),
    ("claude", "GZ-25", "interaction_character"): (False, "answers a different (CMC-minimum/Rubingh) framing; never addresses the required W12/W21 sign asymmetry"),
    ("chatgpt", "GZ-27", "direction_of_correction"): (True, "'smaller' matches gold's core direction"),
    ("claude", "GZ-27", "direction_of_correction"): (True, "'smaller' with the same Perrin F>1 reasoning as gold"),
    ("chatgpt", "GZ-32", "tail_description"): (True, "8-carbon branched 2-ethylhexyl description matches gold"),
    ("gemini", "GZ-32", "tail_description"): (False, "field is null -- no answer given"),
    ("claude", "GZ-32", "tail_description"): (True, "near-verbatim match to gold's 2-ethylhexyl/8-carbon/branch-at-2 description"),
    ("chatgpt", "GZ-35", "valid_hlb_method"): (False, "names Davies; gold requires Griffin (Davies is confirmed structurally invalid for this headgroup)"),
    ("gemini", "GZ-35", "valid_hlb_method"): (False, "names Davies; gold requires Griffin"),
    ("claude", "GZ-35", "valid_hlb_method"): (False, "explicitly states Griffin is NOT valid and recommends Davies -- backwards from gold"),
    ("chatgpt", "GZ-35", "charge_confidence"): (False, "'high' does not reflect gold's required pH-sensitivity hedge ('moderate')"),
    ("gemini", "GZ-35", "charge_confidence"): (False, "'high' does not reflect gold's required pH-sensitivity hedge ('moderate')"),
    ("claude", "GZ-35", "charge_confidence"): (True, "'medium' is a synonym of gold's 'moderate'"),
    ("chatgpt", "GZ-36", "slope_sign"): (False, "literal unfilled '<string>' placeholder, no real content"),
    ("gemini", "GZ-36", "slope_sign"): (False, "'positive' is the opposite of gold's 'negative'"),
    ("claude", "GZ-36", "slope_sign"): (True, "exact match ('negative'); only flagged needs-judgment because gold string is long"),
    ("chatgpt", "GZ-39", "deviation_reason"): (True, "matches gold's pre-CMC-isotherm-doesn't-capture-micellization reasoning"),
    ("gemini", "GZ-39", "deviation_reason"): (True, "matches gold's pre-CMC-isotherm-doesn't-capture-micellization reasoning"),
    ("claude", "GZ-39", "deviation_reason"): (True, "matches gold's pre-CMC-isotherm-doesn't-capture-micellization reasoning"),
    ("chatgpt", "GZ-43", "agreement_assessment"): (True, "matches gold's ~3.5%/strong-agreement core claim"),
    ("gemini", "GZ-43", "agreement_assessment"): (False, "field is null -- no answer given"),
    ("claude", "GZ-43", "agreement_assessment"): (True, "'excellent agreement' matches gold's strong-agreement core claim"),
    ("chatgpt", "GZ-50", "trend"): (True, "matches gold's monotonic-despite-non-monotonic-CMC claim"),
    ("claude", "GZ-50", "trend"): (True, "matches gold's monotonic-despite-non-monotonic-CMC claim essentially verbatim"),
}


def main():
    with open(GOLD_PATH, "r", encoding="utf-8") as f:
        gold_all = json.load(f)
    templates = extract_templates()

    results = {}
    needs_judgment_log = []
    for model, path in SOURCES.items():
        parser = parse_writer_format if model != "claude" else parse_compiled_format
        if not os.path.exists(path):  # source transcript not committed to repo (e.g. claude Sept dir)
            print(f"  SKIP {model}: source file absent from repo: {path}")
            continue
        answers = parser(path)
        q_results = {}
        for n in range(1, 51):
            qid = f"GZ-{n:02d}"
            gold_entry = gold_all.get(qid)
            if not gold_entry:
                continue
            gold = gold_entry["gold"]
            template_keys = templates.get(qid, set())
            graded_keys = (template_keys & set(gold.keys())) - META_KEYS
            model_answer = answers.get(qid)
            if model_answer is None:
                q_results[qid] = {"verdict": "NO_JSON", "fields": {}}
                continue
            field_results = {}
            all_pass = True
            any_judgment = False
            for k in graded_keys:
                gv = gold[k]
                mv = model_answer.get(k) if isinstance(model_answer, dict) else None
                ok, kind = grade_field(mv, gv)
                if kind == "needs-judgment":
                    resolution = JUDGMENT_RESOLUTIONS.get((model, qid, k))
                    if resolution is None:
                        any_judgment = True
                        needs_judgment_log.append((model, qid, k, gv, mv))
                        field_results[k] = {"ok": None, "kind": kind, "gold": gv, "model": mv}
                        continue
                    ok, rationale = resolution
                    field_results[k] = {"ok": ok, "kind": "judgment", "gold": gv, "model": mv, "rationale": rationale}
                    if not ok:
                        all_pass = False
                    continue
                field_results[k] = {"ok": ok, "kind": kind, "gold": gv, "model": mv}
                if not ok:
                    all_pass = False
            verdict = "NEEDS_JUDGMENT" if any_judgment else ("PASS" if all_pass else "FAIL")
            q_results[qid] = {"verdict": verdict, "fields": field_results}
        results[model] = q_results

    print(f"{'Model':<10} {'PASS':>6} {'FAIL':>6} {'NEEDS_JDG':>10} {'NO_JSON':>8} {'Total':>6} {'Pass%(of gradable)':>20}")
    for model, qr in results.items():
        p = sum(1 for v in qr.values() if v["verdict"] == "PASS")
        f = sum(1 for v in qr.values() if v["verdict"] == "FAIL")
        nj = sum(1 for v in qr.values() if v["verdict"] == "NEEDS_JUDGMENT")
        nojson = sum(1 for v in qr.values() if v["verdict"] == "NO_JSON")
        tot = len(qr)
        gradable = p + f
        pct = f"{100*p/gradable:.1f}%" if gradable else "n/a"
        print(f"{model:<10} {p:>6} {f:>6} {nj:>10} {nojson:>8} {tot:>6} {pct:>20}")

    print(f"\n{len(needs_judgment_log)} fields need human/LLM judgment (long free-text gold, no clean string to auto-match):")
    for model, qid, k, gv, mv in needs_judgment_log:
        print(f"\n--- {model} {qid} field '{k}' ---")
        print(f"  grading criteria (gold): {gv[:300]}")
        print(f"  model said: {mv}")

    with open(os.path.join(BENCH, "unaugmented_2026-09-21", "grade_results.json"), "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, default=str)


if __name__ == "__main__":
    main()
