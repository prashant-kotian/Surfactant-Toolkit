"""Grade Antigravity (gemini-3.8-flash-high) GroundZero50 runs, both
conditions. Concatenates the 5 batch files per condition, splits on
'=== GZ## ===' markers, extracts each question's FINAL_JSON, and grades
against gold_answers.json. Reuses grade_unaugmented.py's extractor/
comparator/template logic (same as grade_augmented_deepseek.py
did earlier).
"""
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from grade_unaugmented import extract_last_json_block, extract_templates, grade_field, META_KEYS

HERE = Path(__file__).resolve().parent
GOLD_PATH = HERE / "gold_answers.json"

with open(GOLD_PATH, encoding="utf-8") as f:
    gold_all = json.load(f)
gold_all = {k: v for k, v in gold_all.items() if k.startswith("GZ-") and 1 <= int(k.split("-")[1]) <= 50}
templates = extract_templates()

# Judgment-field resolutions specific to THIS Antigravity run (read fresh,
# not reused from DeepSeek's run -- different model, different answers).
JUDGMENT_RESOLUTIONS = {
    ("GZ-04", "charge_confidence"): (True, "exact match to gold"),
    ("GZ-14", "lateral_interaction_present"): (True, "gold states this is genuinely ambiguous by design; both True (augmented run) and False (unaugmented run) are defensible reads depending on which signal is weighted (nonzero a vs. Langmuir winning BIC)"),
    ("GZ-25", "interaction_character"): (True, "model explicitly states the required asymmetry (opposite-sign W12/W21), matching gold's core requirement"),
    ("GZ-27", "direction_of_correction"): (True, "model says 'smaller' with the same Perrin-factor>1 reasoning as gold"),
    ("GZ-32", "tail_description"): (True, "near-verbatim match to gold's 2-ethylhexyl/8-carbon/C2-branch description"),
    ("GZ-35", "charge_confidence"): (True, "exact match to gold"),
    ("GZ-35", "valid_hlb_method"): (True, "correctly names Griffin and correctly explains Davies' structural invalidity for this headgroup"),
    ("GZ-36", "slope_sign"): (True, "exact match ('negative'); only flagged needs-judgment because gold string is long"),
    ("GZ-39", "deviation_reason"): (True, "matches gold's monomer-saturation/plateau physical reasoning"),
    ("GZ-43", "agreement_assessment"): (True, "matches gold's 'strong agreement, ~3.5%' core claim"),
    ("GZ-50", "trend"): (True, "matches gold's monotonic-despite-non-monotonic-CMC claim essentially verbatim"),
}

# Manual overrides for fields the grader got structurally wrong -- either a
# genuine string_matches false-negative from gold's own parenthetical context
# (same class of issue found and fixed for DeepSeek's grading this session),
# or a gold field type (list) the shared grade_field() doesn't handle at all
# (falls through to its "unknown" branch, which always returns False
# regardless of correctness).
MANUAL_OVERRIDES = {
    ("GZ-01", "isotherm_model", "augmented"): (True, "'langmuir' is an exact match to gold's core model name; '(BIC-selected over Frumkin)' is justification text, not a second required token"),
    ("GZ-03", "structural_family", "augmented"): (True, "'dimeric_gemini_type' is the same classification as gold's 'dimeric/gemini' in this project's own established label convention (see e.g. X-32/GZ-23), just the underscored-identifier form instead of the slash form"),
    ("GZ-23", "valid_hlb_method", "augmented"): (True, "'griffin' exact match to gold's core value"),
    ("GZ-23", "structural_family", "unaugmented"): (True, "'dimeric/gemini' exact match to gold's core value"),
    ("GZ-23", "valid_hlb_method", "unaugmented"): (True, "'griffin' exact match to gold's core value"),
    ("GZ-29", "selected_model", "unaugmented"): (True, "'frumkin' exact match to gold's core value; '(wins by BIC, delta_bic=+4.40)' is justification text"),
    ("GZ-38", "units_used", "unaugmented"): (True, "model said exactly the chart's unit, matching gold's core value; '(as labeled on the chart)' is descriptive context"),
    ("GZ-42", "more_likely_correct", "augmented"): (True, "model correctly identifies Source 2 AND gold's exact numeric range (2.4 to 3.5) -- a genuinely correct answer; the grader's string match failed only because of wrapper phrasing ('Source 2 (Cc in the range ...)' vs gold's number-first phrasing)"),
    ("GZ-42", "more_likely_correct", "unaugmented"): (True, "same genuine correct answer as the augmented run -- Source 2, exact gold numeric range"),
    ("GZ-05", "sub_intervals", "augmented"): (True, "gold is a list, a type grade_field() does not handle at all (falls through to its catch-all False) -- actual values match to 4+ significant figures (-12.5162 vs gold -12.51622692462097, -14.8101 vs -14.810140401192402)"),
    ("GZ-05", "sub_intervals", "unaugmented"): (True, "same grade_field() list-type gap; values match to 3-4 significant figures (-12.516 vs -12.51622692462097, -14.81 vs -14.810140401192402)"),
}


def load_condition(raw_dir: Path) -> dict:
    full_text = ""
    for i in range(1, 6):
        full_text += (raw_dir / f"batch{i}.txt").read_text(encoding="utf-8", errors="replace") + "\n"
    chunks = re.split(r"=== (GZ\d+) ===", full_text)
    blocks = {}
    for i in range(1, len(chunks), 2):
        qid_raw = chunks[i]
        content = chunks[i + 1]
        num = int(qid_raw[2:])
        qid = f"GZ-{num:02d}"
        blocks[qid] = content
    return blocks


def grade_condition(label: str, raw_dir: Path):
    blocks = load_condition(raw_dir)
    expected_ids = {f"GZ-{n:02d}" for n in range(1, 51)}
    missing = sorted(expected_ids - set(blocks.keys()), key=lambda x: int(x.split("-")[1]))

    results = {}
    needs_judgment_log = []
    for qid in sorted(expected_ids, key=lambda x: int(x.split("-")[1])):
        gold_entry = gold_all[qid]
        gold = gold_entry["gold"]
        template_keys = templates.get(qid, set())
        graded_keys = (template_keys & set(gold.keys())) - META_KEYS

        if qid not in blocks:
            results[qid] = {"verdict": "NO_JSON", "fields": {}, "reason": "question ID marker not found"}
            continue
        model_answer = extract_last_json_block(blocks[qid])
        if model_answer is None:
            results[qid] = {"verdict": "NO_JSON", "fields": {}, "reason": "no FINAL_JSON parsed"}
            continue

        field_results = {}
        all_pass = True
        unresolved_judgment = False
        for k in graded_keys:
            gv = gold[k]
            mv = model_answer.get(k) if isinstance(model_answer, dict) else None
            override = MANUAL_OVERRIDES.get((qid, k, label))
            if override is not None:
                ok, rationale = override
                field_results[k] = {"ok": ok, "kind": "manual-override", "gold": gv, "model": mv, "rationale": rationale}
                if not ok:
                    all_pass = False
                continue
            ok, kind = grade_field(mv, gv)
            if kind == "needs-judgment":
                resolution = JUDGMENT_RESOLUTIONS.get((qid, k))
                if resolution is None:
                    unresolved_judgment = True
                    needs_judgment_log.append((qid, k, gv, mv))
                    field_results[k] = {"ok": None, "kind": "needs-judgment", "gold": gv, "model": mv}
                    continue
                else:
                    ok, rationale = resolution
                    field_results[k] = {"ok": ok, "kind": "judgment", "gold": gv, "model": mv, "rationale": rationale}
                    if not ok:
                        all_pass = False
                    continue
            field_results[k] = {"ok": ok, "kind": kind, "gold": gv, "model": mv}
            if not ok:
                all_pass = False

        verdict = "NEEDS_JUDGMENT" if unresolved_judgment else ("PASS" if all_pass else "FAIL")
        results[qid] = {"verdict": verdict, "fields": field_results}

    n_pass = sum(1 for r in results.values() if r["verdict"] == "PASS")
    n_fail = sum(1 for r in results.values() if r["verdict"] == "FAIL")
    n_judgment = sum(1 for r in results.values() if r["verdict"] == "NEEDS_JUDGMENT")
    n_nojson = sum(1 for r in results.values() if r["verdict"] == "NO_JSON")
    print(f"\n=== GroundZero50 {label} (gemini-3.8-flash-high): PASS={n_pass} FAIL={n_fail} "
          f"NEEDS_JUDGMENT={n_judgment} NO_JSON={n_nojson} / 50 ===")
    if missing:
        print(f"  MISSING question IDs: {missing}")
    for qid, r in results.items():
        if r["verdict"] not in ("PASS",):
            bad = [k for k, fv in r.get("fields", {}).items() if fv.get("ok") is False]
            reason = r.get("reason", "")
            print(f"  {qid}: {r['verdict']}  {bad if bad else reason}")

    out_path = HERE / f"antigravity_{label}_grade_results.json"   # the file build_results.py reads
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=1, default=str)

    if needs_judgment_log:
        print(f"\n  {len(needs_judgment_log)} UNRESOLVED judgment fields for {label}:")
        for qid, k, gv, mv in needs_judgment_log:
            print(f"\n  --- {qid} field '{k}' ---")
            print(f"    gold: {str(gv)[:250]}")
            print(f"    model: {str(mv)[:250]}")

    return results, n_pass, needs_judgment_log


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    grade_condition("augmented", HERE / "antigravity_aug_raw")
    # unaugmented condition kept for reference only (2026-10-21): agy had coding/web-search tools, so it was
    # not a true unaugmented baseline; Gemini web is the unaugmented Gemini arm. Not used in comparisons.
    if "--with-unaug" in sys.argv:
        grade_condition("unaugmented", HERE / "antigravity_unaug_raw")
