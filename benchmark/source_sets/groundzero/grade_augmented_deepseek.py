"""Grade the DeepSeek (Better DeepSeek extension + MCP) augmented GroundZero50
(GZ-01..50) run. Reuses grade_unaugmented.py's template-key extraction,
JSON extraction, and field grading so the comparison logic matches what was
already used/validated for this question set.

Usage: python grade_augmented_deepseek.py
"""
import json
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

from grade_unaugmented import extract_last_json_block, extract_templates, grade_field, META_KEYS

HERE = Path(__file__).resolve().parent
GOLD_PATH = HERE / "gold_answers.json"
TEMPLATE_PATH = HERE / "question_bank_augmented.txt"
ANSWERS_DIR = HERE / "augmented_runs_deepseek_2026-10-15" / "groundzero50"

# Manual resolution of the "needs-judgment" free-text fields (long prose gold
# with no clean string to auto-match). Decided by reading each model answer
# against its gold's stated reasoning, not by keyword matching. Where a
# question also has an independent hard-failed field, the judgment call is
# moot -- the question fails regardless -- and is still recorded here for
# transparency/audit rather than left silently unresolved.
JUDGMENT_RESOLUTIONS = {
    ("GZ-04", "charge_confidence"): (True, "model said 'moderate', exact match to gold"),
    ("GZ-14", "lateral_interaction_present"): (True, "gold notes this is genuinely ambiguous by design (a nonzero but Langmuir still wins BIC); model's True reflects the real nonzero magnitude, a defensible read"),
    ("GZ-25", "interaction_character"): (True, "model explicitly flags the asymmetry (opposite-sign W12/W21) gold asks for, not a forced blanket label"),
    ("GZ-27", "direction_of_correction"): (True, "model says 'smaller' with the same Perrin-factor>1 reasoning as gold"),
    ("GZ-32", "tail_description"): (True, "near-verbatim match to gold's 2-ethylhexyl/8-carbon/C2-branch description"),
    ("GZ-35", "valid_hlb_method"): (False, "2026-10-15 retry answer now PREFERS Davies ('Davies group-contribution method (preferred)...Griffin...also acceptable'), contradicting gold's conclusion that Davies is structurally invalid for this class and Griffin is the one valid method -- superseded from the original run's True resolution, which was based on different (now-overwritten) answer text"),
    ("GZ-35", "charge_confidence"): (True, "2026-10-15 retry answer is now cleanly 'moderate', an exact match to gold -- superseded from the original run's False resolution, which was based on a hedged 'high...moderate' answer that this retry replaced"),
    ("GZ-36", "slope_sign"): (True, "exact string match to gold ('negative'); only flagged needs-judgment because gold string is long"),
    ("GZ-39", "deviation_reason"): (True, "model's monomer-saturation/plateau explanation matches gold's physical reasoning"),
    ("GZ-43", "agreement_assessment"): (True, "model's ~3.4% figure and 'consensus value' conclusion matches gold's 'strong agreement, ~3.5%, within normal variation'"),
    ("GZ-50", "trend"): (True, "2026-10-15 retry answer now correctly states 'Monotonically more negative with increasing temperature...no minimum is reproduced', matching gold's core claim -- superseded from the original run's False resolution, which was based on a different (now-overwritten) answer that incorrectly reported a reversal"),
}

KNOWN_LIMITATIONS = {
    "GZ-17": "Better DeepSeek extension markdown-corrupts this DTAB SMILES "
             "([N+](C)(C)C -> [N+](https://c/)(C)C) before the tool call; "
             "model detected the corruption but could not escape it (same "
             "root cause as HARD50's H-30), and degenerated into a repeated-"
             "token loop on retry. Not a toolkit or grading issue.",
    "GZ-49": "No MCP tool exposes the asymmetric-Margules (EOMMM) two-"
             "parameter regression this question needs; the model correctly "
             "recognized the gap and tried to fall back to its sandboxed "
             "code runner, which requires a manual 'Run Code' click the "
             "automated pipeline never provides -- a genuine toolkit-"
             "capability gap, not a DeepSeek/extension bug.",
    "GZ-04": "Better DeepSeek extension truncates this long chiral-ring SMILES "
             "mid-string before the tool call (CCCCCCCC(CC(=O)O)OC(=O)CC(CCCCCCC)"
             "O[C@H]1[C@@H](...)... gets cut off right at '[C@H]', dropping the "
             "entire sugar ring) -- a second, distinct manifestation of the same "
             "extension-level text-scraping unreliability as H-30/GZ-17's "
             "markdown corruption. Retried 3 times total (incl. a fresh page "
             "reload and a dedicated clean-environment retry after fixing an "
             "unrelated MCP SuperAssistant conflict) -- same truncation "
             "reproduced every time; confirmed independent of that conflict, "
             "not fixable via prompting.",
    "GZ-36": "No response was generated at all across 3 attempts (incl. a fresh "
             "page reload and a dedicated clean-environment retry after fixing "
             "an unrelated MCP SuperAssistant conflict) -- same blank-response "
             "symptom every time, so ruled out as the SuperAssistant "
             "interference that explained similar failures elsewhere (e.g. "
             "GZ-26, which DID resolve after that fix). Root cause for GZ-36 "
             "specifically remains unidentified.",
}

# Manual overrides for fields that went through the automatic (non-judgment)
# string_matches path but got a false-negative purely from gold's own
# parenthetical context pushing it over string_matches' 80-char needs-
# judgment threshold or under its 34% token-overlap cutoff, even though the
# model's answer exactly matches gold's actual required value.
MANUAL_OVERRIDES = {
    ("GZ-38", "units_used"): (True, "model said exactly 'µM', an exact match to gold's core "
                              "value -- gold's '(as labeled on the chart)' is descriptive "
                              "context, not a second required token, but it pushed the lenient "
                              "string-match below its 34% overlap threshold (1/3 tokens)."),
    ("GZ-35", "structural_family"): (True, "model's 'alkyl diglycoside (glycosidic surfactant "
                                     "with acetylated sugar headgroup)' is a genuine semantic "
                                     "match to gold's 'glycolipid biosurfactant (disaccharide "
                                     "headgroup)' -- 'di-glycoside' literally means a two-sugar "
                                     "(disaccharide) headgroup, just different vocabulary than "
                                     "gold's phrasing; the lenient token matcher rejected it "
                                     "purely on word choice, not substance."),
}

with open(GOLD_PATH, encoding="utf-8") as f:
    gold_all = json.load(f)
# the file also carries non-question meta keys (_grading_note, _temperature_note)
# and GZ-51..68, which belong to a separate reproducibility_x5 set out of
# scope for this run -- only GZ-01..50 were ever sent to DeepSeek
gold_all = {
    k: v for k, v in gold_all.items()
    if k.startswith("GZ-") and 1 <= int(k.split("-")[1]) <= 50
}

templates = extract_templates()

results = {}
needs_judgment_log = []
for qid, gold_entry in sorted(gold_all.items()):
    num = int(qid.split("-")[1])
    fpath = ANSWERS_DIR / f"GZ{num:02d}.txt"
    gold = gold_entry["gold"]
    template_keys = templates.get(qid, set())
    graded_keys = (template_keys & set(gold.keys())) - META_KEYS

    if not fpath.exists():
        reason = KNOWN_LIMITATIONS.get(qid, "no answer file")
        results[qid] = {"verdict": "NO_JSON", "fields": {}, "reason": reason}
        continue

    text = fpath.read_text(encoding="utf-8", errors="replace")
    model_answer = extract_last_json_block(text)
    if model_answer is None:
        results[qid] = {"verdict": "NO_JSON", "fields": {}, "reason": "unparseable"}
        continue

    field_results = {}
    all_pass = True
    unresolved_judgment = False
    for k in graded_keys:
        gv = gold[k]
        mv = model_answer.get(k) if isinstance(model_answer, dict) else None
        override = MANUAL_OVERRIDES.get((qid, k))
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
            else:
                ok, rationale = resolution
                field_results[k] = {"ok": ok, "kind": "judgment", "gold": gv, "model": mv, "rationale": rationale}
                if not ok:
                    all_pass = False
                continue
        field_results[k] = {"ok": ok, "kind": kind, "gold": gv, "model": mv}
        if kind != "needs-judgment" and not ok:
            all_pass = False
    verdict = "NEEDS_JUDGMENT" if unresolved_judgment else ("PASS" if all_pass else "FAIL")
    results[qid] = {"verdict": verdict, "fields": field_results}

p = sum(1 for v in results.values() if v["verdict"] == "PASS")
f = sum(1 for v in results.values() if v["verdict"] == "FAIL")
nj = sum(1 for v in results.values() if v["verdict"] == "NEEDS_JUDGMENT")
nojson = sum(1 for v in results.values() if v["verdict"] == "NO_JSON")
tot = len(results)
gradable = p + f
pct = f"{100*p/gradable:.1f}%" if gradable else "n/a"

print(f"DeepSeek augmented GroundZero50: PASS={p} FAIL={f} NEEDS_JUDGMENT={nj} NO_JSON={nojson} / {tot}")
print(f"Pass rate (of gradable, excluding needs-judgment/no-json): {pct}")
print()
for qid, r in sorted(results.items()):
    line = f"  {qid}: {r['verdict']}"
    if r["verdict"] == "FAIL":
        bad = [k for k, fv in r["fields"].items() if fv["ok"] is False]
        line += f"  bad fields: {bad}"
    if r["verdict"] == "NO_JSON" and "reason" in r:
        line += f"  ({r['reason']})"
    print(line)

if needs_judgment_log:
    print(f"\n{len(needs_judgment_log)} fields need human/LLM judgment (long free-text gold):")
    for qid, k, gv, mv in needs_judgment_log:
        print(f"\n--- {qid} field '{k}' ---")
        print(f"  grading criteria (gold): {str(gv)[:300]}")
        print(f"  model said: {mv}")

with open(HERE / "deepseek_augmented_grade_results.json", "w", encoding="utf-8") as f:
    json.dump(results, f, indent=1, default=str)
