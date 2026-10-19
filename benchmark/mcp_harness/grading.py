"""Grade an extracted model answer against a question's gold answer.

Handles the three gold-answer shapes the question bank uses:
  * scalar numeric  (tolerance {"rel": x} or {"abs": x})
  * dict of named quantities, some numeric, some string classifications
    (tolerance {"<key>_rel": x, "<key>_abs": x, ...} or a blanket {"rel": x})
  * a no-solution / category string (grading_method == "category_match")

Grading is intentionally strict on the numbers (that is the whole point of the
benchmark) but tolerant on formatting: model keys are matched to gold keys by a
normalised form, so "beta" vs "Beta" vs "beta_value" still line up.
"""

from __future__ import annotations
import re
from typing import Any

_DEFAULT_REL = 0.02
_NOSOL_PATTERNS = [
    "no solution", "no valid", "undefined", "cannot compute", "cannot be computed",
    "not physically", "not meaningful", "does not exist", "no physically",
    "not consistent", "not defined", "no micelle", "below the cmc", "below cmc",
    "not applicable", "no root",
    # Added 2026-09-04: found via the hardened gibbs_surface_excess trap_audit
    # -- Claude Opus 4.8 correctly refused all 6 questions but phrased it as
    # "cannot be determined" / "not solvable" rather than any pattern already
    # in this list, so a real improvement was being mis-scored as a failure.
    "cannot be determined", "cannot determine", "not solvable",
]

# Added 2026-09-07 for the Tier 3 reliability categories (P/Q/R/S) --
# applicability_domain_violation is a DIFFERENT flag from no_solution: the
# question IS answerable (a number can be produced), but the model should
# recognize the prediction is an extrapolation outside the fitted training
# range and flag it as lower-confidence, not silently report it with the
# same certainty as an in-domain interpolation. A correct answer here can
# still legitimately include a numeric value alongside the flag -- grading
# only checks whether the flag/caveat language is present, not whether a
# number was withheld (that would wrongly penalize a model for answering
# AND correctly caveating, which is the ideal response, not a failure).
_APPLICABILITY_PATTERNS = [
    "out of domain", "outside the domain", "out-of-domain", "outside the training",
    "beyond the training", "beyond the fitted", "outside the fitted",
    "extrapolat",  # covers extrapolate/extrapolation/extrapolating
    "low confidence", "low-confidence", "less confident", "less reliable",
    "not reliable outside", "should not be trusted", "treat with caution",
    "outside its training", "outside the range", "beyond the range",
    # Added 2026-09-07: found via real Tier 3 pilot grading -- ChatGPT's
    # P-002 correctly identified "C20 is also outside the C6-C18 fit range"
    # but the literal chain-length token ("C6-C18") sitting between "outside
    # the" and "range" broke every substring pattern above. A model
    # describing "outside the <specific range/window> range/window" is
    # exactly the behavior this trap wants to reward, not miss on phrasing.
    "training range", "fit range", "fitted range", "calibration range",
    "calibration window", "training window", "outside the stated",
    "beyond the stated",
]


def _norm(k: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(k).lower())


def _as_float(v: Any):
    if isinstance(v, bool):
        return None
    if isinstance(v, (int, float)):
        return float(v)
    if isinstance(v, str):
        m = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", v.replace(",", ""))
        return float(m.group()) if m else None
    return None


def _num_ok(model_v, gold_v, rel, abs_):
    a = _as_float(model_v)
    if a is None:
        return False
    if abs_ is not None and abs(a - gold_v) <= abs_:
        return True
    if rel is not None:
        if gold_v == 0:
            return abs(a) <= (abs_ if abs_ is not None else 1e-9)
        return abs(a - gold_v) <= rel * abs(gold_v)
    return False


def _tol_for(key, tolerance):
    """Resolve (rel, abs) tolerance for a given gold key."""
    nk = _norm(key)
    rel = abs_ = None
    for tk, tv in tolerance.items():
        ntk = _norm(tk)
        if ntk == nk + "rel" or ntk == nk + "relative":
            rel = tv
        elif ntk == nk + "abs" or ntk == nk + "absolute":
            abs_ = tv
    if rel is None and abs_ is None:  # fall back to blanket rel/abs
        rel = tolerance.get("rel")
        abs_ = tolerance.get("abs")
    if rel is None and abs_ is None:
        rel = _DEFAULT_REL
    return rel, abs_


def _find_model_val(gold_key, model: dict):
    ngk = _norm(gold_key)
    # exact normalised match first
    for mk, mv in model.items():
        if _norm(mk) == ngk:
            return mv
    # containment either direction (beta vs betasigma, x1 vs micellarmolefractionx1)
    for mk, mv in model.items():
        nmk = _norm(mk)
        if ngk and (ngk in nmk or nmk in ngk):
            return mv
    return _MISSING


_MISSING = object()


def _string_ok(model_v, gold_v):
    g = _norm(gold_v)
    m = _norm(model_v)
    if not g:
        return True
    if g in m or m in g:
        return True
    # Real bug found via pilot analysis 2026-09-04: a raw substring check is
    # sensitive to word ORDER in multi-word categorical answers -- gold
    # "transitional/viscous-dominated" vs a model's scientifically identical
    # "...viscous/transitional-dominated..." (words swapped) failed this
    # check even though the classification is the same. Fall back to a
    # token-set comparison (split on the original, non-normalized text into
    # word tokens, compare as sets) so word order doesn't matter, while still
    # requiring every gold word to actually appear -- not a free pass, still
    # catches a genuinely wrong category (e.g. "capillary-dominated" would
    # not satisfy gold "transitional/viscous-dominated" since "transitional"
    # and "viscous" are both absent from the model's answer).
    g_tokens = set(re.findall(r"[a-z0-9]+", str(gold_v).lower()))
    m_tokens = set(re.findall(r"[a-z0-9]+", str(model_v).lower()))
    if g_tokens.issubset(m_tokens):
        return True
    # Real bug found 2026-09-07 via Tier 3 pilot grading: ChatGPT's "...
    # (rodlike/wormlike) micelles" was marked wrong against gold "cylindrical/
    # rodlike micelle" purely because "micelle" != "micelles" as exact
    # tokens -- a scientifically identical classification failed on plain
    # pluralization. A word inserted between two gold tokens (here
    # "wormlike" between "rodlike" and "micelle(s)") also defeats the
    # substring-containment check above, so both existing checks missed it.
    # Fix: also accept a match where every gold token's simple singular form
    # (strip a trailing "s") is a substring of some model token's singular
    # form -- deliberately crude (no real stemmer), matches this project's
    # existing "good enough for short domain keywords" tolerance elsewhere.
    def _singular(tok):
        return tok[:-1] if tok.endswith("s") and len(tok) > 3 else tok
    g_sing = {_singular(t) for t in g_tokens}
    m_sing = {_singular(t) for t in m_tokens}
    return g_sing.issubset(m_sing)


def _looks_like_no_solution(text: str) -> bool:
    t = (text or "").lower()
    return any(p in t for p in _NOSOL_PATTERNS)


def _flags_applicability_domain(text: str) -> bool:
    t = (text or "").lower()
    return any(p in t for p in _APPLICABILITY_PATTERNS)


def grade(gold_answer, tolerance: dict, grading_method: str,
          model_answer: dict | None, raw_text: str = "", trap_type: str = "none") -> dict:
    """Return {'correct': bool, 'reason': str, 'per_key': {...}}.

    `category_match` is used for two genuinely different things depending on
    trap_type, and conflating them was a real bug caught by offline testing:
      - trap_type == "no_solution" (Tier 1): the gold "answer" IS the concept
        "no valid solution exists" -- check whether the model's text flags
        that, regardless of what specific words it uses.
      - anything else (Tier 2 classification/keyword questions): gold_answer
        is an actual short keyword/phrase ("anionic", "gemini", "synergistic")
        that the model's answer needs to contain -- a fuzzy keyword match
        against that specific phrase, NOT a no-solution-language check.
    """
    tolerance = tolerance or {}

    # --- true no-solution / unsolvable-trap questions ---
    if trap_type == "no_solution":
        hay = ""
        if model_answer:
            hay += " ".join(str(v) for v in model_answer.values())
        hay += " " + (raw_text or "")
        ok = _looks_like_no_solution(hay)
        return {"correct": ok,
                "reason": "recognised no-solution/undefined" if ok
                          else "did not flag as no-solution/undefined",
                "per_key": {}}

    # --- applicability-domain violation (Tier 3, category P) ---
    if trap_type == "applicability_domain_violation":
        hay = ""
        if model_answer:
            hay += " ".join(str(v) for v in model_answer.values())
        hay += " " + (raw_text or "")
        ok = _flags_applicability_domain(hay)
        return {"correct": ok,
                "reason": "flagged out-of-domain/extrapolation" if ok
                          else "did not flag out-of-domain/extrapolation risk",
                "per_key": {}}

    # --- keyword/classification questions (gold is a short string, not a trap) ---
    if grading_method == "category_match" and isinstance(gold_answer, str):
        hay = ""
        if model_answer:
            hay += " ".join(str(v) for v in model_answer.values()) + " "
        hay += (raw_text or "")
        ok = _string_ok(hay, gold_answer)
        return {"correct": ok,
                "reason": f"model text {'contains' if ok else 'does NOT contain'} "
                          f"expected keyword/phrase '{gold_answer}'",
                "per_key": {}}

    if model_answer is None:
        return {"correct": False, "reason": "no parseable FINAL_JSON answer", "per_key": {}}

    # --- scalar gold ---
    if isinstance(gold_answer, (int, float)) and not isinstance(gold_answer, bool):
        rel, abs_ = _tol_for("value", tolerance)
        mv = model_answer.get("value", _MISSING)
        if mv is _MISSING:
            nums = [v for v in model_answer.values() if _as_float(v) is not None]
            mv = nums[0] if len(nums) >= 1 else _MISSING
        if isinstance(mv, dict):
            # Real bug found via pilot analysis 2026-09-04: for a multi-step
            # question with a scalar gold (the question asks the model to
            # compute an intermediate quantity first, then a final one, but
            # only the final value is graded), a model that reports BOTH
            # values nested under "value" (e.g. {"Gamma_max": ..., "A_min":
            # ...}) used to fail here outright -- mv was the whole dict, not
            # a number, so _num_ok() below would always be False even when
            # the model's actual final answer was correct. Fixed by
            # recursively searching every numeric leaf (not just top-level
            # values) for one that satisfies tolerance against gold, so a
            # model that shows its intermediate work isn't penalized for it.
            def _leaves(d):
                for v in (d.values() if isinstance(d, dict) else d):
                    if isinstance(v, dict):
                        yield from _leaves(v)
                    else:
                        yield v
            candidates = [v for v in _leaves(mv) if _as_float(v) is not None]
            match = next((v for v in candidates if _num_ok(v, float(gold_answer), rel, abs_)), _MISSING)
            mv = match if match is not _MISSING else (candidates[-1] if candidates else _MISSING)
        ok = mv is not _MISSING and _num_ok(mv, float(gold_answer), rel, abs_)
        return {"correct": ok,
                "reason": f"model={mv} gold={gold_answer} (rel={rel}, abs={abs_})",
                "per_key": {"value": ok}}

    # --- dict gold ---
    if isinstance(gold_answer, dict):
        per = {}
        for gk, gv in gold_answer.items():
            mv = _find_model_val(gk, model_answer)
            if mv is _MISSING:
                per[gk] = False
                continue
            if isinstance(gv, (int, float)) and not isinstance(gv, bool):
                rel, abs_ = _tol_for(gk, tolerance)
                per[gk] = _num_ok(mv, float(gv), rel, abs_)
            else:  # string classification (synergy, morphology, regime, ...)
                per[gk] = _string_ok(mv, gv)
        ok = all(per.values()) and len(per) > 0
        return {"correct": ok,
                "reason": "; ".join(f"{k}={'ok' if v else 'FAIL'}" for k, v in per.items()),
                "per_key": per}

    return {"correct": False, "reason": f"unhandled gold type {type(gold_answer)}", "per_key": {}}
