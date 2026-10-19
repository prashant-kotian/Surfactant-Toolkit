"""Extract a model's final answer from its free-text response.

Every model (in every condition) is instructed to end with a line:
    FINAL_JSON: {...}
mapping each requested quantity to a numeric value (and any classification to a
string). This standardised, condition-neutral instruction is what makes grading
uniform across models and conditions. Extraction is deliberately forgiving --
a model that gets the physics right but formats the trailing line slightly off
should not be scored wrong for that -- so several fallbacks are attempted.
"""

from __future__ import annotations
import json
import re


def extract_final(text: str) -> dict | None:
    """Return the parsed FINAL_JSON dict, or None if nothing parseable is found."""
    if not text:
        return None

    # 1. explicit FINAL_JSON: {...} marker, last occurrence wins
    for m in reversed(list(re.finditer(r"FINAL_JSON\s*:\s*", text, re.IGNORECASE))):
        obj = _json_after(text, m.end())
        if obj is not None:
            return obj

    # 2. a fenced ```json ... ``` block, last one
    blocks = re.findall(r"```(?:json)?\s*(\{.*?\})\s*```", text, re.DOTALL)
    for blk in reversed(blocks):
        try:
            return json.loads(blk)
        except Exception:
            pass

    # 3. the last balanced {...} object anywhere in the text
    obj = _last_json_object(text)
    if obj is not None:
        return obj

    return None


def _json_after(text: str, start: int) -> dict | None:
    """Parse a JSON object beginning at/after `start` (skips leading whitespace)."""
    i = start
    while i < len(text) and text[i] in " \t\r\n":
        i += 1
    if i >= len(text) or text[i] != "{":
        # maybe a bare value like `FINAL_JSON: 3.14` -> wrap it
        tail = text[start:].strip().splitlines()[0].strip() if text[start:].strip() else ""
        num = _first_number(tail)
        return {"value": num} if num is not None else None
    return _balanced(text, i)


def _last_json_object(text: str) -> dict | None:
    last = None
    for m in re.finditer(r"\{", text):
        obj = _balanced(text, m.start())
        if obj is not None:
            last = obj
    return last


def _balanced(text: str, i: int) -> dict | None:
    depth = 0
    for j in range(i, len(text)):
        if text[j] == "{":
            depth += 1
        elif text[j] == "}":
            depth -= 1
            if depth == 0:
                try:
                    val = json.loads(text[i:j + 1])
                    return val if isinstance(val, dict) else None
                except Exception:
                    return None
    return None


def _first_number(s: str):
    m = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", s)
    return float(m.group()) if m else None
