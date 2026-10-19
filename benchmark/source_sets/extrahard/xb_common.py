"""Shared helpers for EXTRA50 (X-01..X-50) gold construction."""
from __future__ import annotations
import json
import math
import re

R_GAS = 8.314462618        # J/(mol.K) — thermodynamics.R_GAS
WATER_M = 55.5             # mol/L  — thermodynamics.WATER_MOLARITY_M
N_A = 6.02214076e23        # avogadro (exact, CODATA 2019)


# ---------------------------------------------------------------- fields
def num(value: float, tol: float, unit: str = "") -> dict:
    assert tol <= 0.02, f"tol {tol} > 0.02"
    assert math.isfinite(float(value)) and float(value) != 0.0, f"bad gold {value}"
    d = {"type": "num", "tol": float(tol), "value": round(float(value), 10)}
    if unit:
        d["unit"] = unit
    return d


def integer(value: int) -> dict:
    return {"type": "int", "value": int(value)}


def exact(value: str) -> dict:
    return {"type": "str", "value": str(value)}


def boolean(value: bool) -> dict:
    return {"type": "bool", "value": bool(value)}


# ------------------------------------------------------- body assembly
def outputs_block(items: list[tuple[str, str]]) -> str:
    """items: (json_key, human description). Emits the Required outputs
    list + the FINAL_JSON contract + the END OF QUESTION marker, built
    from the SAME key list as the fields dict (key-parity by design)."""
    lines = ["", "Required outputs (report every item, in this order):"]
    for i, (_k, desc) in enumerate(items, 1):
        lines.append(f"{i}. {desc}")
    fj = ", ".join(f'"{k}": <value>' for k, _ in items)
    lines += ["",
              "End with exactly one line of this form (JSON, keys exactly as written):",
              f"FINAL_JSON: {{{fj}}}",
              "END OF QUESTION"]
    return "\n".join(lines)


def qhead(qid: str, title: str) -> str:
    return f"Question {qid} — {title}\n\n"


# ------------------------------------------------------- validations
_FJ = re.compile(r"^FINAL_JSON: \{(.*)\}\s*$", re.M)


def _fj_keys(body: str) -> list[str]:
    m = None
    for m in _FJ.finditer(body):
        pass
    if not m:
        return []
    inner = "{" + m.group(1) + "}"
    inner = inner.replace("<value>", "0").replace("true", "True").replace("false", "False")
    try:
        return list(json.loads(m.group(1) and _pyjson(inner)).keys())
    except Exception:
        # fall back to regex keys
        return re.findall(r'"([A-Za-z0-9_]+)"\s*:', m.group(1))


def _pyjson(s: str) -> str:
    return s


def validate(all_q: dict) -> list[str]:
    errs = []
    banned = [r"\bas above\b", r"\bprevious question", r"\bH-\d\d\b", r"\bGZ-\d\d\b",
              r"\bscale50\b", r"\bFigure\s*\d", r"\bimage\b", r"\bchart\b", r"\bplot below\b", r"attached (file|figure|table|image|plot)",
              r"as shown in"]
    for qid, q in sorted(all_q.items()):
        body = q["body"]
        # 1. key parity
        fk = list(q["fields"].keys())
        bk = _fj_keys(body)
        if fk != bk:
            errs.append(f"{qid}: field keys {fk} != FINAL_JSON keys {bk}")
        # 2. tolerances
        for k, f in q["fields"].items():
            if f["type"] == "num" and f["tol"] > 0.02:
                errs.append(f"{qid}.{k}: tol {f['tol']} > 0.02")
            if f["type"] == "num" and (not math.isfinite(f["value"]) or f["value"] == 0):
                errs.append(f"{qid}.{k}: non-finite/zero gold")
        # 3. self-containment / bans
        for pat in banned:
            if re.search(pat, body, re.I):
                errs.append(f"{qid}: banned pattern {pat!r}")
        # 4. own id only reference
        ids = set(re.findall(r"\bX-\d\d\b", body))
        if ids - {qid}:
            errs.append(f"{qid}: references other question ids {ids - {qid}}")
        # 5. marker
        if not body.rstrip().endswith("END OF QUESTION"):
            errs.append(f"{qid}: missing END OF QUESTION")
        # 6. axes vs fields
        n_arith = sum(1 for f in q["fields"].values() if f["type"] in ("num", "int"))
        n_reason = sum(1 for f in q["fields"].values() if f["type"] in ("str", "bool"))
        if q.get("axes") != {"arith": n_arith, "reason": n_reason}:
            errs.append(f"{qid}: axes {q.get('axes')} != counted {n_arith}/{n_reason}")
        # 7. levers
        if len(q.get("levers", [])) < 2:
            errs.append(f"{qid}: fewer than 2 levers")
        if not (4 <= len(q["fields"]) <= 8):
            errs.append(f"{qid}: {len(q['fields'])} fields outside 4..8")
    # global: no duplicated body text blocks (>60 chars) across questions
    seen = {}
    for qid, q in all_q.items():
        for line in q["body"].split("\n"):
            s = line.strip()
            if s in _BOILERPLATE:
                continue
            if len(s) > 60 and not s.startswith(("Question", "FINAL_JSON", "1.", "2.", "3.", "4.", "5.", "6.", "7.")):
                if s in seen and seen[s] != qid:
                    errs.append(f"{qid}/{seen[s]}: shared body line {s[:70]!r}")
                seen[s] = qid
    return errs


_BOILERPLATE = {
    "Required outputs (report every item, in this order):",
    "End with exactly one line of this form (JSON, keys exactly as written):",
}


def rel_close(a: float, b: float, tol: float) -> bool:
    return abs(a - b) / abs(b) <= tol


def check_near_duplicates(bodies: dict) -> list[str]:
    """Exact normalized-line overlap (>=60 chars) between different prompts."""
    import re as _re
    norm_lines: dict[str, set[str]] = {}
    for k, b in bodies.items():
        lines = set()
        for ln in b.splitlines():
            t = _re.sub(r"\s+", " ", ln).strip().lower()
            if len(t) >= 60 and "keys exactly as written" not in t \
                    and not t.startswith('{"') \
                    and not _re.match(r"^\d+\. ", t):
                lines.add(t)
        norm_lines[k] = lines
    keys = sorted(norm_lines)
    out = []
    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            shared = norm_lines[a] & norm_lines[b]
            if shared:
                ex = sorted(shared)[0]
                out.append(f"{a}/{b}: {len(shared)} shared line(s), e.g. {ex[:70]!r}")
    return out
