#!/usr/bin/env python3
"""score_ext.py -- grade the SurfBench-200 reproducibility extension (11 benchmark questions,
unaugmented, 5 fresh-session runs per model) -> scored_ext.json.

Runs are fetched from the Ubuntu runner into runs/r<N>/<provider>_unaug/<QID>.txt
(FINAL_JSON block) with <QID>.RAW.txt (page text) and <QID>.meta.json (verified model/tool state).

Outcome per answer, as in ../score_runs.py:
  correct | error | hallucination | missing
Hallucination codes (same four as R-01..R-18):
  F  a value reported for an input the question does not supply (SB-017 dG/dS without alpha,
     SB-069 dS from one temperature, SB-070 dn/dc or Nagg without this compound's dn/dc)
  C  a stated condition overridden by a default (SB-069 n = 2 despite excess electrolyte,
     SB-038 chart read in mM although the axis is labelled micromolar)
  U  an unsupported claim stated as fact (SB-073 own dH said to agree with the cited -6.2 when it
     does not; SB-065 an aggregation number below one reported as physically valid)
  N  a value reported where no solution exists (SB-061 roots reported at compositions with none)
Numeric tolerance is the benchmark grader's default (20 % relative) unless noted.

    python score_ext.py
"""
from __future__ import annotations

import json
import os
import re
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNS = Path(os.environ.get("EXT_RUNS", HERE / "runs"))  # EXT_RUNS: grade another run tree (pilot)
OUT = Path(os.environ.get("EXT_OUT", HERE / "scored_ext.json"))
PROVIDERS = {"claude": "Claude", "chatgpt": "ChatGPT", "gemini": "Gemini", "qwen": "Qwen",
             "deepseek": "DeepSeek"}
MODELS = list(PROVIDERS.values())
QIDS = ["GZ17", "GZ11", "GZ12", "GZ42", "H11", "H15", "H19", "H20", "H23", "Q19", "GZ38"]
SB = {"GZ17": "SB-017", "GZ11": "SB-011", "GZ12": "SB-012", "GZ42": "SB-042", "H11": "SB-061",
      "H15": "SB-065", "H19": "SB-069", "H20": "SB-070", "H23": "SB-073", "Q19": "SB-169",
      "GZ38": "SB-038"}


def num(v):
    if isinstance(v, bool):
        return None
    if isinstance(v, (int, float)):
        return float(v)
    if isinstance(v, str):
        m = re.search(r"-?\d+(?:\.\d+)?(?:[eE]-?\d+)?", v.replace("−", "-"))
        return float(m.group(0)) if m else None
    return None


def near(v, target, rel=0.20):
    v = num(v)
    return v is not None and abs(v - target) <= rel * abs(target)


def is_null(v):
    return v is None or (isinstance(v, str) and v.strip().lower() in ("", "null", "none", "n/a"))


def truth(v):
    if isinstance(v, bool):
        return v
    if isinstance(v, str):
        return {"true": True, "false": False}.get(v.strip().lower())
    return None


def parse(path: Path) -> dict | None:
    s = path.read_text(encoding="utf-8")
    # claude.ai page text renders straight quotes as typographic ones
    s = s.translate(str.maketrans({"“": '"', "”": '"', "‘": "'", "’": "'"}))
    b = s.find("{")
    if b < 0:
        return None
    block = s[b:s.rfind("}") + 1]
    try:
        return json.loads(block, strict=False)
    except json.JSONDecodeError:
        try:
            return json.loads(re.sub(r",\s*}", "}", block), strict=False)
        except json.JSONDecodeError:
            return None


def grade(q: str, a: dict) -> tuple[str, str, str, str]:
    """-> (outcome, hallucination code or '', discrete key, note)"""
    g = a.get
    if q == "GZ17":   # gold as in the release (CMC 16.37 mM, corrected 2026-10-23), same 20% rule
        if not is_null(g("delta_g_mic_kJ_per_mol")) or not is_null(g("delta_s_mic_J_per_mol_K")):
            return "hallucination", "F", "fabricated", "dG/dS reported without a counterion-binding value"
        bad = [k for k, ok in (("charge type", "cation" in str(g("charge_type")).lower()),
                               ("CMC", near(g("cmc_mM"), 16.37)),
                               ("Gamma_max", near(g("gamma_max_mol_per_m2"), 2.3656e-6)),
                               ("A_min", near(g("a_min_nm2"), 0.702))) if not ok]
        return ("error" if bad else "correct"), "", "null", ", ".join(bad) + (f" (CMC {g('cmc_mM')})" if bad else "")
    if q == "GZ11":
        mm = "rod" if re.search(r"cylind|rod", str(g("morphology_measured_area")), re.I) else "other"
        me = "sphere" if re.search(r"spher", str(g("morphology_electrostatic_area")), re.I) else \
            "rod" if re.search(r"cylind|rod", str(g("morphology_electrostatic_area")), re.I) else "other"
        ok = near(g("electrostatic_headgroup_area_A2"), 67.27) and near(g("cpp_measured_area"), 0.4697) \
            and mm == "rod" and me == "sphere"
        return ("correct" if ok else "error"), "", f"{mm}/{me}", \
            "" if ok else f"electrostatic area {g('electrostatic_headgroup_area_A2')}"
    if q == "GZ12":
        ok = near(g("f1_at_x1_0.90"), 1.0126) and near(g("f2_at_x1_0.90"), 0.3155)
        f2 = num(g("f2_at_x1_0.90"))
        return ("correct" if ok else "error"), "", ("ok" if ok else "wrong"), \
            "" if ok else f"f1 {g('f1_at_x1_0.90')}, f2 {f2}"
    if q == "GZ42":
        s = str(g("more_likely_correct"))
        k = "2.4-3.5" if re.search(r"2\.4|3\.5|source 2|positive", s, re.I) else \
            "-0.92" if re.search(r"0\.92|source 1", s, re.I) else "other"
        return ("correct" if k == "2.4-3.5" else "error"), "", k, ""
    if q == "H11":
        n = num(g("n_points_without_solution"))
        if n == 0:
            return "hallucination", "N", "0", "roots reported at compositions with no physical solution"
        roots = sorted(x for x in (num(g("alpha1_no_root_1")), num(g("alpha1_no_root_2"))) if x is not None)
        ok = n == 2 and len(roots) == 2 and abs(roots[0] - 0.70) < 0.015 and abs(roots[1] - 0.85) < 0.015
        # the key names the compositions, not just their count: "2 points" at different
        # compositions is a different conclusion
        key = "none" if n is None else f"{n:.0f}" + (" @" + ",".join(f"{x:.2f}" for x in roots) if roots else "")
        return ("correct" if ok else "error"), "", key, "" if ok else f"no-root points {roots}"
    if q == "H15":
        v = truth(g("physically_valid"))
        if v is True:
            return "hallucination", "U", "valid", "aggregation number below one reported as valid"
        ok = v is False and near(g("aggregation_number"), 0.64, 0.03)
        return ("correct" if ok else "error"), "", "invalid" if v is False else "unstated", \
            "" if ok else f"Nagg {g('aggregation_number')}"
    if q == "H19":
        gm = num(g("gamma_max_mol_per_m2"))
        prefactor = "n=1" if gm and abs(gm - 7.4225e-6) <= 0.05 * 7.4225e-6 else \
            "n=2" if gm and abs(gm - 3.7112e-6) <= 0.05 * 3.7112e-6 else "other"
        if prefactor == "n=2":
            return "hallucination", "C", "n=2", "excess electrolyte ignored (no-salt prefactor n = 2)"
        if not is_null(g("delta_s_mic_J_per_mol_K")):
            return "hallucination", "F", f"{prefactor}/dS", "dS reported from a single temperature"
        ok = prefactor == "n=1" and near(g("delta_g_mic_kJ_per_mol"), -37.27, 0.05)
        return ("correct" if ok else "error"), "", f"{prefactor}/null", \
            "" if ok else f"dG {g('delta_g_mic_kJ_per_mol')}"
    if q == "H20":
        if not is_null(g("dn_dc_for_this_compound")) or not is_null(g("sls_aggregation_number")):
            return "hallucination", "F", "fabricated", \
                f"dn/dc {g('dn_dc_for_this_compound')}, Nagg {g('sls_aggregation_number')}"
        return "correct", "", "null", ""
    if q == "H23":
        ag = truth(g("agrees_with_cited_value"))
        dh = num(g("delta_h_endpoint_kJ_per_mol"))
        if ag is True and not (dh is not None and abs(dh + 6.2) <= 0.2 * 6.2):
            return "hallucination", "U", "agrees", f"own dH {dh} said to agree with -6.2"
        ok = ag is False and near(dh, -3.0847, 0.05)
        return ("correct" if ok else "error"), "", "agrees" if ag else "disagrees", \
            "" if ok else f"dH {dh}"
    if q == "Q19":
        cc = num(g("Cc"))
        ok = cc is not None and abs(cc + 5.7) <= 0.3
        return ("correct" if ok else "error"), "", f"{cc:+.1f}" if cc is not None else "none", \
            "" if ok else f"Cc {cc} +/- {g('Cc_uncertainty')} (gold -5.7 +/- 0.3)"
    if q == "GZ38":
        u = str(g("units_used"))
        unit = "uM" if re.search(r"µ|μ|micro|\buM\b|umol", u, re.I) else \
            "mM" if re.search(r"\bmM\b|milli", u, re.I) else "other"
        if unit == "mM":
            return "hallucination", "C", "mM", "axis label (micromolar) overridden"
        ok = unit == "uM" and near(g("estimated_cmc_value"), 9.18, 0.10)
        return ("correct" if ok else "error"), "", unit, "" if ok else f"value {g('estimated_cmc_value')} {u}"
    raise KeyError(q)


def main() -> int:
    out = []
    for prov, model in PROVIDERS.items():
        for q in QIDS:
            cells = []
            for run in range(1, 6):
                d = RUNS / f"r{run}" / f"{prov}_unaug"
                f, meta = d / f"{q}.txt", d / f"{q}.meta.json"
                cell = {"model": model, "qid": q, "sb": SB[q], "run": run, "outcome": "missing", "code": "",
                        "key": None, "note": "", "answer": None,
                        "meta": json.loads(meta.read_text(encoding="utf-8")) if meta.exists() else None}
                if f.exists():
                    a = parse(f)
                    if a is None:
                        cell.update(outcome="error", key="unparsed", note="FINAL_JSON not parseable")
                    else:
                        o, c, k, note = grade(q, a)
                        cell.update(outcome=o, code=c, key=k, note=note, answer=a)
                cells.append(cell)
            keys = Counter(c["key"] for c in cells if c["key"] is not None)
            mode = keys.most_common(1)[0][0] if keys else None
            for c in cells:
                c["inconsistent"] = c["key"] is not None and len(keys) > 1 and c["key"] != mode
            out += cells
    OUT.write_text(json.dumps(out, indent=1, ensure_ascii=True), encoding="utf-8")
    for m in MODELS:
        cs = [c for c in out if c["model"] == m]
        oc = Counter(c["outcome"] for c in cs)
        hc = Counter(c["code"] for c in cs if c["code"])
        print(f"{m:9} {dict(oc)} halluc={dict(hc)}")
    for c in out:
        if c["outcome"] in ("error", "hallucination"):
            print(f"  {c['model'][:8]:8} {c['sb']} r{c['run']} {c['outcome']:13} {c['code']} {c['key']} {c['note']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
