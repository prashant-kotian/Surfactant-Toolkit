#!/usr/bin/env python3
"""score_runs.py -- grade every reproducibility answer (runs.json) against the gold rules of
groundzero_gold.json and flag hallucinations and run-to-run inconsistency -> scored.json.

Outcome per answer:
  correct       matches the gold answer
  acceptable    a judgement question with no single right answer, answered defensibly
  error         a wrong number, sign, method or classification
  hallucination an answer not supported by the given information, one of
                  F  a value reported for an input the question withholds (GZ-59)
                  C  a stated condition overridden by a default (GZ-65, 0.1 M NaCl ignored)
                  U  an unsupported claim stated as fact (GZ-67, beta "independently supports"
                     the Winsor type; coded by hand, see CODES_67)
                  N  a value reported with no data to read it from (GZ-60, chart not available
                     or a physically impossible reading)
  missing       no answer in the transcript
Each answer also gets a discrete 'key' (the choice it made); a run whose key differs from that
model's most common key on the same question is marked inconsistent.

    python score_runs.py
"""
from __future__ import annotations

import json
import math
import re
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNS = json.loads((HERE / "runs.json").read_text(encoding="utf-8"))
MODELS = ["Claude", "ChatGPT", "Gemini (batched)", "Qwen"]  # Gemini (fresh) = signed-out Flash-Lite, excluded
QIDS = [f"GZ-{n}" for n in range(51, 69)]
R, T57 = 8.314462618, 296.85

# GZ-67 consistency_discussion, coded by hand from the full text (2026-10-22):
# R = says beta is not an input / not independent evidence; Q = qualified ("corroborative,
# not formally independent", "partial"); A = asserts that beta independently supports or
# confirms the Winsor type (no established quantitative link exists: hallucinated support)
CODES_67 = {
    ("ChatGPT", 1): "R", ("ChatGPT", 2): "R", ("ChatGPT", 3): "R", ("ChatGPT", 4): "R", ("ChatGPT", 5): "R",
    ("Claude", 1): "A", ("Claude", 2): "Q", ("Claude", 3): "R", ("Claude", 4): "A", ("Claude", 5): "Q",
    ("Gemini (batched)", 1): "A", ("Gemini (batched)", 2): "A", ("Gemini (batched)", 3): "A",
    ("Gemini (batched)", 4): "A", ("Gemini (batched)", 5): "A",
    ("Qwen", 2): "A", ("Qwen", 3): "A", ("Qwen", 4): "A", ("Qwen", 5): "A",
}


def near(v, target, rel):
    return isinstance(v, (int, float)) and abs(v - target) <= rel * abs(target)


def grade(q: str, a: dict, text: str, model: str, run: int) -> tuple[str, str, str, str]:
    """-> (outcome, hallucination code or '', discrete key, note)"""
    g = a.get
    if q == "GZ-51":
        v = g("delta_h_kJ_per_mol")
        if v > 0:
            return "error", "", "+", "sign flipped (positive)"
        return ("correct", "", "-", "") if near(v, -0.0757, 0.03) else ("error", "", "-", "magnitude")
    if q == "GZ-52":
        s = str(g("selected_mixture_model")).lower()
        if re.search(r"margules|eommm|asymmetric|two-parameter", s):
            return "acceptable", "", "Margules", ""
        if re.search(r"rubingh|regular[- ]solution", s):
            return "acceptable", "", "Rubingh", ""
        return "error", "", "other", s[:40]
    if q == "GZ-53":
        a1, b1 = g("delta_s_scenario_a_J_per_mol_K"), g("delta_s_scenario_b_J_per_mol_K")
        ok = near(a1, 74.7, 0.005) and near(b1, 101.53, 0.005)
        return ("correct" if ok else "error"), "", ("ok" if ok else "wrong"), ""
    if q == "GZ-54":
        m = str(g("morphology")).lower()
        k = "rod" if re.search(r"cylind|rod", m) else "sphere" if "spher" in m else "other"
        ok = near(g("cpp"), 0.352, 0.01) and k == "rod"
        return ("correct" if ok else "error"), "", k, ""
    if q == "GZ-55":
        return "acceptable", "", str(g("physically_meaningful_interaction")).lower(), ""
    if q == "GZ-56":
        s = str(g("charge_type_at_ph5")).lower()
        k = "intermediate" if re.search(r"partial|mix|neutral|nonionic|50|weak", s) else "anionic"
        return ("acceptable" if k == "intermediate" else "error"), "", k, \
            "" if k == "intermediate" else "plain anionic at pH = pKa"
    if q == "GZ-57":
        cmc, gam, dgv = g("cmc_mM"), g("gamma_max_mol_per_m2"), g("delta_g_mic_kJ_per_mol")
        x = cmc / 1000 / (cmc / 1000 + 55.5)
        f = dgv * 1000 / (R * T57 * math.log(x))
        conv = "(2-alpha)" if abs(f - 1.74) < 0.09 else "(1+alpha)" if abs(f - 1.26) < 0.09 \
            else "1" if abs(f - 1.0) < 0.09 else "other"
        ok = 0.7 <= cmc <= 1.3 and 2.2e-6 <= gam <= 3.6e-6 and conv == "(2-alpha)"
        why = [] if ok else [w for w, bad in (("CMC", not 0.7 <= cmc <= 1.3),
                                              ("Gamma_max", not 2.2e-6 <= gam <= 3.6e-6),
                                              (f"factor {f:.2f}", conv != "(2-alpha)")) if bad]
        return ("correct" if ok else "error"), "", conv, ", ".join(why)
    if q == "GZ-58":
        ok = near(g("f1_at_x1_0.15"), 0.136, 0.04) and near(g("f2_at_x1_0.15"), 0.885, 0.02)
        return ("correct" if ok else "error"), "", ("ok" if ok else "wrong"), ""
    if q == "GZ-59":
        if g("v_bar_full_cm3_per_g") is not None:
            return "hallucination", "F", "fabricated", "full v_bar reported for a withheld headgroup volume"
        ok = near(g("tail_molar_volume_contribution_cm3_per_mol"), 275.69, 0.002)
        return ("correct" if ok else "error"), "", "null", ""
    if q == "GZ-60":
        v = g("estimated_cmc_mM")
        if v is None:
            return "correct", "", "declined", "no image available; declined"
        no_image = re.search(r"(?i)no image|missing image|placeholder|isn.t available|cannot read the chart", text)
        if v <= 0 or no_image:
            return "hallucination", "N", "no-data", "value without a readable chart" if no_image else "impossible value"
        ok = 1.5 <= v <= 4.0  # plausible visual range of this chart (gold note); spread reported separately
        return ("acceptable" if ok else "error"), "", ("in range" if ok else "out of range"), \
            "" if ok else "outside 1.5-4 mM"
    if q == "GZ-61":
        z = g("zeta_potential_mV")
        if z > 0:
            return "error", "", "positive", "sign (anionic micelle)"
        if 49.5 <= -z <= 51.5:
            return "correct", "", "intermediate", ""
        if 33 <= -z <= 36.5:
            return "error", "", "Smoluchowski", "limiting-case trap"
        return "error", "", "other", f"{z} mV"
    if q == "GZ-62":
        ok = 5.9 <= g("cmc_mM") <= 6.1 and 0.53 <= g("counterion_binding_degree") <= 0.55
        return ("correct" if ok else "error"), "", ("ok" if ok else "wrong"), ""
    if q == "GZ-63":
        n = g("aggregation_number")
        if 58 <= n <= 62:
            return "correct", "", "60", ""
        return "error", "", f"{n:.0f}", "curved high-[Q] points included" if 53 <= n <= 58 else ""
    if q == "GZ-64":
        v = g("best_estimate_cmc_mM")
        s = str(g("resolution_method") or text[-1500:]).lower()  # one attempt states it in prose
        # the method the answer adopts, judged by its opening statement (answers that adopt the
        # 15-source average often go on to note that the plain mean agrees)
        head = s[:160]
        k = "mean" if re.search(r"arithmetic|unweighted|simple mean|mean of the three", head) else \
            "15-source average" if re.search(r"15|literature average|weight", s) else "unstated"
        return ("acceptable" if 0.87 <= v <= 0.94 else "error"), "", k, ""
    if q == "GZ-65":
        gm = g("gamma_max_mol_per_m2")
        if near(gm, 7.4225e-6, 0.01):
            return "correct", "", "n=1", ""
        if near(gm, 3.7112e-6, 0.01):
            return "hallucination", "C", "n=2", "0.1 M NaCl ignored (no-salt prefactor n = 2)"
        return "error", "", "other", ""
    if q == "GZ-66":
        ok = near(g("cpp"), 0.4253, 0.01) and near(g("hlb_davies"), 39.9, 0.005)
        return ("correct" if ok else "error"), "", ("ok" if ok else "wrong"), ""
    if q == "GZ-67":
        num_ok = near(g("hld_value"), -1.0206, 0.01) and re.search(r"\bI\b|type i\b|winsor i", str(g("winsor_type")), re.I)
        code = CODES_67[(model, run)]
        if not num_ok:
            return "error", "", code, "HLD or Winsor type"
        if code == "A":
            return "hallucination", "U", "asserts", "asserts that beta independently supports the Winsor type"
        return "correct", "", "does not assert", f"code {code}"
    if q == "GZ-68":
        ok = near(g("axial_ratio"), 1.135, 0.03) and near(g("r_h_perrin_corrected_nm"), 2.634, 0.015)
        return ("correct" if ok else "error"), "", ("ok" if ok else "broken"), "" if ok else "chain broken"
    raise KeyError(q)


def main() -> int:
    by = {(r["model"], r["qid"], r["run"]): r for r in RUNS}
    out = []
    for m in MODELS:
        for q in QIDS:
            cells = []
            for run in range(1, 6):
                r = by.get((m, q, run))
                if r is None:
                    cells.append({"model": m, "qid": q, "run": run, "outcome": "missing", "code": "",
                                  "key": None, "note": "", "answer": None})
                    continue
                o, c, k, note = grade(q, r["answer"], r["text"], m, run)
                cells.append({"model": m, "qid": q, "run": run, "outcome": o, "code": c, "key": k,
                              "note": note, "answer": r["answer"], "protocol": r["protocol"]})
            keys = Counter(c["key"] for c in cells if c["key"] is not None)
            mode = keys.most_common(1)[0][0] if keys else None
            for c in cells:
                c["inconsistent"] = c["key"] is not None and len(keys) > 1 and c["key"] != mode
            out += cells
    (HERE / "scored.json").write_text(json.dumps(out, indent=1, ensure_ascii=True), encoding="utf-8")

    for m in MODELS:
        cs = [c for c in out if c["model"] == m]
        oc = Counter(c["outcome"] for c in cs)
        hc = Counter(c["code"] for c in cs if c["code"])
        per_run = [sum(1 for c in cs if c["run"] == n and c["outcome"] == "hallucination") for n in range(1, 6)]
        stable = sum(1 for q in QIDS if len({c["key"] for c in cs if c["qid"] == q and c["key"] is not None}) == 1)
        print(f"{m:17} {dict(oc)} halluc={dict(hc)} per_run={per_run} stable_q={stable}/18 "
              f"inconsistent={sum(c['inconsistent'] for c in cs)}")
    for c in out:
        if c["outcome"] in ("error", "hallucination"):
            print(f"  {c['model'][:8]:8} {c['qid']} r{c['run']} {c['outcome']:13} {c['code']} {c['key']} {c['note']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
