#!/usr/bin/env python3
"""build_question_set.py — assemble the single 200-question benchmark.

SurfBench-200 is ONE benchmark. It was authored as four working sets
(groundzero, hard50, extrahard, orig-hard50) purely for construction
convenience; this script merges them into one id space (SB-001..SB-200)
with a per-question classification on two axes:

  type        : conceptual | computational | analytical
  subconcept  : surfactant-science topic key (drives the class-matrix table)

Classification is rule-based over each question's gold metadata and question
text (field names, declared failure axes, task verbs), with a curated
override table (question_set_overrides.py). Deterministic; no network.
Run after the source-set golds exist:  python3 benchmark/build_question_set.py
"""
from __future__ import annotations

import json
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
SS = HERE / "source_sets"
BENCH_NAME = "SurfBench-200"      # single constant to rename the benchmark
BENCH_N = 200

SUBCONCEPT_RULES = [  # first match wins; applied to gold field names + question body
    (r"solubil|extraction", "solubilization"),
    (r"spreading|contact_angle|wetting|work_of_adhesion", "surface_wetting"),
    (r"foam|thin_film|drainage|surface_visc|marangoni|bending|Cc", "foaming_films"),
    (r"debye_length|screening|huckel|zeta|electrophor|poisson|ionic_strength|kappa_a|double_layer|surface_charge|obrien", "colloid_electrostatics"),
    (r"szyszkowski|gibbs|frumkin|gamma_max|gamma_inf|gamma_at|a_min|adsorb|tensiometr|dynamic_surface_tension", "adsorption"),
    (r"cpp|packing|tail.*volume|critical_length|axial|ellipsoid|cylinder|morphology|partial_molar|v_bar", "packing_morphology"),
    (r"n_agg|aggregation_number|micelle.*molar_mass|light.?scattering|sals|\bsls\b|dls|associating cluster|metropolis", "micelle_structure"),
    (r"kinetic|gillespie|diffus|relax|quench|micelliz|stopped|master.?equation|transport|breakthrough", "kinetics_transport"),
    (r"eommm|clint|rubingh|rubingox|regular_solution|activity_coeff|beta|counterion|f1_|f2_|x1_|cmix|c1_c2|mixed|cmc1|cmc2|alpha1", "mixed_micelles"),
    (r"delta_g|delta_h|delta_s|delta_cp|entropy|enthalp|hain.?marks|vant_hoff|thermodynam", "micelle_thermodynamics"),
    (r"hld|emulsion|emulsif|phase_inversion|pit|cloud_point|hlb", "emulsion_hld"),
    (r"cmc|critical.?micelle|conductivity", "cmc_measurement"),
    (r"uncertaint|monte.?carlo|residual|outlier|bootstrap|regress|slope|intercept|tolerance", "data_fitting"),
    (r"qspr|retention|logp|hplc|k_oc", "structure_property"),
    (r"viscos|rheolog|gel", "rheology"),
    (r"safety|toxic|biodegrad|environ|ecotox", "safety_environment"),
    (r"name|identify|classif|charge_type|structural_family|nomenclature|smiles|terminology|define|convention|acronym|speciation", "identification_terminology"),
    (r"conflict|adjudicat|audit|forensic|refut|recompute|cross.?source", "analysis_synthesis"),
    (r"select|choose|recommend|compare", "analysis_synthesis"),
]
ANALYTICAL = re.compile(r"\b(select|choose|recommend|justify|compare|reconcil|audit|determine which|which (model|method|structure)|should|refuse|flag|judge|decide|best|explain why)\b", re.I)
COMPUTATIONAL = re.compile(r"\b(compute|calculate|estimate|determine|fit|derive|regress|propagate|evaluate|round|convert|recompute|count|plot)\b", re.I)
CONCEPTUAL = re.compile(r"\b(name|classify|identify|define|list|cite|what (is|are)|which (surfactant|compound|charge|family)|give the|report the (name|term|charge))\b", re.I)


def question_text(body: str) -> str:
    """Strip reference-only notes / paste markers; keep what the model saw."""
    if ">>> PASTE BELOW >>>" in body:
        body = body.split(">>> PASTE BELOW >>>", 1)[1]
    body = re.sub(r"\[reference only[^\]]*\]", " ", body)
    return body


def classify(body: str, field_names: str, gold_axes: str) -> tuple[str, str]:
    hay = field_names + " " + body[:2600]
    sub = next((name for pat, name in SUBCONCEPT_RULES if re.search(pat, hay, re.I)),
               "general_physchem")
    ax = (gold_axes + " " + body).lower()
    if "reasoning" in ax or "vision" in ax or ANALYTICAL.search(body[:2600]):
        typ = "analytical"
    elif "definition" in ax or (CONCEPTUAL.search(body) and not COMPUTATIONAL.search(body)):
        typ = "conceptual"
    elif COMPUTATIONAL.search(body) or re.search(r"\b(number|value|mM|nm|N_?[av]g)\b", field_names, re.I):
        typ = "computational"
    else:
        typ = "conceptual" if sub in ("identification_terminology", "safety_environment") else "computational"
    return typ, sub


def gz_blocks(path: Path) -> dict:
    txt = path.read_text(encoding="utf-8", errors="replace")
    out = {}
    for m in re.finditer(r"(?m)^GZ-(\d\d)$", txt):
        n = int(m.group(1))
        nxt = re.search(r"(?m)^GZ-\d\d$", txt[m.end():])
        out[n] = txt[m.end(): m.end() + (nxt.start() if nxt else 4000)]
    return out


def h50_blocks(path: Path) -> dict:
    txt = path.read_text(encoding="utf-8", errors="replace")
    ends = [(m.end(), int(m.group(1))) for m in re.finditer(r"END OF QUESTION H-(\d\d)", txt)]
    out = {}
    prev = 0
    for e, n in ends:
        out[n] = txt[prev:e]
        prev = e
    return out


def main() -> int:
    items = []
    gz_gold = json.loads((SS / "groundzero/gold_answers.json").read_text(encoding="utf-8"))
    gz_txt = gz_blocks(SS / "groundzero/question_bank_augmented.txt")
    h_gold = json.loads((SS / "hard50/gold_answers.json").read_text(encoding="utf-8"))
    h_txt = h50_blocks(SS / "hard50/question_bank.txt")
    x_gold = json.loads((SS / "extrahard/gold_answers.json").read_text(encoding="utf-8"))

    def emit(sb, alias, setname, body, fields, axes):
        body = question_text(body or "")
        typ, sub = classify(body, fields, axes)
        items.append(dict(sb=sb, alias=alias, source_set=setname, type=typ, subconcept=sub))

    for n in range(1, 51):
        rec = gz_gold.get(f"GZ-{n:02d}") or {}
        fields = " ".join(k for k, v in (rec.get("gold") or {}).items() if v is not None)
        emit(f"SB-{n:03d}", f"GZ-{n:02d}", "groundzero", gz_txt.get(n, ""),
             fields, str(rec.get("source", ""))[:400])
    for n in range(1, 51):
        rec = h_gold.get(f"H-{n:02d}") or {}
        fields = " ".join((rec.get("fields") or {}).keys())
        axes = f"{rec.get('failure_axis','')} {rec.get('family','')}"
        emit(f"SB-{n + 50:03d}", f"H-{n:02d}", "hard50", h_txt.get(n, ""), fields, axes)
    for n in range(1, 51):
        rec = x_gold.get(f"X-{n:02d}") or {}
        fields = " ".join((rec.get("fields") or {}).keys())
        emit(f"SB-{n + 100:03d}", f"X-{n:02d}", "extrahard", rec.get("body", ""),
             fields, str(rec.get("family", "")) + " " + " ".join(rec.get("levers", []) or []))
    for n in range(1, 51):
        raw = (SS / "orig_hard50" / f"Q{n:02d}.md").read_text(encoding="utf-8", errors="replace")
        m = re.search(r"(?m)^Q\d\d\s*--[\s\S]+$", raw)
        body = m.group(0) if m else raw
        emit(f"SB-{n + 150:03d}", f"Q{n:02d}", "orig_hard50", body, "", "")

    from question_set_overrides import OVERRIDES
    for it in items:
        if it["sb"] in OVERRIDES:
            it["type"], it["subconcept"] = OVERRIDES[it["sb"]]

    out = {"benchmark": BENCH_NAME, "n": len(items),
           "excluded_from_scoring": {
               "groundzero": {"image_input": ["GZ-15", "GZ-36", "GZ-37", "GZ-38", "GZ-39", "GZ-40"],
                               "retired_ambiguous": ["GZ-26"]},
               "hard50": {"image_input": ["H-36", "H-37", "H-38"]},
               "extrahard": {"pre_fix_packet_cells": ["X-19"]},
               "per_arm_corrections": "standard list above; per-arm re-runs, re-grades and exclusions since "
                                      "2026-10-20 are logged in results/results.json cells[*].corrections"},
           "items": items}
    (HERE / "question_set.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    from collections import Counter
    assert len(items) == BENCH_N
    print("types:", dict(Counter(i["type"] for i in items)))
    print("subconcepts:", dict(Counter(i["subconcept"] for i in items).most_common()))
    for it in items:
        if it["subconcept"] == "general_physchem":
            print("FALLBACK:", it["sb"], it["alias"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
