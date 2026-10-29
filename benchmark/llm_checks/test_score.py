#!/usr/bin/env python3
"""Kit self-test (brief-2c/2d key sets + scoring kinds). Fixtures per case x stage:
correct (truth values), paper-repeating (printed values), malformed; plus targeted kind tests:
LC1 interval (sound alternative fit passes, out-of-band fails), LC4 (finite-difference route
passes; printed dH rows fail). All expectations verified, nothing waived."""
import json, subprocess, tempfile, sys
from pathlib import Path
NL = "\n"
HERE = Path(__file__).resolve().parent
CASES = ["LC1","LC2","LC3","LC4","LC5"]
PAPER = {
 "LC1": {"dG_298_kJ_mol": 44.8, "dH_298_kJ_mol": -20.6, "dCp_kJ_mol_K": -2.7, "printed_result_follows": "yes"},
 "LC2": {"x1_P184_P188": 0.40, "beta_P184_P188": 3.68, "x1_P184_P407": 0.39, "beta_P184_P407": 3.23,
         "x1_P188_P407": 0.35, "beta_P188_P407": 1.72, "printed_result_follows": "yes"},
 "LC3": {"x1": 0.84, "beta": 21.33, "printed_result_follows": "yes"},
 "LC4": {"dG_293p15_kJ_mol": -32.4, "dH_303p15_kJ_mol": -5.0, "dH_313p15_kJ_mol": -5.2,
         "dH_323p15_kJ_mol": -6.2, "more_negative_than_CPC": "no", "printed_result_follows": "yes"},
 "LC5": {"A_min_nm2": 0.0723, "CPP": 269.809, "printed_result_follows": "yes"},
}
EXTRA = [  # (case, stage, name, answer dict, min number of FAILs expected)
 ("LC1","S1","altfit", {"dG_298_kJ_mol": -44.8, "dH_298_kJ_mol": -60.0, "dCp_kJ_mol_K": -5.0}, 0),
 ("LC1","S2","outband", {"dG_298_kJ_mol": -44.8, "dH_298_kJ_mol": -100.0, "dCp_kJ_mol_K": 0.0,
                         "printed_result_follows": "cannot_determine"}, 1),
 ("LC4","S1","findiff", {"dG_293p15_kJ_mol": -40.9, "dH_303p15_kJ_mol": 61.0,
                         "dH_313p15_kJ_mol": 117.0, "dH_323p15_kJ_mol": 88.0}, 0),
]
fails = []
with tempfile.TemporaryDirectory() as td:
    adir = Path(td)/"answers"
    def put(run_id, obj):
        r = adir/run_id; r.mkdir(parents=True)
        (r/"answer.txt").write_text("work..." + NL + "FINAL_JSON: " + json.dumps(obj) + NL)
    for case in CASES:
        truth = json.loads((HERE/"truth"/(case + ".json")).read_text())
        for stage in ["S1","S2"]:
            keys = truth["stages"][stage]
            put(case + "-" + stage + "-fix", {k: v["value"] for k, v in keys.items()})
            put(case + "-" + stage + "-paper", {k: PAPER[case].get(k, "n/a") for k in keys})
            r3 = adir/(case + "-" + stage + "-bad"); r3.mkdir(parents=True)
            (r3/"answer.txt").write_text("FINAL_JSON: {broken" + NL)
    for case, stage, name, obj, _ in EXTRA:
        put(case + "-" + stage + "-" + name, obj)
    # brief-28 section 1 parsing normalizations: curly quotes, booleans, leading token, block form
    import json as _j
    def rawobj(case, stage, over=None):
        keys = json.loads((HERE/"truth"/(case + ".json")).read_text())["stages"][stage]
        d = {k: v["value"] for k, v in keys.items()}
        d.update(over or {})
        return d
    o = rawobj("LC3","S1"); txt = "FINAL_JSON: " + _j.dumps(o).replace('"', "\u201c", 1)
    curly = ("{\u201cx1\u201d: \u201cno_solution\u201d, \u201cbeta\u201d: \u201cno_solution\u201d}")
    r = adir/"LC3-S1-raw-curly"; r.mkdir(parents=True)
    (r/"answer.txt").write_text("prose" + NL + "FINAL_JSON: " + curly + NL)
    o2 = rawobj("LC3","S2", {"printed_result_follows": False})   # JSON booleans via dumps: True->true
    r = adir/"LC3-S2-raw-bool"; r.mkdir(parents=True)
    (r/"answer.txt").write_text("FINAL_JSON: " + _j.dumps(o2) + NL)
    r = adir/"LC1-S2-raw-leadtok"; r.mkdir(parents=True)
    lt = rawobj("LC1","S2"); lt["printed_result_follows"] = "cannot_determine - the 95 % interval spans both signs"
    (r/"answer.txt").write_text("FINAL_JSON: " + _j.dumps(lt) + NL)
    r = adir/"LC1-S2-raw-leadno"; r.mkdir(parents=True)
    ln = rawobj("LC1","S2"); ln["printed_result_follows"] = "no -- dG(298) reproduces"
    (r/"answer.txt").write_text("FINAL_JSON: " + _j.dumps(ln) + NL)
    r = adir/"LC2-S2-raw-block"; r.mkdir(parents=True)
    blk = "FINAL_JSON\n\ue056\njson\n" + _j.dumps(rawobj("LC2","S2", {"printed_result_follows": False})) + "\n\nFINAL_JSON: " + _j.dumps(rawobj("LC2","S2", {"printed_result_follows": False})) + NL
    (r/"answer.txt").write_text(blk)
    out = subprocess.run([sys.executable, str(HERE/"score.py"), str(adir), str(Path(td)/"s.tsv")],
                         capture_output=True, text=True, cwd=str(HERE))
    lines = out.stdout.splitlines()
    def find(sub): return [l for l in lines if sub in l]
    for case in CASES:
        for stage in ["S1","S2"]:
            mine = find(case + "-" + stage + "-fix:")
            if len(mine) < 2 or any(" FAIL" in l for l in mine): fails.append("correct " + case + "-" + stage + ": " + str(mine))
            papl = find(case + "-" + stage + "-paper:")
            nfail = sum(1 for l in papl if " FAIL" in l)
            if case == "LC1" and stage == "S1":
                # the paper prints dG(298) without its sign; dH and dCp lie inside the MC intervals
                if nfail != 1 or not any("dG_298_kJ_mol FAIL" in l for l in papl):
                    fails.append("LC1-S1 paper fixture should fail on the unsigned dG only :: " + str(papl))
            elif nfail == 0: fails.append("paper fixture " + case + "-" + stage + " all-PASS")
            if stage == "S2" and not any("printed_result_follows FAIL" in l for l in papl):
                fails.append("paper fixture " + case + "-S2 missing follows-FAIL")
            if case == "LC4" and stage == "S2":
                for kk in ["dH_303p15", "dH_313p15", "dH_323p15"]:
                    if not any(kk in l and "FAIL" in l for l in papl): fails.append("LC4 printed dH " + kk + " not failing")
            if not find(case + "-" + stage + "-bad: parsed=no"): fails.append("malformed " + case + "-" + stage + " not reported")
    for case, stage, name, obj, minf in EXTRA:
        rows = find(case + "-" + stage + "-" + name + ":")
        nf = sum(1 for l in rows if " FAIL" in l)
        if len(rows) < 2: fails.append("extra fixture " + name + " not scored")
        if nf < minf: fails.append("extra fixture " + name + ": expected >=" + str(minf) + " fails, got " + str(nf) + " :: " + str(rows))
    for rid, expect_fail in [("LC3-S1-raw-curly",0),("LC3-S2-raw-bool",0),("LC1-S2-raw-leadtok",0),
                             ("LC1-S2-raw-leadno",1),("LC2-S2-raw-block",0)]:
        rows = find(rid + ":")
        nf = sum(1 for l in rows if " FAIL" in l)
        if len(rows) < 2 or nf != expect_fail:
            fails.append("raw fixture " + rid + ": fails=" + str(nf) + " expected " + str(expect_fail) + " :: " + str(rows[:2]))
print("SELF-TEST:", "PASS" if not fails else "FAIL")
[print("  -", f) for f in fails[:10]]
sys.exit(1 if fails else 0)
