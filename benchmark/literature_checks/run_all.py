"""Run the five literature checks and collect their findings.

    python run_all.py [path/to/toolkit]   -> case1..5_results.json, literature_checks_summary.json
"""
from __future__ import annotations

import importlib
import json

from common import HERE

CASES = ["case1_ctab_vant_hoff", "case2_poloxamer_rubingh", "case3_oleate_pe65_bound",
         "case4_cpc_tx114_gibbs_helmholtz", "case5_dmimcl_des_adsorption"]


def main():
    summary = {}
    for name in CASES:
        r = importlib.import_module(name).main()
        summary[name] = {"doi": r["doi"], "paper": r["paper"], "findings": r["findings"]}
    (HERE / "literature_checks_summary.json").write_text(json.dumps(summary, indent=1), encoding="utf-8")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
