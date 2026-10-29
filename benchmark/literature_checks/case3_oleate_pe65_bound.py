"""Case 3. Sodium oleate + phospholan PE65: a mixed CMC that no pseudo-phase model allows.

V. V. Sergeev, O. V. Cheremisina, A. T. Fedorov, A. A. Gorbacheva, D. A. Balandinsky,
ACS Omega 7 (2022) 3016-3023, doi:10.1021/acsomega.1c06047.
Inputs (text and Table 3): CMC of PE65 33 mg/L, of sodium oleate 330 mg/L, of the mixture 116 mg/L;
alpha1 = 0.4 (fraction of PE65); printed x1 = 0.84, beta = 21.33. Claim examined: "antagonism in
the formation of mixed micelles is observed".

 1. Bound: at the mixed CMC each component's monomer concentration is alpha_i * C12, and in any
    pseudo-phase model (ideal, regular solution, or separate micelles) it cannot exceed that
    component's own CMC, so C12 <= min_i C_i / alpha_i. Mass concentrations suffice for this test.
 2. Toolkit Rubingh solution on the mass basis and on molar bases for PE65 molecular weights
    300-1000 g/mol (oleate 304.4).
 3. The printed (x1, beta) in both Rubingh equations.

    python case3_oleate_pe65_bound.py [path/to/toolkit]   -> case3_results.json
"""
from __future__ import annotations

import math

from common import write
from surfactantkit.mixed_micelle import clint_ideal_cmc, solve_rubingh_x

C1, C2, C12 = 33.0, 330.0, 116.0      # mg/L: PE65, sodium oleate, mixture
A1 = 0.4                              # fraction of PE65
X1_PR, BETA_PR = 0.84, 21.33
MW_OLEATE = 304.4


def main():
    out = {"paper": "Sergeev et al., ACS Omega 7 (2022) 3016", "doi": "10.1021/acsomega.1c06047"}
    out["bound"] = {"pe65_monomer_at_mixed_cmc_mg_L": A1 * C12, "pe65_cmc_mg_L": C1,
                    "max_allowed_mixed_cmc_mg_L": round(min(C1 / A1, C2 / (1 - A1)), 1),
                    "violated": A1 * C12 > C1}
    out["mass_basis"] = {"ideal_cmc_mg_L": round(clint_ideal_cmc(A1, C1, C2), 1),
                         "toolkit_rubingh_x1": solve_rubingh_x(A1, C12, C1, C2)}
    molar = {}
    for mw in (300, 400, 500, 600, 800, 1000):
        a = (A1 / mw) / (A1 / mw + (1 - A1) / MW_OLEATE)
        c12 = C12 * (A1 / mw + (1 - A1) / MW_OLEATE)
        molar[str(mw)] = {"alpha1": round(a, 3), "toolkit_rubingh_x1": solve_rubingh_x(a, c12, C1 / mw, C2 / MW_OLEATE)}
    out["molar_basis_by_pe65_mw"] = molar
    b1 = math.log(A1 * C12 / (X1_PR * C1)) / (1 - X1_PR) ** 2
    b2 = math.log((1 - A1) * C12 / ((1 - X1_PR) * C2)) / X1_PR**2
    out["printed_pair_in_rubingh"] = {"beta_from_eq1": round(b1, 2), "beta_from_eq2": round(b2, 2), "printed_beta": BETA_PR}
    out["findings"] = {"bound_violated": out["bound"]["violated"],
                       "no_rubingh_solution_any_basis": out["mass_basis"]["toolkit_rubingh_x1"] is None
                       and all(v["toolkit_rubingh_x1"] is None for v in molar.values()),
                       "printed_pair_inconsistent": abs(b1 - b2) > 1}
    return write("case3", out)


if __name__ == "__main__":
    import json
    print(json.dumps(main(), indent=1))
