"""Case 2. Poloxamer mixed micelles: sign and size of the Rubingh interaction parameter.

J. Sedlarikova, M. Janalikova, P. Egner, P. Pleva, ACS Omega 9 (2024) 23209-23219,
doi:10.1021/acsomega.3c08917.
Inputs: Table 1 (molecular weights) and Table 2 (CMC with SD, beta, xM); mixtures prepared 1:1 by
weight (Methods). Claims examined: beta = +3.68, +3.23, +1.72, described in the text as "negative
deviations from ideal behavior, which indicates a mild synergic effect"; xM close to x = 0.5.

 1. Reproduce the ideal (Clint) CMC and note that every mixed CMC lies below both pure CMCs.
 2. Toolkit Rubingh solution (solve_rubingh_x, rubingh_beta) for x = 0.5 (the paper's reading) and
    for the molar fraction of a 1:1 mass mixture; Monte Carlo over the printed CMC SDs.
 3. Independent bisection of Rubingh's equation in numpy.

    python case2_poloxamer_rubingh.py [path/to/toolkit]   -> case2_results.json
"""
from __future__ import annotations

import math

import numpy as np

from common import write
from surfactantkit.mixed_micelle import clint_ideal_cmc, rubingh_beta, solve_rubingh_x

MW = {"P184": 2900, "P188": 8400, "P407": 12600}
CMC = {"P184": (0.55, 0.01), "P188": (0.17, 0.05), "P407": (0.10, 0.01)}      # mM, SD
MIX = {("P184", "P188"): ((0.04, 0.01), 3.68, 0.40),
       ("P184", "P407"): ((0.06, 0.01), 3.23, 0.39),
       ("P188", "P407"): ((0.07, 0.01), 1.72, 0.35)}
N_MC = 4000


def np_rubingh(a1, cm, c1, c2):
    f = lambda x: x * x * math.log(a1 * cm / (x * c1)) - (1 - x) ** 2 * math.log((1 - a1) * cm / ((1 - x) * c2))
    lo, hi = 1e-9, 1 - 1e-9
    if f(lo) * f(hi) > 0:
        return None, None
    for _ in range(200):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if f(lo) * f(mid) > 0 else (lo, mid)
    x = (lo + hi) / 2
    return x, math.log(a1 * cm / (x * c1)) / (1 - x) ** 2


def main():
    out = {"paper": "Sedlarikova et al., ACS Omega 9 (2024) 23209", "doi": "10.1021/acsomega.3c08917", "mixtures": {}}
    rng = np.random.default_rng(20261027)
    for (a, b), ((cm, scm), b_pr, x_pr) in MIX.items():
        c1, c2 = CMC[a][0], CMC[b][0]
        a_mass = (1 / MW[a]) / (1 / MW[a] + 1 / MW[b])
        row = {"cmc_mix_mM": cm, "pure_cmc_mM": [c1, c2], "below_both_pure": cm < min(c1, c2),
               "printed": {"beta": b_pr, "xM": x_pr}}
        for label, a1 in (("x_0.5", 0.5), ("one_to_one_by_mass", a_mass)):
            x = solve_rubingh_x(a1, cm, c1, c2)
            beta = rubingh_beta(x, a1, cm, c1)
            xn, bn = np_rubingh(a1, cm, c1, c2)
            bs = []
            for _ in range(N_MC):
                d1, d2, dm = (c1 + CMC[a][1] * rng.standard_normal(), c2 + CMC[b][1] * rng.standard_normal(),
                              cm + scm * rng.standard_normal())
                if min(d1, d2, dm) <= 0:
                    continue
                xx = solve_rubingh_x(a1, dm, d1, d2)
                if xx is not None:
                    bs.append(rubingh_beta(xx, a1, dm, d1))
            row[label] = {"alpha1": round(a1, 3), "ideal_cmc_mM": round(clint_ideal_cmc(a1, c1, c2), 3),
                          "x1": round(x, 3), "beta": round(beta, 2), "numpy_beta": round(bn, 2),
                          "beta_95": np.percentile(bs, [2.5, 97.5]).round(2).tolist()}
        out["mixtures"][f"{a}/{b}"] = row
    out["findings"] = {
        "all_mixed_cmc_below_both_pure": all(r["below_both_pure"] for r in out["mixtures"].values()),
        "recomputed_beta_negative_with_95pct_below_zero": all(r["one_to_one_by_mass"]["beta_95"][1] < 0
                                                              for r in out["mixtures"].values()),
    }
    return write("case2", out)


if __name__ == "__main__":
    import json
    print(json.dumps(main(), indent=1))
