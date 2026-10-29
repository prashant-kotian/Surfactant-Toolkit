#!/usr/bin/env python3
"""compute_p2_p6.py -- recompute the application reference numbers for P2 (sub-CMC census)
and P6 (micellization free-energy conventions) on the released toolkit.

The 2026-10-05 reference results (dev repo origin/main 6fa609a,
application_reference/results_B_C_5_6_7.json) were computed before the ionic
counterion-factor correction (release 01a396c). This script recomputes both
problems through the release toolkit's own functions, compares against the old
reference values, and writes values_p2_p6.json next to itself.

    python compute_p2_p6.py [path/to/Surfactant-Toolkit]
"""
from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
TK = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE.parents[1]
sys.path.insert(0, str(TK / "src"))

from surfactantkit.kinetics import mass_action_equilibrium  # noqa: E402
from surfactantkit.thermodynamics import (  # noqa: E402
    R_GAS, cmc_to_mole_fraction, gibbs_free_energy_micellization, mc_aggregation_distribution)

T = 298.15
OLD = json.loads((HERE / "reference_2026-10-05_results_B_C_5_6_7.json").read_text(encoding="utf-8-sig"))


def census(cmc_M: float, n_max: int = 100) -> list[dict]:
    """Stepwise mass-action census with K = 1/CMC (the 2026-10-05 protocol)."""
    K = 1.0 / cmc_M
    rows = []
    for r in (0.05, 0.10, 0.12, 0.25, 0.50, 1.00, 2.00):
        c = r * cmc_M
        x = mass_action_equilibrium(c, K, n_max)
        mass = [(n + 1) * xn / c for n, xn in enumerate(x)]
        rows.append({"C_over_CMC": r, "C_M": c, "monomer": mass[0], "dimer": mass[1],
                     "trimer": mass[2], "n_ge_4": sum(mass[3:]),
                     "mass_balance_rel_resid": abs(sum(mass) - 1.0)})
    return rows


def compare_census(new: list[dict], old_key: str) -> float:
    old = OLD[old_key]["rows"]
    worst = 0.0
    for a, b in zip(new, old):
        for f, g in (("monomer", "massfrac_monomer"), ("dimer", "massfrac_dimer"),
                     ("trimer", "massfrac_trimer"), ("n_ge_4", "massfrac_n_ge_4")):
            worst = max(worst, abs(a[f] - b[g]))
    return worst


def dg(cmc_M: float, factor: float, molar: bool = False, T: float = T) -> float:
    x = cmc_M if molar else cmc_to_mole_fraction(cmc_M)
    if molar:  # ln of the molar CMC (standard state 1 M), outside the function's mole-fraction guard
        return factor * R_GAS * T * math.log(x) / 1000.0
    return gibbs_free_energy_micellization(x, T, factor)


def main() -> int:
    head = subprocess.run(["git", "-C", str(TK), "rev-parse", "--short", "HEAD"],
                          capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", str(TK), "status", "--porcelain", "--", "src"],
                           capture_output=True, text=True).stdout.strip()
    out: dict = {"toolkit": TK.name, "toolkit_commit": head, "src_dirty": bool(dirty), "T_K": T}

    # ---- P2: sub-CMC census (counterion factor does not enter) ----
    c12e8 = census(1.0e-4)
    sds = census(0.0083)
    mc = mc_aggregation_distribution(seed=0)
    out["P2"] = {
        "C12E8": {"cmc_M": 1.0e-4, "source": "C12E8 CMC 1.0e-4 M (seven-gold set)", "rows": c12e8,
                  "max_abs_diff_vs_2026-10-05": compare_census(c12e8, "CS-B_C12E8")},
        "SDS": {"cmc_M": 0.0083, "source": "SDS cmc0 8.3 mM (Bales 1998); no counterion equilibrium",
                "rows": sds, "max_abs_diff_vs_2026-10-05": compare_census(sds, "CS-B_SDS")},
        "MC": {"seed": 0, "free_monomer_fraction": mc.free_monomer_fraction if hasattr(mc, "free_monomer_fraction") else None,
               "converged": getattr(mc, "converged", None),
               "old_free_monomer_fraction": OLD["CS-B_MC"]["free_monomer_fraction"]},
    }

    # ---- P6: free energy of micellization under each convention ----
    beta_ddgc, cmc_ddgc = 0.74, 6.2e-3          # Bouchal et al. 2016, Table 2 (DDGC)
    alpha_sds, ua, cmc_sds = 0.272, 0.017, 0.0083  # Bales 2001 EPR alpha; Bales 1998 cmc0
    beta_sds = 1.0 - alpha_sds
    ddgc_correct = dg(cmc_ddgc, 1.0 + beta_ddgc, T=298.0)  # Bouchal et al. report 298 K
    ddgc_inverse = dg(cmc_ddgc, 2.0 - beta_ddgc, T=298.0)
    sds_rows = {
        "pseudo-phase, mole fraction, (1 + beta) [correct]": dg(cmc_sds, 1.0 + beta_sds),
        "pseudo-phase, mole fraction, (2 - beta) = (1 + alpha) [inverse]": dg(cmc_sds, 2.0 - beta_sds),
        "mole fraction, factor 1 (counterions ignored)": dg(cmc_sds, 1.0),
        "mole fraction, factor 2 (full binding)": dg(cmc_sds, 2.0),
        "molar standard state, (1 + beta)": dg(cmc_sds, 1.0 + beta_sds, molar=True),
    }
    vals = list(sds_rows.values())
    unc = abs(R_GAS * T * math.log(cmc_to_mole_fraction(cmc_sds)) / 1000.0) * ua
    out["P6"] = {
        "DDGC": {"T_K": 298.0, "cmc_M": cmc_ddgc, "beta": beta_ddgc, "published_kJ": -28.3, "published_unc": 0.9,
                 "correct_(1+beta)": ddgc_correct, "inverse_(2-beta)": ddgc_inverse,
                 "old_reference": OLD["CS-6_conventions"]["DDGC"]["our"]},
        "SDS": {"cmc_M": cmc_sds, "alpha": alpha_sds, "alpha_unc": ua, "rows": sds_rows,
                "unc_from_alpha_kJ": unc, "spread_kJ": max(vals) - min(vals),
                "old_reference_rows": OLD["CS-6_conventions"]["SDS_spread"]},
    }
    (HERE / "values_p2_p6.json").write_text(json.dumps(out, indent=2), encoding="utf-8")

    p2, p6 = out["P2"], out["P6"]
    print(f"toolkit {head} src_dirty={bool(dirty)}")
    print("P2 C12E8 0.12xCMC:", {k: round(v, 4) for k, v in c12e8[2].items() if k in ("monomer", "dimer", "trimer", "n_ge_4")})
    print("P2 max |diff| vs 2026-10-05: C12E8 %.2e, SDS %.2e" % (p2["C12E8"]["max_abs_diff_vs_2026-10-05"], p2["SDS"]["max_abs_diff_vs_2026-10-05"]))
    print("P2 MC free monomer:", p2["MC"]["free_monomer_fraction"], "old", p2["MC"]["old_free_monomer_fraction"], "converged", p2["MC"]["converged"])
    print("P6 DDGC: (1+beta) %.2f | (2-beta) %.2f | published -28.3 | old ref %.2f" % (ddgc_correct, ddgc_inverse, p6["DDGC"]["old_reference"]))
    for k, v in sds_rows.items():
        print("P6 SDS %-62s %.2f" % (k, v))
    print("P6 SDS +/- from alpha: %.2f ; spread %.2f" % (unc, p6["SDS"]["spread_kJ"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
