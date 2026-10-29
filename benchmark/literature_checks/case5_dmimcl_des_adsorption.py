"""Case 5. [C10mim]Cl in amino-acid deep eutectic solvents: surface excess, area and packing parameter.

M. K. Banjare, B. Barman, K. Behera, J. M. Khan, R. K. Banjare, S. Pandey, K. K. Ghosh,
Heliyon 10 (2024) e35598, doi:10.1016/j.heliyon.2024.e35598 (CC BY-NC).
Inputs: Table 2 (CMC by surface tension), Table 3 (gamma_CMC, Gamma_max, A_min, pi_CMC, pC20, CPP),
n = 2 (text), eqs (6)-(8) (CPP with Tanford volume and length of the C10 chain).
Claims examined: Gamma_max 12-25 umol/m2 and A_min 0.065-0.134 nm2; CPP 255-328, from which the
paper infers a change of aggregate shape "from globule to spherical".

 1. Reproduce: A_min = 1/(N_A Gamma_max) from the printed Gamma_max.
 2. Physical floor: the cross-section of one all-trans alkyl chain is about 0.20 nm2 (Tanford volume
    over length, toolkit tanford chain volume and length), so no single-chain surfactant can occupy
    less; CPP = v/(a0 l) must lie between 0 and about 1 for any micelle or bilayer.
 3. Gamma_max recovered from the paper's own printed CMC, C20 (pC20) and pi_CMC: the Szyszkowski
    equation pi = n R T Gamma ln(1 + K c) at pi = 20 (C20) and at the CMC gives two equations in
    Gamma_max and K; checked with the toolkit szyszkowski_surface_tension. Then A_min and CPP with the
    toolkit critical_packing_parameter.

    python case5_dmimcl_des_adsorption.py [path/to/toolkit]   -> case5_results.json
"""
from __future__ import annotations

import math

from common import NA, write
from surfactantkit.adsorption import szyszkowski_surface_tension
from surfactantkit.cpp import critical_packing_parameter

R, T = 8.314462618, 298.15
N_C = 10
V_A3 = 27.4 + 26.9 * N_C           # eq (7)
L_A = 1.5 + 1.265 * N_C            # eq (8)
# system: (CMC mM, gamma_CMC, pi_CMC, pC20, Gamma_max printed umol/m2, A_min printed nm2, CPP printed)
ROWS = {
    "water": (38.7, 38.82, 33.17, 1.92, 22.96, 0.0723, 269.809),
    "DES1 5 wt%": (34.2, 37.71, 33.92, 2.14, 18.47, 0.089, 283.535),
    "DES1 10 wt%": (29.3, 41.68, 29.42, 2.27, 12.32, 0.134, 254.838),
    "DES2 5 wt%": (32.4, 40.82, 31.17, 2.18, 14.21, 0.116, 283.38),
    "DES2 10 wt%": (29.0, 38.45, 33.55, 2.39, 13.89, 0.119, 280.113),
    "DES3 5 wt%": (34.0, 39.21, 30.32, 1.91, 25.47, 0.065, 328.115),
    "DES3 10 wt%": (31.0, 41.68, 27.51, 1.98, 18.67, 0.088, 257.967),
}


def solve_k(cmc, c20, pi_cmc):
    f = lambda k: math.log(1 + k * cmc) / math.log(1 + k * c20) - pi_cmc / 20.0
    lo, hi = 1e-6, 1e5
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        lo, hi = (mid, hi) if f(lo) * f(mid) > 0 else (lo, mid)
    return math.sqrt(lo * hi)


def main():
    out = {"paper": "Banjare et al., Heliyon 10 (2024) e35598", "doi": "10.1016/j.heliyon.2024.e35598",
           "chain": {"v_A3": V_A3, "l_A": L_A, "min_cross_section_nm2": round(V_A3 / L_A / 100, 3)}, "systems": {}}
    for name, (cmc, gc, pc, pc20, g_pr, a_pr, cpp_pr) in ROWS.items():
        c20 = 10 ** (-pc20) * 1000                      # mM
        k = solve_k(cmc, c20, pc)                        # 1/mM
        g = pc * 1e-3 / (2 * R * T * math.log(1 + k * cmc))   # mol/m2, n = 2 as in the paper
        g0 = gc + pc
        check_cmc = szyszkowski_surface_tension(cmc * 1e-3, g0, g, k * 1000, "ionic_no_added_salt", T)
        check_c20 = szyszkowski_surface_tension(c20 * 1e-3, g0, g, k * 1000, "ionic_no_added_salt", T)
        a_nm2 = 1e18 / (NA * g)
        out["systems"][name] = {
            "printed": {"Gamma_umol_m2": g_pr, "A_min_nm2": a_pr, "CPP": cpp_pr},
            "A_from_printed_Gamma_nm2": round(1e18 / (NA * g_pr * 1e-6), 4),
            "recovered": {"Gamma_umol_m2": round(g * 1e6, 2), "A_min_nm2": round(a_nm2, 2),
                          "CPP": round(critical_packing_parameter(V_A3, a_nm2 * 100, L_A), 2),
                          "printed_over_recovered_Gamma": round(g_pr / (g * 1e6), 1),
                          "check_gamma_cmc": round(check_cmc, 2), "check_pi_at_C20": round(g0 - check_c20, 2)},
        }
    s = out["systems"].values()
    out["findings"] = {
        "printed_A_below_chain_cross_section": all(v["printed"]["A_min_nm2"] < out["chain"]["min_cross_section_nm2"] for v in s),
        "printed_CPP_range": [min(v["printed"]["CPP"] for v in s), max(v["printed"]["CPP"] for v in s)],
        "recovered_CPP_range": [min(v["recovered"]["CPP"] for v in s), max(v["recovered"]["CPP"] for v in s)],
        "recovered_all_spherical_le_one_third": all(v["recovered"]["CPP"] <= 0.335 for v in s),
    }
    return write("case5", out)


if __name__ == "__main__":
    import json
    print(json.dumps(main(), indent=1))
