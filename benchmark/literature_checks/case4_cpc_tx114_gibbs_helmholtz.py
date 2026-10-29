"""Case 4. CPC/TX-114 mixed micelles with Reactive Red 2: free energy and enthalpy that cannot both hold.

T. Yaqoob, S. Shaukat, R. Alonaizan, R. Ullah, I. Khan, M. F. Nazar, H. M. Abd Ur Rahman,
Molecules 28 (2023) 3952, doi:10.3390/molecules28093952 (CC BY).
Inputs: Table 4 (T, CMC, deltaG, deltaH, deltaS, alpha) and eqs (5)-(6):
deltaG = (2 - alpha) RT ln X_cmc, deltaH = -2.303 (2 - alpha) RT^2 d(log X_cmc)/dT.
Claims examined: micellization in the mixed system is exothermic (deltaH -5.0 to -7.1 kJ/mol) and
more spontaneous than in CPC alone at every temperature.

 1. Reproduce: deltaG at 303-323 K and the printed deltaH from eq (6) (slope of log X vs T).
 2. deltaG at 293.15 K recomputed (toolkit gibbs_free_energy_micellization) vs printed.
 3. Gibbs-Helmholtz on the recomputed deltaG(T) (toolkit vant_hoff_multi_point_fit, reference
    temperature at each T) and a numpy finite-difference check; and the same with alpha held at
    its mean, which isolates the alpha(T) contribution that eq (6) leaves out.

    python case4_cpc_tx114_gibbs_helmholtz.py [path/to/toolkit]   -> case4_results.json
"""
from __future__ import annotations

import numpy as np

from common import R_KJ, WATER_M, write
from surfactantkit.thermodynamics import gibbs_free_energy_micellization, vant_hoff_multi_point_fit

T = np.array([293.15, 303.15, 313.15, 323.15])
CMC_mM = np.array([0.87, 0.89, 0.92, 0.98])
ALPHA = np.array([0.484, 0.505, 0.332, 0.230])
DG_PR = np.array([-32.4, -41.5, -47.6, -51.9])
DH_PR = np.array([-5.0, -5.2, -6.2, -7.1])
DS_PR = np.array([0.009, 0.120, 0.132, 0.139])       # kJ/(K mol)
DG_CPC_293 = -36.1                                   # Table 3, CPC alone


def gh_dh(dg):
    """deltaH(T_i) from a Kirchhoff fit of dG(T) referenced at each T_i (toolkit)."""
    x = np.exp(dg / (R_KJ * T))
    return np.array([vant_hoff_multi_point_fit(list(x), list(T), 1.0, float(t)).delta_H_kJ_per_mol for t in T])


def main():
    out = {"paper": "Yaqoob et al., Molecules 28 (2023) 3952", "doi": "10.3390/molecules28093952"}
    x = CMC_mM * 1e-3 / WATER_M
    dg = np.array([gibbs_free_energy_micellization(xi, t, 2 - a) for xi, t, a in zip(x, T, ALPHA)])
    out["dG_recomputed"] = dg.round(1).tolist()
    out["dG_printed"] = DG_PR.tolist()
    slope = np.polyfit(T, np.log10(x), 1)[0]
    dh_eq6 = -2.303 * (2 - ALPHA) * R_KJ * T**2 * slope
    out["reproduce_dH_eq6"] = {"dH": dh_eq6.round(1).tolist(), "printed": DH_PR.tolist()}
    out["dS_from_printed_dH_dG"] = ((DH_PR - DG_PR) / T).round(3).tolist()
    out["dS_printed"] = DS_PR.tolist()
    dh_gh = gh_dh(dg)
    dh_fd = -T**2 * np.gradient(dg / T, T)            # numpy: dH = -T^2 d(dG/T)/dT
    dg_const_alpha = np.array([gibbs_free_energy_micellization(xi, t, 2 - ALPHA.mean()) for xi, t in zip(x, T)])
    out["dG_alpha_held_at_mean"] = dg_const_alpha.round(1).tolist()
    out["gibbs_helmholtz_dH"] = {"toolkit": dh_gh.round(0).tolist(), "numpy_finite_difference": dh_fd.round(0).tolist(),
                                 "alpha_held_at_mean_toolkit": gh_dh(dg_const_alpha).round(1).tolist()}
    out["findings"] = {
        "dG_293_recomputed": round(float(dg[0]), 1), "dG_293_printed": float(DG_PR[0]),
        "dS_293_from_printed_dH_dG": round(float((DH_PR[0] - DG_PR[0]) / T[0]), 3), "dS_293_printed": float(DS_PR[0]),
        "mixed_vs_cpc_293_printed_more_spontaneous": bool(DG_PR[0] < DG_CPC_293),
        "mixed_vs_cpc_293_recomputed_more_spontaneous": bool(dg[0] < DG_CPC_293),
        "eq6_reproduces_printed_dH_within_kJ": round(float(np.abs(dh_eq6 - DH_PR).max()), 1),
        "gh_dH_303_323_kJ": dh_gh[1:].round(0).tolist(),
        "gh_dH_positive_303_323": bool((dh_gh[1:] > 0).all() and (dh_fd[1:] > 0).all()),
        "printed_dH_negative_all_T": bool((DH_PR < 0).all()),
        "dG_change_293_323_printed_vs_alpha_held": [round(float(DG_PR[-1] - DG_PR[0]), 1),
                                                    round(float(dg_const_alpha[-1] - dg_const_alpha[0]), 1)],
    }
    return write("case4", out)


if __name__ == "__main__":
    import json
    print(json.dumps(main(), indent=1))
