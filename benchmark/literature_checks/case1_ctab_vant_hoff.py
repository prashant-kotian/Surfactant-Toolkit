"""Case 1. CTAB micellization: van't Hoff vs calorimetric enthalpy.

B. Naskar, Colloids Interfaces 9 (2025) 4, doi:10.3390/colloids9010004 (CC BY).
Inputs: Table 1 (ITC CMC, beta, deltaG_m, van't Hoff and ITC deltaH_m), eq (6)-(7) coefficients,
and the stated uncertainties (Table 1 footnote a: CMC 5 %, beta 8 %, deltaG 3 %, deltaH 6 %).
Claims examined: van't Hoff and ITC enthalpies "showed significant differences"; deltaCp -2.7
(van't Hoff) vs -0.45 kJ/(mol K) (ITC).

 1. Reproduce deltaG_m = (1+beta) RT ln X_cmc (toolkit gibbs_free_energy_micellization).
 2. Evaluate the printed fit coefficients; recover the (B, C) that generate the van't Hoff column.
 3. Kirchhoff fit of the printed deltaG(T) (toolkit vant_hoff_multi_point_fit) and Monte Carlo over
    the paper's 3 % deltaG uncertainty; numpy least squares as an independent check.

    python case1_ctab_vant_hoff.py [path/to/toolkit]   -> case1_results.json
"""
from __future__ import annotations

import numpy as np

from common import R_KJ, WATER_M, write
from surfactantkit.thermodynamics import gibbs_free_energy_micellization, vant_hoff_multi_point_fit

T = np.array([290, 293, 298, 303, 308, 313, 323.0])
CMC_ITC_mM = np.array([1.33, 1.25, 1.18, 1.19, 1.35, 1.44, 1.60])
BETA = np.array([0.71, 0.70, 0.68, 0.67, 0.64, 0.63, 0.59])
DG_PRINTED = np.array([43.9, 44.3, 44.8, 45.2, 44.6, 44.8, 44.6])   # printed without a sign
DH_VH = np.array([2.86, -6.10, -20.6, -34.7, -48.3, -61.4, -86.6])
DH_ITC = np.array([-3.50, -6.02, -8.37, -11.2, -12.0, -15.6, -18.7])
ABC = (-872.37, 12.0, 6909.0)
DCP = {"van_t_hoff": -2.7, "itc": -0.45}
TR, N_MC = 298.15, 4000


def tk_fit(dg):
    """X = exp(dG/RT) with counterion factor 1 carries each point's own (1+beta) exactly."""
    r = vant_hoff_multi_point_fit(list(np.exp(dg / (R_KJ * T))), list(T), 1.0, TR)
    return r.delta_H_kJ_per_mol, r.delta_Cp_J_per_mol_K / 1000.0, r


def np_fit(dg):
    r3 = (T - TR) - T * np.log(T / TR)
    h, _s, cp = np.linalg.lstsq(np.column_stack([np.ones_like(T), -T, r3]), dg, rcond=None)[0]
    return h, cp


def main():
    out = {"paper": "Naskar, Colloids Interfaces 9 (2025) 4", "doi": "10.3390/colloids9010004"}
    dg_tk = np.array([gibbs_free_energy_micellization(c * 1e-3 / WATER_M, t, 1 + b)
                      for c, t, b in zip(CMC_ITC_mM, T, BETA)])
    out["reproduce_dG"] = {"toolkit": dg_tk.round(2).tolist(),
                           "max_abs_diff_kJ": round(float(np.abs(np.abs(dg_tk) - DG_PRINTED).max()), 3),
                           "printed_sign": "positive in Table 1"}
    A, B, C = ABC
    out["printed_coefficients"] = {"dG_kJ_from_ABC": ((A + B / T + C / T**2) * T / 1000).round(1).tolist(),
                                   "dH_kJ_from_BC": ((B + 2 * C / T) / 1000).round(3).tolist()}   # same J-based reading as dG
    s, i = np.polyfit(1 / T, DH_VH, 1)
    out["BC_generating_printed_vH_column"] = {"B": round(float(i), 2), "C": round(float(s / 2), 0)}
    dg = -DG_PRINTED
    h, cp, r = tk_fit(dg)
    hn, cpn = np_fit(dg)
    out["kirchhoff_fit_298"] = {"dH_kJ": round(h, 2), "dCp_kJ_per_mol_K": round(cp, 3), "r2": round(r.r_squared, 3),
                                "rms_dG_kJ": round(float(np.sqrt(np.mean((np.array(r.delta_G_fitted_kJ_per_mol) - dg) ** 2))), 3),
                                "numpy_dH_kJ": round(float(hn), 2), "numpy_dCp": round(float(cpn), 3)}
    out["dG_span_kJ"] = round(float(dg.max() - dg.min()), 2)
    out["stated_dG_uncertainty_kJ"] = round(float(0.03 * DG_PRINTED.mean()), 2)
    rng = np.random.default_rng(20261027)
    draws = np.array([tk_fit(dg * (1 + 0.03 * rng.standard_normal(T.size)))[:2] for _ in range(N_MC)])
    lo_h, hi_h = np.percentile(draws[:, 0], [2.5, 97.5])
    lo_c, hi_c = np.percentile(draws[:, 1], [2.5, 97.5])
    out["monte_carlo_3pct"] = {"dH_298_95": [round(lo_h, 1), round(hi_h, 1)], "dCp_95": [round(lo_c, 2), round(hi_c, 2)]}
    out["findings"] = {
        "itc_dH_298_inside_interval": bool(lo_h <= DH_ITC[2] <= hi_h),
        "itc_dCp_inside_interval": bool(lo_c <= DCP["itc"] <= hi_c),
        "printed_vH_dH_298_inside_interval": bool(lo_h <= DH_VH[2] <= hi_h),
    }
    return write("case1", out)


if __name__ == "__main__":
    import json
    print(json.dumps(main(), indent=1))
