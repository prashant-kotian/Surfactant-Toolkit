#!/usr/bin/env python3
"""compute_p3.py -- P3 validation on the release toolkit: two independent solvers for dynamic
adsorption of a Langmuir surfactant, compared at three concentrations.

  solver A  adsorption.ward_tordai_adsorption (Volterra product integration, diffusion control)
  solver B  dynamic_tension.dynamic_surface_tension_forward (implicit finite-difference diffusion with
            Langmuir adsorption kinetics, k_off = 200 1/s)

Protocol of the 2026-10-05 reference run: D = 1e-9 m2/s, Gamma_inf = 3e-6 mol/m2, K_d = 1e-3 M,
gamma0 = 72 mN/m, t_k = 5 (k/16)^2 s, k = 1..16. Surface tension from the surface excess with the
Szyszkowski (Langmuir) equation of state for both solvers. Writes values_p3.json.

    python compute_p3.py [path/to/Surfactant-Toolkit]
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
TK = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE.parents[1]
sys.path.insert(0, str(TK / "src"))
from surfactantkit.adsorption import ward_tordai_adsorption  # noqa: E402
from surfactantkit.dynamic_tension import dynamic_surface_tension_forward  # noqa: E402

R, T, D, GI, KD, G0, KOFF = 8.314462618, 298.15, 1e-9, 3e-6, 1e-3, 72.0, 200.0


def eos(gamma_excess: float) -> float:
    return G0 + 1000.0 * R * T * GI * math.log(1.0 - gamma_excess / GI)


def main() -> int:
    out = {"protocol": "Langmuir, D=1e-9 m2/s, Gamma_inf=3e-6 mol/m2, Kd=1e-3 M, gamma0=72 mN/m, "
                       "k_off=200 1/s, t_k=5(k/16)^2 s", "per_C0": {}}
    for c0 in (1e-4, 1e-3, 1e-2):
        wt = ward_tordai_adsorption(c0, D, t_end_s=5.0, n_times=16, mode="langmuir",
                                    gamma_inf_mol_per_m2=GI, K_d_mol_per_L=KD)
        pde = dynamic_surface_tension_forward(wt.t_s, c0, D, 1.0 / KD, GI, gamma0_mN_per_m=G0,
                                              temperature_K=T, k_off_per_s=KOFF)
        g_wt = [eos(g) for g in wt.gamma_mol_per_m2]
        rel = [abs(a - b) / b for a, b in zip(pde.surface_excess_mol_per_m2, wt.gamma_mol_per_m2)]
        dg = [abs(a - b) for a, b in zip(pde.gamma_mN_per_m, g_wt)]
        out["per_C0"][f"{c0:g}"] = {
            "worst_Gamma_rel": max(rel), "worst_gamma_mN_m": max(dg),
            "gamma_final_WT": g_wt[-1], "gamma_final_PDE": pde.gamma_mN_per_m[-1],
            "gamma_eq": pde.gamma_equilibrium_mN_per_m, "converged": [wt.converged, pde.converged]}
    (HERE / "values_p3.json").write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")
    for c0, r in out["per_C0"].items():
        print(f"C0 {c0:>6}: worst Gamma {100 * r['worst_Gamma_rel']:.2f} %, worst gamma "
              f"{r['worst_gamma_mN_m']:.3f} mN/m, final {r['gamma_final_WT']:.3f}/{r['gamma_final_PDE']:.3f} "
              f"(eq {r['gamma_eq']:.3f})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
