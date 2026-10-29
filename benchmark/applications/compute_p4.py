#!/usr/bin/env python3
"""compute_p4.py -- P4 reference values on the release toolkit: half-creaming time of the
paraffin-in-water emulsion of Murasiewicz and Illienko (2024), d32 = 235 um.

Inputs as in the 2026-10-05 reference run: paraffin 870 kg/m3 (assumed, typical), water 998 kg/m3,
viscosity 1.0 mPa s, 10 cm column, no coalescence (the paper reports an unchanged drop-size
distribution over a month). Monodisperse, and lognormal distributions of width sigma at the same
d32 (nine equal-volume bins, rescaled so that d32 = 235 um exactly). Writes values_p4.json.

    python compute_p4.py [path/to/Surfactant-Toolkit]
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from statistics import NormalDist

HERE = Path(__file__).resolve().parent
TK = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE.parents[1]
sys.path.insert(0, str(TK / "src"))
from surfactantkit.emulsion import simulate_emulsion_stability  # noqa: E402

D32, RHO_O, RHO_W, ETA, H = 235.0, 870.0, 998.0, 1.0e-3, 0.10


def bins(sigma: float, n: int = 9):
    z = [NormalDist().inv_cdf((k - 0.5) / n) for k in range(1, n + 1)]
    d = [math.exp(sigma * zk) for zk in z]
    w = [1.0 / n] * n
    scale = D32 / (sum(w) / sum(wi / di for wi, di in zip(w, d)))
    return [di * scale for di in d], w


def t50(d, w) -> float:
    tmax = 4.0 * H / min(abs(RHO_O - RHO_W) * 9.80665 * (x * 1e-6) ** 2 / (18 * ETA) for x in d)
    r = simulate_emulsion_stability(d, w, RHO_O, RHO_W, ETA, H, n_times=20001, t_max_s=tmax)
    return r.t50_s


def main() -> int:
    out = {"inputs": {"d32_um": D32, "rho_oil": RHO_O, "rho_water": RHO_W, "viscosity_Pa_s": ETA,
                      "column_m": H, "coalescence_per_s": 0.0}}
    out["mono_t50_s"] = t50([D32], [1.0])
    v = abs(RHO_O - RHO_W) * 9.80665 * (D32 * 1e-6) ** 2 / (18 * ETA)
    assert abs(out["mono_t50_s"] - H / (2 * v)) < 1e-2, out["mono_t50_s"]
    out["poly_t50_s"] = {f"{s:g}": t50(*bins(s)) for s in (0.3, 0.5, 0.8)}
    (HERE / "values_p4.json").write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")
    print(f"monodisperse t50 {out['mono_t50_s']:.2f} s (H/2v {H / (2 * v):.2f} s)")
    for s, t in out["poly_t50_s"].items():
        print(f"lognormal sigma {s}: t50 {t:.2f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
