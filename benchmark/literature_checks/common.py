"""Shared setup for the literature checks: puts the release toolkit on the path and writes results."""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
TK = Path(sys.argv[1]) if len(sys.argv) > 1 and not sys.argv[1].startswith("-") else HERE.parents[1]
sys.path.insert(0, str(TK / "src"))

R_KJ = 8.314462618e-3     # kJ/(mol K)
NA = 6.02214076e23
WATER_M = 55.5            # mol/L, the value these papers use for X_cmc = CMC/55.5


def write(name: str, out: dict) -> dict:
    out["toolkit"] = TK.name
    (HERE / f"{name}_results.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    return out
