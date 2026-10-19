"""Phase 2 MCP wiring: every Phase-2 library function is exposed as a
tool whose name equals the function name (blueprint rule), with the
JSON-callable variants documented in mcp_server. Registration only --
numerics live in the per-worker test files.
"""
import asyncio
import json

import pytest

from surfactantkit.mcp_server import mcp

PHASE2_TOOLS = {
    # O1 / S2
    "mass_action_equilibrium",
    "micelle_relaxation_odes",
    "kmc_micelle_kinetics",
    # E1 / E2
    "obrien_white_numerical",
    "charge_regulation_pb",
    # O3
    "dynamic_surface_tension_forward",
    "dynamic_surface_tension_ode",
    # C3
    "eacn_group_contribution",
    "eacn_linear_mixing",
    # S3
    "brownian_self_assembly",
    "brownian_cmc_scan",
    # I1 / I3 / I4 / I5 / C2
    "mcmc_posterior_fit",
    "nonlinear_fit_with_errors",
    "global_multi_dataset_fit",
    "propagate_uncertainty_chain",
    "cmc_homologous_series",
}


async def _call(tool_name: str, args: dict) -> dict:
    result = await mcp.call_tool(tool_name, args)
    content = result.content[0] if hasattr(result, "content") else result[0]
    return json.loads(content.text)


def test_all_phase2_tools_registered():
    tools = asyncio.run(mcp.list_tools())
    names = {t.name for t in tools}
    missing = PHASE2_TOOLS - names
    assert not missing, f"Phase-2 tools not registered: {sorted(missing)}"


def test_phase2_tool_call_path_matches_library():
    # a cheap end-to-end per wire format: schema keys survive serialization
    out = asyncio.run(_call("mass_action_equilibrium",
                            {"c_total_M": 1e-3, "K_per_M": 5e3,
                             "n_max": 50}))
    assert out["converged"] and out["method"] and out["caveat"]
    # mass balance: sum n * x_n = c_total
    mb = sum((i + 1) * v for i, v in enumerate(out["x_n_M"]))
    assert mb == pytest.approx(1e-3, rel=1e-8)

    out = asyncio.run(_call("eacn_group_contribution",
                            {"smiles": "CCCCCCCCc1ccccc1"}))
    assert out["eacn"] == 8 and out["sources"]
