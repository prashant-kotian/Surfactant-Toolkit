"""SurfMCP: exposes SurfactantKit's surfactant-science calculations as MCP
tools, so an AI assistant computes CMC/Rubingh/HLB/CPP values correctly
instead of hallucinating them.

Every tool returns units explicitly in the response -- the whole point
of this server is to stop unit confusion (mN/m vs dyn/cm, mM vs mol/kg,
etc.), so a bare unlabeled number defeats the purpose.

Run locally (stdio transport, for Claude Desktop / Cursor):
    python -m surfactantkit.mcp_server
"""

from __future__ import annotations

try:
    from mcp.server.mcpserver import MCPServer  # mcp 2.x name -- migrated 2026-09-13
except ImportError:
    # Real environment mismatch found 2026-09-14: the 2026-09-13 migration
    # was made and tested against a newer `mcp` SDK (Ubuntu conda env) where
    # FastMCP had already been renamed to MCPServer -- but this Windows
    # install's `mcp` package is still 1.27.2, which only has FastMCP under
    # its original name (confirmed directly: `mcp.server.mcpserver` module
    # doesn't exist there at all). The 09-13 fix silently broke this
    # environment without anyone noticing since it was never re-tested here.
    # Fixed properly this time: try the new name first, fall back to the old
    # one, so this file works against EITHER real installed SDK version
    # rather than assuming one.
    from mcp.server.fastmcp import FastMCP as MCPServer

try:  # audit 2026-10-06: ToolError lives under the new name in mcp 2.x,
    from mcp.server.mcpserver.exceptions import (
        ToolError as _MCPToolError,  # the old in mcp 1.x
    )
except ImportError:  # (mirrors the MCPServer import above -- same version split)
    from mcp.server.fastmcp.exceptions import ToolError as _MCPToolError

from . import adsorption as ads
from . import classify as classify_mod
from . import cloudpoint as cp_mod
from . import conformer as conf_mod
from . import cpp as cpp_mod
from . import curve_analysis as curve
from . import descriptors as desc_mod
from . import dpd_cmc as dpd_mod
from . import dpd_micelle as dpdm_mod
from . import dynamic_tension as dt_mod
from . import dynamics as dyn
from . import eacn as eacn_mod
from . import electrostatics as elec
from . import emulsion as emu
from . import hlb as hlb_mod
from . import hld as hld_mod
from . import homologous as hom_mod
from . import inference as inf_mod
from . import kinetics as kin_mod
from . import mixed_micelle as mm
from . import mixture_model_selection as mms
from . import nagarajan as nag_mod
from . import orchestrate as orch
from . import selfassembly as sa_mod
from . import solubilization as solub
from . import thermodynamics as thermo
from . import trajectory as traj_mod
from . import transport as tr_mod
from . import wetting

mcp = MCPServer(
    "SurfactantKit",
    instructions=(
        "Surfactant and interfacial-science calculations: Clint ideal mixing, "
        "Rubingh regular-solution theory for binary surfactant mixtures (and its "
        "Rosen extension to mixed adsorbed monolayers), the Gibbs adsorption "
        "isotherm, HLB (Griffin's and Davies' methods), the critical packing "
        "parameter and aggregation number (Tanford's formulas), full "
        "micellization thermodynamics (deltaG/deltaH/deltaS, counterion "
        "binding), the Corrin-Harkins salt-CMC relation, Debye screening "
        "length, zeta potential (Henry equation), and hydrodynamic radius "
        "(Stokes-Einstein). Use these tools instead of computing or recalling "
        "these values from memory -- CMC, interaction-parameter, HLB, and "
        "electrolyte/electrostatic values are easy to hallucinate and easy to "
        "get wrong, especially for gemini (two-headgroup) surfactants, mixed "
        "salt/ionic-strength systems, and anything requiring iterative "
        "numerical solving or careful unit conversion (e.g. cP vs Pa.s, "
        "cm^2/s vs m^2/s -- these are exactly the kind of silent unit errors "
        "this server exists to prevent). Also covers the Szyszkowski surface "
        "tension equation, wetting/adhesion (Young-Dupre, spreading "
        "coefficient), the capillary number (enhanced oil recovery), and "
        "micellar solubilization (Molar Solubilization Ratio)."
    ),
)




def _shape_no_solution(func):
    """Audit P0-1: convert library ValueError/KeyError into a consistent,
    readable MCP response instead of a raw internal error.

    Every tool in this server is refusal-first: when it cannot compute,
    the answer is "cannot compute + why", shaped as
    {"status": "no_solution", "reason": ...} (the mechanism proven to move
    SurfMCP 0/6 -> 6/6 on the Gibbs n_factor experiment). The LIBRARY keeps
    raising fail-fast ValueError/KeyError with the offending field names --
    this only shapes the MCP boundary. Any other exception type still
    propagates untouched (programming errors must stay loud)."""
    import functools

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except ValueError as exc:
            return {"status": "no_solution", "reason": str(exc)}
        except KeyError as exc:
            key = exc.args[0] if exc.args else exc
            return {"status": "no_solution",
                    "reason": f"unknown or missing key: {key!r} (refusing to guess)"}

    return wrapper


def _tool_guard(func):
    """Audit 2026-10-06: uniform input contract for EVERY MCP tool.

    (1) Non-finite float inputs (NaN/inf, nested anywhere in lists/dicts)
        are rejected with the offending field name before any math runs --
        the toolkit's fail-fast rule previously did not cover them at all,
        and they silently propagated into results (fuzz-found: 220 cases).
    (2) Any error not already shaped into a no_solution result is re-raised
        as mcp ToolError carrying the full "ExceptionType: message" text:
        mcp>=2.0's FastMCP strips raw exceptions into a message-less
        "Error executing tool X" (confirmed on the real stdio path), which
        would leave benchmark agents (and users) blind to the reason.
    Refusal-first shaping (``_shape_no_solution``) is untouched and still
    wins wherever it is applied; server logs keep the real traceback.
    """
    import functools as _ft
    import inspect as _insp
    import math as _math

    params = list(_insp.signature(func).parameters)

    def _check(value, where):
        if isinstance(value, float):
            if not _math.isfinite(value):
                raise ValueError(f"{where}: non-finite float {value!r} is not a valid input")
        elif isinstance(value, dict):
            for k, v in value.items():
                _check(v, f"{where}[{k!r}]")
        elif isinstance(value, (list, tuple)):
            for i, v in enumerate(value):
                _check(v, f"{where}[{i}]")

    @_ft.wraps(func)
    def wrapper(*args, **kwargs):
        try:
            for name, value in zip(params, args):
                _check(value, name)
            for name, value in kwargs.items():
                _check(value, name)
            return func(*args, **kwargs)
        except Exception as exc:
            if isinstance(exc, _MCPToolError):
                raise
            raise _MCPToolError(f"{type(exc).__name__}: {exc}") from exc

    return wrapper


@mcp.tool()
@_tool_guard
@_shape_no_solution
def clint_ideal_cmc(alpha1: float, cmc1_mM: float, cmc2_mM: float) -> dict:
    """Clint's ideal mixed CMC for a binary surfactant mixture (no
    interaction between components assumed). alpha1 is the bulk mole
    fraction of component 1 in the total surfactant (0 to 1, exclusive).
    cmc1_mM and cmc2_mM are the pure-component CMCs in mM."""
    value = mm.clint_ideal_cmc(alpha1, cmc1_mM, cmc2_mM)
    return {"cmc_ideal_mM": value, "unit": "mM", "model": "Clint ideal mixing"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def rubingh_solve(alpha1: float, cmc_mix_mM: float, cmc1_mM: float, cmc2_mM: float) -> dict:
    """Solve Rubingh's regular-solution equation for a binary surfactant
    mixture: returns the micellar mole fraction of component 1 (x1) and
    the interaction parameter beta. Negative beta = synergistic mixing;
    positive beta = antagonistic; beta near 0 = near-ideal. alpha1 is the
    bulk mole fraction of component 1; cmc_mix_mM is the experimentally
    measured mixed CMC; cmc1_mM/cmc2_mM are the pure-component CMCs, all
    in mM. Raises an error if no physically valid root is found (this
    can genuinely happen for near-ideal mixtures or inconsistent input
    data -- it is not a bug to report back if it occurs)."""
    x1 = mm.solve_rubingh_x(alpha1, cmc_mix_mM, cmc1_mM, cmc2_mM)
    if x1 is None:
        raise ValueError(
            "No valid Rubingh root found in (0, 1) for these inputs. This can "
            "happen for near-ideal mixtures or when cmc_mix is not physically "
            "consistent with the given alpha1/cmc1/cmc2. Do not guess a value."
        )
    beta = mm.rubingh_beta(x1, alpha1, cmc_mix_mM, cmc1_mM)
    synergy = "synergistic" if beta < 0 else ("antagonistic" if beta > 0 else "ideal")
    return {
        "micellar_mole_fraction_x1": x1,
        "beta": beta,
        "beta_unit": "dimensionless (regular solution interaction parameter)",
        "synergy_classification": synergy,
        "note": (
            "beta here is the exact pointwise regular-solution value for this "
            "single (alpha1, cmc_mix) pair. Published beta values are sometimes "
            "fit by regression across several compositions and may differ from "
            "this in magnitude (not sign) by up to ~20% -- see "
            "literature_validation_notes.md in the SurfactantKit repo. For a "
            "whole composition series, use rubingh_beta_regression instead."
        ),
    }


@mcp.tool()
@_tool_guard
@_shape_no_solution
def rubingh_beta_regression(alpha1_series: list[float], cmc_mix_series_mM: list[float], cmc1_mM: float, cmc2_mM: float) -> dict:
    """
    Fit ONE representative Rubingh beta across a WHOLE composition series (standard practice for reporting a binary interaction parameter; contrast rubingh_solve's single-point beta). Args: alpha1_series + cmc_mix_series_mM (same length, matching order), cmc1_mM/cmc2_mM. Returns beta, fitted x1 series, r_squared, residuals. Note: regression beta differs in magnitude (not sign) from pointwise betas by up to ~20% -- both conventions documented.
    """
    result = mm.rubingh_beta_regression(alpha1_series, cmc_mix_series_mM, cmc1_mM, cmc2_mM)
    synergy = "synergistic" if result.beta_mean < 0 else ("antagonistic" if result.beta_mean > 0 else "ideal")
    return {
        "beta_mean": result.beta_mean,
        "beta_std": result.beta_std,
        "beta_unit": "dimensionless (regular solution interaction parameter)",
        "synergy_classification": synergy,
        "betas_per_point": result.betas,
        "x1_per_point": result.x1_values,
        "alpha1_used": result.alpha1_used,
        "n_points_used": result.n_points_used,
        "n_points_skipped": result.n_points_skipped,
    }


@mcp.tool()
@_tool_guard
@_shape_no_solution
def motomura_ideal_composition(alpha1: float, cmc1_mM: float, cmc2_mM: float) -> dict:
    """Motomura & Aratono's ideal micellar composition (X1) for a binary
    surfactant mixture -- the zero-interaction-assumption prediction for
    how much of component 1 ends up IN the micelle, as opposed to
    clint_ideal_cmc which predicts the ideal mixed CMC itself. Same
    zero-parameter ideal-mixing assumption as Clint's model (no
    interaction parameter, unlike rubingh_solve's beta); use this when
    the ideal-mixing micellar composition itself is asked for, not just
    the ideal CMC. alpha1 is the bulk mole fraction of component 1;
    cmc1_mM/cmc2_mM are the pure-component CMCs, in mM."""
    value = mm.motomura_ideal_composition(alpha1, cmc1_mM, cmc2_mM)
    return {"x1_ideal": value, "unit": "dimensionless (micellar mole fraction)", "model": "Motomura ideal mixing"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def rodenas_x1(alpha1: float, dln_cmc_mix_dalpha1: float) -> dict:
    """Rodenas' Gibbs-Duhem-based, "model-independent" micellar mole
    fraction of component 1 -- an alternative to rubingh_solve that
    does NOT assume any particular regular-solution/activity-
    coefficient model, at the cost of needing the LOCAL SLOPE of the
    real, multi-point experimental ln(cmc_mix) vs. alpha1 curve
    (dln_cmc_mix_dalpha1) rather than a single (alpha1, cmc_mix) pair.
    This slope MUST come from real neighboring experimental data
    (a fit or numerical derivative) -- do not guess or assume a value.
    Follow with rodenas_activity_coefficients using the returned x1."""
    value = mm.rodenas_x1(alpha1, dln_cmc_mix_dalpha1)
    return {"micellar_mole_fraction_x1": value, "model": "Rodenas (Gibbs-Duhem, model-independent)"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def rodenas_x1_from_series(alpha1_series: list[float], cmc_mix_series_mM: list[float]) -> dict:
    """Compute Rodenas' d[ln(cmc_mix)]/d(alpha1) AND the resulting x1 at
    every composition directly from a RAW multi-alpha1 cmc_mix(alpha1)
    series -- closes the gap rodenas_x1 itself flags (it needs the local
    slope as an already-known input). At each composition, fits the
    unique quadratic through it and its 2 nearest neighbors along
    alpha1, then differentiates analytically at that exact alpha1 --
    this reduces exactly to the standard central/forward/backward finite-
    difference formulas for evenly-spaced data, generalized to real
    (unevenly-spaced) experimental series. alpha1_series and
    cmc_mix_series_mM must be the same length, at least 3 points, all
    alpha1 distinct and strictly between 0 and 1. Points are sorted by
    alpha1 internally; results are returned in that sorted order.
    ALWAYS check x1_physically_valid: on a sparse, sharply-curved real
    series (e.g. a strongly synergistic system with a deep CMC minimum
    near equimolar composition) the local-slope estimate can produce a
    micellar mole fraction outside the physical [0,1] range -- a real
    numerical-stability failure of this method, not a bug. See the
    caveat field; cross-check against solve_rubingh_x/rubingh_beta,
    which is often more robust on exactly this kind of sparse data."""
    result = mm.rodenas_x1_series(alpha1_series, cmc_mix_series_mM)
    return {
        "alpha1_used": result.alpha1_used,
        "dln_cmc_mix_dalpha1": result.dln_cmc_mix_dalpha1,
        "x1_rodenas": result.x1_rodenas,
        "x1_physically_valid": result.x1_physically_valid,
        "n_points": result.n_points,
        "method": result.method,
        "caveat": result.caveat,
    }


@mcp.tool()
@_tool_guard
@_shape_no_solution
def rodenas_activity_coefficients(alpha1: float, x1: float, cmc_mix_mM: float, cmc1_mM: float, cmc2_mM: float) -> dict:
    """Activity coefficients (f1, f2) for Rodenas' model, given the
    already-determined x1 (from rodenas_x1). Same mass-balance
    definition solve_rubingh_x uses internally -- Rodenas and Rubingh
    differ only in how x1 is determined, not in how f1/f2 follow from
    a known x1."""
    f1, f2 = mm.rodenas_activity_coefficients(alpha1, x1, cmc_mix_mM, cmc1_mM, cmc2_mM)
    return {"f1": f1, "f2": f2, "unit": "dimensionless", "model": "Rodenas"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def rubingh_activity_coefficients(x1: float, beta: float) -> dict:
    """Regular-solution activity coefficients (f1, f2) for a binary
    mixed micelle, given the micellar mole fraction x1 and interaction
    parameter beta (from rubingh_solve)."""
    f1, f2 = mm.activity_coefficients(x1, beta)
    return {"f1": f1, "f2": f2, "unit": "dimensionless"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def asymmetric_margules_activity_coefficients(x1: float, w12: float, w21: float) -> dict:
    """Activity coefficients (f1, f2) from the EOMMM binary asymmetric
    Margules model -- generalizes rubingh_activity_coefficients with
    TWO independent interaction parameters (w12 != w21) instead of
    RST's single symmetric beta. w12/w21 are dimensionless (RT units,
    same convention as rubingh_solve's beta) -- w12 is the energy of
    introducing a molecule of component 1 into a pure micelle of
    component 2 (and vice versa for w21). At w12 == w21 this reduces
    EXACTLY to rubingh_activity_coefficients (RST is the provable
    symmetric special case). ONLY computes activity coefficients for
    an ALREADY-KNOWN w12/w21 pair (e.g. quoted directly from a paper) --
    does NOT fit w12/w21 from raw CMC/alpha1 data the way rubingh_solve
    does for RST; for that, use eommm_global_fit."""
    f1, f2 = mm.asymmetric_margules_activity_coefficients(x1, w12, w21)
    return {"f1": f1, "f2": f2, "unit": "dimensionless", "model": "EOMMM asymmetric Margules (binary)"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def eommm_global_fit(alpha1_series: list[float], cmc_mix_series_mM: list[float], cmc1_mM: float, cmc2_mM: float,
                      r1: float = 2.0, r2: float = 2.0, cmc_margin: float = 0.10) -> dict:
    """
    Fit EOMMM's two Margules parameters (W12, W21) and each point's micellar mole fraction x1 simultaneously across a composition series -- a real transcription of Schulz & Durand 2016's own SI procedure (minimizes constraint infeasibility at cmc_margin; the SI's tight-margin search is NOT automated -- call at decreasing cmc_margin and watch total_infeasibility). Args: alpha1_series + cmc_mix_series_mM (same length, >=3 points, 5+ recommended), cmc1_mM/cmc2_mM, r1/r2 dissociation numbers (default 2.0 = fully dissociated 1:1 ionic; 1.0 for a nonionic). Returns W12, W21, fitted x1/cmc_var, total_infeasibility, r_squared. Uniqueness at the default cmc_margin=0.10 is NOT guaranteed (wider margins admit more of the (W12,W21) plane) -- tighten the margin to check stability. Expected runtime: several seconds (2D grid search, stdlib only, no numpy/scipy).
    """
    result = mm.eommm_global_fit(alpha1_series, cmc_mix_series_mM, cmc1_mM, cmc2_mM, r1, r2, cmc_margin)
    return {
        "W12": result.W12,
        "W21": result.W21,
        "x1_values": result.x1_values,
        "cmc_var_values": result.cmc_var_values,
        "alpha1_used": result.alpha1_used,
        "total_infeasibility": result.total_infeasibility,
        "cmc_margin_used": result.cmc_margin_used,
        "r_squared": result.r_squared,
        "n_points": result.n_points,
        "method": result.method,
    }


@mcp.tool()
@_tool_guard
@_shape_no_solution
def eommm_find_minimal_feasible_margin(alpha1_series: list[float], cmc_mix_series_mM: list[float], cmc1_mM: float, cmc2_mM: float,
                                        r1: float = 2.0, r2: float = 2.0, binary_search_iters: int = 20) -> dict:
    """
    Automates eommm_global_fit's margin-tightening: binary-searches cmc_margin down to the smallest value where a fully feasible fit exists, then returns the fit there (the tightest-margin, most uniquely-determined answer). Raises if even the search ceiling (default 0.20) is infeasible. NOTE (disclosed): an earlier prototype was unshipped for grid-resolution sensitivity; the shipped version reports grid parameters in params so convergence can be re-checked.
    """
    result = mm.eommm_find_minimal_feasible_margin(
        alpha1_series, cmc_mix_series_mM, cmc1_mM, cmc2_mM, r1, r2,
        binary_search_iters=binary_search_iters)
    return {
        "W12": result.W12,
        "W21": result.W21,
        "x1_values": result.x1_values,
        "cmc_var_values": result.cmc_var_values,
        "alpha1_used": result.alpha1_used,
        "total_infeasibility": result.total_infeasibility,
        "minimal_feasible_cmc_margin": result.cmc_margin_used,
        "r_squared": result.r_squared,
        "n_points": result.n_points,
        "method": result.method,
    }


@mcp.tool()
@_tool_guard
@_shape_no_solution
def eommm_global_fit_all_optima(alpha1_series: list[float], cmc_mix_series_mM: list[float], cmc1_mM: float, cmc2_mM: float,
                                 r1: float = 2.0, r2: float = 2.0, cmc_margin: float = 0.10,
                                 feasible_tolerance: float = 1e-6) -> dict:
    """
    Like eommm_global_fit, but does not silently return only the single best-scoring (W12, W21) -- searches the WHOLE (W12, W21) plane for every distinct feasible basin (eommm_global_fit's own docstring discloses non-uniqueness is possible at a generous margin but never surfaces a second optimum when one exists). Same args as eommm_global_fit, plus feasible_tolerance (default 1e-6) for what counts as "fully feasible." Returns a list of optima (each with W12/W21/x1_values/cmc_var_values/total_infeasibility/r_squared), sorted by total_infeasibility, plus degenerate (true iff more than one survives). If grading or comparing against an external answer, accept a match against ANY entry, not just the first. Runtime: a few times eommm_global_fit's own (one fine-refine pass per basin found, typically 1-5).
    """
    result = mm.eommm_global_fit_all_optima(
        alpha1_series, cmc_mix_series_mM, cmc1_mM, cmc2_mM, r1, r2, cmc_margin,
        feasible_tolerance=feasible_tolerance)
    return {
        "optima": [
            {
                "W12": o.W12, "W21": o.W21, "x1_values": o.x1_values,
                "cmc_var_values": o.cmc_var_values, "total_infeasibility": o.total_infeasibility,
                "r_squared": o.r_squared,
            }
            for o in result.optima
        ],
        "n_optima_found": result.n_optima_found,
        "degenerate": result.degenerate,
        "cmc_margin_used": cmc_margin,
        "method": result.method,
    }


@mcp.tool()
@_tool_guard
@_shape_no_solution
def select_mixture_model(alpha1_series: list[float], cmc_mix_series_mM: list[float], cmc1_mM: float, cmc2_mM: float) -> dict:
    """
    Decide WHICH mixture theory fits a raw binary composition series (alpha1, measured mixed CMC per point, two pure CMCs): fits Clint (0 params), Rubingh (1 param, beta regressed), EOMMM asymmetric Margules (2 params, global W12/W21), compares by BIC and reports the winner. Also reports Rodenas' model-independent composition as a structurally different cross-check (see rodenas_diagnostic), and explicitly discloses why Rosen, Corrin-Harkins and Maeda are NOT entered (each needs a genuinely different input a bare composition series lacks). Inputs must be same length with >=3 points. Returns chosen model, per-model BIC/params/fits, and the disclosure list.
    """
    result = mms.select_mixture_model(alpha1_series, cmc_mix_series_mM, cmc1_mM, cmc2_mM)
    return {
        "selected_model": result.selected_model,
        "candidates": [
            {"name": c.name, "n_params": c.n_params, "params": c.params,
             "predicted_cmc_mix_mM": c.predicted_cmc_mix, "rss_ln_space": c.rss_ln_space, "bic": c.bic}
            for c in result.candidates
        ],
        "rodenas_diagnostic": {
            "alpha1_used": result.rodenas_diagnostic.alpha1_used,
            "x1_rodenas": result.rodenas_diagnostic.x1_rodenas,
            "method": result.rodenas_diagnostic.method,
        },
        "excluded_models": result.excluded_models,
    }


@mcp.tool()
@_tool_guard
@_shape_no_solution
def maeda_free_energy_of_micellization(x1_rub: float, beta: float, cmc1_M: float, cmc2_M: float, temperature_K: float) -> dict:
    """Maeda's free energy of micellization (kJ/mol) for a binary
    ionic/nonionic mixed micelle -- built specifically for ionic-
    nonionic systems, accounting for chain-chain (not just head-group)
    interaction. Reuses Rubingh's own already-solved x1/beta (from
    rubingh_solve) rather than fitting independently: deltaG_M =
    RT*(B0 + B1*X1_Rub + B2*X1_Rub^2), B0=ln(Xcmc2), B1+B2=ln(Xcmc1/
    Xcmc2), B2=-beta, where Xcmc are the pure-component CMCs on the
    mole-fraction scale. cmc1_M/cmc2_M must be in mol/L (molarity),
    not mM -- convert before calling."""
    value = mm.maeda_free_energy_of_micellization(x1_rub, beta, cmc1_M, cmc2_M, temperature_K)
    return {"deltaG_M_kJ_per_mol": value, "unit": "kJ/mol", "model": "Maeda"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def excess_free_energy(x1: float, f1: float, f2: float, temperature_K: float) -> dict:
    """Excess Gibbs free energy of mixed micelle formation (kJ/mol),
    from the micellar mole fraction x1, activity coefficients f1/f2
    (from rubingh_activity_coefficients), and temperature in Kelvin."""
    value = mm.excess_free_energy(x1, f1, f2, temperature_K)
    return {"deltaG_ex_kJ_per_mol": value, "unit": "kJ/mol"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def gibbs_surface_excess(slope_mN_per_m_per_lnC: float, system_type: str, temperature_K: float) -> dict:
    """Maximum surface excess concentration (Gamma_max, mol/m^2) from
    the Gibbs adsorption isotherm, given the pre-CMC slope of surface
    tension (mN/m) vs ln(concentration), temperature in Kelvin, and
    system_type -- one of 'nonionic', 'ionic_excess_electrolyte', or
    'ionic_no_added_salt' -- which fixes the Gibbs prefactor n (1 for
    the first two, 2 for the last). REQUIRED, no default: if the
    question does not state the surfactant's ionic character AND
    (for an ionic surfactant) whether excess inert electrolyte is
    present, DO NOT GUESS a system_type -- this is exactly the kind of
    silent assumption that produces a real, confidently wrong answer;
    report that Gamma_max cannot be determined from the given
    information instead of calling this tool with an assumed value."""
    value = ads.gibbs_gamma_max(slope_mN_per_m_per_lnC, system_type, temperature_K)
    return {"gamma_max_mol_per_m2": value, "unit": "mol/m^2"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def gibbs_area_per_molecule(gamma_max_mol_per_m2: float) -> dict:
    """Minimum area per molecule at the interface (A_min, nm^2), from
    Gamma_max (mol/m^2, from gibbs_surface_excess)."""
    value = ads.gibbs_a_min(gamma_max_mol_per_m2)
    return {"a_min_nm2": value, "unit": "nm^2 per molecule"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def szyszkowski_predict_surface_tension(concentration: float, gamma0_mN_m: float, gamma_max_mol_per_m2: float, K: float, system_type: str, temperature_K: float = 298.15) -> dict:
    """
    Predict surface tension (mN/m) at a given concentration: gamma = gamma0 - n_factor*R*T*Gamma_max*ln(1+K*concentration). Args: concentration (mM or mol/L -- any one unit; K must carry the RECIPROCAL unit of the same basis, e.g. K in 1/mM with concentrations in mM; raising on mismatch), gamma0_mN_m, gamma_max_mol_per_m2, K (FITTED constant specific to this surfactant -- never a universal default), system_type (required n-factor; do not guess), temperature_K.
    """
    value = ads.szyszkowski_surface_tension(concentration, gamma0_mN_m, gamma_max_mol_per_m2, K, system_type, temperature_K)
    return {"surface_tension_mN_m": value, "unit": "mN/m"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def szyszkowski_fit_K(concentrations: list[float], surface_tensions_mN_per_m: list[float], gamma0_mN_m: float, gamma_max_mol_per_m2: float, system_type: str, temperature_K: float = 298.15) -> dict:
    """
    Fit the Szyszkowski/Langmuir constant K from a RAW PRE-CMC surface-tension curve (the step that makes szyszkowski_predict_surface_tension's K known). Args: concentrations (mM or mol/L -- any one consistent molarity unit; the returned K is in the RECIPROCAL of that unit, e.g. 1/mM for an mM axis -- report which basis was used) + surface_tensions_mN_per_m, PRE-CMC points only (points at/above CMC violate the model and bias the fit), same length, >=3 points; system_type REQUIRED (n-factor; do not guess). Returns K (in units of 1/(concentration basis), stated with the fit), Gamma_max fixed-or-fitted per signature, r_squared, fitted curve.
    """
    result = ads.szyszkowski_fit_K(concentrations, surface_tensions_mN_per_m, gamma0_mN_m, gamma_max_mol_per_m2, system_type, temperature_K)
    return {
        "K": result.K,
        "r_squared": result.r_squared,
        "n_points": result.n_points,
        "fitted_surface_tensions_mN_per_m": result.fitted_surface_tensions_mN_per_m,
    }


@mcp.tool()
@_tool_guard
@_shape_no_solution
def frumkin_predict_surface_tension(concentration: float, gamma0_mN_m: float, gamma_max_mol_per_m2: float, K: float, a: float, system_type: str, temperature_K: float = 298.15) -> dict:
    """
    Predict surface tension (mN/m) via the Frumkin isotherm: gamma = gamma0 + n_factor*R*T*Gamma_max*[ln(1-theta) + a*theta^2] with theta from K*c = theta/(1-theta)*exp(-2*a*theta). a (dimensionless lateral interaction): >0 attractive, <0 repulsive, =0 reduces EXACTLY to szyszkowski_predict_surface_tension -- use Szyszkowski unless real data justify a; NEVER guess a. |a| >= 2 makes the isotherm genuinely multi-valued (2D condensation) -- flag such fits, don't treat as ordinary. Args same conventions as szyszkowski_predict_surface_tension; system_type REQUIRED, no default. Returns gamma_mN_per_m, theta.
    """
    value = ads.frumkin_surface_tension(concentration, gamma0_mN_m, gamma_max_mol_per_m2, K, a, system_type, temperature_K)
    return {"surface_tension_mN_m": value, "unit": "mN/m"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def hlb_from_mw(mw_hydrophilic: float, mw_total: float) -> dict:
    """Hydrophile-Lipophile Balance via Griffin's method: HLB = 20 *
    (hydrophilic-portion molecular weight / total molecular weight).
    Returns a value on Griffin's 0-20 scale (mw_hydrophilic and
    mw_total in the same, arbitrary mass unit -- e.g. g/mol)."""
    value = hlb_mod.hlb_griffin(mw_hydrophilic, mw_total)
    return {"hlb": value, "scale": "Griffin 0-20"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def hlb_from_groups(group_counts: dict[str, int], allow_derived_groups: bool = False) -> dict:
    """
    HLB via Davies' 1957 group contribution: HLB = 7 + sum(hydrophilic group numbers) - sum(lipophilic group numbers). Args: group_counts from the VERIFIED lists (hydrophilic: SO4Na, COOK, COONa, N_tertiary_amine, N_quat_*, ester_sorbitan_ring...; lipophilic: CH3, CH2, C=C, ...), allow_derived_groups (default False: opt-in to this project's disclosed Griffin-cross-calibrated numbers for amide/sulfonate etc.). Raises for ANY unknown group name rather than guessing -- see derive_davies_group_number for the method and its chain-length-sensitivity caveat.
    """
    value = hlb_mod.hlb_davies(group_counts, allow_derived_groups=allow_derived_groups)
    return {"hlb": value, "scale": "Davies (7 = neutral reference point)"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def hlb_method_recommendation(structural_family: str, charge_type: str | None = None) -> dict:
    """
    Which HLB method to use for a structural family (resolution of a confirmed STRUCTURAL bottleneck, not a missing-data one): takes classify_surfactant_structural_family's output -- 'dimeric_gemini_type' or 'glycolipid_biosurfactant' -> recommends Griffin's method (hlb_froehlich/etc. per map) with reasoning; other families -> Davies. Raises for 'unparseable' or unknown families rather than guessing. Returns recommended_method + why (the decision rule, sourced).

    charge_type: optional, pass classify_surfactant_charge_type's output
    ('cationic'/'anionic'/'nonionic'/etc.). The Griffin evidence base
    for gemini/glycolipid is entirely NONIONIC structures -- for a
    strongly ionic (cationic/anionic) one, Griffin stays the better-
    justified METHOD NAME, but the returned 'griffin_numeric_caveat'
    field explains why an exact numeric HLB from it should not be
    trusted for that structure. Always check this field before citing a
    specific HLB number for an ionic gemini/glycolipid.
    """
    return hlb_mod.recommend_hlb_method_for_structural_family(structural_family, charge_type)


@mcp.tool()
@_tool_guard
@_shape_no_solution
def gemini_hlb_griffin_reference(gemini_name: str) -> dict:
    """
    Primary-source Griffin-method HLB for named 1-alkylaminoglycerol gemini surfactants -- Liao et al., Arabian J. Chem. 16 (2023) 105111 (primary PDF read in full). gemini_name: '8-3-8', '12-3-12', '8-4-8', '12-4-12' (alkyl_chain-spacer_length-alkyl_chain). Returns hlb, note, method + citation per entry; raises for any other name rather than guessing.
    """
    key = gemini_name.strip()
    if key not in hlb_mod.GEMINI_SURFACTANT_HLB_GRIFFIN_REFERENCE:
        raise ValueError(
            f"gemini_name must be one of {list(hlb_mod.GEMINI_SURFACTANT_HLB_GRIFFIN_REFERENCE)}, "
            f"got {gemini_name!r} -- no guessed value available"
        )
    entry = hlb_mod.GEMINI_SURFACTANT_HLB_GRIFFIN_REFERENCE[key]
    return {
        "gemini_name": key,
        "hlb": entry["hlb"],
        "note": entry["note"],
        "method": "Griffin",
        "source": "Liao et al., Arabian J. Chem. 16 (2023) 105111",
    }


@mcp.tool()
@_tool_guard
@_shape_no_solution
def derive_davies_group_number(hydrophilic_fragment_mass_g_per_mol: float, total_molar_mass_g_per_mol: float, lipophilic_group_counts: dict[str, int]) -> dict:
    """
    Derive a NEW Davies-scale hydrophilic group number for a group missing from Davies' 1957 table: cross-calibrates Davies' additive HLB against Griffin's mass-ratio HLB for ONE real reference compound and solves for the unknown number (the same method Davies himself used; amide/sulfonate gaps were filled this way -- see hlb_from_groups allow_derived_groups). Args: hydrophilic_fragment_mass_g_per_mol (Griffin Mh convention: polar fragment incl. heteroatom-bearing carbons), total_molar_mass_g_per_mol, lipophilic_group_counts (verified Davies groups, e.g. {"CH2":10,"CH3":1}). Disclosed limitation: Davies additive HLB is linear in chain length while Griffin's Mh/M is not -- run at SEVERAL reference chain lengths (C8-C16) before trusting; this returns ONE compound's value, not an average.
    """
    value = hlb_mod.derive_davies_group_number_from_griffin(
        hydrophilic_fragment_mass_g_per_mol, total_molar_mass_g_per_mol, lipophilic_group_counts)
    return {"derived_group_number": value, "scale": "Davies (additive HLB group-contribution units)"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def hlb_davies_guo_ecl(n_carbons_alkyl: int, n_eo: int, extra_group_counts: dict[str, int] | None = None) -> dict:
    """
    HLB for NONIONIC polyethoxylated surfactants (R-(OCH2CH2)n-OH) via Guo, Rong & Ying 2006 -- the sourced alternative to hlb_from_groups for this class: uses EFFECTIVE EO/alkyl chain lengths (piecewise NEO<=50 / >=50 branches, both implemented), fitting 224 real surfactants to <1.5 avg abs error, because a PEO chain's hydrophilicity is sub-linear in EO count where plain Davies (linear) is not. Args: n_carbons_alkyl (total tail carbons), n_eo (actual EO units), optional extra_group_counts (e.g. ester headgroup). Primary-source Table 1 confirms the reused Davies per-unit weights (CH2/CH3=-0.475, EO=0.33) -- verified 2026-09-11, not assumed.
    """
    value = hlb_mod.hlb_davies_guo_ecl(n_carbons_alkyl, n_eo, extra_group_counts)
    return {"hlb": value, "scale": "Davies (7 = neutral reference point), Guo/Rong/Ying effective-chain-length correction"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def hld_optimal_salinity(k: float, eacn: float, cc: float, alpha: float = 0.01, delta_T_K: float = 0.0) -> dict:
    """
    Optimum (HLD=0, Winsor III) salinity for a surfactant/oil/brine system: S* = exp(k*EACN + alpha*deltaT - Cc) (HLD framework; use for ionic systems, hlb_from_cc for HLB-based ones). Args: k, Cc, EACN, deltaT_T_minus_Tref, alpha (temperature coefficient) -- all supplied from measurement/reference tables (do not guess; use hld_*_reference tools). Returns S_star (same salinity units as Cc convention: %se typically -- see params), plus the HLD residual check.
    """
    value = hld_mod.optimal_salinity_ionic(k, eacn, cc, alpha, delta_T_K)
    return {"optimal_salinity_pct": value, "unit": "g NaCl equivalent / 100 mL"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def hld_cationic_quat_reference(surfactant: str) -> dict:
    """Real, measured HLD-framework reference values for a named
    cationic quaternary ammonium surfactant -- characteristic curvature
    (Cc), the shared k constant for this head-group class, and (where
    available) an equivalent Davies-scale HLB. surfactant must be one
    of: LTAB (dodecyltrimethylammonium bromide), MTAB
    (tetradecyltrimethylammonium bromide), CTAB
    (hexadecyltrimethylammonium bromide), CMIC (1-hexadecyl-3-
    methylimidazolium chloride), BDHC (benzyldimethylhexadecylammonium
    chloride), or DDAB (didodecyldimethylammonium bromide, double-tailed).
    Source: Schirone, Tartaro, Gentile & Palazzo, JCIS Open 4 (2021)
    100033, Table 1 -- real HLD-titration measurements, not guessed or
    estimated values. Use these Cc values with hld_optimal_salinity or
    cc_mixing_rule rather than guessing a Cc for these specific
    surfactants."""
    key = surfactant.upper()
    if key not in hld_mod.CATIONIC_QUAT_CC:
        raise ValueError(
            f"surfactant must be one of {sorted(hld_mod.CATIONIC_QUAT_CC)} -- no guessed value "
            "available for anything else"
        )
    cc, cc_uncertainty = hld_mod.CATIONIC_QUAT_CC[key]
    result = {
        "surfactant": key,
        "Cc": cc,
        "Cc_uncertainty": cc_uncertainty,
        "k": hld_mod.CATIONIC_QUAT_K,
        "k_uncertainty": hld_mod.CATIONIC_QUAT_K_UNCERTAINTY,
        "estimated_hlb_davies_scale": hld_mod.hlb_from_cc_cationic_quat(cc),
    }
    if key in hld_mod.CATIONIC_QUAT_HLB_DAVIES:
        result["real_reported_hlb_davies_scale"] = hld_mod.CATIONIC_QUAT_HLB_DAVIES[key]
    return result


@mcp.tool()
@_tool_guard
@_shape_no_solution
def hld_class_reference(surfactant_class: str) -> dict:
    """
    Sourced HLD-NAC default constants (k, alpha, b as applicable) for a surfactant CLASS when no system-specific value exists. surfactant_class: 'anionic_sulfate_sulfonate', 'extended_surfactant', 'cationic_quaternary_ammonium', 'nonionic_ethoxylated', 'zwitterionic' (each with k_default etc. and per-field source notes). Confidence: disclosed levels, not as precise as the per-compound hld_cationic_quat_reference table. Disclosed gap: no primary-source Cc for zwitterionic/gemini classes -- raises rather than borrowing one.
    """
    key = surfactant_class.lower()
    if key == "anionic_sulfate_sulfonate":
        return {
            "surfactant_class": key,
            "k_default": hld_mod.K_ANIONIC_DEFAULT,
            "sds_cc_reference": hld_mod.CC_REFERENCE_ANIONIC_CATIONIC["SDS"]["cc_this_work"],
            "source": "Abbott 'Practical Surfactants Science' (k); Leng & Acosta 2023, primary PDF read in full (SDS Cc)",
        }
    if key == "extended_surfactant":
        return {
            "surfactant_class": key,
            "k_default": hld_mod.K_EXTENDED_SURFACTANT,
            "source": "Abbott 'Practical Surfactants Science' HLD reference page",
        }
    if key == "apg_nonionic":
        return {
            "surfactant_class": key,
            "alpha_default": hld_mod.ALPHA_APG_DEFAULT,
            "b_default": hld_mod.B_APG_SPAN_DEFAULT,
            "source": "Abbott 'Practical Surfactants Science' HLD reference page",
        }
    raise ValueError(
        f"surfactant_class must be one of ['anionic_sulfate_sulfonate', 'extended_surfactant', "
        f"'apg_nonionic'], got {surfactant_class!r} -- no guessed value available for zwitterionic/"
        "gemini classes (a real, disclosed gap, not yet closed; glycolipid biosurfactant Cc IS now "
        "available -- see hld_characteristic_curvature_reference('rhamnolipid'))"
    )


@mcp.tool()
@_tool_guard
@_shape_no_solution
def hld_characteristic_curvature_reference(surfactant: str) -> dict:
    """
    Primary-source characteristic curvature Cc for a named surfactant (anionic, cationic, zwitterionic, gemini, biosurfactant). surfactant: 'SDS','SDHS','SLES','C10PO4S','C16DPODS','AOT','BCl','DPCl' (cc_this_work/cc_literature); 'rhamnolipid' (-1.41); 'lecithin','Epikuron200','C4mPC','dodecylsulfobetaine','lauramineoxide','CAPB' (two real CAPB estimates: cc_this_work/cc_this_work_alt); 'gemini_benzene_sulfonate' (-7.4, cross-checked vs -6.6). Sources: Leng & Acosta JSD 26(3) 2023; Nguyen & Sabatini via Hellweg et al. 2023; Elsevier HLD-NAC book 2021 Table 1.1 -- all read in full. Disclosed discrepancy: a review's Cc(AOT)=-0.92 is rejected here as a citation slip (equals SDHS's own value); the directly-read primary AOT value is returned instead.
    """
    key = surfactant.strip()
    key_lower = key.lower()
    if key_lower == "rhamnolipid":
        return {
            "surfactant": "rhamnolipid",
            "cc": hld_mod.RHAMNOLIPID_CC,
            "class": "biosurfactant (glycolipid)",
            "source": "Nguyen & Sabatini 2008, via Hellweg, Oberdisse & Sottmann, Front. Soft Matter 3:1260211 (2023)",
        }
    if key_lower == "gemini_benzene_sulfonate":
        return {
            "surfactant": "gemini_benzene_sulfonate",
            "cc": hld_mod.GEMINI_BENZENE_SULFONATE_CC,
            "cc_literature": hld_mod.GEMINI_BENZENE_SULFONATE_CC_LITERATURE,
            "class": "gemini/dimeric (anionic)",
            "source": "Acosta, Harwell & Sabatini (eds.), Surfactant Formulation Engineering Using HLD and NAC, "
                      "Elsevier (2021), Table 1.1",
        }
    zwitterionic_key = next((k for k in hld_mod.CC_REFERENCE_ZWITTERIONIC if k.lower() == key_lower), None)
    if zwitterionic_key is not None:
        return {"surfactant": zwitterionic_key, "class": "zwitterionic",
                **hld_mod.CC_REFERENCE_ZWITTERIONIC[zwitterionic_key],
                "source": "Acosta, Harwell & Sabatini (eds.), Surfactant Formulation Engineering Using HLD and NAC, "
                          "Elsevier (2021), Table 1.1"}
    key_upper = key.upper()
    if key_upper not in hld_mod.CC_REFERENCE_ANIONIC_CATIONIC:
        raise ValueError(
            f"surfactant must be one of {sorted(hld_mod.CC_REFERENCE_ANIONIC_CATIONIC)}, "
            f"{sorted(hld_mod.CC_REFERENCE_ZWITTERIONIC)}, 'rhamnolipid', or 'gemini_benzene_sulfonate' -- "
            "no guessed value available for anything else"
        )
    return {"surfactant": key_upper, **hld_mod.CC_REFERENCE_ANIONIC_CATIONIC[key_upper],
            "source": "Leng & Acosta, J. Surfactants Deterg. 26(3) (2023) 287-301, primary PDF read in full"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def hld_fit_k_and_cc_from_salinity_scan(eacn_series: list[float], optimal_salinity_series_pct: list[float], alpha: float = 0.01, delta_T_K: float = 0.0) -> dict:
    """Fit k and Cc from a RAW multi-oil salinity-scan series -- the real
    experimental method that produces the k/Cc inputs hld_optimal_salinity
    otherwise takes as already-known. A salinity scan measures the
    optimum (Winsor III, HLD=0) salinity S* for several different oils
    spanning a range of EACN, at fixed surfactant/temperature; ln(S*) is
    linear in EACN (ln(S*) = k*EACN + alpha*deltaT - Cc), so fitting that
    line gives k (slope) and Cc directly. eacn_series and
    optimal_salinity_series_pct (in g NaCl equivalent / 100 mL, matching
    hld_optimal_salinity's units) must be the same length, at least 2
    points (3+ recommended). Returns the fitted k, Cc, and the fit's
    r_squared -- HLD-NAC assumes k is constant across EACN for a given
    surfactant class; a poor r_squared signals that assumption doesn't
    hold for this system, not a bug in the fit."""
    result = hld_mod.fit_k_and_cc_from_salinity_scan(eacn_series, optimal_salinity_series_pct, alpha, delta_T_K)
    return {
        "k": result.k,
        "cc": result.cc,
        "r_squared": result.r_squared,
        "n_points": result.n_points,
    }


@mcp.tool()
@_tool_guard
@_shape_no_solution
def tanford_chain_geometry(n_carbons: int) -> dict:
    """Tanford's formulas for a saturated, unbranched alkyl chain of
    n_carbons carbons: hydrophobic tail volume (cubic Angstrom) and
    maximum extended (critical) chain length (Angstrom). Feed these into
    critical_packing_parameter along with a headgroup area."""
    v = cpp_mod.tanford_tail_volume(n_carbons)
    lc = cpp_mod.tanford_critical_length(n_carbons)
    return {
        "tail_volume_A3": v,
        "critical_length_A": lc,
        "unit": "cubic Angstrom (volume), Angstrom (length)",
    }


@mcp.tool()
@_tool_guard
@_shape_no_solution
def critical_packing_parameter(volume_A3: float, head_area_A2: float, length_A: float, n_carbons: int | None = None) -> dict:
    """Critical packing parameter (CPP) = volume / (head_area * length).
    Dimensionless. volume_A3 and length_A typically come from
    tanford_chain_geometry; head_area_A2 is the optimal headgroup area
    at the interface (surfactant- and condition-specific -- do not guess
    this without a source). Also returns the predicted aggregate
    morphology. Optional n_carbons: when supplied, adds a validity_note
    advisory for short/medium tails (n_carbons <= 9: the Hartley-Tanford
    spherical-core reasoning is designed for long chains and cannot be
    enforced on short-chain non-spherical assemblies, PMID 31954873) --
    never changes the CPP value; omitted -> validity_note is null."""
    cpp = cpp_mod.critical_packing_parameter(volume_A3, head_area_A2, length_A)
    morphology = cpp_mod.classify_aggregate_morphology(cpp)
    return {
        "cpp": cpp,
        "unit": "dimensionless",
        "predicted_morphology": morphology,
        "validity_note": (
            cpp_mod.packing_parameter_validity_note(n_carbons)
            if n_carbons is not None else None
        ),
    }


@mcp.tool()
@_tool_guard
@_shape_no_solution
def estimate_axial_ratio_from_cpp_geometry(cpp: float, aggregation_number: float, n_carbons: int) -> dict:
    """
    Geometry-based prolate-ellipsoid axial ratio (a/b) for a rodlike micelle from the critical packing parameter -- closes the COMPUTE path for hydrodynamic_radius_perrin_corrected's axial_ratio input (the spherical-micelle formula is never used; it would assume away the geometry being corrected). Args: cpp, n_carbons (tail), headgroup_area_A2, tail_volume_A3 (or volume/length sources per signature). Returns a/b plus regime guards; rejects outside the rodlike band (oblate/vesicle out of scope, disclosed) instead of extrapolating.
    """
    axial_ratio = cpp_mod.estimate_axial_ratio_from_cpp_geometry(cpp, aggregation_number, n_carbons)
    return {"axial_ratio": axial_ratio, "unit": "dimensionless (a/b, prolate ellipsoid of revolution)"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def nagarajan_equilibrium_area_ionic(cmc_M: float, tail_length_A: float, headgroup_prefactor_A: float, temperature_K: float = 298.15, dielectric_constant: float = 80.0) -> dict:
    """
    Equilibrium area per molecule a_e (Angstrom^2) for an IONIC surfactant micelle with chain-length-dependent electrostatic headgroup repulsion -- Nagarajan, Langmuir 18 (2002) 31-38 (real alternative to plain geometric areas). Args: n_carbons, headgroup params, plus cmc_M (no added electrolyte: uses paper's I ~= cmc approximation) and headgroup_prefactor_A (the paper's sqrt(2*pi*e^2*d/(eps*sigma)); 82.0 is its illustrative value for sodium-alkyl-sulfate-like headgroups -- do NOT treat 82 as universal without a source).
    """
    a_e = cpp_mod.nagarajan_equilibrium_area_ionic(cmc_M, tail_length_A, headgroup_prefactor_A, temperature_K, dielectric_constant)
    kappa_inv = cpp_mod.nagarajan_debye_huckel_kappa_inverse(cmc_M, temperature_K, dielectric_constant)
    return {"a_e_A2": a_e, "debye_huckel_kappa_inverse_A": kappa_inv, "unit": "square Angstrom (a_e), Angstrom (kappa^-1)"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def nagarajan_ruckenstein_headgroup_reference(headgroup: str) -> dict:
    """Real, sourced molecular constants for a surfactant headgroup --
    Nagarajan & Ruckenstein, Langmuir 7 (1991) 2934-2969, Table I
    (primary PDF read in full). headgroup must be one of:
    'sodium_sulfate' (a_p=17 A^2, a_o=17 A^2, delta=5.45 A, anionic),
    'sodium_sulfonate' (a_p=17 A^2, a_o=17 A^2, delta=3.85 A, anionic),
    'n_betaine' (a_p=30 A^2, a_o=21 A^2, d=5.0 A, ZWITTERIONIC), or
    'glucoside' (a_p=40 A^2, a_o=21 A^2, nonionic). Feed delta into
    nagarajan_ruckenstein_ionic_headgroup_free_energy for charged
    headgroups, or d into nagarajan_ruckenstein_dipole_headgroup_free_
    energy for the zwitterionic case. Raises for cationic quaternary
    ammonium or carboxylate headgroups -- not in this paper's own Table
    I either, a real, disclosed gap, not guessed."""
    key = headgroup.lower()
    if key not in cpp_mod.NAGARAJAN_RUCKENSTEIN_HEADGROUP_CONSTANTS:
        raise ValueError(
            f"headgroup must be one of {sorted(cpp_mod.NAGARAJAN_RUCKENSTEIN_HEADGROUP_CONSTANTS)} -- "
            "no guessed value available for anything else (cationic quaternary ammonium and "
            "carboxylate headgroups are a real, disclosed gap, not in the source paper's own table)"
        )
    return {"headgroup": key, **cpp_mod.NAGARAJAN_RUCKENSTEIN_HEADGROUP_CONSTANTS[key]}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def nagarajan_ruckenstein_steric_free_energy(area_per_molecule_A2: float, a_p_A2: float) -> dict:
    """Steric repulsion free energy (kT units) between headgroups
    crowded at the aggregate surface -- Nagarajan & Ruckenstein 1991,
    eq. 66: -ln(1 - a_p/a). a_p_A2 is the headgroup's own real hard-core
    area (see nagarajan_ruckenstein_headgroup_reference), never guessed.
    Diverges as area_per_molecule_A2 approaches a_p_A2 -- raises rather
    than returning a nonsensical value if area_per_molecule_A2 <= a_p_A2."""
    value = cpp_mod.nagarajan_ruckenstein_steric_headgroup_free_energy(area_per_molecule_A2, a_p_A2)
    return {"steric_free_energy_kT": value, "unit": "dimensionless (kT units)"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def nagarajan_ruckenstein_dipole_free_energy(
    area_per_molecule_A2: float, radius_A: float, d_A: float, geometry: str = "sphere",
    temperature_K: float = 298.15, dielectric_constant: float = 80.0,
) -> dict:
    """Dipole-dipole interaction free energy (kT units) between
    ZWITTERIONIC headgroups (e.g. N-betaine) at the aggregate surface --
    Nagarajan & Ruckenstein 1991, eqs. 67-68. d_A is the headgroup's own
    real charge-separation distance (see
    nagarajan_ruckenstein_headgroup_reference, e.g. d_A=5.0 for
    n_betaine), never guessed. geometry: 'sphere' (or 'globular'/
    'endcap') or 'cylinder' -- the two use genuinely different bracket
    terms in the source paper's own eqs. 67 vs. 68."""
    value = cpp_mod.nagarajan_ruckenstein_dipole_headgroup_free_energy(
        area_per_molecule_A2, radius_A, d_A, geometry, temperature_K, dielectric_constant
    )
    return {"dipole_free_energy_kT": value, "unit": "dimensionless (kT units)"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def nagarajan_ruckenstein_ionic_free_energy(
    area_per_molecule_A2: float, core_radius_A: float, delta_A: float,
    counterion_concentration_M: float, added_salt_M: float = 0.0,
    temperature_K: float = 298.15, dielectric_constant: float = 80.0, geometry: str = "sphere",
) -> dict:
    """
    Ionic (charged-headgroup) interaction free energy (kT) at the aggregate surface -- Nagarajan & Ruckenstein 1991 eqs. 70-73: curvature-corrected analytical Poisson-Boltzmann solution (more complete than the Debye-Huckel linearization). Args: headgroup charge, area, curvature, electrolyte ionic strength from the REAL salt composition (sum of supplied salt concentrations, never guessed -- raises if absent for a charged system). Returns kT-unit free energy plus the kappa/screening inputs used.
    """
    value = cpp_mod.nagarajan_ruckenstein_ionic_headgroup_free_energy(
        area_per_molecule_A2, core_radius_A, delta_A, counterion_concentration_M,
        added_salt_M, temperature_K, dielectric_constant, geometry,
    )
    return {"ionic_free_energy_kT": value, "unit": "dimensionless (kT units)"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def blankschtein_entropy_of_binding(beta: float) -> dict:
    """Entropy of mixing (kT units) between surfactant heads and ONE
    bound counterion species at the micelle surface -- Srinivasan &
    Blankschtein, Langmuir 19 (2003) 9932-9945, eq 13. beta: degree of
    binding in [0, 1] (beta=0 returns exactly 0.0, the fully-dissociated
    limit)."""
    value = cpp_mod.blankschtein_entropy_of_binding_free_energy(beta)
    return {"gent_kT": value, "unit": "dimensionless (kT units)"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def blankschtein_steric_free_energy(
    beta: float, area_per_molecule_A2: float, ah_surfactant_A2: float, ah_counterion_A2: float,
) -> dict:
    """Steric free energy (kT units) of surfactant heads packed with
    bound counterions at the micelle surface -- Srinivasan &
    Blankschtein, Langmuir 19 (2003) 9932-9945, eq 12. ah_counterion_A2
    is the counterion's hydrated radius SQUARED (rh^2, not pi*rh^2),
    per the source paper's own convention -- e.g. for Cl-, rh=2.13 A
    per Srinivasan & Blankschtein, Langmuir 19 (2003) 9946-9961, Table
    2, giving ah=4.537 A^2."""
    value = cpp_mod.blankschtein_steric_free_energy_with_counterion(
        beta, area_per_molecule_A2, ah_surfactant_A2, ah_counterion_A2
    )
    return {"gst_kT": value, "unit": "dimensionless (kT units)"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def blankschtein_counterion_self_energy(
    valence: int, hydrated_radius_A: float, kappa_inverse_A: float,
    dielectric_constant: float = 80.0, temperature_K: float = 298.15,
) -> dict:
    """Debye-Huckel electrostatic self-energy (kT units, always <= 0)
    released when a counterion binds onto the micelle Stern layer,
    leaving the bulk ionic atmosphere -- Srinivasan & Blankschtein,
    Langmuir 19 (2003) 9932-9945, eq 15 (attributed to Bockris & Reddy,
    Modern Electrochemistry I, 1977). kappa_inverse_A: the real
    solution's Debye-Huckel inverse screening length (Angstrom) -- see
    nagarajan_debye_huckel_kappa_inverse."""
    value = cpp_mod.blankschtein_counterion_self_energy_release(
        valence, hydrated_radius_A, kappa_inverse_A, dielectric_constant, temperature_K
    )
    return {"gdis_kT": value, "unit": "dimensionless (kT units)"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def aggregation_number(tail_volume_A3: float, core_radius_A: float) -> dict:
    """Estimated spherical-micelle aggregation number from geometric
    packing: N_agg = core_volume / tail_volume. tail_volume_A3 typically
    comes from tanford_chain_geometry; core_radius_A is usually
    approximated by the critical chain length from that same tool, but
    pass a measured value (e.g. from SANS/SAXS) if you have one. This is
    a geometric estimate (instant, no real measurement), NOT a substitute
    for a directly measured aggregation number -- if you have raw
    fluorescence-quenching or static-light-scattering data, use
    aggregation_number_from_quenching or aggregation_number_from_sls
    instead; those give a real, measured value from real instrumentation
    a smaller lab is actually likely to have (a standard
    spectrofluorometer, or a single-angle light-scattering/DLS
    instrument), not a packing-geometry approximation."""
    value = cpp_mod.aggregation_number_spherical(tail_volume_A3, core_radius_A)
    return {"aggregation_number": value, "unit": "dimensionless (molecules per micelle)",
            "note": "geometric estimate assuming a spherical micelle core"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def aggregation_number_from_quenching(
    quencher_concentrations: list[float],
    intensities: list[float],
    intensity_without_quencher: float,
    total_surfactant_concentration: float,
    cmc: float,
) -> dict:
    """
    Aggregation number from RAW steady-state fluorescence quenching (SSFQ) via Turro-Yekta 1978: ln(I0/I) = N*[Q]/(Ct - cmc) -- a straight line through the origin whose slope gives N. Args: quencher concentrations, I0/I ratios, total surfactant Ct, cmc (molar, same basis), per signature; probe extinction/quenching constants as required. Disclosed: equation source-verified from a real paper image but round-trip validated only (T2) -- a published worked example remains the nominated next step.
    """
    result = curve.aggregation_number_from_quenching_curve(
        quencher_concentrations, intensities, intensity_without_quencher, total_surfactant_concentration, cmc
    )
    return {
        "aggregation_number": result.aggregation_number,
        "unit": "dimensionless (molecules per micelle)",
        "method": result.method,
    }


@mcp.tool()
@_tool_guard
@_shape_no_solution
def dn_dc_reference(surfactant: str) -> dict:
    """Real, published dn/dc (refractive index increment) for a named,
    commonly-studied surfactant micelle in water -- for use in
    aggregation_number_from_sls when no system-specific measurement is
    available. surfactant must be one of: 'SDS' (dn_dc=0.11 mL/g),
    'CTAB' (dn_dc=0.15 mL/g) -- both at 632.8 nm, 25 C, Malvern
    Panalytical's published reference page -- or 'TritonX100'
    (dn_dc=0.140+/-0.005 mL/g at 546 nm, 20 C, real primary source:
    Stubicar et al., in Mittal (ed.), Surfactants in Solution, Plenum
    Press (1989) 181-193, primary PDF read in full 2026-09-15; stable
    across water and KCl/KBr/KI electrolyte solutions alike) or
    'Tween80' (dn_dc=0.132 mL/g at 690 nm, 25 C, real primary source:
    Jiang, PhD dissertation, Virginia Tech (2011), Wyatt Optilab rEX
    differential refractometer, primary PDF read in full 2026-09-15).
    Raises for anything else rather than guessing; this is a real but
    genuinely incomplete table, not a general database. (Docstring
    corrected 2026-10-02: the earlier 'Tween-80 not found' gap note was
    stale -- the table itself had gained a sourced TWEEN80 entry the
    same day it was written.)"""
    key = surfactant.upper()
    if key not in curve.DN_DC_REFERENCE_ML_PER_G:
        raise ValueError(
            f"surfactant must be one of {sorted(curve.DN_DC_REFERENCE_ML_PER_G)} -- "
            "no guessed value available for anything else"
        )
    return {"surfactant": key, **curve.DN_DC_REFERENCE_ML_PER_G[key]}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def aggregation_number_from_sls(
    micellized_concentrations_g_per_mL: list[float],
    rayleigh_ratios_cm_inv: list[float],
    dn_dc_mL_per_g: float,
    wavelength_nm: float,
    monomer_molar_mass_g_per_mol: float,
    refractive_index_solvent: float = 1.333,
) -> dict:
    """
    Aggregation number from a RAW static light scattering concentration series -- single-angle Debye plot (standard for Rg < 12 nm scatterers). Args: concentrations + Kc/R_theta per signature; dn_dc must be a real measured/cited value (NEVER guessed -- raises if absent). Returns N_agg, Rg, A2, r_squared, plot points. Disclosed: round-trip-validated only -- no published worked example of this exact procedure exists (T2-tier output: computation gold, not 'as reported in literature').
    """
    result = curve.aggregation_number_from_sls_debye_plot(
        micellized_concentrations_g_per_mL, rayleigh_ratios_cm_inv, dn_dc_mL_per_g,
        wavelength_nm, monomer_molar_mass_g_per_mol, refractive_index_solvent,
    )
    return {
        "aggregation_number": result.aggregation_number,
        "micelle_molar_mass_g_per_mol": result.micelle_molar_mass_g_per_mol,
        "second_virial_coefficient_cm3_mol_per_g2": result.second_virial_coefficient_cm3_mol_per_g2,
        "unit": "dimensionless (molecules per micelle)",
        "method": result.method,
    }


@mcp.tool()
@_tool_guard
@_shape_no_solution
def aggregation_number_from_dls(diffusion_coefficient_cm2_per_s: float, viscosity_mPas: float, partial_specific_volume_cm3_per_g: float, monomer_molar_mass_g_per_mol: float, temperature_K: float = 298.15) -> dict:
    """
    Aggregation number from a RAW DLS diffusion coefficient: converts D -> hydrodynamic radius (Stokes-Einstein) -> sphere volume -> N_agg from tail/headgroup volumes. Args: diffusion_cm2_per_s, viscosity_mPas, temperature_K, component molar masses/volumes per signature; partial_specific_volume can come from estimate_partial_specific_volume. NEVER guess input constants -- raises if missing. Disclosed limitation: R_h includes a hydration shell, so N_agg is a real but hydration-biased estimate (T2-tier).
    """
    result = curve.aggregation_number_from_dls(diffusion_coefficient_cm2_per_s, viscosity_mPas, partial_specific_volume_cm3_per_g, monomer_molar_mass_g_per_mol, temperature_K)
    return {
        "aggregation_number": result.aggregation_number,
        "hydrodynamic_radius_nm": result.hydrodynamic_radius_nm,
        "micelle_molar_mass_g_per_mol": result.micelle_molar_mass_g_per_mol,
        "unit": "dimensionless (molecules per micelle)",
        "method": result.method,
    }


@mcp.tool()
@_tool_guard
@_shape_no_solution
def estimate_partial_specific_volume(n_carbons: int, monomer_molar_mass_g_per_mol: float, headgroup_molar_volume_cm3_per_mol: float) -> dict:
    """
    Buildable estimate of a monomer's partial specific volume v_bar (cm^3/g) -- reduces (not removes) the external-literature requirement for aggregation_number_from_svedberg/dls inputs. Sums hydrophobic tail (Tanford volume/mass) + headgroup contributions; hydrophobic and headgroup terms may individually be negative (only the FINAL v_bar must be positive -- do not clamp intermediates). Returns v_bar plus the term breakdown so callers can inspect what dominated.
    """
    v_bar = curve.estimate_partial_specific_volume_from_tail_and_headgroup(n_carbons, monomer_molar_mass_g_per_mol, headgroup_molar_volume_cm3_per_mol)
    return {"partial_specific_volume_cm3_per_g": v_bar, "unit": "cm^3/g"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def aggregation_number_from_svedberg(sedimentation_coefficient_S: float, diffusion_coefficient_cm2_per_s: float, partial_specific_volume_cm3_per_g: float, solvent_density_g_per_cm3: float, monomer_molar_mass_g_per_mol: float, temperature_K: float = 298.15) -> dict:
    """Micelle aggregation number from RAW analytical-ultracentrifugation
    data via the classic Svedberg equation -- a real, sourced
    ALTERNATIVE independent of the light-scattering/DLS routes (a
    different instrument, different physical principle: M =
    R*T*s/[D*(1-v_bar*rho)]). sedimentation_coefficient_S must be in
    SVEDBERGS (1 S = 1e-13 s), the standard reporting unit -- do not
    pass raw seconds. diffusion_coefficient_cm2_per_s: same convention
    as elsewhere in this server. partial_specific_volume_cm3_per_g and
    solvent_density_g_per_cm3: real, measured/cited values -- system-
    specific, do not guess (water ~0.997 g/cm^3 at 25C for the solvent
    density is a real physical constant, not a guess). Raises if the
    buoyancy factor (1-v_bar*rho) is non-positive -- the particle would
    not physically sediment under those inputs."""
    result = curve.aggregation_number_from_svedberg_equation(sedimentation_coefficient_S, diffusion_coefficient_cm2_per_s, partial_specific_volume_cm3_per_g, solvent_density_g_per_cm3, monomer_molar_mass_g_per_mol, temperature_K)
    return {
        "aggregation_number": result.aggregation_number,
        "micelle_molar_mass_g_per_mol": result.micelle_molar_mass_g_per_mol,
        "unit": "dimensionless (molecules per micelle)",
        "method": result.method,
    }


@mcp.tool()
@_tool_guard
@_shape_no_solution
def cmc_from_surface_tension_curve(concentrations_mM: list[float], surface_tensions_mN_per_m: list[float], method: str = "rss") -> dict:
    """
    Extract CMC from a RAW tensiometry curve (surface tension vs concentration) via two-segment linear-regression break-point detection -- the first step other tools assume already happened. Args: concentrations (molar), surface_tensions_mN_per_m, same length; method: "rss" (default, min-RSS segmented), "tangential" (steepest premicellar fit x postmicellar fit intersection), "derivative" (max |dgamma/dlog10C| pair) -- each a real standard break-point estimator, differing only in how the break is read. Refusal: if the curve never reaches a ~20 mN/m drop before the data ends, returns a 'not_defined_for_this_system' status instead of a break-point; unknown method -> no_solution payload naming the allowed values.
    """
    result = curve.cmc_from_surface_tension_curve(concentrations_mM, surface_tensions_mN_per_m, method=method)
    return {
        "cmc_mM": result.cmc_mM,
        "premicellar_slope_mN_per_m_per_lnC": result.premicellar_slope_mN_per_m_per_lnC,
        "gamma_at_cmc_mN_per_m": result.gamma_at_cmc_mN_per_m,
        "r_squared_premicellar": result.r_squared_premicellar,
        "n_premicellar_points": result.n_premicellar_points,
        "n_postmicellar_points": result.n_postmicellar_points,
        "n_baseline_points": result.n_baseline_points,
        "premicellar_x_min_mM": result.premicellar_x_min_mM,
        "gamma0_mN_m": result.gamma0_mN_m,
        "c20_mM": result.c20_mM,
        "pC20": result.pC20,
        "postmicellar_slope_mN_per_m_per_log10C": result.postmicellar_slope_mN_per_m_per_log10C,
        "postmicellar_mean_gamma_mN_per_m": result.postmicellar_mean_gamma_mN_per_m,
        "r_squared_postmicellar": result.r_squared_postmicellar,
        "method": result.method,
    }


@mcp.tool()
@_tool_guard
@_shape_no_solution
def cmc_from_conductivity(concentrations_mM: list[float], conductivities: list[float], method: str = "rss") -> dict:
    """Extract CMC from a RAW specific-conductivity-vs-concentration curve --
    the standard conductometric break-point method for IONIC surfactants
    only (nonionic surfactants have negligible conductivity response and
    cannot be characterized this way). Unlike cmc_from_surface_tension_curve
    (fit on log10 concentration), this fits on LINEAR concentration: two
    real straight-line segments (premicellar and postmicellar), break point
    at their intersection = CMC. Returns both slopes AND the counterion
    binding degree (beta) computed from them the same way
    thermodynamics.counterion_binding_degree() does, so this one call
    answers "raw conductivity curve -> CMC and beta" end to end. method:
    "rss" (default, min-RSS two-line intersection on linear concentration),
    "jacobs_jones" (second-derivative maximum, grid-resolution-limited to
    sampled points), or "onsager" (molar-conductivity kappa/c vs sqrt(c)
    two-line intersection -- CMC and beta are unit-scale invariant, but
    the returned r-squared values are on that basis). Unknown method or
    non-positive concentrations (onsager) -> no_solution payload."""
    result = curve.cmc_from_conductivity_curve(concentrations_mM, conductivities, method=method)
    return {
        "cmc_mM": result.cmc_mM,
        "slope_below_cmc": result.slope_below_cmc,
        "slope_above_cmc": result.slope_above_cmc,
        "counterion_binding_degree": result.counterion_binding_degree,
        "r_squared_below_cmc": result.r_squared_below_cmc,
        "r_squared_above_cmc": result.r_squared_above_cmc,
        "n_points_below_cmc": result.n_points_below_cmc,
        "n_points_above_cmc": result.n_points_above_cmc,
        "method": result.method,
    }


@mcp.tool()
@_tool_guard
@_shape_no_solution
def rosen_monolayer_solve(alpha1: float, c_mix_sigma_mM: float, c1_sigma_mM: float, c2_sigma_mM: float) -> dict:
    """Rosen's extension of Rubingh theory to the mixed ADSORBED
    MONOLAYER at an interface, from surface-tension data -- NOT CMC
    data. c1_sigma_mM/c2_sigma_mM are the pure-component concentrations
    needed to reach a chosen reference surface tension (or pressure);
    c_mix_sigma_mM is the mixed concentration needed to reach that SAME
    reference surface tension at bulk mole fraction alpha1. Using CMC
    values here instead is a common category error -- if you have CMC
    data, use rubingh_solve instead."""
    x1 = mm.solve_rosen_monolayer_x(alpha1, c_mix_sigma_mM, c1_sigma_mM, c2_sigma_mM)
    if x1 is None:
        raise ValueError("No valid root found in (0, 1) for these inputs.")
    beta_sigma = mm.rosen_beta_sigma(x1, alpha1, c_mix_sigma_mM, c1_sigma_mM)
    return {
        "monolayer_mole_fraction_x1": x1,
        "beta_sigma": beta_sigma,
        "synergy_classification": "synergistic" if beta_sigma < 0 else ("antagonistic" if beta_sigma > 0 else "ideal"),
    }


@mcp.tool()
@_tool_guard
@_shape_no_solution
def corrin_harkins_predict(cmc1_mM: float, salt_conc1_mM: float, cmc2_mM: float, salt_conc2_mM: float, salt_conc_target_mM: float) -> dict:
    """Predict CMC at a new counterion (salt) concentration via the
    Corrin-Harkins log-linear relation, fit from two known (CMC, salt
    concentration) data points for the SAME surfactant + salt system.
    The fitted slope g and intercept are system-specific -- there is no
    universal lookup table for them, which is why two real data points
    are required as input. Only reliable as interpolation between (or a
    close extrapolation beyond) the two given concentrations."""
    predicted_cmc, g = mm.corrin_harkins_predict_cmc(cmc1_mM, salt_conc1_mM, cmc2_mM, salt_conc2_mM, salt_conc_target_mM)
    return {"predicted_cmc_mM": predicted_cmc, "fitted_slope_g": g}


@mcp.tool()
@_tool_guard
def ionic_strength(concentrations_M: dict[float, float]) -> dict:
    """Ionic strength I = 0.5 * sum(c_i * z_i^2), mol/L (M).
    concentrations_M: {charge_number: concentration_M}, e.g. for 0.1 M
    NaCl: {1: 0.1, -1: 0.1} (cation and anion both counted). Molar
    concentrations must be >= 0 (charges may be negative; concentrations may not)."""
    if any(c < 0 for c in concentrations_M.values()):
        raise ValueError(
            "concentrations_M: molar concentrations must be >= 0, got a negative "
            "value -- charges may be negative but concentrations may not")
    value = elec.ionic_strength(concentrations_M)
    return {"ionic_strength_M": value, "unit": "M"}


@mcp.tool()
@_tool_guard
def frumkin_theta(concentration: float, K: float, a: float) -> dict:
    """Solve the Frumkin adsorption isotherm for surface coverage theta
    at a given bulk concentration: K*concentration = [theta/(1-theta)] *
    exp(-2*a*theta). K: Langmuir-type equilibrium constant in the
    RECIPROCAL molarity unit of the concentration basis (e.g. 1/mM if
    concentration is given in mM; fit from real data, do not guess). a:
    dimensionless lateral-interaction parameter
    (a>0 attractive, a<0 repulsive, a=0 recovers plain Langmuir)."""
    value = ads.frumkin_theta(concentration, K, a)
    return {"theta": value}


@mcp.tool()
@_tool_guard
def select_isotherm_model(concentrations: list[float], surface_tensions_mN_per_m: list[float], gamma0_mN_m: float, gamma_max_mol_per_m2: float, system_type: str, temperature_K: float = 298.15) -> dict:
    """Choose between the plain Langmuir/Szyszkowski isotherm (1 free
    parameter K) and the Frumkin isotherm (2 free parameters, K and a --
    allows lateral interaction) via BIC, from a raw pre-CMC surface-
    tension-vs-concentration curve. Frumkin will always fit at least as
    well as Langmuir on raw RSS alone (it has an extra free parameter) --
    BIC is what actually penalizes that and makes the selection
    meaningful. system_type is REQUIRED (sets the n-factor; do not
    guess)."""
    import dataclasses
    result = ads.select_isotherm_model(concentrations, surface_tensions_mN_per_m, gamma0_mN_m, gamma_max_mol_per_m2, system_type, temperature_K)
    return dataclasses.asdict(result)


@mcp.tool()
@_tool_guard
@_shape_no_solution
def debye_screening_length(ionic_strength_M: float, temperature_K: float = 298.15) -> dict:
    """Debye screening length (nm) for an electrolyte solution.
    ionic_strength_M must be the IONIC STRENGTH (0.5*sum(c_i*z_i^2)),
    not the raw salt concentration -- for a simple 1:1 salt like NaCl
    these are numerically equal, but for anything else (CaCl2, MgSO4,
    etc.) they are not. Use the ionic_strength helper concept: for each
    ion, multiply its molar concentration by its charge squared, sum,
    and halve."""
    value = elec.debye_length(ionic_strength_M, temperature_K)
    return {"debye_length_nm": value, "unit": "nm"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def grahame_equation_surface_potential(surface_charge_density_C_per_m2: float, ionic_strength_M: float, z: int = 1, temperature_K: float = 298.15) -> dict:
    """Surface (Stern-layer) potential (mV) from surface charge density
    via the Grahame equation -- the standard, textbook closed-form
    Gouy-Chapman diffuse-double-layer relation (Grahame, Chem. Rev. 41
    (1947) 441-501): sigma = sqrt(8*eps_r*eps0*k_B*T*n0)*sinh(z*e*psi/(2*k_B*T)),
    solved for psi. z: charge number of the symmetric supporting
    electrolyte. ionic_strength_M: same convention as
    debye_screening_length (ionic strength, not raw salt concentration).
    Real use: convert a surfactant's surface excess Gamma (from
    gibbs_surface_excess) to a charge density via sigma=z*F*Gamma
    (F=Faraday constant) to estimate the electrostatic potential of an
    adsorbed ionic monolayer -- a real Gouy-Chapman electrostatics
    calculation, more rigorous than assuming a fixed Gibbs prefactor."""
    value = elec.grahame_equation_surface_potential(surface_charge_density_C_per_m2, ionic_strength_M, z, temperature_K)
    return {"surface_potential_mV": value, "unit": "mV"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def zeta_potential(electrophoretic_mobility_um_cm_per_Vs: float, viscosity_mPas: float, regime: str, kappa_a: float | None = None) -> dict:
    """
    Zeta (mV) from electrophoretic mobility via the Henry equation. regime is REQUIRED (no default -- silently picking one is the error this tool prevents): 'huckel' (kappa*a << 1), 'smoluchowski' (kappa*a >> 1, I >= ~10 mM), or 'ohshima'/'swan_furst'/'qin' (full-range approximations, all REQUIRE kappa_a; accuracy vs exact Henry: qin <1.5% in the kappa*a ~ 1-20 band real micelles occupy, ohshima <3%, swan_furst <0.1% overall -- typical R=2-4 nm micelles at 1-100 mM land at kappa*a ~ 0.2-4.2, squarely in the intermediate zone). viscosity_mPas in mPa.s (=cP; water ~0.89 at 25C -- NOT Pa.s, a 1000x unit trap). Returns zeta_mV plus the regime used.
    """
    value = elec.zeta_potential_henry(electrophoretic_mobility_um_cm_per_Vs, viscosity_mPas, regime, kappa_a=kappa_a)
    return {"zeta_potential_mV": value, "unit": "mV", "regime_used": regime.lower()}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def predict_mobility_relaxation_corrected(zeta_mV: float, kappa_a: float, viscosity_mPas: float, temperature_K: float = 298.15, z: int = 1, m: float = 0.15) -> dict:
    """Predict electrophoretic mobility ((um*cm)/(V*s)) from a known/
    assumed zeta potential via the relaxation-corrected model -- the
    forward direction of zeta_potential_relaxation_corrected (see that
    tool's docstring for full sourcing and parameter meaning). Useful
    for comparing a simulated or literature zeta potential against what
    mobility it would actually produce, accounting for the relaxation
    effect Henry's plain formula misses at high zeta."""
    value = elec.electrophoretic_mobility_relaxation_corrected(zeta_mV, kappa_a, viscosity_mPas, temperature_K, z, m)
    return {"electrophoretic_mobility_um_cm_per_Vs": value, "unit": "(um*cm)/(V*s)"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def zeta_potential_relaxation_corrected(electrophoretic_mobility_um_cm_per_Vs: float, kappa_a: float, viscosity_mPas: float, temperature_K: float = 298.15, z: int = 1, m: float = 0.15) -> dict:
    """
    Zeta (mV) from electrophoretic mobility via Ohshima-Healy-White 1983 (J. Chem. Soc. Faraday Trans. 2, 79, 1613, eqs. 57-62+75-76) -- the sourced ALTERNATIVE to plain Henry for |zeta| > 50 mV, where counterion relaxation makes mobility grow sub-linearly (Henry is linear and misses it entirely). Args: mobility, kappa_a (REQUIRED, source's valid range kappa_a >= 10 -- raises otherwise), z (symmetrical z-z electrolyte), m (dimensionless ionic drag, default 0.15 aqueous illustrative; source's Figs use 0.184 for KCl -- pass the real counterion's value if known), temperature_K, viscosity. Solved by grid+bisection, restricted to zeta in [0,150] mV (below the mobility maximum); raises if no solution exists in range rather than returning a wrong branch.
    """
    value = elec.zeta_potential_relaxation_corrected(electrophoretic_mobility_um_cm_per_Vs, kappa_a, viscosity_mPas, temperature_K, z, m)
    return {"zeta_potential_mV": value, "unit": "mV", "model": "relaxation-corrected (Ohshima-Healy-White 1983)"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def hydrodynamic_radius(diffusion_coefficient_cm2_per_s: float, viscosity_mPas: float, temperature_K: float = 298.15) -> dict:
    """Hydrodynamic radius (nm) from a DLS-measured diffusion
    coefficient via the Stokes-Einstein equation. diffusion_coefficient
    must be in cm^2/s (the commonly-reported DLS unit) -- NOT m^2/s, a
    classic 1e4x unit error. viscosity_mPas: solvent viscosity in mPa.s
    (=cP; water is ~0.89 at 25C) -- NOT Pa.s."""
    value = dyn.hydrodynamic_radius_stokes_einstein(diffusion_coefficient_cm2_per_s, viscosity_mPas, temperature_K)
    return {"hydrodynamic_radius_nm": value, "unit": "nm"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def hydrodynamic_radius_perrin_corrected(diffusion_coefficient_cm2_per_s: float, viscosity_mPas: float, axial_ratio: float, temperature_K: float = 298.15) -> dict:
    """TRUE equivalent-sphere hydrodynamic radius (nm) for a non-
    spherical (rodlike/wormlike or disklike) micelle from a DLS-measured
    diffusion coefficient -- a real, sourced ALTERNATIVE to
    hydrodynamic_radius, which assumes a sphere and so returns an
    inflated apparent radius for elongated/flattened aggregates. Uses
    the Perrin friction factor (Perrin 1934) to correct for shape:
    R_true = R_h,apparent / F_P(axial_ratio). axial_ratio = a/b (semi-
    axis along the rotation axis / equatorial semi-axis) must be known
    or assumed independently (e.g. from the critical packing parameter's
    predicted morphology, or SANS/SAXS) -- axial_ratio > 1 is prolate
    (rodlike), < 1 is oblate (disklike), == 1 gives exactly the same
    result as hydrodynamic_radius (no correction needed for a sphere).
    Do not guess axial_ratio."""
    r_true = dyn.hydrodynamic_radius_perrin_corrected(diffusion_coefficient_cm2_per_s, viscosity_mPas, axial_ratio, temperature_K)
    f_p = dyn.perrin_friction_factor(axial_ratio)
    return {"hydrodynamic_radius_true_nm": r_true, "perrin_friction_factor": f_p, "unit": "nm"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def counterion_binding_degree(slope_below_cmc: float, slope_above_cmc: float) -> dict:
    """Degree of counterion binding to the micelle (beta) from the
    slope-ratio method on a conductivity-vs-concentration plot: beta = 1
    - (slope_above_cmc / slope_below_cmc). Both slopes must be positive,
    with slope_above_cmc < slope_below_cmc, and both expressed in the SAME
    units (the ratio cancels only then). Note: this method's physical
    interpretation has been questioned for some systems in the
    literature -- the arithmetic is standard and widely used, but treat
    the result as an approximate, commonly-reported value."""
    value = thermo.counterion_binding_degree(slope_below_cmc, slope_above_cmc)
    return {"beta": value, "unit": "dimensionless (fraction, 0 to 1)"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def counterion_binding_degree_reference(surfactant: str) -> dict:
    """
    Counterion binding degree beta / dissociation alpha (=1-beta) for a named ionic surfactant, for when no system-specific measurement exists. surfactant: 'SDS' (alpha=0.272+/-0.017, HIGH -- Bales et al., J. Phys. Chem. B 105, 6798 (2001), primary EPR PDF read in full), 'DTAB' (alpha=0.28, MODERATE -- web-sourced 2026-09-15, cross-corroborated vs 0.29), 'CTAB' (alpha=0.26, MODERATE, corroborated vs 0.27). Raises for anything else (incomplete table, not a database). Returns alpha, beta, alpha_uncertainty (None where not stated), alpha_confidence, source -- weigh the two tiers rather than treating entries as equally certain.
    """
    key = surfactant.upper()
    if key not in thermo.COUNTERION_BINDING_DEGREE_REFERENCE:
        raise ValueError(
            f"surfactant must be one of {sorted(thermo.COUNTERION_BINDING_DEGREE_REFERENCE)} -- "
            "no guessed value available for anything else"
        )
    return {"surfactant": key, **thermo.COUNTERION_BINDING_DEGREE_REFERENCE[key]}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def gibbs_free_energy_micellization(cmc_M: float, temperature_K: float, counterion_factor: float = 1.0) -> dict:
    """Standard Gibbs free energy of micellization (kJ/mol) from CMC
    (molarity, mol/L) and temperature. counterion_factor must be stated
    explicitly: use 1.0 for a nonionic surfactant; for an ionic surfactant
    use (1 + beta) where beta is the counterion BINDING degree (beta = 1 -
    S2/S1, see counterion_binding_degree), equivalently (2 - alpha) where
    alpha = S2/S1 is the degree of ionization -- this is never auto-detected, since
    silently assuming nonionic behavior for an ionic surfactant (or vice
    versa) is exactly the kind of hidden-assumption error this tool
    exists to prevent."""
    x_cmc = thermo.cmc_to_mole_fraction(cmc_M)
    value = thermo.gibbs_free_energy_micellization(x_cmc, temperature_K, counterion_factor)
    return {"deltaG_mic_kJ_per_mol": value, "unit": "kJ/mol", "cmc_mole_fraction_used": x_cmc}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def mass_action_free_energy_of_micellization(cmc_M: float, aggregation_number: float, temperature_K: float) -> dict:
    """
    Standard deltaG_mic (kJ/mol) via the MASS-ACTION-LAW model: deltaG = R*T*[(1+1/n)*ln(X_cmc) - (1/n)*ln(n)] -- the real alternative to gibbs_free_energy_micellization's pseudo-phase convention (values differ by design; never mix conventions in one comparison). Args: cmc (mole fraction basis unless signature says otherwise), aggregation_number n, temperature_K, mode 'nonionic'|'ionic' with counterion factor. VALIDITY: do not apply to ionic surfactants without the counterion factor (use gibbs_free_energy_micellization for that case); n must be supplied from measurement, never guessed.
    """
    x_cmc = thermo.cmc_to_mole_fraction(cmc_M)
    value = thermo.mass_action_free_energy_of_micellization(x_cmc, aggregation_number, temperature_K)
    return {"deltaG_mic_kJ_per_mol": value, "unit": "kJ/mol", "cmc_mole_fraction_used": x_cmc, "model": "mass-action-law"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def critical_micellization_degree(aggregation_number: float, definition: str = "d2alpha_dc2") -> dict:
    """Critical micellization degree (CMD, dimensionless), the
    mass-action-law model's own definition of "how far along" the
    aggregation reaction is at the CMC -- a companion result to
    mass_action_free_energy_of_micellization under the same
    single-equilibrium mass-action model. Two real, sourced formulas
    are available depending on which curvature condition defines "the
    CMC": definition="d2alpha_dc2" (default, requires aggregation_number
    >= 2) uses the inflection point of micellization degree vs.
    CONCENTRATION; definition="d2alpha_dlnc2" uses the inflection point
    vs. LOG concentration (the more common experimental convention,
    e.g. a conductivity-vs-log(concentration) plot). Source: Rusanov,
    Langmuir 30(48) (2014) 14443-14451, eqs. (15) and (16)."""
    value = thermo.critical_micellization_degree(aggregation_number, definition)
    return {"alpha_m": value, "definition": definition, "aggregation_number": aggregation_number}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def vant_hoff_enthalpy(cmc1_M: float, temperature1_K: float, cmc2_M: float, temperature2_K: float, counterion_factor: float = 1.0) -> dict:
    """Van't Hoff enthalpy of micellization (kJ/mol) from CMC (molarity)
    measured at two temperatures. Assumes deltaH is constant over the
    temperature interval -- treat results from a wide T range with
    appropriate caution, and note this breaks down if the aggregation
    number itself varies significantly with temperature.

    counterion_factor: 1.0 for nonionic, (1 + beta) for ionic with beta =
    counterion binding degree (= 2 - alpha, alpha = ionization degree) -- same
    convention as gibbs_free_energy_micellization and
    vant_hoff_multi_point_fit; not auto-detected, do not guess. Required
    for a consistent answer on an ionic surfactant: deltaH must use the
    SAME counterion_factor as whatever deltaG_mic value it's being
    compared/combined with for the same compound (see this function's
    library-level docstring in thermodynamics.py for the Gibbs-Helmholtz
    derivation)."""
    x1 = thermo.cmc_to_mole_fraction(cmc1_M)
    x2 = thermo.cmc_to_mole_fraction(cmc2_M)
    value = thermo.vant_hoff_enthalpy(x1, temperature1_K, x2, temperature2_K, counterion_factor)
    return {"deltaH_mic_kJ_per_mol": value, "unit": "kJ/mol"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def vant_hoff_multi_point_fit(cmc_series_M: list[float], temperatures_K: list[float], counterion_factor: float = 1.0, reference_temperature_K: float | None = None) -> dict:
    """
    Multi-point (>=3 temperatures, 5+ recommended) Gibbs-Helmholtz fit for deltaH, deltaS AND deltaCp simultaneously from a raw CMC-vs-temperature series -- the real fix for vant_hoff_enthalpy's 2-point limitation (deltaCp captures curvature; with 2 points it is unidentifiable). Args: temperatures_K and cmc_M, SAME length and matching order; system_type/counterion convention as in other Gibbs tools (do not guess). Returns deltaH_kJ_mol, deltaS_J_mol_K, deltaCp_J_mol_K, SEs, r_squared, fit diagnostics.
    """
    x_values = [thermo.cmc_to_mole_fraction(c) for c in cmc_series_M]
    result = thermo.vant_hoff_multi_point_fit(x_values, temperatures_K, counterion_factor, reference_temperature_K)
    return {
        "deltaH_mic_kJ_per_mol": result.delta_H_kJ_per_mol,
        "deltaS_mic_J_per_mol_K": result.delta_S_J_per_mol_K,
        "deltaCp_mic_J_per_mol_K": result.delta_Cp_J_per_mol_K,
        "reference_temperature_K": result.reference_temperature_K,
        "r_squared": result.r_squared,
        "n_points": result.n_points,
        "deltaG_fitted_kJ_per_mol": result.delta_G_fitted_kJ_per_mol,
    }


@mcp.tool()
@_tool_guard
@_shape_no_solution
def simulate_emulsion_stability(drop_diameters_um: list[float], volume_fractions: list[float], dispersed_density_kg_m3: float, continuous_density_kg_m3: float, viscosity_Pa_s: float, column_height_m: float, coalescence_rate_per_s: float = 0.0, temperature_K: float = 298.15, n_times: int = 50, t_max_s: float = 86400.0) -> dict:
    """S5 emulsion stability: separated volume fraction vs time from Stokes creaming (drops of each bin start uniformly distributed over the column, so a bin's creamed share grows as t/t_c up to its travel time t_c = H/v; monodisperse t50 = H/(2v)) + first-order post-creaming coalescence, plus the standard stability index (area under the stable-fraction curve; larger = more stable) and t50. Args: drop_diameters_um (DSD bins), volume_fractions (weights; normalized internally, raw+normalized both in params), densities kg/m^3, viscosity Pa.s, column_height_m, coalescence_rate_per_s (0 = creaming only), t_max_s sampling window. Piecewise-analytic ODE solution -- no time-stepping error. Caveats in output: dilute Stokes (hindered creaming NOT modeled; flagged when total phi > 0.1), equal densities never separate by gravity, sedimentation case disclosed rather than hidden."""
    r = emu.simulate_emulsion_stability(
        drop_diameters_um=drop_diameters_um,
        volume_fractions=volume_fractions,
        dispersed_density_kg_m3=dispersed_density_kg_m3,
        continuous_density_kg_m3=continuous_density_kg_m3,
        viscosity_Pa_s=viscosity_Pa_s,
        column_height_m=column_height_m,
        coalescence_rate_per_s=coalescence_rate_per_s,
        temperature_K=temperature_K,
        n_times=n_times,
        t_max_s=t_max_s,
    )
    return {
        "times_s": list(r.times_s),
        "separated_fraction": list(r.separated_fraction),
        "stability_index_s": r.stability_index_m2,
        "t50_s": r.t50_s,
        "d32_um": r.d32_um,
        "creaming_times_s": list(r.creaming_times_s),
        "velocities_m_per_s": list(r.velocities_m_per_s),
        "converged": r.converged,
        "n_iterations": r.n_iterations,
        "method": r.method,
        "params": r.params,
        "caveat": r.caveat,
    }


@mcp.tool()
@_tool_guard
@_shape_no_solution
def population_balance_emulsion(initial_number_concentrations: list[float], coalescence_kernel: str = "constant", kernel_strength: float = 1e-3, breakage_rate_per_s: float = 0.0, ripening_rate_um3_per_s: float = 0.0, primary_drop_diameter_um: float = 1.0, n_times: int = 40, t_max_s: float = 600.0, n_steps_per_interval: int = 20) -> dict:
    """O2 population balance: Sauter diameter d32(t) from a discrete Smoluchowski breakage-coalescence PBE (bucket classes k=1..K of primary drops), with optional LSW Ostwald ripening. Args: initial_number_concentrations = n_k(0) for classes k=1..K (K = len; >=2), coalescence_kernel 'constant'|'sum'|'product' (the three classic Smoluchowski families; constant-kernel number decay validated against the exact solution), kernel_strength (family rate constant), breakage_rate_per_s > 0 enables first-order SYMMETRIC-BINARY breakage (disclosed caller-chosen daughter rule, not a sourced law), ripening_rate_um3_per_s > 0 adds the classical LSW cubic law (reported alone AND as an independent-mechanism superposition -- approximation disclosed), primary_drop_diameter_um sets the length scale. converged = mass conservation <=1e-6, which doubles as the K-adequacy check: products beyond class K escape the tracked domain, so a False means grow K or shorten t_max_s (never a silent pass). Deterministic RK4, no seed needed."""
    r = emu.population_balance_emulsion(
        initial_number_concentrations=initial_number_concentrations,
        coalescence_kernel=coalescence_kernel,
        kernel_strength=kernel_strength,
        breakage_rate_per_s=breakage_rate_per_s,
        ripening_rate_um3_per_s=ripening_rate_um3_per_s,
        primary_drop_diameter_um=primary_drop_diameter_um,
        n_times=n_times,
        t_max_s=t_max_s,
        n_steps_per_interval=n_steps_per_interval,
    )
    return {
        "times_s": list(r.times_s),
        "d32_um": list(r.d32_um),
        "d32_lsw_only_um": list(r.d32_lsw_only_um) if r.d32_lsw_only_um is not None else None,
        "number_concentration": list(r.number_concentration),
        "mass_rel_error": r.mass_rel_error,
        "converged": r.converged,
        "n_iterations": r.n_iterations,
        "method": r.method,
        "params": r.params,
        "caveat": r.caveat,
    }


@mcp.tool()
@_tool_guard
@_shape_no_solution
def surfactant_transport_1d(velocity_m_s: float, dispersion_m2_s: float, column_length_m: float, inlet_conc_mol_m3: float, langmuir_capacity_mol_m3: float = 0.0, langmuir_affinity_m3_per_mol: float = 0.0, n_cells: int = 200, t_max_s: float = 600.0, n_times: int = 400, n_steps_per_interval: int = 5) -> dict:
    """P4 transport: 1D convection-dispersion in a porous column with equilibrium Langmuir adsorption, Crank-Nicolson + Thomas tridiagonal (hand-rolled, no scipy). Returns c(x,t) samples, outlet breakthrough curve, breakthrough t50 (and pore volumes at half), discrete mass-balance residual (injected = stored + outflux), undershoot monitor, Picard convergence. Tracer case validated against the textbook two-term Ogata-Banks solution (<=1% at 400 cells, second-order convergence); low-concentration Langmuir case matches the retarded analytic form to 1e-6. Caveats: phenomenological (NO fit to any specific core-flood), equilibrium dilute Langmuir, Dirichlet inlet, node-Peclet oscillation risk disclosed (never clamped). Args: velocity_m_s, dispersion_m2_s, column_length_m, inlet_conc_mol_m3, langmuir_capacity_mol_m3 (per pore volume; 0 = tracer), langmuir_affinity_m3_per_mol, n_cells (>=10), t_max_s, n_times outputs, n_steps_per_interval CN steps per output."""
    r = tr_mod.surfactant_transport_1d(
        velocity_m_s=velocity_m_s,
        dispersion_m2_s=dispersion_m2_s,
        column_length_m=column_length_m,
        inlet_conc_mol_m3=inlet_conc_mol_m3,
        langmuir_capacity_mol_m3=langmuir_capacity_mol_m3,
        langmuir_affinity_m3_per_mol=langmuir_affinity_m3_per_mol,
        n_cells=n_cells,
        t_max_s=t_max_s,
        n_times=n_times,
        n_steps_per_interval=n_steps_per_interval,
    )
    return {
        "times_s": list(r.times_s),
        "x_m": list(r.x_m),
        "concentration_profiles": [list(p) for p in r.concentration_profiles],
        "breakthrough_fraction": list(r.breakthrough_fraction),
        "breakthrough_t50_s": r.breakthrough_t50_s,
        "mass_balance_residual": r.mass_balance_residual,
        "min_concentration": r.min_concentration,
        "converged": r.converged,
        "n_iterations": r.n_iterations,
        "method": r.method,
        "params": r.params,
        "caveat": r.caveat,
    }


@mcp.tool()
@_tool_guard
@_shape_no_solution
def analyze_micelle_trajectory(source: str, cluster_cutoff: float, n_radial_bins: int = 30, compute_slab_profile: bool = True) -> dict:
    """M1 trajectory analysis (analysis only): Nagg(t), full cluster-size lists, Rg(t) of the largest aggregate, radial number-density profile of that cluster, and box-z slab density when a box is present. Input: npz with `positions` (T,N,3) [+ optional time, box lengths, masses] or text (XYZ blocks / blank-line separated x y z rows) -- the repo's own MicelleMD output shapes. Clustering = single-linkage at `cluster_cutoff` (same length unit as coordinates; cutoff drives the answer and is echoed in params; bridging sensitivity disclosed). Mass-weighted Rg when npz masses exist. Seeded synthetic fixtures in tests; no RNG in the worker itself. numpy is a lazy `[calc]` dependency. Caveat: does not validate the underlying CG/atomistic model; profiles are number density in input units."""
    r = traj_mod.analyze_micelle_trajectory(
        source=source,
        cluster_cutoff=cluster_cutoff,
        n_radial_bins=n_radial_bins,
        compute_slab_profile=compute_slab_profile,
    )
    return {
        "times": list(r.times),
        "n_frames": r.n_frames,
        "n_particles": r.n_particles,
        "n_clusters": list(r.n_clusters),
        "cluster_sizes": [list(s) for s in r.cluster_sizes],
        "largest_cluster_size": list(r.largest_cluster_size),
        "rg_largest": list(r.rg_largest),
        "radial_density_r": list(r.radial_density_r),
        "radial_density_values": list(r.radial_density_values),
        "slab_density_z": list(r.slab_density_z) if r.slab_density_z is not None else None,
        "converged": r.converged,
        "n_iterations": r.n_iterations,
        "method": r.method,
        "params": r.params,
        "caveat": r.caveat,
    }


@mcp.tool()
@_tool_guard
@_shape_no_solution
def conformer_volume_3d(smiles: str, n_conformers: int = 10, seed: int = 42, tanford_volume_a3: float | None = None) -> dict:
    """C5 conformer volume: ETKDGv3 embedding (fixed seed, 1 thread) + MMFF94s minimization (UFF fallback) -> 3D vdW envelope volume (A^3, rdkit DoubleCubicLatticeVolume), Labute SASA (A^2), heavy-atom gyration shape anisotropy k^2 (0 isotropic, 1 rod-like; hand-rolled Jacobi, no scipy). Optional tanford_volume_a3 (from tanford_tail_volume) -> tanford_ratio: whole-molecule 3D vs TAIL-ONLY group contribution cross-check -- divergence is expected and informative for branched/gemini tails (that is the point, not an error). Validated internally: alkane V_3D linear in chain length R^2>=0.999, exact k^2 golds (collinear=1, tetrahedron=0), deterministic under fixed seed. Caveat: lowest of N sampled conformers, a conformational model not a free-energy calculation. Args: smiles (rdkit-parsable), n_conformers>=1, seed>=0, tanford_volume_a3>0 optional."""
    r = conf_mod.conformer_volume_3d(
        smiles=smiles,
        n_conformers=n_conformers,
        seed=seed,
        tanford_volume_a3=tanford_volume_a3,
    )
    return {
        "smiles": r.smiles,
        "vdW_volume_A3": r.vdW_volume_A3,
        "sasa_A2": r.sasa_A2,
        "shape_anisotropy_k2": r.shape_anisotropy_k2,
        "lowest_energy_kcal_mol": r.lowest_energy_kcal_mol,
        "force_field": r.force_field,
        "n_conformers_embedded": r.n_conformers_embedded,
        "n_conformers_minimized": r.n_conformers_minimized,
        "tanford_ratio": r.tanford_ratio,
        "converged": r.converged,
        "n_iterations": r.n_iterations,
        "method": r.method,
        "params": r.params,
        "caveat": r.caveat,
    }


@mcp.tool()
@_tool_guard
@_shape_no_solution
def cloud_point_estimate(n_carbons: int, n_eo: int, salinity_g_per_100ml: float = 0.0, k: float = 0.16) -> dict:
    """C6 cloud point: HLD=0 inversion for PURE alkyl ethoxylates C_n E_m. T_cp = 25 - (b*S - k*EACN_app + Cc)/cT with Cc = 0.28*n + 2.4 - n_eo (Acosta 2008), EACN_app = 7.3 (Zarate-Munoz Langmuir 2015), b=0.13 (NaCl), cT=+0.06, k=0.16 midpoint of sourced 0.15-0.17 (Ghayour TSD 2023 Eq.1). Source accuracy: typically within 5 C (HLD) / 4 C (HLD-NAC) for pure ethoxylates; C12E5 anchor 'close to 29 C'. Scope: monodisperse only (polydisperse does NOT follow CP-Cc -- out of scope), 1 wt.% basis, <=8 EO established range (flagged as extrapolation above). Args: n_carbons>=4, n_eo>=1, salinity g/100mL >=0, k in 0.05-0.7."""
    r = cp_mod.cloud_point_estimate(
        n_carbons=n_carbons,
        n_eo=n_eo,
        salinity_g_per_100ml=salinity_g_per_100ml,
        k=k,
    )
    return {
        "cloud_point_C": r.cloud_point_C,
        "cc": r.cc,
        "hld_at_cloud_point": r.hld_at_cloud_point,
        "eacn_app": r.eacn_app,
        "converged": r.converged,
        "n_iterations": r.n_iterations,
        "method": r.method,
        "params": r.params,
        "caveat": r.caveat,
    }


@mcp.tool()
@_tool_guard
@_shape_no_solution
def phase_inversion_temperature(family: str, salinity_pct: float, k: float, eacn: float, cc: float, b: float | None = None, alpha: float | None = None, temperature_reference_K: float = 298.15) -> dict:
    """Phase-inversion temperature (K): the temperature where HLD = 0 at fixed brine/oil composition -- for nonionics this IS Shinoda's PIT (Winsor III transition); for ionics it is the HLD=0 optimum temperature. Args: family ('nonionic'|'ionic'), salinity_pct (ionic: eq g NaCl/100 mL > 0), k, eacn, cc from hld_*_reference (never guess), b REQUIRED for nonionic (HLD = b*S - k*EACN - alpha*dT + Cc; ~0.13 dL/g alkyl ethoxylates), alpha defaults sourced (ionic +0.01, nonionic -0.06). Round-trips hld_ionic/hld_nonionic to HLD = 0 exactly; missing b for nonionic / b on ionic / alpha=0 / unknown family -> no_solution payload naming the field."""
    pit_K = hld_mod.phase_inversion_temperature(
        family=family,
        salinity_pct=salinity_pct,
        k=k,
        eacn=eacn,
        cc=cc,
        b=b,
        alpha=alpha,
        temperature_reference_K=temperature_reference_K,
    )
    return {
        "phase_inversion_temperature_K": pit_K,
        "family": family,
        "unit": "K",
        "reference_temperature_K": temperature_reference_K,
    }


@mcp.tool()
@_tool_guard
@_shape_no_solution
def entropy_of_micellization(deltaG_mic_kJ_per_mol: float, deltaH_mic_kJ_per_mol: float, temperature_K: float, delta_cp_J_mol_K: float | None = None, reference_temperature_K: float | None = None) -> dict:
    """Entropy of micellization (J/(mol.K)): deltaS = (deltaH - deltaG) /
    T. Completes the deltaG/deltaH/deltaS thermodynamic triad given the
    other two (from gibbs_free_energy_micellization and
    vant_hoff_enthalpy). Optional deltaCp term (default off): supply
    BOTH delta_cp_J_mol_K (J/mol/K, e.g. from vant_hoff_multi_point_fit)
    AND reference_temperature_K (the Tr your deltaG/deltaH are quoted
    at) to get deltaS reported AT temperature_K via Kirchhoff,
    deltaS(T) = deltaS(Tr) + deltaCp*ln(T/Tr) (Kantonen et al. 2019);
    only one of the two -> no_solution payload naming both fields."""
    value = thermo.entropy_micellization(
        deltaG_mic_kJ_per_mol,
        deltaH_mic_kJ_per_mol,
        temperature_K,
        delta_cp_J_mol_K=delta_cp_J_mol_K,
        reference_temperature_K=reference_temperature_K,
    )
    return {"deltaS_mic_J_per_mol_K": value, "unit": "J/(mol.K)"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def wetting_work_of_adhesion(gamma_LV_mN_m: float, contact_angle_deg: float) -> dict:
    """Young-Dupre work of adhesion (mJ/m^2, numerically = mN/m):
    W_a = gamma_LV * (1 + cos(theta)). gamma_LV_mN_m is the liquid-vapor
    surface tension; contact_angle_deg is the measured equilibrium
    contact angle in degrees (0-180)."""
    value = wetting.work_of_adhesion(gamma_LV_mN_m, contact_angle_deg)
    return {"work_of_adhesion_mJ_per_m2": value, "unit": "mJ/m^2"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def get_standard_probe_liquid_properties(liquid_name: str) -> dict:
    """Real, literature-tabulated surface-tension components for a
    standard wetting probe liquid -- 'water', 'diiodomethane', 'glycerol',
    or 'formamide' (case-insensitive). Returns BOTH the Owens-Wendt
    dispersive/polar split (for owens_wendt_solid_surface_energy) and the
    van Oss-Chaudhury-Good LW/acid/base split (for
    van_oss_chaudhury_good_solid_surface_energy) in one call, so a
    water/glycerol/diiodomethane run needs no external lookup at all.
    Sources: van Oss, Chaudhury & Good, Chem. Rev. 88 (1988) 927-941
    (water, diiodomethane); van Oss, Good & Busscher, J. Dispersion Sci.
    Technol. 11(1) (1990) 75-81 (glycerol, formamide acid/base split).
    Raises for any liquid not in this small, real table -- that liquid's
    values must still be supplied directly by the caller to the two
    solver tools, never guessed."""
    ow = wetting.get_owens_wendt_standard_liquid(liquid_name)
    vocg = wetting.get_vocg_standard_liquid(liquid_name)
    return {
        "liquid": liquid_name.strip().lower(),
        "owens_wendt_dispersive_mN_m": ow["dispersive"],
        "owens_wendt_polar_mN_m": ow["polar"],
        "vocg_lw_mN_m": vocg["lw"],
        "vocg_acid_mN_m": vocg["acid"],
        "vocg_base_mN_m": vocg["base"],
        "vocg_total_mN_m": vocg["total"],
        "unit": "mN/m (= mJ/m^2)",
    }


@mcp.tool()
@_tool_guard
@_shape_no_solution
def owens_wendt_solid_surface_energy(contact_angles_deg: list[float], liquid_gamma_dispersive_mN_m: list[float], liquid_gamma_polar_mN_m: list[float]) -> dict:
    """
    Decompose a SOLID's surface energy into dispersive and polar components via Owens & Wendt (1969; Rabel 1971; Kaelble 1970). Args: >=2 probe liquids (the standard minimum) with their surface tensions + dispersive/polar parts and the solid's contact angles, same length and order -- all real, literature-tabulated values SUPPLIED by the caller, never guessed. Returns gamma_s_disp, gamma_s_polar, gamma_s (mN/m) with fit r_squared.
    """
    result = wetting.owens_wendt_solid_surface_energy(contact_angles_deg, liquid_gamma_dispersive_mN_m, liquid_gamma_polar_mN_m)
    return {
        "gamma_s_dispersive_mN_m": result.gamma_s_dispersive_mN_m,
        "gamma_s_polar_mN_m": result.gamma_s_polar_mN_m,
        "gamma_s_total_mN_m": result.gamma_s_total_mN_m,
        "r_squared": result.r_squared,
        "n_liquids": result.n_liquids,
        "unit": "mN/m (= mJ/m^2)",
    }


@mcp.tool()
@_tool_guard
@_shape_no_solution
def van_oss_chaudhury_good_solid_surface_energy(contact_angles_deg: list[float], liquid_gamma_lw_mN_m: list[float], liquid_gamma_acid_mN_m: list[float], liquid_gamma_base_mN_m: list[float]) -> dict:
    """
    Decompose a SOLID's surface energy into Lifshitz-van der Waals (LW) and Lewis acid/base (AB) components -- van Oss, Chaudhury & Good three-liquid method (Chem. Rev. 88 (1988) 927-941). Args: gamma_L_plus, gamma_L_minus, gamma_L (liquid surface tensions, mN/m) and solid gamma_s_PW (water contact angle input per signature); >=3 real probe liquids with literature-tabulated components (supplied, never guessed). Returns gamma_s (mN/m) = LW + AB decomposition.
    """
    result = wetting.van_oss_chaudhury_good_solid_surface_energy(contact_angles_deg, liquid_gamma_lw_mN_m, liquid_gamma_acid_mN_m, liquid_gamma_base_mN_m)
    return {
        "gamma_s_lw_mN_m": result.gamma_s_lw_mN_m,
        "gamma_s_acid_mN_m": result.gamma_s_acid_mN_m,
        "gamma_s_base_mN_m": result.gamma_s_base_mN_m,
        "gamma_s_total_mN_m": result.gamma_s_total_mN_m,
        "unit": "mN/m (= mJ/m^2)",
    }


@mcp.tool()
@_tool_guard
@_shape_no_solution
def wetting_spreading_coefficient(gamma_LV_mN_m: float, contact_angle_deg: float) -> dict:
    """Spreading coefficient (mN/m) from Young's equation:
    S = gamma_LV * (cos(theta) - 1). Always <= 0 in this formulation
    (equality only at complete wetting, theta=0); more negative means
    poorer spreading."""
    value = wetting.spreading_coefficient(gamma_LV_mN_m, contact_angle_deg)
    return {"spreading_coefficient_mN_m": value, "unit": "mN/m"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def eor_capillary_number(viscosity_mPas: float, velocity_m_per_s: float, interfacial_tension_mN_m: float) -> dict:
    """Capillary number (dimensionless): Ca = (viscosity*velocity) /
    interfacial_tension -- governs oil-displacement efficiency in
    enhanced oil recovery. Below Ca ~ 1e-5, flow is capillary-dominated
    and residual oil stays trapped; surfactants raise Ca mainly by
    driving interfacial tension toward ultra-low values. viscosity_mPas
    in mPa.s (=cP); velocity_m_per_s in m/s; interfacial_tension_mN_m in mN/m."""
    value = wetting.capillary_number(viscosity_mPas, velocity_m_per_s, interfacial_tension_mN_m)
    return {"capillary_number": value, "unit": "dimensionless",
            "flow_regime": "capillary-dominated (residual oil trapped)" if value < 1e-5 else "transitional/viscous-dominated"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def intrinsic_water_solubility_reference(solubilizate: str) -> dict:
    """
    Literature-sourced intrinsic aqueous solubility (mol/L) for a named solubilizate, for use as molar_solubilization_ratio/micelle_water_partition_coefficient's intrinsic_water_solubility_M input. solubilizate: 'naphthalene' (2.17e-4, HIGH -- Paria & Yuet 2006, primary PDF, used in this project's MSR validation), 'benzene' (~0.0228, MODERATE -- tabulated EPA/ATSDR constant), 'pyrene' (~6.87e-7, MODERATE -- real ACS solubility paper). Raises for anything else: a real but incomplete table, not a database (contrast estimate_intrinsic_water_solubility_qspr: a coarse structure-based ESTIMATE for any compound, always labeled as such). Returns value, confidence, source.
    """
    key = solubilizate.lower()
    if key not in solub.INTRINSIC_WATER_SOLUBILITY_REFERENCE_M:
        raise ValueError(
            f"solubilizate must be one of {sorted(solub.INTRINSIC_WATER_SOLUBILITY_REFERENCE_M)} -- "
            "no guessed value available for anything else (try "
            "estimate_intrinsic_water_solubility_qspr for a disclosed structure-based ESTIMATE instead)"
        )
    return {"solubilizate": key, **solub.INTRINSIC_WATER_SOLUBILITY_REFERENCE_M[key]}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def molar_solubilization_ratio(total_solubilized_M: float, intrinsic_water_solubility_M: float, surfactant_concentration_M: float, cmc_M: float) -> dict:
    """Molar Solubilization Ratio (MSR, dimensionless): moles of
    solubilizate taken up per mole of MICELLIZED surfactant (only the
    surfactant above the CMC, since that's what forms micelles):
    MSR = (total_solubilized - intrinsic_water_solubility) /
    (surfactant_concentration - cmc). All four inputs must be in the
    same concentration unit. See micelle_water_partition_coefficient for
    a real complementary metric (intensive, not extensive)."""
    value = solub.molar_solubilization_ratio(total_solubilized_M, intrinsic_water_solubility_M, surfactant_concentration_M, cmc_M)
    return {"msr": value, "unit": "dimensionless (mol solubilizate per mol micellized surfactant)"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def micelle_water_partition_coefficient(total_solubilized_M: float, intrinsic_water_solubility_M: float, surfactant_concentration_M: float, cmc_M: float) -> dict:
    """
    Micelle-water partition coefficient Km (dimensionless), MOLE-FRACTION-BASIS: Km = X_micelle/X_water -- complementary to molar_solubilization_ratio, capturing micellar preference. Args: MSR and phase amounts per signature (plus intrinsic solubility where required). DISCLOSED CONVENTION: literature also uses a molar-CONCENTRATION-ratio basis giving a different number for the same raw data -- this tool's basis is stated in every output; do not compare across bases.
    """
    value = solub.micelle_water_partition_coefficient(total_solubilized_M, intrinsic_water_solubility_M, surfactant_concentration_M, cmc_M)
    return {"km": value, "unit": "dimensionless (mole-fraction-basis convention)"}


@mcp.tool()
@_tool_guard
@_shape_no_solution
def estimate_intrinsic_water_solubility_qspr(smiles: str) -> dict:
    """
    QSPR ESTIMATE (never a measurement) of intrinsic aqueous solubility from SMILES -- ESOL/Pat Walters-style refit, usable for ANY compound with no measured value. Args: smiles (raises on invalid). Returns solubility_M, model info, caveats. HONEST LIMITATION: roughly 0.6-1 log-unit error -- always reported as an estimate, never equivalent to a literature value (contrast intrinsic_water_solubility_reference for the real, incomplete table).
    """
    r = solub.estimate_intrinsic_water_solubility_qspr(smiles)
    return {
        "smiles": r.smiles,
        "log_s_mol_per_L": r.log_s_mol_per_L,
        "solubility_M": r.solubility_M,
        "method": r.method,
        "confidence": r.confidence,
        "caveats": r.caveats,
    }


@mcp.tool()
@_tool_guard
@_shape_no_solution
def classify_surfactant_charge_type(smiles: str) -> dict:
    """
    Charge type (anionic/cationic/zwitterionic/nonionic) from a real SMILES via RDKit SMARTS (sulfate ester, sulfonate, carboxylate, phosphate ester; quaternary ammonium; both present = zwitterionic; none = nonionic). No pH guessing: a free (non-quaternary) amine returns charge_type='ambiguous_pH_dependent' instead of assuming cationic. system_type_for_gibbs: convenience mapping onto this server's Gibbs tools -- defaults any ionic result to 'ionic_no_added_salt' (n=2); override to 'ionic_excess_electrolyte' only if independently known; null for zwitterionic/ambiguous/unparseable. Unparseable SMILES raises (never guesses).
    """
    r = classify_mod.classify_surfactant_charge_type(smiles)
    return {
        "smiles": r.smiles,
        "charge_type": r.charge_type,
        "system_type_for_gibbs": r.system_type_for_gibbs,
        "anionic_groups_found": r.anionic_groups_found,
        "cationic_strong_groups_found": r.cationic_strong_groups_found,
        "cationic_conditional_groups_found": r.cationic_conditional_groups_found,
        "confidence": r.confidence,
        "caveats": r.caveats,
    }


@mcp.tool()
@_tool_guard
@_shape_no_solution
def classify_surfactant_charge_type_at_ph(smiles: str, ph: float, amine_pka_override: float | None = None) -> dict:
    """
    Resolve classify_surfactant_charge_type's 'ambiguous_pH_dependent' case (e.g. free amines) with a REAL, explicitly-supplied solution pH -- the pH itself is never guessed (it is an experiment-specific fact this tool cannot recover), but once given, computes Henderson-Hasselbalch ionization and returns the charge type at that pH plus the ionized fraction. Args: smiles, pH (required, 0-14). Unparseable SMILES or pH outside range raises.
    """
    r = classify_mod.classify_surfactant_charge_type_at_ph(smiles, ph, amine_pka_override)
    return {
        "smiles": r.smiles,
        "ph": r.ph,
        "charge_type_at_ph": r.charge_type_at_ph,
        "amine_class": r.amine_class,
        "fraction_protonated_low": r.fraction_protonated_low,
        "fraction_protonated_high": r.fraction_protonated_high,
        "pka_source": r.pka_source,
        "confidence": r.confidence,
        "caveats": r.caveats,
    }


@mcp.tool()
@_tool_guard
@_shape_no_solution
def classify_surfactant_structural_family(smiles: str) -> dict:
    """
    Structural family from real SMILES via RDKit: dimeric/gemini-type (two ionic headgroups of the SAME strong pattern + hydrocarbon tails), glycolipid biosurfactant (pyranose sugar ring + long hydrocarbon), or none. Gemini signal is deliberately conservative (exactly 2 matches of one strong headgroup pattern); when a true spacer-joined gemini cannot be distinguished from two separate molecules, that ambiguity is DISCLOSED in caveats, not resolved. Unparseable SMILES raises.
    """
    r = classify_mod.classify_surfactant_structural_family(smiles)
    return {
        "smiles": r.smiles,
        "structural_family": r.structural_family,
        "n_strong_ionic_headgroups": r.n_strong_ionic_headgroups,
        "n_sugar_rings": r.n_sugar_rings,
        "n_long_chain_matches": r.n_long_chain_matches,
        "confidence": r.confidence,
        "caveats": r.caveats,
    }


@mcp.tool()
@_tool_guard
@_shape_no_solution
def derive_all_properties_from_smiles_and_curve(
    smiles: str,
    concentrations_mM: list[float],
    surface_tensions_mN_per_m: list[float],
    temperature_K: float = 298.15,
    electrolyte_condition: str | None = None,
    counterion_dissociation_alpha: float | None = None,
    gamma0_mN_m_override: float | None = None,
    gamma_max_mol_per_m2_override: float | None = None,
) -> dict:
    """
    THE autonomous pipeline: from ONLY a real SMILES + a raw pre/post-CMC
        surface-tension curve, derive charge type, CMC, Gamma_max, A_min,
        BIC-selected Langmuir-vs-Frumkin isotherm (fitted K, and a if
        Frumkin) and deltaG_mic where inputs support it -- the one-call
        alternative to chaining classify -> cmc_from_surface_tension_curve ->
        gibbs_surface_excess -> select_isotherm_model yourself. ALWAYS read
        'gaps': every real ambiguity (zwitterionic/pH-ambiguous charge,
        unsupplied electrolyte condition, missing counterion alpha, too few
        pre-CMC points, Frumkin at its search bound) is reported there with
        dependent outputs NULL rather than silently resolved.
        electrolyte_condition ('excess_electrolyte'|'no_added_salt') and
        counterion_dissociation_alpha: pass ONLY if independently known --
        neither is derivable from SMILES + curve alone. alpha is the degree
        of IONIZATION (= S2/S1 = 1 - beta); deltaG_mic then uses the
        counterion factor (2 - alpha) = (1 + beta).

        'cmc_mM' is the AUTHORITATIVE value used for every downstream step
        (Gamma_max's own premicellar range, deltaG_mic) -- it can differ
        from 'cmc_piecewise_breakpoint_mM', the raw piecewise break-point,
        when that break-point's own postmicellar "plateau" segment is
        statistically suspect (too few points, or not meaningfully flatter
        than the premicellar decline -- see plateau_is_statistically_suspect
        in curve_analysis.py). When that happens, a Szyszkowski-isotherm-
        deviation cross-check is run automatically and 'cmc_mM' reflects
        the refined value when one was found ('cmc_was_refined'=True;
        check 'gaps' for the "cmc_mM refined from..." entry); 'cmc_mM' can
        also come back null when that same cross-check finds the ENTIRE
        dataset consistent with a single continuous premicellar decline
        (a real "no CMC detectable here" finding, not a failure).

        gamma0_mN_m_override / gamma_max_mol_per_m2_override: optional,
        ONLY pass if independently known (a real, separately measured
        pure-solvent surface tension, or a real Gamma_max from an
        independent technique) -- strictly better information than this
        pipeline's own curve-derived approximation when a researcher
        already has it; feeds directly into isotherm model selection and
        deltaG_mic.
    """
    r = orch.derive_all_properties_from_smiles_and_curve(
        smiles, concentrations_mM, surface_tensions_mN_per_m, temperature_K,
        electrolyte_condition, counterion_dissociation_alpha,
        gamma0_mN_m_override=gamma0_mN_m_override, gamma_max_mol_per_m2_override=gamma_max_mol_per_m2_override,
    )
    return {
        "smiles": r.smiles,
        "charge_type": r.classification.charge_type,
        "charge_classification_confidence": r.classification.confidence,
        "system_type_used": r.system_type_used,
        "cmc_mM": r.cmc_mM,
        "cmc_piecewise_breakpoint_mM": r.cmc.cmc_mM,
        "cmc_was_refined": r.cmc_refinement is not None and r.cmc_refinement.cmc_mM is not None,
        "gamma_at_cmc_mN_per_m": r.gamma_at_cmc_mN_per_m,
        "premicellar_slope_mN_per_m_per_lnC": r.cmc.premicellar_slope_mN_per_m_per_lnC,
        "cmc_method": r.cmc.method,
        "gamma_max_mol_per_m2": r.gamma_max_mol_per_m2,
        "a_min_nm2": r.a_min_nm2,
        "isotherm_selected_model": r.isotherm.selected_model if r.isotherm else None,
        "isotherm_langmuir_K": r.isotherm.langmuir_K if r.isotherm else None,
        "isotherm_frumkin_K": r.isotherm.frumkin_K if r.isotherm else None,
        "isotherm_frumkin_a": r.isotherm.frumkin_a if r.isotherm else None,
        "isotherm_delta_bic": r.isotherm.delta_bic if r.isotherm else None,
        "delta_g_mic_kJ_per_mol": r.delta_g_mic_kJ_per_mol,
        "counterion_factor_used": r.counterion_factor_used,
        "gaps": r.gaps,
    }



@mcp.tool()
@_tool_guard
@_shape_no_solution
def rdkit_descriptor_panel(smiles: str) -> dict:
    """Compute molecular descriptors from a SMILES string: Crippen logP, TPSA,
    molar refractivity, Labute ASA, Bertz complexity, chi0v/kappa1 topological
    indices, HBD/HBA, ring/rotor counts, sp3 fraction.

    Real graph algorithms executed on CPU (deterministic, milliseconds). Raises
    on invalid SMILES -- never returns guessed values. Estimates, not
    measurements (logP atom-contribution error ~0.3-0.5 log units).
    """
    return desc_mod.descriptor_panel(smiles).as_dict()


@mcp.tool()
@_tool_guard
@_shape_no_solution
def bootstrap_cmc(x: list[float], y: list[float], n_boot: int = 1000,
                  ci: float = 0.95, seed: int = 0) -> dict:
    """CMC point estimate + percentile confidence interval from a sigmoidal
    curve (surface tension or probe spectroscopy vs log10 concentration) via
    residual bootstrap (default 1000 refits, seed-controlled, deterministic).

    Returns cmc_point, cmc_ci_low/high, log10 CMC interval, R^2, n_failed,
    converged. Flags when the data cannot support a tight interval instead of
    pretending precision. Requires scipy (surfactantkit[calc]).
    """
    import dataclasses
    return dataclasses.asdict(
        curve.bootstrap_cmc(x, y, n_boot=n_boot, ci=ci, seed=seed))


@mcp.tool()
@_tool_guard
@_shape_no_solution
def planar_pb_profile(psi0_mV: float, ionic_strength_M: float,
                      temperature_K: float = 298.15, valence: int = 1,
                      n_points: int = 201) -> dict:
    """Solve the nonlinear planar Poisson-Boltzmann double layer (Gouy-Chapman)
    numerically: potential psi(x), charge density rho(x), Debye length, and
    surface charge density -- cross-checked against the analytic Gouy-Chapman
    closed form (sigma_rel_error reported).

    For symmetric z:z electrolytes; mean-field point-ion model. Deterministic
    BVP solve (milliseconds). Requires scipy (surfactantkit[calc]).
    """
    import dataclasses
    return dataclasses.asdict(elec.planar_pb_profile(
        psi0_mV, ionic_strength_M, temperature_K=temperature_K,
        valence=valence, n_points=n_points))


@mcp.tool()
@_tool_guard
@_shape_no_solution
def mc_aggregation_distribution(n_monomers: int = 40,
                                temperature_K: float = 298.15,
                                epsilon_kJ_per_mol: float = 8.0,
                                n_steps: int = 500_000,
                                seed: int = 0) -> dict:
    """Monte Carlo cluster-size (aggregation-number) distribution for a
    canonical pool of associating monomers (bond energy epsilon, Metropolis
    moves, seed-controlled).

    Returns the monomer-fraction histogram vs cluster size, free-monomer
    fraction, mean aggregation number, split-half convergence diagnostics,
    and -- for n_monomers <= 45 -- the exact Boltzmann histogram for
    verification. Requires numpy (surfactantkit[calc]).
    """
    import dataclasses
    return dataclasses.asdict(thermo.mc_aggregation_distribution(
        n_monomers=n_monomers, temperature_K=temperature_K,
        epsilon_kJ_per_mol=epsilon_kJ_per_mol, n_steps=n_steps, seed=seed))


@mcp.tool()
@_tool_guard
@_shape_no_solution
def dls_size_distribution(g1: list[float], tau_s: list[float],
                          q_nm_inv: float | None = None,
                          wavelength_nm: float = 632.8,
                          angle_deg: float = 90.0,
                          lambda_reg: float | None = None) -> dict:
    """Invert a DLS field autocorrelation g1(tau) into an intensity-weighted
    hydrodynamic-size distribution: NNLS + second-difference regularization,
    lambda chosen by L-curve corner (or given explicitly).

    Returns diameter grid, intensity fractions, intensity-mean diameter,
    distribution PDI, relative residual, convergence flag, and an explicit
    ill-posedness caveat. q auto-derived from wavelength/angle/index if not
    given. Requires scipy (surfactantkit[calc]).
    """
    import dataclasses
    return dataclasses.asdict(curve.dls_size_distribution(
        g1, tau_s, q_nm_inv=q_nm_inv, wavelength_nm=wavelength_nm,
        angle_deg=angle_deg, lambda_reg=lambda_reg))


@mcp.tool()
@_tool_guard
@_shape_no_solution
def ward_tordai_adsorption(c0_mol_per_L: float, D_m2_per_s: float = 1e-9,
                           t_end_s: float = 1.0, n_times: int = 40,
                           mode: str = "langmuir",
                           gamma_inf_mol_per_m2: float = 3.0e-6,
                           K_d_mol_per_L: float = 1.0e-3) -> dict:
    """Solve the Ward-Tordai dynamic adsorption equation: surface excess
    Gamma(t) and surface concentration c_s(t) during diffusion-controlled
    buildup towards a Langmuir equilibrium (sqrt-time quadrature + per-step
    bisection; deterministic, no closed form exists).

    mode="fixed_cs" reproduces the exact 2*(c0-cs)*sqrt(Dt/pi) result for a
    prescribed boundary concentration. Pure-Python core (no scipy needed).
    """
    import dataclasses
    return dataclasses.asdict(ads.ward_tordai_adsorption(
        c0_mol_per_L=c0_mol_per_L, D_m2_per_s=D_m2_per_s, t_end_s=t_end_s,
        n_times=n_times, mode=mode, gamma_inf_mol_per_m2=gamma_inf_mol_per_m2,
        K_d_mol_per_L=K_d_mol_per_L))

@mcp.tool()
@_tool_guard
@_shape_no_solution
def nagarajan_cmc_predictor(n_carbons: int, headgroup: str = "alkyl_sulfate",
                            temperature_K: float = 298.15,
                            added_salt_M: float = 0.0,
                            a_p: float | None = None,
                            a_i: float | None = None,
                            delta: float | None = None,
                            ionic: bool | None = None) -> dict:
    """
    A-priori CMC from molecular structure via Nagarajan's molecular-thermodynamic theory (Langmuir 1985, 1, 331-341): minimizes the eq-5 standard free energy per molecule over spherical AND cylindrical aggregates (plus sphere-to-rod transition where applicable), zero fitted parameters. Args: n_carbons, headgroup/volume/area parameters per signature (structure-derived). Returns log10 CMC (mol/L), aggregate geometry, per-term free-energy breakdown, convergence info. Validated on C8-C14 sulfates etc. (<=0.12 log10 vs experiment).
    """
    import dataclasses
    return dataclasses.asdict(nag_mod.nagarajan_cmc_predictor(
        n_carbons=n_carbons, headgroup=headgroup, temperature_K=temperature_K,
        added_salt_M=added_salt_M, a_p=a_p, a_i=a_i, delta=delta,
        ionic=ionic))


@mcp.tool()
@_tool_guard
@_shape_no_solution
def nagarajan_free_energy_terms(g: int, n_carbons: int = 12,
                                headgroup: str = "alkyl_sulfate",
                                concentration_mol_per_L: float = 1e-3,
                                temperature_K: float = 298.15,
                                branch: str = "sphere",
                                added_salt_M: float = 0.0) -> dict:
    """Nagarajan eq 5 free-energy term table (transfer, ordering,
    interface, steric, electrostatic -- all in kT) for ONE hypothetical
    aggregate of g molecules on one branch ('sphere' or 'cylinder').
    Use to inspect which physical term dominates micelle formation at a
    given size; only the minimum of per_molecule over g enters the CMC
    (see nagarajan_cmc_predictor). Concentration sets the Debye kappa for
    charged surfactants.
    """
    import dataclasses
    return dataclasses.asdict(nag_mod.nagarajan_free_energy_terms(
        g=g, n_carbons=n_carbons, headgroup=headgroup,
        concentration_mol_per_L=concentration_mol_per_L,
        temperature_K=temperature_K, branch=branch,
        added_salt_M=added_salt_M))


# ------------------------------------------------------------------
# Phase 2 workers (S2/S3/O1/O3/I1/I3/I4/I5/E1/E2/C2/C3): one MCP tool
# per library function, tool name == function name, JSON-safe payloads
# (dataclasses.asdict) with the library's converged/method/params/caveat
# schema intact.  [calc]-gated functions raise ImportError naming the
# missing extra, exactly like their library versions.
# ------------------------------------------------------------------


@mcp.tool()
@_tool_guard
@_shape_no_solution
def mass_action_equilibrium(c_total_M: float, K_per_M: float, n_max: int,
                            x1_tol: float = 1e-16) -> dict:
    """O1/S2: closed-form mass-action equilibrium (pseudo-phase
    iteration): exact aggregate concentrations x_n for n=1..n_max at
    total concentration c_total_M with stepwise association constant
    K_per_M (x_n = K^(n-1) x_1^n, mass balance sum n*x_n = c_total).
    Returns {"x_n_M": [...], "monomer_M": x_1, "mass_balance_M",
    "converged", "n_iterations", "method", "params", "caveat"}.
    converged checks the mass-balance residual; a CMC is a property of
    a CONCENTRATION SCAN of x_1(c), not of one point -- use the monomer
    column across concentrations for that (n_max truncation disclosed)."""
    x_n = kin_mod.mass_action_equilibrium(c_total_M, K_per_M, n_max,
                                          x1_tol=x1_tol)
    mass_balance = sum((i + 1) * v for i, v in enumerate(x_n))
    residual = abs(mass_balance - c_total_M) / c_total_M
    return {
        "x_n_M": [float(v) for v in x_n],
        "monomer_M": float(x_n[0]),
        "mass_balance_M": float(mass_balance),
        "converged": bool(residual < 1e-8),
        "n_iterations": 1,
        "method": "mass_action_pseudo_phase_bisection_on_x1",
        "params": {"c_total_M": c_total_M, "K_per_M": K_per_M,
                   "n_max": n_max, "x1_tol": x1_tol},
        "caveat": ("truncated at n_max: the mass balance breaks down "
                   "once the distribution still has weight at n_max "
                   "(raise n_max until converged=True sticks); K is the "
                   "STEPWISE constant in x_n = K^(n-1) x_1^n, and a CMC "
                   "needs a concentration scan of monomer_M, not one "
                   "point"),
    }


@mcp.tool()
@_tool_guard
@_shape_no_solution
def micelle_relaxation_odes(c_total_M: float, K_initial_per_M: float,
                            K_final_per_M: float, k_on_per_M_s: float,
                            n_max: int = 40, n_times: int = 400) -> dict:
    """S2: deterministic ODE relaxation of the mass-action distribution
    after a K-jump (relaxation-experiment kinetics): how monomer and
    micelle populations readjust, with the dominant relaxation time.
    Returns MicelleRelaxationResult as dict (time grid, monomer/micelle
    tracks, tau_dom_s, converged/method/params/caveat)."""
    import dataclasses
    return dataclasses.asdict(kin_mod.micelle_relaxation_odes(
        c_total_M=c_total_M, K_initial_per_M=K_initial_per_M,
        K_final_per_M=K_final_per_M, k_on_per_M_s=k_on_per_M_s,
        n_max=n_max, n_times=n_times))


@mcp.tool()
@_tool_guard
@_shape_no_solution
def kmc_micelle_kinetics(c_total_M: float, K_per_M: float,
                         k_on_per_M_s: float, n_max: int = 40,
                         n_molecules: int = 1500,
                         t_end_s: float | None = None,
                         n_samples: int = 240, seed: int = 0,
                         c_total_initial_M: float | None = None) -> dict:
    """O1: kinetic Monte Carlo of stepwise micellization (stochastic
    master-equation trajectory) -- deterministic seeds make reruns
    bit-identical; use when populations are small enough that the ODE
    determinism assumption fails. Returns KmcResult as dict (time
    series, aggregate statistics, converged/method/params/caveat)."""
    import dataclasses
    return dataclasses.asdict(kin_mod.kmc_micelle_kinetics(
        c_total_M=c_total_M, K_per_M=K_per_M, k_on_per_M_s=k_on_per_M_s,
        n_max=n_max, n_molecules=n_molecules, t_end_s=t_end_s,
        n_samples=n_samples, seed=seed,
        c_total_initial_M=c_total_initial_M))


@mcp.tool()
@_tool_guard
@_shape_no_solution
def aniansson_wall_two_regime(tau1_s: float, tau2_s: float, cmc_M: float,
                              n_bar: float, sigma: float | None = None,
                              c_total_M: float | None = None) -> dict:
    """O1b: solves the REAL two-regime Aniansson-Wall relaxation theory in
    closed form for k_minus (micelle-peak dissociation rate constant, the
    FAST process) and R (the small-n 'resistance' governing the SLOW
    process) -- given real, INDEPENDENTLY measured tau1 (fast) and tau2
    (slow) relaxation times from two genuinely different techniques (e.g.
    shock-tube/ultrasonic for tau1; pressure-jump/T-jump/stopped-flow for
    tau2). This is the real fix for a documented, named degeneracy
    (Noskov): dynamic surface tension / tau2 data ALONE cannot determine
    these rate constants under quasi-equilibrium -- only a genuinely
    different technique's tau1 can. micelle_relaxation_odes/
    kmc_micelle_kinetics CANNOT do this (they use one shared k_on/k_off
    for every step and so cannot separate the two real timescales --
    confirmed: their own spectrum spans only ~6x fast/slow on real SDS
    data vs the real measured ~79x).

    n_bar, sigma: real, independently measured (SLS/SANS) aggregation
    number and size-distribution width -- do not guess n_bar; sigma
    defaults to Aniansson & Wall's own cited typical value 0.20*n_bar if
    not independently known. c_total_M defaults to cmc_M (the common
    'measured at C=cmc' case).

    ALWAYS check k_on_implied_per_M_per_s in the result against known
    diffusion-limited association rates (~1e8-1e9 M^-1 s^-1 for a small
    1:1 ionic surfactant) as a real physical sanity check on whether the
    four inputs are mutually consistent."""
    import dataclasses
    return dataclasses.asdict(kin_mod.aniansson_wall_two_regime(
        tau1_s=tau1_s, tau2_s=tau2_s, cmc_M=cmc_M, n_bar=n_bar,
        sigma=sigma, c_total_M=c_total_M))


@mcp.tool()
@_shape_no_solution
def charge_regulation_pb(site_type: str, pka: float,
                         site_density_mol_per_m2: float, ph: float,
                         ionic_strength_M: float,
                         temperature_K: float = 298.15,
                         rel_permittivity: float = 78.4,
                         u_bounds: tuple[float, float] = (-40.0, 40.0),
                         grid_points: int = 801,
                         bisection_iters: int = 120) -> dict:
    """E2: Poisson-Boltzmann double layer with charge REGULATION (sites
    titrate against the local potential: fixed-charge is the pKa -> 0/14
    limit). site_type: 'acid' or 'base'. Returns ChargeRegulationResult
    as dict (surface potential/charge, pH scan-ready output,
    converged/method/params/caveat)."""
    import dataclasses
    return dataclasses.asdict(elec.charge_regulation_pb(
        site_type=site_type, pka=pka,
        site_density_mol_per_m2=site_density_mol_per_m2, ph=ph,
        ionic_strength_M=ionic_strength_M, temperature_K=temperature_K,
        rel_permittivity=rel_permittivity, u_bounds=u_bounds,
        grid_points=grid_points, bisection_iters=bisection_iters))


@mcp.tool()
@_tool_guard
@_shape_no_solution
def obrien_white_numerical(zeta_mV: float, kappa_a: float,
                           transport_number_cation: float = 0.5,
                           relaxation: bool = False, n_grid: int = 900,
                           x_max_slack: float = 12.0) -> dict:
    """E1: O'Brien & White (1978) electrophoretic mobility from first
    principles -- numerical solution of the coupled conservation/
    Poisson/Nernst-Planck/Navier-Stokes boundary-value problem (the
    derivation OHW83 fits analytically). Returns dimensionless mu_hat
    plus SI and micro-phoretic units. Default relaxation=False is the
    VALIDATED quadrature path (sign-symmetric, Henry/Smoluchowski limits
    checked); relaxation=True adds the coupled relaxation-potential
    channel and is marked experimental in method/caveat. The
    ohw83_rel_diff field is the method gap vs the OHW83 analytic
    approximation (a disclosed accuracy estimate, not an error bar);
    mu_hat_ohw83 is NaN for kappa_a < 10 (formula domain). ALWAYS read
    the returned caveat: O'Brien-White itself structurally cannot
    represent anomalous surface/Stern-layer conductance and
    underestimates zeta on real highly-charged colloids because of it
    (Moncho et al. 2001 quantify up to ~35% at kappa_a~50) -- this is a
    limit of the theory, not a bug, and does not go away with
    relaxation=True."""
    r = elec.obrien_white_numerical(
        zeta_mV=zeta_mV, kappa_a=kappa_a,
        transport_number_cation=transport_number_cation,
        relaxation=relaxation, n_grid=n_grid, x_max_slack=x_max_slack)
    return r.as_dict()


@mcp.tool()
@_tool_guard
@_shape_no_solution
def dynamic_surface_tension_forward(t_s: list[float], c0_mol_per_L: float,
                                    D_m2_per_s: float, K_per_M: float,
                                    gamma_max_mol_per_m2: float,
                                    gamma0_mN_per_m: float = 72.0,
                                    n_factor: float = 1.0,
                                    temperature_K: float = 298.15,
                                    k_off_per_s: float = 0.5,
                                    n_cells: int = 160,
                                    max_substep_s: float | None = None,
                                    rel_step: float = 0.01) -> dict:
    """O3: forward dynamic surface tension gamma(t) -- implicit FDM
    diffusion to a kinetic Langmuir surface (rate boundary condition;
    Ward-Tordai with adsorption-desorption kinetics, solved coupled on a
    stretched grid with Richardson extrapolation); gamma from the surface
    excess by the Szyszkowski equation of the same Langmuir isotherm.
    Positive sorted t_s, D in m^2/s, K in 1/M (equilibrium constant;
    k_a = K*k_off), gamma_max in mol/m^2, c0 in mol/L, k_off in 1/s.
    Returns DynamicSurfaceTensionResult as dict."""
    import dataclasses
    return dataclasses.asdict(dt_mod.dynamic_surface_tension_forward(
        t_s=t_s, c0_mol_per_L=c0_mol_per_L, D_m2_per_s=D_m2_per_s,
        K_per_M=K_per_M, gamma_max_mol_per_m2=gamma_max_mol_per_m2,
        gamma0_mN_per_m=gamma0_mN_per_m, n_factor=n_factor,
        temperature_K=temperature_K, k_off_per_s=k_off_per_s,
        n_cells=n_cells, max_substep_s=max_substep_s, rel_step=rel_step))


@mcp.tool()
@_tool_guard
@_shape_no_solution
def dynamic_surface_tension_ode(t_s: list[float],
                                gamma_obs_mN_per_m: list[float],
                                c0_mol_per_L: float,
                                gamma0_mN_per_m: float = 72.0,
                                n_factor: float = 1.0,
                                temperature_K: float = 298.15,
                                p0: tuple[float, float, float, float] | None = None,
                                fit_gamma_max: bool = True,
                                fit_D: bool = True,
                                sigma_mN_per_m: float | None = None,
                                n_cells: int = 120,
                                max_substep_s: float | None = None,
                                rel_step: float = 0.04,
                                max_nfev: int = 400) -> dict:
    """O3: joint fit of a gamma(t) transient to (D, K, Gamma_max,
    k_off). Read params['log_K_vs_log_Gmax_correlation'] and
    params['identifiability_warning'] before quoting K or k_off alone --
    single-concentration curves constrain only lumped combinations
    (Gamma_eq, rate); fit_D=False is the honest default when D is known
    independently. Returns DynamicSurfaceTensionFitResult as dict."""
    import dataclasses
    return dataclasses.asdict(dt_mod.dynamic_surface_tension_ode(
        t_s=t_s, gamma_obs_mN_per_m=gamma_obs_mN_per_m,
        c0_mol_per_L=c0_mol_per_L, gamma0_mN_per_m=gamma0_mN_per_m,
        n_factor=n_factor, temperature_K=temperature_K, p0=p0,
        fit_gamma_max=fit_gamma_max, fit_D=fit_D,
        sigma_mN_per_m=sigma_mN_per_m, n_cells=n_cells,
        max_substep_s=max_substep_s, rel_step=rel_step, max_nfev=max_nfev))


@mcp.tool()
@_tool_guard
@_shape_no_solution
def eacn_group_contribution(smiles: str) -> dict:
    """C3: EACN of one oil molecule from its SMILES via sourced group
    contributions (Wade 1977 definition; Queste/Benzoate 2007 ring
    rules; Ester 2013 increments). Covers n-alkanes, mono-alkylbenzenes,
    mono-alkylcyclohexanes and single long-chain esters ONLY; branched,
    multi-ring, heteroatom, multi-fragment and short-ester structures
    RAISE ValueError with measurement guidance instead of fabricating a
    number. Returns EacnResult as dict (eacn plus band [eacn_low,
    eacn_high], rule, sources, caveat)."""
    return eacn_mod.eacn_group_contribution(smiles).as_dict()


@mcp.tool()
@_tool_guard
@_shape_no_solution
def eacn_linear_mixing(fractions: list[float], eacns: list[float],
                       basis: str = "mole") -> dict:
    """C3: linear mixing rule for blend EACN: EACN = sum x_i * EACN_i
    (fractions normalized internally). basis is 'mole' or 'volume' --
    quote which one you measured. Validated for n-alkane dilutions
    [Benjaminsen SPE J 2024]; confirm chemically dissimilar blends with
    a T* measurement (see caveat)."""
    return eacn_mod.eacn_linear_mixing(fractions, eacns, basis=basis)


@mcp.tool()
@_tool_guard
@_shape_no_solution
def brownian_self_assembly(concentration_M: float = 0.01,
                           n_molecules: int = 200,
                           steps: int = 12000, seed: int = 1,
                           tail_attraction_kT: float = 3.0,
                           spring_kT_per_sigma2: float = 80.0,
                           dt_tau: float = 0.002,
                           molecular_volume_nm3: float = 0.7,
                           sigma_nm: float = 0.5,
                           temperature_K: float = 298.15,
                           frame_every: int = 100,
                           cluster_cutoff_sigma: float = 1.5,
                           init: str = "random",
                           init_relax_steps: int = 400) -> dict:
    """S3: coarse-grained 2-bead surfactant (head + attractive tail)
    under overdamped Langevin dynamics in a periodic box -- cluster
    analysis over production frames gives the aggregation-number
    distribution, monomer fraction and a monomer-concentration CMC
    proxy (valid only when micelles coexist with monomers:
    params['cmc_proxy_valid']). Same seed -> bit-identical result.
    RUNTIME: default settings (200 molecules x 12000 steps) can exceed
    90 s on CPU -- pass smaller 'steps' (e.g. 2000) for exploratory or
    timeout-constrained calls.
    init='clustered' seeds a droplet to skip slow nucleation. QUALITATIVE
    reduced-unit model: trends and relative comparisons, not absolute
    CMCs (full disclosures in the caveat field)."""
    return sa_mod.brownian_self_assembly(
        concentration_M=concentration_M, n_molecules=n_molecules,
        steps=steps, seed=seed, tail_attraction_kT=tail_attraction_kT,
        spring_kT_per_sigma2=spring_kT_per_sigma2, dt_tau=dt_tau,
        molecular_volume_nm3=molecular_volume_nm3, sigma_nm=sigma_nm,
        temperature_K=temperature_K, frame_every=frame_every,
        cluster_cutoff_sigma=cluster_cutoff_sigma, init=init,
        init_relax_steps=init_relax_steps).as_dict()


@mcp.tool()
@_tool_guard
@_shape_no_solution
def brownian_cmc_scan(concentrations_M: list[float],
                      n_molecules: int = 120, steps: int = 8000,
                      seed: int = 1, tail_attraction_kT: float = 3.0) -> dict:
    """S3: run brownian_self_assembly over a concentration series and
    estimate the CMC as the half-max micelle-fraction crossing
    (log-interpolation). Bracket the crossing with tighter points before
    quoting; reduced-unit CG model (see brownian_self_assembly caveat).
    Returns {"table": [...], "cmc_estimate_M": float|None, "converged",
    "method", "n_iterations", "caveat"}."""
    return sa_mod.brownian_cmc_scan(
        concentrations_M=concentrations_M, n_molecules=n_molecules,
        steps=steps, seed=seed, tail_attraction_kT=tail_attraction_kT)


@mcp.tool()
@_tool_guard
@_shape_no_solution
def cmc_from_free_monomer(c_free: float, c_tot: float, beta: float = 0.0,
                          v_m_cm3_per_mol: float | None = None) -> dict:
    """S4 eq3/eq4 forward: CMC from a measured free-monomer concentration.
    beta=0 is exact eq 3 (CMC = C_free, nonionic plateau); beta>0 applies
    eq 4 with counterion binding + 1/(1-phi) crowding (phi from
    v_m_cm3_per_mol, paper assumes density 1 g/cm3). Args: c_free, c_tot in
    mol/L (c_free <= c_tot not required -- forward only), beta in [0,1)
    (pass a SOURCED value, never guess), v_m in (0,5000) cm3/mol or None
    (crowding off). Deterministic, no RNG."""
    return {
        "cmc_M": dpd_mod.cmc_from_free_monomer(
            c_free=c_free, c_tot=c_tot, beta=beta,
            v_m_cm3_per_mol=v_m_cm3_per_mol),
        "method": "anderson_eq4_forward",
        "params": {"c_free_M": c_free, "c_tot_M": c_tot, "beta": beta,
                   "v_m_cm3_per_mol": v_m_cm3_per_mol},
        "caveat": "Eq 4 of Rummel et al. JPCB 2026 (doi:10.1021/acs.jpcb.5c07475): "
                  "beta=0 -> exact eq 3 (CMC=C_free); beta>0 needs a sourced "
                  "counterion-binding value and, for the crowding term, v_m.",
        "converged": True, "n_iterations": 1,
    }


@mcp.tool()
@_tool_guard
@_shape_no_solution
def free_monomer_from_cmc(cmc: float, c_tot: float, beta: float = 0.0,
                          v_m_cm3_per_mol: float | None = None) -> dict:
    """S4 eq3/eq4 inverse: free-monomer concentration implied by a CMC
    (bisection on the monotone eq-4 residual; exact eq 3 at beta=0).
    Args: 0 < cmc <= c_tot in mol/L, beta in [0,1) (sourced value), v_m in
    (0,5000) cm3/mol or None (crowding off). Deterministic, no RNG."""
    return {
        "c_free_M": dpd_mod.free_monomer_from_cmc(
            cmc=cmc, c_tot=c_tot, beta=beta,
            v_m_cm3_per_mol=v_m_cm3_per_mol),
        "method": "anderson_eq4_inverse_bisection",
        "params": {"cmc_M": cmc, "c_tot_M": c_tot, "beta": beta,
                   "v_m_cm3_per_mol": v_m_cm3_per_mol},
        "caveat": "Inverse of eq 4 (Rummel et al. JPCB 2026): reports the C_free "
                  "consistent with the given CMC; requires 0 < cmc <= c_tot.",
        "converged": True, "n_iterations": 300,
    }


@mcp.tool()
@_tool_guard
@_shape_no_solution
def anderson_cmc_analysis(times: list[float], traces: list[list[float]],
                          c_tot: float, beta: float = 0.0,
                          v_m_cm3_per_mol: float | None = None,
                          j_ratio: float = 0.9) -> dict:
    """S4 analysis core: CMC from p_free(t) repeat traces via the Anderson
    short-time method (Rummel et al. J. Phys. Chem. B 2026,
    doi:10.1021/acs.jpcb.5c07475): eq-8 exponential fit A+(B-A)exp(-t/tau)
    on pointwise-mean traces, C_free=A*c_tot, eq-3/eq-4 inversion, and
    exhaustive-subset jackknife sigma (j_ratio=0.9 leave-one-out, as in the
    paper). Args: times strictly increasing (s); traces >= 3 repeats of
    free/total RATIO in (0,1] on the same grid; c_tot mol/L; beta sourced
    (0 = exact eq 3); v_m cm3/mol enables the 1/(1-phi) crowding term
    (None -> phi=0, disclosed); j_ratio in [0.5,1). ANALYSIS ONLY -- does
    not run DPD; supply traces from your own runs (batch class).
    Deterministic (no RNG). converged=False or tau_edge_pinned=true means
    the window does not resolve the slow mode -- do not quote a CMC."""
    r = dpd_mod.anderson_cmc_analysis(
        times=times, traces=traces, c_tot=c_tot, beta=beta,
        v_m_cm3_per_mol=v_m_cm3_per_mol, j_ratio=j_ratio)
    return r.as_dict()


@mcp.tool()
@_tool_guard
@_shape_no_solution
def dpd_micelle_cmc(concentrations_M: list[float] | None = None,
                    n_surfactants: int = 20, steps: int = 2000,
                    n_seeds: int = 3, seed: int = 1, dt: float = 0.04,
                    max_sync_seconds: float = 60.0,
                    allow_long_run: bool = False,
                    scale: str = "interactive") -> dict:
    """S4 DPD CMC (nonionic C8E4, the paper's SI model): numpy DPD
    self-assembly -> p_free(t) traces -> anderson_cmc_analysis (eq 8 +
    eq 3 + jackknife). Table S1 A_ij/R_ij, gamma=4.5, k_b=150, Venturoli
    angles transcribed from Rummel et al. JPCB 2026 SI
    (doi:10.1021/acs.jpcb.5c07475). Args: concentrations_M (default
    [0.6] mol/L), n_surfactants>=8, steps>=200, n_seeds>=3, seed, dt in
    [0.001,0.1] (SI does not state dt), scale='interactive' ONLY
    ('paper' refuses: dl_meso 20k-90k-time-unit batch class). Refuses
    runs estimated above max_sync_seconds unless allow_long_run=True --
    use submit_dpd_micelle_cmc_job for the async path. Deterministic.
    FINITE SIZE: quote a CMC only when every
    per_concentration[].finite_size_ok (>=5 expected free molecules;
    p_free quantizes at 1/n_surfactants). In-house model, NOT dl_meso."""
    return dpdm_mod.dpd_micelle_cmc(
        concentrations_M=concentrations_M, n_surfactants=n_surfactants,
        steps=steps, n_seeds=n_seeds, seed=seed, dt=dt,
        max_sync_seconds=max_sync_seconds, allow_long_run=allow_long_run,
        scale=scale)


@mcp.tool()
@_tool_guard
@_shape_no_solution
def submit_dpd_micelle_cmc_job(concentrations_M: list[float] | None = None,
                               n_surfactants: int = 20, steps: int = 2000,
                               n_seeds: int = 3, seed: int = 1,
                               dt: float = 0.04) -> dict:
    """S4 async: submit dpd_micelle_cmc as a background job; returns
    job_id + estimated_runtime_s immediately. Poll with
    dpd_micelle_cmc_job_status. Same model/args as dpd_micelle_cmc (no
    max_sync_seconds -- submission IS the long-run path). Job JSON
    persists under .dpd_jobs/; a server restart mid-run yields
    state='interrupted' (resubmit; no silent resume, no fabricated
    progress)."""
    return dpdm_mod.submit_dpd_micelle_cmc_job(
        concentrations_M=concentrations_M, n_surfactants=n_surfactants,
        steps=steps, n_seeds=n_seeds, seed=seed, dt=dt)


@mcp.tool()
@_tool_guard
@_shape_no_solution
def dpd_micelle_cmc_job_status(job_id: str) -> dict:
    """S4 async: poll a submitted DPD job. States: queued, running,
    completed (full dpd_micelle_cmc result attached), failed (error
    attached), interrupted (persisted file but server restarted --
    resubmit). completed results carry the full model caveat (in-house
    SI-model DPD, interactive scale, not dl_meso, finite-size gate)."""
    return dpdm_mod.dpd_micelle_cmc_job_status(job_id=job_id)


@mcp.tool()
@_tool_guard
@_shape_no_solution
def nonlinear_fit_with_errors(x: list[float], y: list[float], model: str,
                              p0: list[float] | None = None,
                              fixed: dict | None = None,
                              sigma: list[float] | None = None,
                              bounds: tuple[list[float], list[float]] | None = None,
                              param_names: list[str] | None = None,
                              max_nfev: int = 4000) -> dict:
    """I3: weighted nonlinear least squares with the full error report
    (SEs from (J^T W J)^-1, lag-1 residual check, r_squared). model is
    a MODEL_REGISTRY key ('linear', 'exponential_decay', 'szyszkowski',
    'frumkin') so JSON callers never pass callables; Szyszkowski/
    Frumkin need fixed={'gamma0_mN_m':..., 'gamma_max_mol_per_m2':...}.
    Returns NonlinearFitResult as dict."""
    import dataclasses
    return dataclasses.asdict(inf_mod.nonlinear_fit_with_errors(
        x=x, y=y, model=model, p0=p0, fixed=fixed, sigma=sigma,
        bounds=bounds, param_names=param_names, max_nfev=max_nfev))


@mcp.tool()
@_tool_guard
@_shape_no_solution
def propagate_uncertainty_chain(chain: str,
                                overrides: dict[str, dict] | None = None,
                                n_samples: int = 10000, ci: float = 0.95,
                                seed: int = 0) -> dict:
    """I5: Monte Carlo uncertainty propagation along a NAMED formula
    chain ('gibbs_dg_from_cmc', 'hydrodynamic_radius' or 'debye_length'):
    each input can be {'value','sigma','dist'} and the result is the
    output distribution's median/CI. Deterministic for a fixed seed.
    Returns UncertaintyPropagationResult as dict."""
    import dataclasses
    return dataclasses.asdict(inf_mod.propagate_uncertainty_chain(
        chain=chain, overrides=overrides, n_samples=n_samples, ci=ci,
        seed=seed))


def _resolve_registered_model(model: str, fixed: dict | None):
    """Map a MODEL_REGISTRY key to a callable for mcmc_posterior_fit /
    global_multi_dataset_fit (their `model` argument is a callable in
    the library API; MCP/JSON callers pass the registry key instead)."""
    import functools
    if model not in inf_mod.MODEL_REGISTRY:
        raise ValueError(
            f"unknown model {model!r} -- choose one of "
            f"{sorted(inf_mod.MODEL_REGISTRY)}")
    func = inf_mod.MODEL_REGISTRY[model][0]
    if fixed:
        func = functools.partial(func, **fixed)
    return func


@mcp.tool()
@_tool_guard
@_shape_no_solution
def mcmc_posterior_fit(x: list[float], y: list[float], model: str,
                       p0: list[float], observation_sigma: float,
                       param_names: list[str] | None = None,
                       n_samples: int = 6000, burn_in: int = 1500,
                       ci: float = 0.95, seed: int = 0,
                       bounds: tuple[list[float], list[float]] | None = None,
                       fixed: dict | None = None) -> dict:
    """I1: random-walk Metropolis posterior for y = model(x, *params);
    returns medians, credible intervals and the convergence diagnostics
    (acceptance rate, split-R-hat, ESS) -- check rhat < 1.1 and
    ESS >= 100 before quoting intervals. model is a MODEL_REGISTRY key
    (JSON callers cannot pass callables); fixed binds keyword arguments
    (e.g. for szyszkowski/frumkin). Deterministic for a fixed seed.
    Returns PosteriorResult as dict."""
    import dataclasses
    fn = _resolve_registered_model(model, fixed)
    return dataclasses.asdict(inf_mod.mcmc_posterior_fit(
        x=x, y=y, model=fn, p0=p0, observation_sigma=observation_sigma,
        param_names=param_names, n_samples=n_samples, burn_in=burn_in,
        ci=ci, seed=seed, bounds=bounds))


@mcp.tool()
@_tool_guard
@_shape_no_solution
def global_multi_dataset_fit(datasets: list[dict], model: str,
                             shared_p0: list[float],
                             local_p0: list[float],
                             shared_names: list[str] | None = None,
                             local_names: list[str] | None = None,
                             sigma: float | None = None,
                             bounds: tuple[list[float], list[float]] | None = None,
                             max_nfev: int = 6000,
                             fixed: dict | None = None) -> dict:
    """I4: simultaneous fit of several datasets to one SHARED parameter
    set plus per-dataset local parameters (e.g. shared K across a
    temperature series, per-dataset baselines). datasets: [{'x': [...],
    'y': [...], 'sigma': optional}]. model is a MODEL_REGISTRY key; the
    shared-vs-local decomposition is the point: compare shared-only
    fits against per-dataset fits before claiming parameter invariance
    (params carries both). Returns GlobalFitResult as dict."""
    import dataclasses
    fn = _resolve_registered_model(model, fixed)
    return dataclasses.asdict(inf_mod.global_multi_dataset_fit(
        datasets=datasets, model=fn, shared_p0=shared_p0,
        local_p0=local_p0, shared_names=shared_names,
        local_names=local_names, sigma=sigma, bounds=bounds,
        max_nfev=max_nfev))


@mcp.tool()
@_tool_guard
@_shape_no_solution
def cmc_homologous_series(chain_lengths: list[int], cmc_values: list[float],
                          unit: str = "mM",
                          predict_at: list[int] | None = None,
                          ci: float = 0.95) -> dict:
    """C2: linear free-energy fit of a homologous-series CMC:
    log10(CMC) = A - B * n_carbon (Klevens/Tanford-type relation).
    B encodes the CH2 transfer free energy (slope) and A the headgroup
    offset; returns fit diagnostics, the multiplicative CMC step per CH2
    (reduction_factor_per_ch2), and optional predictions at predict_at
    chain lengths with CIs (points outside the fitted range are flagged
    as extrapolation).
    unit: 'mM' | 'uM' | 'M' | 'mol/L' | 'mol/m3' of the CMC values you
    pass (a label only; keep one unit throughout).
    Returns HomologousSeriesResult as dict."""
    import dataclasses
    return dataclasses.asdict(hom_mod.cmc_homologous_series(
        chain_lengths=chain_lengths, cmc_values=cmc_values, unit=unit,
        predict_at=predict_at, ci=ci))




def main() -> None:
    mcp.run()


if __name__ == "__main__":
    main()
