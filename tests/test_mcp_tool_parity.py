"""MCP tool-name <-> library parity guard (audit TOOLKIT_AUDIT P1-7;
blueprint rule: MCP tool name == function name).

Every ``@mcp.tool`` name must either match a public library function
exactly, or be listed in ``TOOL_LIB_ALIASES`` -- the documented mapping
for the 46 legacy Phase-0 wrappers whose MCP names predate the rule
(alias table: tool -> primary library attribute it calls). Both sides
are re-verified against the source on every run, so

* a NEW tool with a drifted name fails until it matches the library
  name (preferred) or is consciously added to the alias table,
* a renamed/deleted library attribute fails the alias-target check,
* a stale alias entry fails the set-equality check.

Pure static analysis -- no MCP runtime needed, deterministic, fast.
"""
import ast
from pathlib import Path

_SRC = Path(__file__).resolve().parent.parent / "src" / "surfactantkit"

# tool name -> "module.attr" (module = file stem, attr = public function
# or module-level constant the wrapper actually reads). 46 entries: the
# legacy wrappers shipped before the parity rule; new workers must match
# the library name instead of extending this table.
TOOL_LIB_ALIASES: dict[str, str] = {
    "aggregation_number": "cpp.aggregation_number_spherical",
    "aggregation_number_from_quenching": "curve_analysis.aggregation_number_from_quenching_curve",
    "aggregation_number_from_sls": "curve_analysis.aggregation_number_from_sls_debye_plot",
    "aggregation_number_from_svedberg": "curve_analysis.aggregation_number_from_svedberg_equation",
    "blankschtein_counterion_self_energy": "cpp.blankschtein_counterion_self_energy_release",
    "blankschtein_entropy_of_binding": "cpp.blankschtein_entropy_of_binding_free_energy",
    "blankschtein_steric_free_energy": "cpp.blankschtein_steric_free_energy_with_counterion",
    "cmc_from_conductivity": "curve_analysis.cmc_from_conductivity_curve",
    "corrin_harkins_predict": "mixed_micelle.corrin_harkins_predict_cmc",
    "counterion_binding_degree_reference": "thermodynamics.COUNTERION_BINDING_DEGREE_REFERENCE",
    "debye_screening_length": "electrostatics.debye_length",
    "derive_davies_group_number": "hlb.derive_davies_group_number_from_griffin",
    "dn_dc_reference": "curve_analysis.DN_DC_REFERENCE_ML_PER_G",
    "entropy_of_micellization": "thermodynamics.entropy_micellization",
    "eor_capillary_number": "wetting.capillary_number",
    "estimate_partial_specific_volume": "curve_analysis.estimate_partial_specific_volume_from_tail_and_headgroup",
    "frumkin_predict_surface_tension": "adsorption.frumkin_surface_tension",
    "gemini_hlb_griffin_reference": "hlb.GEMINI_SURFACTANT_HLB_GRIFFIN_REFERENCE",
    "get_standard_probe_liquid_properties": "wetting.get_owens_wendt_standard_liquid",
    "gibbs_area_per_molecule": "adsorption.gibbs_a_min",
    "gibbs_surface_excess": "adsorption.gibbs_gamma_max",
    "hlb_from_groups": "hlb.hlb_davies",
    "hlb_from_mw": "hlb.hlb_griffin",
    "hlb_method_recommendation": "hlb.recommend_hlb_method_for_structural_family",
    "hld_cationic_quat_reference": "hld.hlb_from_cc_cationic_quat",
    "hld_characteristic_curvature_reference": "hld.CC_REFERENCE_ANIONIC_CATIONIC",
    "hld_class_reference": "hld.CC_REFERENCE_ANIONIC_CATIONIC",
    "hld_fit_k_and_cc_from_salinity_scan": "hld.fit_k_and_cc_from_salinity_scan",
    "hld_optimal_salinity": "hld.optimal_salinity_ionic",
    "hydrodynamic_radius": "dynamics.hydrodynamic_radius_stokes_einstein",
    "intrinsic_water_solubility_reference": "solubilization.INTRINSIC_WATER_SOLUBILITY_REFERENCE_M",
    "nagarajan_ruckenstein_dipole_free_energy": "cpp.nagarajan_ruckenstein_dipole_headgroup_free_energy",
    "nagarajan_ruckenstein_headgroup_reference": "cpp.NAGARAJAN_RUCKENSTEIN_HEADGROUP_CONSTANTS",
    "nagarajan_ruckenstein_ionic_free_energy": "cpp.nagarajan_ruckenstein_ionic_headgroup_free_energy",
    "nagarajan_ruckenstein_steric_free_energy": "cpp.nagarajan_ruckenstein_steric_headgroup_free_energy",
    "predict_mobility_relaxation_corrected": "electrostatics.electrophoretic_mobility_relaxation_corrected",
    "rdkit_descriptor_panel": "descriptors.descriptor_panel",
    "rodenas_x1_from_series": "mixed_micelle.rodenas_x1_series",
    "rosen_monolayer_solve": "mixed_micelle.solve_rosen_monolayer_x",
    "rubingh_activity_coefficients": "mixed_micelle.activity_coefficients",
    "rubingh_solve": "mixed_micelle.solve_rubingh_x",
    "szyszkowski_predict_surface_tension": "adsorption.szyszkowski_surface_tension",
    "tanford_chain_geometry": "cpp.tanford_tail_volume",
    "wetting_spreading_coefficient": "wetting.spreading_coefficient",
    "wetting_work_of_adhesion": "wetting.work_of_adhesion",
    "zeta_potential": "electrostatics.zeta_potential_henry"
}


def _parse(path: Path) -> ast.Module:
    return ast.parse(path.read_text(encoding="utf-8"))


def _tool_names() -> set[str]:
    tree = _parse(_SRC / "mcp_server.py")
    names = set()
    for node in tree.body:
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        for dec in node.decorator_list:
            if (isinstance(dec, ast.Call) and isinstance(dec.func, ast.Attribute)
                    and isinstance(dec.func.value, ast.Name)
                    and dec.func.value.id == "mcp"):
                names.add(node.name)
                break
    return names


def _module_attrs(tree: ast.Module) -> set[str]:
    """Public module-level attributes: functions, classes, constants --
    including ones nested in module-level try/if blocks, excluding anything
    nested inside another function."""
    attrs: set[str] = set()

    parent: dict[ast.AST, ast.AST] = {}
    for node in ast.walk(tree):
        for child in ast.iter_child_nodes(node):
            parent[child] = node

    def inside_def(node: ast.AST) -> bool:
        # module attrs only: exclude anything nested in a function OR class
        # (class members are not module-level attributes)
        cur = parent.get(node)
        while cur is not None:
            if isinstance(cur, (ast.FunctionDef, ast.AsyncFunctionDef,
                                ast.Lambda, ast.ClassDef)):
                return True
            cur = parent.get(cur)
        return False

    for node in ast.walk(tree):
        if (isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))
                and not node.name.startswith("_") and not inside_def(node)):
            attrs.add(node.name)
        elif isinstance(node, ast.Assign) and not inside_def(node):
            for tgt in node.targets:
                if isinstance(tgt, ast.Name) and not tgt.id.startswith("_"):
                    attrs.add(tgt.id)
        elif (isinstance(node, ast.AnnAssign) and not inside_def(node)
                and isinstance(node.target, ast.Name)
                and not node.target.id.startswith("_")):
            attrs.add(node.target.id)
    return attrs


def _library_attrs() -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    for path in sorted(_SRC.glob("*.py")):
        if path.name in ("mcp_server.py", "__init__.py"):
            continue
        out[path.stem] = _module_attrs(_parse(path))
    return out


def test_every_tool_matches_library_or_documented_alias() -> None:
    tools = _tool_names()
    assert len(tools) >= 100, f"tool discovery broke: only found {len(tools)}"
    lib_fns = set()
    for attrs in _library_attrs().values():
        lib_fns |= attrs
    unmatched = {t for t in tools if t not in lib_fns}
    # exact set equality both ways: no stale aliases, no undocumented drift
    assert unmatched == set(TOOL_LIB_ALIASES), (
        "MCP tool names must match a public library function (preferred) or "
        f"be listed in TOOL_LIB_ALIASES. undocumented={sorted(unmatched - set(TOOL_LIB_ALIASES))} "
        f"stale_aliases={sorted(set(TOOL_LIB_ALIASES) - unmatched)}"
    )


def test_alias_targets_exist() -> None:
    lib = _library_attrs()
    for tool, target in TOOL_LIB_ALIASES.items():
        assert "." in target, f"{tool}: alias target must be 'module.attr', got {target!r}"
        mod, _, attr = target.partition(".")
        assert mod in lib, f"{tool}: no library module {mod!r}"
        assert attr in lib[mod], f"{tool}: {mod}.{attr} no longer exists"


def test_alias_keys_are_tools() -> None:
    tools = _tool_names()
    for tool in TOOL_LIB_ALIASES:
        assert tool in tools, f"alias key {tool!r} is not an MCP tool anymore"
