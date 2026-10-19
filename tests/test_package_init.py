"""Tests for surfactantkit's own top-level __init__.py export surface.

Added 2026-09-15 after a real gap was found during a Tier 1 bottleneck-
resolution completeness audit: __init__.py had not been updated since
before classify.py existed, so a direct `from surfactantkit import X`
user could not reach classify_surfactant_charge_type,
classify_surfactant_structural_family, or any of the 5 new functions
built this session -- despite all of them being real, tested, and
wired into the MCP server. A first manual pass fixed the functions this
session had touched; a second, PROGRAMMATIC pass (AST-walking every
module file for public top-level functions and diffing against
__all__) caught two more real, pre-existing gaps the manual pass missed
entirely (adsorption.py's frumkin_fit_K_and_a and select_isotherm_model,
both real, tested, MCP-wired capabilities from 2026-09-13). That
programmatic check is now test_no_public_function_is_missing_from_dunder_all
below, permanently, per this project's own "close a workflow to its
entirety, and guard against it recurring" standing rule (see
~/.claude/CLAUDE.md) -- a future new public function added to any module
and never added to __all__ now fails loudly instead of silently.
"""

import ast
import os

import surfactantkit


def test_all_names_in_dunder_all_are_real_attributes():
    for name in surfactantkit.__all__:
        assert hasattr(surfactantkit, name), f"'{name}' is listed in __all__ but not actually importable"


def test_dunder_all_has_no_duplicates():
    assert len(surfactantkit.__all__) == len(set(surfactantkit.__all__))


def test_all_names_in_dunder_all_are_callable():
    for name in surfactantkit.__all__:
        assert callable(getattr(surfactantkit, name)), f"'{name}' in __all__ is not callable"


def test_tier0_tier1_bottleneck_fix_functions_are_top_level_importable():
    """Direct regression test for the specific gap found 2026-09-15 --
    every function built during the Tier 0/Tier 1 bottleneck-resolution
    pass must be reachable via `from surfactantkit import X`, not just
    via its own submodule or the MCP server."""
    from surfactantkit import (
        classify_surfactant_charge_type_at_ph,
        estimate_axial_ratio_from_cpp_geometry,
        estimate_intrinsic_water_solubility_qspr,
        estimate_partial_specific_volume_from_tail_and_headgroup,
        get_owens_wendt_standard_liquid,
        get_vocg_standard_liquid,
    )
    assert all(
        callable(f)
        for f in (
            get_vocg_standard_liquid,
            get_owens_wendt_standard_liquid,
            estimate_axial_ratio_from_cpp_geometry,
            estimate_partial_specific_volume_from_tail_and_headgroup,
            classify_surfactant_charge_type_at_ph,
            estimate_intrinsic_water_solubility_qspr,
        )
    )


def test_pre_existing_classify_and_orchestrator_functions_are_top_level_importable():
    """Regression test for the second, wider gap found in the same
    2026-09-15 completeness audit: classify.py's own two pre-existing
    functions, and mixture_model_selection.py/orchestrate.py's entire
    public API, were ALSO never exported from __init__.py -- not new
    this session, but real and now fixed in the same pass."""
    from surfactantkit import (
        classify_surfactant_charge_type,
        classify_surfactant_structural_family,
        derive_all_properties_from_smiles_and_curve,
        rubingh_predict_cmc_mix,
        select_mixture_model,
    )
    assert all(
        callable(f)
        for f in (
            classify_surfactant_charge_type,
            classify_surfactant_structural_family,
            rubingh_predict_cmc_mix,
            select_mixture_model,
            derive_all_properties_from_smiles_and_curve,
        )
    )


def test_no_public_function_is_missing_from_dunder_all():
    """Real, permanent guard against the exact class of bug found
    2026-09-15: AST-walks every module source file for top-level public
    (non-underscore) function definitions and asserts every one of them
    is listed in surfactantkit.__all__. Deliberately excludes
    mcp_server.py (its @mcp.tool()-wrapped functions are a separate,
    intentionally different public surface, exercised by
    test_mcp_server.py's own "expected tools" check) and __init__.py
    itself. Any module deliberately choosing NOT to export a given
    function must add it to _DELIBERATELY_UNEXPORTED below with a real
    reason, not just fall through silently."""
    _DELIBERATELY_UNEXPORTED: dict[str, set[str]] = {
        # Audit P2-9 consolidation (2026-10-02): these are the shared
        # implementations behind the historical private aliases
        # (_bic / _golden_section_minimize / _solve_3x3) that the domain
        # modules import.  They stay private on purpose -- exporting them
        # at package top level would turn an internal numerics helper
        # into public API (and an MCP-parity candidate) it was never
        # meant to be.  Consumers use the domain modules' _-aliases.
        "_optimize": {"golden_section_minimize", "solve_3x3"},
        "_stats": {"bic"},
    }

    src_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "src", "surfactantkit")
    skip_files = {"__init__.py", "mcp_server.py"}
    exported = set(surfactantkit.__all__)

    missing: dict[str, list[str]] = {}
    for fname in sorted(os.listdir(src_dir)):
        if not fname.endswith(".py") or fname in skip_files:
            continue
        module_name = fname[:-3]
        path = os.path.join(src_dir, fname)
        with open(path, encoding="utf-8") as f:
            tree = ast.parse(f.read())
        public_funcs = [
            node.name for node in tree.body
            if isinstance(node, ast.FunctionDef) and not node.name.startswith("_")
        ]
        allowed = _DELIBERATELY_UNEXPORTED.get(module_name, set())
        unexported = [f for f in public_funcs if f not in exported and f not in allowed]
        if unexported:
            missing[module_name] = unexported

    assert missing == {}, (
        f"Public functions defined but not exported in __all__ (fix in __init__.py, "
        f"or add to _DELIBERATELY_UNEXPORTED above with a real reason): {missing}"
    )


def test_version_matches_current_capability_state():
    # bumped 0.3.0 -> 0.4.0 alongside this session's Tier 0/1 additions,
    # matching this project's own established bump-on-feature-addition
    # convention (0.1.0 -> 0.2.0 -> 0.3.0 in prior commits)
    assert surfactantkit.__version__ == "0.4.0"


def test_pyproject_version_matches_init_version():
    """Real, permanent guard against a second real gap found in the same
    2026-09-15 completeness audit: pyproject.toml's own [project].version
    had said "0.1.0" since this repo's very first commit and was NEVER
    bumped alongside __init__.py.__version__ across 0.2.0/0.3.0/0.4.0 --
    meaning `pip show surfactantkit` (or anything reading package
    metadata rather than importing the module) would have reported a
    version 3 real releases behind the actual code. Fixed in the same
    pass; this test keeps the two from drifting apart again."""
    import re

    # Deliberately not using tomllib (Python 3.11+ only; this project
    # supports >=3.10 per pyproject.toml itself) or adding a tomli
    # dependency just for this one check -- a plain regex on the single
    # `version = "..."` line under [project] is enough and keeps this
    # project's real zero-added-dependency design intact.
    pyproject_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "pyproject.toml")
    with open(pyproject_path, encoding="utf-8") as f:
        text = f.read()
    match = re.search(r'^version\s*=\s*"([^"]+)"', text, re.MULTILINE)
    assert match, "could not find a version = \"...\" line in pyproject.toml"
    pyproject_version = match.group(1)
    assert pyproject_version == surfactantkit.__version__, (
        f"pyproject.toml version ({pyproject_version}) and surfactantkit.__version__ "
        f"({surfactantkit.__version__}) have drifted apart -- bump both together"
    )
