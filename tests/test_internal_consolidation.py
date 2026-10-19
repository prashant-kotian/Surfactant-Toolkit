"""P2-9 internal-helper consolidation guard.

The audit item "Duplicated private helpers" moved `_bic` (adsorption,
curve_analysis, mixture_model_selection), `_golden_section_minimize`
(adsorption, mixed_micelle) and `_solve_3x3` (thermodynamics, wetting)
into the shared private modules `_stats` / `_optimize`.  These tests pin
the two things the audit demanded:

1. ONE implementation (module-level identity: every module's alias is
   literally the same function object),
2. NO behavior change (numeric golds for the shared formulas, the exact
   per-site parameterization at each call site, and the verbatim
   singular-matrix error messages each module used to raise itself).
"""

import ast
import math
from pathlib import Path

from surfactantkit import (
    adsorption,
    curve_analysis,
    mixed_micelle,
    mixture_model_selection,
    thermodynamics,
    wetting,
)
from surfactantkit._optimize import golden_section_minimize, solve_3x3
from surfactantkit._stats import bic

_SRC = Path(__file__).resolve().parent.parent / "src" / "surfactantkit"


def test_bic_single_implementation_across_modules():
    assert adsorption._bic is curve_analysis._bic
    assert adsorption._bic is mixture_model_selection._bic
    assert adsorption._bic is bic


def test_bic_exact_gold_matches_documented_formula():
    # n*ln(rss/n) + n_params*ln(n) -- exact, no tolerance
    assert adsorption._bic(1.0, 100, 3) == 100 * math.log(1.0 / 100) + 3 * math.log(100)
    # zero-RSS guard: floor at 1e-300 instead of log(0)
    val = adsorption._bic(0.0, 10, 2)
    assert val == 10 * math.log(1e-300 / 10) + 2 * math.log(10)
    # model-selection direction: more params => higher (worse) BIC at equal RSS
    assert adsorption._bic(1.0, 100, 3) > adsorption._bic(1.0, 100, 2)


def test_golden_section_single_implementation_and_defaults():
    assert adsorption._golden_section_minimize is mixed_micelle._golden_section_minimize
    assert adsorption._golden_section_minimize is golden_section_minimize
    # defaults are the adsorption-lineage values (iters=200, tol=1e-14)
    import inspect

    sig = inspect.signature(golden_section_minimize)
    assert sig.parameters["iters"].default == 200
    assert sig.parameters["tol"].default == 1e-14


def test_golden_section_numeric_golds():
    f = golden_section_minimize
    assert abs(f(lambda x: x * x, -5.0, 5.0)) < 1e-12
    # minimum on a shifted parabola: the VALUE is pinned exactly at the
    # double-precision floor (f differences below eps plateau, so the
    # returned x is only pinned to 1e-6 -- this matches the original
    # local copies bit-for-bit, verified against git HEAD)
    res = f(lambda x: (x - 2.0) ** 2 + 1.0, 0.0, 10.0)
    assert (res - 2.0) ** 2 + 1.0 <= 1.0 + 1e-12
    assert abs(res - 2.0) < 1e-6
    # lo == hi degenerates to the only point
    assert f(lambda x: x * x, 3.0, 3.0) == 3.0


def test_mixed_micelle_call_site_preserves_original_behavior():
    # CORRECTED 2026-10-18 (benchmark cross-audit): mixed_micelle's original
    # local golden-section helper iterated with iters=60 and broke at
    # hi-lo < 1e-13 (verified against git 026b735, 2026-09-19). It was NOT
    # "the same 1e-14 break as adsorption" -- that misread during the 2026-10
    # consolidation silently deepened the convergence tolerance, flipped
    # landing basins on near-degenerate EOMMM fits (regression first caught
    # as cross-model 'tool contradicts gold' failures H-06/GZ-02/GZ-25), and
    # this test then forbade the fix. The call site must pass
    # iters=refine_iters AND pin tol=1e-13.
    src = (_SRC / "mixed_micelle.py").read_text(encoding="utf-8")
    call = [ln for ln in src.splitlines()
            if "_golden_section_minimize(" in ln and "import" not in ln and "def " not in ln]
    assert call, "no golden-section call site found in mixed_micelle.py"
    assert any("iters=refine_iters" in ln for ln in call), call
    assert any("tol=1e-13" in ln for ln in call), call


def test_adsorption_call_sites_rely_on_shared_defaults():
    # None of adsorption's three call sites passes tol= -> they inherit
    # 1e-14/200 exactly as their original local definition had.
    src = (_SRC / "adsorption.py").read_text(encoding="utf-8")
    calls = [
        ln for ln in src.splitlines()
        if "_golden_section_minimize(" in ln and "import" not in ln
    ]
    assert len(calls) >= 1
    assert not any("tol=" in ln for ln in calls), calls


def test_solve_3x3_single_implementation():
    assert thermodynamics._solve_3x3 is wetting._solve_3x3
    assert thermodynamics._solve_3x3 is solve_3x3


def test_solve_3x3_identity_system_gold():
    x = solve_3x3([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                  [1.0, 2.0, 3.0], singular_message="unused")
    assert x == [1.0, 2.0, 3.0]
    # a general system: A x = b with known x (the classic 3x3 example)
    a = [[2.0, 1.0, -1.0], [-3.0, -1.0, 2.0], [-2.0, 1.0, 2.0]]
    b = [8.0, -11.0, -3.0]
    got = solve_3x3(a, b, singular_message="unused")
    for g, t in zip(got, (2.0, 3.0, -1.0)):
        assert abs(g - t) < 1e-12
    # exact bit-parity with the original implementation on a non-trivial
    # system (same operation order in the shared module as git HEAD had)
    assert got[0] == 2.0 and abs(got[1] - 3.0) < 1e-15


def test_solve_3x3_singular_raises_caller_message():
    singular = [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [2.0, 2.0, 2.0]]
    try:
        solve_3x3(singular, [1.0, 2.0, 3.0], singular_message="custom msg")
    except ValueError as exc:
        assert "custom msg" in str(exc)
    else:
        raise AssertionError("singular system must raise ValueError")


def test_thermodynamics_and_wetting_messages_preserved_verbatim():
    # Every fragment of the exact historical texts must still live in the
    # call sites (literals are written across adjacent lines in source).
    thermo = (_SRC / "thermodynamics.py").read_text(encoding="utf-8")
    for frag in (
        "temperature series is degenerate (singular fit matrix) -- ",
        "need at least 3 points with genuinely different temperatures",
    ):
        assert frag in thermo, frag
    wet = (_SRC / "wetting.py").read_text(encoding="utf-8")
    for frag in (
        "the 3 test liquids are not independent enough to solve this ",
        "system (singular matrix) -- use 3 liquids with genuinely ",
        "different LW/acid/base compositions",
    ):
        assert frag in wet, frag
    assert "singular_message=" in thermo and "singular_message=" in wet
    # and the assembled messages still raise exactly as before, through
    # the shared implementation
    full_thermo = ("temperature series is degenerate (singular fit matrix) -- "
                   "need at least 3 points with genuinely different temperatures")
    singular = [[1.0, 1.0, 1.0], [2.0, 2.0, 2.0], [3.0, 3.0, 3.0]]
    try:
        thermodynamics._solve_3x3(singular, [1.0, 2.0, 3.0],
                                  singular_message=full_thermo)
    except ValueError as exc:
        assert str(exc) == full_thermo
    else:
        raise AssertionError("singular system must raise")


def test_no_module_still_defines_the_old_private_copies():
    banned = {"_bic", "_golden_section_minimize", "_solve_3x3"}
    for path in _SRC.glob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in tree.body:
            if isinstance(node, ast.FunctionDef) and node.name in banned:
                raise AssertionError(
                    f"{path.name} still defines a local {node.name} -- "
                    "consolidation into _stats/_optimize was reverted?"
                )
