"""P1-4 machine-readable validation-tier registry guard.

`surfactantkit.validation_tiers.VALIDATION_TIERS` is GENERATED from the
audit's section 1.1 table; this test re-parses that markdown independently
and fails on any drift (new row, changed tier, changed evidence), plus
schema/spot-gold checks so the registry stays usable as a benchmark
filter (tool -> tier -> citation/test evidence).
"""

import re
from pathlib import Path

from surfactantkit.validation_tiers import (
    TIER_DEFINITIONS,
    VALIDATION_TIERS,
    get_tool_tier,
    tools_in_tier,
)

_ROOT = Path(__file__).resolve().parent.parent
_AUDIT = _ROOT / "docs" / "validation_tier_audit.md"


def _parse_audit_table() -> dict[str, tuple[str, str]]:
    """Independent parser: row -> (tool keys, tier, evidence) from 1.1."""
    md = _AUDIT.read_text(encoding="utf-8")
    tbl = md.split("### 1.1 Full technique \u2192 tier table")[1]
    out: dict[str, tuple[str, str]] = {}
    for line in tbl.splitlines():
        if not line.startswith("|"):
            continue
        if line.startswith("|---") or "Technique (tool/function)" in line:
            continue
        cells = [c.strip() for c in line.strip("|").split(" | ")]
        if len(cells) < 3:
            continue
        if not (cells[0].startswith("`") or cells[0].startswith("*")):
            continue
        tier = cells[1].replace("**", "").strip()
        evidence = cells[2]
        for bt in re.findall(r"`([^`]+)`", cells[0]):
            names = ["hld_ionic", "hld_nonionic"] if bt == "hld_ionic/nonionic" else [bt]
            for n in names:
                assert n not in out, f"parser duplicate {n}"
                out[n] = (tier, evidence)
    return out


def test_registry_matches_audit_table_exactly():
    expected = _parse_audit_table()
    assert len(expected) >= 60, f"audit table parsed only {len(expected)} rows"
    assert set(expected) == set(VALIDATION_TIERS), (
        "registry keys drift from audit table: "
        f"only-in-registry={set(VALIDATION_TIERS) - set(expected)} "
        f"only-in-audit={set(expected) - set(VALIDATION_TIERS)}"
    )
    for name, (tier, evidence) in expected.items():
        rec = VALIDATION_TIERS[name]
        assert rec["tier"] == tier, f"{name}: {rec['tier']} != audit {tier}"
        assert rec["evidence"] == evidence, f"{name}: evidence drifted from audit"


def test_schema():
    assert len(VALIDATION_TIERS) >= 60
    allowed = set(TIER_DEFINITIONS)
    for name, rec in VALIDATION_TIERS.items():
        assert rec["tier"] in allowed, f"{name}: unknown tier {rec['tier']!r}"
        assert rec["evidence"].strip(), f"{name}: empty evidence"
        assert set(rec) == {"tier", "evidence"}


def test_spot_golds():
    assert VALIDATION_TIERS["clint_ideal_cmc"]["tier"] == "T1"
    assert "PMC6554738" in VALIDATION_TIERS["clint_ideal_cmc"]["evidence"]
    assert VALIDATION_TIERS["rubingh_beta"]["tier"] == "T1\u2020"
    assert VALIDATION_TIERS["rubingh_beta_regression"]["tier"] == "T2"
    vh = VALIDATION_TIERS["vant_hoff_enthalpy"]["tier"]
    assert vh.startswith("T1") and "T2" in vh  # the split T1/T2 row
    assert VALIDATION_TIERS["entropy_of_micellization"]["tier"] == "dependent"
    assert VALIDATION_TIERS["derive_all_properties_from_smiles_and_curve"]["tier"] == "chained"
    # T4 dead ends are deliberately absent (ship no tool at all)
    assert get_tool_tier("gibbs_n_factor") is None


def test_tools_in_tier_filtering():
    t1 = tools_in_tier("T1")
    assert "clint_ideal_cmc" in t1
    assert "rubingh_beta_regression" not in t1
    both = set(tools_in_tier("T1", "T2"))
    assert set(t1) <= both
    assert "rubingh_beta_regression" in both
    # T2 set is the audit's documented "no literature backing" class
    assert len(tools_in_tier("T2")) >= 20


def test_tools_in_tier_unknown_tier_fails_fast():
    try:
        tools_in_tier("T9")
    except ValueError as exc:
        assert "T9" in str(exc)
    else:
        raise AssertionError("unknown tier must raise")


def test_get_tool_tier_exact_and_wildcard_lookup():
    assert get_tool_tier("clint_ideal_cmc") == "T1"
    # audit-table wildcard rows resolve via fnmatch
    assert get_tool_tier("tanford_chain_length") == VALIDATION_TIERS["tanford_*"]["tier"]
    assert get_tool_tier("anything_reference") == VALIDATION_TIERS["*_reference"]["tier"]
    assert get_tool_tier("derive_davies_group_number_imidazoline") == "T3"
    assert get_tool_tier("no_such_tool") is None
