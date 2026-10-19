"""Tests for surfactantkit.eacn (worker C3).

Gold (sourced rules, see eacn.py citations):
- n-alkanes: EACN = N_C exactly (Wade 1977 definitional);
- octylbenzene 8 (Queste/Benzoate 2007), butylcyclohexane 8 (4+4);
- single long-chain ester: N_C - 11 with band [-10, -12] -> C20 -> 9, [8,10];
- refusals: branched chains, bare rings, ring/ester combos, heteroatoms,
  multi-fragment, short esters, invalid SMILES (never fabricate);
- linear mixing midpoint + HLD feed monotonicity (higher EACN -> higher lnS*).
"""
import pytest

from surfactantkit.eacn import (
    EacnResult,
    eacn_group_contribution,
    eacn_linear_mixing,
)


def test_n_alkanes_exact():
    for smi, want in [("CCCCCC", 6), ("CCCCCCCCCC", 10),
                      ("CCCCCCCCCCCCCCCC", 16)]:
        r = eacn_group_contribution(smi)
        assert isinstance(r, EacnResult)
        assert (r.eacn, r.eacn_low, r.eacn_high) == (want, want, want)
        assert r.rule == "n_alkane"
        assert r.sources and r.caveat and r.method
        d = r.as_dict()
        for k in ("converged", "method", "params", "caveat"):
            assert k in d


def test_alkylbenzene_ring_zero_rule():
    # Queste 2007: benzene ring contribution 0; octylbenzene = 8 exactly
    r = eacn_group_contribution("CCCCCCCCc1ccccc1")
    assert r.eacn == 8 and r.eacn_low == r.eacn_high == 8
    assert r.rule == "alkylbenzene_sidechain"


def test_alkylcyclohexane_additive():
    # butylcyclohexane: side chain 4 + ring 4 = 8 (Queste/Benzoate 2007)
    r = eacn_group_contribution("CCCCC1CCCCC1")
    assert r.eacn == 8 and r.rule.startswith("alkylcyclohexane_sidechain")
    assert (r.eacn_low, r.eacn_high) == (8, 8)


def test_single_ester_band():
    # C20 single ester: N_C - 11 = 9; Ester 2013 increment -10..-12 -> band
    r = eacn_group_contribution("CCCCCCCCCCCCCCCCCC(=O)OCC")
    assert r.eacn == 9
    assert (r.eacn_low, r.eacn_high) == (8, 10)
    assert r.rule == "ester_single_acyclic"
    assert r.eacn_low < r.eacn < r.eacn_high


def test_refusals_never_fabricate():
    cases = {
        "CC(C)C": "branch",                     # branched alkane
        "CCCCCO": "heteroatom",                 # 1-hexanol
        "CCOCC": "heteroatom",                  # ether
        "c1ccccc1": "side chain",               # benzene, no chain
        "C1CCCCC1": "side chain",               # bare cyclohexane
        "CC(C)c1ccccc1": "branch",              # isopropylbenzene
        "Cc1ccccc1C": "substituent",            # o-xylene, 2 chains
        "c1ccc2ccccc2c1": "ring",               # fused rings
        "CCOC(=O)C=C": "long-chain",            # short-chain ester
        "CC(=O)OCC": "long-chain",              # ethyl acetate
        "CC(=O)O.OCC": "fragment",              # multi-fragment
        "not_a_smiles": "SMILES",               # invalid
        "": "non-empty",
    }
    for smi, pat in cases.items():
        with pytest.raises(ValueError) as e:
            eacn_group_contribution(smi)
        assert pat.lower() in str(e.value).lower(), (smi, str(e.value))
    # refusals must tell the user what to do instead
    with pytest.raises(ValueError, match=r"measure|T\*|estimator"):
        eacn_group_contribution("CC(C)C")


def test_unsaturated_chains_refused():
    # EACN = ACN is definitional only for SATURATED n-alkanes (Wade 1977);
    # alkenes/alkynes, non-aromatic ring double bonds and unsaturated
    # esters were never sourced -> refuse, never fabricate (saturation
    # gate added 2026-10-05 after 1-octene slipped through as "n_alkane").
    cases = {
        "CCCCCCC=C": "unsaturated",                       # 1-octene
        "CCCCCCCC/C=C/CCCC": "unsaturated",               # trans-5-decene
        "CCCCCCCCC#C": "unsaturated",                     # 1-decyne
        "C1=CCCCC1": "unsaturated",                       # cyclohexene
        "CCCCCCCCCCCC/C=C/CCCCCCCC(=O)OCC": "unsaturated",  # methyl oleate
        "C=CCc1ccccc1": "unsaturated",                    # allylbenzene
    }
    for smi, pat in cases.items():
        with pytest.raises(ValueError) as e:
            eacn_group_contribution(smi)
        assert pat in str(e.value).lower(), (smi, str(e.value))
    # guidance to measure, not guess, must survive the gate
    with pytest.raises(ValueError, match=r"measure|T\*|estimator"):
        eacn_group_contribution("CCCCCCC=C")
    # saturated controls still pass untouched
    assert eacn_group_contribution("CCCCCCCC").eacn == 8
    assert eacn_group_contribution("CCCCCCCCc1ccccc1").eacn == 8
    r = eacn_group_contribution("CCCCCCCCCCCCCCCCCC(=O)OCC")
    assert r.rule == "ester_single_acyclic" and r.eacn == 9


def test_linear_mixing():
    m = eacn_linear_mixing([0.5, 0.5], [16.0, 8.0])
    assert m["eacn"] == pytest.approx(12.0)
    assert m["basis"] == "mole"
    assert m["converged"] and m["method"] and m["caveat"]
    mv = eacn_linear_mixing([1, 1, 1], [6.0, 10.0, 12.0], basis="volume")
    assert mv["eacn"] == pytest.approx(28.0 / 3.0, rel=1e-6)
    # fractions are normalized (any positive scaling of a blend is the blend)
    mn = eacn_linear_mixing([0.3, 0.3], [16.0, 8.0])
    assert mn["eacn"] == pytest.approx(12.0)
    with pytest.raises(ValueError, match="equal length"):
        eacn_linear_mixing([1.0], [1.0, 2.0])
    with pytest.raises(ValueError, match="non-negative"):
        eacn_linear_mixing([-1.0, 2.0], [1.0, 2.0])
    with pytest.raises(ValueError, match="positive"):
        eacn_linear_mixing([0.0, 0.0], [1.0, 2.0])
    with pytest.raises(ValueError, match="'mole'"):
        eacn_linear_mixing([1.0], [1.0], basis="mass")


def test_hld_feed_monotonicity():
    # higher EACN -> higher optimal salinity (ionic HLD, k>0)
    from surfactantkit.hld import optimal_salinity_ionic
    s8 = optimal_salinity_ionic(k=0.06, eacn=8.0, cc=1.0)
    s10 = optimal_salinity_ionic(k=0.06, eacn=10.0, cc=1.0)
    assert s10 > s8
    # feed a computed EACN straight through
    r = eacn_group_contribution("CCCCCCCCc1ccccc1")
    s_feed = optimal_salinity_ionic(k=0.06, eacn=r.eacn, cc=1.0)
    assert s_feed == pytest.approx(
        optimal_salinity_ionic(k=0.06, eacn=8.0, cc=1.0))
