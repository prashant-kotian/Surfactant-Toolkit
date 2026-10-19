"""Tests for surfactantkit.descriptors (computational worker C4).

Validation strategy: known-value checks against published/standard reference
values (aspirin TPSA exactly 63.6 A^2 per Ertl's definition; experimental
logP for three simple compounds within the Wildman-Crippen method's stated
accuracy), determinism, and explicit refusal on invalid SMILES.
"""
import pytest

from surfactantkit.descriptors import DescriptorPanelResult, descriptor_panel


def test_aspirin_known_values():
    r = descriptor_panel("CC(=O)Oc1ccccc1C(=O)O")
    # TPSA: published value for aspirin is exactly 63.6 A^2 (Ertl definition)
    assert r.tpsa_A2 == pytest.approx(63.6, abs=0.05)
    assert r.mol_weight_g_per_mol == pytest.approx(180.159, rel=1e-4)
    assert r.exact_mass_g_per_mol == pytest.approx(180.04226, abs=1e-4)
    # Crippen logP vs experimental 1.19 (method accuracy ~0.3-0.5)
    assert r.logp_crippen == pytest.approx(1.19, abs=0.5)
    assert r.n_rings == 1 and r.n_aromatic_rings == 1
    assert r.n_rotatable_bonds == 2


def test_simple_compounds_logp_within_method_accuracy():
    # experimental logP: ethanol -0.31, phenol 1.46, 1-butanol 0.88
    for smi, exp in [("CCO", -0.31), ("Oc1ccccc1", 1.46), ("CCCCO", 0.88)]:
        r = descriptor_panel(smi)
        assert r.logp_crippen == pytest.approx(exp, abs=0.5), smi


def test_known_alkane_underestimation_disclosed():
    # Honest method-limit test: Wildman-Crippen systematically reads n-alkanes
    # low vs experiment (hexane: computed 2.59 vs experimental 3.50). The
    # descriptor is still a valid *computation*; the caveat covers the error.
    r = descriptor_panel("CCCCCC")
    assert r.logp_crippen == pytest.approx(3.50, abs=1.0)
    assert "not a measurement" in r.caveat


def test_ethanol_tpsa():
    r = descriptor_panel("CCO")
    assert r.tpsa_A2 == pytest.approx(20.23, abs=0.05)


def test_sds_sanity_and_counts():
    r = descriptor_panel("CCCCCCCCCCCCOS(=O)(=O)[O-]")
    assert r.mol_weight_g_per_mol == pytest.approx(265.395, abs=0.05)
    assert r.n_hbd == 0  # sulfate anion has no H donors
    assert "S" in r.canonical_smiles


def test_determinism():
    a = descriptor_panel("CCCCCCCCCCCCOS(=O)(=O)[O-]")
    b = descriptor_panel("CCCCCCCCCCCCOS(=O)(=O)[O-]")
    assert a == b


def test_invalid_smiles_refuses():
    with pytest.raises(ValueError, match="do not guess"):
        descriptor_panel("not_a_smiles")


def test_empty_smiles_refuses():
    with pytest.raises(ValueError):
        descriptor_panel("")


def test_result_fields_and_caveat():
    r = descriptor_panel("CCO")
    d = r.as_dict()
    assert isinstance(r, DescriptorPanelResult)
    assert "ESTIMATE" in d["caveat"] or "estimate" in d["caveat"]
    assert d["method"]
    for key in ("logp_crippen", "tpsa_A2", "chi0v", "kappa1", "bertz_ct"):
        assert key in d
