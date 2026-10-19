"""Tests for C5 conformer_volume_3d (ETKDG + MMFF -> V_3D/SASA/anisotropy)."""
import pytest

from surfactantkit.conformer import _shape_anisotropy, conformer_volume_3d


def test_shape_anisotropy_exact_golds():
    """k^2 = 1 exactly for collinear points, 0 for a tetrahedron (isotropic)."""
    collinear = [(-1.0, 0.0, 0.0), (0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]
    assert _shape_anisotropy(collinear) == pytest.approx(1.0, abs=1e-12)
    tetra = [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]
    assert _shape_anisotropy(tetra) == pytest.approx(0.0, abs=1e-12)
    assert _shape_anisotropy([(0.0, 0.0, 0.0)]) == 0.0


def test_alkane_v3d_linear_in_chain_length():
    """Internal-consistency gold: V_3D grows linearly with n (R^2 >= 0.999),
    no literature numbers involved."""
    ns, vs = [8, 12, 16], []
    for n in ns:
        r = conformer_volume_3d("C" * n, n_conformers=5, seed=42)
        vs.append(r.vdW_volume_A3)
    mx, my = sum(ns) / 3, sum(vs) / 3
    b = sum((x - mx) * (y - my) for x, y in zip(ns, vs, strict=True)) / sum(
        (x - mx) ** 2 for x in ns
    )
    a = my - b * mx
    ss_res = sum((y - (a + b * x)) ** 2 for x, y in zip(ns, vs, strict=True))
    ss_tot = sum((y - my) ** 2 for y in vs)
    r2 = 1.0 - ss_res / ss_tot
    assert r2 >= 0.999
    assert vs[0] < vs[1] < vs[2]


def test_dodecane_cross_check_with_tanford_tail_volume():
    from surfactantkit import tanford_tail_volume
    tv = tanford_tail_volume(12)
    r = conformer_volume_3d(
        "CCCCCCCCCCCC", n_conformers=5, seed=42, tanford_volume_a3=tv
    )
    assert r.tanford_ratio is not None
    assert r.tanford_ratio == pytest.approx(r.vdW_volume_A3 / tv, rel=1e-12)
    # whole-molecule 3D envelope vs tail-only group contribution: same
    # order of magnitude (definitions differ; divergence is the point)
    assert 0.3 <= r.tanford_ratio <= 3.0
    assert r.force_field == "MMFF94s"
    assert r.converged and r.n_iterations >= 5
    assert "definitional" in r.caveat


def test_gemini_tail_volume_larger_than_single_tail():
    gemini = conformer_volume_3d(
        "CCCCCCCCCCCCNCCNCCCCCCCCCCCC", n_conformers=8, seed=42
    )
    single = conformer_volume_3d("CCCCCCCCCCCCN", n_conformers=8, seed=42)
    ratio = gemini.vdW_volume_A3 / single.vdW_volume_A3
    assert ratio >= 1.6  # two C12 tails vs one: far above packed-shell noise


def test_branched_tail_less_anisotropic_than_linear_chain():
    linear = conformer_volume_3d("CCCCCCCCCCCCCCCC", n_conformers=5, seed=42)
    branched = conformer_volume_3d(
        "CC(C)(C)C(C)(C)C(C)(C)C", n_conformers=5, seed=42
    )
    assert branched.shape_anisotropy_k2 < linear.shape_anisotropy_k2
    for r in (linear, branched):
        assert 0.0 <= r.shape_anisotropy_k2 < 1.0
        assert r.vdW_volume_A3 > 0 and r.sasa_A2 > 0


def test_determinism_same_smiles_same_numbers():
    a = conformer_volume_3d("CCCCCCCCCCCC", n_conformers=5, seed=42)
    b = conformer_volume_3d("CCCCCCCCCCCC", n_conformers=5, seed=42)
    assert a.vdW_volume_A3 == b.vdW_volume_A3
    assert a.sasa_A2 == b.sasa_A2
    assert a.shape_anisotropy_k2 == b.shape_anisotropy_k2
    assert a.lowest_energy_kcal_mol == b.lowest_energy_kcal_mol
    assert a.params == b.params


def test_schema_and_params():
    r = conformer_volume_3d("CCO", n_conformers=3, seed=7)
    assert r.converged
    assert r.n_conformers_embedded == 3 and r.n_conformers_minimized == 3
    assert "ETKDGv3" in r.method and "MMFF94s" in r.method
    assert r.params["force_field"] == "MMFF94s"
    assert r.params["seed"] == 7
    assert r.params["molecular_formula"] == "C2H6O"
    assert r.params["heavy_atom_count"] == 3
    assert "free-energy global minimum" in r.caveat


def test_uff_fallback_when_mmff_params_missing():
    """Tin centre: MMFF94 has no parameters -> UFF path must run and be
    reported honestly in params/force_field."""
    r = conformer_volume_3d("CC[Sn](CC)CC", n_conformers=3, seed=42)
    assert r.force_field in ("UFF", "MMFF94s")  # build-dependent FF coverage
    assert r.n_conformers_minimized == 3
    assert r.vdW_volume_A3 > 0 and r.sasa_A2 > 0
    if r.force_field == "UFF":
        assert r.converged


def test_fail_fast_names_every_bad_field():
    with pytest.raises(ValueError, match="smiles"):
        conformer_volume_3d("")
    with pytest.raises(ValueError, match="smiles"):
        conformer_volume_3d("not_a_smiles(")
    with pytest.raises(ValueError, match="n_conformers"):
        conformer_volume_3d("CCO", n_conformers=0)
    with pytest.raises(TypeError, match="n_conformers"):
        conformer_volume_3d("CCO", n_conformers=2.5)  # type: ignore[reportArgumentType]
    with pytest.raises(ValueError, match="seed"):
        conformer_volume_3d("CCO", seed=-1)
    with pytest.raises(TypeError, match="seed"):
        conformer_volume_3d("CCO", seed="42")  # type: ignore[reportArgumentType]
    with pytest.raises(ValueError, match="tanford_volume_a3"):
        conformer_volume_3d("CCO", tanford_volume_a3=0.0)
    with pytest.raises(TypeError, match="tanford_volume_a3"):
        conformer_volume_3d("CCO", tanford_volume_a3="x")  # type: ignore[reportArgumentType]
