"""Tests for M1 analyze_micelle_trajectory (analysis-only worker)."""
import math

import numpy as np
import pytest

from surfactantkit.trajectory import analyze_micelle_trajectory


def _lattice_blob(n, center, spacing=0.35, jitter=0.03, seed=0):
    """Cubic-lattice aggregate: single-linkage connectivity at cutoff 0.45
    is GUARANTEED (spacing + 2*jitter < cutoff), sizes exactly n."""
    rng = np.random.default_rng(seed)
    side = int(np.ceil(n ** (1.0 / 3.0)))
    grid = np.array([
        [i, j, k]
        for i in range(side)
        for j in range(side)
        for k in range(side)
    ], dtype=float)[:n]
    pts = (grid - (side - 1) / 2.0) * spacing + np.asarray(center)
    return pts + rng.normal(0.0, jitter, size=pts.shape)


def _make_fixture(tmp_path, seed=42, kind="one_cluster"):
    """Seeded synthetic trajectory: deterministic across runs. Aggregates
    and monomers sit on jittered lattices so single-linkage clustering has
    an exactly known answer (no flaky Gaussian outliers, no accidental
    monomer merges)."""
    frames = []
    if kind == "one_cluster":
        for f in range(5):
            # everything inside the 12^3 box with origin 0 (box semantics)
            blob = _lattice_blob(100, (4.0, 6.0, 6.0), seed=100 + f)
            idx = np.arange(40)
            # monomers: clear of the blob (>= 0.6 gap) and of each other
            lattice = np.stack([
                6.0 + 0.7 * (idx % 8),          # 6.0 .. 10.9
                1.0 + 2.0 * ((idx // 8) % 5),   # 1.0 .. 9.0
                np.full(40, 10.0),
            ], axis=1).astype(float)
            rng = np.random.default_rng(seed + f)
            free = lattice + rng.normal(0.0, 0.05, size=lattice.shape)
            frames.append(np.vstack([blob, free]))
    elif kind == "two_clusters":
        for f in range(4):
            a = _lattice_blob(50, (-2.0, 0.0, 0.0), seed=200 + f)
            b = _lattice_blob(50, (2.0, 0.0, 0.0), seed=300 + f)
            idx = np.arange(10)
            lattice = np.stack([
                -5.0 + 1.5 * idx,
                4.0 * np.ones(10),
                np.zeros(10),
            ], axis=1).astype(float)
            rng = np.random.default_rng(seed + f)
            free = lattice + rng.normal(0.0, 0.05, size=lattice.shape)
            frames.append(np.vstack([a, b, free]))
    return np.stack(frames)


def test_seeded_fixture_npz_sizes_and_rg(tmp_path):
    pos = _make_fixture(tmp_path, kind="one_cluster")
    p = tmp_path / "traj.npz"
    np.savez(p, positions=pos, box=np.array([12.0, 12.0, 12.0]))
    r = analyze_micelle_trajectory(str(p), cluster_cutoff=0.45)
    assert r.converged and r.n_iterations == 5 and r.n_frames == 5
    # single compact cluster of 100 + 40 monomers, every frame
    assert set(r.largest_cluster_size) == {100}
    assert set(r.n_clusters) == {41}
    assert all(sz[0] == 100 for sz in r.cluster_sizes)
    # Rg gold: hand-computed unweighted Rg of frame-0 rows 0..99 (the blob)
    sub = pos[0][:100]
    com = sub.mean(axis=0)
    rg_hand = math.sqrt(float(np.mean(np.sum((sub - com) ** 2, axis=1))))
    assert r.rg_largest[0] == pytest.approx(rg_hand, rel=1e-12)
    assert all(g > 0.2 for g in r.rg_largest)
    # radial profile: positive density near core
    assert len(r.radial_density_r) == 30
    # several populated shells; core-adjacent shells are densest region
    assert max(r.radial_density_values) > 0
    assert sum(1 for v in r.radial_density_values if v > 0) >= 5
    # slab profile present (box given): total particle-frame density conserved
    assert r.slab_density_z is not None and len(r.slab_density_z) == 50
    dz = 12.0 / 50
    total = sum(r.slab_density_z) * dz * 12.0 * 12.0
    assert total == pytest.approx(140.0, rel=1e-9)


def test_two_clusters_and_separation(tmp_path):
    pos = _make_fixture(tmp_path, kind="two_clusters")
    p = tmp_path / "traj.npz"
    np.savez(p, positions=pos)
    r = analyze_micelle_trajectory(str(p), cluster_cutoff=0.45)
    assert set(r.largest_cluster_size) == {50}
    assert set(r.n_clusters) == {12}  # 2 clusters + 10 monomers
    assert r.slab_density_z is None
    assert "no box array" in r.caveat


def test_txt_xyz_matches_npz(tmp_path):
    pos = _make_fixture(tmp_path, kind="one_cluster")
    p_npz = tmp_path / "traj.npz"
    np.savez(p_npz, positions=pos)
    lines = []
    for f in pos:
        lines.append(str(len(f)))
        lines.append("comment")
        for x, y, z in f:
            lines.append(f"C {x:.17g} {y:.17g} {z:.17g}")
    p_txt = tmp_path / "traj.xyz"
    p_txt.write_text("\n".join(lines))
    a = analyze_micelle_trajectory(str(p_npz), cluster_cutoff=0.45)
    b = analyze_micelle_trajectory(str(p_txt), cluster_cutoff=0.45)
    assert a.largest_cluster_size == b.largest_cluster_size
    assert a.cluster_sizes == b.cluster_sizes
    assert a.rg_largest == pytest.approx(b.rg_largest, rel=1e-12)
    assert a.n_clusters == b.n_clusters


def test_txt_bare_rows_parse(tmp_path):
    pos = _make_fixture(tmp_path, kind="two_clusters")
    blocks = []
    for f in pos:
        rows = [" ".join(f"{v:.8f}" for v in row) for row in f]
        blocks.append("\n".join(rows))
    p = tmp_path / "traj.txt"
    p.write_text("\n\n".join(blocks))
    r = analyze_micelle_trajectory(str(p), cluster_cutoff=0.45)
    assert set(r.largest_cluster_size) == {50}


def test_mass_weighted_rg_exact_gold(tmp_path):
    """Hand-computed weighted Rg from the same fixture: exact match."""
    pos = _make_fixture(tmp_path, kind="one_cluster")
    masses = np.linspace(0.5, 2.5, pos.shape[1])
    p = tmp_path / "traj.npz"
    np.savez(p, positions=pos, masses=masses)
    r = analyze_micelle_trajectory(str(p), cluster_cutoff=0.45)
    # hand Rg for frame 0: find the 100 largest cluster = first 100 rows
    # (blob built first in fixture); weights equal across blob rows
    sub = pos[0][:100]
    w = masses[:100]
    com = (w[:, None] * sub).sum(axis=0) / w.sum()
    rg_hand = math.sqrt(float((w * np.sum((sub - com) ** 2, axis=1)).sum()) / w.sum())
    assert r.rg_largest[0] == pytest.approx(rg_hand, rel=1e-9)
    assert r.params["mass_weighted_rg"] is True


def test_determinism_same_input_same_output(tmp_path):
    pos = _make_fixture(tmp_path, kind="one_cluster")
    p = tmp_path / "traj.npz"
    np.savez(p, positions=pos)
    a = analyze_micelle_trajectory(str(p), cluster_cutoff=0.45)
    b = analyze_micelle_trajectory(str(p), cluster_cutoff=0.45)
    assert a.largest_cluster_size == b.largest_cluster_size
    assert a.rg_largest == b.rg_largest
    assert a.radial_density_values == b.radial_density_values
    assert a.n_iterations == b.n_iterations


def test_fail_fast_field_names(tmp_path):
    pos = _make_fixture(tmp_path, kind="two_clusters")
    good = tmp_path / "traj.npz"
    np.savez(good, positions=pos)
    with pytest.raises(ValueError, match="cluster_cutoff"):
        analyze_micelle_trajectory(str(good), cluster_cutoff=0.0)
    with pytest.raises(TypeError, match="n_radial_bins"):
        analyze_micelle_trajectory(
            str(good), cluster_cutoff=0.4,
            n_radial_bins="30",  # type: ignore[reportArgumentType]
        )
    with pytest.raises(ValueError, match="n_radial_bins"):
        analyze_micelle_trajectory(str(good), cluster_cutoff=0.4,
                                   n_radial_bins=2)
    with pytest.raises(TypeError, match="compute_slab_profile"):
        analyze_micelle_trajectory(
            str(good), cluster_cutoff=0.4,
            compute_slab_profile=1,  # type: ignore[reportArgumentType]
        )
    # missing positions key
    bad = tmp_path / "bad.npz"
    np.savez(bad, coords=pos)
    with pytest.raises(ValueError, match="positions"):
        analyze_micelle_trajectory(str(bad), cluster_cutoff=0.4)
    # NaN positions
    nan_pos = pos.copy()
    nan_pos[0, 0, 0] = math.nan
    bad2 = tmp_path / "nan.npz"
    np.savez(bad2, positions=nan_pos)
    with pytest.raises(ValueError, match="NaN"):
        analyze_micelle_trajectory(str(bad2), cluster_cutoff=0.4)
    # missing file
    with pytest.raises(ValueError, match="cannot read"):
        analyze_micelle_trajectory(str(tmp_path / "none.npz"),
                                   cluster_cutoff=0.4)
    # ragged txt
    rag = tmp_path / "rag.txt"
    rag.write_text("0 0 0\n1 1 1\n\n0 0 0\n")
    with pytest.raises(ValueError, match="differing particle counts"):
        analyze_micelle_trajectory(str(rag), cluster_cutoff=0.4)


def test_schema_and_caveats(tmp_path):
    pos = _make_fixture(tmp_path, kind="one_cluster")
    p = tmp_path / "traj.npz"
    np.savez(p, positions=pos)
    r = analyze_micelle_trajectory(str(p), cluster_cutoff=0.45)
    assert r.converged and r.n_iterations == r.n_frames
    assert "single-linkage" in r.caveat and "cutoff sensitivity" in r.caveat
    assert "analysis only" in r.caveat
    assert r.method and "single-linkage" in r.method
    assert r.params["cluster_cutoff"] == 0.45
    assert r.params["input_shape"] == [5, 140, 3]
    assert len(r.times) == 5 and r.times[0] == 0.0
