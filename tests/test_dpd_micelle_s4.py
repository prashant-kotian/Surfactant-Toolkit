"""S4 full-worker tests: in-house SI-model DPD (C8E4) + async jobs.

Physics-adjacent claims are NOT asserted here (interactive boxes are
finite-size limited by design -- the worker says so).  What IS pinned:

* force invariants: exact momentum conservation of conservative and
  thermostat pair forces; Venturoli angle forces match a finite
  difference of the potential (the analytic chain rule is verified,
  not trusted);
* the SI model constants themselves (Table S1 spot golds, r0 anchor,
  unit conversion of the box at a given molarity);
* mechanics: schema, deterministic seeding, finite-size gate honesty,
  scale/budget refusals, async submit -> status -> result round trip
  (including the failed-job path), and all input guards.
"""
import math
import time

import numpy as np
import pytest

from surfactantkit.dpd_micelle import (
    _N_PER,
    _SURFACTANT_TYPES,
    _TABLE_S1,
    _bond_angle_arrays,
    _box_for,
    _conservative_forces,
    _pair_tables,
    _place_chains,
    _r0_for,
    _thermostat_forces,
    dpd_micelle_cmc,
    dpd_micelle_cmc_job_status,
    estimate_dpd_runtime_s,
    submit_dpd_micelle_cmc_job,
)

# fast test configuration: N ~ 200 beads, ~2 s per run
TINY = {"concentrations_M": [1.0], "n_surfactants": 8, "steps": 200,
        "n_seeds": 3, "seed": 5, "max_sync_seconds": 600.0}


def _mini_system(n_s=8, c=1.0, seed=5):
    l_box = _box_for(n_s, c)
    n_beads = max(round(3.0 * l_box ** 3), n_s * _N_PER + 8)
    rng = np.random.default_rng(seed)
    x = _place_chains(n_s, n_beads, l_box, rng, np)
    types = np.array(list(_SURFACTANT_TYPES) * n_s
                     + [4] * (n_beads - _N_PER * n_s), dtype=np.int64)
    bonds, angles = _bond_angle_arrays(n_s, np)
    a_mat, r_mat = _pair_tables(np)
    rows, cols = np.triu_indices(n_beads, k=1)
    return x, types, l_box, bonds, angles, a_mat, r_mat, rows, cols


def test_si_model_constants():
    # Table S1 spot golds (transcribed from the paper's Supporting Info)
    table = {(i, j): (a, r) for i, j, a, r in _TABLE_S1}
    assert table[(1, 4)] == (45.0, 1.037)     # CH2CH2 - 2H2O (poor solvent)
    assert table[(4, 4)] == (25.0, 1.0)       # solvent standard
    assert table[(1, 1)] == (22.0, 1.074)     # tail - tail
    assert table[(3, 4)] == (14.5, 0.99)      # head - water (good solvent)
    a_mat, r_mat = _pair_tables(np)
    assert np.allclose(a_mat, a_mat.T) and np.allclose(r_mat, r_mat.T)
    # r0 anchor: SI says 0.39 for CH2CH2-CH2CH2 (combined heavy atoms = 4)
    assert _r0_for(1, 1) == pytest.approx(0.39)
    # linear in combined heavy atoms (SI A.2 stated rule)
    assert _r0_for(0, 1) == pytest.approx(0.39 * 3 / 4)
    # C8E4 composition: 9 beads, 3 CH2CH2 tails
    assert len(_SURFACTANT_TYPES) == 9
    assert _SURFACTANT_TYPES.count(1) == 3


def test_box_unit_conversion_exact():
    c, n_s = 0.6, 20
    l_box = _box_for(n_s, c)
    v_l = l_box ** 3 * (5.64 ** 3) * 1e-27
    assert n_s / (v_l * 6.02214076e23) == pytest.approx(c, rel=1e-12)


def test_force_invariants_momentum_conservation():
    x, types, l_box, bonds, angles, a_mat, r_mat, rows, cols = _mini_system()
    f = _conservative_forces(x, types, l_box, bonds, angles, a_mat, r_mat,
                             rows, cols, np)
    assert np.abs(f.sum(axis=0)).max() < 1e-8
    v = np.random.default_rng(1).normal(size=x.shape)
    ft = _thermostat_forces(x, v, l_box, 1.116, np.random.default_rng(2),
                            0.04, rows, cols, np)
    assert np.abs(ft.sum(axis=0)).max() < 1e-8


def test_angle_force_matches_finite_difference():
    x, types, l_box, _bonds, angles, a_mat, r_mat, rows, cols = _mini_system()
    # isolate angle forces (no bonds) and compare one bead against FD
    f_angle = _conservative_forces(x, types, l_box,
                                   np.zeros((0, 2), dtype=np.int64), angles,
                                   a_mat, r_mat, rows, cols, np)

    def angle_energy(pos):
        e = 0.0
        for i, j, k in angles:
            a = pos[i] - pos[j]
            a -= l_box * np.round(a / l_box)
            b = pos[k] - pos[j]
            b -= l_box * np.round(b / l_box)
            na = np.linalg.norm(a)
            nb = np.linalg.norm(b)
            ct = float(np.clip(np.dot(a, b) / (na * nb), -1.0, 1.0))
            th = math.acos(ct)
            e += 0.5 * 5.0 * (th - math.pi) ** 2
        return e

    h = 1e-6
    bead = int(angles[3, 0])
    fd = np.zeros(3)
    for c in range(3):
        xp = x.copy()
        xp[bead, c] += h
        xm = x.copy()
        xm[bead, c] -= h
        fd[c] = -(angle_energy(xp) - angle_energy(xm)) / (2 * h)
    assert np.allclose(f_angle[bead], fd, atol=1e-4)


def test_tiny_run_schema_disclosures_and_finite_size_gate():
    res = dpd_micelle_cmc(**TINY)
    for k in ("cmc_M", "cmc_sigma_M", "cmc_spread_log10", "per_concentration",
              "converged", "n_iterations", "method", "params", "caveat"):
        assert k in res
    assert res["method"] == "dpd_c8e4_si_model_anderson_pfree_eq3"
    p = res["params"]
    assert p["scale"] == "interactive"
    assert p["gamma"] == 4.5 and p["sigma_length_A"] == 5.64
    assert p["n_submicellar_cutoff"] == 6
    assert res["n_iterations"] == 200 * 3
    cav = res["caveat"]
    assert "doi:10.1021/acs.jpcb.5c07475" in cav
    assert "NOT dl_meso" in cav or "not dl_meso" in cav
    assert "FINITE SIZE" in cav
    pc = res["per_concentration"][0]
    # deterministic trace values in range
    for run in pc["runs"]:
        assert 0.0 <= run["p_free_first"] <= 1.0
        assert 0.0 <= run["p_free_last"] <= 1.0
    # the finite-size gate fires for this tiny box (honest refusal to
    # endorse an unresolvable plateau) -- exact by determinism
    assert pc["finite_size_ok"] is False
    assert "batch" in pc["finite_size_warning"]
    assert res["converged"] is False


def test_determinism_same_seed_bit_identical():
    a = dpd_micelle_cmc(**TINY)
    b = dpd_micelle_cmc(**TINY)
    # wall_seconds/estimate differ; everything physics-facing must match
    for key in ("cmc_M", "converged", "n_iterations", "method"):
        assert a[key] == b[key], key
    assert (a["per_concentration"][0]["runs"]
            == b["per_concentration"][0]["runs"])
    assert (a["per_concentration"][0]["p_free_inf"]
            == b["per_concentration"][0]["p_free_inf"])


def test_guards_and_refusals():
    with pytest.raises(ValueError, match="batch"):
        dpd_micelle_cmc(scale="paper")
    with pytest.raises(ValueError, match="concentrations_M"):
        dpd_micelle_cmc(concentrations_M=[])
    with pytest.raises(ValueError, match="concentrations_M"):
        dpd_micelle_cmc(concentrations_M=[5.0])
    with pytest.raises(ValueError, match="n_surfactants"):
        dpd_micelle_cmc(n_surfactants=4, max_sync_seconds=1e9)
    with pytest.raises(ValueError, match="steps"):
        dpd_micelle_cmc(steps=50, max_sync_seconds=1e9)
    with pytest.raises(ValueError, match="n_seeds"):
        dpd_micelle_cmc(n_seeds=2, max_sync_seconds=1e9)
    with pytest.raises(ValueError, match="dt"):
        dpd_micelle_cmc(dt=0.5, max_sync_seconds=1e9)
    # budget refusal names the async path (no simulation starts)
    with pytest.raises(ValueError, match="submit_dpd_micelle_cmc_job"):
        dpd_micelle_cmc(max_sync_seconds=1e-6)


def test_estimate_is_positive_and_reasonable():
    est = estimate_dpd_runtime_s(steps=200, n_surfactants=8, c_tot_M=1.0,
                                 n_runs=3)
    assert est > 0.0
    assert math.isfinite(est)


def test_async_roundtrip_matches_sync(tmp_path):
    args = dict(TINY)
    args["jobs_dir"] = str(tmp_path)
    sub = submit_dpd_micelle_cmc_job(
        concentrations_M=args["concentrations_M"],
        n_surfactants=args["n_surfactants"], steps=args["steps"],
        n_seeds=args["n_seeds"], seed=args["seed"],
        jobs_dir=args["jobs_dir"])
    assert sub["state"] == "queued" and len(sub["job_id"]) == 12
    # job file exists immediately (persisted on submit)
    assert (tmp_path / f"{sub['job_id']}.json").exists()
    deadline = time.time() + 120.0
    status = None
    while time.time() < deadline:
        status = dpd_micelle_cmc_job_status(sub["job_id"],
                                            jobs_dir=str(tmp_path))
        if status["state"] in ("completed", "failed", "interrupted"):
            break
        time.sleep(0.3)
    assert status is not None and status["state"] == "completed", status
    sync = dpd_micelle_cmc(**TINY)
    async_res = status["result"]
    assert async_res["cmc_M"] == sync["cmc_M"]
    assert (async_res["per_concentration"][0]["runs"]
            == sync["per_concentration"][0]["runs"])


def test_async_failed_job_reports_error(tmp_path):
    sub = submit_dpd_micelle_cmc_job(concentrations_M=[1.0],
                                     n_surfactants=8, steps=200,
                                     n_seeds=1, seed=5,
                                     jobs_dir=str(tmp_path))
    deadline = time.time() + 30.0
    status = None
    while time.time() < deadline:
        status = dpd_micelle_cmc_job_status(sub["job_id"],
                                            jobs_dir=str(tmp_path))
        if status["state"] in ("completed", "failed"):
            break
        time.sleep(0.2)
    assert status["state"] == "failed"
    assert "n_seeds" in (status["error"] or "")


def test_job_status_guards(tmp_path):
    with pytest.raises(ValueError, match="job_id"):
        dpd_micelle_cmc_job_status("")
    with pytest.raises(ValueError, match="unknown job_id"):
        dpd_micelle_cmc_job_status("deadbeef0000", jobs_dir=str(tmp_path))
    # path-like ids are refused before any file is touched
    for bad in ("../../secret", "..\\x", "DEADBEEF0000", "deadbeef00001"):
        with pytest.raises(ValueError, match="12-character hex"):
            dpd_micelle_cmc_job_status(bad, jobs_dir=str(tmp_path))
