"""S4 full worker: DPD CMC prediction for nonionic C8E4 (reduced units).

The published in-silico CMC method of Rummel, Robinson, Warren, Bray &
Anderson, J. Phys. Chem. B 2026 (doi 10.1021/acs.jpcb.5c07475, CC-BY,
PMC12794194) -- simulation side.  The model is the paper's OWN DPD model
from its Supporting Information (Table S1 + sections A.1/A.2, obtained
via Europe PMC supplementary files):

* conservative pair potential U = 1/2 A_ij R_ij (1 - r/R_ij)^2 for
  r <= R_ij (force A_ij (1 - r/R_ij)); the full A_ij/R_ij table for the
  C8E4 + 2H2O bead types is transcribed into _TABLE_S1 (parameters from
  the authors' logP/solubility/density fits, Anderson et al. JCP 2017),
* solvent (2H2O) A = 25, R = 1; DPD length unit = 5.64 A,
* dissipative gamma = 4.5 (sigma = sqrt(2 gamma kT) = 3, kT = 1), the
  standard Groot-Warren thermostat with the SI's cut-off convention
  (max R_ii), velocity Verlet, dt exposed (the SI does not state its
  dt; 0.04 is the conventional stable value),
* harmonic bonds k_b = 150 with r0 proportional to combined heavy-atom
  count (SI anchor: 0.39 DPD units for CH2CH2-CH2CH2, combined 4),
  Venturoli three-body angle theta_0 = 180 deg, k_a = 5,
* surfactant = C8E4 chain [CH3]-[CH2CH2]3-[CH2OCH2]4-[CH2OH]
  (9 beads, 306 g/mol; shorter chain than the paper's C10E6 so
  its higher CMC is finite-size resolvable in interactive boxes
  -- same bead types, same Table S1 parameters), submicellar cut-off n* = 6 molecules (paper
  section 4.2), free-monomer ratio p_free(t) fed to the eq-8 + eq-4
  analysis in this same package (anderson_cmc_analysis).

Scope honesty: an IN-HOUSE numpy reimplementation of the paper's model,
not the dl_meso production engine.  Interactive boxes are small
(N ~ 300-4000 beads, minutes) versus the paper's L = 34-54 boxes with
150-600 surfactants and 20k-90k time-unit runs (scale="paper" refuses
and names the batch-class constraint).  The async job pair
(submit_dpd_micelle_cmc_job / dpd_micelle_cmc_job_status) runs longer
interactive-scale jobs without blocking the server.  Every output
carries these disclosures in `caveat`/`params`.

Deterministic for a fixed seed.  Requires numpy (lazy import).
"""

from __future__ import annotations

import json
import math
import os
import threading
import time
from pathlib import Path

from .dpd_cmc import anderson_cmc_analysis

# ---------------------------------------------------------------------------
# Model constants (SI sections A.1/A.2; C10E6 chemistry)
# ---------------------------------------------------------------------------

# bead type ids: 0 CH3, 1 CH2CH2, 2 CH2OCH2, 3 CH2OH, 4 2H2O
_BEAD_NAMES = ("CH3", "CH2CH2", "CH2OCH2", "CH2OH", "2H2O")
# Table S1 rows restricted to these five types: (i, j, A_ij, R_ij)
_TABLE_S1 = (
    (0, 0, 24.0, 0.955),
    (0, 1, 23.0, 1.0145),
    (0, 2, 28.5, 1.0355),
    (0, 3, 26.0, 0.9675),
    (0, 4, 45.0, 0.9775),
    (1, 1, 22.0, 1.074),
    (1, 2, 28.5, 1.095),
    (1, 3, 26.0, 1.027),
    (1, 4, 45.0, 1.037),
    (2, 2, 25.5, 1.116),
    (2, 3, 25.0, 1.048),
    (2, 4, 24.0, 1.058),
    (3, 3, 14.0, 0.95),
    (3, 4, 14.5, 0.99),
    (4, 4, 25.0, 1.0),
)
_SIGMA_ANGSTROM = 5.64                 # SI: one DPD distance unit
_C8E4_MASS = 306.0                    # g/mol (C8H17(OCH2CH2)4OH)
_WATER_MASS = 18.01528
_N_A = 6.02214076e23
_SIGMA3_L = (_SIGMA_ANGSTROM ** 3) * 1e-27   # sigma^3 in litres
_GAMMA_DP = 4.5                        # SI A.1
_SIGMA_DP = math.sqrt(2.0 * _GAMMA_DP)  # kT = 1
_KB_BOND = 150.0                       # SI A.2
_KA_ANGLE = 5.0                        # SI A.2 (Venturoli, theta0 = 180 deg)
_N_SUBMICELLAR = 6                     # paper 4.2: n* = 6 molecules
_R0_ANCHOR = 0.39                      # SI A.2: combined-heavy-atom anchor
_HEAVY = (1, 2, 3, 2)                  # heavy atoms per bead type
_SURFACTANT_TYPES = (0, 1, 1, 1, 2, 2, 2, 2, 3)  # C8E4, 9 beads
_N_PER = len(_SURFACTANT_TYPES)   # beads per surfactant molecule
_N_TAIL = 3                       # CH2CH2 tail beads per molecule
_N_BEAD_TYPES = 5


def _pair_tables(np):
    """A and R matrices (5x5 symmetric; unused pairs stay inf)."""
    a = np.full((_N_BEAD_TYPES, _N_BEAD_TYPES), np.inf)
    r = np.full((_N_BEAD_TYPES, _N_BEAD_TYPES), np.inf)
    for i, j, aij, rij in _TABLE_S1:
        a[i, j] = a[j, i] = aij
        r[i, j] = r[j, i] = rij
    return a, r


def _r0_for(type_a: int, type_b: int) -> float:
    """SI A.2: nominal bond length proportional to combined heavy atoms,
    anchored at 0.39 DPD units for combined = 4 (CH2CH2-CH2CH2)."""
    return _R0_ANCHOR * (_HEAVY[type_a] + _HEAVY[type_b]) / 4.0


def _box_for(n_surfactants: int, c_tot_M: float) -> float:
    """Exact box edge for n_surfactants at c_tot M (V_L = n/(c N_A))."""
    v_l = n_surfactants / (c_tot_M * _N_A)
    return (v_l / _SIGMA3_L) ** (1.0 / 3.0)


# ---------------------------------------------------------------------------
# Async job registry (in-memory + persisted JSON for restart honesty)
# ---------------------------------------------------------------------------

_JOBS: dict[str, dict] = {}
_JOBS_LOCK = threading.Lock()


def _jobs_dir(jobs_dir: str | None) -> Path:
    return Path(jobs_dir if jobs_dir is not None
                else os.environ.get("SURFACTANTKIT_DPD_JOBS", ".dpd_jobs"))


def _write_job_file(jobs_dir: Path, job_id: str, payload: dict) -> None:
    try:
        jobs_dir.mkdir(parents=True, exist_ok=True)
        (jobs_dir / f"{job_id}.json").write_text(
            json.dumps(payload, default=str))
    except OSError:
        pass  # in-memory state stays authoritative while the process lives


# ---------------------------------------------------------------------------
# Engine internals (numpy passed in; lazily imported by callers)
# ---------------------------------------------------------------------------

def _place_chains(n_surfactants: int, n_beads: int, l_box: float, rng, np):
    """Random-walk C10E6 chains at random COM positions (nominal r0 bond
    lengths) plus uniform solvent, wrapped into the box.  Returns
    (n_beads, 3) positions."""
    x = np.zeros((n_beads, 3))
    coms = rng.uniform(0.0, l_box, size=(n_surfactants, 3))
    for s in range(n_surfactants):
        base = s * _N_PER
        pos = np.zeros(3)
        x[base] = pos
        for k in range(1, _N_PER):
            dirv = rng.normal(size=3)
            dirv /= (np.linalg.norm(dirv) + 1e-300)
            pos = pos + dirv * _r0_for(_SURFACTANT_TYPES[k - 1],
                                       _SURFACTANT_TYPES[k])
            x[base + k] = pos
        x[base:base + _N_PER] += coms[s]
    if n_beads > n_surfactants * _N_PER:
        x[n_surfactants * _N_PER:] = rng.uniform(
            0.0, l_box, size=(n_beads - n_surfactants * _N_PER, 3))
    return np.mod(x, l_box)


def _bond_angle_arrays(n_surfactants: int, np):
    bonds = np.array([(s * _N_PER + k, s * _N_PER + k + 1)
                      for s in range(n_surfactants)
                      for k in range(_N_PER - 1)], dtype=np.int64)
    angles = np.array([(s * _N_PER + k, s * _N_PER + k + 1,
                        s * _N_PER + k + 2)
                       for s in range(n_surfactants)
                       for k in range(_N_PER - 2)], dtype=np.int64)
    return bonds, angles


def _conservative_forces(x, types, l_box, bonds, angles, a_mat, r_mat,
                         rows, cols, np):
    """DPD conservative (Table S1) + harmonic bonds + Venturoli angles.

    Pair terms use precomputed triu indices; every pair contributes
    antisymmetrically and each angle triplet assigns F_j = -(F_i + F_k),
    so total momentum is conserved exactly (up to float roundoff).
    """
    n = x.shape[0]
    dx = x[rows] - x[cols]
    dx -= l_box * np.round(dx / l_box)
    r = np.sqrt(np.einsum("ij,ij->i", dx, dx))
    a_vec = a_mat[types[rows], types[cols]]
    r_vec = r_mat[types[rows], types[cols]]
    inside = r < r_vec
    r_safe = np.where(inside, r, 1.0)
    coeff = np.where(inside, a_vec * (1.0 - r / r_safe), 0.0)
    f = coeff[:, None] * dx
    force = np.zeros_like(x)
    for c in range(3):
        force[:, c] = (np.bincount(rows, weights=f[:, c], minlength=n)
                       - np.bincount(cols, weights=f[:, c], minlength=n))

    # harmonic bonds (kb = 150; r0 by combined heavy atoms)
    if bonds.shape[0]:
        b0, b1 = bonds[:, 0], bonds[:, 1]
        dxb = x[b0] - x[b1]
        dxb -= l_box * np.round(dxb / l_box)
        rb = np.sqrt(np.einsum("ij,ij->i", dxb, dxb))
        r0 = np.array([_r0_for(int(types[i]), int(types[j]))
                       for i, j in zip(b0, b1)])
        rb_safe = np.maximum(rb, 1e-12)
        fb = (_KB_BOND * (rb - r0) / rb_safe)[:, None] * dxb
        for c in range(3):
            force[:, c] += (np.bincount(b0, weights=fb[:, c], minlength=n)
                            - np.bincount(b1, weights=fb[:, c], minlength=n))

    # Venturoli three-body angles (theta0 = 180 deg, ka = 5)
    if angles.shape[0]:
        i_idx, j_idx, k_idx = angles[:, 0], angles[:, 1], angles[:, 2]
        av = x[i_idx] - x[j_idx]
        av -= l_box * np.round(av / l_box)
        bv = x[k_idx] - x[j_idx]
        bv -= l_box * np.round(bv / l_box)
        na = np.sqrt(np.einsum("ij,ij->i", av, av))
        nb = np.sqrt(np.einsum("ij,ij->i", bv, bv))
        na_safe = np.maximum(na, 1e-12)
        nb_safe = np.maximum(nb, 1e-12)
        ua = av / na_safe[:, None]
        ub = bv / nb_safe[:, None]
        cos_t = np.clip(np.einsum("ij,ij->i", ua, ub), -1.0, 1.0)
        theta = np.arccos(cos_t)
        sin_t = np.maximum(np.sqrt(1.0 - cos_t * cos_t), 1e-12)
        # dtheta/da = (cos*ua - ub)/(|a| sin);  F_i = -ka(th-th0) dtheta/da
        fi = ((-_KA_ANGLE * (theta - math.pi) / (na_safe * sin_t))[:, None]
              * (cos_t[:, None] * ua - ub))
        fk = ((-_KA_ANGLE * (theta - math.pi) / (nb_safe * sin_t))[:, None]
              * (cos_t[:, None] * ub - ua))
        for c in range(3):
            force[:, c] += (np.bincount(i_idx, weights=fi[:, c],
                                        minlength=n)
                            + np.bincount(k_idx, weights=fk[:, c],
                                          minlength=n)
                            - np.bincount(j_idx, weights=fi[:, c] + fk[:, c],
                                          minlength=n))
    return force


def _thermostat_forces(x, v, l_box, r_cut, rng, dt, rows, cols, np):
    """Groot-Warren dissipative + random forces, (1 - r/rc)^2 form,
    rc = max R_ii (SI A.1).  Pair-symmetric => momentum conserving."""
    n = x.shape[0]
    dx = x[rows] - x[cols]
    dx -= l_box * np.round(dx / l_box)
    r = np.sqrt(np.einsum("ij,ij->i", dx, dx))
    dv = v[rows] - v[cols]
    w = np.where(r < r_cut, 1.0 - r / r_cut, 0.0)
    w2 = w * w
    r_safe = np.maximum(r, 1e-12)
    rhat = dx / r_safe[:, None]
    v_dot_r = np.einsum("ij,ij->i", dv, rhat)
    fd = (-_GAMMA_DP * w2 * v_dot_r)[:, None] * rhat
    noise = rng.normal(size=rows.shape[0])
    fr = ((_SIGMA_DP * w2 / math.sqrt(dt)) * noise)[:, None] * rhat
    f = fd + fr
    force = np.zeros_like(x)
    for c in range(3):
        force[:, c] = (np.bincount(rows, weights=f[:, c], minlength=n)
                       - np.bincount(cols, weights=f[:, c], minlength=n))
    return force


def _cluster_sizes(tx, l_box, mol_of, n_surfactants, contact, np):
    n_t = tx.shape[0]
    iu, ju = np.triu_indices(n_t, k=1)
    dx = tx[iu] - tx[ju]
    dx -= l_box * np.round(dx / l_box)
    r2 = np.einsum("ij,ij->i", dx, dx)
    parent = list(range(n_surfactants))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for a, b in zip(iu[r2 < contact * contact], ju[r2 < contact * contact]):
        ma, mb = int(mol_of[a]), int(mol_of[b])
        ra, rb = find(ma), find(mb)
        if ra != rb:
            parent[rb] = ra
    sizes: dict[int, int] = {}
    for m in range(n_surfactants):
        root = find(m)
        sizes[root] = sizes.get(root, 0) + 1
    return sorted(sizes.values())


# ---------------------------------------------------------------------------
# Runtime estimate + the worker
# ---------------------------------------------------------------------------

def estimate_dpd_runtime_s(steps: int, n_surfactants: int, c_tot_M: float,
                           n_runs: int) -> float:
    """Wall-clock estimate for n_runs DPD runs: measured conservative-force
    eval cost x 1.8 (thermostat overhead) x steps x runs x 2.5 safety."""
    try:
        import numpy as np
    except ImportError as exc:  # pragma: no cover
        raise ImportError("dpd_micelle_cmc needs numpy (pip install numpy)") \
            from exc
    l_box = _box_for(n_surfactants, c_tot_M)
    n_beads = max(round(3.0 * l_box ** 3), n_surfactants * _N_PER + 8)
    rng = np.random.default_rng(0)
    x = _place_chains(n_surfactants, n_beads, l_box, rng, np)
    types = np.array(list(_SURFACTANT_TYPES) * n_surfactants
                     + [4] * (n_beads - _N_PER * n_surfactants), dtype=np.int64)
    bonds, angles = _bond_angle_arrays(n_surfactants, np)
    a_mat, r_mat = _pair_tables(np)
    rows, cols = np.triu_indices(n_beads, k=1)
    t0 = time.perf_counter()
    for _ in range(3):
        _conservative_forces(x, types, l_box, bonds, angles, a_mat, r_mat,
                             rows, cols, np)
    per_eval = (time.perf_counter() - t0) / 3.0
    return per_eval * 1.8 * steps * n_runs * 2.5


def _run_one_dpd(n_surfactants: int, c_tot_M: float, steps: int, seed: int,
                 dt: float, sample_every: int, np):
    """One deterministic trajectory -> (times, p_free trace, meta)."""
    rng = np.random.default_rng(seed)
    l_box = _box_for(n_surfactants, c_tot_M)
    n_beads = max(round(3.0 * l_box ** 3), n_surfactants * _N_PER + 8)
    x = _place_chains(n_surfactants, n_beads, l_box, rng, np)
    types = np.array(list(_SURFACTANT_TYPES) * n_surfactants
                     + [4] * (n_beads - _N_PER * n_surfactants), dtype=np.int64)
    bonds, angles = _bond_angle_arrays(n_surfactants, np)
    a_mat, r_mat = _pair_tables(np)
    r_cut = float(np.max([v for _, _, _, v in _TABLE_S1]))
    rows, cols = np.triu_indices(n_beads, k=1)
    # tail (CH2CH2) beads for clustering: molecule s owns beads 9s+1..3
    tail_idx = np.array([s * _N_PER + 1 + k
                         for s in range(n_surfactants)
                         for k in range(_N_TAIL)], dtype=np.int64)
    mol_of = np.repeat(np.arange(n_surfactants), _N_TAIL)
    contact = 1.074                      # Table S1 CH2CH2-CH2CH2 cutoff

    v = np.zeros_like(x)
    times: list[float] = []
    p_free: list[float] = []
    for step in range(steps + 1):
        if step % sample_every == 0:
            sizes = _cluster_sizes(x[tail_idx], l_box, mol_of,
                                   n_surfactants, contact, np)
            n_free = sum(s_ for s_ in sizes if s_ <= _N_SUBMICELLAR)
            times.append(float(step * dt))
            p_free.append(n_free / n_surfactants)
        if step == steps:
            break
        f = _conservative_forces(x, types, l_box, bonds, angles,
                                 a_mat, r_mat, rows, cols, np)
        f += _thermostat_forces(x, v, l_box, r_cut, rng, dt, rows, cols, np)
        v = v + 0.5 * f * dt
        x = np.mod(x + v * dt, l_box)
        f = _conservative_forces(x, types, l_box, bonds, angles,
                                 a_mat, r_mat, rows, cols, np)
        f += _thermostat_forces(x, v, l_box, r_cut, rng, dt, rows, cols, np)
        v = v + 0.5 * f * dt
    meta = {
        "seed": int(seed),
        "n_beads": int(n_beads),
        "box_sigma": round(float(l_box), 3),
        "p_free_first": float(p_free[0]),
        "p_free_last": float(p_free[-1]),
    }
    return times, p_free, meta


def dpd_micelle_cmc(
    concentrations_M: list[float] | None = None,
    n_surfactants: int = 20,
    steps: int = 2000,
    n_seeds: int = 3,
    seed: int = 1,
    dt: float = 0.04,
    sample_every: int | None = None,
    max_sync_seconds: float = 60.0,
    allow_long_run: bool = False,
    jobs_dir: str | None = None,
    scale: str = "interactive",
) -> dict:
    """S4: DPD self-assembly CMC for nonionic C8E4 (the paper's SI model).

    Runs n_seeds independent short DPD trajectories per concentration,
    records p_free(t) (free-monomer ratio; aggregates <= n* = 6
    molecules are "free"), and feeds each concentration's traces to
    anderson_cmc_analysis (eq-8 fit + eq 3 + jackknife).  Args:
    concentrations_M (mol/L, default [0.5] -- interactive boxes are small
    and high-wt% by construction), n_surfactants >= 8, steps >= 200,
    n_seeds >= 3 (jackknife), seed, dt (0.04 conventional; SI does not
    state its own), scale='interactive' only ('paper' refuses: dl_meso
    20k-90k-time-unit batch class).  Refuses runs estimated above
    max_sync_seconds unless allow_long_run=True or submitted via
    submit_dpd_micelle_cmc_job.  Deterministic for a fixed seed.
    """
    try:
        import numpy as np
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "dpd_micelle_cmc needs numpy (pip install numpy -- the [calc] "
            "extra ships it)") from exc

    if scale != "interactive":
        raise ValueError(
            f"scale={scale!r} is refused: the production protocol (paper "
            "Table 1: dl_meso engine, L=34-54 boxes of 150-600 surfactants, "
            "60k-90k DPD time units, 10 seeds x 4 concentrations) is a "
            "batch/async-class job for the external dl_meso binary, not "
            "this in-house numpy worker -- run it as a batch job outside "
            "MCP")
    if concentrations_M is None:
        concentrations_M = [0.6]
    if not isinstance(concentrations_M, list) or not concentrations_M:
        raise ValueError("concentrations_M must be a non-empty list of mol/L")
    if any(not (0.0 < c <= 2.0) for c in concentrations_M):
        raise ValueError(
            "concentrations_M entries must be in (0, 2.0] mol/L -- "
            "interactive boxes cannot hold dilute systems (box volume is "
            "fixed by n_s/c)")
    if n_surfactants < 8:
        raise ValueError(f"n_surfactants must be >= 8 (got {n_surfactants}) "
                         "-- smaller boxes cannot form a micelle")
    if steps < 200:
        raise ValueError(f"steps must be >= 200 (got {steps}) -- shorter "
                         "runs cannot record an eq-8 decay window")
    if n_seeds < 3:
        raise ValueError(f"n_seeds must be >= 3 (got {n_seeds}) -- the "
                         "jackknife analysis needs >= 3 repeat traces")
    if not (0.001 <= dt <= 0.1):
        raise ValueError(f"dt must be in [0.001, 0.1] (got {dt}) -- velocity "
                         "Verlet DPD is unstable outside this window")

    n_runs = len(concentrations_M) * n_seeds
    est = estimate_dpd_runtime_s(steps, n_surfactants,
                                 float(concentrations_M[0]), n_runs)
    if est > max_sync_seconds and not allow_long_run:
        raise ValueError(
            f"estimated runtime ~{est:.0f}s exceeds max_sync_seconds="
            f"{max_sync_seconds:g} for {n_runs} runs -- either pass "
            "allow_long_run=True or submit via submit_dpd_micelle_cmc_job "
            "and poll dpd_micelle_cmc_job_status (this is the async path)")

    t_wall0 = time.perf_counter()
    sample_every = sample_every or max(steps // 40, 1)
    per_conc: list[dict] = []
    total_steps_done = 0
    n_beads0 = None

    for ci, c_tot in enumerate(concentrations_M):
        traces_conc: list[list[float]] = []
        times_grid: list[float] = []
        run_meta: list[dict] = []
        for si in range(n_seeds):
            run_seed = seed + 1000 * ci + si
            times_g, p_free_t, meta = _run_one_dpd(
                n_surfactants=n_surfactants, c_tot_M=float(c_tot),
                steps=steps, seed=run_seed, dt=dt,
                sample_every=sample_every, np=np)
            total_steps_done += steps
            if si == 0:
                times_grid = times_g
                n_beads0 = meta["n_beads"]
            traces_conc.append(p_free_t)
            run_meta.append(meta)
        analysis = anderson_cmc_analysis(
            times=times_grid, traces=traces_conc, c_tot=float(c_tot),
            beta=0.0, v_m_cm3_per_mol=_C8E4_MASS,  # density 1 g/cm3 (paper)
            j_ratio=0.9)
        ad = analysis.as_dict()
        per_conc.append({
            "c_tot_M": float(c_tot),
            "cmc_M": ad["cmc_M"],
            "cmc_jackknife_sigma_M": ad["cmc_jackknife_sigma_M"],
            "p_free_inf": ad["p_free_inf"],
            "tau_s": ad["tau_s"],
            "converged": ad["converged"],
            "n_subsets": ad["n_subsets"],
            "runs": run_meta,
            "analysis_params": ad["params"],
        })

    for p in per_conc:
        # finite-size quantization: p_free moves in steps of 1/n_s, so the
        # fitted A implies E[free molecules] = A * n_s.  Below ~5 molecules
        # the box cannot resolve the plateau and the CMC is not quotable.
        a_i = p["p_free_inf"]
        expectancy = (a_i * n_surfactants) if a_i is not None else None
        p["free_molecule_expectancy"] = expectancy
        p["finite_size_ok"] = bool(expectancy is not None and expectancy >= 5.0)
        if not p["finite_size_ok"]:
            p["finite_size_warning"] = (
                "fitted free-monomer expectation "
                f"{expectancy if expectancy is not None else float('nan'):.2f}"
                " molecules < 5: the box cannot resolve the p_free plateau "
                "(1/n_s quantization) -- enlarge the box (batch class) "
                "before quoting a CMC from this run")
    cmcs = [p["cmc_M"] for p in per_conc if p["cmc_M"] is not None]
    converged = bool(cmcs) and len(cmcs) == len(per_conc) and all(
        p["converged"] and p["finite_size_ok"] for p in per_conc)
    cmc_combined = (sum(cmcs) / len(cmcs)) if cmcs else None
    spread = None
    spread_log10 = None
    if len(cmcs) >= 2:
        mean_c = cmc_combined
        spread = math.sqrt(sum((c - mean_c) ** 2 for c in cmcs)
                           / (len(cmcs) - 1))
        if min(cmcs) > 0:
            spread_log10 = math.log10(max(cmcs) / min(cmcs))

    wall = time.perf_counter() - t_wall0
    caveat = (
        "In-house numpy DPD of the paper's OWN SI model (Table S1 "
        "A_ij/R_ij, gamma=4.5, k_b=150, Venturoli angles, r0 anchored at "
        "0.39 for CH2CH2-CH2CH2; Rummel et al. JPCB 2026 "
        "doi:10.1021/acs.jpcb.5c07475 Supporting Information) -- NOT "
        "dl_meso. Interactive scale: small N, high wt%, minutes per run "
        "vs the paper's L=34-54 / 60k-90k-time-unit production protocol "
        "(scale='paper' refuses); finite-size and short-window effects "
        "can shift p_free (the SI shows larger boxes fit better). C10E6 "
        "modeled NONIONIC (beta=0, exact eq 3) -- no SPME "
        "electrostatics; ionic surfactants are out of interactive scope. "
        "dt is exposed (the SI does not state its value). FINITE SIZE: "
        "p_free is quantized at 1/n_s -- quote a CMC only when "
        "per_concentration[].finite_size_ok (>= 5 expected free molecules); "
        "low-CMC systems like C10E6 need L~17+ boxes (batch class). "
        "Compare concentrations before quoting: per_concentration lists "
        "each box's estimate and the combined sigma/spread."
    )
    return {
        "cmc_M": cmc_combined,
        "cmc_sigma_M": spread,
        "cmc_spread_log10": spread_log10,
        "per_concentration": per_conc,
        "converged": converged,
        "n_iterations": total_steps_done,
        "method": "dpd_c8e4_si_model_anderson_pfree_eq3",
        "params": {
            "concentrations_M": [float(c) for c in concentrations_M],
            "n_surfactants": int(n_surfactants),
            "n_beads": int(n_beads0) if n_beads0 else None,
            "steps": int(steps),
            "n_seeds": int(n_seeds),
            "seed": int(seed),
            "dt": float(dt),
            "sample_every": int(sample_every),
            "n_submicellar_cutoff": _N_SUBMICELLAR,
            "wall_seconds": round(wall, 3),
            "estimated_runtime_s": round(est, 1),
            "gamma": _GAMMA_DP,
            "sigma_length_A": _SIGMA_ANGSTROM,
            "scale": "interactive",
            "async_hint": "submit_dpd_micelle_cmc_job for long runs",
        },
        "caveat": caveat,
    }


# ---------------------------------------------------------------------------
# Async submission API
# ---------------------------------------------------------------------------

def submit_dpd_micelle_cmc_job(
    concentrations_M: list[float] | None = None,
    n_surfactants: int = 20,
    steps: int = 2000,
    n_seeds: int = 3,
    seed: int = 1,
    dt: float = 0.04,
    jobs_dir: str | None = None,
) -> dict:
    """S4 async: start a dpd_micelle_cmc run in a background thread and
    return a job_id immediately; poll dpd_micelle_cmc_job_status(job_id).
    Job files persist under jobs_dir (default .dpd_jobs/).  A server
    restart mid-run leaves state='interrupted' (resubmit -- no silent
    resume, no fabricated progress)."""
    import uuid

    job_id = uuid.uuid4().hex[:12]
    payload = {
        "job_id": job_id, "state": "queued", "submitted_at": time.time(),
        "args": {"concentrations_M": concentrations_M,
                 "n_surfactants": n_surfactants, "steps": steps,
                 "n_seeds": n_seeds, "seed": seed, "dt": dt},
        "result": None, "error": None,
    }
    jdir = _jobs_dir(jobs_dir)
    with _JOBS_LOCK:
        _JOBS[job_id] = payload
    _write_job_file(jdir, job_id, payload)

    def _worker():
        with _JOBS_LOCK:
            _JOBS[job_id]["state"] = "running"
            _JOBS[job_id]["started_at"] = time.time()
        _write_job_file(jdir, job_id, dict(_JOBS[job_id]))
        try:
            result = dpd_micelle_cmc(
                concentrations_M=concentrations_M,
                n_surfactants=n_surfactants, steps=steps, n_seeds=n_seeds,
                seed=seed, dt=dt, allow_long_run=True, jobs_dir=jobs_dir)
            with _JOBS_LOCK:
                _JOBS[job_id]["state"] = "completed"
                _JOBS[job_id]["result"] = result
                _JOBS[job_id]["finished_at"] = time.time()
        except Exception as exc:  # noqa: BLE001 -- job boundary captures all
            with _JOBS_LOCK:
                _JOBS[job_id]["state"] = "failed"
                _JOBS[job_id]["error"] = f"{type(exc).__name__}: {exc}"
                _JOBS[job_id]["finished_at"] = time.time()
        _write_job_file(jdir, job_id, dict(_JOBS[job_id]))

    threading.Thread(target=_worker, daemon=True,
                     name=f"dpd-{job_id}").start()
    est = estimate_dpd_runtime_s(
        steps, n_surfactants,
        float(concentrations_M[0] if concentrations_M else 0.5),
        len(concentrations_M or [0.5]) * n_seeds)
    return {
        "job_id": job_id,
        "state": "queued",
        "estimated_runtime_s": round(est, 1),
        "poll": f"dpd_micelle_cmc_job_status(job_id={job_id!r})",
        "jobs_dir": str(jdir),
        "caveat": "Background thread in this server process: a server "
                  "restart interrupts the run (state becomes 'interrupted' "
                  "-- resubmit). The result JSON is persisted to jobs_dir "
                  "on completion.",
    }


def dpd_micelle_cmc_job_status(job_id: str,
                               jobs_dir: str | None = None) -> dict:
    """S4 async: status + result of a submitted DPD job.  States: queued,
    running, completed (full dpd_micelle_cmc result attached), failed
    (error string attached), interrupted (persisted file but no live
    thread -- the server restarted mid-run)."""
    if not job_id or not isinstance(job_id, str):
        raise ValueError(f"job_id must be a non-empty string "
                         f"(got {job_id!r})")
    # ids are uuid4().hex[:12]; anything else could address a file outside
    # the jobs directory
    if len(job_id) != 12 or any(ch not in "0123456789abcdef" for ch in job_id):
        raise ValueError(f"job_id must be the 12-character hex id returned "
                         f"by submit_dpd_micelle_cmc_job (got {job_id!r})")
    with _JOBS_LOCK:
        rec = dict(_JOBS.get(job_id) or {})
    if not rec:
        jdir = _jobs_dir(jobs_dir)
        f = jdir / f"{job_id}.json"
        if f.exists():
            try:
                rec = json.loads(f.read_text(encoding="utf-8"))
                rec["state"] = "interrupted"
                rec["note"] = ("persisted progress found but no live thread "
                               "(server restarted) -- resubmit")
            except (OSError, json.JSONDecodeError) as exc:
                raise ValueError(f"job file for job_id={job_id!r} "
                                 f"unreadable: {exc}") from exc
        else:
            raise ValueError(f"unknown job_id {job_id!r} -- no in-memory "
                             f"job and no file in {jdir}")
    return {
        "job_id": job_id,
        "state": rec.get("state"),
        "result": rec.get("result"),
        "error": rec.get("error"),
        "submitted_at": rec.get("submitted_at"),
        "finished_at": rec.get("finished_at"),
        "note": rec.get("note"),
        "caveat": "completed results carry the full dpd_micelle_cmc caveat "
                  "(in-house SI-model DPD, interactive scale, not dl_meso).",
    }
