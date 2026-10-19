"""S3: coarse-grained Langevin bead-spring self-assembly of surfactants.

Each surfactant molecule = 2 beads (hydrophilic head + hydrophobic tail)
joined by a harmonic spring.  Overdamped Brownian (Langevin) dynamics in a
periodic cube; tail-tail attraction drives micelle formation; cluster
analysis over production frames yields the aggregation-number distribution
and a monomer-fraction based simulated CMC estimate.

This is a QUALITATIVE CG model: reduced units with an explicit time
mapping, no explicit solvent, no electrostatics.  It reproduces trends
(micelles above CMC, monomer plateau at CMC, Nagg growing with
concentration) -- it does not predict absolute CMCs for real surfactants.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class BrownianAssemblyResult:
    """Cluster statistics from a CG Langevin self-assembly run."""
    concentration_M: float
    n_molecules: int
    mean_nagg: float                # micellar (n>=5) mean aggregation number
    largest_cluster_frac: float     # largest cluster / N over frames
    monomer_fraction: float         # molecules alone, time-averaged
    nagg_mean_overall: float        # all molecules incl. monomers
    nagg_p95: float
    micelle_count_mean: float
    n_frames: int
    time_tau: float                 # total simulated time in tau = gamma sigma^2 / kT
    nagg_distribution: dict[str, float]   # cluster size (str for JSON) -> fraction
                                          # of clusters over production frames
    nagg_timeseries: list[float]    # per-frame mean cluster size (stationarity trace)
    converged: bool
    n_iterations: int
    method: str
    params: dict = field(default_factory=dict)
    caveat: str = ""

    def as_dict(self) -> dict:
        return {k: getattr(self, k) for k in
                ("concentration_M", "n_molecules", "mean_nagg",
                 "largest_cluster_frac", "monomer_fraction",
                 "nagg_mean_overall", "nagg_p95", "micelle_count_mean",
                 "n_frames", "time_tau", "nagg_distribution",
                 "nagg_timeseries", "converged", "n_iterations",
                 "method", "params", "caveat")}


def brownian_self_assembly(
    concentration_M: float = 0.01,
    n_molecules: int = 200,
    steps: int = 12000,
    seed: int = 1,
    tail_attraction_kT: float = 3.0,
    spring_kT_per_sigma2: float = 80.0,
    dt_tau: float = 0.002,
    molecular_volume_nm3: float = 0.7,
    sigma_nm: float = 0.5,
    temperature_K: float = 298.15,
    frame_every: int = 100,
    cluster_cutoff_sigma: float = 1.5,
    init: str = "random",
    init_relax_steps: int = 400,
) -> BrownianAssemblyResult:
    """Simulate CG surfactant self-assembly; return Nagg statistics.

    Parameters
    ----------
    concentration_M : total surfactant molarity (sets the box edge from
        molecular_volume_nm3: L = (N * V_mol / conc / N_A)^(1/3)).
    n_molecules : number of 2-bead surfactant molecules (>= 20 for
        meaningful cluster statistics).
    steps : overdamped Langevin steps (>= 500; assembly needs thousands).
    seed : RNG seed -- SAME seed gives bit-identical results (asserted in
        tests), different seed a genuinely different trajectory.
    tail_attraction_kT : depth of the tail-tail LJ attraction in kT
        (0 = ideal-ish gas, no micelles; ~1-1.5 gives micelles above CMC).
    spring_kT_per_sigma2 : head-tail harmonic spring constant k/(kT/sigma^2).
    dt_tau : time step in tau = gamma*sigma^2/kT units.
    molecular_volume_nm3 : excluded volume per molecule (packing).
    cluster_cutoff_sigma : tail-tail contact distance for cluster graph.

    Returns
    -------
    BrownianAssemblyResult with converged (stationarity of mean Nagg over
    frame quintiles), n_iterations (=steps), method, params, caveat.

    Notes
    -----
    - Deterministic given seed (Box-Muller via numpy Generator).
    - Force magnitude capped at 200 kT/sigma (WCA stiffness guard; the
      cap rate is recorded in params.force_cap_fraction).
    - Real-time mapping: tau = gamma sigma^2 / kT ~ 1-10 ns for
      aqueous surfactants (sigma=0.5 nm, gamma from Stokes-Einstein);
      treat kinetics as qualitative.
    """
    try:
        import numpy as np
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "brownian_self_assembly needs numpy (core dependency -- "
            "pip install numpy); not the [calc] extra") from exc

    import math

    if n_molecules < 20:
        raise ValueError("n_molecules must be >= 20 -- smaller boxes cannot "
                         "form a micelle and the Nagg statistics would be "
                         "meaningless (do not report cluster numbers from "
                         "a handful of molecules)")
    if steps < 500:
        raise ValueError("steps must be >= 500 -- shorter runs never "
                         "sample even single-molecule diffusion across "
                         "the box, let alone aggregation")
    if concentration_M <= 0:
        raise ValueError("concentration_M must be positive -- it sets the "
                         "box volume; pass the formulation concentration")
    if not (0.0 <= tail_attraction_kT):
        raise ValueError("tail_attraction_kT must be >= 0 (units of kT)")
    if spring_kT_per_sigma2 <= 0 or dt_tau <= 0:
        raise ValueError("spring_kT_per_sigma2 and dt_tau must be positive")
    if frame_every < 1:
        raise ValueError("frame_every must be >= 1")
    if init_relax_steps < 0:
        raise ValueError("init_relax_steps must be >= 0 (0 only if you "
                         "are certain the start config has no overlaps)")
    if init not in ("random", "clustered"):
        raise ValueError("init must be 'random' or 'clustered' -- "
                         "'clustered' seeds a droplet to avoid the slow "
                         "nucleation barrier; 'random' is unbiased but "
                         "needs many more steps to micellize")

    rng = np.random.default_rng(int(seed))

    # ---- box edge from concentration: L^3 = N * V_mol ; V_mol = V_mol_nm3
    # conc in mol/L = N / (N_A * V_L); V_nm3 = V_L * 1e24
    n_a = 6.02214076e23
    volume_L = n_molecules / (concentration_M * n_a)          # liters
    box_nm = (volume_L * 1e24) ** (1.0 / 3.0)   # 1 L = 1e24 nm^3
    # concentration already fixes the total volume; molecular_volume only
    # informs the packing fraction reported in params.
    packing = n_molecules * molecular_volume_nm3 / (box_nm ** 3)

    # positions: heads and tails paired at distance r0 = sigma
    r0 = 1.0
    if init == "clustered":
        # start all molecules in one loose micelle-like droplet: drastically
        # shortens the (slow) nucleation time; dynamics still equilibrates
        # the droplet and it may dissolve below the CMC -- which is a feature
        # (low-conc runs with this init measure droplet dissolution).
        vecs = rng.normal(size=(n_molecules, 3))
        vecs /= np.linalg.norm(vecs, axis=1, keepdims=True)
        # radius so that n spheres of DIAMETER sigma reach fill ~0.45:
        # V = n*(pi/6)*sigma^3 / 0.45 = (4/3) pi r^3  ->  r = (n/3.6)^(1/3) sigma
        core_r = (n_molecules / 3.6) ** (1.0 / 3.0) * sigma_nm
        rad = core_r * (rng.uniform(0.15, 1.0, size=(n_molecules, 1))
                        ** (1.0 / 3.0))
        c = box_nm / 2.0
        pos = np.empty((n_molecules, 2, 3))
        pos[:, 1, :] = c + vecs * rad                      # tails in core
        pos[:, 0, :] = c + vecs * (rad + r0)               # heads outward
        pos %= box_nm
    elif init == "random":
        pos = rng.uniform(0.0, box_nm, size=(n_molecules, 2, 3))
        directions = rng.normal(size=(n_molecules, 3))
        directions /= np.linalg.norm(directions, axis=1, keepdims=True)
        pos[:, 1, :] = pos[:, 0, :] + directions * r0
        pos %= box_nm
    else:
        raise ValueError("init must be 'random' or 'clustered'")

    sigma = sigma_nm
    # WCA (all pairs) + LJ tail-tail
    eps_tail = tail_attraction_kT           # kT units
    k_spring = spring_kT_per_sigma2         # kT/sigma^2
    gamma = 1.0                             # friction, reduced
    dt = dt_tau
    noise_amp = math.sqrt(2.0 * dt / (gamma ** 2))   # kT=1 units, overdamped
    # overdamped stability: max displacement per step from a capped force
    # must stay << sigma (force_cap*dt <= 0.05 nm = 0.1*sigma), otherwise
    # overlapping beads teleport across the box and dense states explode.
    force_cap = 25.0
    cut = 2.5
    cut2 = cut * cut

    n_beads = 2 * n_molecules
    bead_type = np.zeros(n_beads, dtype=np.int8)   # 0 head, 1 tail
    bead_type[1::2] = 1
    bonded_partner = np.arange(n_beads) ^ 1

    def forces(p):
        flat = p.reshape(-1, 3)
        # --- bonded spring (head-tail within molecule)
        d = flat[1::2] - flat[0::2]
        d -= box_nm * np.round(d / box_nm)
        r = np.linalg.norm(d, axis=1)
        r_safe = np.maximum(r, 1e-9)
        # U = k/2 (r-r0)^2 ; d = r_tail - r_head ; F_head = +k(r-r0) * d/r
        f_spring = k_spring * (r_safe - r0)[:, None] * d / r_safe[:, None]
        # Newton pairs
        f_pair = np.zeros_like(flat)
        f_pair[0::2] += f_spring
        f_pair[1::2] -= f_spring

        # --- nonbonded all-pairs (O(N^2), N beads is small in tests)
        diff = flat[:, None, :] - flat[None, :, :]
        diff -= box_nm * np.round(diff / box_nm)
        d2 = np.einsum('ijk,ijk->ij', diff, diff)
        np.fill_diagonal(d2, np.inf)
        # exclude bonded pairs
        idx = np.arange(n_beads)
        d2[idx, bonded_partner] = np.inf
        mask = d2 < cut2
        inv_r2 = np.where(mask, 1.0 / np.maximum(d2, 1e-12), 0.0)
        s6 = (sigma * sigma) * inv_r2
        s6 = s6 * s6 * s6
        # LJ radial coefficient, positive = repulsive: 48 e (s12 - s6/2)/r^2
        lj_coeff = 48.0 * (s6 * s6 - 0.5 * s6) * inv_r2
        rmin2 = sigma * sigma * 2.0 ** (1.0 / 3.0)
        # WCA core for every nonbonded pair (r < rmin)
        wca_coeff = np.where((d2 < rmin2) & mask, lj_coeff, 0.0)
        # tail-tail: full LJ of depth eps_tail for rmin <= r < cut
        # (at eps=1 this reduces to the LJ continuation of WCA)
        tt = (bead_type[:, None] == 1) & (bead_type[None, :] == 1)
        tail_coeff = np.where(tt & mask & (d2 >= rmin2),
                              eps_tail * lj_coeff, 0.0)
        coeff = wca_coeff + tail_coeff
        # f_i += coeff * diff_ij ; f_j -= ...
        f_pair += np.einsum('ij,ijk->ik', coeff, diff)

        # cap
        norm = np.linalg.norm(f_pair, axis=1)
        too_big = norm > force_cap
        if np.any(too_big):
            f_pair[too_big] *= (force_cap / norm[too_big])[:, None]
        return f_pair.reshape(p.shape), float(np.mean(too_big))

    # ---- deterministic overlap relaxation before dynamics (no noise):
    # random/clustered starts can place two beads on top of each other;
    # the raw cap-limited forces would then blast the config apart in the
    # first step.  Steepest-descent with a small dt resolves overlaps.
    # bounded steepest descent: clamp the per-step displacement so the
    # stiff WCA wall cannot make the integrator oscillate (a bare SD step
    # with dt*F' >> 1 at r << sigma never resolves overlaps).
    max_disp = 0.05 * sigma_nm
    for _ in range(init_relax_steps):
        f, _cap = forces(pos)
        dx = (f / gamma) * (dt / 10.0)
        norms = np.maximum(np.linalg.norm(dx, axis=2, keepdims=True), 1e-300)
        dx = np.where(norms > max_disp, dx * (max_disp / norms), dx)
        pos += dx
        pos %= box_nm

    # ---- production
    half = steps // 2
    frames = []
    cap_samples = []
    for step in range(steps):
        f, cap_frac = forces(pos)
        if step >= half and step % frame_every == 0:   # production only
            cap_samples.append(cap_frac)
        drift = (f / gamma) * dt
        pos += drift + noise_amp * rng.normal(size=pos.shape)
        pos %= box_nm
        if step >= half and (step - half) % frame_every == 0:
            frames.append(pos.copy())

    if len(frames) < 4:
        raise RuntimeError(
            f"not enough production frames ({len(frames)}) -- increase "
            "steps or lower frame_every so at least 4 frames are sampled "
            "after the equilibration half")

    # ---- cluster analysis on tail-tail contact graph (+ bonded)
    def clusters_of(p):
        flat = p.reshape(-1, 3)
        tails = flat[1::2]
        diff = tails[:, None, :] - tails[None, :, :]
        diff -= box_nm * np.round(diff / box_nm)
        d2 = np.einsum('ijk,ijk->ij', diff, diff)
        adj = d2 < (cluster_cutoff_sigma * sigma) ** 2
        np.fill_diagonal(adj, True)
        # BFS
        seen = np.zeros(n_molecules, dtype=bool)
        sizes = []
        for start in range(n_molecules):
            if seen[start]:
                continue
            stack = [start]
            seen[start] = True
            count = 0
            while stack:
                node = stack.pop()
                count += 1
                nbrs = np.nonzero(adj[node] & ~seen)[0]
                for nb in nbrs:
                    seen[nb] = True
                    stack.append(int(nb))
            sizes.append(count)
        return np.array(sizes, dtype=float)

    per_frame_nagg = []
    per_frame_largest = []
    per_frame_monomer = []
    per_frame_micelles = []
    for p in frames:
        sizes = clusters_of(p)
        per_frame_nagg.append(sizes)
        per_frame_largest.append(sizes.max() / n_molecules)
        per_frame_monomer.append(float(np.sum(sizes == 1)) / n_molecules)
        per_frame_micelles.append(float(np.sum(sizes >= 5)))

    micellar_sizes = np.concatenate([s[s >= 5] for s in per_frame_nagg]) \
        if any(np.any(s >= 5) for s in per_frame_nagg) else np.array([])
    mean_nagg = float(np.mean(micellar_sizes)) if micellar_sizes.size else 1.0
    all_sizes = np.concatenate(per_frame_nagg)
    nagg_overall = float(np.mean(all_sizes))
    nagg_p95 = float(np.percentile(all_sizes, 95))

    largest = float(np.mean(per_frame_largest))
    monomer = float(np.mean(per_frame_monomer))
    micelles_mean = float(np.mean(per_frame_micelles))

    # aggregation-number distribution: fraction of clusters of each size
    # over all production frames (the module docstring and MCP tool
    # description both promise this; previously only mean/p95 shipped)
    size_counts: dict[int, int] = {}
    n_clusters_total = 0
    for s in per_frame_nagg:
        for sz in s:
            k = int(sz)
            size_counts[k] = size_counts.get(k, 0) + 1
            n_clusters_total += 1
    nagg_distribution = {
        str(k): size_counts[k] / n_clusters_total
        for k in sorted(size_counts)
    } if n_clusters_total else {}

    # stationarity: quintile drift of micelle-weighted mean cluster size
    frame_means = [float(np.mean(s)) for s in per_frame_nagg]
    n_fr = len(frame_means)
    q = max(1, n_fr // 5)
    first = float(np.mean(frame_means[:q]))
    last = float(np.mean(frame_means[-q:]))
    drift = abs(last - first) / max(first, 1e-12)
    converged = bool(drift < 0.35)      # loose: finite CG runs drift

    monomer_conc = monomer * concentration_M
    cmc_proxy_valid = bool(micelles_mean >= 1.0 and monomer > 0.02)

    return BrownianAssemblyResult(
        concentration_M=float(concentration_M),
        n_molecules=int(n_molecules),
        mean_nagg=mean_nagg,
        largest_cluster_frac=largest,
        monomer_fraction=monomer,
        nagg_mean_overall=nagg_overall,
        nagg_p95=nagg_p95,
        micelle_count_mean=micelles_mean,
        n_frames=len(frames),
        time_tau=float(steps * dt_tau),
        nagg_distribution=nagg_distribution,
        nagg_timeseries=[float(x) for x in frame_means],
        converged=converged,
        n_iterations=int(steps),
        method="overdamped_langevin_2bead_spring_cellfree_n2",
        params={
            "seed": int(seed),
            "tail_attraction_kT": float(tail_attraction_kT),
            "spring_kT_per_sigma2": float(spring_kT_per_sigma2),
            "dt_tau": float(dt_tau),
            "box_nm": float(box_nm),
            "packing_fraction": float(packing),
            "sigma_nm": float(sigma_nm),
            "temperature_K": float(temperature_K),
            "frame_every": int(frame_every),
            "cluster_cutoff_sigma": float(cluster_cutoff_sigma),
            "init": str(init),
            "init_relax_steps": int(init_relax_steps),
            "force_cap_fraction_mean": float(np.mean(cap_samples)),
            "stationarity_drift": float(drift),
            "monomer_concentration_M": float(monomer_conc),
            "cmc_proxy_M": float(monomer_conc) if cmc_proxy_valid else None,
            "cmc_proxy_valid": cmc_proxy_valid,
        },
        caveat=(
            "Coarse-grained 2-bead surfactant, implicit solvent, no "
            "electrostatics: qualitative trends only (micellization above "
            "a CMC-like threshold, Nagg growing with concentration). The "
            "cmc_proxy_M is the measured MONOMER concentration in the "
            "production frames and equals the CMC only when micelles "
            "coexist with monomers (cmc_proxy_valid=True); a real CMC "
            "requires a concentration scan (brownian_cmc_scan) and "
            "extrapolation to the monomer plateau. Time is in reduced "
            "tau = gamma*sigma^2/kT; treat kinetics as order-of-magnitude. "
            "tail_attraction_kT sets the effective critical packing: "
            "quote it with every number. Force cap fraction "
            "(params) should stay small -- if >0.01 the dt/ stiffness "
            "settings are too aggressive."
        ),
    )


def brownian_cmc_scan(
    concentrations_M: list[float],
    n_molecules: int = 120,
    steps: int = 8000,
    seed: int = 1,
    tail_attraction_kT: float = 3.0,
    **kwargs,
) -> dict:
    """Run brownian_self_assembly over a concentration series; estimate CMC
    as the concentration where the micelle fraction crosses 50% of its
    high-conc value (monotone-crossing interpolation).

    Returns dict with per-concentration table, cmc_estimate_M (None if no
    crossing in the scanned window), converged (all runs converged), and
    the caveat that a bracketed scan should be tightened around the
    crossing before quoting a CMC.
    """
    if len(concentrations_M) < 2:
        raise ValueError("concentrations_M needs at least 2 points to scan "
                         "-- a single concentration cannot locate a CMC")
    if any(c <= 0 for c in concentrations_M):
        raise ValueError("all concentrations must be positive")
    rows = []
    for i, c in enumerate(sorted(float(x) for x in concentrations_M)):
        r = brownian_self_assembly(
            concentration_M=c, n_molecules=n_molecules, steps=steps,
            seed=seed + i, tail_attraction_kT=tail_attraction_kT, **kwargs)
        rows.append({"concentration_M": c,
                     "micelle_fraction": r.micelle_count_mean / n_molecules,
                     "monomer_fraction": r.monomer_fraction,
                     "largest_cluster_frac": r.largest_cluster_frac,
                     "mean_nagg": r.mean_nagg,
                     "converged": r.converged})
    # 50% crossing of micelle fraction (normalized to max observed)
    fracs = [r["micelle_fraction"] for r in rows]
    concs = [r["concentration_M"] for r in rows]
    top = max(fracs)
    cmc = None
    if top > 0:
        half = 0.5 * top
        for i in range(len(fracs) - 1):
            a, b = fracs[i], fracs[i + 1]
            if a < half <= b and b > a:
                t = (half - a) / (b - a)
                cmc = concs[i] * (concs[i + 1] / concs[i]) ** t
                break
    return {
        "table": rows,
        "cmc_estimate_M": (float(cmc) if cmc is not None else None),
        "converged": bool(all(r["converged"] for r in rows)),
        "method": "micelle_fraction_half_max_log_interpolation",
        "n_iterations": len(rows),
        "caveat": "Half-max crossing on a coarse scan: bracket the "
                  "crossing with 2-3 tighter concentrations before "
                  "quoting a CMC, and remember this is a reduced-unit "
                  "CG model (see brownian_self_assembly caveat).",
    }
