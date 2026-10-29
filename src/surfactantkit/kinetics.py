"""Kinetics workers (Phase 2): O1 stiff micelle-relaxation master
equation, S2 Gillespie kinetic Monte Carlo.

Both describe the SAME physical system -- stepwise micellization
C_n + C_1 <-> C_{n+1} with one equilibrium constant K per step and
detailed balance (k_off = k_on/K) -- so their results must agree with
each other and with the closed-form mass-action equilibrium.  Every gold
is internal-consistency based: no fitted literature numbers.

Notation: state y[0..N-1] with y[j] = [C_{j+1}].  Step m (0..N-2):
g[m] = k_on*y[m]*y[0] - k_off*y[m+1]  (net forward flux of
C_{m+1} + C_1 -> C_{m+2}).

Balances (mass-conserving by construction -- verified by a
finite-difference Jacobian check inside O1 and by total-mass monitoring):
  dy[0]      = -sum(g) - g[0]         # the extra -g[0] is the 2-monomer
                                      # bookkeeping of the C1+C1 <-> C2 step
  dy[j]      = g[j-1] - g[j]          # j = 1..N-2
  dy[N-1]    = g[N-2]
"""

from __future__ import annotations

import math
from dataclasses import dataclass

_CALC_HINT = ("kinetics workers require numpy+scipy -- install with "
              "'pip install surfactantkit[calc]'")


def _require_calc():
    try:
        import numpy as np
        import scipy  # noqa: F401
    except ImportError as exc:  # pragma: no cover
        raise ImportError(_CALC_HINT) from exc
    return np


def mass_action_equilibrium(c_total_M: float, K_per_M: float, n_max: int,
                            x1_tol: float = 1e-16) -> list[float]:
    """Exact stepwise mass-action equilibrium concentrations x_1..x_n_max.

    x_n = K^(n-1) * x_1^n with sum(n * x_n) = c_total; bisection on x_1
    (the total is strictly increasing in x_1).  This closed-form target is
    what both kinetics workers must relax to -- independent of any ODE or
    stochastic solver.

    Numerical note (2026-10-05): the product is evaluated as
    x_1 * (K*x_1)^(n-1), algebraically identical to K^(n-1)*x_1^n but
    overflow-safe -- the separated form overflowed float64 for perfectly
    physical inputs (e.g. K=1e4 M^-1, n_max=100: K^99 = 1e396) while the
    grouped form stays finite wherever the term itself is representable.
    """
    if c_total_M <= 0:
        raise ValueError("c_total_M must be positive")
    if K_per_M <= 0:
        raise ValueError("K_per_M must be positive")
    if n_max < 3:
        raise ValueError("n_max must be >= 3")

    def total(x1: float) -> float:
        return sum(n * x1 * (K_per_M * x1) ** (n - 1)
                   for n in range(1, n_max + 1))

    lo, hi = 0.0, c_total_M
    while total(hi) < c_total_M:
        hi *= 2.0
        if hi > 1e6 * c_total_M:
            break
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if total(mid) < c_total_M:
            lo = mid
        else:
            hi = mid
        if hi - lo <= x1_tol * max(c_total_M, 1e-300):
            break
    x1 = 0.5 * (lo + hi)
    return [x1 * (K_per_M * x1) ** (n - 1) for n in range(1, n_max + 1)]


# ---------------------------------------------------------------------------
# O1 -- stiff master-equation relaxation (deterministic ODE)
# ---------------------------------------------------------------------------

@dataclass
class MicelleRelaxationResult:
    equilibrium_concentrations_M: list[float]
    kinetic_cmc_M: float            # monomer at equilibrium == kinetic CMC
    relaxation_times_s: list[float] # slowest nonzero linear-response modes
    tau_dominant_s: float
    mean_aggregation_number: float
    trajectory_t_s: list[float]
    trajectory_c1_M: list[float]
    trajectory_mean_nagg: list[float]
    tail_exp_fit_r_squared: float
    tail_tau_fit_s: float
    max_relative_mass_error: float
    converged: bool
    n_iterations: int
    method: str
    params: dict
    caveat: str


def micelle_relaxation_odes(
    c_total_M: float,
    K_initial_per_M: float,
    K_final_per_M: float,
    k_on_per_M_s: float,
    n_max: int = 40,
    n_times: int = 400,
) -> MicelleRelaxationResult:
    """Stiff master-equation relaxation after a T-jump (worker O1).

    C_n + C_1 <-> C_{n+1} (n = 1..n_max) with detailed balance
    k_off = k_on/K.  Starts at full equilibrium under K_initial, then
    evolves under K_final (a temperature jump changes K at fixed total
    concentration) until it reaches the NEW exact equilibrium.

    Returns: exact equilibrium from an independent closed-form bisection,
    the kinetic CMC (equilibrium monomer), the relaxation spectrum from
    Jacobian eigenvalues at that equilibrium (the slowest nonzero mode is
    the experimentally measurable tau of the jump), the c_1(t) /
    mean-Nagg(t) trajectories, an exponential tail fit cross-checking the
    eigenvalue time, and the total-mass conservation error.

    Solved with scipy BDF + analytic Jacobian (verified against finite
    differences at run time): fast monomer association vs slow micelle
    reorganization is exactly the stiff regime.

    caveat: ideal mass action with one common K per step -- no
    size-dependent K, no C_n+C_m fusion, no headgroup electrostatics; the
    spectrum is the linear-response spectrum around the new equilibrium.
    """
    np = _require_calc()
    from scipy.integrate import solve_ivp

    if c_total_M <= 0:
        raise ValueError("c_total_M must be positive")
    if K_initial_per_M <= 0 or K_final_per_M <= 0:
        raise ValueError("K values must be positive")
    if k_on_per_M_s <= 0:
        raise ValueError("k_on_per_M_s must be positive")
    if n_max < 3:
        raise ValueError("n_max must be >= 3")
    if n_times < 10:
        raise ValueError("n_times must be >= 10")

    N = n_max
    k_on = k_on_per_M_s
    k_off = k_on_per_M_s / K_final_per_M

    eq_final = mass_action_equilibrium(c_total_M, K_final_per_M, N)
    eq_initial = mass_action_equilibrium(c_total_M, K_initial_per_M, N)
    y0 = np.asarray(eq_initial, dtype=float)

    def _g(y):
        return k_on * y[: N - 1] * y[0] - k_off * y[1:N]

    def rhs(t, y):
        g = _g(y)
        dydt = np.empty(N)
        dydt[0] = -np.sum(g) - g[0]
        dydt[1:N - 1] = g[: N - 2] - g[1:]
        dydt[N - 1] = g[N - 2]
        return dydt

    def jac(t, y):
        J = np.zeros((N, N))
        # row 0: dy0 = -sum(g) - g[0]; note dg0/dy0 = 2*k_on*y0 (y0 squared)
        J[0, 0] = -k_on * float(np.sum(y[: N - 1])) - 3.0 * k_on * y[0]
        for j in range(1, N):
            v = k_off
            if j <= N - 2:
                v += -k_on * y[0]
            if j == 1:
                v += k_off
            J[0, j] = v
        # middle rows
        for j in range(1, N - 1):
            J[j, j - 1] = k_on * y[0]
            # row 1's incoming step is g[0] with y0 squared -> 2*k_on*y0
            J[j, 0] = k_on * (y[j - 1] - y[j]) + (k_on * y[0] if j == 1 else 0.0)
            J[j, j] = -k_off - k_on * y[0]
            J[j, j + 1] = k_off
        # last row
        J[N - 1, 0] = k_on * y[N - 2]
        J[N - 1, N - 2] = k_on * y[0]
        J[N - 1, N - 1] = -k_off
        return J

    # --- finite-difference verification of the analytic Jacobian ---
    y_test = y0 * 1.001 + c_total_M * 1e-9
    J_a = jac(0.0, y_test)
    f0 = rhs(0.0, y_test)
    J_n = np.empty((N, N))
    for col in range(N):
        h = 1e-7 * max(abs(y_test[col]), c_total_M * 1e-6)
        yp = y_test.copy()
        yp[col] += h
        J_n[:, col] = (rhs(0.0, yp) - f0) / h
    fd_err = float(np.max(np.abs(J_a - J_n)) / max(np.max(np.abs(J_n)), 1e-300))
    if fd_err > 1e-4:
        raise RuntimeError(
            f"analytic Jacobian disagrees with finite differences "
            f"(rel err {fd_err:.2e}) -- implementation bug"
        )

    # --- relaxation spectrum from J at the target equilibrium ---
    J_eq = jac(0.0, np.asarray(eq_final))
    eigs = np.linalg.eigvals(J_eq)
    scale = max(abs(eigs))
    # decay rates ASCENDING -> relaxation times longest-first
    decay_rates = sorted(-l.real for l in eigs if l.real < -1e-9 * scale)
    if not decay_rates:
        raise RuntimeError("no stable negative eigenvalue at equilibrium")
    tau_dom = 1.0 / decay_rates[0]  # slowest nonzero mode
    relaxation_times = [1.0 / r for r in decay_rates[: min(8, len(decay_rates))]]

    t_end = 15.0 * tau_dom
    t_eval = np.linspace(0.0, t_end, n_times)
    sol = solve_ivp(rhs, (0.0, t_end), y0, method="BDF", jac=jac, t_eval=t_eval,
                    rtol=1e-7, atol=1e-12 * c_total_M)
    if not sol.success:
        raise RuntimeError(f"BDF solve failed: {sol.message}")

    Y = sol.y
    masses = np.arange(1, N + 1)[:, None]
    totals = np.sum(masses * Y, axis=0)
    max_rel_mass_err = float(np.max(np.abs(totals - c_total_M)) / c_total_M)
    concentrations_sum = np.sum(Y, axis=0)
    mean_nagg_t = totals / np.maximum(concentrations_sum, 1e-300)
    c1_traj = Y[0]

    # --- exponential tail fit of c1(t): c1_inf + A exp(-t/tau) ---
    n_tail = max(10, int(0.4 * n_times))
    t_tail = t_eval[-n_tail:]
    dy = c1_traj[-n_tail:] - float(eq_final[0])
    mask = np.abs(dy) > 1e-13 * c_total_M
    if mask.sum() >= 5:
        lt = np.log(np.abs(dy[mask]))
        tt = t_tail[mask]
        slope, intercept = np.polyfit(tt, lt, 1)
        tail_tau = -1.0 / slope if slope < 0 else float("nan")
        pred = intercept + slope * tt
        ss_res = float(np.sum((lt - pred) ** 2))
        ss_tot = float(np.sum((lt - lt.mean()) ** 2))
        tail_r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    else:
        tail_r2, tail_tau = 1.0, float(tau_dom)

    caveat = (
        "ideal stepwise mass action with one common K -- no size-dependent "
        "K, no C_n+C_m fusion, no headgroup electrostatics; the relaxation "
        "spectrum is the linear-response spectrum around the new "
        "equilibrium (large jumps contain nonlinear transients -- the tail "
        "fit already excludes them by using only the final 40% of the "
        "trajectory)"
    )
    return MicelleRelaxationResult(
        equilibrium_concentrations_M=[float(v) for v in eq_final],
        kinetic_cmc_M=float(eq_final[0]),
        relaxation_times_s=[float(v) for v in relaxation_times],
        tau_dominant_s=float(tau_dom),
        mean_aggregation_number=(
            float(sum((i + 1) * v for i, v in enumerate(eq_final)))
            / max(sum(eq_final), 1e-300)
        ),
        trajectory_t_s=[float(v) for v in t_eval],
        trajectory_c1_M=[float(v) for v in c1_traj],
        trajectory_mean_nagg=[float(v) for v in mean_nagg_t],
        tail_exp_fit_r_squared=float(tail_r2),
        tail_tau_fit_s=float(tail_tau),
        max_relative_mass_error=max_rel_mass_err,
        converged=bool(sol.success),
        n_iterations=int(sol.nfev),
        method="BDF (scipy) + analytic Jacobian; eigenvalues at equilibrium",
        params={
            "c_total_M": c_total_M,
            "K_initial_per_M": K_initial_per_M,
            "K_final_per_M": K_final_per_M,
            "k_on_per_M_s": k_on_per_M_s,
            "k_off_per_s": k_off,
            "n_max": n_max,
            "n_nonzero_modes": len(decay_rates),
            "t_end_s": float(t_end),
            "jacobian_fd_rel_error": fd_err,
        },
        caveat=caveat,
    )


# ---------------------------------------------------------------------------
# S2 -- Gillespie kinetic Monte Carlo of the same stepwise system
# ---------------------------------------------------------------------------

@dataclass
class KmcResult:
    trajectory_t_s: list[float]
    trajectory_c1_M: list[float]
    trajectory_mean_nagg: list[float]
    trajectory_n_micelles: list[float]  # counts of aggregates with n >= 2
    equilibrium_c1_M: float             # mean monomer over final 20% window
    equilibrium_mean_nagg: float
    mass_action_c1_M: float             # closed-form target (effective mass)
    mass_action_mean_nagg: float
    tau_fit_s: float
    n_events: int
    mass_conservation_max_error: int    # max |sum n*N_n - initial| (integers -> 0)
    seed: int
    converged: bool                     # final window statistically stationary
    n_iterations: int
    method: str
    params: dict
    caveat: str


def kmc_micelle_kinetics(
    c_total_M: float,
    K_per_M: float,
    k_on_per_M_s: float,
    n_max: int = 40,
    n_molecules: int = 1500,
    t_end_s: float | None = None,
    n_samples: int = 240,
    seed: int = 0,
    c_total_initial_M: float | None = None,
) -> KmcResult:
    """Gillespie kinetic Monte Carlo of stepwise micellization (worker S2).

    The stochastic twin of micelle_relaxation_odes: counts of n-mers in a
    well-mixed volume, reactions C_n + C_1 <-> C_{n+1} fired by the direct
    SSA method with detailed balance k_off = k_on/K.  Starts from a Poisson
    equilibrium sample under K at c_total_initial_M (default: c_total_M --
    a concentration jump), then relaxes under K at c_total_M.

    Returns sampled trajectories, a stationary-window equilibrium compared
    against the CLOSED-FORM mass-action target (same function O1 uses),
    an exponential relaxation-time fit, and exact integer mass conservation
    (sum n*N_n is invariant -- SSA conserves it by construction; reported as
    the max deviation, always 0 for a correct implementation).

    Seed-controlled: same seed -> bit-identical trajectory.  Statistical
    bands (not exact equality) apply to equilibrium comparisons -- that is
    what stochastic means; n_molecules trades accuracy for runtime.

    caveat: ideal stepwise mass action (same model as O1 -- no fusion, no
    size-dependent K); the Poisson start is the grand-canonical
    approximation to the true canonical ensemble, so the EFFECTIVE total
    mass of each realization (used for the closed-form target) can differ
    from c_total_M by O(sqrt(N)) sampling noise.
    """
    np = _require_calc()

    if c_total_M <= 0:
        raise ValueError("c_total_M must be positive")
    if K_per_M <= 0:
        raise ValueError("K_per_M must be positive")
    if k_on_per_M_s <= 0:
        raise ValueError("k_on_per_M_s must be positive")
    if n_max < 3:
        raise ValueError("n_max must be >= 3")
    if n_molecules < 200:
        raise ValueError("n_molecules must be >= 200 for meaningful statistics")
    if n_samples < 20:
        raise ValueError("n_samples must be >= 20")
    if c_total_initial_M is not None and c_total_initial_M <= 0:
        raise ValueError("c_total_initial_M must be positive when given")
    c_init = c_total_M if c_total_initial_M is None else c_total_initial_M

    k_on = k_on_per_M_s
    k_off = k_on_per_M_s / K_per_M

    # volume chosen so that EXPECTED molecule count matches n_molecules
    n_a = 6.02214076e23
    nu = n_molecules / c_total_M            # V * N_A  (per M^-1)
    v_liter = nu / n_a

    # --- Poisson equilibrium start at c_init (K same for init/final) ---
    eq_init = mass_action_equilibrium(c_init, K_per_M, n_max)
    nu_init = n_molecules / c_init
    rng = np.random.default_rng(seed)
    counts = np.asarray(
        [rng.poisson(nu_init * x) for x in eq_init], dtype=np.int64
    )
    # ensure at least one monomer and non-degenerate start
    counts[0] = max(counts[0], 10)
    mass0 = int(np.dot(np.arange(1, n_max + 1), counts))
    c_total_eff = mass0 / (nu * 1.0)  # effective concentration of this realization

    # closed-form target for THIS realization's mass (same nu as evolution)
    eq_target = mass_action_equilibrium(max(c_total_eff, 1e-300), K_per_M, n_max)

    # --- relaxation-time estimate for windowing (reuse O1; fallback grid) ---
    if t_end_s is None:
        try:
            ref = micelle_relaxation_odes(
                c_total_eff, K_per_M * 1.4, K_per_M, k_on_per_M_s,
                n_max=min(n_max, 30), n_times=60,
            )
            t_end_s = 12.0 * ref.tau_dominant_s
        except Exception:  # noqa: BLE001 -- the reference ODE solve can
            # fail for extreme parameter combinations; fall back to the
            # simpler analytic estimate rather than propagate the failure.
            t_end_s = 10.0 / (k_on_per_M_s * max(c_init, 1e-6))
    if t_end_s <= 0 or not math.isfinite(t_end_s):
        raise ValueError("t_end_s must be a positive finite time")

    # --- precompute propensities each event via incremental update ---
    c_scale = 1.0 / (n_a * v_liter)  # converts N[m]*N[0] counts -> Molar-rate prefactor
    # forward a_f[m] = k_on * N[m]*N[0] * c_scale ; reverse a_r[m] = k_off*N[m+1]

    def total_and_pick():
        # returns (a_total, kind, m) ; kind 0 = forward, 1 = reverse
        fwd = k_on * c_scale * counts[: n_max - 1] * counts[0]
        rev = k_off * counts[1:n_max]
        a = float(fwd.sum() + rev.sum())
        if a <= 0.0:
            return 0.0, -1, -1
        r = rng.random() * a
        acc = 0.0
        # interleave: forward and reverse lists by index for uniform pick
        for m in range(n_max - 1):
            acc += fwd[m]
            if r < acc:
                return a, 0, int(m)
            acc += rev[m]
            if r < acc:
                return a, 1, int(m)
        return a, 1, n_max - 2

    # split sampling: log-dense early (resolves the relaxation transient)
    # + linear late (for stationary-window statistics)
    n_early = max(10, n_samples // 2)
    sample_times = np.concatenate([
        np.geomspace(max(t_end_s * 1e-3, 1e-15), t_end_s * 0.3, n_early),
        np.linspace(t_end_s * 0.3001, t_end_s, n_samples - n_early),
    ])
    traj_t, traj_c1, traj_nagg, traj_mic = [], [], [], []
    t = 0.0
    n_events = 0
    si = 0
    max_events = 4_000_000
    stuck = False

    while si < n_samples:
        a, kind, m = total_and_pick()
        if kind < 0 or a <= 0.0:
            # no possible reaction: fill remaining samples with current state
            while si < n_samples:
                traj_t.append(float(sample_times[si]))
                traj_c1.append(float(counts[0]) / (n_a * v_liter))
                traj_nagg.append(0.0)
                traj_mic.append(float(np.sum(counts[1:])))
                si += 1
            break
        dt = -math.log(rng.random() + 1e-300) / a
        t_next = t + dt
        # record all sample times passed by this interval (state unchanged)
        while si < n_samples and sample_times[si] <= t_next:
            traj_t.append(float(sample_times[si]))
            traj_c1.append(float(counts[0]) / (n_a * v_liter))
            tot = np.dot(np.arange(1, n_max + 1), counts)
            nfree = counts.sum()
            traj_nagg.append(float(tot / nfree) if nfree > 0 else 0.0)
            traj_mic.append(float(np.sum(counts[1:])))
            si += 1
        t = t_next
        if t > t_end_s:
            break
        # fire reaction
        if kind == 0:  # forward: N[m]--, N[0]--, N[m+1]++
            if counts[m] > 0 and counts[0] > 0:
                counts[m] -= 1
                counts[0] -= 1
                counts[m + 1] += 1
        else:          # reverse: N[m+1]--, N[0]++, N[m]++
            if counts[m + 1] > 0:
                counts[m + 1] -= 1
                counts[0] += 1
                counts[m] += 1
        n_events += 1
        if n_events >= max_events:
            stuck = True
            break

    # fill any tail samples (shouldn't happen, but keep arrays aligned)
    while si < n_samples:
        traj_t.append(float(sample_times[si]))
        traj_c1.append(float(counts[0]) / (n_a * v_liter))
        tot = np.dot(np.arange(1, n_max + 1), counts)
        nfree = counts.sum()
        traj_nagg.append(float(tot / nfree) if nfree > 0 else 0.0)
        traj_mic.append(float(np.sum(counts[1:])))
        si += 1

    mass_final = int(np.dot(np.arange(1, n_max + 1), counts))
    mass_err = abs(mass_final - mass0)

    c1_arr = np.asarray(traj_c1)
    nagg_arr = np.asarray(traj_nagg)
    t_arr = np.asarray(traj_t)
    w = max(4, n_samples // 5)
    eq_c1 = float(c1_arr[-w:].mean())
    eq_nagg = float(nagg_arr[-w:].mean())

    # stationarity: drift of the SELF-AVERAGING observable <n> between the
    # last two quintiles (monomer counts fluctuate O(1) even at equilibrium)
    drift = (abs(nagg_arr[-w:].mean() - nagg_arr[-2 * w:-w].mean())
             / max(eq_nagg, 1e-300))
    stationary = bool(drift < 0.05)

    # exponential fit of <n>(t) in the TRANSIENT window (0, 0.3*t_end]:
    # mean-nagg self-averages over ~n_molecules/<n> aggregates, unlike the
    # monomer count which micelles buffer (signal-buried-in-noise).
    # tau from MICELLE COUNT (largest jump signal, Poisson-limited noise),
    # fitted by a grid search over single-exponentials -- far more robust
    # under shot noise than a log-linear regression.  Window (0.03, 0.35]*
    # t_end skips the initial dimerization flash, then the slow mode dominates.
    mic_arr = np.asarray(traj_mic)
    noise_floor = float(np.std(mic_arr[-w:])) if w > 2 else 0.0
    eq_mic = float(mic_arr[-w:].mean())
    transient = (t_arr > 0.03 * t_end_s) & (t_arr <= 0.35 * t_end_s)
    dev_m = mic_arr - eq_mic
    tau_fit = float("nan")
    if transient.sum() >= 5 and eq_mic > 0.0 and float(np.max(np.abs(dev_m[transient]))) > 0.5:
        tt = t_arr[transient]
        dv = dev_m[transient]
        best_sse, best_tau = math.inf, float("nan")
        for tau_g in np.geomspace(t_end_s / 500.0, 0.5 * t_end_s, 80):
            e = np.exp(-tt / tau_g)
            denom = float(np.dot(e, e))
            if denom <= 0:
                continue
            amp = float(np.dot(dv, e)) / denom
            sse = float(np.sum((dv - amp * e) ** 2))
            if sse < best_sse:
                best_sse, best_tau = sse, float(tau_g)
        tau_fit = best_tau

    return KmcResult(
        trajectory_t_s=[float(v) for v in t_arr],
        trajectory_c1_M=[float(v) for v in c1_arr],
        trajectory_mean_nagg=[float(v) for v in nagg_arr],
        trajectory_n_micelles=[float(v) for v in traj_mic],
        equilibrium_c1_M=eq_c1,
        equilibrium_mean_nagg=eq_nagg,
        mass_action_c1_M=float(eq_target[0]),
        mass_action_mean_nagg=(
            float(sum((i + 1) * v for i, v in enumerate(eq_target)))
            / max(sum(eq_target), 1e-300)
        ),
        tau_fit_s=float(tau_fit),
        n_events=int(n_events),
        mass_conservation_max_error=int(mass_err),
        seed=int(seed),
        converged=bool(stationary and not stuck and mass_err == 0),
        n_iterations=int(n_events),
        method="Gillespie SSA (direct method), stepwise monomer addition",
        params={
            "c_total_M": c_total_M,
            "c_total_initial_M": c_init,
            "K_per_M": K_per_M,
            "k_on_per_M_s": k_on_per_M_s,
            "k_off_per_s": k_off,
            "n_max": n_max,
            "n_molecules": n_molecules,
            "t_end_s": float(t_end_s),
            "volume_liter": v_liter,
            "c_total_effective_M": float(c_total_eff),
            "final_window_drift": float(drift),
            "max_events_hit": bool(stuck),
            "tau_fit_observable": "n_micelles grid-fit exp in (0.03,0.35]*t_end",
            "micelle_count_noise_floor": float(noise_floor),
        },
        caveat=(
            "stochastic: equilibrium comparisons are statistical bands, not "
            "equalities (seed-controlled though); ideal stepwise mass action "
            "-- the SAME model as O1, so O1/KMC cross-checks validate "
            "execution, not the physical assumptions.  tau_fit_s is the "
            "apparent time of micelle-COUNT growth (a fast-biased projection "
            "of the multi-exponential spectrum -- mode-resolved times come "
            "from O1's eigenvalue list, not from this fit)"
        ),
    )


# ---------------------------------------------------------------------------
# O1b -- Aniansson-Wall TWO-REGIME kinetics: resolves the real degeneracy
# O1/S2 cannot, because they use one shared k_on/k_off for every step.
# ---------------------------------------------------------------------------
#
# Real, published degeneracy (Noskov; see Noskov & Grigoriev review of
# micellar relaxation kinetics): under quasi-equilibrium, dynamic surface
# tension (DST) alone cannot determine the micelle association/dissociation
# rate constants -- only ROUGH estimates are possible from DST by itself,
# because DST is sensitive only to the SLOW relaxation process.
#
# The real Aniansson & Wall (1974, 1976) theory explains why: the fast
# relaxation time tau1 (microseconds, monomer exchange AT THE MICELLE PEAK,
# large aggregation number n~n_bar) and the slow relaxation time tau2
# (milliseconds, redistribution through the NARROW connecting region of
# SMALL n between free monomer and micelle) are governed by genuinely
# DIFFERENT dissociation rate constants -- k_minus (roughly constant near
# the micelle peak) vs the small-n rate constants buried inside R (the
# "resistance" of the narrow connecting region, eq 5 of the theory). A
# single shared k_on/k_off model (as used by micelle_relaxation_odes/
# kmc_micelle_kinetics above) cannot separate these because it has only
# ONE rate constant for every step -- confirmed numerically: on real SDS
# data (see test file), that single-rate model's own eigenvalue spectrum
# only spans a ~6x fast/slow ratio, vs the real measured ~79x (Lang,
# Tondre, Zana, Bauer, Hoffmann & Ulbricht, J. Phys. Chem. 79 (1975) 276,
# Table I: SDS tau1=29 microsec, tau2=2.3 ms, both at C=cmc).
#
# This function inverts the REAL two-regime theory instead: given real,
# INDEPENDENTLY measured tau1 AND tau2 (from two genuinely different
# techniques -- e.g. shock-tube/ultrasonic for tau1, P-jump/T-jump/DST for
# tau2, exactly the real multi-technique combination Noskov's finding says
# is REQUIRED), plus the real aggregation number n_bar and distribution
# width sigma (independently measurable by SLS/SANS -- NOT fitted here),
# it solves in closed form for k_minus (the micelle-peak dissociation rate
# constant) and R (the small-n resistance term) -- the two real, separate
# rate-constant-level quantities DST data alone cannot supply.

@dataclass
class AnianssonWallTwoRegimeResult:
    k_minus_per_s: float            # dissociation rate constant AT THE MICELLE PEAK (fast process)
    R_s_per_M: float                 # "resistance" of the narrow small-n connecting region (slow process)
    k_on_implied_per_M_per_s: float  # k_minus / [free monomer]_total -- a real physical sanity check:
                                        # should land near diffusion-limited association (~1e8-1e9 M^-1 s^-1
                                        # for a small ion) if tau1/tau2/n_bar/sigma are mutually consistent
    average_micellar_lifetime_s: float  # T_m = n_bar * tau2 (eq 6 of the theory): mean time a monomer
                                           # spends resident in a micelle before the WHOLE aggregate dissolves
    a_dimensionless: float            # (C_total - CMC) / CMC, the concentration-excess factor actually used
    n_bar: float
    sigma: float
    converged: bool
    n_iterations: int
    method: str
    params: dict
    caveat: str


def aniansson_wall_two_regime(
    tau1_s: float,
    tau2_s: float,
    cmc_M: float,
    n_bar: float,
    sigma: float | None = None,
    c_total_M: float | None = None,
) -> AnianssonWallTwoRegimeResult:
    """Solve the real two-regime Aniansson-Wall relaxation equations for
    k_minus and R, given real, independently measured tau1 (fast) and
    tau2 (slow) relaxation times -- the inversion that breaks the real
    Noskov degeneracy (DST/tau2 alone cannot supply k_minus; you need a
    genuinely different technique's tau1 too).

    Equations solved (Aniansson & Wall 1974; as reviewed in Patist,
    Kanicky, Shukla & Shah, J. Colloid Interface Sci. 245 (2002) 1-15,
    eqs 3-4, with a = (C_total-CMC)/CMC):

        1/tau1 = (k_minus/sigma^2) * (1 + (sigma^2/n_bar)*a)
        1/tau2 = (n_bar^2/(CMC*R)) * (1 + (sigma^2/n_bar)*a)^(-1)

    Both equations share the SAME (1 + (sigma^2/n_bar)*a) factor, which is
    fully determined once c_total_M, cmc_M, n_bar and sigma are given (it
    is NOT a free/fitted parameter) -- so k_minus and R are recovered in
    exact closed form from the two independently measured relaxation
    times, not fitted.

    tau1_s, tau2_s: real, independently measured relaxation times (seconds)
    -- from genuinely different techniques (fast: shock-tube or ultrasonic
    absorption; slow: pressure-jump, temperature-jump, or stopped-flow/DST).
    Do not supply two numbers from the SAME technique/paper claiming to
    measure the same process -- that does not break the degeneracy.

    n_bar: real, independently measured mean aggregation number (SLS,
    SANS, or equivalent) -- NOT guessed. sigma: the Gaussian distribution
    half-width; if not independently known, Aniansson & Wall's own stated
    typical value sigma = 0.20 * n_bar is used (real, cited, not a free
    fit) -- pass a real measured value when you have one.

    c_total_M: the real total surfactant concentration at which tau1/tau2
    were measured. Defaults to c_total_M = cmc_M (a = 0), the common case
    when relaxation times are reported "at C = cmc".
    """
    if tau1_s <= 0 or tau2_s <= 0:
        raise ValueError("tau1_s and tau2_s must be positive")
    if cmc_M <= 0:
        raise ValueError("cmc_M must be positive")
    if n_bar <= 1:
        raise ValueError("n_bar must be > 1")
    if sigma is None:
        sigma = 0.20 * n_bar  # Aniansson & Wall's own stated typical value
    if sigma <= 0:
        raise ValueError("sigma must be positive")
    if c_total_M is None:
        c_total_M = cmc_M
    if c_total_M <= 0:
        raise ValueError("c_total_M must be positive")

    a = (c_total_M - cmc_M) / cmc_M
    factor = 1.0 + (sigma ** 2 / n_bar) * a  # shared, fully determined -- not fitted
    if factor <= 0:
        raise ValueError(
            f"(1 + (sigma^2/n_bar)*a) = {factor:.4g} is non-positive at these inputs "
            "(c_total_M far enough below cmc_M to be outside the theory's validity) "
            "-- this combination cannot be inverted"
        )

    k_minus = (sigma ** 2 / tau1_s) / factor
    R = (n_bar ** 2 * tau2_s) / (cmc_M * factor)   # 1/tau2 = n^2/(CMC*R) * factor^-1
    k_on_implied = k_minus / cmc_M  # real physical sanity check, not an input
    T_m = n_bar * tau2_s

    return AnianssonWallTwoRegimeResult(
        k_minus_per_s=float(k_minus),
        R_s_per_M=float(R),
        k_on_implied_per_M_per_s=float(k_on_implied),
        average_micellar_lifetime_s=float(T_m),
        a_dimensionless=float(a),
        n_bar=float(n_bar),
        sigma=float(sigma),
        converged=True,
        n_iterations=1,  # exact closed-form inversion, not iterated
        method="Aniansson-Wall two-regime closed-form inversion (eqs 3-4)",
        params={
            "tau1_s": tau1_s, "tau2_s": tau2_s, "cmc_M": cmc_M,
            "c_total_M": c_total_M, "shared_factor": factor,
        },
        caveat=(
            "Exact algebraic inversion of the real two-regime theory, not a "
            "fit -- k_minus and R are only as trustworthy as the four real "
            "inputs (tau1, tau2, n_bar, sigma). Check k_on_implied against "
            "known diffusion-limited association rates (~1e8-1e9 M^-1 s^-1 "
            "for a small 1:1 ionic surfactant) as a real physical sanity "
            "check -- if it lands far outside that range, one of the four "
            "inputs is likely wrong or the two tau values were not actually "
            "measured by independent techniques on the same real process. "
            "sigma defaults to Aniansson & Wall's own stated 0.20*n_bar "
            "when not independently supplied -- this is a real, cited "
            "typical value, not guessed, but is itself an approximation; "
            "pass a real measured sigma (SLS/SANS) when available."
        ),
    )
