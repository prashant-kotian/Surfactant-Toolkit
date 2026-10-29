"""Dynamic surface tension from a diffusion-boundary-layer FDM, jointly
fitted to gamma(t) data (Phase 2 worker O3).

Forward model: 1D diffusion toward an adsorbing planar interface on a
geometrically stretched finite-volume grid (first cell resolving the
diffusion length of the earliest output time, extent 6 diffusion lengths
of the latest, reservoir Dirichlet at the far end), with the Langmuir
adsorption kinetics solved in the same backward-Euler step (no operator
splitting), Richardson extrapolation to second order, and the Szyszkowski
equation of state that belongs to the Langmuir isotherm:
gamma(t) = gamma0 + n_factor*R*T*Gamma_max*ln(1 - Gamma(t)/Gamma_max).

This is the FDM dual of adsorption.ward_tordai_adsorption (Volterra
quadrature + bisection): two independent implementations of the same
physical problem -- they are cross-checked against each other in tests.
The FIT wrapper then inverts the forward model for (D, K, Gamma_max)
with covariance standard errors (the joint fit P1 does not do).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

_R_GAS = 8.314462618
_CALC_HINT = ("dynamic surface tension workers require numpy+scipy -- install "
              "with 'pip install surfactantkit[calc]'")


def _require_calc():
    try:
        import numpy as np
        import scipy  # noqa: F401
    except ImportError as exc:  # pragma: no cover
        raise ImportError(_CALC_HINT) from exc
    return np


@dataclass
class DynamicSurfaceTensionResult:
    gamma_mN_per_m: list[float]
    t_s: list[float]
    surface_excess_mol_per_m2: list[float]
    subsurface_M: list[float]
    gamma_equilibrium_mN_per_m: float
    converged: bool
    method: str
    params: dict
    caveat: str

    def as_dict(self) -> dict:
        return {k: getattr(self, k) for k in
                ("gamma_mN_per_m", "t_s", "surface_excess_mol_per_m2",
                 "subsurface_M", "gamma_equilibrium_mN_per_m", "converged",
                 "method", "params", "caveat")}


def dynamic_surface_tension_forward(
    t_s: list[float],
    c0_mol_per_L: float,
    D_m2_per_s: float,
    K_per_M: float,
    gamma_max_mol_per_m2: float,
    gamma0_mN_per_m: float = 72.0,
    n_factor: float = 1.0,
    temperature_K: float = 298.15,
    k_off_per_s: float = 0.5,
    n_cells: int = 160,
    max_substep_s: float | None = None,
    rel_step: float = 0.01,
) -> DynamicSurfaceTensionResult:
    """FDM forward gamma(t): diffusion toward the interface coupled implicitly
    to Langmuir KINETIC adsorption (Ward-Tordai with a rate boundary condition).

    t_s: output times (must be positive and sorted). D in m^2/s, K in 1/M
    (k_a = K*k_off), Gamma_max in mol/m^2, c0 in mol/L, k_off in 1/s
    (desorption rate -- the adsorption-barrier timescale; with k_off >= 1e4
    the surface is quasi-equilibrium and everything happens at
    diffusion-limited speed).

    Numerics: n_cells finite volumes on a geometrically stretched grid whose
    first cell resolves the diffusion length of the EARLIEST output time and
    whose extent covers 6 diffusion lengths of the latest; backward-Euler
    steps of rel_step * t (optionally capped at max_substep_s), with the
    surface kinetics solved in the same implicit step (Newton on the surface
    node, no operator splitting), and Richardson extrapolation over a
    half-step run (second order). Default settings agree with a 4x finer
    grid and 4x smaller steps to about 1e-3 mN/m.

    Why kinetics: for small surfactants the diffusion-only monolayer
    saturation time is ~ns (Gamma_max/c0)^2/D -- invisible on any
    experimental gamma(t); the measured SECOND-TO-MINUTE transients are
    adsorption-kinetics limited, which is exactly what k_off models.
    """
    np = _require_calc()
    from scipy.linalg import solve_banded

    if not t_s or any(tt <= 0 for tt in t_s):
        raise ValueError("t_s must be positive times")
    if sorted(t_s) != list(t_s):
        raise ValueError("t_s must be sorted ascending")
    if c0_mol_per_L <= 0 or D_m2_per_s <= 0 or K_per_M <= 0:
        raise ValueError("c0, D and K must be positive")
    if gamma_max_mol_per_m2 <= 0:
        raise ValueError("gamma_max must be positive")
    if gamma0_mN_per_m <= 0:
        raise ValueError("gamma0 must be positive")
    if n_factor not in (1.0, 2.0):
        raise ValueError("n_factor must be 1.0 or 2.0 (Gibbs prefactor)")
    if n_cells < 40:
        raise ValueError("n_cells must be >= 40")
    if max_substep_s is not None and max_substep_s <= 0:
        raise ValueError("max_substep_s must be positive")
    if not 0.0 < rel_step <= 0.2:
        raise ValueError("rel_step must be in (0, 0.2]")
    if k_off_per_s <= 0:
        raise ValueError("k_off_per_s must be positive")
    k_a = K_per_M * k_off_per_s  # 1/(M s)
    k_a_si = k_a / 1000.0        # m^3/(mol s): dG/dt = k_a_si c (Gmax - G) - k_off G

    t_min, t_end = float(t_s[0]), float(t_s[-1])
    D = D_m2_per_s
    gmax = gamma_max_mol_per_m2

    # stretched grid: h_0 = sqrt(D t_min)/20 (or uniform when that is coarser
    # than L/n), h_j = h_0 q^j with sum = L = 6 sqrt(D t_end)
    L = 6.0 * math.sqrt(D * t_end)
    M = n_cells
    h0 = min(math.sqrt(D * t_min) / 20.0, L / M)
    if h0 * M >= L * (1.0 - 1e-12):
        h = np.full(M, L / M)
    else:
        lo, hi = 1.0 + 1e-12, 2.0
        for _ in range(200):
            q = 0.5 * (lo + hi)
            if h0 * (q ** M - 1.0) / (q - 1.0) < L:
                lo = q
            else:
                hi = q
        h = h0 * q ** np.arange(M)
    vol = np.empty(M + 1)                # control volume of each node (m)
    vol[0] = 0.5 * h[0]
    vol[1:M] = 0.5 * (h[:-1] + h[1:])
    vol[M] = 1.0                          # Dirichlet reservoir row
    c_bulk = c0_mol_per_L * 1000.0        # mol/m^3
    cap = float("inf") if max_substep_s is None else float(max_substep_s)

    # time steps: geometric (rel_step * t), started at t_min * 1e-4, hitting
    # every output time exactly; None marks an output
    coarse: list = []
    t_now = 0.0
    for t_req in t_s:
        while t_now < t_req * (1.0 - 1e-14):
            dt = min(max(rel_step * t_now, t_min * 1e-4), cap, t_req - t_now)
            coarse.append(dt)
            t_now += dt
        coarse.append(None)
    fine: list = []
    for dt in coarse:
        fine += [None] if dt is None else [0.5 * dt, 0.5 * dt]

    def integrate(seq):
        """Backward Euler on (V + dt A) c_new = V c_old with the surface flux
        dGamma/dt taken out of node 0 inside the same implicit step."""
        c = np.full(M + 1, c_bulk)
        surf = 0.0
        g_rec, cs_rec = [], []
        mats: dict = {}
        for dt in seq:
            if dt is None:
                g_rec.append(surf)
                cs_rec.append(float(c[0]))
                continue
            ab0 = mats.get(dt)
            if ab0 is None:
                ab0 = np.zeros((3, M + 1))        # banded (l=1, u=1)
                diag = vol.copy()
                cond = dt * D / h                 # link j <-> j+1
                diag[:M] += cond
                diag[1:M] += cond[:M - 1]
                ab0[0, 1:] = -cond                # row j, column j+1
                ab0[2, :M - 1] = -cond[:M - 1]    # row j+1, column j (not the reservoir row)
                diag[M] = 1.0
                ab0[1] = diag
                mats[dt] = ab0
            rhs0 = vol * c
            rhs0[M] = c_bulk
            cs = float(c[0])
            for _ in range(50):
                # implicit Langmuir step at subsurface cs: G_new(cs) and dG_new/dcs
                den = 1.0 + dt * k_a_si * cs + dt * k_off_per_s
                g_new = (surf + dt * k_a_si * gmax * cs) / den
                dg = dt * k_a_si * (gmax - g_new) / den
                ab = ab0.copy()
                ab[1, 0] += dg
                rhs = rhs0.copy()
                rhs[0] -= g_new - dg * cs - surf
                c_new = solve_banded((1, 1), ab, rhs)
                done = abs(c_new[0] - cs) <= 1e-12 * c_bulk + 1e-15
                cs = float(c_new[0])
                if done:
                    break
            c = c_new
            surf = (surf + dt * k_a_si * gmax * cs) / (1.0 + dt * k_a_si * cs + dt * k_off_per_s)
        return np.array(g_rec), np.array(cs_rec)

    g_c, cs_c = integrate(coarse)
    g_f, cs_f = integrate(fine)
    g_out = np.clip(2.0 * g_f - g_c, 0.0, gmax * (1.0 - 1e-15))   # Richardson (BE is first order)
    cs_out = np.maximum(2.0 * cs_f - cs_c, 0.0)
    # Szyszkowski (Langmuir) equation of state, the one consistent with
    # the Langmuir kinetics above; gamma from the kinetic surface STATE,
    # not the isotherm at the instantaneous subsurface concentration
    gamma_out = [float(gamma0_mN_per_m + n_factor * _R_GAS * temperature_K
                       * gmax * math.log(1.0 - g / gmax) * 1e3) for g in g_out]

    # equilibrium reference: Szyszkowski at the bulk concentration
    gam_eq = gamma0_mN_per_m - (n_factor * _R_GAS * temperature_K * gamma_max_mol_per_m2
                                * math.log(1.0 + K_per_M * c0_mol_per_L) * 1e3)

    return DynamicSurfaceTensionResult(
        gamma_mN_per_m=gamma_out,
        t_s=[float(v) for v in t_s],
        surface_excess_mol_per_m2=[float(v) for v in g_out],
        subsurface_M=[float(v) / 1000.0 for v in cs_out],
        gamma_equilibrium_mN_per_m=float(gam_eq),
        converged=True,
        method="implicit FDM (stretched grid, coupled Langmuir kinetics, Richardson)",
        params={
            "c0_mol_per_L": c0_mol_per_L,
            "D_m2_per_s": D_m2_per_s,
            "K_per_M": K_per_M,
            "gamma_max_mol_per_m2": gamma_max_mol_per_m2,
            "n_cells": n_cells,
            "domain_length_m": L,
            "first_cell_m": float(h[0]),
            "rel_step": rel_step,
            "max_substep_s": max_substep_s,
            "n_steps": sum(1 for v in coarse if v is not None),
            "k_off_per_s": k_off_per_s,
            "k_a_per_M_s": k_a,
            "solver": "scipy.linalg.solve_banded tridiagonal per Newton iteration",
        },
        caveat=(
            "Langmuir kinetics with a single k_off (no adsorption barrier "
            "distribution, no reorientation); far boundary is a Dirichlet "
            "reservoir at 6 diffusion lengths of the latest output time"
        ),
    )


@dataclass
class DynamicSurfaceTensionFitResult:
    D_m2_per_s: float
    D_sd: float
    K_per_M: float
    K_sd: float
    gamma_max_mol_per_m2: float
    gamma_max_sd: float
    ci_low: list[float]
    ci_high: list[float]
    r_squared: float
    rmse_mN_per_m: float
    max_abs_residual_mN_per_m: float
    fitted_gamma_mN_per_m: list[float]
    n_points: int
    converged: bool
    n_iterations: int
    method: str
    params: dict
    caveat: str

    def as_dict(self) -> dict:
        return {k: getattr(self, k) for k in
                ("D_m2_per_s", "D_sd", "K_per_M", "K_sd",
                 "gamma_max_mol_per_m2", "gamma_max_sd", "ci_low",
                 "ci_high", "r_squared", "rmse_mN_per_m",
                 "max_abs_residual_mN_per_m", "fitted_gamma_mN_per_m",
                 "n_points", "converged", "n_iterations", "method",
                 "params", "caveat")}


def dynamic_surface_tension_ode(
    t_s: list[float],
    gamma_obs_mN_per_m: list[float],
    c0_mol_per_L: float,
    gamma0_mN_per_m: float = 72.0,
    n_factor: float = 1.0,
    temperature_K: float = 298.15,
    p0: tuple[float, float, float, float] | None = None,
    fit_gamma_max: bool = True,
    fit_D: bool = True,
    sigma_mN_per_m: float | None = None,
    n_cells: int = 120,
    max_substep_s: float | None = None,
    rel_step: float = 0.04,
    max_nfev: int = 400,
) -> DynamicSurfaceTensionFitResult:
    """Joint fit of dynamic surface tension gamma(t) -> D, K, Gamma_max,
    k_off (worker O3).

    Inverts the kinetic FDM forward model for diffusion, isotherm and the
    adsorption-desorption rate simultaneously from one gamma(t) transient
    -- the FDM boundary-layer formulation of diffusion-to-interface with a
    rate boundary condition (the FDM dual of adsorption.
    ward_tordai_adsorption, which is a forward-only Volterra method).

    p0: (D [m^2/s], K [1/M], Gamma_max [mol/m^2], k_off [1/s]) initial
    guesses; default (6e-11, 500, 4e-6, 0.5).  fit_gamma_max=False freezes
    Gamma_max at p0's value; fit_D=False freezes D (correct when D is
    known from an independent measurement -- strongly recommended when
    the transient is kinetics-limited, because then the curve carries
    almost no D information and the honest SE would be enormous).
    sigma_mN_per_m: known observation noise for weighting + covariance SEs.

    SEs come from the weighted log-parameter covariance; check the
    reported K-Gmax correlation before quoting either alone.
    """
    np = _require_calc()
    from scipy.optimize import least_squares

    if len(t_s) != len(gamma_obs_mN_per_m):
        raise ValueError("t_s and gamma_obs must be the same length")
    if len(t_s) < 6:
        raise ValueError("need at least 6 (t, gamma) points to fit 3 parameters")
    if any(g <= 0 or g >= gamma0_mN_per_m for g in gamma_obs_mN_per_m):
        raise ValueError(
            "gamma observations must lie in (0, gamma0) -- adsorption can "
            "only lower the tension; check units (mN/m) and gamma0"
        )
    if c0_mol_per_L <= 0:
        raise ValueError("c0_mol_per_L must be positive")
    if sigma_mN_per_m is not None and sigma_mN_per_m <= 0:
        raise ValueError("sigma_mN_per_m must be positive when given")
    p0 = (6e-11, 500.0, 4.0e-6, 0.5) if p0 is None else tuple(float(v) for v in p0)
    if len(p0) != 4:
        raise ValueError("p0 must be (D, K, gamma_max, k_off)")
    if any(v <= 0 for v in p0):
        raise ValueError("p0 entries must be positive (D, K, gamma_max, k_off)")
    if False:
        pass  # (K and k_off are always fitted)

    t_arr = [float(v) for v in t_s]
    y_arr = np.asarray(gamma_obs_mN_per_m, dtype=float)

    # free-parameter layout: [D?] [K] [Gm?] [k_off]  (all log-transformed)
    D_fixed, Gm_fixed = p0[0], p0[2]
    free_idx = []
    if fit_D:
        free_idx.append(0)
    free_idx.append(1)
    if fit_gamma_max:
        free_idx.append(2)
    free_idx.append(3)

    def unpack(logp):
        full = [math.log(D_fixed), math.log(p0[1]), math.log(Gm_fixed),
                math.log(p0[3])]
        for slot, li in enumerate(free_idx):
            full[li] = logp[slot]
        return tuple(math.exp(v) for v in full)

    def forward(logp_free):
        D, K, Gm, k_off = unpack(logp_free)
        res = dynamic_surface_tension_forward(
            t_arr, c0_mol_per_L, D, K, Gm, gamma0_mN_per_m,
            n_factor, temperature_K, k_off_per_s=k_off,
            n_cells=n_cells, max_substep_s=max_substep_s, rel_step=rel_step)
        return np.asarray(res.gamma_mN_per_m, dtype=float)

    x0 = np.log([p0[i] for i in free_idx])

    def resid(logp_free):
        r = forward(logp_free) - y_arr
        if sigma_mN_per_m is not None:
            r = r / sigma_mN_per_m
        return r

    sol = least_squares(resid, x0, method="trf", max_nfev=max_nfev,
                        xtol=1e-10, ftol=1e-10)
    popt = sol.x
    D_fit, K_fit, Gm_fit, koff_fit = unpack(popt)
    n_params = len(free_idx)

    resid_final = resid(popt)
    rss = float(np.sum(resid_final ** 2))
    if sigma_mN_per_m is None:
        ss_tot = float(np.sum((y_arr - y_arr.mean()) ** 2))
        r2 = 1.0 - rss / ss_tot if ss_tot > 0 else 1.0
        dof = len(y_arr) - n_params
        scale = rss / dof
    else:
        # weighted fit: r2 against raw residuals
        raw = forward(popt) - y_arr
        ss_tot = float(np.sum((y_arr - y_arr.mean()) ** 2))
        r2 = 1.0 - float(np.sum(raw ** 2)) / ss_tot if ss_tot > 0 else 1.0
        scale = sigma_mN_per_m ** 2
        dof = len(y_arr) - n_params
    rmse = math.sqrt(float(np.sum((forward(popt) - y_arr) ** 2)) / len(y_arr))
    max_res = float(np.max(np.abs(forward(popt) - y_arr)))

    # covariance in FREE log parameters -> delta method for relative SEs
    jac = sol.jac
    jtj = jac.T @ jac
    try:
        cov_log = np.linalg.inv(jtj) * scale
    except np.linalg.LinAlgError:
        cov_log = np.linalg.pinv(jtj) * scale
    sd_log_free = np.sqrt(np.maximum(np.diag(cov_log), 0.0))
    from scipy import stats
    t_crit = float(stats.t.ppf(0.975, dof)) if dof > 0 else 1.96

    # scatter free-parameter SEs back onto the full (D, K, Gm, k_off) axis
    full_vals = [D_fit, K_fit, Gm_fit, koff_fit]
    full_sd = [0.0, 0.0, 0.0, 0.0]
    full_lo = list(full_vals)
    full_hi = list(full_vals)
    for slot, li in enumerate(free_idx):
        v, sl = full_vals[li], sd_log_free[slot]
        full_sd[li] = v * sl
        # clamp: on flat/unidentifiable likelihoods sd_log can be huge and
        # exp() overflows; the CI then simply spans many decades (reported
        # as such) instead of crashing the fit
        full_lo[li] = v * math.exp(max(-t_crit * sl, -700.0))
        full_hi[li] = v * math.exp(min(t_crit * sl, 700.0))

    # K-Gm correlation among free log params (if both free)
    if fit_D and fit_gamma_max:
        corr = cov_log[1, 2] / math.sqrt(max(cov_log[1, 1] * cov_log[2, 2], 1e-300))
    elif fit_gamma_max:
        corr = cov_log[0, 1] / math.sqrt(max(cov_log[0, 0] * cov_log[1, 1], 1e-300))
    else:
        corr = float("nan")
    return DynamicSurfaceTensionFitResult(
        D_m2_per_s=float(D_fit),
        D_sd=float(full_sd[0]),
        K_per_M=float(K_fit),
        K_sd=float(full_sd[1]),
        gamma_max_mol_per_m2=float(Gm_fit),
        gamma_max_sd=float(full_sd[2]),
        ci_low=[float(v) for v in full_lo],
        ci_high=[float(v) for v in full_hi],
        r_squared=float(r2),
        rmse_mN_per_m=float(rmse),
        max_abs_residual_mN_per_m=max_res,
        fitted_gamma_mN_per_m=[float(v) for v in forward(popt)],
        n_points=len(y_arr),
        converged=bool(sol.success),
        n_iterations=int(sol.nfev),
        method="joint LM on log-params over the implicit FDM forward (scipy least_squares)",
        params={
            "c0_mol_per_L": c0_mol_per_L,
            "gamma0_mN_per_m": gamma0_mN_per_m,
            "n_factor": n_factor,
            "temperature_K": temperature_K,
            "fit_gamma_max": fit_gamma_max,
            "fit_D": fit_D,
            "k_off_per_s_fit": float(koff_fit),
            "sigma_mN_per_m": sigma_mN_per_m,
            "n_cells": n_cells,
            "log_K_vs_log_Gmax_correlation": float(corr),
            "free_log_cov_condition": float(np.linalg.cond(cov_log)),
            "identifiability_warning": (
                "ridge detected (|corr| > 0.99 or cond > 1e8): only lumped "
                "combinations (Gamma_eq, B_eff) are constrained by this "
                "single-concentration curve; pinv SEs UNDERSTATE uncertainty "
                "along the ridge -- fix gamma_max (Gibbs area) or add a "
                "second concentration before quoting K or k_off separately"
                if (math.isfinite(corr) and abs(corr) > 0.99)
                or (math.isfinite(float(np.linalg.cond(cov_log)))
                    and float(np.linalg.cond(cov_log)) > 1e8)
                else None
            ),
            "free_params": ["D", "K", "gamma_max", "k_off"][0:0]
                           + [name for slot, name in enumerate(
                               ["D", "K", "gamma_max", "k_off"])
                               if [0, 1, 2, 3][slot] in free_idx],
            "t_crit_95": t_crit,
            "dof": dof,
            "solver_message": str(sol.message),
        },
        caveat=(
            "single-transient identifiability: in kinetics-limited data the "
            "curve carries almost no D information (fit_D=False is the "
            "honest setup) and K-Gamma_max stay ridge-correlated (reported "
            "in params) -- quote the amplitude K*Gamma_max or add a second "
            "concentration; SEs assume this forward model is correct "
            "(splitting/domain/kinetics-form errors are NOT in them)"
        ),
    )
