"""Dynamic surface tension from a diffusion-boundary-layer FDM, jointly
fitted to gamma(t) data (Phase 2 worker O3).

Forward model: 1D diffusion toward an adsorbing planar interface solved
by Crank-Nicolson on a uniform grid sized from the diffusion length
(L = 4*sqrt(D*t_max), reservoir Dirichlet at the far end), with the
Langmuir isotherm applied algebraically at the surface cell after each
diffusion step (operator splitting, disclosed), and
gamma(t) = gamma0 - n_factor*R*T*Gamma(c_s(t)).

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
    max_substep_s: float = 0.1,
) -> DynamicSurfaceTensionResult:
    """FDM forward gamma(t): Crank-Nicolson diffusion + Langmuir KINETIC
    surface (Ward-Tordai with a rate boundary condition).

    t_s: output times (must be positive and sorted); the solver substeps
    internally with dt <= max_substep_s (Crank-Nicolson has no diffusion
    stability limit; this bounds splitting/curvature error instead).
    D in m^2/s, K in 1/M (k_a = K*k_off), Gamma_max in mol/m^2, c0 in
    mol/L, k_off in 1/s (desorption rate -- the adsorption-barrier
    timescale; with k_off >= 1e4 the surface is quasi-equilibrium and
    everything happens at diffusion-limited speed).

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
    if max_substep_s <= 0:
        raise ValueError("max_substep_s must be positive")
    if k_off_per_s <= 0:
        raise ValueError("k_off_per_s must be positive")
    k_a = K_per_M * k_off_per_s  # 1/(M s)

    t_end = float(t_s[-1])
    # domain sized from the FASTEST diffusion expected in the fit window
    L = 4.0 * math.sqrt(D_m2_per_s * t_end) + 1e-12
    M = n_cells
    dx = L / M
    # concentrations internal in mol/m^3, D in m^2/s, Gamma mol/m^2
    c_bulk = c0_mol_per_L * 1000.0
    c = np.full(M + 1, c_bulk, dtype=float)
    c[-1] = c_bulk

    # CN matrices: (I - dt/2 * A) c_new = (I + dt/2 * A) c_old
    # interior + node0 (half-cell: dc0/dt = 2D(c1-c0)/dx^2)
    def build(dt: float):
        ab = np.zeros((3, M + 1))  # banded (l=1,u=1)
        # row 0: half-cell; r = D*dt/(2*dx^2) (D included -- the discrete
        # diffusion number controlling both CN stability-free coefficients)
        r = D_m2_per_s * dt / (2.0 * dx * dx)
        # (I - dt/2 A): diag0 = 1 + 2*r ; super0 = -2*r   [A: dc0 = 2D(c1-c0)/dx2]
        # (I + dt/2 A): rhs diag0' = 1 - 2*r ; super0' = +2*r
        # rows 1..M-1: diag = 1 + 2*r ; sub/super = -r
        # row M: Dirichlet -> identity
        diag = np.ones(M + 1)
        sub = np.zeros(M)
        sup = np.zeros(M)
        diag[0] = 1.0 + 2.0 * r
        sup[0] = -2.0 * r
        for j in range(1, M):
            diag[j] = 1.0 + 2.0 * r
            sub[j - 1] = -r
            sup[j] = -r
        diag[M] = 1.0
        sub[M - 1] = 0.0
        # banded storage rows: ab[0,1:] = super, ab[1,:] = diag, ab[2,:-1] = sub
        ab[0, 1:] = sup
        ab[1, :] = diag
        ab[2, :-1] = sub
        # rhs matrix application done explicitly (symmetric A)
        return ab, r

    def apply_rhs(c_old, r):
        out = np.empty_like(c_old)
        # row0: (1 - 2r) c0 + 2r c1
        out[0] = (1.0 - 2.0 * r) * c_old[0] + 2.0 * r * c_old[1]
        out[1:M] = ((1.0 - 2.0 * r) * c_old[1:M]
                    + r * c_old[0:M - 1] + r * c_old[2:M + 1])
        out[M] = c_old[M]
        return out

    # gamma from the actual SURFACE STATE (kinetic Gamma), NOT the
    # isotherm at instantaneous cs -- during transients they differ by the
    # very kinetics this model exists to describe
    gamma_out: list[float] = []
    g_out: list[float] = []
    cs_out: list[float] = []
    surf_state = 0.0  # Gamma as STATE: starts clean
    t_prev = 0.0
    for t_req in t_s:
        while t_prev < t_req - 1e-15:
            dt = min(max_substep_s, t_req - t_prev)
            ab, r = build(dt)
            rhs = apply_rhs(c, r)
            c_new = solve_banded((1, 1), ab, rhs)
            # Langmuir KINETIC surface step (exact linear ODE over dt with
            # cs held), mass-conserving: the half-cell pays for uptake; the
            # positive cap (dG <= half-cell inventory) keeps the surface from
            # taking what diffusion has not delivered (rarely binding now --
            # kinetics rate-limit uptake, the cap is a guard).
            half_cell = 0.5 * dx
            cs_M = max(c_new[0], 0.0) / 1000.0
            B = k_a * cs_M + k_off_per_s
            g_eq = gamma_max_mol_per_m2 * (k_a * cs_M / B)
            g_new = g_eq + (surf_state - g_eq) * math.exp(-B * dt)
            dG = g_new - surf_state
            if dG > 0.0:
                available = max(c_new[0], 0.0) * half_cell
                dG = min(dG, available)
            else:
                dG = max(dG, -surf_state)
            c_new[0] = max(c_new[0] - dG / half_cell, 0.0)
            surf_state += dG
            c_new[M] = c_bulk
            c = c_new
            t_prev += dt
        gamma_out.append(float(gamma0_mN_per_m - n_factor * _R_GAS
                                * temperature_K * surf_state * 1e3))
        g_out.append(float(surf_state))
        cs_out.append(float(c[0] / 1000.0))

    # equilibrium reference: gamma at the bulk isotherm value
    G_eq = gamma_max_mol_per_m2 * (K_per_M * c0_mol_per_L
                                   / (1.0 + K_per_M * c0_mol_per_L))
    gam_eq = gamma0_mN_per_m - n_factor * _R_GAS * temperature_K * G_eq * 1e3

    return DynamicSurfaceTensionResult(
        gamma_mN_per_m=gamma_out,
        t_s=[float(v) for v in t_s],
        surface_excess_mol_per_m2=g_out,
        subsurface_M=cs_out,
        gamma_equilibrium_mN_per_m=float(gam_eq),
        converged=True,
        method="crank-nicolson + algebraic langmuir split (FDM)",
        params={
            "c0_mol_per_L": c0_mol_per_L,
            "D_m2_per_s": D_m2_per_s,
            "K_per_M": K_per_M,
            "gamma_max_mol_per_m2": gamma_max_mol_per_m2,
            "n_cells": n_cells,
            "domain_length_m": L,
            "max_substep_s": max_substep_s,
            "k_off_per_s": k_off_per_s,
            "k_a_per_M_s": k_a,
            "solver": "scipy.linalg.solve_banded tridiagonal per step",
        },
        caveat=(
            "Langmuir kinetics with a single k_off (no adsorption barrier "
            "distribution, no reorientation); operator splitting is first "
            "order in the surface update (halving max_substep_s changes "
            "early-time gamma by O(dt)); far boundary is a Dirichlet "
            "reservoir at 4 diffusion lengths"
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
    max_substep_s: float = 0.15,
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
            n_cells=n_cells, max_substep_s=max_substep_s)
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
        method="joint LM on log-params over CN-FDM forward (scipy least_squares)",
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
