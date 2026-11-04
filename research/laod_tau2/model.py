"""Reaction-diffusion forward model for LAOD gamma(t) with a micellar reservoir.

Why this and not the two-compartment model in the white paper
-------------------------------------------------------------
The white paper's Code Block 1 uses a FIXED boundary-layer thickness
delta = sqrt(pi*D/omega) and a subsurface ODE
    dc_s/dt = (j_diff - j_ads)/delta,   j_diff = (2D/delta)(c0 - c_s)
That is a one-way filter, not a boundary layer: the bulk reservoir never
loses what the surface gains, so integral surfactant mass is not conserved.
For a freshly created interface, or at the strain amplitudes (eps up to 0.5)
this work is about, that error is first-order, not cosmetic.

This module solves the actual 1D problem with the micelles carried as an
explicit, immobile species, so total surfactant (bulk monomer + micellar +
adsorbed) is conserved *by construction*, not by assertion:

    dc/dt  = D*d2c/dx2 - dM/dt
    dM/dt  = -(M - M_eq(c)) / tau2,      M_eq(c) = C_m*c
    dGamma/dt = D*dc/dx|_0 - Gamma*dlnA/dt
    Gamma  = Gamma_inf*K*c_s/(1 + K*c_s)          (Langmuir, local equilibrium)
    gamma  = gamma0 + n*R*T*Gamma_inf*ln(1 - Gamma/Gamma_inf)

Two reservoir parameters, and both are needed:
    capacity  C_m = d(micellar)/d(monomer)   -- how MUCH surfactant the
                micelles can supply within a cycle
    tau2      -- how FAST they supply it (Aniansson-Wall slow relaxation)
The white-paper prototype had neither; the probe variable "Cm" alone could
not separate "strong reservoir" from "fast reservoir".

Boundary conditions: Dirichlet at the far end (semi-infinite reservoir), and
the surface node carries the Langmuir isotherm plus the dilatational term.

Disclosed limitations -- each is either anchored or bracketed by a test:
  * micellar diffusion neglected. The old justification -- "micelles are much
    slower than monomers, so it matters only at high frequency" -- is WRONG as
    stated (FINDINGS 30.14, audit item 2.6). A smaller D does not imply a
    negligible contribution when the micellar inventory and susceptibility are
    large. For two diffusing pools in local equilibrium the admittance scales as
    sqrt((D + C_m D_m/n)(1 + C_m/n)), so the immobile approximation is low by

        sqrt(1 + (C_m/n)(D_m/D)),        D_m/D ~ 0.25  ->  Fainerman 2008

    which is sqrt(1+(m-1)*0.25): 1.06x at 1.5x CMC, 1.22x at 4x CMC, 1.58x at
    10x CMC. So the neglect is defensible below ~2x CMC and NOT defensible at
    the concentrations where the external data live. Stated as a bounded error
    rather than dismissed.
  * reservoir linearised about equilibrium -- valid for the moderate
    amplitudes a LAOD experiment uses, NOT for a fresh-interface step.
  * no adsorption barrier (interface at local equilibrium). The library's
    dynamic_tension module covers the kinetically limited case; the slow
    process of interest here is the micellar one.
  * planar interface; pendant-drop curvature/depletion is a separate
    correction that lowers E at low omega (documented in the literature).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

R_GAS = 8.314462618


@dataclass
class ModelParams:
    """SI throughout: c mol/m^3, K m^3/mol, Gamma_inf mol/m^2, D m^2/s."""

    c0: float
    K: float
    gamma_inf: float
    D: float
    gamma0: float = 72.0
    n_factor: float = 1.0
    temperature_K: float = 298.15
    frumkin_a: float = 0.0       # lateral interaction; 0 = Langmuir
    capacity: float = 0.0        # C_m = dM_eq/dc, a SUSCEPTIBILITY (dimensionless)
    tau2: float = 1.0            # exchange coefficient, see relaxation_time
    n_agg: float = 1.0           # aggregation number, for the inventory closure
    n_cells: int = 96
    stretch: float = 4.0
    extent_factor: float = 8.0

    # ---- reservoir bookkeeping (FINDINGS 30.1, 30.2) --------------------
    #
    # C_m is dM_eq/dc: a DERIVATIVE, not the ratio M_eq/c. For a monodisperse
    # mass-action closure M_eq = n K_n c^n, so
    #
    #     C_m = n * M_eq / c      =>      M_eq(c) = C_m * c / n_agg.
    #
    # Two distinct errors followed from conflating these:
    #
    #   (a) the inventory was initialised as M0 = C_m*c0, which is n_agg times
    #       too large (5.06 vs 0.0355 mol/m^3 at 1.5x CMC, n=95). Constant
    #       offsets cancel in the linearised dynamics, so phases were unaffected,
    #       but the mass-error denominator was inflated by the same factor and
    #       every positivity diagnostic was meaningless.
    #
    #   (b) the exchange law relaxed M towards C_m*c instead of C_m*c/n_agg, so
    #       the relaxing quantity was the susceptibility times c, not an
    #       inventory. The closed two-pool eigenmode is then
    #
    #           tau_closed = tau2 / (1 + C_m/n_agg)
    #
    #       not tau2. At the old default n_agg = 1 and C_m = 47.5, tau_closed =
    #       tau2/48.5 -- the coefficient we called "tau2 = 4 s" was a 0.08 s
    #       relaxation. Verified to 7 decimals against ODE integration
    #       (FINDINGS 29.1).
    #
    # n_agg defaults to 1.0, which reproduces every pre-2026-11-05 result
    # bit-for-bit. Real use should pass the aggregation number.
    @property
    def micelle_inventory(self) -> float:
        """Equilibrium micellar inventory M_eq(c0) = C_m*c0/n_agg."""
        return self.capacity * self.c0 / self.n_agg

    @property
    def relaxation_time(self) -> float:
        """Closed two-pool relaxation time, tau2/(1 + C_m/n_agg).

        This is the quantity that should be compared with a measured micellar
        relaxation time. `tau2` is the coefficient in the exchange law and is
        NOT itself a physical time unless C_m/n_agg << 1.
        """
        return self.tau2 / (1.0 + self.capacity / self.n_agg)

    def micelle_eq(self, c) -> "object":
        """Equilibrium micellar inventory as a function of local monomer c."""
        return self.capacity * c / self.n_agg

    # ---- equilibrium ---------------------------------------------------
    def _theta_eq(self, c_s: float) -> float:
        """Coverage in equilibrium with subsurface concentration c_s.

        Isotherm (Frumkin, with the sign convention used throughout this
        module: Pi = -nRT Gamma_inf [ln(1-th) + a th^2]):

            K c_s = th/(1-th) * exp(-2 a th)

        so the residual is

            f(th) = ln(th/(1-th)) - 2 a th - ln(K c_s),      df/dth = 1/(th(1-th)) - 2a.

        BUG FIX 2026-11-03 (two separate defects, both previously masked):
          1. theta0 used the plain Langmuir formula K c/(1+K c) even when
             frumkin_a != 0. For a = 4, K c = 0.0066 that returns 0.0065 when
             the correct coverage is 0.892. gamma_s_eq is built on theta0, so
             every Frumkin run started from a near-bare interface, relaxed
             towards equilibrium instead of oscillating about it, and produced
             a flat gamma(t). That was misread (FINDINGS 13.9) as numerical
             stiffness -- it was a sign error in the residual.
          2. The old Newton used df = 1/(th(1-th)) + 2a, i.e. the wrong sign on
             the interaction term. df/dth genuinely goes NEGATIVE past the
             spinodal th(1-th) = 1/(2a), so a plain Newton can walk away from
             the root. We bracket on a grid and bisect instead, which is slower
             but cannot diverge, and we select the HIGH-coverage branch because
             the adsorbed monolayer is the high-coverage solution.

        Above the spinodal (a > ~2) the isotherm is multi-valued and the branch
        choice is a physical assumption, not arithmetic. That is the 2D
        condensation boundary referred to as a bounded model limit.
        """
        if self.frumkin_a == 0.0:
            return self.K * c_s / (1.0 + self.K * c_s)

        lkc = math.log(self.K * c_s)
        a = self.frumkin_a

        def f(th: float) -> float:
            return math.log(th / (1.0 - th)) - 2.0 * a * th - lkc

        # Dense grid: the function is smooth and the roots are well separated.
        n = 4096
        grid = [(i + 0.5) / n for i in range(n)]
        brackets = []
        prev_th, prev_f = grid[0], f(grid[0])
        for th in grid[1:]:
            cur_f = f(th)
            if prev_f == 0.0:
                brackets.append((prev_th, prev_th))
            elif prev_f * cur_f < 0.0:
                brackets.append((prev_th, th))
            prev_th, prev_f = th, cur_f
        if not brackets:
            # No sign change: saturate to whichever end is closest to zero.
            return 1e-6 if abs(f(grid[0])) < abs(f(grid[-1])) else 1.0 - 1e-6

        lo, hi = brackets[-1]                      # high-coverage branch
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if f(lo) * f(mid) <= 0.0:
                hi = mid
            else:
                lo = mid
            if hi - lo < 1e-14:
                break
        return 0.5 * (lo + hi)

    @property
    def theta0(self) -> float:
        return self._theta_eq(self.c0)

    @property
    def gamma_s_eq(self) -> float:
        return self.gamma_inf * self.theta0

    def isotherm(self, c_s: float) -> float:
        """Langmuir if frumkin_a == 0, else the Frumkin isotherm (implicit)."""
        return self.gamma_inf * self._theta_eq(c_s)

    def _cs_of_gamma(self, G: float) -> float:
        """Inverse isotherm: subsurface concentration from surface excess."""
        th = min(max(G / self.gamma_inf, 1e-12), 1.0 - 1e-12)
        return th / (1.0 - th) * math.exp(-2.0 * self.frumkin_a * th) / self.K

    def dgamma_dcs(self, c_s: float) -> float:
        """d gamma / d c_s at the interface, mol/m^2 per (mol/m^3), i.e. metres.

        BUG FIX 2026-11-05 (FINDINGS 30.6). The Langmuir branch returned
        dGamma/dc (mol/m^2 / (mol/m^3) = m) while the Frumkin branch returned
        dPi/dc (mN/m / (mol/m^3)). Those differ by nRT*Gamma_inf, so the two
        branches were returning quantities in different units under one name.
        The audited Langmuir runs were unaffected; any Frumkin extension using
        this would have been silently wrong.

        This function now returns dGamma/dc for BOTH branches. The surface
        tension derivative is exposed separately as dPi_dcs.
        """
        if self.frumkin_a == 0.0:
            return self.gamma_inf * self.K / (1.0 + self.K * c_s) ** 2
        th = self.isotherm(c_s) / self.gamma_inf
        th = min(max(th, 1e-12), 1.0 - 1e-12)
        dth_dcs = th * (1.0 - th) / (
            c_s * (1.0 - 2.0 * self.frumkin_a * th * (1.0 - th))
        )
        return self.gamma_inf * dth_dcs

    def dPi_dcs(self, c_s: float) -> float:
        """dPi/dc_s = -dgamma/dc_s, mN/m per (mol/m^3). Gibbs, either branch."""
        return self.n_factor * R_GAS * self.temperature_K * 1e3 * self.dgamma_dcs(c_s)

    def surface_tension(self, gamma_s):
        import numpy as np

        g = np.minimum(np.asarray(gamma_s, dtype=float), self.gamma_inf * (1.0 - 1e-12))
        th = np.maximum(g / self.gamma_inf, 1e-12)
        rt = self.n_factor * R_GAS * self.temperature_K * self.gamma_inf * 1e3
        return self.gamma0 + rt * (np.log1p(-th) + self.frumkin_a * th * th)


@dataclass
class OscillationResult:
    t: list
    gamma_mN_m: list
    area_rel: list
    gamma_surface_excess: list
    mass_error: float
    n_cycles: int
    mass_residual_abs: float = float("nan")
    mass_endpoint_drift: float = float("nan")
    params: dict = field(default_factory=dict)


def _grid(p: ModelParams, omega_min: float):
    import numpy as np

    L = p.extent_factor * math.sqrt(p.D / omega_min)
    s = (np.exp(p.stretch * (np.arange(p.n_cells) + 1) / p.n_cells) - 1.0) / (math.exp(p.stretch) - 1.0)
    edges = np.concatenate([[0.0], L * s])
    return 0.5 * (edges[1:] + edges[:-1]), np.diff(edges)


def run_oscillation(
    p: ModelParams,
    omega: float,
    eps: float,
    n_cycles: int = 24,
    n_per_cycle: int = 320,
    rtol: float = 1e-7,
    atol: float = 1e-11,
    n_return_cycles: int = 4,
) -> OscillationResult:
    """Integrate the oscillating-area response; return the last `n_return_cycles`.

    n_return_cycles defaults to 4, which reproduces every number published before
    2026-11-03 bit-for-bit. Raising it is free precision: phase noise on the Fourier
    fit falls as 1/sqrt(N), measured 0.046 deg at 4 cycles and 0.026 deg at 12
    (FINDINGS 19.3). Prior to this parameter the function silently returned 4 cycles
    whatever the caller asked for, so the protocol's "retain >= 12 cycles" was never
    honoured by the analysis.
    """
    import numpy as np
    from scipy.integrate import solve_ivp

    centres, widths = _grid(p, omega)
    n_int = p.n_cells - 1                  # integrated bulk cells (last is Dirichlet)
    dV = widths[:n_int]
    h_surf = centres[0]
    h_face = np.diff(centres)              # n_cells-1 faces, h_face[i] between centre i,i+1

    Ginf, K, D = p.gamma_inf, p.K, p.D
    Cm, tau2 = p.capacity, p.tau2
    c10, M0 = p.c0, p.micelle_inventory

    def dlnA_dt(t):
        A = 1.0 + eps * math.sin(omega * t)
        return eps * omega * math.cos(omega * t) / A

    def area(t):
        import numpy as np

        return 1.0 + eps * np.sin(omega * np.asarray(t, dtype=float))

    # Layout: [Gamma, c_0..c_{n-2}, M_0..M_{n-2}].
    #
    # Gamma is the DIFFERENTIAL surface variable and the subsurface
    # concentration cs is ALGEBRAIC through the inverse Langmuir isotherm:
    #     Gamma = Gamma_inf*K*cs/(1+K*cs)   =>   cs = Gamma/(K*(Gamma_inf-Gamma))
    # CORRECTION 2026-11-05 (FINDINGS 30.7): the comment that used to sit here
    # claimed that carrying cs as the differential state makes the interface
    # "infinitely compressible" and collapses E to ~1/E0. That is FALSE. For an
    # invertible isotherm, f'(cs) dcs/dt = j - f(cs) dlnA/dt is the same
    # equation as integrating Gamma, under a change of variable; the two agree
    # to solver tolerance (tests/test_model_equivalence.py). What fails in a
    # naive cs formulation is omitting the f'(cs) storage factor, not the
    # choice of variable. We keep Gamma as the state because it is well
    # conditioned, not because cs is wrong.
    iG = 0
    ib = 1
    iM = 1 + n_int

    G_eq = p.gamma_s_eq
    y0 = np.concatenate([[G_eq], np.full(n_int, c10), np.full(n_int, M0)])

    n_agg = p.n_agg

    def subsurface(G):
        """Inverse isotherm: cs(Gamma), guarded away from full coverage."""
        return p._cs_of_gamma(G)

    def rhs(t, y):
        G = y[iG]
        c = y[ib : ib + n_int]
        M = y[iM : iM + n_int]
        dy = np.zeros_like(y)
        dlnA = dlnA_dt(t)

        cs = subsurface(G)
        j_in = D * (c[0] - cs) / h_surf            # mol/(m^2 s), + = adsorption

        # ---- surface -----------------------------------------------------
        dGdt = j_in - G * dlnA

        # ---- bulk: finite-volume diffusion ------------------------------
        fluxes = np.empty(n_int)
        fluxes[:-1] = D * (c[1:] - c[:-1]) / h_face[: n_int - 1]
        fluxes[-1] = D * (p.c0 - c[-1]) / h_face[n_int - 1]
        dc = np.empty(n_int)
        dc[0] = (fluxes[0] - j_in) / dV[0]
        dc[1:] = (fluxes[1:] - fluxes[:-1]) / dV[1:]

        # ---- micellar exchange (conservative) --------------------------
        if Cm > 0.0:
            # M_eq(c) = Cm*c/n_agg, NOT Cm*c. Relaxing towards Cm*c made the
            # relaxing quantity a susceptibility rather than an inventory and
            # put a spurious n_agg (=95) into the closed eigenmode.
            dM = -((M - Cm * c / n_agg) - (M0 - Cm * c10 / n_agg)) / tau2
        else:
            dM = np.zeros(n_int)
        dc -= dM

        dy[iG] = dGdt
        dy[ib : ib + n_int] = dc
        dy[iM : iM + n_int] = dM
        return dy

    if n_cycles < 5:
        raise ValueError("n_cycles must be >= 5: the last 4 cycles are returned as the steady state")
    T_per = 2.0 * math.pi / omega
    t_end = n_cycles * T_per
    t_eval = np.linspace(t_end - n_return_cycles * T_per, t_end,
                         n_return_cycles * n_per_cycle)

    def total_surfactant(y):
        return float(np.sum(dV * y[1 : 1 + n_int]) + np.sum(dV * y[iM : iM + n_int]) + y[iG])

    sol = solve_ivp(
        rhs, (0.0, t_end), y0, t_eval=t_eval, method="BDF",
        rtol=rtol, atol=atol, max_step=T_per / 160.0,
    )
    if not sol.success:
        raise RuntimeError(f"integration failed: {sol.message}")

    Gs = sol.y[iG]

    # Mass-balance diagnostic. This was previously hard-coded to nan (an
    # unimplemented stub), which was then misread as a numerical failure.
    # Total surfactant is computed the same way the RHS conserves it:
    # integral of bulk monomer + integral of micellar + adsorbed Gamma.
    # A semi-infinite Dirichlet boundary means the domain can exchange with
    # the far field, so this monitors DRIFT of the interior sum, which the
    # conservative exchange term must not change.
    # FINDINGS 30.3: the earlier diagnostic reported the endpoint change of
    # (integral c + integral M + Gamma) and called it a mass error. That sum is
    # NOT conserved by the equations: differentiating it gives
    #
    #     d/dt [ int(c+M) + Gamma ] = J_far(0) - J_far(t) - Gamma*dlnA/dt
    #
    # because (i) the Dirichlet far boundary exchanges surfactant, and (ii) the
    # surface is stretching, so Gamma*dlnA/dt is a genuine sink. A small number
    # from that expression was luck, not conservation. We now report the
    # residual of the actual identity:
    #
    #     residual(t) = [int(c+M) + Gamma](t) - [int(c+M) + Gamma](0)
    #                   + int_0^t Gamma*dlnA/dt dt - int_0^t [J_in - J_far] dt
    #
    # which is zero up to solver tolerance if the discretisation is consistent.
    def _inventory(yv):
        return float(np.sum(yv[1 : 1 + n_int] * dV)
                     + np.sum(yv[iM : iM + n_int] * dV) + yv[iG])

    def _far_flux(yv):
        return float(D * (p.c0 - yv[1 + n_int - 1]) / h_face[n_int - 1])

    def _in_flux(yv):
        return float(D * (yv[1] - p._cs_of_gamma(yv[0])) / h_surf)

    m0 = _inventory(y0)
    m_end = _inventory(sol.y[:, -1])
    # trapezoidal accumulation of the two source terms over the returned window
    ts = sol.t
    src_area = 0.0
    src_flux = 0.0
    for k in range(1, len(ts)):
        dt = ts[k] - ts[k - 1]
        src_area += 0.5 * dt * (sol.y[iG, k] * dlnA_dt(ts[k])
                                + sol.y[iG, k - 1] * dlnA_dt(ts[k - 1]))
        src_flux += 0.5 * dt * ((_in_flux(sol.y[:, k]) - _far_flux(sol.y[:, k]))
                                + (_in_flux(sol.y[:, k - 1]) - _far_flux(sol.y[:, k - 1])))
    delta = m_end - m0
    residual = delta + src_area - src_flux
    # Scale by the surfactant actually EXCHANGED over the window, not by the
    # net integral (which nearly cancels and made the ratio meaningless) and not
    # by the total inventory (which was inflated by n_agg under the old closure).
    exchanged = 0.0
    for k in range(1, len(ts)):
        dt = ts[k] - ts[k - 1]
        j = 0.5 * (abs(_in_flux(sol.y[:, k]) - _far_flux(sol.y[:, k]))
                   + abs(_in_flux(sol.y[:, k - 1]) - _far_flux(sol.y[:, k - 1])))
        g = 0.5 * (abs(sol.y[iG, k] * dlnA_dt(ts[k]))
                   + abs(sol.y[iG, k - 1] * dlnA_dt(ts[k - 1])))
        exchanged += dt * (j + g)
    mass_error = abs(residual) / max(exchanged, 1e-30)
    mass_residual_abs = abs(residual)
    mass_endpoint_drift = abs(delta) / abs(m0) if m0 else float("nan")

    return OscillationResult(
        t=sol.t.tolist(),
        gamma_mN_m=p.surface_tension(Gs).tolist(),
        area_rel=area(sol.t).tolist(),
        gamma_surface_excess=Gs.tolist(),
        mass_error=mass_error,
        mass_residual_abs=mass_residual_abs,
        mass_endpoint_drift=mass_endpoint_drift,
        n_cycles=n_cycles,
        params={"omega": omega, "eps": eps, "capacity": p.capacity, "tau2": p.tau2,
                "n_agg": p.n_agg, "tau_closed": p.relaxation_time},
    )


def total_surfactant_at_equilibrium(p: ModelParams, omega: float) -> float:
    """Reference total (bulk monomer + micellar + adsorbed) at eps = 0."""
    _, widths = _grid(p, omega)
    vol = float(sum(widths[:-1]))
    return p.c0 * vol + p.micelle_inventory * vol + p.gamma_s_eq
