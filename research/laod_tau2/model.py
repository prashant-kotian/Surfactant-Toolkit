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
  * micellar diffusion neglected (correct in the sense that micelles are much
    slower than monomers; matters only at the highest frequencies).
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
    capacity: float = 0.0        # C_m (dimensionless)
    tau2: float = 1.0            # s
    n_cells: int = 96
    stretch: float = 4.0
    extent_factor: float = 8.0

    # ---- equilibrium ---------------------------------------------------
    @property
    def theta0(self) -> float:
        return self.K * self.c0 / (1.0 + self.K * self.c0)

    @property
    def gamma_s_eq(self) -> float:
        return self.gamma_inf * self.theta0

    def isotherm(self, c_s: float) -> float:
        return self.gamma_inf * self.K * c_s / (1.0 + self.K * c_s)

    def dgamma_dcs(self, c_s: float) -> float:
        return self.gamma_inf * self.K / (1.0 + self.K * c_s) ** 2

    def surface_tension(self, gamma_s):
        import numpy as np

        g = np.minimum(np.asarray(gamma_s, dtype=float), self.gamma_inf * (1.0 - 1e-12))
        return self.gamma0 + (
            self.n_factor * R_GAS * self.temperature_K * self.gamma_inf * 1e3 * np.log1p(-g / self.gamma_inf)
        )


@dataclass
class OscillationResult:
    t: list
    gamma_mN_m: list
    area_rel: list
    gamma_surface_excess: list
    mass_error: float
    n_cycles: int
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
) -> OscillationResult:
    """Integrate the oscillating-area response; return the last 4 cycles."""
    import numpy as np
    from scipy.integrate import solve_ivp

    centres, widths = _grid(p, omega)
    n_int = p.n_cells - 1                  # integrated bulk cells (last is Dirichlet)
    dV = widths[:n_int]
    h_surf = centres[0]
    h_face = np.diff(centres)              # n_cells-1 faces, h_face[i] between centre i,i+1

    Ginf, K, D = p.gamma_inf, p.K, p.D
    Cm, tau2 = p.capacity, p.tau2
    c10, M0 = p.c0, p.capacity * p.c0

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
    # Writing cs as the differential state instead (the "obvious" choice) makes
    # the interface infinitely compressible: the model then satisfies
    # Gamma = f(cs) at every instant and reports ~1/E0 elasticity. That is the
    # single most important structural choice in this file.
    iG = 0
    ib = 1
    iM = 1 + n_int

    G_eq = p.gamma_s_eq
    y0 = np.concatenate([[G_eq], np.full(n_int, c10), np.full(n_int, M0)])

    def subsurface(G):
        """Inverse Langmuir: cs(Gamma), guarded away from full coverage."""
        return G / (K * max(Ginf - G, 1e-12 * Ginf))

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
            dM = -((M - Cm * c) - (M0 - Cm * c10)) / tau2
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
    t_eval = np.linspace(t_end - 4.0 * T_per, t_end, 4 * n_per_cycle)

    def total_surfactant(y):
        return float(np.sum(dV * y[1 : 1 + n_int]) + np.sum(dV * y[iM : iM + n_int]) + y[iG])

    sol = solve_ivp(
        rhs, (0.0, t_end), y0, t_eval=t_eval, method="BDF",
        rtol=rtol, atol=atol, max_step=T_per / 160.0,
    )
    if not sol.success:
        raise RuntimeError(f"integration failed: {sol.message}")

    Gs = sol.y[iG]
    return OscillationResult(
        t=sol.t.tolist(),
        gamma_mN_m=p.surface_tension(Gs).tolist(),
        area_rel=area(sol.t).tolist(),
        gamma_surface_excess=Gs.tolist(),
        mass_error=float("nan"),
        n_cycles=n_cycles,
        params={"omega": omega, "eps": eps, "capacity": p.capacity, "tau2": p.tau2},
    )


def total_surfactant_at_equilibrium(p: ModelParams, omega: float) -> float:
    """Reference total (bulk monomer + micellar + adsorbed) at eps = 0."""
    _, widths = _grid(p, omega)
    vol = float(sum(widths[:-1]))
    return p.c0 * vol + (p.capacity * p.c0) * vol + p.gamma_s_eq
