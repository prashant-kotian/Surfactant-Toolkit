"""P4 -- 1D convection-dispersion transport with Langmuir adsorption in a
porous column (Crank-Nicolson). Produces breakthrough curves (EOR /
low-salinity flooding relevance; the `hld` module anchors the formulation
context, not this PDE).

Physics modeled:
  * advection-dispersion: dc/dt = -v dc/dx + D d2c/dx2
  * equilibrium Langmuir adsorption on the solid phase, capacity and
    affinity expressed PER PORE VOLUME (q(c) = qmax*b*c/(1+b*c)), so the
    coupled storage term is d(c + q)/dt
  * Dirichlet inlet (c(0,t) = c_in, Danckwerts corrections not applied --
    disclosed), zero-gradient outlet, zero initial concentration.

Numerics:
  * Crank-Nicolson (theta = 1/2) on a uniform grid, second-order central
    space, hand-rolled Thomas tridiagonal solve (house zero-dependency
    style -- no scipy).
  * Nonlinearity: Picard iteration with a CONSERVATIVE secant storage
    linearization (delta q measured against delta c per node), so the
    aqueous + adsorbed balance is exact up to the spatial quadrature.
  * Output samples at n_times moments (n_steps_per_interval internal CN
    steps each).

Honesty [H]: phenomenological column model. No fit to, or claim against,
any specific core-flood or published breakthrough set. The Langmuir form
is an equilibrium isotherm assumption (same validity limits as the
adsorption module: ideal dilute solution, fixed temperature, no surface
heterogeneity); Crank-Nicolson is A-stable but dispersive node-Peclet
oscillations (v*dx/D >~ 2) are disclosed in params/caveat, not silently
clamped; concentrations are NOT clamped non-negative (undershoot is
reported via min_concentration).
"""

from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class Transport1DResult:
    """Result of `surfactant_transport_1d`."""

    times_s: tuple[float, ...]                  # output sample times
    x_m: tuple[float, ...]                      # grid nodes (0..L)
    concentration_profiles: tuple[tuple[float, ...], ...]  # c(x) per sample
    breakthrough_fraction: tuple[float, ...]    # c_out(t)/c_in at samples
    breakthrough_t50_s: float | None            # first time c_out >= 0.5 c_in
    mass_balance_residual: float                # relative |injected - stored - out|
    min_concentration: float                    # undershoot monitor (should be ~0)
    converged: bool                             # Picard tolerance met at every step
    n_iterations: int                           # total Picard/Thomas solves
    method: str
    params: dict
    caveat: str


def _thomas(lower: list[float], diag: list[float], upper: list[float],
            rhs: list[float]) -> list[float] | None:
    """Solve a tridiagonal system. Returns None if a pivot is ~0."""
    n = len(diag)
    c = [0.0] * n
    d = [0.0] * n
    if abs(diag[0]) < 1e-300:
        return None
    c[0] = upper[0] / diag[0]
    d[0] = rhs[0] / diag[0]
    for i in range(1, n):
        piv = diag[i] - lower[i] * c[i - 1]
        if abs(piv) < 1e-300:
            return None
        if i < n - 1:
            c[i] = upper[i] / piv
        d[i] = (rhs[i] - lower[i] * d[i - 1]) / piv
    out = [0.0] * n
    out[n - 1] = d[n - 1]
    for i in range(n - 2, -1, -1):
        out[i] = d[i] - c[i] * out[i + 1]
    return out


def surfactant_transport_1d(
    velocity_m_s: float,
    dispersion_m2_s: float,
    column_length_m: float,
    inlet_conc_mol_m3: float,
    langmuir_capacity_mol_m3: float = 0.0,
    langmuir_affinity_m3_per_mol: float = 0.0,
    n_cells: int = 200,
    t_max_s: float = 600.0,
    n_times: int = 400,
    n_steps_per_interval: int = 5,
) -> Transport1DResult:
    """P4: Crank-Nicolson breakthrough curves for 1D convection-dispersion
    with equilibrium Langmuir adsorption.

    Args:
        velocity_m_s: Darcy/pore velocity v (m/s), >= 0.
        dispersion_m2_s: hydrodynamic dispersion D (m^2/s), > 0.
        column_length_m: column length L (m), > 0.
        inlet_conc_mol_m3: feed concentration c_in (mol/m^3), > 0.
        langmuir_capacity_mol_m3: qmax per pore volume (mol/m^3), >= 0
            (0 = non-adsorbing tracer).
        langmuir_affinity_m3_per_mol: Langmuir b (m^3/mol), >= 0.
        n_cells: grid intervals (>= 10).
        t_max_s: simulated duration (s), > 0.
        n_times: output samples (>= 2).
        n_steps_per_interval: CN steps per output sample (>= 1).

    Returns:
        Transport1DResult with c(x,t) samples, outlet breakthrough,
        t50, discrete mass-balance residual, Picard diagnostics and the
        standard converged/n_iterations/method/params/caveat schema.
        Fail-fast: ValueError naming every bad field.
    """
    if velocity_m_s < 0:
        raise ValueError("velocity_m_s must be >= 0")
    if dispersion_m2_s <= 0:
        raise ValueError("dispersion_m2_s must be > 0")
    if column_length_m <= 0:
        raise ValueError("column_length_m must be > 0")
    if inlet_conc_mol_m3 <= 0:
        raise ValueError("inlet_conc_mol_m3 must be > 0")
    if langmuir_capacity_mol_m3 < 0:
        raise ValueError("langmuir_capacity_mol_m3 must be >= 0")
    if langmuir_affinity_m3_per_mol < 0:
        raise ValueError("langmuir_affinity_m3_per_mol must be >= 0")
    if n_cells < 10:
        raise ValueError("n_cells must be >= 10")
    if t_max_s <= 0:
        raise ValueError("t_max_s must be > 0")
    if n_times < 2:
        raise ValueError("n_times must be >= 2")
    if n_steps_per_interval < 1:
        raise ValueError("n_steps_per_interval must be >= 1")

    v, D, L, c0 = velocity_m_s, dispersion_m2_s, column_length_m, inlet_conc_mol_m3
    qmax, b = langmuir_capacity_mol_m3, langmuir_affinity_m3_per_mol
    N = n_cells
    dx = L / N
    n_steps = n_times * n_steps_per_interval
    dt = t_max_s / n_steps
    h = 0.5                           # CN theta weight (storage sp carries 1/dt)
    adsorbing = qmax > 0 and b > 0

    def q_of(c: float) -> float:
        return qmax * b * c / (1.0 + b * c)

    c = [0.0] * (N + 1)
    x = [i * dx for i in range(N + 1)]

    out_times: list[float] = []
    out_profiles: list[tuple[float, ...]] = []
    breakthrough: list[float] = []
    out_acc = 0.0       # exact left-Riemann sum of outlet flux v*c_out*dt
    injected = 0.0      # exact sum of inlet flux v*c_in*dt
    n_iter_total = 0
    converged = True
    min_c = 0.0

    for step in range(n_steps):
        c_new = list(c)
        picard_ok = False
        for _k in range(50):
            mfac = [1.0] * (N + 1)
            if adsorbing:
                for i in range(N + 1):
                    dc = c_new[i] - c[i]
                    if abs(dc) > 1e-14:
                        mfac[i] = 1.0 + (q_of(c_new[i]) - q_of(c[i])) / dc
                    else:
                        mfac[i] = 1.0 + qmax * b / (1.0 + b * c[i]) ** 2
            lower = [0.0] * N
            diag = [0.0] * N
            upper = [0.0] * N
            rhs = [0.0] * N
            for i in range(1, N + 1):
                j = i - 1
                if i < N:
                    cm = v / (2.0 * dx) + D / dx ** 2
                    cii = -2.0 * D / dx ** 2
                    cp = -v / (2.0 * dx) + D / dx ** 2
                else:
                    cm = v / (2.0 * dx) + D / dx ** 2
                    cii = -v / (2.0 * dx) - D / dx ** 2
                    cp = 0.0
                sp = mfac[i] / dt
                diag[j] = sp - h * cii
                lower[j] = -h * cm if i > 1 else 0.0
                upper[j] = -h * cp if i < N else 0.0
                c_next = c[i + 1] if i < N else 0.0
                r = sp * c[i] + h * (cm * c[i - 1] + cii * c[i] + cp * c_next)
                if i == 1:
                    r -= (-h * cm) * c0
                rhs[j] = r
            solved = _thomas(lower, diag, upper, rhs)
            if solved is None:
                converged = False
                break
            new_vec = [c0] + solved
            delta = max(abs(new_vec[i] - c_new[i]) for i in range(N + 1))
            c_new = new_vec
            n_iter_total += 1
            if delta <= 1e-10 * max(c0, 1.0):
                picard_ok = True
                break
        if not picard_ok:
            converged = False
        c = c_new

        # exact left-Riemann flux sums each step; undershoot monitor
        injected += v * c0 * dt
        out_acc += v * c[N] * dt
        min_c = min(min_c, min(c))

        if (step + 1) % n_steps_per_interval == 0:
            out_times.append((step + 1) * dt)
            out_profiles.append(tuple(c))
            breakthrough.append(c[N] / c0)
            min_c = min(min_c, min(c))

    # balance closure: injected = outflux integral + inventory AT T
    # (trapezoid in x over the final profile)
    stored = (math.fsum(c) * dx - 0.5 * (c[0] + c[N]) * dx) + (
        math.fsum(q_of(ci) for ci in c) * dx
        - 0.5 * (q_of(c[0]) + q_of(c[N])) * dx
    )
    denom = injected if injected > 0 else 1.0
    residual = abs(injected - stored - out_acc) / denom

    # t50 on the breakthrough samples (linear interpolation)
    t50 = None
    for i in range(1, len(breakthrough)):
        f0, f1 = breakthrough[i - 1], breakthrough[i]
        if f0 < 0.5 <= f1:
            t50 = (out_times[i - 1]
                   + (0.5 - f0) / (f1 - f0) * (out_times[i] - out_times[i - 1]))
            break

    pe_cell = v * dx / D
    caveats = [
        (
            "phenomenological column model: no fit to or claim against any "
            "specific core-flood or published breakthrough set"
        ),
        (
            "equilibrium Langmuir adsorption (dilute, ideal, fixed T); capacity "
            "expressed per pore volume; Dirichlet inlet without Danckwerts "
            "correction, zero-gradient outlet"
        ),
        (
            f"node Peclet v*dx/D = {pe_cell:.3g}"
            + (" (dispersive CN oscillation risk above ~2 -- undershoot reported "
               "in min_concentration, never silently clamped)"
               if pe_cell > 2 else "")
        ),
    ]
    if not adsorbing:
        caveats.append("langmuir_capacity=0 or affinity=0: non-adsorbing tracer")
    if not converged:
        caveats.append("Picard did not reach 1e-10 relative tolerance: grow "
                       "n_steps_per_interval")
    caveats.append("dynamic surface tension / interfacial rheology NOT solved here "
                   "(see ward_tordai_adsorption; the display incident, not transport)")

    params = {
        "v_m_per_s": v, "D_m2_per_s": D, "length_m": L, "c_in_mol_m3": c0,
        "qmax_mol_m3": qmax, "b_m3_per_mol": b,
        "n_cells": N, "dx_m": dx, "t_max_s": t_max_s, "dt_s": dt,
        "n_steps": n_steps, "pe_cell": pe_cell,
        "reduced_velocity_s_inv": v / L,
        "breakthrough_pore_volumes_at_half": (v * t50 / L) if t50 else None,
    }
    return Transport1DResult(
        times_s=tuple(out_times),
        x_m=tuple(x),
        concentration_profiles=tuple(out_profiles),
        breakthrough_fraction=tuple(breakthrough),
        breakthrough_t50_s=t50,
        mass_balance_residual=residual,
        min_concentration=min_c,
        converged=converged,
        n_iterations=n_iter_total,
        method="Crank-Nicolson (theta=1/2) + central space, Thomas tridiagonal, "
               "Picard with secant storage linearization for equilibrium Langmuir",
        params=params,
        caveat="; ".join(caveats),
    )
