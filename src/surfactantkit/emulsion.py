"""Phase-3 emulsion workers: S5 (Stokes creaming + coalescence ensemble)
and O2 (breakage-coalescence population balance + LSW ripening).

All computation is pure stdlib (no scipy/numpy) -- the workers solve
piecewise-analytic ODEs and a discrete Smoluchowski system directly, per
the project's zero-dep core policy. Every result carries the standard
worker schema: converged / n_iterations / method / params / caveat.

S5 physics (standard Stokes creaming, textbook):
    v_settle = |Delta rho| * g * d^2 / (18 * eta)      [m/s]
    t_cream(d) = column_height / v_settle              [s]
    dphi_i/dt = -k_coal * phi_i  for t >= t_cream_i    (linear ODE)
    => phi_i(t) = phi_i0 * exp(-k_coal * (t - t_cream_i)) after creaming
Solved EXACTLY (no time-stepping error); the reported "stability index"
is the area under the stable-fraction curve,
    SI = integral_0^t_max (1 - F_separated(t)) dt,
the standard clarification/creaming stability index (larger = more
stable). Hindered-creaming at high phase fraction is NOT modeled
(dilute Stokes; disclosed when sum(phi) > 0.1).

O2 physics:
  * coalescence: discrete Smoluchowski aggregation
        dn_k/dt = 1/2 sum_{i+j=k} lam(i,j) n_i n_j
                  - n_k sum_j lam(k,j) n_j
    with the three classic kernel families (Smoluchowski 1917):
    constant lam0, sum lam_s*(i+j), product lam_p*i*j -- integrated by
    RK4 with a mass-conservation CONVERGENCE CHECK. Collisions whose
    product would exceed the class-K truncation boundary are handled as
    out-of-domain: reactants are consumed but the product never enters
    the tracked domain (standard PBE truncation). Total tracked volume
    therefore drifts whenever the boundary is active, and
    `mass_rel_error` is the honest adequacy diagnostic -- small only
    when K is large enough for the requested time span (otherwise the
    result reports converged=False: grow K or shorten t_max_s).
  * breakage (optional): first-order, constant rate, SYMMETRIC BINARY
    daughter rule (k -> ceil(k/2) + floor(k/2)) -- a caller-visible
    model choice, not a literature claim; disclosed in method/params.
  * Ostwald ripening (optional): classical LSW cubic law
        d32(t)^3 = d32(0)^3 + K_lsw * t   (Lifshitz-Slyozov-Wagner),
    reported both alone and superposed on the PBE curve (the
    superposition is an independent-mechanism approximation, disclosed).
"""
from __future__ import annotations

import math
from dataclasses import dataclass

G_ACCEL_M_S2 = 9.80665


@dataclass(frozen=True)
class EmulsionStabilityResult:
    times_s: tuple[float, ...]
    separated_fraction: tuple[float, ...]          # volume fraction separated at each t
    stability_index_m2: float  # integral of (1-F) dt  [s]
    t50_s: float | None        # first t with F >= 0.5 (linear interp); None if never
    d32_um: float              # Sauter mean diameter of the feed DSD
    creaming_times_s: tuple[float, ...]  # per-bin time to travel the column
    velocities_m_per_s: tuple[float, ...]
    converged: bool
    n_iterations: int
    method: str
    params: dict
    caveat: str


def _validate_common(drop_diameters_um, volume_fractions):
    if len(drop_diameters_um) != len(volume_fractions):
        raise ValueError(
            "drop_diameters_um and volume_fractions must be the same length "
            f"(got {len(drop_diameters_um)} and {len(volume_fractions)})"
        )
    if not drop_diameters_um:
        raise ValueError("drop_diameters_um must not be empty")
    if any(d <= 0 for d in drop_diameters_um):
        raise ValueError("every drop diameter must be positive (um)")
    if any(w < 0 for w in volume_fractions):
        raise ValueError("volume_fractions must be non-negative")
    if sum(volume_fractions) <= 0:
        raise ValueError("volume_fractions must contain at least one positive entry")


def simulate_emulsion_stability(
    drop_diameters_um: list[float],
    volume_fractions: list[float],
    dispersed_density_kg_m3: float,
    continuous_density_kg_m3: float,
    viscosity_Pa_s: float,
    column_height_m: float,
    coalescence_rate_per_s: float = 0.0,
    temperature_K: float = 298.15,
    n_times: int = 50,
    t_max_s: float = 86400.0,
) -> EmulsionStabilityResult:
    """S5: separated-fraction vs time for an emulsion of given DSD under
    Stokes creaming + post-creaming first-order coalescence, plus the
    standard stability index (area under the stable-fraction curve).

    drop_diameters_um: DSD bin diameters (um), volume_fractions the
    matching volume weights (need not sum to 1; normalized internally,
    the normalization disclosed in params). densities in kg/m^3,
    viscosity in Pa.s, column_height_m the creaming path length,
    coalescence_rate_per_s the FIRST-ORDER film-drainage/coalescence
    rate that applies once a bin has reached the top (0 = creaming
    only). temperature_K is recorded in params (densities/viscosity are
    given directly, not derived from T). n_times samples over
    [0, t_max_s]. Returns per-bin Stokes velocities and creaming times,
    F(t), t50, d32 of the feed, and the schema block. Guard: ValueError
    naming any bad field (lengths, non-positive d/eta/H, negative
    rates, n_times < 2).
    """
    _validate_common(drop_diameters_um, volume_fractions)
    if viscosity_Pa_s <= 0:
        raise ValueError("viscosity_Pa_s must be positive")
    if column_height_m <= 0:
        raise ValueError("column_height_m must be positive")
    if coalescence_rate_per_s < 0:
        raise ValueError("coalescence_rate_per_s must be non-negative")
    if temperature_K <= 0:
        raise ValueError("temperature_K must be positive")
    if n_times < 2:
        raise ValueError("n_times must be at least 2")
    if t_max_s <= 0:
        raise ValueError("t_max_s must be positive")

    d_um = list(drop_diameters_um)
    w = list(volume_fractions)
    w_sum = sum(w)
    w = [x / w_sum for x in w]

    delta_rho = dispersed_density_kg_m3 - continuous_density_kg_m3
    abs_dr = abs(delta_rho)
    direction = "creaming (droplets less dense than continuous; rising)"
    if delta_rho > 0:
        direction = "sedimentation (dispersed denser than continuous; settling, same Stokes magnitude)"
    elif delta_rho == 0:
        direction = "no density drive (never separates by gravity)"

    velocities = []
    creaming_times = []
    for d_um_i, weight in zip(d_um, w):
        d_m = d_um_i * 1e-6
        v = abs_dr * G_ACCEL_M_S2 * d_m * d_m / (18.0 * viscosity_Pa_s)
        velocities.append(v)
        creaming_times.append(column_height_m / v if v > 0 else math.inf)
    # volume-weighted (weight only orders reporting; F uses all bins anyway)
    _ = w

    times = [t_max_s * i / (n_times - 1) for i in range(n_times)]

    def separated_at(t: float) -> float:
        total = 0.0
        for w_i, tc in zip(w, creaming_times):
            if t < tc:
                continue
            if coalescence_rate_per_s > 0:
                total += w_i * (1.0 - math.exp(-coalescence_rate_per_s * (t - tc)))
            else:
                total += w_i
        return total

    F = [separated_at(t) for t in times]
    # trapezoid on the uniform grid for SI of (1 - F)
    dt = times[1] - times[0]
    si = 0.0
    for i in range(n_times - 1):
        si += 0.5 * ((1.0 - F[i]) + (1.0 - F[i + 1])) * dt

    t50 = None
    for i in range(1, n_times):
        if F[i] >= 0.5 > F[i - 1]:
            f0, f1 = F[i - 1], F[i]
            t50 = times[i - 1] + (0.5 - f0) / (f1 - f0) * (times[i] - times[i - 1])
            break

    # true Sauter mean (6V/A): d32 = sum(w) / sum(w/d); monodisperse -> d
    d32 = sum(w) / sum(wi / di for di, wi in zip(d_um, w))

    phi_total = sum(volume_fractions)
    caveats = [
        (
            "dilute Stokes creaming (no hindered-creaming correction); "
            "post-creaming coalescence is first-order with a single rate"
        ),
        f"gravity separation direction: {direction}",
    ]
    if phi_total > 0.1:
        caveats.append(
            f"total volume fraction {phi_total:.3g} > 0.1: hindered creaming "
            "would slow the real system relative to this dilute prediction"
        )
    if coalescence_rate_per_s == 0:
        caveats.append("coalescence_rate_per_s=0: separation = creaming only")

    return EmulsionStabilityResult(
        times_s=tuple(times),
        separated_fraction=tuple(F),
        stability_index_m2=si,
        t50_s=t50,
        d32_um=d32,
        creaming_times_s=tuple(creaming_times),
        velocities_m_per_s=tuple(velocities),
        converged=True,
        n_iterations=n_times,
        method="Stokes creaming (exact per-bin travel time) + analytic "
               "first-order coalescence ODE, trapezoid stability index",
        params={
            "drop_diameters_um": list(drop_diameters_um),
            "volume_fractions_raw": list(volume_fractions),
            "volume_fractions_normalized": w,
            "dispersed_density_kg_m3": dispersed_density_kg_m3,
            "continuous_density_kg_m3": continuous_density_kg_m3,
            "delta_rho_kg_m3": delta_rho,
            "viscosity_Pa_s": viscosity_Pa_s,
            "column_height_m": column_height_m,
            "coalescence_rate_per_s": coalescence_rate_per_s,
            "temperature_K": temperature_K,
            "t_max_s": t_max_s,
            "g_m_s2": G_ACCEL_M_S2,
        },
        caveat="; ".join(caveats),
    )


@dataclass(frozen=True)
class PopulationBalanceResult:
    times_s: tuple[float, ...]
    d32_um: tuple[float, ...]
    d32_lsw_only_um: tuple[float, ...] | None
    number_concentration: tuple[float, ...]
    mass_rel_error: float
    converged: bool
    n_iterations: int
    method: str
    params: dict
    caveat: str


def population_balance_emulsion(
    initial_number_concentrations: list[float],
    coalescence_kernel: str = "constant",
    kernel_strength: float = 1e-3,
    breakage_rate_per_s: float = 0.0,
    ripening_rate_um3_per_s: float = 0.0,
    primary_drop_diameter_um: float = 1.0,
    n_times: int = 40,
    t_max_s: float = 600.0,
    n_steps_per_interval: int = 20,
) -> PopulationBalanceResult:
    """O2: Sauter-diameter evolution d32(t) from a discrete
    breakage-coalescence population balance (bucket classes k = 1..K
    primary drops), with optional LSW Ostwald ripening.

    initial_number_concentrations: n_k(0) for classes k=1..K (any
    units; class k = cluster of k primary drops). coalescence_kernel:
    'constant' (lam0), 'sum' (lam_s*(i+j)), or 'product' (lam_p*i*j) --
    the three classic Smoluchowski kernel families with known analytic
    solution classes (validation uses the constant-kernel exact
    solution n_tot(t) = n0/(1 + lam0*n0*t/2)). kernel_strength is the
    single rate constant of the chosen family (units depend on the
    kernel: per-concentration per-second). breakage_rate_per_s > 0
    switches on first-order SYMMETRIC BINARY breakage (k ->
    ceil(k/2)+floor(k/2)) -- a disclosed model choice, not a sourced
    daughter law. ripening_rate_um3_per_s > 0 adds the classical LSW
    cubic law d32(t)^3 = d32(0)^3 + K*t (reported alone AND superposed
    on the PBE curve; the superposition is an independent-mechanism
    approximation, disclosed). primary_drop_diameter_um sets the length
    scale (d_k = primary * k^(1/3)). RK4 over n_times samples with
    n_steps_per_interval substeps. Truncation: collisions producing a
    cluster beyond class K consume their reactants but create no tracked
    product (standard out-of-domain truncation) -- tracked volume is
    conserved only while the boundary stays inactive, so
    converged = mass conservation relative error <= 1e-6 doubles as the
    "K was large enough for this time span" check. If it fails, grow
    initial_number_concentrations (larger K) or shorten t_max_s.
    d32 = sum(k*n_k) /
    sum(k^(2/3)*n_k) * primary. Guards: ValueError naming bad fields.
    """
    if not initial_number_concentrations:
        raise ValueError("initial_number_concentrations must not be empty")
    if len(initial_number_concentrations) < 2:
        raise ValueError(
            "initial_number_concentrations needs at least 2 classes "
            f"(got {len(initial_number_concentrations)})"
        )
    if any(n < 0 for n in initial_number_concentrations):
        raise ValueError("initial_number_concentrations must be non-negative")
    if sum(initial_number_concentrations) <= 0:
        raise ValueError(
            "initial_number_concentrations must contain at least one positive entry"
        )
    if coalescence_kernel not in ("constant", "sum", "product"):
        raise ValueError(
            "coalescence_kernel must be 'constant', 'sum', or 'product' "
            f"(got {coalescence_kernel!r})"
        )
    if kernel_strength < 0:
        raise ValueError("kernel_strength must be non-negative")
    if breakage_rate_per_s < 0:
        raise ValueError("breakage_rate_per_s must be non-negative")
    if ripening_rate_um3_per_s < 0:
        raise ValueError("ripening_rate_um3_per_s must be non-negative")
    if primary_drop_diameter_um <= 0:
        raise ValueError("primary_drop_diameter_um must be positive")
    if n_times < 2 or n_steps_per_interval < 1:
        raise ValueError("n_times must be >= 2 and n_steps_per_interval >= 1")
    if t_max_s <= 0:
        raise ValueError("t_max_s must be positive")

    K = len(initial_number_concentrations)
    n = [float(x) for x in initial_number_concentrations]
    mass0 = sum(k * n[k - 1] for k in range(1, K + 1))

    def lam(i: int, j: int) -> float:
        if kernel_strength == 0:
            return 0.0
        if coalescence_kernel == "constant":
            return kernel_strength
        if coalescence_kernel == "sum":
            return kernel_strength * (i + j)
        return kernel_strength * i * j

    # Collision accounting: reactant loss is the FULL double sum
    # (including pairs whose product i+j > K escapes the tracked domain);
    # gains only sum in-domain pairs -- standard PBE truncation.
    def rhs(nv: list[float]) -> list[float]:
        gain = [0.0] * (K + 1)
        loss = [0.0] * (K + 1)
        # iterate only populated classes (sparse -> fast; distribution
        # spans maybe tens of buckets even at K=256)
        idx = [k for k in range(1, K + 1) if nv[k - 1] > 0]
        for i in idx:
            li = 0.0
            for j in idx:
                li += lam(i, j) * nv[j - 1]
            loss[i] = nv[i - 1] * li
        for a_i in range(len(idx)):
            i = idx[a_i]
            for b_i in range(a_i, len(idx)):
                j = idx[b_i]
                if i + j > K:
                    continue  # out-of-domain product: loss counted above, no gain
                a = lam(i, j)
                factor = 2.0 if i != j else 1.0
                gain[i + j] += 0.5 * factor * a * nv[i - 1] * nv[j - 1]
        # breakage: first-order symmetric binary
        if breakage_rate_per_s > 0:
            for k in range(2, K + 1):
                nk = nv[k - 1]
                if nk <= 0:
                    continue
                loss[k] += breakage_rate_per_s * nk
                left = (k + 1) // 2
                right = k // 2
                for part in (left, right):
                    # daughters are strictly smaller than k <= K: always in-domain
                    gain[part] += breakage_rate_per_s * nk
        return [gain[k] - loss[k] for k in range(1, K + 1)]

    def d32_of(nv: list[float]) -> float:
        m2 = sum((k ** (2.0 / 3.0)) * nv[k - 1] for k in range(1, K + 1))
        m1 = sum(k * nv[k - 1] for k in range(1, K + 1))
        if m2 <= 0:
            return 0.0
        return (m1 / m2) * primary_drop_diameter_um

    d32_0 = d32_of(n)
    times = [t_max_s * i / (n_times - 1) for i in range(n_times)]
    d32_series = [d32_0]
    ntot_series = [sum(n)]
    h_total = times[-1] / (n_times - 1) / n_steps_per_interval if n_times > 1 else t_max_s
    for step in range(1, n_times):
        dt = times[step] - times[step - 1]
        sub = dt / n_steps_per_interval
        for _ in range(n_steps_per_interval):
            k1 = rhs(n)
            n2 = [a + sub / 2 * b for a, b in zip(n, k1)]
            k2 = rhs(n2)
            n3 = [a + sub / 2 * b for a, b in zip(n, k2)]
            k3 = rhs(n3)
            n4 = [a + sub * b for a, b in zip(n, k3)]
            k4 = rhs(n4)
            n = [max(0.0, a + sub / 6 * (b + 2 * c + 2 * d + e))
                 for a, b, c, d, e in zip(n, k1, k2, k3, k4)]
        d32_series.append(d32_of(n))
        ntot_series.append(sum(n))
    _ = h_total

    mass_end = sum(k * n[k - 1] for k in range(1, K + 1))
    mass_rel_err = abs(mass_end - mass0) / mass0 if mass0 > 0 else 0.0
    converged = mass_rel_err <= 1e-6

    lsw_series = None
    combined_series = None
    if ripening_rate_um3_per_s > 0:
        lsw_series = [
            (d32_0 ** 3 + ripening_rate_um3_per_s * t) ** (1.0 / 3.0)
            for t in times
        ]
        combined_series = [
            (d ** 3 + ripening_rate_um3_per_s * t) ** (1.0 / 3.0)
            for d, t in zip(d32_series, times)
        ]

    caveats = [
        (
            f"{K}-bucket class truncation (products beyond class K escape the "
            "tracked domain; mass_rel_error is the K-adequacy check)"
        ),
        "RK4 on the discrete Smoluchowski system (deterministic, no stochastic term)",
    ]
    if breakage_rate_per_s > 0:
        caveats.append(
            "breakage uses a symmetric-BINARY daughter rule at a constant rate "
            "(caller-chosen model, not a sourced daughter distribution)"
        )
    if ripening_rate_um3_per_s > 0:
        caveats.append(
            "ripening = classical LSW cubic law; any combined series is an "
            "independent-mechanism superposition (disclosed approximation)"
        )

    result_params = {
        "coalescence_kernel": coalescence_kernel,
        "kernel_strength": kernel_strength,
        "breakage_rate_per_s": breakage_rate_per_s,
        "ripening_rate_um3_per_s": ripening_rate_um3_per_s,
        "primary_drop_diameter_um": primary_drop_diameter_um,
        "n_classes": K,
        "t_max_s": t_max_s,
        "n_steps_per_interval": n_steps_per_interval,
        "initial_number_concentrations": list(initial_number_concentrations),
        "d32_final_um": d32_series[-1],
        "d32_combined_final_um": combined_series[-1] if combined_series else None,
    }
    return PopulationBalanceResult(
        times_s=tuple(times),
        d32_um=tuple(d32_series),
        d32_lsw_only_um=tuple(lsw_series) if lsw_series else None,
        number_concentration=tuple(ntot_series),
        mass_rel_error=mass_rel_err,
        converged=converged,
        n_iterations=n_times * n_steps_per_interval,
        method=(
            f"discrete Smoluchowski coalescence ({coalescence_kernel} kernel) "
            " + symmetric-binary breakage, RK4".replace("  ", " ")
            + (" + LSW ripening" if ripening_rate_um3_per_s > 0 else "")
        ),
        params=result_params,
        caveat="; ".join(caveats),
    )
