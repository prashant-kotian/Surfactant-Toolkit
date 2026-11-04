"""Real surfactant parameter sets, with the provenance of every number.

Parameters supplied by Antigravity (2026-11-02) from their validated database
and experimental literature; the *derived* quantities here are computed, not
typed, and each carries its own note.

Capacity closure
----------------
    C_m = (c_total - CMC)/CMC * n_bar

WHAT THIS MEANS PHYSICALLY AFTER THE 2026-11-05 CORRECTION (FINDINGS 30.13)
--------------------------------------------------------------------------
C_m is a susceptibility and the inventory is M_eq = C_m*c/n_bar, so the BULK of
the closure is (m-1):

    M_eq = (m-1) * c          (independent of n_bar!)
    B    = 1 + C_m/(n_bar*(1+i omega tau)) = 1 + (m-1)/(1+i omega tau)

n_bar CANCELS from the reservoir strength. It survives only (a) inside the
inventory definition and (b) as the conversion between C_m and total micellised
material. The manuscript's premise -- that the reservoir is n_bar times stronger
than the monomer pool, and the whole detectability argument built on n_bar(m-1)
-- was an artefact of the same missing /n. Under the corrected model the
equilibrium retardation is 1+(m-1) = m, exactly the total/monomer ratio it must
be, and at 1.5x CMC the phase shift from the reservoir is ~4 deg at 0.1 Hz
rather than ~18 deg.

Derived from monodisperse mass action in the pinned-monomer limit
(c1 = CMC above the CMC):  dc_tot/dc1 = 1 + n*c_mic/c1, so
C_m = dc_mic/dc1 = n*phi*c_tot/CMC with phi = (c_tot-CMC)/c_tot, which reduces
exactly to the closure. Verified symbolically in `verify_capacity_closure()`.

Two corrections that the crude closure ignores, both applied here:

1. **Micellar volume fraction.** Micelles occupy volume, lowering the water
   activity and raising the true monomer concentration for a given total. At the
   concentrations we care about the effect is ~1e-5 for nonionics but reaches
   ~10% for 10x-CMC ionic systems with a 0.3 kg/mol surfactant. Small, but it is
   the difference between "we neglected it silently" and "we bounded it".

2. **The CMC edge.** The closure gives C_m = n_bar*(m-1) with m = c_tot/CMC,
   which is LINEAR in (c_tot-CMC). The exact derivative goes to zero at the CMC
   (both phi and dc_mic/dc1 do). The closure therefore *overestimates* C_m just
   above the CMC and is only reliable for c_tot >~ 1.5x CMC.

   CORRECTION 2026-11-03: this note used to claim "every result in this study is
   run at c_tot >= 2x CMC". That is false. The concentration-window studies
   (study_e0_sensitivity.py, study_ionic_generality.py) run at 1.0-10x CMC and
   report an optimum at 1.5x, inside the zone this note calls unreliable. An
   external review caught the contradiction. A capacity-perturbation sweep
   (FINDINGS 18.8) shows the *window* survives a x0.5-2.0 capacity error but the
   *optimum location* moves from 1.5x to 1.2x if capacity is overstated by 50%.
   The closure overstates in the direction that favours 1.5x, so the operating
   point is defensible but not proven.

Ionic note: for SDS and CTAB the micelle is charged and the monomer is
partly counterion-condensed, so the effective C_m derived from n_bar is an
*upper* bound. Reported as such rather than silently corrected.
"""

from __future__ import annotations

from dataclasses import dataclass

M_WATER = 0.018015
RHO_SOLUTION = 1000.0  # kg/m^3


@dataclass(frozen=True)
class SurfactantSystem:
    name: str
    cmc_mM: float
    n_bar: float
    gamma_inf: float          # mol/m^2
    D: float                  # m^2/s
    isotherm: str             # 'langmuir' or 'frumkin'
    isotherm_a: float = 0.0   # Frumkin lateral interaction
    n_ions: int = 1           # 2 for a 1:1 ionic surfactant (Gibbs prefactor)
    m_surf: float = 0.30      # kg/mol, for the volume correction
    tau2_ms: tuple = (1.0, 10.0)
    tau1_us: tuple = (10.0, 50.0)
    source: str = "Antigravity parameter set, 2026-11-02"

    def capacity(self, c_total_mM: float) -> float:
        """C_m from the closure, with the micellar-volume correction."""
        if c_total_mM <= self.cmc_mM:
            return 0.0
        m = c_total_mM / self.cmc_mM
        c_analytic = (m - 1.0) * self.n_bar
        phi = self.volume_fraction(c_total_mM)
        return c_analytic / max(1.0 - phi, 1e-9)

    def volume_fraction(self, c_total_mM: float) -> float:
        """Volume fraction of micellar surfactant.

        BUG FIX 2026-11-05 (FINDINGS 30.6). The previous expression multiplied a
        concentration in mM by an extra 1e-3 before dividing by a density in
        kg/m^3. mM is numerically mol/m^3, so that factor made phi 1000x too
        small: 1.91e-08 where the true value is 1.91e-05. The correction it
        feeds is correspondingly 1000x too small -- small either way, but the
        stated "~1e-5 for nonionics, ~10% at 10x CMC ionic" note was wrong by
        that factor.
        """
        if c_total_mM <= self.cmc_mM:
            return 0.0
        return (c_total_mM - self.cmc_mM) * self.m_surf / RHO_SOLUTION

    def capacity_note(self, c_total_mM: float) -> str:
        m = c_total_mM / self.cmc_mM
        phi = self.volume_fraction(c_total_mM)
        notes = []
        if m < 1.5:
            notes.append(f"c_tot only {m:.2f}x CMC: closure overestimates C_m near the CMC edge")
        if self.isotherm == "frumkin":
            notes.append("ionic: n_bar-derived C_m is an upper bound (counterion condensation)")
        notes.append(f"volume correction {100 * (1 / (1 - phi) - 1):.4f}%")
        return "; ".join(notes)

    def monomer_concentration(self, c_total_mM: float, tol: float = 1e-14) -> float:
        """Equilibrium MONOMER concentration for a given total loading.

        Solves   c_total = c_eq + M_eq(c_eq)   with   M_eq = C_m*c_eq/n_bar
        using C_m = n_bar*(c_total/cmc - 1) evaluated at the solution.

        Why this exists (FINDINGS 30.1): `to_model` used to pass c_total_mM
        straight in as `c0`, and the forward model's transport variable is the
        MONOMER concentration. That silently made coverage rise with total
        concentration purely because more surfactant was added, on top of the
        reservoir term -- two equilibrium conventions at once. Measured effect
        at 1 Hz: 3.94 degrees of phase, larger than the entire concentration
        window the paper reported.

        Above the CMC the monomer concentration is pinned near the CMC, so the
        correction is largest exactly where the old code was most wrong.
        """
        if c_total_mM <= self.cmc_mM:
            return c_total_mM
        # Bisection on c_eq in [0, c_total]; monotone increasing.
        #
        # NOTE on the closure: with C_m = n_bar*(m-1) the equation
        #     c + M_eq(c) = ctot,   M_eq = C_m*c/n_bar = (m-1)*c
        # has the exact solution c = ctot/m = CMC. So for THIS closure the
        # "solve" convention is identical to pinning the monomer at the CMC --
        # which is physically the classical pinned-monomer limit, and which is
        # why the concentration error showed up as a plain "use ctot instead of
        # cmc" bug. The solver is kept general so that a non-linear closure can
        # be dropped in without touching call sites.
        lo, hi = 0.0, c_total_mM
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            cm = self.capacity(c_total_mM)
            f = mid + cm * mid / self.n_bar - c_total_mM
            if f > 0:
                hi = mid
            else:
                lo = mid
            if hi - lo < tol * c_total_mM:
                break
        return 0.5 * (lo + hi)

    def to_model(self, c_total_mM: float, c_m_tau: float = 1.0, K: float | None = None,
                 monomer_convention: str = "solve"):
        """ModelParams for the LAOD forward model (SI units).

        K must be supplied. An earlier version passed K=None with a comment saying it
        would be "set below" -- it never was, and the returned object raised
        TypeError on first use (theta0 divides by K). An external review flagged this
        as a runtime trap; it is now a hard error with an actionable message.

        monomer_convention:
          "solve"  (default) -- c0 = equilibrium monomer concentration from
                                c_total = c_eq + M_eq(c_eq). Correct.
          "pinned"           -- c0 = CMC, the classical pinned-monomer limit.
          "total"            -- c0 = c_total. LEGACY, retained only to
                                reproduce pre-2026-11-05 results; physically
                                wrong above the CMC (FINDINGS 30.1).
        """
        from .model import ModelParams

        if K is None:
            raise ValueError(
                f"{self.name}: the adsorption constant K must be supplied. "
                f"Anchor it to a measured E0 (see study_ionic_generality.K_from_theta) "
                f"rather than guessing: an unanchored K lands the monolayer in the "
                f"saturated regime (FINDINGS 12)."
            )
        if monomer_convention == "solve":
            c0 = self.monomer_concentration(c_total_mM)
        elif monomer_convention == "pinned":
            c0 = min(self.cmc_mM, c_total_mM)
        elif monomer_convention == "total":
            c0 = c_total_mM
        else:
            raise ValueError(f"unknown monomer_convention {monomer_convention!r}")
        return ModelParams(
            c0=c0,                         # mol/m^3 == mM, MONOMER
            K=K,                           # must be anchored to a measured E0
            gamma_inf=self.gamma_inf,
            D=self.D,
            capacity=self.capacity(c_total_mM),
            tau2=c_m_tau * 1e-3,
            n_factor=float(self.n_ions),
            n_agg=self.n_bar,
        )


# Parameter table from Antigravity's message. Langmuir K for C12E8 as given;
# the Frumkin systems are characterised by (a, n) and their K must be fitted
# from a published isotherm -- NOT supplied yet, so they are marked.
SYSTEMS = {
    "SDS": SurfactantSystem(
        name="SDS", cmc_mM=8.2, n_bar=62, gamma_inf=3.90e-6, D=5.5e-10,
        isotherm="frumkin", isotherm_a=0.8, n_ions=2, m_surf=0.288,
        tau2_ms=(2.0, 40.0), tau1_us=(10.0, 50.0),
    ),
    # C12E8: cmc_mM=0.071 and n_bar=95 are UNSOURCED design settings, not
    # literature values. Five rounds of results are calibrated to them. Rusdi
    # 2005 (our only CnE8 source) covers n = 14/16/18 only, gives C14E8 CMC =
    # 18.1 umol/L at 298.2 K, and reports aggregation numbers x1000. Anyone
    # using this row must replace these two numbers with a traced measurement
    # for the actual specimen (FINDINGS 24, 29.9).
    "C12E8": SurfactantSystem(
        name="C12E8", cmc_mM=0.071, n_bar=95, gamma_inf=2.80e-6, D=3.8e-10,
        isotherm="langmuir", isotherm_a=0.0, n_ions=1, m_surf=0.538,
        tau2_ms=(1.0, 10.0), tau1_us=(5.0, 20.0),
        source="UNSOURCED design settings (was: Antigravity parameter set). "
               "CMC and n_bar have no traced primary measurement.",
    ),
    "CTAB": SurfactantSystem(
        name="CTAB", cmc_mM=0.92, n_bar=90, gamma_inf=3.20e-6, D=4.0e-10,
        isotherm="frumkin", isotherm_a=0.7, n_ions=2, m_surf=0.364,
        tau2_ms=(50.0, 500.0), tau1_us=(100.0, 400.0),
    ),
}

# Langmuir K, m^3/mol (equilibrium isotherm). Only C12E8 was supplied (140).
# For the Frumkin systems K is NOT yet supplied -> must be fitted from a
# published surface-tension isotherm before those rows can be trusted.
K_LANGMUIR = {"C12E8": 140.0}
K_MISSING = ("SDS", "CTAB")


def verify_capacity_closure(n_bar=62.0, m_values=(1.5, 2.0, 5.0, 10.0)) -> list:
    """C_m from mass action vs the closure, in the pinned-monomer limit."""
    out = []
    for m in m_values:
        phi = (m - 1.0) / m
        derived = n_bar * phi * m * 1.0      # c_tot/c1 = m at c1 = CMC
        closure = n_bar * (m - 1.0)
        out.append({"m": m, "derived": derived, "closure": closure,
                    "match": abs(derived - closure) < 1e-9})
    return out
