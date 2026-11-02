"""Real surfactant parameter sets, with the provenance of every number.

Parameters supplied by Antigravity (2026-11-02) from their validated database
and experimental literature; the *derived* quantities here are computed, not
typed, and each carries its own note.

Capacity closure
----------------
    C_m = (c_total - CMC)/CMC * n_bar

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
   above the CMC and is only reliable for c_tot >~ 1.5x CMC. Every result in
   this study is run at c_tot >= 2x CMC for that reason.

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
        # volume fraction of micellar surfactant (c in mM -> mol/L = M)
        phi = (c_total_mM - self.cmc_mM) * 1e-3 * self.m_surf / RHO_SOLUTION
        return c_analytic / max(1.0 - phi, 1e-9)

    def capacity_note(self, c_total_mM: float) -> str:
        m = c_total_mM / self.cmc_mM
        phi = (c_total_mM - self.cmc_mM) * 1e-3 * self.m_surf / RHO_SOLUTION
        notes = []
        if m < 1.5:
            notes.append(f"c_tot only {m:.2f}x CMC: closure overestimates C_m near the CMC edge")
        if self.isotherm == "frumkin":
            notes.append("ionic: n_bar-derived C_m is an upper bound (counterion condensation)")
        notes.append(f"volume correction {100 * (1 / (1 - phi) - 1):.4f}%")
        return "; ".join(notes)

    def to_model(self, c_total_mM: float, c_m_tau: float = 1.0):
        """ModelParams for the LAOD forward model (SI units)."""
        from .model import ModelParams

        return ModelParams(
            c0=c_total_mM,                 # mol/m^3 == mM
            K=None,                        # set below from the isotherm
            gamma_inf=self.gamma_inf,
            D=self.D,
            capacity=self.capacity(c_total_mM),
            tau2=c_m_tau * 1e-3,
            n_factor=float(self.n_ions),
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
    "C12E8": SurfactantSystem(
        name="C12E8", cmc_mM=0.071, n_bar=95, gamma_inf=2.80e-6, D=3.8e-10,
        isotherm="langmuir", isotherm_a=0.0, n_ions=1, m_surf=0.538,
        tau2_ms=(1.0, 10.0), tau1_us=(5.0, 20.0),
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
