"""Czakaj et al. 2021 (Molecules 26:5894, CC BY 4.0) extracted data.

Provenance supplied by Antigravity 2026-11-02: Table 4 from the full-text XML,
Figure 7 from PDF vector paths, amplitude convention confirmed as a ratio of
Fourier AMPLITUDES (not powers). No error bars are drawn on Fig. 7; each point
is a single FFT after Loess smoothing, with an instrumental noise floor of
0.1 mN/m on sigma, giving ~0.025-0.05 on the ratio.

Conditions: T = 295.15 K, RT = 2454 J/mol, pendant drop V0 = 11 uL on a 2 mm
capillary, A0 ~ 21.5 mm^2, area strain amplitude u0 ~ 0.06-0.07.
"""

from __future__ import annotations

TEMPERATURE_K = 295.15
RT_J_PER_MOL = 2454.0
AREA_STRAIN = 0.065          # u0 = dA/A0
CMC_mM = 1.0

# Table 4: limiting elasticity and characteristic diffusion frequency
TABLE_4 = {
    0.2: {"E0_mN_m": 30.1, "nu_D_Hz": 0.0032},
    0.3: {"E0_mN_m": 73.9, "nu_D_Hz": 0.0200},   # peak
    0.5: {"E0_mN_m": 45.3, "nu_D_Hz": 0.0081},
    0.8: {"E0_mN_m": 26.9, "nu_D_Hz": 0.0300},
    1.0: {"E0_mN_m": 12.6, "nu_D_Hz": 0.0600},
}

# Figure 7: ratio of second to first harmonic AMPLITUDES
FIG_7 = {
    0.2: {0.01: 0.0307, 0.10: 0.0387},
    0.5: {0.01: 0.0377, 0.10: 0.0187},
    0.8: {0.01: 0.1247, 0.10: 0.2218},
    1.0: {0.01: 0.3126, 0.10: 0.4327},
    1.5: {0.01: 0.3004, 0.10: 0.3585},
}

NOISE_FLOOR_RATIO = (0.025, 0.050)   # from 0.1 mN/m on sigma_1 = 2-4 mN/m


def lvdt_E(E0_mN_m: float, nu_D_Hz: float, nu_Hz: float):
    """Lucassen-van den Tempel complex modulus, Czakaj's Eq. 1 (frequency in Hz).

        eps = eps0 (1 + xi + i xi)/(1 + 2 xi + 2 xi^2),  xi = sqrt(nu_D/(2 nu))
    """
    import cmath
    import math

    xi = math.sqrt(nu_D_Hz / (2.0 * nu_Hz))
    return E0_mN_m * (1.0 + xi + 1j * xi) / (1.0 + 2.0 * xi + 2.0 * xi * xi)


def check_lvdt_table() -> list[dict]:
    """Consistency check: does Eq. 1 reproduce the paper's own quoted fit?"""
    out = []
    for c, row in TABLE_4.items():
        for nu in (0.01, 0.10):
            E = lvdt_E(row["E0_mN_m"], row["nu_D_Hz"], nu)
            out.append({"c_mM": c, "nu_Hz": nu, "E_real": E.real, "E_imag": E.imag,
                        "lag_deg": abs(__import__("math").degrees(__import__("math").atan2(E.imag, E.real)))})
    return out
