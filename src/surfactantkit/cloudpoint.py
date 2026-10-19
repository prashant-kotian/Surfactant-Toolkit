"""C6 -- cloud-point estimate for pure (monodisperse) alkyl ethoxylates.

Sourcing pass (done before building; every constant traceable):
  * HLD form for nonionic systems (Ghayour, Tenside Surfactants Deterg.
    60(1) 2023, Eq. 1, free access):
        HLD = b*S - k*EACN + Cc + cT*(T - 25 C)
    with b = 0.13 for NaCl (g salt / 100 mL), k in 0.15-0.17 for most
    conventional surfactants, cT ~ +0.06 for ethoxylated nonionics.
  * Cc group contribution for linear alcohol ethoxylates (Acosta 2008,
    as quoted in the same source, Eq. 2):
        Cc = 0.28 * SACN + 2.4 - EON
    (SACN = alkyl tail carbons, EON = ethylene oxide units)
  * Apparent EACN of the ethoxylate micelle core, EACN_app = 7.3
    (Zarate-Munoz, Boza Troncoso & Acosta, Langmuir 2015, 31,
    12000-12008 -- the cloud-point paper itself)
  * Cloud point = temperature at HLD = 0 for the surfactant/water (no
    oil) system:  T_cp = 25 - (b*S - k*EACN_app + Cc) / cT

Honesty [H]: the source paper reports HLD cloud points for PURE
alkyl ethoxylates typically within 5 C (and within 4 C with the HLD-NAC
follow-up); polydisperse/technical-grade ethoxylates do NOT follow the
linear CP-Cc correlation and are out of scope (disclosed, not silently
approximated); the linear CP-Cc observation was established for pure
surfactants with <= 8 EO units; k is a conventional-range constant
(0.15-0.17) -- the default 0.16 is the midpoint of that SOURCED range
and its sensitivity is echoed in params; basis is 1 wt.% surfactant.
No claim against any specific measured cloud point.
"""

from __future__ import annotations

from dataclasses import dataclass

# --- sourced constants (see module docstring for full provenance) ---
B_NACL_DEFAULT = 0.13        # Ghayour 2023 Eq.1 discussion: b = 0.13 for NaCl
K_DEFAULT = 0.16             # midpoint of sourced conventional range 0.15-0.17
CT_ETHOXYLATE_DEFAULT = 0.06  # Ghayour 2023: cT ~ +0.06 for ethoxylated nonionics
EACN_APP_ETHOXYLATE = 7.3     # Zarate-Munoz et al. Langmuir 2015
CC_A = 0.28                   # Acosta 2008 group contribution (Cc = a*SACN + 2.4 - EON)
CC_B = 2.4


@dataclass(frozen=True)
class CloudPointResult:
    """Result of `cloud_point_estimate`."""

    cloud_point_C: float
    cc: float                      # characteristic curvature of the surfactant
    hld_at_cloud_point: float       # should be ~0 (algebraic identity check)
    eacn_app: float
    converged: bool
    n_iterations: int               # closed form: 1 evaluation
    method: str
    params: dict
    caveat: str


def cloud_point_estimate(
    n_carbons: int,
    n_eo: int,
    salinity_g_per_100ml: float = 0.0,
    k: float = K_DEFAULT,
) -> CloudPointResult:
    """C6: HLD-based cloud-point estimate for a PURE (monodisperse) alkyl
    ethoxylate C_n E_m in water.

    Args:
        n_carbons: alkyl tail carbons SACN (integer >= 4).
        n_eo: ethylene oxide units EON (integer >= 1).
        salinity_g_per_100ml: NaCl in g per 100 mL (>= 0; b = 0.13 per
            the sourced NaCl value).
        k: EACN scaling constant (default 0.16 = midpoint of the sourced
            conventional 0.15-0.17 range; pass 0.15/0.17 to sensitivity-
            check).

    Returns:
        CloudPointResult (degrees C) with the standard
        converged/n_iterations/method/params/caveat schema.
        Fail-fast ValueError naming every bad field.
    """
    if isinstance(n_carbons, bool) or not isinstance(n_carbons, int):
        raise TypeError("n_carbons must be an integer")
    if n_carbons < 4:
        raise ValueError("n_carbons must be >= 4")
    if isinstance(n_eo, bool) or not isinstance(n_eo, int):
        raise TypeError("n_eo must be an integer")
    if n_eo < 1:
        raise ValueError("n_eo must be >= 1")
    if not isinstance(salinity_g_per_100ml, (int, float)) or isinstance(
        salinity_g_per_100ml, bool
    ):
        raise TypeError("salinity_g_per_100ml must be a number")
    if salinity_g_per_100ml < 0:
        raise ValueError("salinity_g_per_100ml must be >= 0")
    if not isinstance(k, (int, float)) or isinstance(k, bool):
        raise TypeError("k must be a number")
    if not 0.05 <= k <= 0.7:
        raise ValueError(
            "k must lie in the sourced conventional band 0.05-0.7 "
            "(most conventional surfactants: 0.15-0.17)"
        )

    cc = CC_A * n_carbons + CC_B - n_eo
    b = B_NACL_DEFAULT
    c_t = CT_ETHOXYLATE_DEFAULT
    # HLD = 0 at the cloud point (no oil phase; EACN = EACN_app)
    hld_at_25 = b * salinity_g_per_100ml - k * EACN_APP_ETHOXYLATE + cc
    t_cp = 25.0 - hld_at_25 / c_t
    # identity check: evaluate HLD at the returned temperature
    hld_check = hld_at_25 + c_t * (t_cp - 25.0)

    caveats = [
        (
            "pure (monodisperse) alkyl ethoxylates only: polydisperse/"
            "technical-grade ethoxylates do not follow the CP-Cc "
            "correlation and are out of scope"
        ),
        (
            "source validation: typically within 5 C of experiment (HLD) "
            "and 4 C (HLD-NAC follow-up) for pure ethoxylates; linear "
            "CP-Cc correlation established for <= 8 EO units"
        ),
        (
            f"k={k} (midpoint of sourced conventional 0.15-0.17 unless "
            "caller-overridden; band is a real uncertainty, not a fit "
            "residual)"
        ),
        (
            "1 wt.% surfactant basis; NaCl uses b=0.13; no claim against "
            "any specific measured cloud point"
        ),
    ]
    if n_eo > 8:
        caveats.append(
            "n_eo > 8: beyond the <=8-EO range where the linear CP-Cc "
            "correlation was established (extrapolation)"
        )

    params = {
        "SACN": n_carbons,
        "EON": n_eo,
        "cc": cc,
        "salinity_g_per_100ml": float(salinity_g_per_100ml),
        "b_NaCl": b,
        "k": float(k),
        "cT_per_C": c_t,
        "eacn_app": EACN_APP_ETHOXYLATE,
        "hld_at_cloud_point": hld_check,
        "sources": (
            "HLD form + constants: Ghayour, Tenside Surfactants Deterg. "
            "60(1) 2023 Eq.1 (b=0.13 NaCl, cT~+0.06, k 0.15-0.17); "
            "Cc: Acosta 2008 group contribution (0.28*SACN+2.4-EON); "
            "EACN_app=7.3: Zarate-Munoz et al., Langmuir 2015, 31, 12000"
        ),
    }
    return CloudPointResult(
        cloud_point_C=t_cp,
        cc=cc,
        hld_at_cloud_point=hld_check,
        eacn_app=EACN_APP_ETHOXYLATE,
        converged=abs(hld_check) < 1e-9,
        n_iterations=1,
        method=("HLD=0 inversion (closed form) with Acosta Cc group "
                "contribution and EACN_app=7.3"),
        params=params,
        caveat="; ".join(caveats),
    )
