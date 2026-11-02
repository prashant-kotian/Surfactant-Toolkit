"""Physics sanity bounds for supplied parameter sets.

CORRECTION recorded here: an earlier version of this file claimed
"E0 cannot exceed n*R*T*Gamma_inf". That is WRONG. For the Langmuir/Frumkin
equation of state

    E0 = n R T Gamma_inf * theta/(1-theta)

E0 DIVERGES as theta -> 1, so nRT*Gamma_inf is a prefactor, not a ceiling. The
library's own gibbs_marangoni_elasticity_langmuir docstring warns about exactly
this divergence. What is true instead:

  * a physically sensible E0 for a single-chain surfactant (~20-30 mN/m) pins
    theta to ~0.6-0.8, NOT ~0.95;
  * a supplied isotherm implying theta > 0.9 sits in the regime where the ideal
    equation of state diverges unphysically, so its K is too large.

Both checks below encode THAT, and they are the cheap insurance that caught the
supplied CTAB/SDS isotherm constants.
"""

from __future__ import annotations

import pytest

from research.laod_tau2.systems import SYSTEMS

RT = 8.314462618 * 298.15


def e0_mN_m(n_factor: float, gamma_inf: float, theta: float) -> float:
    rt = n_factor * RT * gamma_inf * 1e3
    return rt * theta / (1.0 - theta)


def theta_for_E0(n_factor: float, gamma_inf: float, E0: float) -> float:
    """Invert E0 = nRT*Ginf*theta/(1-theta) -> theta."""
    r = E0 / (n_factor * RT * gamma_inf * 1e3)
    return r / (1.0 + r)


def test_prefactor_magnitudes_are_sensible():
    """n*R*T*Gamma_inf should be 7-20 mN/m for these surfactants."""
    for s in SYSTEMS.values():
        rt = s.n_ions * RT * s.gamma_inf * 1e3
        assert 6.0 < rt < 22.0, (s.name, rt)


def test_literature_E0_implies_moderate_coverage_not_saturation():
    """E0 ~ 25 mN/m must correspond to theta ~0.55-0.80, not >0.9."""
    for s in SYSTEMS.values():
        th = theta_for_E0(s.n_ions, s.gamma_inf, 25.0)
        assert 0.5 < th < 0.85, (s.name, th)


def test_supplied_ctab_sds_isotherms_sit_in_the_divergent_regime():
    """Recorded defect: the supplied Frumkin K gives theta = 0.94-0.96 and
    E0 = 200-400 mN/m, i.e. deep in the unphysical saturated regime.
    Kept as a test so those numbers cannot be used without the flag being visible.
    """
    cases = {"CTAB": (0.92, 18.5), "SDS": (8.2, 1.25)}
    for name, (cmc, K) in cases.items():
        s = SYSTEMS[name]
        theta = K * cmc / (1.0 + K * cmc)
        assert theta > 0.90, (name, theta)
        assert e0_mN_m(s.n_ions, s.gamma_inf, theta) > 150.0


def test_corrected_isotherm_gives_physical_elasticity():
    """The K that yields E0 = 25 mN/m at the CMC, for reference."""
    expected = {"CTAB": 1.713, "SDS": 0.158}
    for name, K in expected.items():
        s = SYSTEMS[name]
        theta = K * s.cmc_mM / (1.0 + K * s.cmc_mM)
        assert 15.0 < e0_mN_m(s.n_ions, s.gamma_inf, theta) < 35.0, (name, theta)
