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


# --------------------------------------------------------------------------
# Section 13 regressions. These lock in the generality result and the two
# self-corrections so they cannot silently drift back.
# --------------------------------------------------------------------------

def test_e0_anchored_isotherms_are_physical():
    """Inverting E0 -> theta -> K must give moderate coverage for both ionics."""
    from research.laod_tau2.study_ionic_generality import make_params

    for name, expect in (("CTAB", (0.65, 0.80)), ("SDS", (0.55, 0.75))):
        _, info = make_params(name, 1.5, 0.065)
        assert expect[0] < info["theta"] < expect[1], (name, info["theta"])


def test_ctab_is_not_detectable_inside_the_pat_band():
    """The section 13.2 negative: CTAB's I2 stays below the noise floor at 1 Hz."""
    import numpy as np

    from research.laod_tau2 import harmonics as H
    from research.laod_tau2.lae_data import NOISE_FLOOR_RATIO
    from research.laod_tau2.model import run_oscillation
    from research.laod_tau2.study_ionic_generality import make_params

    p, _ = make_params("CTAB", 1.5, 0.065)
    r = run_oscillation(p, 6.28, 0.20, n_cycles=16)
    t = np.asarray(r.t)
    i2 = H.harmonics(t, r.gamma_mN_m, 6.28)[2]
    assert i2 < NOISE_FLOOR_RATIO[1], i2


def test_c12e8_window_optimum_is_robust_to_the_e0_anchor():
    """Section 13.5: the 1.5x CMC optimum must hold across the E0 anchors.

    Reads the recorded sweep rather than re-running it (the sweep takes minutes).
    """
    import json
    from pathlib import Path

    path = Path(__file__).resolve().parents[1] / "e0_sensitivity.json"
    if not path.exists():
        import pytest

        pytest.skip("e0_sensitivity.json not generated")
    data = json.loads(path.read_text())
    assert len(data["rows"]) == 4
    for row in data["rows"]:
        assert row["best_mult"] == 1.5, (row["E0"], row["best_mult"])
        assert 1.5 in row["usable_mult"]
        assert row["best_dlag"] > 15.0


def test_dlag_label_is_per_factor_four_not_per_doubling():
    """Section 13.6: the recorded field spans tau2/2 -> 2*tau2 = 4x, and the
    directly measured per-doubling values differ from it by roughly 1.7x."""
    import json
    from pathlib import Path

    path = Path(__file__).resolve().parents[1] / "e0_sensitivity.json"
    if not path.exists():
        import pytest

        pytest.skip("e0_sensitivity.json not generated")
    data = json.loads(path.read_text())
    row = next(r for r in data["rows"] if abs(r["E0"] - 103.5) < 0.5)
    w = {x["mult"]: x for x in row["window"]}
    per_4x = w[1.5]["dlag_per_doubling"]
    # directly measured per-doubling near the operating point is 7.3-10.2 deg
    assert 15.0 < per_4x < 20.0
    assert abs(per_4x / 2.0 - 8.7) < 3.0


# --------------------------------------------------------------------------
# Section 13.10 regressions: the Frumkin isotherm sign fix and mass balance.
# These exist because a wrong-sign residual silently produced flat gamma(t),
# which was then misdiagnosed as solver stiffness.
# --------------------------------------------------------------------------

def test_frumkin_theta_matches_independent_inversion():
    """theta0 must satisfy the Frumkin isotherm it claims to solve."""
    import math

    from research.laod_tau2.model import ModelParams
    from research.laod_tau2.systems import SYSTEMS

    s = SYSTEMS["CTAB"]
    c = 1.5 * s.cmc_mM
    cases = {0.7: (0.690, 0.724), 2.0: (0.1246, 0.821),
             4.0: (0.004763, 0.892), 7.8: (4.86e-06, 0.940)}
    for a, (K, expect) in cases.items():
        p = ModelParams(c0=c, K=K, gamma_inf=s.gamma_inf, D=s.D, n_factor=2.0,
                        capacity=1.0, frumkin_a=a, temperature_K=298.15)
        th = p.theta0
        assert abs(th - expect) < 2e-3, (a, th, expect)
        residual = math.log(th / (1 - th)) - 2 * a * th - math.log(K * c)
        assert abs(residual) < 1e-9, (a, residual)


def test_langmuir_theta_is_unchanged_by_the_fix():
    """a = 0 must reduce exactly to K c/(1 + K c)."""
    from research.laod_tau2.model import ModelParams
    from research.laod_tau2.systems import SYSTEMS

    s = SYSTEMS["C12E8"]
    c = 1.5 * s.cmc_mM
    p = ModelParams(c0=c, K=140.0, gamma_inf=s.gamma_inf, D=s.D,
                    capacity=0.0, tau2=4.0, temperature_K=298.15)
    assert abs(p.theta0 - 140.0 * c / (1 + 140.0 * c)) < 1e-12


def test_mass_is_conserved_at_the_working_point():
    """Interior surfactant must not drift: the exchange term is conservative."""
    import numpy as np

    from research.laod_tau2.model import ModelParams, run_oscillation
    from research.laod_tau2.systems import SYSTEMS

    s = SYSTEMS["C12E8"]
    c = 1.5 * s.cmc_mM
    p = ModelParams(c0=c, K=140.0, gamma_inf=s.gamma_inf, D=s.D, n_factor=1.0,
                    capacity=s.capacity(c), tau2=4.0, temperature_K=298.15)
    r = run_oscillation(p, 6.28, 0.20, n_cycles=8)
    assert np.isfinite(r.mass_error)
    assert r.mass_error < 1e-2, r.mass_error


# --------------------------------------------------------------------------
# Section 18 regressions, from the external Paper 1 review.
# --------------------------------------------------------------------------

def test_to_model_refuses_to_return_an_unrunnable_object():
    """SurfactantSystem.to_model must raise rather than emit K=None (review #4)."""
    from research.laod_tau2.systems import SYSTEMS

    s = SYSTEMS["C12E8"]
    with pytest.raises(ValueError, match="K must be supplied"):
        s.to_model(1.5 * s.cmc_mM)
    # and must work when K is supplied
    p = s.to_model(1.5 * s.cmc_mM, K=140.0)
    assert p.K == 140.0
    assert 0.0 < p.theta0 < 1.0


def test_relative_se_is_not_divided_by_log():
    """The v3 estimator reported SE(ln tau2)/|ln tau2|, which is not a relative SE.

    Review #6. Guard: at tau2 = 4 s the wrong formula would understate by ln(4).
    """
    import math

    import numpy as np
    from scipy.optimize import least_squares

    # Toy exact-log fit with known noise scale: x = ln(tau2).
    truth = math.log(4.0)
    x = np.linspace(-2, 2, 40)
    y = np.full_like(x, 0.0) + 0.01

    def res(p):
        return np.full_like(x, p[0] - truth) - y

    sol = least_squares(res, x0=[0.0], xtol=1e-12, ftol=1e-12)
    se = float(np.sqrt(np.diag(np.linalg.inv(sol.jac.T @ sol.jac) * 1.0))[0])
    correct = se                      # relative SE of tau2 for small errors
    wrong = se / abs(sol.x[0])        # the withdrawn formula
    assert correct > wrong * 1.3, (correct, wrong)
