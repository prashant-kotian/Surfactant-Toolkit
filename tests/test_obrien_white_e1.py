"""Tests for surfactantkit.electrostatics.obrien_white_numerical (worker E1).

Gold: Huckel limit, monotone ka -> 1 approach, exact sign symmetry,
relaxation gap sign at high zeta (no-relaxation quadrature must EXCEED the
relaxation-inclusive OHW83 approximation), ohw83 cross-check with NaN below
ka=10, bit-identical reruns, and the ka=100 / high-zeta physical sanity that
the earlier always-on Schur implementation failed.
"""
import itertools
import math

import pytest

from surfactantkit.electrostatics import OBrienWhiteResult, obrien_white_numerical


def test_huckel_limit_and_monotone_ka_sweep():
    # zeta = 5 mV: linear regime, microelectrophoretics must recover Huckel
    vals = {}
    for ka in (0.1, 0.5, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0):
        r = obrien_white_numerical(5.0, ka)
        assert isinstance(r, OBrienWhiteResult)
        vals[ka] = r.mu_hat
    # Huckel: mu_hat -> 2/3 as ka -> 0
    assert vals[0.1] == pytest.approx(2.0 / 3.0, rel=0.01)
    assert abs(vals[0.1] - r.mu_hat_huckel) / r.mu_hat_huckel < 0.01
    # monotone increasing in ka, approaching 1 (Henry/Smoluchowski limit)
    ks = sorted(vals)
    seq = [vals[k] for k in ks]
    assert all(b > a for a, b in itertools.pairwise(seq)), seq
    assert vals[300.0] > 0.98
    # gold anchors (fixed solver defaults, n_grid=900)
    assert vals[1.0] == pytest.approx(0.684, abs=0.002)
    assert vals[10.0] == pytest.approx(0.8346, abs=0.003)


def test_exact_sign_symmetry():
    for ka in (1.0, 10.0):
        rp = obrien_white_numerical(50.0, ka)
        rm = obrien_white_numerical(-50.0, ka)
        # exact antisymmetry (carried by the g field, no zeta_sign factor)
        assert rp.mu_hat + rm.mu_hat == pytest.approx(0.0, abs=1e-12)


def test_ka100_high_zeta_physical():
    # regression for the always-on Schur defect: returned -4.63 at ka=100
    r = obrien_white_numerical(100.0, 100.0)
    assert 0.0 < r.mu_hat < 1.0
    assert r.mu_hat == pytest.approx(0.9107, abs=0.005)
    rm = obrien_white_numerical(-100.0, 100.0)
    assert rm.mu_hat == pytest.approx(-r.mu_hat, abs=1e-12)
    # physical value moves monotonically with |zeta| in this regime
    r75 = obrien_white_numerical(75.0, 100.0)
    assert r75.mu_hat > r.mu_hat


def test_relaxation_gap_sign_and_ohw83_crosscheck():
    # ka=10, zeta=100 mV: default (relaxation OFF) quadrature must EXCEED
    # the relaxation-inclusive OHW83 approximation -- sign of the gap is the
    # physics statement, magnitude is disclosed as method error not an
    # error bar (memory: 0.6391 vs 0.5921, +7.9%).
    r = obrien_white_numerical(100.0, 10.0)          # relaxation=False default
    assert r.relaxation_group_gamma_bar is None or r.converged
    assert r.mu_hat > r.mu_hat_ohw83
    rel = (r.mu_hat - r.mu_hat_ohw83) / r.mu_hat_ohw83
    assert 0.0 < rel < 0.25
    assert r.ohw83_rel_diff > 0
    # OHW83 cross-check is NaN below ka=10 (asymptotic formula domain)
    r_lo = obrien_white_numerical(50.0, 3.0)
    assert math.isnan(r_lo.mu_hat_ohw83)


def test_relaxation_channel_runs_and_is_marked_experimental():
    r = obrien_white_numerical(100.0, 10.0, relaxation=True)
    assert r.converged and math.isfinite(r.mu_hat)
    assert r.params["relaxation"] is True
    # coupled channel is experimental: method string must say so
    assert "EXPERIMENTAL" in r.method or "experimental" in r.caveat
    assert "relaxation" in r.method or "relaxation" in str(r.params)


def test_deterministic_and_schema():
    a = obrien_white_numerical(50.0, 10.0)
    b = obrien_white_numerical(50.0, 10.0)
    assert a.mu_hat == b.mu_hat            # bit-identical rerun
    assert a.as_dict() == b.as_dict()
    for k in ("converged", "n_iterations", "method", "params", "caveat"):
        assert k in a.as_dict()
    assert a.converged and a.n_iterations >= 1
    assert a.params["relaxation"] is False  # safe default on the wire
    # units: micro-electrophoretic mobility, um cm / (V s)
    assert a.mu_um_cm_per_Vs == pytest.approx(a.mu_SI * 1e8, rel=1e-10)


def test_error_paths():
    for kwargs, pat in [
        ({"zeta_mV": 0.0, "kappa_a": 10.0}, "zeta"),
        ({"zeta_mV": 50.0, "kappa_a": 0.0}, "kappa"),
        ({"zeta_mV": 50.0, "kappa_a": -1.0}, "kappa"),
        ({"zeta_mV": 50.0, "kappa_a": 10.0, "transport_number_cation": 0.0},
         "transport"),
        ({"zeta_mV": 50.0, "kappa_a": 10.0, "n_grid": 10}, "n_grid"),
    ]:
        with pytest.raises(ValueError) as e:
            obrien_white_numerical(**kwargs)
        assert pat in str(e.value).lower()


def test_caveat_discloses_anomalous_surface_conductance_limit():
    """O'Brien-White 1978 structurally assumes zero ionic conduction
    inside the slipping plane and so underestimates zeta on real,
    highly-charged colloids where anomalous surface conductance is
    present -- Moncho, Martinez-Lopez & Hidalgo-Alvarez, Colloids Surf.
    A 192 (2001) 215-226, quantify this at up to ~35% in zeta at
    kappa_a~50. This is not a solver bug this toolkit can fix (it is a
    structural limit of the theory itself), so the caveat must disclose
    it on every call rather than let a user trust a bare number at face
    value for moderate kappa_a on a real colloid."""
    r = obrien_white_numerical(80.0, 50.0)
    assert "anomalous surface" in r.caveat
    assert "Moncho" in r.caveat
    assert "kappa_a" in r.caveat.lower() or "kappa_a~50" in r.caveat
