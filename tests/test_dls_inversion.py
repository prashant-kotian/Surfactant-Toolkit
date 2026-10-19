"""Tests for surfactantkit.curve_analysis.dls_size_distribution (worker P3).

Golds: (a) scattering-vector geometry vs hand value for 633 nm / 90 deg /
water; (b) monodisperse synthetic correlation recovers the input diameter;
(c) bimodal synthetic recovers both modes with correct mass split.
"""
import math

import numpy as np
import pytest

from surfactantkit.curve_analysis import (
    _q_from_geometry,
    dls_size_distribution,
)

KT = 1.380649e-23 * 298.15
ETA = 0.89e-3
Q = _q_from_geometry(632.8, 90.0, 1.331)


def _gamma(d_nm: float) -> float:
    D_si = KT / (3 * math.pi * ETA * d_nm * 1e-9)
    return D_si * 1e18 * Q * Q


def test_q_geometry():
    # standard 633 nm, 90 deg, water: q ~ 0.0188 nm^-1
    assert Q == pytest.approx(0.0188, abs=0.0004)
    # 60 deg also computable
    q60 = _q_from_geometry(632.8, 60.0, 1.331)
    assert q60 < Q
    assert q60 == pytest.approx(Q * math.sin(math.radians(30)) / math.sin(math.radians(45)),
                                rel=1e-9)


def test_monodisperse_recovers_diameter():
    tau = np.logspace(-4, 0, 40)
    rng = np.random.default_rng(1)
    g = np.exp(-_gamma(50.0) * tau) + rng.normal(0, 0.003, tau.size)
    r = dls_size_distribution(g.tolist(), tau.tolist(), q_nm_inv=Q)
    d = np.array(r.diameter_nm)
    w = np.array(r.intensity_fraction)
    peak = d[int(np.argmax(w))]
    assert peak == pytest.approx(50.0, rel=0.15)
    assert r.dh_intensity_mean_nm == pytest.approx(50.0, rel=0.15)
    assert r.relative_residual < 0.05
    assert r.pdi_from_distribution < 0.1
    assert r.converged


def test_bimodal_recovers_both_modes_and_mass_split():
    tau = np.logspace(-4, 0, 40)
    rng = np.random.default_rng(1)
    g = (0.5 * np.exp(-_gamma(10.0) * tau)
         + 0.5 * np.exp(-_gamma(100.0) * tau)
         + rng.normal(0, 0.003, tau.size))
    r = dls_size_distribution(g.tolist(), tau.tolist(), q_nm_inv=Q)
    d = np.array(r.diameter_nm)
    w = np.array(r.intensity_fraction)
    small, large = d < 30.0, d > 30.0
    # mass split: true intensities 0.5 / 0.5
    assert 0.32 < w[small].sum() < 0.68
    assert 0.32 < w[large].sum() < 0.68
    # each mode has its own local maximum in the right window
    w_s, w_l = w[small].copy(), w[large].copy()
    # small-mode peak may pile onto the lowest resolvable grid node
    # (ill-posed edge effect): accept the resolvable window
    assert 5.0 <= d[small][int(np.argmax(w_s))] < 18.0
    assert 60.0 < d[large][int(np.argmax(w_l))] < 150.0
    # deepest point strictly between the two peaks must be far below them
    i_s = int(np.flatnonzero(small)[np.argmax(w_s)])
    i_l = int(np.flatnonzero(large)[np.argmax(w_l)])
    lo, hi = min(i_s, i_l), max(i_s, i_l)
    valley = w[lo + 1:hi].min()
    assert valley < 0.4 * min(w_s.max(), w_l.max())


def test_deterministic():
    tau = np.logspace(-4, 0, 24)
    g = np.exp(-_gamma(80.0) * tau)
    a = dls_size_distribution(g.tolist(), tau.tolist(), q_nm_inv=Q, n_gamma=40)
    b = dls_size_distribution(g.tolist(), tau.tolist(), q_nm_inv=Q, n_gamma=40)
    assert a == b


def test_explicit_lambda_and_q():
    tau = np.logspace(-4, 0, 30)
    g = np.exp(-_gamma(30.0) * tau)
    r = dls_size_distribution(g.tolist(), tau.tolist(), q_nm_inv=Q,
                              lambda_reg=1e-3)
    assert r.lambda_reg == 1e-3
    d = np.array(r.diameter_nm)
    w = np.array(r.intensity_fraction)
    assert d[int(np.argmax(w))] == pytest.approx(30.0, rel=0.2)


def test_validation_errors():
    tau = np.logspace(-4, 0, 20).tolist()
    g = np.exp(-_gamma(50.0) * np.asarray(tau)).tolist()
    with pytest.raises(ValueError, match="equal length"):
        dls_size_distribution([0.5] * 5, tau)
    with pytest.raises(ValueError, match=">= 8"):
        dls_size_distribution([0.9, 0.8, 0.7], [1e-3, 1e-2, 1e-1])
    tau_bad = list(tau)
    tau_bad[0] = 0.0
    with pytest.raises(ValueError, match="positive and finite"):
        dls_size_distribution([0.5] * len(tau_bad), tau_bad)
    with pytest.raises(ValueError, match="flat"):
        dls_size_distribution([0.5] * 12, np.logspace(-4, 0, 12).tolist())
    with pytest.raises(ValueError, match="n_gamma"):
        dls_size_distribution(g, tau, n_gamma=5)
    with pytest.raises(ValueError, match="lambda_reg"):
        dls_size_distribution(g, tau, lambda_reg=0.0)
    with pytest.raises(ValueError, match="viscosity"):
        dls_size_distribution(g, tau, viscosity_mPa_s=0.0)


def test_result_metadata():
    tau = np.logspace(-4, 0, 30)
    g = np.exp(-_gamma(45.0) * tau)
    r = dls_size_distribution(g.tolist(), tau.tolist(), q_nm_inv=Q)
    assert r.method and r.params["angle_deg"] == 90.0
    assert "ILL-POSED" in r.caveat
    assert r.n_grid == 60
    assert r.lambda_reg > 0
