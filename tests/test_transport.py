"""Tests for P4 surfactant_transport_1d (Crank-Nicolson breakthrough)."""
import itertools
import math

import pytest

from surfactantkit.transport import surfactant_transport_1d


def _ogata_banks(x, t, v, D):
    """Textbook two-term analytic solution of the 1D advection-dispersion
    equation for a continuous step inlet on a semi-infinite domain."""
    return (
        0.5 * math.erfc((x - v * t) / (2.0 * math.sqrt(D * t)))
        + 0.5 * math.exp(v * x / D)
        * math.erfc((x + v * t) / (2.0 * math.sqrt(D * t)))
    )


def _retarded_ob(x, t, v, D, R):
    """Equilibrium-linear retardation form (constant R substitution)."""
    return (
        0.5 * math.erfc((R * x - v * t) / (2.0 * math.sqrt(D * t)))
        + 0.5 * math.exp(v * x / D)
        * math.erfc((R * x + v * t) / (2.0 * math.sqrt(D * t)))
    )


def test_tracer_matches_ogata_banks_and_converges():
    v, D, L, c0 = 1e-3, 5e-7, 0.2, 1.0
    errs = []
    for nx in (100, 400):
        r = surfactant_transport_1d(
            v, D, L, c0, n_cells=nx, t_max_s=100.0,
            n_times=nx, n_steps_per_interval=5,
        )
        e = max(
            abs(p - _ogata_banks(x, 100.0, v, D))
            for p, x in zip(r.concentration_profiles[-1], r.x_m, strict=True)
            if 0.005 < x < 0.19
        )
        errs.append(e)
        # sharp-front trapezoid quadrature limits the MEASURED balance at
        # coarse grids; scheme itself is conservative (sorbing runs hit 1e-12)
        assert r.converged and r.mass_balance_residual <= 1e-2
    # coarse/fine error ratio ~ O(dx^2): at least 2.5x improvement
    assert errs[0] / errs[1] >= 2.5
    assert errs[1] <= 1e-2  # absolute agreement with the textbook solution


def test_low_concentration_langmuir_matches_retarded_analytic():
    """bc << 1: Langmuir -> linear Kd = qmax*b, solution = R-form OB."""
    v, D, L = 1e-3, 5e-7, 0.2
    qmax, b, R = 0.5, 2.0, 1.0 + 0.5 * 2.0
    c0 = 1e-3
    r = surfactant_transport_1d(
        v, D, L, c0, langmuir_capacity_mol_m3=qmax,
        langmuir_affinity_m3_per_mol=b,
        n_cells=200, t_max_s=600.0, n_times=600, n_steps_per_interval=5,
    )
    e = max(
        abs(p / c0 - _retarded_ob(x, 600.0, v, D, R))
        for p, x in zip(r.concentration_profiles[-1], r.x_m, strict=True)
        if x < 0.9 * L
    )
    assert e <= 1e-6
    assert r.converged and r.mass_balance_residual <= 1e-6


def test_langmuir_full_strength_balance_order_and_monotonicity():
    v, D, L, c0 = 1e-3, 5e-7, 0.2, 1.0
    tracer = surfactant_transport_1d(
        v, D, L, c0, n_cells=400, t_max_s=400.0,
        n_times=400, n_steps_per_interval=5,
    )
    sorbing = surfactant_transport_1d(
        v, D, L, c0, langmuir_capacity_mol_m3=0.5,
        langmuir_affinity_m3_per_mol=2.0,
        n_cells=200, t_max_s=800.0, n_times=400, n_steps_per_interval=5,
    )
    # discrete injected = stored + outflux
    assert sorbing.mass_balance_residual <= 1e-6
    assert tracer.mass_balance_residual <= 1e-2
    # breakthrough is monotone and completes
    bf = sorbing.breakthrough_fraction
    assert all(y >= x - 1e-12 for x, y in itertools.pairwise(bf))
    assert bf[-1] == pytest.approx(1.0, abs=1e-6)
    # adsorption delays the front
    assert tracer.breakthrough_t50_s is not None
    assert sorbing.breakthrough_t50_s is not None
    assert sorbing.breakthrough_t50_s > tracer.breakthrough_t50_s
    # no unclamped undershoot for these smooth parameters
    assert sorbing.min_concentration >= -1e-6
    assert sorbing.converged


def test_outlet_t50_near_pore_volume_for_tracer():
    """Plug-flow-ish: t50 ~ L/v (dispersion only shifts it slightly)."""
    v, D, L, c0 = 1e-3, 5e-7, 0.2, 1.0
    r = surfactant_transport_1d(
        v, D, L, c0, n_cells=400, t_max_s=400.0,
        n_times=400, n_steps_per_interval=5,
    )
    assert r.breakthrough_t50_s == pytest.approx(L / v, rel=0.05)


def test_high_cell_peclet_discloses_oscillation_risk():
    v, D, L = 1e-3, 1e-7, 0.2  # Pe_cell with 20 cells = 20
    r = surfactant_transport_1d(
        v, D, L, 1.0, n_cells=20, t_max_s=50.0,
        n_times=10, n_steps_per_interval=10,
    )
    assert "oscillation" in r.caveat
    assert r.params["pe_cell"] > 2


def test_non_adsorbing_caveat_and_schema_fields():
    r = surfactant_transport_1d(
        1e-3, 5e-7, 0.2, 1.0, n_cells=50, t_max_s=50.0,
        n_times=5, n_steps_per_interval=2,
    )
    assert "non-adsorbing tracer" in r.caveat
    assert "phenomenological" in r.caveat
    assert r.n_iterations >= 5 * 2
    assert r.method and "Crank-Nicolson" in r.method
    assert len(r.times_s) == 5
    assert len(r.concentration_profiles) == 5
    assert len(r.concentration_profiles[0]) == len(r.x_m)
    assert "pe_cell" in r.params and "breakthrough_pore_volumes_at_half" in r.params


def test_determinism():
    kw: dict = {"n_cells": 50, "t_max_s": 50.0, "n_times": 5, "n_steps_per_interval": 2}  # untyped on purpose
    a = surfactant_transport_1d(1e-3, 5e-7, 0.2, 1.0, **kw)
    b = surfactant_transport_1d(1e-3, 5e-7, 0.2, 1.0, **kw)
    assert a.concentration_profiles == b.concentration_profiles
    assert a.breakthrough_fraction == b.breakthrough_fraction
    assert a.mass_balance_residual == b.mass_balance_residual


def test_fail_fast_names_every_bad_field():
    good: dict = {"velocity_m_s": 1e-3, "dispersion_m2_s": 5e-7,  # untyped on purpose
                       "column_length_m": 0.2, "inlet_conc_mol_m3": 1.0}
    with pytest.raises(ValueError, match="velocity_m_s"):
        surfactant_transport_1d(**{**good, "velocity_m_s": -1.0})
    with pytest.raises(ValueError, match="dispersion_m2_s"):
        surfactant_transport_1d(**{**good, "dispersion_m2_s": 0.0})
    with pytest.raises(ValueError, match="column_length_m"):
        surfactant_transport_1d(**{**good, "column_length_m": 0.0})
    with pytest.raises(ValueError, match="inlet_conc_mol_m3"):
        surfactant_transport_1d(**{**good, "inlet_conc_mol_m3": -0.5})
    with pytest.raises(ValueError, match="langmuir_capacity_mol_m3"):
        surfactant_transport_1d(**{**good, "langmuir_capacity_mol_m3": -1.0})
    with pytest.raises(ValueError, match="langmuir_affinity_m3_per_mol"):
        surfactant_transport_1d(**{**good, "langmuir_affinity_m3_per_mol": -1.0})
    with pytest.raises(ValueError, match="n_cells"):
        surfactant_transport_1d(**{**good, "n_cells": 5})
    with pytest.raises(ValueError, match="t_max_s"):
        surfactant_transport_1d(**{**good, "t_max_s": 0.0})
    with pytest.raises(ValueError, match="n_times"):
        surfactant_transport_1d(**{**good, "n_times": 1})
    with pytest.raises(ValueError, match="n_steps_per_interval"):
        surfactant_transport_1d(**{**good, "n_steps_per_interval": 0})
