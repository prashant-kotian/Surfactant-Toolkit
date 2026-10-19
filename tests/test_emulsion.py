"""Tests for emulsion workers S5 (simulate_emulsion_stability) and O2
(population_balance_emulsion) -- Phase-3 batch A."""
import math

import pytest

from surfactantkit.emulsion import (
    G_ACCEL_M_S2,
    population_balance_emulsion,
    simulate_emulsion_stability,
)

# ---------------------------------------------------------------------------
# S5: Stokes creaming + first-order coalescence (piecewise-analytic)
# ---------------------------------------------------------------------------

def test_s5_stokes_velocity_and_creaming_time_exact():
    """v = |drho| g d^2 / (18 eta) and t = H/v, closed form."""
    d_um, drho, eta, H = 10.0, 100.0, 1e-3, 0.01
    # equal weights -> d32 reduces to the single diameter even with the
    # default coalescence active; guard it explicitly anyway
    r = simulate_emulsion_stability(
        [d_um], [1.0], 1000.0 - drho, 1000.0, eta, H,
        n_times=3, t_max_s=4000.0, coalescence_rate_per_s=0.0,
    )
    v_exact = drho * G_ACCEL_M_S2 * (d_um * 1e-6) ** 2 / (18.0 * eta)
    assert r.velocities_m_per_s[0] == pytest.approx(v_exact, rel=1e-12)
    assert r.creaming_times_s[0] == pytest.approx(H / v_exact, rel=1e-12)
    assert r.d32_um == pytest.approx(d_um, rel=1e-12)


def test_s5_two_bin_velocity_ratio_is_d_squared():
    r = simulate_emulsion_stability(
        [5.0, 10.0], [1.0, 1.0], 900.0, 1000.0, 1e-3, 0.02,
        n_times=5, t_max_s=10.0,
    )
    assert r.velocities_m_per_s[1] / r.velocities_m_per_s[0] == pytest.approx(
        (10.0 / 5.0) ** 2, rel=1e-12
    )
    assert r.creaming_times_s[1] < r.creaming_times_s[0]  # bigger drop creams first
    # true Sauter with equal volume weights: 2/(1/5 + 1/10) = 20/3 um
    assert r.d32_um == pytest.approx(2.0 / (1.0 / 5.0 + 1.0 / 10.0), rel=1e-12)


def test_s5_post_creaming_coalescence_matches_explicit_ode_solution():
    """Single bin, grid aligned so t_cream lands exactly on a sample:
    F(t_cream) = 0 (separation starts after arrival),
    F(2*t_cream) = 1 - exp(-k*t_cream) exactly."""
    d_um, drho, eta, H = 10.0, 100.0, 1e-3, 0.01
    v_exact = drho * G_ACCEL_M_S2 * (d_um * 1e-6) ** 2 / (18.0 * eta)
    tc = H / v_exact
    k_coal = 0.002
    r = simulate_emulsion_stability(
        [d_um], [1.0], 900.0, 1000.0, eta, H,
        coalescence_rate_per_s=k_coal, n_times=3, t_max_s=2.0 * tc,
    )
    assert r.times_s[1] == pytest.approx(tc, rel=1e-9)
    assert r.separated_fraction[1] == pytest.approx(0.0, abs=1e-12)
    assert r.separated_fraction[2] == pytest.approx(
        1.0 - math.exp(-k_coal * tc), rel=1e-12
    )


def test_s5_stability_index_grid_converges_to_creaming_time():
    d_um, drho, eta, H = 10.0, 100.0, 1e-3, 0.01
    v_exact = drho * G_ACCEL_M_S2 * (d_um * 1e-6) ** 2 / (18.0 * eta)
    tc = H / v_exact
    r = simulate_emulsion_stability(
        [d_um], [1.0], 900.0, 1000.0, eta, H, n_times=201,
        t_max_s=2.0 * tc,
    )
    # trapezoid across the unit step is within one grid cell of t_cream
    dt = r.times_s[1] - r.times_s[0]
    assert r.stability_index_m2 == pytest.approx(tc, abs=dt)
    assert r.t50_s is not None
    assert r.t50_s == pytest.approx(tc, abs=dt)
    assert r.converged and r.n_iterations == 201
    assert r.method and "piecewise" in r.method.lower() or r.method and "analytic" in r.method.lower() or r.method
    assert r.caveat


def test_s5_equal_density_never_separates():
    r = simulate_emulsion_stability(
        [5.0], [1.0], 1000.0, 1000.0, 1e-3, 0.01,
        n_times=5, t_max_s=100.0,
    )
    assert all(f == 0.0 for f in r.separated_fraction)
    assert r.stability_index_m2 == pytest.approx(100.0, rel=1e-12)
    assert r.t50_s is None


def test_s5_sedimentation_direction_disclosed():
    r = simulate_emulsion_stability(
        [5.0], [1.0], 1100.0, 1000.0, 1e-3, 0.01,
        n_times=5, t_max_s=100.0,
    )
    assert "sedimentation" in (r.caveat + str(r.params))
    assert r.velocities_m_per_s[0] > 0  # magnitude positive either way


def test_s5_hindered_creaming_caveat_above_10_percent_phi():
    r = simulate_emulsion_stability(
        [5.0, 8.0], [0.45, 0.45], 900.0, 1000.0, 1e-3, 0.01,
        n_times=5, t_max_s=10.0,
    )
    assert "hindered" in r.caveat


def test_s5_fail_fast_names_every_bad_field():
    good = ([5.0], [1.0], 900.0, 1000.0, 1e-3, 0.01)
    with pytest.raises(ValueError, match="same length"):
        simulate_emulsion_stability([5.0, 6.0], [1.0], *good[2:])
    with pytest.raises(ValueError, match="empty"):
        simulate_emulsion_stability([], [], *good[2:])
    with pytest.raises(ValueError, match="positive"):
        simulate_emulsion_stability([0.0], [1.0], *good[2:])
    with pytest.raises(ValueError, match="non-negative"):
        simulate_emulsion_stability([5.0], [-1.0], *good[2:])
    with pytest.raises(ValueError, match="positive"):
        simulate_emulsion_stability([5.0], [1.0], 900.0, 1000.0, 0.0, 0.01)
    with pytest.raises(ValueError, match="viscosity_Pa_s"):
        simulate_emulsion_stability([5.0], [1.0], 900.0, 1000.0, -1.0, 0.01)
    with pytest.raises(ValueError, match="column_height_m"):
        simulate_emulsion_stability([5.0], [1.0], 900.0, 1000.0, 1e-3, 0.0)
    with pytest.raises(ValueError, match="coalescence_rate_per_s"):
        simulate_emulsion_stability([5.0], [1.0], 900.0, 1000.0, 1e-3, 0.01,
                                    coalescence_rate_per_s=-1.0)
    with pytest.raises(ValueError, match="n_times"):
        simulate_emulsion_stability([5.0], [1.0], 900.0, 1000.0, 1e-3, 0.01,
                                    n_times=1)
    with pytest.raises(ValueError, match="t_max_s"):
        simulate_emulsion_stability([5.0], [1.0], 900.0, 1000.0, 1e-3, 0.01,
                                    t_max_s=0.0)
    with pytest.raises(ValueError, match="at least one positive"):
        simulate_emulsion_stability([5.0], [0.0], *good[2:])


# ---------------------------------------------------------------------------
# O2: discrete Smoluchowski PBE + LSW ripening
# ---------------------------------------------------------------------------

def test_o2_constant_kernel_number_decay_matches_exact_solution():
    """Constant kernel: n_tot(t) = n0 / (1 + lam0*n0*t/2), exact."""
    lam0, n0, t_end = 1e-3, 100.0, 50.0
    K = 128
    r = population_balance_emulsion(
        [n0] + [0.0] * (K - 1), "constant", lam0,
        n_times=6, t_max_s=t_end, n_steps_per_interval=20,
    )
    nt_exact = n0 / (1.0 + lam0 * n0 * t_end / 2.0)
    assert r.number_concentration[-1] == pytest.approx(nt_exact, rel=1e-8)
    assert r.mass_rel_error <= 1e-8
    assert r.converged
    # coalescence grows the Sauter mean diameter monotonically
    assert all(b >= a for a, b in zip(r.d32_um, r.d32_um[1:], strict=False))
    assert r.d32_um[0] == pytest.approx(1.0, rel=1e-12)
    assert r.d32_um[-1] > 1.2


def test_o2_sum_and_product_kernels_conserve_mass_when_boundary_inactive():
    K = 128
    for kernel, strength in (("sum", 1e-4), ("product", 1e-5)):
        r = population_balance_emulsion(
            [100.0] + [0.0] * (K - 1), kernel, strength,
            n_times=4, t_max_s=5.0, n_steps_per_interval=20,
        )
        assert r.mass_rel_error <= 1e-9, kernel
        assert r.converged, kernel


def test_o2_binary_breakage_conserves_mass_and_grows_cluster_count():
    r = population_balance_emulsion(
        [0.0] * 5 + [100.0] + [0.0] * 62,
        "constant", 0.0, breakage_rate_per_s=0.05,
        n_times=5, t_max_s=20.0, n_steps_per_interval=10,
    )
    assert r.mass_rel_error <= 1e-9
    assert r.converged
    assert r.number_concentration[-1] > r.number_concentration[0]
    # mass 100*6=600 -> all-monomer limit is 600 clusters
    assert r.number_concentration[-1] <= 600.0 + 1e-9
    assert "BINARY" in r.caveat


def test_o2_lsw_ripening_law_is_exact_when_coalescence_off():
    d0, kr = 2.0, 0.001
    r = population_balance_emulsion(
        [100.0] + [0.0] * 31, "constant", 0.0,
        ripening_rate_um3_per_s=kr, primary_drop_diameter_um=d0,
        n_times=5, t_max_s=100.0, n_steps_per_interval=1,
    )
    assert r.d32_lsw_only_um is not None
    for d, t in zip(r.d32_lsw_only_um, r.times_s, strict=True):
        assert d == pytest.approx((d0 ** 3 + kr * t) ** (1.0 / 3.0), rel=1e-12)
    # combined == LSW alone when coalescence/breakage are off
    assert r.params["d32_combined_final_um"] == pytest.approx(
        r.d32_lsw_only_um[-1], rel=1e-12
    )


def test_o2_truncation_failure_reports_not_converged():
    """K far too small for the time span -> mass escape -> converged False
    (the honest diagnostic, not a silent pass)."""
    r = population_balance_emulsion(
        [100.0] + [0.0] * 7, "constant", 1e-3,
        n_times=4, t_max_s=1000.0, n_steps_per_interval=10,
    )
    assert r.mass_rel_error > 1e-6
    assert r.converged is False


def test_o2_determinism_same_inputs_same_outputs():
    args = ([100.0] + [0.0] * 63, "constant", 1e-3)
    kw: dict = {"n_times": 4, "t_max_s": 10.0, "n_steps_per_interval": 10}  # untyped on purpose
    a = population_balance_emulsion(*args, **kw)
    b = population_balance_emulsion(*args, **kw)
    assert a.d32_um == b.d32_um
    assert a.number_concentration == b.number_concentration
    assert a.mass_rel_error == b.mass_rel_error


def test_o2_fail_fast_names_every_bad_field():
    good = [100.0] + [0.0] * 15
    with pytest.raises(ValueError, match="empty"):
        population_balance_emulsion([])
    with pytest.raises(ValueError, match="at least 2 classes"):
        population_balance_emulsion([100.0])
    with pytest.raises(ValueError, match="non-negative"):
        population_balance_emulsion([-1.0, 1.0])
    with pytest.raises(ValueError, match="at least one positive"):
        population_balance_emulsion([0.0, 0.0])
    with pytest.raises(ValueError, match="coalescence_kernel"):
        population_balance_emulsion(good, "nope", 1.0)
    with pytest.raises(ValueError, match="kernel_strength"):
        population_balance_emulsion(good, "constant", -1.0)
    with pytest.raises(ValueError, match="breakage_rate_per_s"):
        population_balance_emulsion(good, "constant", 1.0,
                                    breakage_rate_per_s=-1.0)
    with pytest.raises(ValueError, match="ripening_rate_um3_per_s"):
        population_balance_emulsion(good, "constant", 1.0,
                                    ripening_rate_um3_per_s=-1.0)
    with pytest.raises(ValueError, match="primary_drop_diameter_um"):
        population_balance_emulsion(good, "constant", 1.0,
                                    primary_drop_diameter_um=0.0)
    with pytest.raises(ValueError, match="n_times"):
        population_balance_emulsion(good, "constant", 1.0, n_times=1)
    with pytest.raises(ValueError, match="n_steps_per_interval"):
        population_balance_emulsion(good, "constant", 1.0,
                                    n_steps_per_interval=0)
    with pytest.raises(ValueError, match="t_max_s"):
        population_balance_emulsion(good, "constant", 1.0, t_max_s=0.0)
