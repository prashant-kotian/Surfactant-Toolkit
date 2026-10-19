"""Tests for C6 cloud_point_estimate (HLD-based, fully sourced)."""

import pytest

from surfactantkit.cloudpoint import (
    B_NACL_DEFAULT,
    CT_ETHOXYLATE_DEFAULT,
    EACN_APP_ETHOXYLATE,
    K_DEFAULT,
    cloud_point_estimate,
)


def test_closed_form_gold_recomputation():
    """Exact arithmetic gold: T = 25 - (bS - k*7.3 + Cc)/cT."""
    for n, eo in [(12, 5), (12, 6), (10, 6), (8, 6)]:
        r = cloud_point_estimate(n, eo)
        cc = 0.28 * n + 2.4 - eo
        expected = 25.0 - (0.0 - K_DEFAULT * 7.3 + cc) / 0.06
        assert r.cloud_point_C == pytest.approx(expected, rel=1e-12)
        assert r.cc == pytest.approx(cc, rel=1e-12)
        assert abs(r.hld_at_cloud_point) < 1e-9
        assert r.converged and r.n_iterations == 1


def test_c12e5_within_source_paper_accuracy_window():
    """Zarate-Munoz 2015: C12E5 cloud point 'close to 29 C' and HLD
    predictions typically within 5 C -> our estimate must land inside
    that window (25-34 C)."""
    r = cloud_point_estimate(12, 5)
    assert 25.0 <= r.cloud_point_C <= 34.0
    assert abs(r.cloud_point_C - 29.0) <= 5.0


def test_trends_more_eo_higher_cp_more_tail_lower_cp():
    """Hydration trend: CP rises with EO, falls with tail length."""
    cp_e5 = cloud_point_estimate(12, 5).cloud_point_C
    cp_e6 = cloud_point_estimate(12, 6).cloud_point_C
    cp_e8 = cloud_point_estimate(12, 8).cloud_point_C
    assert cp_e5 < cp_e6 < cp_e8
    cp_c8 = cloud_point_estimate(8, 6).cloud_point_C
    cp_c12 = cloud_point_estimate(12, 6).cloud_point_C
    assert cp_c8 > cp_c12


def test_nacl_depresses_cloud_point():
    """The source paper validated HLD prediction of NaCl CP depression:
    adding salt must lower the estimate (b*S > 0 term)."""
    pure = cloud_point_estimate(12, 6)
    salty = cloud_point_estimate(12, 6, salinity_g_per_100ml=1.0)
    assert salty.cloud_point_C < pure.cloud_point_C
    # exact closed-form shift: dT = -b*S/cT
    shift = B_NACL_DEFAULT * 1.0 / CT_ETHOXYLATE_DEFAULT
    assert pure.cloud_point_C - salty.cloud_point_C == pytest.approx(
        shift, rel=1e-12
    )


def test_k_sensitivity_band_and_constant_echo():
    r_lo = cloud_point_estimate(12, 6, k=0.15)
    r_hi = cloud_point_estimate(12, 6, k=0.17)
    assert r_lo.cloud_point_C < r_hi.cloud_point_C
    # full sourced band spans only a few K for C12E6
    assert r_hi.cloud_point_C - r_lo.cloud_point_C < 5.0
    r = cloud_point_estimate(12, 6)
    assert r.params["k"] == K_DEFAULT
    assert r.params["eacn_app"] == EACN_APP_ETHOXYLATE == 7.3
    assert r.params["b_NaCl"] == B_NACL_DEFAULT == 0.13
    assert "Zarate" in r.params["sources"] and "Ghayour" in r.params["sources"]
    assert "Acosta 2008" in r.params["sources"]


def test_caveats_and_extrapolation_flag():
    r = cloud_point_estimate(12, 6)
    assert "polydisperse" in r.caveat and "out of scope" in r.caveat
    assert "within 5 C" in r.caveat
    assert "extrapolation" not in r.caveat
    r12 = cloud_point_estimate(12, 12)
    assert "extrapolation" in r12.caveat


def test_fail_fast_names_every_bad_field():
    with pytest.raises(TypeError, match="n_carbons"):
        cloud_point_estimate("12", 6)  # type: ignore[reportArgumentType]
    with pytest.raises(ValueError, match="n_carbons"):
        cloud_point_estimate(2, 6)
    with pytest.raises(TypeError, match="n_eo"):
        cloud_point_estimate(12, 6.0)  # type: ignore[reportArgumentType]
    with pytest.raises(ValueError, match="n_eo"):
        cloud_point_estimate(12, 0)
    with pytest.raises(ValueError, match="salinity_g_per_100ml"):
        cloud_point_estimate(12, 6, salinity_g_per_100ml=-0.1)
    with pytest.raises(TypeError, match="salinity_g_per_100ml"):
        cloud_point_estimate(12, 6, salinity_g_per_100ml="1")  # type: ignore[reportArgumentType]
    with pytest.raises(ValueError, match="k"):
        cloud_point_estimate(12, 6, k=0.01)
    with pytest.raises(TypeError, match="k"):
        cloud_point_estimate(12, 6, k="0.16")  # type: ignore[reportArgumentType]


def test_determinism():
    a = cloud_point_estimate(12, 6, salinity_g_per_100ml=0.5)
    b = cloud_point_estimate(12, 6, salinity_g_per_100ml=0.5)
    assert a.cloud_point_C == b.cloud_point_C
    assert a.params == b.params
