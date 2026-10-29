"""Regression tests for the 2026-10-30 surface-tension CMC fixes, found on a
46-isotherm catanionic dataset (amido-amine cationics + sodium lauroyl
glutamate) generated around known CMCs:

* rss/intersection: a 3-segment model whose "baseline" is steeper than its
  "decline" is no longer accepted (C12-Et+SLG x=0.6 read 1.68 mM against
  0.406 mM built in).
* method="intersection": regression-line intersection, continuous.
* method="tangential": the tangent is the steepest run of consecutive points
  below the plateau, not a prefix from the lowest concentration (3-10x low
  on curves with a knee near the solvent value).
* _fit_szyszkowski: the b search range scales with the data (b was pinned
  at 1e4 per mM for CMC ~0.01 mM).
* cmc_via_szyszkowski_plateau_intersection: Szyszkowski nonlinear curve x
  plateau line.
"""

import math

import pytest

from surfactantkit.curve_analysis import (
    _fit_szyszkowski,
    cmc_from_surface_tension_curve,
    cmc_via_szyszkowski_plateau_intersection,
)

# Means of three runs, C12-Et + SLG at x(cationic) = 0.6 and C14-s3OH-C14 + SLG at x = 0.6;
# each series was generated around CMC = C / (C/CMC) = 0.406017 mM and 0.028208 mM.
C12ET_X06_C = [0.00324, 0.004628, 0.006612, 0.009445, 0.013493, 0.019276, 0.027537, 0.039339, 0.056198,
               0.080283, 0.11469, 0.163843, 0.234061, 0.334373, 0.477676, 0.682394, 0.974848, 1.39264,
               1.989486, 2.842123, 4.060176]
C12ET_X06_G = [70.72, 68.25, 65.813, 63.413, 60.71, 58.31, 55.747, 53.217, 50.823, 48.357, 45.967, 43.493,
               40.767, 38.497, 37.36, 37.0, 37.12, 36.817, 36.907, 36.963, 37.17]
C14OH_X06_C = [0.000225, 0.000322, 0.000459, 0.000656, 0.000937, 0.001339, 0.001913, 0.002733, 0.003904,
               0.005578, 0.007968, 0.011383, 0.016262, 0.023231, 0.033187, 0.04741, 0.067729, 0.096755,
               0.138222, 0.19746, 0.282085]
C14OH_X06_G = [61.237, 59.373, 57.657, 55.79, 54.123, 52.14, 50.52, 48.727, 46.88, 45.103, 43.353, 41.493,
               39.813, 37.92, 37.25, 36.91, 37.07, 37.07, 37.013, 37.043, 36.97]


def _logspace(lo, hi, n):
    return [10 ** (math.log10(lo) + i * (math.log10(hi) - math.log10(lo)) / (n - 1)) for i in range(n)]


@pytest.mark.parametrize("c, g, built", [(C12ET_X06_C, C12ET_X06_G, 0.406017),
                                         (C14OH_X06_C, C14OH_X06_G, 0.028208)])
def test_curved_decline_is_not_split_into_a_steep_baseline(c, g, built):
    """Before the fix the 3-segment search called the real decline "baseline" and won the BIC:
    CMC +315% and +40%. Both now land within 5% of the CMC the series was built around."""
    for method in ("rss", "intersection"):
        r = cmc_from_surface_tension_curve(c, g, method=method)
        assert r.cmc_mM == pytest.approx(built, rel=0.05), method
        assert r.n_postmicellar_points >= 6


def test_accepted_baseline_is_always_flatter_than_the_decline():
    c, g = C12ET_X06_C, C12ET_X06_G
    r = cmc_from_surface_tension_curve(c, g)
    if r.n_baseline_points:
        lc = [math.log10(x) for x in c[: r.n_baseline_points]]
        base_slope = (g[r.n_baseline_points - 1] - g[0]) / (lc[-1] - lc[0])
        assert abs(base_slope) < abs(r.premicellar_slope_mN_per_m_per_log10C)


def test_intersection_hits_true_cmc_exactly_between_grid_points():
    """gamma = 50 - 10 log10 C up to CMC = 7 mM (not a grid point), flat above: the regression lines
    meet at exactly 7 mM, while rss can only report a midpoint between two measured points."""
    c = [0.5, 1, 2, 3, 4, 5, 6, 8, 10, 12, 14, 16, 20]
    plateau = 50 - 10 * math.log10(7)
    g = [50 - 10 * math.log10(x) if x <= 7 else plateau for x in c]
    r = cmc_from_surface_tension_curve(c, g, method="intersection")
    assert r.cmc_mM == pytest.approx(7.0, rel=1e-9)
    assert r.gamma_at_cmc_mN_per_m == pytest.approx(plateau, rel=1e-9)
    assert r.method.startswith("regression intersection")
    assert cmc_from_surface_tension_curve(c, g).cmc_mM != pytest.approx(7.0, rel=1e-3)


def test_tangential_ignores_the_knee_near_the_solvent_value():
    """Flat at 72 mN/m, then gamma = 72 - 15 log10(C/0.01) to CMC = 0.5 mM, flat above. The old prefix
    search fitted from the lowest concentration through the knee and read the CMC far too low; the
    steepest consecutive run is the decline itself, so the tangent meets the plateau exactly at 0.5."""
    c = _logspace(0.001, 5, 25)
    g = [72.0 if x <= 0.01 else (72 - 15 * math.log10(x / 0.01) if x <= 0.5 else 72 - 15 * math.log10(50))
         for x in c]
    r = cmc_from_surface_tension_curve(c, g, method="tangential")
    assert r.cmc_mM == pytest.approx(0.5, rel=1e-9)
    assert r.premicellar_slope_mN_per_m_per_log10C == pytest.approx(-15.0, rel=1e-9)


def test_fit_szyszkowski_recovers_a_large_b_for_a_low_cmc_surfactant():
    """b = 2e5 per mM (CMC ~0.01 mM): the old fixed 1e-4..1e4 range pinned b at 1e4."""
    c = _logspace(1e-5, 5e-3, 15)
    g = [72.0 - 4.0 * math.log(2e5 * x + 1.0) for x in c]
    g0, a, b, r2 = _fit_szyszkowski(c, g)
    assert b == pytest.approx(2e5, rel=1e-5)
    assert a == pytest.approx(4.0, rel=1e-6)
    assert g0 == pytest.approx(72.0, rel=1e-8)
    assert r2 == pytest.approx(1.0, abs=1e-10)


def test_szyszkowski_plateau_intersection_recovers_exact_cmc():
    """Exact Szyszkowski curve below CMC = 0.01 mM, flat at gamma(CMC) above: the fitted curve meets
    the plateau line at exactly 0.01 mM, and gamma_CMC is that plateau value."""
    cmc = 0.01
    gam = lambda x: 72.0 - 4.0 * math.log(2e5 * x + 1.0)
    c = _logspace(1e-5, 1e-1, 25)
    g = [gam(x) if x <= cmc else gam(cmc) for x in c]
    r = cmc_via_szyszkowski_plateau_intersection(c, g, system_type="nonionic")
    assert r.cmc_mM == pytest.approx(cmc, rel=1e-6)
    assert r.gamma_at_cmc_mN_per_m == pytest.approx(gam(cmc), rel=1e-9)
    assert r.b_fitted_per_mM == pytest.approx(2e5, rel=1e-5)
    assert r.gamma_max_from_fit_mol_per_m2 == pytest.approx(4.0 / (8.314462618 * 298.15 * 1000.0), rel=1e-6)


def test_szyszkowski_plateau_intersection_on_measured_curves():
    for c, g, built in ((C12ET_X06_C, C12ET_X06_G, 0.406017), (C14OH_X06_C, C14OH_X06_G, 0.028208)):
        r = cmc_via_szyszkowski_plateau_intersection(c, g)
        assert r.cmc_mM == pytest.approx(built, rel=0.06)
        assert r.gamma_max_from_fit_mol_per_m2 is None


def test_szyszkowski_plateau_intersection_reports_none_without_a_crossing():
    c = _logspace(0.01, 10, 12)
    g = [72.0 - 4.0 * math.log(50 * x + 1.0) for x in c]   # no plateau: keeps falling
    r = cmc_via_szyszkowski_plateau_intersection(c, g)
    assert r.cmc_mM is None and r.reason


def test_divalent_ionic_no_added_salt_uses_gibbs_prefactor_three():
    from surfactantkit.adsorption import gibbs_gamma_max
    one = gibbs_gamma_max(-10.0, "nonionic", 298.15)
    assert gibbs_gamma_max(-10.0, "divalent_ionic_no_added_salt", 298.15) == pytest.approx(one / 3.0, rel=1e-12)


def test_szyszkowski_plateau_intersection_rejects_unknown_system_type():
    with pytest.raises(ValueError, match="system_type"):
        cmc_via_szyszkowski_plateau_intersection(C12ET_X06_C, C12ET_X06_G, system_type="cationic-ish")
