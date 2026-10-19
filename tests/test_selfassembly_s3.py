"""Tests for surfactantkit.selfassembly (worker S3).

Gold (reduced-unit 2-bead Langevin model, seed-fixed):
- same seed -> bit-identical summary, different seed -> different trajectory;
- concentration trend (the CMC-like statement): at eps_tail=4 kT,
  c=0.15 M micellizes (monomer fraction < 1, micelle_count_mean >= 1,
  largest cluster grows) while c=0.002 M stays monomeric (monomer >= 0.95,
  no micelles), both from unbiased random starts;
- clusterd init: droplet survives above the CMC-like threshold;
- schema + refusals for degenerate setups.

Dynamics claim is qualitative (see result caveat): trends, not absolute
CMCs. Runtime ~10 s/run, keep the number of runs small.
"""
import pytest

from surfactantkit.selfassembly import (
    BrownianAssemblyResult,
    brownian_cmc_scan,
    brownian_self_assembly,
)

FAST = {"n_molecules": 60, "steps": 8000, "tail_attraction_kT": 4.0}


def test_determinism_and_schema():
    a = brownian_self_assembly(concentration_M=0.15, seed=7, **FAST)
    b = brownian_self_assembly(concentration_M=0.15, seed=7, **FAST)
    assert isinstance(a, BrownianAssemblyResult)
    assert a.as_dict() == b.as_dict()          # bit-identical same seed
    c = brownian_self_assembly(concentration_M=0.15, seed=8, **FAST)
    assert (c.mean_nagg, c.monomer_fraction) != (a.mean_nagg, a.monomer_fraction)
    d = a.as_dict()
    for k in ("converged", "n_iterations", "method", "params", "caveat",
              "nagg_distribution", "nagg_timeseries"):
        assert k in d
    # distribution promised by the module docstring + MCP description:
    # fractions of clusters sum to 1, JSON-string keys parse back to sizes
    dist = d["nagg_distribution"]
    assert dist and abs(sum(dist.values()) - 1.0) < 1e-12
    assert all(int(k) >= 1 for k in dist)
    # timeseries = per-frame mean cluster size (stationarity trace)
    assert len(d["nagg_timeseries"]) == d["n_frames"]
    assert d["n_iterations"] == 8000
    assert d["params"]["seed"] == 7
    assert d["params"]["force_cap_fraction_mean"] < 0.5
    assert "qualitative" in d["caveat"] or "reduced" in d["caveat"]


def test_concentration_trend_micelles_vs_monomers():
    # the physics S3 must show: above the CMC-like threshold molecules
    # aggregate, below it they stay monomeric
    hi = brownian_self_assembly(concentration_M=0.15, seed=1, **FAST)
    lo = brownian_self_assembly(concentration_M=0.002, seed=1, **FAST)
    # high concentration: micelles exist
    assert hi.micelle_count_mean >= 1.0
    assert hi.monomer_fraction < 0.9
    assert hi.largest_cluster_frac > 0.05
    # low concentration: essentially all monomers, no micelles
    assert lo.monomer_fraction >= 0.95
    assert lo.micelle_count_mean == 0.0
    # ordering
    assert hi.monomer_fraction < lo.monomer_fraction
    assert hi.mean_nagg > lo.mean_nagg
    assert hi.largest_cluster_frac > lo.largest_cluster_frac
    # cmc proxy only claimed when micelles coexist with monomers
    assert hi.params["cmc_proxy_valid"] is True
    assert 0 < hi.params["cmc_proxy_M"] <= 0.15
    assert lo.params["cmc_proxy_valid"] is False
    assert lo.params["cmc_proxy_M"] is None


def test_clustered_init_survives_above_threshold():
    r = brownian_self_assembly(concentration_M=0.15, seed=2,
                               init="clustered", **FAST)
    assert r.params["init"] == "clustered"
    # droplet survives: majority of molecules still in a large aggregate
    assert r.largest_cluster_frac > 0.3
    assert r.monomer_fraction < 0.6
    assert r.converged


def test_time_and_frames_reporting():
    r = brownian_self_assembly(concentration_M=0.15, seed=1,
                               steps=2000, n_molecules=40,
                               tail_attraction_kT=4.0, frame_every=200)
    assert r.n_frames >= 4
    assert r.time_tau == pytest.approx(2000 * 0.002)


def test_error_paths():
    calls = [
        ({"n_molecules": 10}, ">= 20"),
        ({"steps": 100}, ">= 500"),
        ({"concentration_M": 0.0}, "positive"),
        ({"concentration_M": -1.0}, "positive"),
        ({"tail_attraction_kT": -0.5}, ">= 0"),
        ({"spring_kT_per_sigma2": 0.0}, "positive"),
        ({"dt_tau": 0.0}, "positive"),
        ({"frame_every": 0}, ">= 1"),
        ({"init": "spherical"}, "'random' or 'clustered'"),
        ({"init_relax_steps": -1}, ">= 0"),
    ]
    for kwargs, pat in calls:
        with pytest.raises(ValueError) as e:
            brownian_self_assembly(**kwargs)
        assert pat in str(e.value), (kwargs, str(e.value))


def test_cmc_scan_series():
    out = brownian_cmc_scan([0.002, 0.15], n_molecules=60, steps=6000,
                            seed=5, tail_attraction_kT=4.0)
    assert len(out["table"]) == 2
    assert out["method"] and out["caveat"]
    # micelle fraction increases with concentration across the window
    frac = [row["micelle_fraction"] for row in out["table"]]
    assert frac[1] > frac[0]
    # crossing estimate either None (window too narrow) or in-window
    if out["cmc_estimate_M"] is not None:
        assert 0.002 <= out["cmc_estimate_M"] <= 0.15
    with pytest.raises(ValueError, match="at least 2"):
        brownian_cmc_scan([0.1])
    with pytest.raises(ValueError, match="positive"):
        brownian_cmc_scan([0.1, -1.0])
