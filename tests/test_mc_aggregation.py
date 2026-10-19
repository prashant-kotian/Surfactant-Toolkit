"""Tests for surfactantkit.thermodynamics.mc_aggregation_distribution (S1).

Gold: exact Boltzmann histogram over integer partitions (labeled cluster
model).  MC marginals must match it within sampling error; chain must be
seed-reproducible and flag its own convergence against the statistical
noise floor.
"""
import math

import numpy as np
import pytest

from surfactantkit.thermodynamics import (
    _exact_aggregation_marginals,
    _integer_partitions_of_n,
    mc_aggregation_distribution,
)


def test_integer_partition_generator_counts():
    # p(1)=1, p(5)=7, p(10)=42
    assert sum(1 for _ in _integer_partitions_of_n(1)) == 1
    assert sum(1 for _ in _integer_partitions_of_n(5)) == 7
    assert sum(1 for _ in _integer_partitions_of_n(10)) == 42


def test_exact_marginals_normalized_and_balanced():
    ex = np.array(_exact_aggregation_marginals(8, 1.6))
    assert ex.sum() == pytest.approx(1.0, abs=1e-9)
    # independent check: fraction of free monomers = Z_(N-1)/Z_N via
    # recurrence Z_n = sum_k C(n-1,k-1) exp(beta(k-1)) Z_{n-k}
    beta = 1.6
    Z = [1.0, 1.0]
    for n in range(2, 9):
        Z.append(sum(
            math.comb(n - 1, k - 1) * np.exp(beta * (k - 1)) * Z[n - k]
            for k in range(1, n + 1)))
    assert ex[1] == pytest.approx(Z[7] / Z[8], rel=1e-9)


def test_mc_matches_exact_default_run():
    r = mc_aggregation_distribution(seed=4)
    assert r.exact_computed
    assert r.exact_max_abs_deviation is not None
    assert r.exact_max_abs_deviation < 0.02
    assert r.converged
    assert r.free_monomer_fraction > 0


def test_mc_deterministic_for_seed():
    a = mc_aggregation_distribution(n_monomers=12, epsilon_kJ_per_mol=6.0,
                                    n_steps=60_000, seed=9)
    b = mc_aggregation_distribution(n_monomers=12, epsilon_kJ_per_mol=6.0,
                                    n_steps=60_000, seed=9)
    assert a.histogram_n_monomer_fraction == b.histogram_n_monomer_fraction
    c = mc_aggregation_distribution(n_monomers=12, epsilon_kJ_per_mol=6.0,
                                    n_steps=60_000, seed=10)
    assert c.histogram_n_monomer_fraction != a.histogram_n_monomer_fraction


def test_stronger_binding_shifts_distribution_up():
    weak = mc_aggregation_distribution(n_monomers=25, epsilon_kJ_per_mol=4.0,
                                       n_steps=250_000, seed=3)
    strong = mc_aggregation_distribution(n_monomers=25, epsilon_kJ_per_mol=12.0,
                                         n_steps=250_000, seed=3)
    assert strong.mean_cluster_size_monomer_weighted > \
        weak.mean_cluster_size_monomer_weighted * 1.5
    assert strong.free_monomer_fraction < weak.free_monomer_fraction


def test_histogram_sums_to_one():
    r = mc_aggregation_distribution(n_monomers=20, n_steps=100_000, seed=2)
    total = sum(f for _, f in r.histogram_n_monomer_fraction)
    assert total == pytest.approx(1.0, abs=1e-9)
    assert all(n_ >= 1 for n_, _ in r.histogram_n_monomer_fraction)
    assert r.largest_cluster_observed <= 20


def test_validation_errors():
    with pytest.raises(ValueError, match="n_monomers"):
        mc_aggregation_distribution(n_monomers=1)
    with pytest.raises(ValueError, match="epsilon"):
        mc_aggregation_distribution(epsilon_kJ_per_mol=0.0)
    with pytest.raises(ValueError, match="n_steps"):
        mc_aggregation_distribution(n_steps=10)
    with pytest.raises(ValueError, match="temperature"):
        mc_aggregation_distribution(temperature_K=-5.0)
    with pytest.raises(ValueError, match="burn_in"):
        mc_aggregation_distribution(burn_in_fraction=0.7)
    with pytest.raises(ValueError, match="thin"):
        mc_aggregation_distribution(thin=0)


def test_result_metadata():
    r = mc_aggregation_distribution(n_monomers=10, n_steps=50_000, seed=1)
    assert r.method and r.params["seed"] == 1
    assert "bond-energy cluster model" in r.caveat
    assert r.n_samples > 500
    assert r.split_half_noise_floor > 0
