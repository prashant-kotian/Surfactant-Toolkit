"""Real numeric validation for electrostatics.charge_regulation_pb (E2).

Like kinetics.py, this worker had no dedicated test anywhere in the repo
before this file -- only a generic MCP tool-registration smoke test --
despite ROADMAP.md claiming "ψ0 band, θ range, pH monotone [done]". This
file checks real physical identities: the weak-coupling (vanishing site
density) limit must recover the bare Henderson-Hasselbalch fraction, the
strong-coupling regime must show genuine self-limiting (theta below the
bare HH value), the site/double-layer charges must self-consistently
balance, and the sign/monotonicity behavior must match the stated acid/
base convention.
"""

import pytest

from surfactantkit.electrostatics import charge_regulation_pb, debye_length


def bare_hh_ionized_fraction(pka: float, ph: float) -> float:
    return 1.0 / (1.0 + 10.0 ** (pka - ph))


def test_weak_coupling_limit_recovers_bare_henderson_hasselbalch():
    """At vanishing site density there's nothing to self-limit: the
    regulated surface potential/pH must collapse to the bulk values and
    theta must match the bare (unregulated) Henderson-Hasselbalch
    fraction."""
    pka, ph = 5.0, 7.0
    r = charge_regulation_pb("acid", pka=pka, site_density_mol_per_m2=1e-12,
                              ph=ph, ionic_strength_M=0.1)
    assert abs(r.psi0_mV) < 1e-2
    assert r.surface_pH == pytest.approx(ph, abs=1e-3)
    assert r.theta_charged == pytest.approx(bare_hh_ionized_fraction(pka, ph), rel=1e-4)


def test_strong_coupling_shows_real_self_limiting():
    """At realistic saturated-monolayer site density and a large pKa-pH
    mismatch, charge regulation must suppress ionization well below the
    bare (fixed-sigma) Henderson-Hasselbalch prediction -- this
    suppression is the entire physical point of the model."""
    pka, ph = 5.0, 9.0
    bare = bare_hh_ionized_fraction(pka, ph)
    r = charge_regulation_pb("acid", pka=pka, site_density_mol_per_m2=3e-6,
                              ph=ph, ionic_strength_M=0.01)
    assert r.converged is True
    assert bare > 0.999  # bare model says "essentially fully ionized"
    assert r.theta_charged < 0.95  # regulation pulls it down substantially
    assert r.theta_charged < bare


def test_site_and_double_layer_charge_self_consistent():
    r = charge_regulation_pb("acid", pka=5.0, site_density_mol_per_m2=3e-6,
                              ph=9.0, ionic_strength_M=0.01)
    assert r.sigma_sites_C_per_m2 == pytest.approx(r.sigma_double_layer_C_per_m2, rel=1e-8)
    assert r.params["relative_sigma_imbalance"] < 1e-8


def test_acid_site_gives_nonpositive_potential_base_site_nonnegative():
    acid = charge_regulation_pb("acid", pka=5.0, site_density_mol_per_m2=3e-6,
                                  ph=9.0, ionic_strength_M=0.01)
    base = charge_regulation_pb("base", pka=9.0, site_density_mol_per_m2=3e-6,
                                  ph=5.0, ionic_strength_M=0.01)
    assert acid.psi0_mV <= 1e-9
    assert base.psi0_mV >= -1e-9
    # acid's (pKa-pH) and base's (pH-pKa) are the SAME exponent (5-9 ==
    # 5-9) here, so these two cases are a true mirror pair -> same |psi0|.
    assert abs(acid.psi0_mV) == pytest.approx(abs(base.psi0_mV), rel=1e-9)


def test_theta_monotone_nondecreasing_with_ph_for_acid_site():
    pka = 7.0
    phs = [3.0, 5.0, 7.0, 9.0, 11.0]
    thetas = [
        charge_regulation_pb("acid", pka=pka, site_density_mol_per_m2=2e-6,
                               ph=p, ionic_strength_M=0.05).theta_charged
        for p in phs
    ]
    assert all(thetas[i] <= thetas[i + 1] + 1e-12 for i in range(len(thetas) - 1))


def test_debye_length_matches_direct_formula():
    r = charge_regulation_pb("acid", pka=5.0, site_density_mol_per_m2=3e-6,
                               ph=9.0, ionic_strength_M=0.01, temperature_K=298.15)
    assert r.debye_length_nm == pytest.approx(debye_length(0.01, 298.15), rel=1e-12)


@pytest.mark.parametrize("bad_kwargs", [
    {"site_type": "amphoteric", "pka": 5.0, "site_density_mol_per_m2": 3e-6, "ph": 7.0, "ionic_strength_M": 0.1},
    {"site_type": "acid", "pka": float("nan"), "site_density_mol_per_m2": 3e-6, "ph": 7.0, "ionic_strength_M": 0.1},
    {"site_type": "acid", "pka": 5.0, "site_density_mol_per_m2": 0.0, "ph": 7.0, "ionic_strength_M": 0.1},
    {"site_type": "acid", "pka": 5.0, "site_density_mol_per_m2": 3e-6, "ph": float("inf"), "ionic_strength_M": 0.1},
    {"site_type": "acid", "pka": 5.0, "site_density_mol_per_m2": 3e-6, "ph": 7.0, "ionic_strength_M": 0.0},
    {"site_type": "acid", "pka": 5.0, "site_density_mol_per_m2": 3e-6, "ph": 7.0, "ionic_strength_M": 0.1, "temperature_K": 0.0},
    {"site_type": "acid", "pka": 5.0, "site_density_mol_per_m2": 3e-6, "ph": 7.0, "ionic_strength_M": 0.1, "rel_permittivity": 0.0},
    {"site_type": "acid", "pka": 5.0, "site_density_mol_per_m2": 3e-6, "ph": 7.0, "ionic_strength_M": 0.1, "u_bounds": (5.0, -5.0)},
])
def test_error_paths(bad_kwargs):
    with pytest.raises(ValueError):
        charge_regulation_pb(**bad_kwargs)
