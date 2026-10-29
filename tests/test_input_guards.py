"""Refusal of physically invalid inputs at the public functions: negative
absolute temperatures, solvent concentrations, ion concentrations and
salinities, and mole fractions outside [0, 1]. A valid call must still
return the textbook value."""

import math

import pytest

from surfactantkit.adsorption import gibbs_gamma_max
from surfactantkit.electrostatics import ionic_strength
from surfactantkit.hld import cc_mixing_rule, hld_nonionic
from surfactantkit.mixed_micelle import activity_coefficients, excess_free_energy
from surfactantkit.thermodynamics import cmc_to_mole_fraction


def test_gibbs_gamma_max_rejects_nonpositive_temperature():
    with pytest.raises(ValueError):
        gibbs_gamma_max(-10.0, "nonionic", -298.15)
    with pytest.raises(ValueError):
        gibbs_gamma_max(-10.0, "nonionic", 0.0)
    assert gibbs_gamma_max(-10.0, "ionic_no_added_salt", 298.15) == pytest.approx(
        0.010 / (2 * 8.314462618 * 298.15))


def test_cmc_to_mole_fraction_rejects_nonpositive_water_molarity():
    with pytest.raises(ValueError):
        cmc_to_mole_fraction(0.0083, -55.5)
    assert cmc_to_mole_fraction(0.0083) == pytest.approx(0.0083 / 55.5083)


def test_ionic_strength_rejects_negative_concentrations():
    with pytest.raises(ValueError):
        ionic_strength({1: -0.1, -1: -0.1})
    assert ionic_strength({1: 0.1, -1: 0.1}) == pytest.approx(0.1)


def test_hld_nonionic_rejects_negative_salinity():
    with pytest.raises(ValueError):
        hld_nonionic(-1.0, 0.13, 0.16, 10.0, 1.0)
    assert hld_nonionic(0.0, 0.13, 0.16, 10.0, 1.0) == pytest.approx(-0.6)


def test_cc_mixing_rule_rejects_fractions_outside_unit_interval():
    with pytest.raises(ValueError):
        cc_mixing_rule([-1.0, 3.0], [-1.0, 2.0])
    assert cc_mixing_rule([-1.0, 3.0], [0.5, 0.5]) == pytest.approx(1.0)


def test_regular_solution_functions_reject_x1_outside_unit_interval():
    with pytest.raises(ValueError):
        activity_coefficients(2.0, -2.0)
    with pytest.raises(ValueError):
        excess_free_energy(1.5, 0.5, 0.5, 298.15)
    with pytest.raises(ValueError):
        excess_free_energy(0.5, 0.5, 0.5, -298.15)
    f1, f2 = activity_coefficients(0.5, -2.0)
    assert f1 == pytest.approx(math.exp(-0.5)) and f2 == pytest.approx(math.exp(-0.5))
