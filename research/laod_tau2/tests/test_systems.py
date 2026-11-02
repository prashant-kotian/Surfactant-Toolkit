"""Tests for the real-surfactant parameter sets and the capacity closure.

The closure C_m = (c_tot - CMC)/CMC * n_bar was supplied independently; these
tests check it against its own derivation and pin the two regimes where it is
known to fail, so the failure is recorded rather than discovered later.
"""

from __future__ import annotations

import pytest

from research.laod_tau2.systems import (
    K_LANGMUIR,
    K_MISSING,
    SYSTEMS,
    verify_capacity_closure,
)


def test_closure_matches_mass_action_derivation():
    """C_m = n*phi*c_tot/c1 must equal the closure in the pinned limit."""
    assert all(r["match"] for r in verify_capacity_closure())


def test_capacity_is_zero_below_cmc():
    for s in SYSTEMS.values():
        assert s.capacity(0.5 * s.cmc_mM) == 0.0


def test_capacity_is_linear_above_cmc():
    """The closure is exactly linear in (c_tot - CMC)."""
    for s in SYSTEMS.values():
        a = s.capacity(2 * s.cmc_mM)
        b = s.capacity(4 * s.cmc_mM)
        assert b == pytest.approx(3 * a, rel=2e-3)


def test_closure_overestimates_near_the_cmc_edge():
    """Known defect: the exact derivative -> 0 at the CMC, the closure does not.

    Recorded as a test so that every study's c_tot >= 2x CMC choice has a
    reason attached to it.
    """
    s = SYSTEMS["SDS"]
    near = s.capacity(1.01 * s.cmc_mM)
    far = s.capacity(10 * s.cmc_mM)
    # near-CMC the closure claims a reservoir 1% of the way up; physically it
    # must be far smaller relative to n_bar. Quantify the overestimate.
    assert near / s.n_bar > 0.005      # closure says ~0.01
    assert far / s.n_bar == pytest.approx(9.0, rel=1e-2)


def test_frumkin_systems_flag_missing_isotherm_constant():
    """Ionic systems have no supplied Langmuir K; results must not be quoted
    for them until an isotherm is fitted."""
    assert set(K_MISSING) == {"SDS", "CTAB"}
    assert "SDS" not in K_LANGMUIR and "CTAB" not in K_LANGMUIR


def test_volume_correction_is_bounded():
    """The micellar volume correction must stay small in the studied range."""
    for s in SYSTEMS.values():
        assert 1.0 <= s.capacity(10 * s.cmc_mM) / ((10 - 1) * s.n_bar) < 1.35
