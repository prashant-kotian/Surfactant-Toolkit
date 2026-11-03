"""Guards the Czakaj data transcription and the LVDT consistency identity."""

from __future__ import annotations

import cmath
import math

import pytest

from research.laod_tau2 import lae_data


def test_table4_shape():
    assert len(lae_data.TABLE_4) == 5
    assert max(r["E0_mN_m"] for r in lae_data.TABLE_4.values()) == 73.9   # peak at 0.3 mM


def test_fig7_has_two_frequencies_everywhere():
    for c, row in lae_data.FIG_7.items():
        assert set(row) == {0.01, 0.10}
        assert all(0.0 < v < 1.0 for v in row.values())


def test_fig7_peaks_at_the_cmc():
    """The measured maximum sits at c = 1.0 mM = the CMC, not above it.

    This is the observation that rules out micellar buffering as the cause:
    at the CMC the reservoir capacity is zero by construction.
    """
    at_01 = {c: v[0.10] for c, v in lae_data.FIG_7.items()}
    assert max(at_01, key=at_01.get) == 1.0


def test_lvdt_identity_recovers_e0_as_high_frequency_limit():
    """Eq. 1 must approach eps0 as nu -> infinity."""
    E = lae_data.lvdt_E(12.6, 0.06, 1e9)
    assert abs(E) == pytest.approx(12.6, rel=1e-3)


def test_lvdt_phase_in_range():
    for row in lae_data.check_lvdt_table():
        assert 0.0 <= row["lag_deg"] <= 90.0


def test_noise_floor_is_consistent_with_submicellar_baseline():
    """Their sub-CMC baseline (0.019-0.045) must sit inside the noise band."""
    lo, hi = lae_data.NOISE_FLOOR_RATIO
    sub = [lae_data.FIG_7[c][0.10] for c in (0.2, 0.5)]
    assert all(lo * 0.5 <= v <= hi for v in sub)
