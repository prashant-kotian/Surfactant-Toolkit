"""S4 analysis-core tests (Anderson short-time CMC, Rummel et al. JPCB 2026).

Golds:
- eq 3 exactness (beta=0 -> CMC = C_free bit-identical, phi-independent);
- eq 4 forward/inverse round-trip + independent hand evaluation of the
  published equation;
- eq-8 profile fit recovers A/B/tau of an exactly generated trace to
  machine precision;
- jackknife sigma reproduces the closed-form leave-one-out value on
  traces whose subset means are analytic (exact exponential family);
- beta=0 gives sigma_CMC = c_tot * sigma_A exactly;
- all input guards name their fields; flat traces -> converged=False.
"""
import math

import pytest

from surfactantkit.dpd_cmc import (
    AndersonCMCResult,
    anderson_cmc_analysis,
    cmc_from_free_monomer,
    free_monomer_from_cmc,
)

TIMES = [i * 5.0 for i in range(60)]
A_TRUE, B_TRUE, TAU_TRUE = 0.02, 0.5, 120.0


def _exp_trace(a: float, b: float = B_TRUE, tau: float = TAU_TRUE):
    return [a + (b - a) * math.exp(-t / tau) for t in TIMES]


# ---------------- eq 3 / eq 4 ----------------

def test_eq3_beta0_is_exact_c_free():
    assert cmc_from_free_monomer(0.0025, 0.05, beta=0.0) == 0.0025
    # phi cannot matter at beta=0 (eq 4 reduces to eq 3)
    assert (cmc_from_free_monomer(0.0025, 0.05, beta=0.0, v_m_cm3_per_mol=400.0)
            == 0.0025)


def test_eq3_inverse_beta0():
    x = free_monomer_from_cmc(0.004, 0.08, beta=0.0)
    assert abs(x - 0.004) < 1e-12


def test_eq4_forward_matches_hand_equation():
    c_free, c_tot, beta, v_m = 0.003, 0.06, 0.6, 350.0
    phi = v_m * c_tot / 1000.0
    inner = beta * c_free + (1.0 - beta) * c_tot / (1.0 - phi)
    expected = math.exp(
        (math.log(c_free) + beta * math.log(inner)) / (1.0 + beta))
    got = cmc_from_free_monomer(c_free, c_tot, beta=beta, v_m_cm3_per_mol=v_m)
    assert abs(got - expected) <= 1e-15 * expected


def test_eq4_round_trip_forward_inverse():
    c_free, c_tot, beta, v_m = 0.003, 0.06, 0.6, 350.0
    cmc = cmc_from_free_monomer(c_free, c_tot, beta=beta, v_m_cm3_per_mol=v_m)
    back = free_monomer_from_cmc(cmc, c_tot, beta=beta, v_m_cm3_per_mol=v_m)
    assert abs(back - c_free) / c_free < 1e-9


def test_eq4_crowding_increases_cmc_at_fixed_beta():
    beta, v_m = 0.6, 500.0
    with_phi = cmc_from_free_monomer(0.002, 0.1, beta=beta, v_m_cm3_per_mol=v_m)
    without = cmc_from_free_monomer(0.002, 0.1, beta=beta, v_m_cm3_per_mol=None)
    assert with_phi > without


def test_eq4_guard_messages_name_fields():
    with pytest.raises(ValueError, match="c_free"):
        cmc_from_free_monomer(0.0, 0.05)
    with pytest.raises(ValueError, match="c_tot"):
        cmc_from_free_monomer(0.002, -1.0)
    with pytest.raises(ValueError, match="beta"):
        cmc_from_free_monomer(0.002, 0.05, beta=1.2)
    with pytest.raises(ValueError, match="v_m_cm3_per_mol"):
        cmc_from_free_monomer(0.002, 0.05, beta=0.5, v_m_cm3_per_mol=9000.0)
    with pytest.raises(ValueError, match="phi"):
        cmc_from_free_monomer(0.002, 100.0, beta=0.5, v_m_cm3_per_mol=4000.0)
    with pytest.raises(ValueError, match="cmc"):
        free_monomer_from_cmc(0.5, 0.1)
    with pytest.raises(ValueError, match="cmc"):
        free_monomer_from_cmc(-1.0, 0.1)


# ---------------- eq-8 fit + jackknife analysis ----------------

def test_analysis_recovers_exact_trace_parameters():
    traces = [_exp_trace(A_TRUE) for _ in range(10)]
    r = anderson_cmc_analysis(TIMES, traces, c_tot=0.1, beta=0.0)
    assert isinstance(r, AndersonCMCResult)
    assert r.converged
    assert abs(r.p_free_inf - A_TRUE) < 1e-6
    assert abs(r.p_free_0 - B_TRUE) < 1e-6
    assert abs(r.tau_s - TAU_TRUE) < 1e-3
    assert r.n_subsets == 10          # C(10, 9) leave-one-out, j=0.9
    assert len(r.subset_cmc_M) == 10
    assert r.params["tau_edge_pinned"] is False
    assert r.n_iterations == 11       # full set + 10 subsets


def test_beta0_cmc_is_a_times_c_tot_and_sigma_scales_exactly():
    traces = [_exp_trace(A_TRUE) for _ in range(10)]
    r = anderson_cmc_analysis(TIMES, traces, c_tot=0.1, beta=0.0)
    assert abs(r.cmc_M - r.p_free_inf * 0.1) <= 1e-15 * r.cmc_M
    # beta=0: CMC = A * c_tot linearly -> jackknife sigma transforms exactly
    assert r.cmc_jackknife_sigma_M == 0.0          # identical traces
    assert r.p_free_inf_jackknife_sigma == 0.0


def test_jackknife_sigma_matches_closed_form_leave_one_out():
    # traces from the exact exponential family with distinct A_i:
    # every subset mean trace is again exact, so subset A = mean of the
    # member A_i and the jackknife has a hand-computable value.
    a_vals = [A_TRUE * (1.0 + 0.05 * (i - 4.5) / 4.5) for i in range(10)]
    traces = [_exp_trace(a) for a in a_vals]
    r = anderson_cmc_analysis(TIMES, traces, c_tot=0.1, beta=0.0)
    subset_a = [(sum(a_vals) - a_vals[i]) / 9 for i in range(10)]
    mean = sum(subset_a) / 10
    expected = math.sqrt(9 / 10 * sum((x - mean) ** 2 for x in subset_a))
    assert abs(r.p_free_inf_jackknife_sigma - expected) < 1e-6
    assert r.cmc_jackknife_sigma_M > 0.0
    # beta=0 exact linear propagation of sigma through CMC = A * c_tot
    assert abs(r.cmc_jackknife_sigma_M - 0.1 * r.p_free_inf_jackknife_sigma) < 1e-12


def test_determinism_same_input_same_output():
    traces = [_exp_trace(A_TRUE * (1 + 0.01 * i)) for i in range(10)]
    a = anderson_cmc_analysis(TIMES, traces, c_tot=0.05, beta=0.3,
                              v_m_cm3_per_mol=350.0)
    b = anderson_cmc_analysis(TIMES, traces, c_tot=0.05, beta=0.3,
                              v_m_cm3_per_mol=350.0)
    assert a.as_dict() == b.as_dict()
    assert a.converged and b.converged


def test_schema_keys_and_ionic_path():
    traces = [_exp_trace(A_TRUE) for _ in range(10)]
    r = anderson_cmc_analysis(TIMES, traces, c_tot=0.1, beta=0.6,
                              v_m_cm3_per_mol=350.0, j_ratio=0.9)
    d = r.as_dict()
    for k in ("converged", "n_iterations", "method", "params", "caveat",
              "cmc_M", "cmc_jackknife_sigma_M", "p_free_inf", "tau_s"):
        assert k in d
    assert d["method"] == "anderson_pfree_exp_profile_fit_jackknife_eq4"
    assert "j_ratio" in d["caveat"] and "doi:10.1021/acs.jpcb.5c07475" in d["caveat"]
    # ionic eq 4 gives a CMC above the naive plateau reading
    assert r.cmc_M > r.p_free_inf * 0.1
    assert r.cmc_jackknife_sigma_M == 0.0   # identical traces


def test_flat_traces_do_not_converge():
    traces = [[0.5] * len(TIMES) for _ in range(10)]
    r = anderson_cmc_analysis(TIMES, traces, c_tot=0.1, beta=0.0)
    assert not r.converged
    assert r.cmc_M is None
    assert r.cmc_jackknife_sigma_M is None


def test_analysis_guards():
    ok_traces = [_exp_trace(A_TRUE) for _ in range(10)]
    with pytest.raises(ValueError, match="traces"):
        anderson_cmc_analysis(TIMES, ok_traces[:2], c_tot=0.1)
    with pytest.raises(ValueError, match="times"):
        anderson_cmc_analysis([1.0, 1.0, 2.0, 3.0], ok_traces, c_tot=0.1)
    with pytest.raises(ValueError, match=r"traces\[0\]"):
        anderson_cmc_analysis(TIMES, [t[:10] for t in ok_traces], c_tot=0.1)
    bad = [_exp_trace(A_TRUE) for _ in range(10)]
    bad[0][5] = 2.5
    with pytest.raises(ValueError, match="RATIO"):
        anderson_cmc_analysis(TIMES, bad, c_tot=0.1)
    with pytest.raises(ValueError, match="c_tot"):
        anderson_cmc_analysis(TIMES, ok_traces, c_tot=0.0)
    with pytest.raises(ValueError, match="beta"):
        anderson_cmc_analysis(TIMES, ok_traces, c_tot=0.1, beta=-0.1)
    with pytest.raises(ValueError, match="j_ratio"):
        anderson_cmc_analysis(TIMES, ok_traces, c_tot=0.1, j_ratio=0.3)
    with pytest.raises(ValueError, match="max_subsets"):
        anderson_cmc_analysis(TIMES, ok_traces, c_tot=0.1, j_ratio=0.9,
                              max_subsets=3)
    with pytest.raises(ValueError, match="j_ratio"):
        # n=3, j=0.5 -> 1 trace per subset: below the >=2 floor
        anderson_cmc_analysis(TIMES, ok_traces[:3], c_tot=0.1, j_ratio=0.5)
