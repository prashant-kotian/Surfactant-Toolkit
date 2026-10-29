"""Each literature check reproduces what the paper got right and reaches its stated finding."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import case1_ctab_vant_hoff as c1  # noqa: E402
import case2_poloxamer_rubingh as c2  # noqa: E402
import case3_oleate_pe65_bound as c3  # noqa: E402
import case4_cpc_tx114_gibbs_helmholtz as c4  # noqa: E402
import case5_dmimcl_des_adsorption as c5  # noqa: E402


def test_case1_reproduces_dG_and_vant_hoff_is_undetermined():
    r = c1.main()
    assert r["reproduce_dG"]["max_abs_diff_kJ"] < 0.1
    assert r["kirchhoff_fit_298"]["dH_kJ"] == r["kirchhoff_fit_298"]["numpy_dH_kJ"]
    assert all(r["findings"].values())


def test_case2_beta_is_negative():
    r = c2.main()
    assert all(r["findings"].values())
    for m in r["mixtures"].values():
        assert m["printed"]["beta"] > 0 > m["one_to_one_by_mass"]["beta"]
        assert abs(m["x_0.5"]["beta"] - m["x_0.5"]["numpy_beta"]) < 0.02


def test_case3_mixed_cmc_violates_the_pseudo_phase_bound():
    assert all(c3.main()["findings"].values())


def test_case4_dG_and_dH_cannot_both_hold():
    f = c4.main()["findings"]
    assert f["gh_dH_positive_303_323"] and f["printed_dH_negative_all_T"]
    assert f["mixed_vs_cpc_293_recomputed_more_spontaneous"] and not f["mixed_vs_cpc_293_printed_more_spontaneous"]


def test_case5_recovered_packing_is_spherical():
    r = c5.main()
    assert r["findings"]["printed_A_below_chain_cross_section"]
    assert r["findings"]["recovered_all_spherical_le_one_third"]
    for v in r["systems"].values():
        assert abs(v["A_from_printed_Gamma_nm2"] - v["printed"]["A_min_nm2"]) < 0.005
        assert abs(v["recovered"]["check_pi_at_C20"] - 20.0) < 0.01
