"""Tests for surfactantkit.nagarajan (worker C1) -- a-priori CMC
prediction from the Nagarajan molecular-thermodynamic free-energy model
(Langmuir 1985, 1, 331-341, eqs 5-7 + appendix geometry + Table I).

Golds:
  (a) seven experimental CMCs (the source's own validation class) -- all
      within 0.15 log10 units (worst observed at implementation time 0.09);
  (b) Klevens slope of the sulfate series;
  (c) term-table identities recomputed independently in the test;
  (d) sphere/cylinder geometry continuity at the sphere limit;
  (e) Debye kappa vs the source's closed-form 25 C expression;
  (f) salt lowers the CMC, self-consistent kappa iteration converges,
      determinism, and validation errors.
"""
import math

import pytest

from surfactantkit.nagarajan import (
    HEADGROUPS,
    NagarajanCMCResult,
    _kappa_per_m,
    nagarajan_cmc_predictor,
    nagarajan_free_energy_terms,
)

# (label, n_carbons, headgroup, experimental CMC in mol/L)
SEVEN = [
    ("sodium octyl sulfate", 8, "alkyl_sulfate", 0.13),
    ("sodium decyl sulfate", 10, "alkyl_sulfate", 0.040),
    ("sodium dodecyl sulfate", 12, "alkyl_sulfate", 0.0082),
    ("sodium tetradecyl sulfate", 14, "alkyl_sulfate", 0.0021),
    ("C10PO", 10, "dialkyl_phosphine_oxide", 0.0037),
    ("C10SO", 10, "alkyl_methyl_sulfoxide", 0.0016),
    ("C12E8", 12, "polyoxyethylene", 0.00010),
]

N_A = 6.02214076e23
K_B = 1.380649e-23
E_CHARGE = 1.602176634e-19
EPS0 = 8.8541878128e-12


def test_all_seven_experimental_cmcs_within_015_log_units():
    worst = 0.0
    for label, nc, hg, exp in SEVEN:
        r = nagarajan_cmc_predictor(nc, hg)
        err = abs(math.log10(r.cmc_mol_per_L / exp))
        worst = max(worst, err)
        assert err <= 0.15, f"{label}: log10 pred/exp = {err:.3f}"
    assert worst <= 0.15


def test_klevens_slope_of_sulfate_series():
    logs = [math.log10(nagarajan_cmc_predictor(nc, "alkyl_sulfate").cmc_mol_per_L)
            for nc in (8, 10, 12, 14)]
    slope = (logs[-1] - logs[0]) / 6.0
    # empirical Klevens slope ~ -0.28..-0.32 log10 CMC per CH2
    assert -0.40 <= slope <= -0.20, slope


def test_kappa_iteration_converges_for_all_seven():
    for label, nc, hg, _ in SEVEN:
        r = nagarajan_cmc_predictor(nc, hg)
        assert r.converged, f"{label}: not converged in {r.n_iterations}"
        assert r.n_iterations <= 20
        assert r.monomer_conc_iterated_M > 0
        assert r.kappa_nm_inv > 0


def test_results_deterministic():
    a = nagarajan_cmc_predictor(12, "alkyl_sulfate")
    b = nagarajan_cmc_predictor(12, "alkyl_sulfate")
    assert a.as_dict() == b.as_dict()


def test_term_table_recomputed_independently():
    """Literal re-derivation of eq 5 sphere terms at one point (n=12, g=40)."""
    nc, g = 12, 40
    t = nagarajan_free_energy_terms(
        g, nc, "alkyl_sulfate", concentration_mol_per_L=1e-3)
    kT = K_B * 298.15
    v0 = (27 * (nc - 1) + 54) * 1e-30
    l0 = (1.50 + 1.26 * nc) * 1e-10
    r0 = (3 * g * v0 / (4 * math.pi)) ** (1 / 3)
    assert r0 <= l0
    a = 4 * math.pi * r0 * r0 / g

    transfer = -g * (2.05 + 1.49 * nc)
    ordering = g * (-0.5 + 0.24 * nc)
    interface = g * 0.05 * (a - 17e-20) / kT          # a0 = min(17, 21) A2
    steric = -g * math.log(1 - 17e-20 / a)
    # source's 25 C kappa closed form (SI conversion): sqrt(C)/3.08e-8 cm^-1
    kappa = (math.sqrt(1e-3) / 3.08e-8) * 100.0
    q2_over_4pe0 = E_CHARGE ** 2 / (4 * math.pi * EPS0 * 80.0)
    brac = (1 + kappa * 1e-10) / (1 + kappa * 1e-10 + kappa * (r0 + 2.3e-10))
    elec = g * g * 0.46 * q2_over_4pe0 / (2 * kT * (r0 + 2.3e-10)) * brac

    assert t.transfer == pytest.approx(transfer, rel=1e-12)
    assert t.ordering == pytest.approx(ordering, rel=1e-12)
    assert t.interface == pytest.approx(interface, rel=1e-10)
    assert t.steric == pytest.approx(steric, rel=1e-12)
    # module uses general-SI kappa (~1% from the paper form): allow 3%
    assert t.electrostatic == pytest.approx(elec, rel=0.03)
    assert t.total == pytest.approx(
        t.transfer + t.ordering + t.interface + t.steric + t.electrostatic,
        rel=1e-12)
    assert t.per_molecule == pytest.approx(t.total / g, rel=1e-12)


def test_nonionic_has_zero_electrostatics():
    t = nagarajan_free_energy_terms(
        40, 12, "polyoxyethylene", concentration_mol_per_L=1e-3)
    assert t.electrostatic == 0.0
    t_ion = nagarajan_free_energy_terms(
        40, 12, "alkyl_sulfate", concentration_mol_per_L=1e-3)
    assert t_ion.electrostatic > 0


def test_sphere_cylinder_geometry_continuous_at_sphere_limit():
    nc = 12
    v0 = (27 * (nc - 1) + 54) * 1e-30
    l0 = (1.50 + 1.26 * nc) * 1e-10
    g_s = int(4 * math.pi * l0 ** 3 / (3 * v0))
    # sphere just inside the limit
    sp = nagarajan_free_energy_terms(g_s - 1, nc, "alkyl_sulfate",
                                     branch="sphere",
                                     concentration_mol_per_L=1e-3)
    # cylinder just outside (L tiny -> r0 = l0, a ~ 4 pi l0^2 / g)
    cy = nagarajan_free_energy_terms(g_s + 1, nc, "alkyl_sulfate",
                                     branch="cylinder",
                                     concentration_mol_per_L=1e-3)
    assert sp.core_radius_nm == pytest.approx(l0 * 1e9, rel=0.02)
    assert cy.core_radius_nm == pytest.approx(sp.core_radius_nm, rel=0.05)
    # areas match to the small g-difference
    assert cy.area_per_molecule_A2 == pytest.approx(
        sp.area_per_molecule_A2, rel=0.05)


def test_kappa_matches_source_closed_form_at_25C():
    for conc in (1e-4, 1e-3, 1e-2, 0.1):
        paper = (math.sqrt(conc) / 3.08e-8) * 100.0  # cm^-1 -> m^-1
        assert _kappa_per_m(conc, 298.15) == pytest.approx(paper, rel=0.02)


def test_added_salt_lowers_cmc():
    plain = nagarajan_cmc_predictor(12, "alkyl_sulfate")
    salty = nagarajan_cmc_predictor(12, "alkyl_sulfate", added_salt_M=0.1)
    assert salty.cmc_mol_per_L < 0.5 * plain.cmc_mol_per_L
    assert salty.kappa_nm_inv > plain.kappa_nm_inv


def test_delta_g_in_physical_range():
    for label, nc, hg, _ in SEVEN:
        r = nagarajan_cmc_predictor(nc, hg)
        assert -25.0 < r.delta_g_mic_kT_per_molecule < -4.0, label


def test_result_metadata_complete():
    r = nagarajan_cmc_predictor(10, "alkyl_methyl_sulfoxide")
    assert isinstance(r, NagarajanCMCResult)
    assert r.method.startswith("Nagarajan")
    assert r.caveat
    assert r.params["a_p_A2"] == 37.0
    assert r.branch_optimum in ("sphere", "cylinder")
    assert r.terms_at_optimum["method"].startswith("Nagarajan")


def test_custom_headgroup_requires_all_fields():
    with pytest.raises(ValueError, match="a_p, a_i, delta, ionic"):
        nagarajan_cmc_predictor(12, "custom", a_p=20.0)
    r = nagarajan_cmc_predictor(12, "custom", a_p=17.0, a_i=1.0,
                                delta=2.3, ionic=True)
    # custom sulfates should behave like the preset
    ref = nagarajan_cmc_predictor(12, "alkyl_sulfate")
    assert r.cmc_mol_per_L == pytest.approx(ref.cmc_mol_per_L, rel=1e-9)


def test_validation_errors():
    with pytest.raises(ValueError, match="headgroup"):
        nagarajan_cmc_predictor(12, "nonexistent_headgroup")
    with pytest.raises(ValueError, match="n_carbons"):
        nagarajan_cmc_predictor(3, "alkyl_sulfate")
    with pytest.raises(ValueError, match="added_salt"):
        nagarajan_cmc_predictor(12, "alkyl_sulfate", added_salt_M=-0.1)
    with pytest.raises(ValueError, match="g must be >= 2"):
        nagarajan_free_energy_terms(1, 12, "alkyl_sulfate")
    with pytest.raises(ValueError, match="branch"):
        nagarajan_free_energy_terms(40, 12, "alkyl_sulfate", branch="torus")


def test_table_i_headgroup_presets():
    assert HEADGROUPS["alkyl_sulfate"] == (17.0, 1.0, 2.3, True)
    assert HEADGROUPS["potassium_carboxylate"] == (11.0, 1.35, 1.0, True)
    assert HEADGROUPS["polyoxyethylene"] == (36.0, 0.0, 0.0, False)
    # nonionics are non-ionic by construction
    for hg, (_, _, _, ionic) in HEADGROUPS.items():
        if hg in ("dialkyl_phosphine_oxide", "alkyl_methyl_sulfoxide",
                  "polyoxyethylene"):
            assert ionic is False
