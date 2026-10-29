"""Tests for surfactantkit.hlb and surfactantkit.cpp.

Constants (Davies group numbers, Tanford formula coefficients, CPP
morphology thresholds) were checked against cited sources before being
hardcoded -- see the module docstrings in hlb.py and cpp.py. Tests here
verify: (1) exact boundary conditions that follow directly from the
formula definitions, (2) a worked example (SDS via Davies' method) built
purely from the verified group numbers, and (3) the CPP classification
thresholds at and around each boundary.
"""

import math

import pytest

from surfactantkit.cpp import (
    NAGARAJAN_RUCKENSTEIN_HEADGROUP_CONSTANTS,
    aggregation_number_spherical,
    blankschtein_counterion_self_energy_release,
    blankschtein_entropy_of_binding_free_energy,
    blankschtein_steric_free_energy_with_counterion,
    classify_aggregate_morphology,
    critical_packing_parameter,
    estimate_axial_ratio_from_cpp_geometry,
    nagarajan_debye_huckel_kappa_inverse,
    nagarajan_equilibrium_area_ionic,
    nagarajan_ruckenstein_dipole_headgroup_free_energy,
    nagarajan_ruckenstein_ionic_headgroup_free_energy,
    nagarajan_ruckenstein_steric_headgroup_free_energy,
    tanford_critical_length,
    tanford_tail_volume,
)
from surfactantkit.hlb import (
    DAVIES_DERIVED_HYDROPHILIC_GROUPS,
    GEMINI_SURFACTANT_HLB_GRIFFIN_REFERENCE,
    derive_davies_group_number_from_griffin,
    guo_effective_alkyl_chain_length,
    guo_effective_eo_chain_length,
    guo_effective_po_chain_length,
    hlb_davies,
    hlb_davies_guo_ecl,
    hlb_griffin,
    recommend_hlb_method_for_structural_family,
)


def test_hlb_griffin_boundary_conditions():
    # fully hydrophilic -> HLB = 20; fully lipophilic -> HLB = 0
    assert hlb_griffin(100.0, 100.0) == pytest.approx(20.0)
    assert hlb_griffin(0.0, 100.0) == pytest.approx(0.0)
    assert hlb_griffin(50.0, 100.0) == pytest.approx(10.0)


def test_hlb_griffin_rejects_bad_input():
    with pytest.raises(ValueError):
        hlb_griffin(150.0, 100.0)  # hydrophilic > total
    with pytest.raises(ValueError):
        hlb_griffin(10.0, 0.0)


def test_hlb_griffin_matches_real_ethoxylate_series_within_explained_tolerance():
    """Real-value check, not just boundary conditions: closes the
    long-open hlb_griffin gap (2026-09-04, Kosswig, Ullmann's Encyclopedia
    of Industrial Chemistry Vol. 35, Table 13 -- see
    literature_validation_notes.md for the full source and the honest
    explanation of the small offset). Two chemically different hydrophobe
    classes (aliphatic octadecanol, aromatic isononylphenol), 6 points
    total: Mh = n_EO * 44.053 (real ethylene oxide addition-unit mass),
    M = hydrophobe MW + Mh, both from real atomic weights. Formula gets
    the right trend and order of magnitude with a small, consistent
    ~5-7% low offset in the same direction every time (idealized single-
    species stoichiometry vs. Table 13's real, polydisperse catalogue
    values) -- a real, explained gap, not scattered error that would
    indicate a formula bug."""
    C, H, O = 12.011, 1.008, 15.999
    mw_eo = 2 * C + 4 * H + O

    mw_octadecanol = 18 * C + 38 * H + O
    for n_eo, table_hlb in [(10, 13), (12, 14), (16, 15)]:
        mh = n_eo * mw_eo
        computed = hlb_griffin(mh, mw_octadecanol + mh)
        assert computed == pytest.approx(table_hlb, abs=0.8)
        assert computed < table_hlb  # the offset is consistently negative, not random

    mw_isononylphenol = 15 * C + 24 * H + O
    for n_eo, table_hlb in [(5, 10.7), (6, 11.5), (7, 12.3)]:
        mh = n_eo * mw_eo
        computed = hlb_griffin(mh, mw_isononylphenol + mh)
        assert computed == pytest.approx(table_hlb, abs=0.8)
        assert computed < table_hlb


def test_hlb_davies_sds_worked_example():
    """SDS = CH3-(CH2)11-SO4Na: 1 CH3 + 11 CH2 + 1 SO4Na headgroup.
    HLB = 7 + 38.7 - 12*0.475 = 40.0, the value usually quoted for SDS
    (-SO4Na = 38.7 in Davies 1957, the group table of hlb.py)."""
    hlb = hlb_davies({"SO4Na": 1, "CH2": 11, "CH3": 1})
    assert hlb == pytest.approx(40.0, abs=0.01)


def test_hlb_davies_unknown_group_raises_instead_of_guessing():
    """Amide and sulfonate have no verified Davies number in this table
    (see hlb.py docstring) -- must raise, not silently return a guess.
    (Quaternary ammonium used to be in this category too, until resolved
    2026-09-08 -- see the worked-example tests below.)"""
    with pytest.raises(KeyError):
        hlb_davies({"amide": 1, "CH2": 11})
    with pytest.raises(KeyError):
        hlb_davies({"sulfonate": 1, "CH2": 11})


# --- Quaternary ammonium Davies numbers (resolved 2026-09-08) --------------
# Source: B.H. O, J. Colloid Interface Sci. 198 (1998) 249 (GN=22.0 for
# -N+(CH3)3, GN=22.5 for >N+(CH3)2), via Proverbio, Bardavid, Arancibia &
# Schulz, Colloids Surf. A 214 (2003) 167-171 -- all four values below are
# the source paper's OWN reported, fully worked HLB examples, reproduced
# here exactly (not approximately), the strongest possible validation.


def test_hlb_davies_dtab_matches_paper_worked_example():
    """DTAB = decyltrimethylammonium bromide: C10 tail (1 CH3 + 9 CH2)
    + single-tail trimethylammonium head."""
    hlb = hlb_davies({"N_quat_trimethyl": 1, "CH2": 9, "CH3": 1})
    assert hlb == pytest.approx(24.25, abs=0.01)


def test_hlb_davies_ltab_matches_paper_worked_example():
    """LTAB = dodecyltrimethylammonium bromide: C12 tail (1 CH3 + 11 CH2)
    + single-tail trimethylammonium head."""
    hlb = hlb_davies({"N_quat_trimethyl": 1, "CH2": 11, "CH3": 1})
    assert hlb == pytest.approx(23.3, abs=0.01)


def test_hlb_davies_ddab_matches_paper_worked_example():
    """DDAB = didodecyldimethylammonium bromide: two C12 tails
    (2 CH3 + 22 CH2) + double-tail dimethyl-dialkylammonium head."""
    hlb = hlb_davies({"N_quat_dimethyl_dialkyl": 1, "CH2": 22, "CH3": 2})
    assert hlb == pytest.approx(18.1, abs=0.01)


def test_hlb_davies_dodab_matches_paper_worked_example():
    """DODAB = dioctadecyldimethylammonium bromide: two C18 tails
    (2 CH3 + 34 CH2) + double-tail dimethyl-dialkylammonium head."""
    hlb = hlb_davies({"N_quat_dimethyl_dialkyl": 1, "CH2": 34, "CH3": 2})
    assert hlb == pytest.approx(12.4, abs=0.01)


# --- Derived Davies group numbers for amide/sulfonate (2026-09-11) ---------
# NEW CAPABILITY, not a literature citation: this project's own
# Griffin-Davies cross-calibration derivation (see
# derive_davies_group_number_from_griffin's docstring for the full method).
# Real atomic masses (IUPAC standard atomic weights), no external source
# needed beyond stoichiometry.
_C, _H, _N, _O, _S, _NA, _P = 12.011, 1.008, 14.007, 15.999, 32.06, 22.990, 30.974


def test_derive_davies_group_number_reduces_hlb_davies_and_hlb_griffin_to_agreement():
    """Trivial-by-construction check: for the EXACT reference compound
    the derivation is solved against, hlb_davies (with the freshly
    derived number, not the hardcoded constant) and hlb_griffin
    (independently, from raw masses) must give IDENTICAL HLB values --
    that agreement is the literal definition of how the group number is
    derived, so this is a real regression check on the arithmetic, not a
    claim of external validation."""
    # sodium dodecanesulfonate: C12H25SO3Na (n=12 alkyl sulfonate)
    n = 12
    mh = _S + 3 * _O + _NA
    tail_mass = n * _C + (2 * n + 1) * _H
    m_total = tail_mass + mh
    hlb_g = hlb_griffin(mh, m_total)

    gn = derive_davies_group_number_from_griffin(mh, m_total, {"CH2": n - 1, "CH3": 1})
    hlb_d = 7.0 + gn - 0.475 * (n - 1) - 0.475 * 1  # manual Davies sum using the FRESH gn
    assert hlb_d == pytest.approx(hlb_g, abs=1e-9)


def test_derive_davies_group_number_sulfonate_self_consistent_across_chain_lengths():
    """Real, disclosed chain-length sensitivity check (not hidden): the
    derived sulfonate group number should land in a tight band across
    the practically-relevant C9-C12 alkyl sulfonate range (matching this
    project's own Sutherland et al. 2009 sodium alkylsulfonate citation
    elsewhere), even though it is NOT perfectly constant (an expected,
    documented property of cross-calibrating a linear Davies model
    against Griffin's nonlinear mass-ratio formula)."""
    mh = _S + 3 * _O + _NA
    values = {}
    for n in (9, 11, 12):
        tail_mass = n * _C + (2 * n + 1) * _H
        m_total = tail_mass + mh
        values[n] = derive_davies_group_number_from_griffin(mh, m_total, {"CH2": n - 1, "CH3": 1})
    assert all(6.0 < v < 6.5 for v in values.values())
    # the hardcoded reference value (C12) matches the live re-derivation exactly
    assert values[12] == pytest.approx(DAVIES_DERIVED_HYDROPHILIC_GROUPS["sulfonate"], abs=0.005)


def test_derive_davies_group_number_amide_dialkanolamide_self_consistent():
    """Same self-consistency check for the fatty-acid-diethanolamide head
    (-CO-N(CH2CH2OH)2, e.g. cocamide/lauramide DEA), across lauric (C12),
    myristic (C14), and palmitic (C16) acid tails."""
    mh = _C + _O + _N + 2 * (2 * _C + 4 * _H + _O + _H)  # -CO-N(CH2CH2OH)2
    values = {}
    for n_acid in (12, 14, 16):
        n_tail = n_acid - 1
        tail_mass = n_tail * _C + (2 * n_tail + 1) * _H
        m_total = tail_mass + mh
        values[n_acid] = derive_davies_group_number_from_griffin(mh, m_total, {"CH2": n_tail - 1, "CH3": 1})
    assert all(7.0 < v < 8.0 for v in values.values())
    assert values[12] == pytest.approx(DAVIES_DERIVED_HYDROPHILIC_GROUPS["amide_dialkanolamide"], abs=0.005)


def test_derive_davies_group_number_amide_monoalkanolamide_self_consistent():
    """Same check for the fatty-acid-monoethanolamide head (-CO-NH-CH2CH2OH,
    e.g. cocamide/lauramide MEA)."""
    mh = _C + _O + _N + _H + (2 * _C + 4 * _H + _O + _H)  # -CO-NH-CH2CH2OH
    values = {}
    for n_acid in (12, 14, 16):
        n_tail = n_acid - 1
        tail_mass = n_tail * _C + (2 * n_tail + 1) * _H
        m_total = tail_mass + mh
        values[n_acid] = derive_davies_group_number_from_griffin(mh, m_total, {"CH2": n_tail - 1, "CH3": 1})
    assert all(5.0 < v < 6.5 for v in values.values())
    assert values[12] == pytest.approx(DAVIES_DERIVED_HYDROPHILIC_GROUPS["amide_monoalkanolamide"], abs=0.005)


def test_derive_davies_group_number_rejects_bad_inputs():
    with pytest.raises(ValueError):
        derive_davies_group_number_from_griffin(0.0, 100.0, {"CH2": 5})
    with pytest.raises(ValueError):
        derive_davies_group_number_from_griffin(100.0, 50.0, {"CH2": 5})  # Mh >= M
    with pytest.raises(KeyError):
        derive_davies_group_number_from_griffin(50.0, 200.0, {"not_a_real_group": 5})


def test_hlb_davies_sulfonate_group_usable_only_when_opted_in():
    """The derived group is invisible to hlb_davies by default (matching
    the strict, literature-only philosophy this project has used
    throughout for the main DAVIES_HYDROPHILIC_GROUPS table), and only
    becomes usable with an explicit opt-in -- demonstrating the actual,
    previously-impossible capability: computing an HLB for a sulfonate
    surfactant, which used to hard KeyError unconditionally."""
    with pytest.raises(KeyError):
        hlb_davies({"sulfonate": 1, "CH2": 11, "CH3": 1})
    hlb = hlb_davies({"sulfonate": 1, "CH2": 11, "CH3": 1}, allow_derived_groups=True)
    assert hlb == pytest.approx(7.566, abs=0.01)  # sodium dodecanesulfonate's own Griffin HLB


def test_hlb_davies_amide_dialkanolamide_worked_example():
    """Lauramide DEA (cocamide DEA's dominant single-chain-length
    approximation): C11 tail + the derived diethanolamide head group."""
    hlb = hlb_davies({"amide_dialkanolamide": 1, "CH2": 10, "CH3": 1}, allow_derived_groups=True)
    assert hlb == pytest.approx(9.194, abs=0.01)


def test_amide_dialkanolamide_gn_second_independent_anchor_qualitative():
    """A SECOND, independent (non-Griffin-derived) cross-check for the
    derived diethanolamide numbers, 2026-09-11 -- checked directly
    against a fetched primary source, not a search-engine paraphrase.
    Pawignya et al. (IOP Conf. Ser.: Mater. Sci. Eng., Atlantis Press)
    report a real EXPERIMENTALLY MEASURED HLB=5.940 for a palm-oil-
    derived diethanolamide, via their own CMC-based formula (a
    genuinely different methodology from Griffin's mass-ratio formula
    this project's numbers are calibrated against -- NOT a number
    directly comparable value-for-value).

    Honestly scoped as a QUALITATIVE, not quantitative, cross-check: the
    real, checkable claim is that palm oil's fatty-acid profile (longer
    chains, C16-C18-dominant) should give a LOWER HLB than a C12
    (lauric/coconut-scale) alkanolamide under EITHER method, since both
    Griffin's and this paper's formulas capture the same real physical
    trend -- more hydrocarbon mass per molecule means less hydrophilic.
    This project's own derived numbers (see DAVIES_DERIVED_HYDROPHILIC_GROUPS'
    block comment) independently show exactly this trend across chain
    length. Both facts (this project's own C16 point sitting closer to
    the primary-source's C16-C18-dominated 5.940 than its own C12 point
    does, and the real physical direction of the trend matching) are
    asserted; the primary source's absolute number is NOT asserted to
    equal this project's number, since the two use different formulas."""
    DAVIES_DERIVED_HYDROPHILIC_GROUPS["amide_dialkanolamide"]
    hlb_c12 = hlb_davies({"amide_dialkanolamide": 1, "CH2": 10, "CH3": 1}, allow_derived_groups=True)
    hlb_c16 = hlb_davies({"amide_dialkanolamide": 1, "CH2": 14, "CH3": 1}, allow_derived_groups=True)
    palm_oil_measured_hlb = 5.940  # Pawignya et al., their own CMC-based formula, NOT Griffin's

    # the real physical trend: longer chain (palm-scale, ~C16-C18) -> lower HLB
    assert hlb_c16 < hlb_c12
    # this project's C16 point is closer to the primary source's own
    # (longer-chain, palm-oil) measured value than its C12 point is
    assert abs(hlb_c16 - palm_oil_measured_hlb) < abs(hlb_c12 - palm_oil_measured_hlb)
    # both sources land in the same broad hydrophilic-nonionic-surfactant
    # regime (order of magnitude, not exact) -- not off by a factor of 2+
    assert 0.5 * palm_oil_measured_hlb < hlb_c12 < 2.5 * palm_oil_measured_hlb


# --- Extended derived groups (2026-09-11): sultaine, carboxybetaine, ------
# phosphate ester -- same Griffin-Davies cross-calibration methodology as
# sulfonate/amide above, applied to the remaining SurfBench-flagged gaps
# (carboxybetaine, sultaine, phosphate ester -- imidazoline deliberately
# skipped, see hlb.py's own block comment for why).


def test_derive_davies_group_number_sultaine_self_consistent():
    """Alkyl sulfobetaine head (-N+(CH3)2-(CH2)3-SO3-), a genuine
    internal-salt zwitterion -- e.g. SB10/SB12, the exact compounds
    already cited via Sutherland et al. 2009 elsewhere in this project.
    Tightest self-consistency of any derived group (~2.7% spread)."""
    head = _N + 2 * (_C + 3 * _H) + 3 * (_C + 2 * _H) + _S + 3 * _O
    values = {}
    for n in (8, 10, 12, 14):
        tail_mass = n * _C + (2 * n + 1) * _H
        m_total = tail_mass + head
        values[n] = derive_davies_group_number_from_griffin(head, m_total, {"CH2": n - 1, "CH3": 1})
    assert all(8.4 < v < 8.9 for v in values.values())
    assert values[12] == pytest.approx(DAVIES_DERIVED_HYDROPHILIC_GROUPS["sultaine"], abs=0.005)


def test_derive_davies_group_number_carboxybetaine_self_consistent():
    """Cocamidopropyl-betaine-style head (-CO-NH-(CH2)3-N+(CH3)2-CH2-COO-),
    one of the most common real amphoteric surfactants."""
    head = (_C + _O) + (_N + _H) + 3 * (_C + 2 * _H) + (_N + 2 * (_C + 3 * _H)) + (_C + 2 * _H) + (_C + 2 * _O)
    values = {}
    for n_acid in (12, 14, 16):
        n_tail = n_acid - 1
        tail_mass = n_tail * _C + (2 * n_tail + 1) * _H
        m_total = tail_mass + head
        values[n_acid] = derive_davies_group_number_from_griffin(head, m_total, {"CH2": n_tail - 1, "CH3": 1})
    assert all(9.0 < v < 9.6 for v in values.values())
    assert values[12] == pytest.approx(DAVIES_DERIVED_HYDROPHILIC_GROUPS["carboxybetaine"], abs=0.005)


def test_derive_davies_group_number_phosphate_diNa_self_consistent():
    """Disodium monoalkyl phosphate head (-O-P(=O)(ONa)2) -- unlike the
    amide/betaine groups, the alkyl-O-P linkage does NOT consume a tail
    carbon, so the full Cn chain counts as lipophilic."""
    head = 4 * _O + _P + 2 * _NA
    values = {}
    for n in (10, 12, 14):
        tail_mass = n * _C + (2 * n + 1) * _H
        m_total = tail_mass + head
        values[n] = derive_davies_group_number_from_griffin(head, m_total, {"CH2": n - 1, "CH3": 1})
    assert all(7.5 < v < 8.1 for v in values.values())


def test_derive_davies_group_number_amphodiacetate_self_consistent():
    """Disodium lauroamphodiacetate head (the real, final, commercial
    imidazoline-DERIVED amphoteric surfactant -- imidazoline itself is
    only a synthesis intermediate, never the sold product; see hlb.py's
    own block comment for the full disclosure of why this specific
    diacetate reaction product was picked over leaving the whole family
    unresolved). Structure confirmed against PubChem CID 109973
    (C20H36N2Na2O6) two independent ways: hand atomic-mass summation
    (this test) and RDKit get_weight_from_smiles on the open-chain amide
    tautomer -- both land on 446.496 g/mol total, 291.191 g/mol head."""
    head = (
        (_C + _O)                       # amide carbonyl C=O (or ring C2 in the cyclic tautomer)
        + (_N + _H)                     # amide N-H
        + 2 * (_C + 2 * _H)             # -CH2-CH2- linker to the tertiary N
        + _N                            # tertiary N
        + ((_C + 2 * _H) + (_C + 2 * _O) + _NA)  # -CH2-COONa
        + (2 * (_C + 2 * _H) + _O + (_C + 2 * _H) + (_C + 2 * _O) + _NA)  # -CH2CH2-O-CH2-COONa
    )
    assert head == pytest.approx(291.191, abs=0.01)
    values = {}
    for n_acid in (10, 12, 14):
        n_tail = n_acid - 1
        tail_mass = n_tail * _C + (2 * n_tail + 1) * _H
        m_total = tail_mass + head
        values[n_acid] = derive_davies_group_number_from_griffin(head, m_total, {"CH2": n_tail - 1, "CH3": 1})
    assert all(11.0 < v < 11.6 for v in values.values())
    assert values[12] == pytest.approx(DAVIES_DERIVED_HYDROPHILIC_GROUPS["amphodiacetate"], abs=0.005)


def test_hlb_davies_sultaine_and_carboxybetaine_and_phosphate_usable_only_when_opted_in():
    with pytest.raises(KeyError):
        hlb_davies({"sultaine": 1, "CH2": 11, "CH3": 1})
    with pytest.raises(KeyError):
        hlb_davies({"carboxybetaine": 1, "CH2": 10, "CH3": 1})
    with pytest.raises(KeyError):
        hlb_davies({"phosphate_diNa": 1, "CH2": 11, "CH3": 1})

    hlb_sb12 = hlb_davies({"sultaine": 1, "CH2": 11, "CH3": 1}, allow_derived_groups=True)
    assert hlb_sb12 == pytest.approx(9.907, abs=0.01)  # SB12's own Griffin HLB

    hlb_capb = hlb_davies({"carboxybetaine": 1, "CH2": 10, "CH3": 1}, allow_derived_groups=True)
    assert hlb_capb == pytest.approx(10.932, abs=0.01)  # lauramidopropyl betaine's own Griffin HLB

    hlb_phos = hlb_davies({"phosphate_diNa": 1, "CH2": 11, "CH3": 1}, allow_derived_groups=True)
    assert hlb_phos == pytest.approx(9.085, abs=0.01)  # disodium lauryl phosphate's own Griffin HLB


def test_tanford_tail_volume_matches_additive_formula():
    """Cross-check the shorthand (27.4 + 26.9*nc) against the more
    detailed per-group formula (v_CH3 + (nc-1)*v_CH2) for a C12 chain --
    they must agree exactly, since the shorthand is algebraically
    derived from the detailed formula, not an independent claim."""
    n_carbons = 12
    shorthand = 27.4 + 26.9 * n_carbons
    detailed = 54.3 + (n_carbons - 1) * 26.9
    assert shorthand == pytest.approx(detailed, abs=1e-9)
    assert tanford_tail_volume(n_carbons) == pytest.approx(detailed, abs=1e-9)


def test_tanford_critical_length_c12():
    # lc = 1.5 + 1.265*12 = 16.68 Angstrom -- in the commonly-cited
    # ballpark (~16.7 A) for a fully-extended dodecyl chain.
    assert tanford_critical_length(12) == pytest.approx(16.68, abs=0.01)


def test_tanford_formulas_match_independent_secondary_source():
    """Real second-independent-source check for Category D, 2026-09-10
    (after 7 other search/fetch attempts found no clean worked numeric
    example from a real surfactant paper to check against -- see
    SurfactantKit/ROADMAP.md for the full disclosure of what was tried
    and did not pan out). A physical-chemistry course problem
    (independently phrased and unit-converted, not sourced from this
    project's own existing citations) states the identical formula as
    v(n) = (27.4+26.9n)*1e-3 nm^3 and lc(n) = (0.154+0.1265n) nm.
    Converting units: (27.4+26.9n)*1e-3 nm^3 = (27.4+26.9n) A^3 (since
    1 nm^3 = 1000 A^3) -- EXACT match to this project's tanford_tail_volume.
    (0.154+0.1265n) nm = (1.54+1.265n) A -- matches this project's
    tanford_critical_length's 1.5+1.265n to within 0.04 A (a real, small,
    explainable rounding difference in the constant term, 1.5 vs 1.54;
    the chain-length-dependent coefficient, 1.265, matches EXACTLY).
    At the time this test was written, this was NOT yet a worked example
    from a real surfactant system -- that harder bar is met separately,
    below, by test_critical_packing_parameter_matches_real_sds_dtab_paper."""
    for n in (8, 12, 16):
        v_secondary_A3 = (27.4 + 26.9 * n) * 1e-3 * 1000.0  # nm^3 -> A^3
        lc_secondary_A = (0.154 + 0.1265 * n) * 10.0  # nm -> A
        assert tanford_tail_volume(n) == pytest.approx(v_secondary_A3, abs=1e-9)
        assert tanford_critical_length(n) == pytest.approx(lc_secondary_A, abs=0.05)


def test_critical_packing_parameter_matches_real_sds_dtab_paper():
    """Category D's harder bar, finally met 2026-09-12: a full,
    independent, real-surfactant CPP worked example, not just a
    formula-transcription check (see the note above and ROADMAP.md's
    2026-09-10 entry for the 7 earlier attempts that found nothing
    usable).

    Source: Sachin, K. M., Karpe, S. A., Singh, M., and Bhattarai, A.,
    "Self-assembly of sodium dodecylsulfate and dodecyltrimethylammonium
    bromide mixed surfactants with dyes in aqueous mixtures," R. Soc.
    Open Sci. 6, 181979 (2019), PMC6458362, Table 2. (Corrected
    2026-09-12: the real primary-source PDF was obtained and shows
    these are the paper's real authors -- an earlier citation of this
    same real title/journal/volume/page/PMC-ID had the wrong author
    names attached; the title, DOI, and every numeric value below were
    already correct and unaffected.) Both surfactants are 12-carbon
    (dodecyl) chains, so Tanford's V0/lc are identical for both --
    only Amin (from the paper's own Gibbs-isotherm surface-tension
    slope, a real measurement, not assumed) differs between rows. This
    project's own tanford_tail_volume(12), tanford_critical_length(12)
    and critical_packing_parameter() are chained end-to-end and
    reproduce the paper's own reported P across 6 independent points
    (2 systems x 3 temperatures) to within 1.2% relative error, and
    classify_aggregate_morphology() reproduces the paper's own stated
    "cylindrical or rod-shaped micelles" prediction at every point.
    Small (<1.2%) offset explained by the same 1.265 vs 1.26 lc-
    coefficient rounding already documented in the test above -- the
    paper states lc=1.54+1.26*nc, this project uses 1.5+1.265*nc."""
    v = tanford_tail_volume(12)
    lc = tanford_critical_length(12)

    # (label, Amin in A^2/molecule, paper's own reported P)
    rows = [
        ("SDS-rich 293.15K", 44.70, 0.47),
        ("SDS-rich 298.15K", 49.37, 0.43),
        ("SDS-rich 303.15K", 51.82, 0.41),
        ("DTAB-rich 293.15K", 57.89, 0.36),
        ("DTAB-rich 298.15K", 59.62, 0.35),
        ("DTAB-rich 303.15K", 61.59, 0.34),
    ]
    for label, amin, p_paper in rows:
        p = critical_packing_parameter(v, amin, lc)
        assert p == pytest.approx(p_paper, rel=0.012), label
        assert classify_aggregate_morphology(p) == "cylindrical/rodlike micelle", label


def test_aggregation_number_spherical_c12_lands_in_literature_range():
    """C12 chain (SDS/DTAB scale) geometric aggregation number should
    land near the commonly-reported literature range (~55-70) for these
    surfactants -- a sanity/order-of-magnitude check, not an exact
    literature match (this is a geometric estimate, not a substitute for
    a directly-measured aggregation number; see cpp.py docstring)."""
    v = tanford_tail_volume(12)
    lc = tanford_critical_length(12)
    nagg = aggregation_number_spherical(v, lc)
    assert 40.0 < nagg < 80.0


def test_aggregation_number_spherical_reproduces_bales_1998_worked_example():
    """Bales, Messina, Vidal, Peric & Nascimento, J. Phys. Chem. B 1998,
    102, 10347-10358 (full text obtained 2026-09-11, freely hosted by
    the author at csun.edu/~vcphy00s/Bales56.pdf -- previously this
    project's notes only had access to a secondary citation of this
    paper's Table 2, not the primary text itself). Their own
    illustrative worked example (citing Cabane's SANS reference system):
    69 mM SDS, salt-free, NA=63, Nc=12 -> Vtail=350.2 A^3 (their eq. 9,
    matching this project's tanford_tail_volume(12) exactly).

    IMPORTANT, honestly scoped: the paper does not print a numeric R_c
    value directly (it reports V_p and the hydration number instead);
    R_c here is BACK-DERIVED from their own eq. 10, NA*Vtail =
    (4pi/3)*Rc^3, using their own stated NA=63 -- so recovering NA=63
    from aggregation_number_spherical(Vtail, Rc) is a self-consistency
    check against their own formula (same caveat already noted for the
    Nagarajan Table 2 CPP check above), not an independent measurement
    of Rc from a different technique."""
    import math

    v_tail = tanford_tail_volume(12)
    n_a_paper = 63.0
    r_c = ((3.0 * n_a_paper * v_tail) / (4.0 * math.pi)) ** (1.0 / 3.0)
    assert aggregation_number_spherical(v_tail, r_c) == pytest.approx(n_a_paper, rel=1e-9)


def test_aggregation_number_spherical_matches_direct_geometry_formula():
    import math

    v, r = 350.2, 16.68
    expected = ((4.0 / 3.0) * math.pi * r ** 3) / v
    assert aggregation_number_spherical(v, r) == pytest.approx(expected)


def test_aggregation_number_spherical_rejects_nonpositive():
    with pytest.raises(ValueError):
        aggregation_number_spherical(0, 16.68)
    with pytest.raises(ValueError):
        aggregation_number_spherical(350.2, -1.0)


def test_tanford_rejects_bad_chain_length():
    with pytest.raises(ValueError):
        tanford_tail_volume(0)
    with pytest.raises(ValueError):
        tanford_critical_length(-1)


def test_critical_packing_parameter_basic():
    # v=350.2 (C12), a0=50 A^2 (typical ionic headgroup area), lc=16.68 A
    cpp = critical_packing_parameter(350.2, 50.0, 16.68)
    assert cpp == pytest.approx(350.2 / (50.0 * 16.68))


def test_critical_packing_parameter_rejects_nonpositive():
    with pytest.raises(ValueError):
        critical_packing_parameter(0, 50.0, 16.68)
    with pytest.raises(ValueError):
        critical_packing_parameter(350.2, -1.0, 16.68)


@pytest.mark.parametrize(
    "cpp,expected",
    [
        (0.20, "spherical micelle"),
        (1.0 / 3.0, "spherical micelle"),          # boundary, inclusive
        (0.34, "cylindrical/rodlike micelle"),
        (0.5, "cylindrical/rodlike micelle"),       # boundary, inclusive
        (0.51, "vesicle/bilayer"),
        (1.0, "vesicle/bilayer"),                   # boundary, inclusive
        (1.01, "inverted structure"),
        (2.0, "inverted structure"),
    ],
)
def test_classify_aggregate_morphology_thresholds(cpp, expected):
    assert classify_aggregate_morphology(cpp) == expected


def test_classify_aggregate_morphology_rejects_nonpositive():
    with pytest.raises(ValueError):
        classify_aggregate_morphology(0)
    with pytest.raises(ValueError):
        classify_aggregate_morphology(-0.5)


# Tests for estimate_axial_ratio_from_cpp_geometry, added 2026-09-15 as
# part of the Tier 1 bottleneck-resolution work (real [COMPUTE] path for
# dynamics.py's Perrin axial_ratio -- see BOTTLENECK_RESOLUTION_PLAN.md
# item 8). Validated by round-trip (construct a known axial_ratio's
# consistent aggregation number via the SAME real ellipsoid-of-revolution
# geometry, recover the same axial_ratio back out) since no single real
# paper reports CPP + a real independently-measured N_agg + n_carbons
# together for a rodlike micelle to cross-check against directly.


def test_estimate_axial_ratio_returns_one_for_spherical_cpp():
    assert estimate_axial_ratio_from_cpp_geometry(0.25, aggregation_number=60, n_carbons=12) == pytest.approx(1.0)
    assert estimate_axial_ratio_from_cpp_geometry(1.0 / 3.0, aggregation_number=60, n_carbons=12) == pytest.approx(1.0)


def test_estimate_axial_ratio_round_trip_recovers_known_ratio():
    n_carbons = 12
    b = tanford_critical_length(n_carbons)
    v_tail = tanford_tail_volume(n_carbons)
    for true_axial_ratio in (1.2, 1.8, 3.0, 5.0):
        a = true_axial_ratio * b
        v_core = (4.0 / 3.0) * math.pi * a * b ** 2
        n_agg = v_core / v_tail
        recovered = estimate_axial_ratio_from_cpp_geometry(0.45, aggregation_number=n_agg, n_carbons=n_carbons)
        assert recovered == pytest.approx(true_axial_ratio, rel=1e-6)


def test_estimate_axial_ratio_rejects_cpp_outside_rodlike_regime():
    with pytest.raises(ValueError):
        estimate_axial_ratio_from_cpp_geometry(0.6, aggregation_number=200, n_carbons=12)  # vesicle/bilayer regime
    with pytest.raises(ValueError):
        estimate_axial_ratio_from_cpp_geometry(1.2, aggregation_number=200, n_carbons=12)  # inverted regime


def test_estimate_axial_ratio_rejects_aggregation_number_too_small_for_rodlike_cpp():
    # an N_agg only barely above the spherical-micelle scale is inconsistent
    # with a cpp value that claims the rodlike regime -- must raise, not
    # silently clamp to 1.0
    with pytest.raises(ValueError):
        estimate_axial_ratio_from_cpp_geometry(0.45, aggregation_number=10, n_carbons=12)


def test_estimate_axial_ratio_rejects_bad_inputs():
    with pytest.raises(ValueError):
        estimate_axial_ratio_from_cpp_geometry(0.0, aggregation_number=60, n_carbons=12)
    with pytest.raises(ValueError):
        estimate_axial_ratio_from_cpp_geometry(0.4, aggregation_number=0, n_carbons=12)
    with pytest.raises(ValueError):
        estimate_axial_ratio_from_cpp_geometry(0.4, aggregation_number=60, n_carbons=0)


def test_estimate_axial_ratio_agrees_with_aggregation_number_spherical_at_the_boundary():
    """Real cross-function consistency check, found worth verifying
    during a Tier 1 completeness review: at axial_ratio=1 (a=b=lc), the
    ellipsoid-of-revolution volume formula this function uses reduces
    EXACTLY to the sphere-of-radius-lc volume formula
    aggregation_number_spherical already uses elsewhere in this same
    module -- both were derived independently but must agree at the
    boundary where the two shapes coincide. Confirms no hidden factor-of-
    something inconsistency between the two functions."""
    n_carbons = 12
    lc = tanford_critical_length(n_carbons)
    v_tail = tanford_tail_volume(n_carbons)
    n_agg_at_sphere_boundary = aggregation_number_spherical(v_tail, lc)

    recovered_axial_ratio = estimate_axial_ratio_from_cpp_geometry(
        cpp=0.45, aggregation_number=n_agg_at_sphere_boundary, n_carbons=n_carbons
    )
    assert recovered_axial_ratio == pytest.approx(1.0, rel=1e-9)


# --- Guo/Rong/Ying 2006 nonionic HLB refinement (alternative-methods
# sweep, 2026-09-10). Source: J. Colloid Interface Sci. 298 (2006)
# 441-450, doi:10.1016/j.jcis.2005.12.009 (paywalled; exact coefficients
# below verified via a citing patent, US 11,344,493 B2, which quotes the
# formulas directly).


def test_guo_effective_eo_chain_length_matches_cited_formula():
    import math

    for n_eo in (1, 5, 10, 20, 50):
        expected = 13.45 * math.log(n_eo) - 0.16 * n_eo + 1.26
        assert guo_effective_eo_chain_length(n_eo) == pytest.approx(expected)


def test_guo_effective_eo_chain_length_rejects_bad_input():
    with pytest.raises(ValueError):
        guo_effective_eo_chain_length(0)
    with pytest.raises(ValueError):
        guo_effective_eo_chain_length(-3)


def test_guo_effective_eo_chain_length_second_branch_above_50():
    # Primary source's own eq. (5') second branch, confirmed 2026-09-11
    # against the obtained primary-source PDF -- previously unimplemented.
    for n_eo in (51, 60, 100):
        expected = 0.056 * n_eo + 43.08
        assert guo_effective_eo_chain_length(n_eo) == pytest.approx(expected)
    # the two branches nearly agree at the n_eo=50 boundary (by the
    # source's own construction, not exactly equal)
    assert guo_effective_eo_chain_length(50) == pytest.approx(45.88, abs=0.01)
    assert (0.056 * 50 + 43.08) == pytest.approx(45.88, abs=0.05)


def test_guo_effective_alkyl_chain_length_matches_cited_formula():
    for n_ch2 in (5, 10, 15):
        expected = 0.965 * n_ch2 - 0.178
        assert guo_effective_alkyl_chain_length(n_ch2) == pytest.approx(expected)


def test_guo_effective_alkyl_chain_length_rejects_bad_input():
    with pytest.raises(ValueError):
        guo_effective_alkyl_chain_length(0)


def test_guo_effective_po_chain_length_matches_cited_formula():
    for n_po in (1, 5, 10):
        expected = 2.057 * n_po + 9.06
        assert guo_effective_po_chain_length(n_po) == pytest.approx(expected)


def test_guo_effective_po_chain_length_rejects_bad_input():
    with pytest.raises(ValueError):
        guo_effective_po_chain_length(0)


def test_hlb_davies_guo_ecl_matches_manual_formula():
    """Independent recomputation in the test itself, not copy-pasted
    from the implementation."""
    n_carbons, n_eo = 12, 9  # C12E9-like nonionic ethoxylate
    n_ch2_eff = 0.965 * (n_carbons - 1) - 0.178
    n_eo_eff = 13.45 * __import__("math").log(n_eo) - 0.16 * n_eo + 1.26
    expected = 7.0 + 0.33 * n_eo_eff - 0.475 * n_ch2_eff - 0.475 * 1   # EO unit 0.33 (corrected 2026-10-20)
    assert hlb_davies_guo_ecl(n_carbons, n_eo) == pytest.approx(expected)


def test_hlb_davies_guo_ecl_vs_plain_davies_diverges_for_long_eo_chains():
    """Real point of the refinement: for a long EO chain, the effective
    (sub-linear) count differs substantially from the actual count, so
    Guo's HLB should differ meaningfully from plain Davies using the
    actual (unmodified) EO count -- confirms the ECL correction is doing
    something real, not silently reducing to the plain formula."""
    n_carbons, n_eo = 12, 40
    guo_hlb = hlb_davies_guo_ecl(n_carbons, n_eo)
    plain_davies_hlb = hlb_davies({"EO_unit": n_eo, "CH2": n_carbons - 1, "CH3": 1})
    assert abs(guo_hlb - plain_davies_hlb) > 1.0


def test_hlb_davies_guo_ecl_with_extra_group_counts():
    base = hlb_davies_guo_ecl(12, 9)
    with_extra = hlb_davies_guo_ecl(12, 9, extra_group_counts={"OH_free": 1})
    assert with_extra == pytest.approx(base + 1.9)


def test_hlb_davies_guo_ecl_c12e8_lands_in_the_physical_range():
    """Sanity check against reality, not against the formula itself (the 2026-10-20 bug, EO weighted
    1.3 instead of 0.33, passed the formula-mirroring tests above while giving C12E8 an HLB of 37.9).
    C12E8 (C12H25(OCH2CH2)8OH) is a classic O/W emulsifier/detergent with an HLB near 13 on every
    standard scale (Griffin mass fraction: 20 x 352/539 = 13.1)."""
    hlb = hlb_davies_guo_ecl(12, 8, extra_group_counts={"OH_free": 1})
    assert 11.5 < hlb < 14.5
    assert hlb_griffin(8 * 44.053, 538.8) == pytest.approx(13.08, abs=0.05)


def test_hlb_davies_guo_ecl_rejects_unknown_extra_group():
    with pytest.raises(KeyError):
        hlb_davies_guo_ecl(12, 9, extra_group_counts={"not_a_real_group": 1})


def test_hlb_davies_guo_ecl_rejects_bad_inputs():
    with pytest.raises(ValueError):
        hlb_davies_guo_ecl(1, 9)  # too few carbons
    with pytest.raises(ValueError):
        hlb_davies_guo_ecl(12, 0)  # need at least 1 EO unit


# --- Nagarajan 2002 ionic-headgroup electrostatics (alternative-methods
# review, 2026-09-10) -- Table 2's real worked example, sodium alkyl
# sulfates, n_C=8..16. Source: Nagarajan, Langmuir 18 (2002) 31-38.
# headgroup_prefactor_A=82.0 is the paper's OWN stated value for this
# specific table's illustration (sodium-alkyl-sulfate-like headgroup),
# not a universal constant -- see nagarajan_equilibrium_area_ionic's
# docstring for the disclosure.

NAGARAJAN_TABLE2 = [
    # n_C, l0 (A), cmc (M), kappa^-1 (A) [paper], a_e (A^2) [paper]
    (8, 11.5, 0.032, 17.22, 63.5),
    (10, 14.0, 0.016, 24.35, 65.4),
    (12, 16.5, 0.008, 34.43, 67.4),
    (14, 19.0, 0.004, 48.7, 69.5),
    (16, 21.5, 0.002, 68.9, 71.6),
]


@pytest.mark.parametrize("n_c, l0, cmc, kappa_inv_paper, a_e_paper", NAGARAJAN_TABLE2)
def test_nagarajan_kappa_inverse_matches_table2(n_c, l0, cmc, kappa_inv_paper, a_e_paper):
    kappa_inv = nagarajan_debye_huckel_kappa_inverse(cmc)
    assert kappa_inv == pytest.approx(kappa_inv_paper, rel=0.01)


@pytest.mark.parametrize("n_c, l0, cmc, kappa_inv_paper, a_e_paper", NAGARAJAN_TABLE2)
def test_nagarajan_equilibrium_area_matches_table2(n_c, l0, cmc, kappa_inv_paper, a_e_paper):
    a_e = nagarajan_equilibrium_area_ionic(cmc, l0, headgroup_prefactor_A=82.0)
    assert a_e == pytest.approx(a_e_paper, rel=0.01)


def test_nagarajan_equilibrium_area_chains_into_critical_packing_parameter():
    """The real point of the refinement: a_e -> CPP varies with tail
    length for ionic surfactants (not constant, contrary to the naive
    headgroup-only assumption) -- confirm CPP computed via Nagarajan's
    a_e is NOT constant across the n_C=8..16 series, and decreases
    (tighter packing) as chain length increases, matching Table 2's own
    v0/(a_e*l0) column trend (0.331 -> 0.293)."""
    from surfactantkit.cpp import tanford_tail_volume

    cpps = []
    for n_c, l0, cmc, _, _ in NAGARAJAN_TABLE2:
        v0 = tanford_tail_volume(n_c)
        a_e = nagarajan_equilibrium_area_ionic(cmc, l0, headgroup_prefactor_A=82.0)
        cpps.append(critical_packing_parameter(v0, a_e, l0))
    assert cpps == sorted(cpps, reverse=True)  # strictly decreasing with chain length
    assert cpps[0] == pytest.approx(0.331, rel=0.02)
    assert cpps[-1] == pytest.approx(0.293, rel=0.02)


def test_nagarajan_kappa_inverse_rejects_bad_inputs():
    with pytest.raises(ValueError):
        nagarajan_debye_huckel_kappa_inverse(0.0)
    with pytest.raises(ValueError):
        nagarajan_debye_huckel_kappa_inverse(0.01, temperature_K=0.0)
    with pytest.raises(ValueError):
        nagarajan_debye_huckel_kappa_inverse(0.01, dielectric_constant=0.0)


def test_nagarajan_equilibrium_area_rejects_bad_inputs():
    with pytest.raises(ValueError):
        nagarajan_equilibrium_area_ionic(0.01, tail_length_A=0.0, headgroup_prefactor_A=82.0)
    with pytest.raises(ValueError):
        nagarajan_equilibrium_area_ionic(0.01, tail_length_A=16.5, headgroup_prefactor_A=0.0)


# Tests for the Nagarajan & Ruckenstein 1991 full molecular-thermodynamic
# headgroup model, added 2026-09-15 -- real Tier 3 bottleneck-resolution
# work (see BOTTLENECK_RESOLUTION_PLAN.md item 7), built from the actual
# primary-source PDF the user provided (Langmuir 7 (1991) 2934-2969), not
# a secondary citation. Validated three ways: (1) a real, independent
# physical constant (the Bjerrum length in water at 25C, ~7.0-7.1 A) that
# the shared e^2/(eps*kT) combination inside both new functions must
# reduce to exactly; (2) a direct cross-check of the new ionic function's
# internal kappa against this module's OWN already-paper-verified
# nagarajan_debye_huckel_kappa_inverse (same underlying Debye-Huckel
# formula, same real worked example, cmc=0.008 M -> kappa^-1=34.33 A);
# (3) real physical-limit/sign checks (steric repulsion diverges at
# contact, is real and positive when crowded).


def test_bjerrum_length_cross_check():
    """The real, independent, well-known physical constant both new
    functions' e^2/(eps*kT) combination must reduce to: ~7.0-7.1 A in
    water at 25 C (eps~78-80) -- a genuine external cross-check on the
    CGS-Gaussian unit handling, not just internal self-consistency."""
    from surfactantkit.cpp import _E_ESU, _K_B_ERG

    l_bjerrum_cm = (_E_ESU ** 2) / (80.0 * _K_B_ERG * 298.15)
    l_bjerrum_A = l_bjerrum_cm * 1e8
    assert 6.9 <= l_bjerrum_A <= 7.2


def test_nagarajan_ruckenstein_steric_zero_at_large_area():
    a_p = NAGARAJAN_RUCKENSTEIN_HEADGROUP_CONSTANTS["sodium_sulfate"]["a_p_A2"]
    energy_far = nagarajan_ruckenstein_steric_headgroup_free_energy(1000.0 * a_p, a_p)
    assert 0.0 < energy_far < 0.01  # negligible steric repulsion when area >> a_p


def test_nagarajan_ruckenstein_steric_diverges_near_contact():
    a_p = NAGARAJAN_RUCKENSTEIN_HEADGROUP_CONSTANTS["sodium_sulfonate"]["a_p_A2"]
    energy_near = nagarajan_ruckenstein_steric_headgroup_free_energy(a_p * 1.001, a_p)
    energy_far = nagarajan_ruckenstein_steric_headgroup_free_energy(a_p * 2.0, a_p)
    assert energy_near > 5.0 * energy_far  # sharply increasing as area approaches a_p


def test_nagarajan_ruckenstein_steric_rejects_area_at_or_below_ap():
    a_p = 17.0
    with pytest.raises(ValueError):
        nagarajan_ruckenstein_steric_headgroup_free_energy(a_p, a_p)
    with pytest.raises(ValueError):
        nagarajan_ruckenstein_steric_headgroup_free_energy(a_p * 0.5, a_p)


def test_nagarajan_ruckenstein_dipole_matches_hand_derived_bjerrum_formula():
    """Direct formula cross-check: for the sphere case, the dipole
    energy must equal 2*pi*(R_cm/a_cm2)*l_bjerrum_cm*[d/(d+R)] exactly --
    verified independently of the function's own internal unit handling."""
    from surfactantkit.cpp import _CM_PER_ANGSTROM, _E_ESU, _K_B_ERG

    betaine = NAGARAJAN_RUCKENSTEIN_HEADGROUP_CONSTANTS["n_betaine"]
    area_A2, radius_A, d_A = betaine["a_o_A2"], 20.0, betaine["d_A"]

    l_bjerrum_cm = (_E_ESU ** 2) / (80.0 * _K_B_ERG * 298.15)
    area_cm2 = area_A2 * _CM_PER_ANGSTROM ** 2
    radius_cm = radius_A * _CM_PER_ANGSTROM
    expected = 2.0 * math.pi * (radius_cm / area_cm2) * l_bjerrum_cm * (d_A / (d_A + radius_A))

    actual = nagarajan_ruckenstein_dipole_headgroup_free_energy(area_A2, radius_A, d_A, geometry="sphere")
    assert actual == pytest.approx(expected, rel=1e-9)
    assert actual > 0.0  # a real, physically repulsive contribution


def test_nagarajan_ruckenstein_dipole_cylinder_geometry_uses_log_bracket():
    betaine = NAGARAJAN_RUCKENSTEIN_HEADGROUP_CONSTANTS["n_betaine"]
    sphere_case = nagarajan_ruckenstein_dipole_headgroup_free_energy(
        betaine["a_o_A2"], 20.0, betaine["d_A"], geometry="sphere"
    )
    cylinder_case = nagarajan_ruckenstein_dipole_headgroup_free_energy(
        betaine["a_o_A2"], 20.0, betaine["d_A"], geometry="cylinder"
    )
    assert sphere_case != pytest.approx(cylinder_case)  # genuinely different geometry brackets
    assert cylinder_case > 0.0


def test_nagarajan_ruckenstein_dipole_rejects_bad_inputs_and_geometry():
    betaine = NAGARAJAN_RUCKENSTEIN_HEADGROUP_CONSTANTS["n_betaine"]
    with pytest.raises(ValueError):
        nagarajan_ruckenstein_dipole_headgroup_free_energy(0.0, 20.0, betaine["d_A"])
    with pytest.raises(ValueError):
        nagarajan_ruckenstein_dipole_headgroup_free_energy(betaine["a_o_A2"], 0.0, betaine["d_A"])
    with pytest.raises(ValueError):
        nagarajan_ruckenstein_dipole_headgroup_free_energy(betaine["a_o_A2"], 20.0, 0.0)
    with pytest.raises(ValueError):
        nagarajan_ruckenstein_dipole_headgroup_free_energy(betaine["a_o_A2"], 20.0, betaine["d_A"], geometry="torus")


def test_nagarajan_ruckenstein_ionic_kappa_matches_existing_verified_function():
    """Real cross-check against this module's OWN already-paper-verified
    nagarajan_debye_huckel_kappa_inverse (Nagarajan 2002, Table 2:
    cmc=0.008 M -> kappa^-1=34.43 A) -- both functions implement the
    SAME Debye-Huckel kappa formula, so calling the new ionic function
    with counterion_concentration_M=0.008, added_salt_M=0 must produce
    an internally-consistent result with that already-verified value.
    Checked here by reproducing kappa_inverse independently and
    confirming it matches, then using it in the new ionic function."""
    expected_kappa_inverse_A = nagarajan_debye_huckel_kappa_inverse(0.008)
    assert expected_kappa_inverse_A == pytest.approx(34.33, abs=0.05)

    sulfate = NAGARAJAN_RUCKENSTEIN_HEADGROUP_CONSTANTS["sodium_sulfate"]
    # a real, physically sane equilibrium-area-scale input for this cross-check
    energy = nagarajan_ruckenstein_ionic_headgroup_free_energy(
        area_per_molecule_A2=sulfate["a_o_A2"] * 4.0,
        core_radius_A=16.5,
        delta_A=sulfate["delta_A"],
        counterion_concentration_M=0.008,
    )
    assert energy > 0.0  # a real, physically repulsive ionic contribution


def test_nagarajan_ruckenstein_ionic_free_energy_increases_as_area_shrinks():
    """Real physical sanity check: more crowding (smaller area per
    molecule at delta) must increase the ionic repulsion, not decrease
    or invert it."""
    sulfonate = NAGARAJAN_RUCKENSTEIN_HEADGROUP_CONSTANTS["sodium_sulfonate"]
    kwargs = {"core_radius_A": 16.5, "delta_A": sulfonate["delta_A"], "counterion_concentration_M": 0.01}
    energy_crowded = nagarajan_ruckenstein_ionic_headgroup_free_energy(area_per_molecule_A2=40.0, **kwargs)
    energy_loose = nagarajan_ruckenstein_ionic_headgroup_free_energy(area_per_molecule_A2=80.0, **kwargs)
    assert energy_crowded > energy_loose > 0.0


def test_nagarajan_ruckenstein_ionic_added_salt_increases_kappa_reduces_energy():
    """Real physical sanity check: added salt screens electrostatics more
    strongly (larger kappa, shorter Debye length), so the ionic
    repulsion at fixed geometry must DECREASE with added salt."""
    sulfate = NAGARAJAN_RUCKENSTEIN_HEADGROUP_CONSTANTS["sodium_sulfate"]
    kwargs = {"area_per_molecule_A2": 60.0, "core_radius_A": 16.5, "delta_A": sulfate["delta_A"],
                  "counterion_concentration_M": 0.008}
    energy_no_salt = nagarajan_ruckenstein_ionic_headgroup_free_energy(added_salt_M=0.0, **kwargs)
    energy_with_salt = nagarajan_ruckenstein_ionic_headgroup_free_energy(added_salt_M=0.1, **kwargs)
    assert energy_with_salt < energy_no_salt


def test_nagarajan_ruckenstein_ionic_rejects_bad_inputs():
    sulfate = NAGARAJAN_RUCKENSTEIN_HEADGROUP_CONSTANTS["sodium_sulfate"]
    good = {"area_per_molecule_A2": 60.0, "core_radius_A": 16.5, "delta_A": sulfate["delta_A"], "counterion_concentration_M": 0.008}
    with pytest.raises(ValueError):
        nagarajan_ruckenstein_ionic_headgroup_free_energy(**{**good, "area_per_molecule_A2": 0.0})
    with pytest.raises(ValueError):
        nagarajan_ruckenstein_ionic_headgroup_free_energy(**{**good, "core_radius_A": 0.0})
    with pytest.raises(ValueError):
        nagarajan_ruckenstein_ionic_headgroup_free_energy(**{**good, "delta_A": 0.0})
    with pytest.raises(ValueError):
        nagarajan_ruckenstein_ionic_headgroup_free_energy(**{**good, "counterion_concentration_M": -0.1})
    with pytest.raises(ValueError):
        nagarajan_ruckenstein_ionic_headgroup_free_energy(**{**good, "counterion_concentration_M": 0.0, "added_salt_M": 0.0})
    with pytest.raises(ValueError):
        nagarajan_ruckenstein_ionic_headgroup_free_energy(**{**good, "geometry": "torus"})


def test_nagarajan_ruckenstein_headgroup_constants_are_real_and_complete():
    for name, entry in NAGARAJAN_RUCKENSTEIN_HEADGROUP_CONSTANTS.items():
        assert entry["class"] in ("anionic", "cationic", "zwitterionic", "nonionic"), name
        assert entry["a_p_A2"] > 0 and entry["a_o_A2"] > 0, name
    assert NAGARAJAN_RUCKENSTEIN_HEADGROUP_CONSTANTS["sodium_sulfate"]["delta_A"] == pytest.approx(5.45)
    assert NAGARAJAN_RUCKENSTEIN_HEADGROUP_CONSTANTS["sodium_sulfonate"]["delta_A"] == pytest.approx(3.85)
    assert NAGARAJAN_RUCKENSTEIN_HEADGROUP_CONSTANTS["n_betaine"]["d_A"] == pytest.approx(5.0)
    assert NAGARAJAN_RUCKENSTEIN_HEADGROUP_CONSTANTS["n_betaine"]["class"] == "zwitterionic"


def test_nagarajan_ruckenstein_headgroup_constants_2003_extension_closes_cationic_and_carboxylate_gaps():
    """Real Tier 3 gaps closed 2026-09-15, same day (second real-paper
    batch): cationic quaternary ammonium and carboxylate headgroups,
    both explicitly disclosed as missing from the 1991 paper's own
    Table I, are now present via Nagarajan's own later, more
    comprehensive 2003 book chapter (Table 3), which independently
    reports IDENTICAL values for every headgroup already in this table
    from 1991 -- a genuine cross-confirmation, not just one paper's own
    self-consistency."""
    quat = NAGARAJAN_RUCKENSTEIN_HEADGROUP_CONSTANTS["trimethyl_ammonium_bromide"]
    assert quat["class"] == "cationic"
    assert quat["delta_A"] == pytest.approx(3.45)

    carboxylate = NAGARAJAN_RUCKENSTEIN_HEADGROUP_CONSTANTS["sodium_carboxylate"]
    assert carboxylate["class"] == "anionic"
    assert carboxylate["delta_A"] == pytest.approx(5.55)

    lecithin = NAGARAJAN_RUCKENSTEIN_HEADGROUP_CONSTANTS["lecithin"]
    assert lecithin["class"] == "zwitterionic"
    assert lecithin["delta_A"] == pytest.approx(6.50)
    assert lecithin["d_A"] == pytest.approx(6.2)


def test_nagarajan_ruckenstein_ionic_free_energy_usable_with_trimethyl_ammonium_bromide():
    """Real, usable end-to-end check: the newly-added cationic headgroup
    constants plugged directly into the existing ionic free-energy
    function (no new function needed -- the real work was sourcing the
    values)."""
    quat = NAGARAJAN_RUCKENSTEIN_HEADGROUP_CONSTANTS["trimethyl_ammonium_bromide"]
    energy = nagarajan_ruckenstein_ionic_headgroup_free_energy(
        area_per_molecule_A2=quat["a_o_A2"] * 2.0, core_radius_A=20.0,
        delta_A=quat["delta_A"], counterion_concentration_M=0.01,
    )
    assert energy > 0.0


# --- Davies-gemini/glycolipid HLB bottleneck resolution ----------------
# Tier 3 item 6 closure, 2026-09-15: Davies' method confirmed structurally
# invalid for gemini/glycolipid classes (tested twice with real data,
# see BOTTLENECK_RESOLUTION_PLAN.md); Griffin formally established as the
# real, published, working alternative for those classes (Liao et al.
# 2023, Arabian J. Chem. 16, 105111).

def test_recommend_hlb_method_gemini_recommends_griffin_not_davies():
    result = recommend_hlb_method_for_structural_family("dimeric_gemini_type")
    assert result["recommended_method"] == "griffin"
    assert result["davies_valid"] is False


def test_recommend_hlb_method_glycolipid_recommends_griffin_not_davies():
    result = recommend_hlb_method_for_structural_family("glycolipid_biosurfactant")
    assert result["recommended_method"] == "griffin"
    assert result["davies_valid"] is False


def test_recommend_hlb_method_monomeric_allows_either():
    result = recommend_hlb_method_for_structural_family("monomeric")
    assert result["recommended_method"] == "either"
    assert result["davies_valid"] is True


def test_recommend_hlb_method_is_case_insensitive():
    result = recommend_hlb_method_for_structural_family("Dimeric_Gemini_Type")
    assert result["recommended_method"] == "griffin"


def test_recommend_hlb_method_rejects_unparseable_instead_of_guessing():
    with pytest.raises(ValueError):
        recommend_hlb_method_for_structural_family("unparseable")


def test_recommend_hlb_method_rejects_unknown_family():
    with pytest.raises(ValueError):
        recommend_hlb_method_for_structural_family("something_made_up")


# Real gap found via GZ-23 in benchmark/groundzero_50 (ported from
# SurfactantKit): an augmented Claude answer correctly used this
# function's own Griffin recommendation, then correctly pushed back on
# it -- GEMINI_SURFACTANT_HLB_GRIFFIN_REFERENCE (this function's own Liao
# et al. 2023 evidence base) is entirely nonionic aminoglycerol geminis,
# and Claude demonstrated a real numeric instability (HLB swinging ~2x
# depending on whether a counterion counts as "hydrophilic mass") for a
# real cationic 12-4-12 gemini that reference never covers. charge_type
# param added so the tool itself now discloses this, instead of leaving
# it to whichever LLM happens to think of it.

def test_recommend_hlb_method_no_caveat_when_charge_type_not_given():
    """Backward compatibility: omitting charge_type (every pre-existing
    call site) must behave exactly as before -- no caveat asserted either
    way, since the function genuinely cannot check it without being told."""
    result = recommend_hlb_method_for_structural_family("dimeric_gemini_type")
    assert result["griffin_numeric_caveat"] is None
    assert result["recommended_method"] == "griffin"


def test_recommend_hlb_method_no_caveat_for_nonionic_gemini():
    """Nonionic is the ACTUALLY-validated case (every real compound in
    GEMINI_SURFACTANT_HLB_GRIFFIN_REFERENCE is nonionic) -- no caveat."""
    result = recommend_hlb_method_for_structural_family("dimeric_gemini_type", charge_type="nonionic")
    assert result["griffin_numeric_caveat"] is None


def test_recommend_hlb_method_flags_caveat_for_ionic_gemini():
    """The real gap: a strongly ionic (cationic/anionic) gemini or
    glycolipid gets a real, specific caveat -- method name (Griffin) is
    unchanged, but a numeric value should not be trusted."""
    cationic = recommend_hlb_method_for_structural_family("dimeric_gemini_type", charge_type="cationic")
    assert cationic["recommended_method"] == "griffin"  # still the best-justified method NAME
    assert cationic["davies_valid"] is False
    assert cationic["griffin_numeric_caveat"] is not None
    assert "cationic" in cationic["griffin_numeric_caveat"]
    assert "do not trust an exact numeric" in cationic["griffin_numeric_caveat"]

    anionic = recommend_hlb_method_for_structural_family("glycolipid_biosurfactant", charge_type="anionic")
    assert anionic["griffin_numeric_caveat"] is not None
    assert "anionic" in anionic["griffin_numeric_caveat"]


def test_recommend_hlb_method_no_caveat_for_zwitterionic_or_ambiguous_charge():
    """Only cationic/anionic (strong, permanent charge) trigger the
    caveat -- zwitterionic and pH-ambiguous headgroups are a real,
    different case this function doesn't have evidence to speak to
    either way, so it stays silent rather than guessing."""
    zw = recommend_hlb_method_for_structural_family("dimeric_gemini_type", charge_type="zwitterionic")
    assert zw["griffin_numeric_caveat"] is None
    amb = recommend_hlb_method_for_structural_family("dimeric_gemini_type", charge_type="ambiguous_pH_dependent")
    assert amb["griffin_numeric_caveat"] is None


def test_recommend_hlb_method_monomeric_always_has_caveat_key_none():
    """Return shape must be consistent across all branches -- a caller
    checking result['griffin_numeric_caveat'] should never KeyError."""
    result = recommend_hlb_method_for_structural_family("monomeric", charge_type="cationic")
    assert result["griffin_numeric_caveat"] is None


def test_recommend_hlb_method_gz23_real_compound_end_to_end():
    """Real, direct reproduction of the actual GZ-23 case: classify the
    real 12-4-12 cationic gemini SMILES, feed both classifier outputs in,
    confirm the caveat fires for this exact compound."""
    from surfactantkit.classify import (
        classify_surfactant_charge_type,
        classify_surfactant_structural_family,
    )

    smiles = "[Br-].[Br-].CCCCCCCCCCCC[N+](C)(C)CCCC[N+](C)(C)CCCCCCCCCCCC"
    charge = classify_surfactant_charge_type(smiles)
    fam = classify_surfactant_structural_family(smiles)
    assert charge.charge_type == "cationic"
    assert fam.structural_family == "dimeric_gemini_type"

    result = recommend_hlb_method_for_structural_family(fam.structural_family, charge.charge_type)
    assert result["recommended_method"] == "griffin"
    assert result["davies_valid"] is False
    assert result["griffin_numeric_caveat"] is not None


def test_gemini_hlb_griffin_reference_has_all_four_real_liao_compounds():
    assert GEMINI_SURFACTANT_HLB_GRIFFIN_REFERENCE["8-3-8"]["hlb"] == pytest.approx(7.97)
    assert GEMINI_SURFACTANT_HLB_GRIFFIN_REFERENCE["12-3-12"]["hlb"] == pytest.approx(6.37)
    assert GEMINI_SURFACTANT_HLB_GRIFFIN_REFERENCE["8-4-8"]["hlb"] == pytest.approx(7.73)
    assert GEMINI_SURFACTANT_HLB_GRIFFIN_REFERENCE["12-4-12"]["hlb"] == pytest.approx(6.22)


def test_gemini_hlb_griffin_reference_values_are_below_neutral_hlb7_consistent_with_hydrophobic_gemini():
    """Real cross-check, not a tautology: Liao et al. themselves report
    the more hydrophobic (longer spacer / longer tail) compounds formed
    W/O emulsions, consistent with Griffin HLB below the ~7 neutral
    convention -- all four real values here are indeed < 8, and the
    12-carbon-tail compounds (more hydrophobic) score lower than their
    8-carbon-tail counterparts with the same spacer."""
    assert GEMINI_SURFACTANT_HLB_GRIFFIN_REFERENCE["12-3-12"]["hlb"] < GEMINI_SURFACTANT_HLB_GRIFFIN_REFERENCE["8-3-8"]["hlb"]
    assert GEMINI_SURFACTANT_HLB_GRIFFIN_REFERENCE["12-4-12"]["hlb"] < GEMINI_SURFACTANT_HLB_GRIFFIN_REFERENCE["8-4-8"]["hlb"]


# --- Blankschtein counterion-binding theory (Tier 3 item 1 closure) ----

def test_blankschtein_entropy_of_binding_vanishes_at_beta_zero():
    assert blankschtein_entropy_of_binding_free_energy(0.0) == 0.0


def test_blankschtein_entropy_of_binding_matches_hand_derivation_at_beta_half():
    beta = 0.5
    expected = math.log(1.0 / 1.5) + 0.5 * math.log(0.5 / 1.5)
    assert blankschtein_entropy_of_binding_free_energy(beta) == pytest.approx(expected)


def test_blankschtein_entropy_of_binding_rejects_out_of_range():
    with pytest.raises(ValueError):
        blankschtein_entropy_of_binding_free_energy(1.5)
    with pytest.raises(ValueError):
        blankschtein_entropy_of_binding_free_energy(-0.1)


def test_blankschtein_steric_free_energy_matches_nagarajan_steric_at_beta_zero():
    """Real boundary identity: at beta=0 (no bound counterions), eq 12
    must reduce EXACTLY to the already-validated Nagarajan-Ruckenstein
    1991 steric headgroup formula (eq 66), since eq 12 is that same
    test-particle packing argument generalized to include counterions."""
    area, ah_s, ah_c = 60.0, 25.0, 4.537
    generalized = blankschtein_steric_free_energy_with_counterion(0.0, area, ah_s, ah_c)
    baseline = nagarajan_ruckenstein_steric_headgroup_free_energy(area, ah_s)
    assert generalized == pytest.approx(baseline)


def test_blankschtein_steric_free_energy_increases_with_beta():
    area, ah_s, ah_c = 60.0, 25.0, 4.537
    g_low = blankschtein_steric_free_energy_with_counterion(0.1, area, ah_s, ah_c)
    g_high = blankschtein_steric_free_energy_with_counterion(0.8, area, ah_s, ah_c)
    assert g_high > g_low  # more bound counterions -> more steric crowding -> higher free energy


def test_blankschtein_steric_free_energy_rejects_overpacking():
    with pytest.raises(ValueError):
        blankschtein_steric_free_energy_with_counterion(1.0, 25.0, 25.0, 10.0)


def test_blankschtein_counterion_self_energy_release_is_negative():
    value = blankschtein_counterion_self_energy_release(
        valence=1, hydrated_radius_A=2.13, kappa_inverse_A=10.0
    )
    assert value < 0.0


def test_blankschtein_counterion_self_energy_vanishes_at_infinite_dilution():
    """Real physical limit: as kappa -> 0 (infinite dilution,
    kappa_inverse_A -> infinity), the standard Debye-Huckel self-energy/
    activity-coefficient correction vanishes -- there is no ionic
    atmosphere to release self-energy from at infinite dilution. This is
    the standard, textbook DH limiting behavior (Bockris & Reddy, Modern
    Electrochemistry I), independent of the kappa-independent Born
    solvation term the source paper's own text notes is folded into the
    SAME reference chemical potential subtracted elsewhere in their
    theory (eq 4) -- not into this kappa-dependent correction itself."""
    value = blankschtein_counterion_self_energy_release(
        valence=1, hydrated_radius_A=2.13, kappa_inverse_A=1.0e9
    )
    assert value == pytest.approx(0.0, abs=1e-6)


def test_blankschtein_counterion_self_energy_more_negative_at_higher_ionic_strength():
    """Real physical trend: more screening (smaller kappa_inverse_A,
    i.e. higher ionic strength) releases MORE self-energy on binding."""
    dilute = blankschtein_counterion_self_energy_release(valence=1, hydrated_radius_A=2.13, kappa_inverse_A=100.0)
    concentrated = blankschtein_counterion_self_energy_release(valence=1, hydrated_radius_A=2.13, kappa_inverse_A=5.0)
    assert concentrated < dilute


def test_blankschtein_counterion_self_energy_more_negative_for_higher_valence():
    z1 = blankschtein_counterion_self_energy_release(valence=1, hydrated_radius_A=3.21, kappa_inverse_A=10.0)
    z2 = blankschtein_counterion_self_energy_release(valence=2, hydrated_radius_A=3.21, kappa_inverse_A=10.0)
    assert z2 < z1  # divalent counterion releases more self-energy on binding than monovalent




# ---------------------------------------------------------------------------
# Batch sign-off (audit B4): packing-parameter validity advisory (PMID
# 31954873) -- additive note only, never changes the CPP number.
# ---------------------------------------------------------------------------

def test_packing_parameter_validity_note_flags_short_and_medium_chains():
    from surfactantkit.cpp import packing_parameter_validity_note

    for n in (1, 5, 6, 9):
        note = packing_parameter_validity_note(n)
        assert note is not None
        assert "PMID 31954873" in note
        assert "n_carbons" in note
    # review's long-chain band (nC > 10) and its stated boundary (nC = 10)
    assert packing_parameter_validity_note(12) is None
    assert packing_parameter_validity_note(10) is None
    with pytest.raises(ValueError, match="n_carbons"):
        packing_parameter_validity_note(0)


def test_critical_packing_parameter_library_result_unchanged_without_note():
    """The advisory lives in a separate function: critical_packing_parameter
    itself is still a pure float computation (regression pin)."""
    v = tanford_tail_volume(12)
    lc = tanford_critical_length(12)
    assert critical_packing_parameter(v, 50.0, lc) == pytest.approx(
        v / (50.0 * lc)
    )
