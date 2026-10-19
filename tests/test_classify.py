"""Tests for surfactantkit.classify.classify_surfactant_charge_type -- the
SMILES-based charge-type classifier that closes a real gap: every function
elsewhere in this library that needs a surfactant's ionic character
(gibbs_gamma_max, szyszkowski_fit_K, ...) required it as an explicit
human-supplied string, with no way to determine it from the surfactant's
own structure. Validated against real, well-known surfactant SMILES
spanning every charge class, including the deliberately ambiguous
(free-amine, no stated pH) case this function must NOT silently resolve.
"""

import pytest

from surfactantkit.classify import (
    classify_surfactant_charge_type,
    classify_surfactant_charge_type_at_ph,
    classify_surfactant_structural_family,
)


def test_sds_anionic_sulfate():
    r = classify_surfactant_charge_type("CCCCCCCCCCCCOS(=O)(=O)[O-].[Na+]")
    assert r.charge_type == "anionic"
    assert r.system_type_for_gibbs == "ionic_no_added_salt"
    assert r.confidence == "high"


def test_sds_acid_form_still_anionic():
    r = classify_surfactant_charge_type("CCCCCCCCCCCCOS(=O)(=O)O")
    assert r.charge_type == "anionic"


def test_sdbs_anionic_sulfonate():
    r = classify_surfactant_charge_type("CCCCCCCCCCCCc1ccc(cc1)S(=O)(=O)[O-].[Na+]")
    assert r.charge_type == "anionic"


def test_sodium_oleate_anionic_carboxylate_moderate_confidence():
    r = classify_surfactant_charge_type("CCCCCCCC/C=C\\CCCCCCCC(=O)[O-].[Na+]")
    assert r.charge_type == "anionic"
    assert r.confidence == "moderate"  # carboxylate is pH-sensitive near its own pKa, disclosed


def test_ctab_cationic_quaternary_ammonium():
    r = classify_surfactant_charge_type("CCCCCCCCCCCCCCCC[N+](C)(C)C.[Br-]")
    assert r.charge_type == "cationic"
    assert r.system_type_for_gibbs == "ionic_no_added_salt"
    assert r.confidence == "high"


def test_dtab_cationic():
    r = classify_surfactant_charge_type("CCCCCCCCCCCC[N+](C)(C)C.[Br-]")
    assert r.charge_type == "cationic"


def test_capb_zwitterionic():
    r = classify_surfactant_charge_type("CCCCCCCCCCCCC(=O)NCCC[N+](C)(C)CC(=O)[O-]")
    assert r.charge_type == "zwitterionic"
    assert r.system_type_for_gibbs is None  # deliberately not defaulted -- real open question


def test_pyridinium_ionic_liquid_surfactant_cationic():
    """Real gap found and fixed 2026-09-18 while building GZ-26: the
    original cationic patterns only matched sp3 [NX4+], silently missing
    a real N-alkylated aromatic ring nitrogen -- [C12mpy][Br] (Fu et al.
    2019) was misclassified as nonionic before this fix."""
    r = classify_surfactant_charge_type("CCCCCCCCCCCC[N+]1=CC=CC(=C1)C.[Br-]")
    assert r.charge_type == "cationic"
    assert r.cationic_strong_groups_found == ["aromatic_quaternary_nitrogen"]


def test_imidazolium_ionic_liquid_surfactant_cationic():
    """Same real fix, second real compound: CMIC (1-hexadecyl-3-
    methylimidazolium chloride) is already used elsewhere in this
    project's own HLD-NAC reference data (CATIONIC_QUAT_CC) -- it was
    ALSO silently misclassified as nonionic before this fix."""
    r = classify_surfactant_charge_type("CCCCCCCCCCCCCCCC[n+]1ccn(C)c1.[Cl-]")
    assert r.charge_type == "cationic"


def test_plain_pyridine_and_imidazole_not_cationic():
    """Negative control -- neither ring nitrogen carries an alkyl
    substituent here, so neither is permanently charged; the new pattern
    must not over-match plain neutral heterocycles."""
    assert classify_surfactant_charge_type("c1ccncc1").charge_type == "nonionic"
    assert classify_surfactant_charge_type("c1cnc[nH]1").charge_type == "nonionic"


def test_dhpc_phospholipid_zwitterionic_via_phosphodiester():
    """Real bug found and fixed 2026-09-15 while building Paper 3 ground-
    zero benchmark questions: DHPC (1,2-diheptanoyl-sn-glycero-3-
    phosphocholine, real zwitterionic phospholipid, PubChem CID 181610)
    has a phosphoDIester (C-O-P(=O)(O-)-O-C), which the original
    phosphate_ester pattern (monoester only) silently missed entirely --
    this molecule was misclassified as plain 'cationic' (only its
    quaternary ammonium was detected) instead of 'zwitterionic'."""
    r = classify_surfactant_charge_type(
        "CCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCC")
    assert r.charge_type == "zwitterionic"
    assert "phosphate_diester" in r.anionic_groups_found


def test_c12e8_nonionic():
    r = classify_surfactant_charge_type("CCCCCCCCCCCCOCCOCCOCCOCCOCCOCCOCCO")
    assert r.charge_type == "nonionic"
    assert r.system_type_for_gibbs == "nonionic"
    assert r.confidence == "moderate"  # absence-based claim, disclosed as such


def test_free_amine_is_ambiguous_not_guessed():
    """A free (non-quaternary) amine's charge state genuinely depends on
    solution pH, which this function is never given -- must NOT silently
    default to cationic."""
    r = classify_surfactant_charge_type("CCCCCCCCCCCCN")
    assert r.charge_type == "ambiguous_pH_dependent"
    assert r.system_type_for_gibbs is None


# --- classify_surfactant_charge_type_at_ph (Tier 1 bottleneck fix #14,
# 2026-09-15 -- see BOTTLENECK_RESOLUTION_PLAN.md). Dodecylamine
# (CCCCCCCCCCCCN, a real primary alkylamine surfactant) is the same
# fixture as test_free_amine_is_ambiguous_not_guessed above -- this
# closes the follow-on gap: even when a caller DOES supply a real pH,
# the base classifier had no way to use it.


def test_ph_conditional_dodecylamine_low_ph_is_cationic():
    # well below the primary-amine pKaH range (10.4-10.8) -> mostly protonated
    r = classify_surfactant_charge_type_at_ph("CCCCCCCCCCCCN", ph=4.0)
    assert r.charge_type_at_ph == "cationic"
    assert r.fraction_protonated_low > 0.99
    assert r.fraction_protonated_high > 0.99


def test_ph_conditional_dodecylamine_high_ph_is_nonionic():
    # well above the pKaH range -> mostly deprotonated (neutral amine)
    r = classify_surfactant_charge_type_at_ph("CCCCCCCCCCCCN", ph=13.0)
    assert r.charge_type_at_ph == "nonionic"
    assert r.fraction_protonated_low < 0.01
    assert r.fraction_protonated_high < 0.01


def test_ph_conditional_dodecylamine_near_pka_is_still_ambiguous():
    # right in the middle of the pKaH range -> genuinely mixed population
    r = classify_surfactant_charge_type_at_ph("CCCCCCCCCCCCN", ph=10.6)
    assert r.charge_type_at_ph == "still_ambiguous"
    assert 0.0 < r.fraction_protonated_low < r.fraction_protonated_high < 1.0


def test_ph_conditional_with_caller_supplied_exact_pka_is_single_point():
    r = classify_surfactant_charge_type_at_ph("CCCCCCCCCCCCN", ph=9.0, amine_pka_override=10.6)
    assert r.fraction_protonated_low == pytest.approx(r.fraction_protonated_high)
    assert r.confidence == "high"  # exact caller pKa, not a literature range
    assert r.charge_type_at_ph == "cationic"


def test_ph_conditional_delegates_for_pH_independent_charge_type():
    # SDS is anionic regardless of pH -- must delegate, not force an amine analysis
    r = classify_surfactant_charge_type_at_ph("CCCCCCCCCCCCOS(=O)(=O)[O-].[Na+]", ph=7.0)
    assert r.charge_type_at_ph == "anionic"
    assert r.amine_class is None


def test_ph_conditional_multiple_amine_classes_disclosed_not_silently_picked():
    """Real completeness gap found and fixed 2026-09-15 during a Tier 1
    review: a molecule with more than one distinct free-amine class (here,
    DMAPA, 3-dimethylaminopropylamine, CAS 109-55-7 -- a real industrial
    amine, precursor to cocamidopropyl betaine synthesis) must not have
    this function silently resolve using only the first-listed class
    without saying so."""
    dmapa = "NCCCN(C)C"  # primary amine one end, tertiary amine the other
    r = classify_surfactant_charge_type_at_ph(dmapa, ph=4.0)
    assert len(r.caveats) >= 1
    assert any("Multiple distinct amine classes" in c for c in r.caveats)
    assert r.confidence == "low"


def test_ph_conditional_rejects_bad_ph():
    with pytest.raises(ValueError):
        classify_surfactant_charge_type_at_ph("CCCCCCCCCCCCN", ph=-1.0)
    with pytest.raises(ValueError):
        classify_surfactant_charge_type_at_ph("CCCCCCCCCCCCN", ph=15.0)


def test_unparseable_smiles_reported_not_raised():
    r = classify_surfactant_charge_type("not a smiles $$$")
    assert r.charge_type == "unparseable"
    assert r.system_type_for_gibbs is None
    assert r.confidence == "none"


# --- classify_surfactant_structural_family -------------------------------
# Real, sourced structures, not invented: gemini G6/12-4-12 are the same
# literature-verified compounds already used in SurfQSPR's own
# GEMINI_CANDIDATE_STRUCTURES (Azum et al. 2022). Monorhamnolipid and
# sophorolipid SMILES pulled directly from PubChem (CID 162246, CID
# 11856871) via live lookup before writing this test, not recalled from
# training -- see classify.py's own module comment for full provenance.

GEMINI_G6 = "[Br-].[Br-].CCCCCCCCCCCCCCCC[N+](C)(C)CCCCCC[N+](C)(C)CCCCCCCCCCCCCCCC"
GEMINI_12_4_12 = "[Br-].[Br-].CCCCCCCCCCCC[N+](C)(C)CCCC[N+](C)(C)CCCCCCCCCCCC"
MONORHAMNOLIPID = "CCCCCCCC(CC(=O)O)OC(=O)CC(CCCCCCC)O[C@H]1[C@@H]([C@@H]([C@H]([C@@H](O1)C)O)O)O"
SOPHOROLIPID = "CC(CCCCCC/C=C/CCCCCCCC(=O)O)O[C@H]1[C@@H]([C@H]([C@@H]([C@H](O1)COC(=O)C)O)O)O[C@H]2[C@@H]([C@H]([C@@H]([C@H](O2)COC(=O)C)O)O)O"
GLUCOSE_ALONE = "C(C1C(C(C(C(O1)O)O)O)O)O"


def test_gemini_g6_detected_as_dimeric():
    r = classify_surfactant_structural_family(GEMINI_G6)
    assert r.structural_family == "dimeric_gemini_type"
    assert r.n_strong_ionic_headgroups == 2
    assert r.confidence == "high"
    assert r.caveats  # the bolaform-ambiguity disclosure must be present


def test_gemini_12_4_12_detected_as_dimeric():
    r = classify_surfactant_structural_family(GEMINI_12_4_12)
    assert r.structural_family == "dimeric_gemini_type"
    assert r.n_strong_ionic_headgroups == 2


def test_dtab_monomeric_not_dimeric():
    r = classify_surfactant_structural_family("CCCCCCCCCCCC[N+](C)(C)C.[Br-]")
    assert r.structural_family == "monomeric"
    assert r.n_strong_ionic_headgroups == 1


def test_ctab_monomeric_not_dimeric():
    r = classify_surfactant_structural_family("CCCCCCCCCCCCCCCC[N+](C)(C)C.[Br-]")
    assert r.structural_family == "monomeric"
    assert r.n_strong_ionic_headgroups == 1


def test_sds_monomeric_no_sugar_ring():
    r = classify_surfactant_structural_family("CCCCCCCCCCCCOS(=O)(=O)[O-].[Na+]")
    assert r.structural_family == "monomeric"
    assert r.n_sugar_rings == 0


def test_monorhamnolipid_detected_as_glycolipid_biosurfactant():
    r = classify_surfactant_structural_family(MONORHAMNOLIPID)
    assert r.structural_family == "glycolipid_biosurfactant"
    assert r.n_sugar_rings == 1
    assert r.n_long_chain_matches >= 1
    assert r.confidence == "high"


def test_sophorolipid_detected_as_glycolipid_biosurfactant():
    r = classify_surfactant_structural_family(SOPHOROLIPID)
    assert r.structural_family == "glycolipid_biosurfactant"
    assert r.n_sugar_rings == 2  # sophorose is a disaccharide -- two rings, unlike rhamnolipid's one


def test_bare_glucose_not_a_biosurfactant():
    """A sugar with no lipid tail is not a surfactant at all -- the long-
    chain requirement must correctly exclude it, not just check for a
    sugar ring alone."""
    r = classify_surfactant_structural_family(GLUCOSE_ALONE)
    assert r.structural_family == "monomeric"
    assert r.n_sugar_rings == 1
    assert r.n_long_chain_matches == 0


def test_structural_family_unparseable_smiles_reported_not_raised():
    r = classify_surfactant_structural_family("not a smiles $$$")
    assert r.structural_family == "unparseable"
    assert r.confidence == "none"
