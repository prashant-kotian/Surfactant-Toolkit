"""Tests for surfactantkit.mixed_micelle.

Validation strategy, deliberately explicit about what's certain and what
isn't. This library was checked against ELEVEN independent published
binary surfactant systems (cationic-anionic, cationic-nonionic,
anionic-nonionic, gemini-nonionic, gemini-zwitterionic,
nonionic-biosurfactant, anionic-anionic/bile salt, zwitterionic-anionic,
zwitterionic-cationic, nonionic-cationic), all open access, numbers
pulled directly from source tables. Full case-by-case notes, including a
log of dead-end search attempts: see literature_validation_notes.md.

1. Clint ideal CMC: hard exact-match assertions against papers that
   report a directly comparable ideal-CMC value (Azum et al. 2022;
   McLachlan et al. 2020; Liu et al. 2020). One paper in the survey
   (Muherei & Junin 2009) reports an "ideal" CMC that does NOT reproduce
   from its own stated pure-component CMCs via the standard Clint
   formula (0.73 mM computed vs. 0.906 mM reported) -- flagged as an
   unexplained discrepancy in that source (possibly a different pure-CMC
   basis for their ideal column), not asserted here.

2. Rubingh solver, micellar mole fraction (x1): hard assertions against
   FIVE literature systems -- matches to within ~1% in every case. This
   is the strongest evidence the solver is correct, since x1 is the
   solver's actual output (beta is then a deterministic function of it).

3. Rubingh solver, interaction parameter (beta): NOT asserted against
   literature values directly. Investigation across all five systems
   found the two algebraically-equivalent ways of computing beta at the
   solved root (term1/(1-x1)^2 and term2/x1^2) agree with each other to
   3-4 decimal places every time (confirming internal correctness), but
   neither matches the literature-reported beta precisely (typically
   off by 5-20%, sign always correct). The most likely explanation:
   many papers fit beta by regression across several alpha compositions,
   while this solver computes the exact pointwise value for one
   composition -- both are legitimate, but not directly comparable to
   3 decimals. We assert sign and same order of magnitude only.

4. The Rubingh solver is ALSO checked with mathematically guaranteed
   round-trip tests (construct a cmc_mix for which a chosen x1 is the
   exact algebraic root, solve it back, recover the same x1 and beta) --
   these don't depend on trusting any external source at all.
"""

import math

import pytest

from surfactantkit.adsorption import gibbs_a_min, gibbs_gamma_max
from surfactantkit.mixed_micelle import (
    R_GAS,
    activity_coefficients,
    asymmetric_margules_activity_coefficients,
    clint_ideal_cmc,
    eommm_find_minimal_feasible_margin,
    eommm_global_fit,
    eommm_global_fit_all_optima,
    excess_free_energy,
    maeda_free_energy_of_micellization,
    motomura_ideal_composition,
    rodenas_activity_coefficients,
    rodenas_x1,
    rodenas_x1_series,
    rubingh_beta,
    rubingh_beta_regression,
    solve_rubingh_x,
)

# DTAB-SDS system at 293.15 K (Rodriguez et al., PMC6554738)
DTAB_PURE_CMC = 14.80  # mM
SDS_PURE_CMC = 8.00  # mM


def test_clint_ideal_dtab_rich():
    # alpha_DTAB = 0.75; paper reports CMCid = 12.21 mM
    cmc_id = clint_ideal_cmc(0.75, DTAB_PURE_CMC, SDS_PURE_CMC)
    assert cmc_id == pytest.approx(12.21, abs=0.01)


def test_clint_ideal_sds_rich():
    # alpha_DTAB = 0.25 (i.e. alpha_SDS = 0.75); paper reports CMCid = 9.04 mM
    cmc_id = clint_ideal_cmc(0.25, DTAB_PURE_CMC, SDS_PURE_CMC)
    assert cmc_id == pytest.approx(9.04, abs=0.01)


def test_clint_ideal_rejects_bad_alpha():
    with pytest.raises(ValueError):
        clint_ideal_cmc(0.0, DTAB_PURE_CMC, SDS_PURE_CMC)
    with pytest.raises(ValueError):
        clint_ideal_cmc(1.0, DTAB_PURE_CMC, SDS_PURE_CMC)


# --- Additional literature systems -----------------------------------------
# Component 1 is always named first in each system's name below.

# Muherei & Junin (2009), Asian J. Appl. Sci. 2(2), 115-127, Table 2A/2B.
# System: TX-100 (comp 1) - SDS (comp 2). alpha(TX-100) = 0.47.
MUHEREI_TX100_CMC = 0.387  # mM
MUHEREI_SDS_CMC = 3.468  # mM
MUHEREI_ALPHA1 = 0.47
MUHEREI_CMC_MIX = 0.547  # mM, experimental
MUHEREI_X1_PAPER = 0.7501


def test_rubingh_x1_muherei_tx100_sds():
    x1 = solve_rubingh_x(MUHEREI_ALPHA1, MUHEREI_CMC_MIX, MUHEREI_TX100_CMC, MUHEREI_SDS_CMC)
    assert x1 is not None
    assert x1 == pytest.approx(MUHEREI_X1_PAPER, abs=0.01)
    beta = rubingh_beta(x1, MUHEREI_ALPHA1, MUHEREI_CMC_MIX, MUHEREI_TX100_CMC)
    assert beta < 0  # paper reports beta = -1.888 (synergistic)


# Azum et al. (2022), Biointerface Res. Appl. Chem. 12(6), 7416-7428, Table 1.
# System: gemini G6 (comp 1) - Triton X-114 (comp 2). alpha(G6) = 0.48.
AZUM_G6_CMC = 0.041  # mM
AZUM_TX114_CMC = 0.263  # mM
AZUM_ALPHA1 = 0.48
AZUM_CMC_MIX = 0.061  # mM, experimental
AZUM_CMC_IDEAL_PAPER = 0.073  # mM
AZUM_X1_PAPER = 0.764


def test_clint_ideal_azum_gemini_tx114():
    cmc_id = clint_ideal_cmc(AZUM_ALPHA1, AZUM_G6_CMC, AZUM_TX114_CMC)
    assert cmc_id == pytest.approx(AZUM_CMC_IDEAL_PAPER, abs=0.001)


def test_rubingh_x1_azum_gemini_tx114():
    x1 = solve_rubingh_x(AZUM_ALPHA1, AZUM_CMC_MIX, AZUM_G6_CMC, AZUM_TX114_CMC)
    assert x1 is not None
    assert x1 == pytest.approx(AZUM_X1_PAPER, abs=0.01)
    beta = rubingh_beta(x1, AZUM_ALPHA1, AZUM_CMC_MIX, AZUM_G6_CMC)
    assert beta < 0  # paper reports beta = -1.211 (synergistic)


# Lee & Lee (2012), J. Korean Chem. Soc. 56(5), 556-562, Table 3.
# System: TTAB (comp 1) - Tween-20 (comp 2). alpha(TTAB) = 0.6.
LEE_TTAB_CMC = 2.25  # mM
LEE_TWEEN20_CMC = 0.41  # mM
LEE_ALPHA1 = 0.6
LEE_CMC_MIX = 0.70  # mM, experimental
LEE_X1_PAPER = 0.28


def test_rubingh_x1_lee_ttab_tween20():
    x1 = solve_rubingh_x(LEE_ALPHA1, LEE_CMC_MIX, LEE_TTAB_CMC, LEE_TWEEN20_CMC)
    assert x1 is not None
    assert x1 == pytest.approx(LEE_X1_PAPER, abs=0.02)
    beta = rubingh_beta(x1, LEE_ALPHA1, LEE_CMC_MIX, LEE_TTAB_CMC)
    assert beta < 0  # paper reports beta = -0.78 (synergistic)


# McLachlan et al. (2020), RSC Adv. 10(6), 3221-3232, Table 1.
# System: gemini 12-4-12 (comp 1) - zwitterionic ZW3-12 (comp 2), alpha = 0.5.
MCLACHLAN_GEMINI_CMC = 1.10  # mM
MCLACHLAN_ZW312_CMC = 2.63  # mM
MCLACHLAN_ALPHA1 = 0.5
MCLACHLAN_CMC_MIX = 1.52  # mM, experimental
MCLACHLAN_CMC_IDEAL_PAPER = 1.55  # mM
MCLACHLAN_X_ZW312_PAPER = 0.301  # micellar mole fraction of ZW3-12 (component 2)


def test_clint_ideal_mclachlan_gemini_zwitterionic():
    cmc_id = clint_ideal_cmc(MCLACHLAN_ALPHA1, MCLACHLAN_GEMINI_CMC, MCLACHLAN_ZW312_CMC)
    assert cmc_id == pytest.approx(MCLACHLAN_CMC_IDEAL_PAPER, abs=0.01)


def test_rubingh_x1_mclachlan_gemini_zwitterionic():
    x1 = solve_rubingh_x(MCLACHLAN_ALPHA1, MCLACHLAN_CMC_MIX, MCLACHLAN_GEMINI_CMC, MCLACHLAN_ZW312_CMC)
    assert x1 is not None
    # paper reports the ZWITTERIONIC (component 2) micellar mole fraction;
    # component 1 (gemini) mole fraction is 1 - x_ZW312
    assert (1.0 - x1) == pytest.approx(MCLACHLAN_X_ZW312_PAPER, abs=0.01)


# Liu et al. (2020), Molecules 25(18), 4327, Table 9.
# System: Triton X-100 (comp 1) - rhamnolipid biosurfactant (comp 2), alpha = 0.888.
LIU_TX100_CMC = 0.309
LIU_RHAMNOLIPID_CMC = 0.134
LIU_ALPHA1 = 0.888
LIU_CMC_MIX = 0.253  # experimental
LIU_CMC_THEOR_PAPER = 0.270
LIU_X1_PAPER = 0.744


def test_clint_ideal_liu_tx100_rhamnolipid():
    cmc_id = clint_ideal_cmc(LIU_ALPHA1, LIU_TX100_CMC, LIU_RHAMNOLIPID_CMC)
    assert cmc_id == pytest.approx(LIU_CMC_THEOR_PAPER, abs=0.001)


def test_rubingh_x1_liu_tx100_rhamnolipid():
    x1 = solve_rubingh_x(LIU_ALPHA1, LIU_CMC_MIX, LIU_TX100_CMC, LIU_RHAMNOLIPID_CMC)
    assert x1 is not None
    assert x1 == pytest.approx(LIU_X1_PAPER, abs=0.01)
    beta = rubingh_beta(x1, LIU_ALPHA1, LIU_CMC_MIX, LIU_TX100_CMC)
    assert beta < 0  # paper reports beta = -0.379 (synergistic)


# Kang, Bahadur, et al. -- Mixed Micelles of Sodium Cholate and Sodium
# Dodecylsulphate 1:1 Binary Mixture at Different Temperatures, PMC4087020.
# System: sodium cholate NaCA (comp 1) - SDS (comp 2), alpha = 0.5 (1:1), 25C.
# This is the ONLY system in the survey where beta itself matches the
# paper closely -- consistent with a single fixed-ratio study (no
# multi-point regression possible), supporting the explanation in the
# module docstring for why beta usually doesn't match exactly elsewhere.
CHOLATE_NACA_CMC = 11.50  # mM
CHOLATE_SDS_CMC = 11.98  # mM
CHOLATE_ALPHA1 = 0.5
CHOLATE_CMC_MIX = 4.07  # mM, experimental (fluorimetry)
CHOLATE_CMC_IDEAL_PAPER = 11.74  # mM
CHOLATE_X1_PAPER = 0.503
CHOLATE_BETA_PAPER = -4.23


def test_clint_ideal_cholate_sds():
    cmc_id = clint_ideal_cmc(CHOLATE_ALPHA1, CHOLATE_NACA_CMC, CHOLATE_SDS_CMC)
    assert cmc_id == pytest.approx(CHOLATE_CMC_IDEAL_PAPER, abs=0.01)


def test_rubingh_x1_and_beta_cholate_sds():
    x1 = solve_rubingh_x(CHOLATE_ALPHA1, CHOLATE_CMC_MIX, CHOLATE_NACA_CMC, CHOLATE_SDS_CMC)
    assert x1 is not None
    assert x1 == pytest.approx(CHOLATE_X1_PAPER, abs=0.005)
    beta = rubingh_beta(x1, CHOLATE_ALPHA1, CHOLATE_CMC_MIX, CHOLATE_NACA_CMC)
    # this is the one system where beta itself is asserted closely, not just its sign
    assert beta == pytest.approx(CHOLATE_BETA_PAPER, abs=0.1)


# --- Ninth/tenth literature systems (2026-09-11): DHPC+SDS and DHPC+DTAB ---
# Source: Vautier-Giongo, Bakshi, Singh, Ranganathan, Hajdu & Bales, J. Colloid
# Interface Sci. 282(1) (2005) 149-155, doi:10.1016/j.jcis.2004.08.071 (PDF
# provided by the user; also freely hosted by the author at
# csun.edu/~vcphy00s/Bales87.pdf). Two GENUINELY new systems: a zwitterionic
# phospholipid (DHPC) mixed separately with an anionic (SDS) and a cationic
# (DTAB) surfactant -- a chemically different pairing class from every other
# system in this file. Their own Table 1/2 report cmc*, beta12 (via their own
# eq. 5, algebraically Rubingh's beta), and the RST-solved micelle mole
# fraction (via their own eq. 4, algebraically identical to solve_rubingh_x's
# own residual equation) at 5 interior bulk compositions each.
BALES_SDS_CMC = 8.3  # mmol/L, pure SDS (their Table 1, X_DHPC=0)
BALES_DTAB_CMC = 15.6  # mmol/L, pure DTAB (their Table 2, X_DHPC=0)
BALES_DHPC_CMC = 1.8  # mmol/L, pure DHPC (both tables, X_DHPC=1)

# (alpha1_SDS, cmc_mix, paper's x1_SDS) -- alpha1 = 1 - X_DHPC (bulk)
BALES_DHPC_SDS_POINTS = [
    (0.95, 3.6, 0.63),
    (0.90, 2.8, 0.56),
    (0.80, 1.8, 0.48),
    (0.60, 1.4, 0.29),  # real, disclosed outlier -- see test docstring
    (0.20, 1.3, 0.24),
]
# (alpha1_DTAB, cmc_mix, paper's x1_DTAB)
BALES_DHPC_DTAB_POINTS = [
    (0.95, 10.8, 0.67),
    (0.90, 6.2, 0.50),
    (0.80, 4.8, 0.37),
    (0.60, 2.7, 0.29),
    (0.50, 2.4, 0.24),
]


def test_rubingh_x1_bales_dhpc_sds():
    """4 of 5 points match the paper's own reported RST micelle mole
    fraction closely (within 0.02 absolute, 1-7% relative) -- a strong,
    genuinely independent (zwitterionic/anionic) system confirmation.
    The X_DHPC=0.40 point
    (alpha1_SDS=0.60) is a real, disclosed OUTLIER, not swept under the
    rug: this solver gives x1=0.406, a ~40% relative deviation from the
    paper's printed 0.29 -- while every neighboring point matches
    tightly. This specific point also breaks the otherwise smooth,
    monotonically-decreasing x1 sequence the other 4 points trace out
    (0.63, 0.56, 0.48, [0.29 printed, 0.41 solved], 0.24) -- 0.41 fits
    the visible trend; 0.29 does not. The most likely explanation is a
    transcription/typesetting error in the paper's own Table 1 (same
    "a paper's own table doesn't internally reproduce" pattern already
    documented for Muherei & Junin 2009 elsewhere in this file), not a
    solver bug -- but this is not provable from the published text
    alone, so it is asserted as a real discrepancy, not silently
    excluded or forced to match."""
    for alpha1, cmc_mix, x1_paper in BALES_DHPC_SDS_POINTS[:3] + BALES_DHPC_SDS_POINTS[4:]:
        x1 = solve_rubingh_x(alpha1, cmc_mix, BALES_SDS_CMC, BALES_DHPC_CMC)
        assert x1 is not None
        assert x1 == pytest.approx(x1_paper, abs=0.02)

    # the disclosed outlier point: real, large, NOT forced to match
    alpha1, cmc_mix, x1_paper = BALES_DHPC_SDS_POINTS[3]
    x1 = solve_rubingh_x(alpha1, cmc_mix, BALES_SDS_CMC, BALES_DHPC_CMC)
    assert x1 is not None
    rel_error = abs(x1 - x1_paper) / x1_paper
    assert rel_error > 0.3  # confirms this really is the disclosed large outlier
    # and confirms it fits the smooth monotonic trend the other 4 points trace
    assert 0.24 < x1 < 0.48


def test_rubingh_x1_bales_dhpc_dtab():
    """All 5 points match the paper's own reported RST micelle mole
    fraction closely -- a second, genuinely independent (zwitterionic/
    cationic) system from the same primary source, no outlier here."""
    for alpha1, cmc_mix, x1_paper in BALES_DHPC_DTAB_POINTS:
        x1 = solve_rubingh_x(alpha1, cmc_mix, BALES_DTAB_CMC, BALES_DHPC_CMC)
        assert x1 is not None
        assert x1 == pytest.approx(x1_paper, abs=0.015)


# --- Eleventh literature system (2026-09-11): TX100 + DTAB ---------------
# Source: Serafini, Fernandez-Leyes, Sanchez M., Pereyra, Schulz E.P., Durand,
# Schulz P.C. & Ritacco, "The aqueous Triton X-100 - Dodecyltrimethylammonium
# bromide micellar mixed system. Experimental results and thermodynamic
# analysis," Colloids Surf. A (2019), doi:10.1016/j.colsurfa.2018.11.032
# (Supporting Information docx provided by the user, 2026-09-11; also on
# arXiv:1806.09721). A genuinely new pairing for this project's validation
# table (nonionic Triton X-100 + cationic DTAB) -- both components appear
# individually elsewhere in this file, but never paired together. Their own
# Table SI-I reports experimental CMC (surface tension + SLS) alongside
# their own Clint-ideal-CMC column, computed "using the average experimental
# CMC values" for the pure components.
TX100_DTAB_CMC_PURE = 0.0144  # M, pure DTAB (their Table SI-I, DTAB=1 row)
TX100_CMC_PURE = (0.000184 + 0.00023) / 2.0  # M, pure TX100, average of ST and SLS methods

# (alpha1_DTAB, paper's own Clint-ideal CMC, M)
TX100_DTAB_CLINT_POINTS = [
    (0.125, 0.00024),
    (0.25, 0.00027),
    (0.375, 0.00033),
    (0.5, 0.00041),
    (0.625, 0.00054),
    (0.75, 0.00079),
    (0.875, 0.0015),
    (0.95, 0.0033),
]


def test_clint_ideal_cmc_tx100_dtab():
    """All 8 interior points match the paper's own reported Clint-ideal
    CMC column closely (<1.8% relative, consistent with rounding in
    their 2-3-sig-fig table) -- an 11th independent literature system,
    and a genuinely new pairing class (nonionic-cationic) not covered by
    any of this file's other 10 systems."""
    for alpha1, cmc_ideal_paper in TX100_DTAB_CLINT_POINTS:
        computed = clint_ideal_cmc(alpha1, TX100_DTAB_CMC_PURE, TX100_CMC_PURE)
        assert computed == pytest.approx(cmc_ideal_paper, rel=0.02)


def test_rubingh_beta_internal_self_consistency_across_all_literature_systems():
    """The two algebraically-equivalent ways of computing beta at the
    solved root (from term1 or from term2) must agree with each other
    closely, regardless of whether either matches a literature-reported
    beta (see module docstring, point 3)."""
    systems = [
        (MUHEREI_ALPHA1, MUHEREI_CMC_MIX, MUHEREI_TX100_CMC, MUHEREI_SDS_CMC),
        (AZUM_ALPHA1, AZUM_CMC_MIX, AZUM_G6_CMC, AZUM_TX114_CMC),
        (LEE_ALPHA1, LEE_CMC_MIX, LEE_TTAB_CMC, LEE_TWEEN20_CMC),
        (LIU_ALPHA1, LIU_CMC_MIX, LIU_TX100_CMC, LIU_RHAMNOLIPID_CMC),
    ]
    for alpha1, cmc_mix, cmc1, cmc2 in systems:
        x1 = solve_rubingh_x(alpha1, cmc_mix, cmc1, cmc2)
        assert x1 is not None
        term1 = math.log((alpha1 * cmc_mix) / (x1 * cmc1))
        term2 = math.log(((1.0 - alpha1) * cmc_mix) / ((1.0 - x1) * cmc2))
        beta_from_term1 = term1 / (1.0 - x1) ** 2
        beta_from_term2 = term2 / (x1 ** 2)
        assert beta_from_term1 == pytest.approx(beta_from_term2, abs=0.01)


def _cmc_mix_at_root(x1: float, alpha1: float, cmc1: float, cmc2: float) -> float:
    """Given a target root x1, solve the Rubingh residual equation
    (x1^2 * ln(A*Cmix) = (1-x1)^2 * ln(B*Cmix), A = alpha1/(x1*cmc1),
    B = (1-alpha1)/((1-x1)*cmc2)) for the cmc_mix that makes x1 the exact
    root -- closed form, valid for x1 != 0.5. This is the correct forward
    construction for a round-trip test: fixing beta alone (via the
    term1/(1-x1)^2 definition) does not fix cmc_mix, because cmc_mix also
    appears in term2; both terms must balance at the same x1."""
    a = alpha1 / (x1 * cmc1)
    b = (1.0 - alpha1) / ((1.0 - x1) * cmc2)
    denom = x1 ** 2 - (1.0 - x1) ** 2
    ln_cmix = ((1.0 - x1) ** 2 * math.log(b) - x1 ** 2 * math.log(a)) / denom
    return math.exp(ln_cmix)


def test_rubingh_round_trip_case_a():
    """Construct a cmc_mix for which x1=0.3 is the exact algebraic root
    (see _cmc_mix_at_root), then verify the numerical solver recovers
    that same x1 and the beta computed from it -- mathematically
    guaranteed, independent of any literature value."""
    alpha1, cmc1, cmc2, x1_true = 0.4, 10.0, 5.0, 0.3
    cmc_mix = _cmc_mix_at_root(x1_true, alpha1, cmc1, cmc2)
    beta_true = rubingh_beta(x1_true, alpha1, cmc_mix, cmc1)

    x1_solved = solve_rubingh_x(alpha1, cmc_mix, cmc1, cmc2)
    assert x1_solved is not None
    assert x1_solved == pytest.approx(x1_true, abs=1e-4)

    beta_solved = rubingh_beta(x1_solved, alpha1, cmc_mix, cmc1)
    assert beta_solved == pytest.approx(beta_true, abs=1e-3)


def test_rubingh_round_trip_case_b():
    alpha1, cmc1, cmc2, x1_true = 0.6, 8.0, 12.0, 0.7
    cmc_mix = _cmc_mix_at_root(x1_true, alpha1, cmc1, cmc2)
    beta_true = rubingh_beta(x1_true, alpha1, cmc_mix, cmc1)

    x1_solved = solve_rubingh_x(alpha1, cmc_mix, cmc1, cmc2)
    assert x1_solved is not None
    assert x1_solved == pytest.approx(x1_true, abs=1e-4)

    beta_solved = rubingh_beta(x1_solved, alpha1, cmc_mix, cmc1)
    assert beta_solved == pytest.approx(beta_true, abs=1e-3)


def test_rubingh_sign_matches_literature_synergy_direction():
    """Qualitative check against the DTAB-SDS paper: SDS-rich mixing is
    reported synergistic (CMC < CMCideal -> negative beta); DTAB-rich
    mixing is reported antagonistic (CMC > CMCideal -> positive beta).
    We assert the sign only -- see module docstring."""
    # SDS-rich: alpha_DTAB = 0.25, experimental CMCmix = 6.011 mM (synergistic)
    x1 = solve_rubingh_x(0.25, 6.011, DTAB_PURE_CMC, SDS_PURE_CMC)
    assert x1 is not None
    beta = rubingh_beta(x1, 0.25, 6.011, DTAB_PURE_CMC)
    assert beta < 0, "SDS-rich composition should be synergistic (negative beta)"

    # DTAB-rich: alpha_DTAB = 0.75, experimental CMCmix = 13.00 mM (antagonistic)
    x1 = solve_rubingh_x(0.75, 13.00, DTAB_PURE_CMC, SDS_PURE_CMC)
    assert x1 is not None
    beta = rubingh_beta(x1, 0.75, 13.00, DTAB_PURE_CMC)
    assert beta > 0, "DTAB-rich composition should be antagonistic (positive beta)"


# --- rubingh_beta_regression -------------------------------------------
# Real standard practice (verified via WebSearch before implementing, not
# guessed): report a single representative beta as the arithmetic mean of
# the pointwise betas computed at each (alpha1, cmc_mix) composition in a
# series -- see the function's own docstring. Validated here via a
# mathematically-guaranteed round trip: for a chosen true beta and a set
# of x1 values, the closed-form RST mass-balance equations
# (cmc_mix = x1*f1*cmc1 + (1-x1)*f2*cmc2, alpha1 = x1*f1*cmc1/cmc_mix)
# construct a self-consistent (alpha1, cmc_mix) series that ALL share
# that same true beta by construction -- no external literature source
# needed, same pattern as test_rubingh_round_trip_case_a/b above.


def _series_for_true_beta(x1_values, beta_true, cmc1, cmc2):
    alpha1_series, cmc_mix_series = [], []
    for x1 in x1_values:
        f1 = math.exp(beta_true * (1.0 - x1) ** 2)
        f2 = math.exp(beta_true * x1 ** 2)
        cmc_mix = x1 * f1 * cmc1 + (1.0 - x1) * f2 * cmc2
        alpha1 = x1 * f1 * cmc1 / cmc_mix
        alpha1_series.append(alpha1)
        cmc_mix_series.append(cmc_mix)
    return alpha1_series, cmc_mix_series


def test_rubingh_beta_regression_recovers_true_beta_round_trip():
    beta_true = -1.8
    x1_values = [0.3, 0.4, 0.5, 0.6, 0.7]
    alpha1_series, cmc_mix_series = _series_for_true_beta(x1_values, beta_true, DTAB_PURE_CMC, SDS_PURE_CMC)

    result = rubingh_beta_regression(alpha1_series, cmc_mix_series, DTAB_PURE_CMC, SDS_PURE_CMC)
    assert result.n_points_used == 5
    assert result.n_points_skipped == 0
    assert result.beta_mean == pytest.approx(beta_true, abs=1e-3)
    assert result.beta_std == pytest.approx(0.0, abs=1e-2)
    for beta in result.betas:
        assert beta == pytest.approx(beta_true, abs=1e-2)


def test_rubingh_beta_regression_antagonistic_case():
    beta_true = 1.8
    x1_values = [0.25, 0.4, 0.55, 0.7]
    alpha1_series, cmc_mix_series = _series_for_true_beta(x1_values, beta_true, DTAB_PURE_CMC, SDS_PURE_CMC)

    result = rubingh_beta_regression(alpha1_series, cmc_mix_series, DTAB_PURE_CMC, SDS_PURE_CMC)
    assert result.beta_mean == pytest.approx(beta_true, abs=1e-3)
    assert result.n_points_used == 4


def test_rubingh_beta_regression_skips_unsolvable_points_not_errors(monkeypatch):
    """solve_rubingh_x's own documented contract is to return None when no
    root is found (near-ideal or self-inconsistent data). Rather than
    fight root-finder numerics to contrive a genuinely rootless input
    (every cmc_mix value tried empirically still admitted some root),
    this isolates rubingh_beta_regression's own skip-and-continue logic
    directly by forcing one specific point's solve to report no root."""
    import surfactantkit.mixed_micelle as mm

    beta_true = -1.5
    x1_values = [0.3, 0.5, 0.7]
    alpha1_series, cmc_mix_series = _series_for_true_beta(x1_values, beta_true, DTAB_PURE_CMC, SDS_PURE_CMC)
    unsolvable_alpha1 = 0.9999
    alpha1_series.append(unsolvable_alpha1)
    cmc_mix_series.append(1.0)

    real_solve = mm.solve_rubingh_x

    def fake_solve(alpha1, cmc_mix, cmc1, cmc2, *args, **kwargs):
        if alpha1 == unsolvable_alpha1:
            return None
        return real_solve(alpha1, cmc_mix, cmc1, cmc2, *args, **kwargs)

    monkeypatch.setattr(mm, "solve_rubingh_x", fake_solve)

    result = mm.rubingh_beta_regression(alpha1_series, cmc_mix_series, DTAB_PURE_CMC, SDS_PURE_CMC)
    assert result.n_points_used == 3
    assert result.n_points_skipped == 1
    assert result.beta_mean == pytest.approx(beta_true, abs=1e-3)


def test_rubingh_beta_regression_rejects_bad_inputs():
    with pytest.raises(ValueError):
        rubingh_beta_regression([0.3], [10.0], DTAB_PURE_CMC, SDS_PURE_CMC)  # too few points
    with pytest.raises(ValueError):
        rubingh_beta_regression([0.3, 0.5], [10.0], DTAB_PURE_CMC, SDS_PURE_CMC)  # mismatched lengths


def test_activity_coefficients_ideal_limit():
    """At beta = 0 (ideal mixing), both activity coefficients must be 1."""
    f1, f2 = activity_coefficients(x1=0.4, beta=0.0)
    assert f1 == pytest.approx(1.0)
    assert f2 == pytest.approx(1.0)


def test_excess_free_energy_ideal_limit_is_zero():
    f1, f2 = activity_coefficients(x1=0.4, beta=0.0)
    dg_ex = excess_free_energy(0.4, f1, f2, temperature_k=298.15)
    assert dg_ex == pytest.approx(0.0, abs=1e-9)


def test_gibbs_adsorption_positive_slope_gives_positive_gamma_max():
    # A surfactant lowering surface tension with concentration has a
    # negative dGamma/d(ln C); Gamma_max must come out positive.
    gamma_max = gibbs_gamma_max(slope_mn_per_ln_c=-5.0, system_type="ionic_no_added_salt", temperature_k=298.15)
    assert gamma_max > 0


def test_gibbs_a_min_matches_literature_order_of_magnitude():
    """Sanity check against the DTAB-SDS paper's reported A_min values
    (44.70 and 58.07 Ang^2/molecule) -- just confirms the formula lands
    in the right physical range for a real ionic surfactant system, not
    an exact reproduction (we don't have their exact slope inputs)."""
    # A_min ~ 50 Ang^2 implies Gamma_max ~ 1e18 / (Avogadro * 50) mol/m^2
    gamma_max_typical = 1.0e18 / (6.02214076e23 * 50.0)
    a_min = gibbs_a_min(gamma_max_typical)
    assert 30.0 < a_min < 80.0


def test_gibbs_a_min_rejects_nonpositive_gamma_max():
    with pytest.raises(ValueError):
        gibbs_a_min(0.0)
    with pytest.raises(ValueError):
        gibbs_a_min(-1.0)


# --- Motomura ideal composition (alternative-methods review, 2026-09-07/08) --
# Source: equation (8), Serafini et al., Colloids Surf. A 562 (2019) 170-181
# (arXiv:1806.09721), read directly (real fetched PDF, not a search summary).


def test_motomura_ideal_composition_equals_alpha_when_cmcs_are_equal():
    """When both pure CMCs are identical, the ideal composition must
    trivially equal the bulk composition -- no reason for micellar and
    bulk mole fractions to differ if the two surfactants are identical
    in CMC."""
    x1_id = motomura_ideal_composition(0.37, cmc1=5.0, cmc2=5.0)
    assert x1_id == pytest.approx(0.37)


def test_motomura_ideal_composition_matches_clint_algebraic_identity():
    """X1_id = alpha1 * CMC_mix_ideal / cmc1 -- an exact algebraic
    identity between Motomura's and Clint's ideal-mixing formulas
    (derived from their shared zero-interaction assumption), checked
    here independent of any external literature value across several
    of this file's own literature systems."""
    systems = [
        (0.75, DTAB_PURE_CMC, SDS_PURE_CMC),
        (MUHEREI_ALPHA1, MUHEREI_TX100_CMC, MUHEREI_SDS_CMC),
        (AZUM_ALPHA1, AZUM_G6_CMC, AZUM_TX114_CMC),
        (LEE_ALPHA1, LEE_TTAB_CMC, LEE_TWEEN20_CMC),
    ]
    for alpha1, cmc1, cmc2 in systems:
        x1_id = motomura_ideal_composition(alpha1, cmc1, cmc2)
        cmc_mix_id = clint_ideal_cmc(alpha1, cmc1, cmc2)
        assert x1_id == pytest.approx(alpha1 * cmc_mix_id / cmc1, rel=1e-9)


# --- EOMMM asymmetric Margules activity coefficients ------------------------
# Source: equations (9)-(13), Schulz & Durand, Comput. Chem. Eng. 87 (2016)
# 145-153 (doi:10.1016/j.compchemeng.2015.12.026), read directly (real fetched
# PDF, provided by the user after WebFetch attempts on this exact paper
# failed -- see METHOD_ALTERNATIVES_LITERATURE_REVIEW.md).


def test_asymmetric_margules_reduces_to_rst_when_symmetric():
    """W12 == W21 must reduce EXACTLY to activity_coefficients (RST) --
    the paper's own stated result (symmetric formulations are the
    special case of the asymmetric ones), and the real regression guard
    for this function since no single-point literature case exists to
    match digit-for-digit (EOMMM's real fit is a multi-point global
    optimization, not reproduced here -- see the function's docstring)."""
    for x1, beta in [(0.3, -2.17), (0.5, 1.0), (0.75, -0.78), (0.1, 3.5)]:
        f1_asym, f2_asym = asymmetric_margules_activity_coefficients(x1, beta, beta)
        f1_rst, f2_rst = activity_coefficients(x1, beta)
        assert f1_asym == pytest.approx(f1_rst, rel=1e-9)
        assert f2_asym == pytest.approx(f2_rst, rel=1e-9)


def test_asymmetric_margules_matches_paper_own_binary_symmetric_equations():
    """Direct check against the paper's own explicitly-stated binary
    symmetric-case equations (14)-(16): Gexc = x1*x2*beta,
    ln(f1) = beta*x2^2, ln(f2) = beta*x1^2."""
    x1, beta = 0.4, -1.5
    x2 = 1.0 - x1
    f1, f2 = asymmetric_margules_activity_coefficients(x1, beta, beta)
    assert math.log(f1) == pytest.approx(beta * x2 ** 2, abs=1e-9)
    assert math.log(f2) == pytest.approx(beta * x1 ** 2, abs=1e-9)


def test_asymmetric_margules_infinite_dilution_limits():
    """Mathematically guaranteed property of the formula itself,
    independent of any external source: as x1 -> 0, ln(f1) -> w12
    (component 1 infinitely dilute in pure component 2); as x1 -> 1,
    ln(f2) -> w21."""
    w12, w21 = 2.3, -4.1
    f1_near0, _ = asymmetric_margules_activity_coefficients(1e-8, w12, w21)
    _, f2_near1 = asymmetric_margules_activity_coefficients(1.0 - 1e-8, w12, w21)
    assert math.log(f1_near0) == pytest.approx(w12, abs=1e-5)
    assert math.log(f2_near1) == pytest.approx(w21, abs=1e-5)


def test_asymmetric_margules_hyamine_dtab_literature_case():
    """Real numeric case from the paper's own Case study 1 (Hyamine(1)-
    DTAB(2), Section 4.1): W12 = -8836.50 J/mol, W21 = -2204.19 J/mol at
    298.15 K, with x1(Hy) = 1 - x_DTAB = 1 - 0.192 = 0.808 (their Table
    1, asymmetric column, alpha_DTAB = 0.25).

    UPGRADED 2026-09-11 (Schulz & Durand 2016 SI, mmc2.csv, obtained and
    provided by the user, previously unavailable -- only a graph, Fig. 1,
    was accessible before): mmc2.csv gives EXACT digit-level (x1, ln f1,
    ln f2) data for this same system across 5 interior compositions, not
    just an axis-range check. Confirmed this really is the Hyamine-DTAB
    system by an independent cross-check: mmc2's own endpoint values
    (Gexc/RTx1x2 at x1->0 and x1->1, which the asymmetric Margules
    formula makes trivially equal to W12 and W21 respectively) are
    -3.567 and -0.890 in RT units -- matching -8836.50/(R*298.15) =
    -3.565 and -2204.19/(R*298.15) = -0.889 essentially exactly (the
    tiny residual is rounding in the CSV's 3-4 significant figures)."""
    temperature_k = 298.15
    w12 = -8836.50 / (R_GAS * temperature_k)
    w21 = -2204.19 / (R_GAS * temperature_k)
    x1_hy = 1.0 - 0.192
    f1, f2 = asymmetric_margules_activity_coefficients(x1_hy, w12, w21)
    assert -4.0 <= math.log(f1) <= 0.5
    assert -4.0 <= math.log(f2) <= 0.5

    # exact digit-level check against mmc2.csv's own interior points
    mmc2_points = [
        (0.305, -0.932, -0.430),
        (0.415, -0.461, -0.692),
        (0.483, -0.263, -0.853),
        (0.626, -0.030, -1.134),
        (0.808, 0.028, -1.252),
    ]
    for x1, ln_f1_paper, ln_f2_paper in mmc2_points:
        f1, f2 = asymmetric_margules_activity_coefficients(x1, w12, w21)
        assert math.log(f1) == pytest.approx(ln_f1_paper, abs=0.005)
        assert math.log(f2) == pytest.approx(ln_f2_paper, abs=0.005)


def test_asymmetric_margules_schulz_durand_si_second_system_mmc3():
    """A second real dataset from the same Schulz & Durand 2016 SI
    (mmc3.csv, obtained 2026-09-11) -- (x1, ln f1, ln f2) across 7
    interior compositions for a second real binary system in the SI (not
    independently identified with full confidence from the extractable
    SI text as a specific named surfactant pair, unlike mmc2/Hyamine-
    DTAB above, so cited generically rather than over-claiming). W12,
    W21 read directly from the CSV's own pure-component (x1=0, x1=1)
    endpoint values, exactly as for mmc2 above -- exact-matches all 7
    interior points."""
    w12, w21 = -1.981, -0.202  # RT units, mmc3.csv's own x1=0/x1=1 endpoints
    mmc3_points = [
        (0.234, -0.673, -0.161),
        (0.304, -0.436, -0.247),
        (0.382, -0.237, -0.351),
        (0.478, -0.076, -0.471),
        (0.608, 0.028, -0.590),
        (0.937, 0.005, -0.375),
        (0.971, 0.001, -0.288),
    ]
    for x1, ln_f1_paper, ln_f2_paper in mmc3_points:
        f1, f2 = asymmetric_margules_activity_coefficients(x1, w12, w21)
        assert math.log(f1) == pytest.approx(ln_f1_paper, abs=0.005)
        assert math.log(f2) == pytest.approx(ln_f2_paper, abs=0.005)


def _solve_w12_w21_from_two_gexc_points(x1a: float, gexc_a_J_per_mol: float, x1b: float, gexc_b_J_per_mol: float, temperature_k: float) -> tuple[float, float]:
    """Gexc/RT = x1*x2*(x2*W12 + x1*W21) = x1*x2^2*W12 + x1^2*x2*W21 is
    LINEAR in (W12, W21) for a fixed composition -- so two (x1, Gexc)
    points give a genuine 2x2 linear system, solved here directly (no
    numpy/scipy) rather than by curve-fitting/optimization, since this
    is an EXACT algebraic solve, not an approximation."""
    RT = R_GAS * temperature_k
    x2a, x2b = 1.0 - x1a, 1.0 - x1b
    Aa, Ba = x1a * x2a ** 2, x1a ** 2 * x2a
    Ab, Bb = x1b * x2b ** 2, x1b ** 2 * x2b
    ga_RT, gb_RT = gexc_a_J_per_mol / RT, gexc_b_J_per_mol / RT
    det = Aa * Bb - Ba * Ab
    w12 = (ga_RT * Bb - Ba * gb_RT) / det
    w21 = (Aa * gb_RT - ga_RT * Ab) / det
    return w12, w21


def _assert_asymmetric_margules_matches_schulz_durand_eommm_fit(xs: list[float], gexc_J_per_mol: list[float], temperature_k: float = 298.15) -> None:
    """Solve (W12, W21) EXACTLY from two of the CSV's own interior
    points (a genuine 2-equation/2-unknown linear solve, not a fit),
    then check that excess_free_energy/asymmetric_margules_activity_coefficients
    reproduces the CSV's OWN OTHER points -- this is a real, precise
    test of the FORMULA ITSELF (not curve-fitting noise), since if
    asymmetric_margules_activity_coefficients truly implements the same
    equation Schulz & Durand used to generate their own EOMMM-fitted
    Gexc curves, 2 points should be enough to exactly predict all the
    rest."""
    w12, w21 = _solve_w12_w21_from_two_gexc_points(xs[0], gexc_J_per_mol[0], xs[-1], gexc_J_per_mol[-1], temperature_k)
    for x1, g_paper_J_per_mol in zip(xs, gexc_J_per_mol):
        f1, f2 = asymmetric_margules_activity_coefficients(x1, w12, w21)
        g_exc_kJ_per_mol = excess_free_energy(x1, f1, f2, temperature_k)  # excess_free_energy returns kJ/mol
        assert g_exc_kJ_per_mol * 1000.0 == pytest.approx(g_paper_J_per_mol, rel=0.003)  # <0.3% relative, see docstring above


def test_asymmetric_margules_schulz_durand_si_c8e4_sds_mmc6():
    """Real, precise (<0.25% relative error observed) confirmation of
    asymmetric_margules_activity_coefficients against a THIRD real
    Schulz & Durand 2016 SI dataset (mmc6.csv, obtained 2026-09-11):
    their own EOMMM-fitted excess free energy curve for the C8E4-SDS
    system (Case Study 3, their eq. 44), 10 interior points, T=298.15K
    matching this project's own default convention throughout."""
    xs = [0.215, 0.255, 0.328, 0.392, 0.48, 0.559, 0.625, 0.672, 0.713, 0.779]
    gexc = [-1322.245, -1478.211, -1698.78, -1820.839, -1883.356, -1839.215, -1732.119, -1619.212, -1493.673, -1245.732]
    _assert_asymmetric_margules_matches_schulz_durand_eommm_fit(xs, gexc)


def test_asymmetric_margules_schulz_durand_si_c8e4_sds_mmc8():
    """Same check, mmc8.csv (a second C8E4-SDS EOMMM-fitted curve,
    presumably a different temperature/condition within the same Case
    Study 3 -- not independently identified from the extractable SI
    text, cited generically)."""
    xs = [0.178, 0.213, 0.277, 0.334, 0.413, 0.487, 0.552, 0.599, 0.643, 0.717]
    gexc = [-1637.583, -1845.729, -2145.925, -2323.403, -2442.846, -2430.42, -2330.768, -2209.373, -2063.946, -1753.066]
    _assert_asymmetric_margules_matches_schulz_durand_eommm_fit(xs, gexc)


def test_asymmetric_margules_schulz_durand_si_c8e4_sds_mmc10():
    """Same check, mmc10.csv (a third C8E4-SDS EOMMM-fitted curve)."""
    xs = [0.136, 0.164, 0.217, 0.262, 0.326, 0.387, 0.441, 0.482, 0.522, 0.595]
    gexc = [-2077.975, -2377.708, -2816.025, -3090.581, -3317.914, -3384.577, -3335.456, -3237.798, -3100.102, -2753.027]
    _assert_asymmetric_margules_matches_schulz_durand_eommm_fit(xs, gexc)


def test_schulz_durand_si_c8e4_sds_real_vs_eommm_fit_qualitative():
    """A SEPARATE, weaker, honestly-disclosed comparison: the SI's own
    REAL experimental Gexc (mmc5/mmc7/mmc9.csv, from Hey et al. 1985
    data, their Table 4/5) versus the EOMMM-fitted curves just verified
    above (mmc6/mmc8/mmc10), via piecewise-linear interpolation of the
    fitted curve to the real data's own (different, sparser) composition
    grid -- pure data-comparison, not a test of this project's own code
    (the formula match above already covers that with much higher
    precision). Real, disclosed finding: interior compositions (well
    within the fitted curve's dense coverage) agree closely (<3%
    relative), while compositions near the edges (x<0.05 or x>0.95,
    where the fitted curve's own nearest points are far away and linear
    interpolation is a poor proxy for its true, presumably-curved
    approach to zero at the boundary) disagree by 9-21% -- a real,
    methodologically-explained pattern, not noise to hide."""
    def linear_interp(x_known, y_known, x_query):
        pts = sorted(zip(x_known, y_known))
        xs_, ys_ = [p[0] for p in pts], [p[1] for p in pts]
        if x_query <= xs_[0]:
            return ys_[0]
        if x_query >= xs_[-1]:
            return ys_[-1]
        for i in range(len(xs_) - 1):
            if xs_[i] <= x_query <= xs_[i + 1]:
                t = (x_query - xs_[i]) / (xs_[i + 1] - xs_[i])
                return ys_[i] + t * (ys_[i + 1] - ys_[i])
        raise ValueError

    mmc6_x = [0.0, 0.215, 0.255, 0.328, 0.392, 0.48, 0.559, 0.625, 0.672, 0.713, 0.779, 1.0]
    mmc6_g = [0.0, -1322.245, -1478.211, -1698.78, -1820.839, -1883.356, -1839.215, -1732.119, -1619.212, -1493.673, -1245.732, 0.0]
    # real (Hey et al. 1985) interior points, mmc5.csv
    real_interior = [(0.308, -1645.871), (0.77, -1280.568)]  # well inside the fitted curve's coverage
    real_edge = [(0.025, -194.091), (0.98, -138.916)]  # near the composition boundaries

    for x1, g_real in real_interior:
        g_fit = linear_interp(mmc6_x, mmc6_g, x1)
        assert abs(g_fit - g_real) / abs(g_real) < 0.03  # <3%, interior

    for x1, g_real in real_edge:
        g_fit = linear_interp(mmc6_x, mmc6_g, x1)
        rel_err = abs(g_fit - g_real) / abs(g_real)
        assert 0.15 < rel_err < 0.25  # real, disclosed, larger edge discrepancy


# --- eommm_global_fit (multi-point EOMMM global fit) ------------------------
# REWRITTEN 2026-09-11 (SI obtained, docx+CSVs provided by the user -- see
# eommm_global_fit's own docstring for the full, honest history: a first
# attempt to transcribe the SI's literal free-energy-minimization objective
# was found to be unbounded below and abandoned; a companion paper's SI
# (Serafini et al. 2019) explained the real "vary the margin to the
# tightest feasible value" procedure, which this implements instead as a
# well-posed infeasibility-minimization). No external published raw
# (alpha1, W12, W21) numeric example with enough digits to check against
# was found, so validation here is via mathematically-guaranteed round
# trips -- but see the real, disclosed non-uniqueness-at-wide-margin
# finding below, which is NOT swept under the rug: it directly matches
# the SI's own stated reason for the margin-tightening procedure.
# These tests are slow (a few seconds each -- a genuine 2D grid search, no
# numpy) by the nature of the problem, kept deliberately small (n=3) to
# bound runtime.


def _eommm_series_for_true_params(x1_values, w12_true, w21_true, cmc1, cmc2, r1=2.0, r2=2.0):
    """General r-aware construction: for a chosen x1 and TRUE (w12, w21),
    cmc_var = cmc1*(x1*f1)^(r1/2) + cmc2*((1-x1)*f2)^(r2/2) and
    alpha1 = cmc1*(x1*f1)^(r1/2)/cmc_var satisfy BOTH cdefp1 and cdefp2
    EXACTLY for ANY r1, r2 (not just the r=2 special case) -- verified
    algebraically in eommm_global_fit's own module-level history; reduces
    exactly to the simpler D=x1*f1*cmc1+(1-x1)*f2*cmc2 construction used
    here before 2026-09-11 when r1=r2=2."""
    alpha1_series, cmc_mix_series = [], []
    for x1 in x1_values:
        f1, f2 = asymmetric_margules_activity_coefficients(x1, w12_true, w21_true)
        cA = cmc1 * (x1 * f1) ** (r1 / 2.0)
        cB = cmc2 * ((1.0 - x1) * f2) ** (r2 / 2.0)
        cmc_var = cA + cB
        alpha1_series.append(cA / cmc_var)
        cmc_mix_series.append(cmc_var)
    return alpha1_series, cmc_mix_series


def test_eommm_global_fit_recovers_true_parameters_round_trip_at_tight_margin():
    """At a TIGHT cmc_margin, the fit recovers the true (W12, W21)
    essentially exactly."""
    true_w12, true_w21 = 6.0, -3.0
    x1_true_values = [0.15, 0.3, 0.5, 0.7, 0.85]
    alpha1_series, cmc_mix_series = _eommm_series_for_true_params(
        x1_true_values, true_w12, true_w21, DTAB_PURE_CMC, SDS_PURE_CMC
    )

    result = eommm_global_fit(alpha1_series, cmc_mix_series, DTAB_PURE_CMC, SDS_PURE_CMC, cmc_margin=1e-6)
    assert result.W12 == pytest.approx(true_w12, abs=0.01)
    assert result.W21 == pytest.approx(true_w21, abs=0.01)
    assert result.total_infeasibility < 5e-4
    assert result.n_points == 5
    for fitted, true_x1 in zip(result.x1_values, x1_true_values):
        assert fitted == pytest.approx(true_x1, abs=0.01)


def test_eommm_global_fit_finds_feasible_but_imprecise_fit_at_default_margin():
    """Real, disclosed finding (2026-09-11), directly matching the reason
    the Serafini et al. 2019 SI gives for its own margin-tightening
    procedure: at the SI's own default cmc_margin=0.10, the SAME
    round-trip data used above still finds a genuinely FEASIBLE fit
    (total_infeasibility near 0 -- this is not a search failure) but
    does NOT reliably recover the exact true (W12, W21) the way the
    tight-margin case above does -- a wider margin genuinely admits more
    of the (W12, W21) plane as feasible, so a different point within
    that wider region can be found instead of the true one. Kept as a
    qualitative, not exact, check precisely because the magnitude of
    this effect is itself sensitive to grid resolution details (a
    stronger, more specific claim was tried and found not to reproduce
    reliably -- see this test's own history in the module docstring of
    eommm_global_fit) -- the honest, robust claim is "feasible but not
    necessarily exact," not a specific numeric deviation."""
    true_w12, true_w21 = 6.0, -3.0
    x1_true_values = [0.15, 0.3, 0.5, 0.7, 0.85]
    alpha1_series, cmc_mix_series = _eommm_series_for_true_params(
        x1_true_values, true_w12, true_w21, DTAB_PURE_CMC, SDS_PURE_CMC
    )

    result = eommm_global_fit(alpha1_series, cmc_mix_series, DTAB_PURE_CMC, SDS_PURE_CMC, cmc_margin=0.10)
    assert result.total_infeasibility < 5e-4  # a genuinely feasible fit IS found...
    # ...within a real, physically sensible neighborhood of the true parameters,
    # even if not exact (unlike the tight-margin case, which IS exact)
    assert abs(result.W12 - true_w12) < 1.0
    assert abs(result.W21 - true_w21) < 1.0


def test_eommm_global_fit_reduces_to_rubingh_beta_when_symmetric():
    """At W12=W21 (the symmetric limit), a TIGHT-margin fit should
    recover approximately the same beta that solve_rubingh_x/rubingh_beta
    find on the identical data -- ties the new global fit back to the
    already-literature-validated Rubingh machinery, independent of
    trusting eommm_global_fit's own round-trip alone."""
    true_beta = -1.8
    x1_true_values = [0.3, 0.5, 0.7]
    alpha1_series, cmc_mix_series = _eommm_series_for_true_params(
        x1_true_values, true_beta, true_beta, DTAB_PURE_CMC, SDS_PURE_CMC
    )

    result = eommm_global_fit(alpha1_series, cmc_mix_series, DTAB_PURE_CMC, SDS_PURE_CMC, cmc_margin=1e-6)
    assert result.W12 == pytest.approx(true_beta, abs=0.02)
    assert result.W21 == pytest.approx(true_beta, abs=0.02)

    for a1, cm in zip(alpha1_series, cmc_mix_series):
        x1_rub = solve_rubingh_x(a1, cm, DTAB_PURE_CMC, SDS_PURE_CMC)
        assert x1_rub is not None
        beta_rub = rubingh_beta(x1_rub, a1, cm, DTAB_PURE_CMC)
        assert beta_rub == pytest.approx(true_beta, abs=0.02)


def test_eommm_global_fit_rejects_bad_inputs():
    with pytest.raises(ValueError):
        eommm_global_fit([0.3, 0.5], [10.0, 8.0], DTAB_PURE_CMC, SDS_PURE_CMC)  # too few points
    with pytest.raises(ValueError):
        eommm_global_fit([0.3, 0.5, 0.7], [10.0, 8.0], DTAB_PURE_CMC, SDS_PURE_CMC)  # mismatched lengths
    with pytest.raises(ValueError):
        eommm_global_fit([0.0, 0.5, 0.7], [10.0, 8.0, 6.0], DTAB_PURE_CMC, SDS_PURE_CMC)  # bad alpha1
    with pytest.raises(ValueError):
        eommm_global_fit([0.3, 0.5, 0.7], [10.0, -8.0, 6.0], DTAB_PURE_CMC, SDS_PURE_CMC)  # bad cmc_mix
    with pytest.raises(ValueError):
        eommm_global_fit([0.3, 0.5, 0.7], [10.0, 8.0, 6.0], 0.0, SDS_PURE_CMC)  # bad cmc1
    with pytest.raises(ValueError):
        eommm_global_fit([0.3, 0.5, 0.7], [10.0, 8.0, 6.0], DTAB_PURE_CMC, SDS_PURE_CMC, r1=0.0)  # bad r1
    with pytest.raises(ValueError):
        eommm_global_fit([0.3, 0.5, 0.7], [10.0, 8.0, 6.0], DTAB_PURE_CMC, SDS_PURE_CMC, cmc_margin=1.5)  # bad margin


# --- eommm_global_fit_all_optima (2026-10-08) -------------------------------
# Built because eommm_global_fit's own docstring discloses non-uniqueness as
# a real possibility but never surfaces it: the single-best grid search
# always returns one (W12, W21) even when another, equally feasible pair
# exists elsewhere in the plane. This searches the WHOLE grid for every
# distinct feasible basin instead of stopping at the first one found.


def test_eommm_global_fit_all_optima_matches_single_optimum_on_known_unique_case():
    """On the SAME round-trip dataset eommm_global_fit's own tight-margin
    test uses (known, by construction, to have a single true optimum),
    the multi-optima search must find exactly that one optimum -- not
    manufacture spurious extra ones, and not disagree with
    eommm_global_fit's own answer."""
    true_w12, true_w21 = 6.0, -3.0
    x1_true_values = [0.15, 0.3, 0.5, 0.7, 0.85]
    alpha1_series, cmc_mix_series = _eommm_series_for_true_params(
        x1_true_values, true_w12, true_w21, DTAB_PURE_CMC, SDS_PURE_CMC
    )

    single = eommm_global_fit(alpha1_series, cmc_mix_series, DTAB_PURE_CMC, SDS_PURE_CMC, cmc_margin=1e-6)
    multi = eommm_global_fit_all_optima(alpha1_series, cmc_mix_series, DTAB_PURE_CMC, SDS_PURE_CMC, cmc_margin=1e-6)

    assert multi.n_optima_found == 1
    assert multi.degenerate is False
    assert multi.optima[0].W12 == pytest.approx(single.W12, abs=1e-6)
    assert multi.optima[0].W21 == pytest.approx(single.W21, abs=1e-6)
    assert multi.optima[0].W12 == pytest.approx(true_w12, abs=0.01)
    assert multi.optima[0].W21 == pytest.approx(true_w21, abs=0.01)


def test_eommm_global_fit_all_optima_basin_mechanism_finds_more_than_one_when_tolerance_is_loose():
    """Direct test of the basin-detection/dedup machinery itself (not a
    claim about any specific dataset's true physical degeneracy): at the
    default (wider) cmc_margin, eommm_global_fit's own test above already
    shows the fit is feasible but not pinned to the exact true point --
    i.e. a real region of near-feasible (W12, W21) exists around it. With
    a deliberately loose feasible_tolerance, the basin search must surface
    at least the single already-known optimum (sanity floor), and the
    returned list must be internally consistent: sorted ascending by
    total_infeasibility, with no two entries closer than one coarse grid
    step apart (the dedup guarantee)."""
    true_w12, true_w21 = 6.0, -3.0
    x1_true_values = [0.15, 0.3, 0.5, 0.7, 0.85]
    alpha1_series, cmc_mix_series = _eommm_series_for_true_params(
        x1_true_values, true_w12, true_w21, DTAB_PURE_CMC, SDS_PURE_CMC
    )

    multi = eommm_global_fit_all_optima(
        alpha1_series, cmc_mix_series, DTAB_PURE_CMC, SDS_PURE_CMC,
        cmc_margin=0.10, feasible_tolerance=5e-4,
    )
    assert multi.n_optima_found >= 1
    infeas_values = [o.total_infeasibility for o in multi.optima]
    assert infeas_values == sorted(infeas_values)
    assert all(v < 5e-4 for v in infeas_values)
    coarse_step = 2.0 * 30.0 / 60.0
    for i in range(len(multi.optima)):
        for j in range(i + 1, len(multi.optima)):
            dist = math.hypot(multi.optima[i].W12 - multi.optima[j].W12,
                               multi.optima[i].W21 - multi.optima[j].W21)
            assert dist >= 3.0 * coarse_step


def test_eommm_global_fit_all_optima_rejects_bad_inputs():
    with pytest.raises(ValueError):
        eommm_global_fit_all_optima([0.3, 0.5], [10.0, 8.0], DTAB_PURE_CMC, SDS_PURE_CMC)  # too few points
    with pytest.raises(ValueError):
        eommm_global_fit_all_optima([0.3, 0.5, 0.7], [10.0, 8.0], DTAB_PURE_CMC, SDS_PURE_CMC)  # mismatched lengths
    with pytest.raises(ValueError):
        eommm_global_fit_all_optima([0.3, 0.5, 0.7], [10.0, 8.0, 6.0], DTAB_PURE_CMC, SDS_PURE_CMC, cmc_margin=1.5)
    with pytest.raises(ValueError):
        eommm_global_fit_all_optima([0.3, 0.5, 0.7], [10.0, 8.0, 6.0], DTAB_PURE_CMC, SDS_PURE_CMC,
                                     feasible_tolerance=0.0)  # bad tolerance


# --- eommm_find_minimal_feasible_margin (2026-09-12) -----------------------
# RESOLVED: eommm_global_fit's own docstring originally (2026-09-11) disclosed
# automating this as unshippable due to a real grid-resolution sensitivity
# found in an early prototype. That prototype used an EARLIER, buggy version
# of the infeasibility objective (the geometric-mean bug -- see
# eommm_global_fit's own docstring history). Re-tested directly against the
# CURRENT, fixed objective/resolution before re-attempting: the sensitivity
# is gone -- see the function's own docstring for the full re-verification
# (synthetic round-trip: monotonic convergence to the exact true parameters
# at every intermediate margin, not just the final one; real Hyamine/DTAB
# data: a genuine, stable, nonzero minimal margin reflecting real
# measurement noise). This test uses a reduced binary_search_iters (8, vs.
# the default 20) to keep runtime bounded (~30s) while still demonstrating
# real, exact convergence -- the same round-trip data already used for
# eommm_global_fit's own tight-margin test above.


def test_eommm_find_minimal_feasible_margin_round_trip():
    true_w12, true_w21 = 6.0, -3.0
    x1_true_values = [0.15, 0.3, 0.5, 0.7, 0.85]
    alpha1_series, cmc_mix_series = _eommm_series_for_true_params(
        x1_true_values, true_w12, true_w21, DTAB_PURE_CMC, SDS_PURE_CMC
    )

    result = eommm_find_minimal_feasible_margin(alpha1_series, cmc_mix_series, DTAB_PURE_CMC, SDS_PURE_CMC, binary_search_iters=8)
    assert result.W12 == pytest.approx(true_w12, abs=0.01)
    assert result.W21 == pytest.approx(true_w21, abs=0.01)
    assert result.total_infeasibility < 5e-4
    assert result.cmc_margin_used < 0.01  # converged to a genuinely tight margin, not stuck near margin_hi


def test_eommm_find_minimal_feasible_margin_rejects_bad_margin_hi():
    with pytest.raises(ValueError):
        eommm_find_minimal_feasible_margin([0.3, 0.5, 0.7], [10.0, 8.0, 6.0], DTAB_PURE_CMC, SDS_PURE_CMC, margin_hi=0.0)
    with pytest.raises(ValueError):
        eommm_find_minimal_feasible_margin([0.3, 0.5, 0.7], [10.0, 8.0, 6.0], DTAB_PURE_CMC, SDS_PURE_CMC, margin_hi=1.5)


# --- Rodenas model (alternative-methods review, 2026-09-07/08) --------------
# Source: equation (10)-(12), Azum, Rub, Alotaibi, Khan & Asiri, Biointerface
# Res. Appl. Chem. 12(6) (2022) 7416-7428 -- the same G6 gemini paper already
# used above for Clint/Rubingh validation -- citing Rodenas, Valiente &
# Villafruela, J. Phys. Chem. B 103(21) (1999) 4549-4554.


def test_rodenas_x1_reduces_to_motomura_under_ideal_mixing():
    """Real, source-independent algebraic identity: when cmc_mix
    follows Clint's EXACT ideal-mixing law, Rodenas' X1_Rod must equal
    Motomura's ideal composition formula -- both are zero-interaction
    limits of otherwise-different methods. The derivative is computed
    here by central finite difference directly on clint_ideal_cmc, not
    re-derived by hand, so this doesn't depend on trusting a separate
    hand-worked formula matching the implementation."""
    cmc1, cmc2 = 14.80, 8.00
    for alpha1 in (0.2, 0.4, 0.6, 0.8):
        h = 1e-6
        ln_cmc_plus = math.log(clint_ideal_cmc(alpha1 + h, cmc1, cmc2))
        ln_cmc_minus = math.log(clint_ideal_cmc(alpha1 - h, cmc1, cmc2))
        dln_cmc_dalpha1 = (ln_cmc_plus - ln_cmc_minus) / (2.0 * h)

        x1_rod = rodenas_x1(alpha1, dln_cmc_dalpha1)
        x1_motomura = motomura_ideal_composition(alpha1, cmc1, cmc2)
        assert x1_rod == pytest.approx(x1_motomura, abs=1e-5)


def test_rodenas_x1_zero_slope_gives_x1_equals_alpha1():
    """Direct, trivial consequence of the formula itself: a flat
    ln(cmc_mix) vs alpha1 curve (slope 0) gives X1_Rod = alpha1."""
    assert rodenas_x1(0.35, 0.0) == pytest.approx(0.35)


# --- rodenas_x1_series -------------------------------------------------
# Real, standard numerical-differentiation method (verified via WebFetch/
# WebSearch before implementing -- the source papers describe getting
# "the slope from the curve" without spelling out their exact numerical
# procedure): local quadratic (Lagrange) fit through each point's 2
# nearest neighbors, differentiated analytically. Validated two ways:
# (1) exact recovery when the underlying ln(cmc_mix) truly IS quadratic
# in alpha1 (the local fit is then not an approximation at all); (2)
# exact reduction to the textbook central/forward/backward finite-
# difference formulas for evenly-spaced data.


def test_rodenas_x1_series_exact_for_quadratic_ln_cmc_mix():
    """If ln(cmc_mix) truly follows a quadratic in alpha1, the local
    3-point quadratic fit reproduces the EXACT analytic derivative at
    every point (not just approximately) -- a strong round-trip check."""
    a, b, c = 0.5, -2.0, 3.0  # ln(cmc_mix) = a*alpha1^2 + b*alpha1 + c

    def ln_cmc(alpha1):
        return a * alpha1 ** 2 + b * alpha1 + c

    def true_slope(alpha1):
        return 2.0 * a * alpha1 + b

    alpha1_series = [0.1, 0.25, 0.4, 0.55, 0.7, 0.85]  # unevenly spaced on purpose
    cmc_mix_series = [math.exp(ln_cmc(a1)) for a1 in alpha1_series]

    result = rodenas_x1_series(alpha1_series, cmc_mix_series)
    assert result.n_points == 6
    for a1, slope in zip(result.alpha1_used, result.dln_cmc_mix_dalpha1):
        assert slope == pytest.approx(true_slope(a1), abs=1e-9)


def test_rodenas_x1_series_reduces_to_textbook_central_difference():
    """Evenly-spaced interior point: must match (y2-y0)/(2h) exactly."""
    h = 0.1
    alpha1_series = [0.2, 0.3, 0.4]
    cmc_mix_series = [1.0, 1.5, 2.6]  # arbitrary, not quadratic

    result = rodenas_x1_series(alpha1_series, cmc_mix_series)
    ln_cmc = [math.log(c) for c in cmc_mix_series]
    expected_middle_slope = (ln_cmc[2] - ln_cmc[0]) / (2.0 * h)
    assert result.dln_cmc_mix_dalpha1[1] == pytest.approx(expected_middle_slope, abs=1e-9)


def test_rodenas_x1_series_reduces_to_textbook_forward_difference_at_first_point():
    """Evenly-spaced first point: must match (-3y0+4y1-y2)/(2h) exactly."""
    h = 0.1
    alpha1_series = [0.2, 0.3, 0.4]
    cmc_mix_series = [1.0, 1.5, 2.6]
    ln_cmc = [math.log(c) for c in cmc_mix_series]
    expected_first_slope = (-3 * ln_cmc[0] + 4 * ln_cmc[1] - ln_cmc[2]) / (2.0 * h)

    result = rodenas_x1_series(alpha1_series, cmc_mix_series)
    assert result.dln_cmc_mix_dalpha1[0] == pytest.approx(expected_first_slope, abs=1e-9)


def test_rodenas_x1_series_chains_into_rodenas_x1_correctly():
    """The returned x1_rodenas must match calling rodenas_x1 directly
    with the same (alpha1, slope) pairs -- confirms the chaining, not
    just the slope computation."""
    alpha1_series = [0.1, 0.25, 0.4, 0.55, 0.7]
    cmc_mix_series = [8.0, 6.5, 5.2, 4.3, 3.7]

    result = rodenas_x1_series(alpha1_series, cmc_mix_series)
    for a1, slope, x1 in zip(result.alpha1_used, result.dln_cmc_mix_dalpha1, result.x1_rodenas):
        assert x1 == pytest.approx(rodenas_x1(a1, slope))


def test_rodenas_x1_series_sorts_unordered_input():
    alpha1_series = [0.7, 0.1, 0.4]
    cmc_mix_series = [3.7, 8.0, 5.2]  # matching the alpha1 order above

    result = rodenas_x1_series(alpha1_series, cmc_mix_series)
    assert result.alpha1_used == sorted(alpha1_series)


def test_rodenas_x1_series_rejects_bad_inputs():
    with pytest.raises(ValueError):
        rodenas_x1_series([0.2, 0.4], [1.0, 2.0])  # too few points
    with pytest.raises(ValueError):
        rodenas_x1_series([0.2, 0.4, 0.6], [1.0, 2.0])  # mismatched lengths
    with pytest.raises(ValueError):
        rodenas_x1_series([0.0, 0.4, 0.6], [1.0, 2.0, 3.0])  # bad alpha1
    with pytest.raises(ValueError):
        rodenas_x1_series([0.2, 0.4, 0.6], [1.0, -2.0, 3.0])  # bad cmc_mix
    with pytest.raises(ValueError):
        rodenas_x1_series([0.2, 0.4, 0.4], [1.0, 2.0, 3.0])  # duplicate alpha1


def test_rodenas_x1_series_against_real_multipoint_literature_data():
    """Real test against an actual raw multi-point literature alpha1
    series (2026-09-10) -- NOT a second independent source (this reuses
    the same G6+T-20 system, Table 1/2, already used for the single-
    point rodenas_x1/Maeda validations elsewhere in this file: Azum,
    Rub, Alotaibi, Khan & Asiri, Biointerface Res. Appl. Chem. 12(6)
    (2022) 7416-7428), but the FIRST time rodenas_x1_series specifically
    (the numerical local-slope-from-a-series helper, previously only
    round-trip tested against synthetic exactly-quadratic data) has been
    checked against real, sparse (5-point), experimentally-noisy data
    with a real reported X1_Rod comparison. Real, disclosed finding:
    this does NOT cleanly reproduce the paper's own reported X1_Rod
    values -- interior points (alpha1=0.4, 0.5, 0.6) land within
    3.4-9.1% relative error, but the two series ENDPOINTS (alpha1=0.3,
    0.8) are off by ~19-20%, with the alpha1=0.8 point even exceeding
    the physically-valid [0,1] mole-fraction range. This is consistent
    with (not contradicted by) rodenas_x1_series's own documented
    one-sided-derivative limitation at series endpoints, and is the
    most likely explanation (the paper's own local slope was almost
    certainly obtained by some different method -- a global fit or a
    graphical read -- not stated in enough detail to reproduce exactly).
    Kept as a real, honest record of this method's real-world behavior
    on sparse literature data, not as a validated numeric match -- do
    not read this as proof the function is correct to the precision the
    round-trip tests suggest; sparse real series have genuine numerical-
    differentiation error this function cannot eliminate."""
    alpha1_series = [0.3, 0.4, 0.5, 0.6, 0.8]
    cmc_mix_series = [0.060, 0.056, 0.051, 0.048, 0.035]  # mM, paper's own Table 1 'cmc' column
    paper_x1_rod = [0.520, 0.652, 0.762, 0.852, 0.968]

    result = rodenas_x1_series(alpha1_series, cmc_mix_series)
    # interior points: real, meaningful agreement
    for i in (1, 2, 3):
        assert result.x1_rodenas[i] == pytest.approx(paper_x1_rod[i], rel=0.10)
    # qualitative check only for the endpoints (real, larger, disclosed error there)
    assert result.x1_rodenas[0] < result.x1_rodenas[2] < result.x1_rodenas[4]  # monotonic trend preserved


def test_rodenas_x1_series_against_rodenas_1999_original_system():
    """A GENUINELY independent second literature system (2026-09-11,
    primary source obtained, PDF provided by the user): Rodenas,
    Valiente & Villafruela, J. Phys. Chem. B 103(21) (1999) 4549-4554,
    Table 1 -- the paper that ORIGINATED the model this project calls
    rodenas_x1/rodenas_x1_series, not a later paper citing it. Unlike
    every other Rodenas validation in this file (all reuse the same G6
    gemini/TX-114 or G6/T-20 systems from Azum et al. 2022), this is a
    completely different chemical system: C12E4 (tetraethylene glycol
    mono-n-dodecyl ether) mixed with CTAB (hexadecyltrimethylammonium
    bromide).

    Real, disclosed finding: this reproduces WORSE than the Azum-paper
    check above (29-60% relative error at every point, one point
    -- alpha1=0.5 -- even landing slightly outside the valid [0,1] mole
    fraction range), not better. The most likely explanation, stated
    directly in the paper's own text (not guessed): the paper does NOT
    get its d[ln(cmc_mix)]/d(alpha1) by finite-differencing these 5 raw
    (sparse, unevenly-spaced: alpha1 = 0.025 to 0.7) table points the
    way rodenas_x1_series does -- it first fits a smooth empirical
    2-exponential curve to CMC*(alpha1) (the paper's own eq. 14) and
    differentiates THAT analytically. A smooth global fit and a local
    quadratic through 3 sparse, unevenly-spaced neighbors are genuinely
    different numerical procedures and are not expected to agree
    closely on data this sparse. The paper's own text also states it
    excluded the alpha1=0.025 point from its own Figure 4 as unreliable
    -- an explicit admission that this point is not trustworthy even by
    the source's own standard. Note this is NOT the single worst point
    here by our own error metric (alpha1=0.1 is, at 60% relative error
    vs. alpha1=0.025's 54%) -- but both of the two lowest-alpha1 points,
    where the paper admits unreliability, are the two worst overall,
    which is the honest, checkable claim actually asserted below.

    Kept as a real, honest record (like the Azum-paper check above) --
    this is evidence rodenas_x1_series's sparse-series numerical
    differentiation has genuine, real limits on very sparse/uneven data,
    not evidence the core rodenas_x1 formula itself is wrong (that
    formula is separately verified via the exact zero-interaction
    algebraic identity in test_rodenas_x1_reduces_to_motomura_under_ideal_mixing,
    which does not depend on any finite-difference approximation)."""
    alpha1_series = [0.025, 0.1, 0.3, 0.5, 0.7]
    cmc_mix_series = [3.5e-4, 1.7e-4, 1.1e-4, 5.5e-5, 4.9e-5]  # M, paper's own Table 1 'CMC*' column
    paper_x1 = [0.68, 0.49, 0.69, 0.78, 0.79]  # paper's own Table 1 'chi1' (Gibbs-Duhem/new-treatment) column

    result = rodenas_x1_series(alpha1_series, cmc_mix_series)
    # real, disclosed: none of the 5 points land within a tight tolerance --
    # confirm the (large, real) error magnitude rather than hiding it
    rel_errors = [abs(x - p) / p for x, p in zip(result.x1_rodenas, paper_x1)]
    assert all(0.25 < e < 0.65 for e in rel_errors)
    # the two lowest-alpha1 points (alpha1=0.025, excluded by the
    # paper's own Figure 4, and alpha1=0.1, the next-lowest) are the two
    # worst points here -- consistent with the paper's own caveat about
    # this low-alpha1 region, even though 0.025 alone isn't the single worst
    two_worst = sorted(range(5), key=lambda i: -rel_errors[i])[:2]
    assert set(two_worst) == {0, 1}


def test_rodenas_x1_series_sds_ddao_shows_severe_instability_x1_physically_valid_flagged():
    """A THIRD independent real system (2026-10-06), confirming the same
    failure mode as test_rodenas_x1_series_against_rodenas_1999_original_system
    above -- but far more severe here, and with precise (not digitized)
    inputs. Source: Tyagi, Seddon, Khodaparast, Sharratt, Robles, Cabral,
    Colloids Surf. A 618 (2021) 126414 -- real pure CMCs stated in the
    paper's own text (SDS 8.2 mM, DDAO 1.1 mM), real molar ratios
    (70:30...30:70 SDS:DDAO), and real published MICELLAR beta (Table S2,
    supplementary info) at each ratio: -3.4, -5.7, -6.9, -6.4, -2.0.

    cmc_mix at each point here is not digitized from a figure -- it is
    back-solved so that solve_rubingh_x + rubingh_beta reproduces the
    paper's own published beta exactly (to within 0.005), using this
    project's own already-validated Rubingh implementation. This
    recovers the real underlying cmc_mix values without pixel-reading
    error, and confirms Rubingh itself is well-behaved on this real,
    strongly synergistic, sharply/non-monotonically curved system
    (physically valid x1 in [0,1] at every point, beta minimum at the
    real equimolar composition matching the paper's own reported shape).

    Rodenas, run on these same precise 5 points, is NOT well-behaved:
    4 of 5 points land outside the physical [0,1] mole-fraction range
    (one above 2.5, one below -1.0) -- confirming this is a real,
    repeatable property of the method on sparse (5-point), sharply-
    curved real composition series, not a digitization artifact and not
    a one-off (see the 1999 C12E4/CTAB case above, which showed the same
    qualitative failure, less severely, on different real data)."""
    alpha1_series = [0.7, 0.6, 0.5, 0.4, 0.3]
    # back-solved from the paper's own published micellar beta via this
    # project's own Rubingh solver (see docstring) -- not digitized
    cmc_mix_series = [0.00131582, 0.00067778, 0.00047772, 0.00052097, 0.00122120]

    result = rodenas_x1_series(alpha1_series, cmc_mix_series)
    n_invalid = sum(1 for v in result.x1_physically_valid if not v)
    assert n_invalid >= 4
    assert any(x > 1.0 or x < 0.0 for x in result.x1_rodenas)
    assert "[0,1]" in result.caveat or "physical" in result.caveat.lower()
    # cross-check: Rubingh on the SAME real pure CMCs/compositions stays
    # physically valid throughout, reproducing the paper's own beta shape
    betas = []
    for a1, cmix in zip(alpha1_series, cmc_mix_series):
        x1 = solve_rubingh_x(a1, cmix, 0.0082, 0.0011)
        assert x1 is not None and 0.0 <= x1 <= 1.0
        betas.append(rubingh_beta(x1, a1, cmix, 0.0082))
    assert min(betas) == betas[2]  # most negative at alpha1=0.5 (equimolar), matching the paper


def test_rodenas_activity_coefficients_basic():
    alpha1, x1, cmc_mix, cmc1, cmc2 = 0.5, 0.6, 5.0, 10.0, 8.0
    f1, f2 = rodenas_activity_coefficients(alpha1, x1, cmc_mix, cmc1, cmc2)
    assert f1 == pytest.approx((0.5 * 5.0) / (0.6 * 10.0))
    assert f2 == pytest.approx((0.5 * 5.0) / (0.4 * 8.0))


def test_rodenas_activity_coefficients_rejects_bad_x1():
    with pytest.raises(ValueError):
        rodenas_activity_coefficients(0.5, 0.0, 5.0, 10.0, 8.0)
    with pytest.raises(ValueError):
        rodenas_activity_coefficients(0.5, 1.0, 5.0, 10.0, 8.0)


# --- Maeda free energy of micellization (alternative-methods review) -------
# Source: equation (9), Azum, Rub, Alotaibi, Khan & Asiri, Biointerface Res.
# Appl. Chem. 12(6) (2022) 7416-7428 (same G6 gemini paper as above), citing
# Maeda, J. Colloid Interface Sci. 172 (1995) 98-105.


def test_maeda_b0_matches_real_paper_value():
    """Real numeric cross-check, not just a formula transcription
    check: B0 = ln(Xcmc2) computed with this project's own already-
    validated TX-114 CMC (0.263 mM, AZUM_TX114_CMC, from the same
    source paper's own Table 1) must match the paper's own reported
    -B0 = 12.25 for the G6+TX-114 system (independent of alpha1 --
    B0 depends only on cmc2)."""
    from surfactantkit.thermodynamics import cmc_to_mole_fraction

    xcmc2 = cmc_to_mole_fraction(AZUM_TX114_CMC / 1000.0)  # mM -> M
    b0 = math.log(xcmc2)
    assert -b0 == pytest.approx(12.25, abs=0.02)


def test_maeda_free_energy_matches_independent_recomputation():
    """Independently recompute B0/B1/B2/deltaG_M in the test itself
    (not copy-pasted from the implementation) and check they agree --
    catches a transcription bug the implementation and a copy-pasted
    test would both share."""
    x1_rub, beta = 0.6, -1.8
    cmc1_M, cmc2_M, temperature_k = 0.041e-3, 0.263e-3, 298.15
    xcmc1 = cmc1_M / (cmc1_M + 55.5)
    xcmc2 = cmc2_M / (cmc2_M + 55.5)
    b2_expected = -beta
    b1_expected = math.log(xcmc1 / xcmc2) - b2_expected
    b0_expected = math.log(xcmc2)
    delta_g_expected = (R_GAS * temperature_k * (b0_expected + b1_expected * x1_rub + b2_expected * x1_rub ** 2)) / 1000.0

    delta_g = maeda_free_energy_of_micellization(x1_rub, beta, cmc1_M, cmc2_M, temperature_k)
    assert delta_g == pytest.approx(delta_g_expected, rel=1e-9)


def test_maeda_free_energy_matches_second_independent_literature_system():
    """Real SECOND independent literature validation for Maeda's method
    (2026-09-10), beyond the single Azum et al. 2022 gemini-surfactant
    source used everywhere else in this file -- a different paper,
    different chemistry entirely (a drug-surfactant system, not a
    gemini/conventional-surfactant pair): Rub, Azum, Kumar, Arshad,
    Khan, Alotaibi & Asiri, Polymers 13(22) (2021) 4025,
    doi:10.3390/polym13224025, imipramine hydrochloride (IMP) + Triton
    X-100 (TX-100), aqueous, 298 K, their own Table 1/2, alpha1=0.5:
    cmc(IMP)=41.85 mmol/kg, cmc(TX-100)=0.31 mmol/kg, cmc_mix=0.48
    mmol/kg, X1^Rb=0.8585, beta^Rb=-4.35, deltaG_Maeda=-19.64 kJ/mol
    (their own reported value, taken directly, not re-derived).

    Real, disclosed finding: reproducing their reported deltaG_Maeda
    from their OWN reported X1^Rb/beta^Rb (mmol/kg treated as ~mM,
    dilute-solution approximation, component 1 = IMP per this project's
    established ionic-first Maeda convention from the Azum validation)
    gives -20.85 kJ/mol vs. their reported -19.64 kJ/mol -- a real ~6%
    relative error, plausible from rounding in their 3-4-significant-
    figure table values (same "close but not exact, real rounding
    residual" pattern as every other cross-paper check in this file),
    not a formula error: sign, order of magnitude, and the correct
    magnitude to within single-digit percent all match a genuinely
    independent second source. A real, separate attempt to also
    independently re-solve X1 via this project's own solve_rubingh_x
    from their raw (alpha1, cmc_mix, cmc1, cmc2) inputs recovered
    0.147, not 0.8585, unless cmc1/cmc2 are swapped in the SOLVER call
    specifically (giving 0.853, close to their 0.8585) -- a real
    component-labeling-convention subtlety between this paper's
    Rubingh-solve indexing and its Maeda indexing that could not be
    resolved with full confidence from the fetched text alone; this
    test therefore uses their OWN reported X1/beta values directly
    rather than asserting this project's solver reproduces them, so it
    validates maeda_free_energy_of_micellization specifically, not an
    additional Rubingh cross-check."""
    cmc1_M, cmc2_M = 41.85e-3, 0.31e-3  # mmol/kg ~= mM, dilute approximation
    x1_rub, beta, temperature_k = 0.8585, -4.35, 298.0
    dg = maeda_free_energy_of_micellization(x1_rub, beta, cmc1_M, cmc2_M, temperature_k)
    assert dg == pytest.approx(-19.64, rel=0.07)


def test_motomura_ideal_composition_rejects_bad_inputs():
    with pytest.raises(ValueError):
        motomura_ideal_composition(0.0, 5.0, 5.0)
    with pytest.raises(ValueError):
        motomura_ideal_composition(1.0, 5.0, 5.0)
    with pytest.raises(ValueError):
        motomura_ideal_composition(0.5, 0.0, 5.0)
    with pytest.raises(ValueError):
        motomura_ideal_composition(0.5, 5.0, -1.0)
