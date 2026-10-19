"""Standard thermodynamics of micellization: Gibbs free energy, van't
Hoff enthalpy, entropy, and counterion binding degree.

Oddly absent from the rest of this library until now, despite being
central to the field -- these are the ΔG°mic/ΔH°mic/ΔS°mic triad that
essentially every experimental surfactant paper reports.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from ._optimize import solve_3x3 as _solve_3x3

R_GAS = 8.314462618  # J/(mol.K)
WATER_MOLARITY_M = 55.5  # mol/L, standard dilute-aqueous-solution approximation


def cmc_to_mole_fraction(cmc_M: float, water_molarity_M: float = WATER_MOLARITY_M) -> float:
    """Convert a CMC in molarity (mol/L) to the mole-fraction scale used
    in the standard Gibbs free energy of micellization formula:
    X_cmc = CMC / (CMC + water_molarity) ~= CMC / water_molarity for
    dilute solutions (the usual approximation, water_molarity ~= 55.5
    mol/L for water at room temperature)."""
    if cmc_M <= 0:
        raise ValueError("cmc_M must be positive")
    return cmc_M / (cmc_M + water_molarity_M)


def counterion_binding_degree(slope_below_cmc: float, slope_above_cmc: float) -> float:
    """Degree of counterion binding to the micelle, beta, from the
    slope-ratio method applied to a conductivity-vs-concentration plot:
    beta = 1 - (slope_above_cmc / slope_below_cmc).

    Both slopes should be positive (conductivity increases with
    surfactant concentration in both regimes, just at different rates).
    Note: the physical interpretation of this slope-ratio method has
    been questioned for some systems in the literature (the arithmetic
    is standard and widely used, but treat beta from this method as an
    approximate, commonly-reported value, not an unimpeachable one --
    see literature_validation_notes.md).
    """
    if slope_below_cmc <= 0 or slope_above_cmc <= 0:
        raise ValueError("both slopes must be positive")
    if slope_above_cmc >= slope_below_cmc:
        raise ValueError(
            "slope_above_cmc should be less than slope_below_cmc (conductivity "
            "rises more slowly above the CMC as counterions associate with "
            "micelles) -- check which slope is which"
        )
    return 1.0 - (slope_above_cmc / slope_below_cmc)


# Tier 2 bottleneck-resolution addition, 2026-09-15 (see benchmark/
# audit findings, 2026-09 item 1): a real
# reference table of counterion binding degree (beta) / dissociation
# degree (alpha = 1-beta) for the most-studied ionic surfactants,
# closing this bottleneck for compounds where a real measured value
# already exists in the literature. Different entries carry genuinely
# different confidence levels, disclosed per entry rather than presented
# as uniformly precise:
#
# SDS: alpha=0.272 +/- 0.017 (beta=0.728), Bales, Messina, Vidal & Peric,
# J. Phys. Chem. B 105 (2001) 6798-6804 -- the HIGHEST-confidence entry
# here: a real, novel EPR method, independently cross-validated in the
# paper itself against multiple other literature techniques, and the
# exact PDF was read in full by this project (provided by the user
# 2026-09-11, already used elsewhere -- see gibbs_free_energy_
# micellization's own docstring).
#
# DTAB and CTAB: found via real, cross-corroborated web search results
# this session (2026-09-15) -- each value independently agreed with a
# SECOND, separately-cited literature value in the same search result
# (DTAB: 0.28, cross-checked against a cited 0.29; CTAB: 0.26, cross-
# checked against a cited 0.27), both at 298.15 K via the same slope-
# ratio conductivity method this module's own counterion_binding_degree
# implements. The primary papers themselves were not independently
# fetched this session (search-result-level corroboration only, NOT the
# same standard as SDS's directly-read primary PDF) -- disclosed
# explicitly via the per-entry alpha_confidence field, not silently
# presented as equally certain.
COUNTERION_BINDING_DEGREE_REFERENCE = {
    "SDS": {"alpha": 0.272, "alpha_uncertainty": 0.017, "beta": 0.728,
            "source": "Bales, Messina, Vidal & Peric, J. Phys. Chem. B 105 (2001) 6798-6804 (EPR method, primary PDF read in full)",
            "alpha_confidence": "high"},
    "DTAB": {"alpha": 0.28, "alpha_uncertainty": None, "beta": 0.72,
             "source": "web search 2026-09-15, cross-corroborated against a second cited literature value (0.29); primary paper not independently fetched",
             "alpha_confidence": "moderate"},
    "CTAB": {"alpha": 0.26, "alpha_uncertainty": None, "beta": 0.74,
             "source": "web search 2026-09-15, cross-corroborated against a second cited literature value (0.27); primary paper not independently fetched",
             "alpha_confidence": "moderate"},
    "AOT": {"alpha": 0.61, "alpha_uncertainty": None, "beta": 0.39,
            "source": "Thapa, Ray, Dey, Sultana, Aswal & Ismail, RSC Adv. 5 (2015) 45956-45964, Table 1 "
                      "(primary PDF read in full 2026-09-15) -- Corrin-Harkins plot, NaCl medium, below the "
                      "critical salt concentration c* where AOT's own real counterion-binding transition occurs; "
                      "HIGH confidence, real primary-source measurement, same underlying method (ionic-slope-based "
                      "binding constant) as SDS's own entry above",
            "alpha_confidence": "high"},
}


def gibbs_free_energy_micellization(cmc_mole_fraction: float, temperature_K: float, counterion_factor: float = 1.0) -> float:
    """Standard Gibbs free energy of micellization, kJ/mol:
    deltaG_mic = counterion_factor * R * T * ln(X_cmc).

    counterion_factor: explicit, not guessed. Use 1.0 for a nonionic
    surfactant. For an ionic surfactant with counterion binding degree
    beta (see counterion_binding_degree), use (2 - beta) -- equivalently,
    (1 + alpha) if what's in hand is the counterion DISSOCIATION/
    ionization degree alpha (= 1 - beta) instead, e.g. as directly
    reported by Bales, Messina, Vidal & Peric, J. Phys. Chem. B 105
    (2001) 6798-6804 (PDF provided by the user, 2026-09-11), which gives
    a real, precisely-measured alpha = 0.272 +/- 0.017 for SDS via a
    novel EPR method, cross-validated there against multiple independent
    literature techniques -- i.e. counterion_factor = 1.272 for SDS
    specifically, a real sourced value rather than an arbitrary
    illustrative one. The (2-beta) = (1+alpha) identity (alpha=1-beta)
    is checked via its own physical limits, not just algebra: fully
    dissociated (alpha=1, beta=0) behaves like two independent particles
    losing translational entropy together, counterion_factor -> 2 (both
    forms agree); fully bound (alpha=0, beta=1, counterion moves as one
    package with the surfactant) behaves like a single nonionic-style
    particle, counterion_factor -> 1 (both forms agree too) -- see
    tests/test_thermodynamics.py for the worked check. This is NOT
    auto-detected from the inputs, because silently assuming nonionic
    behavior for an ionic surfactant (or vice versa) is exactly the kind
    of hidden-assumption error this library exists to prevent.
    """
    if not (0.0 < cmc_mole_fraction < 1.0):
        raise ValueError("cmc_mole_fraction must be strictly between 0 and 1")
    if temperature_K <= 0:
        raise ValueError("temperature_K must be positive")
    return (counterion_factor * R_GAS * temperature_K * math.log(cmc_mole_fraction)) / 1000.0


def mass_action_free_energy_of_micellization(cmc_mole_fraction: float, aggregation_number: float, temperature_K: float) -> float:
    """Standard Gibbs free energy of micellization, kJ/mol, via the
    MASS-ACTION-LAW model (single equilibrium n*S <-> M_n, K =
    [M_n]/[S]^n) -- a real, sourced ALTERNATIVE to
    gibbs_free_energy_micellization's PSEUDO-PHASE-SEPARATION
    convention, for NONIONIC surfactants only (see below):

    deltaG_mic = R*T * [(1 + 1/n)*ln(X_cmc) - (1/n)*ln(n)]

    Verified via an open-access source (Boneva-Aroca et al., RSC
    Advances 13 (2023) 9387, "Refined definition of the critical
    micelle concentration and application to alkyl maltosides used in
    membrane protein research," their eq. 21, quoted directly, not
    guessed or reconstructed), AND directly confirmed 2026-09-11 against
    the primary source itself (Rusanov, Langmuir 30(48) (2014) 14443-14451,
    PDF provided by the user, previously flagged paywalled): Rusanov's
    own eq. (18), cm = 1/Kj = K^(1-n), is the same single-equilibrium
    mass-action-law CMC definition this formula is built on -- his
    abstract states it directly ("CMC = K^(1-n)"). This is the SAME
    well-known classical result (also widely cited elsewhere, e.g.
    Moroi's "Micelles" monograph and Rosen's "Surfactants and
    Interfacial Phenomena"), now confirmed against both an open-access
    restatement and the primary source itself.

    The pseudo-phase-separation model this project otherwise uses
    (gibbs_free_energy_micellization) is the n -> infinity LIMIT of
    this formula: as n grows, (1+1/n) -> 1 and (1/n)*ln(n) -> 0, so
    deltaG_mic -> R*T*ln(X_cmc) exactly (verified in
    tests/test_thermodynamics.py -- confirms this is a real
    generalization, not a competing/inconsistent formula). For small-
    to-moderate aggregation numbers (n below roughly 50-100, per this
    project's own literature review), the correction terms are NOT
    negligible and the pseudo-phase result can be a real, quantifiable
    over- or under-estimate -- this function makes that correction
    available directly, matching the "let a researcher choose the
    method suited to their system" goal.

    NONIONIC SURFACTANTS ONLY: this is the simple single-equilibrium
    mass-action result with no counterion association step. Extending
    to ionic surfactants requires a materially different mass-action
    treatment (an explicit counterion-binding equilibrium alongside
    monomer aggregation) that was not verified this pass -- do not
    apply this function to an ionic surfactant's counterion_factor
    convention; use gibbs_free_energy_micellization with the
    appropriate counterion_factor for that case instead.
    """
    if not (0.0 < cmc_mole_fraction < 1.0):
        raise ValueError("cmc_mole_fraction must be strictly between 0 and 1")
    if aggregation_number <= 1:
        raise ValueError("aggregation_number must be greater than 1")
    if temperature_K <= 0:
        raise ValueError("temperature_K must be positive")
    n = aggregation_number
    g_per_RT = (1.0 + 1.0 / n) * math.log(cmc_mole_fraction) - (1.0 / n) * math.log(n)
    return (R_GAS * temperature_K * g_per_RT) / 1000.0


def critical_micellization_degree(aggregation_number: float, definition: str = "d2alpha_dc2") -> float:
    """Critical micellization degree (CMD), alpha_m: the dimensionless
    fraction n*K*c^(n-1) reaches at the CMC under the single-component
    mass-action-law model (see mass_action_free_energy_of_micellization
    for the same model's free-energy side), n(Kc)^(n-1) = alpha/(1-alpha)^n
    (Rusanov's eq. 17) -- a real, NEW capability (2026-09-11, added on
    request), not a validation of something already implemented.

    Two different formulas are given by the primary source depending on
    which curvature condition defines "the CMC" (both are legitimate;
    they simply answer "sharpest change in WHAT, as a function of
    WHAT" differently):

    definition="d2alpha_dc2" (default): CMC defined by d^2(alpha)/dc^2 = 0
    (sharpest change in micellization degree vs. CONCENTRATION):
        alpha_m = (sqrt(n/2) - 1) / (n - 1)

    definition="d2alpha_dlnc2": CMC defined by d^2(alpha)/d(ln c)^2 = 0
    (sharpest change vs. LOG concentration -- the more common
    experimental convention, e.g. the inflection point of a
    conductivity-vs-log(concentration) plot):
        alpha_m = 1 / (sqrt(n) + 1)

    Source: Rusanov, "The Mass-Action-Law Theory of Micellization
    Revisited," Langmuir 30(48) (2014) 14443-14451 (PDF provided by the
    user 2026-09-11), eqs. (15) and (16) respectively. Both formulas
    verified exactly against the paper's own worked numeric example
    (n=100): eq. 15 gives alpha_m=0.061, eq. 16 gives alpha_m=0.091 --
    both reproduced here to 3 decimal places (see
    tests/test_thermodynamics.py), not just approximately.

    Real, checkable consistency property (verified in the test suite,
    not just asserted): both formulas satisfy alpha_m -> 0 as n -> infinity
    (a sharper, more "abrupt" transition for larger aggregates, the same
    qualitative behavior the pseudo-phase-separation limit assumes),
    and definition="d2alpha_dlnc2" always gives a strictly larger
    alpha_m than definition="d2alpha_dc2" for the same n (verified
    algebraically for n>1, not just at the one worked example).

    aggregation_number must be > 1 (n=1 is not a micelle)."""
    if aggregation_number <= 1:
        raise ValueError("aggregation_number must be greater than 1")
    n = aggregation_number
    if definition == "d2alpha_dc2":
        alpha_m = (math.sqrt(n / 2.0) - 1.0) / (n - 1.0)
        if alpha_m < 0:
            raise ValueError(
                f"definition='d2alpha_dc2' gives a negative (unphysical) alpha_m for n={n} "
                "(this formula requires n >= 2 -- use definition='d2alpha_dlnc2' for small n instead)"
            )
        return alpha_m
    if definition == "d2alpha_dlnc2":
        return 1.0 / (math.sqrt(n) + 1.0)
    raise ValueError('definition must be "d2alpha_dc2" or "d2alpha_dlnc2"')


def vant_hoff_enthalpy(cmc1_mole_fraction: float, temperature1_K: float, cmc2_mole_fraction: float, temperature2_K: float, counterion_factor: float = 1.0) -> float:
    """Van't Hoff enthalpy of micellization, kJ/mol, from CMC (mole
    fraction scale) measured at two temperatures:
    deltaH_mic = counterion_factor * R * (ln(X_cmc,2) - ln(X_cmc,1)) / (1/T2 - 1/T1).

    Real sign bug found and fixed 2026-09-04: this function previously
    divided by (1/T1 - 1/T2) instead of (1/T2 - 1/T1), giving the exact
    opposite sign of the true van't Hoff enthalpy. Caught via two
    independent, convergent checks: (1) an unaugmented Claude Opus 4.8
    pilot answer for a real SurfBench question matched textbook Gibbs-
    Helmholtz-derivation and Le Chatelier physical reasoning exactly
    (CMC decreasing with rising T => endothermic, positive dH) while
    this function's old output had the opposite sign; (2) re-examining
    this project's own existing literature_validation_notes.md entry for
    this exact function against Fu et al. 2019 (RSC Advances) showed the
    old two-point estimates (+6.446, +0.462 kJ/mol) had the opposite sign
    from the paper's own local-derivative value (-4.472 kJ/mol) -- at the
    time misread as pure two-point-vs-polynomial magnitude noise, but the
    sign mismatch was real and is this same bug. Full re-derivation from
    the Gibbs-Helmholtz equation (d(dG/T)/dT = -dH/T^2, with
    dG = R*T*ln(X_cmc)) gives d(ln X_cmc) = (dH/R) d(1/T), i.e. the
    (1/T2 - 1/T1) denominator used now, not the old (1/T1 - 1/T2).

    counterion_factor: explicit, not guessed -- same convention as
    gibbs_free_energy_micellization's own counterion_factor (1.0 for
    nonionic, (2-beta)=(1+alpha) for ionic). Added 2026-09-23, real gap
    found via GZ-26 in benchmark/groundzero_50: this function
    previously had NO counterion_factor parameter at all (unlike its
    siblings gibbs_free_energy_micellization and vant_hoff_multi_point_fit),
    silently always returning the nonionic-assumption value even for an
    ionic surfactant. An augmented Claude Opus 4.8 answer for a real ionic
    ([C12mpy][Br]) two-point van't Hoff question caught this: it used
    gibbs_free_energy_micellization(counterion_factor=1.732) for deltaG
    (consistent with this project's own literature_validation_notes.md
    entry for the same compound), then manually multiplied this
    function's raw output by the same 1.732 for deltaH, since applying
    the SAME counterion_factor to both is not optional but required by
    Gibbs-Helmholtz self-consistency: if deltaG/T = counterion_factor * R
    * ln(X_cmc), then deltaH = -T^2 * d(deltaG/T)/dT = counterion_factor
    * R * d(ln X_cmc)/d(1/T) -- the identical factor falls straight out
    of differentiating the SAME deltaG formula this library already uses.
    Verified: with counterion_factor=1.0 (the default), every existing
    call site and test is completely unaffected -- this is purely
    additive, matching the sibling functions' pattern exactly.

    Assumes deltaH is constant over the T1-T2 interval -- this
    assumption weakens for wide temperature ranges or when the
    aggregation number itself varies significantly with temperature;
    treat results from a wide T range with appropriate caution.
    """
    if temperature1_K <= 0 or temperature2_K <= 0:
        raise ValueError("temperatures must be positive")
    if temperature1_K == temperature2_K:
        raise ValueError("temperature1_K and temperature2_K must differ")
    if not (0.0 < cmc1_mole_fraction < 1.0) or not (0.0 < cmc2_mole_fraction < 1.0):
        raise ValueError("CMC mole fractions must be strictly between 0 and 1")
    delta_ln_cmc = math.log(cmc2_mole_fraction) - math.log(cmc1_mole_fraction)
    delta_inv_T = (1.0 / temperature2_K) - (1.0 / temperature1_K)
    return (counterion_factor * R_GAS * delta_ln_cmc / delta_inv_T) / 1000.0


@dataclass
class VantHoffMultiPointFitResult:
    delta_H_kJ_per_mol: float  # at reference_temperature_K
    delta_S_J_per_mol_K: float  # at reference_temperature_K
    delta_Cp_J_per_mol_K: float  # assumed constant over the fitted range
    reference_temperature_K: float
    r_squared: float
    n_points: int
    delta_G_fitted_kJ_per_mol: list[float]  # model-predicted deltaG at each input T, same order as input


def vant_hoff_multi_point_fit(
    cmc_mole_fractions: list[float],
    temperatures_K: list[float],
    counterion_factor: float = 1.0,
    reference_temperature_K: float | None = None,
) -> VantHoffMultiPointFitResult:
    """Full multi-point (3+ temperatures, 5+ recommended) Gibbs-Helmholtz
    fit for deltaH, deltaS, AND deltaCp simultaneously -- the real fix
    this project's own benchmark/METHOD_ALTERNATIVES_LITERATURE_REVIEW.md
    already recommended for vant_hoff_enthalpy's strictly-2-point
    limitation, per Kantonen, Henriksen & Gilson, BBA Gen. Subj. 1862
    (2018) 692-704 (doi:10.1016/j.bbagen.2017.11.020, their Eq. 4):

        deltaG(Tj) = deltaH(Tr) - Tj*deltaS(Tr)
                     + deltaCp(Tr) * [(Tj-Tr) - Tj*ln(Tj/Tr)]

    where deltaG(Tj) = counterion_factor * R * Tj * ln(X_cmc(Tj)) (same
    formula and counterion_factor convention as
    gibbs_free_energy_micellization), and Tr is a reference temperature.

    The source paper describes fitting deltaH(Ti)/deltaS(Ti)/deltaCp(Ti)
    by "nonlinear optimization" -- but for a FIXED reference temperature
    Tr, this equation is actually LINEAR in the three unknowns
    (deltaH, deltaS, deltaCp): the Tj-dependence lives entirely in the
    three regressor coefficients [1, -Tj, (Tj-Tr)-Tj*ln(Tj/Tr)], not in
    the parameters themselves. This implementation exploits that and
    solves the exact ordinary-least-squares closed form (via Gaussian
    elimination on the 3x3 normal-equations system) instead of iterative
    nonlinear optimization -- mathematically equivalent for a fixed Tr,
    verified here via exact round-trip recovery of synthetic data (see
    tests/test_thermodynamics.py), and avoids adding a nonlinear-solver
    dependency to a project that is zero-external-dependency by design.

    reference_temperature_K: defaults to the mean of temperatures_K,
    which keeps the fit well-conditioned (reduces the numerical
    correlation between the fitted deltaH and deltaCp) -- deltaH and
    deltaS come out "at" this reference temperature (a real, unavoidable
    consequence of deltaCp being nonzero: enthalpy and entropy genuinely
    differ at different reference temperatures, unlike deltaCp itself or
    the fitted deltaG(T) curve, both of which are reference-independent).
    Override to get deltaH/deltaS reported at a specific literature-
    comparable temperature (e.g. 298.15 K) instead.

    Needs at least 3 points (the fit has 3 free parameters) with
    genuinely different temperatures; 5+ is recommended for deltaCp to
    be reliably resolved from real (noisy) data, per the source paper's
    own experimental design (they use 5 temperature points).
    """
    if len(cmc_mole_fractions) != len(temperatures_K):
        raise ValueError("cmc_mole_fractions and temperatures_K must be the same length")
    if len(cmc_mole_fractions) < 3:
        raise ValueError("need at least 3 temperature points to fit deltaH, deltaS, and deltaCp simultaneously")
    if any(not (0.0 < x < 1.0) for x in cmc_mole_fractions):
        raise ValueError("all CMC mole fractions must be strictly between 0 and 1")
    if any(t <= 0 for t in temperatures_K):
        raise ValueError("all temperatures must be positive")

    n = len(temperatures_K)
    tr = reference_temperature_K if reference_temperature_K is not None else sum(temperatures_K) / n

    delta_g = [
        (counterion_factor * R_GAS * t * math.log(x)) / 1000.0  # kJ/mol, same as gibbs_free_energy_micellization
        for x, t in zip(cmc_mole_fractions, temperatures_K)
    ]

    # Regressors per point: r1=1 (coeff deltaH), r2=-Tj (coeff deltaS),
    # r3=(Tj-Tr)-Tj*ln(Tj/Tr) (coeff deltaCp); all in kJ-consistent units
    # (deltaS, deltaCp come out in kJ/(mol.K) here, converted to J after).
    r3 = [(t - tr) - t * math.log(t / tr) for t in temperatures_K]

    m00 = float(n)
    m01 = -sum(temperatures_K)
    m02 = sum(r3)
    m11 = sum(t * t for t in temperatures_K)
    m12 = -sum(t * g for t, g in zip(temperatures_K, r3))
    m22 = sum(g * g for g in r3)
    matrix = [[m00, m01, m02], [m01, m11, m12], [m02, m12, m22]]

    b0 = sum(delta_g)
    b1 = -sum(t * g for t, g in zip(temperatures_K, delta_g))
    b2 = sum(g * dg for g, dg in zip(r3, delta_g))
    vector = [b0, b1, b2]

    delta_h_kJ, delta_s_kJ, delta_cp_kJ = _solve_3x3(
            matrix, vector,
            singular_message=(
                "temperature series is degenerate (singular fit matrix) -- "
                "need at least 3 points with genuinely different temperatures"
            ),
        )

    fitted = [delta_h_kJ - t * delta_s_kJ + delta_cp_kJ * g for t, g in zip(temperatures_K, r3)]
    mean_g = sum(delta_g) / n
    ss_res = sum((y - f) ** 2 for y, f in zip(delta_g, fitted))
    ss_tot = sum((y - mean_g) ** 2 for y in delta_g)
    r_squared = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 1.0

    return VantHoffMultiPointFitResult(
        delta_H_kJ_per_mol=delta_h_kJ,
        delta_S_J_per_mol_K=delta_s_kJ * 1000.0,
        delta_Cp_J_per_mol_K=delta_cp_kJ * 1000.0,
        reference_temperature_K=tr,
        r_squared=r_squared,
        n_points=n,
        delta_G_fitted_kJ_per_mol=fitted,
    )


def entropy_micellization(
    delta_g_mic_kJ_per_mol: float,
    delta_h_mic_kJ_per_mol: float,
    temperature_K: float,
    delta_cp_J_mol_K: float | None = None,
    reference_temperature_K: float | None = None,
) -> float:
    """Entropy of micellization, J/(mol.K): deltaS = (deltaH - deltaG) / T.
    Completes the deltaG/deltaH/deltaS triad given the other two.

    Optional deltaCp term (default off -- omit both optional arguments
    for the unchanged classic result): when delta_cp_J_mol_K AND
    reference_temperature_K are BOTH supplied, deltaG/deltaH are read as
    their values AT reference_temperature_K and the returned deltaS is
    reported AT temperature_K via the Kirchhoff relation (the exact
    form Kantonen, Henriksen & Gilson, BBA Gen. Subj. 1862 (2019) uses
    for their deltaG(T) expression, already cited by
    vant_hoff_multi_point_fit):

        deltaS(T) = deltaS(Tr) + deltaCp * ln(T / Tr)
        deltaS(Tr) = (deltaH(Tr) - deltaG(Tr)) / Tr

    Supplying only one of the two optional arguments -> ValueError
    naming both fields (never half-correct a reference silently);
    reference_temperature_K <= 0 or temperature_K <= 0 -> ValueError.
    With reference_temperature_K == temperature_K the result is
    identically the classic formula (ln(1) = 0)."""
    if temperature_K <= 0:
        raise ValueError("temperature_K must be positive")
    if delta_cp_J_mol_K is None and reference_temperature_K is None:
        delta_s_kJ = (delta_h_mic_kJ_per_mol - delta_g_mic_kJ_per_mol) / temperature_K
        return delta_s_kJ * 1000.0  # kJ/(mol.K) -> J/(mol.K)
    if delta_cp_J_mol_K is None or reference_temperature_K is None:
        raise ValueError(
            "delta_cp_J_mol_K and reference_temperature_K must be supplied "
            "together (deltaS(T) = deltaS(Tr) + deltaCp*ln(T/Tr) needs both "
            "the deltaCp and the reference the deltaG/deltaH are quoted at) "
            "-- or omit both for the classic deltaS = (deltaH - deltaG)/T"
        )
    if reference_temperature_K <= 0:
        raise ValueError("reference_temperature_K must be positive")
    delta_s_ref_J = (
        (delta_h_mic_kJ_per_mol - delta_g_mic_kJ_per_mol) / reference_temperature_K
    ) * 1000.0
    return delta_s_ref_J + delta_cp_J_mol_K * math.log(temperature_K / reference_temperature_K)


# ---------------------------------------------------------------------------
# Monte Carlo aggregation-number distribution (computational worker S1).
# Equilibrium cluster-size statistics for a canonical pool of associating
# monomers: Metropolis sampling of labeled partitions with bond energy
# epsilon, validated against exact integer-partition enumeration. Sampling a
# Boltzmann-weighted partition space is pure computation -- text-only models
# either refuse or invent histograms.
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class MCAggregationResult:
    """Monte Carlo cluster-size distribution for a canonical associating pool.

    ``histogram_n_monomer_fraction`` pairs (n, fraction of all monomers that
    reside in clusters of size exactly n); fractions sum to 1.
    """
    n_monomers: int
    temperature_K: float
    epsilon_kJ_per_mol: float
    beta_epsilon: float
    n_steps: int
    burn_in_fraction: float
    thin: int
    n_samples: int
    seed: int
    histogram_n_monomer_fraction: tuple[tuple[int, float], ...]
    free_monomer_fraction: float
    mean_cluster_size_monomer_weighted: float
    largest_cluster_observed: int
    split_half_l1: float
    split_half_noise_floor: float
    converged: bool
    exact_computed: bool
    exact_max_abs_deviation: float | None
    method: str
    params: dict
    caveat: str


def _integer_partitions_of_n(n: int):
    """Yield each integer partition of n as {size: multiplicity} dicts."""
    parts = [0] * (n + 1)

    def rec(remaining: int, mx: int):
        if remaining == 0:
            # NOTE: must scan ALL sizes up to n, not just this frame's mx --
            # outer frames contribute the larger parts
            yield {s: parts[s] for s in range(1, n + 1) if parts[s]}
            return
        for s in range(min(remaining, mx), 0, -1):
            parts[s] += 1
            yield from rec(remaining - s, s)
            parts[s] -= 1

    yield from rec(n, n)


def _exact_aggregation_marginals(n_monomers: int, beta_eps: float):
    """Exact monomer-fraction histogram over all integer partitions.

    Labeled-monomer cluster model: every set partition with cluster sizes
    {n_k} carries statistical weight  N! / prod_s((s!)^m_s m_s!)  times the
    Boltzmann factor exp(-beta E), E = -eps * (N - n_clusters).
    """
    import math as _math
    log_w = []
    margs = []  # per partition: list over n of fraction of monomers in size-n
    for m in _integer_partitions_of_n(n_monomers):
        n_clusters = sum(m.values())
        # log of N! / prod_s (s!)^m_s m_s!
        lw = _math.lgamma(n_monomers + 1)
        for s, ms in m.items():
            lw -= ms * _math.lgamma(s + 1) + _math.lgamma(ms + 1)
        lw += beta_eps * (n_monomers - n_clusters)  # -beta E = +beta eps bonds
        log_w.append(lw)
        row = [0.0] * (n_monomers + 1)
        for s, ms in m.items():
            row[s] = (ms * s) / n_monomers
        margs.append(row)
    mx = max(log_w)
    w = [_math.exp(v - mx) for v in log_w]
    z = sum(w)
    hist = [0.0] * (n_monomers + 1)
    for wi, row in zip(w, margs):
        for n in range(1, n_monomers + 1):
            hist[n] += wi * row[n] / z
    return hist


def mc_aggregation_distribution(
    n_monomers: int = 40,
    temperature_K: float = 298.15,
    epsilon_kJ_per_mol: float = 8.0,
    n_steps: int = 500_000,
    burn_in_fraction: float = 0.1,
    thin: int = 10,
    seed: int = 0,
    exact_if_possible: bool = True,
) -> MCAggregationResult:
    """Equilibrium cluster-size distribution via canonical Metropolis MC.

    Model: ``n_monomers`` labeled monomers in a box; an aggregate of size n
    has energy E(n) = -epsilon*(n-1) (n-1 bonds).  Moves: detach a random
    monomer to free, attach it to a random other cluster, or transfer it
    between clusters (Metropolis acceptance on the bond energy).  The
    histogram is the fraction of monomers residing in clusters of each size.

    When ``n_monomers <= 45`` and ``exact_if_possible``, the exact Boltzmann
    histogram over all integer partitions is computed alongside and the
    Monte Carlo deviation is reported in ``exact_max_abs_deviation`` --
    verification of the simulation against its own exact solution.

    Seed-controlled and reproducible.  ``split_half_l1`` compares the first
    and second half of the chain (``converged=False`` if poorly mixed).

    Requires numpy (surfactantkit[calc]).
    """
    import numpy as np

    if n_monomers < 2:
        raise ValueError("n_monomers must be >= 2")
    if temperature_K <= 0:
        raise ValueError("temperature_K must be positive")
    if epsilon_kJ_per_mol <= 0:
        raise ValueError(
            "epsilon_kJ_per_mol must be positive (binding energy per bond); "
            "the associating-cluster model is undefined without attraction"
        )
    if n_steps < 1000:
        raise ValueError("n_steps must be >= 1000 for usable statistics")
    if not 0.0 <= burn_in_fraction < 0.5:
        raise ValueError("burn_in_fraction must be in [0, 0.5)")
    if thin < 1:
        raise ValueError("thin must be a positive integer")

    rng = np.random.default_rng(seed)
    n = int(n_monomers)
    beta_eps = (epsilon_kJ_per_mol * 1000.0) / (8.314462618 * temperature_K)

    # state: parent[i] = cluster id of monomer i; clusters tracked as dicts
    # id -> set(labels).  All initially free (each monomer its own cluster).
    parent = list(range(n))
    members = {i: [i] for i in range(n)}

    def cluster_of(i: int) -> int:
        return parent[i]

    burn = int(n_steps * burn_in_fraction)
    hist_a = np.zeros(n + 1, dtype=np.int64)
    hist_b = np.zeros(n + 1, dtype=np.int64)
    half = burn + (n_steps - burn) // 2
    n_samples = 0
    largest = 1

    for step in range(n_steps):
        i = int(rng.integers(n))
        ci = cluster_of(i)
        roll = float(rng.random())
        if roll < 0.4:
            # detach i from its cluster (if cluster size > 1)
            if len(members[ci]) <= 1:
                continue  # already free: no move
            members[ci].remove(i)
            new_id = n + step  # fresh id for the liberated monomer
            members[new_id] = [i]
            parent[i] = new_id
            # Detailed balance with asymmetric proposals: the reverse
            # (attach) move must also choose its target cluster uniformly
            # among the K-1 clusters other than the new singleton, so
            # alpha_detach = min(1, exp(-beta*eps) / (K-1)) where K is the
            # cluster count after the detach.  (The reverse factor cancels
            # for attach and transfer moves -- see below -- but not here.)
            k_after = len(members)
            ratio = float(np.exp(-beta_eps)) / max(1, k_after - 1)
            if float(rng.random()) >= min(1.0, ratio):
                # reject: undo
                del members[new_id]
                members[ci].append(i)
                parent[i] = ci
        elif roll < 0.8:
            # attach free monomer i to a random other cluster
            if len(members[ci]) > 1:
                continue  # not free: this branch is the free->cluster move
            if len(members) <= 1:
                continue  # no target
            targets = [cid for cid in members if cid != ci]
            tgt = targets[int(rng.integers(len(targets)))]
            # E decreases by epsilon: always accepted
            members[tgt].append(i)
            del members[ci]
            parent[i] = tgt
        else:
            # transfer i between two distinct multi-monomer clusters (dE=0)
            if len(members[ci]) <= 1:
                continue
            others = [cid for cid in members if cid != ci and len(members[cid]) >= 1]
            if not others:
                continue
            tgt = others[int(rng.integers(len(others)))]
            members[ci].remove(i)
            members[tgt].append(i)
            parent[i] = tgt

        # maintain cluster bookkeeping for emptied clusters
        dead = [cid for cid, mem in members.items() if not mem]
        for cid in dead:
            del members[cid]

        if step >= burn and (step - burn) % thin == 0:
            n_samples += 1
            sizes = [len(mem) for mem in members.values()]
            largest = max(largest, max(sizes))
            sz_hist = np.zeros(n + 1, dtype=np.int64)
            for mem in members.values():
                sz_hist[len(mem)] += len(mem)
            target = hist_a if step < half else hist_b
            target += sz_hist

    if n_samples < 100:
        raise RuntimeError(
            f"too few samples collected ({n_samples}) -- increase n_steps"
        )

    total_a = hist_a.sum()
    hist_f = hist_a / total_a
    total_b = hist_b.sum()
    hist_fb = hist_b / total_b
    split_l1 = float(np.abs(hist_f - hist_fb).sum())
    # Expected noise floor of the two-half L1 distance under perfect mixing:
    # each bin's half-to-half difference is ~ Normal(0, 2*p*(1-p)/n_half), and
    # E|Z| = sqrt(2/pi)*sigma.  A converged chain should land near this floor.
    n_half_eff = max(1.0, total_a / n)
    p_pool = (hist_a + hist_b) / (total_a + total_b)
    noise_floor = float(
        sum(math.sqrt(2.0 / math.pi) *
            math.sqrt(max(0.0, 2.0 * p * (1.0 - p) / n_half_eff))
            for p in p_pool if p > 0.0)
    )

    free_fraction = float(hist_f[1])
    mean_n = float(sum(n_ * hist_f[n_] for n_ in range(1, n + 1)))

    exact_hist = None
    max_dev = None
    if exact_if_possible and n <= 45:
        eh = _exact_aggregation_marginals(n, beta_eps)
        exact_hist = eh
        max_dev = float(np.max(np.abs(hist_f - np.asarray(eh))))

    pairs = tuple((n_, float(hist_f[n_]))
                  for n_ in range(1, n + 1) if hist_f[n_] > 0.0)

    return MCAggregationResult(
        n_monomers=n,
        temperature_K=temperature_K,
        epsilon_kJ_per_mol=epsilon_kJ_per_mol,
        beta_epsilon=float(beta_eps),
        n_steps=int(n_steps),
        burn_in_fraction=float(burn_in_fraction),
        thin=int(thin),
        n_samples=int(n_samples),
        seed=int(seed),
        histogram_n_monomer_fraction=pairs,
        free_monomer_fraction=free_fraction,
        mean_cluster_size_monomer_weighted=mean_n,
        largest_cluster_observed=int(largest),
        split_half_l1=split_l1,
        split_half_noise_floor=noise_floor,
        converged=bool(split_l1 <= 2.5 * noise_floor),
        exact_computed=exact_hist is not None,
        exact_max_abs_deviation=max_dev,
        method="canonical Metropolis over labeled cluster partitions "
               "(detach/attach/transfer moves, bond energy epsilon) "
               + ("+ exact integer-partition enumeration"
                  if exact_hist is not None else "+ exact enumeration skipped "
                                                  "(N > 45)"),
        params={"n_monomers": n, "temperature_K": temperature_K,
                "epsilon_kJ_per_mol": epsilon_kJ_per_mol,
                "n_steps": n_steps, "seed": seed, "thin": thin,
                "burn_in_fraction": burn_in_fraction},
        caveat=(
            "Equilibrium statistics of a minimal bond-energy cluster model "
            "(no shape, no electrostatics, no solvent): it demonstrates "
            "cooperative size distributions and their temperature/epsilon "
            "dependence, not quantitative CMCs for a specific surfactant. "
            "Exact enumeration available only for N <= 45."
        ),
    )
