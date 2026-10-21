"""The actual point of this library: given ONLY a surfactant's real SMILES
and a raw pre-and-post-CMC surface-tension-vs-concentration curve -- the
same two things a real researcher has in hand before doing any analysis,
and nothing else (no stated ionic character, no named isotherm model, no
pre-fit parameters) -- autonomously derive every property this codebase
can legitimately derive, the way a competent researcher would work through
it by hand: classify the surfactant, find the CMC, decide whether the
pre-CMC adsorption is plain Langmuir or genuinely non-ideal (Frumkin), fit
whichever is actually justified, and compute the standard geometric and
thermodynamic quantities that follow.

Every other function in this library is deliberately atomic and refuses to
guess a single missing input (system_type, isotherm model, counterion
binding) rather than silently assuming one -- that discipline is a real
strength for a human calling one function with full context, but it means
NONE of those functions alone can answer "what can we learn from this
surfactant" without a human first making several judgment calls. This
module is that human's judgment, made explicit and inspectable: every
place a real ambiguity exists (zwitterionic/pH-dependent charge, unknown
electrolyte condition, no independent counterion-binding measurement) is
carried through to the final result as an explicit, named gap rather than
resolved by a hidden default -- the same "do not guess" discipline this
library already applies function-by-function, just applied once across
the whole pipeline instead of separately at every step.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from .adsorption import (
    IsothermModelSelectionResult,
    gibbs_a_min,
    gibbs_gamma_max,
    select_isotherm_model,
)
from .classify import ChargeClassificationResult, classify_surfactant_charge_type
from .curve_analysis import (
    CmcFromCurveResult,
    SzyszkowskiDeviationCmcResult,
    cmc_from_surface_tension_curve,
    cmc_via_szyszkowski_deviation,
    plateau_is_statistically_suspect,
)
from .thermodynamics import cmc_to_mole_fraction, gibbs_free_energy_micellization


@dataclass
class OrchestratorResult:
    smiles: str
    classification: ChargeClassificationResult
    system_type_used: str | None  # what was actually passed to Gibbs-prefactor-dependent steps, or None if undetermined
    cmc: CmcFromCurveResult  # RAW piecewise-linear result, unmodified even when cmc_mM below is a refined value -- see cmc_refinement
    cmc_mM: float | None  # the AUTHORITATIVE cmc used for every downstream step -- equals cmc.cmc_mM unless cmc_refinement below found a real refinement, or None when cmc_refinement.no_break_found_in_range (a real "no CMC detectable in this data" finding, not a computation failure)
    gamma_at_cmc_mN_per_m: float | None  # surface tension AT cmc_mM above (not cmc.gamma_at_cmc_mN_per_m's fitted value when refined) -- the real measured tension at the refined concentration when it is itself one of the input points, else the raw piecewise value; None exactly when cmc_mM is None
    cmc_refinement: SzyszkowskiDeviationCmcResult | None  # set only when cmc.cmc_mM's own plateau segment was statistically suspect (see plateau_is_statistically_suspect) and a Szyszkowski-deviation cross-check was run
    gamma_max_mol_per_m2: float | None
    a_min_nm2: float | None
    isotherm: IsothermModelSelectionResult | None
    delta_g_mic_kJ_per_mol: float | None
    counterion_factor_used: float | None
    gaps: list[str] = field(default_factory=list)  # named, honest list of everything NOT determined and why


def derive_all_properties_from_smiles_and_curve(
    smiles: str,
    concentrations_mM: list[float],
    surface_tensions_mN_per_m: list[float],
    temperature_K: float = 298.15,
    electrolyte_condition: str | None = None,
    counterion_dissociation_alpha: float | None = None,
    isotherm_a_bounds: tuple[float, float] = (-4.0, 4.0),
    gamma0_mN_m_override: float | None = None,
    gamma_max_mol_per_m2_override: float | None = None,
) -> OrchestratorResult:
    """Full (SMILES, raw curve) -> all-derivable-properties pipeline.

    gamma0_mN_m_override / gamma_max_mol_per_m2_override: optional, ONLY
    pass these if independently known (e.g. a real, separately measured
    pure-solvent surface tension, or a real Gamma_max from an independent
    technique such as neutron reflectivity or radiotracer measurement).
    By default this pipeline DERIVES both from the curve's own data
    (gamma0 from the detected baseline/lowest-concentration point,
    Gamma_max from the premicellar log-linear slope) -- real, standard
    practice, but a real independently-measured value is strictly better
    information when a researcher already has it, and silently ignoring
    it in favor of a same-curve-derived approximation is exactly the
    kind of "discard real given information" gap this project's own
    discipline exists to catch. When given,
    gamma0_mN_m_override/gamma_max_mol_per_m2_override feed directly
    into isotherm model selection (Langmuir vs Frumkin) and deltaG_mic
    instead of the curve-derived values -- a real, sometimes material
    difference. A_min_nm2/Gamma_max_mol_per_m2 in the returned result
    still ALSO report the curve-derived value (for comparison/audit),
    but the isotherm fit and deltaG_mic use the override when supplied.

    electrolyte_condition: optional override, one of 'excess_electrolyte'
    or 'no_added_salt' -- ONLY pass this if it is independently known
    about the real solution being characterized (e.g. stated in the
    experiment's own description). A bare SMILES carries no information
    about what else is dissolved in the water, so classify.py's own
    default ('ionic_no_added_salt', i.e. Gibbs prefactor n=2) is used
    for any ionic surfactant when this is left as None -- surfaced as an
    explicit gap in the result rather than silently assumed to be
    correct, since a real excess-electrolyte system would need n=1
    instead and get a genuinely different Gamma_max.

    counterion_dissociation_alpha: optional, ONLY pass this if
    independently measured (e.g. via the conductivity slope-ratio method,
    counterion_binding_degree, or a literature/EPR value like Bales et
    al.'s alpha=0.272 for SDS -- see gibbs_free_energy_micellization's
    own docstring). It is the degree of IONIZATION (S2/S1 = 1 - beta), and
    deltaG_mic uses the counterion factor 2 - alpha (= 1 + beta; 1.728 for
    SDS). Without it, deltaG_mic for an IONIC surfactant is
    left undetermined (reported as a named gap) rather than guessed --
    this is real information a surface-tension curve alone cannot supply,
    the same reasoning gibbs_free_energy_micellization's own docstring
    already applies.

    isotherm_a_bounds: passed through to select_isotherm_model.
    """
    gaps: list[str] = []

    classification = classify_surfactant_charge_type(smiles)

    if classification.charge_type == "unparseable":
        system_type: str | None = None
        gaps.append(f"SMILES could not be parsed ('{smiles}') -- charge type, Gamma_max, A_min, "
                     "isotherm model, and deltaG_mic are all undetermined until a valid SMILES is supplied.")
    elif classification.charge_type in ("zwitterionic", "ambiguous_pH_dependent"):
        system_type = None
        gaps.append(f"charge_type='{classification.charge_type}' -- {classification.caveats[0] if classification.caveats else 'ambiguous from structure alone'} "
                     "Gamma_max, A_min, and the isotherm fit all require a known Gibbs prefactor, which cannot "
                     "be assigned here without resolving this ambiguity (e.g. independently confirming solution pH, "
                     "or the specific zwitterionic-system literature precedent for its prefactor convention).")
    else:
        system_type = classification.system_type_for_gibbs  # 'nonionic' or 'ionic_no_added_salt'
        if electrolyte_condition == "excess_electrolyte" and system_type == "ionic_no_added_salt":
            system_type = "ionic_excess_electrolyte"
        elif electrolyte_condition == "no_added_salt":
            pass  # already the default for ionic; no-op for nonionic
        elif electrolyte_condition is None and system_type == "ionic_no_added_salt":
            gaps.append("electrolyte_condition not supplied for an ionic surfactant -- defaulted to "
                        "'ionic_no_added_salt' (Gibbs prefactor n=2), the correct choice ONLY if the real "
                        "solution truly has no significant excess inert electrolyte; if it does, re-run with "
                        "electrolyte_condition='excess_electrolyte' (n=1) for a materially different Gamma_max.")
        if classification.caveats:
            gaps.extend(f"classification note ({classification.charge_type}): {c}" for c in classification.caveats)

    cmc = cmc_from_surface_tension_curve(concentrations_mM, surface_tensions_mN_per_m)
    cmc_mM = cmc.cmc_mM
    gamma_at_cmc_mN_per_m = cmc.gamma_at_cmc_mN_per_m
    cmc_refinement: SzyszkowskiDeviationCmcResult | None = None

    gamma_max = a_min = None
    isotherm: IsothermModelSelectionResult | None = None
    if system_type is not None:
        if plateau_is_statistically_suspect(cmc):
            cmc_refinement = cmc_via_szyszkowski_deviation(concentrations_mM, surface_tensions_mN_per_m, cmc, system_type, temperature_K)
            if cmc_refinement.cmc_mM is not None:
                cmc_mM = cmc_refinement.cmc_mM
                # cmc_via_szyszkowski_deviation always returns a real measured concentration
                # (the first point past the premicellar range that deviates), so this is a
                # direct lookup, never an extrapolation.
                gamma_at_cmc_mN_per_m = next(
                    g for c, g in zip(concentrations_mM, surface_tensions_mN_per_m) if c == cmc_mM
                )
                gaps.append(
                    f"cmc_mM refined from {cmc.cmc_mM:.4g} mM (piecewise break-point) to {cmc_mM:.4g} mM "
                    f"(Szyszkowski-isotherm deviation cross-check): the piecewise algorithm's own "
                    f"postmicellar segment was statistically suspect (too few points and/or not "
                    f"meaningfully flatter than the premicellar decline -- see "
                    f"plateau_is_statistically_suspect) -- {cmc_refinement.reason}."
                )
            elif cmc_refinement.no_break_found_in_range:
                # Real, distinct finding from "refined but still uncertain" above: the
                # deviation check didn't just fail to IMPROVE on the piecewise value, it
                # found the ENTIRE dataset consistent with a single continuous premicellar
                # isotherm -- i.e. real evidence there is NO detectable CMC in this
                # concentration range at all. Falling back to the untrustworthy piecewise
                # cmc_mM here would keep wrongly truncating the premicellar window. cmc_mM
                # becomes None and every given point is treated as premicellar instead of
                # trusting an artifact boundary.
                cmc_mM = None
                gamma_at_cmc_mN_per_m = None
                gaps.append(
                    f"No CMC could be reliably detected anywhere in this concentration range: the "
                    f"piecewise algorithm's own {cmc.cmc_mM:.4g} mM is a statistically suspect "
                    f"artifact (see plateau_is_statistically_suspect), and the Szyszkowski-deviation "
                    f"cross-check found the ENTIRE dataset consistent with a single continuous "
                    f"premicellar decline ({cmc_refinement.reason}) -- not just an inconclusive check, "
                    f"real evidence against a break existing in this range at all. Treating every given "
                    f"point as premicellar for Gamma_max/isotherm-selection purposes below, rather than "
                    f"truncating to a subset bounded by an artifact CMC."
                )
            else:
                gaps.append(
                    f"cmc_mM={cmc.cmc_mM:.4g} mM comes from a statistically suspect piecewise plateau "
                    f"segment (see plateau_is_statistically_suspect), and the Szyszkowski-deviation "
                    f"cross-check could not locate an alternative break either ({cmc_refinement.reason}) "
                    f"-- treat this CMC value with real caution."
                )

        gamma_max = gibbs_gamma_max(cmc.premicellar_slope_mN_per_m_per_lnC, system_type, temperature_K)
        a_min = gibbs_a_min(gamma_max)

        gamma0_for_isotherm = cmc.gamma0_mN_m
        gamma_max_for_isotherm = gamma_max
        if gamma0_mN_m_override is not None:
            gamma0_for_isotherm = gamma0_mN_m_override
            gaps.append(f"gamma0 for isotherm fitting/deltaG_mic uses the supplied override "
                        f"{gamma0_mN_m_override:.4g} mN/m, not the curve-derived {cmc.gamma0_mN_m:.4g} mN/m.")
        if gamma_max_mol_per_m2_override is not None:
            gamma_max_for_isotherm = gamma_max_mol_per_m2_override
            gaps.append(f"Gamma_max for isotherm fitting uses the supplied override "
                        f"{gamma_max_mol_per_m2_override:.4g} mol/m^2, not the curve-derived "
                        f"{gamma_max:.4g} mol/m^2 reported above (A_min above is still the curve-derived value).")

        if cmc_mM is None:
            premicellar_pairs = sorted(zip(concentrations_mM, surface_tensions_mN_per_m))
        else:
            premicellar_pairs = sorted(
                (c, g) for c, g in zip(concentrations_mM, surface_tensions_mN_per_m)
                if cmc.premicellar_x_min_mM <= c < cmc_mM
            )
        if len(premicellar_pairs) >= 4:
            pre_c = [p[0] for p in premicellar_pairs]
            pre_g = [p[1] for p in premicellar_pairs]
            isotherm = select_isotherm_model(pre_c, pre_g, gamma0_for_isotherm, gamma_max_for_isotherm, system_type, temperature_K, a_bounds=isotherm_a_bounds)
            gaps.extend(f"isotherm fit note: {w}" for w in isotherm.warnings)
        else:
            gaps.append(f"only {len(premicellar_pairs)} pre-CMC data point(s) available (need >=4) -- "
                        "Langmuir-vs-Frumkin isotherm model selection skipped; Gamma_max/A_min above come "
                        "from the curve's own premicellar slope directly and do not depend on this.")
    else:
        gaps.append("Gamma_max, A_min, and isotherm model selection all skipped: system_type undetermined (see above).")

    delta_g_mic = None
    counterion_factor_used = None
    if cmc_mM is None:
        if system_type is not None:
            gaps.append("deltaG_mic skipped: no CMC could be reliably determined in this concentration "
                        "range (see the earlier gap explaining why) -- there is no real cmc_mM to convert "
                        "to a mole fraction.")
    elif system_type == "nonionic":
        counterion_factor_used = 1.0
        delta_g_mic = gibbs_free_energy_micellization(cmc_to_mole_fraction(cmc_mM / 1000.0), temperature_K, counterion_factor_used)
    elif system_type in ("ionic_no_added_salt", "ionic_excess_electrolyte"):
        if counterion_dissociation_alpha is not None:
            # alpha = degree of IONIZATION (S2/S1); pseudo-phase factor = 2 - alpha = 1 + beta.
            # Was 1 + alpha until 2026-10-20 (inverse convention; see gibbs_free_energy_micellization).
            counterion_factor_used = 2.0 - counterion_dissociation_alpha
            delta_g_mic = gibbs_free_energy_micellization(cmc_to_mole_fraction(cmc_mM / 1000.0), temperature_K, counterion_factor_used)
        else:
            gaps.append("deltaG_mic skipped for this ionic surfactant: no independent counterion-binding "
                        "measurement supplied (counterion_dissociation_alpha). This is real information a "
                        "surface-tension curve and SMILES alone cannot provide -- it requires an independent "
                        "measurement (e.g. conductivity slope-ratio method, see counterion_binding_degree) or "
                        "a literature value for this specific surfactant. Do not assume alpha; report the gap.")
    elif system_type is not None:
        gaps.append(f"deltaG_mic skipped: unrecognized system_type '{system_type}'.")

    return OrchestratorResult(
        smiles=smiles, classification=classification, system_type_used=system_type, cmc=cmc,
        cmc_mM=cmc_mM, gamma_at_cmc_mN_per_m=gamma_at_cmc_mN_per_m, cmc_refinement=cmc_refinement,
        gamma_max_mol_per_m2=gamma_max, a_min_nm2=a_min, isotherm=isotherm,
        delta_g_mic_kJ_per_mol=delta_g_mic, counterion_factor_used=counterion_factor_used, gaps=gaps,
    )
