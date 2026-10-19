"""Micellar solubilization capacity.

Originally implemented only the Molar Solubilization Ratio (MSR), which
has one clean, unambiguous definition in the literature -- the
micelle-water partition coefficient (Km) was deliberately left out,
because a check of the literature found its exact formula is study-
specific (some define it on a mole fraction basis, others on a molar-
concentration-ratio basis, with no single dominant convention), and
implementing one silently risked reproducing the wrong convention for a
given paper.

UPDATE, 2026-09-10 (alternative-methods sweep): re-checked, and the
ambiguity is real and still unresolved -- but Km is also a real,
commonly co-reported COMPLEMENTARY metric to MSR (this project's own
literature review: apparent solubilization capacity depends heavily on
WHERE the solute sits, palisade/interface vs. core -- MSR alone doesn't
capture this). Rather than continue omitting it entirely, this module
now implements the MOLE-FRACTION-BASIS convention specifically (X_mic/
X_w, one of the two real conventions found, not invented) with the
ambiguity explicitly disclosed in micelle_water_partition_coefficient's
own docstring -- letting a researcher use it with eyes open, matching
this project's "let a researcher choose the method suited to their
system" goal, rather than a blanket omission. Same discipline as the
missing Davies HLB groups and the unimplemented Krafft point for
anything that remains genuinely unresolved.
"""

from __future__ import annotations

from dataclasses import dataclass, field

WATER_MOLARITY_M = 55.5  # mol/L, same dilute-aqueous-solution convention as thermodynamics.py

# Tier 2 bottleneck-resolution addition, 2026-09-15 (see benchmark/
# audit findings, 2026-09 item 13): real
# intrinsic aqueous solubility values for common solubilizates used in
# real micellar-solubilization studies, closing this bottleneck for
# compounds molar_solubilization_ratio/micelle_water_partition_
# coefficient's own intrinsic_water_solubility_M input would otherwise
# need externally supplied for every use.
#
# naphthalene: 2.17e-4 mol/L (27.8 mg/L), the SAME real, already-used
# value from this project's own literature_validation_notes.md Round 5
# (Paria & Yuet, Ind. Eng. Chem. Res. 45(10) (2006) 3552-3558, the
# paper's OWN measured value, stated there as consistent with an
# independent literature value of 2.44e-4 mol/L) -- HIGH confidence, a
# real primary-source value already validated end-to-end against this
# project's own molar_solubilization_ratio function.
#
# benzene and pyrene: found via real web search 2026-09-15, both cross-
# checked against multiple independent citing sources in the same
# search (benzene ~1.8 g/L / ~1780 mg/L at 25C, a very widely tabulated
# constant, e.g. EPA/ATSDR-compiled; pyrene 0.139 mg/L at 25C, a
# specific experimental value found via a real ACS solubility-
# measurement paper's own reported table) -- MODERATE confidence: real,
# cross-checked search results, but neither primary source's own table
# was independently re-fetched and read in full this session, unlike
# naphthalene's already-established primary-source status.
INTRINSIC_WATER_SOLUBILITY_REFERENCE_M = {
    "naphthalene": {"solubility_M": 2.17e-4, "mw_g_per_mol": 128.17,
                     "source": "Paria & Yuet 2006, Ind. Eng. Chem. Res. 45(10) 3552-3558 (primary PDF read in full, already used in this project's own molar_solubilization_ratio validation)",
                     "confidence": "high"},
    "benzene": {"solubility_M": 1780e-3 / 78.11, "mw_g_per_mol": 78.11,
                "source": "web search 2026-09-15, cross-corroborated across multiple citing sources (~1.8 g/L, EPA/ATSDR-compiled); primary source not independently fetched",
                "confidence": "moderate"},
    "pyrene": {"solubility_M": 0.139e-3 / 202.25, "mw_g_per_mol": 202.25,
               "source": "web search 2026-09-15, real experimental value (0.139 mg/L at 25C) via a real ACS solubility-measurement paper; primary table not independently fetched",
               "confidence": "moderate"},
}


@dataclass
class QsprSolubilityEstimate:
    smiles: str
    log_s_mol_per_L: float
    solubility_M: float
    method: str = "ESOL (Delaney 2004), RDKit-refit coefficients (Walters)"
    confidence: str = "low-moderate"
    caveats: list[str] = field(default_factory=list)


def estimate_intrinsic_water_solubility_qspr(smiles: str) -> QsprSolubilityEstimate:
    """Real, buildable [COMPUTE]/QSPR estimate of a solubilizate's
    intrinsic aqueous solubility from structure alone -- closes part of
    the intrinsic_water_solubility_M bottleneck for
    molar_solubilization_ratio/micelle_water_partition_coefficient (see
    audit findings, 2026-09 item 13),
    for compounds no real measured value has been sourced for yet.

    Method: ESOL (Delaney, J. Chem. Inf. Comput. Sci. 2004, 44, 1000-1005)
    -- a real, published, widely-used linear QSPR fit on 2874 measured
    aqueous solubilities:

        log10(S, mol/L) = a - b*clogP - c*MW + d*RB - e*AP

    Coefficients used here are Pat Walters' RDKit-specific REFIT of
    Delaney's original model (practicalcheminformatics.blogspot.com,
    "Predicting Aqueous Solubility -- It's Harder Than It Looks", 2018;
    code: github.com/PatWalters/solubility), not Delaney's original
    published coefficients -- deliberately, because Delaney's own
    coefficients were fit against his own (Molconn-Z-derived) clogP
    values, which are NOT numerically identical to RDKit's Crippen.MolLogP
    used here; using RDKit's own descriptors with Delaney's original
    coefficients would silently mix two different logP implementations.
    Walters' refit was specifically calibrated for RDKit's own
    Descriptors.MolWt, Crippen.MolLogP, Lipinski.NumRotatableBonds, and
    an aromatic-atom-fraction proportion -- live-verified 2026-09-15
    before being hardcoded here.

    HONEST, PROMINENT LIMITATION: this is a coarse QSPR ESTIMATE, not a
    measurement -- Delaney's own validation reports roughly 0.6-1 log-
    unit RMSE (i.e. this can be off by up to an order of magnitude or
    more for an individual compound), and is real, published estimation
    error, not a defect of this implementation. Use a real measured
    value (or this project's own mined literature values) whenever one
    is available; treat this only as a fallback for compounds nothing
    else has been sourced for, always disclosed as an estimate, never
    silently substituted for a real measurement.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Crippen, Descriptors, Lipinski
    except ImportError as e:  # pragma: no cover
        raise ImportError(
            "estimate_intrinsic_water_solubility_qspr requires rdkit (pip install rdkit) "
            "-- the ESOL QSPR model needs real molecular descriptors (clogP, MW, rotatable "
            "bonds, aromatic proportion), not hand-rolled string parsing."
        ) from e

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse this SMILES: '{smiles}'. Check for typos/invalid valence.")

    mw = Descriptors.MolWt(mol)
    clogp = Crippen.MolLogP(mol)
    rotors = Lipinski.NumRotatableBonds(mol)
    aromatic_pattern = Chem.MolFromSmarts("a")
    n_aromatic = len(mol.GetSubstructMatches(aromatic_pattern))
    n_heavy = mol.GetNumAtoms()
    aromatic_proportion = (n_aromatic / n_heavy) if n_heavy > 0 else 0.0

    log_s = (
        0.26121066137801696
        - 0.7416739523408995 * clogp
        - 0.0066138847738667125 * mw
        + 0.003451545565957996 * rotors
        - 0.42624840441316975 * aromatic_proportion
    )
    solubility_M = 10.0 ** log_s

    return QsprSolubilityEstimate(
        smiles=smiles,
        log_s_mol_per_L=log_s,
        solubility_M=solubility_M,
        caveats=[
            ("This is a QSPR ESTIMATE (ESOL, RDKit-refit coefficients), not a measurement -- "
            "real published validation error is roughly 0.6-1 log unit (up to an order of "
            "magnitude or more for an individual compound). Prefer a real measured/cited "
            "value when one is available for this specific solubilizate."),
        ],
    )


def molar_solubilization_ratio(total_solubilized_M: float, intrinsic_water_solubility_M: float, surfactant_concentration_M: float, cmc_M: float) -> float:
    """Molar Solubilization Ratio (MSR, dimensionless): moles of
    solubilizate taken up per mole of MICELLIZED surfactant (i.e. only
    the surfactant above the CMC, which is what actually forms micelles):
    MSR = (total_solubilized - intrinsic_water_solubility) / (surfactant_concentration - cmc).

    All four inputs must be in the same concentration unit (e.g. all mM
    or all mol/L). surfactant_concentration must exceed cmc -- below the
    CMC there are no micelles to solubilize anything into.
    """
    if surfactant_concentration_M <= cmc_M:
        raise ValueError("surfactant_concentration_M must exceed cmc_M (no micelles form below the CMC)")
    if total_solubilized_M < intrinsic_water_solubility_M:
        raise ValueError("total_solubilized_M cannot be less than intrinsic_water_solubility_M")
    if intrinsic_water_solubility_M < 0:
        raise ValueError("intrinsic_water_solubility_M must be non-negative")
    return (total_solubilized_M - intrinsic_water_solubility_M) / (surfactant_concentration_M - cmc_M)


def micelle_water_partition_coefficient(
    total_solubilized_M: float,
    intrinsic_water_solubility_M: float,
    surfactant_concentration_M: float,
    cmc_M: float,
    water_molarity_M: float = WATER_MOLARITY_M,
) -> float:
    """Micelle-water partition coefficient (Km, dimensionless), MOLE-
    FRACTION-BASIS convention: Km = X_solubilizate,micelle /
    X_solubilizate,water, a real, commonly co-reported complementary
    metric to molar_solubilization_ratio -- captures how STRONGLY the
    solubilizate favors the micellar pseudophase over the bulk aqueous
    phase (an intensive, concentration-normalized quantity), whereas MSR
    is an extensive uptake-per-mole-of-surfactant quantity. A high Km
    with a low MSR (or vice versa) is real, meaningful information MSR
    alone does not capture (this project's own literature review: Ianiro
    et al., Langmuir 2019, PMC6448116 -- apparent solubilization
    capacity depends heavily on WHERE the solute sits, palisade/
    interface vs. core).

    IMPORTANT, disclosed honestly: this is ONE of (at least) two real
    conventions found in the literature for Km -- some papers instead
    use a molar-CONCENTRATION-RATIO basis (Km = C_solubilizate,micelle /
    C_solubilizate,water, no mole-fraction normalization), which gives a
    numerically DIFFERENT value from this mole-fraction-basis
    calculation for the same raw data. There is no single dominant
    convention -- verify which one a given paper you are comparing
    against actually uses before treating numbers as directly
    comparable; do not assume they match.

    Same four raw-data inputs as molar_solubilization_ratio (see that
    function for their meaning), reused directly for consistency:

        X_micelle = (total_solubilized - intrinsic_water_solubility)
                    / [(total_solubilized - intrinsic_water_solubility)
                       + (surfactant_concentration - cmc)]
        X_water = intrinsic_water_solubility / water_molarity_M
        Km = X_micelle / X_water

    water_molarity_M: dilute-water-molarity approximation (default 55.5
    mol/L, same convention as thermodynamics.cmc_to_mole_fraction).
    intrinsic_water_solubility_M is taken as the actual aqueous-phase
    solubilizate concentration at equilibrium (valid when excess
    undissolved solute is present, saturating the water phase at its own
    intrinsic solubility -- the same assumption molar_solubilization_ratio
    already makes for this same input).
    """
    if surfactant_concentration_M <= cmc_M:
        raise ValueError("surfactant_concentration_M must exceed cmc_M (no micelles form below the CMC)")
    if total_solubilized_M < intrinsic_water_solubility_M:
        raise ValueError("total_solubilized_M cannot be less than intrinsic_water_solubility_M")
    if intrinsic_water_solubility_M <= 0:
        raise ValueError("intrinsic_water_solubility_M must be positive (X_water is undefined at zero)")
    if water_molarity_M <= 0:
        raise ValueError("water_molarity_M must be positive")

    n_solubilized_micelle = total_solubilized_M - intrinsic_water_solubility_M
    n_surfactant_micellized = surfactant_concentration_M - cmc_M
    x_micelle = n_solubilized_micelle / (n_solubilized_micelle + n_surfactant_micellized)
    x_water = intrinsic_water_solubility_M / water_molarity_M
    return x_micelle / x_water
