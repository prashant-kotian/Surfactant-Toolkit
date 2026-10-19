"""Structure-based surfactant charge-type classification from SMILES.

This closes a real, previously-missing capability: every function in this
project that needs to know whether a surfactant is ionic or nonionic
(gibbs_gamma_max, szyszkowski_fit_K, hld_ionic/hld_nonionic, ...) has always
required that as an explicit, human-supplied string parameter -- there was
no way to determine it from the surfactant's own structure. A researcher
handing over a real dataset AND the surfactant's real SMILES should not
also have to separately tell the toolkit "this is ionic" -- the toolkit
should be able to work that out itself, the same way a real chemist reading
a structure would.

Method: RDKit SMARTS substructure matching for the functional groups that
actually determine surfactant charge behavior in water -- sulfate esters,
sulfonates, carboxylates/carboxylic acids, phosphate MONOesters AND
phosphate DIesters (real, genuine gap found and fixed 2026-09-15: the
original pattern only matched a monoester, C-O-P(=O)(O-)(O-) -- a real
phospholipid/phosphatidylcholine-type zwitterionic surfactant like DHPC has
a phosphoDIester instead, C-O-P(=O)(O-)-O-C, which the monoester-only
pattern silently missed entirely, misclassifying an otherwise-correctly-
detected quaternary-ammonium-bearing phospholipid as plain "cationic"
instead of "zwitterionic" -- caught while building real Paper 3 ground-
zero benchmark questions against DHPC, not hypothesized), quaternary
AROMATIC nitrogen -- pyridinium/imidazolium-type ionic-liquid headgroups
(real, second genuine gap found and fixed 2026-09-18: the original
cationic patterns only matched sp3 [NX4+], silently missing a real
N-alkylated aromatic ring nitrogen like pyridinium's [n+] entirely --
misclassifying real pyridinium (e.g. the Fu et al. 2019 ionic-liquid
surfactant [C12mpy][Br]) AND real imidazolium (e.g. CMIC, 1-hexadecyl-
3-methylimidazolium, already used elsewhere in this project's own HLD-NAC
reference data) as "nonionic" -- caught the same way, while building GZ-26.
Verified the new pattern correctly rejects plain neutral pyridine/imidazole,
which have no alkyl substituent on the ring nitrogen and are not
permanently charged) (anionic,
strong-to-moderate acids, treated as ionized at normal surfactant testing
pH regardless of how the SMILES happens to write the protonation state,
since that's how these groups actually behave in real aqueous
formulations), quaternary ammonium (permanently cationic, pH-independent),
and free primary/secondary/tertiary amines (only conditionally cationic,
protonated below their pKa -- this ambiguity is disclosed explicitly rather
than resolved by guessing, matching this project's standing discipline of
refusing to silently assume unstated information).

A molecule matching both an anionic and a permanently-cationic (quaternary
ammonium) group is classified zwitterionic. A molecule matching neither is
classified nonionic -- but this is a claim by absence, so the low-confidence
case (a completely unfamiliar backbone with no recognized group at all) is
flagged, not silently trusted.
"""
from __future__ import annotations

from dataclasses import dataclass, field

try:
    from rdkit import Chem
except ImportError as e:  # pragma: no cover
    raise ImportError(
        "classify.py requires rdkit (pip install rdkit) -- this is the one "
        "module in SurfactantKit with an external dependency, since correct "
        "functional-group detection from arbitrary real SMILES genuinely "
        "needs a real cheminformatics substructure-matching engine, not "
        "hand-rolled string parsing (which would be fragile against atom "
        "ordering/branching in real SMILES)."
    ) from e


# --- SMARTS patterns, one per real functional group that determines charge
# behavior in water. Each is a (name, smarts, charge_contribution, strength)
# tuple. "strength": "strong" = permanently/always ionized at normal aqueous
# pH (sulfate, sulfonate, quaternary ammonium); "moderate" = ionized at
# typical surfactant-testing pH (>~pH 6) but not at low pH (carboxylic
# acid/carboxylate, pKa ~4-5); "conditional" = only ionized below the
# group's own pKa, genuinely pH-dependent, disclosed as an open question
# rather than resolved (free amines, pKa ~9-11 so USUALLY protonated at
# neutral pH, but this project does not silently assume that).
_ANIONIC_PATTERNS = [
    ("sulfate_ester", "[#6]-[OX2]-[SX4](=[OX1])(=[OX1])-[OX1,OX2H0-1,OX2H1]", "strong"),
    ("sulfonate", "[#6]-[SX4](=[OX1])(=[OX1])-[OX1,OX2H0-1,OX2H1]", "strong"),
    ("phosphate_ester", "[#6]-[OX2]-[PX4](=[OX1])([OX1,OX2H,OX2H0-1])[OX1,OX2H,OX2H0-1]", "strong"),
    ("phosphate_diester", "[#6]-[OX2]-[PX4](=[OX1])([OX1,OX2H,OX2H0-1])[OX2]-[#6]", "strong"),
    ("carboxylate", "[CX3](=[OX1])[OX1H0-1,OX2H1]", "moderate"),
]
_CATIONIC_STRONG_PATTERNS = [
    ("quaternary_ammonium", "[#6][NX4+]([#6])([#6])[#6]", "strong"),
    ("quaternary_ammonium_ring", "[#6][NX4+;R]", "strong"),
    ("aromatic_quaternary_nitrogen", "[#6][n+;R]", "strong"),
]
_CATIONIC_CONDITIONAL_PATTERNS = [
    ("primary_amine", "[NX3H2;!$(NC=O);!$(N=*)]", "conditional"),
    ("secondary_amine", "[NX3H1;!$(NC=O);!$(N=*);!R]", "conditional"),
    ("tertiary_amine", "[NX3H0;!$(NC=O);!$(N=*);!$([N+])]", "conditional"),
]

_COMPILED_ANIONIC = [(name, Chem.MolFromSmarts(sm), strength) for name, sm, strength in _ANIONIC_PATTERNS]
_COMPILED_CATIONIC_STRONG = [(name, Chem.MolFromSmarts(sm), strength) for name, sm, strength in _CATIONIC_STRONG_PATTERNS]
_COMPILED_CATIONIC_COND = [(name, Chem.MolFromSmarts(sm), strength) for name, sm, strength in _CATIONIC_CONDITIONAL_PATTERNS]


@dataclass
class ChargeClassificationResult:
    smiles: str
    charge_type: str  # "anionic" | "cationic" | "zwitterionic" | "nonionic" | "ambiguous_pH_dependent" | "unparseable"
    system_type_for_gibbs: str | None  # "ionic_no_added_salt" | "ionic_excess_electrolyte" | "nonionic" | None if undetermined
    anionic_groups_found: list[str] = field(default_factory=list)
    cationic_strong_groups_found: list[str] = field(default_factory=list)
    cationic_conditional_groups_found: list[str] = field(default_factory=list)
    confidence: str = ""  # "high" | "moderate" | "low" -- honestly stated, not hidden
    caveats: list[str] = field(default_factory=list)


def classify_surfactant_charge_type(smiles: str) -> ChargeClassificationResult:
    """Determine a surfactant's charge type (anionic/cationic/zwitterionic/
    nonionic) from its real SMILES structure via functional-group matching.

    Does NOT guess when the structure is genuinely ambiguous (e.g. a free
    amine with no quaternary nitrogen, whose charge state depends on a pH
    this function is never given) -- returns "ambiguous_pH_dependent" with
    the caveat spelled out, rather than silently picking cationic or
    nonionic, matching this project's standing "do not guess system_type"
    discipline already applied throughout adsorption.py/hld.py.

    system_type_for_gibbs is provided as a convenience default mapping onto
    the exact string this project's other functions (gibbs_gamma_max,
    szyszkowski_fit_K, ...) already expect -- but note it defaults ionic
    surfactants to "ionic_no_added_salt" (Gibbs prefactor n=2), since a bare
    SMILES carries no information about whether excess inert electrolyte is
    present in the actual solution being characterized; the caller must
    override this to "ionic_excess_electrolyte" if that is independently
    known to be true for their real system.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ChargeClassificationResult(
            smiles=smiles, charge_type="unparseable", system_type_for_gibbs=None,
            confidence="none",
            caveats=[f"RDKit could not parse this SMILES: '{smiles}'. Check for typos/invalid valence."],
        )

    anionic_found = [name for name, patt, _ in _COMPILED_ANIONIC if patt is not None and mol.HasSubstructMatch(patt)]
    cationic_strong_found = [name for name, patt, _ in _COMPILED_CATIONIC_STRONG if patt is not None and mol.HasSubstructMatch(patt)]
    cationic_cond_found = [name for name, patt, _ in _COMPILED_CATIONIC_COND if patt is not None and mol.HasSubstructMatch(patt)]

    is_anionic = len(anionic_found) > 0
    is_cationic_strong = len(cationic_strong_found) > 0
    is_cationic_conditional_only = (not is_cationic_strong) and len(cationic_cond_found) > 0

    caveats: list[str] = []

    if is_anionic and is_cationic_strong:
        return ChargeClassificationResult(
            smiles=smiles, charge_type="zwitterionic", system_type_for_gibbs=None,
            anionic_groups_found=anionic_found, cationic_strong_groups_found=cationic_strong_found,
            cationic_conditional_groups_found=cationic_cond_found, confidence="high",
            caveats=[("Zwitterionic surfactants are net-neutral overall; the Gibbs prefactor "
                     "convention for zwitterionics is a separate real question this function "
                     "does not resolve -- do not default to either the ionic or nonionic "
                     "prefactor without checking the specific system's own literature "
                     "precedent (e.g. this project's own sulfobetaine/betaine validation "
                     "entries).")],
        )

    if is_anionic and is_cationic_conditional_only:
        return ChargeClassificationResult(
            smiles=smiles, charge_type="ambiguous_pH_dependent", system_type_for_gibbs=None,
            anionic_groups_found=anionic_found, cationic_conditional_groups_found=cationic_cond_found,
            confidence="low",
            caveats=[(f"Contains both a real anionic group ({anionic_found}) and a free amine "
                     f"({cationic_cond_found}) with no stated pH -- this could be zwitterionic "
                     "(amino-acid-type), purely anionic (if the amine is unprotonated at the "
                     "real solution pH), or something else. Cannot be resolved from structure "
                     "alone; report the pH-dependence rather than guessing.")],
        )

    if is_anionic:
        strengths = {name: s for name, _, s in _ANIONIC_PATTERNS if name in anionic_found}
        if all(strengths.get(n) == "moderate" for n in anionic_found):
            caveats.append("Only a carboxylate/carboxylic acid group found (pKa ~4-5) -- "
                            "this is ionized (anionic) at typical surfactant-testing pH (>~6) "
                            "but would behave as a nonionic (protonated) acid at low pH. "
                            "Confirm the real solution pH if operating near the pKa.")
        return ChargeClassificationResult(
            smiles=smiles, charge_type="anionic", system_type_for_gibbs="ionic_no_added_salt",
            anionic_groups_found=anionic_found, confidence="high" if "moderate" not in strengths.values() else "moderate",
            caveats=caveats,
        )

    if is_cationic_strong:
        return ChargeClassificationResult(
            smiles=smiles, charge_type="cationic", system_type_for_gibbs="ionic_no_added_salt",
            cationic_strong_groups_found=cationic_strong_found, confidence="high", caveats=caveats,
        )

    if is_cationic_conditional_only:
        return ChargeClassificationResult(
            smiles=smiles, charge_type="ambiguous_pH_dependent", system_type_for_gibbs=None,
            cationic_conditional_groups_found=cationic_cond_found, confidence="low",
            caveats=[(f"Only a free amine ({cationic_cond_found}, typical pKa ~9-11) found, no "
                     "quaternary (permanently charged) nitrogen. This is USUALLY protonated "
                     "(cationic) at neutral/acidic pH, but this function does not assume that "
                     "-- report as pH-dependent rather than silently picking cationic, matching "
                     "this project's standing 'do not guess system_type' discipline.")],
        )

    # No ionizable/charged group matched at all.
    return ChargeClassificationResult(
        smiles=smiles, charge_type="nonionic", system_type_for_gibbs="nonionic",
        confidence="moderate",
        caveats=[("No anionic or cationic functional group was matched by this function's "
                 "known patterns. This is a claim by absence, not a positive structural "
                 "confirmation of a nonionic headgroup (e.g. PEG/glycoside) -- correct for "
                 "all real nonionic surfactants tested so far, but if this SMILES contains an "
                 "unusual/unrecognized ionizable group this function's pattern list does not "
                 "cover, it would be silently missed. Flagged as moderate, not high, "
                 "confidence for exactly this reason.")],
    )


# --- pH-conditional resolution of the "ambiguous_pH_dependent" case ------
#
# Real, buildable [COMPUTE] closure of part of the pH-ambiguous-charge
# bottleneck (see benchmark/paper3_groundzero/BOTTLENECK_RESOLUTION_PLAN.md
# item 14): classify_surfactant_charge_type above correctly refuses to
# resolve a free-amine surfactant's charge type when solution pH is
# unstated -- that refusal stays correct and is NOT removed here. What WAS
# a real, separate gap: even when a caller DOES supply a real solution pH,
# classify_surfactant_charge_type had no way to use it. This closes that
# gap via the standard Henderson-Hasselbalch relation applied to literature
# pKaH ranges for the three free-amine classes this module already detects.
#
# pKaH ranges (aqueous, conjugate acid of the free amine) sourced from
# real small-molecule reference values, live-verified 2026-09-15:
# methylamine (primary) pKaH=10.64, dimethylamine (secondary) pKaH=10.73,
# trimethylamine (tertiary) pKaH=9.79. Reported here as RANGES, not a
# single false-precise number, and explicitly disclosed as small-molecule
# analogues -- a long-chain alkylamine's pKa is well known to sit close to
# its small-molecule analogue (alkyl chain length beyond ~C3 has little
# further inductive effect), typically within a few tenths of a pH unit,
# but a real system-specific measured pKa should be preferred via the
# amine_pka_override parameter when precision matters.
_AMINE_PKAH_RANGES = {
    "primary_amine": (10.4, 10.8),
    "secondary_amine": (10.5, 11.0),
    "tertiary_amine": (9.5, 10.0),
}


def _henderson_hasselbalch_fraction_protonated(ph: float, pka: float) -> float:
    """Fraction of a basic (amine) group protonated (charged) at a given
    pH: f = 1 / (1 + 10^(pH - pKaH)). Local helper, no external state."""
    return 1.0 / (1.0 + 10.0 ** (ph - pka))


@dataclass
class PHConditionalChargeResult:
    smiles: str
    ph: float
    charge_type_at_ph: str  # "anionic" | "cationic" | "zwitterionic" | "nonionic" | "still_ambiguous"
    amine_class: str | None
    fraction_protonated_low: float | None
    fraction_protonated_high: float | None
    pka_source: str = ""
    confidence: str = ""
    caveats: list[str] = field(default_factory=list)


def classify_surfactant_charge_type_at_ph(
    smiles: str, ph: float, amine_pka_override: float | None = None
) -> PHConditionalChargeResult:
    """Resolve classify_surfactant_charge_type's "ambiguous_pH_dependent"
    case using a REAL, explicitly-supplied solution pH -- the pH itself is
    still never guessed (must be a real caller-supplied value, an
    experiment-specific fact this function cannot recover any other way),
    but once given, this computes real ionization behavior via Henderson-
    Hasselbalch instead of leaving a usable input unused.

    Only meaningful for structures classify_surfactant_charge_type itself
    reports as "ambiguous_pH_dependent" (a free primary/secondary/tertiary
    amine, alone or alongside a real anionic group) -- for any other
    charge_type, delegates straight back to the pH-independent result
    (a permanently-charged or nonionic group does not need pH information
    at all, and this function does not pretend otherwise).

    amine_pka_override: a real, system-specific measured pKa for THIS
    exact compound, if independently known -- when supplied, used exactly
    (single-point fraction_protonated, low==high). When not supplied,
    falls back to this module's own literature pKaH RANGE for the
    detected amine class (see _AMINE_PKAH_RANGES above) and reports a
    fraction_protonated RANGE (low, high) rather than a single false-
    precise number -- both bounds computed via Henderson-Hasselbalch.

    Still returns "still_ambiguous" (not a forced pick) when the computed
    protonation fraction straddles a genuinely ambiguous middle ground
    (0.3-0.7, i.e. within about 0.5 pH units of pKaH) -- a real physical
    statement that the population is genuinely mixed at that pH, not a
    resolvable single answer.
    """
    if ph < 0.0 or ph > 14.0:
        raise ValueError("ph must be between 0 and 14")
    base = classify_surfactant_charge_type(smiles)
    if base.charge_type != "ambiguous_pH_dependent":
        return PHConditionalChargeResult(
            smiles=smiles, ph=ph, charge_type_at_ph=base.charge_type,
            amine_class=None, fraction_protonated_low=None, fraction_protonated_high=None,
            pka_source="not applicable (charge type is pH-independent)",
            confidence=base.confidence,
            caveats=base.caveats + [("This structure's charge type does not depend on pH; "
                                     "the supplied ph value was not needed.")],
        )

    amine_names = base.cationic_conditional_groups_found
    if not amine_names:
        return PHConditionalChargeResult(
            smiles=smiles, ph=ph, charge_type_at_ph="still_ambiguous",
            amine_class=None, fraction_protonated_low=None, fraction_protonated_high=None,
            confidence="low",
            caveats=[("Base classification reported ambiguous_pH_dependent but no amine group "
                     "was recorded -- cannot compute a pKaH-based resolution.")],
        )
    amine_class = amine_names[0]  # the detected free-amine substitution class
    multi_amine_caveat: list[str] = []
    if len(amine_names) > 1:
        # A real, honest completeness gap, disclosed rather than silently
        # resolved: the base classifier's substructure match only reports
        # WHICH amine classes are present, not how many separate nitrogen
        # atoms of each -- a real polyamine surfactant can genuinely carry
        # more than one distinct amine class (e.g. one primary + one
        # tertiary nitrogen), each with its own pKaH and its own real
        # contribution to the molecule's net charge state at a given pH.
        # This function resolves using only the FIRST-listed class below
        # and does not attempt a combined/summed multi-site answer (that
        # would require per-site protonation states, not a single
        # charge_type label) -- lowered confidence and an explicit
        # caveat flag this rather than silently picking one site.
        multi_amine_caveat.append(
            f"Multiple distinct amine classes were detected ({amine_names}), not just one -- "
            f"this resolution uses only the first-listed class ('{amine_class}') and its own "
            "pKaH range/override. A real polyamine surfactant can have each site at a genuinely "
            "different protonation state at the same pH; treat this result as characterizing "
            "only that one site, not a verified net-charge answer for the whole molecule."
        )

    if amine_pka_override is not None:
        pka_low = pka_high = amine_pka_override
        pka_source = f"caller-supplied real pKa={amine_pka_override:.2f} for this specific compound"
    else:
        pka_low, pka_high = _AMINE_PKAH_RANGES.get(amine_class, (9.0, 11.0))
        pka_source = (f"literature pKaH range for {amine_class} class ({pka_low:.1f}-{pka_high:.1f}), "
                      "small-molecule analogue-based -- see this module's own docstring")

    # f (fraction protonated) is monotonically increasing in pKaH at fixed pH,
    # so pka_low gives the lower fraction-protonated bound and pka_high the upper.
    frac_low = _henderson_hasselbalch_fraction_protonated(ph, pka_low)
    frac_high = _henderson_hasselbalch_fraction_protonated(ph, pka_high)

    is_anionic_too = len(base.anionic_groups_found) > 0
    mostly_protonated = frac_low >= 0.7
    mostly_deprotonated = frac_high <= 0.3
    caveats = [f"Amine class: {amine_class}. pKaH source: {pka_source}."] + multi_amine_caveat

    if not (mostly_protonated or mostly_deprotonated):
        return PHConditionalChargeResult(
            smiles=smiles, ph=ph, charge_type_at_ph="still_ambiguous",
            amine_class=amine_class, fraction_protonated_low=frac_low, fraction_protonated_high=frac_high,
            pka_source=pka_source, confidence="low",
            caveats=caveats + [(f"At pH={ph}, the amine's computed protonation fraction spans "
                                f"{frac_low:.2f}-{frac_high:.2f} -- genuinely mixed/ambiguous population "
                                "this close to pKaH, not a resolvable single charge state.")],
        )

    if mostly_protonated:
        resolved = "zwitterionic" if is_anionic_too else "cationic"
        caveats.append(f"At pH={ph}, computed protonation fraction {frac_low:.2f}-{frac_high:.2f} "
                        "(mostly protonated/charged).")
    else:
        resolved = "anionic" if is_anionic_too else "nonionic"
        caveats.append(f"At pH={ph}, computed protonation fraction {frac_low:.2f}-{frac_high:.2f} "
                        "(mostly deprotonated/neutral amine).")

    base_confidence = "moderate" if amine_pka_override is None else "high"
    final_confidence = "low" if multi_amine_caveat else base_confidence

    return PHConditionalChargeResult(
        smiles=smiles, ph=ph, charge_type_at_ph=resolved,
        amine_class=amine_class, fraction_protonated_low=frac_low, fraction_protonated_high=frac_high,
        pka_source=pka_source, confidence=final_confidence,
        caveats=caveats,
    )


# --- structural-family classification (orthogonal to charge type) --------
#
# Real, previously-missing capability flagged by TOOLKIT_CAPABILITY_AUDIT.md:
# classify_surfactant_charge_type only determines charge type (anionic/
# cationic/zwitterionic/nonionic) -- it says nothing about whether a
# surfactant is a gemini/dimeric (two headgroups, one molecule), a glycolipid
# biosurfactant (rhamnolipid/sophorolipid-type), or an ordinary single-
# headgroup ("monomeric") surfactant. This closes that gap using the same
# real RDKit SMARTS substructure-matching discipline as the charge-type
# classifier above, verified against real, sourced structures before being
# treated as correct (not guessed):
#   - gemini/dimeric: G6 and 12-4-12 (real, literature-verified bis-
#     quaternary-ammonium gemini surfactants already used elsewhere in this
#     project -- see SurfQSPR's GEMINI_CANDIDATE_STRUCTURES) each show the
#     SAME strong ionic headgroup pattern matching exactly TWICE, vs exactly
#     ONCE for monomeric DTAB/CTAB -- confirmed directly, not assumed.
#   - glycolipid biosurfactant: real monorhamnolipid (PubChem CID 162246,
#     CAS 37134-61-5) and sophorolipid (PubChem CID 11856871) both match a
#     real pyranose-sugar-ring pattern combined with a long hydrocarbon
#     chain; every plain surfactant already in this project's own test
#     suite (SDS, DTAB, CTAB, an ethoxylate) shows zero sugar-ring matches,
#     and bare glucose (a sugar with NO lipid tail, not a surfactant at all)
#     is correctly excluded by also requiring a long chain.
#
# Honest limitation, disclosed rather than silently resolved: this
# heuristic (exactly 2 matching strong ionic headgroups in one molecule)
# cannot structurally distinguish a true gemini (two separate tail+head
# units joined near the headgroups by a short spacer) from a bolaform
# surfactant (a single long hydrophobic backbone with one headgroup at
# EACH end, no separate spacer/tail pair) -- both real structural classes
# produce the same "2 matching headgroups" signal. Bolaform surfactants are
# a real but comparatively rare structural class; this function reports
# "dimeric (gemini-type)" and discloses the bolaform ambiguity explicitly
# in caveats rather than silently picking one, matching this module's
# standing "do not guess, disclose the gap" discipline.

_SUGAR_RING_PATTERN = "[C;R1]1([OX2,OH])[C;R1]([OX2,OH])[C;R1][C;R1]([OX2,OH])[C;R1][O;R1]1"
_LONG_CHAIN_PATTERN = "[CH2][CH2][CH2][CH2][CH2][CH2]"  # >=6 contiguous CH2 -- crude but real long-tail signal

_COMPILED_SUGAR_RING = Chem.MolFromSmarts(_SUGAR_RING_PATTERN)
_COMPILED_LONG_CHAIN = Chem.MolFromSmarts(_LONG_CHAIN_PATTERN)


@dataclass
class StructuralFamilyResult:
    smiles: str
    structural_family: str  # "dimeric_gemini_type" | "glycolipid_biosurfactant" | "monomeric" | "unparseable"
    n_strong_ionic_headgroups: int  # count of the SAME strong charge-type pattern found, deduplicated across types
    n_sugar_rings: int
    n_long_chain_matches: int
    confidence: str = ""
    caveats: list[str] = field(default_factory=list)


def classify_surfactant_structural_family(smiles: str) -> StructuralFamilyResult:
    """Determine a surfactant's structural family (dimeric/gemini-type,
    glycolipid biosurfactant, or ordinary monomeric) from its real SMILES
    structure -- orthogonal to classify_surfactant_charge_type's charge-type
    axis (a gemini surfactant can be cationic, anionic, or nonionic; a
    glycolipid biosurfactant is typically anionic via its own fatty-acid
    carboxylate, already caught by classify_surfactant_charge_type
    separately). See this module's own comment block above for the real,
    sourced structures this was verified against before shipping.

    Counts the largest number of matches among the SAME strong ionic
    headgroup pattern (quaternary ammonium, quaternary ammonium in a ring,
    sulfate ester, sulfonate, phosphate ester -- the "strong" patterns
    already used by classify_surfactant_charge_type, reused here rather
    than duplicated) -- exactly 2 matches of the SAME pattern is this
    function's real, verified signal for a dimeric/gemini-type structure.
    A count of 3+ is flagged as low-confidence (a real oligomeric/trimeric
    surfactant, or a false-positive multi-match on an unusual structure)
    rather than silently forced into "dimeric".

    Separately checks for a real pyranose-sugar-ring pattern combined with
    a long hydrocarbon chain -- the real structural signature of a
    glycolipid biosurfactant (rhamnolipid/sophorolipid-type), distinct from
    a bare sugar (no lipid tail, not a surfactant) or an ordinary
    surfactant (no sugar ring at all).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return StructuralFamilyResult(
            smiles=smiles, structural_family="unparseable",
            n_strong_ionic_headgroups=0, n_sugar_rings=0, n_long_chain_matches=0,
            confidence="none",
            caveats=[f"RDKit could not parse this SMILES: '{smiles}'. Check for typos/invalid valence."],
        )

    strong_headgroup_counts = []
    for name, smarts, _strength in _ANIONIC_PATTERNS + _CATIONIC_STRONG_PATTERNS:
        if _strength_of(name) != "strong":
            continue
        patt = Chem.MolFromSmarts(smarts)
        if patt is not None:
            strong_headgroup_counts.append(len(mol.GetSubstructMatches(patt)))
    max_headgroup_count = max(strong_headgroup_counts) if strong_headgroup_counts else 0

    n_sugar_rings = len(mol.GetSubstructMatches(_COMPILED_SUGAR_RING)) if _COMPILED_SUGAR_RING is not None else 0
    n_long_chain = len(mol.GetSubstructMatches(_COMPILED_LONG_CHAIN)) if _COMPILED_LONG_CHAIN is not None else 0

    is_glycolipid = n_sugar_rings >= 1 and n_long_chain >= 1

    if is_glycolipid:
        return StructuralFamilyResult(
            smiles=smiles, structural_family="glycolipid_biosurfactant",
            n_strong_ionic_headgroups=max_headgroup_count, n_sugar_rings=n_sugar_rings,
            n_long_chain_matches=n_long_chain, confidence="high",
            caveats=[],
        )

    if max_headgroup_count == 2:
        return StructuralFamilyResult(
            smiles=smiles, structural_family="dimeric_gemini_type",
            n_strong_ionic_headgroups=2, n_sugar_rings=n_sugar_rings,
            n_long_chain_matches=n_long_chain, confidence="high",
            caveats=[("This structural signal (2 matching strong ionic headgroups in one molecule) "
                     "cannot distinguish a true gemini (two separate tail+head units joined near the "
                     "headgroups by a short spacer) from a bolaform surfactant (a single long "
                     "hydrophobic backbone with one headgroup at EACH end, no separate spacer/tail "
                     "pair) -- both produce the same signal. Bolaform surfactants are real but "
                     "comparatively rare; confirm against the actual connectivity if this distinction "
                     "matters for your use case.")],
        )

    if max_headgroup_count >= 3:
        return StructuralFamilyResult(
            smiles=smiles, structural_family="monomeric",
            n_strong_ionic_headgroups=max_headgroup_count, n_sugar_rings=n_sugar_rings,
            n_long_chain_matches=n_long_chain, confidence="low",
            caveats=[(f"{max_headgroup_count} matches of the same strong ionic headgroup pattern were "
                     "found -- not confidently classified as dimeric/gemini-type (which this function "
                     "only asserts at exactly 2) or as a genuinely different oligomeric structural "
                     "class; reported as monomeric by default but flagged low-confidence rather than "
                     "silently picked.")],
        )

    return StructuralFamilyResult(
        smiles=smiles, structural_family="monomeric",
        n_strong_ionic_headgroups=max_headgroup_count, n_sugar_rings=n_sugar_rings,
        n_long_chain_matches=n_long_chain, confidence="high" if max_headgroup_count <= 1 else "moderate",
        caveats=[],
    )


def _strength_of(pattern_name: str) -> str:
    for name, _smarts, strength in _ANIONIC_PATTERNS + _CATIONIC_STRONG_PATTERNS:
        if name == pattern_name:
            return strength
    return ""
