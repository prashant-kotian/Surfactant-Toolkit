"""EACN (Equivalent Alkane Carbon Number) from molecular structure.

Group-contribution rules with an explicit, cited rule table -- every rule
below is traceable to a published source, anything else is REFUSED rather
than guessed (a hallucinated EACN silently poisons every HLD calculation
downstream: optimal salinity, fish diagrams, phase inversion).

SOURCING PASS (2026-10-01, web):
  [Wade1977]  Wade et al. 1977 introduced EACN; by DEFINITION the EACN of
              an n-alkane equals its carbon number (ACN).  Reaffirmed in
              ACS Omega 2022 (doi 10.1021/acsomega.2c05334): "for the
              n-alkane series, the EACN is, by definition, equal to the
              carbon number ... denoted ACN".
  [Benjaminsen] "Reliable Equivalent Alkane Carbon Number Determination..."
              SPE J 2024: linear lnS* vs ACN interpolation; pure C6..C16
              have EACN = ACN; the linear mixing rule for blends is
              sufficient ("these data do not support a nonlinear mixing
              rule").
  [Benzoate2007] The EACN scale revisited (Queste et al. 2007) + HLD
              review quotes: "an alkyl benzene gave the same interfacial
              tension as an alkane having the same number of carbon atoms
              as the alkyl side chain" (benzene ring contributes 0), and
              "a cyclohexyl ring was approximately equivalent to four
              carbons"; the same sources assign butylcyclohexane EACN = 8
              (= 4 side-chain + 4 ring) and octylbenzene EACN = 8
              (= 8 side-chain + 0 ring).
  [Ester2013]  "Classification of ester oils according to their EACN"
              (J Colloid Interface Sci 2013) + COSMO-RS review quote:
              "an ester group decreases the EACN by about 10-12 units with
              respect to the total number of carbon atoms"; ester EACN
              also depends on ester POSITION and chain split, so we report
              a band, not a point truth.

Everything not covered by these four rules (branched alkanes, ethers,
ketones, alcohols, nitriles, polyaromatics, ring+ester combinations,
multi-ring systems, UNSATURATED chains -- alkenes/alkynes and
unsaturated esters, since EACN = carbon number is definitional only for
SATURATED n-alkanes and Ester2013's -10..-12 increment is sourced for
saturated chains, ...) raises ValueError with guidance: measure T* via
fish diagrams or use a COSMO-RS/ML estimator instead of guessing.

SATURATION GATE (added 2026-10-05, found while building the CS-A
application demo): the original implementation applied the n-alkane and
single-ester rules without checking bond saturation, so e.g. 1-octene
was silently returned as "n_alkane, EACN = 8" and methyl oleate as an
ester EACN -- a documented-rule violation that would quietly poison HLD
predictions. Any carbon-carbon double/triple bond outside an aromatic
ring is now refused on every path (the carbonyl C=O of the ester group
itself is NOT such a bond and stays allowed).
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class EacnResult:
    """Group-contribution EACN with band, rule provenance and caveats."""
    eacn: float
    eacn_low: float
    eacn_high: float
    rule: str                     # which table row fired
    n_carbons: int
    breakdown: list[str] = field(default_factory=list)
    sources: list[str] = field(default_factory=list)
    method: str = "group_contribution_rule_table"
    converged: bool = True
    n_iterations: int = 1
    params: dict = field(default_factory=dict)
    caveat: str = ""

    def as_dict(self) -> dict:
        return {
            "eacn": self.eacn, "eacn_low": self.eacn_low,
            "eacn_high": self.eacn_high, "rule": self.rule,
            "n_carbons": self.n_carbons, "breakdown": self.breakdown,
            "sources": self.sources, "method": self.method,
            "converged": self.converged, "n_iterations": self.n_iterations,
            "params": self.params, "caveat": self.caveat,
        }


_REFUSAL = (
    "no sourced group-contribution rule for this scaffold -- EACN rules "
    "cover ONLY (1) n-alkanes, (2) alkylbenzenes with one unbranched "
    "saturated side chain, (3) alkylcyclohexanes with one unbranched "
    "saturated side chain, (4) single-ester acyclic saturated chains. "
    "For {reason} do not guess: measure T* (fish diagram / Winsor III "
    "scan) against an n-alkane ladder, or use a COSMO-RS or ML EACN "
    "estimator. Refusing is the correct behavior here -- a fabricated "
    "EACN would corrupt optimal-salinity predictions downstream."
)


def _refuse(reason: str) -> ValueError:
    return ValueError(_REFUSAL.format(reason=reason))


def eacn_group_contribution(smiles: str) -> EacnResult:
    """EACN of an oil molecule from its SMILES via sourced group increments.

    Rules (see module docstring for full citations):
      n-alkane (acyclic, saturated, unbranched, C-only):
          EACN = carbon number, exact            [Wade1977 definitional]
      alkylbenzene (1 benzene ring + 1 unbranched saturated chain):
          EACN = side-chain carbons, exact       [Benzoate2007; octylbenzene=8]
      alkylcyclohexane (1 cyclohexane ring + 1 unbranched saturated chain):
          EACN = side-chain carbons + 4, exact   [Benzoate2007; butylcyclohexane=8]
      single ester, acyclic saturated chain:
          EACN = N_C - 11, band [N_C-12, N_C-10]  [Ester2013; COSMO-RS review]
      everything else: ValueError (measured-T* / COSMO-RS guidance).

    Blend EACNs: use eacn_linear_mixing() -- linear mole/volume mixing is
    experimentally sufficient [Benjaminsen].
    """
    try:
        from rdkit import Chem
    except ImportError as exc:  # pragma: no cover - env dependent
        raise ImportError(
            "eacn_group_contribution needs rdkit for SMILES parsing -- "
            "install the toolkit's core dependencies (pip install rdkit); "
            "it is a core dep, not the [calc] extra"
        ) from exc

    if not isinstance(smiles, str) or not smiles.strip():
        raise ValueError("smiles must be a non-empty SMILES string -- "
                         "pass the oil molecule, e.g. 'CCCCCCCC' for n-octane")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"invalid SMILES {smiles!r} -- rdkit could not parse "
                         "it; check the structure string, do not invent one")

    n_c = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
    if n_c == 0:
        raise _refuse("the molecule contains no carbon")

    ri = mol.GetRingInfo()
    n_rings = ri.NumRings()
    n_hetero = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() not in (1, 6))
    frags = len(Chem.GetMolFrags(mol))

    # ester count: C(=O)O-C  (acid/ester carbonyl with alkoxy oxygen)
    ester_matches = mol.GetSubstructMatches(
        Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CX4]"))
    n_ester = len(ester_matches)

    # branch points: ACYCLIC carbon bonded to >=3 other ACYCLIC carbons.
    # Ring neighbours do not count (the chain attachment point on a ring is
    # not a branch), and ring atoms themselves are never branch points.
    branch_points = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6 or atom.IsInRing():
            continue
        # >=3 carbon neighbours TOTAL (ring atoms included): the chain
        # attachment carbon on a ring has only 2 (ring + chain), while an
        # isopropyl CH has 3 (ring-or-chain + 2 methyls) -> flagged
        c_neighbors = [nb for nb in atom.GetNeighbors()
                       if nb.GetAtomicNum() == 6]
        if len(c_neighbors) >= 3:
            branch_points += 1
    # ring substitution sites: ring carbons carrying an acyclic carbon
    ring_substitution_sites = 0
    for atom in mol.GetAtoms():
        if atom.IsInRing() and any(
                nb.GetAtomicNum() == 6 and not nb.IsInRing()
                for nb in atom.GetNeighbors()):
            ring_substitution_sites += 1

    aromatic_ring_carbons = sum(
        1 for a in mol.GetAtoms() if a.IsInRing() and a.GetIsAromatic())

    # Saturation gate (see module docstring): a C-C double/triple bond
    # OUTSIDE an aromatic ring falls outside every sourced rule. The
    # ester's own carbonyl C=O (C-O bond) is allowed; kekulized/aromatic
    # ring bonds are allowed (Benzene/Queste rules cover those rings).
    c_unsaturated = any(
        b.GetBeginAtom().GetAtomicNum() == 6
        and b.GetEndAtom().GetAtomicNum() == 6
        and b.GetBondTypeAsDouble() > 1.0
        and not (b.IsInRing() and b.GetIsAromatic())
        for b in mol.GetBonds()
    )

    params = {"smiles": smiles, "n_carbons": n_c, "n_rings": n_rings,
              "n_ester": n_ester, "branch_points": branch_points,
              "n_fragments": frags, "n_heteroatoms": n_hetero,
              "c_unsaturated": c_unsaturated}

    def result(eacn, low, high, rule, breakdown, sources, caveat):
        return EacnResult(
            eacn=float(eacn), eacn_low=float(low), eacn_high=float(high),
            rule=rule, n_carbons=n_c, breakdown=breakdown, sources=sources,
            params=dict(params), caveat=caveat)

    if frags != 1:
        raise _refuse("the structure has multiple disconnected fragments")

    # ---- single-ester rule (checked before ring refusals below? no: ester
    # must be acyclic per Ester2013, so ring+ester is refused)
    if n_ester > 0:
        if n_ester > 1:
            raise _refuse("multiple ester groups (only the single-ester "
                          "increment -10..-12 per ester is sourced, and "
                          "polyesters were not validated)")
        if n_rings > 0:
            raise _refuse("ester + ring combination (the -10..-12 ester "
                          "increment is sourced only for acyclic chains)")
        if branch_points > 0:
            raise _refuse("branched ester (ester increments are sourced for "
                          "linear chains only; position matters)")
        if n_c < 12:
            raise _refuse("a short-chain ester (the -10..-12 increment is "
                          "sourced for long-chain oil tails, N_C >= 12; "
                          "small esters' EACNs vary strongly with chain "
                          "split and need a measured T*)")
        if c_unsaturated:
            raise _refuse("an UNSATURATED ester (C=C in the chain; "
                          "Ester2013's -10..-12 increment is sourced for "
                          "SATURATED single-ester chains only -- a double "
                          "bond's own EACN contribution was never sourced "
                          "here)")
        e = n_c - 11
        return result(
            e, n_c - 12, n_c - 10, "ester_single_acyclic",
            [f"total carbons N_C = {n_c}", "ester increment -11 (band -10..-12)"],
            ["Ester2013: J Colloid Interface Sci, ester oils EACN classification",
             ("COSMO-RS review: 'an ester group decreases the EACN by about "
             "10-12 units with respect to the total number of carbon atoms'")],
            "Ester EACN depends on ester POSITION and acid/alcohol chain "
            "split (Ester2013 shows both matter): the -11 increment with a "
            "2-unit band is the sourced central estimate, not a point "
            "truth. Validate against a T* measurement before quoting to "
            "3 significant figures.")

    # ---- hydrocarbon rules: no heteroatoms allowed
    if n_hetero > 0:
        raise _refuse(f"a heteroatom-containing oil (n_heteroatoms={n_hetero}; "
                      "only hydrocarbons and single esters are sourced)")

    if c_unsaturated:
        raise _refuse("an UNSATURATED hydrocarbon (alkene/alkyne, or a "
                      "double bond in a non-aromatic ring: the EACN = "
                      "carbon-number identity is DEFINITIONAL only for "
                      "SATURATED n-alkanes, and no unsaturated increment "
                      "is sourced for the alkylbenzene/alkylcyclohexane "
                      "rules either)")

    if n_rings == 0:
        # acyclic hydrocarbon: must be unbranched alkane
        if branch_points > 0:
            raise _refuse("a BRANCHED alkane (branched-alkane EACNs were "
                          "assigned only case-by-case by T* measurement -- "
                          "Queste 2007 -- with no table we can cite)")
        return result(
            n_c, n_c, n_c, "n_alkane",
            [f"carbon number = {n_c} (EACN = ACN by definition)"],
            [("Wade1977 (EACN definition); ACS Omega 2022: 'for the n-alkane "
             "series, the EACN is, by definition, equal to the carbon number'")],
            "Exact by definition for unbranched n-alkanes. Do not extend to "
            "iso-alkanes: their EACN must be measured (T* fish diagram).")

    if n_rings == 1 and aromatic_ring_carbons == 6:
        # alkylbenzene: ring = 0, EACN = unbranched side chain carbons
        if branch_points > 0:
            raise _refuse("a branched side chain on the alkylbenzene (only "
                          "unbranched chains are sourced)")
        if ring_substitution_sites != 1:
            raise _refuse("a benzene ring with != 1 alkyl substituent "
                          "(the side-chain rule is sourced for "
                          "mono-alkylbenzenes only)")
        side = n_c - 6
        if side < 1:
            raise _refuse("benzene with no alkyl side chain (the "
                          "side-chain rule needs >=1 chain carbon)")
        return result(
            side, side, side, "alkylbenzene_sidechain",
            [f"side chain carbons = {side}", "benzene ring contribution = 0"],
            [("Benzoate2007/HLD review: 'an alkyl benzene gave the same "
             "interfacial tension as an alkane having the same number of "
             "carbon atoms as the alkyl side chain'; octylbenzene EACN = 8")],
            "Ring-contributes-zero holds for mono-alkylbenzenes with "
            "unbranched saturated side chains; poly-substituted or "
            "polyaromatic oils need measurement.")

    if n_rings == 1 and aromatic_ring_carbons == 0:
        ring_size = len(ri.AtomRings()[0])
        if ring_size == 6:
            # saturated 6-ring = cyclohexane: ring contribution +4
            if branch_points > 0:
                raise _refuse("a branched side chain on the "
                              "alkylcyclohexane (only unbranched chains "
                              "are sourced)")
            if ring_substitution_sites != 1:
                raise _refuse("a cyclohexane ring with != 1 alkyl "
                              "substituent (mono-alkyl only)")
            side = n_c - 6
            if side < 1:
                raise _refuse("cyclohexane with no alkyl side chain (the "
                              "side-chain + ring-4 rule needs >=1 chain "
                              "carbon; bare cyclohexane's EACN should be "
                              "measured)")
            return result(
                side + 4, side + 4, side + 4, "alkylcyclohexane_sidechain_plus4",
                [f"side chain carbons = {side}",
                 "cyclohexyl ring contribution = 4"],
                [("Benzoate2007/HLD review: 'a cyclohexyl ring was "
                 "approximately equivalent to four carbons'; "
                 "butylcyclohexane EACN = 8 (= 4 + 4)")],
                "Ring = +4 comes from the butylcyclohexane = 8 assignment "
                "plus the four-carbon rule; treat other ring sizes as "
                "uncovered (refused).")
        raise _refuse(f"a {ring_size}-membered non-aromatic ring (only the "
                      "cyclohexane +4 increment is sourced)")

    raise _refuse("multiple or partially aromatic rings (only mono-"
                  "benzene and mono-cyclohexane are sourced)")


def eacn_linear_mixing(fractions: list[float],
                       eacns: list[float],
                       basis: str = "mole") -> dict:
    """Linear EACN mixing rule for oil blends: EACN = sum x_i * EACN_i.

    Linear mole- or volume-fraction mixing is experimentally sufficient
    for alkane dilutions [Benjaminsen, SPE J 2024: "the linear mixing
    rule ... is sufficient to estimate the diluted oil EACN; these data
    do not support a nonlinear mixing rule"].
    """
    if len(fractions) != len(eacns) or not fractions:
        raise ValueError("fractions and eacns must be non-empty lists of "
                         "equal length")
    if basis not in ("mole", "volume"):
        raise ValueError("basis must be 'mole' or 'volume' -- pick the "
                         "fraction type you actually measured")
    if any(f < 0 for f in fractions):
        raise ValueError("fractions must be non-negative")
    total = float(sum(fractions))
    if total <= 0:
        raise ValueError("fractions must sum to a positive value")
    e = sum(f * v for f, v in zip(fractions, eacns)) / total
    return {
        "eacn": float(e),
        "basis": basis,
        "method": "linear_volume_or_mole_weighting",
        "converged": True,
        "n_iterations": 1,
        "params": {"fractions": list(map(float, fractions)),
                   "eacns": list(map(float, eacns))},
        "caveat": "Linear mixing validated for n-alkane dilutions "
                  "[Benjaminsen 2024]; for chemically dissimilar blends "
                  "(e.g. aromatics into paraffins) confirm with a T* "
                  "measurement before quoting.",
    }
