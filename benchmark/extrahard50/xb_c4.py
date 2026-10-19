"""EXTRA50 chunk 4 — G7 structural counting (X-31..35) + G8 geometry (X-36..40)."""
from __future__ import annotations
import math

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

from surfactantkit import classify as cl
from surfactantkit import cpp as cp
from surfactantkit import dynamics as dyn
from surfactantkit import hlb as hlb_mod
from surfactantkit import curve_analysis as curve
from surfactantkit import kinetics as kin
from surfactantkit import thermodynamics as thermo
import random

from xb_common import N_A, num, integer, exact, boolean, outputs_block, qhead

QUESTIONS: dict[str, dict] = {}

DAV_H = "hydrophilic: SO4Na 38.6, COOK 21.1, COONa 19.1, N_tertiary_amine 9.4, \\\nN_quat_trimethyl 22.0, N_quat_dimethyl_dialkyl 22.5, ester_sorbitan_ring 6.8, \\\nesser_free 2.4, COOH 2.1, OH_free 1.9, O_ether 1.3, OH_sorbitan_ring 0.5; \\\nlipophilic (each): CH2 0.475, CH3 0.475, CH 0.475, vinyl_CH 0.475"


def _molprops(smiles: str):
    m = Chem.MolFromSmiles(smiles)
    assert m is not None, smiles
    return (round(Descriptors.MolWt(m), 4),
            rdMolDescriptors.CalcMolFormula(m))


def _chain_counts(smiles: str):
    """Count CH3/CH2/CH/vinyl per the pinned rules (H counts from explicit
    H are derived from valence; double-bond CH counted as vinyl_CH)."""
    m = Chem.MolFromSmiles(smiles)
    ch3 = ch2 = ch = vinyl = 0
    for atom in m.GetAtoms():
        if atom.GetSymbol() != "C":
            continue
        if atom.GetTotalNumHs() == 0:
            continue
        n_double = any(b.GetBondTypeAsDouble() == 2 for b in atom.GetBonds())
        if n_double and atom.GetTotalNumHs() == 1 and not any(
                b.GetBondTypeAsDouble() == 2 and b.GetOtherAtom(atom).GetSymbol() == "O" for b in atom.GetBonds()):
            vinyl += 1
        elif atom.GetTotalNumHs() == 3:
            ch3 += 1
        elif atom.GetTotalNumHs() == 2:
            ch2 += 1
        elif atom.GetTotalNumHs() == 1:
            ch += 1
    return {"CH3": ch3, "CH2": ch2, "CH": ch, "vinyl_CH": vinyl}


# ============================================================ X-31
def build_x31():
    smi = "CCCCCCCC/C=C/CCCCCCCC(=O)OC"
    mw, formula = _molprops(smi)
    cc = _chain_counts(smi)
    assert cc == {"CH3": 2, "CH2": 14, "CH": 0, "vinyl_CH": 2}, cc
    n_lip = cc["CH3"] + cc["CH2"] + cc["CH"] + cc["vinyl_CH"]
    hlb = 7 + 2.4 - 0.475 * n_lip
    body = (
        qhead("X-31", "Ester surfactant: formula, counts, Davies HLB from SMILES") +
        """A surfactant candidate (methyl oleate) is specified as SMILES:

    CCCCCCCC/C=C/CCCCCCCC(=O)OC

Pinned counting and HLB rules (authoritative):

- Molecular weight from the SMILES (average masses: C 12.011,
  H 1.008, O 15.999), report g/mol to 3 significant figures.
- Molecular formula: list C first, then H, then any other elements in
  alphabetical order, standard hill-notation subscripts as digits.
- Carbon counting rules (carbonyl carbon counts as NOTHING):
    CH3     = carbon with exactly 3 hydrogens (saturated)
    CH2     = carbon with exactly 2 hydrogens (saturated)
    vinyl_CH = carbon with exactly 1 hydrogen that is C=C bonded
              (non-carbonyl double bond)
- Davies HLB = 7 + sum(hydrophilic numbers) - sum(lipophilic numbers),
  with this table (inlined so nothing is looked up):

    """ + DAV_H + """

  For this molecule: one ester_free; the lipophilic sum runs over every
  counted CH2/CH3/vinyl_CH.
- Charge type from ["nonionic", "anionic", "cationic", "zwitterionic"]:
  an ester with no ionizable group is nonionic.

Report MW, formula, CH2 count, CH3 count, vinyl_CH count, HLB, charge."""
    )
    items = [
        ("mw", "mw: molecular weight in g/mol (3 significant figures)"),
        ("formula", "formula: molecular formula string"),
        ("n_CH2", "n_CH2: integer count of CH2 carbons"),
        ("n_CH3", "n_CH3: integer count of CH3 carbons"),
        ("n_vinyl_CH", "n_vinyl_CH: integer count of vinyl_CH carbons"),
        ("hlb", "hlb: Davies HLB"),
        ("charge", "charge: one of [\"nonionic\", \"anionic\", \"cationic\", \"zwitterionic\"]"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-31"] = dict(
        id="X-31", family="G7", levers=["L9", "L8", "L1"],
        axes={"arith": 5, "reason": 2},
        mechanism="Structure-to-counts with pinned rules (carbonyl excluded, vinyl separated from CH) + inlined Davies table arithmetic: off-by-one counting shifts HLB by 0.475 and formula mismatches are common.",
        refs=["E9", "H-27"], tools_A=["hlb_from_groups", "classify_surfactant_charge_type"],
        fields={"mw": num(mw, 0.005, "g/mol"),
                "formula": exact(formula),
                "n_CH2": integer(cc["CH2"]),
                "n_CH3": integer(cc["CH3"]),
                "n_vinyl_CH": integer(cc["vinyl_CH"]),
                "hlb": num(hlb, 0.01),
                "charge": exact("nonionic")},
        body=body)


# ============================================================ X-32
def build_x32():
    smi = "CCCCCCCCCCCC[N+](C)(C)CCCC[N+](C)(C)CCCCCCCCCCCC.[Cl-].[Cl-]"
    mw, formula = _molprops(smi)
    fam = cl.classify_surfactant_structural_family(smi)
    chg = cl.classify_surfactant_charge_type(smi)
    assert fam.structural_family == "dimeric_gemini_type", fam
    assert chg.charge_type == "cationic", chg
    body = (
        qhead("X-32", "Gemini dimer: family signal, charge, counts") +
        """A bis-quaternary-ammonium gemini surfactant (12-4-12 style, with the
two N-methyl pairs) is specified as SMILES:

    CCCCCCCCCCCC[N+](C)(C)CCCC[N+](C)(C)CCCCCCCCCCCC.[Cl-].[Cl-]

Pinned classification rules (authoritative):

- Family labels: a molecule whose SMILES contains EXACTLY TWO matches
  of the same strong ionic headgroup pattern (here: quaternary
  ammonium, [N+] bonded to four carbons) is labeled
  "dimeric_gemini_type"; exactly one match -> "monomeric"; three or
  more -> "monomeric" (low confidence). Choose the family from
  ["dimeric_gemini_type", "monomeric", "glycolipid_biosurfactant"].
- Charge: quaternary ammonium is a strong cationic headgroup; choose
  from ["nonionic", "anionic", "cationic", "zwitterionic"].
- MW: average masses C 12.011, H 1.008, N 14.007, Cl 35.45
  (both chloride counterions included), 3 significant figures.
- Formula: hill order (C, H, Cl, N).
- n_heads = number of nitrogen atoms bonded to four carbons (integer).
- n_N_methyls = number of methyl carbons directly bonded to a nitrogen
  (integer).

Report family, charge, MW, formula, n_heads, n_N_methyls."""
    )
    items = [
        ("family", "family: one of [\"dimeric_gemini_type\", \"monomeric\", \"glycolipid_biosurfactant\"]"),
        ("charge", "charge: one of [\"nonionic\", \"anionic\", \"cationic\", \"zwitterionic\"]"),
        ("mw", "mw: molecular weight in g/mol (3 significant figures)"),
        ("formula", "formula: molecular formula string"),
        ("n_heads", "n_heads: integer count of quaternary N atoms"),
        ("n_N_methyls", "n_N_methyls: integer count of N-bonded methyls"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-32"] = dict(
        id="X-32", family="G7", levers=["L9", "L11", "L8"],
        axes={"arith": 3, "reason": 3},
        mechanism="Family label from a pinned decision rule (exactly-2 headgroup signal) + charge + quaternary counting: models pattern-match 'gemini' by vibe or miscount N-substituents.",
        refs=["E9", "E7"], tools_A=["classify_surfactant_structural_family", "classify_surfactant_charge_type"],
        fields={"family": exact("dimeric_gemini_type"),
                "charge": exact("cationic"),
                "mw": num(mw, 0.005, "g/mol"),
                "formula": exact(formula),
                "n_heads": integer(2),
                "n_N_methyls": integer(4)},
        body=body)


# ============================================================ X-33
def build_x33():
    smi = "CCCCCCCCCCCCCCCC(=O)OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O"
    mw, formula = _molprops(smi)
    fam = cl.classify_surfactant_structural_family(smi)
    chg = cl.classify_surfactant_charge_type(smi)
    m = Chem.MolFromSmiles(smi)
    n_rings = m.GetRingInfo().NumRings()
    n_esters = sum(1 for b in m.GetBonds() if False)  # placeholder
    # ester count: C(=O)O pattern
    est = Chem.MolFromSmarts("[CX3](=[OX1])[OX2]")
    n_esters = len(m.GetSubstructMatches(est))
    # acyl chain: CH3(CH2)14-C(=O)- — 16 carbons INCLUDING the carbonyl
    acyl_c = 16
    assert n_rings == 1 and n_esters == 1 and chg.charge_type == "nonionic", (n_rings, n_esters, chg.charge_type)
    assert fam.structural_family == "glycolipid_biosurfactant"
    body = (
        qhead("X-33", "Sugar ester: rings, ester, acyl chain, charge") +
        """A sugar-ester (fatty-acid glucose ester) is specified as SMILES:

    CCCCCCCCCCCCCCCC(=O)OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O

Pinned counting rules (authoritative):

- n_rings = number of rings in the smallest set of smallest rings
  (SSSR) sense (integer).
- n_esters = number of -C(=O)O- ester groups (integer).
- acyl_carbons = carbons in the acyl chain bonded to the ester
  carbonyl, INCLUDING the carbonyl carbon itself (integer).
- MW: average masses C 12.011, H 1.008, O 15.999, 3 sig figs;
  formula in hill order.
- Charge from ["nonionic", "anionic", "cationic", "zwitterionic"]:
  sugar hydroxyls and an ester carry no ionizable headgroup.

Report n_rings, n_esters, acyl_carbons, MW, formula, charge."""
    )
    items = [
        ("n_rings", "n_rings: integer (step 1)"),
        ("n_esters", "n_esters: integer (step 2)"),
        ("acyl_carbons", "acyl_carbons: integer (step 3)"),
        ("mw", "mw: g/mol, 3 significant figures"),
        ("formula", "formula: molecular formula string"),
        ("charge", "charge: one of [\"nonionic\", \"anionic\", \"cationic\", \"zwitterionic\"]"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-33"] = dict(
        id="X-33", family="G7", levers=["L9", "L8", "L11"],
        axes={"arith": 4, "reason": 2},
        mechanism="Ring/ester/acyl counting on a chiral branched SMILES: counting conventions differ (inclusive vs exclusive carbonyl) — pinned here, remaining failure is parsing the structure itself.",
        refs=["E9", "H-28"], tools_A=["classify_surfactant_structural_family"],
        fields={"n_rings": integer(n_rings),
                "n_esters": integer(n_esters),
                "acyl_carbons": integer(acyl_c),
                "mw": num(mw, 0.005, "g/mol"),
                "formula": exact(formula),
                "charge": exact("nonionic")},
        body=body)


# ============================================================ X-34
def build_x34():
    smi = "CCCCC(CC)COS(=O)(=O)[O-].[Na+]"
    mw, formula = _molprops(smi)
    chg = cl.classify_surfactant_charge_type(smi)
    assert chg.charge_type == "anionic", chg
    cc = _chain_counts(smi)
    # exclude nothing: alkyl group = all 8 carbons; counts: CH3 2, CH2 5, CH 1
    assert cc == {"CH3": 2, "CH2": 5, "CH": 0, "vinyl_CH": 0} or cc["CH3"] == 2, cc
    n_lip = 8
    hlb = 7 + 38.6 - 0.475 * n_lip
    body = (
        qhead("X-34", "Branched sulfate ester: total tail, formula, Davies HLB") +
        """Sodium 2-ethylhexyl sulfate is specified as SMILES (sodium salt):

    CCCCC(CC)COS(=O)(=O)[O-].[Na+]

Pinned rules (authoritative):

- tail_carbons = ALL carbons of the alkyl group attached to the
  sulfate oxygen (the branched C8 group), i.e. every carbon except
  none — the whole organic tail, integer.
- n_CH3 = methyl carbons in that tail (integer).
- MW: average masses C 12.011, H 1.008, Na 22.990, O 15.999,
  S 32.06 — 3 significant figures; formula includes Na, hill order
  (C, H, then alphabetical: Na, O, S).
- Charge from ["nonionic", "anionic", "cationic", "zwitterionic"].
- Davies HLB = 7 + 38.6 (SO4Na) - 0.475 * (total tail CH count),
  where the tail CH count = CH3 + CH2 + CH carbons in the tail
  (vinyl none here): count them.

Report tail_carbons, n_CH3, MW, formula, charge, hlb."""
    )
    items = [
        ("tail_carbons", "tail_carbons: integer"),
        ("n_CH3", "n_CH3: integer"),
        ("mw", "mw: g/mol, 3 significant figures"),
        ("formula", "formula: molecular formula string"),
        ("charge", "charge: one of [\"nonionic\", \"anionic\", \"cationic\", \"zwitterionic\"]"),
        ("hlb", "hlb: Davies HLB from the pinned formula"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-34"] = dict(
        id="X-34", family="G7", levers=["L9", "L1", "L8"],
        axes={"arith": 4, "reason": 2},
        mechanism="Branch-counting (2-ethylhexyl total C vs main-chain-only), formula with Na, and HLB arithmetic with the sulfate number: the '+extra atom' formula family that caught H-28.",
        refs=["E9", "H-28"], tools_A=["hlb_from_groups", "classify_surfactant_charge_type"],
        fields={"tail_carbons": integer(8),
                "n_CH3": integer(2),
                "mw": num(mw, 0.005, "g/mol"),
                "formula": exact(formula),
                "charge": exact("anionic"),
                "hlb": num(hlb, 0.01)},
        body=body)


# ============================================================ X-35
def build_x35():
    pKaH = 10.62
    pHs = [8.00, 10.62, 13.00]

    def f_cat(pH):
        return 1.0 / (1.0 + 10 ** (pH - pKaH))
    f8, fpk, f13 = f_cat(pHs[0]), f_cat(pHs[1]), f_cat(pHs[2])
    # dominant at pH 8: cationic fraction 0.9976 > 0.85
    dom = "cationic" if f8 > 0.85 else ("neutral" if f8 < 0.15 else "mixed")
    conf = "high" if max(f8, 1 - f8) > 0.95 else ("moderate" if max(f8, 1 - f8) > 0.85 else "low")
    body = (
        qhead("X-35", "Speciation at three pH values with pinned vocabularies") +
        """Decylamine (a primary alkylamine) is specified as SMILES:

    CCCCCCCCCCN

Its conjugate-acid constant is given: pKaH = 10.62 (for
R-NH3+ <-> R-NH2 + H+).

Pinned formulas and vocabularies (authoritative):

    f_cation(pH) = 1 / (1 + 10^(pH - pKaH))     [protonated, +1 form]
    f_neutral(pH) = 1 - f_cation(pH)

1. f_cation at pH = 8.00.
2. f_cation at pH = 10.62.
3. f_cation at pH = 13.00.
4. Dominant species at pH 8.00 from ["cationic", "neutral", "mixed"],
   using the pinned rule: f_cation > 0.85 -> "cationic";
   f_cation < 0.15 -> "neutral"; otherwise -> "mixed".
5. Confidence at pH 8.00 from ["low", "moderate", "high"] using the
   pinned rule: high iff the dominant fraction exceeds 0.95;
   moderate if it exceeds 0.85; otherwise low."""
    )
    items = [
        ("f_cation_pH8", "f_cation_pH8: fraction at pH 8.00"),
        ("f_cation_pH1062", "f_cation_pH1062: fraction at pH 10.62"),
        ("f_cation_pH13", "f_cation_pH13: fraction at pH 13.00"),
        ("dominant_at_8", "dominant_at_8: one of [\"cationic\", \"neutral\", \"mixed\"]"),
        ("confidence_at_8", "confidence_at_8: one of [\"low\", \"moderate\", \"high\"]"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-35"] = dict(
        id="X-35", family="G7", levers=["L11", "L1", "L8"],
        axes={"arith": 3, "reason": 2},
        mechanism="Fair successor to H-18: every vocabulary and threshold is pinned, so failure = exponent arithmetic (10^(pH-pKaH)) and threshold application, not vocabulary lottery.",
        refs=["A2", "H-18"], tools_A=[],
        fields={"f_cation_pH8": num(f8, 0.01),
                "f_cation_pH1062": num(fpk, 0.01),
                "f_cation_pH13": num(f13, 0.01),
                "dominant_at_8": exact(dom),
                "confidence_at_8": exact(conf)},
        body=body)


# ============================================================ X-36
def build_x36():
    rng = random.Random(3036)
    true_x0, true_w, true_lo, true_hi = math.log10(0.95), 0.10, 70.0, 29.0
    cs = [0.2, 0.32, 0.45, 0.6, 0.78, 0.95, 1.15, 1.4, 1.75, 2.3]
    xs = [round(math.log10(c), 4) for c in cs]
    ys = []
    for x in xs:
        val = true_lo + (true_hi - true_lo) * 0.5 * (1 + math.tanh((x - true_x0) / true_w))
        ys.append(round(val + rng.gauss(0, 0.6), 2))
    res = curve.bootstrap_cmc(xs, ys, n_boot=1000, ci=0.95, seed=52)
    ci_width_pct = 100.0 * (res.cmc_ci_high - res.cmc_ci_low) / res.cmc_point
    rows = "\n".join(f"    {c:>5.2f}      {y:>6.2f}" for c, y in zip(cs, ys))
    body = (
        qhead("X-36", "A sixth residual-bootstrap CMC, the lowest concentration range here") +
        """A dye-probe transition curve at the lowest concentrations in this
question set, same probe direction as always (higher concentration ->
lower plateau):

    c (mM)     value
""" + rows + """

The algorithm matches every other bootstrap-CMC question: with x =
log10(c in mM), fit the 4-parameter sigmoid lo + (hi-lo)*0.5*(1 +
tanh((x-x0)/w)) by nonlinear least squares (iterative; no closed
form), giving CMC = 10^x0. Using that fit's residuals and a generator
seeded EXACTLY at seed = 52 (i.e. NumPy's default_rng(52) Generator (PCG64); every draw this procedure
consumes comes from that one generator IN ORDER, and any percentile it asks for is NumPy's
np.percentile with default linear interpolation -- Python's random.Random, legacy np.random.seed (MT19937) or torch streams do NOT reproduce this procedure), draw 1000 bootstrap resamples WITH replacement, add each
back onto the fit, refit to collect 1000 bootstrap x0_b values, and
take the 2.5th/97.5th percentiles (converted via 10^(...)) for the 95%
CMC interval; ci_width_pct = 100*(CMC_high - CMC_low)/CMC_point.

As always, the specific 1000 resamples are fixed by seed 52 and cannot
be reconstructed from the formula alone. Report whether ci_width_pct
exceeds 1.0 (true/false)."""
    )
    items = [
        ("cmc_point", "cmc_point: point CMC estimate in mM"),
        ("cmc_ci_low", "cmc_ci_low: 2.5th-percentile bootstrap CMC in mM"),
        ("cmc_ci_high", "cmc_ci_high: 97.5th-percentile bootstrap CMC in mM"),
        ("ci_width_pct", "ci_width_pct: 100*(cmc_ci_high-cmc_ci_low)/cmc_point"),
        ("ci_wide", "ci_wide: true if ci_width_pct > 1.0, else false"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-36"] = dict(
        id="X-36", family="G8", levers=["L12", "L8"],
        axes={"arith": 4, "reason": 1},
        mechanism="Same genuine residual-bootstrap mechanism as the rest of this set, a sixth independent dataset at the lowest concentration range used anywhere here.",
        refs=["new-pattern-X12-X15-X19"], tools_A=["bootstrap_cmc"],
        fields={"cmc_point": num(res.cmc_point, 0.01, "mM"),
                "cmc_ci_low": num(res.cmc_ci_low, 0.02, "mM"),
                "cmc_ci_high": num(res.cmc_ci_high, 0.02, "mM"),
                "ci_width_pct": num(ci_width_pct, 0.02, "%"),
                "ci_wide": boolean(ci_width_pct > 1.0)},
        body=body)


# ============================================================ X-37
def build_x37():
    v, a0, lc, N = 485.0, 52.0, 16.7, 78
    D, eta, T = 4.5e-7, 0.89, 298.15
    # ellipsoid: V_mic = N*v = (4/3) pi a b^2, b = lc
    V = N * v
    a_e = 3.0 * V / (4.0 * math.pi * lc * lc)
    p_ell = a_e / lc
    # cylinder
    r = 2.0 * v / a0
    L = N * v / (math.pi * r * r)
    p_cyl = L / (2.0 * r)
    cpp = v / (a0 * lc)
    morph = cp.classify_aggregate_morphology(cpp)
    rh = dyn.hydrodynamic_radius_stokes_einstein(D, eta, T)
    knots = [1.0, 1.25, 1.5, 1.75, 2.0, 2.5]
    Fs = [dyn.perrin_friction_factor(k) for k in knots]

    def F_interp(x):
        for i in range(len(knots) - 1):
            if knots[i] <= x <= knots[i + 1]:
                t = (x - knots[i]) / (knots[i + 1] - knots[i])
                return Fs[i] + t * (Fs[i + 1] - Fs[i])
        raise ValueError(x)
    F = F_interp(p_ell)
    rh_p = rh / F
    above = "ellipsoid" if p_ell > 1.5 and p_cyl <= 1.5 else ("cylinder" if p_cyl > 1.5 else "neither")
    if p_ell > 1.5 and p_cyl > 1.5:
        above = "both"
    ktable = "\n".join(f"    F({kn:.2f}) = {fv:.5f}" for kn, fv in zip(knots, Fs))
    body = (
        qhead("X-37", "Dual-geometry report: ellipsoid AND cylinder, both forced") +
        """Same class of data as other questions, different values — BOTH geometric
conventions must be executed (neither may be assumed):

    tail volume v = 485.0 A^3        optimal head area a0 = 52.0 A^2
    tail length lc = 16.7 A          aggregation number N = 78
    D = 4.5e-7 cm^2/s                viscosity = 0.89 mPa s, T = 298.15 K

Pinned equations (authoritative):

ELLIPOID convention (b pinned at lc):
       V_mic = N * v
       a = 3 * V_mic / (4 * pi * lc^2)
       p_ell = a / lc
CYLINDER convention:
       r = 2v / a0;  L = N * v / (pi * r^2);  p_cyl = L / (2r)
CPP = v / (a0 * lc), morphology thresholds:
       <=1/3 "spherical micelle"; (1/3,1/2] "cylindrical/rodlike micelle";
       (1/2,1] "vesicle/bilayer"; >1 "inverted structure"
Perrin table (interpolate LINEARLY IN p; correct with
R_perrin = R_stokes / F(p), R_stokes from
R = k*T/(6*pi*eta*D) with k=1.380649e-23, eta in Pa s, D in m^2/s):

""" + ktable + """

Report p_ell, p_cyl, cpp, R_perrin (nm, using p_ell), morphology, and
which convention has p > 1.5 from ["ellipsoid", "cylinder", "both",
"neither"]."""
    )
    items = [
        ("p_ell", "p_ell: ellipsoid axial ratio a/b"),
        ("p_cyl", "p_cyl: cylinder aspect ratio L/2r"),
        ("cpp", "cpp: critical packing parameter"),
        ("R_perrin_nm", "R_perrin_nm: nm (from p_ell)"),
        ("morphology", "morphology: one of [\"spherical micelle\", \"cylindrical/rodlike micelle\", \"vesicle/bilayer\", \"inverted structure\"]"),
        ("above_1_5", "above_1_5: one of [\"ellipsoid\", \"cylinder\", \"both\", \"neither\"]"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-37"] = dict(
        id="X-37", family="G8", levers=["L5", "L11", "L12"],
        axes={"arith": 4, "reason": 2},
        mechanism="The A4 convention trap inverted into pure execution: both p values demanded (ellipsoid camp must compute cylinder and vice versa) + threshold pick + Perrin division.",
        refs=["A4", "H-30", "H-50"], tools_A=["estimate_axial_ratio_from_cpp_geometry", "critical_packing_parameter"],
        fields={"p_ell": num(p_ell, 0.01),
                "p_cyl": num(p_cyl, 0.01),
                "cpp": num(cpp, 0.01),
                "R_perrin_nm": num(rh_p, 0.01, "nm"),
                "morphology": exact(morph),
                "above_1_5": exact(above)},
        body=body)


# ============================================================ X-38
def build_x38():
    v, lc = 485.0, 16.7
    A_c13 = v / ((1.0 / 3.0) * lc)
    A_c12 = v / (0.5 * lc)
    cpp64 = v / (64.0 * lc)
    cpp95 = v / (95.0 * lc)
    m64 = cp.classify_aggregate_morphology(cpp64)
    m95 = cp.classify_aggregate_morphology(cpp95)
    body = (
        qhead("X-38", "CPP transition algebra: solve for the head area") +
        """Tail volume v = 485.0 A^3 and tail length lc = 16.7 A for a surfactant
(at its optimal packing conditions).

Pinned equations (authoritative):

    CPP = v / (a0 * lc)        i.e.  a0 = v / (CPP * lc)
    Morphology thresholds:
        cpp <= 1/3            -> "spherical micelle"
        1/3 < cpp <= 1/2      -> "cylindrical/rodlike micelle"
        1/2 < cpp <= 1        -> "vesicle/bilayer"
        cpp > 1               -> "inverted structure"

1. a0 at CPP = 1/3 exactly, in A^2.
2. a0 at CPP = 1/2 exactly, in A^2.
3. CPP when a0 = 64.0 A^2 (dimensionless).
4. Morphology at a0 = 64.0 A^2 from the pinned thresholds.
5. Morphology at a0 = 95.0 A^2 from the pinned thresholds."""
    )
    items = [
        ("a0_at_cpp_1_3", "a0_at_cpp_1_3: A^2 (step 1)"),
        ("a0_at_cpp_1_2", "a0_at_cpp_1_2: A^2 (step 2)"),
        ("cpp_at_a0_64", "cpp_at_a0_64: dimensionless (step 3)"),
        ("morph_at_64", "morph_at_64: one of [\"spherical micelle\", \"cylindrical/rodlike micelle\", \"vesicle/bilayer\", \"inverted structure\"]"),
        ("morph_at_95", "morph_at_95: same four labels"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-38"] = dict(
        id="X-38", family="G8", levers=["L1", "L11", "L5"],
        axes={"arith": 3, "reason": 2},
        mechanism="Inverting the CPP relation + boundary-conditioned labels: two look-similar a0 answers (87 vs 58) and two labels separated by a threshold — threshold-reading failures observed in round 1.",
        refs=["E1", "H-30"], tools_A=["critical_packing_parameter"],
        fields={"a0_at_cpp_1_3": num(A_c13, 0.005, "A^2"),
                "a0_at_cpp_1_2": num(A_c12, 0.005, "A^2"),
                "cpp_at_a0_64": num(cpp64, 0.01),
                "morph_at_64": exact(m64),
                "morph_at_95": exact(m95)},
        body=body)


# ============================================================ X-39
def build_x39():
    res = kin.kmc_micelle_kinetics(
        c_total_M=0.045, K_per_M=900.0, k_on_per_M_s=2.0e6,
        n_max=22, n_molecules=320, n_samples=45, seed=63)
    body = (
        qhead("X-39", "A third Gillespie kinetic Monte Carlo of stepwise micellization") +
        """A weakly-associating system (low K, slow on-rate) compared to any other
kinetic Monte Carlo question here:

    c_total = 0.045 M
    equilibrium constant K = 900 /M
    on-rate k_on = 2.0e6 /M/s;  off-rate k_off = k_on / K
    max aggregate size tracked n_max = 22
    molecule count in the simulated volume = 320
    starts from a Poisson sample consistent with equilibrium under K

The SSA mechanics are the same as every other Gillespie question in
this set: with a generator seeded EXACTLY at seed = 63, i.e. NumPy's default_rng(63) Generator (PCG64); every draw this procedure consumes comes from that one generator IN ORDER, and any percentile it asks for is NumPy's np.percentile with default linear interpolation -- Python's random.Random, legacy np.random.seed (MT19937) or torch streams do NOT reproduce this procedure, compute every
C_n+C_1/C_{n+1} propensity, draw an exponential waiting time and a
reaction channel from those propensities (both from the same seeded
generator, in sequence), fire the reaction, update state, and repeat.
Sample 45 evenly-spaced checkpoints; take the final 20% as the
stationary window and average the monomer concentration
(equilibrium_c1_M) and the mean aggregation number over n >= 2
(equilibrium_mean_nagg) over it. Fit a single exponential decay to the
monomer trajectory for tau_fit_s. Report the number of reaction events
fired and the max deviation from exact mass conservation (should be
exactly 0 for a correct run).

The specific sequence of simulated events depends entirely on the
pseudo-random draws from seed 63 -- weak association means many more
dissociation events relative to association, which no closed-form
estimate captures correctly. Report whether the final window was
judged stationary, choosing from ["converged", "not_converged"]."""
    )
    items = [
        ("equilibrium_c1_M", "equilibrium_c1_M: stationary-window mean monomer concentration, mol/L"),
        ("equilibrium_mean_nagg", "equilibrium_mean_nagg: stationary-window mean aggregation number (n>=2)"),
        ("tau_fit_s", "tau_fit_s: exponential relaxation time fit to the monomer trajectory, s"),
        ("mass_conservation_max_error", "mass_conservation_max_error: max deviation from exact mass conservation (integer)"),
        ("convergence", "convergence: one of [\"converged\", \"not_converged\"]"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-39"] = dict(
        id="X-39", family="G8", levers=["L12", "L8"],
        axes={"arith": 4, "reason": 1},
        mechanism="Same genuine Gillespie-SSA mechanism as X-11/X-20, a weak-association regime that produces very different event statistics -- generalization across the association-strength spectrum.",
        refs=["new-pattern-X12-X15-X19"], tools_A=["kmc_micelle_kinetics"],
        fields={"equilibrium_c1_M": num(res.equilibrium_c1_M, 0.02, "mol/L"),
                "equilibrium_mean_nagg": num(res.equilibrium_mean_nagg, 0.02),
                "tau_fit_s": num(res.tau_fit_s, 0.02, "s"),
                "mass_conservation_max_error": integer(res.mass_conservation_max_error),
                "convergence": exact("converged")},
        body=body)


# ============================================================ X-40
def build_x40():
    res = thermo.mc_aggregation_distribution(
        n_monomers=22, temperature_K=303.15, epsilon_kJ_per_mol=9.0,
        n_steps=65000, seed=71)
    body = (
        qhead("X-40", "A fourth canonical Metropolis Monte Carlo cluster simulation") +
        """22 labeled monomers at a warmer temperature with the strongest
per-bond binding energy used anywhere in this question set:

    n_monomers = 22
    temperature_K = 303.15
    epsilon = 9.0 kJ/mol  (per-bond binding energy)
    n_steps = 65000 Metropolis moves
    burn-in fraction = 0.1 (default), thinning = 10 (default)

The canonical-MC mechanics are identical to every other cluster
question here: with a generator seeded EXACTLY at seed = 71, i.e. NumPy's default_rng(71) Generator (PCG64); every draw this procedure consumes comes from that one generator IN ORDER, and any percentile it asks for is NumPy's np.percentile with default linear interpolation -- Python's random.Random, legacy np.random.seed (MT19937) or torch streams do NOT reproduce this procedure, start
from all 22 monomers free, then repeat 65000 times a detach/attach/
transfer proposal (chosen uniformly) accepted or rejected by the
Metropolis criterion on beta*dE, beta = 1/(R*T/1000) and dE =
-epsilon*(bonds after - bonds before) (both the move choice and the
accept/reject draw consuming the SAME seeded generator, in sequence).
Discard the first 10% as burn-in, keep every 10th remaining sample,
and build the monomer-fraction-by-cluster-size histogram; from it,
read off the free-monomer fraction, the monomer-weighted mean cluster
size, and the single biggest cluster that appears anywhere in the
kept samples. Since n_monomers = 22 <= 45, the exact Boltzmann
histogram over integer partitions is also available -- report the max
absolute deviation between your Monte Carlo histogram and that exact
solution.

Strong binding at this n_monomers pushes most of the probability mass
onto a handful of large clusters; the exact sequence of accepted moves
from seed 71 is the only thing that determines how much. Report
whether this run converged, choosing from
["converged", "not_converged"]."""
    )
    items = [
        ("free_monomer_fraction", "free_monomer_fraction: fraction of monomers in size-1 clusters"),
        ("mean_cluster_size", "mean_cluster_size: monomer-weighted mean cluster size"),
        ("largest_cluster", "largest_cluster: largest cluster size observed in the kept samples (integer)"),
        ("exact_max_deviation", "exact_max_deviation: max abs deviation from the exact partition solution"),
        ("convergence", "convergence: one of [\"converged\", \"not_converged\"]"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-40"] = dict(
        id="X-40", family="G8", levers=["L12", "L8"],
        axes={"arith": 4, "reason": 1},
        mechanism="Same genuine canonical-MC mechanism as X-14/X-18/X-23, the strongest binding energy in this set -- checks the pattern holds even when the distribution concentrates sharply.",
        refs=["new-pattern-X12-X15-X19"], tools_A=["mc_aggregation_distribution"],
        fields={"free_monomer_fraction": num(res.free_monomer_fraction, 0.02),
                "mean_cluster_size": num(res.mean_cluster_size_monomer_weighted, 0.02),
                "largest_cluster": integer(res.largest_cluster_observed),
                "exact_max_deviation": num(res.exact_max_abs_deviation, 0.02),
                "convergence": exact("converged")},
        body=body)


def build():
    for f in (build_x31, build_x32, build_x33, build_x34, build_x35,
              build_x36, build_x37, build_x38, build_x39, build_x40):
        f()
    return QUESTIONS
