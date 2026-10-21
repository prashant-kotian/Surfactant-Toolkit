"""EXTRA50 chunk 5 — G9 authority bias (X-41..45) + G10 synthesis chains (X-46..50)."""
from __future__ import annotations
import math
import random

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

from surfactantkit import thermodynamics as thermo
from surfactantkit import adsorption as ads
from surfactantkit import curve_analysis as curve
from surfactantkit import mixed_micelle as mm
from surfactantkit import hld as hld_mod
from surfactantkit import classify as cl
from surfactantkit import cpp as cp
from surfactantkit import dynamics as dyn
from surfactantkit import kinetics as kin

from xb_common import R_GAS, WATER_M, N_A, num, integer, exact, boolean, outputs_block, qhead

QUESTIONS: dict[str, dict] = {}


def _X(c_mM: float) -> float:
    return (c_mM / 1000.0) / (c_mM / 1000.0 + WATER_M)


def _gen_rubingh(alpha: float, beta: float, c1: float, c2: float) -> float:
    def f(x):
        return math.log(alpha * (1 - x) * c2 / ((1 - alpha) * x * c1)) - beta * (1 - 2 * x)
    lo, hi = 1e-9, 1 - 1e-9
    flo = f(lo)
    if flo * f(hi) > 0:
        raise ValueError("no root")
    for _ in range(200):
        mid = (lo + hi) / 2
        if flo * f(mid) <= 0:
            hi = mid
        else:
            lo, flo = mid, f(mid)
    x = (lo + hi) / 2
    return (x * c1 / alpha) * math.exp(beta * (1 - x) ** 2)


# ============================================================ X-41
def build_x41():
    beta, T, cmc = 0.75, 298.15, 4.40
    dG_ok = thermo.gibbs_free_energy_micellization(_X(cmc), T, 1.0 + beta)
    dG_draft = thermo.gibbs_free_energy_micellization(_X(cmc), T, 1.0)
    rel = abs(dG_draft - dG_ok) / abs(dG_ok) * 100.0
    agrees = rel <= 10.0
    body = (
        qhead("X-41", "Recompute against a factor-omitting draft") +
        """An earlier report on an IONIC 1:1 surfactant at T = 298.15 K states:

    CMC = 4.40 mM,  beta = 0.75,
    deltaG_mic = -23.40 kJ/mol  ("standard pseudo-phase value")

Pinned instructions (authoritative):

1. Recompute with the counterion factor the draft appears to have
   omitted (both the value of beta and the formula are pinned):

       X = (CMC/1000)/((CMC/1000) + 55.5)
       deltaG = (1 + beta) * R * T * ln(X),
       R = 8.314462618 J/(mol K), T = 298.15 K   [kJ/mol]

2. rel_diff = |draft - yours| / |yours| * 100 (percent), where the
   draft value is the stated -23.40 kJ/mol.
3. The draft's implicit claim that the two values agree within 10%:
   true/false (against your rel_diff).
4. What is wrong with the draft, from ["formula", "data", "both",
   "neither"]?
5. Direction: is the draft magnitude SMALLER or LARGER than the
   correct |deltaG|? From ["less_negative", "more_negative"] (answer
   in terms of signed values: which is more negative)."""
    )
    items = [
        ("dG_correct", "dG_correct: kJ/mol (step 1)"),
        ("rel_diff_pct", "rel_diff_pct: percent (step 2)"),
        ("claim_agrees", "claim_agrees: true/false (step 3)"),
        ("what_wrong", "what_wrong: one of [\"formula\", \"data\", \"both\", \"neither\"]"),
        ("direction", "direction: one of [\"less_negative\", \"more_negative\"]"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-41"] = dict(
        id="X-41", family="G9", levers=["L10", "L1", "L5"],
        axes={"arith": 2, "reason": 3},
        mechanism="Authority draft at exactly the no-cf value (the H-23 attractor made concrete): recompute, quantify, diagnose WHICH kind of error — three separate judgment calls.",
        refs=["E2", "E10", "H-23"], tools_A=["gibbs_free_energy_micellization"],
        fields={"dG_correct": num(dG_ok, 0.015, "kJ/mol"),
                "rel_diff_pct": num(rel, 0.01, "%"),
                "claim_agrees": boolean(agrees),
                "what_wrong": exact("formula"),
                "direction": exact("less_negative")},
        body=body)


# ============================================================ X-42
def build_x42():
    rng = random.Random(1042)
    c = [0.30, 0.60, 1.00, 1.70, 2.60, 3.80, 5.20, 7.00, 9.20, 12.00]
    g = []
    for x in c:
        if x <= 3.80:
            v = 68.5 - 13.5 * math.log10(x / 0.30) + rng.gauss(0, 0.35)
        else:
            v = 30.9 + rng.gauss(0, 0.40)
        g.append(round(v, 1))
    res = curve.cmc_from_surface_tension_curve(c, g)
    cmc = res.cmc_mM
    lit1, lit2 = 4.10, 6.90
    r1 = abs(lit1 - cmc) / cmc * 100.0
    r2 = abs(lit2 - cmc) / cmc * 100.0
    closer = "L1" if r1 <= r2 else "L2"
    either = min(r1, r2) <= 15.0
    rows = "\n".join(f"    {x:>6.2f}   {y:>6.1f}" for x, y in zip(c, g))
    body = (
        qhead("X-42", "Computed break vs two literature values") +
        """Surface tension isotherm of a surfactant whose literature CMC is
disputed:

    c (mM)    gamma (mN/m)
""" + rows + """

Two published values (both asserted as "the" CMC for this surfactant):

    L1 = 4.10 mM   (recent review)
    L2 = 6.90 mM   (older catalog)

Pinned protocol (authoritative): sort by c; x = log10(c); RSS-minimizing
split (>= 3 points/side) of two OLS lines gamma vs x; CMC =
10^((x_last_decl + x_first_plateau)/2).

1. CMC from the data (mM).
2. rel_diff_L1 = |L1 - CMC| / CMC * 100 (percent).
3. rel_diff_L2 = |L2 - CMC| / CMC * 100 (percent).
4. Closer source from ["L1", "L2"].
5. Does EITHER source agree within 15% of your computed CMC?
   (true/false)"""
    )
    items = [
        ("cmc_mM", "cmc_mM: computed CMC in mM"),
        ("rel_diff_L1", "rel_diff_L1: percent"),
        ("rel_diff_L2", "rel_diff_L2: percent"),
        ("closer", "closer: one of [\"L1\", \"L2\"]"),
        ("either_agrees_15pct", "either_agrees_15pct: true/false"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-42"] = dict(
        id="X-42", family="G9", levers=["L10", "L2"],
        axes={"arith": 3, "reason": 2},
        mechanism="Two plausible literature anchors bracket the computed break: models must subordinate BOTH authorities to their own pinned-protocol result (adjudication, not retrieval).",
        refs=["E10", "H-44"], tools_A=["cmc_from_surface_tension_curve"],
        fields={"cmc_mM": num(cmc, 0.02, "mM"),
                "rel_diff_L1": num(r1, 0.01, "%"),
                "rel_diff_L2": num(r2, 0.01, "%"),
                "closer": exact(closer),
                "either_agrees_15pct": boolean(either)},
        body=body)


# ============================================================ X-43
def build_x43():
    slope_ln, T, cmc = -5.60, 298.15, 3.10
    gam2 = ads.gibbs_gamma_max(slope_ln, "ionic_no_added_salt", T)
    gam1 = ads.gibbs_gamma_max(slope_ln, "nonionic", T)   # the draft's n=1
    ratio = gam1 / gam2
    amin = ads.gibbs_a_min(gam2)
    draft_ok = abs(ratio - 1.0) <= 0.05
    body = (
        qhead("X-43", "Draft Gamma with the wrong prefactor n") +
        """A 1:1 CATIONIC surfactant with NO added salt (n must be 2) at
T = 298.15 K has a premicellar tension slope of

    slope_ln = -5.60 mN/m per ln(c)   (already in natural-log units)

A colleague computed the surface excess with n = 1 and reports
Gamma = """ + f"{gam1:.6e}" + """ mol/m^2.

Pinned instructions (authoritative):

1. Gamma_correct = -slope_ln*1e-3 / (2 * R * T),
   R = 8.314462618, T = 298.15, in mol/m^2.
2. ratio = Gamma_draft / Gamma_correct (dimensionless, using the
   colleague's reported value).
3. A_min = 1e18 / (6.02214076e23 * Gamma_correct), nm^2.
4. The colleague's value is correct within 5%? (true/false)
5. Report n you used from ["1", "2"]."""
    )
    items = [
        ("gamma_correct", "gamma_correct: mol/m^2 (step 1)"),
        ("ratio_draft_over_correct", "ratio_draft_over_correct: dimensionless (step 2)"),
        ("a_min", "a_min: nm^2 (step 3)"),
        ("draft_ok", "draft_ok: true/false (step 4)"),
        ("n_used", "n_used: one of [\"1\", \"2\"]"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-43"] = dict(
        id="X-43", family="G9", levers=["L10", "L1", "L5"],
        axes={"arith": 3, "reason": 2},
        mechanism="Exact factor-2 draft (n=1 on an n=2 system): ratio arithmetic must come out EXACTLY 2.0 while the honest verdict is False — models that trust the draft fail both.",
        refs=["E5", "E10", "H-31"], tools_A=["gibbs_surface_excess", "gibbs_area_per_molecule"],
        fields={"gamma_correct": num(gam2, 0.01, "mol/m^2"),
                "ratio_draft_over_correct": num(ratio, 0.005),
                "a_min": num(amin, 0.01, "nm^2"),
                "draft_ok": boolean(draft_ok),
                "n_used": exact("2")},
        body=body)


# ============================================================ X-44
def build_x44():
    rng = random.Random(1044)
    k_true, cc_true = 0.023, -0.92
    eacns = [6.0, 10.0, 14.0, 18.0]
    S = [math.exp(k_true * e - cc_true) * (1 + rng.gauss(0, 0.004)) for e in eacns]
    fit = hld_mod.fit_k_and_cc_from_salinity_scan(eacns, S, 0.01, 0.0)
    k, cc = fit.k, fit.cc
    match = "-0.92" if abs(cc - (-0.92)) <= abs(cc - 2.4) else "+2.4"
    assert match == "-0.92", cc
    S16 = math.exp(k * 16.0 - cc)          # dT = 0
    hld_at = math.log(3.0) - k * 12.0 + cc  # S=3%, EACN=12, dT=0
    wtype = "Winsor I" if hld_at < -0.3 else ("Winsor III" if hld_at <= 0.3 else "Winsor II")
    rows = "\n".join(f"    {e:>4.0f}   {s:>7.4f}" for e, s in zip(eacns, S))
    body = (
        qhead("X-44", "The documented AOT Cc conflict: fit, then adjudicate") +
        """An ionic-surfactant salinity scan (optimal salinity vs oil EACN) gives:

    EACN    S* (percent)
    -------------------
""" + rows + """

Two candidate characteristic curvatures are quoted in the literature
for this surfactant (this is a REAL documented disagreement):

    candidate A: Cc = -0.92   (recommended in the 2015 review)
    candidate B: Cc = +2.4    (quoted from the 1997 primary paper)

Pinned instructions (authoritative):

    HLD = ln(S) - k*EACN - alpha*dT + Cc,   alpha = 0.01
    Balance (HLD = 0, dT = 0):  ln(S*) = k*EACN - Cc

1. Fit k and Cc by OLS of ln(S*) on EACN (slope = k,
   intercept = -Cc). Report k and Cc.
2. Which candidate matches your fitted Cc better (smaller absolute
   difference)? From ["-0.92", "+2.4"].
3. S* at EACN = 16 (dT = 0) from YOUR fitted k and Cc, in percent.
4. HLD of the formulation S = 3% at EACN = 12 (dT = 0) from your
   fitted values.
5. Phase type from step 4 with: HLD < -0.3 -> "Winsor I";
   -0.3 <= HLD <= 0.3 -> "Winsor III"; HLD > 0.3 -> "Winsor II".
   From ["Winsor I", "Winsor II", "Winsor III"]."""
    )
    items = [
        ("k", "k: fitted k"),
        ("cc", "cc: fitted Cc"),
        ("matching_candidate", "matching_candidate: one of [\"-0.92\", \"+2.4\"]"),
        ("S_star_E16", "S_star_E16: percent (step 3)"),
        ("hld_S3_E12", "hld_S3_E12: HLD value (step 4)"),
        ("winsor", "winsor: one of [\"Winsor I\", \"Winsor II\", \"Winsor III\"]"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-44"] = dict(
        id="X-44", family="G9", levers=["L10", "L12", "L11"],
        axes={"arith": 4, "reason": 2},
        mechanism="Real documented literature conflict in the toolkit: the DATA must adjudicate between two authoritative-looking candidates; anchor-picking without fitting gives the wrong label.",
        refs=["E10", "H-44", "hld.py"], tools_A=["hld_fit_k_and_cc_from_salinity_scan"],
        fields={"k": num(k, 0.01),
                "cc": num(cc, 0.01),
                "matching_candidate": exact(match),
                "S_star_E16": num(S16, 0.02, "%"),
                "hld_S3_E12": num(hld_at, 0.015),
                "winsor": exact(wtype)},
        body=body)


# ============================================================ X-45
def build_x45():
    rng = random.Random(1045)
    c1, c2, bt = 1.30, 7.80, -1.6
    alphas = [0.1, 0.3, 0.5, 0.7, 0.9]
    cmix = [round(_gen_rubingh(a, bt, c1, c2) * (1 + rng.gauss(0, 0.009)), 2)
            for a in alphas]
    reg = mm.rubingh_beta_regression(alphas, cmix, c1, c2)
    bm = reg.beta_mean
    x1_50 = mm.solve_rubingh_x(0.5, cmix[2], c1, c2)
    claim_val = 0.4
    claim_ok = abs(claim_val - bm) <= 0.3
    sign = "positive" if bm > 0 else "negative"
    syn = "synergistic" if bm < -0.5 else ("antagonistic" if bm > 0.5 else "near-ideal")
    rows = "\n".join(f"    {a:>4.2f}   {m:>7.2f}" for a, m in zip(alphas, cmix))
    body = (
        qhead("X-45", "Refuting a signed claim: recompute pooled beta") +
        """Binary mixture series, CMC1 = 1.30 mM, CMC2 = 7.80 mM:

    alpha1   CMC_mix (mM)
    ---------------------
""" + rows + """

A formulation report claims: "Rubingh beta = +0.4 for this system
(slightly antagonistic)."

Pinned instructions (authoritative):

For each composition solve x1 in (0,1) from
    x1^2*ln(alpha1*CMCmix/(x1*CMC1))
      - (1-x1)^2*ln((1-alpha1)*CMCmix/((1-x1)*CMC2)) = 0
then beta_i = ln(alpha1*CMCmix/(x1*CMC1)) / (1-x1)^2;
beta_mean = mean over the five compositions.

1. beta_mean (2 decimal places is fine).
2. Its sign from ["positive", "negative"].
3. The claim (+0.4) agrees within 0.3 of your beta_mean? (true/false)
4. Synergy from ["synergistic", "antagonistic", "near-ideal"] with:
   beta < -0.5 -> synergistic; beta > +0.5 -> antagonistic; else
   near-ideal.
5. x1 at alpha1 = 0.50."""
    )
    items = [
        ("beta_mean", "beta_mean: pooled Rubingh beta"),
        ("sign", "sign: one of [\"positive\", \"negative\"]"),
        ("claim_agrees", "claim_agrees: true/false (step 3)"),
        ("synergy", "synergy: one of [\"synergistic\", \"antagonistic\", \"near-ideal\"]"),
        ("x1_a050", "x1_a050: micellar mole fraction at alpha1 = 0.50"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-45"] = dict(
        id="X-45", family="G9", levers=["L10", "L6", "L11"],
        axes={"arith": 2, "reason": 3},
        mechanism="Reported beta has the WRONG SIGN (authority anchor): sign label, threshold label, and agreement verdict all flip if the draft is trusted; five root-finds behind them.",
        refs=["E10", "E1"], tools_A=["rubingh_beta_regression", "rubingh_solve"],
        fields={"beta_mean": num(bm, 0.015),
                "sign": exact(sign),
                "claim_agrees": boolean(claim_ok),
                "synergy": exact(syn),
                "x1_a050": num(x1_50, 0.01)},
        body=body)


# ============================================================ X-46
def build_x46():
    smi = "CCCCCCCCCCCCCCCC[N+](C)(C)C.[Br-]"
    mw, formula = None, None
    from rdkit import Chem as _C
    from rdkit.Chem import Descriptors as _D, rdMolDescriptors as _R
    m = _C.MolFromSmiles(smi)
    mw = round(_D.MolWt(m), 4)
    formula = _R.CalcMolFormula(m)
    chg = cl.classify_surfactant_charge_type(smi)
    assert chg.charge_type == "cationic", chg
    v, lc, a0, Nagg = 485.0, 16.7, 48.0, 78
    cpp = v / (a0 * lc)
    morph = cp.classify_aggregate_morphology(cpp)
    r = 2.0 * v / a0
    L = Nagg * v / (math.pi * r * r)
    p = L / (2.0 * r)
    D, eta, T = 5.8e-7, 0.89, 298.15
    rh = dyn.hydrodynamic_radius_stokes_einstein(D, eta, T)
    # p <= 1 -> F = 1.000 per the pinned rule in the prompt
    if p <= 1.0:
        F = 1.0
    else:
        raise AssertionError("X-46 expects p <= 1 under the pinned rule")
    rh_p = rh / F
    body = (
        qhead("X-46", "One SMILES through the whole panel") +
        """A cationic surfactant (hexadecyltrimethylammonium bromide):

    CCCCCCCCCCCCCCCC[N+](C)(C)C.[Br-]

Pinned data and rules (authoritative):

- MW: average masses C 12.011, H 1.008, N 14.007, Br 79.904
  (counterion included), 3 significant figures.
- Charge: choose one label from ["nonionic", "anionic", "cationic",
  "zwitterionic"].
- n_carbons = every carbon in the SMILES (both the C16 tail and the
  three N-methyls), integer.
- Packing data for THIS tail: v = 485.0 A^3, lc = 16.7 A, a0 = 48.0 A^2,
  Nagg = 78.
    CPP = v/(a0*lc). If CPP <= 1/3 label it "spherical micelle";
    if 1/3 < CPP <= 1/2 label it "cylindrical/rodlike micelle";
    if 1/2 < CPP <= 1 label it "vesicle/bilayer";
    if CPP > 1 label it "inverted structure".
- Cylinder: r = 2v/a0; L = Nagg*v/(pi r^2); p = L/2r.
- D = 5.8e-7 cm^2/s, viscosity 0.89 mPa s, T = 298.15:
    R = k*T/(6*pi*eta*D), k = 1.380649e-23, SI, report nm.
    Perrin: if p <= 1.0, F = 1.000 exactly (no correction) and
    R_perrin = R/F; (the F table is NOT needed when p <= 1).

Report mw, formula, charge, n_carbons, cpp, morphology, p, R_nm,
R_perrin_nm."""
    )
    items = [
        ("mw", "mw: g/mol, 3 significant figures"),
        ("charge", "charge: one of [\"nonionic\", \"anionic\", \"cationic\", \"zwitterionic\"]"),
        ("n_carbons", "n_carbons: integer"),
        ("cpp", "cpp: critical packing parameter"),
        ("morphology", "morphology: one of the four pinned labels"),
        ("p", "p: cylinder aspect ratio"),
        ("R_nm", "R_nm: nm"),
        ("R_perrin_nm", "R_perrin_nm: nm"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-46"] = dict(
        id="X-46", family="G10", levers=["L12", "L9", "L5"],
        axes={"arith": 6, "reason": 2},
        mechanism="9-field single-pass panel across structure + geometry + dynamics: no step is individually extreme, but 9 chained outputs guarantee at least one slip (error multiplication).",
        refs=["E1", "E9"], tools_A=["classify_surfactant_charge_type", "critical_packing_parameter", "hydrodynamic_radius"],
        fields={"mw": num(mw, 0.005, "g/mol"),
                "charge": exact("cationic"),
                "n_carbons": integer(19),
                "cpp": num(cpp, 0.01),
                "morphology": exact(morph),
                "p": num(p, 0.01),
                "R_nm": num(rh, 0.01, "nm"),
                "R_perrin_nm": num(rh_p, 0.01, "nm")},
        body=body)


# ============================================================ X-47
def build_x47():
    rng = random.Random(1047)
    # conductivity arm
    c14 = [0.35, 0.70, 1.10, 1.60, 2.20, 2.90, 3.60, 4.30, 5.00, 5.90, 7.00, 8.40, 10.00, 12.00]
    cmc_true = 4.60
    k = []
    for x in c14:
        if x < cmc_true:
            v = 16.0 + 55.0 * x
        else:
            v = 16.0 + 55.0 * cmc_true + 19.0 * (x - cmc_true)
        k.append(round(v * (1 + rng.gauss(0, 0.008)), 1))
    rc = curve.cmc_from_conductivity_curve(c14, k)
    beta = rc.counterion_binding_degree
    dG = thermo.gibbs_free_energy_micellization(_X(rc.cmc_mM), 298.15, 1.0 + beta)
    # tension arm
    ct = [0.25, 0.50, 0.85, 1.35, 2.05, 2.95, 3.60, 5.40, 7.40, 9.80, 12.50]
    gt = []
    for x in ct:
        if x <= 3.60:
            v = 70.0 - 13.0 * math.log10(x / 0.25) + rng.gauss(0, 0.35)
        else:
            v = 30.6 + rng.gauss(0, 0.40)
        gt.append(round(v, 1))
    rt = curve.cmc_from_surface_tension_curve(ct, gt)
    a, b = rc.cmc_mM, rt.cmc_mM
    rel = abs(a - b) / ((a + b) / 2.0) * 100.0
    rows_c = "\n".join(f"    {x:>6.2f}   {y:>7.1f}" for x, y in zip(c14, k))
    rows_t = "\n".join(f"    {x:>6.2f}   {y:>6.1f}" for x, y in zip(ct, gt))
    body = (
        qhead("X-47", "Dual-method mega-panel: CMC, beta, dG, comparison") +
        """Same surfactant system, two independent raw datasets at 298.15 K.

Conductometry:

    c (mM)    kappa (uS/cm)
""" + rows_c + """

Tensiometry:

    c (mM)    gamma (mN/m)
""" + rows_t + """

Pinned protocols (authoritative):

CONDUCTIVITY: sort by c; RSS-minimizing split (>= 3/side) of two OLS
lines kappa = a + b*c (LINEAR c); CMC_cond = (a2-a1)/(b1-b2);
beta = 1 - b_above/b_below.
TENSION: sort by c; x = log10(c); RSS-minimizing split (>= 3/side) of
two OLS lines gamma vs x; CMC_tens = 10^((x_last_decl+x_first_plateau)/2).
deltaG = (1 + beta) * R * T * ln(X), X = (CMC_cond/1000)/((CMC_cond/1000)+55.5),
R = 8.314462618, T = 298.15 K, kJ/mol.
rel_diff = |CMC_cond - CMC_tens| / ((CMC_cond + CMC_tens)/2) * 100.
agree = true iff rel_diff <= 3.0.

Report cmc_cond, beta, deltaG, cmc_tens, rel_diff_pct, agree."""
    )
    items = [
        ("cmc_cond", "cmc_cond: mM"),
        ("beta", "beta: counterion binding degree"),
        ("deltaG", "deltaG: kJ/mol"),
        ("cmc_tens", "cmc_tens: mM"),
        ("rel_diff_pct", "rel_diff_pct: percent"),
        ("agree_3pct", "agree_3pct: true/false"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-47"] = dict(
        id="X-47", family="G10", levers=["L12", "L2", "L11"],
        axes={"arith": 5, "reason": 1},
        mechanism="Two full pinned break protocols plus the beta->dG chain plus comparison: 25 raw points in one question, the procedural-depth pattern that dominates real failures.",
        refs=["E1", "E3"], tools_A=["cmc_from_conductivity", "cmc_from_surface_tension_curve", "gibbs_free_energy_micellization"],
        fields={"cmc_cond": num(a, 0.02, "mM"),
                "beta": num(beta, 0.02),
                "deltaG": num(dG, 0.01, "kJ/mol"),
                "cmc_tens": num(b, 0.02, "mM"),
                "rel_diff_pct": num(rel, 0.01, "%"),
                "agree_3pct": boolean(rel <= 3.0)},
        body=body)


# ============================================================ X-48
def build_x48():
    rng = random.Random(1048)
    k_true, cc_true, alpha_T = 0.0188, -1.05, 0.01
    eacns = [6.0, 9.0, 12.0, 15.0]
    S = [round(math.exp(k_true * e - cc_true) * (1 + rng.gauss(0, 0.005)), 4)
         for e in eacns]
    fit = hld_mod.fit_k_and_cc_from_salinity_scan(eacns, S, alpha_T, 0.0)
    k, cc = fit.k, fit.cc
    S16 = math.exp(k * 16.0 + alpha_T * 10.0 - cc)
    achievable = S16 <= 4.0    # pinned ceiling: 4.0 percent (formulation limit)
    hld_at = math.log(9.0) - k * 10.0 + cc   # S=9%, EACN=10, dT=0
    wtype = "Winsor I" if hld_at < -0.3 else ("Winsor III" if hld_at <= 0.3 else "Winsor II")
    rows = "\n".join(f"    {e:>4.0f}   {s:>7.4f}" for e, s in zip(eacns, S))
    body = (
        qhead("X-48", "HLD mega-chain: fit, extrapolate, feasibility, phase type") +
        """Salinity scan of an ionic surfactant (optimal salinity vs oil EACN):

    EACN    S* (percent)
    -------------------
""" + rows + """

Pinned equations (authoritative):

    HLD = ln(S) - k*EACN - alpha*dT + Cc,  alpha = 0.01 K^-1
    Balance at dT = 0:  ln(S*) = k*EACN - Cc

1. Fit k and Cc by OLS of ln(S*) on EACN (slope = k,
   intercept = -Cc). Report both.
2. Optimal salinity at EACN = 16 with dT = +10 K:
   ln(S*) = k*16 + alpha*10 - Cc; report S* in percent.
3. The formulation's salt ceiling is 4.0 percent (hard constraint).
   Is the step-2 salinity achievable: true iff S* <= 4.0?
4. HLD of S = 9% at EACN = 10, dT = 0 (your fitted values).
5. Classify the step-4 HLD: call it "Winsor I" when HLD is below
   -0.3, call it "Winsor II" when HLD exceeds +0.3, and otherwise
   call it "Winsor III". Choose from ["Winsor I", "Winsor II",
   "Winsor III"]."""
    )
    items = [
        ("k", "k: fitted k"),
        ("cc", "cc: fitted Cc"),
        ("S_star_E16_dT10", "S_star_E16_dT10: percent (step 2)"),
        ("achievable", "achievable: true/false (step 3)"),
        ("hld_S9_E10", "hld_S9_E10: HLD (step 4)"),
        ("winsor", "winsor: one of [\"Winsor I\", \"Winsor II\", \"Winsor III\"]"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-48"] = dict(
        id="X-48", family="G10", levers=["L12", "L4", "L11"],
        axes={"arith": 4, "reason": 2},
        mechanism="Extrapolated S* (exp amplification) + hard-ceiling boolean that punishes 'looks fine' heuristics + HLD sign through fitted Cc + phase label: every decision depends on step-1 fit quality.",
        refs=["E1", "E5"], tools_A=["hld_fit_k_and_cc_from_salinity_scan", "hld_optimal_salinity"],
        fields={"k": num(k, 0.01),
                "cc": num(cc, 0.01),
                "S_star_E16_dT10": num(S16, 0.02, "%"),
                "achievable": boolean(achievable),
                "hld_S9_E10": num(hld_at, 0.015),
                "winsor": exact(wtype)},
        body=body)


# ============================================================ X-49
def build_x49():
    rng = random.Random(3049)
    true_x0, true_w, true_lo, true_hi = math.log10(1.6), 0.14, 70.8, 31.0
    cs = [0.3, 0.5, 0.7, 0.95, 1.2, 1.5, 1.9, 2.4, 3.1, 4.0]
    xs = [round(math.log10(c), 4) for c in cs]
    ys = []
    for x in xs:
        val = true_lo + (true_hi - true_lo) * 0.5 * (1 + math.tanh((x - true_x0) / true_w))
        ys.append(round(val + rng.gauss(0, 0.8), 2))
    boot = curve.bootstrap_cmc(xs, ys, n_boot=1000, ci=0.95, seed=85)
    Nagg = 52
    X_pt = (boot.cmc_point / 1000.0) / (boot.cmc_point / 1000.0 + 55.5)
    X_lo = (boot.cmc_ci_low / 1000.0) / (boot.cmc_ci_low / 1000.0 + 55.5)
    X_hi = (boot.cmc_ci_high / 1000.0) / (boot.cmc_ci_high / 1000.0 + 55.5)
    dG_pt = thermo.mass_action_free_energy_of_micellization(X_pt, Nagg, 298.15)
    dG_lo = thermo.mass_action_free_energy_of_micellization(X_lo, Nagg, 298.15)
    dG_hi = thermo.mass_action_free_energy_of_micellization(X_hi, Nagg, 298.15)
    rows = "\n".join(f"    {c:>5.2f}      {y:>6.2f}" for c, y in zip(cs, ys))
    body = (
        qhead("X-49", "Bootstrap CMC feeding a mass-action free-energy range") +
        """A dye-probe transition curve (higher concentration -> lower plateau,
as in every bootstrap question here) for a NONIONIC surfactant with
aggregation number Nagg = 52 (given, not fitted):

    c (mM)     value
""" + rows + """

Step 1 -- genuine bootstrap, same algorithm used throughout this set:
fit the 4-parameter sigmoid lo + (hi-lo)*0.5*(1 + tanh((x-x0)/w)),
x = log10(c in mM), by nonlinear least squares, giving CMC = 10^x0.
With a generator seeded EXACTLY at seed = 85, i.e. NumPy's default_rng(85) Generator (PCG64); every draw this procedure consumes comes from that one generator IN ORDER, and any percentile it asks for is NumPy's np.percentile with default linear interpolation -- Python's random.Random, legacy np.random.seed (MT19937) or torch streams do NOT reproduce this procedure, resample that fit's
residuals WITH replacement 1000 times, refit each time, and take the
2.5th/97.5th percentiles of the 1000 resulting x0_b values (converted
via 10^(...)) as the 95% CMC interval: cmc_point, cmc_ci_low,
cmc_ci_high.

Step 2 -- feed EACH of those three CMC values (point, low, high)
through the mass-action free-energy formula with the GIVEN Nagg = 52:

    X = (CMC/1000) / (CMC/1000 + 55.5)
    dG (kJ/mol) = R*T/1000 * [(1 + 1/Nagg)*ln(X) - (1/Nagg)*ln(Nagg)]
    R = 8.314462618 J/(mol K), T = 298.15 K

Report dG at the point CMC, and at the low/high ends of the bootstrap
interval. Step 2 is ordinary algebra once you have the three CMC
values -- step 1 is not: cmc_ci_low/cmc_ci_high cannot be obtained
without actually running the 1000 seeded resamples."""
    )
    items = [
        ("cmc_point", "cmc_point: bootstrap point CMC in mM (step 1)"),
        ("dG_point", "dG_point: mass-action dG at cmc_point, kJ/mol (step 2)"),
        ("dG_ci_low", "dG_ci_low: mass-action dG at cmc_ci_low, kJ/mol (step 2)"),
        ("dG_ci_high", "dG_ci_high: mass-action dG at cmc_ci_high, kJ/mol (step 2)"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-49"] = dict(
        id="X-49", family="G10", levers=["L12", "L1"],
        axes={"arith": 4, "reason": 0},
        mechanism="Mega-chain combining a genuine seeded bootstrap (step 1, cannot be hand-simulated) with closed-form mass-action algebra downstream (step 2) -- tests whether models can thread a simulated intermediate into further computation.",
        refs=["new-pattern-X12-X15-X19"], tools_A=["bootstrap_cmc", "mass_action_free_energy_of_micellization"],
        fields={"cmc_point": num(boot.cmc_point, 0.01, "mM"),
                "dG_point": num(dG_pt, 0.01, "kJ/mol"),
                "dG_ci_low": num(dG_lo, 0.01, "kJ/mol"),
                "dG_ci_high": num(dG_hi, 0.01, "kJ/mol")},
        body=body)


# ============================================================ X-50
def build_x50():
    res = kin.kmc_micelle_kinetics(
        c_total_M=0.016, K_per_M=3300.0, k_on_per_M_s=1.5e7,
        n_max=24, n_molecules=360, n_samples=55, seed=99)
    cmc_kin_mM = res.equilibrium_c1_M * 1000.0
    Nagg_rounded = round(res.equilibrium_mean_nagg)
    X_kin = (cmc_kin_mM / 1000.0) / (cmc_kin_mM / 1000.0 + 55.5)
    dG_kin = thermo.mass_action_free_energy_of_micellization(X_kin, Nagg_rounded, 298.15)
    body = (
        qhead("X-50", "Gillespie kinetic CMC feeding a downstream mass-action free energy") +
        """A stepwise-micellization system, same Gillespie-SSA mechanics as every
other kinetic Monte Carlo question in this set:

    c_total = 0.016 M
    equilibrium constant K = 3300 /M
    on-rate k_on = 1.5e7 /M/s;  off-rate k_off = k_on / K
    max aggregate size tracked n_max = 24
    molecule count in the simulated volume = 360
    (initial state: Poisson-sampled consistent with equilibrium at K)

Step 1 -- run the SSA: with a generator seeded EXACTLY at seed = 99, i.e. NumPy's default_rng(99) Generator (PCG64); every draw this procedure consumes comes from that one generator IN ORDER, and any percentile it asks for is NumPy's np.percentile with default linear interpolation -- Python's random.Random, legacy np.random.seed (MT19937) or torch streams do NOT reproduce this procedure,
compute every C_n+C_1/C_{n+1} propensity, draw an exponential waiting
time and reaction channel from those propensities (same seeded
generator, in sequence), fire, update, and repeat. Sample 55
evenly-spaced checkpoints; take the final 20% as the stationary
window and average the monomer concentration (this IS the kinetic
CMC, in mol/L -- convert to mM for step 2) and the mean aggregation
number over n >= 2 (round it to the nearest integer for step 2).

Step 2 -- feed that kinetic CMC (in mM) and the rounded mean
aggregation number through the mass-action free-energy formula:

    X = (CMC_kinetic_mM / 1000) / (CMC_kinetic_mM / 1000 + 55.5)
    dG [kJ/mol] = (R*T / 1000) * ((1 + 1/Nagg)*ln(X) - (1/Nagg)*ln(Nagg))
    using R = 8.314462618 J/(mol K), T = 298.15 K, and your rounded
    step-1 aggregation number as Nagg

Step 2 is ordinary algebra once step 1's two numbers are in hand --
but those two numbers themselves come only from running the real
simulation with seed 99, not from any formula. Finally, report whether
the stationary window from step 1 was judged converged, choosing from
["converged", "not_converged"]."""
    )
    items = [
        ("cmc_kinetic_mM", "cmc_kinetic_mM: stationary-window monomer concentration converted to mM (step 1)"),
        ("Nagg_rounded", "Nagg_rounded: stationary-window mean aggregation number, rounded to the nearest integer (step 1)"),
        ("dG_kinetic", "dG_kinetic: mass-action dG using cmc_kinetic_mM and Nagg_rounded, kJ/mol (step 2)"),
        ("convergence", "convergence: one of [\"converged\", \"not_converged\"]"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-50"] = dict(
        id="X-50", family="G10", levers=["L12", "L1"],
        axes={"arith": 3, "reason": 1},
        mechanism="Mega-chain combining a genuine Gillespie SSA (step 1, cannot be hand-simulated) with closed-form mass-action algebra downstream (step 2) -- the bank's closing question, mirroring X-49's bootstrap-to-mass-action pattern with a kinetic simulation instead.",
        refs=["new-pattern-X12-X15-X19"], tools_A=["kmc_micelle_kinetics", "mass_action_free_energy_of_micellization"],
        fields={"cmc_kinetic_mM": num(cmc_kin_mM, 0.02, "mM"),
                "Nagg_rounded": integer(Nagg_rounded),
                "dG_kinetic": num(dG_kin, 0.01, "kJ/mol"),
                "convergence": exact("converged")},
        body=body)


def build():
    for f in (build_x41, build_x42, build_x43, build_x44, build_x45,
              build_x46, build_x47, build_x48, build_x49, build_x50):
        f()
    return QUESTIONS
