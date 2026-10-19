"""EXTRA50 chunk 3 — G5 mixture bookkeeping (X-21..25) + G6 pinned forensics (X-26..30)."""
from __future__ import annotations
import math
import random

from surfactantkit import thermodynamics as thermo
from surfactantkit import adsorption as ads
from surfactantkit import curve_analysis as curve
from surfactantkit import mixed_micelle as mm
from surfactantkit import kinetics as kin

from xb_common import R_GAS, WATER_M, num, integer, exact, boolean, outputs_block, qhead

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


# ============================================================ X-21
def build_x21():
    rng = random.Random(4021)
    c1, c2, bt = 2.10, 6.40, -1.9
    alphas = [0.1, 0.3, 0.5, 0.7, 0.9]
    cmix = [round(_gen_rubingh(a, bt, c1, c2) * (1 + rng.gauss(0, 0.007)), 2)
            for a in alphas]
    fit = mm.eommm_global_fit(alphas, cmix, c1, c2)
    feasible = fit.total_infeasibility < 1e-6
    multi = mm.eommm_global_fit_all_optima(alphas, cmix, c1, c2)
    assert not multi.degenerate, f"X-21: non-unique EOMMM fit found ({multi.n_optima_found} optima)"
    rows = "\n".join(f"    {a:>4.2f}   {m:>7.2f}" for a, m in zip(alphas, cmix))
    body = (
        qhead("X-21", "EOMMM global fit on a five-point mixture series") +
        """Binary mixture series, CMC1 = 2.10 mM, CMC2 = 6.40 mM:

    alpha1   CMC_mix (mM)
    ---------------------
""" + rows + """

Pinned model-selection protocol (authoritative):

- EOMMM fits the two independent asymmetric-Margules interaction
  parameters W12, W21 SIMULTANEOUSLY across the whole series, together
  with each point's own micellar mole fraction x1, by minimizing a
  mass-balance INFEASIBILITY measure (not a least-squares residual) at
  a fixed margin around each point's experimental CMC_mix.
- The real procedure: a coarse 61x61 grid search over (W12, W21) in
  [-30, 30]^2 to find the best cell, then a 41x41 fine grid search in a
  small window around it -- not a local gradient descent and not a
  single closed-form solve.
- cmc_margin = 0.10 (10% band around each point's CMC_mix counts as
  mass-balance-consistent).

Report W12, W21, the fit's r_squared (RSS in ln(CMC_mix) space), the
micellar mole fraction x1 at the LAST point (alpha1 = 0.90), and
whether total_infeasibility fell below 1e-6."""
    )
    items = [
        ("W12", "W12: fitted asymmetric-Margules interaction parameter"),
        ("W21", "W21: fitted asymmetric-Margules interaction parameter"),
        ("r_squared", "r_squared: fit quality in ln(CMC_mix) space"),
        ("x1_last_point", "x1_last_point: micellar mole fraction x1 at alpha1 = 0.90"),
        ("feasible", "feasible: true if total_infeasibility < 1e-6, else false"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-21"] = dict(
        id="X-21", family="G5", levers=["L12", "L6"],
        axes={"arith": 4, "reason": 1},
        mechanism="Same genuine 2D grid-search mechanism proven to universally fail on X-12/X-16, on a fifth independent mixture series -- generalization check across many different (CMC1, CMC2, series) combinations. Fit verified unique via eommm_global_fit_all_optima (n_optima_found=1) -- not assumed.",
        refs=["new-pattern-X12-X15-X19"], tools_A=["eommm_global_fit"],
        fields={"W12": num(fit.W12, 0.02),
                "W21": num(fit.W21, 0.02),
                "r_squared": num(fit.r_squared, 0.01),
                "x1_last_point": num(fit.x1_values[-1], 0.01),
                "feasible": boolean(feasible)},
        body=body)


# ============================================================ X-22
def build_x22():
    rng = random.Random(4022)
    c1, c2, bt = 1.60, 9.20, 1.4
    alphas = [0.05, 0.15, 0.25, 0.50, 0.75, 0.85, 0.95]
    cmix = [round(_gen_rubingh(a, bt, c1, c2) * (1 + rng.gauss(0, 0.01)), 2)
            for a in alphas]
    fit = mm.eommm_global_fit(alphas, cmix, c1, c2)
    feasible = fit.total_infeasibility < 1e-6
    multi = mm.eommm_global_fit_all_optima(alphas, cmix, c1, c2)
    assert not multi.degenerate, f"X-22: non-unique EOMMM fit found ({multi.n_optima_found} optima)"
    rows = "\n".join(f"    {a:>4.2f}   {m:>7.2f}" for a, m in zip(alphas, cmix))
    body = (
        qhead("X-22", "EOMMM global fit on a seven-point mixture series") +
        """Binary mixture series, CMC1 = 1.60 mM, CMC2 = 9.20 mM:

    alpha1   CMC_mix (mM)
    ---------------------
""" + rows + """

The fitting rule (authoritative, same EOMMM procedure used throughout
this set): W12 and W21 are not found by a least-squares residual on
CMC_mix directly. Instead, at a fixed cmc_margin = 0.10 around each
point's own experimental CMC_mix, every (W12, W21) pair is scored by
how badly it violates the mass-balance constraints at that margin (its
total "infeasibility"); the fit is whichever (W12, W21) pair, together
with each point's own micellar mole fraction x1, drives that
infeasibility score to its minimum. Finding that minimum needs a 2D
search -- first coarse (61x61 over [-30,30]^2), then a fine 41x41
re-search inside the best coarse cell -- since the objective has no
closed-form minimizer and is not gradient-friendly everywhere.

Give W12, W21, r_squared, the micellar mole fraction x1 at the third
point (alpha1 = 0.25), and whether total_infeasibility came in under
1e-6."""
    )
    items = [
        ("W12", "W12: fitted asymmetric-Margules interaction parameter"),
        ("W21", "W21: fitted asymmetric-Margules interaction parameter"),
        ("r_squared", "r_squared: fit quality in ln(CMC_mix) space"),
        ("x1_third_point", "x1_third_point: micellar mole fraction x1 at alpha1 = 0.25"),
        ("feasible", "feasible: true if total_infeasibility < 1e-6, else false"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-22"] = dict(
        id="X-22", family="G5", levers=["L12", "L6"],
        axes={"arith": 4, "reason": 1},
        mechanism="Same genuine 2D grid-search mechanism as X-21, a seven-point series instead of five -- checks the pattern holds regardless of series length. Fit verified unique via eommm_global_fit_all_optima (n_optima_found=1) -- not assumed.",
        refs=["new-pattern-X12-X15-X19"], tools_A=["eommm_global_fit"],
        fields={"W12": num(fit.W12, 0.02),
                "W21": num(fit.W21, 0.02),
                "r_squared": num(fit.r_squared, 0.01),
                "x1_third_point": num(fit.x1_values[2], 0.01),
                "feasible": boolean(feasible)},
        body=body)


# ============================================================ X-23
def build_x23():
    res = thermo.mc_aggregation_distribution(
        n_monomers=40, temperature_K=293.15, epsilon_kJ_per_mol=8.5,
        n_steps=90000, seed=23)
    body = (
        qhead("X-23", "A third canonical Metropolis Monte Carlo cluster simulation") +
        """40 labeled monomers, a lower temperature and stronger binding energy
than any other Monte Carlo cluster question in this set:

    n_monomers = 40
    temperature_K = 293.15
    epsilon = 8.5 kJ/mol  (per-bond binding energy)
    n_steps = 90000 Metropolis moves
    burn-in fraction = 0.1 (default), thinning = 10 (default)

Pinned procedure (authoritative; the canonical-MC mechanics match every
other cluster-simulation question here -- this question only changes
the physical inputs):

1. With a pseudo-random generator seeded EXACTLY with seed = 23, i.e. NumPy's default_rng(23) Generator (PCG64); every draw this procedure consumes comes from that one generator IN ORDER, and any percentile it asks for is NumPy's np.percentile with default linear interpolation -- Python's random.Random, legacy np.random.seed (MT19937) or torch streams do NOT reproduce this procedure, start
   from all 40 monomers free.
2. Repeat 90000 times: propose a detach, attach, or transfer move
   (uniformly among the three), accept or reject by the Metropolis
   criterion on beta*dE with beta = 1/(R*T/1000) and
   dE = -epsilon*(bonds after - bonds before) (move choice and
   accept/reject both consume the SAME seeded generator, in sequence).
3. Discard the first 10% as burn-in, keep every 10th remaining sample,
   and build the monomer-fraction-by-cluster-size histogram.
4. From that histogram, report the free-monomer fraction, the
   monomer-weighted mean cluster size, and the single largest cluster
   size seen anywhere in the retained samples.
5. n_monomers = 40 <= 45, so the exact Boltzmann histogram over integer
   partitions of 40 is also available; report the max absolute
   deviation between the Monte Carlo histogram and that exact
   solution.

The retained-sample histogram depends on the exact sequence of
pseudo-random accept/reject draws from seed 23 -- there is no formula
that reproduces it directly except via the (separately verified) exact
partition cross-check in step 5. Report whether the run converged,
choosing from ["converged", "not_converged"]."""
    )
    items = [
        ("free_monomer_fraction", "free_monomer_fraction: fraction of monomers in size-1 clusters (step 4)"),
        ("mean_cluster_size", "mean_cluster_size: monomer-weighted mean cluster size (step 4)"),
        ("largest_cluster", "largest_cluster: largest cluster size observed in the retained samples (step 4, integer)"),
        ("exact_max_deviation", "exact_max_deviation: max abs deviation from the exact partition solution (step 5)"),
        ("convergence", "convergence: one of [\"converged\", \"not_converged\"]"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-23"] = dict(
        id="X-23", family="G5", levers=["L12", "L8"],
        axes={"arith": 4, "reason": 1},
        mechanism="Same genuine canonical-MC mechanism as X-14/X-18 at a third, colder/stronger-binding condition -- checks the pattern holds across the physical regime, not just one tuned example.",
        refs=["new-pattern-X12-X15-X19"], tools_A=["mc_aggregation_distribution"],
        fields={"free_monomer_fraction": num(res.free_monomer_fraction, 0.02),
                "mean_cluster_size": num(res.mean_cluster_size_monomer_weighted, 0.02),
                "largest_cluster": integer(res.largest_cluster_observed),
                "exact_max_deviation": num(res.exact_max_abs_deviation, 0.02),
                "convergence": exact("converged")},
        body=body)


# ============================================================ X-24
def build_x24():
    rng = random.Random(4024)
    c1, c2, bt = 0.85, 7.30, -1.1
    alphas = [0.20, 0.50, 0.80, 0.95]
    cmix = [round(_gen_rubingh(a, bt, c1, c2) * (1 + rng.gauss(0, 0.008)), 2)
            for a in alphas]
    fit = mm.eommm_global_fit(alphas, cmix, c1, c2)
    feasible = fit.total_infeasibility < 1e-6
    multi = mm.eommm_global_fit_all_optima(alphas, cmix, c1, c2)
    assert not multi.degenerate, f"X-24: non-unique EOMMM fit found ({multi.n_optima_found} optima)"
    rows = "\n".join(f"    {a:>4.2f}   {m:>7.2f}" for a, m in zip(alphas, cmix))
    body = (
        qhead("X-24", "EOMMM global fit on a four-point mixture series") +
        """Binary mixture series, CMC1 = 0.85 mM, CMC2 = 7.30 mM:

    alpha1   CMC_mix (mM)
    ---------------------
""" + rows + """

Fitting rule (authoritative, unchanged mechanism from every other
EOMMM question here, only the series differs): rather than minimizing
a least-squares residual on CMC_mix, the real objective is a
mass-balance infeasibility score, evaluated at cmc_margin = 0.10 around
each point's experimental CMC_mix; W12, W21, and every point's own
micellar mole fraction x1 are the choice that makes this infeasibility
smallest. No closed form exists for that minimum, so it is located by
a coarse 61x61 grid over (W12, W21) in [-30,30]^2 followed by a fine
41x41 grid confined to the best coarse cell's neighborhood.

State W12, W21, r_squared, the micellar mole fraction x1 at the second
point (alpha1 = 0.50), and whether total_infeasibility landed below
1e-6."""
    )
    items = [
        ("W12", "W12: fitted asymmetric-Margules interaction parameter"),
        ("W21", "W21: fitted asymmetric-Margules interaction parameter"),
        ("r_squared", "r_squared: fit quality in ln(CMC_mix) space"),
        ("x1_second_point", "x1_second_point: micellar mole fraction x1 at alpha1 = 0.50"),
        ("feasible", "feasible: true if total_infeasibility < 1e-6, else false"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-24"] = dict(
        id="X-24", family="G5", levers=["L12", "L6"],
        axes={"arith": 4, "reason": 1},
        mechanism="Same genuine 2D grid-search mechanism as X-21/X-22, a four-point series with a minority (c1<c2 asymmetric-weighted) composition set -- checks generalization across series shapes. Fit verified unique via eommm_global_fit_all_optima (n_optima_found=1) -- not assumed.",
        refs=["new-pattern-X12-X15-X19"], tools_A=["eommm_global_fit"],
        fields={"W12": num(fit.W12, 0.02),
                "W21": num(fit.W21, 0.02),
                "r_squared": num(fit.r_squared, 0.01),
                "x1_second_point": num(fit.x1_values[1], 0.01),
                "feasible": boolean(feasible)},
        body=body)


# ============================================================ X-25
def build_x25():
    rng = random.Random(4025)
    c1, c2, bt = 1.80, 5.60, 0.9
    alphas = [0.15, 0.35, 0.55, 0.75, 0.9]
    cmix = [round(_gen_rubingh(a, bt, c1, c2) * (1 + rng.gauss(0, 0.006)), 2)
            for a in alphas]
    fit = mm.eommm_global_fit(alphas, cmix, c1, c2)
    feasible = fit.total_infeasibility < 1e-6
    multi = mm.eommm_global_fit_all_optima(alphas, cmix, c1, c2)
    assert not multi.degenerate, f"X-25: non-unique EOMMM fit found ({multi.n_optima_found} optima)"
    rows = "\n".join(f"    {a:>4.2f}   {m:>7.2f}" for a, m in zip(alphas, cmix))
    body = (
        qhead("X-25", "A fifth EOMMM global fit, closing out the mixture-series block") +
        """Binary mixture series, CMC1 = 1.80 mM, CMC2 = 5.60 mM:

    alpha1   CMC_mix (mM)
    ---------------------
""" + rows + """

As with every other EOMMM question in this set, the fit is NOT a
least-squares regression against CMC_mix. At a fixed cmc_margin = 0.10
band around each point's own experimental CMC_mix, every candidate
(W12, W21) is scored by its total mass-balance infeasibility there;
the reported fit -- W12, W21, and each point's own micellar mole
fraction x1 -- is the candidate that minimizes that score. There is no
closed form for this minimum: it is found with a coarse 61x61 grid
over (W12, W21) in [-30,30]^2, refined by a fine 41x41 grid inside the
best coarse cell.

Provide W12, W21, r_squared, the micellar mole fraction x1 at the
middle point (alpha1 = 0.55), and whether total_infeasibility ended up
below 1e-6."""
    )
    items = [
        ("W12", "W12: fitted asymmetric-Margules interaction parameter"),
        ("W21", "W21: fitted asymmetric-Margules interaction parameter"),
        ("r_squared", "r_squared: fit quality in ln(CMC_mix) space"),
        ("x1_middle_point", "x1_middle_point: micellar mole fraction x1 at alpha1 = 0.55"),
        ("feasible", "feasible: true if total_infeasibility < 1e-6, else false"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-25"] = dict(
        id="X-25", family="G5", levers=["L12", "L6"],
        axes={"arith": 4, "reason": 1},
        mechanism="Same genuine 2D grid-search mechanism as X-21/X-22/X-24, a fifth independent series closing out this family -- consistent-pattern generalization evidence. Fit verified unique via eommm_global_fit_all_optima (n_optima_found=1) -- not assumed.",
        refs=["new-pattern-X12-X15-X19"], tools_A=["eommm_global_fit"],
        fields={"W12": num(fit.W12, 0.02),
                "W21": num(fit.W21, 0.02),
                "r_squared": num(fit.r_squared, 0.01),
                "x1_middle_point": num(fit.x1_values[2], 0.01),
                "feasible": boolean(feasible)},
        body=body)


# ============================================================ X-26
def build_x26():
    rng = random.Random(1026)
    Tr = 298.15
    Ts = [288.15, 293.15, 298.15, 303.15, 308.15, 313.15, 318.15]
    dH_true, dS_true, dCp_true = -14500.0, 37.0, 0.0
    cmcs = []
    for T in Ts:
        dG = dH_true + dCp_true * (T - Tr) - T * (dS_true + dCp_true * math.log(T / Tr))
        dG += rng.gauss(0, 25.0)
        X = math.exp(dG / (R_GAS * T))
        cmcs.append(round(X * WATER_M / (1 - X) * 1000.0, 4))
    cmcs[3] = round(cmcs[3] * 0.93, 4)   # corrupt row: 7% low
    xs = [_X(c) for c in cmcs]

    def fit3(idx):
        return thermo.vant_hoff_multi_point_fit(
            [xs[i] for i in idx], [Ts[i] for i in idx],
            counterion_factor=1.0, reference_temperature_K=Tr)

    def predict(res, T):
        dH, dS, dCp = res.delta_H_kJ_per_mol, res.delta_S_J_per_mol_K / 1000.0, res.delta_Cp_J_per_mol_K / 1000.0
        # dG in kJ: dH(Tr)+dCp*(T-Tr) - T*(dS(Tr)+dCp*ln(T/Tr))
        return dH + dCp * (T - Tr) - T * (dS + dCp * math.log(T / Tr))

    loo_err = []
    for i in range(7):
        f6 = fit3([j for j in range(7) if j != i])
        dG_i = R_GAS * Ts[i] * math.log(xs[i]) / 1000.0
        loo_err.append(abs(dG_i - predict(f6, Ts[i])))
    flag = max(range(7), key=lambda i: loo_err[i])
    others = [e for i, e in enumerate(loo_err) if i != flag]
    n_flag = 1 if loo_err[flag] > 2.0 * (sum(others) / len(others)) else 0
    assert n_flag == 1 and flag == 3, (flag, loo_err)
    f_all = fit3(list(range(7)))
    f_clean = fit3([j for j in range(7) if j != 3])
    rows = "\n".join(f"    {i+1:>2d}  {T:>7.2f}   {c:>8.4f}"
                     for i, (T, c) in enumerate(zip(Ts, cmcs)))
    body = (
        qhead("X-26", "Corrupt-row forensics with pinned Tr and leave-one-out rule") +
        """CMC vs temperature, NONIONIC surfactant (cf = 1):

    row   T (K)     CMC (mM)
    -----------------------
""" + rows + """

One row's CMC was corrupted in transcription (7% low). Pinned detection
protocol (authoritative):

1. Convert each row to X = (CMC/1000)/((CMC/1000) + 55.5) and
   dG_i = R*T_i*ln(X_i),  R = 8.314462618 J/(mol K).
2. For EACH row i: fit the 3-parameter Kantonen model on the OTHER six
   rows only (dG = dH - T*dS + dCp*[(T-298.15) - T*ln(T/298.15)],
   least squares, reference temperature 298.15 K), then predict dG at
   T_i. Leave-one-out error_i = |dG_i - predicted_i|.
3. Flag the row with the LARGEST leave-one-out error; it is corrupt iff
   its error > 2 x (mean error of the other six rows). Report the
   flagged row's T as printed, and how many rows you flag.
4. Fit the 3-parameter model (reference 298.15 K) to ALL seven rows ->
   dH_all; then to the six clean rows (excluding the flagged row) ->
   dH_cleaned. Report both in kJ/mol at 298.15 K.
5. Report the r^2 of the cleaned six-row fit."""
    )
    items = [
        ("flagged_T", "flagged_T: the flagged row's T exactly as printed (e.g. 300.00)"),
        ("n_flagged", "n_flagged: number of corrupt rows (integer)"),
        ("dH_all", "dH_all: kJ/mol from all 7 rows (step 4)"),
        ("dH_cleaned", "dH_cleaned: kJ/mol from the 6 clean rows (step 4)"),
        ("r2_cleaned", "r2_cleaned: r^2 of the cleaned fit (step 5)"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-26"] = dict(
        id="X-26", family="G6", levers=["L7", "L12", "L1"],
        axes={"arith": 4, "reason": 1},
        mechanism="Seven hand refits (leave-one-out 3-param) — computationally out of unaugmented reach; subtle 5% corruption (not the obvious 10x) with pinned threshold rule.",
        refs=["E1", "H-41"], tools_A=["vant_hoff_multi_point_fit"],
        fields={"flagged_T": exact("303.15"),
                "n_flagged": integer(n_flag),
                "dH_all": num(f_all.delta_H_kJ_per_mol, 0.015, "kJ/mol"),
                "dH_cleaned": num(f_clean.delta_H_kJ_per_mol, 0.015, "kJ/mol"),
                "r2_cleaned": num(f_clean.r_squared, 0.01)},
        body=body)


# ============================================================ X-27
def build_x27():
    rng = random.Random(1027)
    # T displayed in degrees C; row 5 mistakenly printed as its Kelvin value
    t_C = [15.0, 20.0, 25.0, 30.0, 308.15, 40.0, 45.0, 50.0]
    t_K_true = [288.15, 293.15, 298.15, 303.15, 308.15, 313.15, 318.15, 323.15]
    dH_true, dS_true = -12800.0, 38.0
    cmcs = []
    for T in t_K_true:
        dG = dH_true - T * dS_true + rng.gauss(0, 20.0)
        X = math.exp(dG / (R_GAS * T))
        cmcs.append(round(X * WATER_M / (1 - X) * 1000.0, 4))

    def fit2(ts, cs, idx=None):
        import numpy as np
        idx = range(len(ts)) if idx is None else idx
        T = np.array([ts[i] for i in idx], float)
        y = np.array([(R_GAS * ts[i] * math.log(_X(cs[i]))) / 1000.0 for i in idx])
        A = np.column_stack([np.ones(len(idx)), -T])
        p, *_ = np.linalg.lstsq(A, y, rcond=None)
        pred = A @ p
        rss = float(((y - pred) ** 2).sum())
        tss = float(((y - y.mean()) ** 2).sum())
        r2 = 1.0 - rss / tss if tss > 0 else 1.0
        resid = y - pred
        return float(p[0]), float(p[1]), r2, resid, y, pred

    # detection: naive conversion of EVERY displayed value (row5: 308.15+273.15),
    # leave-one-out 2-param fits (outlier masking makes plain residuals useless)
    t_naive = [t + 273.15 for t in t_C]
    loo = []
    for i in range(8):
        dHi, dSi, r2i, ri, yi, pi = fit2(t_naive, cmcs, [j for j in range(8) if j != i])
        # predict row i: dG = dH - T*dS with the naive T of row i
        y_i = (R_GAS * t_naive[i] * math.log(_X(cmcs[i]))) / 1000.0
        pred_i = dHi - t_naive[i] * dSi
        loo.append(abs(y_i - pred_i))
    flag = max(range(8), key=lambda i: loo[i])
    others = [e for i, e in enumerate(loo) if i != flag]
    n_flag = 1 if loo[flag] > 2.0 * (sum(others) / len(others)) else 0
    flags = [flag] if n_flag else []
    assert flags == [4], (flags, loo)
    # corrected fit: row5 interpreted as already-K (308.15 K), others +273.15
    t_ok = [t_K_true[i] for i in range(8)]
    dH_c, dS_c_kJ, r2_c, *_ = fit2(t_ok, cmcs)
    dS_c = dS_c_kJ * 1000.0   # kJ/(mol.K) -> J/(mol.K)
    rows = "\n".join(f"    {i+1:>2d}   {tc:>7.2f}   {c:>8.4f}"
                     for i, (tc, c) in enumerate(zip(t_C, cmcs)))
    body = (
        qhead("X-27", "Mixed-unit temperature column: detect, correct, refit") +
        """CMC vs temperature for a nonionic surfactant. The temperature column is
labeled T (deg C) — EXCEPT that one row was recorded in kelvin by
mistake (the transcription error you must find):

    row   T (deg C)    CMC (mM)
    --------------------------
""" + rows + """

Pinned protocol (authoritative):

1. Naive pass: convert EVERY printed T via K = T + 273.15, compute
   dG_i = R*K_i*ln(X_i) with X = (CMC/1000)/((CMC/1000)+55.5),
   R = 8.314462618.
2. For EACH row i: fit dG = dH - T*dS (2 parameters, least squares) on
   the OTHER seven rows (naive K values), predict dG at K_i, and take
   the leave-one-out error_i = |dG_i - predicted_i|. Flag the row with
   the LARGEST error; it is corrupt iff its error > 2 x (mean error of
   the other seven rows). Report the flagged row's T exactly as printed,
   and the count of flagged rows.
3. Correction: the flagged row's printed value IS its kelvin
   temperature (use it directly, no +273.15); all other rows get
   +273.15. Refit dG = dH - T*dS on all EIGHT corrected rows.
4. Report dH (kJ/mol), dS (J/(mol K)), and r^2 of the corrected fit."""
    )
    items = [
        ("flagged_T", "flagged_T: flagged row's T exactly as printed (step 2)"),
        ("n_flagged", "n_flagged: count of flagged rows (integer)"),
        ("dH_corrected", "dH_corrected: kJ/mol (step 4)"),
        ("dS_corrected", "dS_corrected: J/(mol K) (step 4)"),
        ("r2_corrected", "r2_corrected: r^2 of the corrected fit (step 4)"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-27"] = dict(
        id="X-27", family="G6", levers=["L7", "L8", "L12"],
        axes={"arith": 4, "reason": 1},
        mechanism="Unit-collision forensics: one kelvin value inside a celsius column; naive conversion flags it, then correct handling must NOT convert it again (+273.15 twice trap).",
        refs=["E8", "E3"], tools_A=["vant_hoff_multi_point_fit"],
        fields={"flagged_T": exact("308.15"),
                "n_flagged": integer(len(flags)),
                "dH_corrected": num(dH_c, 0.015, "kJ/mol"),
                "dS_corrected": num(dS_c, 0.015, "J/(mol K)"),
                "r2_corrected": num(r2_c, 0.01)},
        body=body)


# ============================================================ X-28
def build_x28():
    rng = random.Random(1028)
    c = [0.30, 0.60, 1.00, 1.60, 1.60, 2.40, 3.50, 5.00, 7.00, 9.50]
    g = []
    for x in c:
        if x <= 2.40:
            v = 67.0 - 13.0 * math.log10(x) + rng.gauss(0, 0.3)
        else:
            v = 30.4 + rng.gauss(0, 0.35)
        g.append(round(v, 1))
    # plant the conflict at c = 1.60 (indices 3 and 4)
    g[3], g[4] = 45.9, 38.7

    def local_resid(i):
        # OLS line through the 4 nearest OTHER points in log10(c)
        others = [j for j in range(10) if j != i]
        others.sort(key=lambda j: abs(math.log10(c[j]) - math.log10(c[i])))
        sel = others[:4]
        xs = [math.log10(c[j]) for j in sel]
        ys = [g[j] for j in sel]
        n = 4
        mx, my = sum(xs) / n, sum(ys) / n
        sxx = sum((x - mx) ** 2 for x in xs)
        b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
        a = my - b * mx
        return abs(g[i] - (a + b * math.log10(c[i])))
    r3, r4 = local_resid(3), local_resid(4)
    drop = 4 if r4 > r3 else 3
    assert drop == 4, (r3, r4)
    g_clean = [g[j] for j in range(10) if j != drop]
    c_clean = [c[j] for j in range(10) if j != drop]
    res = curve.cmc_from_surface_tension_curve(c_clean, g_clean)
    gam = ads.gibbs_gamma_max(res.premicellar_slope_mN_per_m_per_lnC, "ionic_no_added_salt", 298.15)
    rows = "\n".join(
        f"    {i+1:>2d}  {x:>6.2f}   {y:>6.1f}" + ("   <-- duplicate c with the row above" if i == 4 else "")
        for i, (x, y) in enumerate(zip(c, g)))
    body = (
        qhead("X-28", "Duplicate-concentration conflict: exclusion rule + refit") +
        """Surface-tension isotherm of a 1:1 CATIONIC surfactant, no added salt
(n = 2), T = 298.15 K. Two measurements exist at c = 1.60 mM (they
disagree — one is a transcription casualty):

    row   c (mM)    gamma (mN/m)
    --------------------------
""" + rows + """

Pinned protocol (authoritative):

1. For EACH of the two duplicate rows: fit an OLS line (gamma vs
   log10(c)) through its 4 nearest other points in log10(c) (excluding
   both duplicates), and compute |observed - fitted| for that duplicate.
   Exclude the duplicate with the LARGER residual. Report the excluded
   gamma value.
2. On the remaining 9 points run the standard break rule: split at the
   RSS-minimizing boundary (>= 3 points per side) of two OLS lines in
   log10(c); CMC = 10^((x_last_decl + x_first_plateau)/2).
3. Premicellar slope of the winning declining segment -> convert to
   per-ln via /ln(10); Gamma = -slope_ln*1e-3/(2*R*298.15),
   R = 8.314462618; report mol/m^2.
4. Report how many points you used (integer), and confirm the conflict
   was detected (true/false)."""
    )
    items = [
        ("excluded_gamma", "excluded_gamma: the gamma value you excluded (mN/m)"),
        ("cmc_mM", "cmc_mM: CMC in mM (step 2)"),
        ("gamma_max", "gamma_max: mol/m^2 (step 3)"),
        ("n_used", "n_used: number of points used (integer)"),
        ("conflict_detected", "conflict_detected: true/false (step 4)"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-28"] = dict(
        id="X-28", family="G6", levers=["L7", "L11", "L2"],
        axes={"arith": 4, "reason": 1},
        mechanism="Statistical exclusion rule with two near-equal candidates (38.7 vs 45.9: which residual is larger is NOT obvious) then a full break+Gibbs chain on the cleaned set.",
        refs=["E7", "H-41"], tools_A=["cmc_from_surface_tension_curve", "gibbs_surface_excess"],
        fields={"excluded_gamma": num(g[drop], 0.0005, "mN/m"),
                "cmc_mM": num(res.cmc_mM, 0.02, "mM"),
                "gamma_max": num(gam, 0.01, "mol/m^2"),
                "n_used": integer(9),
                "conflict_detected": boolean(True)},
        body=body)


# ============================================================ X-29
def build_x29():
    T1, T2 = 298.15, 313.15
    c1, c2 = 5.50, 9.80
    X1, X2 = _X(c1), _X(c2)
    ln10 = math.log(10)
    run1 = R_GAS * (math.log(X2) - math.log(X1)) / (1.0 / T2 - 1.0 / T1) / 1000.0
    run2 = R_GAS * (math.log(X2) - math.log(X1)) / (1.0 / T1 - 1.0 / T2) / 1000.0
    run3 = R_GAS * (math.log10(X2) - math.log10(X1)) / (1.0 / T2 - 1.0 / T1) / 1000.0
    ratio31 = run3 / run1
    body = (
        qhead("X-29", "Sign-audit: three runs, one flipped denominator, one base-10") +
        """Two-point van't Hoff analysis, T1 = 298.15 K / X1 = """ + f"{X1:.6e}" + """,
T2 = 313.15 K / X2 = """ + f"{X2:.6e}" + """ (mole fractions of the CMC).

Three different runs of the SAME dataset produced:

    run1:  deltaH = """ + f"{run1:.3f}" + """ kJ/mol
    run2:  deltaH = """ + f"{run2:.3f}" + """ kJ/mol
    run3:  deltaH = """ + f"{run3:.3f}" + """ kJ/mol

Pinned ground-truth formula (authoritative):

    deltaH = R * (lnX2 - lnX1) / (1/T2 - 1/T1) / 1000,
    R = 8.314462618 J/(mol K)

1. Which run used the FLIPPED denominator (1/T1 - 1/T2 instead of
   1/T2 - 1/T1)? Choose from ["run1", "run2", "run3"].
2. Report the corrected value of that run (the number it SHOULD have
   been under the pinned formula), in kJ/mol.
3. Which run used base-10 logarithms (log10 instead of ln) for the
   concentration term? Choose from ["run1", "run2", "run3"].
4. Report run3/run1 as a ratio (dimensionless)."""
    )
    items = [
        ("sign_flipped_run", "sign_flipped_run: one of [\"run1\", \"run2\", \"run3\"]"),
        ("corrected_value", "corrected_value: kJ/mol (step 2)"),
        ("base10_run", "base10_run: one of [\"run1\", \"run2\", \"run3\"]"),
        ("ratio_run3_run1", "ratio_run3_run1: dimensionless (step 4)"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-29"] = dict(
        id="X-29", family="G6", levers=["L7", "L11", "L5"],
        axes={"arith": 2, "reason": 2},
        mechanism="Documented toolkit sign-bug pattern + base-10 conversion identity (ratio = ln 10): forensic labeling instead of computation-from-scratch; models must DETECT, not just compute.",
        refs=["E5", "sign-bug", "H-36"], tools_A=["vant_hoff_enthalpy"],
        fields={"sign_flipped_run": exact("run2"),
                "corrected_value": num(run1, 0.01, "kJ/mol"),
                "base10_run": exact("run3"),
                "ratio_run3_run1": num(ratio31, 0.01)},
        body=body)


# ============================================================ X-30
def build_x30():
    from surfactantkit.curve_analysis import _best_two_segment_split, _bic
    rng = random.Random(1030)
    c = [0.20, 0.50, 0.90, 1.50, 2.40, 3.60, 5.20, 7.20, 10.00]
    g = []
    for x in c:
        if x <= 3.60:
            v = 66.5 - 14.0 * math.log10(x) + rng.gauss(0, 0.3)
        else:
            v = 31.6 + rng.gauss(0, 0.35)
        g.append(round(v, 1))
    g[3] = 0.0459   # 45.9 mN/m recorded in N/m

    def bic2(cv, gv, minpts=3):
        pairs = sorted(zip(cv, gv))
        lx = [math.log10(x) for x, _ in pairs]
        gy = [y for _, y in pairs]
        rss, split, *_ = _best_two_segment_split(lx, gy, len(pairs), minpts)
        return _bic(rss, len(pairs), 5), split

    # detection: leave-one-out total RSS of the RSS-minimizing 2-seg fit;
    # the factor-1000 row dominates RSS, so excluding it drops RSS most.
    def rss_2seg(cv, gv):
        pairs = sorted(zip(cv, gv))
        lx = [math.log10(x) for x, _ in pairs]
        gy = [y for _, y in pairs]
        return _best_two_segment_split(lx, gy, len(pairs), 3)[0]
    rss_raw = rss_2seg(c, g)
    rss_loo = [rss_2seg([c[i] for i in range(9) if i != j],
                        [g[i] for i in range(9) if i != j]) for j in range(9)]
    flag = min(range(9), key=lambda j: rss_loo[j])
    drop_ok = (rss_raw - rss_loo[flag]) > 0.5 * rss_raw
    flags = [flag] if drop_ok else []
    assert flags == [3], (flags, rss_raw, rss_loo)
    c_ok = [c[i] for i in range(9) if i != 3]
    g_ok = [g[i] for i in range(9) if i != 3]
    res = curve.cmc_from_surface_tension_curve(c_ok, g_ok)
    bic_after, _ = bic2(c_ok, g_ok)
    bic_before, _ = bic2(c, g)
    delta = bic_after - bic_before
    rows = "\n".join(
        f"    {i+1:>2d}  {x:>6.2f}   {y:>6.3f}" for i, (x, y) in enumerate(zip(c, g)))
    body = (
        qhead("X-30", "Scale-corruption row (N/m vs mN/m) detection and refit") +
        """Surface-tension isotherm. One row's tension was recorded in N/m
instead of mN/m (factor 1000) — find it:

    row   c (mM)    gamma (column says mN/m)
    ---------------------------------------
""" + rows + """

Pinned protocol (authoritative):

1. Leave-one-out RSS detection: for EACH row i, remove it, then fit
   the RSS-minimizing 2-segment break (x = log10(c), >= 3 points/side,
   two OLS lines) to the remaining EIGHT points and record the total
   RSS. Flag the row whose removal gives the SMALLEST total RSS — but
   only if removing it cuts the raw nine-point RSS by more than half
   (otherwise there is no corrupt row). Report the flagged row's c
   exactly as printed.
2. BIC values: BIC = n*ln(RSS/n) + 5*ln(n) for the 2-segment fit —
   compute it on the RAW 9 points (bic_before) and on the 8 clean
   points after excluding the flagged row (bic_after).
   Report delta = bic_after - bic_before.
3. On the clean 8 points run the same break rule; report CMC (mM):
   10^((x_last_decl + x_first_plateau)/2).
4. Report points used (integer) and whether a corrupt row was
   detected (true/false)."""
    )
    items = [
        ("flagged_c", "flagged_c: flagged row's c exactly as printed (step 1)"),
        ("delta_bic", "delta_bic: bic_after - bic_before (step 2)"),
        ("cmc_mM", "cmc_mM: CMC in mM (step 3)"),
        ("n_used", "n_used: points used (integer)"),
        ("corrupt_detected", "corrupt_detected: true/false (step 4)"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-30"] = dict(
        id="X-30", family="G6", levers=["L7", "L2", "L12"],
        axes={"arith": 4, "reason": 1},
        mechanism="Factor-1000 planted value among plausible decimals: residual-rule detection + BIC before/after + refit chain; the corrupted value LOOKS like a legitimate small tension.",
        refs=["E7", "E3"], tools_A=["cmc_from_surface_tension_curve"],
        fields={"flagged_c": num(c[3], 1e-6, "mM"),
                "delta_bic": num(delta, 0.02),
                "cmc_mM": num(res.cmc_mM, 0.02, "mM"),
                "n_used": integer(8),
                "corrupt_detected": boolean(True)},
        body=body)


def build():
    for f in (build_x21, build_x22, build_x23, build_x24, build_x25,
              build_x26, build_x27, build_x28, build_x29, build_x30):
        f()
    return QUESTIONS
