"""EXTRA50 chunk 1 — G1 genuine-computation execution (X-01..05, rewritten
2026-10-08: bootstrap/Monte-Carlo questions whose gold requires an actual
seeded simulation, not a closed-form formula) + G2 break estimation (X-06..10).
All golds computed by SurfactantKit library functions from DISPLAYED inputs."""
from __future__ import annotations
import math
import random

from surfactantkit import thermodynamics as thermo
from surfactantkit import adsorption as ads
from surfactantkit import curve_analysis as curve
from surfactantkit import inference as infer
from surfactantkit import kinetics as kin

from xb_common import R_GAS, WATER_M, N_A, num, integer, exact, boolean, outputs_block, qhead

QUESTIONS: dict[str, dict] = {}


def _X(c_mM: float) -> float:
    return (c_mM / 1000.0) / (c_mM / 1000.0 + WATER_M)


# ============================================================ X-01
def build_x01():
    rng = random.Random(2001)
    true_x0, true_w, true_lo, true_hi = math.log10(3.2), 0.18, 68.0, 29.5
    cs = [0.5, 0.8, 1.2, 1.8, 2.4, 3.0, 3.6, 4.8, 6.0, 7.5]
    xs = [round(math.log10(c), 4) for c in cs]
    ys = []
    for x in xs:
        val = true_lo + (true_hi - true_lo) * 0.5 * (1 + math.tanh((x - true_x0) / true_w))
        ys.append(round(val + rng.gauss(0, 1.1), 2))
    res = curve.bootstrap_cmc(xs, ys, n_boot=1000, ci=0.95, seed=0)
    ci_width_pct = 100.0 * (res.cmc_ci_high - res.cmc_ci_low) / res.cmc_point
    rows = "\n".join(f"    {c:>5.2f}      {y:>6.2f}" for c, y in zip(cs, ys))
    body = (
        qhead("X-01", "Residual-bootstrap CMC: point estimate and 95% interval") +
        """A dye-probe transition curve for an ionic surfactant (property value vs.
concentration; higher concentration -> lower plateau):

    c (mM)     value
""" + rows + """

Pinned procedure (authoritative -- every step is exact, nothing is left
to convention):

1. Let x = log10(c in mM). Fit the 4-parameter sigmoidal transition
       value(x) = lo + (hi - lo) * 0.5 * (1 + tanh((x - x0) / w))
   to the 10 points above by nonlinear least squares (there is no closed
   form for this fit; it requires an iterative solver). CMC = 10^x0.
2. Compute the fit's own residuals r_i = value_i - fitted_i at each of
   the 10 points.
3. Using a pseudo-random generator seeded EXACTLY with seed = 0, i.e. NumPy's default_rng(0) Generator (PCG64); every draw this procedure consumes comes from that one generator IN ORDER, and any percentile it asks for is NumPy's np.percentile with default linear interpolation -- Python's random.Random, legacy np.random.seed (MT19937) or torch streams do NOT reproduce this procedure, draw
   1000 bootstrap resamples: for each of the 1000 repeats, resample the
   10 residuals WITH replacement (some residuals may appear more than
   once in a repeat, others not at all), add the resampled residuals
   back onto the original fitted curve, and refit the same 4-parameter
   sigmoid to that resampled curve to get a bootstrap x0_b.
4. From the 1000 bootstrap x0_b values (the full discard rule: a refit is discarded if and only if its least-squares call
   raises/fails or its fitted transition width w lies outside the open interval (1e-4, 10*(xmax-xmin))
   on the log10 axis; every other refit is kept), take the 2.5th and 97.5th
   percentiles and convert back via 10^(...) to get the 95% confidence
   interval on the CMC.
5. ci_width_pct = 100 * (CMC_high - CMC_low) / CMC_point.

This is a literal simulation procedure, not a formula you can evaluate
by hand: steps 3-4 require actually constructing a pseudo-random
generator from the stated seed and running 1000 real nonlinear refits,
not estimating a plausible-looking interval. Report whether all 1000
resamples converged to a usable fit, choosing from
["converged", "not_converged"]."""
    )
    items = [
        ("cmc_point", "cmc_point: point CMC estimate in mM (step 1)"),
        ("cmc_ci_low", "cmc_ci_low: 2.5th-percentile bootstrap CMC in mM (step 4)"),
        ("cmc_ci_high", "cmc_ci_high: 97.5th-percentile bootstrap CMC in mM (step 4)"),
        ("ci_width_pct", "ci_width_pct: 100*(cmc_ci_high-cmc_ci_low)/cmc_point (step 5)"),
        ("convergence", "convergence: one of [\"converged\", \"not_converged\"]"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-01"] = dict(
        id="X-01", family="G1", levers=["L12", "L8"],
        axes={"arith": 4, "reason": 1},
        mechanism="Gold requires an actual seeded residual-bootstrap (1000 real nonlinear refits), not a closed-form formula -- the full procedure is stated, but text-reasoning cannot execute 1000 refits, so a model must fabricate a plausible-looking interval rather than compute the real one.",
        refs=["new-pattern-X12-X15-X19"], tools_A=["bootstrap_cmc"],
        fields={"cmc_point": num(res.cmc_point, 0.01, "mM"),
                "cmc_ci_low": num(res.cmc_ci_low, 0.02, "mM"),
                "cmc_ci_high": num(res.cmc_ci_high, 0.02, "mM"),
                "ci_width_pct": num(ci_width_pct, 0.02, "%"),
                "convergence": exact("converged")},
        body=body)


# ============================================================ X-02
def build_x02():
    rng = random.Random(2002)
    true_x0, true_w, true_lo, true_hi = math.log10(5.6), 0.22, 0.42, 1.95
    cs = [1.0, 1.6, 2.4, 3.2, 4.0, 5.0, 6.2, 7.6, 9.5, 12.0, 15.0]
    xs = [round(math.log10(c), 4) for c in cs]
    ys = []
    for x in xs:
        val = true_lo + (true_hi - true_lo) * 0.5 * (1 + math.tanh((x - true_x0) / true_w))
        ys.append(round(val + rng.gauss(0, 0.05), 3))
    res = curve.bootstrap_cmc(xs, ys, n_boot=1000, ci=0.90, seed=3)
    ci_width_pct = 100.0 * (res.cmc_ci_high - res.cmc_ci_low) / res.cmc_point
    wide = ci_width_pct > 5.0
    rows = "\n".join(f"    {c:>5.2f}      {y:>6.3f}" for c, y in zip(cs, ys))
    body = (
        qhead("X-02", "Residual-bootstrap CMC on a rising (conductivity-like) curve") +
        """A molar-conductivity-like property for an ionic surfactant (property
value vs. concentration; higher concentration -> HIGHER plateau, the
opposite direction from a tension/dye-probe curve):

    c (mM)     value
""" + rows + """

Pinned procedure (authoritative -- identical algorithm to a standard
dye-probe fit; the only difference is this curve rises, not falls):

1. Let x = log10(c in mM). Fit the 4-parameter sigmoidal transition
       value(x) = lo + (hi - lo) * 0.5 * (1 + tanh((x - x0) / w))
   to the 11 points above by nonlinear least squares (no closed form;
   an iterative solver is required). CMC = 10^x0.
2. Compute the fit's own residuals at each of the 11 points.
3. Using a pseudo-random generator seeded EXACTLY with seed = 3, i.e. NumPy's default_rng(3) Generator (PCG64); every draw this procedure consumes comes from that one generator IN ORDER, and any percentile it asks for is NumPy's np.percentile with default linear interpolation -- Python's random.Random, legacy np.random.seed (MT19937) or torch streams do NOT reproduce this procedure, draw
   1000 bootstrap resamples of the 11 residuals WITH replacement, add
   each resample back onto the fitted curve, and refit to get a
   bootstrap x0_b each time.
4. From the 1000 bootstrap x0_b values, take the 5th and 95th
   percentiles (this is a 90% interval, not 95%) and convert via
   10^(...) to get the CMC confidence interval.
5. ci_width_pct = 100 * (CMC_high - CMC_low) / CMC_point.

As in any such bootstrap, steps 3-4 require actually running 1000 real
refits from the stated seed -- there is no shortcut formula for a
percentile bootstrap interval. Report whether ci_width_pct exceeds 5.0
(a plain true/false on the number you computed in step 5)."""
    )
    items = [
        ("cmc_point", "cmc_point: point CMC estimate in mM (step 1)"),
        ("cmc_ci_low", "cmc_ci_low: 5th-percentile bootstrap CMC in mM (step 4)"),
        ("cmc_ci_high", "cmc_ci_high: 95th-percentile bootstrap CMC in mM (step 4)"),
        ("ci_width_pct", "ci_width_pct: 100*(cmc_ci_high-cmc_ci_low)/cmc_point (step 5)"),
        ("ci_wide", "ci_wide: true if ci_width_pct > 5.0, else false"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-02"] = dict(
        id="X-02", family="G1", levers=["L12", "L8"],
        axes={"arith": 4, "reason": 1},
        mechanism="Same genuine-bootstrap pattern as X-01 on a rising curve, different seed/CI level/dataset -- tests whether a model's fabricated interval generalizes or was tuned to one example.",
        refs=["new-pattern-X12-X15-X19"], tools_A=["bootstrap_cmc"],
        fields={"cmc_point": num(res.cmc_point, 0.01, "mM"),
                "cmc_ci_low": num(res.cmc_ci_low, 0.02, "mM"),
                "cmc_ci_high": num(res.cmc_ci_high, 0.02, "mM"),
                "ci_width_pct": num(ci_width_pct, 0.02, "%"),
                "ci_wide": boolean(wide)},
        body=body)


# ============================================================ X-03
def build_x03():
    res = infer.propagate_uncertainty_chain(
        "gibbs_dg_from_cmc",
        overrides={"cmc_mole_fraction": {"value": 0.0085, "sigma": 0.0006},
                   "aggregation_number": {"value": 45, "sigma": 4},
                   "temperature_K": {"value": 303.15, "sigma": 0.3}},
        n_samples=10000, ci=0.95, seed=7)
    skew = (res.output_mean - res.output_median) / res.output_sd
    label = "right-skewed" if skew > 0.05 else ("left-skewed" if skew < -0.05 else "symmetric")
    body = (
        qhead("X-03", "Monte Carlo error propagation through the mass-action free energy") +
        """A researcher's three measured inputs, each with a stated 1-sigma
uncertainty, feed the NONIONIC mass-action free-energy-of-micellization
formula:

    cmc_mole_fraction = 0.0085, sigma = 0.0006 (normal)
    aggregation_number = 45,    sigma = 4       (normal)
    temperature_K       = 303.15, sigma = 0.3   (normal)

    deltaG_mic (kJ/mol) = R*T/1000 * [(1 + 1/n)*ln(X) - (1/n)*ln(n)]
    R = 8.314462618 J/(mol K), X = cmc_mole_fraction, n = aggregation_number

Pinned Monte Carlo procedure (authoritative; this is a real propagation,
not a linearized/first-order error bar):

1. Using a pseudo-random generator seeded EXACTLY with seed = 7, i.e. NumPy's default_rng(7) Generator (PCG64); every draw this procedure consumes comes from that one generator IN ORDER, and any percentile it asks for is NumPy's np.percentile with default linear interpolation -- Python's random.Random, legacy np.random.seed (MT19937) or torch streams do NOT reproduce this procedure, draw
   10000 independent samples. In each sample, draw X ~ Normal(0.0085,
   0.0006), n ~ Normal(45, 4), T ~ Normal(303.15, 0.3) (all three drawn
   independently, same sample index), then evaluate deltaG_mic with the
   formula above using that sample's (X, n, T).
2. Over the 10000 deltaG_mic values collected this way, report the
   mean, the sample standard deviation (ddof=1), and the 2.5th/97.5th
   percentile values (a 95% interval).
3. Let skew = (mean - median) / sd for those 10000 values. Call the
   result "right-skewed" when skew exceeds +0.05, "left-skewed" when it
   is below -0.05, and "symmetric" for anything in between.

A linearized (first-order Taylor) error propagation does NOT reproduce
this distribution's mean, sd, or interval exactly -- only an actual
10000-draw simulation from the stated seed does, since the formula is
nonlinear in X and n."""
    )
    items = [
        ("output_mean", "output_mean: mean of the 10000 deltaG_mic samples, kJ/mol"),
        ("output_sd", "output_sd: sample standard deviation (ddof=1), kJ/mol"),
        ("ci_low", "ci_low: 2.5th percentile, kJ/mol"),
        ("ci_high", "ci_high: 97.5th percentile, kJ/mol"),
        ("skew_label", "skew_label: one of [\"left-skewed\", \"symmetric\", \"right-skewed\"] (step 3)"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-03"] = dict(
        id="X-03", family="G1", levers=["L12", "L1"],
        axes={"arith": 4, "reason": 1},
        mechanism="Gold is a seeded 10000-draw Monte Carlo propagation through a nonlinear formula; the formula is fully disclosed but the output distribution's exact mean/sd/interval cannot be derived analytically, only simulated.",
        refs=["new-pattern-X12-X15-X19"], tools_A=["propagate_uncertainty_chain"],
        fields={"output_mean": num(res.output_mean, 0.01, "kJ/mol"),
                "output_sd": num(res.output_sd, 0.02, "kJ/mol"),
                "ci_low": num(res.ci_low, 0.02, "kJ/mol"),
                "ci_high": num(res.ci_high, 0.02, "kJ/mol"),
                "skew_label": exact(label)},
        body=body)


# ============================================================ X-04
def build_x04():
    res = infer.propagate_uncertainty_chain(
        "debye_length",
        overrides={"ionic_strength_M": {"value": 0.015, "sigma": 0.004},
                   "temperature_K": {"value": 298.15, "sigma": 0.2}},
        n_samples=10000, ci=0.95, seed=11)
    skew = (res.output_mean - res.output_median) / res.output_sd
    label = "right-skewed" if skew > 0.05 else ("left-skewed" if skew < -0.05 else "symmetric")
    body = (
        qhead("X-04", "Monte Carlo error propagation through the Debye screening length") +
        """Two measured inputs, each with a stated 1-sigma uncertainty, feed the
Debye screening-length formula:

    ionic_strength_M = 0.015, sigma = 0.004  (normal)
    temperature_K    = 298.15, sigma = 0.2   (normal)

    lambda_D (nm) = sqrt(eps_r*eps0*k_B*T / (2*N_A*I_SI*e^2)) * 1e9,
        I_SI = ionic_strength_M * 1000   (mol/L -> mol/m^3)
    eps_r = 78.4 (water, 25C), eps0 = 8.8541878128e-12 F/m,
    k_B = 1.380649e-23 J/K, N_A = 6.02214076e23 /mol,
    e = 1.602176634e-19 C

Run this as a real Monte Carlo propagation (authoritative; do not
substitute a linearized/first-order error bar):

1. Using a pseudo-random generator seeded EXACTLY with seed = 11, i.e. NumPy's default_rng(11) Generator (PCG64); every draw this procedure consumes comes from that one generator IN ORDER, and any percentile it asks for is NumPy's np.percentile with default linear interpolation -- Python's random.Random, legacy np.random.seed (MT19937) or torch streams do NOT reproduce this procedure, draw
   10000 independent samples. In each sample, draw I ~ Normal(0.015,
   0.004) and T ~ Normal(298.15, 0.2) (independently, same sample
   index), then evaluate lambda_D with the formula above using that
   sample's (I, T).
2. From the 10000 lambda_D values this produces, give the mean, the
   sample standard deviation (ddof=1), and the values at the 2.5th and
   97.5th percentiles (together, a 95% interval).
3. Define skew = (mean - median) / sd across those same 10000 values,
   and label the distribution "right-skewed" for skew above +0.05,
   "left-skewed" for skew below -0.05, "symmetric" otherwise.

lambda_D is proportional to 1/sqrt(I) -- a CONVEX function of I for
I > 0 -- so a symmetric input distribution on I does NOT propagate to
a symmetric output distribution; only an actual 10000-draw simulation
from the stated seed reproduces the true (skewed) distribution."""
    )
    items = [
        ("output_mean", "output_mean: mean of the 10000 lambda_D samples, nm"),
        ("output_sd", "output_sd: sample standard deviation (ddof=1), nm"),
        ("ci_low", "ci_low: 2.5th percentile, nm"),
        ("ci_high", "ci_high: 97.5th percentile, nm"),
        ("skew_label", "skew_label: one of [\"left-skewed\", \"symmetric\", \"right-skewed\"] (step 3)"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-04"] = dict(
        id="X-04", family="G1", levers=["L12", "L1"],
        axes={"arith": 4, "reason": 1},
        mechanism="Gold is a seeded 10000-draw Monte Carlo propagation through a convex (1/sqrt(I)) formula -- the true output is right-skewed; a model doing linearized symmetric error propagation gets a systematically wrong answer that cannot pass the real tolerance.",
        refs=["new-pattern-X12-X15-X19"], tools_A=["propagate_uncertainty_chain"],
        fields={"output_mean": num(res.output_mean, 0.01, "nm"),
                "output_sd": num(res.output_sd, 0.02, "nm"),
                "ci_low": num(res.ci_low, 0.02, "nm"),
                "ci_high": num(res.ci_high, 0.02, "nm"),
                "skew_label": exact(label)},
        body=body)


# ============================================================ X-05
def build_x05():
    res = infer.propagate_uncertainty_chain(
        "hydrodynamic_radius",
        overrides={"diffusion_coefficient_cm2_per_s": {"value": 6.5e-7, "sigma": 1.2e-7},
                   "viscosity_mPas": {"value": 0.89, "sigma": 0.015},
                   "temperature_K": {"value": 298.15, "sigma": 0.3}},
        n_samples=10000, ci=0.95, seed=19)
    skew = (res.output_mean - res.output_median) / res.output_sd
    label = "right-skewed" if skew > 0.05 else ("left-skewed" if skew < -0.05 else "symmetric")
    body = (
        qhead("X-05", "Monte Carlo error propagation through the Stokes-Einstein radius") +
        """Three measured inputs, each with a stated 1-sigma uncertainty, feed the
Stokes-Einstein hydrodynamic-radius formula:

    diffusion_coefficient_cm2_per_s = 6.5e-7, sigma = 1.2e-7  (normal)
    viscosity_mPas                  = 0.89,   sigma = 0.015   (normal)
    temperature_K                   = 298.15, sigma = 0.3     (normal)

    R_h (nm) = k_B*T / (6*pi*eta_SI*D_SI) * 1e9, k_B = 1.380649e-23 J/K,
    where D_SI is diffusion_coefficient_cm2_per_s rescaled from cm^2/s
    to m^2/s (x1e-4) and eta_SI is viscosity_mPas rescaled from mPa.s
    to Pa.s (x1e-3)

Treat this as a genuine Monte Carlo propagation (authoritative; a
linearized/first-order error bar is not an acceptable substitute):

1. Using a pseudo-random generator seeded EXACTLY with seed = 19, i.e. NumPy's default_rng(19) Generator (PCG64); every draw this procedure consumes comes from that one generator IN ORDER, and any percentile it asks for is NumPy's np.percentile with default linear interpolation -- Python's random.Random, legacy np.random.seed (MT19937) or torch streams do NOT reproduce this procedure, draw
   10000 independent samples. In each sample, draw D ~ Normal(6.5e-7,
   1.2e-7), eta ~ Normal(0.89, 0.015), T ~ Normal(298.15, 0.3)
   (independently, same sample index), then evaluate R_h with the
   formula above using that sample's (D, eta, T).
2. Across the 10000 resulting R_h values, state the mean, the sample
   standard deviation (ddof=1), and the 2.5th- and 97.5th-percentile
   values (a 95% interval).
3. With skew = (mean - median) / sd computed over those 10000 values,
   name the distribution "right-skewed" once skew passes +0.05,
   "left-skewed" once it passes -0.05, and "symmetric" in between.

R_h is proportional to 1/D -- a CONVEX function of D for D > 0 -- so a
symmetric input distribution on D does NOT propagate to a symmetric
output distribution; only an actual 10000-draw simulation from the
stated seed reproduces the true (skewed) distribution."""
    )
    items = [
        ("output_mean", "output_mean: mean of the 10000 R_h samples, nm"),
        ("output_sd", "output_sd: sample standard deviation (ddof=1), nm"),
        ("ci_low", "ci_low: 2.5th percentile, nm"),
        ("ci_high", "ci_high: 97.5th percentile, nm"),
        ("skew_label", "skew_label: one of [\"left-skewed\", \"symmetric\", \"right-skewed\"] (step 3)"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-05"] = dict(
        id="X-05", family="G1", levers=["L12", "L1"],
        axes={"arith": 4, "reason": 1},
        mechanism="Gold is a seeded 10000-draw Monte Carlo propagation through a convex (1/D) formula -- the true output is right-skewed, same genuine-simulation pattern as X-04 on a different chain/domain (hydrodynamics vs. electrostatics).",
        refs=["new-pattern-X12-X15-X19"], tools_A=["propagate_uncertainty_chain"],
        fields={"output_mean": num(res.output_mean, 0.01, "nm"),
                "output_sd": num(res.output_sd, 0.02, "nm"),
                "ci_low": num(res.ci_low, 0.02, "nm"),
                "ci_high": num(res.ci_high, 0.02, "nm"),
                "skew_label": exact(label)},
        body=body)


# ============================================================ X-06
def build_x06():
    rng = random.Random(3006)
    true_x0, true_w, true_lo, true_hi = math.log10(4.5), 0.20, 69.0, 33.5
    cs = [0.6, 1.0, 1.5, 2.2, 3.0, 4.0, 5.2, 6.8, 8.8, 11.0]
    xs = [round(math.log10(c), 4) for c in cs]
    ys = []
    for x in xs:
        val = true_lo + (true_hi - true_lo) * 0.5 * (1 + math.tanh((x - true_x0) / true_w))
        ys.append(round(val + rng.gauss(0, 1.0), 2))
    res = curve.bootstrap_cmc(xs, ys, n_boot=1000, ci=0.95, seed=14)
    ci_width_pct = 100.0 * (res.cmc_ci_high - res.cmc_ci_low) / res.cmc_point
    rows = "\n".join(f"    {c:>5.2f}      {y:>6.2f}" for c, y in zip(cs, ys))
    body = (
        qhead("X-06", "A fourth residual-bootstrap CMC, a broader concentration span") +
        """A dye-probe transition curve spanning a wider concentration range than
the other bootstrap questions here (higher concentration -> lower
plateau, as always for this probe type):

    c (mM)     value
""" + rows + """

The fitting rule is unchanged from any other bootstrap question in this
set, restated for self-containment: fit the 4-parameter sigmoid
value(x) = lo + (hi - lo)*0.5*(1 + tanh((x - x0)/w)) with x = log10(c
in mM) by nonlinear least squares (iterative; no closed form), giving
CMC = 10^x0. Take that fit's residuals, then with a generator seeded
EXACTLY at seed = 14, using NumPy's default_rng(seed) PCG64 bit generator as the single generator consumed in a fixed order (the order the reference implementation draws its resample indices), and percentiles computed with np.percentile's default linear interpolation -- Python's random.Random, legacy np.random.seed (MT19937) or torch streams do NOT reproduce this procedure, draw 1000 bootstrap resamples of the residuals
WITH replacement, add each resample back onto the fit, and refit to
collect 1000 bootstrap x0_b values; the 2.5th/97.5th percentiles of
those (converted via 10^(...)) give the 95% CMC interval, and
ci_width_pct = 100*(CMC_high - CMC_low)/CMC_point.

None of this is computable from the formula alone -- the specific 1000
resamples depend on the stated seed, and only running them reproduces
the interval. Report whether ci_width_pct exceeds 6.0 (true/false)."""
    )
    items = [
        ("cmc_point", "cmc_point: point CMC estimate in mM"),
        ("cmc_ci_low", "cmc_ci_low: 2.5th-percentile bootstrap CMC in mM"),
        ("cmc_ci_high", "cmc_ci_high: 97.5th-percentile bootstrap CMC in mM"),
        ("ci_width_pct", "ci_width_pct: 100*(cmc_ci_high-cmc_ci_low)/cmc_point"),
        ("ci_wide", "ci_wide: true if ci_width_pct > 6.0, else false"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-06"] = dict(
        id="X-06", family="G2", levers=["L12", "L8"],
        axes={"arith": 4, "reason": 1},
        mechanism="Same genuine residual-bootstrap mechanism as X-01/X-02/X-17, a fourth independent dataset/seed on a wider concentration span.",
        refs=["new-pattern-X12-X15-X19"], tools_A=["bootstrap_cmc"],
        fields={"cmc_point": num(res.cmc_point, 0.01, "mM"),
                "cmc_ci_low": num(res.cmc_ci_low, 0.02, "mM"),
                "cmc_ci_high": num(res.cmc_ci_high, 0.02, "mM"),
                "ci_width_pct": num(ci_width_pct, 0.02, "%"),
                "ci_wide": boolean(ci_width_pct > 6.0)},
        body=body)


# ============================================================ X-07
def build_x07():
    res = infer.propagate_uncertainty_chain(
        "gibbs_dg_from_cmc",
        overrides={"cmc_mole_fraction": {"value": 0.0062, "sigma": 0.0004},
                   "aggregation_number": {"value": 38, "sigma": 3},
                   "temperature_K": {"value": 298.15, "sigma": 0.25}},
        n_samples=10000, ci=0.95, seed=33)
    skew = (res.output_mean - res.output_median) / res.output_sd
    label = "right-skewed" if skew > 0.05 else ("left-skewed" if skew < -0.05 else "symmetric")
    body = (
        qhead("X-07", "A second Monte Carlo propagation through the mass-action free energy") +
        """A different nonionic system's measured inputs, each with a stated
1-sigma uncertainty, feed the same mass-action free-energy formula used
elsewhere in this set:

    cmc_mole_fraction = 0.0062, sigma = 0.0004 (normal)
    aggregation_number = 38,    sigma = 3       (normal)
    temperature_K       = 298.15, sigma = 0.25  (normal)

    deltaG_mic [kJ/mol] = (R*T / 1000) * ((1 + 1/n)*ln(X) - (1/n)*ln(n)),
    using R = 8.314462618 J/(mol K); X stands for cmc_mole_fraction and n
    for aggregation_number in that expression.

With a pseudo-random generator seeded EXACTLY at seed = 33, i.e. NumPy's default_rng(33) Generator (PCG64); every draw this procedure consumes comes from that one generator IN ORDER, and any percentile it asks for is NumPy's np.percentile with default linear interpolation -- Python's random.Random, legacy np.random.seed (MT19937) or torch streams do NOT reproduce this procedure, draw 10000
independent samples -- X ~ Normal(0.0062, 0.0004), n ~ Normal(38, 3),
T ~ Normal(298.15, 0.25), same sample index across all three -- and
evaluate deltaG_mic at each sample. From the resulting 10000 values,
state the mean, the sample standard deviation (ddof=1), and the
2.5th/97.5th percentiles (a 95% interval). Taking skew = (mean -
median)/sd over those 10000 values, call the distribution
"right-skewed" above skew = +0.05, "left-skewed" below -0.05, and
"symmetric" otherwise.

A first-order (linearized) error bar does not reproduce this
distribution's true mean/sd/interval -- the formula is nonlinear in X
and n, so only the actual 10000-draw simulation from seed 33 gets it
right."""
    )
    items = [
        ("output_mean", "output_mean: mean of the 10000 deltaG_mic samples, kJ/mol"),
        ("output_sd", "output_sd: sample standard deviation (ddof=1), kJ/mol"),
        ("ci_low", "ci_low: 2.5th percentile, kJ/mol"),
        ("ci_high", "ci_high: 97.5th percentile, kJ/mol"),
        ("skew_label", "skew_label: one of [\"left-skewed\", \"symmetric\", \"right-skewed\"]"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-07"] = dict(
        id="X-07", family="G2", levers=["L12", "L1"],
        axes={"arith": 4, "reason": 1},
        mechanism="Same genuine 10000-draw Monte Carlo propagation as X-03, a second independent system/seed through the same nonlinear mass-action formula.",
        refs=["new-pattern-X12-X15-X19"], tools_A=["propagate_uncertainty_chain"],
        fields={"output_mean": num(res.output_mean, 0.01, "kJ/mol"),
                "output_sd": num(res.output_sd, 0.02, "kJ/mol"),
                "ci_low": num(res.ci_low, 0.02, "kJ/mol"),
                "ci_high": num(res.ci_high, 0.02, "kJ/mol"),
                "skew_label": exact(label)},
        body=body)


# ============================================================ X-08
def build_x08():
    res = infer.propagate_uncertainty_chain(
        "hydrodynamic_radius",
        overrides={"diffusion_coefficient_cm2_per_s": {"value": 4.2e-7, "sigma": 0.9e-7},
                   "viscosity_mPas": {"value": 1.02, "sigma": 0.02},
                   "temperature_K": {"value": 298.15, "sigma": 0.25}},
        n_samples=10000, ci=0.95, seed=41)
    skew = (res.output_mean - res.output_median) / res.output_sd
    label = "right-skewed" if skew > 0.05 else ("left-skewed" if skew < -0.05 else "symmetric")
    body = (
        qhead("X-08", "A second Monte Carlo propagation through the Stokes-Einstein radius") +
        """A more viscous solvent and a smaller, noisier diffusion coefficient
than any other hydrodynamic-radius question here, each input carrying
a stated 1-sigma uncertainty:

    diffusion_coefficient_cm2_per_s = 4.2e-7, sigma = 0.9e-7  (normal)
    viscosity_mPas                  = 1.02,   sigma = 0.02    (normal)
    temperature_K                   = 298.15, sigma = 0.25    (normal)

    R_h (nm) = k_B*T / (6*pi*eta_SI*D_SI) * 1e9,
        D_SI = diffusion_coefficient_cm2_per_s * 1e-4   (cm^2/s -> m^2/s)
        eta_SI = viscosity_mPas * 1e-3                   (mPa.s -> Pa.s)
    k_B = 1.380649e-23 J/K

With a pseudo-random generator seeded EXACTLY at seed = 41, i.e. NumPy's default_rng(41) Generator (PCG64); every draw this procedure consumes comes from that one generator IN ORDER, and any percentile it asks for is NumPy's np.percentile with default linear interpolation -- Python's random.Random, legacy np.random.seed (MT19937) or torch streams do NOT reproduce this procedure, draw 10000
independent samples -- D ~ Normal(4.2e-7, 0.9e-7), eta ~ Normal(1.02,
0.02), T ~ Normal(298.15, 0.25), same sample index for all three -- and
evaluate R_h at each. Across those 10000 values, give the mean, sample
standard deviation (ddof=1), and the 2.5th/97.5th percentiles. With
skew = (mean - median)/sd over the same 10000 values, call the result
"right-skewed" above skew = +0.05, "left-skewed" below -0.05, else
"symmetric".

R_h is proportional to 1/D, a convex relationship for D > 0, so this
distribution is NOT symmetric even though both D and eta are drawn
from symmetric normals -- a linearized error bar misses this, and only
the real 10000-draw simulation from seed 41 captures it."""
    )
    items = [
        ("output_mean", "output_mean: mean of the 10000 R_h samples, nm"),
        ("output_sd", "output_sd: sample standard deviation (ddof=1), nm"),
        ("ci_low", "ci_low: 2.5th percentile, nm"),
        ("ci_high", "ci_high: 97.5th percentile, nm"),
        ("skew_label", "skew_label: one of [\"left-skewed\", \"symmetric\", \"right-skewed\"]"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-08"] = dict(
        id="X-08", family="G2", levers=["L12", "L1"],
        axes={"arith": 4, "reason": 1},
        mechanism="Same genuine 10000-draw Monte Carlo propagation as X-05, a second independent condition through the same convex Stokes-Einstein formula.",
        refs=["new-pattern-X12-X15-X19"], tools_A=["propagate_uncertainty_chain"],
        fields={"output_mean": num(res.output_mean, 0.01, "nm"),
                "output_sd": num(res.output_sd, 0.02, "nm"),
                "ci_low": num(res.ci_low, 0.02, "nm"),
                "ci_high": num(res.ci_high, 0.02, "nm"),
                "skew_label": exact(label)},
        body=body)


# ============================================================ X-09
def build_x09():
    rng = random.Random(3009)
    true_x0, true_w, true_lo, true_hi = math.log10(1.4), 0.12, 71.2, 30.5
    cs = [0.25, 0.4, 0.6, 0.85, 1.1, 1.4, 1.8, 2.3, 3.0, 4.0]
    xs = [round(math.log10(c), 4) for c in cs]
    ys = []
    for x in xs:
        val = true_lo + (true_hi - true_lo) * 0.5 * (1 + math.tanh((x - true_x0) / true_w))
        ys.append(round(val + rng.gauss(0, 0.7), 2))
    res = curve.bootstrap_cmc(xs, ys, n_boot=1000, ci=0.95, seed=27)
    ci_width_pct = 100.0 * (res.cmc_ci_high - res.cmc_ci_low) / res.cmc_point
    rows = "\n".join(f"    {c:>5.2f}      {y:>6.2f}" for c, y in zip(cs, ys))
    body = (
        qhead("X-09", "A fifth residual-bootstrap CMC, the tightest concentration spacing here") +
        """A dye-probe transition curve with closely-spaced low concentrations,
same probe direction as every other such question (higher
concentration -> lower plateau):

    c (mM)     value
""" + rows + """

The algorithm is the one used throughout this set: fit the
4-parameter sigmoid value(x) = lo + (hi - lo)*0.5*(1 + tanh((x -
x0)/w)) with x = log10(c in mM) by nonlinear least squares, giving
CMC = 10^x0. Take the residuals of that fit; with a generator seeded
EXACTLY at seed = 27, using NumPy's default_rng(seed) PCG64 bit generator as the single generator consumed in a fixed order (the order the reference implementation draws its resample indices), and percentiles computed with np.percentile's default linear interpolation -- Python's random.Random, legacy np.random.seed (MT19937) or torch streams do NOT reproduce this procedure, resample those residuals WITH replacement 1000
times, add each resample back onto the fit, and refit each time to
build 1000 bootstrap x0_b values. Their 2.5th/97.5th percentiles,
converted via 10^(...), give the 95% CMC interval, and
ci_width_pct = 100*(CMC_high - CMC_low)/CMC_point.

Nothing here reduces to a formula you can evaluate directly -- the
1000 specific resamples are fixed by seed 27, and only actually running
them reproduces the interval. Report whether ci_width_pct exceeds 2.0
(true/false)."""
    )
    items = [
        ("cmc_point", "cmc_point: point CMC estimate in mM"),
        ("cmc_ci_low", "cmc_ci_low: 2.5th-percentile bootstrap CMC in mM"),
        ("cmc_ci_high", "cmc_ci_high: 97.5th-percentile bootstrap CMC in mM"),
        ("ci_width_pct", "ci_width_pct: 100*(cmc_ci_high-cmc_ci_low)/cmc_point"),
        ("ci_wide", "ci_wide: true if ci_width_pct > 2.0, else false"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-09"] = dict(
        id="X-09", family="G2", levers=["L12", "L8"],
        axes={"arith": 4, "reason": 1},
        mechanism="Same genuine residual-bootstrap mechanism as every other bootstrap question here, a fifth independent dataset with the tightest spacing in the set.",
        refs=["new-pattern-X12-X15-X19"], tools_A=["bootstrap_cmc"],
        fields={"cmc_point": num(res.cmc_point, 0.01, "mM"),
                "cmc_ci_low": num(res.cmc_ci_low, 0.02, "mM"),
                "cmc_ci_high": num(res.cmc_ci_high, 0.02, "mM"),
                "ci_width_pct": num(ci_width_pct, 0.02, "%"),
                "ci_wide": boolean(ci_width_pct > 2.0)},
        body=body)


# ============================================================ X-10
def build_x10():
    res = kin.micelle_relaxation_odes(
        c_total_M=0.022, K_initial_per_M=6000.0, K_final_per_M=2200.0,
        k_on_per_M_s=4.0e7, n_max=26)
    body = (
        qhead("X-10", "A second stiff master-equation relaxation, K decreasing this time") +
        """A micellar system at c_total = 0.022 M undergoes a temperature jump
where K governing every C_n + C_1 <-> C_{n+1} step (n = 1..26) DROPS
from K_initial = 6000 /M to K_final = 2200 /M at t = 0 (the opposite
direction from any other relaxation question here), with k_on =
4.0e7 /M/s fixed (k_off = k_on/K, detailed balance).

The method is identical to any other relaxation-ODE question in this
set: build the 26-state master equation using K_final and k_on for
every rate, initialize at the EXACT equilibrium under K_initial (via
closed-form bisection on the single-K mass-action equilibrium), then
integrate forward with a stiff solver (this mixes a fast monomer
timescale with slow reorganization -- a naive fixed-step integrator
will not converge; a BDF-type method with the analytic Jacobian is
required) until the new equilibrium under K_final is reached.

The new equilibrium's monomer concentration is the new kinetic CMC
(kinetic_cmc_M); its mean aggregation number is
mean_aggregation_number. Linearizing around that new equilibrium and
taking the Jacobian's eigenvalues gives the relaxation spectrum; the
slowest nonzero mode's time constant is tau_dominant_s. A correct
stiff solve keeps the integrator's own max relative mass-conservation
error below 1e-8 throughout.

tau_dominant_s in particular has no closed form here -- it comes from
an eigenvalue decomposition of a 26x26 Jacobian evaluated at a
numerically-integrated equilibrium. Report whether the integration
converged, choosing from ["converged", "not_converged"]."""
    )
    items = [
        ("kinetic_cmc_M", "kinetic_cmc_M: new equilibrium monomer concentration, mol/L"),
        ("tau_dominant_s", "tau_dominant_s: slowest relaxation time constant, s"),
        ("mean_aggregation_number", "mean_aggregation_number: new equilibrium mean Nagg"),
        ("mass_conserved", "mass_conserved: true if the max relative mass-conservation error is below 1e-8, else false"),
        ("convergence", "convergence: one of [\"converged\", \"not_converged\"]"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-10"] = dict(
        id="X-10", family="G2", levers=["L12", "L1"],
        axes={"arith": 3, "reason": 2},
        mechanism="Same genuine stiff-ODE mechanism as X-13, with K decreasing instead of increasing and a different molecule-count/rate regime -- checks the pattern is not direction-specific.",
        refs=["new-pattern-X12-X15-X19"], tools_A=["micelle_relaxation_odes"],
        fields={"kinetic_cmc_M": num(res.kinetic_cmc_M, 0.01, "mol/L"),
                "tau_dominant_s": num(res.tau_dominant_s, 0.02, "s"),
                "mean_aggregation_number": num(res.mean_aggregation_number, 0.01),
                "mass_conserved": boolean(res.max_relative_mass_error < 1e-8),
                "convergence": exact("converged")},
        body=body)


def build():
    for f in (build_x01, build_x02, build_x03, build_x04, build_x05,
              build_x06, build_x07, build_x08, build_x09, build_x10):
        f()
    return QUESTIONS
