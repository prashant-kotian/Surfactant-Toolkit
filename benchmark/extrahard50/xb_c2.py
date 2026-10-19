"""EXTRA50 chunk 2 — G3 model races (X-11..15) + G4 transformed coordinates (X-16..20)."""
from __future__ import annotations
import math
import random

from surfactantkit import thermodynamics as thermo
from surfactantkit import adsorption as ads
from surfactantkit import curve_analysis as curve
from surfactantkit import mixed_micelle as mm
from surfactantkit import mixture_model_selection as mms
from surfactantkit import hld as hld_mod
from surfactantkit import electrostatics as elec
from surfactantkit import kinetics as kin

from xb_common import R_GAS, WATER_M, num, integer, exact, boolean, outputs_block, qhead

QUESTIONS: dict[str, dict] = {}


def _X(c_mM: float) -> float:
    return (c_mM / 1000.0) / (c_mM / 1000.0 + WATER_M)


def _gen_rubingh(alpha: float, beta: float, c1: float, c2: float) -> float:
    """Forward Rubingh: given alpha and true beta, return Cmix (same unit as c1,c2).
    Solve ln[alpha*(1-x)*c2/((1-alpha)*x*c1)] = beta*(1-2x) for x in (0,1),
    then Cmix = (x*c1/alpha)*exp(beta*(1-x)^2)."""
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


# ============================================================ X-11
def build_x11():
    res = kin.kmc_micelle_kinetics(
        c_total_M=0.012, K_per_M=4200.0, k_on_per_M_s=8.0e6,
        n_max=25, n_molecules=400, n_samples=60, seed=31)
    body = (
        qhead("X-11", "Gillespie kinetic Monte Carlo of stepwise micellization") +
        """A well-mixed micellization system is modeled at the molecular level:
stepwise addition/loss of single monomers, C_n + C_1 <-> C_{n+1}, fired
by the direct Gillespie stochastic simulation algorithm (SSA):

    c_total = 0.012 M
    equilibrium constant K = 4200 /M (same K governs every n)
    on-rate k_on = 8.0e6 /M/s;  off-rate k_off = k_on / K (detailed balance)
    max aggregate size tracked n_max = 25
    molecule count in the simulated volume = 400
    starts from a Poisson sample at equilibrium under K, then relaxes

Pinned procedure (authoritative -- this is a real stochastic simulation,
not a formula):

1. Using a pseudo-random generator seeded EXACTLY with seed = 31, i.e. NumPy's default_rng(31) Generator (PCG64); every draw this procedure consumes comes from that one generator IN ORDER, and any percentile it asks for is NumPy's np.percentile with default linear interpolation -- Python's random.Random, legacy np.random.seed (MT19937) or torch streams do NOT reproduce this procedure,
   construct the Poisson-equilibrium initial state under K at
   c_total = 0.012 M for 400 molecules in a volume sized so the expected
   molecule count equals 400.
2. Run the direct SSA: at each step, compute all C_n+C_1 forward
   (association) and C_{n+1} backward (dissociation) propensities,
   draw an exponential waiting time and a reaction channel from those
   propensities (both draws from the SAME seeded generator, consumed in
   sequence), fire the chosen reaction, update state, repeat.
3. Sample the running state at 60 evenly-spaced checkpoints in time.
   Take the FINAL 20% of checkpoints as the stationary window and
   average the monomer concentration over it (equilibrium_c1_M) and the
   mean aggregation number over aggregates with n >= 2
   (equilibrium_mean_nagg).
4. Fit a single exponential decay to the monomer trajectory to get a
   relaxation time tau_fit_s.
5. Exact integer mass conservation: sum(n * count_n) over all n must be
   EXACTLY invariant across every step (SSA conserves molecule count by
   construction) -- report the number of simulated reaction events
   fired and the max observed deviation from exact conservation
   (should be exactly 0 for a correct run).

There is no closed-form shortcut for steps 1-4: the exact sequence of
reaction events depends on the specific pseudo-random draws from the
stated seed, and only running the actual SSA reproduces the stated
quantities to the required precision. Report whether the final window
was judged statistically stationary, choosing from
["converged", "not_converged"]."""
    )
    items = [
        ("equilibrium_c1_M", "equilibrium_c1_M: stationary-window mean monomer concentration, mol/L"),
        ("equilibrium_mean_nagg", "equilibrium_mean_nagg: stationary-window mean aggregation number (n>=2)"),
        ("tau_fit_s", "tau_fit_s: exponential relaxation time fit to the monomer trajectory, s"),
        ("mass_conservation_max_error", "mass_conservation_max_error: max deviation from exact mass conservation (integer, step 5)"),
        ("convergence", "convergence: one of [\"converged\", \"not_converged\"]"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-11"] = dict(
        id="X-11", family="G3", levers=["L12", "L8"],
        axes={"arith": 4, "reason": 1},
        mechanism="Gold is a seeded Gillespie SSA trajectory -- the exact reaction sequence depends on specific pseudo-random draws; no closed-form or mental simulation reproduces equilibrium_c1_M/tau_fit_s to the required precision.",
        refs=["new-pattern-X12-X15-X19"], tools_A=["kmc_micelle_kinetics"],
        fields={"equilibrium_c1_M": num(res.equilibrium_c1_M, 0.02, "mol/L"),
                "equilibrium_mean_nagg": num(res.equilibrium_mean_nagg, 0.02),
                "tau_fit_s": num(res.tau_fit_s, 0.02, "s"),
                "mass_conservation_max_error": integer(res.mass_conservation_max_error),
                "convergence": exact("converged")},
        body=body)


# ============================================================ X-12
def build_x12():
    rng = random.Random(1012)
    c1, c2, beta_true = 1.05, 8.60, -2.1
    alphas = [0.1, 0.3, 0.5, 0.7, 0.9]
    cmix = []
    for a in alphas:
        v = _gen_rubingh(a, beta_true, c1, c2)
        cmix.append(round(v * (1 + rng.gauss(0, 0.006)), 2))
    sel = mms.select_mixture_model(alphas, cmix, c1, c2)
    cand = {c.name: c for c in sel.candidates}
    beta_reg = mm.rubingh_beta_regression(alphas, cmix, c1, c2)
    bic_rub = cand["rubingh_regular_solution"].bic
    bic_eom = cand["eommm_asymmetric_margules"].bic
    # verify (not assume) this dataset's EOMMM fit is unique at the default
    # margin select_mixture_model uses -- eommm_global_fit's own docstring
    # discloses non-uniqueness is possible for SOME datasets; check this one.
    multi = mm.eommm_global_fit_all_optima(alphas, cmix, c1, c2)
    assert not multi.degenerate, f"X-12: non-unique EOMMM fit found ({multi.n_optima_found} optima) -- grading must list all, not assume one"
    rows = "\n".join(f"    {a:>4.2f}   {m:>7.2f}" for a, m in zip(alphas, cmix))
    body = (
        qhead("X-12", "Three-way mixture-model race (Clint / Rubingh / EOMMM)") +
        """Mixed-micelle CMC series for a binary pair with pure-component CMCs
CMC1 = 1.05 mM, CMC2 = 8.60 mM:

    alpha1   CMC_mix (mM)
    ---------------------
""" + rows + """

Pinned model-selection protocol (authoritative):

- Candidate models and their parameter counts:
      clint_ideal                k = 0
      rubingh_regular_solution   k = 1   (the beta parameter)
      eommm_asymmetric_margules  k = 2   (W12, W21)
- For EACH candidate: fit it to the 5 points (least squares of residuals
  in LN of CMC_mix), get RSS in that ln space, then
      BIC = n * ln(RSS/n) + k * ln(n),  n = 5
- Select the lowest BIC among the three.
- Independently of the race, report the pointwise-averaged Rubingh beta
  for the series:

      x1 solved from:  x1^2 * ln(alpha1*CMCmix/(x1*CMC1))
                     - (1-x1)^2 * ln((1-alpha1)*CMCmix/((1-x1)*CMC2)) = 0
      beta_i = ln(alpha1*CMCmix/(x1*CMC1)) / (1-x1)^2
      beta_mean = mean of beta_i over the points that solve

Report BIC of rubingh and of eommm, delta = BIC_eommm - BIC_rubingh,
the selected model id (choose from ["clint_ideal",
"rubingh_regular_solution", "eommm_asymmetric_margules"]), and
beta_mean."""
    )
    items = [
        ("bic_rubingh", "bic_rubingh: BIC of rubingh_regular_solution"),
        ("bic_eommm", "bic_eommm: BIC of eommm_asymmetric_margules"),
        ("delta_bic", "delta_bic: BIC_eommm - BIC_rubingh"),
        ("selected", "selected: one of [\"clint_ideal\", \"rubingh_regular_solution\", \"eommm_asymmetric_margules\"]"),
        ("beta_mean", "beta_mean: pointwise-averaged Rubingh beta (dimensionless)"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-12"] = dict(
        id="X-12", family="G3", levers=["L3", "L6", "L11"],
        axes={"arith": 4, "reason": 1},
        mechanism="Hand-fitting three models (incl. implicit Rubingh root-find) in ln space + BIC + separate beta regression: execution explosion; round-1 BIC races all failed numerically. EOMMM fit verified unique for this dataset via eommm_global_fit_all_optima (n_optima_found=1) -- not assumed.",
        refs=["E4", "E1", "H-06"], tools_A=["select_mixture_model", "rubingh_beta_regression"],
        fields={"bic_rubingh": num(bic_rub, 0.02),
                "bic_eommm": num(bic_eom, 0.02),
                "delta_bic": num(bic_eom - bic_rub, 0.02),
                "selected": exact(sel.selected_model),
                "beta_mean": num(beta_reg.beta_mean, 0.015)},
        body=body)


# ============================================================ X-13
def build_x13():
    res = kin.micelle_relaxation_odes(
        c_total_M=0.018, K_initial_per_M=2500.0, K_final_per_M=9500.0,
        k_on_per_M_s=3.0e7, n_max=28)
    body = (
        qhead("X-13", "Stiff master-equation relaxation after a temperature jump") +
        """A micellar system at total concentration c_total = 0.018 M undergoes a
sudden temperature (T-)jump: the equilibrium constant K governing EVERY
stepwise association C_n + C_1 <-> C_{n+1} (n = 1..28) switches from
K_initial = 2500 /M to K_final = 9500 /M at t = 0, with on-rate
k_on = 3.0e7 /M/s (off-rate k_off = k_on/K, detailed balance) held fixed
throughout.

Pinned procedure (authoritative -- a real stiff ODE integration, not a
formula):

1. Build the 28-state master equation dy/dt = A(y)*y for the aggregate
   concentrations y_1..y_28 (y_n = concentration of n-mers), using
   K_final and k_on for every forward/backward rate.
2. Set the initial condition y(0) to the EXACT equilibrium distribution
   under K_initial at c_total = 0.018 M (found by closed-form bisection
   on the single-K mass-action equilibrium).
3. Integrate the stiff ODE system forward in time (this system mixes a
   fast monomer-association timescale with a slow micelle-reorganization
   timescale -- a stiff BDF-type integrator with the analytic Jacobian
   is required; a naive fixed-step integrator will not converge) until
   it reaches the NEW equilibrium under K_final.
4. The new equilibrium's monomer concentration IS the new kinetic CMC
   (kinetic_cmc_M). Its mean aggregation number is
   mean_aggregation_number.
5. Linearize the ODE system around the new equilibrium and take the
   Jacobian's eigenvalues; the slowest nonzero mode's time constant is
   tau_dominant_s -- the experimentally measurable relaxation time of
   the jump.
6. Total mass sum(n*y_n) must be conserved exactly by the true ODE
   solution; a correct stiff solve keeps the integrator's own max
   relative mass-conservation error below 1e-8 throughout.

There is no closed-form shortcut for steps 1-5: tau_dominant_s comes
from an eigenvalue decomposition of a 28x28 Jacobian evaluated at a
numerically-integrated equilibrium, not from any single formula.
Report whether the integration was judged converged, choosing from
["converged", "not_converged"]."""
    )
    items = [
        ("kinetic_cmc_M", "kinetic_cmc_M: new equilibrium monomer concentration, mol/L (step 4)"),
        ("tau_dominant_s", "tau_dominant_s: slowest relaxation time constant, s (step 5)"),
        ("mean_aggregation_number", "mean_aggregation_number: new equilibrium mean Nagg (step 4)"),
        ("mass_conserved", "mass_conserved: true if the max relative mass-conservation error is below 1e-8, else false (step 6)"),
        ("convergence", "convergence: one of [\"converged\", \"not_converged\"]"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-13"] = dict(
        id="X-13", family="G3", levers=["L12", "L1"],
        axes={"arith": 3, "reason": 2},
        mechanism="Gold requires a real stiff-ODE integration (28-state master equation, BDF + analytic Jacobian) and an eigenvalue decomposition for tau_dominant_s -- deterministic but not hand-derivable from any closed-form expression.",
        refs=["new-pattern-X12-X15-X19"], tools_A=["micelle_relaxation_odes"],
        fields={"kinetic_cmc_M": num(res.kinetic_cmc_M, 0.01, "mol/L"),
                "tau_dominant_s": num(res.tau_dominant_s, 0.02, "s"),
                "mean_aggregation_number": num(res.mean_aggregation_number, 0.01),
                "mass_conserved": boolean(res.max_relative_mass_error < 1e-8),
                "convergence": exact("converged")},
        body=body)


# ============================================================ X-14
def build_x14():
    res = thermo.mc_aggregation_distribution(
        n_monomers=35, temperature_K=298.15, epsilon_kJ_per_mol=7.0,
        n_steps=80000, seed=44)
    body = (
        qhead("X-14", "Canonical Metropolis Monte Carlo of an associating cluster system") +
        """35 labeled monomers sit in a box; a cluster of size n has energy
E(n) = -epsilon*(n-1) (one bond per additional monomer):

    n_monomers = 35
    temperature_K = 298.15
    epsilon = 7.0 kJ/mol  (per-bond binding energy)
    n_steps = 80000 Metropolis moves
    burn-in fraction = 0.1 (default), thinning = 10 (default)

Pinned procedure (authoritative -- a real canonical Monte Carlo
simulation, not a formula):

1. Using a pseudo-random generator seeded EXACTLY with seed = 44, i.e. NumPy's default_rng(44) Generator (PCG64); every draw this procedure consumes comes from that one generator IN ORDER, and any percentile it asks for is NumPy's np.percentile with default linear interpolation -- Python's random.Random, legacy np.random.seed (MT19937) or torch streams do NOT reproduce this procedure, start
   from all 35 monomers free (each its own cluster of size 1).
2. Repeat 80000 times: propose one of three moves chosen uniformly at
   random -- detach a random monomer to become free, attach a random
   free-or-clustered monomer to a random other cluster, or transfer a
   random monomer directly between two clusters -- then accept or
   reject the move with the Metropolis criterion on the bond-energy
   change beta*dE, beta = 1/(R*T/1000) in units consistent with
   epsilon in kJ/mol (each accept/reject draw also consumes the SAME
   seeded generator, in sequence with the move proposals).
3. Discard the first 10% of steps as burn-in, then keep every 10th
   remaining step. From the retained samples, build the histogram of
   the fraction of monomers residing in clusters of each size n.
4. From that histogram, report: the fraction of monomers that are free
   monomers (n=1), the monomer-weighted mean cluster size
   (sum of n * (fraction of monomers in size-n clusters)), and the
   largest cluster size observed anywhere in the retained samples.
5. Since n_monomers = 35 <= 45, the EXACT Boltzmann histogram over all
   integer partitions of 35 is also computable in closed form; report
   the maximum absolute deviation between your Monte Carlo histogram
   and this exact solution (a genuine verification step, not a free
   parameter).

There is no shortcut for steps 2-4: the specific sequence of accepted
and rejected moves depends on the exact pseudo-random draws from the
stated seed, and only running the actual 80000-step chain reproduces
these statistics. Report whether the first and second half of the
chain agreed well enough to call the run converged, choosing from
["converged", "not_converged"]."""
    )
    items = [
        ("free_monomer_fraction", "free_monomer_fraction: fraction of monomers in size-1 clusters (step 4)"),
        ("mean_cluster_size", "mean_cluster_size: monomer-weighted mean cluster size (step 4)"),
        ("largest_cluster", "largest_cluster: largest cluster size observed in the retained samples (step 4, integer)"),
        ("exact_max_deviation", "exact_max_deviation: max abs deviation from the exact partition solution (step 5)"),
        ("convergence", "convergence: one of [\"converged\", \"not_converged\"]"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-14"] = dict(
        id="X-14", family="G3", levers=["L12", "L8"],
        axes={"arith": 4, "reason": 1},
        mechanism="Gold is a seeded 80000-step canonical Metropolis MC chain over cluster configurations -- the retained-sample histogram depends on the exact accept/reject sequence, not derivable analytically except via the (separately verified) exact-partition cross-check.",
        refs=["new-pattern-X12-X15-X19"], tools_A=["mc_aggregation_distribution"],
        fields={"free_monomer_fraction": num(res.free_monomer_fraction, 0.02),
                "mean_cluster_size": num(res.mean_cluster_size_monomer_weighted, 0.02),
                "largest_cluster": integer(res.largest_cluster_observed),
                "exact_max_deviation": num(res.exact_max_abs_deviation, 0.02),
                "convergence": exact("converged")},
        body=body)


# ============================================================ X-15
def build_x15():
    rng = random.Random(1015)
    c1, c2, beta_true = 1.05, 8.60, -1.5
    alphas = [0.1, 0.3, 0.5, 0.7, 0.9]
    cmix = [round(_gen_rubingh(a, beta_true, c1, c2) * (1 + rng.gauss(0, 0.015)), 2)
            for a in alphas]
    margin = mm.eommm_find_minimal_feasible_margin(alphas, cmix, c1, c2)
    # verify (not assume) uniqueness at the minimal feasible margin just found
    multi = mm.eommm_global_fit_all_optima(alphas, cmix, c1, c2, cmc_margin=margin.cmc_margin_used,
                                            feasible_tolerance=1e-3)
    assert not multi.degenerate, f"X-15: non-unique EOMMM fit found at the minimal margin ({multi.n_optima_found} optima)"
    rod = mm.rodenas_x1_series(alphas, cmix)
    i = list(rod.alpha1_used).index(0.5)
    x1_rod = rod.x1_rodenas[i]
    x1_rub = mm.solve_rubingh_x(0.5, cmix[2], c1, c2)
    agree = abs(x1_rod - x1_rub) <= 0.03
    feas = margin.cmc_margin_used <= 0.01
    rows = "\n".join(f"    {a:>4.2f}   {m:>7.2f}" for a, m in zip(alphas, cmix))
    body = (
        qhead("X-15", "EOMMM feasibility margin + x1 method race") +
        """Binary mixture series, CMC1 = 1.05 mM, CMC2 = 8.60 mM:

    alpha1   CMC_mix (mM)
    ---------------------
""" + rows + """

Pinned instructions (authoritative):

1. Fit the asymmetric-Margules EOMMM model (two interaction parameters
   W12, W21) by least squares, allowing BOTH CMCs to relax by a common
   multiplicative margin; find the SMALLEST margin (as a fraction of the
   CMC values, e.g. 0.06 = 6%) for which a fully feasible fit exists
   (total infeasibility <= 0.001). Report that minimal margin.
   (This is a search — report the converged value to ~2 significant
   figures; the intended answer is a small nonzero number.)
2. feasible_within_1pct = (margin <= 0.01), true/false.
3. Rodenas micellar mole fraction from the local slope of ln(CMC_mix)
   vs alpha1 at alpha1 = 0.50, using the central finite difference
   over the points alpha1 = 0.30, 0.50, 0.70:

       x1_Rod = -(1 - alpha1)*alpha1 * dln(CMCmix)/d(alpha1) + alpha1

4. Rubingh x1 at alpha1 = 0.50: solve
       x1^2*ln(alpha1*CMCmix/(x1*CMC1)) - (1-x1)^2*
       ln((1-alpha1)*CMCmix/((1-x1)*CMC2)) = 0
   for x1 in (0, 1).
5. agree = |x1_Rod - x1_Rub| <= 0.03, true/false."""
    )
    items = [
        ("margin", "margin: minimal feasible CMC margin (fraction, step 1)"),
        ("feasible_within_1pct", "feasible_within_1pct: true/false (step 2)"),
        ("x1_rodenas", "x1_rodenas: Rodenas x1 at alpha1 = 0.50 (step 3)"),
        ("x1_rubingh", "x1_rubingh: Rubingh x1 at alpha1 = 0.50 (step 4)"),
        ("agree_003", "agree_003: true/false (step 5)"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-15"] = dict(
        id="X-15", family="G3", levers=["L3", "L6", "L12"],
        axes={"arith": 3, "reason": 2},
        mechanism="Two implicit root-finds + margin search + finite-difference derivative: computationally out of unaugmented reach (abandonment), and the two x1 methods historically disagree in models' hands. EOMMM fit verified unique at the minimal feasible margin via eommm_global_fit_all_optima (n_optima_found=1) -- not assumed.",
        refs=["E1", "E4", "H-10"], tools_A=["eommm_find_minimal_feasible_margin", "rodenas_x1_from_series", "rubingh_solve"],
        fields={"margin": num(margin.cmc_margin_used, 0.02),
                "feasible_within_1pct": boolean(feas),
                "x1_rodenas": num(x1_rod, 0.01),
                "x1_rubingh": num(x1_rub, 0.01),
                "agree_003": boolean(agree)},
        body=body)


# ============================================================ X-16
def build_x16():
    rng = random.Random(4016)
    alphas = [0.15, 0.35, 0.55, 0.75, 0.95]
    c1, c2, beta_true = 1.20, 7.40, 2.6

    def _gen(a, beta, cc1, cc2):
        def f(x):
            return math.log(a * (1 - x) * cc2 / ((1 - a) * x * cc1)) - beta * (1 - 2 * x)
        lo, hi = 1e-9, 1 - 1e-9
        flo = f(lo)
        for _ in range(200):
            mid = (lo + hi) / 2
            fm = f(mid)
            if (flo < 0) == (fm < 0):
                lo, flo = mid, fm
            else:
                hi = mid
        x = (lo + hi) / 2
        return (x * cc1 / a) * math.exp(beta * (1 - x) ** 2)
    cmix = [round(_gen(a, beta_true, c1, c2) * (1 + rng.gauss(0, 0.006)), 2) for a in alphas]
    fit = mm.eommm_global_fit(alphas, cmix, c1, c2)
    feasible = fit.total_infeasibility < 1e-6
    multi = mm.eommm_global_fit_all_optima(alphas, cmix, c1, c2)
    assert not multi.degenerate, f"X-16: non-unique EOMMM fit found ({multi.n_optima_found} optima)"
    rows = "\n".join(f"    {a:>4.2f}   {m:>6.2f}" for a, m in zip(alphas, cmix))
    body = (
        qhead("X-16", "EOMMM asymmetric-Margules global fit: W12/W21 by 2D grid search") +
        """A binary surfactant pair has pure-component CMCs CMC1 = 1.20 mM,
CMC2 = 7.40 mM, and the mixed-micelle series below:

    alpha1   CMC_mix (mM)
    ---------------------
""" + rows + """

Pinned model-selection protocol (authoritative). Unlike an ordinary
regression, EOMMM does not minimize a least-squares residual against
CMC_mix at all: at a fixed cmc_margin = 0.10 band around each point's
own experimental CMC_mix, it scores every candidate (W12, W21) by how
badly that pair -- together with each point's own micellar mole
fraction x1 -- violates the underlying mass-balance constraints there
(its "infeasibility"), and reports whichever candidate drives that
score closest to zero. That minimum has no closed form and is not
gradient-friendly everywhere, so it is located by a 2D search: first a
coarse 61x61 pass over (W12, W21) in [-30, 30]^2, then a fine 41x41
pass confined to the best coarse cell's neighborhood.

Report W12, W21, the fit's r_squared (RSS in ln(CMC_mix) space against
the fitted cmc_var values versus the data), and the micellar mole
fraction x1 at the FIRST point (alpha1 = 0.15). Finally, report whether
the fit's total_infeasibility fell below 1e-6 (fully feasible at the
stated margin), choosing true or false."""
    )
    items = [
        ("W12", "W12: fitted asymmetric-Margules interaction parameter"),
        ("W21", "W21: fitted asymmetric-Margules interaction parameter"),
        ("r_squared", "r_squared: fit quality in ln(CMC_mix) space"),
        ("x1_first_point", "x1_first_point: micellar mole fraction x1 at alpha1 = 0.15"),
        ("feasible", "feasible: true if total_infeasibility < 1e-6, else false"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-16"] = dict(
        id="X-16", family="G4", levers=["L12", "L6"],
        axes={"arith": 4, "reason": 1},
        mechanism="Gold is a genuine 2D (coarse+fine) grid search over a non-closed-form infeasibility objective -- the same mechanism proven to universally fail on X-12, applied here directly as the primary task rather than inside a BIC race. Fit verified unique for this dataset via eommm_global_fit_all_optima (n_optima_found=1) -- not assumed.",
        refs=["new-pattern-X12-X15-X19"], tools_A=["eommm_global_fit"],
        fields={"W12": num(fit.W12, 0.02),
                "W21": num(fit.W21, 0.02),
                "r_squared": num(fit.r_squared, 0.01),
                "x1_first_point": num(fit.x1_values[0], 0.01),
                "feasible": boolean(feasible)},
        body=body)


# ============================================================ X-17
def build_x17():
    rng = random.Random(2017)
    true_x0, true_w, true_lo, true_hi = math.log10(2.1), 0.15, 70.5, 32.0
    cs = [0.35, 0.55, 0.85, 1.2, 1.6, 2.0, 2.5, 3.2, 4.0, 5.2]
    xs = [round(math.log10(c), 4) for c in cs]
    ys = []
    for x in xs:
        val = true_lo + (true_hi - true_lo) * 0.5 * (1 + math.tanh((x - true_x0) / true_w))
        ys.append(round(val + rng.gauss(0, 0.9), 2))
    res = curve.bootstrap_cmc(xs, ys, n_boot=1000, ci=0.95, seed=6)
    ci_width_pct = 100.0 * (res.cmc_ci_high - res.cmc_ci_low) / res.cmc_point
    rows = "\n".join(f"    {c:>5.2f}      {y:>6.2f}" for c, y in zip(cs, ys))
    body = (
        qhead("X-17", "A third residual-bootstrap CMC, tighter noise and a different seed") +
        """A dye-probe transition curve for a different ionic surfactant (same
direction as a standard tension/dye-probe curve: higher concentration
-> lower plateau):

    c (mM)     value
""" + rows + """

Pinned procedure (authoritative; the algorithm matches every other
residual-bootstrap CMC question in this set, only the data/seed move):

1. Let x = log10(c in mM). Fit the 4-parameter sigmoidal transition
       value(x) = lo + (hi - lo) * 0.5 * (1 + tanh((x - x0) / w))
   to the 10 points above by nonlinear least squares (iterative solver
   required; no closed form). CMC = 10^x0.
2. Take this fit's own residuals at the 10 points.
3. With a pseudo-random generator seeded EXACTLY with seed = 6, i.e. NumPy's default_rng(6) Generator (PCG64); every draw this procedure consumes comes from that one generator IN ORDER, and any percentile it asks for is NumPy's np.percentile with default linear interpolation -- Python's random.Random, legacy np.random.seed (MT19937) or torch streams do NOT reproduce this procedure, draw
   1000 bootstrap resamples of the residuals WITH replacement, add each
   resample back onto the fitted curve, and refit to get a bootstrap
   x0_b each time.
4. From the 1000 bootstrap x0_b values, take the 2.5th and 97.5th
   percentiles and convert via 10^(...) for the 95% CMC interval.
5. ci_width_pct = 100 * (CMC_high - CMC_low) / CMC_point.

As always, steps 3-4 require actually running 1000 real nonlinear
refits from the stated seed -- there is no formula shortcut for a
percentile bootstrap interval. Report whether ci_width_pct exceeds 3.0
(true/false on the number from step 5)."""
    )
    items = [
        ("cmc_point", "cmc_point: point CMC estimate in mM (step 1)"),
        ("cmc_ci_low", "cmc_ci_low: 2.5th-percentile bootstrap CMC in mM (step 4)"),
        ("cmc_ci_high", "cmc_ci_high: 97.5th-percentile bootstrap CMC in mM (step 4)"),
        ("ci_width_pct", "ci_width_pct: 100*(cmc_ci_high-cmc_ci_low)/cmc_point (step 5)"),
        ("ci_wide", "ci_wide: true if ci_width_pct > 3.0, else false"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-17"] = dict(
        id="X-17", family="G4", levers=["L12", "L8"],
        axes={"arith": 4, "reason": 1},
        mechanism="Same genuine-bootstrap mechanism as X-01/X-02 on a third independent dataset/seed -- tests whether the pattern generalizes rather than a model having memorized one specific interval.",
        refs=["new-pattern-X12-X15-X19"], tools_A=["bootstrap_cmc"],
        fields={"cmc_point": num(res.cmc_point, 0.01, "mM"),
                "cmc_ci_low": num(res.cmc_ci_low, 0.02, "mM"),
                "cmc_ci_high": num(res.cmc_ci_high, 0.02, "mM"),
                "ci_width_pct": num(ci_width_pct, 0.02, "%"),
                "ci_wide": boolean(ci_width_pct > 3.0)},
        body=body)


# ============================================================ X-18
def build_x18():
    res = thermo.mc_aggregation_distribution(
        n_monomers=28, temperature_K=308.15, epsilon_kJ_per_mol=6.2,
        n_steps=70000, seed=55)
    body = (
        qhead("X-18", "A second canonical Metropolis Monte Carlo cluster simulation") +
        """28 labeled monomers, a different temperature and bond strength from
any other Monte Carlo question in this set:

    n_monomers = 28
    temperature_K = 308.15
    epsilon = 6.2 kJ/mol  (per-bond binding energy)
    n_steps = 70000 Metropolis moves
    burn-in fraction = 0.1 (default), thinning = 10 (default)

Pinned procedure (authoritative; this matches the canonical-MC
algorithm used throughout this set -- only the physical inputs move):

1. With a pseudo-random generator seeded EXACTLY with seed = 55, i.e. NumPy's default_rng(55) Generator (PCG64); every draw this procedure consumes comes from that one generator IN ORDER, and any percentile it asks for is NumPy's np.percentile with default linear interpolation -- Python's random.Random, legacy np.random.seed (MT19937) or torch streams do NOT reproduce this procedure, start
   from all 28 monomers free.
2. Repeat 70000 times: propose a detach, attach, or transfer move
   (uniformly at random among the three), accept or reject by the
   Metropolis criterion on beta*dE with beta = 1/(R*T/1000) and
   dE = -epsilon*(bonds after - bonds before) (both the move choice and
   the accept/reject draw consume the SAME seeded generator, in
   sequence).
3. Discard the first 10% as burn-in, keep every 10th remaining sample,
   and build the monomer-fraction-by-cluster-size histogram.
4. Report the free-monomer fraction, the monomer-weighted mean cluster
   size, and the largest cluster size observed in the retained samples.
5. n_monomers = 28 <= 45, so the exact Boltzmann histogram over integer
   partitions is also available; report the max absolute deviation
   between the Monte Carlo histogram and this exact solution.

The specific sequence of accepted/rejected moves depends on the exact
pseudo-random draws from seed 55 -- there is no analytic shortcut for
steps 2-4. Report whether this run's first and second half agreed well
enough to call it converged, choosing from
["converged", "not_converged"]."""
    )
    items = [
        ("free_monomer_fraction", "free_monomer_fraction: fraction of monomers in size-1 clusters (step 4)"),
        ("mean_cluster_size", "mean_cluster_size: monomer-weighted mean cluster size (step 4)"),
        ("largest_cluster", "largest_cluster: largest cluster size observed in the retained samples (step 4, integer)"),
        ("exact_max_deviation", "exact_max_deviation: max abs deviation from the exact partition solution (step 5)"),
        ("convergence", "convergence: one of [\"converged\", \"not_converged\"]"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-18"] = dict(
        id="X-18", family="G4", levers=["L12", "L8"],
        axes={"arith": 4, "reason": 1},
        mechanism="Same genuine canonical-MC mechanism as X-14 with different n_monomers/T/epsilon/seed -- tests generalization of the simulated-answer requirement across physical conditions.",
        refs=["new-pattern-X12-X15-X19"], tools_A=["mc_aggregation_distribution"],
        fields={"free_monomer_fraction": num(res.free_monomer_fraction, 0.02),
                "mean_cluster_size": num(res.mean_cluster_size_monomer_weighted, 0.02),
                "largest_cluster": integer(res.largest_cluster_observed),
                "exact_max_deviation": num(res.exact_max_abs_deviation, 0.02),
                "convergence": exact("converged")},
        body=body)


# ============================================================ X-19
def build_x19():
    from surfactantkit.electrostatics import EPS0, WATER_REL_PERMITTIVITY_25C
    mu, eta = 4.10, 0.89
    z_plain = elec.zeta_potential_henry(mu, eta, "smoluchowski")
    ka = 14.5
    z_corr = elec.zeta_potential_relaxation_corrected(mu, ka, eta, 298.15, 1, 0.15)
    shift = abs(z_corr - z_plain) / abs(z_plain) * 100.0
    larger = "plain" if abs(z_plain) >= abs(z_corr) else "corrected"
    mu_SI = mu * 1e-8
    body = (
        qhead("X-19", "Zeta potential: unit chain + relaxation correction") +
        """Electrophoretic mobility of a colloidal dispersion (unit as commonly
reported):

    mu = 4.10 (um*cm)/(V*s),   viscosity = 0.89 mPa*s

Pinned formulas (authoritative):

    mobility_SI = mu * 1e-8                      [m^2/(V s)]
    eta_SI = viscosity * 1e-3                    [Pa s]
    Henry:  zeta = 3 * eta_SI * mobility_SI /
                   (2 * eps_r * eps0 * f)
    with eps_r = """ + f"{WATER_REL_PERMITTIVITY_25C}" + """, eps0 = 8.8541878128e-12 F/m,
    and for the smoluchowski regime f = 1.5; report zeta in mV.

    Then a relaxation-corrected determination on the SAME mobility with
    kappa*a = 14.5, z = 1, m = 0.15, T = 298.15 K, using the iterative
    Henry-equation correction (report its zeta in mV).

1. mobility_SI in m^2/(V s).
2. zeta_plain in mV (smoluchowski, formula above).
3. zeta_corr in mV (relaxation-corrected value).
4. shift = |zeta_corr - zeta_plain| / |zeta_plain| * 100 (percent).
5. Which determination is larger in magnitude? ["plain", "corrected"]."""
    )
    items = [
        ("mu_SI", "mu_SI: m^2/(V s) (step 1)"),
        ("zeta_plain_mV", "zeta_plain_mV: mV (step 2)"),
        ("zeta_corr_mV", "zeta_corr_mV: mV (step 3)"),
        ("shift_pct", "shift_pct: percent (step 4)"),
        ("larger", "larger: one of [\"plain\", \"corrected\"]"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-19"] = dict(
        id="X-19", family="G4", levers=["L8", "L1", "L12"],
        axes={"arith": 4, "reason": 1},
        mechanism="Unit chain (um cm -> m^2) + Henry factor bookkeeping + iterative correction value models cannot reach unaided: silent x1e-8 misses and f=1.5 dropped are the expected slips.",
        refs=["E8", "E1"], tools_A=["zeta_potential", "zeta_potential_relaxation_corrected"],
        fields={"mu_SI": num(mu_SI, 0.005, "m^2/(V s)"),
                "zeta_plain_mV": num(z_plain, 0.01, "mV"),
                "zeta_corr_mV": num(z_corr, 0.01, "mV"),
                "shift_pct": num(shift, 0.01, "%"),
                "larger": exact(larger)},
        body=body)


# ============================================================ X-20
def build_x20():
    res = kin.kmc_micelle_kinetics(
        c_total_M=0.030, K_per_M=1800.0, k_on_per_M_s=5.0e6,
        n_max=25, n_molecules=350, n_samples=50, seed=78)
    body = (
        qhead("X-20", "A second Gillespie kinetic Monte Carlo of stepwise micellization") +
        """A different micellization system, different K/concentration/molecule
count/seed from any other kinetic Monte Carlo question in this set:

    c_total = 0.030 M
    equilibrium constant K = 1800 /M (same K governs every n)
    on-rate k_on = 5.0e6 /M/s;  off-rate k_off = k_on / K
    max aggregate size tracked n_max = 25
    molecule count in the simulated volume = 350
    (initial state: a Poisson sample consistent with equilibrium under K)

Pinned procedure (authoritative; the SSA mechanics are unchanged from
every other kinetic Monte Carlo question here -- only the physical
inputs are new):

1. With a pseudo-random generator seeded EXACTLY with seed = 78, i.e. NumPy's default_rng(78) Generator (PCG64); every draw this procedure consumes comes from that one generator IN ORDER, and any percentile it asks for is NumPy's np.percentile with default linear interpolation -- Python's random.Random, legacy np.random.seed (MT19937) or torch streams do NOT reproduce this procedure,
   construct the Poisson-equilibrium initial state under K at
   c_total = 0.030 M for 350 molecules in a volume sized so the
   expected molecule count equals 350.
2. Run the direct SSA: compute every C_n+C_1/C_{n+1} propensity, draw
   an exponential waiting time and a reaction channel from those
   propensities (both from the SAME seeded generator, in sequence),
   fire the reaction, update state, repeat.
3. Sample the running state at 50 evenly-spaced checkpoints. Take the
   final 20% as the stationary window; average the monomer
   concentration (equilibrium_c1_M) and the mean aggregation number
   over n >= 2 (equilibrium_mean_nagg) over that window.
4. Fit a single exponential decay to the monomer trajectory for
   tau_fit_s.
5. Report the number of reaction events fired and the max deviation
   from exact mass conservation (sum n*count_n invariant by
   construction -- should be exactly 0 for a correct run).

As in any Gillespie SSA, the exact event sequence depends on the
specific pseudo-random draws from seed 78 -- only running the real
simulation reproduces equilibrium_c1_M/tau_fit_s to the required
precision. Report whether the final window was judged stationary,
choosing from ["converged", "not_converged"]."""
    )
    items = [
        ("equilibrium_c1_M", "equilibrium_c1_M: stationary-window mean monomer concentration, mol/L"),
        ("equilibrium_mean_nagg", "equilibrium_mean_nagg: stationary-window mean aggregation number (n>=2)"),
        ("tau_fit_s", "tau_fit_s: exponential relaxation time fit to the monomer trajectory, s"),
        ("mass_conservation_max_error", "mass_conservation_max_error: max deviation from exact mass conservation (integer, step 5)"),
        ("convergence", "convergence: one of [\"converged\", \"not_converged\"]"),
    ]
    body += "\n" + outputs_block(items)
    QUESTIONS["X-20"] = dict(
        id="X-20", family="G4", levers=["L12", "L8"],
        axes={"arith": 4, "reason": 1},
        mechanism="Same genuine Gillespie-SSA mechanism as X-11 with different K/concentration/seed -- tests generalization rather than a model having memorized one specific trajectory's statistics.",
        refs=["new-pattern-X12-X15-X19"], tools_A=["kmc_micelle_kinetics"],
        fields={"equilibrium_c1_M": num(res.equilibrium_c1_M, 0.02, "mol/L"),
                "equilibrium_mean_nagg": num(res.equilibrium_mean_nagg, 0.02),
                "tau_fit_s": num(res.tau_fit_s, 0.02, "s"),
                "mass_conservation_max_error": integer(res.mass_conservation_max_error),
                "convergence": exact("converged")},
        body=body)


def build():
    for f in (build_x11, build_x12, build_x13, build_x14, build_x15,
              build_x16, build_x17, build_x18, build_x19, build_x20):
        f()
    return QUESTIONS
