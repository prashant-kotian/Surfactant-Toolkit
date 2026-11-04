"""RECOVERY UNCERTAINTY under the optimized design. The decisive computation.

Design (per Antigravity's brief): C12E8, c = 1.5x CMC, capacity 47.5,
tau2 = 4.0 s, eps = 0.20, phase lag as the observable.

A single frequency cannot give an uncertainty (zero degrees of freedom), so
this runs a SWEEP across the usable window and reports three independent
estimates:

  1. Monte-Carlo on the phase lag: inject 0.05 mN/m tension noise into
     gamma(t), re-demodulate, measure the scatter of the recovered lag.
  2. Delta method: relSE = sigma_lag / |d lag / d ln tau2|, which needs no fit.
  3. Joint least-squares with covariance: 1-parameter (capacity from the
     equilibrium closure) and 2-parameter (capacity free), reporting the
     tau2/capacity correlation -- the degeneracy that decides whether the
     2-parameter version is honest.
"""
import json, math
import numpy as np
from research.laod_tau2 import harmonics as H
from research.laod_tau2.model import ModelParams, run_oscillation
from research.laod_tau2.systems import SYSTEMS

S = SYSTEMS["C12E8"]; K = 140.0
C_TOT = 1.5 * S.cmc_mM
CAP = S.capacity(C_TOT)
TAU2 = 4.0
EPS = 0.20
SIGMA_G = 0.05                 # mN/m
OMEGAS = (3.0, 4.0, 5.0, 6.28)  # the usable part of the PAT window

# ---------------------------------------------------------------------------
# SAMPLING RATE (audit item 4.3, FINDINGS 30.15)
#
# run_oscillation defaults to 320 samples per cycle. At the frequencies used
# here (omega = 3-6.28 rad/s, i.e. 0.48-1.0 Hz) that is 153-320 samples/second,
# while the instrument the paper proposes to use acquires 25-50 frames/second.
# White noise was then attached independently to every solver sample, so the
# synthetic experiment contained 2.5-3.6x more independent measurements than the
# real one. Precision figures computed that way are optimistic by up to
# sqrt(3.58) in the mean, i.e. by a factor of ~1.6-1.9.
#
# We now decimate to a declared instrument frame rate before adding noise, so
# the Monte Carlo is conditioned on the acquisition the paper actually proposes.
# ---------------------------------------------------------------------------
INSTRUMENT_FPS = 25.0          # PAT-1M, conservative end of 25-50 fps


def _acquire(t, y, fps=INSTRUMENT_FPS):
    """Resample onto the instrument's uniform frame clock."""
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    n = max(8, int(round((t[-1] - t[0]) * fps)) + 1)
    tt = np.linspace(t[0], t[-1], n)
    return tt, np.interp(tt, t, y)


def model_lag(tau2, cap, om, eps=EPS):
    p = ModelParams(c0=C_TOT, K=K, gamma_inf=S.gamma_inf, D=S.D, gamma0=72.0,
                    n_factor=1.0, capacity=cap, tau2=tau2)
    r = run_oscillation(p, om, eps, n_cycles=16)
    t, g = _acquire(r.t, r.gamma_mN_m)
    a = np.interp(t, np.asarray(r.t), np.asarray(r.area_rel))
    return H.lag_deg(t, g, a, om)

def noisy_lag(tau2, cap, om, seed):
    p = ModelParams(c0=C_TOT, K=K, gamma_inf=S.gamma_inf, D=S.D, gamma0=72.0,
                    n_factor=1.0, capacity=cap, tau2=tau2)
    r = run_oscillation(p, om, EPS, n_cycles=16)
    t, g = _acquire(r.t, r.gamma_mN_m)
    a = np.interp(t, np.asarray(r.t), np.asarray(r.area_rel))
    g = H.add_noise(t, g, SIGMA_G, seed=seed)
    return H.lag_deg(t, g, a, om)

print("="*78)
print("RECOVERY UNCERTAINTY, C12E8 at 1.5xCMC (capacity 47.5), tau2 = 4.0 s")
print(f"acquired at {INSTRUMENT_FPS:g} fps (instrument-realistic), not 320/cycle")
print(f"tension noise {SIGMA_G} mN/m, eps = {EPS}, window {OMEGAS} rad/s")
print("="*78)

# ---- 1. Monte-Carlo scatter of the lag -----------------------------------
print("\n[1] Monte-Carlo: scatter of the recovered phase lag (60 realisations)")
mc = {}
for om in OMEGAS:
    # (c) independent stream per frequency: offset by omega-index, not reused
    vals = [noisy_lag(TAU2, CAP, om, seed=1000*int(round(om*100)) + i) for i in range(60)]
    mc[om] = dict(mean=float(np.mean(vals)), sd=float(np.std(vals, ddof=1)))
    print(f"  omega={om:5.2f}: lag = {mc[om]['mean']:7.2f} +/- {mc[om]['sd']:.3f} deg")

# ---- 2. Delta method on a single point -----------------------------------
print("\n[2] Delta method (no fit needed): relSE = sigma_lag / |d lag/d ln tau2|")
delta = {}
for om in OMEGAS:
    l1 = model_lag(TAU2*0.8, CAP, om); l2 = model_lag(TAU2*1.25, CAP, om)
    slope = abs(l2 - l1) / math.log(1.25/0.8)      # deg per e-fold of tau2
    relse = mc[om]["sd"] / slope if slope > 0 else float("inf")
    delta[om] = dict(slope_deg_per_efold=float(slope), relSE=float(relse))
    print(f"  omega={om:5.2f}: dlag/dln tau2 = {slope:6.2f} deg  ->  relSE = {relse:6.3f}")

# ---- 3. Joint least squares over the sweep -------------------------------
from scipy.optimize import least_squares
obs = [model_lag(TAU2, CAP, om) for om in OMEGAS]

def fit(one_param):
    x0 = [math.log(TAU2)] if one_param else [math.log(TAU2), math.log(CAP)]
    lo = [math.log(0.05)] if one_param else [math.log(0.05), math.log(1.0)]
    hi = [math.log(200.)] if one_param else [math.log(200.), math.log(2000.)]
    dof = max(len(OMEGAS) - len(x0), 1)
    # FINDINGS 30.16 (audit item 4.4). Three defects fixed here:
    #
    #  (a) The fit was UNWEIGHTED although the Monte-Carlo phase SD varies by a
    #      factor of 1.6 across the sweep (0.212 deg at omega=3 vs 0.133 at
    #      omega=6.28). A scalar noise estimate across unequal-SD frequencies
    #      misweights the fit. We now weight by the measured per-frequency SD.
    #  (b) Both parameters are fitted in LOG coordinates, so
    #      sqrt(cov[i,i]) is ALREADY an approximate RELATIVE uncertainty. The
    #      capacity uncertainty was then divided by capacity a second time.
    #  (c) Seeds were reused across frequencies, which makes the residuals
    #      artificially correlated between omega values. Each frequency now gets
    #      an independent stream, and we say so.
    sd = np.array([mc[o]["sd"] for o in OMEGAS])
    w = 1.0 / sd                                        # (a)
    def res_w(x):
        return [wi * (model_lag(math.exp(x[0]), CAP if one_param else math.exp(x[1]), om) - o)
                for wi, om, o in zip(w, OMEGAS, obs)]

    sol = least_squares(res_w, x0=x0, bounds=(lo, hi),
                        xtol=1e-10, ftol=1e-10, diff_step=1e-3)
    chi2 = float(np.sum(np.square(sol.fun)))
    scale = max(chi2 / dof, 1.0)      # 1.0 when the model is adequate and sigma known
    cov = np.linalg.inv(sol.jac.T @ sol.jac) * scale
    se = np.sqrt(np.diag(cov))        # (b) already relative, log coordinates
    out = dict(tau2=float(math.exp(sol.x[0])), relSE=float(se[0]),
               chi2=chi2, dof=dof)
    if not one_param:
        out["capacity"] = float(math.exp(sol.x[1]))
        out["capacity_relSE"] = float(se[1])        # (b) no second division
        out["corr"] = float(cov[0,1]/math.sqrt(cov[0,0]*cov[1,1]))
    return out

print("\n[3] Joint least squares over the sweep")
f1 = fit(True);  print(f"  1-parameter (capacity fixed): tau2 = {f1['tau2']:.3f} s, relSE = {f1['relSE']:.3f}")
f2 = fit(False); print(f"  2-parameter (capacity free): tau2 = {f2['tau2']:.3f} s, relSE = {f2['relSE']:.3f}, "
                      f"capacity = {f2['capacity']:.1f} (relSE {f2['capacity_relSE']:.2f}), "
                      f"corr = {f2['corr']:.2f}")

res = dict(design=dict(system="C12E8", c_tot_mM=C_TOT, capacity=CAP, tau2=TAU2,
                       eps=EPS, sigma_g=SIGMA_G, omegas=list(OMEGAS)),
           monte_carlo=mc, delta={str(k): v for k, v in delta.items()},
           joint_1p=f1, joint_2p=f2)
json.dump(res, open("research/laod_tau2/recovery_uncertainty.json", "w"), indent=2)
print("\nwrote research/laod_tau2/recovery_uncertainty.json")
