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

def model_lag(tau2, cap, om, eps=EPS):
    p = ModelParams(c0=C_TOT, K=K, gamma_inf=S.gamma_inf, D=S.D, gamma0=72.0,
                    n_factor=1.0, capacity=cap, tau2=tau2)
    r = run_oscillation(p, om, eps, n_cycles=16)
    t = np.asarray(r.t)
    return H.lag_deg(t, r.gamma_mN_m, r.area_rel, om)

def noisy_lag(tau2, cap, om, seed):
    p = ModelParams(c0=C_TOT, K=K, gamma_inf=S.gamma_inf, D=S.D, gamma0=72.0,
                    n_factor=1.0, capacity=cap, tau2=tau2)
    r = run_oscillation(p, om, EPS, n_cycles=16)
    t = np.asarray(r.t)
    g = H.add_noise(t, r.gamma_mN_m, SIGMA_G, seed=seed)
    return H.lag_deg(t, g, r.area_rel, om)

print("="*78)
print("RECOVERY UNCERTAINTY, C12E8 at 1.5xCMC (capacity 47.5), tau2 = 4.0 s")
print(f"tension noise {SIGMA_G} mN/m, eps = {EPS}, window {OMEGAS} rad/s")
print("="*78)

# ---- 1. Monte-Carlo scatter of the lag -----------------------------------
print("\n[1] Monte-Carlo: scatter of the recovered phase lag (60 realisations)")
mc = {}
for om in OMEGAS:
    vals = [noisy_lag(TAU2, CAP, om, seed=100+ i) for i in range(60)]
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
    def res(x):
        tau2 = math.exp(x[0])
        cap = CAP if one_param else math.exp(x[1])
        return [model_lag(tau2, cap, om) - o for om, o in zip(OMEGAS, obs)]
    x0 = [math.log(TAU2)] if one_param else [math.log(TAU2), math.log(CAP)]
    lo = [math.log(0.05)] if one_param else [math.log(0.05), math.log(1.0)]
    hi = [math.log(200.)] if one_param else [math.log(200.), math.log(2000.)]
    sol = least_squares(res, x0=x0, bounds=(lo, hi), xtol=1e-10, ftol=1e-10, diff_step=1e-3)
    dof = max(len(OMEGAS) - len(x0), 1)
    # scale by the MC-measured lag noise
    s2 = max(2*sol.cost/dof, np.mean([mc[o]['sd'] for o in OMEGAS])**2)
    cov = np.linalg.inv(sol.jac.T @ sol.jac) * s2
    se = np.sqrt(np.diag(cov))
    out = dict(tau2=float(math.exp(sol.x[0])), relSE=float(se[0]))
    if not one_param:
        out["capacity"] = float(math.exp(sol.x[1]))
        out["capacity_relSE"] = float(se[1]/math.exp(sol.x[1]))
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
