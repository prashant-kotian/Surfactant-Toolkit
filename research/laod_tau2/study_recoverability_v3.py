"""Final recoverability map: WHERE in (tau2, capacity) space does LAOD work?

Pass 1 failed because it used capacity=0.2 (reservoir 25x too weak) and a
single tau2 range. Pass 2 showed capacity~5 gives sub-20% tau2 recovery but
degraded at tau2 >~ 0.3 s. This pass maps the boundary honestly, so the paper
can state its domain of validity instead of implying universality.
"""
import json, math
import numpy as np
from research.laod_tau2 import harmonics as H
from research.laod_tau2.model import ModelParams, run_oscillation

BASE = dict(c0=1.0, K=5.0, gamma_inf=3.5e-6, D=4.0e-10, gamma0=72.0)
GRID = [(w, 0.20) for w in (0.3, 1.0, 3.0, 6.28)]   # <= 1 Hz drop tensiometer
SCALE = np.array([0.005, 0.005, 1.0])

def feat(tau2, cap, om, eps=0.20, sigma=0.0, seed=0):
    p = ModelParams(**BASE, tau2=tau2, capacity=cap)
    r = run_oscillation(p, om, eps, n_cycles=14)
    t = np.asarray(r.t)
    g = H.add_noise(t, r.gamma_mN_m, sigma, seed) if sigma else np.asarray(r.gamma_mN_m)
    hg = H.harmonics(t, g, om)
    return dict(I2=hg[2]/hg[1], I3=hg[3]/hg[1], lag=H.lag_deg(t, g, r.area_rel, om))

def vector(tau2, cap, sigma=0.0, seed=0):
    return [feat(tau2, cap, om, sigma=sigma, seed=seed+i) for i,(om,_) in enumerate(GRID)]

def fit(obs, cap_known=None):
    from scipy.optimize import least_squares
    def res(x):
        tau2 = math.exp(x[0])
        cap = cap_known if cap_known is not None else math.exp(x[1])
        q = vector(tau2, cap); r=[]
        for a,b in zip(obs,q):
            r += [(b['I2']-a['I2'])/SCALE[0],(b['I3']-a['I3'])/SCALE[1],(b['lag']-a['lag'])/SCALE[2]]
        return np.asarray(r)
    x0 = [math.log(0.1)] if cap_known is not None else [math.log(0.1), math.log(5.0)]
    lo = [math.log(1e-4)] if cap_known is not None else [math.log(1e-4), math.log(1e-2)]
    hi = [math.log(50.)] if cap_known is not None else [math.log(10.), math.log(50.)]
    sol = least_squares(res, x0=x0, bounds=(lo,hi), xtol=1e-8, ftol=1e-8,
                        diff_step=1e-2, max_nfev=40)
    dof = max(len(sol.fun)-len(x0),1); s2 = 2*sol.cost/dof
    try:
        cov = np.linalg.inv(sol.jac.T@sol.jac)*max(s2,1.0); se=np.sqrt(np.diag(cov))
        rel = float(se[0]/abs(sol.x[0]))
    except np.linalg.LinAlgError:
        rel = float('nan')
    return dict(tau2=float(math.exp(sol.x[0])), rel_se=rel)

def main():
    out={}
    print("TAU2 RECOVERY MAP  (noise 0.05 mN/m, window 0.3-6.28 rad/s, eps=0.20)")
    print("two protocols: capacity fitted jointly, and capacity known from")
    print("equilibrium data (CMC, aggregation number -> capacity)\n")
    print(f"{'cap':>5} {'tau2_t':>8} | {'2-par fit':>10} {'relSE':>6} | {'1-par fit':>10} {'relSE':>6} {'verdict':>10}")
    for cap in (5.0, 1.4):
        for tau2 in (1e-3, 1e-2, 1e-1, 3e-1, 1.0):
            obs = vector(tau2, cap, sigma=0.05, seed=11)
            f2 = fit(obs); f1 = fit(obs, cap_known=cap)
            ok = f1['rel_se'] < 0.30
            out[f"c{cap}|t{tau2}"] = dict(two_param=f2, one_param=f1)
            print(f"{cap:5.2f} {tau2:8.0e} | {f2['tau2']:10.2e} {f2['rel_se']:6.2f} | "
                  f"{f1['tau2']:10.2e} {f1['rel_se']:6.2f} {'USABLE' if ok else 'not usable':>10}")
    json.dump(out, open('research/laod_tau2/recoverability_map.json','w'), indent=2)
    n_ok = sum(1 for v in out.values() if v['one_param']['rel_se'] < 0.30)
    print(f"\nusable cells (1-parameter, rel SE < 30%): {n_ok}/{len(out)}")
    print("wrote research/laod_tau2/recoverability_map.json")

if __name__ == "__main__": main()
