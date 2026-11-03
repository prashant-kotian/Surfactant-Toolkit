"""Does the usable tau2 window follow the instrument bandwidth?

Claim under test: LAOD recovers tau2 only when omega*tau2 ~ 1 falls INSIDE
the instrument's frequency window, so the required bandwidth is set by tau2:
    omega_max >~ 1/tau2   (and omega_min open enough to see the low-frequency branch)
If that rule holds, it is the paper's practical deliverable: an instrument
specification per micellar system, rather than a blanket "we can measure it".
"""
import json, math
import numpy as np
from research.laod_tau2 import harmonics as H
from research.laod_tau2.model import ModelParams, run_oscillation

BASE = dict(c0=1.0, K=5.0, gamma_inf=3.5e-6, D=4.0e-10, gamma0=72.0)
SCALE = np.array([0.005, 0.005, 1.0])

def feat(tau2, cap, om, eps=0.20, sigma=0.0, seed=0):
    p = ModelParams(**BASE, tau2=tau2, capacity=cap)
    r = run_oscillation(p, om, eps, n_cycles=10)
    t = np.asarray(r.t)
    g = H.add_noise(t, r.gamma_mN_m, sigma, seed) if sigma else np.asarray(r.gamma_mN_m)
    hg = H.harmonics(t, g, om)
    return dict(I2=hg[2]/hg[1], I3=hg[3]/hg[1], lag=H.lag_deg(t, g, r.area_rel, om))

def fit_1p(tau2_true, cap, grid, sigma=0.05, seed=11):
    from scipy.optimize import least_squares
    obs = [feat(tau2_true, cap, om, sigma=sigma, seed=seed+i) for i,(om,_) in enumerate(grid)]
    def res(x):
        tau2 = math.exp(x[0]); r=[]
        for (om,_),a in zip(grid, obs):
            b = feat(tau2, cap, om)
            r += [(b['I2']-a['I2'])/SCALE[0],(b['I3']-a['I3'])/SCALE[1],(b['lag']-a['lag'])/SCALE[2]]
        return np.asarray(r)
    sol = least_squares(res, x0=[math.log(max(tau2_true,1e-4))], bounds=([math.log(1e-5)],[math.log(50.)]),
                        xtol=1e-7, ftol=1e-7, diff_step=1e-2, max_nfev=30)
    dof = max(len(sol.fun)-1,1); s2 = 2*sol.cost/dof
    try:
        cov = np.linalg.inv(sol.jac.T@sol.jac)*max(s2,1.0); se = math.sqrt(cov[0,0])
        rel = se/abs(sol.x[0])
    except np.linalg.LinAlgError:
        rel = float('nan')
    return math.exp(sol.x[0])/tau2_true, rel

WINDOWS = {
    "1 Hz drop (2-1000 mHz)":      [0.01, 0.03, 0.1, 0.3, 1.0, 6.28],
    "10 Hz CPT":                    [0.1, 1.0, 6.28, 20.0, 62.8],
    "100 Hz microtensiometer":      [1.0, 10.0, 62.8, 200.0, 628.0],
}
TAUS = [1e-3, 1e-2, 1e-1]

out = {}
print("BANDWIDTH RULE TEST (capacity=5, sigma=0.05 mN/m, eps=0.20)")
print(f"{'window':28} {'tau2':>8} {'fit/true':>9} {'relSE':>7} {'omega*tau2 cover':>18}")
for name, ws in WINDOWS.items():
    grid = [(w, 0.20) for w in ws]
    for tau2 in TAUS:
        try:
            ratio, rel = fit_1p(tau2, 5.0, grid)
        except Exception as e:
            ratio, rel = float('nan'), float('nan')
        cover = f"{min(ws)*tau2:.2g}-{max(ws)*tau2:.2g}"
        out[f"{name}|{tau2}"] = dict(ratio=ratio, rel_se=rel, window=ws)
        print(f"{name:28} {tau2:8.0e} {ratio:9.2f} {rel:7.2f} {cover:>18}")
json.dump(out, open('research/laod_tau2/bandwidth_rule.json','w'), indent=2)
print("\nwrote research/laod_tau2/bandwidth_rule.json")
