"""Recoverability, second pass: realistic micellar capacity.

Pass 1 used capacity = 0.2 (a 25x-too-weak reservoir) and failed. The
leverage diagnostic showed identifiability is controlled by the capacity,
not by coverage: capacity ~5 (micelles holding most of the surfactant, the
normal situation well above the CMC) raises the tau2 leverage from ~0.9x to
~13x the noise floor. This pass inverts in that regime.
"""
import json, math, sys
import numpy as np
from . import harmonics as H
from .model import ModelParams, run_oscillation

BASE = dict(c0=1.0, K=5.0, gamma_inf=3.5e-6, D=4.0e-10, gamma0=72.0)
GRID = [(w, 0.20) for w in (0.3, 1.0, 3.0, 2*math.pi*1.0)]
SCALE = np.array([0.005, 0.005, 1.0])   # floors on I2/I1, I3/I1, lag[deg]

def obs_at(tau2, cap, sigma, seed, grid=GRID):
    out=[]
    for i,(om,eps) in enumerate(grid):
        p=ModelParams(**BASE,tau2=tau2,capacity=cap)
        r=run_oscillation(p,om,eps,n_cycles=14)
        t=np.asarray(r.t); g=H.add_noise(t,r.gamma_mN_m,sigma,seed+i)
        hg=H.harmonics(t,g,om)
        out.append(dict(I2=hg[2]/hg[1],I3=hg[3]/hg[1],lag=H.lag_deg(t,g,r.area_rel,om),I1=hg[1]))
    return out

def pred(tau2, cap, grid=GRID):
    out=[]
    for (om,eps) in grid:
        p=ModelParams(**BASE,tau2=tau2,capacity=cap)
        r=run_oscillation(p,om,eps,n_cycles=14)
        t=np.asarray(r.t); hg=H.harmonics(t,r.gamma_mN_m,om)
        out.append(dict(I2=hg[2]/hg[1],I3=hg[3]/hg[1],lag=H.lag_deg(t,r.gamma_mN_m,r.area_rel,om),I1=hg[1]))
    return out

def fit(o):
    from scipy.optimize import least_squares
    def res(x):
        tau2,cap=math.exp(x[0]),math.exp(x[1])
        q=pred(tau2,cap); r=[]
        for a,b in zip(o,q):
            r += [(b['I2']-a['I2'])/SCALE[0],(b['I3']-a['I3'])/SCALE[1],(b['lag']-a['lag'])/SCALE[2]]
        return np.asarray(r)
    sol=least_squares(res,x0=[math.log(0.1),math.log(5.0)],
                      bounds=([math.log(1e-4),math.log(1e-2)],[math.log(20.),math.log(50.)]),
                      xtol=1e-8,ftol=1e-8,diff_step=1e-2,max_nfev=60)
    dof=max(len(sol.fun)-2,1); s2=2*sol.cost/dof
    try:
        cov=np.linalg.inv(sol.jac.T@sol.jac)*max(s2,1.0); se=np.sqrt(np.diag(cov))
        rel=se/np.abs(sol.x)
    except np.linalg.LinAlgError:
        rel=np.array([np.nan,np.nan])
    return dict(tau2=float(math.exp(sol.x[0])),cap=float(math.exp(sol.x[1])),
                rel_se=rel.tolist(),cost=float(sol.cost))

def main():
    res={}
    print(f"{'sigma':>6} {'tau2_t':>8} {'cap_t':>6} {'tau2_f':>9} {'relSE':>7} {'cap_f':>7} {'relSE':>7} {'I1[mN/m]':>9}")
    for sigma in (0.05,0.01):
        for tau2,cap in [(1e-2,5.0),(1e-1,5.0),(1.0,5.0),(1e-1,1.4)]:
            o=obs_at(tau2,cap,sigma,11); f=fit(o)
            res[f"s{sigma}|t{tau2}|c{cap}"]=f
            print(f"{sigma:6.2f} {tau2:8.0e} {cap:6.2f} {f['tau2']:9.3e} {f['rel_se'][0]:7.2f} "
                  f"{f['cap']:7.2f} {f['rel_se'][1]:7.2f} {np.mean([d['I1'] for d in o]):9.3f}")
    json.dump(res,open('research/laod_tau2/recoverability_v2.json','w'),indent=2)

if __name__=="__main__": main()
