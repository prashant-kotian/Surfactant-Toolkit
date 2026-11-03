"""Is the harmonic signal even MEASURABLE at real capacity?

The feasibility study used capacity 5 and found I1 ~ 0.05-0.23 mN/m --
already at the edge of the 0.05 mN/m noise floor. Real systems have
capacity 62-855, and a large capacity means the micelles buffer the
interface: Gamma barely moves, so gamma barely moves, so the harmonic
AMPLITUDES shrink even if the RATIOS stay informative.

A ratio is useless if its numerator is below the noise. This script
therefore reports ABSOLUTE amplitudes, not ratios, against the noise floor.
"""
import math, json
import numpy as np
from research.laod_tau2 import harmonics as H
from research.laod_tau2.model import ModelParams, run_oscillation
from research.laod_tau2.systems import SYSTEMS

NOISE = 0.05          # mN/m, pendant-drop resolution
OMEGAS = (0.3, 1.0, 3.0, 6.28)

def probe(K, sysname, c_tot, tau2_ms, om, eps=0.20):
    s = SYSTEMS[sysname]
    p = ModelParams(c0=c_tot, K=K, gamma_inf=s.gamma_inf, D=s.D,
                    gamma0=72.0, n_factor=float(s.n_ions),
                    capacity=s.capacity(c_tot), tau2=tau2_ms*1e-3)
    r = run_oscillation(p, om, eps, n_cycles=14)
    t = np.asarray(r.t); hg = H.harmonics(t, r.gamma_mN_m, om)
    return dict(I1=hg[1], I2=hg[2], I3=hg[3],
                I2_over_I1=hg[2]/hg[1] if hg[1] else 0.0,
                I3_over_I1=hg[3]/hg[1] if hg[1] else 0.0,
                lag=H.lag_deg(t, r.gamma_mN_m, r.area_rel, om),
                capacity=p.capacity, theta=p.theta0)

print("MEASURABILITY at REAL capacity (eps=0.20, noise floor 0.05 mN/m)")
print("I2, I3 are ABSOLUTE amplitudes in mN/m. If I2 < 0.05 the second")
print("harmonic is not observable, whatever its ratio to I1.\n")
print(f"{'system':7} {'c_tot':>6} {'Cm':>7} {'tau2':>7} {'w':>5} {'I1':>8} {'I2':>9} {'I2/I1':>7} {'I2>noise':>9}")
rows=[]
# only C12E8 has a supplied Langmuir K, so it is the only runnable system
for c_mult in (2, 10):
    ct = c_mult*SYSTEMS['C12E8'].cmc_mM
    for tau2 in (1.0, 10.0):
        for om in OMEGAS:
            d = probe(140.0, 'C12E8', ct, tau2, om)
            rows.append(dict(system='C12E8', c_tot=ct, tau2=tau2, omega=om, **d))
            print(f"{'C12E8':7} {ct:6.3f} {d['capacity']:7.1f} {tau2:7.1f} {om:5.2f} "
                  f"{d['I1']:8.5f} {d['I2']:9.6f} {d['I2_over_I1']:7.4f} "
                  f"{'YES' if d['I2']>NOISE else 'no':>9}")
json.dump(rows, open('research/laod_tau2/real_amplitude.json','w'), indent=2)
print("\nwrote research/laod_tau2/real_amplitude.json")
