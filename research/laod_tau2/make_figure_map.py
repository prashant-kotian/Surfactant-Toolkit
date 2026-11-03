"""Figure 1: the (tau2, instrument) recoverability map for real systems.

The claim it must make unambiguously:
  an instrument resolves tau2 only when omega_max*tau2 >= ~0.3 (usable) or
  >= ~3 (precise), so each instrument covers a tau2 range - and no single
  instrument covers the real systems.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from research.laod_tau2.systems import SYSTEMS

OMEGA = {"Pendant drop / bubble (<=1 Hz)": 2*np.pi*1.0,
         "Capillary-pressure microtensiometer (100 Hz)": 2*np.pi*100.0}
USABLE, PRECISE = 0.3, 3.0

fig, ax = plt.subplots(figsize=(9.4, 5.2))
x = np.logspace(-4.5, 1.2, 900)

rows = []
for name, om in OMEGA.items():
    rows.append((name, "usable  (relSE < 30%)", om, USABLE, "#2b6cb0", 0.30))
    rows.append((name, "precise (relSE < 10%)", om, PRECISE, "#2b6cb0", 0.60))

ypos = []
y = 5.9
for name, quality, om, thr, colour, alpha in rows:
    left = thr / om
    ax.fill_between([left, 30], y - 0.28, y + 0.28, color=colour, alpha=alpha, lw=0)
    ax.text(6.3e-5, y, f"{name}\n{quality}", va="center", ha="left", fontsize=8,
            color="#1a202c")
    label = f"{left*1e3:.1f} ms" if left < 0.01 else f"{left*1e3:.0f} ms"
    ax.text(left*1.15, y + 0.0, label, va="center", fontsize=8, color="white",
            fontweight="bold")
    ypos.append(y)
    y -= 0.78

sysinfo = [("SDS, no salt", 0.002, 0.040, "#c53030"),
           ("SDS, 0.1 M NaCl", 1.0, 5.0, "#9b2c2c"),
           ("C12E8", 0.001, 0.010, "#2c5282"),
           ("CTAB, pure", 0.050, 0.500, "#6b46c1"),
           ("CTAB, with NaBr", 0.500, 5.0, "#44337a")]
y -= 0.30
for label, lo, hi, col in sysinfo:
    ax.plot([lo, hi], [y, y], color=col, lw=7, solid_capstyle="butt", alpha=0.9)
    ax.text(hi*1.18, y, label, va="center", fontsize=9, color=col)
    y -= 0.62

ax.set_xscale("log")
ax.set_xlim(1e-4, 2e1)
ax.set_ylim(y + 0.35, 6.6)
ax.set_yticks(ypos)
ax.set_yticklabels(["", "", "", ""])
ax.set_xlabel(r"micellar slow relaxation time $\tau_2$  (s)")
ax.set_title("Recoverable $\\tau_2$ is set by instrument bandwidth, not by the surfactant\n"
             "filled band = instrument can resolve it;  coloured bar = where the real system sits",
             fontsize=10.5)
ax.spines[["top", "right", "left"]].set_visible(False)
ax.grid(axis="x", alpha=0.15)
ax.axvline(1e-3, color="0.75", ls=":", lw=0.8)
ax.axvline(1e0, color="0.75", ls=":", lw=0.8)
for xv, lab in ((1e-3, "1 ms"), (1e0, "1 s")):
    ax.text(xv, y + 0.28, lab, ha="center", va="top", fontsize=7.5, color="0.5")
fig.tight_layout()
fig.savefig("research/laod_tau2/figures/fig1_recoverability_map.png", dpi=200)
print("wrote fig1")
