# Identifiability of interfacial dilatational relaxation mechanisms

Exploratory code for the candidate in `NEXT_PAPER_GAP_ANALYSIS.md` §2.

All four mechanism models are written through one object, the surface relaxation
admittance B(omega), which makes them comparable and makes the 45 degree
generalisation derivable rather than numerical.

| file | what it does |
|---|---|
| `models.py` | B(omega) for four mechanisms; the observable; the 45 + |arg B|/2 cap |
| `fit_data.py` | AICc comparison and profile likelihood on the C14E8 data; nuisance scan |
| `design_map.py` | which (tau, band, C_m) combinations are detectable at 3 sigma |
| `amplitude_test.py` | phase-only vs (|E|, phase); the h/sqrt(D) degeneracy |
| `tests/` | 11 tests, each pinning a step against an independent reference |

## Conventions that have bitten this project before

* `lag_deg` returns |delta|, delta = -arg E. Pinned by test 1 against
  `laod_tau2.linear_response`, which was itself verified against a PDE solver.
* E0 never enters the phase. Do not fit it to phase-only data; it will wander
  to any value and the fit will still look perfect.
* h and D are degenerate through h/sqrt(D), structurally, at any precision.
  Never report both from one experiment.

## Reading order

`FINDINGS_IDENTIFIABILITY.md` first, then the code.

## Caveat

The 67.5 degree class bound was WITHDRAWN in round 29 -- it was a scan
artefact. The ceiling for one relaxation is
45 + atan(C/(2 sqrt(C+1)))/2, which rises to 90 deg, and Fainerman 2008 already
states that observation in words. The surviving contribution is the closed-form
inversion (ceiling -> lower bound on C), not the bound itself.
