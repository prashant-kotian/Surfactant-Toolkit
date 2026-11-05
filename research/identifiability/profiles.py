"""Identifiability analysis: profile likelihood, and what it says that a
standard error cannot.

The point of this module is narrow and deliberate. Fitting a model and
reporting a standard error answers "how precisely is this parameter pinned,
assuming the model is right". It does not answer:

  (a) is the parameter bounded at all, or only bounded above / below?
  (b) is the model distinguishable from a rival mechanism?
  (c) does the answer survive the nuisance parameters being wrong?

Profile likelihood answers all three. Standard errors, covariance matrices and
R^2 answer none of them, which is why the interfacial rheology literature --
which almost universally reports the first kind of quantity -- can assign
mechanisms that its own data do not support.

Method
------
For a parameter of interest, hold it at a sequence of values and maximise the
likelihood over everything else. The resulting curve of best-achievable
chi-square is the profile. If it rises above the threshold on BOTH sides the
parameter is bounded; one side means it is only bounded from one direction;
neither means it is unidentifiable.

Thresholds use the chi-square quantile at 1 degree of freedom (3.841 for 95%)
because one parameter is being constrained at a time. This is the standard
profile-likelihood construction; it is not a heuristic.

References for the method itself, not for its application here:
  Raue et al., Bioinformatics 2009, 25, 1923 -- profile likelihood for
  practical identifiability, including the one-sided case.
  Cox & Snell / standard asymptotic likelihood theory for the threshold.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares, minimize_scalar

CHI2_95_1DOF = 3.841458820694124


def residuals(model, omega, data, theta, kind="phase"):
    return model.observables(omega, theta, kind) - data


def fit(model, omega, data, p0=None, kind="phase", n_starts=12, seed=0):
    """Multi-start bounded least squares.

    Multi-start is not decoration: the profile scans need the *global*
    optimum at each fixed value, and single-start gradient descent will
    silently return a local one, which flattens the profile and manufactures
    false identifiability. That is a failure mode I have already been bitten
    by once in this project.
    """
    rng = np.random.default_rng(seed)
    lo = np.array([model.lb[k] for k in model.params], float)
    hi = np.array([model.ub[k] for k in model.params], float)
    # work in log space for time constants, linear for E0
    lo_t = np.array([np.log10(lo[i]) if k in model.log_params else lo[i]
                     for i, k in enumerate(model.params)])
    hi_t = np.array([np.log10(hi[i]) if k in model.log_params else hi[i]
                     for i, k in enumerate(model.params)])

    def to_theta(t):
        return np.array([10.0 ** t[i] if k in model.log_params else t[i]
                         for i, k in enumerate(model.params)])

    def res_t(t):
        return residuals(model, omega, data, to_theta(t), kind)

    starts = []
    if p0 is not None:
        starts.append(np.array([np.log10(p0[i]) if k in model.log_params
                                else p0[i] for i, k in enumerate(model.params)]))
    for _ in range(n_starts):
        starts.append(lo_t + rng.random(len(lo_t)) * (hi_t - lo_t))

    best = None
    for s in starts:
        try:
            r = least_squares(res_t, np.clip(s, lo_t, hi_t),
                              bounds=(lo_t, hi_t), method="trf",
                              xtol=1e-12, ftol=1e-12, gtol=1e-12,
                              max_nfev=4000)
        except Exception:
            continue
        sse = float(np.sum(r.fun ** 2))
        if best is None or sse < best[0]:
            best = (sse, to_theta(r.x))
    return best


def profile(model, omega, data, name, grid, kind="phase", **fitkw):
    """Profile chi-square over one parameter.

    Returns (grid, chi2, best_other_params). chi2 is measured relative to the
    global minimum found on the grid, so profile[0] == 0 by construction.
    """
    idx = model.params.index(name)
    out_chi2, out_theta = [], []
    for g in grid:
        def res_rest(rest_t, g=g):
            theta = np.empty(len(model.params))
            for j, k in enumerate(model.params):
                if j == idx:
                    theta[j] = g
                else:
                    theta[j] = (10.0 ** rest_t[j] if k in model.log_params
                                else rest_t[j])
            return residuals(model, omega, data, theta, kind)

        others = [k for j, k in enumerate(model.params) if j != idx]
        lo = np.array([np.log10(model.lb[k]) if k in model.log_params
                       else model.lb[k] for k in others], float)
        hi = np.array([np.log10(model.ub[k]) if k in model.log_params
                       else model.ub[k] for k in others], float)
        if len(others) == 0:
            out_chi2.append(float(np.sum(res_rest(np.array([])) ** 2)))
            out_theta.append(np.array([g]))
            continue

        rng = np.random.default_rng(1)
        best = None
        for _ in range(8):
            s = lo + rng.random(len(lo)) * (hi - lo)
            try:
                r = least_squares(res_rest, s, bounds=(lo, hi),
                                  method="trf", xtol=1e-12, ftol=1e-12,
                                  gtol=1e-12, max_nfev=3000)
            except Exception:
                continue
            sse = float(np.sum(r.fun ** 2))
            if best is None or sse < best[0]:
                best = (sse, r.x)
        if best is None:
            out_chi2.append(np.inf)
            out_theta.append(np.full(len(others), np.nan))
        else:
            out_chi2.append(best[0])
            out_theta.append(best[1])
    chi2 = np.array(out_chi2)
    return grid, chi2, out_theta


def classify(grid, chi2, thresh=CHI2_95_1DOF):
    """Describe the profile shape.

    'bounded'     both tails exceed the threshold: an interval exists
    'one_sided'   only one tail exceeds: the parameter has a limit, not a value
    'flat'        neither: the data do not constrain it at all
    Also reports whether the best value sits on a bound, which is the
    signature of a parameter being pushed to a physical limit rather than
    measured.
    """
    d = chi2 - chi2.min()
    inside = grid[d <= thresh]
    lo_edge = d[0] <= thresh
    hi_edge = d[-1] <= thresh
    if len(inside) == 0:
        return {"shape": "resolved_but_narrow", "interval": None,
                "lower_open": False, "upper_open": False}
    lower_open = bool(lo_edge)
    upper_open = bool(hi_edge)
    if lower_open and upper_open:
        shape = "flat"
    elif lower_open or upper_open:
        shape = "one_sided"
    else:
        shape = "bounded"
    return {"shape": shape, "interval": (float(inside.min()), float(inside.max())),
            "lower_open": lower_open, "upper_open": upper_open}


def aic(sse, n, k):
    """AIC for Gaussian errors with unknown variance (constants dropped)."""
    return n * np.log(max(sse, 1e-300) / n) + 2 * k


def aicc(sse, n, k):
    """Small-sample correction, needed because n is 6 here, not 600."""
    a = aic(sse, n, k)
    if n - k - 1 <= 0:
        return np.inf
    return a + 2 * k * (k + 1) / (n - k - 1)


def akaike_weights(scores):
    s = np.asarray(scores, float)
    s = s - s.min()
    w = np.exp(-0.5 * s)
    return w / w.sum()
