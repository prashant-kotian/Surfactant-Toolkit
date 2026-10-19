"""Raw experimental curve analysis -- the missing first step every other
module in this library assumes has already happened. Every function in
adsorption.py, mixed_micelle.py, etc. takes an already-known CMC or slope as
input; nothing in this library previously took a raw surface-tension-vs-
concentration table and extracted CMC from it. Built 2026-09-07 specifically
to close that gap, after a real methodology correction: the benchmark's
original plan tested models on pre-parameterized textbook questions, but the
actual research workflow starts from raw measured data, and that's the gap
this module (and the benchmark built on top of it) targets.

Method: piecewise-linear break-point detection on gamma vs log10(C) -- the
same real method described in Shah/Das/Bhattarai 2025 (Heliyon,
PMC11835642): "A graph is plotted with gamma and log C... gamma shows a
breakpoint -- this point is called CMC." Standard, not invented: fit every
possible split of the sorted (log C, gamma) points into a declining
premicellar line and a flat-ish post-CMC segment, and take the split that
minimizes total residual sum of squares across both segments.

UPGRADED 2026-09-12 after a real researcher-submitted dataset exposed a
genuine failure mode of the plain two-segment version: real tensiometry
data often has a flat "lag" baseline at low concentration (surface tension
near pure solvent, below where surfactant adsorption becomes significant)
BEFORE the real declining Gibbs-adsorption-linear region starts. The plain
two-segment model can't represent three regimes (flat baseline, decline,
flat plateau) and, verified numerically on the real submitted data, picks
the WRONG breakpoint (baseline-to-decline instead of decline-to-plateau)
because forcing one line through [decline+plateau] combined can score a
deceptively lower combined RSS than forcing one line through [baseline+
decline] combined -- even though the latter split is the physically
correct one. Fixed by also trying a three-segment (baseline, decline,
plateau) fit and choosing between the two- and three-segment models via
BIC (Bayesian Information Criterion, the standard statistic for this exact
model-selection problem in segmented/piecewise regression -- penalizes the
three-segment model's extra parameters so it is only selected when a real
third regime is present, not merely because more segments can always fit
at least as well). Verified this does not regress the existing literature
validation: the AOT pilot dataset (Shah/Das/Bhattarai 2025, no baseline
lag) still correctly selects the two-segment model and reports the same
CMC as before.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from ._stats import bic as _bic
from .cpp import tanford_tail_volume


def _linreg(xs: list[float], ys: list[float]) -> tuple[float, float, float]:
    """Ordinary least squares. Returns (slope, intercept, residual_sum_sq)."""
    n = len(xs)
    mean_x = sum(xs) / n
    mean_y = sum(ys) / n
    sxx = sum((x - mean_x) ** 2 for x in xs)
    sxy = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    if sxx == 0:
        slope = 0.0
    else:
        slope = sxy / sxx
    intercept = mean_y - slope * mean_x
    rss = sum((y - (slope * x + intercept)) ** 2 for x, y in zip(xs, ys))
    return slope, intercept, rss


@dataclass
class CmcFromCurveResult:
    cmc_mM: float
    premicellar_slope_mN_per_m_per_lnC: float  # dGamma/d(ln C), the form
                                                 # gibbs_gamma_max() expects --
                                                 # always the DECLINING segment's
                                                 # slope, whether or not a flat
                                                 # baseline was also detected
    premicellar_slope_mN_per_m_per_log10C: float  # as directly fit, for reference
    gamma_at_cmc_mN_per_m: float
    r_squared_premicellar: float  # of the declining segment specifically
    n_premicellar_points: int  # points in the declining segment specifically
    n_postmicellar_points: int
    premicellar_x_min_mM: float  # lowest concentration actually in the declining
                                   # segment -- NOT necessarily the dataset's lowest
                                   # concentration when a flat baseline was detected
    postmicellar_slope_mN_per_m_per_log10C: float = 0.0  # the plateau segment's own
                                                            # fitted slope -- near 0
                                                            # for a real flat plateau;
                                                            # a real diagnostic, not
                                                            # assumed to be exactly 0
    postmicellar_intercept_mN_per_m: float = 0.0
    postmicellar_mean_gamma_mN_per_m: float = 0.0  # the plateau's mean measured
                                                      # surface tension -- the simplest,
                                                      # most directly citable "postmicellar
                                                      # plateau value"
    r_squared_postmicellar: float = 0.0
    n_baseline_points: int = 0  # 0 when no flat pre-onset baseline was detected/used
    gamma0_mN_m: float = 0.0  # pure-solvent-like surface tension: the detected flat
                                # baseline's own mean if one was found, else the
                                # dataset's own lowest-concentration point (the
                                # standard practical proxy for a true C=0 blank) --
                                # real measured data, never a looked-up constant
    c20_mM: float | None = None  # concentration needed to reduce gamma0 by 20 mN/m
                                    # (the standard surfactant efficiency parameter) --
                                    # None if the premicellar decline never reaches a
                                    # 20 mN/m drop before the CMC (a real, disclosed
                                    # "not defined for this system" case, not an error)
    pC20: float | None = None  # -log10(c20_mM/1000) on the molar scale; None with c20_mM
    method: str = "two-segment linear regression break-point (min total RSS)"


def _best_two_segment_split(log_c: list[float], gamma: list[float], n: int, min_points_per_segment: int):
    """Returns (total_rss, split, pre_slope, pre_intercept, pre_rss)
    minimizing combined RSS across exactly 2 segments."""
    best = None
    for split in range(min_points_per_segment, n - min_points_per_segment + 1):
        pre_slope, pre_intercept, pre_rss = _linreg(log_c[:split], gamma[:split])
        _, _, post_rss = _linreg(log_c[split:], gamma[split:])
        total_rss = pre_rss + post_rss
        if best is None or total_rss < best[0]:
            best = (total_rss, split, pre_slope, pre_intercept, pre_rss)
    return best


def _best_three_segment_split(log_c: list[float], gamma: list[float], n: int, min_points_per_segment: int):
    """Returns (total_rss, split1, split2, decline_slope, decline_intercept,
    decline_rss) minimizing combined RSS across exactly 3 segments
    (baseline, decline, plateau). None if there aren't enough points for
    3 segments of min_points_per_segment each."""
    if n < 3 * min_points_per_segment:
        return None
    best = None
    for split1 in range(min_points_per_segment, n - 2 * min_points_per_segment + 1):
        for split2 in range(split1 + min_points_per_segment, n - min_points_per_segment + 1):
            _, _, base_rss = _linreg(log_c[:split1], gamma[:split1])
            decline_slope, decline_intercept, decline_rss = _linreg(log_c[split1:split2], gamma[split1:split2])
            _, _, plateau_rss = _linreg(log_c[split2:], gamma[split2:])
            total_rss = base_rss + decline_rss + plateau_rss
            if best is None or total_rss < best[0]:
                best = (total_rss, split1, split2, decline_slope, decline_intercept, decline_rss)
    return best


def cmc_from_surface_tension_curve(
    concentrations_mM: list[float],
    surface_tensions_mN_per_m: list[float],
    min_points_per_segment: int = 3,
    method: str = "rss",
) -> CmcFromCurveResult:
    """Extracts CMC from a raw gamma-vs-concentration tensiometry curve by
    finding the concentration that best splits the data into a declining
    premicellar line (gamma vs log10 C, real Gibbs-adsorption-linear region)
    and a flat post-micellar segment, minimizing combined residual sum of
    squares across every valid split point.

    Real, standard method -- not a guess: this is exactly what
    Shah/Das/Bhattarai 2025 describe doing by eye/software ("the breakpoint
    -- this point is called CMC"), just made a deterministic, repeatable
    algorithm instead of a by-eye read.

    Also tries a 3-segment model (flat baseline, decline, flat plateau) for
    data with a flat pre-onset lag region at low concentration, and picks
    between the 2- and 3-segment models via BIC -- see this module's own
    docstring for the real dataset that made this necessary. The returned
    premicellar_slope/r_squared/n_premicellar_points always describe the
    DECLINING segment specifically, whether or not a baseline was found.

    Also returns C20/pC20, the standard surfactant efficiency parameter
    (added 2026-09-12): C20 is the concentration required to reduce gamma0
    (the pure-solvent-like baseline surface tension, itself derived from
    THIS curve's own data -- the detected flat baseline's mean, or the
    lowest-concentration point if no baseline was found) by exactly
    20 mN/m; pC20 = -log10(C20 in mol/L). Confirmed real, standard
    definition against an open-access primary source (a cardanol-
    surfactant paper reporting pC20=4.03 alongside CMC/gamma_CMC from the
    same curve, C20 << CMC as physically expected). c20_mM/pC20 are None
    when the premicellar decline never reaches a 20 mN/m drop before the
    CMC -- a real "not defined for this system" case, not extrapolated
    past where the fit is valid.

    Real bug fixed 2026-09-13, found via the new autonomous orchestrator on
    a real 8-point SDS dataset: _best_three_segment_split requires
    n >= 3*min_points_per_segment, so with the default min_points_per_segment=3
    the 3-segment model was silently UNREACHABLE for any dataset under 9
    points -- which includes most real tensiometry curves (typically
    6-10 points) and every MCP-exposed call (mcp_server.py's tool wrapper
    does not expose this parameter at all, so every tool-augmented run was
    silently stuck on 2-segment-only for such datasets). This produced a
    spurious early breakpoint on the SDS case (CMC=4.59 mM, A_min=111 A^2
    -- physically implausible) purely because the correct 3-segment fit
    was never attempted.

    Simply lowering the default to 2 outright was tried and REJECTED: it
    regressed the AOT literature validation case (a 9-point dataset with
    no real baseline lag), because min_points_per_segment=2 lets the
    3-segment search exploit a trivial 2-point segment (a line through 2
    points always has RSS=0) to spuriously beat the 2-segment model's BIC
    even when no real third regime exists -- verified directly (AOT case
    went from correctly selecting 2-segment to spuriously selecting a
    5-point "baseline" that doesn't exist in that data).

    Fixed properly instead: min_points_per_segment stays 3 by default
    (preserving the AOT case's protection against the 2-point-segment
    exploit whenever 3-segment is already reachable at that value), but
    the 3-segment attempt specifically relaxes its OWN minimum -- down to
    a floor of 2, and only as far as n actually requires
    (max(2, n // 3)) -- when n is too small to reach 3-segment at the
    caller's min_points_per_segment at all. The 2-segment search always
    keeps the caller's original min_points_per_segment, untouched, since
    that search was never the unreachable one. Verified this recovers the
    correct SDS result (CMC=9.18 mM, matching the value obtained by
    manually passing min_points_per_segment=2, the pre-fix workaround)
    while leaving the AOT case's result exactly as it was before this fix.

    method selects the break-point estimator: "rss" (default -- the
    two-/three-segment BIC selection described above, unchanged),
    "tangential" (intersection of the steepest premicellar fit with the
    postmicellar fit -- Rosen's graphical tangent method, deterministic),
    or "derivative" (midpoint of the steepest |dgamma/dlog10C| adjacent
    pair). All three return the same result fields; only the `method`
    string and the break estimate differ. Invalid method -> ValueError
    naming the allowed values.
    """
    if len(concentrations_mM) != len(surface_tensions_mN_per_m):
        raise ValueError("concentrations and surface tensions must be the same length")
    if len(concentrations_mM) < 2 * min_points_per_segment:
        raise ValueError(
            f"need at least {2*min_points_per_segment} points "
            f"({min_points_per_segment} per segment) to fit a break-point"
        )
    method_choice = method
    if method_choice not in ("rss", "derivative", "tangential"):
        raise ValueError(
            "method must be 'rss', 'derivative', or 'tangential' "
            f"(got {method_choice!r})"
        )
    tangential_cmc_log_c = None
    derivative_cmc_log_c = None

    pairs = sorted(zip(concentrations_mM, surface_tensions_mN_per_m))
    log_c = [math.log10(c) for c, _ in pairs]
    gamma = [g for _, g in pairs]
    n = len(pairs)

    if method_choice == "rss":
        two_seg = _best_two_segment_split(log_c, gamma, n, min_points_per_segment)
        assert two_seg is not None  # guaranteed by the n >= 2*min_points_per_segment check above
        # 3-segment gets its own, possibly-relaxed minimum -- see this
        # function's own docstring for why: only relax when the caller's
        # min_points_per_segment makes 3-segment entirely unreachable for
        # this n, and only down to the smallest floor (2) that still makes
        # it reachable, never below what n actually requires.
        three_seg_min_pts = min_points_per_segment
        if n < 3 * three_seg_min_pts:
            three_seg_min_pts = max(2, n // 3)
        three_seg = _best_three_segment_split(log_c, gamma, n, three_seg_min_pts)

        use_three_seg = False
        if three_seg is not None:
            bic_two = _bic(two_seg[0], n, n_params=2 * 2 + 1)
            bic_three = _bic(three_seg[0], n, n_params=2 * 3 + 2)
            use_three_seg = bic_three < bic_two

        if use_three_seg and three_seg is not None:
            _, split1, split, decline_slope, decline_intercept, decline_rss = three_seg
            n_baseline_points = split1
            method = "three-segment linear regression break-point (baseline + decline + plateau, BIC-selected)"
        else:
            _, split, decline_slope, decline_intercept, decline_rss = two_seg
            split1 = 0
            n_baseline_points = 0
            method = "two-segment linear regression break-point (min total RSS)"
    elif method_choice == "tangential":
        # Tangential-intersection break (Rosen's graphical tangent method,
        # made deterministic): scan every valid pre-line [0:j] and keep the
        # split whose PRE line is steepest (largest |slope|, ties -> last,
        # i.e. the end of the steepest region where the break sits); on
        # ideal segmented data that pre-line IS the premicellar segment.
        # The CMC is the algebraic intersection of that pre-line with the
        # post-line [j:n], read back in log10-concentration space.
        best_j = -1
        pre_slope_t = pre_int_t = pre_rss_t = 0.0
        for j in range(min_points_per_segment, n - min_points_per_segment + 1):
            s, i_, rss_ = _linreg(log_c[:j], gamma[:j])
            # 1e-9 slope tolerance: float noise must not stop a genuine
            # tie from advancing to the LAST (break-side) steepest fit
            if best_j < 0 or abs(s) >= abs(pre_slope_t) - 1e-9:
                best_j, pre_slope_t, pre_int_t, pre_rss_t = j, s, i_, rss_
        if best_j < 0:
            raise ValueError(
                "no valid pre/post split range for method='tangential' -- "
                "need at least 2*min_points_per_segment points"
            )
        split = best_j
        split1 = 0
        decline_slope, decline_intercept, decline_rss = pre_slope_t, pre_int_t, pre_rss_t
        n_baseline_points = 0
        post_slope_t, post_int_t, _ = _linreg(log_c[split:], gamma[split:])
        denom = pre_slope_t - post_slope_t
        if abs(denom) < 1e-12:
            raise ValueError(
                "premicellar and postmicellar tangential fits are parallel (no slope "
                "break) -- this data does not show a surface-tension break point; "
                "check the input rather than trusting this fit"
            )
        tangential_cmc_log_c = (post_int_t - pre_int_t) / denom
        method = "tangential intersection (steepest premicellar fit x postmicellar fit, log10C)"
    else:  # method_choice == "derivative"
        # Maximum first-derivative break: the steepest adjacent-point pair
        # of |dgamma/dlog10C|, ties -> LAST (end of the steepest region).
        # Restricted to pairs leaving >=2 points on each side so the shared
        # reporting regressions stay well-defined.
        best_i = None
        best_mag = -1.0
        for i in range(1, n - 3):
            dx = log_c[i + 1] - log_c[i]
            if dx == 0.0:
                continue
            mag = abs((gamma[i + 1] - gamma[i]) / dx)
            if mag >= best_mag - 1e-9:
                best_i, best_mag = i, mag
        if best_i is None:
            raise ValueError(
                "no adjacent pair with distinct concentrations -- cannot take the "
                "first-derivative break-point (method='derivative'); check the input"
            )
        derivative_cmc_log_c = (log_c[best_i] + log_c[best_i + 1]) / 2.0
        split = best_i + 2
        split1 = 0
        n_baseline_points = 0
        decline_slope, decline_intercept, decline_rss = _linreg(
            log_c[:split], gamma[:split]
        )
        method = ("maximum first-derivative |dgamma/dlog10C| break-point "
                  "(midpoint of steepest adjacent pair)")

    _decline_x, decline_y = log_c[split1:split], gamma[split1:split]
    mean_y = sum(decline_y) / len(decline_y)
    ss_tot = sum((y - mean_y) ** 2 for y in decline_y)
    r_squared = 1.0 - (decline_rss / ss_tot) if ss_tot > 0 else 1.0

    # CMC = concentration at the boundary between the declining segment and
    # the postmicellar (plateau) segment -- the midpoint (in log-concentration
    # space) between the last declining point and the first postmicellar
    # point, matching how a break-point is actually read off a real plot
    # (the intersection region between the two branches, not either branch's
    # own last/first raw data point).
    if tangential_cmc_log_c is not None:
        cmc_log_c = tangential_cmc_log_c
    elif derivative_cmc_log_c is not None:
        cmc_log_c = derivative_cmc_log_c
    else:
        cmc_log_c = (log_c[split - 1] + log_c[split]) / 2.0
    cmc_mM = 10 ** cmc_log_c
    gamma_at_cmc = decline_slope * cmc_log_c + decline_intercept

    # d(gamma)/d(ln C) = d(gamma)/d(log10 C) / ln(10) -- the form
    # gibbs_gamma_max() actually consumes.
    slope_per_lnC = decline_slope / math.log(10)

    # Postmicellar (plateau) segment's own regression -- computed once here
    # rather than threaded through the split-search helpers, since it's not
    # needed during the search itself, only for reporting.
    post_x, post_y = log_c[split:], gamma[split:]
    post_slope, post_intercept, post_rss = _linreg(post_x, post_y)
    post_mean_y = sum(post_y) / len(post_y)
    post_ss_tot = sum((y - post_mean_y) ** 2 for y in post_y)
    r_squared_post = 1.0 - (post_rss / post_ss_tot) if post_ss_tot > 0 else 1.0

    # gamma0 (pure-solvent-like surface tension): the detected flat baseline's own
    # mean if one was found, else the dataset's own lowest-concentration point --
    # real measured data, never a looked-up pure-water constant.
    if n_baseline_points > 0:
        gamma0 = sum(gamma[:split1]) / split1
    else:
        gamma0 = gamma[0]

    # C20 (surfactant efficiency parameter): the concentration at which the
    # premicellar decline line has reduced gamma0 by exactly 20 mN/m -- real,
    # standard definition (Rosen's "Surfactants and Interfacial Phenomena";
    # confirmed 2026-09-12 against an open-access primary source, Table 1 of a
    # real cardanol-surfactant paper reporting pC20 alongside CMC/gamma_CMC from
    # the same curve). Only defined if the premicellar decline actually reaches
    # a 20 mN/m drop before the CMC -- NOT extrapolated past that point.
    c20_mM = None
    pc20 = None
    if (gamma0 - gamma_at_cmc) >= 20.0:
        target_gamma = gamma0 - 20.0
        log_c20 = (target_gamma - decline_intercept) / decline_slope
        c20_mM = 10 ** log_c20
        pc20 = -math.log10(c20_mM / 1000.0)

    return CmcFromCurveResult(
        cmc_mM=cmc_mM,
        premicellar_slope_mN_per_m_per_lnC=slope_per_lnC,
        premicellar_slope_mN_per_m_per_log10C=decline_slope,
        gamma_at_cmc_mN_per_m=gamma_at_cmc,
        r_squared_premicellar=r_squared,
        n_premicellar_points=split - split1,
        n_postmicellar_points=n - split,
        premicellar_x_min_mM=10 ** log_c[split1],
        postmicellar_slope_mN_per_m_per_log10C=post_slope,
        postmicellar_intercept_mN_per_m=post_intercept,
        postmicellar_mean_gamma_mN_per_m=post_mean_y,
        r_squared_postmicellar=r_squared_post,
        n_baseline_points=n_baseline_points,
        gamma0_mN_m=gamma0,
        c20_mM=c20_mM,
        pC20=pc20,
        method=method,
    )


def plateau_is_statistically_suspect(cmc_result: CmcFromCurveResult, ratio_threshold: float = 0.2) -> bool:
    """Real, disclosed limitation of cmc_from_surface_tension_curve's own
    break-point placement, found via GZ-17 in benchmark/groundzero_50
    (ported from SurfactantKit): the postmicellar "plateau" segment can
    have as few as 2 points (an unavoidable consequence of allowing the
    3-segment model to engage at all for the 6-10 point datasets real
    tensiometry curves typically have), and any 2 points fit a line with
    a perfect R^2=1.0 by mathematical construction, REGARDLESS of
    whether that segment is actually flat. This flags exactly that
    situation so a caller (see cmc_via_szyszkowski_deviation below)
    knows to cross-check rather than trust cmc_mM at face value.

    Calibrated against SurfactantKit's own already-validated real cases,
    not guessed: SDS's genuinely flat plateau has
    |postmicellar_slope|/|premicellar_slope| = 0.048; DTAB's and CTAB's
    non-flat plateaus (both later found to give a wrong cmc_mM) have
    0.48 and 0.43 respectively -- a wide, clean gap, so 0.2 sits
    comfortably in between rather than being finely tuned to one case.
    Also flags outright whenever the plateau has fewer than 3 points,
    since 2 points carry zero residual degrees of freedom and cannot
    statistically demonstrate flatness at all, independent of the ratio.
    """
    if cmc_result.n_postmicellar_points < 3:
        return True
    if cmc_result.premicellar_slope_mN_per_m_per_log10C == 0:
        return False
    ratio = abs(cmc_result.postmicellar_slope_mN_per_m_per_log10C) / abs(cmc_result.premicellar_slope_mN_per_m_per_log10C)
    return ratio > ratio_threshold


@dataclass
class SzyszkowskiDeviationCmcResult:
    cmc_mM: float | None  # None when no real deviation point could be located (see reason)
    gamma0_fitted_mN_m: float
    a_fitted: float  # = n_factor * R * T * Gamma_max, the Szyszkowski slope-scale parameter
    b_fitted: float  # = K, the Szyszkowski adsorption equilibrium constant
    gamma_max_from_fit_mol_per_m2: float  # a_fitted converted back to Gamma_max, given system_type -- a real, independent cross-check against gibbs_gamma_max's own (log-linear-slope-based) value
    fit_r_squared: float
    deviation_threshold_mN_per_m: float
    n_points_used_in_fit: int
    reason: str
    no_break_found_in_range: bool = False  # True ONLY when the fit succeeded AND every point
                                              # beyond the premicellar range still tracked it (a
                                              # real "no CMC detectable in this data" finding) --
                                              # distinct from cmc_mM=None due to too few points to
                                              # even fit (n_points_used_in_fit < 3), where nothing
                                              # can honestly be concluded either way. A caller
                                              # (e.g. orchestrate.py) can check this field directly
                                              # instead of string-matching `reason`.


def _fit_szyszkowski(concs: list[float], gammas: list[float], gamma0: float | None = None) -> tuple[float, float, float, float]:
    """Shared fitting engine behind cmc_via_szyszkowski_deviation: fit
    gamma = gamma0 - a*ln(b*C+1) to (concs, gammas) via real iterative
    optimization, not a single first guess -- for ANY candidate b,
    gamma0 and a (or just a, when gamma0 is pinned) are the EXACT
    closed-form OLS solution (this module's own _linreg) for the
    linear-in-x=ln(b*C+1) transformed problem; only b itself needs a
    search, done via the shared _golden_section_minimize helper over
    ln(b), minimizing the resulting RSS globally over the whole search
    range -- a real joint optimization of all free parameters, not "one
    fit is considered final" by construction.

    gamma0: pass a fixed value to fit ONLY a and b with gamma0 pinned.
    Leave None to fit gamma0 freely as well (3 free parameters total).

    Returns (gamma0_used_or_fitted, a, b, r_squared).
    """
    from ._optimize import golden_section_minimize as _golden_section_minimize

    if gamma0 is None:
        def rss_for_ln_b(ln_b: float) -> float:
            b = math.exp(ln_b)
            xs = [math.log(b * c + 1.0) for c in concs]
            _, _, rss = _linreg(xs, gammas)
            return rss

        ln_b_best = _golden_section_minimize(rss_for_ln_b, math.log(1e-4), math.log(1e4))
        b_fit = math.exp(ln_b_best)
        xs = [math.log(b_fit * c + 1.0) for c in concs]
        slope, intercept, rss = _linreg(xs, gammas)
        a_fit, gamma0_fit = -slope, intercept
    else:
        def rss_for_ln_b(ln_b: float) -> float:
            b = math.exp(ln_b)
            xs = [math.log(b * c + 1.0) for c in concs]
            ys = [gamma0 - g for g in gammas]
            sxy = sum(x * y for x, y in zip(xs, ys))
            sxx = sum(x * x for x in xs)
            a = sxy / sxx if sxx > 0 else 0.0
            return sum((y - a * x) ** 2 for x, y in zip(xs, ys))

        ln_b_best = _golden_section_minimize(rss_for_ln_b, math.log(1e-4), math.log(1e4))
        b_fit = math.exp(ln_b_best)
        xs = [math.log(b_fit * c + 1.0) for c in concs]
        ys = [gamma0 - g for g in gammas]
        sxy = sum(x * y for x, y in zip(xs, ys))
        sxx = sum(x * x for x in xs)
        a_fit = sxy / sxx
        gamma0_fit = gamma0
        rss = sum((y - a_fit * x) ** 2 for x, y in zip(xs, ys))

    mean_gamma = sum(gammas) / len(gammas)
    ss_tot = sum((g - mean_gamma) ** 2 for g in gammas)
    r_squared = 1.0 - (rss / ss_tot) if ss_tot > 0 else 1.0
    return gamma0_fit, a_fit, b_fit, r_squared


def cmc_via_szyszkowski_deviation(
    concentrations_mM: list[float],
    surface_tensions_mN_per_m: list[float],
    cmc_result: CmcFromCurveResult,
    system_type: str,
    temperature_K: float = 298.15,
    r_gas: float = 8.314462618,
) -> SzyszkowskiDeviationCmcResult:
    """Cross-check/fallback CMC estimate for when
    plateau_is_statistically_suspect(cmc_result) is True -- i.e. when the
    piecewise-linear break-point's own postmicellar segment is too small
    or not flat enough to trust (ported from SurfactantKit). Fits the
    real Szyszkowski isotherm gamma(C) = gamma0 - a*ln(b*C + 1) through
    the RELIABLE premicellar (baseline+decline) data, extrapolates it,
    and finds the CMC as the concentration where the REAL measured curve
    first departs meaningfully from that smooth model -- rather than
    trying to fit a second, separate "plateau" line from data too sparse
    to prove it's flat.

    Fits gamma0, a AND b jointly and FREELY from the raw data -- a real
    2-parameter nonlinear fit. Fits against exactly the same
    baseline+premicellar points cmc_from_surface_tension_curve itself
    already identified (its n_baseline_points/n_premicellar_points, so
    both methods agree on which points are "the reliable decline"), then
    walks every point ABOVE that range in increasing concentration order
    and reports the first one whose measured surface tension differs
    from the fitted isotherm's own prediction by more than 3x the fit's
    own RMS residual (floored at 0.5 mN/m, a real typical tensiometry
    precision, so a suspiciously perfect fit doesn't produce an
    unrealistically tight threshold) -- a self-calibrating threshold
    from the fit's own noise level, not an arbitrary constant.

    Returns cmc_mM=None (with a stated reason) rather than guessing when
    fewer than 3 baseline+premicellar points are available to fit
    gamma0/a/b meaningfully, or when no point in the given data actually
    deviates from the fitted isotherm (the model tracks the real curve
    throughout the given concentration range -- a real "cannot locate a
    break from this data" case, not an error to paper over).
    """
    from .adsorption import N_FACTOR_BY_SYSTEM_TYPE

    if system_type not in N_FACTOR_BY_SYSTEM_TYPE:
        raise ValueError(
            f"system_type must be one of {sorted(N_FACTOR_BY_SYSTEM_TYPE)}, got {system_type!r} -- "
            "do not guess; if the surfactant's ionic character or electrolyte condition isn't "
            "stated, this deviation check cannot be run"
        )

    pairs = sorted(zip(concentrations_mM, surface_tensions_mN_per_m))
    n_fit = cmc_result.n_baseline_points + cmc_result.n_premicellar_points
    fit_pairs = pairs[:n_fit]
    beyond_pairs = pairs[n_fit:]

    if len(fit_pairs) < 3:
        return SzyszkowskiDeviationCmcResult(
            cmc_mM=None, gamma0_fitted_mN_m=0.0, a_fitted=0.0, b_fitted=0.0,
            gamma_max_from_fit_mol_per_m2=0.0, fit_r_squared=0.0, deviation_threshold_mN_per_m=0.0,
            n_points_used_in_fit=len(fit_pairs),
            reason="fewer than 3 baseline+premicellar points available -- cannot fit gamma0/a/b meaningfully",
        )

    fit_concs = [c for c, _ in fit_pairs]
    fit_gammas = [g for _, g in fit_pairs]

    gamma0_fit, a_fit, b_fit, r_squared = _fit_szyszkowski(fit_concs, fit_gammas, gamma0=None)

    n_factor = N_FACTOR_BY_SYSTEM_TYPE[system_type]
    gamma_max_from_fit = a_fit / (n_factor * r_gas * temperature_K * 1000.0)

    residuals = [g - (gamma0_fit - a_fit * math.log(b_fit * c + 1.0)) for c, g in fit_pairs]
    rms = (sum(r * r for r in residuals) / len(residuals)) ** 0.5
    threshold = max(3.0 * rms, 0.5)

    for c, g in beyond_pairs:
        predicted = gamma0_fit - a_fit * math.log(b_fit * c + 1.0)
        if abs(g - predicted) > threshold:
            return SzyszkowskiDeviationCmcResult(
                cmc_mM=c, gamma0_fitted_mN_m=gamma0_fit, a_fitted=a_fit, b_fitted=b_fit,
                gamma_max_from_fit_mol_per_m2=gamma_max_from_fit, fit_r_squared=r_squared,
                deviation_threshold_mN_per_m=threshold, n_points_used_in_fit=len(fit_pairs),
                reason=f"first point past the baseline+premicellar range where measured tension departs "
                       f"from the fitted Szyszkowski isotherm by more than {threshold:.2f} mN/m",
            )

    return SzyszkowskiDeviationCmcResult(
        cmc_mM=None, gamma0_fitted_mN_m=gamma0_fit, a_fitted=a_fit, b_fitted=b_fit,
        gamma_max_from_fit_mol_per_m2=gamma_max_from_fit, fit_r_squared=r_squared,
        deviation_threshold_mN_per_m=threshold, n_points_used_in_fit=len(fit_pairs),
        reason="no point beyond the baseline+premicellar range deviates from the fitted isotherm -- "
               "cannot locate a break from this data",
        no_break_found_in_range=True,
    )


@dataclass
class CmcFromConductivityResult:
    cmc_mM: float
    slope_below_cmc: float  # d(conductivity)/d(concentration), premicellar segment
    slope_above_cmc: float  # d(conductivity)/d(concentration), postmicellar segment
    counterion_binding_degree: float  # 1 - slope_above/slope_below, the same
                                       # formula thermodynamics.counterion_binding_degree()
                                       # computes -- provided directly here so this
                                       # function alone answers "what's beta from
                                       # this raw conductivity curve"
    r_squared_below_cmc: float
    r_squared_above_cmc: float
    n_points_below_cmc: int
    n_points_above_cmc: int
    method: str = "two-segment linear regression, linear-concentration break-point (min total RSS)"


def cmc_from_conductivity_curve(
    concentrations_mM: list[float],
    conductivities: list[float],
    min_points_per_segment: int = 3,
    method: str = "rss",
) -> CmcFromConductivityResult:
    """Extracts CMC from a raw specific-conductivity-vs-concentration curve
    -- the standard conductometric break-point method for IONIC surfactants
    (nonionic surfactants have negligible conductivity response and cannot
    be characterized this way -- a real, disclosed method limitation, not
    an oversight).

    Unlike cmc_from_surface_tension_curve (which fits gamma vs LOG10
    concentration, per the Gibbs adsorption isotherm's own log-concentration
    dependence), conductometry is fit on LINEAR concentration -- specific
    conductivity increases linearly with surfactant concentration in both
    the premicellar (free ion/monomer) and postmicellar (micelle-dominated)
    regimes, just with two different real slopes, because a micelle is a
    much larger, slower-diffusing, less efficient charge carrier than a
    free monomer/counterion. The CMC is the concentration where the two
    fitted lines intersect -- both slopes are real, physically meaningful
    quantities in their own right (not just a means to locate the CMC),
    which is why this function returns both directly, in exactly the form
    thermodynamics.counterion_binding_degree() consumes, AND the resulting
    beta itself computed the same way that function does (so this one
    function answers the whole "raw conductivity curve -> CMC and beta"
    question end to end).

    Real, standard method -- not invented: segmented/piecewise linear
    regression with the break point taken at the two fitted lines'
    intersection is the textbook conductometric CMC method (see e.g. the
    two-linear-segments description in real surfactant characterization
    literature), same "fit every valid split, take the one minimizing
    total residual sum of squares" discipline as this module's
    cmc_from_surface_tension_curve, adapted for the linear (not log)
    x-axis this specific technique uses.

    method selects the break-point estimator: "rss" (default -- the
    two-segment min-RSS intersection above, unchanged), "jacobs_jones"
    (second-derivative maximum, grid-resolution-limited to the sampled
    concentrations; operator-free second-derivative CMC determination,
    Langmuir 2020, 10.1021/acs.langmuir.0c00420), or "onsager"
    (molar-conductivity kappa/c vs sqrt(concentration) two-line
    intersection -- CMC and beta are invariant to the kappa/c unit
    scale, but the returned r-squared values are on that basis). The
    slope-order guard applies to all three. Invalid method -> ValueError
    naming the allowed values.
    """
    if len(concentrations_mM) != len(conductivities):
        raise ValueError("concentrations and conductivities must be the same length")
    if len(concentrations_mM) < 2 * min_points_per_segment:
        raise ValueError(
            f"need at least {2*min_points_per_segment} points "
            f"({min_points_per_segment} per segment) to fit a break-point"
        )
    method_choice = method
    if method_choice not in ("rss", "jacobs_jones", "onsager"):
        raise ValueError(
            "method must be 'rss', 'jacobs_jones', or 'onsager' "
            f"(got {method_choice!r})"
        )

    pairs = sorted(zip(concentrations_mM, conductivities))
    conc = [c for c, _ in pairs]
    cond = [k for _, k in pairs]
    n = len(pairs)

    if method_choice == "rss":
        best = None
        for split in range(min_points_per_segment, n - min_points_per_segment + 1):
            pre_x, pre_y = conc[:split], cond[:split]
            post_x, post_y = conc[split:], cond[split:]
            pre_slope, pre_intercept, pre_rss = _linreg(pre_x, pre_y)
            post_slope, post_intercept, post_rss = _linreg(post_x, post_y)
            total_rss = pre_rss + post_rss
            if best is None or total_rss < best[0]:
                best = (total_rss, split, pre_slope, pre_intercept, pre_rss, post_slope, post_intercept, post_rss)

        if best is None:
            raise ValueError(
                "two-segment break-point fit needs at least 2 points on each "
                "segment: no valid split found for this concentration grid")
        _, split, pre_slope, pre_intercept, pre_rss, post_slope, post_intercept, post_rss = best
        cmc_mM = None  # set after the slope-order guard, via intersection
        method_label = "two-segment linear regression, linear-concentration break-point (min total RSS)"
    elif method_choice == "jacobs_jones":
        # Second-derivative family: the break is where the slope of the
        # conductivity curve changes most -- discrete d2(kappa)/dc2 as the
        # adjacent-slope jump |s_next - s_prev|; CMC = the sampled
        # concentration at the maximum (grid-resolution-limited), per the
        # operator-free second-derivative CMC determination (Langmuir 2020,
        # 10.1021/acs.langmuir.0c00420; JCIS 2011, 10.1016/j.jcis.2011.12.037
        # lineage). Strict > keeps the FIRST maximum on plateaus.
        best_i = None
        best_mag = -1.0
        for i in range(1, n - 2):
            if conc[i] == conc[i - 1] or conc[i + 1] == conc[i]:
                continue
            s_prev = (cond[i] - cond[i - 1]) / (conc[i] - conc[i - 1])
            s_next = (cond[i + 1] - cond[i]) / (conc[i + 1] - conc[i])
            mag = abs(s_next - s_prev)
            if mag > best_mag:
                best_i, best_mag = i, mag
        if best_i is None:
            raise ValueError(
                "no adjacent pair with distinct concentrations -- cannot take the "
                "second-derivative break-point (method='jacobs_jones'); check the input"
            )
        split = best_i + 1
        pre_x, pre_y = conc[:split], cond[:split]
        post_x, post_y = conc[split:], cond[split:]
        pre_slope, pre_intercept, pre_rss = _linreg(pre_x, pre_y)
        post_slope, post_intercept, post_rss = _linreg(post_x, post_y)
        cmc_mM = conc[best_i]
        method_label = ("second-derivative maximum (discrete slope-jump d2kappa/dc2, "
                        "grid-resolution-limited to sampled concentrations)")
    else:  # method_choice == "onsager"
        # Onsager/Kohlrausch limiting-slope form: molar-conductivity basis
        # y = kappa/c (proportional to Lambda; the absolute unit scale
        # cancels out of both CMC and beta) vs x = sqrt(concentration);
        # two-line min-RSS split, CMC = (intersection x)^2. Real standard
        # practice: plot Lambda vs sqrt(c), the two linear regions below
        # and above the break, intersection read off as the CMC.
        if any(c <= 0 for c in conc):
            raise ValueError(
                "concentrations_mM must all be positive for method='onsager' "
                "(molar conductivity kappa/c is undefined at zero concentration)"
            )
        xs = [math.sqrt(c) for c in conc]
        ys = [k / c for k, c in zip(cond, conc)]
        best = None
        for split in range(min_points_per_segment, n - min_points_per_segment + 1):
            x_pre_s, x_pre_i, x_pre_rss = _linreg(xs[:split], ys[:split])
            x_post_s, x_post_i, x_post_rss = _linreg(xs[split:], ys[split:])
            total_rss = x_pre_rss + x_post_rss
            if best is None or total_rss < best[0]:
                best = (total_rss, split, x_pre_s, x_pre_i, x_pre_rss, x_post_s, x_post_i, x_post_rss)

        if best is None:
            raise ValueError(
                "sqrt(c)-basis break-point fit needs at least 2 points on each "
                "segment: no valid split found for this concentration grid")
        _, split, x_pre_s, x_pre_i, _, x_post_s, x_post_i, _ = best
        x_denom = x_pre_s - x_post_s
        if abs(x_denom) < 1e-12:
            raise ValueError(
                "sqrt(c)-basis premicellar and postmicellar fits are parallel "
                "(no break) -- this data does not show a conductometric break "
                "point; check the input rather than trusting this fit"
            )
        x_int = (x_post_i - x_pre_i) / x_denom
        cmc_mM = x_int * x_int
        # Slope-order guard, beta and r-squared all stay on the LINEAR
        # kappa-vs-c basis (the basis thermodynamics.counterion_binding_degree
        # consumes) -- only the CMC location comes from the sqrt(c)
        # intersection, so beta means the same thing for every method.
        pre_slope, pre_intercept, pre_rss = _linreg(conc[:split], cond[:split])
        post_slope, post_intercept, post_rss = _linreg(conc[split:], cond[split:])
        method_label = ("molar-conductivity (kappa/c) vs sqrt(concentration) two-line "
                        "intersection (Onsager/Kohlrausch limiting-slope form)")

    if pre_slope <= post_slope:
        raise ValueError(
            "fitted premicellar slope is not greater than the postmicellar slope -- this "
            "contradicts the expected conductometric behavior (micelles are less efficient "
            "charge carriers than free monomer/counterions) and likely means this data does "
            "not show a real conductometric break point (e.g. a nonionic surfactant, or data "
            "that doesn't actually span the CMC) -- check the input rather than trusting this fit"
        )

    if cmc_mM is None:
        # CMC = the actual intersection of the two fitted lines, the standard
        # convention for this method (both segments are genuinely linear here,
        # unlike surface tension's curved transition region, so the precise
        # algebraic intersection is the right convention, not a data-point
        # midpoint).
        cmc_mM = (post_intercept - pre_intercept) / (pre_slope - post_slope)

    pre_y = cond[:split]
    mean_pre = sum(pre_y) / len(pre_y)
    ss_tot_pre = sum((y - mean_pre) ** 2 for y in pre_y)
    r_squared_pre = 1.0 - (pre_rss / ss_tot_pre) if ss_tot_pre > 0 else 1.0

    post_y = cond[split:]
    mean_post = sum(post_y) / len(post_y)
    ss_tot_post = sum((y - mean_post) ** 2 for y in post_y)
    r_squared_post = 1.0 - (post_rss / ss_tot_post) if ss_tot_post > 0 else 1.0

    beta = 1.0 - (post_slope / pre_slope)

    return CmcFromConductivityResult(
        cmc_mM=cmc_mM,
        slope_below_cmc=pre_slope,
        slope_above_cmc=post_slope,
        counterion_binding_degree=beta,
        r_squared_below_cmc=r_squared_pre,
        r_squared_above_cmc=r_squared_post,
        n_points_below_cmc=split,
        n_points_above_cmc=n - split,
        method=method_label,
    )


@dataclass
class AggregationNumberQuenchingResult:
    aggregation_number: float
    slope_per_Q: float  # fitted N/(Ct - cmc), from ln(I0/I) = slope * [Q]
    method: str = "steady-state fluorescence quenching (Turro-Yekta, 1978)"


def aggregation_number_from_quenching_curve(
    quencher_concentrations: list[float],
    intensities: list[float],
    intensity_without_quencher: float,
    total_surfactant_concentration: float,
    cmc: float,
) -> AggregationNumberQuenchingResult:
    """Micelle aggregation number from a raw steady-state fluorescence
    quenching (SSFQ) series -- the Turro-Yekta method (1978), verified
    directly from a real fetched source (equation image, quoting exactly):

    ln(I / I0) = -N*[Q] / (Ct - cmc)

    i.e. ln(I0/I) = N*[Q] / (Ct - cmc), a straight line through the
    origin when ln(I0/I) is plotted against quencher concentration [Q].

    Deliberately the SIMPLEST real aggregation-number method to run:
    needs only a standard steady-state spectrofluorometer (not the
    time-resolved/TCSPC instrumentation TRFQ requires, not SANS, not an
    ultracentrifuge) -- a fluorescent probe (commonly pyrene) solubilized
    in the micelles, quenched by a compound that partitions almost
    entirely into the micellar phase (e.g. cetylpyridinium chloride,
    benzophenone), measured at a series of quencher concentrations.

    quencher_concentrations, intensities: paired raw data -- the probe's
    fluorescence intensity (arbitrary units, e.g. the first vibronic
    peak height for pyrene) at each quencher concentration [Q] (same
    concentration unit as total_surfactant_concentration/cmc).
    intensity_without_quencher: I0, the probe's intensity at [Q]=0
    (measured separately, the same units as `intensities`).
    total_surfactant_concentration, cmc: same concentration unit as
    quencher_concentrations -- Ct must be well above cmc for a
    meaningful (Ct - cmc) micellized-surfactant denominator.

    The slope is fit forced through the origin (least-squares with zero
    intercept: slope = sum(x*y)/sum(x^2)), matching the theoretical
    model exactly (ln(I0/I0)=0 at [Q]=0 identically, not just
    approximately) rather than fitting a free intercept that could
    absorb real experimental baseline noise into a spurious offset.

    Source: Turro & Yekta, J. Am. Chem. Soc. 100 (1978) 5951-5952 (the
    original method); equation verified directly from a real fetched
    secondary source (Bhattacharjee & Mahapatra-style DDM/SAXS
    comparison paper, SCIRP, image-extracted equation), cross-confirmed
    by two independent web summaries of the same equation before
    fetching the image directly.
    """
    if len(quencher_concentrations) != len(intensities):
        raise ValueError("quencher_concentrations and intensities must be the same length")
    if len(quencher_concentrations) < 2:
        raise ValueError("need at least 2 (quencher concentration, intensity) points to fit a slope")
    if intensity_without_quencher <= 0:
        raise ValueError("intensity_without_quencher must be positive")
    if any(i <= 0 for i in intensities):
        raise ValueError("all intensities must be positive")
    if total_surfactant_concentration <= cmc:
        raise ValueError("total_surfactant_concentration must exceed cmc (need a real micellized fraction)")

    xs = quencher_concentrations
    ys = [math.log(intensity_without_quencher / i) for i in intensities]
    sum_xy = sum(x * y for x, y in zip(xs, ys))
    sum_xx = sum(x * x for x in xs)
    if sum_xx == 0:
        raise ValueError("quencher_concentrations cannot all be zero")
    slope = sum_xy / sum_xx

    aggregation_number = slope * (total_surfactant_concentration - cmc)
    return AggregationNumberQuenchingResult(aggregation_number=aggregation_number, slope_per_Q=slope)


AVOGADRO_NUMBER = 6.02214076e23  # /mol, exact SI 2019

# Tier 2 bottleneck-resolution addition, 2026-09-15 (see benchmark/
# audit findings, 2026-09 item 10): real,
# published dn/dc (refractive index increment) values for common
# surfactant micelles in water, closing this bottleneck for compounds
# where a real measured value already exists -- for
# aggregation_number_from_sls_debye_plot's own dn_dc_mL_per_g input.
# Source: Malvern Panalytical, "Refractive index increment dn/dc
# values" (malvernpanalytical.com/en/learn/knowledge-center/insights/
# refractive-index-increment-dndc-values), live-fetched 2026-09-15,
# 632.8 nm (He-Ne red laser), 25 C, water as solvent -- itself drawing
# on the standard compiled reference (Theisen, Johann, Deacon & Harding,
# "Refractive Increment Data-Book for Polymer and Biomolecular
# Scientists," Nottingham University Press, 2000, ISBN 1-897676-29-8).
# REVISED 2026-09-15 (same day): Triton X-100's real dn/dc was located
# in a real primary source the user provided directly (Stubicar,
# Matejas, Zipper & Wilfing, in K.L. Mittal (ed.), "Surfactants in
# Solution," Plenum Press, New York (1989) 181-193 -- their own
# instrument-calibration section states it explicitly, stable across
# both pure water and KCl/KBr/KI electrolyte solutions), closing that
# part of the earlier-disclosed gap. Tween-80 CLOSED same day: Jiang,
# F., PhD dissertation, Virginia Tech (2011), "Effects of the Non-ionic
# Surfactant Tween 80 on the Enzymatic Hydrolysis of Model Cellulose and
# Lignocellulosic Substrates" -- states directly, in its own SPR
# (surface plasmon resonance) methods section, a real dn/dc measured
# with a Wyatt Optilab rEX differential refractometer, primary PDF read
# in full 2026-09-15.
DN_DC_REFERENCE_ML_PER_G = {
    "SDS": {"dn_dc": 0.11, "solvent": "water", "wavelength_nm": 632.8, "temperature_C": 25.0,
            "source": "Malvern Panalytical dn/dc reference page, live-fetched 2026-09-15"},
    "CTAB": {"dn_dc": 0.15, "solvent": "water", "wavelength_nm": 632.8, "temperature_C": 25.0,
             "source": "Malvern Panalytical dn/dc reference page, live-fetched 2026-09-15"},
    "TRITONX100": {"dn_dc": 0.140, "dn_dc_uncertainty": 0.005, "solvent": "water (also stable in KCl/KBr/KI)",
                    "wavelength_nm": 546.0, "temperature_C": 20.0,
                    "source": "Stubicar, Matejas, Zipper & Wilfing, in Mittal (ed.), Surfactants in Solution, "
                              "Plenum Press (1989) 181-193, primary PDF read in full 2026-09-15"},
    "TWEEN80": {"dn_dc": 0.132, "solvent": "water", "wavelength_nm": 690.0, "temperature_C": 25.0,
                "source": "Jiang, PhD dissertation, Virginia Tech (2011), Wyatt Optilab rEX differential "
                          "refractometer, primary PDF read in full 2026-09-15"},
}


@dataclass
class AggregationNumberSLSResult:
    micelle_molar_mass_g_per_mol: float
    aggregation_number: float
    second_virial_coefficient_cm3_mol_per_g2: float
    method: str = "static light scattering, single-angle Debye plot"


def aggregation_number_from_sls_debye_plot(
    micellized_concentrations_g_per_mL: list[float],
    rayleigh_ratios_cm_inv: list[float],
    dn_dc_mL_per_g: float,
    wavelength_nm: float,
    monomer_molar_mass_g_per_mol: float,
    refractive_index_solvent: float = 1.333,
) -> AggregationNumberSLSResult:
    """Micelle aggregation number from a raw static light scattering
    (SLS) concentration series -- the single-angle Debye plot method,
    real and standard for small (Rg < 12 nm) scatterers like surfactant
    micelles (a multi-angle instrument, needed for larger macromolecules,
    is NOT required here) -- more instrumentation than SSFQ
    (aggregation_number_from_quenching_curve) but real, standard, and
    achievable on the same class of light-scattering/DLS instrument
    already used elsewhere in this line of work for micelle sizing.

    Debye equation (verified directly against a real fetched primary
    technical source, Brookhaven Instruments' "SLS FAQ: The Debye
    Plot," which reduces the general Zimm equation to this form
    explicitly when Rg^2*q^2/3 << 1, the small-scatterer limit):

    K*c / DeltaR = 1/Mw + 2*A2*c

    K = 4*pi^2*n^2*(dn/dc)^2 / (NA * lambda0^4)   (optical constant;
    n = solvent refractive index, dn/dc = refractive index increment,
    lambda0 = VACUUM wavelength, NA = Avogadro's number -- this specific
    4*pi^2/vacuum-wavelength form cross-confirmed across three
    independent real sources: Malvern, Anton Paar, and Brookhaven
    technical documentation, all agreeing).

    micellized_concentrations_g_per_mL: concentration of MICELLIZED
    surfactant only (i.e. total surfactant mass concentration minus the
    mass-concentration-equivalent of the cmc) at each measurement point
    -- using total surfactant concentration instead would include free
    monomer that contributes no excess scattering above the solvent
    baseline, biasing the fitted molar mass. Same "must be a real,
    already-corrected input, not guessed" discipline as elsewhere in
    this module.
    rayleigh_ratios_cm_inv: the excess Rayleigh ratio (DeltaR, solvent-
    subtracted, units of inverse cm) at each corresponding concentration
    -- a real instrument-calibrated quantity, not raw photon counts.
    dn_dc_mL_per_g: refractive index increment for this surfactant in
    this solvent at this wavelength -- system-specific, must be a real
    measured or cited value, never guessed (see Zhao, Brown & Schuck,
    Biophys. J. 100 (2011) 2309-2317 for how variable this can be even
    within one macromolecule class).
    wavelength_nm: VACUUM wavelength of the laser used, in nm.
    monomer_molar_mass_g_per_mol: pure-surfactant-monomer molar mass --
    aggregation_number = micelle_molar_mass / monomer_molar_mass.
    refractive_index_solvent: default 1.333 (water, visible light,
    ~20-25 C) -- override for other solvents/temperatures/wavelengths.

    Fits Mw and the second virial coefficient A2 by ordinary least
    squares on Kc/DeltaR vs. c (the standard Debye-plot linear fit),
    reusing this module's own _linreg helper.
    """
    if len(micellized_concentrations_g_per_mL) != len(rayleigh_ratios_cm_inv):
        raise ValueError("micellized_concentrations_g_per_mL and rayleigh_ratios_cm_inv must be the same length")
    if len(micellized_concentrations_g_per_mL) < 2:
        raise ValueError("need at least 2 (concentration, Rayleigh ratio) points to fit the Debye plot")
    if any(c <= 0 for c in micellized_concentrations_g_per_mL):
        raise ValueError("all micellized concentrations must be positive")
    if any(r <= 0 for r in rayleigh_ratios_cm_inv):
        raise ValueError("all Rayleigh ratios must be positive")
    if dn_dc_mL_per_g == 0:
        raise ValueError("dn_dc_mL_per_g must be nonzero -- a real, measured/cited value, not a guess")
    if monomer_molar_mass_g_per_mol <= 0:
        raise ValueError("monomer_molar_mass_g_per_mol must be positive")

    wavelength_cm = wavelength_nm * 1e-7
    optical_constant_K = (
        4.0 * math.pi ** 2 * refractive_index_solvent ** 2 * dn_dc_mL_per_g ** 2
        / (AVOGADRO_NUMBER * wavelength_cm ** 4)
    )

    xs = micellized_concentrations_g_per_mL
    ys = [optical_constant_K * c / r for c, r in zip(micellized_concentrations_g_per_mL, rayleigh_ratios_cm_inv)]
    slope, intercept, _ = _linreg(xs, ys)

    if intercept <= 0:
        raise ValueError(
            "fitted Kc/DeltaR intercept is non-positive -- this gives a non-physical (infinite or negative) "
            "molar mass; check the input data (Rayleigh ratios, dn/dc, wavelength) rather than trusting this fit"
        )
    micelle_molar_mass = 1.0 / intercept
    aggregation_number = micelle_molar_mass / monomer_molar_mass_g_per_mol
    second_virial_coefficient = slope / 2.0

    return AggregationNumberSLSResult(
        micelle_molar_mass_g_per_mol=micelle_molar_mass,
        aggregation_number=aggregation_number,
        second_virial_coefficient_cm3_mol_per_g2=second_virial_coefficient,
    )


@dataclass
class AggregationNumberDLSResult:
    hydrodynamic_radius_nm: float
    micelle_molar_mass_g_per_mol: float
    aggregation_number: float
    method: str = "dynamic light scattering, hydrodynamic-volume/partial-specific-volume route"


def aggregation_number_from_dls(
    diffusion_coefficient_cm2_per_s: float,
    viscosity_mPas: float,
    partial_specific_volume_cm3_per_g: float,
    monomer_molar_mass_g_per_mol: float,
    temperature_K: float = 298.15,
) -> AggregationNumberDLSResult:
    """Micelle aggregation number from a raw DLS-measured diffusion
    coefficient -- an alternative-methods-sweep addition (2026-09-10),
    real and standard (the same instrument class already used elsewhere
    in this line of work for micelle sizing): converts the hydrodynamic
    radius (dynamics.hydrodynamic_radius_stokes_einstein, reused
    directly here, not reimplemented) to a hydrodynamic VOLUME, then to
    a molar mass via the micelle's partial specific volume, then to an
    aggregation number via the monomer molar mass.

        R_h = k_B*T / (6*pi*eta*D)                 [Stokes-Einstein]
        V_h = (4/3)*pi*R_h^3                        [per-particle hydrodynamic volume]
        M_micelle = V_h * N_A / v_bar                [v_bar = partial specific volume]
        N_agg = M_micelle / M_monomer

    Real, disclosed limitation (standard in the DLS/micelle-sizing
    literature, not specific to this implementation): R_h includes a
    thin bound-water hydration shell, so the hydrodynamic volume is a
    real, known OVERESTIMATE of the anhydrous molecular volume implied
    by partial_specific_volume_cm3_per_g -- this gives an aggregation
    number that is a real, usable estimate but typically biased slightly
    HIGH relative to a hydration-corrected value; treat as an estimate
    of the same character as aggregation_number_spherical's geometric
    estimate, not a hydration-corrected exact value. Same units and
    inputs as dynamics.hydrodynamic_radius_stokes_einstein for the first
    two arguments -- chain directly from that function's own docstring
    guidance on diffusion coefficient units (cm^2/s, not m^2/s).

    partial_specific_volume_cm3_per_g: real, measured or cited value
    for this surfactant (e.g. from densitometry) -- system-specific, do
    not guess (same discipline as dn/dc in aggregation_number_from_sls_debye_plot).
    """
    if diffusion_coefficient_cm2_per_s <= 0:
        raise ValueError("diffusion_coefficient_cm2_per_s must be positive")
    if viscosity_mPas <= 0:
        raise ValueError("viscosity_mPas must be positive")
    if partial_specific_volume_cm3_per_g <= 0:
        raise ValueError("partial_specific_volume_cm3_per_g must be positive")
    if monomer_molar_mass_g_per_mol <= 0:
        raise ValueError("monomer_molar_mass_g_per_mol must be positive")
    if temperature_K <= 0:
        raise ValueError("temperature_K must be positive")

    k_B = 1.380649e-23  # J/K, exact SI 2019
    eta_SI = viscosity_mPas * 1e-3  # mPa.s -> Pa.s
    D_SI = diffusion_coefficient_cm2_per_s * 1e-4  # cm^2/s -> m^2/s
    r_h_m = (k_B * temperature_K) / (6.0 * math.pi * eta_SI * D_SI)
    r_h_nm = r_h_m * 1e9

    r_h_cm = r_h_m * 100.0  # m -> cm, to get volume directly in cm^3
    v_h_cm3 = (4.0 / 3.0) * math.pi * r_h_cm ** 3
    micelle_molar_mass = v_h_cm3 * AVOGADRO_NUMBER / partial_specific_volume_cm3_per_g
    aggregation_number = micelle_molar_mass / monomer_molar_mass_g_per_mol

    return AggregationNumberDLSResult(
        hydrodynamic_radius_nm=r_h_nm,
        micelle_molar_mass_g_per_mol=micelle_molar_mass,
        aggregation_number=aggregation_number,
    )


def estimate_partial_specific_volume_from_tail_and_headgroup(
    n_carbons: int,
    monomer_molar_mass_g_per_mol: float,
    headgroup_molar_volume_cm3_per_mol: float,
) -> float:
    """Real, buildable [COMPUTE] estimate of the WHOLE monomer's partial
    specific volume v_bar (cm^3/g) -- the input
    aggregation_number_from_svedberg_equation and aggregation_number_from_dls
    both require as an external, fully-formed, system-specific value (see
    audit findings, 2026-09 item 11).
    This REDUCES, not eliminates, that external requirement: the
    hydrophobic TAIL's molar-volume contribution is computed exactly
    from cpp.py's already-verified tanford_tail_volume (an EXACT unit
    conversion, cubic-Angstrom/molecule -> cm^3/mol via Avogadro's
    number, of an already-sourced constant -- not a new empirical fit),
    leaving only the HEADGROUP's own molar volume as the remaining real,
    external input.

        V_tail_molar_cm3 = tanford_tail_volume(n_carbons) * N_A * 1e-24
        v_bar = (V_tail_molar_cm3 + headgroup_molar_volume_cm3_per_mol) / monomer_molar_mass_g_per_mol

    headgroup_molar_volume_cm3_per_mol: the real, literature-sourced (or
    independently measured) partial molar volume of the headgroup PLUS
    counterion in water (cm^3/mol) -- deliberately NOT computed here from
    geometric/RDKit volume, because ionic headgroups undergo real
    ELECTROSTRICTION: local water structuring around a charged group
    measurably shrinks its effective partial molar volume below its
    simple geometric size, and can even make it NEGATIVE for small,
    highly-charged ions (e.g. Na+ in water is well documented near
    -1 to -2 cm^3/mol in standard ionic partial-molar-volume tables).
    This value must still come from a real ionic/apparent-molar-volume
    source -- do not derive it from geometry, and do not assume it must
    be positive (the function only requires the FINAL v_bar to be
    positive, not the headgroup term on its own).

    No specific real headgroup value is hardcoded anywhere in this
    library -- unlike VOCG_STANDARD_LIQUIDS in wetting.py, partial molar
    volumes of specific surfactant headgroups were not verified against a
    primary source with enough confidence this session to hardcode; the
    real, exact part shipped here is the tail-volume unit conversion
    alone. Validated by round trip only (construct a known v_bar's
    consistent headgroup_molar_volume_cm3_per_mol via this same formula,
    recover it back exactly) -- see tests/test_curve_analysis.py.
    """
    if n_carbons < 1:
        raise ValueError("n_carbons must be at least 1")
    if monomer_molar_mass_g_per_mol <= 0:
        raise ValueError("monomer_molar_mass_g_per_mol must be positive")
    v_tail_A3 = tanford_tail_volume(n_carbons)
    v_tail_molar_cm3 = v_tail_A3 * AVOGADRO_NUMBER * 1e-24  # A^3/molecule -> cm^3/mol
    v_bar = (v_tail_molar_cm3 + headgroup_molar_volume_cm3_per_mol) / monomer_molar_mass_g_per_mol
    if v_bar <= 0:
        raise ValueError(
            f"computed v_bar={v_bar:.4f} cm^3/g is non-positive -- check "
            "headgroup_molar_volume_cm3_per_mol and monomer_molar_mass_g_per_mol; "
            "a real partial specific volume for a surfactant monomer must be positive "
            "even when the headgroup's own contribution is negative (electrostriction)"
        )
    return v_bar


@dataclass
class AggregationNumberSvedbergResult:
    micelle_molar_mass_g_per_mol: float
    aggregation_number: float
    method: str = "analytical ultracentrifugation, Svedberg equation"


def aggregation_number_from_svedberg_equation(
    sedimentation_coefficient_S: float,
    diffusion_coefficient_cm2_per_s: float,
    partial_specific_volume_cm3_per_g: float,
    solvent_density_g_per_cm3: float,
    monomer_molar_mass_g_per_mol: float,
    temperature_K: float = 298.15,
) -> AggregationNumberSvedbergResult:
    """Micelle aggregation number from raw analytical-ultracentrifugation
    data (sedimentation velocity) -- an alternative-methods-sweep
    addition (2026-09-10), the classic, standard route (Svedberg
    equation) to a macromolecule's molar mass, real and independent of
    the light-scattering/DLS-based routes elsewhere in this module
    (different instrument, different physical principle -- a genuine
    cross-check when both are available on the same system).

        M = R*T*s / [D*(1 - v_bar*rho)]
        N_agg = M / M_monomer

    s: sedimentation coefficient in SVEDBERGS (S, 1 S = 1e-13 s) -- the
    standard reporting unit, converted internally; do not pass raw
    seconds here. D: translational diffusion coefficient in cm^2/s (same
    convention as elsewhere in this module/dynamics.py). v_bar: partial
    specific volume (cm^3/g, real measured/cited value, system-specific
    -- do not guess; same discipline as aggregation_number_from_dls).
    rho: solvent density (g/cm^3, ~0.997 for water at 25 C). The
    (1 - v_bar*rho) "buoyancy factor" must be positive -- a real,
    physical requirement (a particle denser than its buoyant
    displacement sediments; one that isn't, doesn't; this function
    raises rather than silently returning a nonsensical negative or
    infinite mass if the inputs violate it).
    """
    if diffusion_coefficient_cm2_per_s <= 0:
        raise ValueError("diffusion_coefficient_cm2_per_s must be positive")
    if partial_specific_volume_cm3_per_g <= 0:
        raise ValueError("partial_specific_volume_cm3_per_g must be positive")
    if solvent_density_g_per_cm3 <= 0:
        raise ValueError("solvent_density_g_per_cm3 must be positive")
    if monomer_molar_mass_g_per_mol <= 0:
        raise ValueError("monomer_molar_mass_g_per_mol must be positive")
    if temperature_K <= 0:
        raise ValueError("temperature_K must be positive")
    buoyancy_factor = 1.0 - partial_specific_volume_cm3_per_g * solvent_density_g_per_cm3
    if buoyancy_factor <= 0:
        raise ValueError(
            "buoyancy factor (1 - v_bar*rho) is non-positive -- the particle would not sediment under "
            "these inputs (check partial_specific_volume_cm3_per_g and solvent_density_g_per_cm3); "
            "the Svedberg equation gives a non-physical result and is not applied"
        )

    R_GAS = 8.314462618  # J/(mol.K)
    s_seconds = sedimentation_coefficient_S * 1e-13
    D_SI = diffusion_coefficient_cm2_per_s * 1e-4  # cm^2/s -> m^2/s
    micelle_molar_mass_SI = (R_GAS * temperature_K * s_seconds) / (D_SI * buoyancy_factor)  # kg/mol
    micelle_molar_mass = micelle_molar_mass_SI * 1000.0  # kg/mol -> g/mol
    aggregation_number = micelle_molar_mass / monomer_molar_mass_g_per_mol

    return AggregationNumberSvedbergResult(
        micelle_molar_mass_g_per_mol=micelle_molar_mass,
        aggregation_number=aggregation_number,
    )




# ---------------------------------------------------------------------------
# Bootstrap confidence intervals for sigmoidal CMC curves (computational
# worker I2).  A point CMC without an uncertainty interval overstates what a
# 10-20 point curve supports; this worker resamples residuals and refits the
# sigmoid thousands of times -- infeasible by hand, impossible by guessing.
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class BootstrapCMCResult:
    """Point estimates and percentile confidence intervals from a residual
    bootstrap over a 4-parameter sigmoid CMC model."""
    cmc_point: float
    cmc_ci_low: float
    cmc_ci_high: float
    log10_cmc_point: float
    log10_cmc_ci_low: float
    log10_cmc_ci_high: float
    transition_property_point: float  # fitted curve value at the midpoint
    transition_property_ci_low: float
    transition_property_ci_high: float
    ci_level: float
    n_boot_requested: int
    n_boot_used: int
    n_failed: int
    converged: bool
    r_squared: float
    method: str
    params: dict
    caveat: str


def _sigmoid4(x, x0, w, lo, hi):
    """4-parameter logistic (Eissler-type) transition:
    value = lo + (hi-lo)*0.5*(1+tanh((x-x0)/w)).  For surface tension,
    lo/hi are the pre/post plateaus and x0 (the midpoint) is log10(cmc)."""
    import numpy as np
    return lo + (hi - lo) * 0.5 * (1.0 + np.tanh((x - x0) / w))


def _fit_sigmoid4(x, y, p0=None):
    """Least-squares fit with plateau-derived initialization.

    Returns (popt tuple (x0,w,lo,hi), r_squared).  Raises ValueError when no
    transition exists (degenerate curve).
    """
    import numpy as np
    from scipy.optimize import least_squares

    n = len(x)
    order = np.argsort(x)
    xs, ys = x[order], y[order]
    k = max(2, n // 4)
    lo0, hi0 = float(ys[:k].mean()), float(ys[-k:].mean())
    if abs(hi0 - lo0) < 1e-12 * max(1.0, abs(ys).max()):
        raise ValueError(
            "curve shows no pre/post transition (plateau levels identical) -- "
            "a CMC cannot be bootstrapped from a flat curve"
        )
    dy = np.diff(ys)
    dx = np.diff(xs)
    slopes = dy / np.where(np.abs(dx) < 1e-300, 1e-300, dx)
    i = int(np.argmax(np.abs(slopes)))
    x00 = float(0.5 * (xs[i] + xs[i + 1]))
    span = float(xs[-1] - xs[0]) or 1.0
    w0 = span / 8.0
    if p0 is None:
        p0 = (x00, w0, lo0, hi0)
    x_lo, x_hi = float(xs[0] - span), float(xs[-1] + span)

    def resid(p):
        return _sigmoid4(x, *p) - y

    sol = least_squares(
        resid, p0,
        bounds=([x_lo, 1e-4, -np.inf, -np.inf], [x_hi, span * 10.0, np.inf, np.inf]),
        max_nfev=2000,
    )
    popt = tuple(sol.x)
    ss_res = float(np.sum(sol.fun ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return popt, r2


def bootstrap_cmc(
    x: list[float],
    y: list[float],
    n_boot: int = 1000,
    ci: float = 0.95,
    seed: int = 0,
) -> BootstrapCMCResult:
    """Point CMC + percentile confidence interval via residual bootstrap.

    Fits a 4-parameter sigmoid (Eissler-type) to y vs x (x on a log scale,
    e.g. log10 concentration), then resamples residuals with replacement and
    refits ``n_boot`` times.  ``cmc = 10**x0`` where x0 is the transition
    midpoint -- the standard sigmoidal CMC estimator for surface tension or
    spectroscopic probe curves.  Rising curves (e.g. conductivity-like) are
    handled by allowing an inverted transition.

    Catches the two failure modes naive CMC reporting hides: interval width
    (sparse/noisy curves give wide CIs -- check them) and fit failures
    (``n_failed``; >10% marks ``converged=False``).

    Requires scipy (surfactantkit[calc]).
    """
    import numpy as np

    if n_boot < 50:
        raise ValueError("n_boot must be >= 50 for a usable percentile interval")
    if not 0.5 < ci < 1.0:
        raise ValueError("ci must be a confidence level strictly between 0.5 and 1.0")
    xa = np.asarray(x, dtype=float)
    ya = np.asarray(y, dtype=float)
    if xa.ndim != 1 or ya.ndim != 1 or len(xa) != len(ya):
        raise ValueError("x and y must be 1-D and of equal length")
    if len(xa) < 8:
        raise ValueError(
            f"need >= 8 points to bootstrap a 4-parameter sigmoid "
            f"(got {len(xa)}) -- collect more concentrations"
        )
    if np.unique(xa).size < 6:
        raise ValueError("x must have >= 6 distinct values")
    if not np.all(np.isfinite(xa)) or not np.all(np.isfinite(ya)):
        raise ValueError("x and y must be finite (NaN/inf rejected -- clean the data)")

    popt, r2 = _fit_sigmoid4(xa, ya)
    x0, w, lo, hi = popt
    y_hat = _sigmoid4(xa, *popt)
    resid = ya - y_hat

    rng = np.random.default_rng(seed)
    boot_x0, boot_mid = [], []
    n_failed = 0
    for _ in range(n_boot):
        y_star = y_hat + rng.choice(resid, size=len(resid), replace=True)
        try:
            p, _ = _fit_sigmoid4(xa, y_star, p0=popt)
            # reject non-physical widths (solver collapsed onto bounds)
            if not (1e-4 < p[1] < (xa.max() - xa.min()) * 10 + 1e-9):
                n_failed += 1
                continue
            boot_x0.append(p[0])
            boot_mid.append(0.5 * (p[2] + p[3]))
        except Exception:  # noqa: BLE001 -- a bad bootstrap resample can
            # fail the fit via any numerical exception (non-convergence,
            # singular matrix, etc.); count as failed and keep resampling.
            n_failed += 1

    n_used = len(boot_x0)
    if n_used < int(0.9 * n_boot):
        # still return the interval, but flag it
        converged = False
    else:
        converged = True
    if n_used == 0:
        raise RuntimeError(
            "all bootstrap refits failed -- the curve is likely too noisy or "
            "too short to support any CMC uncertainty estimate"
        )
    a = (1.0 - ci) / 2.0
    lo_x0, hi_x0 = np.percentile(boot_x0, [100 * a, 100 * (1 - a)])
    lo_mid, hi_mid = np.percentile(boot_mid, [100 * a, 100 * (1 - a)])
    mid_point = 0.5 * (lo + hi)

    return BootstrapCMCResult(
        cmc_point=float(10.0 ** x0),
        cmc_ci_low=float(10.0 ** lo_x0),
        cmc_ci_high=float(10.0 ** hi_x0),
        log10_cmc_point=float(x0),
        log10_cmc_ci_low=float(lo_x0),
        log10_cmc_ci_high=float(hi_x0),
        transition_property_point=float(mid_point),
        transition_property_ci_low=float(lo_mid),
        transition_property_ci_high=float(hi_mid),
        ci_level=float(ci),
        n_boot_requested=int(n_boot),
        n_boot_used=int(n_used),
        n_failed=int(n_failed),
        converged=bool(converged),
        r_squared=float(r2),
        method="4-parameter sigmoid (Eissler-type) + residual bootstrap "
               "percentile interval",
        params={"x0_log10_conc": float(x0), "width_w": float(w),
                "plateau_low": float(lo), "plateau_high": float(hi)},
        caveat=(
            "Percentile bootstrap interval quantifies sampling uncertainty of "
            "the fitted midpoint only; it does not include model-form error "
            "(a true break-point two-line CMC may differ from the sigmoid "
            "midpoint) and assumes independent, homoscedastic residuals."
        ),
    )




# ---------------------------------------------------------------------------
# DLS size-distribution inversion (computational worker P3).  The ill-posed
# Laplace inversion g1(tau) = sum_j c_j exp(-Gamma_j tau) dlnGamma with
# non-negativity + L-curve regularization cannot be done by hand; text-only
# models either refuse or fabricate "z-average" numbers.
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class DLSInversionResult:
    """Regularized non-negative inversion of a DLS field autocorrelation
    into a hydrodynamic-size distribution."""
    diameter_nm: tuple[float, ...]
    intensity_fraction: tuple[float, ...]
    dh_intensity_mean_nm: float
    dh_z_avg_nm: float
    pdi_from_distribution: float
    q_nm_inv: float
    lambda_reg: float
    residual_norm: float
    relative_residual: float
    n_grid: int
    converged: bool
    method: str
    params: dict
    caveat: str


def _q_from_geometry(wavelength_nm: float, angle_deg: float,
                     refractive_index: float) -> float:
    """Scattering vector magnitude q = (4 pi n / lambda0) sin(theta/2) [nm^-1]."""
    import math as _m
    return (4.0 * _m.pi * refractive_index / wavelength_nm) * _m.sin(
        _m.radians(angle_deg) / 2.0)


def dls_size_distribution(
    g1: list[float],
    tau_s: list[float],
    q_nm_inv: float | None = None,
    wavelength_nm: float = 632.8,
    angle_deg: float = 90.0,
    refractive_index: float = 1.331,
    viscosity_mPa_s: float = 0.89,
    temperature_K: float = 298.15,
    n_gamma: int = 60,
    lambda_reg: float | None = None,
) -> DLSInversionResult:
    """Invert a DLS field autocorrelation into a size distribution.

    Model: g1(tau) = sum_j c_j exp(-Gamma_j tau) dln(Gamma) on a log-spaced
    Gamma grid (Gamma = D q^2), coefficients c_j >= 0.  Solution minimizes
    ||K c - g1||^2 + lambda ||D2 c||^2 (second-difference smoothing) with
    non-negativity (scipy NNLS); lambda is chosen by the L-curve corner
    unless given.  D is converted to hydrodynamic diameter by
    Stokes-Einstein d = kT/(3 pi eta D).

    Reported distribution is INTENSITY-weighted (as DLS always is);
    ``dh_z_avg_nm`` is the intensity-weighted mean diameter.

    Requires scipy (surfactantkit[calc]).
    """
    import numpy as np
    from scipy.optimize import nnls

    g = np.asarray(g1, dtype=float)
    tau = np.asarray(tau_s, dtype=float)
    if g.ndim != 1 or tau.ndim != 1 or g.size != tau.size:
        raise ValueError("g1 and tau_s must be 1-D and of equal length")
    if g.size < 8:
        raise ValueError("need >= 8 correlation points to invert")
    if np.any(tau <= 0) or not np.all(np.isfinite(tau)) or not np.all(np.isfinite(g)):
        raise ValueError("tau_s must be positive and finite; g1 finite (clean the data)")
    if np.allclose(g, g[0]):
        raise ValueError("g1 is flat -- no relaxation to invert")
    if n_gamma < 20:
        raise ValueError("n_gamma must be >= 20")
    if viscosity_mPa_s <= 0 or temperature_K <= 0 or refractive_index <= 0:
        raise ValueError("viscosity, temperature, refractive_index must be positive")
    order = np.argsort(tau)
    tau = tau[order]
    g = g[order]

    if q_nm_inv is None:
        q_nm_inv = _q_from_geometry(wavelength_nm, angle_deg, refractive_index)
    if q_nm_inv <= 0:
        raise ValueError("q_nm_inv must be positive")

    g_min = 0.3 / float(tau.max())
    g_max = 3.0 / float(tau.min())
    gamma = np.logspace(np.log10(g_min), np.log10(g_max), n_gamma)

    K = np.exp(-np.outer(tau, gamma))          # (n_tau, n_gamma)
    # second-difference smoothness operator
    D2 = np.zeros((n_gamma - 2, n_gamma))
    for i in range(n_gamma - 2):
        D2[i, i:i + 3] = (1.0, -2.0, 1.0)

    def solve(lam: float):
        A = np.vstack((K, np.sqrt(lam) * D2))
        b = np.concatenate((g, np.zeros(n_gamma - 2)))
        c, _ = nnls(A, b)
        return c

    if lambda_reg is None:
        # L-curve: max curvature of (log resid, log solution-norm)
        lams = np.logspace(-8, 2, 36)
        pts, sols = [], []
        for lam in lams:
            c = solve(lam)
            res = float(np.linalg.norm(K @ c - g))
            sm = float(np.linalg.norm(D2 @ c))
            if res <= 0 or sm <= 0:
                continue
            pts.append((np.log(res), np.log(sm)))
            sols.append(c)
        if len(pts) < 8:
            raise RuntimeError("L-curve sweep degenerate -- check input data")
        P = np.asarray(pts)
        # discrete curvature of parametric curve
        curv = []
        for i in range(1, len(P) - 1):
            p0, p1, p2 = P[i - 1], P[i], P[i + 1]
            u, v = p1 - p0, p2 - p1
            cross2d = abs(u[0] * v[1] - u[1] * v[0])  # 2-D cross by hand
            den = (np.linalg.norm(u) * np.linalg.norm(v)) ** 1.5
            curv.append(cross2d / den if den > 0 else 0.0)
        k = int(np.argmax(curv)) + 1
        lam_used = float(lams[k]) if k < len(lams) else float(lams[len(curv)])
        c = solve(lam_used)
    else:
        if lambda_reg <= 0:
            raise ValueError("lambda_reg must be positive when given")
        lam_used = float(lambda_reg)
        c = solve(lam_used)

    # intensity distribution vs diameter (Stokes-Einstein)
    kT = 1.380649e-23 * temperature_K
    eta = viscosity_mPa_s * 1e-3
    # D (m^2/s) for each Gamma: Gamma = D q^2 with q in nm^-1 -> q^2 in nm^-2
    D_nm2_s = gamma / (q_nm_inv ** 2)
    D_si = D_nm2_s * 1e-18
    d_m = kT / (3.0 * np.pi * eta * D_si)
    diam = d_m * 1e9  # nm

    weights = c / c.sum() if c.sum() > 0 else c
    # sort by diameter ascending
    idx = np.argsort(diam)
    diam = diam[idx]
    weights = weights[idx]

    dh_mean = float(np.sum(weights * diam))
    #                 intensity distribution (weights already sum to 1)
    var = float(np.sum(weights * (diam - dh_mean) ** 2))
    pdi = var / (dh_mean ** 2) if dh_mean > 0 else float("nan")

    rel_res = float(np.linalg.norm(K @ c - g) / np.linalg.norm(g))

    return DLSInversionResult(
        diameter_nm=tuple(float(v) for v in diam),
        intensity_fraction=tuple(float(v) for v in weights),
        dh_intensity_mean_nm=dh_mean,
        dh_z_avg_nm=dh_mean,
        pdi_from_distribution=float(pdi),
        q_nm_inv=float(q_nm_inv),
        lambda_reg=lam_used,
        residual_norm=float(np.linalg.norm(K @ c - g)),
        relative_residual=rel_res,
        n_grid=int(n_gamma),
        converged=bool(rel_res < 0.15),
        method="non-negative least squares with second-difference Tikhonov "
               "smoothing; lambda from L-curve corner; Stokes-Einstein sizing",
        params={"wavelength_nm": wavelength_nm, "angle_deg": angle_deg,
                "refractive_index": refractive_index,
                "viscosity_mPa_s": viscosity_mPa_s,
                "temperature_K": temperature_K, "n_gamma": n_gamma,
                "lambda_reg_requested": lambda_reg},
        caveat=(
            "ILL-POSED PROBLEM: resolution is limited by the tau range and "
            "noise; narrow peaks are broaderened and sizes outside the "
            "Gamma-grid range are invisible. Distribution is intensity-"
            "weighted (large particles dominate); do not read it as a number "
            "distribution without Mie/Rayleigh conversion. PDI here is the "
            "relative variance of the recovered distribution, not the "
            "ISO 22412 cumulants PDI."
        ),
    )

if __name__ == "__main__":
    # Real pilot dataset: AOT in water, reconstructed from Shah/Das/Bhattarai
    # 2025 (Heliyon, PMC11835642) Figure 1 -- 6 premicellar points computed
    # from the paper's own published regression line (y = -16.4x - 16.5,
    # r^2=0.994, x = log10[AOT mol/L]), 3 post-CMC plateau points estimated
    # by eye from the same figure's flat region. Paper's own reported CMC
    # for AOT in water: 2.51 mM.
    concentrations = [0.01585, 0.03981, 0.10000, 0.25119, 0.63096, 1.58489,
                       3.981, 6.31, 10.0]
    tensions = [62.22, 55.66, 49.10, 42.54, 35.98, 29.42,
                27.5, 27.8, 28.0]
    result = cmc_from_surface_tension_curve(concentrations, tensions)
    print(f"Extracted CMC: {result.cmc_mM:.3f} mM  (paper's reported value: 2.51 mM)")
    print(f"Premicellar slope: {result.premicellar_slope_mN_per_m_per_lnC:.3f} mN/m per ln(C)")
    print(f"Gamma at CMC: {result.gamma_at_cmc_mN_per_m:.2f} mN/m")
    print(f"R^2 (premicellar fit): {result.r_squared_premicellar:.4f}")
    print(f"Split: {result.n_premicellar_points} premicellar / {result.n_postmicellar_points} post-CMC points")
