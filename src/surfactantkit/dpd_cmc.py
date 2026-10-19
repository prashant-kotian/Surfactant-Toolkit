"""S4 analysis core: Anderson short-time CMC prediction from p_free traces.

Implements the analysis pipeline of Rummel, Robinson, Warren, Bray &
Anderson, J. Phys. Chem. B 2026 (doi 10.1021/acs.jpcb.5c07475, CC-BY,
PMC12794194), the part that is interactive-Python feasible without the
dl_meso production DPD engine:

* eq 3 (nonionic): C_free ~= CMC above the CMC (plateau reading),
* eq 4 (ionic, charge-phase): log C_free + beta*log(beta*C_free +
  (1-beta)*C_tot/(1-phi)) ~= (1+beta)*log CMC, with phi = V_m*C_tot the
  surfactant packing fraction (paper: V_m at assumed density 1 g/cm3);
  beta=0 reduces exactly to eq 3,
* eq 8: p_free(t) = A + (B-A) exp(-t/tau) exponential fit of free-
  monomer ratio traces (Aniansson-Wall/Becker-Doring motivated),
* jackknife resampling with j_ratio (0.9 = leave-one-out over the repeat
  traces, exhaustively enumerated as in the paper) for the variance of A
  and of the derived CMC.

Everything is pure stdlib + the repo's own golden-section minimizer --
no numpy, no scipy, no RNG (fully deterministic for fixed inputs).
The production DPD simulation worker (`dpd_micelle_cmc`) is a separate,
batch-class capability; this module only ANALYZES traces you already have
(yours or from any simulator).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from itertools import combinations, pairwise

from ._optimize import golden_section_minimize


@dataclass
class AndersonCMCResult:
    """Short-time p_free analysis output (blueprint schema + values)."""
    cmc_M: float | None
    cmc_jackknife_mean_M: float | None
    cmc_jackknife_sigma_M: float | None
    p_free_inf: float | None          # eq 8 A (full data set)
    p_free_0: float | None            # eq 8 B (full data set)
    tau_s: float | None               # eq 8 tau (full data set)
    p_free_inf_jackknife_sigma: float | None
    n_subsets: int
    subset_cmc_M: list[float]
    converged: bool
    n_iterations: int
    method: str
    params: dict = field(default_factory=dict)
    caveat: str = ""

    def as_dict(self) -> dict:
        return {k: getattr(self, k) for k in
                ("cmc_M", "cmc_jackknife_mean_M", "cmc_jackknife_sigma_M",
                 "p_free_inf", "p_free_0", "tau_s",
                 "p_free_inf_jackknife_sigma", "n_subsets", "subset_cmc_M",
                 "converged", "n_iterations", "method", "params", "caveat")}


def _phi_from_v_m(c_tot: float, v_m_cm3_per_mol: float | None) -> float:
    """phi = V_m * C_tot with V_m in L/mol (v_m cm3/mol / 1000), C_tot in M."""
    if v_m_cm3_per_mol is None:
        return 0.0
    return v_m_cm3_per_mol * c_tot / 1000.0


def _eq4_lhs_minus_rhs(c_free: float, c_tot: float, beta: float,
                       phi: float) -> float:
    """log C_free + beta*log(beta*C_free + (1-beta)*C_tot/(1-phi))
    - (1+beta)*log CMC  evaluated at CMC := cmc_from_free_monomer; i.e. the
    eq-4 residual (0 at the solution). Base of log is irrelevant (same base
    on both sides)."""
    crowding = (1.0 - phi) if phi < 1.0 else 1e-300
    inner = beta * c_free + (1.0 - beta) * c_tot / crowding
    log_cmc = (math.log(c_free) + beta * math.log(inner)) / (1.0 + beta)
    return math.log(c_free) + beta * math.log(inner) - (1.0 + beta) * log_cmc


def cmc_from_free_monomer(c_free: float, c_tot: float, beta: float = 0.0,
                          v_m_cm3_per_mol: float | None = None) -> float:
    """Eq 3/eq 4 forward: CMC from a measured free-monomer concentration.

    Parameters
    ----------
    c_free : free (unassociated) surfactant concentration, mol/L
    c_tot  : total surfactant concentration in the same box, mol/L
    beta   : counterion-binding degree (0 = nonionic -> exact eq 3)
    v_m_cm3_per_mol : surfactant molar volume, cm3/mol (the paper assumes
        density 1 g/cm3 when deriving V_m). None disables the 1/(1-phi)
        crowding factor (phi = 0).

    Returns CMC in mol/L. All concentrations must be > 0, 0 <= beta < 1,
    0 < phi < 1 -- ValueError names the offending field.
    """
    if not (c_free > 0.0):
        raise ValueError(f"c_free must be > 0 (got {c_free!r}) -- log-space eq 4 "
                         "is undefined at or below zero")
    if not (c_tot > 0.0):
        raise ValueError(f"c_tot must be > 0 (got {c_tot!r})")
    if not (0.0 <= beta < 1.0):
        raise ValueError(f"beta must be in [0, 1) (got {beta!r})")
    phi = _phi_from_v_m(c_tot, v_m_cm3_per_mol)
    if v_m_cm3_per_mol is not None and not (0.0 < v_m_cm3_per_mol < 5000.0):
        raise ValueError(f"v_m_cm3_per_mol must be in (0, 5000) cm3/mol "
                         f"(got {v_m_cm3_per_mol!r})")
    if not (0.0 <= phi < 0.95):
        raise ValueError(f"packing fraction phi={phi:.4f} from "
                         f"v_m_cm3_per_mol={v_m_cm3_per_mol!r} x c_tot="
                         f"{c_tot!r} must be in [0, 0.95) -- the 1/(1-phi) "
                         "crowding factor diverges as phi -> 1")
    if beta == 0.0:
        # exact eq 3: CMC = C_free, independent of phi (bit-identical, no
        # exp(log()) round-trip error)
        return c_free
    crowding = 1.0 - phi
    inner = beta * c_free + (1.0 - beta) * c_tot / crowding
    log_cmc = (math.log(c_free) + beta * math.log(inner)) / (1.0 + beta)
    return math.exp(log_cmc)


def free_monomer_from_cmc(cmc: float, c_tot: float, beta: float = 0.0,
                          v_m_cm3_per_mol: float | None = None) -> float:
    """Eq 3/eq 4 inverse: free-monomer concentration implied by a CMC.

    Solves the monotonically increasing eq-4 residual for C_free in
    (0, c_tot] by bisection (300 iterations, machine-precision stop).
    Same field guards as cmc_from_free_monomer, plus 0 < cmc <= c_tot.
    """
    if not (cmc > 0.0):
        raise ValueError(f"cmc must be > 0 (got {cmc!r})")
    if not (c_tot > 0.0):
        raise ValueError(f"c_tot must be > 0 (got {c_tot!r})")
    if cmc > c_tot:
        raise ValueError(f"cmc ({cmc!r}) cannot exceed c_tot ({c_tot!r})")
    # validate beta/v_m via the forward function at c_free = c_tot
    cmc_from_free_monomer(c_tot, c_tot, beta=beta, v_m_cm3_per_mol=v_m_cm3_per_mol)
    phi = _phi_from_v_m(c_tot, v_m_cm3_per_mol)

    def residual(x: float) -> float:
        # 0 at the solution; increasing in x; equals (1+beta)*(log x - log cmc)
        # at beta=0 (exact eq 3 -> x = cmc)
        crowding = 1.0 - phi
        inner = beta * x + (1.0 - beta) * c_tot / crowding
        return math.log(x) + beta * math.log(inner) - (1.0 + beta) * math.log(cmc)

    lo, hi = 1e-300, c_tot
    if residual(hi) < 0.0:  # pragma: no cover -- guarded by cmc <= c_tot
        raise ValueError(f"cmc ({cmc!r}) is above the eq-4 maximum for "
                         f"c_tot ({c_tot!r})")
    for _ in range(300):
        mid = math.sqrt(lo * hi)          # geometric bisection (positive domain)
        if residual(mid) < 0.0:
            lo = mid
        else:
            hi = mid
        if hi / lo - 1.0 < 1e-14:
            break
    return math.sqrt(lo * hi)


def _linfit_ab(times: list[float], ys: list[float], tau: float) -> tuple[float, float]:
    """Closed-form OLS for y = A*(1-E) + B*E with E = exp(-t/tau)."""
    e = [math.exp(-t / tau) for t in times]
    u = [1.0 - ei for ei in e]            # coefficient of A
    suu = sum(ui * ui for ui in u)
    suv = sum(ui * ei for ui, ei in zip(u, e))
    svv = sum(ei * ei for ei in e)
    suy = sum(ui * yi for ui, yi in zip(u, ys))
    svy = sum(ei * yi for ei, yi in zip(e, ys))
    det = suu * svv - suv * suv
    if abs(det) < 1e-300:
        # degenerate (tau far from any scale): A/B unidentifiable here
        return float(sum(ys) / len(ys)), float(sum(ys) / len(ys))
    a = (suy * svv - svy * suv) / det
    b = (svy * suu - suy * suv) / det
    return a, b


def _fit_pfree_exp(times: list[float], ys: list[float]) -> tuple[float, float, float, float, float, bool, bool]:
    """Eq-8 fit p = A + (B-A) exp(-t/tau) by tau-profile + golden section.

    Returns (A, B, tau, sse, sse_over_sst, ok, tau_edge_pinned).
    """
    n = len(times)
    ybar = sum(ys) / n
    sst = sum((y - ybar) ** 2 for y in ys)
    if sst <= 1e-30:
        # flat trace: tau unidentifiable; A = mean, converged False
        return ybar, ybar, float("nan"), 0.0, 0.0, False, False

    def sse_at(tau: float) -> float:
        a, b = _linfit_ab(times, ys, tau)
        e = [math.exp(-t / tau) for t in times]
        return sum((a * (1.0 - ei) + b * ei - y) ** 2
                   for ei, y in zip(e, ys))

    span = times[-1] - times[0]
    dt_min = min(b - a for a, b in pairwise(times))
    lo_bound = max(dt_min / 100.0, span * 1e-6)
    hi_bound = span * 100.0

    # coarse log grid, then golden-section refine in the best bracket
    n_grid = 64
    logs = [math.log(lo_bound) + i * (math.log(hi_bound) - math.log(lo_bound))
            / (n_grid - 1) for i in range(n_grid)]
    sses = [sse_at(math.exp(lt)) for lt in logs]
    k = min(range(n_grid), key=lambda i: sses[i])
    g_lo = math.exp(logs[max(k - 1, 0)])
    g_hi = math.exp(logs[min(k + 1, n_grid - 1)])
    tau = golden_section_minimize(sse_at, g_lo, g_hi, iters=200, tol=1e-14)
    a, b = _linfit_ab(times, ys, tau)
    sse = sse_at(tau)
    ok = sse <= sst + 1e-30 and math.isfinite(sse)
    edge = (k in (0, n_grid - 1))
    return a, b, tau, sse, sse / sst, ok, edge


def _jackknife_sigma(estimates: list[float]) -> float:
    """Standard jackknife sigma: sqrt((m-1)/m * sum (t_i - mean)^2)."""
    m = len(estimates)
    if m < 2:
        return 0.0
    mean = sum(estimates) / m
    return math.sqrt((m - 1) / m * sum((t - mean) ** 2 for t in estimates))


def _comb_count(n: int, k: int) -> int:
    return math.comb(n, k)


def anderson_cmc_analysis(
    times: list[float],
    traces: list[list[float]],
    c_tot: float,
    beta: float = 0.0,
    v_m_cm3_per_mol: float | None = None,
    j_ratio: float = 0.9,
    max_subsets: int = 500,
) -> AndersonCMCResult:
    """S4 analysis core: eq-3/eq-4 CMC from p_free(t) repeat traces.

    Parameters
    ----------
    times  : strictly increasing time points, same grid for all traces (s)
    traces : repeat traces of p_free = free/total surfactant RATIO
             (>= 3 repeats; the paper uses 10 x 20-40k DPD time units)
    c_tot  : total surfactant concentration of the runs, mol/L
    beta   : counterion binding (0 -> eq 3 nonionic; ionic surfactants
             typically 0.4-0.8 -- pass a sourced value, never guess)
    v_m_cm3_per_mol : molar volume for the 1/(1-phi) crowding factor
             (paper assumes density 1 g/cm3); None -> phi = 0 (disclosed)
    j_ratio: jackknife ratio (paper: 0.9 leave-one-out; 0.7/0.8 also
             reported in the paper as qualitatively similar)
    max_subsets: hard cap on exhaustively enumerated subsets (paper
             enumerates; we refuse rather than silently sample)

    Returns AndersonCMCResult (blueprint schema: converged,
    n_iterations, method, params, caveat).
    """
    # ---- input guards (fail fast, field names) ----
    if not isinstance(traces, list) or len(traces) < 3:
        raise ValueError(f"traces must be a list of >= 3 repeat traces "
                         f"(got {len(traces) if isinstance(traces, list) else type(traces).__name__}) "
                         "-- jackknife resampling needs repeats (the paper uses 10)")
    if not all(isinstance(t, (list, tuple)) for t in traces):
        raise ValueError("traces must be a list of numeric lists (one per repeat)")
    n_pts = len(times)
    if n_pts < 4:
        raise ValueError(f"times needs >= 4 points for an eq-8 fit (got {n_pts})")
    t = [float(x) for x in times]
    if any(b <= a for a, b in pairwise(t)):
        raise ValueError("times must be strictly increasing")
    yss: list[list[float]] = []
    for i, tr in enumerate(traces):
        if len(tr) != n_pts:
            raise ValueError(f"traces[{i}] has {len(tr)} points but times has {n_pts}")
        y = [float(x) for x in tr]
        if any(not (0.0 <= x <= 1.0 + 1e-9) for x in y):
            raise ValueError(f"traces[{i}] contains values outside [0, 1] -- "
                             "p_free must be a free/total RATIO, not a "
                             "concentration (0 is a valid finite-box reading "
                             "meaning every molecule sits in an aggregate "
                             "larger than n*)")
        yss.append(y)
    if not (c_tot > 0.0):
        raise ValueError(f"c_tot must be > 0 (got {c_tot!r})")
    if not (0.0 <= beta < 1.0):
        raise ValueError(f"beta must be in [0, 1) (got {beta!r})")
    if not (0.5 <= j_ratio < 1.0):
        raise ValueError(f"j_ratio must be in [0.5, 1) (got {j_ratio!r})")

    n_rep = len(yss)
    k_sub = math.floor(j_ratio * n_rep + 1e-9)
    if k_sub < 2:
        raise ValueError(f"j_ratio={j_ratio} with n_traces={n_rep} leaves "
                         f"{k_sub} traces per subset -- need >= 2; raise "
                         "j_ratio or add repeats")
    if k_sub >= n_rep:
        k_sub = n_rep - 1
    n_subsets = _comb_count(n_rep, k_sub)
    if n_subsets > max_subsets:
        raise ValueError(f"j_ratio={j_ratio} with n_traces={n_rep} enumerates "
                         f"{n_subsets} subsets > max_subsets={max_subsets} -- "
                         "the paper enumerates exhaustively; we refuse rather "
                         "than silently subsample (raise max_subsets or use "
                         "j_ratio closer to 1)")

    # phi guard via the forward equation (validates v_m too)
    phi = _phi_from_v_m(c_tot, v_m_cm3_per_mol)
    cmc_from_free_monomer(c_tot, c_tot, beta=beta, v_m_cm3_per_mol=v_m_cm3_per_mol)

    def mean_trace(members: list[int]) -> list[float]:
        n = len(members)
        return [sum(yss[m][j] for m in members) / n for j in range(n_pts)]

    def fit_and_cmc(members: list[int]) -> tuple[float | None, float | None, tuple | None]:
        a, b, tau, sse, ratio, ok, edge = _fit_pfree_exp(t, mean_trace(members))
        if not ok or a <= 0.0 or a > 1.0 + 1e-9:
            return None, None, (a, b, tau, sse, ratio, ok, edge)
        c_free = a * c_tot
        cmc = cmc_from_free_monomer(c_free, c_tot, beta=beta,
                                    v_m_cm3_per_mol=v_m_cm3_per_mol)
        return cmc, c_free, (a, b, tau, sse, ratio, ok, edge)

    # full data set (paper: average over the complete data set)
    cmc_full, _cf, fit_full = fit_and_cmc(list(range(n_rep)))

    # jackknife subsets (exhaustive combinations of size k_sub)
    subset_cmc: list[float] = []
    subset_a: list[float] = []
    all_ok = fit_full is not None and fit_full[5]
    for combo in combinations(range(n_rep), k_sub):
        cmc_i, _cf_i, fit_i = fit_and_cmc(list(combo))
        all_ok = all_ok and fit_i is not None and fit_i[5] and cmc_i is not None
        if cmc_i is not None and fit_i is not None and fit_i[0] > 0.0:
            subset_cmc.append(cmc_i)
            subset_a.append(fit_i[0])
    n_iter = 1 + n_subsets   # one eq-8 fit per subset plus the full set

    sigma_cmc = (_jackknife_sigma(subset_cmc)
                 if all_ok and len(subset_cmc) == n_subsets else None)
    sigma_a = (_jackknife_sigma(subset_a)
               if all_ok and len(subset_a) == n_subsets else None)
    mean_cmc = (sum(subset_cmc) / len(subset_cmc)
                if subset_cmc and len(subset_cmc) == n_subsets else None)

    converged = bool(all_ok and cmc_full is not None)
    a_full = fit_full[0] if fit_full else None
    b_full = fit_full[1] if fit_full else None
    tau_full = fit_full[2] if fit_full else None
    edge = bool(fit_full[6]) if fit_full else True

    caveat = (
        "Anderson short-time analysis (Rummel et al., J. Phys. Chem. B 2026, "
        "doi:10.1021/acs.jpcb.5c07475, CC-BY): eq-8 p_free(t)=A+(B-A)exp(-t/tau) "
        "fit to pointwise-mean traces, C_free = A*c_tot, CMC via eq 4 (beta>0) "
        "or exact eq 3 (beta=0); jackknife over exhaustive C(n, k) subsets at "
        "j_ratio (paper j=0.9 leave-one-out). Inputs are YOUR traces (this tool "
        "does not run DPD). phi="
        + (f"{phi:.4f} from v_m={v_m_cm3_per_mol} cm3/mol (paper assumes "
           "density 1 g/cm3)" if v_m_cm3_per_mol is not None else
           "0 -- v_m not given, crowding 1/(1-phi) term DISABLED")
        + "; beta is a user input (use a sourced value, never guess one). "
        "A tau pinned at the search bound (params.tau_edge_pinned) or "
        "converged=False means the time window does not resolve the slow mode "
        "-- lengthen the traces before quoting a CMC. No RNG: deterministic."
    )
    return AndersonCMCResult(
        cmc_M=cmc_full,
        cmc_jackknife_mean_M=mean_cmc,
        cmc_jackknife_sigma_M=sigma_cmc,
        p_free_inf=a_full,
        p_free_0=b_full,
        tau_s=(tau_full if tau_full is not None and math.isfinite(tau_full) else None),
        p_free_inf_jackknife_sigma=sigma_a,
        n_subsets=n_subsets,
        subset_cmc_M=[float(x) for x in subset_cmc],
        converged=converged,
        n_iterations=n_iter,
        method="anderson_pfree_exp_profile_fit_jackknife_eq4",
        params={
            "c_tot_M": float(c_tot),
            "beta": float(beta),
            "v_m_cm3_per_mol": (float(v_m_cm3_per_mol)
                                if v_m_cm3_per_mol is not None else None),
            "phi": float(phi),
            "j_ratio": float(j_ratio),
            "n_traces": int(n_rep),
            "n_time_points": int(n_pts),
            "subset_size": int(k_sub),
            "sse": (float(fit_full[3]) if fit_full else None),
            "sse_over_sst": (float(fit_full[4]) if fit_full else None),
            "tau_edge_pinned": edge,
            "time_span_s": float(t[-1] - t[0]),
        },
        caveat=caveat,
    )
