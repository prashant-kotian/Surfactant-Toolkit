"""Nagarajan molecular-thermodynamic CMC predictor (computational worker C1).

Implements, from primary source Nagarajan, R. "Molecular Theory for Mixed
Micelles", Langmuir 1985, 1, 331-341 (eqs 5-7, Table I, Appendix -- the
single-component free energy recast from Nagarajan & Ruckenstein, JCIS 1976
and 1979), the complete decomposition of the standard free energy of
micelle formation:

    (mu_g - g mu_1)/kT = g[transfer] + g[ordering]
                         + g*sigma*(a - a0)/kT          (core-water interface)
                         - g*ln(1 - a_p/a)              (steric heads)
                         + [electrostatics]             (charged heads only)

with parameters fixed by the source (sigma = 50 dyne/cm, beta = 0.46,
epsilon = 80, l0 = 1.50 + 1.26 nc A, v0 = 27 A3/CH2 + 54 A3/CH3, Table I
headgroup constants) and the appendix geometry (spheres up to r0 = l0,
then cylinders with hemispherical ends; the appendix's three
electrostatic cases, including the Bessel-function infinite-cylinder
expression).

CMC condition (eq 12 of the 2003 restatement): X_CMC = exp(min_g f(g))
with f(g) = per-molecule free energy, minimized over sphere and cylinder
branches; monomer concentration from C = X * 55.55 M feeds back into the
Debye length (iterated to self-consistency).

A priori: no fitted parameters beyond the constants published in the
source.  Validated in tests against seven experimental CMCs (sodium
C8-C14 alkyl sulfates, C10PO, C10SO, C12E8) to within 0.15 log10 units.

Theory caveats: point-ion electrostatics scaled by beta; eqs 6-7
calibrated at 25 C; no added-electrolyte-specific terms beyond kappa;
g_opt is the CMC-relevant free-energy minimum, NOT a prediction of the
equilibrium aggregate size.
"""
from __future__ import annotations

import math
from dataclasses import asdict, dataclass

# --- physical constants (SI unless noted) -------------------------------
K_B = 1.380649e-23          # J/K
E_CHARGE = 1.602176634e-19  # C
EPS0 = 8.8541878128e-12     # F/m
N_A = 6.02214076e23
EPS_R = 80.0                # Nagarajan's water dielectric constant
SIGMA_INTERFACE = 0.050     # 50 dyne/cm = 50 mN/m core-water interfacial tension
BETA_ELEC = 0.46            # Debye-Huckel correction constant (paper: 0.46)
WATER_MOLAR = 55.55         # mol/L water at 25 C (source: C = 55.55 X)

# --- Table I of Langmuir 1985, 1, 331-341 -------------------------------
# headgroup: (a_p [A2], a_i counterion radius [A], delta charge plane [A], ionic)
HEADGROUPS: dict[str, tuple[float, float, float, bool]] = {
    "alkyl_sulfate": (17.0, 1.0, 2.3, True),        # SOS/SDeS/SDS/STS
    "potassium_carboxylate": (11.0, 1.35, 1.0, True),  # KC8 ... KC14
    "dialkyl_phosphine_oxide": (44.0, 0.0, 0.0, False),  # C10PO (nonionic)
    "alkyl_methyl_sulfoxide": (37.0, 0.0, 0.0, False),   # C10SO (nonionic)
    "polyoxyethylene": (36.0, 0.0, 0.0, False),      # C12E8 (nonionic)
}


def _tail_v0(n_carbons: float) -> float:
    """Tail molecular volume (m^3): 27 A3 per CH2 + 54 A3 for CH3."""
    return (27.0 * (n_carbons - 1.0) + 54.0) * 1e-30


def _tail_l0(n_carbons: float) -> float:
    """Extended tail length (m): l0 = 1.50 + 1.26 nc (Angstrom)."""
    return (1.50 + 1.26 * n_carbons) * 1e-10


def _kappa_per_m(conc_mol_L: float, temperature_K: float) -> float:
    """Debye reciprocal length (1/m) for 1:1 electrolyte, general SI.

    Agrees with the source's kappa = sqrt(C)/(3.08e-8 cm^-1) at 25 C to ~1%.
    """
    if conc_mol_L <= 0:
        return 1e-6  # effectively field-free; never called with 0 by the CMC loop
    i_m3 = conc_mol_L * 1000.0  # mol/m^3
    kappa_sq = (2.0 * N_A * E_CHARGE ** 2 * i_m3 /
                (EPS0 * EPS_R * K_B * temperature_K))
    return math.sqrt(kappa_sq)


def _sphere_geometry(g: int, n_carbons: float):
    """(r0, a) for a spherical micelle of g molecules; None if r0 > l0."""
    v0, l0 = _tail_v0(n_carbons), _tail_l0(n_carbons)
    r0 = (3.0 * g * v0 / (4.0 * math.pi)) ** (1.0 / 3.0)
    if r0 > l0:
        return None
    a = 4.0 * math.pi * r0 * r0 / g
    return r0, a


def _cylinder_geometry(g: int, n_carbons: float):
    """(r0, L, a) for a cylinder with hemispherical ends; None if invalid."""
    v0, l0 = _tail_v0(n_carbons), _tail_l0(n_carbons)
    r0 = l0
    L = (g * v0 - 4.0 * math.pi * r0 ** 3 / 3.0) / (math.pi * r0 ** 2)
    if L <= 0:
        return None
    a = (4.0 * math.pi * r0 ** 2 + 2.0 * math.pi * r0 * L) / g
    return r0, L, a


@dataclass(frozen=True)
class FreeEnergyTerms:
    """eq 5 decomposition for one aggregate size (units: kT)."""
    g: int
    branch: str
    transfer: float
    ordering: float
    interface: float
    steric: float
    electrostatic: float
    total: float
    per_molecule: float
    core_radius_nm: float | None
    area_per_molecule_A2: float | None
    cylinder_length_nm: float | None
    method: str
    caveat: str

    def as_dict(self) -> dict:
        return asdict(self)


def _terms_at_g(g, n_carbons, a_p, a_i, delta, ionic,
                kappa, temperature_K, branch) -> FreeEnergyTerms | None:
    """Evaluate all six eq5 terms for one g on one branch (kT units)."""
    kT = K_B * temperature_K
    if branch == "sphere":
        geo = _sphere_geometry(g, n_carbons)
        if geo is None:
            return None
        r0, a = geo
        L = None
    else:
        geo = _cylinder_geometry(g, n_carbons)
        if geo is None:
            return None
        r0, L, a = geo
    a_p_m2 = a_p * 1e-20
    if a <= a_p_m2 * 1.0001:
        return None  # steric log undefined: packing below head area

    a0 = min(a_p, 21.0) * 1e-20  # shielded core area: min(a_p, 21 A2)
    t_tr = -g * (2.05 + 1.49 * n_carbons)                     # eq 6
    t_or = g * (-0.5 + 0.24 * n_carbons)                      # eq 7
    t_int = g * SIGMA_INTERFACE * (a - a0) / kT
    t_ste = -g * math.log(1.0 - a_p_m2 / a)

    t_el = 0.0
    if ionic:
        q2_over_4pe0 = E_CHARGE ** 2 / (4.0 * math.pi * EPS0 * EPS_R)
        if branch == "sphere":
            r = r0 + delta * 1e-10
            pref = g * g * BETA_ELEC * q2_over_4pe0 / (2.0 * kT * r)
            ka_i = kappa * a_i * 1e-10
            brac = (1.0 + ka_i) / (1.0 + ka_i + kappa * r)
            t_el = pref * brac
        else:
            r = r0 + delta * 1e-10
            if L / (2.0 * r0) > 3.0:
                # appendix rule (b): sphere for hemispherical ends,
                # infinite cylinder for the body
                from scipy.special import k0, k1
                g_ends = (4.0 * math.pi * r0 ** 3 / 3.0) / _tail_v0(n_carbons)
                g_body = g - g_ends
                ka_i = kappa * a_i * 1e-10
                pref_s = g_ends * g_ends * BETA_ELEC * q2_over_4pe0 / (2.0 * kT * r)
                e_s = pref_s * (1.0 + ka_i) / (1.0 + ka_i + kappa * r)
                kr = kappa * r
                pref_c = g_body * g_body * BETA_ELEC * q2_over_4pe0 / (kT * L)
                e_c = pref_c * (k0(kr) / (kr * k1(kr))
                                + math.log(1.0 + a_i * 1e-10 / r))
                t_el = e_s + e_c
            else:
                # appendix rule (c): equivalent sphere of equal area
                from scipy.special import k0
                r_eqv = math.sqrt(r * r + r * L / 2.0)
                pref = g * g * BETA_ELEC * q2_over_4pe0 / (2.0 * kT * r_eqv)
                ka_i = kappa * a_i * 1e-10
                t_el = pref * (1.0 + ka_i) / (1.0 + ka_i + kappa * r_eqv)

    total = t_tr + t_or + t_int + t_ste + t_el
    return FreeEnergyTerms(
        g=g,
        branch=branch,
        transfer=t_tr,
        ordering=t_or,
        interface=t_int,
        steric=t_ste,
        electrostatic=t_el,
        total=total,
        per_molecule=total / g,
        # sphere branch and cylinder-body branch share the same core radius
        # r0 by construction (hemispherical ends on the cylinder): not a bug,
        # but the old identical-branch ternary (RUF034) obscured that.
        core_radius_nm=r0 * 1e9,
        area_per_molecule_A2=a * 1e20,
        cylinder_length_nm=(L * 1e9 if L is not None else None),
        method="Nagarajan Langmuir 1985 eq 5 (terms in kT)",
        caveat="Terms are the standard free-energy decomposition for a "
               "hypothetical aggregate of size g; only their minimum over g "
               "enters the CMC.",
    )


def nagarajan_free_energy_terms(
    g: int,
    n_carbons: int = 12,
    headgroup: str = "alkyl_sulfate",
    *,
    concentration_mol_per_L: float = 1e-3,
    temperature_K: float = 298.15,
    branch: str = "sphere",
    a_p: float | None = None,
    a_i: float | None = None,
    delta: float | None = None,
    ionic: bool | None = None,
    added_salt_M: float = 0.0,
) -> FreeEnergyTerms:
    """Full eq 5 free-energy term table for one aggregation number g.

    Useful for inspecting which physical term dominates (transfer vs
    interface vs steric vs electrostatics) at a given aggregate size.
    """
    if g < 2:
        raise ValueError("g must be >= 2 (dimers are the smallest aggregates)")
    if n_carbons < 4:
        raise ValueError("n_carbons must be >= 4")
    if branch not in ("sphere", "cylinder"):
        raise ValueError("branch must be 'sphere' or 'cylinder'")
    ap, ai, dl, io = _resolve_headgroup(headgroup, a_p, a_i, delta, ionic)
    kappa = _kappa_per_m(concentration_mol_per_L + added_salt_M, temperature_K)
    res = _terms_at_g(g, n_carbons, ap, ai, dl, io, kappa, temperature_K, branch)
    if res is None:
        raise ValueError(
            f"g={g} is not realizable as {branch} for C{n_carbons} with "
            f"a_p={ap} A2 (check sphere limit r0 <= l0 or packing a > a_p)"
        )
    return res


def _resolve_headgroup(headgroup, a_p, a_i, delta, ionic):
    if a_p is None:
        if headgroup not in HEADGROUPS:
            raise ValueError(
                f"unknown headgroup {headgroup!r}; valid: "
                f"{sorted(HEADGROUPS)} -- or pass a_p/a_i/delta/ionic "
                "explicitly for a custom headgroup"
            )
        ap, ai, dl, io = HEADGROUPS[headgroup]
    else:
        if None in (a_i, delta, ionic):
            raise ValueError(
                "custom headgroups need all of a_p, a_i, delta, ionic"
            )
        ap, ai, dl, io = float(a_p), float(a_i), float(delta), bool(ionic)
    if ap <= 0 or dl < 0 or ai < 0:
        raise ValueError("a_p must be > 0; a_i, delta must be >= 0")
    return ap, ai, dl, io


@dataclass(frozen=True)
class NagarajanCMCResult:
    """A-priori CMC prediction from the Nagarajan free-energy minimization."""
    n_carbons: int
    headgroup: str
    cmc_mol_per_L: float
    cmc_mole_fraction: float
    log10_cmc_mol_per_L: float
    delta_g_mic_kT_per_molecule: float
    delta_g_mic_kJ_per_mol: float
    terms_at_optimum: dict
    g_optimum: int
    branch_optimum: str
    largest_sphere_g: int
    kappa_nm_inv: float
    monomer_conc_iterated_M: float
    n_iterations: int
    converged: bool
    method: str
    params: dict
    caveat: str

    def as_dict(self) -> dict:
        return asdict(self)


def nagarajan_cmc_predictor(
    n_carbons: int,
    headgroup: str = "alkyl_sulfate",
    *,
    temperature_K: float = 298.15,
    added_salt_M: float = 0.0,
    n_scan_max: int = 20000,
    a_p: float | None = None,
    a_i: float | None = None,
    delta: float | None = None,
    ionic: bool | None = None,
) -> NagarajanCMCResult:
    """Predict the critical micelle concentration from molecular structure.

    Minimizes the Nagarajan eq 5 standard free energy per molecule over
    sphere (r0 <= l0) and cylinder (hemispherical ends) aggregates, then
    applies the phase-separation CMC condition X_CMC = exp(min f) with
    self-consistent Debye screening (monomer concentration feeds kappa).
    Returns the CMC, the term-by-term free-energy breakdown at the
    optimum, and full provenance.

    No fitted parameters: sigma=50 dyne/cm, beta=0.46, Table I headgroup
    constants, eqs 6-7 exactly as published (25 C calibration -- see
    caveat). Requires scipy for charged cylinders (Bessel K0/K1;
    surfactantkit[calc]); nonionic surfactants run without it.
    """
    if not (4 <= n_carbons <= 30):
        raise ValueError("n_carbons must be in [4, 30] (model validity range)")
    if temperature_K <= 0:
        raise ValueError("temperature_K must be positive")
    if added_salt_M < 0:
        raise ValueError("added_salt_M must be >= 0")
    if n_scan_max < 100:
        raise ValueError("n_scan_max must be >= 100")
    ap, ai, dl, io = _resolve_headgroup(headgroup, a_p, a_i, delta, ionic)

    v0, l0 = _tail_v0(n_carbons), _tail_l0(n_carbons)
    g_s = int(4.0 * math.pi * l0 ** 3 / (3.0 * v0))

    def best_f(c1_guess_M: float):
        kappa = _kappa_per_m(c1_guess_M + added_salt_M, temperature_K)
        best = None
        for g in range(2, max(g_s, 3)):
            res = _terms_at_g(g, n_carbons, ap, ai, dl, io, kappa,
                              temperature_K, "sphere")
            if res is not None and (best is None or res.per_molecule < best.per_molecule):
                best = res
        for g in range(g_s + 1, n_scan_max + 1):
            res = _terms_at_g(g, n_carbons, ap, ai, dl, io, kappa,
                              temperature_K, "cylinder")
            if res is None:
                # packing fell below a_p: cylinder branch closed
                break
            if best is None or res.per_molecule < best.per_molecule:
                best = res
            elif g > g_s + 50 and res.per_molecule > best.per_molecule + 1e-6:
                break
        return best

    c1 = 1e-3  # initial monomer guess for kappa
    best = None
    n_iter = 0
    converged = False
    for n_iter in range(1, 21):
        best = best_f(c1)
        if best is None:
            raise ValueError(
                "Nagarajan model found no physically valid aggregate for "
                "these parameters (every trial aggregation number is steric-"
                "blocked or undefined): check headgroup area / tail volume")
        new_c1 = math.exp(best.per_molecule) * WATER_MOLAR
        if abs(math.log(new_c1 / c1)) <= 1e-9:
            c1 = new_c1
            converged = True
            break
        # log-space damping: the kappa-CMC fixed point can oscillate for
        # charged surfactants (higher CMC -> higher kappa -> weaker
        # electrostatics -> higher CMC ...); geometric averaging is stable
        c1 = new_c1 if n_iter == 1 else math.sqrt(c1 * new_c1)

    x_cmc = math.exp(best.per_molecule)
    cmc = x_cmc * WATER_MOLAR
    kappa = _kappa_per_m(c1 + added_salt_M, temperature_K)

    return NagarajanCMCResult(
        n_carbons=n_carbons,
        headgroup=headgroup,
        cmc_mol_per_L=cmc,
        cmc_mole_fraction=x_cmc,
        log10_cmc_mol_per_L=math.log10(cmc) if cmc > 0 else float("-inf"),
        delta_g_mic_kT_per_molecule=best.per_molecule,
        delta_g_mic_kJ_per_mol=best.per_molecule * K_B * temperature_K * N_A / 1000.0,
        terms_at_optimum=best.as_dict(),
        g_optimum=best.g,
        branch_optimum=best.branch,
        largest_sphere_g=g_s,
        kappa_nm_inv=kappa * 1e-9,
        monomer_conc_iterated_M=c1,
        n_iterations=n_iter,
        converged=converged,
        method="Nagarajan eq 5 free-energy minimization over sphere+cylinder "
               "aggregates; CMC = exp(min_g f) with self-consistent Debye "
               "screening (Langmuir 1985, 1, 331-341; appendix geometry)",
        params={
            "n_carbons": n_carbons, "headgroup": headgroup,
            "temperature_K": temperature_K, "added_salt_M": added_salt_M,
            "a_p_A2": ap, "a_i_A": ai, "delta_A": dl, "ionic": io,
            "sigma_mN_m": SIGMA_INTERFACE * 1000.0, "beta_elec": BETA_ELEC,
            "eps_r": EPS_R, "n_scan_max": n_scan_max,
        },
        caveat=(
            "A-priori molecular-thermodynamic estimate: point-ion "
            "electrostatics scaled by beta=0.46, eqs 6-7 calibrated at 25 C "
            "(temperature dependence of transfer/ordering terms is not "
            "modeled), no specific-ion chemistry. g_optimum is where the "
            "per-molecule free energy is minimized (it sets the CMC) -- it "
            "is NOT a prediction of the equilibrium micelle size."
        ),
    )
