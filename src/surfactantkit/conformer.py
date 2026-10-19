"""C5 -- 3D molecular volume/shape from ETKDG conformers + force-field
minimization: cross-check (and where branched/gemini tails break the
group contribution) replace `tanford_tail_volume`.

Pipeline (rdkit, deterministic):
  1. SMILES -> explicit-H mol, ETKDGv3 embedding of `n_conformers`
     conformers with a FIXED random seed (and 1 thread)
  2. MMFF94s minimization of every conformer (UFF fallback when MMFF
     parameters are missing; whichever ran is reported in params)
  3. Lowest-energy conformer selected -> outputs:
     * vdW envelope volume (rdkit DoubleCubicLatticeVolume, GetVDWVolume,
       probe independent) in A^3
     * Labute approximate SASA in A^2 (CalcLabuteASA)
     * relative shape anisotropy k^2 of the HEAVY-ATOM gyration tensor
       (k^2 = 1 - 3(l1l2+l2l3+l3l1)/(l1+l2+l3)^2; 0 = sphere-like,
       1 = rod-like; eigenvalues by hand-rolled Jacobi rotation, house
       no-scipy style)
     * optional ratio against a caller-supplied `tanford_tail_volume`
       value (the caller computes it from the existing worker so both
       numbers stay honest products of their own definitions)

Honesty [R]/[H]: V_3D is a WHOLE-MOLECULE property of one embedded
conformer; tanford_tail_volume is a TAIL-ONLY group contribution for
linear chains -- the ratio is a definitional cross-check, not a
validation, and branched/gemini tails are exactly where the two
definitions are expected to diverge (that divergence is the point of
the worker, not an error). ETKDG/MMFF are conformational models, not
free-energy calculations; lowest of N sampled conformers is reported.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class ConformerVolumeResult:
    """Result of `conformer_volume_3d`."""

    smiles: str
    vdW_volume_A3: float            # rdkit DCLV envelope volume
    sasa_A2: float                  # Labute approximate SASA
    shape_anisotropy_k2: float      # heavy-atom gyration anisotropy, [0, 1)
    lowest_energy_kcal_mol: float | None
    force_field: str                # "MMFF94s" | "UFF" | "none"
    n_conformers_embedded: int
    n_conformers_minimized: int
    tanford_ratio: float | None     # V_3D / caller's tanford tail volume
    converged: bool                 # embedding + minimization both succeeded
    n_iterations: int               # effort counter: 1/conformer if a single pass converged, 500 if the cap was hit
    method: str
    params: dict
    caveat: str


def _sym3_eigenvalues(a: list[list[float]]) -> list[float]:
    """Eigenvalues of a symmetric 3x3 matrix via cyclic Jacobi rotations
    (hand-rolled: house zero-dependency style). Returns ascending."""
    m = [row[:] for row in a]
    for _sweep in range(50):
        off = math.sqrt(m[0][1] ** 2 + m[0][2] ** 2 + m[1][2] ** 2)
        if off < 1e-15:
            break
        # largest off-diagonal pivot
        pairs = [(abs(m[0][1]), 0, 1), (abs(m[0][2]), 0, 2), (abs(m[1][2]), 1, 2)]
        _, p, q = max(pairs)
        if abs(m[p][q]) < 1e-18:
            break
        theta = (m[q][q] - m[p][p]) / (2.0 * m[p][q])
        t = math.copysign(1.0, theta) / (abs(theta) + math.sqrt(theta * theta + 1.0))
        c = 1.0 / math.sqrt(t * t + 1.0)
        s = t * c
        # apply rotation
        mpp, mqq, mpq = m[p][p], m[q][q], m[p][q]
        m[p][p] = c * c * mpp - 2.0 * s * c * mpq + s * s * mqq
        m[q][q] = s * s * mpp + 2.0 * s * c * mpq + c * c * mqq
        m[p][q] = m[q][p] = 0.0
        for r in range(3):
            if r in (p, q):
                continue
            mrp, mrq = m[r][p], m[r][q]
            m[r][p] = m[p][r] = c * mrp - s * mrq
            m[r][q] = m[q][r] = s * mrp + c * mrq
    eigs = [m[0][0], m[1][1], m[2][2]]
    eigs.sort()
    return eigs


def _shape_anisotropy(coords) -> float:
    """Relative shape anisotropy k^2 of the gyration tensor of N x 3
    coordinates (no masses): 0 isotropic, 1 collinear."""
    n = len(coords)
    if n < 2:
        return 0.0
    cx = sum(p[0] for p in coords) / n
    cy = sum(p[1] for p in coords) / n
    cz = sum(p[2] for p in coords) / n
    xx = yy = zz = xy = xz = yz = 0.0
    for px, py, pz in coords:
        dx, dy, dz = px - cx, py - cy, pz - cz
        xx += dx * dx
        yy += dy * dy
        zz += dz * dz
        xy += dx * dy
        xz += dx * dz
        yz += dy * dz
    g = [[xx / n, xy / n, xz / n],
         [xy / n, yy / n, yz / n],
         [xz / n, yz / n, zz / n]]
    l1, l2, l3 = _sym3_eigenvalues(g)
    s = l1 + l2 + l3
    if s <= 0:
        return 0.0
    k2 = 1.0 - 3.0 * (l1 * l2 + l2 * l3 + l3 * l1) / (s * s)
    # numerical guard into [0, 1)
    if k2 < 0 and k2 > -1e-12:
        k2 = 0.0
    return k2


def conformer_volume_3d(
    smiles: str,
    n_conformers: int = 10,
    seed: int = 42,
    tanford_volume_a3: float | None = None,
) -> ConformerVolumeResult:
    """C5: ETKDG conformer generation + MMFF/UFF minimization -> 3D vdW
    volume, SASA, shape anisotropy (cross-check for tanford_tail_volume
    on branched/gemini tails).

    Args:
        smiles: non-empty SMILES string (valid; ValueError names it).
        n_conformers: conformers to embed before energy ranking (>= 1).
        seed: ETKDG random seed (>= 0; fixed seed = reproducible output).
        tanford_volume_a3: optional caller-supplied value from
            `tanford_tail_volume` (A^3) -> tanford_ratio is returned.

    Returns:
        ConformerVolumeResult with the standard
        converged/n_iterations/method/params/caveat schema.
        Fail-fast ValueError naming every bad field.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import rdDistGeom, rdForceFieldHelpers, rdMolDescriptors
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "rdkit is required for conformer_volume_3d (surfactantkit core dep)"
        ) from exc

    if not isinstance(smiles, str) or not smiles.strip():
        raise ValueError("smiles must be a non-empty string")
    if isinstance(n_conformers, bool) or not isinstance(n_conformers, int):
        raise TypeError("n_conformers must be an integer")
    if n_conformers < 1:
        raise ValueError("n_conformers must be >= 1")
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise TypeError("seed must be an integer")
    if seed < 0:
        raise ValueError("seed must be >= 0 (fixed seed keeps output reproducible)")
    if tanford_volume_a3 is not None:
        if not isinstance(tanford_volume_a3, (int, float)) or (
            isinstance(tanford_volume_a3, bool)
        ):
            raise TypeError("tanford_volume_a3 must be a number")
        if tanford_volume_a3 <= 0:
            raise ValueError("tanford_volume_a3 must be > 0")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"smiles {smiles!r} is not parseable by rdkit")
    mol = Chem.AddHs(mol)
    params = rdDistGeom.ETKDGv3()
    params.randomSeed = seed
    params.numThreads = 1  # deterministic across machines
    cids = list(rdDistGeom.EmbedMultipleConfs(mol, numConfs=n_conformers, params=params))
    if not cids:
        raise ValueError(
            f"ETKDG embedding failed for smiles {smiles!r}: "
            "try another seed or a valid structure"
        )

    force_field = "none"
    energies: dict[int, float] = {}
    n_min_iters = 0
    n_minimized = 0
    converged_min = False
    mmff = rdForceFieldHelpers.MMFFGetMoleculeProperties(mol)
    if mmff is not None:
        force_field = "MMFF94s"
        for cid in cids:
            try:
                ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol, mmff, confId=cid)
                if ff is not None:
                    it = ff.Minimize(maxIts=500)
                    n_min_iters += 500 if it else 1
                    n_minimized += 1
                    energies[cid] = ff.CalcEnergy()
            except (ValueError, RuntimeError):
                continue
        converged_min = n_minimized == len(cids)
    else:
        # UFF: no separate properties object in current rdkit builds --
        # the force-field helper builds UFF params itself
        force_field = "UFF"
        for cid in cids:
            try:
                ff = rdForceFieldHelpers.UFFGetMoleculeForceField(
                    mol, confId=cid
                )
                if ff is not None:
                    it = ff.Minimize(maxIts=500)
                    n_min_iters += 500 if it else 1
                    n_minimized += 1
                    energies[cid] = ff.CalcEnergy()
            except (ValueError, RuntimeError, AttributeError):
                continue
        converged_min = n_minimized == len(cids)

    if energies:
        best = min(energies, key=lambda k: energies[k])
        e_low = energies[best]
    else:
        best = cids[0]
        e_low = None

    # promote the selected conformer to id 0 so downstream 3D descriptors
    # (DCLV, Labute ASA) unambiguously use it
    sel = Chem.Conformer(mol.GetConformer(best))
    mol.RemoveAllConformers()
    mol.AddConformer(sel, assignId=True)
    best = 0
    # DCLV via runtime dict lookup: class signature inference varies
    # across rdkit builds; runtime access is version-stable (rdkit >= 2020.09)
    dclv = rdMolDescriptors.__dict__.get("DoubleCubicLatticeVolume")
    if dclv is None:
        raise ImportError(
            "rdkit DoubleCubicLatticeVolume missing -- need rdkit >= 2020.09"
        )
    dvl = dclv(mol, confId=0)  # type: ignore[operator]  # runtime-stable API
    vdw = float(dvl.GetVDWVolume())
    sasa = float(rdMolDescriptors.CalcLabuteASA(mol))
    heavy = [
        (mol.GetConformer(best).GetAtomPosition(i))
        for i in range(mol.GetNumAtoms())
        if mol.GetAtomWithIdx(i).GetAtomicNum() > 1
    ]
    coords = [(p.x, p.y, p.z) for p in heavy]
    k2 = _shape_anisotropy(coords)

    ratio = (vdw / float(tanford_volume_a3)
             if tanford_volume_a3 is not None else None)

    caveats = [
        (
            "whole-molecule property of the lowest-energy ETKDG/MMFF "
            "conformer (lowest of N sampled, not a free-energy global minimum)"
        ),
        (
            "shape anisotropy k^2 from the heavy-atom gyration tensor "
            "(0 isotropic, 1 rod-like)"
        ),
    ]
    if tanford_volume_a3 is not None:
        caveats.append(
            "tanford_ratio compares whole-molecule 3D vdW volume against the "
            "TAIL-ONLY linear-chain group contribution: definitional "
            "cross-check (divergence expected and informative for branched/"
            "gemini tails, not an error)"
        )
    if force_field == "none":
        caveats.append(
            "no force field available for this structure: ETKDG "
            "coordinates used unminimized (converged=False)"
        )
    elif not converged_min:
        caveats.append(
            "some conformers did not converge in 500 minimization iterations"
        )

    result_params = {
        "n_conformers": n_conformers,
        "seed": seed,
        "force_field": force_field,
        "selected_conf_id": int(best),
        "n_conformers_embedded": len(cids),
        "n_conformers_minimized": n_minimized,
        "heavy_atom_count": len(coords),
        "tanford_volume_a3": (float(tanford_volume_a3)
                              if tanford_volume_a3 is not None else None),
        "molecular_formula": rdMolDescriptors.CalcMolFormula(mol),
    }
    return ConformerVolumeResult(
        smiles=smiles,
        vdW_volume_A3=vdw,
        sasa_A2=sasa,
        shape_anisotropy_k2=k2,
        lowest_energy_kcal_mol=e_low,
        force_field=force_field,
        n_conformers_embedded=len(cids),
        n_conformers_minimized=n_minimized,
        tanford_ratio=ratio,
        converged=bool(converged_min and energies),
        n_iterations=n_min_iters,
        method=(f"ETKDGv3 (seed={seed}, 1 thread) + {force_field} minimization "
                "(<=500 its/conf) -> DCLV vdW volume, Labute SASA, "
                "gyration-tensor k^2 (hand Jacobi)"),
        params=result_params,
        caveat="; ".join(caveats),
    )
