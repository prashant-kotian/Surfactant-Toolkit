"""M1 -- analysis-only worker over existing micelle trajectories.

Starts from the repo's own output formats (per the Phase-3 plan):
  * `.npz`  -- key `positions` float array shaped (T, N, 3); optional
    `time` (T,), `box` (3,) or (T, 3) as box LENGTHS with origin at 0,
    `masses` (N,)
  * `.txt`/`.dat`/`.xyz` -- XYZ frames (natoms line + comment line +
    element x y z rows) or blank-line-separated `x y z` rows

Computes per frame: cluster count and full cluster-size list (single-
linkage connected components under a distance cutoff), largest-cluster
size N_agg(t), radius of gyration Rg(t) of that cluster (mass-weighted
if `masses` given, unweighted otherwise), radial number-density profile
of the largest cluster about its own centre, and -- when `box` exists --
the slab density profile along z.

Honesty: single-linkage clustering at a caller-supplied cutoff (the
cutoff drives every answer; it is echoed in params, and bridging
sensitivity is disclosed); analysis only -- does NOT validate the
underlying CG/atomistic model; density profiles are number density in
the input length unit (no unit conversion guessed).

Dependency: numpy is imported lazily (house `[calc]` policy).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

try:  # numpy only inside the worker: keep package import light
    import numpy as _np
except ImportError:  # pragma: no cover - exercised only without the extra
    _np = None


@dataclass(frozen=True)
class TrajectoryAnalysisResult:
    """Result of `analyze_micelle_trajectory`."""

    times: tuple[float, ...]             # frame index if no time axis given
    n_frames: int
    n_particles: int
    n_clusters: tuple[int, ...]              # per frame
    cluster_sizes: tuple[tuple[int, ...], ...]  # per frame, descending
    largest_cluster_size: tuple[int, ...]   # N_agg(t)
    rg_largest: tuple[float, ...]           # Rg(t) of largest cluster
    radial_density_r: tuple[float, ...]     # bin centres
    radial_density_values: tuple[float, ...]  # number density per shell volume
    slab_density_z: tuple[float, ...] | None   # None when no box given
    converged: bool                   # input parsed and all frames analyzable
    n_iterations: int                 # frames analyzed
    method: str
    params: dict
    caveat: str


def _require_numpy():
    if _np is None:
        raise ImportError(
            "numpy is required for analyze_micelle_trajectory; "
            "install the calc extra: pip install surfactantkit[calc]"
        )
    return _np


def _xyz_row(toks: list[str]) -> list[float]:
    """Parse one XYZ coordinate line: element-prefixed (e x y z) or bare
    numeric (x y z [..]). Use explicit x y z columns, never positional
    guessing from the end."""
    if len(toks) >= 4:
        try:
            float(toks[0])
            # bare numeric row: x y z are the first three
            return [float(toks[0]), float(toks[1]), float(toks[2])]
        except ValueError:
            # element then x y z
            return [float(toks[1]), float(toks[2]), float(toks[3])]
    if len(toks) == 3:
        return [float(toks[0]), float(toks[1]), float(toks[2])]
    raise ValueError("coordinate line needs >= 3 columns")


def _parse_text_frames(lines, np):
    """Segment a text file into frames: XYZ blocks and/or blank-line
    separated bare x y z row groups. Returns (F, N, 3) or None."""
    frames: list[list[list[float]]] = []
    cur: list[list[float]] = []
    i = 0
    n = len(lines)
    saw_any = False

    def flush():
        nonlocal cur
        if cur:
            frames.append(cur)
            cur = []

    while i < n:
        line = lines[i].strip()
        if not line:
            flush()
            i += 1
            continue
        parts = line.split()
        if len(parts) == 1 and parts[0].isdigit() and not cur:
            n_at = int(parts[0])
            if i + 1 >= n:
                raise ValueError("XYZ header truncated before comment line")
            i += 2
            block = []
            for _ in range(n_at):
                if i >= n:
                    raise ValueError(
                        "XYZ frame truncated (expected natoms coordinate lines)"
                    )
                try:
                    block.append(_xyz_row(lines[i].split()))
                except ValueError:
                    raise ValueError(
                        f"unparseable XYZ coordinate line: {lines[i]!r}"
                    ) from None
                i += 1
            frames.append(block)
            saw_any = True
            continue
        try:
            cur.append(_xyz_row(parts))
        except ValueError:
            raise ValueError(f"unparseable trajectory line: {line!r}") from None
        saw_any = True
        i += 1
    flush()
    if not saw_any or not frames:
        return None
    try:
        arr = np.asarray(frames, dtype=float)
    except ValueError:
        raise ValueError(
            "text trajectory frames have differing particle counts "
            "(pad frames or use npz)"
        ) from None
    if arr.ndim != 3:
        return None
    return arr


def _load_trajectory(source, np):
    """Return (positions (T,N,3), times|None, box|None, masses|None)."""
    if str(source).lower().endswith(".npz"):
        try:
            data = np.load(source, allow_pickle=False)
        except OSError as exc:
            raise ValueError(f"cannot read npz file {source!r}: {exc}") from exc
        if "positions" not in data:
            raise ValueError(
                "npz must contain a 'positions' array of shape (T, N, 3)"
            )
        pos = np.asarray(data["positions"], dtype=float)
        times = (np.asarray(data["time"], dtype=float)
                 if "time" in data else None)
        box = (np.asarray(data["box"], dtype=float)
               if "box" in data else None)
        masses = (np.asarray(data["masses"], dtype=float)
                  if "masses" in data else None)
    else:
        try:
            with open(source) as fh:
                lines = fh.read().splitlines()
        except OSError as exc:
            raise ValueError(f"cannot read trajectory {source!r}: {exc}") from exc
        parsed = _parse_text_frames(lines, np)
        if parsed is None:
            raise ValueError(
                "no coordinate frames found (XYZ blocks or blank-line "
                "separated x y z rows expected)"
            )
        pos, times, box, masses = parsed, None, None, None
    if pos.ndim != 3 or pos.shape[2] != 3:
        raise ValueError(
            f"positions must have shape (T, N, 3), got {tuple(pos.shape)}"
        )
    if pos.shape[0] == 0 or pos.shape[1] == 0:
        raise ValueError("trajectory has zero frames or zero particles")
    if not np.isfinite(pos).all():
        raise ValueError("positions contain NaN or inf")
    if times is not None and len(times) != pos.shape[0]:
        raise ValueError("time axis length does not match number of frames")
    if masses is not None and len(masses) != pos.shape[1]:
        raise ValueError("masses length does not match particle count")
    if box is not None and box.shape[-1] != 3:
        raise ValueError("box must be (3,) or (T, 3)")
    return pos, times, box, masses


def _labels_single_linkage(xyz, cutoff, np):
    """Connected components (single linkage) from a boolean adjacency."""
    n = xyz.shape[0]
    d2 = np.sum((xyz[:, None, :] - xyz[None, :, :]) ** 2, axis=-1)
    parent = list(range(n))

    def find(a: int) -> int:
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    ii, jj = np.nonzero(np.triu(d2 < cutoff * cutoff, k=1))
    for a, b in zip(ii.tolist(), jj.tolist(), strict=True):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra
    root_to_label: dict[int, int] = {}
    labels = [0] * n
    for k in range(n):
        rt = find(k)
        if rt not in root_to_label:
            root_to_label[rt] = len(root_to_label)
        labels[k] = root_to_label[rt]
    return labels, len(root_to_label)


def analyze_micelle_trajectory(
    source,
    cluster_cutoff: float,
    n_radial_bins: int = 30,
    compute_slab_profile: bool = True,
) -> TrajectoryAnalysisResult:
    """M1: Nagg(t), cluster sizes, Rg(t) and density profiles from an
    existing CG/atomistic micelle trajectory (analysis only).

    Args:
        source: path to `.npz` (keys: positions (T,N,3); optional time,
            box as lengths with origin 0, masses) or text trajectory
            (XYZ blocks or blank-line separated x y z rows).
        cluster_cutoff: neighbour distance for single-linkage clustering,
            in the SAME length unit as the coordinates (the answer is
            cutoff-sensitive: echoed in params/caveat).
        n_radial_bins: bins for the largest-cluster radial number-density
            profile (integer >= 4).
        compute_slab_profile: also return the box-z slab density profile
            when a `box` is present (bool).

    Returns:
        TrajectoryAnalysisResult with the standard
        converged/n_iterations/method/params/caveat schema.
        Fail-fast ValueError naming every bad field.
    """
    np = _require_numpy()
    if isinstance(cluster_cutoff, bool) or not isinstance(
        cluster_cutoff, (int, float)
    ):
        raise TypeError("cluster_cutoff must be a real number")
    if not math.isfinite(float(cluster_cutoff)) or cluster_cutoff <= 0:
        raise ValueError("cluster_cutoff must be a finite number > 0")
    if isinstance(n_radial_bins, bool) or not isinstance(n_radial_bins, int):
        raise TypeError("n_radial_bins must be an integer >= 4")
    if n_radial_bins < 4:
        raise ValueError("n_radial_bins must be an integer >= 4")
    if not isinstance(compute_slab_profile, bool):
        raise TypeError("compute_slab_profile must be a bool")

    pos, times, box, masses = _load_trajectory(source, np)
    n_frames, n_particles = int(pos.shape[0]), int(pos.shape[1])

    cluster_count: list[int] = []
    sizes_all: list[tuple[int, ...]] = []
    largest: list[int] = []
    rg_series: list[float] = []
    rad_hist = np.zeros(n_radial_bins, dtype=float)
    rad_frames = 0
    slab_hist = np.zeros(50, dtype=float)
    slab_frames = 0
    slab_lo = slab_hi = 0.0
    ok = True
    has_box = box is not None
    mass_weighted = masses is not None
    last_centres = [0.0] * n_radial_bins

    for f in range(n_frames):
        xyz = pos[f]
        labels, n_cl = _labels_single_linkage(xyz, cluster_cutoff, np)
        counts: dict[int, int] = {}
        for lb in labels:
            counts[lb] = counts.get(lb, 0) + 1
        order = sorted(counts.values(), reverse=True)
        big = max(counts, key=lambda k: counts[k])
        members = np.nonzero(np.asarray(labels) == big)[0]
        cluster_count.append(n_cl)
        sizes_all.append(tuple(order))
        largest.append(order[0])
        try:
            sub = xyz[members]
            if masses is not None:
                w = masses[members]
                wsum = float(w.sum())
                if wsum <= 0:
                    raise ValueError("largest-cluster masses sum to zero")
                com = (w[:, None] * sub).sum(axis=0) / wsum
                rg = float(np.sqrt(
                    float((w * np.sum((sub - com) ** 2, axis=1)).sum()) / wsum
                ))
            else:
                com = sub.mean(axis=0)
                rg = float(np.sqrt(float(np.mean(
                    np.sum((sub - com) ** 2, axis=1)
                ))))
            rad = np.sqrt(np.sum((sub - com) ** 2, axis=1))
            r_max = float(rad.max())
            if r_max > 0:
                hist, edges = np.histogram(
                    rad, bins=n_radial_bins, range=(0.0, r_max)
                )
                shell = (4.0 / 3.0 * math.pi
                         * (edges[1:] ** 3 - edges[:-1] ** 3))
                with np.errstate(divide="ignore", invalid="ignore"):
                    dens = np.where(shell > 0, hist / shell, 0.0)
                rad_hist = rad_hist + dens
                rad_frames += 1
                last_centres = ((edges[:-1] + edges[1:]) / 2).tolist()
            rg_series.append(rg)
            if compute_slab_profile and has_box:
                b = box if getattr(box, "ndim", 1) == 1 else box[f]
                slab_lo, slab_hi = 0.0, float(b[2])
                histz, _ = np.histogram(
                    xyz[:, 2], bins=50, range=(slab_lo, slab_hi)
                )
                slab_hist = slab_hist + histz
                slab_frames += 1
        except (ValueError, FloatingPointError):
            ok = False
        # note: cluster bookkeeping above always succeeds

    n_rg = max(rad_frames, 1)
    radial_vals = (rad_hist / n_rg).tolist()
    if times is None:
        times_t = tuple(float(f) for f in range(n_frames))
    else:
        times_t = tuple(float(t) for t in times)

    slab: list[float] | None = None
    if compute_slab_profile and has_box and slab_frames > 0:
        dz = (slab_hi - slab_lo) / 50
        b0 = box if getattr(box, "ndim", 1) == 1 else box[0]
        dx_area = float(b0[0]) * float(b0[1])
        area = dx_area if dx_area > 0 else 1.0
        slab = [float(h) / (slab_frames * dz * area)
                for h in slab_hist.tolist()]
    else:
        slab = None

    caveat = (
        f"single-linkage clustering at caller cutoff {cluster_cutoff} "
        "(same unit as coordinates): bridging can fuse distinct aggregates; "
        "cutoff sensitivity is inherent; analysis only -- does not validate "
        "the underlying CG/atomistic model; radial/slab profiles are number "
        "density in input units (no unit conversion guessed)"
    )
    if not has_box:
        caveat += "; no box array: slab density profile unavailable"
    if not mass_weighted:
        caveat += ("; unweighted (geometric) Rg "
                   "(provide npz masses for mass-weighted)")

    params = {
        "cluster_cutoff": float(cluster_cutoff),
        "n_radial_bins": n_radial_bins,
        "compute_slab_profile": compute_slab_profile,
        "has_box": bool(has_box),
        "mass_weighted_rg": bool(mass_weighted),
        "input_shape": [n_frames, n_particles, 3],
        "mean_largest_cluster": (sum(largest) / len(largest)) if largest else 0.0,
    }
    return TrajectoryAnalysisResult(
        times=times_t,
        n_frames=n_frames,
        n_particles=n_particles,
        n_clusters=tuple(cluster_count),
        cluster_sizes=tuple(sizes_all),
        largest_cluster_size=tuple(largest),
        rg_largest=tuple(rg_series),
        radial_density_r=tuple(last_centres),
        radial_density_values=tuple(radial_vals),
        slab_density_z=tuple(slab) if slab is not None else None,
        converged=bool(ok and rad_frames + 0 >= 0),
        n_iterations=n_frames,
        method=("single-linkage distance clustering per frame + gyration "
                "Rg (mass-weighted when masses given) + shell-binned "
                "number density"),
        params=params,
        caveat=caveat,
    )
