"""Shared internal numerics helpers (audit TOOLKIT_AUDIT_2026-09-30
P2-9 consolidation).

* `golden_section_minimize` -- previously duplicated in adsorption.py
  (iters=200, tol=1e-14) and mixed_micelle.py (iters=60, tol=1e-13);
  both variants are preserved as explicit parameters at each call site
  so behavior is unchanged.
* `solve_3x3` -- previously duplicated in thermodynamics.py and
  wetting.py with module-specific error messages; the message is now a
  required parameter so each caller keeps its exact wording.

Private module -- not part of the public API, not an MCP tool.
"""

from __future__ import annotations

import math


def golden_section_minimize(f, lo: float, hi: float, iters: int = 200,
                            tol: float = 1e-14) -> float:
    """Golden-section search for the minimizer of a UNIMODAL function f
    on [lo, hi]. Pure-Python (no scipy) -- SurfactantKit is
    zero-external-dependency by design. `iters` and `tol` expose the
    two historical call-site variants exactly (see module docstring).
    """
    gr = (math.sqrt(5.0) - 1.0) / 2.0
    c = hi - gr * (hi - lo)
    d = lo + gr * (hi - lo)
    fc, fd = f(c), f(d)
    for _ in range(iters):
        if fc < fd:
            hi, d, fd = d, c, fc
            c = hi - gr * (hi - lo)
            fc = f(c)
        else:
            lo, c, fc = c, d, fd
            d = lo + gr * (hi - lo)
            fd = f(d)
        if hi - lo < tol:
            break
    return (lo + hi) / 2.0


def solve_3x3(matrix: list[list[float]], vector: list[float],
              singular_message: str) -> list[float]:
    """Gaussian elimination with partial pivoting for a 3x3 linear
    system. Small, pure-Python -- avoids numpy for a single fixed-size
    solve. Raises ValueError(singular_message) on a degenerate pivot so
    each caller keeps its own field-specific wording.
    """
    a = [row[:] + [vector[i]] for i, row in enumerate(matrix)]
    n = 3
    for col in range(n):
        pivot_row = max(range(col, n), key=lambda r: abs(a[r][col]))
        if abs(a[pivot_row][col]) < 1e-15:
            raise ValueError(singular_message)
        a[col], a[pivot_row] = a[pivot_row], a[col]
        for r in range(col + 1, n):
            factor = a[r][col] / a[col][col]
            for c in range(col, n + 1):
                a[r][c] -= factor * a[col][c]
    # back substitution
    x = [0.0] * n
    for r in range(n - 1, -1, -1):
        s = a[r][n] - sum(a[r][c] * x[c] for c in range(r + 1, n))
        x[r] = s / a[r][r]
    return x
