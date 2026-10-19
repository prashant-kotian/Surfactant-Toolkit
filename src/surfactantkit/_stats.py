"""Shared internal statistics helpers (audit TOOLKIT_AUDIT_2026-09-30
P2-9 consolidation; previously duplicated as private `_bic` copies in
adsorption.py, curve_analysis.py and mixture_model_selection.py).

Private module -- not part of the public API, not an MCP tool.
"""

from __future__ import annotations

import math


def bic(total_rss: float, n: int, n_params: int) -> float:
    """Bayesian Information Criterion for a least-squares fit,
    n*ln(RSS/n) + n_params*ln(n) -- lower is better; the n_params
    penalty is what makes this a genuine model-comparison metric rather
    than "whichever model has more free parameters always wins".
    Numerically identical to the three former private copies (guard
    against log(0) on a mathematically-perfect fit via 1e-300 floor).
    """
    safe_rss = total_rss if total_rss > 0.0 else 1e-300
    return n * math.log(safe_rss / n) + n_params * math.log(n)
