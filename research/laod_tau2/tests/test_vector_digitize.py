"""End-to-end test of the vector-PDF extractor against a known function.

Guards the claim that matters for data mining: extracting from a vector PDF is
~30x more accurate than clicking on pixels (0.06% vs 1-3%), so the result is
reproducible from the PDF alone.
"""

from __future__ import annotations

import warnings

import numpy as np
import pytest

pytest.importorskip("matplotlib")
pytest.importorskip("pymupdf")

from research.laod_tau2.vector_digitize import validate_synthetic


@pytest.fixture(scope="module")
def synthetic_pdf(tmp_path_factory):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t = np.logspace(-3, 0, 300)
    y = 72.0 - 25.0 * (1.0 - np.exp(-np.sqrt(t / 0.02)))
    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    ax.semilogx(t, y, "-", color="#cc2222", lw=1.0)
    ax.set_xlabel("t (s)")
    ax.set_ylabel("gamma (mN/m)")
    path = tmp_path_factory.mktemp("fig") / "vector.pdf"
    fig.savefig(path, format="pdf")
    plt.close(fig)
    return str(path), (lambda x: 72.0 - 25.0 * (1.0 - np.exp(-np.sqrt(x / 0.02))))


def test_vector_extraction_is_sub_percent(synthetic_pdf):
    path, truth = synthetic_pdf
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = validate_synthetic(path, truth)
    assert r["n_points"] > 20
    assert r["median_rel_error"] < 2e-3      # 0.2%
    assert r["max_rel_error"] < 5e-3


def test_calibration_recovers_log_ticks(synthetic_pdf):
    path, truth = synthetic_pdf
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = validate_synthetic(path, truth)
    decades = [e for _, e in r["calibration"]["x_ticks"]]
    assert decades == [-3.0, -2.0, -1.0, 0.0]
