"""Benchmark golds that the toolkit's orchestrator computes must equal what the
toolkit reports to a caller, i.e. the AUTHORITATIVE OrchestratorResult.cmc_mM,
not the raw piecewise r.cmc.cmc_mM.

GZ-17 (SB-017) shipped with the raw piecewise 8.185 mM, the geometric mean of
the 4.09/16.37 mM gap, although the toolkit refines that suspect plateau to
16.37 mM (the dev-repo correction of 2026-09-23 was lost when the release gold
was assembled; restored 2026-10-23). GZ-16 is the negative control: its plateau
is genuinely flat, so raw and authoritative values agree.
"""

import json
from pathlib import Path

import pytest

from surfactantkit.orchestrate import derive_all_properties_from_smiles_and_curve

GOLD = json.loads((Path(__file__).resolve().parent.parent / "benchmark" / "source_sets" / "groundzero"
                   / "gold_answers.json").read_text(encoding="utf-8"))

# the same PendantProp curves as tests/test_orchestrate.py
CASES = {
    "GZ-16": ("CCCCCCCCCCCCOS(=O)(=O)[O-].[Na+]", 296.15,
              [0.40577109375, 0.8115421875, 1.623084375, 3.24616875, 6.4923375, 12.984675, 25.96935, 51.9387],
              [69.85933451975876, 67.9370105943862, 62.892649575815454, 54.57301518905874,
               42.87619733622359, 37.91402789831731, 37.44295519808575, 36.79102755846683]),
    "GZ-17": ("CCCCCCCCCCCC[N+](C)(C)C.[Br-]", 296.35,
              [0.25578125, 0.5115625, 1.023125, 2.04625, 4.0925, 16.37, 32.74],
              [70.26907673, 68.20699668, 65.39702756, 57.50386199, 49.23642286, 32.34750118, 36.22382546]),
}


@pytest.mark.parametrize("qid", sorted(CASES))
def test_groundzero_cmc_gold_is_the_authoritative_orchestrator_value(qid):
    smiles, temperature_K, concs, gammas = CASES[qid]
    r = derive_all_properties_from_smiles_and_curve(smiles, concs, gammas, temperature_K=temperature_K)
    gold = GOLD[qid]["gold"]

    assert gold["cmc_mM"] == pytest.approx(r.cmc_mM, rel=1e-3)
    assert gold["gamma_max_mol_per_m2"] == pytest.approx(r.gamma_max_mol_per_m2, rel=1e-3)
    assert gold["a_min_nm2"] == pytest.approx(r.a_min_nm2, rel=1e-3)
