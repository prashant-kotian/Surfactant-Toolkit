"""Regression tests for benchmark/build_figure_data.py's F8 failure-anatomy inputs."""
import importlib.util
import json
from pathlib import Path

import pytest

_SPEC = importlib.util.spec_from_file_location(
    "build_figure_data", Path(__file__).resolve().parents[1] / "benchmark" / "build_figure_data.py")
bfd = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(bfd)


def _field(ok, gold, model):
    return {"ok": ok, "gold": gold, "model": model}


def test_multi_model_json_does_not_cross_attribute(tmp_path):
    # Same qid under several models: each model must get its own record, not the
    # last one walked (the bug that labelled failed Gemini answers "all_pass").
    p = tmp_path / "multi.json"
    p.write_text(json.dumps({
        "chatgpt": {"GZ-01": {"verdict": "PASS", "fields": {"x": _field(True, 1.0, 1.0)}}},
        "gemini": {"GZ-01": {"verdict": "FAIL", "fields": {"x": _field(False, 1.0, 5.0)}}},
        "claude": {"GZ-01": {"verdict": "PASS", "fields": {"x": _field(True, 1.0, 1.0)}}},
    }), encoding="utf-8")
    assert bfd.classify(bfd.normalize_detail(p, "Gemini")["GZ-01"]) == "magnitude_error"
    assert bfd.normalize_detail(p, "ChatGPT")["GZ-01"][0][0] is True


def test_gemini_cli_resolves_to_antigravity_key(tmp_path):
    p = tmp_path / "ag.json"
    p.write_text(json.dumps({
        "Antigravity": {"H-01": {"verdict": "FAIL", "fields": {"s": _field(False, "a", "b")}}},
        "Claude": {"H-01": {"verdict": "PASS", "fields": {"s": _field(True, "a", "a")}}},
    }), encoding="utf-8")
    assert bfd.classify(bfd.normalize_detail(p, "Gemini-CLI")["H-01"]) == "selection_or_convention"


def test_bare_bool_fields_are_count_only(tmp_path):
    p = tmp_path / "bools.json"
    p.write_text(json.dumps({"X-01": {"pass": False, "fields": {"a": False, "b": True}}}),
                 encoding="utf-8")
    assert bfd.classify(bfd.normalize_detail(p, "Gemini-CLI")["X-01"]) == "unclassified"


def test_failed_record_without_fields_is_no_answer():
    assert bfd.classify([]) == "no_answer"


def test_failed_record_with_only_passing_fields_fails_loudly():
    with pytest.raises(SystemExit):
        bfd.classify([(True, 1.0, 1.0)])
