"""Tests for the superset scenario registry (config.py), loaded directly."""

import importlib.util
import os
import sys
import types

import pytest

_prospective_dir = os.path.join(
    os.path.dirname(__file__), "..", "dynamic_characterization", "prospective"
)
sys.modules.setdefault(
    "dynamic_characterization", types.ModuleType("dynamic_characterization")
)
_pkg = types.ModuleType("dynamic_characterization.prospective")
_pkg.__path__ = [_prospective_dir]
sys.modules["dynamic_characterization.prospective"] = _pkg
_spec = importlib.util.spec_from_file_location(
    "dynamic_characterization.prospective.config",
    os.path.join(_prospective_dir, "config.py"),
)
config = importlib.util.module_from_spec(_spec)
sys.modules["dynamic_characterization.prospective.config"] = config
_spec.loader.exec_module(config)


@pytest.fixture(autouse=True)
def _reset():
    config.reset_scenario()
    yield
    config.reset_scenario()


def test_prospective_scenarios_count():
    assert len(config.available_scenarios("prospective")) == 18


def test_fair_native_scenarios_present():
    fair_native = set(config.available_fair_scenarios())
    assert ("SSP3", "7.0") in fair_native
    assert ("SSP5", "3.4-over") in fair_native
    assert len(fair_native) == 8


def test_available_scenarios_fair_lists_only_dual_iam_scenarios():
    fair = set(config.available_scenarios("fair"))
    # the 4 dual-support IAM scenarios
    assert fair == {
        ("GCAM4", "SSP4", "6.0"),
        ("IMAGE", "SSP1", "2.6"),
        ("MESSAGE", "SSP2", "4.5"),
        ("REMIND", "SSP5", "8.5"),
    }
    # FAIR-native entries are NOT in the IAM registry
    assert ("FAIR", "SSP3", "7.0") not in fair


def test_dual_support_scenario_metrics():
    entry = config.SCENARIO_REGISTRY[("IMAGE", "SSP1", "2.6")]
    assert entry["metrics"] == frozenset({"prospective", "fair"})
    assert entry["fair_marker"] == "ssp126"


def test_set_fair_scenario_marker_and_support():
    config.set_fair_scenario("SSP3", "7.0")
    assert config.scenario_supports("fair") is True
    assert config.scenario_supports("prospective") is False
    assert config.current_fair_marker() == "ssp370"
    assert config.get_scenario()["iam"] is None


def test_set_fair_scenario_rejects_unknown():
    with pytest.raises(ValueError, match="Invalid FAIR scenario"):
        config.set_fair_scenario("SSP3", "4.5")  # not a FAIR marker


def test_scenario_supports_and_marker():
    config.set_scenario("MESSAGE", "SSP2", "4.5")
    assert config.scenario_supports("fair") is True
    assert config.scenario_supports("prospective") is True
    assert config.current_fair_marker() == "ssp245"


def test_prospective_only_scenario_has_no_fair_marker():
    config.set_scenario("AIM", "SSP3", "4.5")
    assert config.scenario_supports("fair") is False
    with pytest.raises(ValueError, match="fair"):
        config.current_fair_marker()


def test_metric_family():
    assert config.metric_family("pGWP") == "prospective"
    assert config.metric_family("fair_temperature") == "fair"
