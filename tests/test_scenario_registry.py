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


def test_fair_native_markers_present():
    fair = set(config.available_scenarios("fair"))
    assert ("FAIR", "SSP3", "7.0") in fair
    assert ("FAIR", "SSP5", "3.4-over") in fair
    # the 4 dual-support prospective scenarios are also fair-capable
    assert ("IMAGE", "SSP1", "2.6") in fair
    assert ("MESSAGE", "SSP2", "4.5") in fair


def test_dual_support_scenario_metrics():
    entry = config.SCENARIO_REGISTRY[("IMAGE", "SSP1", "2.6")]
    assert entry["metrics"] == frozenset({"prospective", "fair"})
    assert entry["fair_marker"] == "ssp126"


def test_fair_native_is_fair_only():
    entry = config.SCENARIO_REGISTRY[("FAIR", "SSP3", "7.0")]
    assert entry["metrics"] == frozenset({"fair"})
    assert entry["fair_marker"] == "ssp370"


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
