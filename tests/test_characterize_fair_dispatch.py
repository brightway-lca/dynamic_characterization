# tests/test_characterize_fair_dispatch.py
import pandas as pd
import pytest

import dynamic_characterization as dc
import dynamic_characterization.prospective as prospective
from dynamic_characterization import dynamic_characterization as dcmod


@pytest.fixture(autouse=True)
def _reset():
    prospective.reset_scenario()
    yield
    prospective.reset_scenario()


def _inv():
    return pd.DataFrame(
        {
            "date": pd.to_datetime(["2030-01-01"]),
            "amount": [10.0],
            "flow": [1],
            "activity": ["a"],
        }
    )


def test_invalid_metric_rejected():
    with pytest.raises(ValueError, match="Metric must be one of"):
        dc.characterize(_inv(), metric="nonsense")


def test_fair_metric_under_prospective_only_scenario_raises():
    prospective.set_scenario("AIM", "SSP3", "4.5")  # prospective-only
    with pytest.raises(ValueError, match="fair"):
        dc.characterize(_inv(), metric="fair_temperature")


def test_fair_metric_dispatches_to_core(monkeypatch):
    prospective.set_fair_scenario("SSP2", "4.5")
    called = {}

    def fake(df, output="radiative_forcing", quantiles=(50.0,),
             time_horizon=None, workers=None):
        called["output"] = output
        return pd.DataFrame(
            columns=["date", "amount", "flow", "activity", "quantile"]
        )

    # patch the lazily-imported symbol used inside characterize
    import dynamic_characterization.fair.core as core
    monkeypatch.setattr(core, "characterize_with_fair", fake)
    dc.characterize(_inv(), metric="fair_radiative_forcing")
    assert called["output"] == "radiative_forcing"


def test_fair_metric_without_any_scenario_raises():
    # No set_scenario call (autouse fixture resets it).
    with pytest.raises(ValueError, match="no scenario is set"):
        dc.characterize(_inv(), metric="fair_temperature")


def test_add_fair_flow_names_is_best_effort():
    # No bw2 project/biosphere needed: helper must not raise, returns a frame.
    df = _inv()
    out = dcmod._add_fair_flow_names(df)
    assert "flow" in out.columns  # never crashes; may or may not add flow_name
    assert len(out) == len(df)  # frame is not lost
