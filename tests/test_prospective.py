"""
Prospective characterization (scenario data, configuration and the public API of
``prospective.radiative_forcing``, ``prospective.agwp`` and ``prospective.agtp``).

Since version 2.0 the public functions compute the corrected Barbosa Watanabe & Cherubini (2026)
method; ``test_paper_watanabe2026.py`` checks that method against the paper.
"""

import warnings
from collections import namedtuple

import numpy as np
import pandas as pd
import pytest

from dynamic_characterization import characterize
from dynamic_characterization.prospective import agtp, agwp, config, data_loader, watanabe2026
from dynamic_characterization.prospective import radiative_forcing as rf


@pytest.fixture(autouse=True)
def reset_scenario():
    """Reset scenario before each test."""
    config.reset_scenario()
    yield
    config.reset_scenario()


def test_load_irf_ch4():
    """IRF CH4 should be a 101-element array with first value 1.0."""
    irf = data_loader.load_irf_ch4()
    assert len(irf) == 101
    assert irf[0] == pytest.approx(1.0)
    # After 1 year, should decay by exp(-1/11.8)
    assert irf[1] == pytest.approx(0.918746, rel=1e-4)


def test_load_irf_co2():
    """IRF CO2 should have 4 RCP scenarios."""
    irf = data_loader.load_irf_co2()
    assert "RCP26" in irf
    assert "RCP45" in irf
    assert "RCP60" in irf
    assert "RCP85" in irf
    assert len(irf["RCP26"]) == 101
    assert irf["RCP26"][0] == pytest.approx(1.0)


def test_load_irf_n2o():
    """IRF N2O should be a 100-element array with first value 1.0."""
    irf = data_loader.load_irf_n2o()
    assert len(irf) == 100
    assert irf[0] == pytest.approx(1.0)


def test_load_re_ch4_image():
    """RE CH4 for IMAGE should have SSP1 with multiple RCPs."""
    re = data_loader.load_re_ch4("IMAGE")
    assert "SSP1" in re
    assert "2.6" in re["SSP1"]
    assert "_years" in re
    # Years should span 2020-2150
    assert re["_years"][0] == 2020
    assert re["_years"][-1] == 2150


def test_load_re_co2_gcam4():
    """RE CO2 for GCAM4 should have SSP4 scenarios."""
    re = data_loader.load_re_co2("GCAM4")
    assert "SSP4" in re
    assert "4.5" in re["SSP4"]


def test_load_re_n2o_aim():
    """RE N2O for AIM should have SSP scenarios with RCPs."""
    re = data_loader.load_re_n2o("AIM")
    assert "SSP3" in re
    assert "4.5" in re["SSP3"]
    assert "_years" in re
    # Years should span 2020-2150
    assert re["_years"][0] == 2020
    assert re["_years"][-1] == 2150


def test_load_re_invalid_iam():
    """Loading RE with invalid IAM should raise ValueError."""
    with pytest.raises(ValueError, match="Unknown IAM"):
        data_loader.load_re_ch4("INVALID")


def test_set_scenario():
    """Setting scenario should store configuration."""
    config.set_scenario(iam="IMAGE", ssp="SSP1", rcp="2.6")
    scenario = config.get_scenario()
    assert scenario["iam"] == "IMAGE"
    assert scenario["ssp"] == "SSP1"
    assert scenario["rcp"] == "2.6"


def test_set_invalid_scenario():
    """Setting invalid scenario combination should raise."""
    with pytest.raises(ValueError, match="Invalid scenario"):
        config.set_scenario(iam="IMAGE", ssp="SSP3", rcp="2.6")


def test_get_scenario_without_set():
    """Getting scenario without setting should raise RuntimeError."""
    with pytest.raises(RuntimeError, match="No scenario set"):
        config.get_scenario()


def test_get_scenario_without_set_message_is_actionable():
    """The error has to be usable by someone who only knows the calling package."""
    with pytest.raises(RuntimeError) as excinfo:
        config.get_scenario()

    message = str(excinfo.value)

    # A copy-pasteable fix, with the full import path: users reaching this via
    # bw_timex have never imported dynamic_characterization themselves.
    assert "from dynamic_characterization.prospective import set_scenario" in message
    assert 'set_scenario(iam="IMAGE", ssp="SSP1", rcp="2.6")' in message

    # Which calculations need it, and which scenarios exist.
    assert "pGWP" in message
    assert "IMAGE-SSP1" in message
    assert "REMIND-SSP5" in message


def test_reset_scenario():
    """Reset should clear the current scenario."""
    config.set_scenario(iam="IMAGE", ssp="SSP1", rcp="2.6")
    config.reset_scenario()
    with pytest.raises(RuntimeError, match="No scenario set"):
        config.get_scenario()


def test_valid_scenarios():
    """Valid scenarios should match paper's IAM-SSP-RCP combinations."""
    scenarios = config.VALID_SCENARIOS
    # IMAGE only has SSP1
    assert ("IMAGE", "SSP1", "2.6") in scenarios
    assert ("IMAGE", "SSP1", "4.5") in scenarios
    assert ("IMAGE", "SSP3", "2.6") not in scenarios
    # AIM only has SSP3
    assert ("AIM", "SSP3", "4.5") in scenarios


# --- the public API, now on the corrected Watanabe & Cherubini (2026) method -----------------

Row = namedtuple("Row", ["date", "amount", "flow", "activity"])
SCENARIO = ("MESSAGE", "SSP2", "4.5")


def row(year=2030, amount=1.0):
    return Row(date=pd.Timestamp(f"{year}-01-01"), amount=amount, flow=1, activity=2)


@pytest.mark.parametrize("gas", ["CH4", "N2O"])
@pytest.mark.parametrize("horizon", [20, 50, 100])
def test_agwp_and_agtp_give_the_watanabe2026_metrics(gas, horizon):
    config.set_scenario(*SCENARIO)
    agwp_gas = getattr(agwp, f"agwp_{gas.lower()}")(2030, horizon)
    agtp_gas = getattr(agtp, f"agtp_{gas.lower()}")(2030, horizon)
    assert agwp_gas / agwp.agwp_co2(2030, horizon) == pytest.approx(watanabe2026.pgwp(gas, 2030, horizon, *SCENARIO), rel=1e-12)
    assert agtp_gas / agtp.agtp_co2(2030, horizon) == pytest.approx(watanabe2026.pgtp(gas, 2030, horizon, *SCENARIO), rel=1e-12)


def test_indirect_effects_are_included():
    """Before 2.0 these were direct-only (CH4 pGWP100 ~18-20)."""
    config.set_scenario(*SCENARIO)
    assert 28 < agwp.agwp_ch4(2030) / agwp.agwp_co2(2030) < 33
    assert 230 < agwp.agwp_n2o(2030) / agwp.agwp_co2(2030) < 290


@pytest.mark.parametrize("name", ["co2", "ch4", "n2o"])
def test_characterization_functions_add_up_to_agwp(name):
    config.set_scenario(*SCENARIO)
    result = getattr(rf, f"characterize_{name}")(row(amount=3.0), 100)
    assert result.amount.shape == (100,)
    assert result.amount.sum() == pytest.approx(3.0 * getattr(agwp, f"agwp_{name}")(2030, 100), rel=1e-12)
    cumulative = getattr(rf, f"characterize_{name}")(row(amount=3.0), 100, cumulative=True).amount
    assert np.allclose(np.cumsum(result.amount), cumulative)


def test_characterize_co2_uptake():
    config.set_scenario(*SCENARIO)
    assert np.array_equal(rf.characterize_co2_uptake(row(), 50).amount, -rf.characterize_co2(row(), 50).amount)


def test_characterize_date_array():
    config.set_scenario(*SCENARIO)
    result = rf.characterize_ch4(row(2040), 30)
    assert result.date[0] == np.datetime64("2040-01-01T00:00:00")
    # years are mean Gregorian years (timedelta64[Y]), as for all functions of the package
    assert result.date[-1] == np.datetime64("2040-01-01T00:00:00") + np.timedelta64(29, "Y").astype("timedelta64[s]")
    assert len(result.date) == len(result.amount) == 30


def test_time_varying_re():
    config.set_scenario("REMIND", "SSP5", "8.5")
    fixed = rf.characterize_n2o(row(), 100).amount.sum()
    varying = rf.characterize_n2o(row(), 100, time_varying_re=True).amount.sum()
    assert varying != pytest.approx(fixed, rel=1e-3)
    series = watanabe2026.scenario_series("N2O", "agwp", "REMIND", "SSP5", "8.5", 2030, True)
    assert varying == pytest.approx(series[99], rel=1e-12)
    assert agwp.agwp_n2o(2030, 100, time_varying_re=True) == pytest.approx(series[99], rel=1e-12)


def test_time_varying_re_equals_fixed_for_constant_re():
    re = 1.3e-5
    irf = data_loader.load_irf_co2()["RCP45"]
    assert np.array_equal(
        watanabe2026.agwp_co2_series(re, irf), watanabe2026.agwp_co2_series(np.full(130, re), irf)
    )


def test_missing_re_years_take_the_last_available_value():
    """REMIND CH4 RE ends in 2100; later years must not give NaN."""
    config.set_scenario("REMIND", "SSP5", "8.5")
    assert np.isfinite(agwp.agwp_ch4(2120))
    assert np.isfinite(agwp.agwp_ch4(2030, time_varying_re=True))


def test_horizons_beyond_100_years_are_truncated():
    config.set_scenario(*SCENARIO)
    with pytest.warns(UserWarning, match="truncated"):
        value = agwp.agwp_co2(2030, 150)
    assert value == agwp.agwp_co2(2030, 100)
    result = rf.characterize_ch4(row(), 150)
    assert result.amount.shape == (100,)


@pytest.mark.parametrize("metric", ["pGWP", "pGTP"])
def test_public_functions_through_characterize(metric):
    """Custom characterization functions from prospective.radiative_forcing work with pGWP and pGTP."""
    config.set_scenario(*SCENARIO)
    df = pd.DataFrame(
        {
            "date": pd.Series(["2030-01-01"] * 3, dtype="datetime64[s]"),
            "amount": pd.Series([2.0, 2.0, 2.0], dtype="float64"),
            "flow": pd.Series([1, 2, 3], dtype="int"),
            "activity": pd.Series([1, 1, 1], dtype="int"),
        }
    )
    functions = {1: rf.characterize_co2, 2: rf.characterize_ch4, 3: rf.characterize_n2o}
    result = characterize(df, metric=metric, characterization_functions=functions, time_horizon=100)
    by_flow = result.groupby("flow").amount.sum() / 2.0
    port = watanabe2026.pgwp if metric == "pGWP" else watanabe2026.pgtp
    assert by_flow[1] == pytest.approx(1.0, rel=1e-12)
    assert by_flow[2] == pytest.approx(port("CH4", 2030, 100, *SCENARIO), rel=1e-9)
    assert by_flow[3] == pytest.approx(port("N2O", 2030, 100, *SCENARIO), rel=1e-9)


def test_emission_year_bounds_follow_data():
    """Emission years are clamped to the range covered by the RE data, not 2030-2100."""
    config.set_scenario(iam="IMAGE", ssp="SSP1", rcp="2.6")

    years = data_loader.load_re_co2("IMAGE")["_years"]
    assert (int(years[0]), int(years[-1])) == (2020, 2150)

    # Years the SI tables don't report but the data covers: no warning, distinct results
    for year in (2020, 2025, 2120):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            agwp.agwp_co2(emission_year=year, time_horizon=100)

    assert agwp.agwp_co2(emission_year=2020) != agwp.agwp_co2(emission_year=2030)

    # Outside the data: clamped with a warning, equal to the boundary year
    with pytest.warns(UserWarning, match="clamping to 2020"):
        assert agwp.agwp_co2(emission_year=2010) == agwp.agwp_co2(emission_year=2020)

    with pytest.warns(UserWarning, match="clamping to 2150"):
        assert agwp.agwp_co2(emission_year=2200) == agwp.agwp_co2(emission_year=2150)


def test_valid_scenarios_match_data_files():
    """VALID_SCENARIOS lists exactly the IAM-SSP-RCP combinations present in the data."""
    loaders = {
        "CO2": data_loader.load_re_co2,
        "CH4": data_loader.load_re_ch4,
        "N2O": data_loader.load_re_n2o,
    }
    iams = {iam for iam, _, _ in config.VALID_SCENARIOS}

    for gas, loader in loaders.items():
        in_data = {
            (iam, ssp, rcp)
            for iam in iams
            for ssp, rcps in loader(iam).items()
            if ssp != "_years"
            for rcp in rcps
        }
        assert in_data == config.VALID_SCENARIOS, (
            f"{gas} RE data and VALID_SCENARIOS disagree: "
            f"only in data {sorted(in_data - config.VALID_SCENARIOS)}, "
            f"only in config {sorted(config.VALID_SCENARIOS - in_data)}"
        )
