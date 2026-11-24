"""Prospective metrics with custom characterization functions (IPCC AR6 or user-made)."""

import numpy as np
import pandas as pd
import pytest

from dynamic_characterization import characterize
from dynamic_characterization import ar6_prospective
from dynamic_characterization.ipcc_ar6 import gases
from dynamic_characterization.ipcc_ar6.radiative_forcing import (
    characterize_co,
    create_generic_characterization_function,
    make_characterization_function,
)
from dynamic_characterization.prospective import set_scenario, watanabe2026
from dynamic_characterization.prospective.radiative_forcing import (
    characterize_ch4 as prospective_characterize_ch4,
)

# Flow ids: 1 is covered by the Watanabe module, 2 and 3 only by IPCC AR6.
CH4_FLOW = 1
CO_FLOW = 2
GENERIC_FLOW = 3


def dynamic_inventory() -> pd.DataFrame:
    return pd.DataFrame(
        data={
            "date": pd.Series(
                data=["2030-01-01", "2030-01-01", "2030-01-01"],
                dtype="datetime64[s]",
            ),
            "amount": pd.Series(data=[1.0, 1.0, 1.0], dtype="float64"),
            "flow": pd.Series(data=[CH4_FLOW, CO_FLOW, GENERIC_FLOW], dtype="int"),
            "activity": pd.Series(data=[10, 10, 10], dtype="int"),
        }
    )


def mixed_characterization_functions() -> dict:
    # a user-made function from a (stand-in) cumulative decay series
    decay_series = np.linspace(0, 1e-13, 100)
    return {
        CH4_FLOW: prospective_characterize_ch4,
        CO_FLOW: characterize_co,
        GENERIC_FLOW: create_generic_characterization_function(decay_series),
    }


@pytest.mark.parametrize("time_varying_re", [False, True])
@pytest.mark.parametrize("metric", ["pGWP", "pGTP", "prospective_radiative_forcing"])
def test_prospective_metrics_with_ipcc_fallback_functions(metric, time_varying_re):
    """IPCC fallback functions don't take `time_varying_re` and must not get it."""
    set_scenario(iam="IMAGE", ssp="SSP1", rcp="2.6")

    characterized = characterize(
        dynamic_inventory(),
        metric=metric,
        characterization_functions=mixed_characterization_functions(),
        time_horizon=100,
        time_varying_re=time_varying_re,
    )

    # All three flows are characterized, none silently dropped.
    assert set(characterized.flow.unique()) == {CH4_FLOW, CO_FLOW, GENERIC_FLOW}


# --- pGTP of custom functions --------------------------------------------------------------

SCENARIO = ("MESSAGE", "SSP2", "4.5")
DATE = pd.Timestamp("2030-01-01")


def _pgtp(functions, amount=1.0, time_horizon=100):
    set_scenario(*SCENARIO)
    flows = list(functions)
    df = pd.DataFrame(
        {
            "date": pd.Series(["2030-01-01"] * len(flows), dtype="datetime64[s]"),
            "amount": pd.Series([amount] * len(flows), dtype="float64"),
            "flow": pd.Series(flows, dtype="int"),
            "activity": pd.Series([1] * len(flows), dtype="int"),
        }
    )
    result = characterize(df, metric="pGTP", characterization_functions=functions, time_horizon=time_horizon)
    return {flow: result.amount[result.flow == flow].sum() for flow in flows}


def _paper_co2_as_custom_function():
    """A user-made function with the forcing of the scenario CO2 (cumulative AGWP per year)."""
    series = watanabe2026.scenario_series("CO2", "agwp", *SCENARIO, 2030)
    return create_generic_characterization_function(np.asarray(series))


def test_pgtp_is_linear_in_amount():
    functions = {1: characterize_co, 2: _paper_co2_as_custom_function()}
    single = _pgtp(functions, amount=1.0)
    double = _pgtp(functions, amount=2.0)
    for flow in functions:
        assert double[flow] == pytest.approx(2 * single[flow], rel=1e-12)


@pytest.mark.parametrize("time_horizon", [20, 50, 100])
def test_custom_function_uses_the_same_temperature_model_as_the_reference(time_horizon):
    """A custom function with the scenario CO2 forcing gets pGTP ~ 1."""
    result = _pgtp({1: _paper_co2_as_custom_function()}, time_horizon=time_horizon)
    assert result[1] == pytest.approx(1.0, rel=0.02)


def test_pgtp20_of_an_ar6_halocarbon_is_close_to_its_gtp20():
    """
    AR6 functions passed as custom functions use their own AGTP against the scenario CO2. For CF4
    at H = 20 that is close to the GTP20 of Hodnebrog et al. 2020 (6346, hodnebrog20.csv; AR6
    Table 7.SM.7 gives no GTP20).
    """
    result = _pgtp({1: make_characterization_function(gases.halocarbon("75-73-0"))}, time_horizon=20)
    assert result[1] == pytest.approx(6346, rel=0.10)


def test_pgtp_uses_own_agtp_if_available():
    def characterize_with_agtp(series, period=100, cumulative=False):
        return characterize_co(series, period, cumulative=cumulative)

    own_agtp = 1e-15  # K/kg
    characterize_with_agtp.agtp = lambda time_horizon: own_agtp
    result = _pgtp({1: characterize_with_agtp}, time_horizon=100)
    expected = own_agtp / ar6_prospective.co2_function("watanabe2026").agtp_for(DATE, 100)
    assert result[1] == pytest.approx(expected, rel=1e-12)
