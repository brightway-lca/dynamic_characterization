"""Prospective metrics with IPCC fallback functions for non-Watanabe GHGs.

Note: test_prospective.py replaces `dynamic_characterization` in sys.modules with a
stub at import time, which breaks importing the real package afterwards. This module
therefore has to sort before it alphabetically.
"""

import json
import os

import numpy as np
import pandas as pd
import pytest

from dynamic_characterization import characterize
from dynamic_characterization.ipcc_ar6.radiative_forcing import (
    characterize_co,
    characterize_co2,
    create_generic_characterization_function,
)
from dynamic_characterization.prospective import agtp, set_scenario
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
    # A generic function built from a decay series, as done for GHGs that are
    # only in decay_multipliers.json (here: a stand-in decay series).
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


# --- pGTP of IPCC fallback functions -------------------------------------------------

DECAY_MULTIPLIERS = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "dynamic_characterization",
    "ipcc_ar6",
    "data",
    "decay_multipliers.json",
)


def _halocarbon_function(cas_key):
    with open(DECAY_MULTIPLIERS) as f:
        return create_generic_characterization_function(np.array(json.load(f)[cas_key]))


def _pgtp(functions, amount=1.0, time_horizon=100):
    set_scenario(iam="MESSAGE", ssp="SSP2", rcp="4.5")
    flows = list(functions)
    df = pd.DataFrame(
        {
            "date": pd.Series(["2030-01-01"] * len(flows), dtype="datetime64[s]"),
            "amount": pd.Series([amount] * len(flows), dtype="float64"),
            "flow": pd.Series(flows, dtype="int"),
            "activity": pd.Series([1] * len(flows), dtype="int"),
        }
    )
    result = characterize(
        df, metric="pGTP", characterization_functions=functions, time_horizon=time_horizon
    )
    return {flow: result.amount[result.flow == flow].sum() for flow in flows}


def test_pgtp_fallback_is_linear_in_amount():
    """The amount used to be applied twice, so doubling the emission quadrupled the pGTP."""
    functions = {1: characterize_co, 2: _halocarbon_function("000075-73-0")}
    single = _pgtp(functions, amount=1.0)
    double = _pgtp(functions, amount=2.0)
    for flow in functions:
        assert double[flow] == pytest.approx(2 * single[flow], rel=1e-12)


@pytest.mark.parametrize("time_horizon", [20, 100])
def test_pgtp_fallback_co_is_scaled_co2(time_horizon):
    """characterize_co is the IPCC CO2 response times M_CO2/M_CO, so its pGTP is too."""
    result = _pgtp({1: characterize_co, 2: characterize_co2}, time_horizon=time_horizon)
    assert result[1] / result[2] == pytest.approx(44.01 / 28.01, rel=1e-12)


@pytest.mark.parametrize(
    "cas_key, ar6_gtp20",
    [
        ("000075-73-0", 6346),  # CF4, Hodnebrog et al. 2020 / AR6 Table 7.SM.7
        ("000075-69-4", 7931),  # CFC-11
    ],
)
def test_pgtp20_fallback_close_to_ar6_gtp20(cas_key, ar6_gtp20):
    """
    At H = 20 the scenario CO2 response barely differs from AR6, so the fallback pGTP must be
    close to the published GTP20. (At H = 100 it is ~25 % lower: the scenario CO2 response and
    the missing climate-carbon feedback of the legacy functions matter there.)
    """
    result = _pgtp({1: _halocarbon_function(cas_key)}, time_horizon=20)
    assert result[1] == pytest.approx(ar6_gtp20, rel=0.05)


def test_pgtp_fallback_uses_own_agtp_if_available():
    def characterize_with_agtp(series, period=100, cumulative=False):
        return characterize_co(series, period, cumulative=cumulative)

    own_agtp = 1e-15  # K/kg
    characterize_with_agtp.agtp = lambda time_horizon: own_agtp
    result = _pgtp({1: characterize_with_agtp}, time_horizon=100)
    expected = own_agtp / agtp.agtp_co2(emission_year=2030, time_horizon=100)
    assert result[1] == pytest.approx(expected, rel=1e-12)
