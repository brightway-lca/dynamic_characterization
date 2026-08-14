"""Prospective metrics with IPCC fallback functions for non-Watanabe GHGs.

Note: test_prospective.py replaces `dynamic_characterization` in sys.modules with a
stub at import time, which breaks importing the real package afterwards. This module
therefore has to sort before it alphabetically.
"""

import numpy as np
import pandas as pd
import pytest

from dynamic_characterization import characterize
from dynamic_characterization.ipcc_ar6.radiative_forcing import (
    characterize_co,
    create_generic_characterization_function,
)
from dynamic_characterization.prospective import set_scenario
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
