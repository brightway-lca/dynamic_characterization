"""The prospective scenario can be scoped to one call instead of a session.

Named `test_characterize_scenario_scoping` rather than
`test_prospective_scenario_scoping` so it sorts alphabetically before
`test_prospective.py`, which replaces `dynamic_characterization` and
`dynamic_characterization.prospective` in `sys.modules` with stubs at import
time (see the note at the top of `test_characterize_prospective_fallback.py`).
A module that imports the real package must be collected before that happens.
"""

import warnings

import numpy as np
import pandas as pd
import pytest

from dynamic_characterization import characterize
from dynamic_characterization.prospective import (
    get_scenario,
    reset_scenario,
    scenario_context,
    set_scenario,
)
from dynamic_characterization.prospective.radiative_forcing import (
    characterize_ch4 as prospective_characterize_ch4,
    characterize_co2 as prospective_characterize_co2,
    _get_year_index,
    _reset_bound_warnings,
)

IMAGE = {"iam": "IMAGE", "ssp": "SSP1", "rcp": "2.6"}
MESSAGE = {"iam": "MESSAGE", "ssp": "SSP2", "rcp": "4.5"}

# Flow ids, matching the pattern used in test_characterize_prospective_fallback.py:
# these are plain ints standing in for biosphere flow ids, paired with explicit
# `characterization_functions` dicts so the tests don't need a real bw2data
# project, biosphere database, or LCIA method.
CO2_FLOW = 1
CH4_FLOW = 2


@pytest.fixture(autouse=True)
def clean_scenario():
    reset_scenario()
    _reset_bound_warnings()
    yield
    reset_scenario()
    _reset_bound_warnings()


def test_context_applies_scenario_inside_block():
    with scenario_context(IMAGE):
        assert get_scenario() == IMAGE


def test_context_restores_previous_scenario():
    set_scenario(**MESSAGE)
    with scenario_context(IMAGE):
        assert get_scenario() == IMAGE
    assert get_scenario() == MESSAGE


def test_context_restores_unset_state():
    with scenario_context(IMAGE):
        pass
    with pytest.raises(RuntimeError):
        get_scenario()


def test_context_restores_on_exception():
    set_scenario(**MESSAGE)
    with pytest.raises(ValueError):
        with scenario_context(IMAGE):
            raise ValueError("boom")
    assert get_scenario() == MESSAGE


def test_context_with_none_is_a_no_op():
    set_scenario(**MESSAGE)
    with scenario_context(None):
        assert get_scenario() == MESSAGE
    assert get_scenario() == MESSAGE


def test_context_validates_the_scenario():
    with pytest.raises(ValueError):
        with scenario_context({"iam": "REMIND", "ssp": "SSP2", "rcp": "2.6"}):
            pass


def _inventory(flow_id):
    return pd.DataFrame(
        {
            "date": pd.to_datetime(["2030-01-01"]),
            "amount": [1.0],
            "flow": [flow_id],
            "activity": [0],
        }
    )


def _co2_functions():
    return {CO2_FLOW: prospective_characterize_co2}


def _ch4_functions():
    return {CH4_FLOW: prospective_characterize_ch4}


def test_characterize_scenario_argument_matches_session_scenario():
    set_scenario(**IMAGE)
    from_session = characterize(
        _inventory(CO2_FLOW),
        metric="pGWP",
        characterization_functions=_co2_functions(),
    )
    reset_scenario()

    from_argument = characterize(
        _inventory(CO2_FLOW),
        metric="pGWP",
        characterization_functions=_co2_functions(),
        scenario=IMAGE,
    )

    pd.testing.assert_frame_equal(from_session, from_argument)


def test_characterize_scenario_argument_does_not_leak():
    set_scenario(**MESSAGE)
    characterize(
        _inventory(CO2_FLOW),
        metric="pGWP",
        characterization_functions=_co2_functions(),
        scenario=IMAGE,
    )
    assert get_scenario() == MESSAGE


def test_two_scenarios_in_one_process_differ():
    low = characterize(
        _inventory(CH4_FLOW),
        metric="prospective_radiative_forcing",
        characterization_functions=_ch4_functions(),
        scenario=IMAGE,
    )
    high = characterize(
        _inventory(CH4_FLOW),
        metric="prospective_radiative_forcing",
        characterization_functions=_ch4_functions(),
        scenario={"iam": "IMAGE", "ssp": "SSP1", "rcp": "8.5"},
    )
    assert low["amount"].sum() != high["amount"].sum()


def test_pgtp_identity_dispatch_survives_scoping():
    """CO2's pGTP is 1.0 by definition - proves the AGTP branch still matched."""
    result = characterize(
        _inventory(CO2_FLOW),
        metric="pGTP",
        characterization_functions=_co2_functions(),
        scenario=IMAGE,
    )
    assert result["amount"].sum() == pytest.approx(1.0, rel=1e-6)


def test_prospective_metric_without_any_scenario_still_raises():
    with pytest.raises(RuntimeError):
        characterize(
            _inventory(CO2_FLOW),
            metric="pGWP",
            characterization_functions=_co2_functions(),
        )


def test_out_of_bounds_year_warns_once_per_bound():
    years = np.arange(2020, 2151)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        for _ in range(50):
            _get_year_index(1990, years)
    assert len(caught) == 1
    assert "clamping to 2020" in str(caught[0].message)


def test_both_bounds_warn_separately():
    years = np.arange(2020, 2151)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        for _ in range(10):
            _get_year_index(1990, years)
            _get_year_index(2200, years)
    assert len(caught) == 2


def test_characterize_warns_on_each_call_not_just_the_first():
    """Regression test: `_WARNED_BOUNDS` used to be module-global and never
    reset between `characterize()` calls, so in a multi-call sweep (e.g. one
    `characterize()` per row of a `TimexLCA.compare()`) only the first call
    ever warned about a clamped out-of-bounds emission year - later calls
    with the same out-of-bounds year clamped silently. `characterize()` now
    resets the dedup set at the start of its scoped body (see
    `_reset_bound_warnings`'s call site in
    `dynamic_characterization.dynamic_characterization.characterize`), so
    the dedup scope is "once per bound per call", not "once per bound per
    process".
    """
    early_inventory = pd.DataFrame(
        {
            "date": pd.to_datetime(["1990-01-01"]),
            "amount": [1.0],
            "flow": [CO2_FLOW],
            "activity": [0],
        }
    )
    for _ in range(2):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            characterize(
                early_inventory,
                metric="prospective_radiative_forcing",
                characterization_functions=_co2_functions(),
                scenario=IMAGE,
            )
        clamp_warnings = [w for w in caught if "clamping" in str(w.message)]
        assert len(clamp_warnings) == 1
