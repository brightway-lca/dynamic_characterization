"""
Prospective absolute global warming potentials (AGWP) of CO2, CH4 and N2O, in W m-2 yr kg-1.

Since version 2.0 these are computed with :mod:`.watanabe2026` (Barbosa Watanabe & Cherubini
2026, with two errors of the published code corrected), for the scenario of
``prospective.set_scenario()``. CH4 and N2O include their indirect effects and a
climate-carbon feedback. The signatures are unchanged; horizons are limited to 100 years.
"""

import warnings

from ._common import (  # noqa: F401  (re-exported for backwards compatibility)
    CONST_CH4,
    CONST_CO2,
    CONST_N2O,
    M_AIR,
    M_ATMOSPHERE,
    M_CH4,
    M_CO2,
    M_N2O,
    PPB_TO_KG_CH4,
    PPB_TO_KG_CO2,
    PPB_TO_KG_N2O,
    _get_year_index,
    _rcp_to_irf_key,
)
from .config import get_scenario
from . import watanabe2026

MAX_YEARS = 100


def _value(gas: str, metric: str, emission_year: int, time_horizon: int, time_varying_re: bool) -> float:
    if time_horizon <= 0:
        return 0.0
    if time_horizon > MAX_YEARS:
        warnings.warn(f"Time horizon {time_horizon} > {MAX_YEARS} years, truncated to {MAX_YEARS}.")
        time_horizon = MAX_YEARS
    scenario = get_scenario()
    series = watanabe2026.scenario_series(
        gas, metric, scenario["iam"], scenario["ssp"], scenario["rcp"], int(emission_year), bool(time_varying_re)
    )
    return float(series[int(time_horizon) - 1])


def agwp_co2(emission_year: int, time_horizon: int = 100, time_varying_re: bool = False) -> float:
    """AGWP of 1 kg CO2 emitted in ``emission_year`` (2020-2150, clamped outside)."""
    return _value("CO2", "agwp", emission_year, time_horizon, time_varying_re)


def agwp_ch4(emission_year: int, time_horizon: int = 100, time_varying_re: bool = False) -> float:
    """AGWP of 1 kg CH4, incl. ozone, stratospheric water vapour and climate-carbon feedback."""
    return _value("CH4", "agwp", emission_year, time_horizon, time_varying_re)


def agwp_n2o(emission_year: int, time_horizon: int = 100, time_varying_re: bool = False) -> float:
    """AGWP of 1 kg N2O, incl. ozone, its CH4 effect and climate-carbon feedback."""
    return _value("N2O", "agwp", emission_year, time_horizon, time_varying_re)
