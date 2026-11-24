"""
Prospective radiative forcing of CO2, CH4 and N2O (Barbosa Watanabe & Cherubini 2026).

Since version 2.0 these functions are thin wrappers around :mod:`.watanabe2026` (the paper's
method, with two errors of its published code corrected): CH4 and N2O include their indirect
effects (ozone, stratospheric water vapour, the CH4 effect of N2O) and a climate-carbon
feedback. The signatures are unchanged.

The scenario must be set with ``prospective.set_scenario()``. As for all characterization
functions of this package, the value at ``date + k`` years is the forcing integrated over the
year starting then (W yr m-2 per kg), so the annual values add up to the AGWP of the period.
Periods beyond 100 years are truncated to 100 (the length of the scenario CO2 impulse response).

Reference: Barbosa Watanabe, M.D., Cherubini, F. (2026), Environ. Sci. Technol. 60, 3202-3215,
https://doi.org/10.1021/acs.est.5c12391
"""

from dynamic_characterization.classes import CharacterizedRow

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


def _function(gas: str, uptake: bool = False):
    from dynamic_characterization import ar6_prospective

    if gas == "CO2":
        return ar6_prospective.co2_function("watanabe2026", uptake=uptake)
    return ar6_prospective.paper_function(gas)


def characterize_co2(series, period: int = 100, cumulative: bool = False, time_varying_re: bool = False) -> CharacterizedRow:
    """
    Radiative forcing of a CO2 emission (``series``: a row with date, amount, flow, activity),
    with the scenario CO2 radiative efficiency and impulse response.

    Parameters
    ----------
    series : namedtuple
        Row of the dynamic inventory.
    period : int
        Number of years (at most 100).
    cumulative : bool
        Cumulative instead of annual forcing.
    time_varying_re : bool
        Let the radiative efficiency follow the scenario over the response; by default that of
        the emission year is used throughout, as in the paper.
    """
    return _function("CO2")(series, period, cumulative=cumulative, time_varying_re=time_varying_re)


def characterize_co2_uptake(series, period: int = 100, cumulative: bool = False, time_varying_re: bool = False) -> CharacterizedRow:
    """As :func:`characterize_co2`, with a negative sign for CO2 uptake."""
    return _function("CO2", uptake=True)(series, period, cumulative=cumulative, time_varying_re=time_varying_re)


def characterize_ch4(series, period: int = 100, cumulative: bool = False, time_varying_re: bool = False) -> CharacterizedRow:
    """
    Radiative forcing of a CH4 emission, incl. ozone, stratospheric water vapour and the
    climate-carbon feedback. Parameters as for :func:`characterize_co2`.
    """
    return _function("CH4")(series, period, cumulative=cumulative, time_varying_re=time_varying_re)


def characterize_n2o(series, period: int = 100, cumulative: bool = False, time_varying_re: bool = False) -> CharacterizedRow:
    """
    Radiative forcing of an N2O emission, incl. ozone, its effect on CH4 and the climate-carbon
    feedback. Parameters as for :func:`characterize_co2`.
    """
    return _function("N2O")(series, period, cumulative=cumulative, time_varying_re=time_varying_re)


def _delegate_metrics(wrapper, gas: str, uptake: bool = False) -> None:
    """Give a wrapper the ``agwp_for`` / ``agtp_for`` of the function it wraps (used by pGTP)."""
    wrapper.agwp_for = lambda date, time_horizon, time_varying_re=False: _function(gas, uptake).agwp_for(
        date, time_horizon, time_varying_re=time_varying_re
    )
    wrapper.agtp_for = lambda date, time_horizon, time_varying_re=False: _function(gas, uptake).agtp_for(
        date, time_horizon, time_varying_re=time_varying_re
    )


_delegate_metrics(characterize_co2, "CO2")
_delegate_metrics(characterize_co2_uptake, "CO2", uptake=True)
_delegate_metrics(characterize_ch4, "CH4")
_delegate_metrics(characterize_n2o, "N2O")
