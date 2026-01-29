"""
Watanabe prospective characterization functions.

Returns radiative forcing time series matching the existing ipcc_ar6 API,
but using scenario-based radiative efficiencies from Barbosa Watanabe et al. (2026).

Reference: https://doi.org/10.1021/acs.est.5b01118
"""

import warnings

import numpy as np

from dynamic_characterization.classes import CharacterizedRow

from .config import get_scenario
from .data_loader import (
    load_irf_ch4,
    load_irf_co2,
    load_irf_n2o,
    load_re_ch4,
    load_re_co2,
    load_re_n2o,
)

# Constants for unit conversion
# Convert RE from W/m^2/ppb to W/m^2/kg
M_AIR = 28.97  # g/mol (average molar mass of dry air)
M_ATMOSPHERE = 5.13252e18  # kg (total mass of atmosphere)

# Molar masses (g/mol)
M_CO2 = 44.01
M_CH4 = 16.04
M_N2O = 44.01

# Conversion constants: ppb to kg
PPB_TO_KG_CO2 = 1e-9 * M_ATMOSPHERE / M_AIR * M_CO2  # kg per ppb
PPB_TO_KG_CH4 = 1e-9 * M_ATMOSPHERE / M_AIR * M_CH4  # kg per ppb
PPB_TO_KG_N2O = 1e-9 * M_ATMOSPHERE / M_AIR * M_N2O  # kg per ppb

CONST_CO2 = 1.0 / PPB_TO_KG_CO2  # ppb/kg (to convert RE from per-ppb to per-kg)
CONST_CH4 = 1.0 / PPB_TO_KG_CH4  # ppb/kg
CONST_N2O = 1.0 / PPB_TO_KG_N2O  # ppb/kg


def _get_year_index(emission_year: int, years: np.ndarray) -> int:
    """
    Get index for emission year in RE data, with clamping and warnings.

    RE data spans 2020-2150. Years outside 2030-2100 are clamped.
    """
    min_year, max_year = 2030, 2100

    if emission_year < min_year:
        warnings.warn(
            f"Emission year {emission_year} < {min_year}, clamping to {min_year}"
        )
        emission_year = min_year
    elif emission_year > max_year:
        warnings.warn(
            f"Emission year {emission_year} > {max_year}, clamping to {max_year}"
        )
        emission_year = max_year

    # Find index in years array
    idx = np.searchsorted(years, emission_year)
    return idx


def _rcp_to_irf_key(rcp: str) -> str:
    """Convert RCP string (e.g., '2.6') to IRF dict key (e.g., 'RCP26')."""
    return f"RCP{rcp.replace('.', '')}"


def characterize_co2(
    series,
    period: int = 100,
    cumulative: bool = False,
    time_varying_re: bool = False,
) -> CharacterizedRow:
    """
    Calculate radiative forcing time series for 1 kg CO2 emission.

    Uses scenario-based radiative efficiencies from Barbosa Watanabe et al. (2026).
    Scenario must be set via prospective.set_scenario() before calling.

    Parameters
    ----------
    series : namedtuple
        Row from dynamic inventory with date, amount, flow, activity
    period : int
        Time horizon in years (default: 100)
    cumulative : bool
        If True, return cumulative radiative forcing;
        If False, return marginal (yearly) forcing (default)
    time_varying_re : bool
        If True, use RE that evolves over the decay period.
        If False, use fixed RE from emission year (IPCC standard, default).

    Returns
    -------
    CharacterizedRow
        namedtuple with date, amount, flow, activity arrays.
        Amount is in W*yr/m^2/kg CO2.
    """
    scenario = get_scenario()
    iam, ssp, rcp = scenario["iam"], scenario["ssp"], scenario["rcp"]

    # Load data
    re_data = load_re_co2(iam)
    irf_data = load_irf_co2()

    years = re_data["_years"]
    re_series = re_data[ssp][rcp]  # W/m^2/ppb
    irf_series = irf_data[_rcp_to_irf_key(rcp)]

    # Get emission year from series date
    date_beginning = series.date.to_numpy()
    emission_year = int(str(date_beginning)[:4])
    year_idx = _get_year_index(emission_year, years)

    # Limit time horizon to available IRF data
    max_years = min(period, len(irf_series))

    # Create date array
    dates_characterized = date_beginning + np.arange(
        start=0, stop=max_years, dtype="timedelta64[Y]"
    ).astype("timedelta64[s]")

    # Calculate cumulative radiative forcing at each time step
    # AGWP(t) = integral_0^t RE(t') * IRF(t') dt'
    # Start from t=1 so forcing[0]=0 (no time elapsed = no forcing yet)
    forcing = np.zeros(max_years, dtype="float64")
    cumulative_forcing = 0.0

    for t in range(1, max_years):
        irf_t = irf_series[t]

        if time_varying_re:
            # RE evolves: use RE at emission_year + t
            re_idx = min(year_idx + t, len(re_series) - 1)
            re_t = re_series[re_idx]
        else:
            # Fixed RE from emission year
            re_t = re_series[year_idx]

        # Contribution: RE * IRF * conversion * dt (dt=1 year)
        delta_forcing = re_t * irf_t * CONST_CO2
        cumulative_forcing += delta_forcing
        forcing[t] = cumulative_forcing

    # Scale by emission amount
    forcing = forcing * series.amount

    if not cumulative:
        # Convert to marginal (yearly) forcing
        forcing = np.diff(forcing, prepend=0)

    return CharacterizedRow(
        date=np.array(dates_characterized, dtype="datetime64[s]"),
        amount=forcing,
        flow=series.flow,
        activity=series.activity,
    )


def characterize_co2_uptake(
    series,
    period: int = 100,
    cumulative: bool = False,
    time_varying_re: bool = False,
) -> CharacterizedRow:
    """
    Calculate radiative forcing time series for CO2 uptake (negative forcing).

    Same as characterize_co2 but with negative sign for uptake.
    """
    result = characterize_co2(series, period, cumulative, time_varying_re)
    return CharacterizedRow(
        date=result.date,
        amount=-result.amount,
        flow=result.flow,
        activity=result.activity,
    )


def characterize_ch4(
    series,
    period: int = 100,
    cumulative: bool = False,
    time_varying_re: bool = False,
) -> CharacterizedRow:
    """
    Calculate radiative forcing time series for 1 kg CH4 emission.

    Uses scenario-based radiative efficiencies from Barbosa Watanabe et al. (2026).
    Scenario must be set via prospective.set_scenario() before calling.

    This includes only direct effects. For indirect effects (ozone, water vapor),
    multiply the result by the appropriate factor (~1.43 for GWP100).

    Parameters
    ----------
    series : namedtuple
        Row from dynamic inventory with date, amount, flow, activity
    period : int
        Time horizon in years (default: 100)
    cumulative : bool
        If True, return cumulative radiative forcing;
        If False, return marginal (yearly) forcing (default)
    time_varying_re : bool
        If True, use RE that evolves over the decay period.
        If False, use fixed RE from emission year (IPCC standard, default).

    Returns
    -------
    CharacterizedRow
        namedtuple with date, amount, flow, activity arrays.
        Amount is in W*yr/m^2/kg CH4.
    """
    scenario = get_scenario()
    iam, ssp, rcp = scenario["iam"], scenario["ssp"], scenario["rcp"]

    # Load data
    re_data = load_re_ch4(iam)
    irf_series = load_irf_ch4()

    years = re_data["_years"]
    re_series = re_data[ssp][rcp]  # W/m^2/ppb

    # Get emission year from series date
    date_beginning = series.date.to_numpy()
    emission_year = int(str(date_beginning)[:4])
    year_idx = _get_year_index(emission_year, years)

    # Limit time horizon to available IRF data
    max_years = min(period, len(irf_series))

    # Create date array
    dates_characterized = date_beginning + np.arange(
        start=0, stop=max_years, dtype="timedelta64[Y]"
    ).astype("timedelta64[s]")

    # Calculate cumulative radiative forcing at each time step
    # Start from t=1 so forcing[0]=0 (no time elapsed = no forcing yet)
    forcing = np.zeros(max_years, dtype="float64")
    cumulative_forcing = 0.0

    for t in range(1, max_years):
        irf_t = irf_series[t]

        if time_varying_re:
            re_idx = min(year_idx + t, len(re_series) - 1)
            re_t = re_series[re_idx]
        else:
            re_t = re_series[year_idx]

        delta_forcing = re_t * irf_t * CONST_CH4
        cumulative_forcing += delta_forcing
        forcing[t] = cumulative_forcing

    # Scale by emission amount
    forcing = forcing * series.amount

    if not cumulative:
        forcing = np.diff(forcing, prepend=0)

    return CharacterizedRow(
        date=np.array(dates_characterized, dtype="datetime64[s]"),
        amount=forcing,
        flow=series.flow,
        activity=series.activity,
    )


def characterize_n2o(
    series,
    period: int = 100,
    cumulative: bool = False,
    time_varying_re: bool = False,
) -> CharacterizedRow:
    """
    Calculate radiative forcing time series for 1 kg N2O emission.

    Uses scenario-based radiative efficiencies from Barbosa Watanabe et al. (2026).
    Scenario must be set via prospective.set_scenario() before calling.

    Parameters
    ----------
    series : namedtuple
        Row from dynamic inventory with date, amount, flow, activity
    period : int
        Time horizon in years (default: 100)
    cumulative : bool
        If True, return cumulative radiative forcing;
        If False, return marginal (yearly) forcing (default)
    time_varying_re : bool
        If True, use RE that evolves over the decay period.
        If False, use fixed RE from emission year (IPCC standard, default).

    Returns
    -------
    CharacterizedRow
        namedtuple with date, amount, flow, activity arrays.
        Amount is in W*yr/m^2/kg N2O.
    """
    scenario = get_scenario()
    iam, ssp, rcp = scenario["iam"], scenario["ssp"], scenario["rcp"]

    # Load data
    re_data = load_re_n2o(iam)
    irf_series = load_irf_n2o()

    years = re_data["_years"]
    re_series = re_data[ssp][rcp]  # W/m^2/ppb

    # Get emission year from series date
    date_beginning = series.date.to_numpy()
    emission_year = int(str(date_beginning)[:4])
    year_idx = _get_year_index(emission_year, years)

    # Limit time horizon to available IRF data
    max_years = min(period, len(irf_series))

    # Create date array
    dates_characterized = date_beginning + np.arange(
        start=0, stop=max_years, dtype="timedelta64[Y]"
    ).astype("timedelta64[s]")

    # Calculate cumulative radiative forcing at each time step
    # Start from t=1 so forcing[0]=0 (no time elapsed = no forcing yet)
    forcing = np.zeros(max_years, dtype="float64")
    cumulative_forcing = 0.0

    for t in range(1, max_years):
        irf_t = irf_series[t]

        if time_varying_re:
            re_idx = min(year_idx + t, len(re_series) - 1)
            re_t = re_series[re_idx]
        else:
            re_t = re_series[year_idx]

        delta_forcing = re_t * irf_t * CONST_N2O
        cumulative_forcing += delta_forcing
        forcing[t] = cumulative_forcing

    # Scale by emission amount
    forcing = forcing * series.amount

    if not cumulative:
        forcing = np.diff(forcing, prepend=0)

    return CharacterizedRow(
        date=np.array(dates_characterized, dtype="datetime64[s]"),
        amount=forcing,
        flow=series.flow,
        activity=series.activity,
    )
