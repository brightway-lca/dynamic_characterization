"""
Absolute Global Warming Potential (AGWP) calculations.

Based on Barbosa Watanabe et al. (2026) equations using prospective RE and IRF.

AGWP = integral_0^t RE(t') * IRF(t') dt'

For CO2, both RE and IRF vary by scenario.
For CH4/N2O, RE varies by scenario but IRF uses fixed lifetimes.
"""

import warnings

import numpy as np

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
# All RE values in Watanabe SI files are in W/m^2/ppb (including CO2)
# 1 ppb gas = 1e-9 * (M_ATMOSPHERE / M_AIR) * M_gas in kg
M_AIR = 28.97  # g/mol (average molar mass of dry air)
M_ATMOSPHERE = 5.13252e18  # kg (total mass of atmosphere)

# Molar masses (g/mol)
M_CO2 = 44.01
M_CH4 = 16.04
M_N2O = 44.01

# Conversion constants: ppb to kg
# 1 ppb CO2 = 1e-9 * (M_ATMOSPHERE/M_AIR) * M_CO2 kg = 7.79e9 kg
# CONST converts W/m^2/ppb to W/m^2/kg: multiply by ppb/kg
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


def agwp_co2(
    emission_year: int,
    time_horizon: int = 100,
    time_varying_re: bool = False,
) -> float:
    """
    Calculate AGWP for 1 kg CO2.

    Parameters
    ----------
    emission_year : int
        Year of emission (2030-2100, clamped if outside)
    time_horizon : int
        Integration period in years
    time_varying_re : bool
        If True, use RE that evolves over the decay period.
        If False, use fixed RE from emission year (IPCC standard).

    Returns
    -------
    float
        AGWP in W*yr/m^2/kg
    """
    scenario = get_scenario()
    iam, ssp, rcp = scenario["iam"], scenario["ssp"], scenario["rcp"]

    # Load data
    re_data = load_re_co2(iam)
    irf_data = load_irf_co2()

    years = re_data["_years"]
    re_series = re_data[ssp][rcp]  # W/m^2/ppb
    irf_series = irf_data[_rcp_to_irf_key(rcp)]

    year_idx = _get_year_index(emission_year, years)

    # Limit time horizon to available IRF data
    max_years = min(time_horizon, len(irf_series))

    agwp = 0.0
    for t in range(max_years):
        irf_t = irf_series[t]

        if time_varying_re:
            # RE evolves: use RE at emission_year + t
            re_idx = min(year_idx + t, len(re_series) - 1)
            re_t = re_series[re_idx]
        else:
            # Fixed RE from emission year
            re_t = re_series[year_idx]

        # AGWP contribution: RE * IRF * conversion * dt (dt=1 year)
        agwp += re_t * irf_t * CONST_CO2

    return agwp


def agwp_ch4(
    emission_year: int,
    time_horizon: int = 100,
    time_varying_re: bool = False,
) -> float:
    """
    Calculate AGWP for 1 kg CH4.

    Parameters
    ----------
    emission_year : int
        Year of emission (2030-2100, clamped if outside)
    time_horizon : int
        Integration period in years
    time_varying_re : bool
        If True, use RE that evolves over the decay period.

    Returns
    -------
    float
        AGWP in W*yr/m^2/kg
    """
    scenario = get_scenario()
    iam, ssp, rcp = scenario["iam"], scenario["ssp"], scenario["rcp"]

    # Load data
    re_data = load_re_ch4(iam)
    irf_series = load_irf_ch4()

    years = re_data["_years"]
    re_series = re_data[ssp][rcp]  # W/m^2/ppb

    year_idx = _get_year_index(emission_year, years)

    max_years = min(time_horizon, len(irf_series))

    agwp = 0.0
    for t in range(max_years):
        irf_t = irf_series[t]

        if time_varying_re:
            re_idx = min(year_idx + t, len(re_series) - 1)
            re_t = re_series[re_idx]
        else:
            re_t = re_series[year_idx]

        agwp += re_t * irf_t * CONST_CH4

    return agwp


def agwp_n2o(
    emission_year: int,
    time_horizon: int = 100,
    time_varying_re: bool = False,
) -> float:
    """
    Calculate AGWP for 1 kg N2O.

    Parameters
    ----------
    emission_year : int
        Year of emission (2030-2100, clamped if outside)
    time_horizon : int
        Integration period in years
    time_varying_re : bool
        If True, use RE that evolves over the decay period.

    Returns
    -------
    float
        AGWP in W*yr/m^2/kg
    """
    scenario = get_scenario()
    iam, ssp, rcp = scenario["iam"], scenario["ssp"], scenario["rcp"]

    # Load data
    re_data = load_re_n2o(iam)
    irf_series = load_irf_n2o()

    years = re_data["_years"]
    re_series = re_data[ssp][rcp]  # W/m^2/ppb

    year_idx = _get_year_index(emission_year, years)

    max_years = min(time_horizon, len(irf_series))

    agwp = 0.0
    for t in range(max_years):
        irf_t = irf_series[t]

        if time_varying_re:
            re_idx = min(year_idx + t, len(re_series) - 1)
            re_t = re_series[re_idx]
        else:
            re_t = re_series[year_idx]

        agwp += re_t * irf_t * CONST_N2O

    return agwp
