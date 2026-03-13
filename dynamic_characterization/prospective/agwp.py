"""
Absolute Global Warming Potential (AGWP) calculations.

Based on Watanabe et al. (2026) equations using prospective RE and IRF.

AGWP = integral_0^t RE(t') * IRF(t') dt'

For CO2, both RE and IRF vary by scenario.
For CH4/N2O, RE varies by scenario but IRF uses fixed lifetimes.
"""

import warnings

import numpy as np

from dynamic_characterization.prospective.config import get_scenario
from dynamic_characterization.prospective.data_loader import (
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

CH4_INDIRECT_RE = 0.00014 + 0.00004


def _check_indirect_effects_mode(mode: str) -> str:
    """Validate the selected mode for indirect effects."""
    valid_modes = {"none", "no_carbon_cycle", "all"}
    if mode not in valid_modes:
        raise ValueError(
            f"Invalid indirect effect mode '{mode}'. Valid modes: {sorted(valid_modes)}"
        )
    return mode


def _get_re_value(
    re_series: np.ndarray,
    year_idx: int,
    step: int,
    time_varying_re: bool,
) -> float:
    """Configuration-aware radiate efficiency getter."""
    if time_varying_re:
        # RE evolves: use RE at emission_year + t
        re_idx = min(year_idx + step, len(re_series) - 1)
        re_t = re_series[re_idx]
    else:
        re_t = re_series[year_idx]
    return re_t


def _agwp_ch4_cumulative(
    emission_year: int,
    time_horizon: int = 100,
    time_varying_re: bool = False,
    indirect_effects: str = "all",
) -> np.ndarray:
    """
    Return cumulative CH4 AGWP values for each year of the selected horizon.

    Based on Watanabe's Code #4 CH4 AGWP function.

    Parameters
    ----------
    emission_year : int
        Year of emission (2030-2100, clamped if outside)
    time_horizon : int
        Integration period in years
    time_varying_re : bool
        If True, use RE that evolves over the decay period.
    indirect_effects : str
        Indirect effects to include. Defaults to "all".

    Returns
    -------
    np.ndarray
        Cumulative AGWP in W*yr/m^2/kg
    """
    indirect_effects = _check_indirect_effects_mode(indirect_effects)

    scenario = get_scenario()
    iam, ssp, rcp = scenario["iam"], scenario["ssp"], scenario["rcp"]

    # Load data
    re_data_ch4 = load_re_ch4(iam)
    re_data_co2 = load_re_co2(iam)
    irf_series = load_irf_ch4()

    years = re_data_ch4["_years"]
    re_series_ch4 = re_data_ch4[ssp][rcp]
    re_series_co2 = re_data_co2[ssp][rcp]

    year_idx = _get_year_index(emission_year, years)
    max_years = min(time_horizon, len(irf_series))

    if max_years <= 0:
        # Year 0 -> no effect
        return np.array([0], dtype="float64")

    total_rf = np.zeros(max_years, dtype="float64")

    temp_s = np.zeros(max_years, dtype="float64")
    temp_d = np.zeros(max_years, dtype="float64")
    c1 = np.zeros(max_years, dtype="float64")
    c2 = np.zeros(max_years, dtype="float64")
    c3 = np.zeros(max_years, dtype="float64")

    for t in range(max_years):
        re_ch4 = _get_re_value(re_series_ch4, year_idx, t, time_varying_re)
        re_co2 = _get_re_value(re_series_co2, year_idx, t, time_varying_re)

        # Apply indirect effects if requested
        re_ch4 = (
            re_ch4 + CH4_INDIRECT_RE
            if indirect_effects in {"no_carbon_cycle", "all"}
            else re_ch4
        )
        factor1 = re_ch4 * irf_series[t] * CONST_CH4
        factor1_model = factor1 * 1e12

        factor2 = 0.0
        if indirect_effects == "all" and t > 0:
            temp_s_prev = temp_s[t - 1] if t > 0 else 0
            temp_d_prev = temp_d[t - 1] if t > 0 else 0
            c1_prev = c1[t - 1] if t > 0 else 0

            temp_s[t] = (
                temp_s_prev
                + factor1_model / 7.7
                - 1.31 * temp_s_prev / 7.7
                - 0.88 * (temp_s_prev - temp_d_prev) / 7.7
            )
            temp_d[t] = temp_d_prev + (0.88 / 1.03) * (temp_s_prev - temp_d_prev) / 147

            c1[t] = temp_s[t] * 11.06 * 0.6368 + c1_prev * np.exp(-1 / 2.376)
            c2[t] = temp_s[t] * 11.06 * 0.3322 + c1_prev * np.exp(-1 / 30.14)
            c3[t] = temp_s[t] * 11.06 * 0.031 + c1_prev * np.exp(-1 / 490.1)

            factor2 = (c1[t] + c2[t] + c3[t]) * re_co2 * CONST_CO2

        total_rf[t] = factor1 + factor2

    cumulative_agwp = np.zeros(max_years, dtype="float64")
    for t in range(1, max_years):
        cumulative_agwp[t] = (
            cumulative_agwp[t - 1] + (total_rf[t - 1] + total_rf[t]) / 2
        )

    return cumulative_agwp


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

        re_t = _get_re_value(re_series, year_idx, t, time_varying_re)
        # AGWP contribution: RE * IRF * conversion * dt (dt=1 year)
        agwp += re_t * irf_t * CONST_CO2

    return agwp


def agwp_ch4(
    emission_year: int,
    time_horizon: int = 100,
    time_varying_re: bool = False,
    indirect_effects: str = "all",
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
    indirect_effects : str
        Indirect effects to include. Defaults to "all".

    Returns
    -------
    float
        AGWP in W*yr/m^2/kg
    """
    cumulative_agwp = _agwp_ch4_cumulative(
        emission_year=emission_year,
        time_horizon=time_horizon,
        time_varying_re=time_varying_re,
        indirect_effects=indirect_effects,
    )
    return float(cumulative_agwp[-1])


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

        re_t = _get_re_value(re_series, year_idx, t, time_varying_re)
        agwp += re_t * irf_t * CONST_N2O

    return agwp
