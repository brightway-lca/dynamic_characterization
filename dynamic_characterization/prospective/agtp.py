"""
Absolute Global Temperature Potential (AGTP) calculations.

Based on Watanabe et al. (2026) using two-layer energy balance model.

AGTP = integral_0^t RF(t') * R(t-t') dt'

where R is the temperature impulse response function with surface and deep ocean components.
"""

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
from .agwp import (
    CONST_CH4,
    CONST_CO2,
    CONST_N2O,
    _get_year_index,
    _rcp_to_irf_key,
)

# Temperature response parameters from Watanabe SI code
# Two-layer model: surface layer and deep ocean
C_S = 7.7  # Heat capacity surface layer (W*yr/m^2/K)
C_D = 147  # Heat capacity deep ocean (W*yr/m^2/K)
LAMBDA = 1.31  # Climate feedback parameter magnitude (W/m^2/K)
GAMMA = 0.88  # Heat exchange coefficient (W/m^2/K)
EFFICACY = 1.03  # Efficacy factor


def _temperature_response(t: float) -> float:
    """
    Calculate temperature impulse response R(t) at time t.

    Uses two-box model from Geoffroy et al. (2013).

    The two-layer energy balance model has eigenvalues that give two
    characteristic timescales: a fast response (~3-8 years) for the
    surface/mixed layer and a slow response (~200-400 years) for the
    deep ocean.

    Parameters
    ----------
    t : float
        Time in years since forcing

    Returns
    -------
    float
        Temperature response in K per (W/m^2)
    """
    # Eigenvalues of the two-layer system
    # From the coupled ODEs: C_S * dT_S/dt = F - LAMBDA*T_S - GAMMA*(T_S - T_D)
    #                        C_D * dT_D/dt = GAMMA*(T_S - T_D)
    b = (LAMBDA + GAMMA) / C_S + GAMMA / C_D
    c = (LAMBDA * GAMMA) / (C_S * C_D)

    discriminant = b * b - 4 * c
    if discriminant < 0:
        return 0.0

    sqrt_disc = np.sqrt(discriminant)
    lambda1 = (b + sqrt_disc) / 2  # Fast mode (surface)
    lambda2 = (b - sqrt_disc) / 2  # Slow mode (deep ocean)

    # Time constants
    tau1 = 1 / lambda1  # ~3-4 years
    tau2 = 1 / lambda2  # ~280 years

    # Coefficients from eigenvector analysis
    q1 = (1 / C_S) * (LAMBDA + GAMMA - C_S * lambda2) / (lambda1 - lambda2)
    q2 = (1 / C_S) * (C_S * lambda1 - LAMBDA - GAMMA) / (lambda1 - lambda2)

    return EFFICACY * (q1 * np.exp(-t / tau1) + q2 * np.exp(-t / tau2))


def agtp_co2(
    emission_year: int,
    time_horizon: int = 100,
    time_varying_re: bool = False,
) -> float:
    """
    Calculate AGTP for 1 kg CO2.

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
        AGTP in K/kg
    """
    scenario = get_scenario()
    iam, ssp, rcp = scenario["iam"], scenario["ssp"], scenario["rcp"]

    re_data = load_re_co2(iam)
    irf_data = load_irf_co2()

    years = re_data["_years"]
    re_series = re_data[ssp][rcp]
    irf_series = irf_data[_rcp_to_irf_key(rcp)]

    year_idx = _get_year_index(emission_year, years)

    max_years = min(time_horizon, len(irf_series))

    agtp = 0.0
    for t_prime in range(max_years):
        irf_t = irf_series[t_prime]

        if time_varying_re:
            re_idx = min(year_idx + t_prime, len(re_series) - 1)
            re_t = re_series[re_idx]
        else:
            re_t = re_series[year_idx]

        rf_t = re_t * irf_t * CONST_CO2

        # Temperature response R(time_horizon - t')
        delta_t = time_horizon - t_prime - 1
        if delta_t >= 0:
            r_t = _temperature_response(delta_t)
            agtp += rf_t * r_t

    return agtp


def agtp_ch4(
    emission_year: int,
    time_horizon: int = 100,
    time_varying_re: bool = False,
) -> float:
    """
    Calculate AGTP for 1 kg CH4.

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
        AGTP in K/kg
    """
    scenario = get_scenario()
    iam, ssp, rcp = scenario["iam"], scenario["ssp"], scenario["rcp"]

    re_data = load_re_ch4(iam)
    irf_series = load_irf_ch4()

    years = re_data["_years"]
    re_series = re_data[ssp][rcp]

    year_idx = _get_year_index(emission_year, years)

    max_years = min(time_horizon, len(irf_series))

    agtp = 0.0
    for t_prime in range(max_years):
        irf_t = irf_series[t_prime]

        if time_varying_re:
            re_idx = min(year_idx + t_prime, len(re_series) - 1)
            re_t = re_series[re_idx]
        else:
            re_t = re_series[year_idx]

        rf_t = re_t * irf_t * CONST_CH4

        delta_t = time_horizon - t_prime - 1
        if delta_t >= 0:
            r_t = _temperature_response(delta_t)
            agtp += rf_t * r_t

    return agtp


def agtp_n2o(
    emission_year: int,
    time_horizon: int = 100,
    time_varying_re: bool = False,
) -> float:
    """
    Calculate AGTP for 1 kg N2O.

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
        AGTP in K/kg
    """
    scenario = get_scenario()
    iam, ssp, rcp = scenario["iam"], scenario["ssp"], scenario["rcp"]

    re_data = load_re_n2o(iam)
    irf_series = load_irf_n2o()

    years = re_data["_years"]
    re_series = re_data[ssp][rcp]

    year_idx = _get_year_index(emission_year, years)

    max_years = min(time_horizon, len(irf_series))

    agtp = 0.0
    for t_prime in range(max_years):
        irf_t = irf_series[t_prime]

        if time_varying_re:
            re_idx = min(year_idx + t_prime, len(re_series) - 1)
            re_t = re_series[re_idx]
        else:
            re_t = re_series[year_idx]

        rf_t = re_t * irf_t * CONST_N2O

        delta_t = time_horizon - t_prime - 1
        if delta_t >= 0:
            r_t = _temperature_response(delta_t)
            agtp += rf_t * r_t

    return agtp
