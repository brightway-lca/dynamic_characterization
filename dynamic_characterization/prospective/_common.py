"""Constants and helpers shared by the prospective modules (Watanabe et al. 2026 data)."""

import warnings

import numpy as np

# Constants for unit conversion: RE from W/m^2/ppb to W/m^2/kg
M_AIR = 28.97  # g/mol (average molar mass of dry air)
M_ATMOSPHERE = 5.13252e18  # kg (total mass of atmosphere)

# Molar masses (g/mol)
M_CO2 = 44.01
M_CH4 = 16.04
M_N2O = 44.01

# kg per ppb
PPB_TO_KG_CO2 = 1e-9 * M_ATMOSPHERE / M_AIR * M_CO2
PPB_TO_KG_CH4 = 1e-9 * M_ATMOSPHERE / M_AIR * M_CH4
PPB_TO_KG_N2O = 1e-9 * M_ATMOSPHERE / M_AIR * M_N2O

# ppb per kg (to convert RE from per-ppb to per-kg)
CONST_CO2 = 1.0 / PPB_TO_KG_CO2
CONST_CH4 = 1.0 / PPB_TO_KG_CH4
CONST_N2O = 1.0 / PPB_TO_KG_N2O


def _get_year_index(emission_year: int, years: np.ndarray) -> int:
    """
    Index of the emission year in the RE data, clamped (with a warning) to the years the data
    covers (2020-2150).
    """
    min_year, max_year = int(years[0]), int(years[-1])
    if emission_year < min_year:
        warnings.warn(f"Emission year {emission_year} < {min_year}, clamping to {min_year}")
        emission_year = min_year
    elif emission_year > max_year:
        warnings.warn(f"Emission year {emission_year} > {max_year}, clamping to {max_year}")
        emission_year = max_year
    return int(np.searchsorted(years, emission_year))


def _rcp_to_irf_key(rcp: str) -> str:
    """Convert an RCP string (e.g. '2.6') to the CO2 IRF key (e.g. 'RCP26')."""
    return f"RCP{rcp.replace('.', '')}"
