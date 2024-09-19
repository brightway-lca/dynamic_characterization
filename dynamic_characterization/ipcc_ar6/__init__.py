"""
Dynamic characterization functions from based on IPCC AR6.

Relevant scientific publications:
Schivley2015: Relevant scientific publication on the numerical calculation of CRF: https://doi.org/10.1021/acs.est.5b01118
Forster2023: Numerical values from IPCC AR6 Chapter 7 (Table 7.15): https://doi.org/10.1017/9781009157896.009
"""

__all__ = (
    "__version__",
    "characterize_co2",
    "characterize_co2_uptake",
    "characterize_co",
    "characterize_ch4",
    "characterize_n2o",
    "create_generic_characterization_function",
    # Add functions and variables you want exposed in `dynamic_characterization.` namespace here
)

__version__ = "0.0.1"

from .radiative_forcing import (
    characterize_ch4,
    characterize_co,
    characterize_co2,
    characterize_co2_uptake,
    characterize_n2o,
    create_generic_characterization_function,
)
