"""
Prospective characterization factors based on Barbosa Watanabe et al. (2026).

Implements scenario-based characterization using time-varying radiative
efficiencies from IAM-SSP-RCP scenarios.

Reference: https://doi.org/10.1021/acs.est.5b01118
"""

from . import agtp
from . import agwp
from .config import (
    VALID_SCENARIOS,
    get_scenario,
    reset_scenario,
    set_scenario,
)
from .data_loader import (
    load_irf_ch4,
    load_irf_co2,
    load_irf_n2o,
    load_re_co2,
    load_re_ch4,
    load_re_n2o,
)
from .radiative_forcing import (
    characterize_ch4,
    characterize_co2,
    characterize_co2_uptake,
    characterize_n2o,
)

__all__ = [
    "agtp",
    "agwp",
    "VALID_SCENARIOS",
    "get_scenario",
    "reset_scenario",
    "set_scenario",
    "load_irf_ch4",
    "load_irf_co2",
    "load_irf_n2o",
    "load_re_co2",
    "load_re_ch4",
    "load_re_n2o",
    "characterize_ch4",
    "characterize_co2",
    "characterize_co2_uptake",
    "characterize_n2o",
]
