"""
Prospective characterization based on Barbosa Watanabe & Cherubini (2026).

Scenario-based characterization of CO2, CH4 and N2O for IAM-SSP-RCP scenarios: the radiative
efficiencies and, for CO2, the impulse response function depend on the scenario. CH4 and N2O
include their indirect effects and a climate-carbon feedback (see ``watanabe2026``, the paper's
method with two errors of its published code corrected).

Reference: https://doi.org/10.1021/acs.est.5c12391
"""

from . import agtp
from . import agwp
from . import watanabe2026
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
    "watanabe2026",
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
