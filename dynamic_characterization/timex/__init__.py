"""
Dynamic characterization functions from the bw_timex package (https://github.com/brightway-lca/bw_timex).
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

from .radiative_forcing import characterize_co2
from .radiative_forcing import characterize_co2_uptake
from .radiative_forcing import characterize_co
from .radiative_forcing import characterize_ch4
from .radiative_forcing import characterize_n2o
from .radiative_forcing import create_generic_characterization_function