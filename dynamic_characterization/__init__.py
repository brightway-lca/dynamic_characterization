"""
Collection of dynamic characterization functions for life cycle inventories with temporal information.
"""

__all__ = (
    "__version__",
    "dynamic_characterization.characterize_dynamic_inventory",
    "original_temporalis_functions",
    "ipcc_ar6",
    # Add functions and variables you want exposed in `dynamic_characterization.` namespace here
)

__version__ = "0.0.4"

from .dynamic_characterization import characterize
from . import original_temporalis_functions
from . import ipcc_ar6