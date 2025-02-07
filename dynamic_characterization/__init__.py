"""
Collection of dynamic characterization functions for life cycle inventories with temporal information.
"""

__all__ = (
    "__version__",
    "dynamic_characterization.characterize",
    "dynamic_characterization.create_characterization_functions_from_method",
    "original_temporalis_functions",
    "ipcc_ar6",
    # Add functions and variables you want exposed in `dynamic_characterization.` namespace here
)

__version__ = "1.1.0"

from . import ipcc_ar6, original_temporalis_functions
from .dynamic_characterization import (
    characterize,
    create_characterization_functions_from_method,
)
