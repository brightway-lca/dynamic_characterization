"""
Collection of dynamic characterization functions for life cycle inventories with temporal information.

Provides:
- ipcc_ar6: IPCC AR6-based characterization functions for CO2, CH4, N2O, etc.
- prospective: Prospective characterization factors from Watanabe et al. (2026)
- original_temporalis_functions: Legacy functions from bw_temporalis

For prospective metrics (pGWP, pGTP), set the scenario first:
    import dynamic_characterization.prospective as prospective
    prospective.set_scenario(iam="IMAGE", ssp="SSP1", rcp="2.6")
"""

__all__ = (
    "__version__",
    "characterize",
    "create_characterization_functions_from_method",
    "original_temporalis_functions",
    "ipcc_ar6",
    "prospective",
)

__version__ = "1.3.1"

from . import ipcc_ar6, original_temporalis_functions, prospective
from .dynamic_characterization import (
    characterize,
    create_characterization_functions_from_method,
)
