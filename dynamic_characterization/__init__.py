"""
Collection of dynamic characterization functions for life cycle inventories with temporal information.

Provides:
- ipcc_ar6: characterization functions following the IPCC AR6 emission metrics (CO2, CH4, N2O,
  halocarbons) and the indirect effects of CO, VOC, H2 and NOx
- prospective: prospective characterization after Barbosa Watanabe & Cherubini (2026)
- original_temporalis_functions: Legacy functions from bw_temporalis

For prospective metrics (pGWP, pGTP), set the scenario first:
    import dynamic_characterization.prospective as prospective
    prospective.set_scenario(iam="IMAGE", ssp="SSP1", rcp="2.6")
"""

__all__ = (
    "__version__",
    "characterize",
    "create_characterization_functions_from_method",
    "compare_with_static",
    "CharacterizationOptions",
    "original_temporalis_functions",
    "ipcc_ar6",
    "prospective",
)

__version__ = "2.0.0"

from . import ipcc_ar6, original_temporalis_functions, prospective
from .dynamic_characterization import (
    characterize,
    compare_with_static,
    create_characterization_functions_from_method,
)
from .options import CharacterizationOptions
