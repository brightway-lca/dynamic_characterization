"""
Dynamic characterization functions from the bw_temporalis package (https://github.com/brightway-lca/bw_temporalis).
"""

__all__ = (
    "__version__",
    "characterize_co2",
    "characterize_methane",
    # Add functions and variables you want exposed in `dynamic_characterization.` namespace here
)

__version__ = "0.0.1dev1"

from . import characterize_co2
from . import characterize_methane
