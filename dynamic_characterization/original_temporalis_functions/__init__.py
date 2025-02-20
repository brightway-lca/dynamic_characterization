"""
Dynamic characterization functions originally introduced by the bw_temporalis package (https://github.com/brightway-lca/bw_temporalis).
"""

__all__ = (
    "characterize_co2",
    "characterize_methane",
    # Add functions and variables you want exposed in `dynamic_characterization.` namespace here
)

from .radiative_forcing import characterize_co2, characterize_methane
