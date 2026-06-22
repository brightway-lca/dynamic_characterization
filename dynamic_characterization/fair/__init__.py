"""
FAIR climate-model dynamic LCIA path (optional).

Requires the optional `fair` dependency:
    pip install dynamic_characterization[fair]

Importing this submodule never imports `fair`; the dependency is only
loaded when a FAIR characterization is actually run.
"""

FAIR_IMPORT_ERROR_MSG = (
    "The `fair` climate model is required for FAIR metrics. Install it with "
    "`pip install dynamic_characterization[fair]`."
)

__all__ = ["FAIR_IMPORT_ERROR_MSG"]
