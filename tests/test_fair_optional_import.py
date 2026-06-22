"""The fair submodule must import even when `fair` is not installed."""

import importlib


def test_fair_submodule_imports_without_fair():
    mod = importlib.import_module("dynamic_characterization.fair")
    assert hasattr(mod, "FAIR_IMPORT_ERROR_MSG")
    assert "dynamic_characterization[fair]" in mod.FAIR_IMPORT_ERROR_MSG
