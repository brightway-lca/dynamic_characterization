"""The fair submodule must import even when `fair` is not installed."""

import importlib
import sys


def test_fair_submodule_imports_without_fair():
    # Other test modules load `dynamic_characterization` sub-packages in
    # isolation by installing lightweight stubs into ``sys.modules`` at
    # collection time. Those stubs can shadow the real package here (test
    # order is randomized), so this test — which validates the *real*
    # package contract — purges them and imports the genuine modules.
    for name in list(sys.modules):
        if name == "dynamic_characterization" or name.startswith(
            "dynamic_characterization."
        ):
            del sys.modules[name]

    mod = importlib.import_module("dynamic_characterization.fair")
    assert hasattr(mod, "FAIR_IMPORT_ERROR_MSG")
    assert "dynamic_characterization[fair]" in mod.FAIR_IMPORT_ERROR_MSG
