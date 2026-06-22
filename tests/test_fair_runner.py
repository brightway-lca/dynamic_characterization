import importlib.util
import os
import sys
import types

import pytest

_fair_dir = os.path.join(
    os.path.dirname(__file__), "..", "dynamic_characterization", "fair"
)
_dc = sys.modules.setdefault(
    "dynamic_characterization", types.ModuleType("dynamic_characterization")
)
_pkg = types.ModuleType("dynamic_characterization.fair")
_pkg.__path__ = [_fair_dir]
_pkg.FAIR_IMPORT_ERROR_MSG = "install dynamic_characterization[fair]"
sys.modules["dynamic_characterization.fair"] = _pkg
_spec = importlib.util.spec_from_file_location(
    "dynamic_characterization.fair.runner",
    os.path.join(_fair_dir, "runner.py"),
)
runner = importlib.util.module_from_spec(_spec)
sys.modules["dynamic_characterization.fair.runner"] = runner
_spec.loader.exec_module(runner)


def test_require_fair_raises_when_absent():
    if importlib.util.find_spec("fair") is not None:
        pytest.skip("fair installed")
    with pytest.raises(ImportError, match=r"dynamic_characterization\[fair\]"):
        runner.require_fair()


def test_find_calibration_files_returns_tuple():
    pytest.importorskip("fair")
    params, props = runner.find_calibration_files()
    # Either both found or both None, but always a 2-tuple
    assert (params is None) == (props is None)
    assert isinstance((params, props), tuple)


def test_run_fair_rejects_bad_output():
    import pytest

    with pytest.raises(ValueError, match="output must be"):
        runner.run_fair(
            "ssp245", None, __import__("numpy").array([2030, 2031]), output="bogus"
        )
