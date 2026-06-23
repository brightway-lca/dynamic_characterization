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
    result = runner.find_calibration_files()
    # Always a 2-tuple; each entry is a path or None depending on the install.
    assert isinstance(result, tuple)
    assert len(result) == 2
    for entry in result:
        assert entry is None or isinstance(entry, str)


def test_get_calibration_files_returns_existing_paths():
    pytest.importorskip("fair")
    # Downloads the calibration 1.4.1 pair on first use (cached afterwards).
    params, props = runner.get_calibration_files()
    assert os.path.exists(params)
    assert os.path.exists(props)


def test_unit_factor_co2_vs_other():
    assert runner._unit_factor("CO2 FFI") == 1e-12
    assert runner._unit_factor("CO2 AFOLU") == 1e-12
    assert runner._unit_factor("CH4") == 1e-9
    assert runner._unit_factor("Sulfur") == 1e-9


def test_run_fair_rejects_bad_output():
    import pytest

    with pytest.raises(ValueError, match="output must be"):
        runner.run_fair(
            "ssp245", None, __import__("numpy").array([2030, 2031]), output="bogus"
        )
