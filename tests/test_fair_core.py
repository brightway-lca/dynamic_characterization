# tests/test_fair_core.py
import importlib.util
import os
import sys
import types

import numpy as np
import pandas as pd
import pytest

_fair_dir = os.path.join(
    os.path.dirname(__file__), "..", "dynamic_characterization", "fair"
)
_prospective_dir = os.path.join(
    os.path.dirname(__file__), "..", "dynamic_characterization", "prospective"
)

# ---------------------------------------------------------------------------
# Minimal package shims so core's relative imports resolve without bw2data.
# ---------------------------------------------------------------------------

# 1. dynamic_characterization stub with __path__ so sub-packages are findable.
_pkg_dir = os.path.join(os.path.dirname(__file__), "..", "dynamic_characterization")
_dc = sys.modules.get("dynamic_characterization") or types.ModuleType(
    "dynamic_characterization"
)
_dc.__path__ = [_pkg_dir]
sys.modules["dynamic_characterization"] = _dc

# 2. dynamic_characterization.fair package stub.
_fair_pkg = types.ModuleType("dynamic_characterization.fair")
_fair_pkg.__path__ = [_fair_dir]
_fair_pkg.FAIR_IMPORT_ERROR_MSG = "install dynamic_characterization[fair]"
sys.modules["dynamic_characterization.fair"] = _fair_pkg

# 3. dynamic_characterization.prospective package stub.
_prospective_pkg = sys.modules.get(
    "dynamic_characterization.prospective"
) or types.ModuleType("dynamic_characterization.prospective")
_prospective_pkg.__path__ = [_prospective_dir]
_prospective_pkg.__package__ = "dynamic_characterization.prospective"
sys.modules["dynamic_characterization.prospective"] = _prospective_pkg


def _load_fair(modname, filename):
    """Load a module from the fair subpackage directory."""
    spec = importlib.util.spec_from_file_location(
        f"dynamic_characterization.fair.{modname}",
        os.path.join(_fair_dir, filename),
    )
    mod = importlib.util.module_from_spec(spec)
    sys.modules[f"dynamic_characterization.fair.{modname}"] = mod
    spec.loader.exec_module(mod)
    return mod


def _load_prospective(modname, filename):
    """Load a module from the prospective subpackage directory."""
    spec = importlib.util.spec_from_file_location(
        f"dynamic_characterization.prospective.{modname}",
        os.path.join(_prospective_dir, filename),
    )
    mod = importlib.util.module_from_spec(spec)
    sys.modules[f"dynamic_characterization.prospective.{modname}"] = mod
    spec.loader.exec_module(mod)
    return mod


# Load prospective modules that core.py imports via relative import.
if "dynamic_characterization.prospective.config" not in sys.modules:
    _load_prospective("config", "config.py")

# Load fair submodules.
_load_fair("allocation", "allocation.py")
_load_fair("species_map", "species_map.py")
_load_fair("runner", "runner.py")
core = _load_fair("core", "core.py")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _inv():
    # NOTE: flow_name column is required for species resolution.
    # The brief's _inv() omits flow_name; without it, str(flow_id) is passed to
    # resolve_species, which cannot match "carbon dioxide" from id "1". Task 7
    # guarantees the real inventory includes flow names; we add them here so the
    # helper test is meaningful.
    return pd.DataFrame(
        {
            "date": pd.to_datetime(["2030-01-01", "2030-01-01", "2031-01-01"]),
            "amount": [10.0, 2.0, -5.0],
            "flow": [1, 2, 1],
            "activity": ["a", "b", "a"],
            "flow_name": [
                "Carbon dioxide, fossil",
                "Methane, fossil",
                "Carbon dioxide, fossil",
            ],
        }
    )


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_emissions_by_species_aggregates_and_signs():
    years = np.array([2030, 2031])
    by_species, _ = core._inventory_emissions_by_species(_inv(), years)
    # flow 1 = CO2 (sign +1): 10 in 2030, -5 in 2031
    assert "CO2 FFI" in by_species
    np.testing.assert_allclose(by_species["CO2 FFI"], [10.0, -5.0])


def test_characterize_with_fair_runs(monkeypatch):
    pytest.importorskip("fair")
    import dynamic_characterization.prospective.config as cfg  # noqa

    # use a fair-capable scenario
    core.config.set_scenario("FAIR", "SSP2", "4.5")
    out = core.characterize_with_fair(
        _inv(), output="radiative_forcing", quantiles=(50.0,)
    )
    assert list(out.columns) == ["date", "amount", "flow", "activity", "quantile"]
    assert set(out["quantile"].unique()) == {50.0}
