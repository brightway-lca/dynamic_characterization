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
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(autouse=True)
def _reset_scenario():
    """Reset scenario state before and after each test."""
    try:
        core.config.reset_scenario()
    except Exception:
        pass
    yield
    try:
        core.config.reset_scenario()
    except Exception:
        pass


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


def test_zero_row_inventory_returns_empty():
    empty = pd.DataFrame(
        {"date": pd.to_datetime([]), "amount": [], "flow": [], "activity": []}
    )
    out = core.characterize_with_fair(empty)
    assert list(out.columns) == ["date", "amount", "flow", "activity", "quantile"]
    assert len(out) == 0


def test_unmapped_flows_warn_and_return_empty():
    df = pd.DataFrame(
        {
            "date": pd.to_datetime(["2030-01-01"]),
            "amount": [1.0],
            "flow": [999],
            "activity": ["a"],
            "flow_name": ["Some unmappable flow"],
        }
    )
    out = core.characterize_with_fair(df)
    assert list(out.columns) == ["date", "amount", "flow", "activity", "quantile"]
    assert len(out) == 0


def test_characterize_with_fair_runs(monkeypatch):
    pytest.importorskip("fair")
    import dynamic_characterization.prospective.config as cfg  # noqa

    # use a fair-capable scenario
    core.config.set_fair_scenario("SSP2", "4.5")
    out = core.characterize_with_fair(
        _inv(), output="radiative_forcing", quantiles=(50.0,)
    )
    assert list(out.columns) == ["date", "amount", "flow", "activity", "quantile"]
    assert set(out["quantile"].unique()) == {50.0}


# ---------------------------------------------------------------------------
# Vectorized aggregation / attribution
# ---------------------------------------------------------------------------


def _naive_aggregate(df, years):
    """Row-by-row reference for _aggregate_inventory."""
    year_index = {int(y): i for i, y in enumerate(years)}
    per_pair = {}
    for _, row in df.iterrows():
        species, sign = core.species_map.resolve_species(str(row["flow_name"]))
        if species is None:
            continue
        yi = year_index.get(int(pd.Timestamp(row["date"]).year))
        if yi is None:
            continue
        key = (row["flow"], row["activity"])
        signed = per_pair.setdefault(key, (species, np.zeros(len(years))))[1]
        signed[yi] += sign * float(row["amount"])
    return per_pair


def test_aggregate_inventory_matches_row_loop():
    df = pd.DataFrame(
        {
            "date": pd.to_datetime(
                ["2030-01-01", "2030-06-01", "2031-01-01", "2030-01-01", "2031-01-01"]
            ),
            "amount": [10.0, 4.0, -5.0, 2.0, 7.0],
            "flow": [1, 1, 1, 2, 3],
            "activity": ["a", "a", "a", "b", "a"],
            "flow_name": [
                "Carbon dioxide, fossil",
                "Carbon dioxide, fossil",
                "Carbon dioxide, fossil",
                "Methane, fossil",
                "Some unmappable flow",
            ],
        }
    )
    years = np.array([2030, 2031])
    species_names, species_of_pair, flow_ids, activity_ids, yearly = (
        core._aggregate_inventory(df, years)
    )
    expected = _naive_aggregate(df, years)

    assert len(flow_ids) == len(expected)  # the unmappable flow is dropped
    for i, (flow, activity) in enumerate(zip(flow_ids, activity_ids)):
        species, signed = expected[(flow, activity)]
        assert species_names[species_of_pair[i]] == species
        np.testing.assert_allclose(yearly[i], signed)


def test_uptake_flows_get_mirrored_quantiles(monkeypatch):
    """A flow with negative cumulative emissions flips the quantile order."""
    quantiles = (2.5, 50.0, 97.5)
    # Three-config ensemble with a clearly ordered response, in both years.
    perturbed = np.array([[1.0, 2.0], [2.0, 4.0], [3.0, 6.0]])

    def fake_runs(marker, years_, perturbations, max_batch=None):
        zero = np.zeros_like(perturbed)
        return [
            {
                "radiative_forcing": zero if not perturbation else perturbed,
                "temperature": zero if not perturbation else perturbed,
                "scale": 1.0,
            }
            for perturbation in perturbations
        ]

    monkeypatch.setattr(core.runner, "run_perturbations", fake_runs)
    monkeypatch.setattr(core.runner, "require_fair", lambda: None)
    core.config.set_fair_scenario("SSP2", "4.5")

    df = pd.DataFrame(
        {
            "date": pd.to_datetime(["2030-01-01", "2030-01-01"]),
            "amount": [3.0, 1.0],
            "flow": [1, 2],
            "activity": ["release", "uptake"],
            # Both are biogenic CO2 (same FAIR species), one released, one
            # taken up - so they share a per-kg response and mirror each other.
            "flow_name": ["Carbon dioxide, non-fossil", "Carbon dioxide, in air"],
        }
    )
    out = core.characterize_with_fair(
        df, output="radiative_forcing", quantiles=quantiles, time_horizon=1
    )

    # Net species emission is 3 - 1 = 2 kg, constant over both years.
    per_kg = perturbed / 2.0
    for year_index, year in enumerate((2030, 2031)):
        rows = out[out["date"] == np.datetime64(f"{year}-01-01")]
        release = rows[rows["activity"] == "release"].set_index("quantile")
        uptake = rows[rows["activity"] == "uptake"].set_index("quantile")
        for quantile in quantiles:
            expected_release = (
                np.percentile(per_kg[:, year_index], quantile) * 3.0
            )
            # Cumulative uptake is negative, so the mirrored quantile applies.
            expected_uptake = (
                np.percentile(per_kg[:, year_index], 100.0 - quantile) * -1.0
            )
            np.testing.assert_allclose(
                release.loc[quantile, "amount"], expected_release
            )
            np.testing.assert_allclose(
                uptake.loc[quantile, "amount"], expected_uptake
            )
