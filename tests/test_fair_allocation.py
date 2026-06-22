import importlib.util
import os
import sys
import types

import numpy as np

_fair_dir = os.path.join(
    os.path.dirname(__file__), "..", "dynamic_characterization", "fair"
)
sys.modules.setdefault(
    "dynamic_characterization", types.ModuleType("dynamic_characterization")
)
_pkg = types.ModuleType("dynamic_characterization.fair")
_pkg.__path__ = [_fair_dir]
sys.modules["dynamic_characterization.fair"] = _pkg
_spec = importlib.util.spec_from_file_location(
    "dynamic_characterization.fair.allocation",
    os.path.join(_fair_dir, "allocation.py"),
)
allocation = importlib.util.module_from_spec(_spec)
sys.modules["dynamic_characterization.fair.allocation"] = allocation
_spec.loader.exec_module(allocation)


def test_per_kg_handles_zero_denominator():
    out = allocation.per_kg_response(np.array([1.0, 2.0]), np.array([0.0, 4.0]))
    assert out[0] == 0.0
    assert out[1] == 0.5


def test_allocation_sums_to_total_response():
    # Two flows, 3 years. Sum over flows of allocated == species response.
    flow_cum = np.array([[1.0, 2.0, 3.0], [1.0, 1.0, 1.0]])  # (2 flows, 3 yr)
    species_cum = flow_cum.sum(axis=0)                        # (3,)
    response = np.array([0.0, 10.0, 30.0])
    per_kg = allocation.per_kg_response(response, species_cum)
    alloc = allocation.allocate_to_flows(flow_cum, per_kg)
    np.testing.assert_allclose(alloc.sum(axis=0), response)


def test_safe_nanpercentile_all_nan_slice_is_zero():
    stack = np.array([[np.nan, 1.0], [np.nan, 3.0]])  # (2 configs, 2 yr)
    out = allocation.safe_nanpercentile(stack, [50.0], axis=0)
    assert out[0, 0] == 0.0
    assert out[0, 1] == 2.0
