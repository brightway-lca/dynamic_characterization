"""The legacy IPCC AR6 characterization functions must stay unchanged from v1.4.3.

The reference arrays are written by ``dev/make_legacy_golden.py`` (on macOS arm64).

Dates and flows must match exactly. Amounts must match to floating-point round-off: NumPy's
``exp``/``log`` differ in the last bits between platforms and SIMD implementations (the
snapshot is bit-identical on macOS, but not on Linux x86_64), and ``np.diff`` and the GWP ratio
carry that through. A relative tolerance of 1e-12 of each array's scale is far below any real
change of the functions.

The generator is loaded at import time: test_prospective.py replaces `dynamic_characterization`
in sys.modules with a stub when it is collected, so the real package must be imported before.
"""

import importlib.util
import os

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GOLDEN = os.path.join(ROOT, "tests", "data", "legacy_golden_v1.4.3.npz")

RTOL = 1e-12


def _load_generator():
    spec = importlib.util.spec_from_file_location(
        "make_legacy_golden", os.path.join(ROOT, "dev", "make_legacy_golden.py")
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_GENERATOR = _load_generator()


@pytest.fixture(scope="module")
def current():
    return _GENERATOR.compute()


@pytest.fixture(scope="module")
def golden():
    with np.load(GOLDEN) as data:
        return {key: data[key] for key in data.files}


def _aligned(data):
    """
    Sort the rows of each characterize() result by (flow, date).

    characterize() sorts by (date, amount); with round-off differences in the amounts, rows
    with nearly equal amounts could swap places. (flow, date) is unique per row here.
    """
    data = dict(data)
    for prefix in {key.rsplit("__", 1)[0] for key in data if key.startswith("characterize__")}:
        order = np.lexsort((data[prefix + "__date"], data[prefix + "__flow"]))
        for field in ("amount", "date", "flow"):
            data[f"{prefix}__{field}"] = data[f"{prefix}__{field}"][order]
    return data


def test_same_quantities(current, golden):
    assert set(current) == set(golden)


def test_legacy_output_unchanged(current, golden):
    current, golden = _aligned(current), _aligned(golden)
    changed = []
    for key, expected in golden.items():
        actual = current[key]
        if actual.shape != expected.shape:
            changed.append(f"{key} (shape {actual.shape} != {expected.shape})")
        elif key.endswith("__amount"):
            scale = np.abs(expected).max() if expected.size else 0.0
            deviation = np.abs(actual - expected).max() / scale if scale else np.abs(actual).max()
            if deviation > RTOL:
                changed.append(f"{key} (max deviation {deviation:.1e} of max |value|)")
        elif not np.array_equal(actual, expected):
            changed.append(key)
    assert not changed, "legacy output changed for: " + ", ".join(changed)
