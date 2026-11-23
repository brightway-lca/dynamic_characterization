"""The legacy IPCC AR6 characterization functions must stay bit-identical to v1.4.3.

The reference arrays are written by ``dev/make_legacy_golden.py``.

The generator is loaded at import time: test_prospective.py replaces `dynamic_characterization`
in sys.modules with a stub when it is collected, so the real package must be imported before.
"""

import importlib.util
import os

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GOLDEN = os.path.join(ROOT, "tests", "data", "legacy_golden_v1.4.3.npz")


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


def test_same_quantities(current, golden):
    assert set(current) == set(golden)


def test_bit_identical(current, golden):
    mismatched = [
        key for key in golden if not np.array_equal(current[key], golden[key])
    ]
    assert not mismatched, f"legacy output changed for: {mismatched}"
