"""Write the golden snapshot of the legacy (v1.4.3) IPCC AR6 characterization functions.

The snapshot pins the numbers of the default ("legacy") characterization functions, so that
later changes (bug fixes, new parameter sets) can be checked to leave them bit-identical.
Regenerate only on purpose, from the commit whose numbers should become the reference:

    python dev/make_legacy_golden.py

Output: tests/data/legacy_golden_v1.4.3.npz
"""

import json
import os
from collections import namedtuple

import numpy as np
import pandas as pd

from dynamic_characterization import characterize
from dynamic_characterization.ipcc_ar6.radiative_forcing import (
    characterize_ch4,
    characterize_co,
    characterize_co2,
    characterize_co2_uptake,
    characterize_n2o,
    create_generic_characterization_function,
)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "tests", "data", "legacy_golden_v1.4.3.npz")
DECAY_MULTIPLIERS = os.path.join(
    ROOT, "dynamic_characterization", "ipcc_ar6", "data", "decay_multipliers.json"
)

PERIODS = (20, 100, 500)
# CFC-11, CF4, HFC-23, HFC-134a, SF6 (keys as stored in decay_multipliers.json)
HALOCARBON_KEYS = ("000075-69-4", "000075-73-0", "000075-46-7", "000811-97-2", "002551-62-4")

Row = namedtuple("Row", ["date", "amount", "flow", "activity"])


def named_functions():
    with open(DECAY_MULTIPLIERS, encoding="utf-8") as f:
        decay_multipliers = json.load(f)
    functions = {
        "co2": characterize_co2,
        "co2_uptake": characterize_co2_uptake,
        "ch4": characterize_ch4,
        "n2o": characterize_n2o,
        "co": characterize_co,
    }
    for key in HALOCARBON_KEYS:
        functions[f"generic_{key}"] = create_generic_characterization_function(
            np.array(decay_multipliers[key])
        )
    return functions


def inventory(flow_ids):
    """Two emissions per flow, different years and amounts."""
    n = len(flow_ids)
    return pd.DataFrame(
        {
            "date": pd.Series(["2024-01-01"] * n + ["2030-07-01"] * n, dtype="datetime64[s]"),
            "amount": pd.Series([1.0 + i for i in range(n)] + [3.5] * n, dtype="float64"),
            "flow": pd.Series(list(flow_ids) * 2, dtype="int"),
            "activity": pd.Series([1] * n + [2] * n, dtype="int"),
        }
    )


def compute():
    """Return {name: array} for every pinned quantity."""
    functions = named_functions()
    row = Row(date=pd.Timestamp("2024-01-01"), amount=2.5, flow=1, activity=1)
    out = {}
    for name, func in functions.items():
        for period in PERIODS:
            for cumulative in (False, True):
                result = func(row, period, cumulative=cumulative)
                tag = f"{name}__p{period}__{'cum' if cumulative else 'marg'}"
                out[tag + "__amount"] = np.asarray(result.amount, dtype="float64")
                out[tag + "__date"] = np.asarray(result.date, dtype="datetime64[s]").astype("int64")

    # end-to-end through characterize(), with custom functions (no bw2data project needed)
    flow_functions = {i + 1: func for i, func in enumerate(functions.values())}
    df = inventory(flow_functions)
    for metric in ("radiative_forcing", "GWP"):
        for fixed in (False, True):
            result = characterize(
                df,
                metric=metric,
                characterization_functions=flow_functions,
                time_horizon=100,
                fixed_time_horizon=fixed,
                time_horizon_start=pd.Timestamp("2020-01-01"),
            )
            tag = f"characterize__{metric}__fixed{int(fixed)}"
            out[tag + "__amount"] = result["amount"].to_numpy(dtype="float64")
            out[tag + "__date"] = result["date"].to_numpy(dtype="datetime64[s]").astype("int64")
            out[tag + "__flow"] = result["flow"].to_numpy(dtype="int64")
    return out


if __name__ == "__main__":
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    data = compute()
    np.savez_compressed(OUT, **data)
    print(f"wrote {len(data)} arrays to {OUT}")
