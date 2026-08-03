"""Tests for the default characterization functions built from a static LCIA method."""

import bw2data as bd
import numpy as np
import pandas as pd
import pytest
from bw2data.tests import bw2test

from dynamic_characterization import characterize
from dynamic_characterization.dynamic_characterization import (
    clear_characterization_function_cache,
    create_characterization_functions_from_method,
)
from dynamic_characterization.ipcc_ar6.radiative_forcing import (
    characterize_ch4,
    characterize_co2,
    characterize_co2_uptake,
)

METHOD = ("test", "climate change", "GWP100")

BIOFLOWS = {
    "co2_fossil_high_stacks": {
        "name": "Carbon dioxide, fossil",
        "categories": ("air", "non-urban air or from high stacks"),
        "type": "emission",
        "CAS number": "000124-38-9",
    },
    "co2_fossil_unspecified": {
        "name": "Carbon dioxide, fossil",
        "categories": ("air",),
        "type": "emission",
        "CAS number": "000124-38-9",
    },
    "co2_fossil_urban": {
        "name": "Carbon dioxide, fossil",
        "categories": ("air", "urban air close to ground"),
        "type": "emission",
        "CAS number": "000124-38-9",
    },
    "co2_fossil_stratosphere": {
        "name": "Carbon dioxide, fossil",
        "categories": ("air", "lower stratosphere + upper troposphere"),
        "type": "emission",
        "CAS number": "000124-38-9",
    },
    "co2_non_fossil": {
        "name": "Carbon dioxide, non-fossil",
        "categories": ("air", "non-urban air or from high stacks"),
        "type": "emission",
        "CAS number": "000124-38-9",
    },
    "co2_from_soil_or_biomass_stock": {
        "name": "Carbon dioxide, from soil or biomass stock",
        "categories": ("air",),
        "type": "emission",
        "CAS number": "000124-38-9",
    },
    "co2_to_soil": {
        "name": "Carbon dioxide, to soil or biomass stock",
        "categories": ("soil",),
        "type": "emission",
        "CAS number": "000124-38-9",
    },
    "co2_in_air": {
        "name": "Carbon dioxide, in air",
        "categories": ("natural resource", "in air"),
        "type": "natural resource",
        "CAS number": "000124-38-9",
    },
    "ch4_fossil": {
        "name": "Methane, fossil",
        "categories": ("air", "non-urban air or from high stacks"),
        "type": "emission",
        "CAS number": "000074-82-8",
    },
}

# static CFs, roughly EF v3.1 GWP100; uptake flows are negative
CHARACTERIZATION_FACTORS = {
    "co2_fossil_high_stacks": 1.0,
    "co2_fossil_unspecified": 1.0,
    "co2_fossil_urban": 1.0,
    "co2_fossil_stratosphere": 1.0,
    "co2_non_fossil": 1.0,
    "co2_from_soil_or_biomass_stock": 1.0,
    "co2_to_soil": -1.0,
    "co2_in_air": -1.0,
    "ch4_fossil": 29.8,
}


def setup_project():
    """Write a minimal biosphere database and static method into the current project."""
    clear_characterization_function_cache()
    biosphere = bd.Database(bd.config.biosphere)
    biosphere.write(
        {
            (bd.config.biosphere, code): data | {"code": code}
            for code, data in BIOFLOWS.items()
        }
    )
    bd.Method(METHOD).write(
        [
            ((bd.config.biosphere, code), factor)
            for code, factor in CHARACTERIZATION_FACTORS.items()
        ]
    )
    return {code: biosphere.get(code).id for code in BIOFLOWS}


@pytest.mark.parametrize(
    "code",
    [
        "co2_fossil_high_stacks",
        "co2_fossil_unspecified",
        "co2_fossil_urban",
        "co2_fossil_stratosphere",
        "co2_non_fossil",
        "co2_from_soil_or_biomass_stock",
    ],
)
@bw2test
def test_co2_emission_flows_get_co2_function(code):
    ids = setup_project()
    functions = create_characterization_functions_from_method(METHOD)
    assert functions[ids[code]] is characterize_co2


@bw2test
def test_co2_emission_flows_characterized_without_uptake():
    """Disabling uptake must not drop the ordinary CO2 emission flows."""
    ids = setup_project()
    functions = create_characterization_functions_from_method(
        METHOD, characterize_uptake=False
    )
    assert functions[ids["co2_fossil_high_stacks"]] is characterize_co2
    assert functions[ids["co2_from_soil_or_biomass_stock"]] is characterize_co2


@bw2test
def test_uptake_flows_use_uptake_function():
    ids = setup_project()
    functions = create_characterization_functions_from_method(METHOD)
    assert functions[ids["co2_in_air"]] is characterize_co2_uptake
    assert functions[ids["co2_to_soil"]] is characterize_co2_uptake


@bw2test
def test_uptake_flows_absent_when_uptake_disabled():
    ids = setup_project()
    functions = create_characterization_functions_from_method(
        METHOD, characterize_uptake=False
    )
    assert ids["co2_in_air"] not in functions
    assert ids["co2_to_soil"] not in functions


@bw2test
def test_all_method_flows_are_characterized():
    ids = setup_project()
    functions = create_characterization_functions_from_method(METHOD)
    assert set(functions) == set(ids.values())
    assert functions[ids["ch4_fossil"]] is characterize_ch4


@bw2test
def test_dynamic_gwp100_matches_static_score():
    """
    End-to-end: with fixed_time_horizon=False, dynamic GWP100 of an inventory
    must match the static GWP100 score of the same inventory within a few percent.
    """
    ids = setup_project()

    amounts = {
        "co2_fossil_high_stacks": 17864.7,
        "co2_fossil_unspecified": 2322.3,
        "co2_fossil_urban": 1287.9,
        "co2_fossil_stratosphere": 1.9,
        "co2_from_soil_or_biomass_stock": 24.7,
        "co2_non_fossil": 100.0,
        "co2_to_soil": 50.0,
        "co2_in_air": 200.0,
        "ch4_fossil": 30.0,
    }

    dynamic_inventory_df = pd.DataFrame(
        {
            "date": pd.Series(
                ["2024-01-01"] * len(amounts), dtype="datetime64[s]"
            ),
            "amount": pd.Series(list(amounts.values()), dtype="float64"),
            "flow": pd.Series([ids[code] for code in amounts], dtype="int"),
            "activity": pd.Series([1] * len(amounts), dtype="int"),
        }
    )

    static_score = sum(
        amount * CHARACTERIZATION_FACTORS[code] for code, amount in amounts.items()
    )

    characterized = characterize(
        dynamic_inventory_df,
        metric="GWP",
        base_lcia_method=METHOD,
        time_horizon=100,
        fixed_time_horizon=False,
    )

    assert np.isclose(characterized.amount.sum(), static_score, rtol=0.05)
