"""Default characterization functions with characterization_options={"parameter_set": "ar6"}."""

import warnings

import bw2data as bd
import numpy as np
import pandas as pd
import pytest
from bw2data.tests import bw2test

import dynamic_characterization.dynamic_characterization as dc
from dynamic_characterization import CharacterizationOptions, characterize, compare_with_static
from dynamic_characterization.dynamic_characterization import (
    clear_characterization_function_cache,
    create_characterization_functions_from_method,
)
from dynamic_characterization.ipcc_ar6.radiative_forcing import (
    characterize_co2,
    characterize_co2_ar6,
    characterize_co2_uptake_ar6,
)

AR6 = {"parameter_set": "ar6"}
METHOD = ("test", "IPCC 2021", "GWP100")

# name, categories, type, CAS, static GWP100 (ecoinvent IPCC 2021 / AR6 Table 7.15 and 7.SM.7)
BIOFLOWS = {
    "co2_fossil": ("Carbon dioxide, fossil", ("air",), "emission", "000124-38-9", 1.0),
    "co2_non_fossil": ("Carbon dioxide, non-fossil", ("air",), "emission", "000124-38-9", 1.0),
    "co2_in_air": ("Carbon dioxide, in air", ("natural resource", "in air"), "natural resource", "000124-38-9", -1.0),
    "co2_to_soil": ("Carbon dioxide, to soil or biomass stock", ("soil",), "emission", "000124-38-9", -1.0),
    "ch4_fossil": ("Methane, fossil", ("air",), "emission", "000074-82-8", 29.8),
    "ch4_non_fossil": ("Methane, non-fossil", ("air",), "emission", "000074-82-8", 27.0),
    "ch4_soil": ("Methane, from soil or biomass stock", ("air",), "emission", "000074-82-8", 29.8),
    "ch4_unspecified": ("Methane", ("air",), "emission", "74-82-8", 29.8),
    "n2o": ("Dinitrogen monoxide", ("air",), "emission", "010024-97-2", 273.0),
    "cf4": ("Tetrafluoromethane", ("air",), "emission", "000075-73-0", 7380.0),
    "dichloromethane": ("Methane, dichloro-, HCC-30", ("air",), "emission", "000075-09-2", 11.2),
    "vinylidene_fluoride": ("Vinylidene fluoride", ("air",), "emission", "75-38-7", 0.052),
    "co_fossil": ("Carbon monoxide, fossil", ("air",), "emission", "000630-08-0", 4.0624),
    "water": ("Water", ("air",), "emission", "007732-18-5", 0.0),
}
AR6_MATCHED = [code for code in BIOFLOWS if code not in ("co_fossil", "water")]


def setup_project():
    clear_characterization_function_cache()
    biosphere = bd.Database(bd.config.biosphere)
    biosphere.write(
        {
            (bd.config.biosphere, code): {
                "code": code,
                "name": name,
                "categories": categories,
                "type": type_,
                "CAS number": cas,
            }
            for code, (name, categories, type_, cas, _) in BIOFLOWS.items()
        }
    )
    bd.Method(METHOD).write(
        [((bd.config.biosphere, code), cf) for code, (*_, cf) in BIOFLOWS.items() if cf]
    )
    return {code: biosphere.get(code).id for code in BIOFLOWS}


@bw2test
def test_compare_with_static_matches_ar6_factors():
    ids = setup_project()
    table = compare_with_static(METHOD, 100, AR6).set_index("flow")
    for code in AR6_MATCHED:
        static = BIOFLOWS[code][-1]
        assert table.loc[ids[code], "dynamic_gwp"] == pytest.approx(static, rel=0.01, abs=1e-3), code


@bw2test
def test_compare_with_static_legacy_reproduces_known_gap():
    ids = setup_project()
    table = compare_with_static(METHOD, 100, {"parameter_set": "legacy"}).set_index("flow")
    assert table.loc[ids["ch4_fossil"], "dynamic_gwp"] == pytest.approx(26.71, abs=0.01)
    assert table.loc[ids["n2o"], "dynamic_gwp"] == pytest.approx(263.6, abs=0.1)


@bw2test
def test_function_assignment():
    ids = setup_project()
    functions = create_characterization_functions_from_method(METHOD, characterization_options=AR6)
    assert functions[ids["co2_fossil"]] is characterize_co2_ar6
    assert functions[ids["co2_non_fossil"]] is characterize_co2_ar6
    assert functions[ids["co2_in_air"]] is characterize_co2_uptake_ar6
    assert functions[ids["co2_to_soil"]] is characterize_co2_uptake_ar6
    assert functions[ids["ch4_soil"]] is functions[ids["ch4_fossil"]]
    assert functions[ids["ch4_unspecified"]] is functions[ids["ch4_fossil"]]
    assert functions[ids["ch4_non_fossil"]] is not functions[ids["ch4_fossil"]]
    assert functions[ids["dichloromethane"]].spec.name == "Methylene chloride"  # CAS 75-09-2, not CH4
    assert ids["water"] not in functions


@bw2test
def test_uptake_can_be_disabled():
    ids = setup_project()
    functions = create_characterization_functions_from_method(
        METHOD, characterize_uptake=False, characterization_options=AR6
    )
    assert ids["co2_in_air"] not in functions
    assert ids["co2_to_soil"] not in functions
    assert functions[ids["co2_fossil"]] is characterize_co2_ar6


@bw2test
def test_unspecified_origin_and_overrides():
    ids = setup_project()
    options = dict(
        AR6,
        unspecified_origin="non-fossil",
        flow_overrides={(bd.config.biosphere, "ch4_soil"): {"origin": "non-fossil"}, ids["water"]: {"species": None}},
    )
    table = compare_with_static(METHOD, 100, options).set_index("flow")
    assert table.loc[ids["ch4_unspecified"], "dynamic_gwp"] == pytest.approx(27.0, rel=0.01)
    assert table.loc[ids["ch4_soil"], "dynamic_gwp"] == pytest.approx(27.0, rel=0.01)
    assert table.loc[ids["ch4_fossil"], "dynamic_gwp"] == pytest.approx(29.8, rel=0.01)


@bw2test
def test_options_change_the_functions():
    ids = setup_project()
    plain = create_characterization_functions_from_method(METHOD, characterization_options=AR6)
    no_ccf = create_characterization_functions_from_method(
        METHOD, characterization_options=dict(AR6, climate_carbon_feedback=False)
    )
    legacy = create_characterization_functions_from_method(
        METHOD, characterization_options={"parameter_set": "legacy"}
    )
    assert plain[ids["n2o"]] is not no_ccf[ids["n2o"]]
    assert plain[ids["n2o"]].agwp(100) > no_ccf[ids["n2o"]].agwp(100)
    assert legacy[ids["co2_fossil"]] is characterize_co2


def inventory(ids):
    amounts = {"co2_fossil": 1000.0, "co2_in_air": 300.0, "ch4_fossil": 3.0, "ch4_non_fossil": 2.0, "n2o": 0.5, "cf4": 0.001}
    return amounts, pd.DataFrame(
        {
            "date": pd.Series(["2024-01-01"] * len(amounts), dtype="datetime64[s]"),
            "amount": pd.Series(list(amounts.values()), dtype="float64"),
            "flow": pd.Series([ids[code] for code in amounts], dtype="int"),
            "activity": pd.Series([1] * len(amounts), dtype="int"),
        }
    )


@bw2test
def test_dynamic_gwp100_matches_static_ar6_score():
    ids = setup_project()
    amounts, df = inventory(ids)
    static = sum(amount * BIOFLOWS[code][-1] for code, amount in amounts.items())
    result = characterize(df, metric="GWP", base_lcia_method=METHOD, time_horizon=100, characterization_options=AR6)
    assert result.amount.sum() == pytest.approx(static, rel=0.005)


@bw2test
def test_vectorized_radiative_forcing_equals_row_loop():
    ids = setup_project()
    _, df = inventory(ids)
    fast = characterize(df, base_lcia_method=METHOD, time_horizon=50, characterization_options=AR6)
    # a fixed time horizon starting at the emission date gives the same horizon but takes the row loop
    loop = characterize(
        df,
        base_lcia_method=METHOD,
        time_horizon=50,
        fixed_time_horizon=True,
        time_horizon_start=pd.Timestamp("2024-01-01"),
        characterization_options=AR6,
    )
    pd.testing.assert_frame_equal(
        fast.sort_values(["flow", "date"]).reset_index(drop=True),
        loop.sort_values(["flow", "date"]).reset_index(drop=True),
    )


@bw2test
def test_prospective_metrics_are_not_yet_available_with_ar6():
    setup_project()
    _, df = inventory(setup_project())
    with pytest.raises(NotImplementedError):
        characterize(df, metric="pGWP", base_lcia_method=METHOD, characterization_options=AR6)


@bw2test
def test_future_warning_once_without_options(monkeypatch):
    setup_project()
    monkeypatch.setattr(dc, "_DEFAULT_PARAMETER_SET_WARNED", False)
    with pytest.warns(FutureWarning, match="parameter_set"):
        create_characterization_functions_from_method(METHOD)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        clear_characterization_function_cache()
        create_characterization_functions_from_method(METHOD)  # only once per process
        create_characterization_functions_from_method(METHOD, characterization_options={"parameter_set": "legacy"})


def test_option_validation():
    with pytest.raises(ValueError, match="Unknown characterization option"):
        CharacterizationOptions.from_value({"parameterset": "ar6"})
    with pytest.raises(ValueError):
        CharacterizationOptions.from_value({"parameter_set": "ar7"})
    with pytest.raises(ValueError):
        CharacterizationOptions.from_value({"ch4_oxidation_yield": 1.5})
    with pytest.raises(NotImplementedError):
        CharacterizationOptions.from_value({"parameter_set": "fair_irf"})
    assert CharacterizationOptions.from_value(AR6) == CharacterizationOptions(parameter_set="ar6")
    assert hash(CharacterizationOptions.from_value(dict(AR6, flow_overrides={1: {"origin": "fossil"}})))
