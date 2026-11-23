"""
Reproduction of Barbosa Watanabe & Cherubini (2026), doi:10.1021/acs.est.5c12391:
``prospective.watanabe2026`` and the default "watanabe2026" method of the "ar6" parameter set
with the prospective metrics.

Imports must stay at module level and this module must sort before test_prospective.py, which
replaces `dynamic_characterization` in sys.modules when it is collected.
"""

import os
import re
import warnings

import bw2data as bd
import numpy as np
import pandas as pd
import pytest
from bw2data.tests import bw2test

from dynamic_characterization import ar6_prospective, characterize
from dynamic_characterization.dynamic_characterization import (
    clear_characterization_function_cache,
    create_characterization_functions_from_method,
)
from dynamic_characterization.ipcc_ar6 import ar6_core, gases
from dynamic_characterization.prospective import set_scenario, watanabe2026

DATA = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "dynamic_characterization", "prospective", "data"
)
# the SI tables are rounded to one decimal; N2O deviates by up to 3e-5 relative (0.06 absolute),
# the size of the float differences between copies of the input files
TOLERANCE = 0.06


def _si_table(number):
    df = pd.read_excel(os.path.join(DATA, f"es5c12391_si_{number}.xlsx"))
    df = df.rename(columns={df.columns[0]: "Scenario"})
    cells = []
    for _, row in df.iterrows():
        scenario = re.match(r"^([A-Z0-9]+)(?:-[A-Z]+)?_-_?(SSP\d)_-_([\d.]+)$", row["Scenario"].strip()).groups()
        for year in df.columns[1:]:
            if not pd.isna(row[year]):
                cells.append((scenario, int(year), float(row[year])))
    return cells


@pytest.mark.parametrize(
    "number, gas, metric",
    [("001", "CH4", "pgtp"), ("002", "N2O", "pgtp"), ("003", "CH4", "pgwp"), ("004", "N2O", "pgwp")],
)
def test_si_tables_100_year_horizon(number, gas, metric):
    """Every cell of es5c12391_si_001-004.xlsx (18 scenarios x 2030-2050)."""
    cells = _si_table(number)
    assert len(cells) == 378
    deviations = [abs(getattr(watanabe2026, metric)(gas, year, 100, *scenario) - value) for scenario, year, value in cells]
    assert max(deviations) <= TOLERANCE


# SI Tables S.4, S.5, S.8 and S.9 (in es5c12391_si_028.pdf only): GWP20 and GTP50
PDF_TABLES = [
    # gas, metric, horizon, scenario, {year: value}
    ("CH4", "pgwp", 20, ("AIM", "SSP3", "4.5"), {2030: 91.1, 2040: 95.9, 2050: 102.3}),
    ("CH4", "pgwp", 20, ("IMAGE", "SSP1", "8.5"), {2030: 92.7, 2040: 98.0, 2050: 103.8}),
    ("CH4", "pgwp", 20, ("MESSAGE", "SSP2", "2.6"), {2030: 92.0, 2040: 97.0, 2050: 102.1}),
    ("N2O", "pgwp", 20, ("AIM", "SSP3", "4.5"), {2030: 280.8, 2040: 297.1, 2050: 308.7}),
    ("N2O", "pgwp", 20, ("IMAGE", "SSP1", "8.5"), {2030: 278.1, 2040: 292.3, 2050: 306.3}),
    ("N2O", "pgwp", 20, ("MESSAGE", "SSP2", "2.6"), {2030: 273.1, 2040: 279.3, 2050: 280.4}),
    ("CH4", "pgtp", 50, ("AIM", "SSP3", "4.5"), {2030: 10.5, 2040: 11.0, 2050: 11.7}),
    ("CH4", "pgtp", 50, ("IMAGE", "SSP1", "8.5"), {2030: 10.7, 2040: 11.3, 2050: 11.9}),
    ("CH4", "pgtp", 50, ("MESSAGE", "SSP2", "2.6"), {2030: 10.6, 2040: 11.2, 2050: 11.7}),
    ("N2O", "pgtp", 50, ("AIM", "SSP3", "4.5"), {2030: 285.6, 2040: 301.1, 2050: 311.8}),
    ("N2O", "pgtp", 50, ("IMAGE", "SSP1", "8.5"), {2030: 282.9, 2040: 296.4, 2050: 309.5}),
    ("N2O", "pgtp", 50, ("MESSAGE", "SSP2", "2.6"), {2030: 278.1, 2040: 283.7, 2050: 284.4}),
]


@pytest.mark.parametrize("gas, metric, horizon, scenario, values", PDF_TABLES)
def test_si_tables_other_horizons(gas, metric, horizon, scenario, values):
    for year, value in values.items():
        assert getattr(watanabe2026, metric)(gas, year, horizon, *scenario) == pytest.approx(value, abs=TOLERANCE)


def test_co2_irf_selection_of_the_paper_code():
    """The paper's code picks the CO2 impulse response by the scenario's position."""
    assert watanabe2026.co2_irf_key("MESSAGE", "SSP2", "8.5") == "RCP85"
    assert watanabe2026.co2_irf_key("IMAGE", "SSP1", "8.5") == "RCP60"
    assert watanabe2026.co2_irf_key("AIM", "SSP3", "4.5") == "RCP26"
    assert watanabe2026.co2_irf_key("IMAGE", "SSP1", "8.5", paper_code=False) == "RCP85"


def test_corrected_variant_differs():
    paper = watanabe2026.pgtp("CH4", 2030, 100, "MESSAGE", "SSP2", "4.5")
    corrected = watanabe2026.pgtp("CH4", 2030, 100, "MESSAGE", "SSP2", "4.5", paper_code=False)
    assert corrected > 1.2 * paper  # feedback boxes decaying from their own values
    # the CO2 response is the same for MESSAGE (all four RCPs), so CO2 itself does not change
    a = watanabe2026.scenario_series("CO2", "agwp", "MESSAGE", "SSP2", "4.5", 2030, True)
    b = watanabe2026.scenario_series("CO2", "agwp", "MESSAGE", "SSP2", "4.5", 2030, False)
    assert np.array_equal(a, b)


def test_tabulated_co2_response_equals_paper_co2():
    """Gases outside the paper use the same CO2 response, so the pGWP reference is consistent."""
    response = ar6_prospective.scenario_co2_response("IMAGE", "SSP1", "8.5", 2040, method="watanabe2026")
    paper = watanabe2026.scenario_series("CO2", "agwp", "IMAGE", "SSP1", "8.5", 2040)
    for horizon in (1, 20, 50, 100):
        assert float(ar6_core.agwp(gases.co2(), float(horizon), co2_response=response)) == pytest.approx(
            paper[horizon - 1], rel=1e-9
        )


def test_years_outside_the_data_are_clamped():
    with pytest.warns(UserWarning, match="clamping"):
        early = watanabe2026.pgwp("CH4", 2010, 100, "MESSAGE", "SSP2", "4.5")
    assert early == watanabe2026.pgwp("CH4", 2020, 100, "MESSAGE", "SSP2", "4.5")


# --- through characterize() with the "ar6" parameter set ------------------------------------

AR6 = {"parameter_set": "ar6"}
METHOD = ("test", "IPCC 2021", "prospective watanabe2026")
BIOFLOWS = {
    "co2": ("Carbon dioxide, fossil", ("air",), "emission", "000124-38-9"),
    "co2_in_air": ("Carbon dioxide, in air", ("natural resource", "in air"), "natural resource", "000124-38-9"),
    "ch4": ("Methane, fossil", ("air",), "emission", "000074-82-8"),
    "ch4_non_fossil": ("Methane, non-fossil", ("air",), "emission", "000074-82-8"),
    "n2o": ("Dinitrogen monoxide", ("air",), "emission", "010024-97-2"),
    "cf4": ("Tetrafluoromethane", ("air",), "emission", "000075-73-0"),
}


def setup_project():
    clear_characterization_function_cache()
    biosphere = bd.Database(bd.config.biosphere)
    biosphere.write(
        {
            (bd.config.biosphere, code): {"code": code, "name": name, "categories": cats, "type": type_, "CAS number": cas}
            for code, (name, cats, type_, cas) in BIOFLOWS.items()
        }
    )
    bd.Method(METHOD).write([((bd.config.biosphere, code), 1.0) for code in BIOFLOWS])
    return {code: biosphere.get(code).id for code in BIOFLOWS}


def inventory(ids, year):
    return pd.DataFrame(
        {
            "date": pd.Series([f"{year}-01-01"] * len(ids), dtype="datetime64[s]"),
            "amount": pd.Series([3.0] * len(ids), dtype="float64"),
            "flow": pd.Series(list(ids.values()), dtype="int"),
            "activity": pd.Series([1] * len(ids), dtype="int"),
        }
    )


@bw2test
@pytest.mark.parametrize("metric, table, gas", [("pGWP", "003", "CH4"), ("pGWP", "004", "N2O"), ("pGTP", "001", "CH4"), ("pGTP", "002", "N2O")])
def test_characterize_reproduces_the_si_tables(metric, table, gas):
    ids = setup_project()
    cells = {(sc, year): value for sc, year, value in _si_table(table)}
    for scenario, year in ((("MESSAGE", "SSP2", "4.5"), 2030), (("IMAGE", "SSP1", "8.5"), 2045), (("AIM", "SSP3", "6.0"), 2050)):
        set_scenario(*scenario)
        result = characterize(inventory(ids, year), metric=metric, base_lcia_method=METHOD, time_horizon=100, characterization_options=AR6)
        by_flow = result.groupby("flow").amount.sum() / 3.0
        code = gas.lower()
        assert by_flow[ids[code]] == pytest.approx(cells[(scenario, year)], abs=TOLERANCE)
        assert by_flow[ids["co2"]] == pytest.approx(1.0, rel=1e-12)
        assert by_flow[ids["co2_in_air"]] == pytest.approx(-1.0, rel=1e-12)
        if gas == "CH4":  # the paper does not distinguish fossil and non-fossil CH4
            assert by_flow[ids["ch4_non_fossil"]] == pytest.approx(by_flow[ids["ch4"]])


@bw2test
def test_prospective_radiative_forcing_adds_up_to_paper_agwp():
    ids = setup_project()
    set_scenario("REMIND", "SSP5", "8.5")
    result = characterize(
        inventory(ids, 2040), metric="prospective_radiative_forcing", base_lcia_method=METHOD,
        time_horizon=100, characterization_options=AR6,
    )
    by_flow = result.groupby("flow").amount.sum()
    for gas in ("CO2", "CH4", "N2O"):
        paper = watanabe2026.scenario_series(gas, "agwp", "REMIND", "SSP5", "8.5", 2040)
        assert by_flow[ids[gas.lower()]] == pytest.approx(3.0 * paper[99], rel=1e-12)


@bw2test
def test_other_gases_use_the_paper_co2_reference():
    ids = setup_project()
    set_scenario("IMAGE", "SSP1", "8.5")
    functions = create_characterization_functions_from_method(METHOD, use_prospective=True, characterization_options=AR6)
    result = characterize(inventory(ids, 2030), metric="pGWP", base_lcia_method=METHOD, characterization_options=AR6)
    cf4 = result.groupby("flow").amount.sum()[ids["cf4"]] / 3.0
    paper_co2 = watanabe2026.scenario_series("CO2", "agwp", "IMAGE", "SSP1", "8.5", 2030)[99]
    assert cf4 == pytest.approx(functions[ids["cf4"]].agwp_for(pd.Timestamp("2030-01-01"), 100) / paper_co2, rel=1e-9)
    assert 6000 < cf4 < 9000


@bw2test
def test_methods_and_paper_code_options():
    ids = setup_project()
    set_scenario("IMAGE", "SSP1", "8.5")
    values = {}
    for label, options in {
        "paper": AR6,
        "corrected": dict(AR6, prospective_paper_code=False),
        "ar6": dict(AR6, prospective_method="ar6"),
    }.items():
        result = characterize(inventory(ids, 2030), metric="pGWP", base_lcia_method=METHOD, characterization_options=options)
        values[label] = result.groupby("flow").amount.sum()[ids["ch4"]] / 3.0
    assert values["paper"] == pytest.approx(28.2, abs=TOLERANCE)  # SI table 003, IMAGE-SSP1-8.5, 2030
    assert len({round(v, 3) for v in values.values()}) == 3
