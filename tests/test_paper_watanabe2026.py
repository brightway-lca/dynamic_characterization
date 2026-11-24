"""
Barbosa Watanabe & Cherubini (2026), doi:10.1021/acs.est.5c12391, with two errors of the
published code corrected: ``prospective.watanabe2026`` and the default "watanabe2026" method of
the "ar6" parameter set with the prospective metrics.

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
from dynamic_characterization.options import CharacterizationOptions
from dynamic_characterization.prospective import set_scenario, watanabe2026
from dynamic_characterization.prospective.data_loader import load_irf_co2, load_re_co2

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


UNAFFECTED = {("GCAM4", "SSP4"), ("MESSAGE", "SSP2"), ("REMIND", "SSP5")}  # all four RCPs: CO2 IRF selection unaffected


def test_matches_the_si_tables_where_the_corrections_do_not_apply():
    """
    N2O pGTP (Code #9) has no feedback-box error, and GCAM4, MESSAGE and REMIND are not affected
    by the CO2 impulse-response selection: there the corrected port equals the SI table
    (es5c12391_si_002.xlsx), cell by cell.
    """
    cells = [(sc, year, value) for sc, year, value in _si_table("002") if sc[:2] in UNAFFECTED]
    assert len(cells) == 12 * 21
    for scenario, year, value in cells:
        assert watanabe2026.pgtp("N2O", year, 100, *scenario) == pytest.approx(value, abs=TOLERANCE)


@pytest.mark.parametrize("year, value", [(2030, 278.1), (2040, 283.7), (2050, 284.4)])
def test_matches_si_table_s9_where_the_corrections_do_not_apply(year, value):
    """N2O pGTP50, MESSAGE-SSP2-2.6 (SI Table S.9)."""
    assert watanabe2026.pgtp("N2O", year, 50, "MESSAGE", "SSP2", "2.6") == pytest.approx(value, abs=TOLERANCE)


@pytest.mark.parametrize("number, metric", [("003", "pgwp"), ("001", "pgtp")])
def test_corrected_feedback_raises_ch4_above_the_si_tables(number, metric):
    """
    In the published CH4 code, boxes 2 and 3 decay from box 1 (lifetime 2.4 yr), which loses most
    of the slow feedback CO2. Corrected, CH4 pGWP100 and pGTP100 are higher than the tables in
    every scenario that the CO2 impulse-response selection does not affect.
    """
    cells = [(sc, year, value) for sc, year, value in _si_table(number) if sc[:2] in UNAFFECTED]
    ratios = [getattr(watanabe2026, metric)("CH4", year, 100, *sc) / value for sc, year, value in cells]
    assert min(ratios) > 1.0
    assert max(ratios) < (1.10 if metric == "pgwp" else 1.40)


def test_feedback_boxes_decay_from_their_own_values():
    boxes = watanabe2026._ccf_boxes(0.0, (1.0, 2.0, 3.0))
    expected = tuple(c * np.exp(-1 / tau) for c, tau in zip((1.0, 2.0, 3.0), (2.376, 30.14, 490.1)))
    assert np.allclose(boxes, expected, rtol=1e-15)


def test_co2_impulse_response_of_the_scenarios_rcp():
    """IMAGE-SSP1-8.5 uses the RCP8.5 response (the published code used RCP6.0)."""
    assert watanabe2026.co2_irf_key("8.5") == "RCP85"
    re = load_re_co2("IMAGE")
    re_2040 = float(re["SSP1"]["8.5"][list(re["_years"]).index(2040)])
    expected = watanabe2026.agwp_co2_series(re_2040, load_irf_co2()["RCP85"])
    assert np.array_equal(watanabe2026.scenario_series("CO2", "agwp", "IMAGE", "SSP1", "8.5", 2040), expected)


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
@pytest.mark.parametrize("metric, gas", [("pGWP", "CH4"), ("pGWP", "N2O"), ("pGTP", "CH4"), ("pGTP", "N2O")])
def test_characterize_gives_the_watanabe2026_metrics(metric, gas):
    ids = setup_project()
    port = watanabe2026.pgwp if metric == "pGWP" else watanabe2026.pgtp
    for scenario, year in ((("MESSAGE", "SSP2", "4.5"), 2030), (("IMAGE", "SSP1", "8.5"), 2045), (("AIM", "SSP3", "6.0"), 2050)):
        set_scenario(*scenario)
        result = characterize(inventory(ids, year), metric=metric, base_lcia_method=METHOD, time_horizon=100, characterization_options=AR6)
        by_flow = result.groupby("flow").amount.sum() / 3.0
        assert by_flow[ids[gas.lower()]] == pytest.approx(port(gas, year, 100, *scenario), rel=1e-9)
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
def test_the_two_methods_differ():
    ids = setup_project()
    set_scenario("IMAGE", "SSP1", "8.5")
    values = {}
    for method in ("watanabe2026", "ar6"):
        options = dict(AR6, prospective_method=method)
        result = characterize(inventory(ids, 2030), metric="pGWP", base_lcia_method=METHOD, characterization_options=options)
        values[method] = result.groupby("flow").amount.sum()[ids["ch4"]] / 3.0
    assert values["watanabe2026"] == pytest.approx(watanabe2026.pgwp("CH4", 2030, 100, "IMAGE", "SSP1", "8.5"), rel=1e-9)
    assert abs(values["watanabe2026"] / values["ar6"] - 1) > 0.02


def test_paper_code_option_is_gone():
    with pytest.raises(ValueError, match="Unknown characterization option"):
        CharacterizationOptions.from_value({"parameter_set": "ar6", "prospective_paper_code": True})
