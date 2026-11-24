"""The "ar6" parameter set with a prospective scenario, method "ar6" (AR6 recipe with scenario inputs).

The default method "watanabe2026" is tested in test_paper_watanabe2026.py.

"""

import bw2data as bd
import numpy as np
import pandas as pd
import pytest
from bw2data.tests import bw2test

from dynamic_characterization import ar6_prospective, characterize
from dynamic_characterization.ar6_dispatch import FlowClass, function_for
from dynamic_characterization.dynamic_characterization import (
    clear_characterization_function_cache,
    create_characterization_functions_from_method,
)
from dynamic_characterization.ipcc_ar6 import ar6_core, gases
from dynamic_characterization.options import CharacterizationOptions
from dynamic_characterization.prospective import set_scenario

CCF = ar6_core.ccf_parameters()
AR6 = {"parameter_set": "ar6", "prospective_method": "ar6"}
OPTIONS = CharacterizationOptions(parameter_set="ar6", prospective_method="ar6")
DATE = pd.Timestamp("2030-01-01")


def prospective_function(species, origin=None, options=OPTIONS):
    return function_for(FlowClass(species, origin), options, prospective=True)


# --- the tabulated CO2 response reproduces the analytic AR6 one ---------------------------

AR6_SAMPLED = ar6_core.TabulatedCO2Response(
    "AR6, sampled annually", tuple(float(v) for v in ar6_core.radiative_forcing(gases.co2(), np.arange(101.0)))
)


@pytest.mark.parametrize(
    "spec",
    [
        gases.co2(),
        gases.ch4("fossil", ccf=CCF),
        gases.ch4("non-fossil", ccf=CCF),
        gases.n2o(ccf=CCF),
        gases.halocarbon("75-73-0", ccf=CCF),
        gases.co("fossil"),
    ],
    ids=lambda spec: spec.name,
)
@pytest.mark.parametrize("horizon", [20.0, 50.0, 100.0])
def test_tabulated_co2_response_matches_analytic(spec, horizon):
    """Only the CO2 response differs between the two modes, so with the AR6 response sampled
    at annual steps (like the scenario data) the results must agree with the analytic ones."""
    assert float(ar6_core.agwp(spec, horizon, co2_response=AR6_SAMPLED)) == pytest.approx(
        float(ar6_core.agwp(spec, horizon)), rel=5e-4
    )
    assert float(ar6_core.agtp(spec, horizon, co2_response=AR6_SAMPLED)) == pytest.approx(
        float(ar6_core.agtp(spec, horizon)), rel=5e-4
    )


def test_tabulated_response_rejects_longer_horizons():
    with pytest.raises(ValueError, match="covers 100 years"):
        ar6_core.agwp(gases.co2(), 101.0, co2_response=AR6_SAMPLED)


def test_scenario_radiative_efficiency_gets_ar6_adjustments():
    """With the AR6 (unadjusted) radiative efficiencies as input, the builders give the AR6 gases."""
    p = ar6_core.load_parameters()
    sarf_ch4 = p["ch4"]["radiative_efficiency_direct_W_m2_ppb"] / (1 + p["ch4"]["tropospheric_adjustment"])
    sarf_n2o = p["n2o"]["radiative_efficiency_direct_W_m2_ppb"] / (1 + p["n2o"]["tropospheric_adjustment"])
    for horizon in (20.0, 100.0):
        assert float(ar6_core.agwp(gases.ch4("fossil", ccf=CCF, direct_radiative_efficiency=sarf_ch4), horizon)) == pytest.approx(
            float(ar6_core.agwp(gases.ch4("fossil", ccf=CCF), horizon)), rel=1e-12
        )
        assert float(
            ar6_core.agwp(gases.n2o(CCF, direct_radiative_efficiency=sarf_n2o, ch4_direct_radiative_efficiency=sarf_ch4), horizon)
        ) == pytest.approx(float(ar6_core.agwp(gases.n2o(CCF), horizon)), rel=1e-12)


# --- scenario functions ---------------------------------------------------------------------


@pytest.fixture
def scenario():
    set_scenario(iam="MESSAGE", ssp="SSP2", rcp="4.5")
    return ("MESSAGE", "SSP2", "4.5")


def test_ch4_and_n2o_follow_the_scenario(scenario):
    response = ar6_prospective.scenario_co2_response(*scenario, 2030)
    re_ch4 = ar6_prospective.scenario_radiative_efficiency("CH4", *scenario, 2030)
    re_n2o = ar6_prospective.scenario_radiative_efficiency("N2O", *scenario, 2030)
    ch4 = gases.ch4("fossil", ccf=CCF, direct_radiative_efficiency=re_ch4)
    n2o = gases.n2o(CCF, direct_radiative_efficiency=re_n2o, ch4_direct_radiative_efficiency=re_ch4)
    assert prospective_function("CH4", "fossil").agwp_for(DATE, 100) == pytest.approx(
        float(ar6_core.agwp(ch4, 100.0, co2_response=response)), rel=1e-12
    )
    assert prospective_function("N2O").agwp_for(DATE, 100) == pytest.approx(
        float(ar6_core.agwp(n2o, 100.0, co2_response=response)), rel=1e-12
    )


def test_gases_without_co2_terms_do_not_depend_on_the_scenario(scenario):
    options = CharacterizationOptions(parameter_set="ar6", prospective_method="ar6", climate_carbon_feedback=False)
    cf4 = prospective_function("75-73-0", options=options)
    assert cf4.agwp_for(DATE, 100) == pytest.approx(float(ar6_core.agwp(gases.halocarbon("75-73-0"), 100.0)), rel=1e-12)


def test_near_term_values_are_close_to_ar6():
    """In a low-emission scenario around 2030 the background is close to AR6's (2019)."""
    set_scenario(iam="IMAGE", ssp="SSP1", rcp="2.6")
    co2 = ar6_prospective.co2_function("ar6")
    ch4 = prospective_function("CH4", "fossil").agwp_for(DATE, 100) / co2.agwp_for(DATE, 100)
    n2o = prospective_function("N2O").agwp_for(DATE, 100) / co2.agwp_for(DATE, 100)
    assert ch4 == pytest.approx(29.8, rel=0.05)
    assert n2o == pytest.approx(273, rel=0.05)


def test_functions_are_memoized():
    assert prospective_function("CH4", "fossil") is prospective_function("CH4", "fossil")
    assert prospective_function("CH4", "fossil") is not prospective_function("CH4", "non-fossil")


# --- through characterize() ---------------------------------------------------------------

METHOD = ("test", "IPCC 2021", "prospective")
BIOFLOWS = {
    "co2": ("Carbon dioxide, fossil", ("air",), "emission", "000124-38-9"),
    "co2_in_air": ("Carbon dioxide, in air", ("natural resource", "in air"), "natural resource", "000124-38-9"),
    "ch4": ("Methane, fossil", ("air",), "emission", "000074-82-8"),
    "n2o": ("Dinitrogen monoxide", ("air",), "emission", "010024-97-2"),
    "cf4": ("Tetrafluoromethane", ("air",), "emission", "000075-73-0"),
    "h2": ("Hydrogen", ("air",), "emission", "001333-74-0"),
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


def inventory(ids, year=2030):
    return pd.DataFrame(
        {
            "date": pd.Series([f"{year}-01-01"] * len(ids), dtype="datetime64[s]"),
            "amount": pd.Series([2.0] * len(ids), dtype="float64"),
            "flow": pd.Series(list(ids.values()), dtype="int"),
            "activity": pd.Series([1] * len(ids), dtype="int"),
        }
    )


@bw2test
@pytest.mark.parametrize("metric, method", [("pGWP", "agwp_for"), ("pGTP", "agtp_for")])
def test_prospective_metrics(metric, method):
    ids = setup_project()
    set_scenario(iam="MESSAGE", ssp="SSP2", rcp="4.5")
    result = characterize(inventory(ids), metric=metric, base_lcia_method=METHOD, time_horizon=100, characterization_options=AR6)
    by_flow = result.groupby("flow").amount.sum()
    assert by_flow[ids["co2"]] == pytest.approx(2.0, rel=1e-12)
    assert by_flow[ids["co2_in_air"]] == pytest.approx(-2.0, rel=1e-12)
    functions = create_characterization_functions_from_method(METHOD, use_prospective=True, characterization_options=AR6)
    reference = getattr(ar6_prospective.co2_function("ar6"), method)(DATE, 100)
    for code in ("ch4", "n2o", "cf4", "h2"):
        expected = 2.0 * getattr(functions[ids[code]], method)(DATE, 100) / reference
        assert by_flow[ids[code]] == pytest.approx(expected, rel=1e-9), code


@bw2test
def test_prospective_radiative_forcing_adds_up_to_agwp():
    ids = setup_project()
    set_scenario(iam="REMIND", ssp="SSP5", rcp="8.5")
    result = characterize(
        inventory(ids, year=2050), metric="prospective_radiative_forcing", base_lcia_method=METHOD,
        time_horizon=100, characterization_options=AR6,
    )
    functions = create_characterization_functions_from_method(METHOD, use_prospective=True, characterization_options=AR6)
    by_flow = result.groupby("flow").amount.sum()
    for code in ("co2", "ch4", "n2o", "h2"):
        assert by_flow[ids[code]] == pytest.approx(2.0 * functions[ids[code]].agwp_for(pd.Timestamp("2050-01-01"), 100), rel=1e-9)


@bw2test
def test_results_depend_on_scenario_and_year():
    ids = setup_project()
    values = []
    for scenario, year in ((("IMAGE", "SSP1", "2.6"), 2030), (("REMIND", "SSP5", "8.5"), 2030), (("REMIND", "SSP5", "8.5"), 2070)):
        set_scenario(*scenario)
        result = characterize(inventory(ids, year), metric="pGWP", base_lcia_method=METHOD, characterization_options=AR6)
        values.append(result.groupby("flow").amount.sum()[ids["n2o"]])
    assert len(set(np.round(values, 6))) == 3


@bw2test
def test_fallback_to_ipcc_false_keeps_only_scenario_gases():
    ids = setup_project()
    functions = create_characterization_functions_from_method(
        METHOD, use_prospective=True, fallback_to_ipcc=False, characterization_options=AR6
    )
    assert set(functions) == {ids["co2"], ids["co2_in_air"], ids["ch4"], ids["n2o"]}


@bw2test
def test_time_varying_re_and_long_horizons():
    """The "ar6" method ignores time_varying_re (with a warning); horizons > 100 yr are truncated."""
    ids = setup_project()
    set_scenario(iam="MESSAGE", ssp="SSP2", rcp="4.5")
    plain = characterize(inventory(ids), metric="pGWP", base_lcia_method=METHOD, characterization_options=AR6)
    ar6_prospective._WARNED.clear()
    with pytest.warns(UserWarning, match="time_varying_re"):
        varying = characterize(
            inventory(ids), metric="pGWP", base_lcia_method=METHOD, time_varying_re=True, characterization_options=AR6
        )
    pd.testing.assert_frame_equal(plain, varying)
    with pytest.warns(UserWarning, match="truncated"):
        long = characterize(inventory(ids), metric="pGWP", base_lcia_method=METHOD, time_horizon=150, characterization_options=AR6)
    pd.testing.assert_frame_equal(plain, long)
