"""Indirect forcing of CO, VOC, H2 and NOx: the parameter sets reproduce their sources."""

import bw2data as bd
import numpy as np
import pandas as pd
import pytest
from bw2data.tests import bw2test

from dynamic_characterization import CharacterizationOptions, characterize, compare_with_static
from dynamic_characterization.dynamic_characterization import (
    clear_characterization_function_cache,
    create_characterization_functions_from_method,
)
from dynamic_characterization.ipcc_ar6 import ar6_core, gases

CO2 = gases.co2()
DATA = gases.load_precursor_parameters()


def _agwp_co2(published, horizon):
    reference = published["agwp_co2"]
    return float(ar6_core.agwp(CO2, horizon)) if reference == "ar6" else reference[str(horizon)]


def _gwp(spec, published, horizon):
    return float(ar6_core.agwp(spec, float(horizon))) / _agwp_co2(published, horizon)


def _rounding(value):
    """Half a unit in the last digit given (the sources give 1-3 significant digits)."""
    text = repr(float(value)).rstrip("0").rstrip(".")
    decimals = len(text.split(".")[1]) if "." in text else 0
    return 0.5 * 10.0**-decimals


PUBLISHED_TOTALS = [
    (species, name)
    for species, sets in DATA["sets"].items()
    for name, p in sets.items()
    if p["published"] and not (species == "NOx")  # NOx: see test_fry2012_nox
]


@pytest.mark.parametrize("species, name", PUBLISHED_TOTALS, ids=[f"{s}-{n}" for s, n in PUBLISHED_TOTALS])
def test_parameter_set_reproduces_published_total(species, name):
    published = DATA["sets"][species][name]["published"]
    spec = gases.precursor(species, name)
    n_components = len(DATA["sets"][species][name]["components"])
    for horizon, value in zip(published["horizons"], published["total"]):
        computed = _gwp(spec, published, horizon)
        # Fry et al. 2012 round each component to one decimal, so the total carries their rounding
        rounding = n_components * _rounding(value) if "components" in published else _rounding(value)
        tolerance = max(0.03 * abs(value), rounding)
        tolerance = max(tolerance, 0.05 * abs(value)) if name == "sand2023" else tolerance
        assert computed == pytest.approx(value, abs=tolerance), (horizon, computed, value)


COMPONENT_GROUPS = {
    # published component -> parameter-set components it is made of
    "O3": ("O3, short-term", "O3, long-term (CH4-induced)", "O3, CH4-induced", "O3, tropospheric"),
    "CH4": ("CH4",),
    "sulfate": ("sulfate",),
    "stratospheric H2O": ("stratospheric H2O", "stratospheric H2O, CH4-induced"),
}

PUBLISHED_COMPONENTS = [
    (species, name, group)
    for species, sets in DATA["sets"].items()
    for name, p in sets.items()
    if p["published"] and "components" in p["published"]
    for group in p["published"]["components"]
    if group != "sulfate"  # see test_fry2012_sulfate
]


@pytest.mark.parametrize(
    "species, name, group", PUBLISHED_COMPONENTS, ids=[f"{s}-{n}-{g}" for s, n, g in PUBLISHED_COMPONENTS]
)
def test_parameter_set_reproduces_published_components(species, name, group):
    published = DATA["sets"][species][name]["published"]
    specs = [spec for comp, spec in gases.precursor_component_specs(species, name) if comp in COMPONENT_GROUPS[group]]
    assert specs, group
    for horizon, value in zip(published["horizons"], published["components"][group]):
        computed = sum(_gwp(spec, published, horizon) for spec in specs)
        tolerance = max(0.03 * abs(value), _rounding(value))
        assert computed == pytest.approx(value, abs=tolerance), (horizon, computed, value)


@pytest.mark.parametrize("species", ["CO", "VOC"])
def test_fry2012_sulfate(species):
    """
    Fry et al. 2012 Table S2 sulfate includes a long-term response that is not given
    separately; the short-term sulfate forcing of Table S1 differs from it by less than 2 % of
    the total GWP (CO: -0.33 vs -0.4 at H = 20).
    """
    published = DATA["sets"][species]["fry2012_4regions"]["published"]
    (spec,) = [spec for comp, spec in gases.precursor_component_specs(species) if comp == "sulfate"]
    for horizon, value, total in zip(published["horizons"], published["components"]["sulfate"], published["total"]):
        assert abs(_gwp(spec, published, horizon) - value) < 0.02 * abs(total)


def test_fry2012_nox():
    """
    O3 and CH4 match Fry et al. 2012 Table S2; the sulfate forcing (from Table S1) gives
    -4.0 / -1.1 instead of the -3.1 / -0.8 of Table S2, which includes a long-term sulfate
    response not given separately. The total stays well within the published uncertainty.
    """
    published = DATA["sets"]["NOx"]["fry2012_4regions"]["published"]
    spec = gases.precursor("NOx", "fry2012_4regions")
    for horizon, value, uncertainty in zip(published["horizons"], published["total"], (26.1, 10.2)):
        assert abs(_gwp(spec, published, horizon) - value) < 0.1 * uncertainty


def test_mean_of_three_co_sets():
    sets = ("derwent2001", "berntsen2005_uio", "berntsen2005_lmdz")
    mean = np.mean([float(ar6_core.agwp(gases.precursor("CO", s), 1e4)) for s in sets])
    assert float(ar6_core.agwp(gases.precursor("CO", "fuglestvedt2010_mean3"), 1e4)) == pytest.approx(mean, rel=1e-3)


@pytest.mark.parametrize("precursor_lifetime", [0.0, 1e-6, 2.0])
def test_precursor_component_integrates_to_steady_state_forcing(precursor_lifetime):
    spec = ar6_core.GasSpec("x", (ar6_core.precursor_component(3e-13, 11.8, precursor_lifetime),))
    assert float(ar6_core.agwp(spec, 1e4)) == pytest.approx(3e-13, rel=1e-9)
    assert float(ar6_core.radiative_forcing(spec, 0.0)) == pytest.approx(0.0 if precursor_lifetime else 3e-13 / 11.8, abs=1e-25)


def test_precursor_lifetime_limit():
    """A tiny precursor lifetime gives the same metrics as none."""
    a = ar6_core.GasSpec("a", (ar6_core.precursor_component(3e-13, 11.8, 0.0),))
    b = ar6_core.GasSpec("b", (ar6_core.precursor_component(3e-13, 11.8, 1e-6),))
    for horizon in (1.0, 20.0, 100.0):
        assert float(ar6_core.agwp(b, horizon)) == pytest.approx(float(ar6_core.agwp(a, horizon)), rel=1e-5)
        assert float(ar6_core.agtp(b, horizon)) == pytest.approx(float(ar6_core.agtp(a, horizon)), rel=1e-5)


def test_precursor_component_rejects_bad_time_constants():
    with pytest.raises(ValueError):
        ar6_core.precursor_component(1e-13, 0.0, 0.0)
    with pytest.raises(ValueError):
        ar6_core.precursor_component(1e-13, 2.0, 2.0)


def test_signs():
    assert float(ar6_core.agwp(gases.nox("NO"), 100.0)) < 0
    assert float(ar6_core.agwp(gases.nox("NO", nitrate=True), 100.0)) < float(ar6_core.agwp(gases.nox("NO"), 100.0))
    for spec in (gases.co("non-fossil"), gases.voc(), gases.h2()):
        assert float(ar6_core.agwp(spec, 100.0)) > 0


def test_fossil_co_adds_oxidation_co2():
    oxidation = float(ar6_core.agwp(gases.co_oxidation(), 100.0))
    fossil = float(ar6_core.agwp(gases.co("fossil"), 100.0))
    non_fossil = float(ar6_core.agwp(gases.co("non-fossil"), 100.0))
    assert fossil - non_fossil == pytest.approx(oxidation, rel=1e-12)
    assert float(ar6_core.agwp(gases.co("non-fossil", biogenic_oxidation_co2="explicit"), 100.0)) == pytest.approx(fossil)


def test_mass_basis():
    per_kg_n = float(ar6_core.agwp(gases.precursor("NOx"), 100.0))
    assert float(ar6_core.agwp(gases.nox("NO"), 100.0)) == pytest.approx(per_kg_n * 14.007 / 30.006)
    assert float(ar6_core.agwp(gases.nox("NO2"), 100.0)) == pytest.approx(per_kg_n * 14.007 / 46.006)
    per_kg_c = float(ar6_core.agwp(gases.precursor("VOC"), 100.0))
    assert float(ar6_core.agwp(gases.voc(), 100.0)) == pytest.approx(per_kg_c * 0.845)
    assert float(ar6_core.agwp(gases.voc("collins2002"), 100.0)) == pytest.approx(
        float(ar6_core.agwp(gases.precursor("VOC", "collins2002"), 100.0))
    )


def test_option_validation():
    with pytest.raises(ValueError):
        CharacterizationOptions.from_value({"include_precursors": ["CO", "NOx"]})
    with pytest.raises(ValueError):
        CharacterizationOptions.from_value({"precursor_parameter_sets": {"CO": "nonexistent"}})
    with pytest.raises(ValueError):
        CharacterizationOptions.from_value({"nox_nitrate": True})
    options = CharacterizationOptions.from_value({"precursor_parameter_sets": {"H2": "sand2023"}, "include_precursors": {"H2", "CO"}})
    assert options.precursor_parameter_set("H2") == "sand2023"
    assert options.include_precursors == ("CO", "H2")
    hash(options)


# --- flows of the ecoinvent IPCC 2021 "incl. SLCFs" methods ---------------------------------

AR6 = {"parameter_set": "ar6"}
METHOD = ("test", "IPCC 2021", "incl. SLCFs")
BIOFLOWS = {
    # name, CAS, static GWP100 of ecoinvent 3.12 IPCC 2021 incl. SLCFs
    "co_fossil": ("Carbon monoxide, fossil", "000630-08-0", 4.0624),
    "co_soil": ("Carbon monoxide, from soil or biomass stock", "000630-08-0", 4.0624),
    "co_non_fossil": ("Carbon monoxide, non-fossil", "000630-08-0", 2.491),
    "h2": ("Hydrogen", "001333-74-0", 11.6),
    "no": ("Nitric oxide", "010102-43-9", -10.8),
    "nox": ("Nitrogen oxides", "011104-93-1", 0.0),
    "voc": ("VOC, volatile organic compounds", None, 4.5),
    "nmvoc": ("NMVOC, non-methane volatile organic compounds", None, 0.0),
    "hydrogen_sulfide": ("Hydrogen sulfide", "007783-06-4", 0.0),
}


def setup_project():
    clear_characterization_function_cache()
    biosphere = bd.Database(bd.config.biosphere)
    biosphere.write(
        {
            (bd.config.biosphere, code): {"code": code, "name": name, "categories": ("air",), "type": "emission"}
            | ({"CAS number": cas} if cas else {})
            for code, (name, cas, _) in BIOFLOWS.items()
        }
    )
    bd.Method(METHOD).write([((bd.config.biosphere, code), cf or 1.0) for code, (*_, cf) in BIOFLOWS.items()])
    return {code: biosphere.get(code).id for code in BIOFLOWS}


@bw2test
def test_default_precursor_functions():
    ids = setup_project()
    table = compare_with_static(METHOD, 100, AR6).set_index("flow")
    gwp = table["dynamic_gwp"]
    assert gwp[ids["co_non_fossil"]] == pytest.approx(1.63, abs=0.01)  # Fry et al. 2012 with AR6 CO2
    assert gwp[ids["co_fossil"]] == pytest.approx(1.63 + 1.57, abs=0.02)
    assert gwp[ids["co_soil"]] == gwp[ids["co_fossil"]]
    assert gwp[ids["h2"]] == pytest.approx(11.37, abs=0.02)  # Warwick et al. 2023
    assert gwp[ids["voc"]] == pytest.approx(4.5, rel=0.01)  # Fry et al. 2012, 0.845 kg C / kg VOC
    assert gwp[ids["nmvoc"]] == gwp[ids["voc"]]
    for code in ("no", "nox", "hydrogen_sulfide"):
        assert np.isnan(gwp[ids[code]]), code  # NOx is opt-in


@bw2test
def test_nox_opt_in_and_parameter_sets():
    ids = setup_project()
    options = dict(AR6, include_nox=True, precursor_parameter_sets={"H2": "sand2023", "CO": "derwent2001"})
    gwp = compare_with_static(METHOD, 100, options).set_index("flow")["dynamic_gwp"]
    assert gwp[ids["no"]] == pytest.approx(-8.59 * 14.007 / 30.006 * 8.56 / 8.947, rel=0.01)
    assert gwp[ids["nox"]] == pytest.approx(gwp[ids["no"]] * 30.006 / 46.006)
    assert gwp[ids["h2"]] == pytest.approx(11.73, abs=0.02)
    assert gwp[ids["co_non_fossil"]] == pytest.approx(1.90, abs=0.01)


@bw2test
def test_excluded_precursors():
    ids = setup_project()
    functions = create_characterization_functions_from_method(
        METHOD, characterization_options=dict(AR6, include_precursors=[])
    )
    assert ids["h2"] not in functions and ids["voc"] not in functions
    # CO falls back to its oxidation CO2 only, for every origin
    assert functions[ids["co_non_fossil"]] is functions[ids["co_fossil"]]
    assert functions[ids["co_fossil"]].agwp(100) == pytest.approx(float(ar6_core.agwp(gases.co_oxidation(), 100.0)))


@bw2test
def test_precursor_ccf_increases_gwp():
    ids = setup_project()
    plain = compare_with_static(METHOD, 100, AR6).set_index("flow")["dynamic_gwp"]
    with_ccf = compare_with_static(METHOD, 100, dict(AR6, precursor_ccf=True)).set_index("flow")["dynamic_gwp"]
    assert 1.03 < with_ccf[ids["h2"]] / plain[ids["h2"]] < 1.08


@bw2test
def test_characterize_with_precursors():
    ids = setup_project()
    df = pd.DataFrame(
        {
            "date": pd.Series(["2024-01-01"] * 3, dtype="datetime64[s]"),
            "amount": pd.Series([10.0, 2.0, 5.0], dtype="float64"),
            "flow": pd.Series([ids["co_fossil"], ids["h2"], ids["voc"]], dtype="int"),
            "activity": pd.Series([1, 1, 1], dtype="int"),
        }
    )
    rf = characterize(df, base_lcia_method=METHOD, time_horizon=100, characterization_options=AR6)
    assert set(rf.flow) == {ids["co_fossil"], ids["h2"], ids["voc"]}
    gwp = characterize(df, metric="GWP", base_lcia_method=METHOD, time_horizon=100, characterization_options=AR6)
    assert gwp.amount.sum() == pytest.approx(10 * 3.20 + 2 * 11.37 + 5 * 4.50, rel=0.01)
