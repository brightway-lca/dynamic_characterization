"""The "ar6" pulse responses reproduce the IPCC AR6 emission metrics (Tables 7.15 and 7.SM.7)."""

import csv
import importlib.util
import os

import numpy as np
import pandas as pd
import pytest

from dynamic_characterization.flows import normalize_cas
from dynamic_characterization.ipcc_ar6 import ar6_core, gases
from dynamic_characterization.ipcc_ar6.radiative_forcing import (
    characterize_co2_ar6,
    make_characterization_function,
)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AR6_TABLE = os.path.join(ROOT, "tests", "data", "ar6_metrics_supplement_cleaned.csv")


def _load_dev_script(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(ROOT, "dev", f"{name}.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


BUILD_AR6_DATA = _load_dev_script("build_ar6_data")

CCF = ar6_core.ccf_parameters("gasser2017_mean")
CO2 = gases.co2()
GWP_HORIZONS = (20.0, 100.0, 500.0)
GTP_HORIZONS = (50.0, 100.0)


def gwp(spec, horizon):
    return float(ar6_core.agwp(spec, horizon) / ar6_core.agwp(CO2, horizon))


def gtp(spec, horizon):
    return float(ar6_core.agtp(spec, horizon) / ar6_core.agtp(CO2, horizon))


def published_tolerance(value):
    """Table 7.SM.7 is rounded to 3 significant figures, small values to 3 decimals."""
    return max(0.01 * abs(value), 6e-4)


def _ar6_rows():
    with open(AR6_TABLE, newline="") as f:
        rows = list(csv.DictReader(f))
    for row in rows:
        row["CAS"] = normalize_cas(row["CAS"].strip('="'))
    return rows


AR6_ROWS = _ar6_rows()
SF6_CAS = "2551-62-4"
HALOCARBON_ROWS = [row for row in AR6_ROWS if row["CAS"] is not None and row["Lifetime (yr)"]]
# SF6 follows AR6 Table 7.15 (lifetime 1000 yr), not Table 7.SM.7 (3200 yr), see test_sf6_matches_table_7_15
HALOCARBON_ROWS_SM7 = [row for row in HALOCARBON_ROWS if row["CAS"] != SF6_CAS]


def test_parameter_file_is_up_to_date():
    path = os.path.join(ROOT, "dynamic_characterization", "ipcc_ar6", "data", "ar6_parameters.json")
    with open(path) as f:
        assert f.read() == BUILD_AR6_DATA.render(), "run python dev/build_ar6_data.py"


def test_co2_agwp_and_agtp():
    co2_row = AR6_ROWS[0]
    assert co2_row["Formula"] == "CO2"
    for horizon in GWP_HORIZONS:
        published = float(co2_row[f"AGWP{int(horizon)} (W m-2 yr kg-1)"])
        assert ar6_core.agwp(CO2, horizon) == pytest.approx(published, rel=0.005)
    for horizon in GTP_HORIZONS:
        published = float(co2_row[f"AGTP{int(horizon)} (K kg-1)"])
        assert ar6_core.agtp(CO2, horizon) == pytest.approx(published, rel=0.005)
    # AR6 Table 7.SM.7 / Section 7.6.1: AGWP20 and AGWP100 of CO2
    assert ar6_core.agwp(CO2, 20.0) == pytest.approx(2.434e-14, rel=1e-3)
    assert ar6_core.agwp(CO2, 100.0) == pytest.approx(8.947e-14, rel=1e-3)


@pytest.mark.parametrize(
    "formula, spec",
    [
        ("CH4", gases.ch4("non-fossil", "neutral", ccf=CCF)),  # the plain CH4 pulse of Table 7.SM.7
        ("N2O", gases.n2o(ccf=CCF)),
    ],
)
def test_ch4_and_n2o_match_table_7sm7(formula, spec):
    row = next(r for r in AR6_ROWS if r["Formula"] == formula)
    for horizon in GWP_HORIZONS:
        published = float(row[f"GWP{int(horizon)}"])
        assert gwp(spec, horizon) == pytest.approx(published, abs=published_tolerance(published))
    for horizon in GTP_HORIZONS:
        published = float(row[f"GTP{int(horizon)}"])
        assert gtp(spec, horizon) == pytest.approx(published, abs=published_tolerance(published))


@pytest.mark.parametrize(
    "origin, published",
    [
        # AR6 WG1 Table 7.15: GWP20, GWP100, GWP500
        ("fossil", (82.5, 29.8, 10.0)),
        ("non-fossil", (79.7, 27.0, 7.2)),
    ],
)
def test_ch4_fossil_and_non_fossil_match_table_7_15(origin, published):
    spec = gases.ch4(origin, ccf=CCF)
    for horizon, value in zip(GWP_HORIZONS, published):
        assert gwp(spec, horizon) == pytest.approx(value, abs=max(0.01 * value, 0.05))


def test_n2o_matches_table_7_15():
    for horizon, value in zip(GWP_HORIZONS, (273, 273, 130)):
        assert gwp(gases.n2o(ccf=CCF), horizon) == pytest.approx(value, rel=0.01)


def test_every_ar6_halocarbon_is_available():
    missing = [row["Name"] for row in HALOCARBON_ROWS if gases.halocarbon(row["CAS"]) is None]
    assert not missing


@pytest.mark.parametrize("row", HALOCARBON_ROWS_SM7, ids=[r["Name"] for r in HALOCARBON_ROWS_SM7])
def test_halocarbons_match_table_7sm7(row):
    spec = gases.halocarbon(row["CAS"], ccf=CCF)
    for horizon in GWP_HORIZONS:
        published = float(row[f"GWP{int(horizon)}"])
        assert gwp(spec, horizon) == pytest.approx(published, abs=published_tolerance(published))
    for horizon in GTP_HORIZONS:
        published = float(row[f"GTP{int(horizon)}"])
        assert gtp(spec, horizon) == pytest.approx(published, abs=published_tolerance(published))


def test_sf6_matches_table_7_15():
    """AR6 Table 7.15 and the static IPCC 2021 methods: SF6 GWP100 24 300 (lifetime 1000 yr)."""
    assert gwp(gases.halocarbon(SF6_CAS, ccf=CCF), 100.0) == pytest.approx(24300, rel=0.005)


def test_ccf_parameter_sets_are_ordered():
    values = [gwp(gases.ch4(ccf=ar6_core.ccf_parameters(name)), 100.0) for name in ("gasser2017_lower", "gasser2017_mean", "gasser2017_upper")]
    no_ccf = gwp(gases.ch4(ccf=None), 100.0)
    assert no_ccf < values[0] < values[1] < values[2]


def test_ccf_grows_relative_to_direct_forcing():
    """The CCF tail: small at first, larger than the direct CH4 forcing after several decades."""
    with_ccf = gases.ch4("non-fossil", "neutral", ccf=CCF)
    without = gases.ch4("non-fossil", "neutral", ccf=None)
    t = np.array([10.0, 50.0, 100.0])
    share = ar6_core.radiative_forcing(with_ccf, t) / ar6_core.radiative_forcing(without, t) - 1
    assert 0.0 < share[0] < 0.1
    assert share[1] > 0.2
    assert share[2] > 1.0


def test_oxidation_co2_is_the_convolution_of_source_and_co2_response():
    """co2_from_decay (analytic) equals a numerical convolution of the source with the CO2 RF."""
    dt = 0.001
    t = np.arange(0, 100 + dt / 2, dt)
    lifetime = 9.7
    source = np.exp(-t / lifetime) / lifetime
    rf_numeric = np.convolve(source, ar6_core.radiative_forcing(CO2, t))[: t.size] * dt
    rf_analytic = ar6_core.radiative_forcing(ar6_core.co2_from_decay(1.0, lifetime), t)
    assert np.allclose(rf_numeric[1000::1000], rf_analytic[1000::1000], rtol=2e-3)


def test_agtp_is_the_convolution_of_rf_and_temperature_response():
    dt = 0.001
    t = np.arange(0, 100 + dt / 2, dt)
    spec = ar6_core.exponential_spec("test", 1e-13, 12.0)
    q, d = (np.asarray(v) for v in (ar6_core.load_parameters()["temperature_response"][k] for k in "qd"))
    response = sum(qj / dj * np.exp(-t / dj) for qj, dj in zip(q, d))
    numeric = np.convolve(ar6_core.radiative_forcing(spec, t), response)[: t.size] * dt
    assert numeric[-1] == pytest.approx(float(ar6_core.agtp(spec, 100.0)), rel=2e-3)


# --- characterization functions built from the specs ---------------------------------------


class Row:
    def __init__(self, amount=1.0):
        self.date = pd.Timestamp("2024-01-01")
        self.amount = amount
        self.flow = 1
        self.activity = 2


@pytest.mark.parametrize("period", [1, 20, 100, 500])
def test_annual_forcing_adds_up_to_agwp(period):
    func = make_characterization_function(gases.ch4(ccf=CCF))
    result = func(Row(amount=2.0), period)
    assert result.amount.shape == (period,)
    assert result.amount.sum() == pytest.approx(2.0 * func.agwp(period), rel=1e-12)
    cumulative = func(Row(amount=2.0), period, cumulative=True).amount
    assert np.allclose(np.cumsum(result.amount), cumulative)


def test_first_year_is_stamped_at_the_emission_date():
    result = characterize_co2_ar6(Row(), 3)
    assert result.date[0] == np.datetime64("2024-01-01T00:00:00")
    assert result.amount[0] == pytest.approx(float(ar6_core.agwp(CO2, 1.0)))


def test_functions_are_memoized_per_spec():
    assert make_characterization_function(gases.n2o(ccf=CCF)) is make_characterization_function(gases.n2o(ccf=CCF))
    assert make_characterization_function(gases.n2o(ccf=CCF)) is not make_characterization_function(gases.n2o())
