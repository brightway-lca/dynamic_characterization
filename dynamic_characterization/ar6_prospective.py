"""
The ``"ar6"`` parameter set with a prospective scenario (metrics pGWP, pGTP and
prospective_radiative_forcing).

The pulse responses of ``ipcc_ar6`` are evaluated with the scenario of
``prospective.set_scenario`` (Watanabe et al. 2026) for the year of each emission:

- every CO2 term (emitted CO2, CO2 uptake, oxidation CO2 of fossil CH4 and CO, the CO2 of the
  climate-carbon feedback) follows the scenario's CO2 radiative efficiency and impulse response;
- CH4 and N2O use the scenario's direct radiative efficiencies, with the AR6 tropospheric
  adjustments and indirect terms (ozone, stratospheric water vapour, the CH4 effect of N2O);
- all other gases (halocarbons, precursors) keep their AR6 / literature forcing.

The temperature response is the AR6 one, for gases and CO2 alike. The scenario CO2 impulse
response covers 100 years, so do these metrics.
"""

from dataclasses import dataclass
from functools import lru_cache
from typing import Optional, Tuple

import numpy as np

from dynamic_characterization.classes import CharacterizedRow
from dynamic_characterization.ipcc_ar6 import ar6_core, gases
from dynamic_characterization.prospective.config import get_scenario
from dynamic_characterization.prospective.data_loader import (
    load_irf_co2,
    load_re_ch4,
    load_re_co2,
    load_re_n2o,
)
from dynamic_characterization.prospective.radiative_forcing import _get_year_index, _rcp_to_irf_key

_RE_LOADERS = {"CO2": load_re_co2, "CH4": load_re_ch4, "N2O": load_re_n2o}


def _emission_year(date) -> int:
    return int(str(np.datetime64(date, "s"))[:4])


@lru_cache(maxsize=None)
def scenario_radiative_efficiency(gas: str, iam: str, ssp: str, rcp: str, year: int) -> float:
    """Direct radiative efficiency [W m-2 ppb-1] of CO2, CH4 or N2O in a scenario and year."""
    data = _RE_LOADERS[gas](iam)
    series = data[ssp][rcp]
    return float(series[min(_get_year_index(year, data["_years"]), len(series) - 1)])


@lru_cache(maxsize=None)
def scenario_co2_response(iam: str, ssp: str, rcp: str, year: int) -> ar6_core.TabulatedCO2Response:
    """RF of 1 kg CO2 emitted in ``year`` of the scenario, from its CO2 RE and impulse response."""
    per_kg = scenario_radiative_efficiency("CO2", iam, ssp, rcp, year) / ar6_core.kg_per_ppb(
        ar6_core.load_parameters()["constants"]["M_co2"]
    )
    irf = load_irf_co2()[_rcp_to_irf_key(rcp)]
    return ar6_core.TabulatedCO2Response(
        name=f"{iam}-{ssp}-{rcp}, {year}", rf=tuple(per_kg * float(value) for value in irf)
    )


@dataclass(frozen=True)
class Recipe:
    """
    How to build the pulse response of a flow for a scenario and year: ``kind`` "CH4" and "N2O"
    take the scenario's radiative efficiencies; ``kind`` "fixed" uses ``spec`` as it is (its CO2
    terms still follow the scenario).
    """

    kind: str
    options: Tuple[Tuple[str, object], ...] = ()
    spec: Optional[ar6_core.GasSpec] = None

    def spec_for(self, iam: str, ssp: str, rcp: str, year: int) -> ar6_core.GasSpec:
        options = dict(self.options)
        if self.kind == "CH4":
            return gases.ch4(
                direct_radiative_efficiency=scenario_radiative_efficiency("CH4", iam, ssp, rcp, year),
                **options,
            )
        if self.kind == "N2O":
            return gases.n2o(
                direct_radiative_efficiency=scenario_radiative_efficiency("N2O", iam, ssp, rcp, year),
                ch4_direct_radiative_efficiency=scenario_radiative_efficiency("CH4", iam, ssp, rcp, year),
                **options,
            )
        return self.spec


@lru_cache(maxsize=None)
def _spec_and_response(recipe: Recipe, iam: str, ssp: str, rcp: str, year: int):
    return recipe.spec_for(iam, ssp, rcp, year), scenario_co2_response(iam, ssp, rcp, year)


def _check_period(period: float, response: ar6_core.TabulatedCO2Response) -> None:
    if period > response.max_years:
        raise ValueError(
            f"parameter_set='ar6' with a prospective scenario covers time horizons up to "
            f"{response.max_years} years (the length of the scenario CO2 impulse response), "
            f"not {period}."
        )


@lru_cache(maxsize=None)
def make_prospective_function(recipe: Recipe):
    """
    Characterization function for a recipe, with the legacy signature ``(series, period=100,
    cumulative=False)``. The scenario is read from ``prospective.get_scenario()`` and the year
    from the emission date. As the other ``"ar6"`` functions, the value at ``date + k`` years is
    the forcing of the year starting then. ``.agwp_for(date, H)`` and ``.agtp_for(date, H)`` give
    the per-kg metrics.
    """

    def resolve(date):
        scenario = get_scenario()
        return _spec_and_response(recipe, scenario["iam"], scenario["ssp"], scenario["rcp"], _emission_year(date))

    def characterize_ar6_prospective(series, period: int = 100, cumulative: bool = False) -> CharacterizedRow:
        spec, response = resolve(series.date)
        _check_period(period, response)
        date_beginning: np.datetime64 = series.date.to_numpy()
        dates = date_beginning + np.arange(start=0, stop=period, dtype="timedelta64[Y]").astype("timedelta64[s]")
        forcing = np.array(series.amount * ar6_core.annual_cumulative(spec, period, response), dtype="float64")
        if not cumulative:
            forcing = np.diff(forcing, prepend=0)
        return CharacterizedRow(
            date=np.array(dates, dtype="datetime64[s]"), amount=forcing, flow=series.flow, activity=series.activity
        )

    def agwp_for(date, time_horizon) -> float:
        spec, response = resolve(date)
        _check_period(time_horizon, response)
        return float(ar6_core.agwp(spec, float(time_horizon), co2_response=response))

    def agtp_for(date, time_horizon) -> float:
        spec, response = resolve(date)
        _check_period(time_horizon, response)
        return float(ar6_core.agtp(spec, float(time_horizon), co2_response=response))

    name = recipe.spec.name if recipe.spec is not None else recipe.kind
    characterize_ar6_prospective.__name__ = f"characterize_ar6_prospective[{name}]"
    characterize_ar6_prospective.__qualname__ = characterize_ar6_prospective.__name__
    characterize_ar6_prospective.recipe = recipe
    characterize_ar6_prospective.agwp_for = agwp_for
    characterize_ar6_prospective.agtp_for = agtp_for
    return characterize_ar6_prospective


def fixed(spec: ar6_core.GasSpec):
    """Prospective function for a pulse response that does not depend on the scenario (except CO2)."""
    return make_prospective_function(Recipe(kind="fixed", spec=spec))


characterize_co2_ar6_prospective = fixed(gases.co2())
characterize_co2_uptake_ar6_prospective = fixed(gases.co2_uptake())
