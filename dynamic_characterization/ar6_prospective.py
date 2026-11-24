"""
The ``"ar6"`` parameter set with a prospective scenario (metrics pGWP, pGTP and
prospective_radiative_forcing), for the scenario of ``prospective.set_scenario`` and the year of
each emission. Two methods (option ``prospective_method``):

``"watanabe2026"`` (default)
    CO2, CH4 and N2O as in Barbosa Watanabe & Cherubini (2026), with two errors of the published
    code corrected (see ``prospective.watanabe2026``): scenario RE and impulse responses, AR6 ozone
    and stratospheric water vapour terms, the CH4 effect of N2O, a climate-carbon feedback, and an
    annual two-layer temperature model.
    All other gases (halocarbons, CO, VOC, H2, NOx) keep their AR6 / literature forcing, with their
    CO2 terms following the same scenario CO2 response and their AGTP from the same temperature
    scheme, so that the CO2 reference is consistent.

``"ar6"``
    The AR6 recipe with scenario inputs: CO2 terms follow the scenario CO2 response, CH4 and N2O
    use the scenario direct RE with the AR6 tropospheric adjustments and indirect terms, and the
    temperature response is the AR6 one.

The scenario CO2 impulse response covers 100 years, so do these metrics: longer horizons are
truncated to 100 years, with a warning. ``time_varying_re=True`` lets the radiative
efficiencies of CO2, CH4 and N2O follow the scenario over the response (``"watanabe2026"`` only;
the ``"ar6"`` method uses those of the emission year and warns).
"""

from dataclasses import dataclass
from functools import lru_cache
from typing import Optional, Tuple

import warnings

import numpy as np

from dynamic_characterization.classes import CharacterizedRow
from dynamic_characterization.ipcc_ar6 import ar6_core, gases
from dynamic_characterization.prospective import watanabe2026
from dynamic_characterization.prospective.config import get_scenario
from dynamic_characterization.prospective.data_loader import (
    load_irf_co2,
    load_re_ch4,
    load_re_co2,
    load_re_n2o,
)
from dynamic_characterization.prospective._common import CONST_CO2, _get_year_index

METHODS = ("watanabe2026", "ar6")
MAX_YEARS = 100
_RE_LOADERS = {"CO2": load_re_co2, "CH4": load_re_ch4, "N2O": load_re_n2o}


def _emission_year(date) -> int:
    return int(str(np.datetime64(date, "s"))[:4])


def _current_scenario() -> Tuple[str, str, str]:
    scenario = get_scenario()
    return scenario["iam"], scenario["ssp"], scenario["rcp"]


@lru_cache(maxsize=None)
def scenario_radiative_efficiency(gas: str, iam: str, ssp: str, rcp: str, year: int) -> float:
    """Direct radiative efficiency [W m-2 ppb-1] of CO2, CH4 or N2O in a scenario and year."""
    data = _RE_LOADERS[gas](iam)
    series = data[ssp][rcp]
    return float(series[min(_get_year_index(year, data["_years"]), len(series) - 1)])


@lru_cache(maxsize=None)
def scenario_co2_response(iam: str, ssp: str, rcp: str, year: int, method: str = "ar6") -> ar6_core.TabulatedCO2Response:
    """
    RF of 1 kg CO2 emitted in ``year`` of the scenario: CO2 RE of the emission year x impulse
    response of the scenario's RCP. With ``method="watanabe2026"`` the ppb-to-kg conversion is
    that of the paper's code, so that its AGWP equals the CO2 AGWP of ``prospective.watanabe2026``.
    """
    re_co2 = scenario_radiative_efficiency("CO2", iam, ssp, rcp, year)
    if method == "watanabe2026":
        per_kg = re_co2 * CONST_CO2
    else:
        per_kg = re_co2 / ar6_core.kg_per_ppb(ar6_core.load_parameters()["constants"]["M_co2"])
    irf = load_irf_co2()[watanabe2026.co2_irf_key(rcp)]
    return ar6_core.TabulatedCO2Response(name=f"{iam}-{ssp}-{rcp}, {year}", rf=tuple(per_kg * float(v) for v in irf))


@dataclass(frozen=True)
class Recipe:
    """
    How to build the response of a flow for a scenario and year.

    ``kind``: "CH4" / "N2O" (scenario radiative efficiencies), "fixed" (``spec`` as it is, with
    its CO2 terms following the scenario), or "paper" (``gas`` "CO2", "CH4" or "N2O" from
    ``prospective.watanabe2026``, multiplied by ``sign``).
    """

    kind: str
    method: str = "ar6"
    options: Tuple[Tuple[str, object], ...] = ()
    spec: Optional[ar6_core.GasSpec] = None
    gas: Optional[str] = None
    sign: float = 1.0

    def spec_for(self, iam: str, ssp: str, rcp: str, year: int) -> ar6_core.GasSpec:
        options = dict(self.options)
        if self.kind == "CH4":
            return gases.ch4(direct_radiative_efficiency=scenario_radiative_efficiency("CH4", iam, ssp, rcp, year), **options)
        if self.kind == "N2O":
            return gases.n2o(
                direct_radiative_efficiency=scenario_radiative_efficiency("N2O", iam, ssp, rcp, year),
                ch4_direct_radiative_efficiency=scenario_radiative_efficiency("CH4", iam, ssp, rcp, year),
                **options,
            )
        return self.spec

    @property
    def label(self) -> str:
        if self.kind == "paper":
            return f"{self.gas} (Watanabe & Cherubini 2026)" + (" uptake" if self.sign < 0 else "")
        return self.spec.name if self.spec is not None else self.kind


_WARNED = set()


def _warn_once(key: str, message: str) -> None:
    if key not in _WARNED:
        _WARNED.add(key)
        warnings.warn(message, UserWarning, stacklevel=3)


def _truncated(period) -> int:
    """The horizon actually characterized: at most 100 years (with a one-time warning)."""
    if period > MAX_YEARS:
        _warn_once(
            "horizon",
            f"Prospective characterization covers {MAX_YEARS} years (the length of the scenario CO2 "
            f"impulse response); time horizons beyond that are truncated to {MAX_YEARS} years.",
        )
        return MAX_YEARS
    return int(period)


@lru_cache(maxsize=None)
def _annual(recipe: Recipe, iam: str, ssp: str, rcp: str, year: int, time_varying_re: bool = False):
    """(cumulative AGWP at the end of years 1..100, AGTP at years 1..100), per kg."""
    if recipe.kind == "paper":
        agwp = recipe.sign * watanabe2026.scenario_series(recipe.gas, "agwp", iam, ssp, rcp, year, time_varying_re)[:MAX_YEARS]
        agtp = recipe.sign * watanabe2026.scenario_series(recipe.gas, "agtp", iam, ssp, rcp, year, time_varying_re)[:MAX_YEARS]
        return np.asarray(agwp), np.asarray(agtp)
    if time_varying_re:
        if recipe.method == "ar6":
            _warn_once(
                "time_varying_re",
                "time_varying_re is not supported by prospective_method='ar6'; the radiative "
                "efficiencies of the emission year are used.",
            )
        time_varying_re = False
    spec = recipe.spec_for(iam, ssp, rcp, year)
    response = scenario_co2_response(iam, ssp, rcp, year, recipe.method)
    agwp = ar6_core.annual_cumulative(spec, MAX_YEARS, response)
    if recipe.method == "watanabe2026":
        # AGTP with the paper's temperature scheme (Code #7), applied to the annual RF
        rf = ar6_core.radiative_forcing(spec, np.arange(MAX_YEARS, dtype="float64"), co2_response=response)
        agtp = watanabe2026.agtp_from_rf(rf)
    else:
        agtp = ar6_core.agtp(spec, np.arange(1, MAX_YEARS + 1, dtype="float64"), co2_response=response)
    return np.asarray(agwp), np.asarray(agtp)


@lru_cache(maxsize=None)
def make_prospective_function(recipe: Recipe):
    """
    Characterization function for a recipe, with the legacy signature ``(series, period=100,
    cumulative=False)``. The scenario is read from ``prospective.get_scenario()`` and the year
    from the emission date. As for the other ``"ar6"`` functions, the value at ``date + k`` years
    is the forcing of the year starting then; the annual values add up to the AGWP of the
    horizon. ``.agwp_for(date, H)`` and ``.agtp_for(date, H)`` give the per-kg metrics.
    """

    def annual(date, time_varying_re=False):
        return _annual(recipe, *_current_scenario(), _emission_year(date), bool(time_varying_re))

    def characterize_ar6_prospective(
        series, period: int = 100, cumulative: bool = False, time_varying_re: bool = False
    ) -> CharacterizedRow:
        period = _truncated(period)
        agwp, _ = annual(series.date, time_varying_re)
        date_beginning: np.datetime64 = series.date.to_numpy()
        dates = date_beginning + np.arange(start=0, stop=period, dtype="timedelta64[Y]").astype("timedelta64[s]")
        forcing = np.array(series.amount * agwp[:period], dtype="float64")
        if not cumulative:
            forcing = np.diff(forcing, prepend=0)
        return CharacterizedRow(
            date=np.array(dates, dtype="datetime64[s]"), amount=forcing, flow=series.flow, activity=series.activity
        )

    def agwp_for(date, time_horizon, time_varying_re: bool = False) -> float:
        horizon = _truncated(time_horizon)
        return float(annual(date, time_varying_re)[0][horizon - 1]) if horizon >= 1 else 0.0

    def agtp_for(date, time_horizon, time_varying_re: bool = False) -> float:
        horizon = _truncated(time_horizon)
        return float(annual(date, time_varying_re)[1][horizon - 1]) if horizon >= 1 else 0.0

    characterize_ar6_prospective.__name__ = f"characterize_ar6_prospective[{recipe.label}, {recipe.method}]"
    characterize_ar6_prospective.__qualname__ = characterize_ar6_prospective.__name__
    characterize_ar6_prospective.recipe = recipe
    characterize_ar6_prospective.agwp_for = agwp_for
    characterize_ar6_prospective.agtp_for = agtp_for
    return characterize_ar6_prospective


def fixed(spec: ar6_core.GasSpec, method: str = "ar6"):
    """Prospective function for a pulse response that does not depend on the scenario (except CO2)."""
    return make_prospective_function(Recipe(kind="fixed", spec=spec, method=method))


def co2_function(method: str = "watanabe2026", uptake: bool = False):
    """The CO2 (or CO2 uptake) function of a method; the CO2 one is also the pGWP/pGTP reference."""
    if method == "watanabe2026":
        return make_prospective_function(
            Recipe(kind="paper", method=method, gas="CO2", sign=-1.0 if uptake else 1.0)
        )
    return fixed(gases.co2_uptake() if uptake else gases.co2(), method="ar6")


def paper_function(gas: str):
    """CH4 or N2O as in Watanabe & Cherubini (2026), corrected."""
    return make_prospective_function(Recipe(kind="paper", method="watanabe2026", gas=gas))
