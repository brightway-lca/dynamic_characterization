"""
Dynamic characterization functions based on the IPCC AR6 emission metrics (WG1 Ch. 7, Tables 7.15
and 7.SM.7), built from the pulse responses of :mod:`.ar6_core` and :mod:`.gases`.

All functions have the signature ``(series, period=100, cumulative=False) -> CharacterizedRow``.
The value stamped at ``date + k`` years is the forcing integrated over the year starting then
[W yr m-2 per unit emission], so the annual values add up to the AGWP of the period.

The public names of earlier versions are kept:

- ``characterize_co2`` / ``characterize_co2_uptake``: CO2 (AR6 Table 7.SM.7 CO2 response);
- ``characterize_ch4``: fossil CH4 incl. ozone, stratospheric water vapour, oxidation CO2 and the
  climate-carbon feedback (GWP100 29.8);
- ``characterize_n2o``: N2O incl. its ozone and CH4 effects and the climate-carbon feedback (273);
- ``characterize_co``: fossil CO incl. its ozone, CH4 and sulfate effects (Fry et al. 2012) and
  oxidation CO2 (GWP100 3.2);
- ``create_generic_characterization_function``: a function from a user-supplied decay series.

Before version 2.0 these were simplified functions without climate-carbon feedback, without the
oxidation CO2 of fossil CH4, with CO as oxidation only, and with forcing summed to H - 1.
"""

from functools import lru_cache

import numpy as np

from dynamic_characterization.classes import CharacterizedRow

from . import ar6_core, gases


def create_generic_characterization_function(decay_series) -> CharacterizedRow:
    """
    Creates a characterization function for a GHG based on a decay series, by calling the nested method `characterize_generic()`.

    Parameters
    ----------
    decay_series : np.ndarray
        Cumulative forcing per unit emission [W yr m-2] at the years 0, 1, ... after the
        emission (the forcing of year k is ``decay_series[k] - decay_series[k - 1]``).

    Returns
    -------
    A function called `characterize_generic`, which in turn returns a TimeSeries dataframe that contains the forcing of the emission of the row over the given period based on the decay series of that biosphere flow.

    """

    def characterize_generic(
        series,
        period: int = 100,
        cumulative=False,
    ) -> CharacterizedRow:
        """
        Characterization function built from a user-supplied cumulative decay series.

        Parameters
        ----------
        series : array-like
            A single row of the dynamic inventory dataframe.
        period : int, optional
            Time period for calculation (number of years), by default 100
        cumulative : bool,
            cumulative impact

        Returns
        -------
        A TimeSeries dataframe that contains the forcing of the point emission from the row for each year in the given period.
          date: datetime64[s]
          amount: float (forcing at this timestep)
          flow: str
          activity: str

        See also
        --------
        Joos2013: Relevant scientific publication on CRF: https://doi.org/10.5194/acp-13-2793-2013
        Schivley2015: Relevant scientific publication on the numerical calculation of CRF: https://doi.org/10.1021/acs.est.5b01118
        Forster2023: Updated numerical values from IPCC AR6 Chapter 7 (Table 7.15): https://doi.org/10.1017/9781009157896.009

        """

        date_beginning: np.datetime64 = series.date.to_numpy()

        dates_characterized: np.ndarray = date_beginning + np.arange(
            start=0, stop=period, dtype="timedelta64[Y]"
        ).astype("timedelta64[s]")

        decay_multipliers = decay_series[:period]

        forcing = np.array(series.amount * decay_multipliers, dtype="float64")

        if not cumulative:
            forcing = np.diff(forcing, prepend=0)

        return CharacterizedRow(
            date=np.array(dates_characterized, dtype="datetime64[s]"),
            amount=forcing,
            flow=series.flow,
            activity=series.activity,
        )

    return characterize_generic


# --- characterization functions of the AR6 pulse responses -----------------------------------


@lru_cache(maxsize=None)
def make_characterization_function(spec):
    """
    Characterization function for a pulse response from ``ar6_core`` (a ``GasSpec``).

    The returned function has the same signature and output as the legacy functions in this
    module, ``(series, period=100, cumulative=False) -> CharacterizedRow``, with one
    difference in the time convention: the value stamped at ``date + k`` years is the forcing
    integrated over the year *starting* then, so the annual values add up to ``AGWP(period)``
    (the legacy functions add up to ``AGWP(period - 1)``).

    The function also exposes ``.spec``, ``.agwp(time_horizon)`` and ``.agtp(time_horizon)``
    (per kg). It is memoized per spec, so the same spec always gives the same function object.
    """
    def characterize_ar6(series, period: int = 100, cumulative: bool = False) -> CharacterizedRow:
        date_beginning: np.datetime64 = series.date.to_numpy()
        dates_characterized: np.ndarray = date_beginning + np.arange(
            start=0, stop=period, dtype="timedelta64[Y]"
        ).astype("timedelta64[s]")

        forcing = np.array(series.amount * ar6_core.annual_cumulative(spec, period), dtype="float64")
        if not cumulative:
            forcing = np.diff(forcing, prepend=0)

        return CharacterizedRow(
            date=np.array(dates_characterized, dtype="datetime64[s]"),
            amount=forcing,
            flow=series.flow,
            activity=series.activity,
        )

    characterize_ar6.__name__ = f"characterize_ar6[{spec.name}]"
    characterize_ar6.__qualname__ = characterize_ar6.__name__
    characterize_ar6.spec = spec
    characterize_ar6.agwp = lambda time_horizon: float(ar6_core.agwp(spec, float(time_horizon)))
    characterize_ar6.agtp = lambda time_horizon: float(ar6_core.agtp(spec, float(time_horizon)))
    return characterize_ar6


# CO2 with the AR6 impulse response: the GWP reference, and CO2 taken up from the atmosphere
characterize_co2_ar6 = make_characterization_function(gases.co2())
characterize_co2_uptake_ar6 = make_characterization_function(gases.co2_uptake())


# --- the public names of earlier versions ----------------------------------------------------

characterize_co2 = characterize_co2_ar6
characterize_co2_uptake = characterize_co2_uptake_ar6
characterize_ch4 = make_characterization_function(gases.ch4("fossil", ccf=ar6_core.ccf_parameters()))
characterize_n2o = make_characterization_function(gases.n2o(ccf=ar6_core.ccf_parameters()))
characterize_co = make_characterization_function(gases.co("fossil"))
